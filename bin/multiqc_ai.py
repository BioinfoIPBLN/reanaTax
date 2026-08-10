#!/usr/bin/env python3
"""
multiqc_ai.py -- add per-section AI summaries to an existing MultiQC report.

Ported from reanalyzerGSE (scripts/multiqc_ai.py). There, the script also ran
MultiQC itself; here MultiQC is its own Nextflow process, so this script only
does the part that MultiQC cannot do on its own: discover every section in
multiqc_report.html, send that section's exported data file to an
OpenAI-compatible chat endpoint, and inject a pre-baked summary inline under the
section heading. It works for EVERY MultiQC module, not just the ones this
pipeline happens to run.

MultiQC's own AI feature (one global report summary) is orthogonal and is
enabled with --multiqc_ai_builtin, which passes the flags to MultiQC directly.
This script still runs afterwards, and its final job is to scrub the endpoint
that the builtin feature bakes into the report.

The LLM is any OpenAI-compatible /v1/chat/completions endpoint (a local vLLM /
Ollama / llama.cpp server, or api.openai.com, etc.). No endpoint or model is
baked in -- supply them at runtime via CLI flags or env vars:
    LLM_ENDPOINT   OpenAI-compatible /v1/chat/completions URL   (or --llm-endpoint)
    LLM_MODEL      model name                                   (or --llm-model)
    LLM_API_KEY    Bearer token; "dummy" for servers that ignore it (or --llm-api-key)
"""
import argparse
import os
import re
import sys

# Import the sibling shared module regardless of how this script is invoked
# (PATH, absolute path, symlink): put its own directory on sys.path first.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import llm_common  # noqa: E402

# ------------------------------------------------------- LLM config defaults
# Precedence for every parameter:  CLI flag  >  env var  >  fallback below.
# NOTE: endpoint and model have NO baked-in fallback on purpose -- no site-
# specific host or model name is hardcoded anywhere. They must be supplied at
# runtime (flag or env var); main() enforces this.
DEF_ENDPOINT = os.environ.get("LLM_ENDPOINT")
DEF_MODEL = os.environ.get("LLM_MODEL")
DEF_API_KEY = os.environ.get("LLM_API_KEY", "dummy")
DEF_CONTEXT_WINDOW = int(os.environ.get("LLM_CONTEXT_WINDOW", "128000"))
DEF_TIMEOUT = int(os.environ.get("LLM_TIMEOUT", "300"))

DEF_SYS_PROMPT = "\n".join(
    [
        "You are an expert in bioinformatics, sequencing, and QC reports.",
        "The report you are annotating comes from a metagenomic pipeline: reads are trimmed,",
        "aligned against a host genome to remove host sequence, and whatever does not align is",
        "classified taxonomically. A HIGH host alignment rate is expected and is not a QC failure;",
        "what matters is whether enough non-host reads are left to classify.",
        "You are given the data for a single section of a MultiQC report (a plot, table, or heatmap).",
        "Produce 1-2 concise bullet points summarising the key observations and any QC issues.",
        "If nothing is concerning, say so in one bullet.",
        "Use markdown. Highlight severity with directives like :span[39.2%]{.text-red}, .text-orange, .text-yellow, .text-green.",
        "Highlight sample names with :sample[name]{.text-red} etc.",
        "Use 4 spaces to indent nested lists. Do not add headers.",
        "Base every statement only on the numbers provided; never invent samples or values.",
        "Always respond in English.",
    ]
)


# --------------------------------------------------------- LLM call helper
def ask_llm(cfg, section_id, data):
    """One section summary. Sequential + text-only via llm_common. A timeout is
    re-raised (the caller then marks this and the remaining sections as timed
    out); any other error is swallowed so one bad section never aborts the pass."""
    data = llm_common.cap_text(data, cfg.context_window)
    messages = [
        {"role": "system", "content": cfg.sys_prompt},
        {"role": "user", "content": f"Section: {section_id}\n\nData:\n{data}"},
    ]
    try:
        text, usage = llm_common.chat_completion(
            cfg.endpoint, cfg.model, cfg.api_key, messages, timeout=cfg.timeout
        )
    except llm_common.LLMTimeout:
        raise
    except Exception as e:
        llm_common.log(f"[MultiQC AI] error on {section_id}: {e}")
        return "_(no AI response)_", {}
    return (text or "_(no AI response)_"), usage


# ------------------------------------------ discover sections FROM the HTML
def discover_sections(html):
    """Return the set of MultiQC section anchor ids actually present in the
    report. Generic: works for every module, not just the ones we ship."""
    ids = set(re.findall(r'id="mqc-section-wrapper-([^"]+)"', html))
    ids |= set(re.findall(r'id="([^"]+)_ai_summary_response"', html))
    return ids


def data_file_for_section(data_dir, sec_id):
    """Best-effort map a MultiQC HTML section id to its exported data file.
    The file name often differs from the section id in case (featureCounts vs
    featurecounts), punctuation (-/_) and plural/singular (assignment vs
    assignments); some sections also carry a generic '<module>-section-N' id
    that encodes no file name at all. Match tolerantly."""
    if sec_id == "general_stats_table":
        p = os.path.join(data_dir, "multiqc_general_stats.txt")
        return p if os.path.exists(p) else None
    try:
        files = [f for f in os.listdir(data_dir) if f.lower().endswith(".txt")]
    except OSError:
        return None

    def first_match(patterns):
        for pat in patterns:
            rx = re.compile(pat + r"$", re.IGNORECASE)
            hits = sorted(f for f in files if rx.match(f))
            if hits:
                return os.path.join(data_dir, hits[0])
        return None

    # Case-insensitive, punctuation- (-/_) and plural-tolerant spellings of the id.
    variants = set()
    for base in (sec_id, sec_id.replace("-", "_"), sec_id.replace("_", "-")):
        variants.add(base)
        variants.add(base[:-1] if base.endswith("s") else base + "s")
    precise = []
    for v in sorted(variants):
        e = re.escape(v)
        precise += [rf"{e}[-_]plot.*\.txt", rf"{e}[-_]table\.txt", rf"{e}.*-heatmap\.txt", rf"{e}\.txt"]
    hit = first_match(precise)
    if hit:
        return hit

    # Check key distinguishing descriptors in sec_id (e.g. length, duplication, adapter)
    # to avoid mis-mapping when the module has multiple distinct plot files.
    key_descriptors = ["length", "duplication", "adapter", "overrepresented", "gc"]
    sec_low = sec_id.lower()
    required_terms = [d for d in key_descriptors if d in sec_low]
    if "n_content" in sec_low or "per_base_n" in sec_low:
        required_terms.append("n_content")

    # Fuzzy token fallback for multi-token / hyphenated ids
    skip = (
        "multiqc_citations",
        "multiqc_sources",
        "multiqc_software_versions",
        "multiqc_data_sources",
        "section_prompts",
        "llms-full",
        "multiqc_general_stats",
    )
    toks = [t for t in re.split(r"[-_]", sec_id.lower()) if t and t != "section" and not t.isdigit()]
    module = toks[0] if toks else ""
    if not toks or module == "multiqc":
        return None

    need = 2 if len(toks) >= 2 else 1
    best = None
    for f in files:
        low = f.lower()
        if module not in low or low.rsplit(".", 1)[0] in skip:
            continue
        if any(req not in low for req in required_terms):
            continue
        n = sum(1 for t in toks if t in low)
        if n < need:
            continue
        if re.search(r"(plot|heatmap)", low):
            kind = 3
        elif "table" in low:
            kind = 2
        elif low.startswith("multiqc_" + module):
            kind = 1
        else:
            kind = 0
        key = (n, kind, -len(low))
        if best is None or key > best[0]:
            best = (key, f)
    return os.path.join(data_dir, best[1]) if best else None


# --------------------------------------------------- inject one summary
def inject(html, sec_id, summary_html):
    # 1) MultiQC's own pre-created empty AI div (present when its AI feature is on)
    empty = rf'<div class="ai-summary-response" id="{re.escape(sec_id)}_ai_summary_response"[^>]*></div>'
    if re.search(empty, html):
        html = re.sub(
            empty,
            f'<div class="ai-summary-response" id="{sec_id}_ai_summary_response" '
            f'style="margin-bottom:-5px;">{summary_html}</div>',
            html,
            count=1,
        )
        # its wrapper ships hidden; reveal it since we pre-baked the answer
        html = re.sub(
            rf'(id="{re.escape(sec_id)}_ai_summary_wrapper"[^>]*style=")display:\s*none;',
            r"\1display: block;",
            html,
            count=1,
        )
        return html, True
    # 2) fallback: inject right after the section wrapper opening tag
    wrap = rf'(<div [^>]*id="mqc-section-wrapper-{re.escape(sec_id)}"[^>]*>)'
    if re.search(wrap, html):
        block = (
            f'\\1\n<div class="ai-summary-response" style="margin:10px 0;padding:8px 12px;'
            f'border-left:3px solid #4a90d9;background:#f5f9fd;">{summary_html}</div>'
        )
        return re.sub(wrap, block, html, count=1), True
    return html, False


# --------------------------------------------------------------- per-section
def per_section(cfg, html_path, data_dir):
    with open(html_path, encoding="utf-8", errors="surrogateescape") as fh:
        html = fh.read()
    sec_ids = discover_sections(html)
    llm_common.log(f"[MultiQC AI] {len(sec_ids)} sections found in report")

    mappings_log = ["Section ID -> Matched Data File Mappings:", "=" * 60]
    for sec_id in sorted(sec_ids):
        f = data_file_for_section(data_dir, sec_id)
        status = os.path.basename(f) if f else "SKIPPED (no matching data file)"
        mappings_log.append(f"  - {sec_id:<45} -> {status}")
    llm_common.log("[MultiQC AI] Section -> Data File Mappings:\n" + "\n".join(mappings_log[2:]))
    with open(os.path.join(data_dir, "section_mappings.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(mappings_log))

    prompts_log = [f"Per-section AI prompts -- model: {cfg.model}", ""]
    n_injected = 0
    timed_out = False  # once the LLM times out we stop calling and mark the rest
    for sec_id in sorted(sec_ids):
        f = data_file_for_section(data_dir, sec_id)
        if not f:
            llm_common.log(f"[MultiQC AI] no data file for {sec_id}; skipping")
            continue
        if timed_out:
            # The LLM already timed out: every remaining box shows the
            # placeholder rather than a stale or missing summary.
            html, _ = inject(html, sec_id, llm_common.timeout_box())
            continue
        with open(f, encoding="utf-8", errors="replace") as fh:
            data = fh.read()
        try:
            text, usage = ask_llm(cfg, sec_id, data)
        except llm_common.LLMTimeout:
            llm_common.log(f"[MultiQC AI] TIMEOUT on {sec_id}; marking this and remaining sections")
            timed_out = True
            html, _ = inject(html, sec_id, llm_common.timeout_box())
            continue
        # Refuse an implausibly long (degraded / hallucinated) answer: omit the box.
        if llm_common.is_degraded(text):
            llm_common.log(f"[MultiQC AI] {sec_id}: response too long ({len(text)} chars); omitting")
            continue
        html, ok = inject(
            html, sec_id, llm_common.md_to_html(text) + "\n" + llm_common.ai_footer(cfg.model, usage, data_filename=f)
        )
        n_injected += 1 if ok else 0
        llm_common.log(
            f"[MultiQC AI] {sec_id}: {'injected' if ok else 'NO ANCHOR'} "
            f"({usage.get('duration_s', 0):.1f}s, {usage.get('total_tokens', 0)} tok)"
        )
        prompts_log += ["=" * 76, f"Section: {sec_id}", "=" * 76, "", "[USER DATA FILE]", os.path.basename(f), ""]
    with open(os.path.join(data_dir, "section_prompts.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(prompts_log))
    with open(html_path, "w", encoding="utf-8", errors="surrogateescape") as fh:
        fh.write(html)
    llm_common.log(f"[MultiQC AI] wrote {html_path} ({n_injected} section summaries injected)")


# --------------------------------------------------------------------- main
class Config:
    """Resolved LLM parameters threaded through the run."""

    def __init__(self, a):
        self.endpoint = a.llm_endpoint
        self.model = a.llm_model
        self.api_key = a.llm_api_key
        self.context_window = a.llm_context_window
        self.timeout = a.llm_timeout
        # system prompt: --sys-prompt-file wins over --sys-prompt over default
        if a.sys_prompt_file:
            with open(a.sys_prompt_file, encoding="utf-8") as fh:
                self.sys_prompt = fh.read()
        else:
            self.sys_prompt = a.sys_prompt


def main():
    ap = argparse.ArgumentParser(
        description="Add per-section AI summaries to an existing MultiQC report.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--report", required=True, help="multiqc_report.html to annotate in place")
    ap.add_argument("--data-dir", required=True, help="the matching multiqc_data directory")
    ap.add_argument(
        "--redact-only",
        action="store_true",
        help="skip the LLM pass; only scrub the endpoint/key from the report (used when "
        "MultiQC's builtin AI ran but per-section annotation was not requested)",
    )

    g = ap.add_argument_group("LLM parameters (supply endpoint/model at runtime)")
    g.add_argument(
        "--llm-endpoint",
        default=DEF_ENDPOINT,
        help="REQUIRED: OpenAI-compatible /v1/chat/completions URL (env LLM_ENDPOINT)",
    )
    g.add_argument("--llm-model", default=DEF_MODEL, help="REQUIRED: model name (env LLM_MODEL)")
    g.add_argument(
        "--llm-api-key", default=DEF_API_KEY, help="Bearer token; 'dummy' for local servers (env LLM_API_KEY)"
    )
    g.add_argument(
        "--llm-context-window",
        type=int,
        default=DEF_CONTEXT_WINDOW,
        help="context window in tokens (env LLM_CONTEXT_WINDOW)",
    )
    g.add_argument(
        "--llm-timeout", type=int, default=DEF_TIMEOUT, help="per-request timeout, seconds (env LLM_TIMEOUT)"
    )
    g.add_argument("--sys-prompt", default=DEF_SYS_PROMPT, help="system prompt for per-section calls")
    g.add_argument("--sys-prompt-file", help="read system prompt from a file (wins over --sys-prompt)")

    a = ap.parse_args()
    cfg = Config(a)

    # Redaction must happen whatever else does or does not, so it is never
    # conditional on the LLM being reachable.
    targets = [a.report, a.data_dir]
    try:
        if not a.redact_only:
            if not cfg.endpoint or not cfg.model:
                llm_common.log(
                    "[MultiQC AI] LLM endpoint/model not configured; leaving the report unannotated."
                )
            else:
                per_section(cfg, a.report, a.data_dir)
    except Exception as e:
        # An unannotated report is an acceptable outcome; a failed pipeline is not.
        llm_common.log(f"[MultiQC AI] error during AI pass: {e}; leaving the report as MultiQC wrote it.")
    finally:
        # per_section writes only the model (kept) + QC data, never the endpoint,
        # but MultiQC's builtin AI bakes the endpoint into the HTML, so scrub
        # unconditionally -- and last, so nothing written above can survive.
        llm_common.redact_tree(targets, endpoint=cfg.endpoint, api_key=cfg.api_key, label="MultiQC AI")
    return 0


if __name__ == "__main__":
    sys.exit(main())
