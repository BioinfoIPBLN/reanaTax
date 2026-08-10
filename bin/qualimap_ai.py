#!/usr/bin/env python3
"""
qualimap_ai.py -- parse Qualimap BamQC reports and add AI summaries via an
OpenAI-compatible endpoint.

Ported from reanalyzerGSE (scripts/qualimap_ai.py), trimmed to what this
pipeline produces: one BamQC report per sample for the host alignment. The
multi-sample branch of the original is dropped because nf-core has no
qualimap/multibamqc module; the RNA-seq branch is kept because a Qualimap
report directory is recognised by its contents, not by how it was made.

Every `qualimapReport.html` under --analysis-dir is annotated in place:
  - one AI summary box under the "Summary" heading, built from
    genome_results.txt (BamQC) or rnaseq_qc_results.txt (RNA-seq);
  - optionally (--sections) one box per plot, built from the matching table in
    raw_data_qualimapReport/. That is one extra LLM call per plot per sample,
    so it is off by default.

The LLM is any OpenAI-compatible /v1/chat/completions endpoint (local vLLM,
Ollama, OpenAI, ...). No endpoint or model is baked in -- supply them at runtime
via CLI flags or env vars (LLM_ENDPOINT, LLM_MODEL, LLM_API_KEY).
"""
import argparse
import os
import re
import sys
from html.parser import HTMLParser

# Import the sibling shared module regardless of how this script is invoked.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import llm_common  # noqa: E402

# Config defaults (CLI flag > env var > fallback)
DEF_ENDPOINT = os.environ.get("LLM_ENDPOINT")
DEF_MODEL = os.environ.get("LLM_MODEL")
DEF_API_KEY = os.environ.get("LLM_API_KEY", "dummy")
DEF_CONTEXT_WINDOW = int(os.environ.get("LLM_CONTEXT_WINDOW", "128000"))
DEF_TIMEOUT = int(os.environ.get("LLM_TIMEOUT", "300"))

# The BAM being described is the HOST alignment of a metagenomic library, so the
# prompt has to say that: otherwise the model reads a low genome-fraction or a
# low mapping rate as a failure, when here it is the point of the experiment.
CONTEXT = (
    "This BAM comes from a metagenomic pipeline: reads were aligned against the HOST genome "
    "purely to remove host sequence, and only the reads that did NOT align are carried forward "
    "for taxonomic classification. So a low genome fraction, uneven coverage and a low mapping "
    "rate are all expected and are not QC failures by themselves - they mean the library is not "
    "dominated by host. Comment on what the numbers say about the host fraction and about "
    "library quality (duplication, mapping quality, GC), not on whether the alignment 'worked'."
)

SYS_PROMPT_BAMQC = "\n".join(
    [
        "You are an expert bioinformatician evaluating a Qualimap BAM QC report for a single sample.",
        CONTEXT,
        "Produce 2-3 concise bullet points summarizing key alignment metrics (host mapping rate,",
        "duplication rate, mean coverage over the covered part of the genome, mapping quality).",
        "Use standard plain markdown list items (- bullet point).",
        "Base every statement only on the numbers provided; never invent values.",
        "State clearly if quality is satisfactory or if there are flags.",
        "Always respond in English.",
    ]
)

SYS_PROMPT_RNASEQ = "\n".join(
    [
        "You are an expert bioinformatician evaluating a Qualimap RNA-Seq QC report for a single sample.",
        "Produce 2-3 concise bullet points summarizing key observations.",
        "Cover total mapped reads & mapping rate, genomic origin breakdown (Exonic vs Intronic/Intergenic %),",
        "and 5'-3' transcript coverage bias.",
        "Use standard plain markdown list items (- bullet point).",
        "Base every statement only on the numbers provided; never invent values.",
        "State clearly if quality is satisfactory or if there are contamination/degradation flags.",
        "Always respond in English.",
    ]
)

SYS_PROMPT_SECTION = "\n".join(
    [
        "You are an expert bioinformatician evaluating a specific plot/section of a Qualimap report.",
        CONTEXT,
        "Produce 1-2 concise bullet points summarizing key observations from the provided plot raw data.",
        "Use markdown. Use **bold** for important findings.",
        "Base every statement only on the numbers provided; never invent values.",
        "Always respond in English.",
    ]
)

SECTION_MAP = {
    # BAM QC raw data files -> HTML section heading names
    "coverage_across_reference.txt": "Coverage across reference",
    "coverage_histogram.txt": "Coverage Histogram",
    "duplication_rate_histogram.txt": "Duplication Rate Histogram",
    "genome_fraction_coverage.txt": "Genome Fraction Coverage",
    "homopolymer_indels.txt": "Homopolymer Indels",
    "mapped_reads_clipping_profile.txt": "Mapped Reads Clipping Profile",
    "mapped_reads_gc-content_distribution.txt": "Mapped Reads GC-content Distribution",
    "mapped_reads_nucleotide_content.txt": "Mapped Reads Nucleotide Content",
    "mapping_quality_across_reference.txt": "Mapping Quality Across Reference",
    "mapping_quality_histogram.txt": "Mapping Quality Histogram",
    # RNA-Seq QC raw data files -> HTML section heading names
    "coverage_profile_along_genes_(total).txt": "Coverage Profile Along Genes (Total)",
    "coverage_profile_along_genes_(low).txt": "Coverage Profile Along Genes (Low)",
    "coverage_profile_along_genes_(high).txt": "Coverage Profile Along Genes (High)",
}


# ------------------------------------------------------------ HTML table parser
class HTMLTableTextExtractor(HTMLParser):
    """Fallback data source: flatten every <table> in the report to TSV text when
    the raw results file is missing."""

    def __init__(self):
        super().__init__()
        self.in_table = False
        self.in_cell = False
        self.current_table = []
        self.current_row = []
        self.current_cell = []
        self.tables = []

    def handle_starttag(self, tag, attrs):
        if tag == "table":
            self.in_table = True
            self.current_table = []
        elif tag in ("td", "th") and self.in_table:
            self.in_cell = True
            self.current_cell = []
        elif tag == "tr" and self.in_table:
            self.current_row = []

    def handle_endtag(self, tag):
        if tag == "table" and self.in_table:
            self.in_table = False
            if self.current_table:
                self.tables.append(self.current_table)
        elif tag in ("td", "th") and self.in_cell:
            self.in_cell = False
            self.current_row.append(" ".join("".join(self.current_cell).split()))
        elif tag == "tr" and self.in_table:
            if self.current_row:
                self.current_table.append(self.current_row)

    def handle_data(self, data):
        if self.in_cell:
            self.current_cell.append(data)


def extract_tables_as_text(html_content):
    parser = HTMLTableTextExtractor()
    parser.feed(html_content)
    lines = []
    for table in parser.tables:
        for row in table:
            lines.append("\t".join(row))
        lines.append("")
    return "\n".join(lines)


def read_report(report_dir, html_path, results_name):
    """Return (html, data_text). Prefer Qualimap's plain-text results file; fall
    back to the tables scraped out of the HTML itself."""
    with open(html_path, encoding="utf-8", errors="replace") as fh:
        html_content = fh.read()
    txt_path = os.path.join(report_dir, results_name)
    if os.path.exists(txt_path):
        with open(txt_path, encoding="utf-8", errors="replace") as fh:
            return html_content, fh.read()
    return html_content, extract_tables_as_text(html_content)


# ---------------------------------------------------------------- HTML injection
def inject_ai_box(html_content, summary_body_html, footer_html):
    # Ensure a UTF-8 charset meta tag is declared in <head> so browsers do not
    # fall back to Windows-1252 (which turns the robot emoji into mojibake).
    if not re.search(r"<meta\s+charset=", html_content, re.IGNORECASE) and not re.search(
        r'http-equiv=["\']content-type["\']', html_content, re.IGNORECASE
    ):
        if re.search(r"<head[^>]*>", html_content, re.IGNORECASE):
            html_content = re.sub(
                r"(<head[^>]*>)", r'\1\n\t<meta charset="utf-8">', html_content, count=1, flags=re.IGNORECASE
            )
        elif re.search(r"<html[^>]*>", html_content, re.IGNORECASE):
            html_content = re.sub(
                r"(<html[^>]*>)",
                r'\1\n<head><meta charset="utf-8"></head>',
                html_content,
                count=1,
                flags=re.IGNORECASE,
            )

    box_html = (
        f'\n<div class="ai-summary-box" style="background-color: #f8f9fa; border-left: 4px solid #007bff; '
        f'padding: 12px 16px; margin: 15px 0; border-radius: 4px; font-family: sans-serif;">\n'
        f'  <h4 style="margin-top:0; margin-bottom: 8px; color: #007bff; font-weight: 600;">\U0001f916 AI Summary</h4>\n'
        f'  <div class="ai-summary-content" style="font-size: 0.95em; color: #333;">\n{summary_body_html}\n  </div>\n'
        f"  {footer_html}\n"
        f"</div>\n"
    )

    if 'class="ai-summary-box"' in html_content:  # already injected
        return html_content

    pattern = r"(<h2[^>]*>\s*Summary\s*<.*?</h2>)"
    if re.search(pattern, html_content, re.IGNORECASE):
        return re.sub(pattern, r"\1" + box_html, html_content, count=1, flags=re.IGNORECASE)
    if '<div class="content">' in html_content:
        return html_content.replace('<div class="content">', '<div class="content">' + box_html, 1)
    return html_content + box_html


def inject_section_ai_box(html_content, heading_name, summary_body_html, footer_html):
    box_html = (
        f'\n<div class="ai-summary-box-section" style="background-color: #f8f9fa; border-left: 3px solid #17a2b8; '
        f'padding: 8px 12px; margin: 10px 0; border-radius: 4px; font-family: sans-serif; font-size: 0.9em;">\n'
        f'  <div class="ai-summary-content" style="color: #333;">\n{summary_body_html}\n  </div>\n'
        f"  {footer_html}\n"
        f"</div>\n"
    )
    pattern = r"(<h2[^>]*>\s*" + re.escape(heading_name) + r"[^>]*>.*?</h2>)"
    match = re.search(pattern, html_content, re.IGNORECASE)
    if match:
        snippet_after = html_content[match.end() : match.end() + 250]
        if "ai-summary-box" in snippet_after:
            return html_content
        return re.sub(pattern, r"\1" + box_html, html_content, count=1, flags=re.IGNORECASE)
    return html_content


# ----------------------------------------------------------------- LLM helper
def ask_llm(cfg, sys_prompt, data_text):
    data_text = llm_common.cap_text(data_text, cfg.llm_context_window)
    messages = [
        {"role": "system", "content": sys_prompt},
        {"role": "user", "content": f"Data:\n{data_text}"},
    ]
    try:
        text, usage = llm_common.chat_completion(
            cfg.llm_endpoint, cfg.llm_model, cfg.llm_api_key, messages, timeout=cfg.llm_timeout
        )
    except llm_common.LLMTimeout:
        raise
    except Exception as e:
        llm_common.log(f"[Qualimap AI] error: {e}")
        return "_(no AI response)_", {}
    return (text or "_(no AI response)_"), usage


def process_report(cfg, html_path, report_type):
    report_dir = os.path.dirname(html_path)
    sample_name = os.path.basename(report_dir) or "sample"
    llm_common.log(f"[Qualimap AI] Annotating {report_type} report for {sample_name}")

    if report_type == "rnaseq":
        parsed_name, sys_prompt = "rnaseq_qc_results.txt", SYS_PROMPT_RNASEQ
    else:
        parsed_name, sys_prompt = "genome_results.txt", SYS_PROMPT_BAMQC
    html_content, text_data = read_report(report_dir, html_path, parsed_name)

    if not text_data.strip():
        llm_common.log(f"[Qualimap AI] Warning: no text data for {sample_name}; skipping")
        return

    text, usage = ask_llm(cfg, f"{sys_prompt}\nSample: {sample_name}", text_data)
    if llm_common.is_degraded(text):
        llm_common.log(f"[Qualimap AI] Degraded response dropped for {sample_name}")
    else:
        html_content = inject_ai_box(
            html_content, llm_common.md_to_html(text), llm_common.ai_footer(cfg.llm_model, usage, parsed_name)
        )

    # Per-plot boxes: one extra LLM call per plot per sample, hence opt-in.
    raw_data_dir = os.path.join(report_dir, "raw_data_qualimapReport")
    if cfg.sections and os.path.isdir(raw_data_dir):
        for fname, heading in SECTION_MAP.items():
            fpath = os.path.join(raw_data_dir, fname)
            if not os.path.exists(fpath):
                continue
            try:
                with open(fpath, encoding="utf-8", errors="replace") as fh:
                    raw_text = fh.read().strip()
                if not raw_text:
                    continue
                # Cap large raw histograms (coverage_histogram.txt is often 400 kB+)
                lines = [line for line in raw_text.splitlines() if line.strip()]
                if len(lines) > 100:
                    raw_text = "\n".join(
                        lines[:50] + ["\n... [middle histogram data points omitted] ...\n"] + lines[-50:]
                    )
                sec_text, sec_usage = ask_llm(
                    cfg, f"{SYS_PROMPT_SECTION}\nSection: {heading}\nSample: {sample_name}", raw_text
                )
                if not llm_common.is_degraded(sec_text):
                    html_content = inject_section_ai_box(
                        html_content,
                        heading,
                        llm_common.md_to_html(sec_text),
                        llm_common.ai_footer(cfg.llm_model, sec_usage, fname),
                    )
            except llm_common.LLMTimeout:
                raise
            except Exception as e:
                llm_common.log(f"[Qualimap AI] Section {fname} error for {sample_name}: {e}")

    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(html_content)
    llm_common.log(f"[Qualimap AI] Done {sample_name}")


def discover_reports(analysis_dir):
    """Find every Qualimap report under `analysis_dir` and classify it by the
    results file that sits next to it."""
    reports = []
    for root, _dirs, files in os.walk(analysis_dir):
        if "qualimapReport.html" not in files:
            continue
        full = os.path.join(root, "qualimapReport.html")
        if "rnaseq_qc_results.txt" in files or "rnaseq" in root.lower():
            reports.append((full, "rnaseq"))
        else:
            reports.append((full, "bamqc"))
    return sorted(reports)


def main():
    p = argparse.ArgumentParser(description="Add AI summaries to Qualimap reports.")
    p.add_argument("--analysis-dir", default=".", help="root directory containing Qualimap report directories")
    p.add_argument(
        "--sections",
        action="store_true",
        help="also summarise every individual plot (one extra LLM call per plot per sample)",
    )
    p.add_argument("--llm-endpoint", default=DEF_ENDPOINT, help="OpenAI-compatible URL (env LLM_ENDPOINT)")
    p.add_argument("--llm-model", default=DEF_MODEL, help="model name (env LLM_MODEL)")
    p.add_argument("--llm-api-key", default=DEF_API_KEY, help="API key (env LLM_API_KEY)")
    p.add_argument("--llm-context-window", type=int, default=DEF_CONTEXT_WINDOW)
    p.add_argument("--llm-timeout", type=int, default=DEF_TIMEOUT)
    cfg = p.parse_args()

    if not cfg.llm_endpoint or not cfg.llm_model:
        llm_common.log("[Qualimap AI] LLM endpoint/model not configured; leaving reports unannotated.")
        return 0

    analysis_dir = os.path.abspath(cfg.analysis_dir)
    reports = discover_reports(analysis_dir)
    llm_common.log(f"[Qualimap AI] {len(reports)} Qualimap report(s) found under {analysis_dir}")

    for html_path, report_type in reports:
        try:
            process_report(cfg, html_path, report_type)
        except llm_common.LLMTimeout:
            # One timeout means the endpoint is not keeping up; annotating the
            # remaining samples would just repeat the wait.
            llm_common.log(f"[Qualimap AI] TIMEOUT on {html_path}; stopping the AI pass here.")
            break
        except Exception as e:
            llm_common.log(f"[Qualimap AI] error on {html_path}: {e}; leaving that report unannotated.")

    # Unconditional, and last: no endpoint/key may survive in a published report.
    llm_common.redact_tree(analysis_dir, endpoint=cfg.llm_endpoint, api_key=cfg.llm_api_key, label="Qualimap AI")
    return 0


if __name__ == "__main__":
    sys.exit(main())
