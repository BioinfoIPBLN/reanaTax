#!/usr/bin/env python3
"""
llm_insight.py - generate a short, human-readable AI insight for one or more
pipeline artifacts (the combined Bracken abundance table, the combined Kraken2
report, ...) and write it next to them so the reports can embed it.

Ported from reanalyzerGSE (scripts/llm_insight.py); the DGE/enrichment tasks of
the original are replaced by the tasks that make sense here, and the optional
--mqc-out writes a MultiQC custom-content section so the insight lands at the
top of the MultiQC report instead of in a file nobody opens.

    LLM_ENDPOINT          OpenAI-compatible /v1/chat/completions URL  (or --endpoint)
    LLM_MODEL             model name                                  (or --model)
    LLM_API_KEY           Bearer token ("dummy" for servers that ignore it)
    LLM_CONTEXT_WINDOW    approx token budget, used to cap input size

Usage:
    llm_insight.py --input <file> [<file> ...] --task <taxonomy|generic> \
                   [--title "..."] [--out <file>] [--mqc-out <file>] [--max-rows N]

With several --input files the script runs a sequential MAP-REDUCE: it summarises
each file in its own LLM call (one at a time), then makes a final call that
synthesises those summaries into a single narrative. Every call is text-only and
serialised (only one query hits the LLM at a time).

It is opt-in and non-fatal by construction: with no endpoint/model configured it
exits 0 without writing anything, and any runtime error (network, HTTP, bad
response) is also non-fatal, so the reports simply omit the insight box rather
than the pipeline failing. Exit code 42 means the LLM timed out.
"""
import argparse
import html as html_mod
import os
import sys
import time

# Import the sibling shared module regardless of how this script is invoked.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import llm_common  # noqa: E402

# One shared instruction block, plus a task-specific framing. All prompts insist
# on brevity and on not inventing anything that is not in the supplied data.
COMMON = (
    "You are a careful bioinformatics assistant. Base every statement only "
    "on the data provided; never invent taxa, numbers or samples. "
    "Write plain prose (no markdown headings, tables or bullet points), a few "
    "short sentences at most. Your text is shown to researchers under an "
    "'AI summary (verify)' label, so hedge biological interpretation. "
    "Always respond in English."
)

TASKS = {
    "taxonomy": (
        "Given a combined taxonomic abundance table (taxa as rows, samples as columns, from "
        "Kraken2/Bracken applied to the NON-HOST fraction of the libraries), summarise in 3-5 "
        "sentences: which taxa dominate overall, how consistent the profiles are across samples, "
        "any sample that stands out, and whether anything looks like a likely reagent/index "
        "contaminant or a database artefact rather than real signal. Name only taxa that appear "
        "in the table. Make clear this is a hypothesis to be verified, not a conclusion."
    ),
    "generic": "Summarise the following bioinformatics output concisely and factually.",
}

# For map-reduce over SEVERAL tables: summarise each on its own (MAP), then
# synthesise all summaries (REDUCE). Both stay text-only and SEQUENTIAL.
MAP = {
    "taxonomy": (
        "Summarise this SINGLE taxonomic abundance table in 1-2 sentences: name the most abundant "
        "taxa and say how uniform the profile is across the samples present. Use only the rows "
        "given. If the table is empty or nothing was classified, say so briefly."
    ),
}
REDUCE = {
    "taxonomy": (
        "You are given short summaries of several taxonomic abundance tables produced from the same "
        "set of samples (e.g. a Kraken2 report and a Bracken re-estimation). Synthesise them into "
        "3-5 sentences describing the dominant taxa, the consistency across samples, and anything "
        "that warrants a closer look. Make clear this is a hypothesis to be verified, not a "
        "conclusion. Do not introduce taxa that are not present in the summaries."
    ),
}


def _map_prompt(task):
    return MAP.get(task, TASKS[task])


def _reduce_prompt(task):
    return REDUCE.get(task, "Synthesise the following per-item summaries into a single concise, factual paragraph.")


def _row_weight(line):
    """Total of every numeric field in a row. Used to rank taxa by abundance
    without having to know which column layout produced the table (Bracken's
    <sample>_num/_frac pairs and combine_kreports' <sample>_all/_lvl pairs both
    rank sensibly under this)."""
    total = 0.0
    for field in line.split("\t")[1:]:
        try:
            total += float(field)
        except ValueError:
            continue
    return total


def resolve_rows(raw, task, max_rows):
    """Return a bounded slice of the artifact to send. Taxonomic tables are
    ranked by abundance first, so if the context window has to trim anything it
    is the long tail of near-zero taxa rather than the dominant organisms."""
    lines = raw.splitlines()
    if not lines:
        return raw
    if task == "taxonomy":
        # Both table formats put their header(s) first: Bracken a single plain
        # header row, combine_kreports a run of '#'-prefixed lines.
        head, body = [], []
        for i, line in enumerate(lines):
            if not body and (line.startswith("#") or i == 0):
                head.append(line)
            elif line.strip():
                body.append(line)
        body.sort(key=_row_weight, reverse=True)
        sel = body if max_rows <= 0 else body[:max_rows]
        note = (
            f"# Data note: {len(body)} taxa in the table; showing the {len(sel)} most abundant."
            if len(sel) < len(body)
            else f"# Data note: all {len(body)} taxa in the table are shown."
        )
        return "\n".join([note] + head + sel)
    return "\n".join(lines if max_rows <= 0 else lines[: max_rows + 1])


def _read_capped(path, task, max_rows, context_window):
    """Read one artifact, keep the rows we care about (still sorted so the most
    abundant survive if the context window trims)."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        raw = fh.read()
    return llm_common.cap_text(resolve_rows(raw, task, max_rows), context_window, floor=2000)


def _call_tps(u):
    """Per-call tokens-per-second from one usage dict (0 if no timing)."""
    sec = float((u or {}).get("duration_s", 0.0) or 0.0)
    return (llm_common.answer_tokens(u) / sec) if sec > 0 else 0.0, sec


def _summarize_one(a, path):
    """Single artifact -> one query -> (summary text, tps, seconds)."""
    data = _read_capped(path, a.task, a.max_rows, a.context_window)
    user = f"{TASKS[a.task]}\n\n" + (f"Context: {a.title}\n\n" if a.title else "") + f"Data:\n{data}"
    text, u = llm_common.chat_completion(
        a.endpoint,
        a.model,
        a.api_key,
        [{"role": "system", "content": COMMON + "\n\n" + TASKS[a.task]}, {"role": "user", "content": user}],
        timeout=a.timeout,
    )
    tps, sec = _call_tps(u)
    return (text or "").strip(), tps, sec


def _map_reduce(a, paths):
    """MANY tables -> one sequential query PER table, then one query that
    synthesises the per-table summaries into a single narrative. The loop is
    serial by construction: the next query starts only after the previous
    returns. Returns (narrative, avg per-call tps, total model seconds)."""
    summaries = []
    tps_calls = []
    seconds = 0.0
    for i, path in enumerate(paths, 1):
        data = _read_capped(path, a.task, a.max_rows, a.context_window)
        name = os.path.basename(path)
        if len([line for line in data.splitlines() if line.strip() and not line.startswith("#")]) <= 1:
            llm_common.log(f"[llm_insight] skip {i}/{len(paths)}: {name} (no rows)")
            continue
        user = f"{_map_prompt(a.task)}\n\nTable: {name}\n\nData:\n{data}"
        text, u = llm_common.chat_completion(
            a.endpoint,
            a.model,
            a.api_key,
            [{"role": "system", "content": COMMON + "\n\n" + _map_prompt(a.task)}, {"role": "user", "content": user}],
            timeout=a.timeout,
        )
        tps, sec = _call_tps(u)
        seconds += sec
        if tps > 0:
            tps_calls.append(tps)
        # A single per-table summary should be short; drop a runaway one so it
        # cannot poison the reduce step.
        text = "" if llm_common.is_degraded(text) else (text or "").strip()
        llm_common.log(
            f"[llm_insight] mapped {i}/{len(paths)}: {name} ({tps:.0f} TPS, {sec:.0f}s)"
            f"{'' if text else ' (empty)'}"
        )
        if text:
            summaries.append((name, text))
    avg_tps = (sum(tps_calls) / len(tps_calls)) if tps_calls else 0.0
    if not summaries:
        return "", avg_tps, seconds
    joined = "\n".join(f"- {n}: {t}" for n, t in summaries)
    user = (
        f"{_reduce_prompt(a.task)}\n\n"
        + (f"Context: {a.title}\n\n" if a.title else "")
        + f"Per-table summaries:\n{joined}"
    )
    text, u = llm_common.chat_completion(
        a.endpoint,
        a.model,
        a.api_key,
        [{"role": "system", "content": COMMON + "\n\n" + _reduce_prompt(a.task)}, {"role": "user", "content": user}],
        timeout=a.timeout,
    )
    tps, sec = _call_tps(u)
    seconds += sec
    if tps > 0:
        tps_calls.append(tps)
    return (text or "").strip(), (sum(tps_calls) / len(tps_calls)) if tps_calls else 0.0, seconds


def _mqc_description(inputs):
    """The one-line blurb MultiQC prints under the section heading. Escaped
    because it is embedded in the HTML comment MultiQC parses as YAML."""
    tables = ", ".join(os.path.basename(f) for f in inputs)
    return html_mod.escape(f"{llm_common.AI_DISCLAIMER}. Built from: {tables}.")


def write_mqc(path, section_id, section_name, description, body_html):
    """MultiQC custom content: an HTML fragment whose leading comment carries the
    section metadata. MultiQC picks up any *_mqc.html handed to it."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(
            "<!--\n"
            f"id: '{section_id}'\n"
            f"section_name: '{section_name}'\n"
            f"description: '{description}'\n"
            "-->\n"
            f"{body_html}\n"
        )


def main():
    p = argparse.ArgumentParser(description="Write a short AI insight for one or more pipeline tables.")
    p.add_argument("--input", required=True, nargs="+", help="one or more tables; multiple => sequential map-reduce")
    p.add_argument("--task", default="generic", choices=list(TASKS))
    p.add_argument("--title", default="")
    p.add_argument("--out", default=None, help="markdown output (default: <first input>.ai_insight.md)")
    p.add_argument("--mqc-out", default=None, help="also write a MultiQC custom-content *_mqc.html section here")
    p.add_argument("--mqc-section-name", default="AI summary: taxonomic profile")
    p.add_argument("--max-rows", type=int, default=50, help="rows kept per table, ranked by abundance (0 = all)")
    p.add_argument("--endpoint", default=os.environ.get("LLM_ENDPOINT"))
    p.add_argument("--model", default=os.environ.get("LLM_MODEL"))
    p.add_argument("--api-key", default=os.environ.get("LLM_API_KEY", "dummy"))
    p.add_argument("--context-window", type=int, default=int(os.environ.get("LLM_CONTEXT_WINDOW") or 128000))
    p.add_argument("--timeout", type=int, default=int(os.environ.get("LLM_TIMEOUT") or 300))
    a = p.parse_args()

    # Opt-in: silently do nothing if no LLM is configured.
    if not a.endpoint or not a.model:
        llm_common.log("[llm_insight] no LLM endpoint/model configured; skipping.")
        return 0
    inputs = [f for f in a.input if os.path.isfile(f)]
    if not inputs:
        llm_common.log(f"[llm_insight] no input file found ({', '.join(a.input)}); skipping.")
        return 0

    out = a.out or (inputs[0] + ".ai_insight.md")
    try:
        t0 = time.time()
        # One artifact -> single query; many -> sequential map-reduce (all serialised).
        text, tps, seconds = _summarize_one(a, inputs[0]) if len(inputs) == 1 else _map_reduce(a, inputs)
        if not text:
            llm_common.log(f"[llm_insight] empty response for {out}; nothing written.")
            return 0
        # Guard against LLM degradation / hallucination: a wildly long "answer"
        # (garbage filler) is refused rather than shown.
        if llm_common.is_degraded(text):
            llm_common.log(
                f"[llm_insight] response for {out} exceeds {llm_common.max_answer_chars()} chars "
                f"(likely degraded); nothing written."
            )
            return 0
        stamp = time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime())
        # Model name is kept (useful provenance); the endpoint is never included.
        header = (
            f"**AI summary** — {llm_common.AI_DISCLAIMER} — Generated by `{a.model}` on {stamp}, "
            f"{tps:.0f} TPS, {seconds:.0f} seconds"
        )
        with open(out, "w", encoding="utf-8") as fh:
            fh.write(header + "\n\n" + text.strip() + "\n")
        if a.mqc_out:
            body = llm_common.md_to_html(text) + llm_common.ai_footer(a.model, {}, tps=tps, seconds=seconds)
            write_mqc(a.mqc_out, "reanatax-ai-insight", a.mqc_section_name, _mqc_description(inputs), body)
        llm_common.log(
            f"[llm_insight] wrote {out} ({time.time() - t0:.1f}s, {len(inputs)} table(s), "
            f"{tps:.0f} TPS avg, {seconds:.0f}s total)"
        )
    except llm_common.LLMTimeout as e:
        # Say it was tried and timed out rather than silently omitting the box,
        # so the user knows there is something to retry.
        llm_common.log(f"[llm_insight] TIMEOUT for {out}: {e}")
        with open(out, "w", encoding="utf-8") as fh:
            fh.write(f"**AI summary** — {llm_common.AI_TIMEOUT_MSG}.\n")
        if a.mqc_out:
            write_mqc(
                a.mqc_out,
                "reanatax-ai-insight",
                a.mqc_section_name,
                _mqc_description(inputs),
                llm_common.timeout_box(),
            )
        return 0
    except Exception as e:
        llm_common.log(f"[llm_insight] error for {out}: {e}; skipping.")
        return 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
