//
// Post-hoc AI annotation of the HTML reports: per-section summaries in the
// MultiQC report and a summary box in every Qualimap report.
//
// This subworkflow also owns the pipeline's rule that NO TWO LLM CALLS EVER
// OVERLAP. The LLM this talks to is usually a single self-hosted server, and
// hitting it from several tasks at once either queues badly or fails outright.
// Serialisation has two halves:
//
//   1. Within a process: `maxForks = 1` on the AI processes (conf/base.config),
//      so the per-sample QUALIMAP_AI tasks run one after another.
//   2. Between processes: an explicit gate. QUALIMAP_AI's input is combined
//      with a value channel that only resolves once MULTIQC_AI has finished,
//      so the two can never be in flight together.
//
// The rest of the chain is already ordered by the data: LLM_INSIGHT feeds
// MULTIQC, which feeds MULTIQC_AI.
//
// bin/llm_common.py additionally serialises with an in-process lock and a
// flock, but that only covers one host, so it is a backstop and not the
// mechanism relied on here.
//

include { MULTIQC_AI  } from '../../../modules/local/multiqc/ai/main'
include { QUALIMAP_AI } from '../../../modules/local/qualimap/ai/main'

workflow AI_ANNOTATE_REPORTS {

    take:
    ch_multiqc_report // channel: [ val(meta), path(html) ]
    ch_multiqc_data // channel: [ val(meta), path(multiqc_data) ]
    ch_qualimap // channel: [ val(meta), path(qualimap_results_dir) ]
    llm_endpoint // string
    llm_model // string
    llm_api_key // string
    annotate_multiqc // boolean: add per-section summaries (otherwise redact only)
    annotate_qualimap // boolean

    main:

    //
    // MODULE: MultiQC report.
    //
    // This runs whenever an endpoint is configured, even when per-section
    // annotation was not requested, because it is also the pipeline's redaction
    // point: the run's parameter summary - and, with --multiqc_ai_builtin,
    // MultiQC's own AI section - would otherwise carry the endpoint URL into a
    // published, shareable HTML file.
    //
    MULTIQC_AI(
        ch_multiqc_report.join(ch_multiqc_data),
        llm_endpoint,
        llm_model,
        llm_api_key,
        annotate_multiqc,
    )

    //
    // The serialisation gate: resolves as soon as MULTIQC_AI is done, or
    // immediately if MULTIQC_AI never ran (--skip_multiqc). `concat` is what
    // makes the fallback safe - it is only reached once the first channel has
    // closed, so a real MULTIQC_AI result always wins the race with it.
    //
    def ch_gate = MULTIQC_AI.out.report
        .map { _meta, _report -> true }
        .concat(channel.of(true))
        .first()

    //
    // MODULE: Qualimap reports, one task per sample, one at a time.
    //
    def ch_qualimap_ai = channel.empty()

    if (annotate_qualimap) {
        QUALIMAP_AI(
            ch_qualimap.combine(ch_gate).map { meta, results, _gate -> [meta, results] },
            llm_endpoint,
            llm_model,
            llm_api_key,
        )
        ch_qualimap_ai = QUALIMAP_AI.out.results
    }

    emit:
    multiqc_report = MULTIQC_AI.out.report // channel: [ val(meta), path(html) ]
    multiqc_data = MULTIQC_AI.out.data // channel: [ val(meta), path(multiqc_data) ]
    qualimap = ch_qualimap_ai // channel: [ val(meta), path(qualimap_results_dir) ]
}
