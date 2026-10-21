#!/usr/bin/env nextflow
// ^ "Shebang" line: lets the file run directly as ./main.nf. Optional, but conventional.

/*
 * ==========================================================================
 *  RNA-seq DE Results Cleaning Pipeline (Proof of Concept)
 * --------------------------------------------------------------------------
 *  Takes a CSV of differential expression results (e.g. DESeq2 output),
 *  removes rows with missing statistics, classifies significance, and
 *  produces a cleaned table plus a volcano plot.
 *
 *  Usage:
 *      nextflow run main.nf --input data/de_results.csv
 *      nextflow run main.nf --input data/de_results.csv --padj 0.01 --lfc 1.5
 * ==========================================================================
 */

// NOTE: DSL2 is the default since Nextflow 22.12, so the old
// `nextflow.enable.dsl = 2` line is no longer needed.


/*
 * --------------------------------------------------------------------------
 *  PARAMETERS
 * --------------------------------------------------------------------------
 *  Default parameter values (input, outdir, padj, lfc) now live in
 *  nextflow.config, so ALL settings are in one place. They are still
 *  available here as params.input, params.outdir, etc., and can still
 *  be overridden on the command line, e.g. --padj 0.01
 */


/*
 * --------------------------------------------------------------------------
 *  PROCESS: CLEAN_DE_RESULTS
 * --------------------------------------------------------------------------
 *  A `process` is one step of the pipeline. Nextflow runs each process in
 *  its own isolated folder under `work/`, so tasks never overwrite each
 *  other's files. Process names are UPPER_CASE by convention.
 */
process CLEAN_DE_RESULTS {

    // `tag` sets a label shown in the terminal log for each task.
    // Here it shows the input file name without its extension.
    tag "${de_results.baseName}"

    // `container` names the Docker image this process runs inside.
    // rocker/tidyverse is a public, widely used R image that already
    // contains dplyr, ggplot2 and readr, so we don't need to build our own.
    // The tag is PINNED to an exact R version (never use `latest`): for
    // past R releases, Rocker installs CRAN packages from a frozen
    // snapshot date, so everyone gets the same package versions.
    // This directive only takes effect when Docker is enabled
    // (see `docker.enabled` in nextflow.config).
    container 'rocker/tidyverse:4.4.3'

    // `publishDir` copies the declared outputs from the hidden work/
    // folder into a human-readable folder once the task succeeds.
    // mode: 'copy' makes real copies. The default is symlinks, which break
    // if work/ is deleted.
    publishDir params.outdir, mode: 'copy'

    // ---- INPUT block --------------------------------------------------
    // `path` tells Nextflow this value is a FILE. Nextflow stages
    // (symlinks) it into the task's work folder, and `de_results` becomes
    // a variable holding its file name.
    input:
    path de_results

    // ---- OUTPUT block -------------------------------------------------
    // Declares which files the script MUST create. If one is missing,
    // the task fails, which acts as a built-in sanity check.
    // `emit:` gives each output a name so later processes can use it as
    // CLEAN_DE_RESULTS.out.table, etc.
    output:
    path "cleaned_de_results.csv", emit: table
    path "significant_genes.csv",  emit: sig
    path "volcano_plot.pdf",       emit: plot

    // ---- SCRIPT block -------------------------------------------------
    // The triple-quoted string is a Bash command run inside the task folder.
    // ${de_results} is filled in by Nextflow with the staged input file.
    // Any script in the project's bin/ folder is automatically put on
    // $PATH, so we can call it by name with no absolute path needed.
    // (The file must be executable: chmod +x bin/clean_rna_seq.R)
    // The thresholds come from `params`, so they can be changed on the
    // command line (e.g. --padj 0.01 --lfc 1.5) without editing any code.
    script:
    """
    clean_rna_seq.R ${de_results} ${params.padj} ${params.lfc}
    """
}


/*
 * --------------------------------------------------------------------------
 *  WORKFLOW
 * --------------------------------------------------------------------------
 *  The unnamed `workflow` block is the entry point: it says how data
 *  flows into and between processes.
 */
workflow {

    // Fail early with a clear message if the user forgot --input.
    if (!params.input) {
        error "Please provide an input CSV, e.g. --input data/de_results.csv"
    }

    // ---- INPUT CHANNEL ------------------------------------------------
    // A `channel` is a queue that passes data into processes.
    // `channel.fromPath()` turns a path (or glob such as 'data/*.csv')
    // into a channel that emits one file per item.
    // `checkIfExists: true` stops the run right away if the file is
    // missing, instead of failing later inside a task.
    de_results_ch = channel.fromPath(params.input, checkIfExists: true)

    // ---- RUN THE PROCESS ----------------------------------------------
    // Calling a process like a function connects the channel to its input.
    // Nextflow starts one task for each item the channel emits.
    CLEAN_DE_RESULTS(de_results_ch)
}
