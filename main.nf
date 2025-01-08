#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process SEQKIT {
    container "community.wave.seqera.io/library/seqkit:2.9.0--e0e29e1f5c28842a"

    publishDir "${params.outdir}/seqkit/", mode: 'copy'

    input:
    path reads

    output:
    path "reads_ids", emit: reads_ids

    script:
    def prefix = reads.baseName
    """
    seqkit seq -ni ${reads} > reads_ids
    """
}

process ALIGN_BATCH {
    container "community.wave.seqera.io/library/minimap2_python_tqdm_zstandard:d323a957c053da50"

    publishDir "${params.outdir}/minimap2/", mode: 'copy'

    input:
    path reads
    path read_ids
    path (batch_py, stageAs: 'tmp/batch.py')

    output:
    path "batches", emit: batches

    script:
    def prefix = reads.baseName
    """
    mkdir batches
    chmod +x tmp/batch.py
    minimap2 -K8g -cx ava-ont -k25 -w17 -e200 -r150 -m2500 -z200 -f 0.005 -t${task.cpus} --dual=yes $reads $reads | tmp/batch.py $read_ids - batches
    """
}

process ERROR_CORRECTION {
    container "public.ecr.aws/l0c7p3y3/lbcb-sci/herro:2025-01-07"

    publishDir "${params.outdir}/corrected/", mode: 'copy'

    input:
    path reads
    path (batches, stageAs: 'batches/*')
    path model

    output:
    path "*", emit: correct

    script:
    def prefix = reads.baseName
    """
    herro inference --read-alns batches/ -t 3 -d 0,1,2,3 -m $model -b 32 $reads ${prefix}_herro_output.fasta
    """
}

workflow {
    ch_reads = Channel.fromPath(params.reads)
    ch_batch_py = Channel.fromPath(params.batch_py)
    ch_model = Channel.fromPath(params.model)
    ch_batches = Channel.fromPath(params.batches)

    SEQKIT(reads_ch)

    ALIGN_BATCH(
        reads_ch,
        SEQKIT.out.reads_ids,
        batch_py
    )

    // ERROR_CORRECTION(
    //    ch_reads,
    //    ch_batches,
    //    ch_model
    //)
}
