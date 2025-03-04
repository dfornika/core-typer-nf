#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { hash_files }            from './modules/hash_files.nf'
include { fastp }                 from './modules/fastp.nf'
include { core_typer }            from './modules/core_typer.nf'
include { pipeline_provenance }   from './modules/provenance.nf'
include { collect_provenance }    from './modules/provenance.nf'


if (params.profile){
    println("Profile should have a single dash: -profile")
    System.exit(1)
}

workflow {

    ch_workflow_metadata = Channel.value([
	workflow.sessionId,
	workflow.runName,
	workflow.manifest.name,
	workflow.manifest.version,
	workflow.start,
    ])

    if (params.samplesheet_input != 'NO_FILE') {
	ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], [it['R1'], it['R2']]] }
    } else {
	ch_fastq = Channel.fromFilePairs( params.fastq_illumina_search_path, flat: true ).map{ it -> [it[0].split('_')[0], [it[1], it[2]]] }.unique{ it -> it[0] }
    }

    ch_scheme = Channel.fromPath( "${params.scheme}")
    
    main:
    ch_sample_ids = ch_fastq.map{ it -> it[0] }

    hash_files(ch_fastq.combine(Channel.of("fastq-input")))

    fastp(ch_fastq)

    trimmed_reads = fastp.out.trimmed_reads
    
    core_typer(trimmed_reads.combine(ch_scheme))

    if (params.collect_outputs) {
	fastp.out.csv.map{ it -> it[1] }.collectFile(name: params.collected_outputs_prefix + "_fastp.csv", storeDir: params.outdir, keepHeader: true, sort: { it -> it.readLines()[1].split(',')[0] })
    }

    // Collect Provenance
    // The basic idea is to build up a channel with the following structure:
    // [sample_id, [provenance_file_1.yml, provenance_file_2.yml, provenance_file_3.yml...]]
    // At each step, we add another provenance file to the list using the << operator...
    // ...and then concatenate them all together in the 'collect_provenance' process.
    ch_provenance = ch_sample_ids

    ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)

    ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it -> [it[0], [it[1]]] }
    
    ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(fastp.out.provenance).map{ it ->      [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(core_typer.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

    collect_provenance(ch_provenance)
}
