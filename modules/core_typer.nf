
process core_typer {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_*.csv", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), val(scheme)

    output:
    tuple val(sample_id), path("${sample_id}_allele_calls.csv"),  emit: allele_calls
    tuple val(sample_id), path("${sample_id}_core-typer_qc.csv"), emit: qc
    tuple val(sample_id), path("${sample_id}_core_typer_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: core_typer\\n"                >> ${sample_id}_core_typer_provenance.yml
    printf -- "  tools:\\n"                                  >> ${sample_id}_core_typer_provenance.yml
    printf -- "    - tool_name: core-typer\\n"               >> ${sample_id}_core_typer_provenance.yml
    printf -- "      tool_version: \$(core-typer -v 2>&1 | cut -d '-' -f 2)\\n" >> ${sample_id}_core_typer_provenance.yml
    printf -- "      parameters:\\n"                         >> ${sample_id}_core_typer_provenance.yml
    printf -- "        - parameter: --min-identity\\n"       >> ${sample_id}_core_typer_provenance.yml
    printf -- "          value: ${params.min_identity}\\n"   >> ${sample_id}_core_typer_provenance.yml
    printf -- "        - parameter: --min-coverage\\n"       >> ${sample_id}_core_typer_provenance.yml
    printf -- "          value: ${params.min_coverage}\\n"   >> ${sample_id}_core_typer_provenance.yml
    printf -- "        - parameter: --scheme\\n"             >> ${sample_id}_core_typer_provenance.yml
    printf -- "          value: ${params.scheme}\\n"         >> ${sample_id}_core_typer_provenance.yml

    
    core-typer \
        --threads ${task.cpus} \
	--R1 ${reads[0]} \
	--R2 ${reads[1]} \
	--scheme ${scheme} \
	--tmpdir tmp \
	--outdir core-typer-output \
	--no-cleanup

    cp core-typer-output/allele_calls.csv ${sample_id}_allele_calls.csv
    cp core-typer-output/qc.csv ${sample_id}_core-typer_qc.csv
    """
}
