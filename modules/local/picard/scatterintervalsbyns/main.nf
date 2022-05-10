process PICARD_SCATTERINTERVALSBYNS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::picard=2.26.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.9--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.9--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("*.interval_list"), emit: interval_list
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard ScatterIntervalsByNs] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard \\
        -Xmx${avail_mem}g \\
        ScatterIntervalsByNs  \\
        $args \\
        REFERENCE=$fasta \\
        OUTPUT_TYPE=ACGT \\
        OUTPUT=${prefix}.interval_list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard ScatterIntervalsByNs --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
