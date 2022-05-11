//
// Check input samplesheet and get read channels
//

include { READSETS_CHECK } from '../../modules/local/readsets_check'

workflow INPUT_CHECK {
    take:
    readsets // file: /path/to/samplesheet.csv

    main:
    READSETS_CHECK( readsets ).csv
    .splitCsv( header:true, sep:',')
    .branch {
        bams: it.BAM
        fastqs: it.FASTQ1
            return create_fastq_channel(it)
    }.set { reads }

    fastqs = with_readsets_per_sample(reads.fastqs)

    emit:
    reads = fastqs                   // channel: [ val(meta), [ reads ] ]
    versions = READSETS_CHECK.out.versions // channel: [ versions.yml ]
}

def with_readsets_per_sample(chan) {
    def readsets_per_sample = chan
    | map { meta, reads -> [meta.sample, reads] }
    | groupTuple
    | map { sample, reads -> [sample, [readsets_per_sample:reads.size()]] }

    chan
    | map { meta, reads -> [meta.sample, meta, reads] }
    | combine ( readsets_per_sample, by: 0 )
    | map { sample, meta1, reads, meta2 -> [meta1 + meta2, reads] }
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id            = row.Readset
    meta.sample        = row.Sample
    meta.library       = row.Library
    meta.runtype       = row.RunType
    meta.run           = row.Run
    meta.lane          = row.Lane
    meta.adapter1      = row.Adapter1
    meta.adapter2      = row.Adapter2
    meta.qualityoffset = row.QualityOffset
    meta.bed           = row.BED
    meta.single_end    = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.FASTQ1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.FASTQ1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.FASTQ1) ] ]
    } else {
        if (!file(row.FASTQ2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.FASTQ2}"
        }
        fastq_meta = [ meta, [ file(row.FASTQ1), file(row.FASTQ2) ] ]
    }
    return fastq_meta
}
