include { GUNZIP as GUNZIP_FASTA                                } from '../../modules/nf-core/modules/gunzip/main'
include { UNTAR as UNTAR_BWA_INDEX                              } from '../../modules/nf-core/modules/untar/main'
include { PICARD_CREATESEQUENCEDICTIONARY                       } from '../../modules/nf-core/modules/picard/createsequencedictionary/main'
include { BWA_INDEX                                             } from '../../modules/nf-core/modules/bwa/index/main'
include { SAMTOOLS_FAIDX                                        } from '../../modules/nf-core/modules/samtools/faidx/main'

include { PICARD_SCATTERINTERVALSBYNS       } from '../../modules/local/picard/scatterintervalsbyns/main'
include { GATK4_CONCAT_INTERVALLIST         } from '../../modules/local/gatk4/concatintervallist/main'

workflow PREPARE_GENOME {
    take:
    prepare_tool_indices

    main:
    ch_versions = Channel.empty()
    if (params.fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], params.fasta ] ).gunzip.map { it[1] }.first()
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = file(params.fasta)
    }

    // Uncompress bwa index or generate from scratch if required
    ch_bwa_index = Channel.empty()
    if ('bwa' in prepare_tool_indices) {
        if (params.bwa) {
            if (params.bwa.endsWith('.tar.gz')) {
                ch_bwa_index = UNTAR_BWA_INDEX ( params.bwa ).untar
                ch_versions = ch_versions.mix(UNTAR_BWA_INDEX.out.versions)
            } else {
                ch_bwa_index = file(params.bwa)
            }
        } else {
            ch_bwa_index = BWA_INDEX ( ch_fasta ).index
            ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
        }
    }

    ch_dict_index = Channel.empty()
    if ('dict' in prepare_tool_indices) {
        if (params.dict_index) {
            ch_dict_index = file(params.dict_index)
        } else {
            Channel.from(ch_fasta)
            | map { [[id:it.baseName], it]}
            | PICARD_CREATESEQUENCEDICTIONARY
            ch_dict_index = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict.map { meta, path -> path }
            ch_versions = ch_versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions)
        }
    }

    ch_fai_index = Channel.empty()
    if ('fai' in prepare_tool_indices) {
        if (params.fai_index) {
            ch_fai_index = file(params.fai_index)
        } else {
            Channel.from(ch_fasta)
            | map { [[id:it.baseName], it]}
            | SAMTOOLS_FAIDX
            ch_fai_index = SAMTOOLS_FAIDX.out.fai.map { meta, path -> path }
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        }
    }

    ch_interval_list = Channel.empty()
    if ('intervals' in prepare_tool_indices) {
        if (params.intervals) {
            ch_interval_list = Channel.fromPath(params.intervals)
            if (params.scatter_count) {
                Channel.fromPath(params.intervals)
            }
        } else {
            Channel.from(ch_fasta)
            | map { [[id:it.baseName], it]}
            | join( SAMTOOLS_FAIDX.out.fai )
            | join( PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict )
            | PICARD_SCATTERINTERVALSBYNS

            ch_interval_list = PICARD_SCATTERINTERVALSBYNS.out.interval_list.map { meta, path -> path }
            ch_versions = ch_versions.mix(PICARD_SCATTERINTERVALSBYNS.out.versions)
        }
        ch_interval_list.toList() | map { [[id:"all"], it]} | GATK4_CONCAT_INTERVALLIST
        ch_interval_list_concat = GATK4_CONCAT_INTERVALLIST.out.interval_list.map { meta, path -> path }
    }

    emit:
    fasta                = ch_fasta                       //    path: genome.fasta (singleton channel)
    fai                  = ch_fai_index.first()           //    path: genome.fasta.fai (singleton channel)
    dict                 = ch_dict_index.first()          //    path: genome.fasta.dict (singleton channel)
    bwa_index            = ch_bwa_index                   //    path: bwa/index/
    interval_list        = ch_interval_list               //    path: genome.interval_list
    interval_list_concat = ch_interval_list_concat        //    path: all.interval_list (singleton channel)
    // -- //
    versions = ch_versions.ifEmpty(null)                  // channel: [ versions.yml ]
}