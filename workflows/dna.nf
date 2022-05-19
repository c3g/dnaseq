/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowDna.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.readsets, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.readsets)   { ch_input = file(params.readsets) } else { exit 1, 'Input readset not specified!' }

// Check alignment parameters
def prepareToolIndices  = ['dict', 'fai', 'intervals']
if (!params.skip_alignment) { prepareToolIndices << params.aligner }



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTP                                      } from '../modules/local/fastp/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK            } from '../subworkflows/local/input_check'
include { PREPARE_GENOME         } from '../subworkflows/local/prepare_genome'
include { BWA_MEM                } from '../modules/local/bwa/mem/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                                          } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                      } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { SAMTOOLS_INDEX                                   } from '../modules/nf-core/modules/samtools/index/main'
include { PICARD_MERGESAMFILES                             } from '../modules/nf-core/modules/picard/mergesamfiles/main'
include { PICARD_SORTSAM as PICARD_SORT_BEFORE_DEDUP       } from '../modules/nf-core/modules/picard/sortsam/main'
include { PICARD_MARKDUPLICATES                            } from '../modules/nf-core/modules/picard/markduplicates/main'
include { PICARD_SORTSAM as PICARD_SORT_AFTER_DEDUP        } from '../modules/nf-core/modules/picard/sortsam/main'
include { GATK4_BASERECALIBRATOR                           } from '../modules/nf-core/modules/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR                                  } from '../modules/nf-core/modules/gatk4/applybqsr/main'
include { GATK4_HAPLOTYPECALLER                            } from '../modules/nf-core/modules/gatk4/haplotypecaller/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FOR_HAPLOTYPING } from '../modules/nf-core/modules/samtools/index/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow DNA {

    ch_versions = Channel.empty()

    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    INPUT_CHECK(ch_input)

    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    PREPARE_GENOME( prepareToolIndices )

    // Module: Run Fastp
    FASTP(INPUT_CHECK.out.reads, false, false)
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    // Module: Run bwa mem
    BWA_MEM (FASTP.out.reads, PREPARE_GENOME.out.bwa_index, false).bam
    | map { meta, bam -> [ [id:meta.sample], bam] }
    | groupTuple()
    | PICARD_MERGESAMFILES // BAM files are merged per sample

    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())
    ch_versions = ch_versions.mix(PICARD_MERGESAMFILES.out.versions.first())


    PICARD_SORT_BEFORE_DEDUP (PICARD_MERGESAMFILES.out.bam, 'queryname').bam | PICARD_MARKDUPLICATES
    PICARD_SORT_AFTER_DEDUP (PICARD_MARKDUPLICATES.out.bam, 'coordinate')

    ch_versions = ch_versions.mix(PICARD_SORT_BEFORE_DEDUP.out.versions.first())
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())

    PICARD_SORT_AFTER_DEDUP.out.bam | SAMTOOLS_INDEX

    ch_recal_intervals = Channel.empty()
    if (params.intervals) {
        Channel.fromPath(params.intervals)
        | splitText()
        ch_recal_intervals.mix(Channel.fromPath(params.intervals))
    }

    // Base recalibration
    if(!params.skip_recalibration && params.recalibration_vcfs) {

        ch_bam_bai = PICARD_SORT_AFTER_DEDUP.out.bam.join( SAMTOOLS_INDEX.out.bai )
        ch_bam_bai_intervals = ch_bam_bai.combine( PREPARE_GENOME.out.interval_list_concat )

        ArrayList<String> vcfs = params.recalibration_vcfs.split(",")
        ArrayList<String> tbis = vcfs.collect{ it + ".tbi"}
        ch_recalibration_vcfs = Channel.fromPath(vcfs).toList()
        ch_recalibration_tbis = Channel.fromPath(tbis).toList()

        GATK4_BASERECALIBRATOR (
            ch_bam_bai_intervals,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.fai,
            PREPARE_GENOME.out.dict,
            ch_recalibration_vcfs,
            ch_recalibration_tbis
        )

        ch_bam_bai_recal = ch_bam_bai.join( GATK4_BASERECALIBRATOR.out.table )
        ch_bam_bai_recal_intervals = ch_bam_bai_recal.combine( PREPARE_GENOME.out.interval_list_concat )

        GATK4_APPLYBQSR (
            ch_bam_bai_recal_intervals,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.fai,
            PREPARE_GENOME.out.dict
        )

        ch_bam_for_haplotypecaller = GATK4_APPLYBQSR.out.bam

    } else if (!params.recalibration_vcfs) {
        log.warn "No recalibration VCFs provided - skipping GATK recalibration step. Use --skip_recalibration to supress this warning."
    } else {
        ch_bam_for_haplotypecaller = PICARD_SORT_AFTER_DEDUP.out.bam
    }

    ch_bam_for_haplotypecaller
    | join( SAMTOOLS_INDEX_FOR_HAPLOTYPING( ch_bam_for_haplotypecaller ).bai )
    | combine( PREPARE_GENOME.out.interval_list )
    | set { ch_haplotypecaller_inputs }

    ch_dbsnp_vcfs = Channel.fromPath( params.dbsnp ).flatten().first()
    ch_dbsnp_tbis = Channel.fromPath( params.dbsnp + ".idx" ).flatten().first()

    GATK4_HAPLOTYPECALLER (
        ch_haplotypecaller_inputs,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.dict,
        ch_dbsnp_vcfs,
        ch_dbsnp_tbis
    )

    GATK4_HAPLOTYPECALLER.out.vcf
    | join( GATK4_HAPLOTYPECALLER.out.tbi )
    | groupTuple
    | view

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // MODULE: MultiQC
    workflow_summary    = WorkflowDna.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map{meta, json -> json})
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.metrics.collectFile(name: 'markduplicates.metrics.txt'))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/