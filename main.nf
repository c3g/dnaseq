#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    robsyme/dna
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/robsyme/dna
    Website: https://github.com/robsyme/dna
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta      = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.bwa        = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.hapmap     = WorkflowMain.getGenomeAttribute(params, 'hapmap')
params.hapmap_res = WorkflowMain.getGenomeAttribute(params, 'hapmap_res')
params.omni       = WorkflowMain.getGenomeAttribute(params, 'omni')
params.omni_res   = WorkflowMain.getGenomeAttribute(params, 'omni_res')
params.onekg      = WorkflowMain.getGenomeAttribute(params, 'onekg')
params.onekg_res  = WorkflowMain.getGenomeAttribute(params, 'onekg_res')
params.dbsnp      = WorkflowMain.getGenomeAttribute(params, 'dbsnp')
params.dbsnp_res  = WorkflowMain.getGenomeAttribute(params, 'dbsnp_res')


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DNA } from './workflows/dna'

workflow GENFLOW_DNA {
    DNA ()
}

workflow {
    GENFLOW_DNA ()
}
