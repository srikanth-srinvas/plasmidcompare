/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { TRIMGALORE             } from '../modules/local/trim_galore/main'
include { PLASMIDSPADES          } from '../modules/local/plasmidspades/main'
include { MOBRECON               } from '../modules/local/mobrecon/main'
include { QUAST                  } from '../modules/nf-core/quast/main'
include { PROKKA                 } from '../modules/nf-core/prokka/main'
include { DFAST                  } from '../modules/local/dfast/main'
include { BLASTN                 } from '../modules/local/blastn/main'
include { BRIG                   } from '../modules/local/brig/main'
include { MASH                   } from '../modules/local/mash/main'
include { FASTTREE               } from '../modules/local/fasttree/main'
include { RESFINDER              } from '../modules/local/resfinder/main'
include { CARD                   } from '../modules/local/card/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_plasmidcompare_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PLASMIDCOMPARE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // MODULE: Run FastQC
    FASTQC (
        samplesheet: ch_samplesheet,
        output: ch_fastqc
    )

    // MODULE: Run Trim Galore! (optional)
    if (params.run_trim_galore) {
        TRIMGALORE (
            input: ch_fastqc,
            output: ch_trimmed
        )
    } else {
        ch_trimmed = ch_fastqc
    }

    // MODULE: Plasmid Assembly
    if (params.assembled_genomes) {
        MOBRECON (
            input: ch_trimmed,
            output: ch_assembled
        )
    } else {
        PLASMIDSPADES (
            input: ch_trimmed,
            output: ch_assembled
        )
    }

    // MODULE: Calculate assembly statistics using QUAST
    QUAST (
        input: ch_assembled,
        output: ch_quast
    )

    // MODULE: Annotation
    if (params.use_prokka) {
        PROKKA (
            input: ch_assembled,
            output: ch_annotated
        )
    } else {
        DFAST (
            input: ch_assembled,
            output: ch_annotated
        )
    }

    // MODULE: Comparison and Visualization
    BLASTN (
        input: ch_annotated,
        output: ch_blastn
    )

    BRIG (
        input: ch_blastn,
        output: ch_brig
    )

    MASH (
        input: ch_annotated,
        output: ch_mash
    )

    if (params.run_phylogenetic_tree) {
        FASTTREE (
            input: ch_mash,
            output: ch_phylogenetic_tree
        )
    }

    // MODULE: Antibiotic Resistance Gene (ARG) Detection
    RESFINDER (
        input: ch_annotated,
        output: ch_resfinder
    )

    if (params.run_card) {
        CARD (
            input: ch_annotated,
            output: ch_card
        )
    }

    // MODULE: MultiQC
    MULTIQC (
        input: ch_multiqc_files,
        output: ch_multiqc_report
    )

    // Collect software versions
    softwareVersionsToYAML(ch_versions)

    // Generate methods description text
    methodsDescriptionText()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
