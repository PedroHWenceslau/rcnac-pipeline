process ANALYSE_GENES {

    container 'rcnac_pipeline'

    publishDir 'results', mode: 'copy'

    input:
    path transcript_data

    output:
    path "Analise_completa_genes_estresse.pdf"

    script:
    """
    Rscript scripts/microarray_analysis.R ${transcript_data} .
    """
}

workflow {

    transcript_data = file("data/Transcripts_ID_RcNAC.csv")

    ANALYSE_GENES(transcript_data)

}
