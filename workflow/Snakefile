rule all:
    input:
        "results/Analise_completa_genes_estresse.pdf"


rule microarray_analysis:

    input:
        "data/Transcripts_ID_RcNAC.csv"

    output:
        "results/Analise_completa_genes_estresse.pdf"

    container:
        "rcnac_pipeline"

    shell:
        """
        Rscript scripts/microarray_analysis.R {input} results
        """
