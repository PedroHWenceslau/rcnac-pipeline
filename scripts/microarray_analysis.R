################################################
# Microarray analysis pipeline
# RcNAC stress response
################################################

library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(VennDiagram)
library(RColorBrewer)

################################################
# 1) Input
################################################

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
outdir <- args[2]

dir.create(outdir, showWarnings = FALSE)

################################################
# 2) Load data
################################################

df <- read.csv(input_file, row.names = 1)

df <- as.data.frame(
  lapply(df, function(x) as.numeric(gsub(",", ".", x)))
)

df <- df[, colSums(is.na(df)) < nrow(df)]

expr <- as.matrix(df)

expr <- log2(expr + 1)

################################################
# 3) Define groups
################################################

group <- factor(gsub("\\.\\d+$", "", colnames(expr)))

################################################
# 4) Differential expression
################################################

design <- model.matrix(~0 + group)

colnames(design) <- levels(group)

fit <- lmFit(expr, design)

contrast.matrix <- makeContrasts(
  
  X6h20_vs_Dry = X6h20 - Dry,
  X6h25_vs_Dry = X6h25 - Dry,
  X6h35_vs_Dry = X6h35 - Dry,
  
  RP20_vs_Dry = RP20 - Dry,
  RP25_vs_Dry = RP25 - Dry,
  RP35_vs_Dry = RP35 - Dry,
  
  R2_20_vs_Dry = R2_20 - Dry,
  R2_25_vs_Dry = R2_25 - Dry,
  R2_35_vs_Dry = R2_35 - Dry,
  
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

################################################
# 5) Extract results
################################################

results_list <- list()

for (contrast in colnames(contrast.matrix)) {
  
  tab <- topTable(
    fit2,
    coef = contrast,
    adjust = "fdr",
    number = Inf
  )
  
  results_list[[contrast]] <- tab
  
}

################################################
# 6) Combine results
################################################

all_results <- bind_rows(
  lapply(names(results_list), function(x){
    
    tab <- results_list[[x]]
    
    tab$Gene <- rownames(tab)
    tab$Contrast <- x
    
    tab
    
  })
)

################################################
# 7) Identify DEGs
################################################

deg_table <- all_results %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

deg_genes <- unique(deg_table$Gene)

deg_matrix <- expr[deg_genes, ]

################################################
# 8) Heatmap DEGs
################################################

heatmap_data <- t(scale(t(deg_matrix)))

heatmap_plot <- pheatmap(
  heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  fontsize_row = 6,
  silent = TRUE
)

################################################
# 9) Expression summary
################################################

expr_df <- as.data.frame(t(expr))

expr_df$Condition <- group

expr_long <- expr_df %>%
  pivot_longer(
    cols = -Condition,
    names_to = "Gene",
    values_to = "Expression"
  )

expr_long <- expr_long %>%
  filter(Gene %in% deg_genes) %>%
  group_by(Gene, Condition) %>%
  summarise(
    Mean = mean(Expression),
    SE = sd(Expression)/sqrt(n()),
    .groups = "drop"
  )

################################################
# 10) Ranking genes
################################################

deg_rank <- deg_table %>%
  group_by(Gene) %>%
  summarise(
    Mean_logFC = mean(abs(logFC)),
    Min_FDR = min(adj.P.Val),
    Score = Mean_logFC * -log10(Min_FDR)
  ) %>%
  arrange(desc(Score))

write.csv(
  deg_rank,
  file.path(outdir,"ranking_genes_candidatos_qPCR.csv"),
  row.names = FALSE
)

################################################
# 11) Top genes
################################################

top_genes <- head(deg_rank$Gene,10)

write.table(
  top_genes,
  file.path(outdir,"top_genes_qpcr.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

################################################
# 12) Plot expression
################################################

expr_long_top <- expr_long %>%
  filter(Gene %in% top_genes)

bar_plot <- ggplot(expr_long_top,
                   aes(x=Condition,y=Mean,fill=Condition)) +
  
  geom_bar(stat="identity") +
  
  geom_errorbar(
    aes(ymin=Mean-SE,ymax=Mean+SE),
    width=0.2
  ) +
  
  facet_wrap(~Gene,scales="free_y") +
  
  theme_minimal()

################################################
# 13) Heatmap top genes
################################################

expr_heatmap <- expr_long %>%
  select(Gene,Condition,Mean) %>%
  pivot_wider(
    names_from = Condition,
    values_from = Mean
  )

expr_heatmap <- as.data.frame(expr_heatmap)

rownames(expr_heatmap) <- expr_heatmap$Gene

expr_heatmap$Gene <- NULL

ordem_condicoes <- c(
  "Dry","X6h20","X6h25","X6h35",
  "RP20","RP25","RP35",
  "R2_20","R2_25","R2_35"
)

top_matrix <- expr_heatmap[top_genes,]

top_matrix <- top_matrix[,ordem_condicoes]

heatmap_top <- pheatmap(
  t(scale(t(top_matrix))),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  silent = TRUE
)

################################################
# 14) Save figure
################################################

pdf(file.path(outdir,"Analise_completa_genes_estresse.pdf"),
    width=14,height=10)

grid.arrange(
  heatmap_plot$gtable,
  heatmap_top$gtable,
  bar_plot,
  ncol=2
)

dev.off()