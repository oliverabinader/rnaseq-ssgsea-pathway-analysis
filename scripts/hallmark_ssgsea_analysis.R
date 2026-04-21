# ===============================
# Hallmark ssGSEA Pathway Analysis
# ===============================

library(dplyr)
library(tibble)
library(msigdbr)
library(GSVA)
library(pheatmap)
library(openxlsx)

# -------------------------------
# Input files (EDIT THESE PATHS)
# -------------------------------
tpm_file <- "/path/to/non-logged-scaled-TPM.allsamples.tsv"
metadata_file <- "/path/to/metadata.tsv"

# -------------------------------
# Load expression data
# -------------------------------
TPM <- read.delim(tpm_file, sep = ",") %>%
  column_to_rownames("GeneSymbol") %>%
  select(-Geneid)

# Log transform
TPM <- log(TPM + 1)
TPM <- as.matrix(TPM)

# -------------------------------
# Load gene sets (Hallmark)
# -------------------------------
gs.MH <- msigdbr(species = "Mus musculus", category = "H") %>%
  select(gs_name, gene_symbol) %>%
  group_by(gs_name) %>%
  summarise(all.genes = list(unique(gene_symbol))) %>%
  deframe()

# -------------------------------
# Run ssGSEA
# -------------------------------
ssgsea_param <- ssgseaParam(TPM, gs.MH)
ssgsea_scores <- gsva(ssgsea_param)

# Clean pathway names
rownames(ssgsea_scores) <- gsub("^HALLMARK_", "", rownames(ssgsea_scores))

# Scale scores (Z-score)
ssgsea_scaled <- t(scale(t(ssgsea_scores)))
#rownames(ssgsea_scaled) <- gsub("^HALLMARK_", "", rownames(ssgsea_scaled))

# -------------------------------
# Load metadata
# -------------------------------
metadata <- read.delim(metadata_file)

ann_col <- metadata %>%
  select(sample, condition1, condition2) %>%
  column_to_rownames("sample") %>%
  #mutate(condition2 = gsub("_", " ", condition2)) %>%
  rename(Group = condition1, Model = condition2)

# -------------------------------
# Plot heatmap
# -------------------------------
p <- pheatmap(
  ssgsea_scaled,
  annotation_col = ann_col,
  main = "Hallmark ssGSEA",
  show_colnames = TRUE
)

# -------------------------------
# Save outputs
# -------------------------------
dir.create("/path/to/results", showWarnings = FALSE)

png("results/Hallmark_ssGSEA.png", width = 1200, height = 1000, res = 150)
grid::grid.newpage()
grid::grid.draw(p$gtable)
dev.off()

wb <- createWorkbook()
addWorksheet(wb, "ssGSEA_Raw")
addWorksheet(wb, "ssGSEA_Scaled")

writeDataTable(wb, 1, as.data.frame(ssgsea_scores), rowNames = TRUE)
writeDataTable(wb, 2, as.data.frame(ssgsea_scaled), rowNames = TRUE)

saveWorkbook(wb, "results/Hallmark_ssGSEA.xlsx", overwrite = TRUE)

cat("Analysis complete.\n")
