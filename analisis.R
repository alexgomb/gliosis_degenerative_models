set.seed(42)
dir.create("resultados", showWarnings = FALSE)

library(GEOquery)
library(Matrix)
library(Seurat)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(EnhancedVolcano)
library(patchwork)

# Define a vector with sample IDs
gsm_ids <- c(
  "GSM8728964", 
  "GSM8728965", 
  "GSM8728966", 
  "GSM8728967", 
  "GSM8728968", 
  "GSM3854512"  
)

# Download metadata for each GSM without quantitative expression series
gsm_list <- lapply(gsm_ids, function(id) {
    getGEO(id, GSEMatrix = FALSE)
})

# Create a folder for each GSM and download supplementary files
for (gsm in gsm_ids) {
  dir.create(paste0(gsm, "_data"), showWarnings = FALSE)
  getGEOSuppFiles(gsm, makeDirectory = FALSE, baseDir = paste0(gsm, "_data"))
}

# Verify contents of all folders
for (gsm in gsm_ids) {
  cat("Files in", gsm, "_data:\n")
  print(list.files(paste0(gsm, "_data")))
  cat("\n")
}

# Define IDs and folders containing files for each GSM
gsm_ids <- c(
  "GSM8728964", 
  "GSM8728965", 
  "GSM8728966", 
  "GSM8728967", 
  "GSM8728968", 
  "GSM3854512"
)

# Construct folder names for each GSM
carpetas <- paste0(gsm_ids, "_data")

# Prepare an empty list to store matrices
lista_expr <- vector("list", length(gsm_ids))
names(lista_expr) <- gsm_ids

# Loop to read each sample
for (i in seq_along(gsm_ids)) {
  gsm   <- gsm_ids[i]
  dir   <- carpetas[i]
  
  # Identify filenames within the folder
  ficheros <- list.files(dir, full.names = TRUE)
  
  # Read barcodes 
  fn_barcodes <- ficheros[grep("_barcodes.tsv.gz$", ficheros)]
  barcodes    <- read.delim(fn_barcodes, header = FALSE, stringsAsFactors = FALSE)
  
  # Read features or genes 
  fn_features <- ficheros[grep("_features.tsv.gz$|_genes.tsv.gz$", ficheros)]
  features    <- read.delim(fn_features, header = FALSE, stringsAsFactors = FALSE)
  
  # Read the sparse matrix
  fn_matrix <- ficheros[grep("_matrix.mtx.gz$", ficheros)]
  expr_mat  <- readMM(fn_matrix)
  
  # Assign row and column names
  rownames(expr_mat) <- features$V2
  colnames(expr_mat) <- barcodes$V1
  
  # Save the matrix in the list
  lista_expr[[gsm]] <- expr_mat
  
  # Show a brief summary of dimensions
  cat("->", gsm, ": number of genes =", nrow(expr_mat),
      ", number of cells =", ncol(expr_mat), "\n")
  
  # Free temporary variables
  rm(expr_mat, barcodes, features)
  gc()
}

# Make gene names unique for each sample
for (gsm in names(lista_expr)) {
  genes_orig      <- rownames(lista_expr[[gsm]])
  genes_unicos    <- make.unique(genes_orig)
  rownames(lista_expr[[gsm]]) <- genes_unicos
  
  rm(genes_orig, genes_unicos)
}
gc()

# Verify duplicated genes in each sample
for (gsm in names(lista_expr)) {
  n_dup <- sum(duplicated(rownames(lista_expr[[gsm]])))
  cat("->", gsm, ": duplicated genes =", n_dup, "\n")
  
  rm(n_dup)
}
gc()

# Initialize an empty list to store Seurat objects
lista_seurat <- vector("list", length(lista_expr))
names(lista_seurat) <- names(lista_expr)

# Loop to create each Seurat object
for (gsm in names(lista_expr)) {
  lista_seurat[[gsm]] <- CreateSeuratObject(
    counts  = lista_expr[[gsm]],
    project = gsm
  )
  cat("-> Created Seurat for", gsm, ":",
      nrow(lista_seurat[[gsm]]), "genes x",
      ncol(lista_seurat[[gsm]]), "cells\n")
}

# Remove the original list of matrices and force garbage collection to manage RAM limit
rm(lista_expr)
gc()

for (gsm in names(lista_seurat)) {
  seu <- lista_seurat[[gsm]]
  
  # Calculate percent.mt if it does not exist
  if (!"percent.mt" %in% colnames(seu@meta.data)) {
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
    lista_seurat[[gsm]] <- seu  # save the updated object
  }
  
  # Print the summary for that sample
  cat(">>>", gsm, "<<<\n")
  cat("Summary of nFeature_RNA:\n")
  print(summary(seu@meta.data$nFeature_RNA))
  cat("Summary of percent.mt:\n")
  print(summary(seu@meta.data$percent.mt))
  cat("\n-------------------------\n\n")
  
  rm(seu)
}
gc()

# Apply filtering to each Seurat object and display final cell count
for (gsm in names(lista_seurat)) {
  seu <- lista_seurat[[gsm]]
  
  # Filter objects based on transcripts and mitochondrial RNA fraction
  seu_filtrado <- subset(
    seu,
    subset = nFeature_RNA > 200 & percent.mt < 10
  )
  
  # Save the filtered Seurat object in the list
  lista_seurat[[gsm]] <- seu_filtrado
  
  # Print the number of remaining cells
  cat("->", gsm, ": cells after QC =", ncol(seu_filtrado), "\n")
  
  rm(seu, seu_filtrado)
}
gc()

# Extract the vector of present genes in each filtered sample using GetAssayData
lista_genes_por_muestra <- lapply(lista_seurat, function(s) {
  counts_mat <- GetAssayData(s, layer = "counts")
  rownames(counts_mat)
})

# Calculate the intersection to extract common genes
genes_comunes_QC <- Reduce(intersect, lista_genes_por_muestra)

# Print the number of common genes
length(genes_comunes_QC)

rm(lista_genes_por_muestra)
gc()

# Create empty matrix: rows = common genes, columns = each sample
pseudo_bulk_mat_QC <- matrix(
  0,
  nrow = length(genes_comunes_QC),
  ncol = length(lista_seurat),
  dimnames = list(genes_comunes_QC, names(lista_seurat))
)

# Traverse each sample and aggregate gene transcript counts
for (gsm in names(lista_seurat)) {
  # Extract the counts matrix from the Seurat object
  counts_qc <- GetAssayData(lista_seurat[[gsm]], layer = "counts")
  
  # Select only the common genes
  submat    <- counts_qc[genes_comunes_QC, , drop = FALSE]
  
  # Compute row sums and register in the corresponding pseudo-bulk column
  pseudo_bulk_mat_QC[, gsm] <- Matrix::rowSums(submat)
  
  rm(counts_qc, submat)
}
gc()

# Convert the structural matrix format to data.frame
pseudo_bulk_df_QC <- as.data.frame(pseudo_bulk_mat_QC)

# Print dimensions and dataframe sample
cat("Dimensions of pseudo-bulk after QC:", dim(pseudo_bulk_df_QC), "\n\n")
head(pseudo_bulk_df_QC)

# Discard large redundant objects
rm(pseudo_bulk_mat_QC, lista_seurat)
gc()

# Install and load DESeq2 if necessary
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

# Create a data.frame with the structural sample information
colData <- data.frame(
  sample     = colnames(pseudo_bulk_df_QC),
  condition  = c(
    "rd1", "rd1", "rd1",   # GSM8728964, 8965, 8966 correspond to rd1
    "rd10","rd10",         # GSM8728967, 8968 correspond to rd10
    "wt"                  # GSM3854512 corresponds to C57 WT
  ),
  row.names = colnames(pseudo_bulk_df_QC),
  stringsAsFactors = FALSE
)
colData

# Construct the colData data.frame properly formatted
colData <- data.frame(
  condition = c(
    "rd1", "rd1", "rd1",   
    "rd10","rd10",         
    "wt"                   
  ),
  row.names = colnames(pseudo_bulk_df_QC),
  stringsAsFactors = FALSE
)

# Verify dimensional parameters and content
dim(colData)
colData

# Create the DESeqDataSet architecture utilizing the count matrix and colData
dds <- DESeqDataSetFromMatrix(
  countData = pseudo_bulk_df_QC,
  colData   = colData,
  design    = ~ condition
)

# Render class overview
dds

# Execute the core differential expression analysis routine
dds <- DESeq(dds)

# Extract differential outputs across contrast architectures
res_rd1_vs_wt   <- results(dds, contrast = c("condition", "rd1",  "wt"))
res_rd10_vs_wt  <- results(dds, contrast = c("condition", "rd10", "wt"))
res_rd1_vs_rd10 <- results(dds, contrast = c("condition", "rd1",  "rd10"))

# Auxiliary function calculating Total, Up, Down, Base ratios 
resumen_contraste_pct <- function(res_obj) {
  # Coerce to structural data.frame
  df <- as.data.frame(res_obj)
  # Filter underlying genes possessing mean counts >= 6
  df_filt <- df[!is.na(df$padj) & df$baseMean >= 6, ]
  
  total   <- nrow(df_filt)
  up      <- sum(df_filt$padj < 0.1 & df_filt$log2FoldChange >  0, na.rm = TRUE)
  down    <- sum(df_filt$padj < 0.1 & df_filt$log2FoldChange <  0, na.rm = TRUE)
  pct_up   <- up    / total * 100
  pct_down <- down  / total * 100
  
  return(c(
    Total   = total,
    Up      = up,
    Down    = down,
    PctUp   = round(pct_up, 2),
    PctDown = round(pct_down, 2)
  ))
}

# Consolidate output ratio framework
tabla_resumen_pct <- rbind(
  `rd1_vs_wt`   = resumen_contraste_pct(res_rd1_vs_wt),
  `rd10_vs_wt`  = resumen_contraste_pct(res_rd10_vs_wt),
  `rd1_vs_rd10` = resumen_contraste_pct(res_rd1_vs_rd10)
)

# Output summary layout
tabla_resumen_pct

# Process adaptive shrinkage for rd1 vs wt 
resLFC_rd1_vs_wt <- lfcShrink(
  dds,
  contrast = c("condition", "rd1", "wt"),
  type     = "ashr"
)

# Process adaptive shrinkage for rd10 vs wt
resLFC_rd10_vs_wt <- lfcShrink(
  dds,
  contrast = c("condition", "rd10", "wt"),
  type     = "ashr"
)

# Process adaptive shrinkage for rd1 vs rd10 
resLFC_rd1_vs_rd10 <- lfcShrink(
  dds,
  contrast = c("condition", "rd1", "rd10"),
  type     = "ashr"
)

# Declare comprehensive gliosis markers
genes_gliosis <- c(
  "Gfap", "Vim", "S100b", "Aldh1l1",
  "Rlbp1", "Glul", "Slc1a3",
  "C1qa", "C1qb", "Trem2", "Aif1"
)

# Coerce derived structural data elements
tablaLFC_rd1  <- as.data.frame(resLFC_rd1_vs_wt)
tablaLFC_rd10 <- as.data.frame(resLFC_rd10_vs_wt)

# Filter markers parallel defined profiles
tabla_gliosis_rd1  <- tablaLFC_rd1[rownames(tablaLFC_rd1) %in% genes_gliosis, ]
tabla_gliosis_rd10 <- tablaLFC_rd10[rownames(tablaLFC_rd10) %in% genes_gliosis, ]

# Output resultant data structures
head(tabla_gliosis_rd1)
head(tabla_gliosis_rd10)

# Coerce filtered responses
df_rd1_vs_wt  <- as.data.frame(resLFC_rd1_vs_wt)
df_rd10_vs_wt <- as.data.frame(resLFC_rd10_vs_wt)

# Filter parameters enforcing padj < 0.1 criteria coupled with proportional baseline increases
genes_up_rd1  <- rownames(df_rd1_vs_wt)[df_rd1_vs_wt$padj < 0.1 & df_rd1_vs_wt$log2FoldChange >  0]
genes_up_rd10 <- rownames(df_rd10_vs_wt)[df_rd10_vs_wt$padj < 0.1 & df_rd10_vs_wt$log2FoldChange >  0]

entrez_rd1 <- bitr(
  genes_up_rd1,
  fromType   = "SYMBOL",
  toType     = "ENTREZID",
  OrgDb      = org.Mm.eg.db
)$ENTREZID

entrez_rd10 <- bitr(
  genes_up_rd10,
  fromType   = "SYMBOL",
  toType     = "ENTREZID",
  OrgDb      = org.Mm.eg.db
)$ENTREZID

# Perform GO Biological Process statistical derivation representing rd1
ego_rd1 <- enrichGO(
  gene          = entrez_rd1,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE   # Transpose internal ENTREZID configurations to standard notation
)

# Perform GO derivation representing rd10
ego_rd10 <- enrichGO(
  gene          = entrez_rd10,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Subset high priority integrated routes
top10_rd1  <- as.data.frame(ego_rd1)[1:10,  ]
top10_rd10 <- as.data.frame(ego_rd10)[1:10, ]

# Output derivations
cat("=== Top 10 GO BP in rd1 vs WT ===\n")
print(top10_rd1[, c("ID", "Description", "GeneRatio", "p.adjust")])
cat("\n=== Top 10 GO BP in rd10 vs WT ===\n")
print(top10_rd10[, c("ID", "Description", "GeneRatio", "p.adjust")])

# Map explicit representations rd1 vs WT
df_volc_rd1 <- as.data.frame(resLFC_rd1_vs_wt)
df_volc_rd1$gene <- rownames(df_volc_rd1)

# Map explicit representations rd10 vs WT
df_volc_rd10 <- as.data.frame(resLFC_rd10_vs_wt)
df_volc_rd10$gene <- rownames(df_volc_rd10)

p_rd1 <- EnhancedVolcano(
  df_volc_rd1,
  lab            = df_volc_rd1$gene,
  x              = "log2FoldChange",
  y              = "padj",
  xlab           = bquote(~Log[2]~ "(Fold Change)"),
  ylab           = bquote(~-Log[10]~ "(adj P-value)"),
  pCutoff        = 0.05,
  FCcutoff       = 1,
  pointSize      = 2.0,
  labSize        = 3.0,
  title          = "rd1 vs WT",
  subtitle       = "padj < 0.05, |log2FC| > 1",
  caption        = "Upregulated targets represented in red",
  drawConnectors = TRUE,
  max.overlaps   = 30,
  legendPosition = "bottom"
)

p_rd10 <- EnhancedVolcano(
  df_volc_rd10,
  lab            = df_volc_rd10$gene,
  x              = "log2FoldChange",
  y              = "padj",
  xlab           = bquote(~Log[2]~ "(Fold Change)"),
  ylab           = bquote(~-Log[10]~ "(adj P-value)"),
  pCutoff        = 0.05,
  FCcutoff       = 1,
  pointSize      = 2.0,
  labSize        = 3.0,
  title          = "rd10 vs WT",
  subtitle       = "padj < 0.05, |log2FC| > 1",
  caption        = "Upregulated targets represented in red",
  drawConnectors = TRUE,
  max.overlaps   = 30,
  legendPosition = "bottom"
)

volcano_combinado <- p_rd1 + p_rd10 + plot_layout(ncol = 2)
print(volcano_combinado)

# Export structural datasets in CSV format
write.csv(tabla_resumen_pct, file = "resultados/tabla_resumen_pct.csv")
write.csv(tabla_gliosis_rd1, file = "resultados/tabla_gliosis_rd1.csv")
write.csv(tabla_gliosis_rd10, file = "resultados/tabla_gliosis_rd10.csv")
write.csv(top10_rd1, file = "resultados/top10_rd1.csv")
write.csv(top10_rd10, file = "resultados/top10_rd10.csv")

# Export consolidated high-resolution volcano visualization
ggsave("resultados/volcano_plot_combinado.png", plot = volcano_combinado, width = 14, height = 7, dpi = 300)
