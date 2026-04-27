# ============================================================
# Bulk RNA-seq Analysis (Tau ASO iPSC neurons)
# ============================================================
# - Merge RSEM outputs
# - Run DESeq2
# - Annotate genes (Ensembl → Symbol/Entrez)
# ============================================================

setwd("rsem_output")

# ============================
# Load libraries
# ============================
library(data.table)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ============================
# 1. Build count matrix
# ============================
files <- list.files(pattern = "SRR.*\\.genes.results$", full.names = TRUE)

counts <- fread(files[1], select = c("gene_id", "expected_count"))
colnames(counts)[2] <- gsub(".genes.results", "", basename(files[1]))

for (f in files[-1]) {
  tmp <- fread(f, select = c("gene_id", "expected_count"))
  colnames(tmp)[2] <- gsub(".genes.results", "", basename(f))
  counts <- merge(counts, tmp, by = "gene_id")
}

fwrite(counts, "RSEM_gene_count_matrix.tsv", sep = "\t")

# ============================
# 2. Prepare DESeq2 input
# ============================
counts <- as.data.frame(counts)
rownames(counts) <- counts$gene_id
counts$gene_id <- NULL

# Metadata
sample_info <- data.frame(
  row.names = colnames(counts),
  condition = c(
    "Tau.ASO", "Tau.ASO", "Tau.ASO",
    "ASO.NT", "ASO.NT", "ASO.NT"
  )
)

# Convert to integer counts
counts_int <- round(counts)

# ============================
# 3. Run DESeq2
# ============================
dds <- DESeqDataSetFromMatrix(
  countData = counts_int,
  colData = sample_info,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "Tau.ASO", "ASO.NT"))
res <- res[order(res$padj), ]

write.csv(as.data.frame(res), "DESeq2_ASO_vs_Control.csv")

# ============================
# 4. QC Plots
# ============================
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")

plotMA(res, ylim = c(-5, 5))

# ============================
# 5. Gene annotation
# ============================
res_annot <- as.data.frame(res)
res_annot$gene_id <- rownames(res_annot)

res_annot$symbol <- mapIds(
  org.Hs.eg.db,
  keys = res_annot$gene_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

res_annot$entrez <- mapIds(
  org.Hs.eg.db,
  keys = res_annot$gene_id,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

write.csv(
  res_annot,
  file = "DESeq2_results_annotated.csv",
  row.names = FALSE
)

cat("Analysis complete.\n")
