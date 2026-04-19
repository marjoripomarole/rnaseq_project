# scripts/02_deseq2.R
library(DESeq2)
library(ggplot2)
library(pheatmap)

# ── Load counts ──
raw <- read.delim("results/counts/gene_counts.txt", comment.char="#")
counts <- raw[, 7:ncol(raw)]
rownames(counts) <- raw$Geneid

# Clean column names
colnames(counts) <- gsub(".*aligned\\.|\\_Aligned.*", "", colnames(counts))

# ── Sample metadata ──
coldata <- data.frame(
  condition = factor(c("untreated","dex","untreated","dex",
                       "untreated","dex","untreated","dex")),
  row.names = colnames(counts)
)

# ── DESeq2 ──
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","dex","untreated"), alpha=0.05)
res <- res[order(res$padj), ]

# ── Summary ──
summary(res)
sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
cat("\nSignificant DE genes:", nrow(sig), "\n")
cat("Upregulated:", sum(sig$log2FoldChange > 0), "\n")
cat("Downregulated:", sum(sig$log2FoldChange < 0), "\n")

# ── Validate: check for known dex-responsive genes ──
known_genes <- c("CRISPLD2","DUSP1","KLF15","PER1","TSC22D3")
cat("\nValidation - known dex-responsive genes:\n")
for (g in known_genes) {
  idx <- grep(g, rownames(res))
  if (length(idx) > 0) {
    r <- res[idx[1], ]
    cat(sprintf("  %s: log2FC=%.2f, padj=%.2e %s\n",
        g, r$log2FoldChange, r$padj,
        ifelse(!is.na(r$padj) && r$padj < 0.05, "✓ SIGNIFICANT", "")))
  }
}

# ── Save results ──
write.csv(as.data.frame(res), "results/de/all_results.csv")
write.csv(as.data.frame(sig), "results/de/significant_genes.csv")
