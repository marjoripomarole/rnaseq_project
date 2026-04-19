# ── Volcano plot ──
df <- as.data.frame(res)
df$sig <- !is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 1

pdf("results/figures/volcano_plot.pdf", width=8, height=6)
ggplot(df, aes(log2FoldChange, -log10(padj), color=sig)) +
  geom_point(size=0.5, alpha=0.5) +
  scale_color_manual(values=c("grey70","red3"), labels=c("NS","Significant")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50") +
  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="grey50") +
  labs(title="Dexamethasone vs Untreated — Airway Smooth Muscle",
       subtitle=paste(nrow(sig), "significant DE genes (|log2FC| > 1, padj < 0.05)"),
       x="Log2 Fold Change", y="-Log10 Adjusted P-value", color="") +
  theme_minimal(base_size=12)
dev.off()

# ── PCA ──
vsd <- vst(dds)
pdf("results/figures/pca_plot.pdf", width=6, height=5)
plotPCA(vsd, intgroup="condition") +
  theme_minimal() +
  labs(title="PCA — Samples should cluster by condition")
dev.off()

# ── Heatmap of top 30 DE genes ──
top30 <- head(rownames(res[order(res$padj),]), 30)
mat <- assay(vsd)[top30, ]

pdf("results/figures/heatmap_top30.pdf", width=8, height=10)
pheatmap(mat, scale="row", annotation_col=coldata,
         color=colorRampPalette(c("navy","white","firebrick3"))(100),
         fontsize_row=8, main="Top 30 DE Genes")
dev.off()

cat("Figures saved to results/figures/\n")
