
# 0) Install / Load Required Packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
pkgs <- c(
  "GenomicRanges", "rtracklayer", "Rsamtools", "GenomicAlignments",
  "SummarizedExperiment", "DESeq2", "biomaRt", "S4Vectors",
  "ggplot2", "ggrepel"
)
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
  library(pkg, character.only = TRUE)
}

# ─────────────────────────────────────────────────────────────────────────────
# 1) Import master peaks (BED) + annotated peaks table
# ─────────────────────────────────────────────────────────────────────────────
master_bed_path <- "/c4/home/kraval/fix_AH/bam_files/macs/bed_tests/master_peaks.bed" #CHANGE WHERE YOUR master_peaks.bed file is located
annot_bed_path  <- "/c4/home/kraval/fix_AH/bam_files/macs/bed_tests/master_peaks.annotated.bed" #CHANGE WHERE THE ANNOTATED FILE IS (bedtools closest output)

master_peaks <- rtracklayer::import(master_bed_path)
annotated_peaks <- read.delim(
  annot_bed_path,
  header = FALSE,
  stringsAsFactors = FALSE
)
colnames(annotated_peaks) <- c(
  "chr", "start", "end",
  "match_chr", "match_start", "match_end",
  "geneSymbol", "distance"
)
annotated_peaks$peak_id <- paste0("Peak_", seq_len(nrow(annotated_peaks)))
if (is.null(names(master_peaks))) {
  names(master_peaks) <- paste0("Peak_", seq_along(master_peaks))
}

# ─────────────────────────────────────────────────────────────────────────────
# 2) Convert ENSG IDs in geneSymbol → HGNC using biomaRt (hg38)
# ─────────────────────────────────────────────────────────────────────────────
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
all_ids <- unique(annotated_peaks$geneSymbol)
ensembl_like <- grep("^ENSG\\d{11}", all_ids, value = TRUE)
if (length(ensembl_like) > 0) {
  conversion <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters    = "ensembl_gene_id",
    values     = ensembl_like,
    mart       = ensembl
  )
  colnames(conversion) <- c("geneSymbol_old", "geneSymbol_new")
  annotated_peaks <- merge(
    annotated_peaks, conversion,
    by.x = "geneSymbol", by.y = "geneSymbol_old",
    all.x = TRUE, sort = FALSE
  )
  annotated_peaks$geneSymbol <- ifelse(
    is.na(annotated_peaks$geneSymbol_new),
    annotated_peaks$geneSymbol,
    annotated_peaks$geneSymbol_new
  )
  annotated_peaks$geneSymbol_new <- NULL
}

# ─────────────────────────────────────────────────────────────────────────────
# 3) Force-label any peak overlapping the TERT locus as “TERT”
#    (hg38 Chr 5:1,259,000–1,298,000)
# ─────────────────────────────────────────────────────────────────────────────
tertlocus <- GRanges("5", IRanges(start = 1259000, end = 1298000))
hits_tert_all <- subsetByOverlaps(master_peaks, tertlocus)
if (length(hits_tert_all) > 0) {
  tert_ids_all <- names(hits_tert_all)
  annotated_peaks$geneSymbol[
    annotated_peaks$peak_id %in% tert_ids_all
  ] <- "TERT"
  message("Forced TERT labeling on peaks: ", paste(tert_ids_all, collapse = ", "))
} else {
  message("No peaks overlap the TERT locus.")
}
gene_map <- setNames(annotated_peaks$geneSymbol, annotated_peaks$peak_id)

# ─────────────────────────────────────────────────────────────────────────────
# 4) Define & load all eight BAM files into a named BamFileList
# ─────────────────────────────────────────────────────────────────────────────
bam_files <- c(
  "/c4/home/kraval/fix_AH/bam_files/SF7996_Coff_sgScrambled_H3K9me3_1_S194_L008.rmdup.bam",
  "/c4/home/kraval/fix_AH/bam_files/SF7996_Coff_sgScrambled_H3K9me3_2_S195_L008.rmdup.bam",
  "/c4/home/kraval/fix_AH/bam_files/SF7996_Coff_sgTERT_pooledCi_cloneA2_H3K9me3_1_S198_L008.rmdup.bam",
  "/c4/home/kraval/fix_AH/bam_files/SF7996_Coff_sgTERT_pooledCi_cloneA2_H3K9me3_2_S199_L008.rmdup.bam",
  "/c4/home/kraval/fix_AH/bam_files/SF7996_Coff_sgTERT_pooledCi_cloneA4_H3K9me3_1_S200_L008.rmdup.bam",
  "/c4/home/kraval/fix_AH/bam_files/SF7996_Coff_sgTERT_pooledCi_cloneA4_H3K9me3_2_S201_L008.rmdup.bam",
  "/c4/home/kraval/fix_AH/bam_files/SF7996_Coff_sgTERT_pooledCi_polycolonal_H3K9me3_1_S196_L008.rmdup.bam",
  "/c4/home/kraval/fix_AH/bam_files/SF7996_Coff_sgTERT_pooledCi_polycolonal_H3K9me3_2_S197_L008.rmdup.bam"
)
bam_names <- c(
  "sgScramb_1", "sgScramb_2",
  "A2_1",       "A2_2",
  "A4_1",       "A4_2",
  "Poly_1",     "Poly_2"
)
bam_list <- BamFileList(bam_files)
names(bam_list) <- bam_names

# ─────────────────────────────────────────────────────────────────────────────
# 5) Count overlaps: build SummarizedExperiment 'se'
# ─────────────────────────────────────────────────────────────────────────────
se <- summarizeOverlaps(
  features      = master_peaks,
  reads         = bam_list,
  mode          = "Union",
  singleEnd     = FALSE,
  ignore.strand = TRUE,
  fragments     = TRUE
)
colnames(se) <- bam_names

# ─────────────────────────────────────────────────────────────────────────────
# 6) Filter out peaks with total reads < 50 (across all 8 samples)
# ─────────────────────────────────────────────────────────────────────────────
tot_counts_pre  <- rowSums(assay(se))
keep_mask       <- tot_counts_pre >= 50
se_filt         <- se[keep_mask, ]
message("Number of peaks before filtering: ", nrow(se))
message("Number of peaks after filtering (≥50 reads): ", nrow(se_filt))
tot_counts_post <- rowSums(assay(se_filt))
if (!all(tot_counts_post >= 50)) {
  bad <- which(tot_counts_post < 50)
  warning("Some peaks in se_filt still <50 reads: ",
          paste(rownames(se_filt)[bad], collapse = ", "))
} else {
  message("All ", nrow(se_filt), " peaks in se_filt have ≥50 reads.")
}

# ─────────────────────────────────────────────────────────────────────────────
# 7) Build DESeq2 dataset on 'se_filt' and run DESeq()
# ─────────────────────────────────────────────────────────────────────────────
sample_info <- DataFrame(
  condition = factor(c(
    "sgScramb","sgScramb",
    "A2","A2",
    "A4","A4",
    "Poly","Poly"
  )),
  row.names  = colnames(se_filt)
)
colData(se_filt) <- sample_info
dds <- DESeqDataSet(se_filt, design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "sgScramb")
dds <- DESeq(dds)

# ─────────────────────────────────────────────────────────────────────────────
# 8) Extract DESeq2 results for four contrasts
# ─────────────────────────────────────────────────────────────────────────────
res_A2   <- results(dds, contrast = c("condition","A2","sgScramb"))
res_A4   <- results(dds, contrast = c("condition","A4","sgScramb"))
res_Poly <- results(dds, contrast = c("condition","Poly","sgScramb"))

# Combined A2+A4: average log2FC, take max p-value, then BH correction
log2_A2   <- res_A2$log2FoldChange
log2_A4   <- res_A4$log2FoldChange
log2_comb <- rowMeans(cbind(log2_A2, log2_A4), na.rm = FALSE)
p_comb    <- pmax(res_A2$pvalue, res_A4$pvalue, na.rm = TRUE)
res_A2A4           <- res_A2
res_A2A4$log2FoldChange <- log2_comb
res_A2A4$pvalue         <- p_comb
res_A2A4$padj           <- p.adjust(p_comb, method = "BH")

# ─────────────────────────────────────────────────────────────────────────────
# 9) Attach 'gene' column to each result (so TERT appears)
# ─────────────────────────────────────────────────────────────────────────────
res_A2$gene   <- gene_map[rownames(res_A2)]
res_A2$gene[is.na(res_A2$gene)] <- ""
res_A4$gene   <- gene_map[rownames(res_A4)]
res_A4$gene[is.na(res_A4$gene)] <- ""
res_Poly$gene <- gene_map[rownames(res_Poly)]
res_Poly$gene[is.na(res_Poly$gene)] <- ""
res_A2A4$gene <- gene_map[rownames(res_A2A4)]
res_A2A4$gene[is.na(res_A2A4$gene)] <- ""

# ─────────────────────────────────────────────────────────────────────────────
# 10) Write out CSVs for each contrast (filtered)
#       Columns: peak_id, gene, log2FC, pvalue, negLog10_pval, sig_label, total_counts
# ─────────────────────────────────────────────────────────────────────────────
make_volcano_csv_filtered <- function(res, se_obj, output_csv) {
  df <- as.data.frame(res)
  df$peak_id       <- rownames(df)
  df$negLog10_pval <- -log10(df$pvalue + 1e-8)
  df$sig_label     <- (df$pvalue < 0.01 & abs(df$log2FoldChange) > 1)
  if (!"gene" %in% colnames(df)) {
    df$gene <- ""
  } else {
    df$gene[is.na(df$gene)] <- ""
  }
  total_counts_vec <- rowSums(assay(se_obj))
  df$total_counts <- total_counts_vec[df$peak_id]
  df_out <- df[ , c(
    "peak_id", "gene", "log2FoldChange", "pvalue", 
    "negLog10_pval", "sig_label", "total_counts"
  ) ]
  write.csv(df_out, file = output_csv, row.names = FALSE)
  message("Wrote filtered CSV: ", output_csv, "  [rows = ", nrow(df_out), "]")
  invisible(df_out)
}

csv_A2   <- make_volcano_csv_filtered(res_A2,   se_filt, "volcano_A2_filtered.csv")
csv_A4   <- make_volcano_csv_filtered(res_A4,   se_filt, "volcano_A4_filtered.csv")
csv_Poly <- make_volcano_csv_filtered(res_Poly, se_filt, "volcano_Poly_filtered.csv")
csv_A2A4 <- make_volcano_csv_filtered(res_A2A4, se_filt, "volcano_A2A4_filtered.csv")

# ─────────────────────────────────────────────────────────────────────────────
# 11) Define volcano-plot function that always labels TERT (blue) + sig peaks (red)
# ─────────────────────────────────────────────────────────────────────────────
plot_filtered_with_TERT <- function(res, title_text, file_name = NULL) {
  df <- as.data.frame(res)
  df$peak_id       <- rownames(df)
  df$negLog10_pval <- -log10(df$pvalue + 1e-8)
  df$label_std     <- ifelse(df$pvalue < 0.01 & abs(df$log2FoldChange) > 1, df$gene, "")
  df$label_TERT    <- ifelse(df$gene == "TERT", df$gene, "")
  df$color_group   <- ifelse(
    df$gene == "TERT",                                     "TERT",
    ifelse(df$pvalue < 0.01 & abs(df$log2FoldChange) > 1,  "sig", "ns")
  )
  df$color_group <- factor(df$color_group, levels = c("ns","sig","TERT"))
  p <- ggplot(df, aes(x = log2FoldChange, y = negLog10_pval)) +
    geom_point(aes(color = color_group), alpha = 0.6, size = 1.8) +
    scale_color_manual(
      values = c("ns" = "grey70", "sig" = "red", "TERT" = "blue")
    ) +
    geom_text_repel(
      data = subset(df, label_std != ""),
      aes(label = label_std),
      color        = "black",
      size         = 3,
      max.overlaps = 20
    ) +
    geom_text_repel(
      data = subset(df, label_TERT != ""),
      aes(label = label_TERT),
      color        = "blue",
      size         = 3,
      max.overlaps = 20
    ) +
    theme_minimal() +
    labs(
      title = paste0(title_text, " (filtered ≥50 reads)"),
      x     = "log₂ Fold Change",
      y     = "-log₁₀(p-value)"
    ) +
    theme(
      legend.position = "none",
      plot.title     = element_text(hjust = 0.5)
    )
  if (!is.null(file_name)) {
    ggsave(file_name, plot = p, width = 7, height = 6)
  }
  return(p)
}

# ─────────────────────────────────────────────────────────────────────────────
# 12) Generate and save volcano plots for each contrast
# ─────────────────────────────────────────────────────────────────────────────
pA2_filtered_TERT   <- plot_filtered_with_TERT(
  res_A2,   "A2 vs sgScramb",        "volcano_A2_filtered_with_TERT.png"
)
pA4_filtered_TERT   <- plot_filtered_with_TERT(
  res_A4,   "A4 vs sgScramb",        "volcano_A4_filtered_with_TERT.png"
)
pPoly_filtered_TERT <- plot_filtered_with_TERT(
  res_Poly, "Polyclonal vs sgScramb","volcano_Poly_filtered_with_TERT.png"
)
pA2A4_filtered_TERT <- plot_filtered_with_TERT(
  res_A2A4, "A2+A4 vs sgScramb",     "volcano_A2A4_filtered_with_TERT.png"
)

# (Optional) Display all four plots in your R plotting window
print(pA2_filtered_TERT)
print(pA4_filtered_TERT)
print(pPoly_filtered_TERT)
print(pA2A4_filtered_TERT)

message("Script complete: CSVs written and volcano plots saved.")
