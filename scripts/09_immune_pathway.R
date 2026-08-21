#!/usr/bin/env Rscript
# 09_immune_pathway.R —— Hallmark 通路评分 + 免疫标记相关性（真实实现）
# 基因集来源：msigdbr 7.5.1（本地内置 MSigDB Hallmark，无需联网）。28 免疫细胞反卷积需
#   Charoentong 基因集（未随数据提供），本步用经典免疫标记基因做相关性替代（见 deviation_log 待办）。

suppressPackageStartupMessages({
  library(optparse)
  library(GSVA)
  library(msigdbr)
  library(ggpubr)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "09_immune_pathway.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--expr"), type = "character",
              default = file.path(DATA_PROC, "expr_human_main.rds"), help = "GSE115857 表达矩阵"),
  make_option(c("--meta"), type = "character",
              default = file.path(DATA_PROC, "meta_human_main.csv"), help = "GSE115857 表型"),
  make_option(c("--subtype"), type = "character",
              default = file.path(DATA_PROC, "subtype_assignment.csv"), help = "08 分型（仅 CKD）"),
  make_option(c("--sig"), type = "character",
              default = file.path(DATA_PROC, "progression_signature.rds"), help = "06 signature"),
  make_option(c("--outdir"), type = "character", default = DATA_PROC, help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

main <- function() {
  ensure_dirs()
  expr <- readRDS(opt$expr)
  meta <- read.csv(opt$meta, stringsAsFactors = FALSE)
  sig <- readRDS(opt$sig)
  hall <- msigdbr(species = "Homo sapiens", category = "H")
  want <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
            "HALLMARK_TGF_BETA_SIGNALING",
            "HALLMARK_INFLAMMATORY_RESPONSE",
            "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
            "HALLMARK_IL6_JAK_STAT3_SIGNALING",
            "HALLMARK_HYPOXIA",
            "HALLMARK_APICAL_JUNCTION",
            "HALLMARK_COMPLEMENT")
  sets <- split(hall$gene_symbol[hall$gs_name %in% want], hall$gs_name[hall$gs_name %in% want])
  score <- gsva(ssgseaParam(exprData = as.matrix(expr), geneSets = sets))
  safe_write(file.path(opt$outdir, "hallmark_scores.csv"), function(p) {
    out <- data.frame(sample = colnames(score), t(score), check.names = FALSE)
    write_csv_utf8(out, p)
  })

  meta$disease <- ifelse(meta$status == "IgAN patient", "IgAN", "control")
  ckd <- meta$disease == "IgAN"
  emt <- as.numeric(score["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", ])
  names(emt) <- colnames(score)
  tgf <- as.numeric(score["HALLMARK_TGF_BETA_SIGNALING", ])
  names(tgf) <- colnames(score)
  log_msg("pathway", paste0("EMT IgAN vs control p=",
                            signif(wilcox.test(emt[ckd], emt[!ckd])$p.value, 3),
                            " | TGFB p=", signif(wilcox.test(tgf[ckd], tgf[!ckd])$p.value, 3)))

  sub <- read.csv(opt$subtype, stringsAsFactors = FALSE)
  common <- intersect(sub$sample, colnames(expr))
  df <- data.frame(subtype = factor(sub$subtype[match(common, sub$sample)]),
                   EMT = emt[common])
  p <- ggboxplot(df, x = "subtype", y = "EMT", color = "subtype", palette = CB_PALETTE,
                 add = "jitter", xlab = NULL, ylab = "EMT pathway score") +
    stat_compare_means(method = "wilcox.test") + pub_theme()
  save_pub_fig(p, file.path(FIG_DIR, "fig_subtype_emt"))

  markers <- intersect(c("PTPRC", "CD8A", "CD4", "FOXP3", "CD68", "ITGAX",
                         "NKG7", "CD19", "MS4A1", "LYZ", "CD163"), rownames(expr))
  cor_mat <- cor(t(expr[sig, ckd, drop = FALSE]), t(expr[markers, ckd, drop = FALSE]),
                 method = "spearman")
  safe_write(file.path(opt$outdir, "sig_immune_marker_cor.csv"), function(p) {
    write_csv_utf8(as.data.frame(cor_mat), p)
  })
  hp <- file.path(FIG_DIR, "fig_sig_immune_heatmap.pdf")
  hp_png <- file.path(FIG_DIR, "fig_sig_immune_heatmap.png")
  if (file.exists(hp)) unlink(hp)
  if (file.exists(hp_png)) unlink(hp_png)
  pdf(hp, width = 5, height = 4)
    pheatmap(cor_mat, color = colorRampPalette(c("#0072B2", "white", "#E69F00"))(50),
             fontsize = 6, cluster_cols = FALSE)
  dev.off()
  png(hp_png, width = 5, height = 4, units = "in", res = 300)
  pheatmap(cor_mat, color = colorRampPalette(c("#0072B2", "white", "#E69F00"))(50),
           fontsize = 6, cluster_cols = FALSE)
  dev.off()
  log_msg("immune", "DONE")
}

if (!interactive()) main()
