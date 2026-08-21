#!/usr/bin/env Rscript
# 15_immune_deconvolution.R —— 免疫细胞反卷积（Bindea 24 签名，ssGSEA）
# 来源说明（D20）：Charoentong 2017 的 28 细胞基因集（Table S2）在本网络不可达
# （Cell Press/ScienceDirect 403）；采用同框架的 **Bindea et al. Immunity 2013 24 免疫签名**，
# 基因集取自 IOBR 包 `signature_collection`（data/reference/signature_collection_IOBR.rda）。
# 输出：IgAN vs 对照各细胞评分差异（wilcox+BH）、C1/C2 差异、签名基因-免疫评分相关。

suppressPackageStartupMessages({
  library(optparse)
  library(GSVA)
  library(ggpubr)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "15_immune_deconvolution.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--expr"), type = "character",
              default = file.path(DATA_PROC, "expr_human_main.rds"), help = "主队列表达矩阵"),
  make_option(c("--meta"), type = "character",
              default = file.path(DATA_PROC, "meta_human_main.csv"), help = "主队列表型"),
  make_option(c("--subtype"), type = "character",
              default = file.path(DATA_PROC, "subtype_assignment.csv"), help = "08 分型"),
  make_option(c("--sig"), type = "character",
              default = file.path(DATA_PROC, "progression_signature.rds"), help = "06 signature"),
  make_option(c("--rda"), type = "character",
              default = file.path(DATA_REF, "signature_collection_IOBR.rda"), help = "IOBR 签名集合"),
  make_option(c("--outdir"), type = "character", default = DATA_PROC, help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

main <- function() {
  ensure_dirs()
  expr <- readRDS(opt$expr)
  meta <- read.csv(opt$meta, stringsAsFactors = FALSE)
  meta$disease <- ifelse(meta$status == "IgAN patient", "IgAN", "control")
  sig <- readRDS(opt$sig)

  e <- new.env()
  load(opt$rda, envir = e)
  sc <- e$signature_collection
  bindea <- grep("_Bindea_et_al$", names(sc), value = TRUE)
  tissue <- c("SW480_cancer_cells_Bindea_et_al", "Normal_mucosa_Bindea_et_al",
              "Blood_vessels_Bindea_et_al", "Lymph_vessels_Bindea_et_al")
  bindea <- setdiff(bindea, tissue)
  sets <- lapply(bindea, function(nm) intersect(sc[[nm]], rownames(expr)))
  names(sets) <- sub("_Bindea_et_al$", "", bindea)
  keep <- lengths(sets) >= 5
  sets <- sets[keep]
  log_msg("immune", paste0("Bindea signatures used = ", length(sets)))

  score <- gsva(ssgseaParam(exprData = as.matrix(expr), geneSets = sets))
  safe_write(file.path(opt$outdir, "immune_bindea_scores.csv"), function(p) {
    out <- data.frame(sample = colnames(score), t(score), check.names = FALSE)
    write_csv_utf8(out, p)
  })

  res <- do.call(rbind, lapply(rownames(score), function(nm) {
    x <- as.numeric(score[nm, ])
    data.frame(celltype = nm,
               IgAN = mean(x[meta$disease == "IgAN"]),
               control = mean(x[meta$disease == "control"]),
               p = wilcox.test(x[meta$disease == "IgAN"], x[meta$disease == "control"])$p.value)
  }))
  res$padj <- p.adjust(res$p, method = "BH")
  res <- res[order(res$p), ]
  safe_write(file.path(opt$outdir, "immune_deconvolution_results.csv"), function(p) {
    write_csv_utf8(res, p)
  })
  log_msg("immune", paste0("top: ", paste(head(res$celltype, 5), collapse = ", ")))

  sub <- read.csv(opt$subtype, stringsAsFactors = FALSE)
  common <- intersect(sub$sample, colnames(expr))
  df <- data.frame(subtype = factor(sub$subtype[match(common, sub$sample)]),
                   Macrophages = as.numeric(score["Macrophages", common]))
  p <- ggboxplot(df, x = "subtype", y = "Macrophages", color = "subtype",
                 palette = CB_PALETTE, add = "jitter", xlab = NULL,
                 ylab = "Macrophage score") +
    stat_compare_means(method = "wilcox.test") + pub_theme()
  save_pub_fig(p, file.path(FIG_DIR, "fig_subtype_macrophage"))

  # 热图：两组各细胞评分均值
  hm <- sapply(rownames(score), function(nm) {
    c(IgAN = mean(as.numeric(score[nm, meta$disease == "IgAN"])),
      control = mean(as.numeric(score[nm, meta$disease == "control"])))
  })
  hp <- file.path(FIG_DIR, "fig_immune_heatmap.pdf")
  hp_png <- file.path(FIG_DIR, "fig_immune_heatmap.png")
  for (f in c(hp, hp_png)) if (file.exists(f)) unlink(f)
  pdf(hp, width = 4.5, height = 7)
  pheatmap(hm, color = colorRampPalette(c("#0072B2", "white", "#E69F00"))(50),
           fontsize = 6, cluster_cols = FALSE)
  dev.off()
  png(hp_png, width = 4.5, height = 7, units = "in", res = 300)
  pheatmap(hm, color = colorRampPalette(c("#0072B2", "white", "#E69F00"))(50),
           fontsize = 6, cluster_cols = FALSE)
  dev.off()

  # 签名基因与免疫评分 Spearman 相关
  cor_mat <- cor(t(expr[sig, , drop = FALSE]), t(score), method = "spearman")
  safe_write(file.path(opt$outdir, "sig_immune_cell_cor.csv"), function(p) {
    write_csv_utf8(as.data.frame(cor_mat), p)
  })
  log_msg("immune", "DONE")
}

if (!interactive()) main()
