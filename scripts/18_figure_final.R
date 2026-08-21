#!/usr/bin/env Rscript
# 18_figure_final.R —— 手稿终版作图（修复命名重叠、空图、热图比例、P值格式、UMAP字体）
# 输出全部替换 figures/ 下同名文件；绘图规范统一由 00_config.R 提供。

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(pheatmap)
  library(glmnet)
  library(pROC)
  library(rms)
  library(Seurat)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "18_figure_final.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--outdir"), type = "character", default = DATA_PROC, help = "数据目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

pval_lab <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.0001) return("p < 0.0001")
  sprintf("p = %.4f", p)
}

box_plot <- function(df, xcol, ycol, xlab = NULL, ylab, pval) {
  df$x <- df[[xcol]]
  df$y <- df[[ycol]]
  ggplot(df, aes(x, y, color = x)) +
    geom_boxplot(outlier.size = 0.5, width = 0.55) +
    geom_jitter(width = 0.12, size = 0.6, alpha = 0.7) +
    scale_color_manual(values = CB_PALETTE, guide = "none") +
    labs(x = xlab, y = ylab) +
    annotate("text", x = 1.5, y = max(df$y, na.rm = TRUE) * 1.06, size = 2.5,
             label = pval_lab(pval)) +
    pub_theme()
}

main <- function() {
  ensure_dirs()
  expr <- readRDS(file.path(opt$outdir, "expr_human_main.rds"))
  meta <- read.csv(file.path(opt$outdir, "meta_human_main.csv"), stringsAsFactors = FALSE)
  meta$disease <- ifelse(meta$status == "IgAN patient", "IgAN", "control")
  score <- read.csv(file.path(opt$outdir, "ssgsea_mech_scores.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
  sub <- read.csv(file.path(opt$outdir, "subtype_assignment.csv"), stringsAsFactors = FALSE)
  hallmark <- read.csv(file.path(opt$outdir, "hallmark_scores.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  imm28 <- read.csv(file.path(opt$outdir, "immune28_charoentong_scores.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
  sig <- readRDS(file.path(opt$outdir, "progression_signature.rds"))

  log_msg("fig1", "mechanical score boxplot")
  p1 <- box_plot(data.frame(group = factor(meta$disease, levels = c("control", "IgAN")),
                            value = score$MRG_up),
                 "group", "value", ylab = "ssGSEA score (MRG up)",
                 pval = wilcox.test(score$MRG_up ~ meta$disease)$p.value)
  save_pub_fig(p1, file.path(FIG_DIR, "fig_mrg_score_igag"))

  log_msg("fig2a", "ROC curve (recompute with seed 42)")
  uni <- read.csv(file.path(opt$outdir, "univariate_igag.csv"), stringsAsFactors = FALSE)
  keep <- uni$gene %in% rownames(expr)
  x <- t(expr[uni$gene[keep], meta$sample, drop = FALSE])
  y <- factor(meta$disease, levels = c("control", "IgAN"))
  set.seed(42)
  cv <- cv.glmnet(x, y, family = "binomial", nfolds = 10)
  pr <- as.numeric(predict(cv, newx = x, s = "lambda.min", type = "response"))
  ro <- roc(y, pr, quiet = TRUE)
  roc_df <- data.frame(spec = 1 - ro$specificities, sens = ro$sensitivities)
  p2 <- ggplot(roc_df, aes(spec, sens)) + geom_line(color = CB_PALETTE[1], linewidth = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    annotate("text", x = 0.72, y = 0.24, size = 2.5,
             label = paste0("AUC = ", sprintf("%.3f", auc(ro)))) +
    labs(x = "1 - specificity", y = "sensitivity", title = "Mechanical signature (IgAN vs control)") +
    pub_theme()
  save_pub_fig(p2, file.path(FIG_DIR, "fig_roc_igag_signature"))

  log_msg("fig2b", "nomogram (wide layout)")
  idx <- meta$disease %in% c("control", "IgAN")
  dat <- data.frame(IgAN = as.integer(meta$disease[idx] == "IgAN"),
                    t(expr[intersect(sig, rownames(expr)), meta$sample[idx], drop = FALSE]))
  full <- glm(IgAN ~ ., data = dat, family = binomial())
  step <- MASS::stepAIC(full, direction = "both", trace = 0)
  assign("dd", datadist(dat), envir = .GlobalEnv)
  options(datadist = "dd")
  lrm_fit <- lrm(IgAN ~ ., data = dat, x = TRUE, y = TRUE)
  nomo <- nomogram(lrm_fit, fun = plogis, funlabel = "IgAN probability",
                   fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9))
  pdf_path <- file.path(FIG_DIR, "fig_nomogram_igag.pdf")
  png_path <- file.path(FIG_DIR, "fig_nomogram_igag.png")
  for (f in c(pdf_path, png_path)) if (file.exists(f)) unlink(f)
  pdf(pdf_path, width = 10, height = 5)
  par(cex = 0.85, cex.axis = 0.8, mar = c(5, 4, 3, 2))
  plot(nomo)
  dev.off()
  png(png_path, width = 10, height = 5, units = "in", res = 300)
  par(cex = 0.85, cex.axis = 0.8, mar = c(5, 4, 3, 2))
  plot(nomo)
  dev.off()

  log_msg("fig3a", "subtype mechanical score")
  cls <- setNames(sub$subtype, sub$sample)
  d3a <- data.frame(subtype = factor(paste0("C", cls[meta$sample[meta$disease == "IgAN"]])),
                    value = score$MRG_up[meta$disease == "IgAN"])
  p3a <- box_plot(d3a, "subtype", "value", ylab = "ssGSEA score (MRG up)",
                  pval = wilcox.test(value ~ subtype, data = d3a)$p.value)
  save_pub_fig(p3a, file.path(FIG_DIR, "fig_subtype_mech_score"))

  log_msg("fig3b", "EMT subtype boxplot (name-index fix)")
  emt <- as.numeric(hallmark[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]])
  names(emt) <- hallmark$sample
  ig <- meta$sample[meta$disease == "IgAN"]
  d3b <- data.frame(subtype = factor(paste0("C", cls[ig])),
                    value = emt[ig])
  p3b <- box_plot(d3b, "subtype", "value", ylab = "EMT pathway score",
                  pval = wilcox.test(value ~ subtype, data = d3b)$p.value)
  save_pub_fig(p3b, file.path(FIG_DIR, "fig_subtype_emt"))

  log_msg("fig3c", "immune28 heatmap (wide layout)")
  cells28 <- setdiff(colnames(imm28), "sample")
  hm <- sapply(cells28, function(nm) {
    c(IgAN = mean(as.numeric(imm28[[nm]][meta$disease == "IgAN"])),
      control = mean(as.numeric(imm28[[nm]][meta$disease == "control"])))
  })
  hp <- file.path(FIG_DIR, "fig_immune28_heatmap.pdf")
  hp_png <- file.path(FIG_DIR, "fig_immune28_heatmap.png")
  for (f in c(hp, hp_png)) if (file.exists(f)) unlink(f)
  pdf(hp, width = 11, height = 2.6)
  pheatmap(hm, color = colorRampPalette(c("#0072B2", "white", "#E69F00"))(50),
           fontsize = 7, cluster_cols = FALSE, cluster_rows = FALSE,
           angle_col = 45, border_color = NA)
  dev.off()
  png(hp_png, width = 11, height = 2.6, units = "in", res = 300)
  pheatmap(hm, color = colorRampPalette(c("#0072B2", "white", "#E69F00"))(50),
           fontsize = 7, cluster_cols = FALSE, cluster_rows = FALSE,
           angle_col = 45, border_color = NA)
  dev.off()

  log_msg("fig3d", "central memory CD8 subtype boxplot")
  ct <- "Central memory CD8 T cell"
  xct <- as.numeric(imm28[[ct]])
  names(xct) <- imm28$sample
  d3d <- data.frame(subtype = factor(paste0("C", cls[ig])), value = xct[ig])
  p3d <- box_plot(d3d, "subtype", "value", ylab = paste0(ct, " score"),
                  pval = wilcox.test(value ~ subtype, data = d3d)$p.value)
  save_pub_fig(p3d, file.path(FIG_DIR, "fig_subtype_Central_memory_CD8_T_cell"))

  log_msg("fig4a/b", "UMAP and MRG violin")
  obj <- readRDS(file.path(opt$outdir, "uuo_seurat.rds"))
  p4a <- DimPlot(obj, group.by = "celltype", label = TRUE, repel = TRUE,
                 label.size = 2.5, label.box = FALSE) +
    scale_color_manual(values = rep(CB_PALETTE, 2)) + pub_theme()
  save_pub_fig(p4a, file.path(FIG_DIR, "fig_uuo_umap_celltype"))

  d4b <- obj@meta.data
  p4b <- ggplot(d4b, aes(group, MRG1, color = group)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.size = 0.2) +
    scale_color_manual(values = CB_PALETTE, guide = "none") +
    labs(x = NULL, y = "Mechanical core module score") +
    annotate("text", x = 1.5, y = max(d4b$MRG1, na.rm = TRUE) * 1.08, size = 2.5,
             label = pval_lab(wilcox.test(MRG1 ~ group, data = d4b)$p.value)) +
    pub_theme()
  save_pub_fig(p4b, file.path(FIG_DIR, "fig_uuo_mrg_score"))

  log_msg("fig", "DONE")
}

if (!interactive()) main()
