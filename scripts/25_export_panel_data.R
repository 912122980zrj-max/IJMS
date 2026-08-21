#!/usr/bin/env Rscript
# 25_export_panel_data.R —— 组图重绘所需面板底层数据导出（只读，不改动结果文件）
# 输出到 data/processed/panel_export/：ROC 曲线、列线图系数与校准、UMAP 坐标。
# 字体规范：PLOS 要求图内文字 8-12pt，故组图由 24_composite_figures.py 以 8pt+ 重绘，
# 本脚本只负责把 R 端不可直接由 CSV 得到的量导出为 CSV。

suppressPackageStartupMessages({
  library(optparse)
  library(glmnet)
  library(pROC)
  library(rms)
  library(MASS)
  library(Seurat)
})

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
this_script <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else getwd()
project_root <- dirname(dirname(normalizePath(this_script)))

option_list <- list(
  make_option(c("--root"), type = "character", default = project_root, help = "项目根目录"),
  make_option(c("--outdir"), type = "character",
              default = file.path(project_root, "data", "processed", "panel_export"),
              help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

main <- function() {
  root <- normalizePath(opt$root)
  dp <- file.path(root, "data", "processed")
  out <- normalizePath(opt$outdir)
  if (!dir.exists(out)) dir.create(out, recursive = TRUE)

  expr <- readRDS(file.path(dp, "expr_human_main.rds"))
  meta <- read.csv(file.path(dp, "meta_human_main.csv"), stringsAsFactors = FALSE)
  meta$disease <- ifelse(meta$status == "IgAN patient", "IgAN", "control")
  uni <- read.csv(file.path(dp, "univariate_igag.csv"), stringsAsFactors = FALSE)
  sig <- readRDS(file.path(dp, "progression_signature.rds"))

  # ---- ROC 曲线（与 06/18 一致：top40 单因素基因 → cv.glmnet λ.min） ----
  keep <- uni$gene %in% rownames(expr)
  x <- t(expr[uni$gene[keep], meta$sample, drop = FALSE])
  y <- factor(meta$disease, levels = c("control", "IgAN"))
  set.seed(42)
  cv <- cv.glmnet(x, y, family = "binomial", nfolds = 10)
  pr <- as.numeric(predict(cv, newx = x, s = "lambda.min", type = "response"))
  ro <- roc(y, pr, quiet = TRUE)
  write.csv(data.frame(spec = 1 - ro$specificities, sens = ro$sensitivities),
            file.path(out, "roc_curve.csv"), row.names = FALSE)
  message("[panel-export] roc_curve.csv (AUC=", round(auc(ro), 3), ")")

  # ---- 列线图：逐步 logistic（按 Methods 排除 FSGS 的口径，85 样本） ----
  sel <- meta$status %in% c("IgAN patient", "minimal change disease",
                            "Membranous glomerulonephritis", "Living donor")
  dat <- data.frame(IgAN = as.integer(meta$status == "IgAN patient")[sel],
                    t(expr[intersect(sig, rownames(expr)), meta$sample[sel], drop = FALSE]))
  full <- glm(IgAN ~ ., data = dat, family = binomial())
  step <- MASS::stepAIC(full, direction = "both", trace = 0)
  terms <- names(coef(step))[-1]
  coefs <- data.frame(term = c("Intercept", terms),
                      coef = unname(coef(step)))
  ranges <- do.call(rbind, lapply(terms, function(g) {
    data.frame(term = g, min = min(dat[[g]]), max = max(dat[[g]]))
  }))
  write.csv(coefs, file.path(out, "nomogram_coefs.csv"), row.names = FALSE)
  write.csv(ranges, file.path(out, "nomogram_ranges.csv"), row.names = FALSE)

  assign("dd", datadist(dat), envir = .GlobalEnv)
  options(datadist = "dd")
  lf <- lrm(IgAN ~ ., data = dat, x = TRUE, y = TRUE)
  pred <- predict(lf, type = "fitted")
  cc <- val.prob(pred, dat$IgAN, pl = FALSE)
  write.csv(data.frame(pred = pred, obs = dat$IgAN),
            file.path(out, "calibration.csv"), row.names = FALSE)
  message("[panel-export] nomogram terms=", paste(terms, collapse = "+"),
          " C=", round(unname(cc["C (ROC)"]), 3),
          " slope=", round(unname(cc["Slope"]), 3))

  # ---- 主队列样本级面板数据（Fig 1B/3A-3C/3D 共用） ----
  score <- read.csv(file.path(dp, "ssgsea_mech_scores.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
  hallmark <- read.csv(file.path(dp, "hallmark_scores.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  sub <- read.csv(file.path(dp, "subtype_assignment.csv"), stringsAsFactors = FALSE)
  imm <- read.csv(file.path(dp, "immune28_charoentong_scores.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
  panel_df <- meta[, c("sample", "status", "disease")]
  panel_df <- merge(panel_df, score[, c("sample", "MRG_up")], by = "sample")
  panel_df <- merge(panel_df, hallmark[, c("sample", "HALLMARK_TGF_BETA_SIGNALING",
                                           "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")],
                    by = "sample")
  panel_df <- merge(panel_df, sub, by = "sample", all.x = TRUE)
  panel_df <- merge(panel_df, imm, by = "sample")
  write.csv(panel_df, file.path(out, "main_panel_df.csv"), row.names = FALSE)

  # ---- UMAP 坐标（uuo_seurat.rds） ----
  obj <- readRDS(file.path(dp, "uuo_seurat.rds"))
  emb <- as.data.frame(Seurat::Embeddings(obj, reduction = "umap"))
  names(emb) <- c("umap_1", "umap_2")
  umap_df <- data.frame(cell = rownames(emb), emb,
                        group = obj$group, celltype = obj$celltype, MRG1 = obj$MRG1)
  write.csv(umap_df, file.path(out, "umap.csv"), row.names = FALSE)
  message("[panel-export] umap.csv (n=", nrow(umap_df), ")")
  message("[panel-export] DONE")
}

if (!interactive()) main()
