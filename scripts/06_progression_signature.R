#!/usr/bin/env Rscript
# 06_progression_signature.R —— 疾病/纤维化轴机械基因 signature（真实实现）
# 数据现实：无生存/随访（见 deviation_log D1/D2），终点 = CKD vs normal（GSE66494，53 vs 8）。
# 流程（沿用 IPF 论文漏斗）：单因素 logistic + BH → LASSO（10 折 lambda.min）→ RF（重要性 top15）
#   → 交集核心 signature → ROC（整体 + 5 折 CV 交叉验证）。

suppressPackageStartupMessages({
  library(optparse)
  library(glmnet)
  library(randomForest)
  library(pROC)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "06_progression_signature.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--expr"), type = "character",
              default = file.path(DATA_PROC, "expr_human_main.rds"), help = "GSE115857 表达矩阵"),
  make_option(c("--meta"), type = "character",
              default = file.path(DATA_PROC, "meta_human_main.csv"), help = "GSE115857 表型"),
  make_option(c("--core"), type = "character",
              default = file.path(DATA_PROC, "wgcna_core.rds"), help = "05 输出的核心机械基因"),
  make_option(c("--outdir"), type = "character", default = DATA_PROC, help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

main <- function() {
  ensure_dirs()
  expr <- readRDS(opt$expr)
  meta <- read.csv(opt$meta, stringsAsFactors = FALSE)
  meta$disease <- ifelse(meta$status == "IgAN patient", "IgAN", "control")
  meta <- meta[meta$status %in% c("IgAN patient", "minimal change disease",
                                  "Membranous glomerulonephritis", "Living donor"), ]
  expr <- expr[, meta$sample]
  group <- factor(meta$disease, levels = c("control", "IgAN"))
  core <- readRDS(opt$core)
  genes <- intersect(core, rownames(expr))
  if (length(genes) == 0L) stop("无候选基因：请检查 05 输出的核心基因集")
  log_msg("screen", paste0("candidates = ", length(genes)))

  uni <- do.call(rbind, lapply(genes, function(g) {
    fit <- glm(group ~ x, data = data.frame(group = group, x = as.numeric(expr[g, ])),
               family = binomial())
    s <- summary(fit)$coefficients
    data.frame(gene = g, p = s[2, 4], or = exp(s[2, 1]))
  }))
  uni$padj <- p.adjust(uni$p, method = "BH")
  uni <- uni[order(uni$padj), ]
  # D18 QC：剔除 OR>20/<0.05（近完全分离）与 P 后缀假基因（如 OR2A20P，易为注释/交叉杂交伪影）
  drop <- uni$or > 20 | uni$or < 0.05 | grepl("P$", uni$gene)
  log_msg("qc", paste0("dropped extreme-OR/pseudogene candidates = ", sum(drop)))
  uni <- uni[!drop, ]
  top <- head(uni, 40)
  safe_write(file.path(opt$outdir, "univariate_igag.csv"), function(p) write_csv_utf8(top, p))
  log_msg("screen", paste0("univariate top40 (BH<0.05 n=", sum(top$padj < 0.05), ")"))

  x <- t(expr[top$gene, , drop = FALSE])

  set.seed(42)
  cv <- cv.glmnet(x, group, family = "binomial", nfolds = 10)
  co <- coef(cv, s = "lambda.min")
  lasso_genes <- rownames(co)[as.numeric(co) != 0][-1]
  log_msg("lasso", paste0("LASSO genes = ", length(lasso_genes)))

  set.seed(42)
  rf <- randomForest(x = x, y = group, ntree = 1000, importance = TRUE)
  rf_genes <- names(head(sort(rf$importance[, 1], decreasing = TRUE), 15))

  sig <- intersect(lasso_genes, rf_genes)
  log_msg("core", paste0("signature = ", paste(sig, collapse = ", ")))
  safe_write(file.path(opt$outdir, "progression_signature.rds"), function(p) saveRDS(sig, p))

  score <- as.numeric(predict(cv, newx = x, s = "lambda.min", type = "response"))
  roc_all <- roc(group, score, quiet = TRUE)
  log_msg("roc", paste0("overall AUC = ", round(auc(roc_all), 3)))

  set.seed(42)
  folds <- sample(rep(1:5, length.out = length(group)))
  cv_auc <- sapply(1:5, function(k) {
    tr <- folds != k
    cvm <- cv.glmnet(x[tr, , drop = FALSE], group[tr], family = "binomial", nfolds = 10)
    sc <- as.numeric(predict(cvm, newx = x[!tr, , drop = FALSE],
                             s = "lambda.min", type = "response"))
    as.numeric(auc(roc(group[!tr], sc, quiet = TRUE)))
  })
  log_msg("roc", paste0("5-fold CV AUC = ", round(mean(cv_auc), 3),
                        " ± ", round(sd(cv_auc), 3)))

  roc_df <- data.frame(spec = 1 - roc_all$specificities, sens = roc_all$sensitivities)
  p <- ggplot(roc_df, aes(spec, sens)) + geom_line(color = CB_PALETTE[1], linewidth = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    annotate("text", x = 0.75, y = 0.25, size = 2.5,
             label = paste0("AUC = ", round(auc(roc_all), 3))) +
    labs(x = "1 - specificity", y = "sensitivity", title = "Mechanical signature: IgAN vs control") +
    pub_theme()
  save_pub_fig(p, file.path(FIG_DIR, "fig_roc_igag_signature"))
  safe_write(file.path(opt$outdir, "signature_roc.csv"), function(p) {
    write_csv_utf8(data.frame(fold = c(0, 1:5),
                              AUC = c(as.numeric(auc(roc_all)), cv_auc)), p)
  })
  log_msg("signature", "DONE")
}

if (!interactive()) main()
