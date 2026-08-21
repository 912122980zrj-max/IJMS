#!/usr/bin/env Rscript
# 54_export_signature_scores.R —— 导出与 Fig2a 口径一致的签名预测值与逐基因 ROC
# 复现 06_progression_signature.R / 25_export_panel_data.R 的 LASSO(lambda.min) 预测，
# 口径统一为 85 样本（IgAN 55 + 对照 30，排除 FSGS），并额外导出逐样本签名得分，
# 供 Python 端用 bootstrap 给 Fig2a 的签名 ROC 加置信带与 AUC 95% CI。
#
# 输出（data/processed/panel_export/）:
#   signature_scores.csv   sample,score,y       （LASSO 线性预测 + 结局）
#   per_gene_roc_85.csv    gene,fpr,tpr        （NDNF/PCDHB7/RRAGB 各基因，85 集）
#   per_gene_auc_85.csv    gene,AUC,direction,n_iga,n_ctrl

suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
})

root <- "E:/sheng xin/ObstructiveNephropathy_MRG"
dp <- file.path(root, "data", "processed")
out <- file.path(dp, "panel_export")
if (!dir.exists(out)) dir.create(out, recursive = TRUE)

expr <- readRDS(file.path(dp, "expr_human_main.rds"))
meta <- read.csv(file.path(dp, "meta_human_main.csv"), stringsAsFactors = FALSE)
uni <- read.csv(file.path(dp, "univariate_igag.csv"), stringsAsFactors = FALSE)

# ---- 口径：85 样本（排除 FSGS），列线图/校准同源 ----
sel <- meta$status %in% c("IgAN patient", "minimal change disease",
                          "Membranous glomerulonephritis", "Living donor")
meta <- meta[sel, , drop = FALSE]
meta$disease <- ifelse(meta$status == "IgAN patient", "IgAN", "control")
expr <- expr[, meta$sample, drop = FALSE]
group <- factor(meta$disease, levels = c("control", "IgAN"))

keep <- uni$gene %in% rownames(expr)
x <- t(expr[uni$gene[keep], , drop = FALSE])
y <- group

# ---- LASSO（与 06/25 一致：top40 单因素基因 → cv.glmnet lambda.min）----
set.seed(42)
cv <- cv.glmnet(x, y, family = "binomial", nfolds = 10)
score <- as.numeric(predict(cv, newx = x, s = "lambda.min", type = "response"))
ro <- roc(y, score, quiet = TRUE)
cat(sprintf("[sig] overall AUC = %.4f  (n = %d, pos = %d, neg = %d)\n",
            as.numeric(auc(ro)), length(y),
            sum(y == "IgAN"), sum(y == "control")))

# 逐样本得分
score_df <- data.frame(sample = meta$sample,
                       score = score,
                       y = as.integer(y == "IgAN"))
write.csv(score_df, file.path(out, "signature_scores.csv"), row.names = FALSE)

# ---- 逐基因 ROC（85 集）----
genes <- c("NDNF", "PCDHB7", "RRAGB")
rows_auc <- list(); rows_roc <- list()
for (g in genes) {
  v <- as.numeric(expr[g, ])
  a <- as.numeric(auc(roc(y, v, quiet = TRUE)))
  dirn <- if (a >= 0.5) "up" else "down"
  aa <- if (dirn == "up") a else 1 - a
  # 为画图统一方向：disease 上调 gene 用原值，下调用负值，使 AUC>=0.5
  plotv <- if (dirn == "up") v else -v
  r <- roc(y, plotv, quiet = TRUE)
  rows_roc[[g]] <- data.frame(gene = g, fpr = 1 - r$specificities, tpr = r$sensitivities)
  rows_auc[[g]] <- data.frame(gene = g, AUC = aa, direction = dirn,
                              n_iga = sum(y == "IgAN"),
                              n_ctrl = sum(y == "control"))
}
per_auc <- do.call(rbind, rows_auc)
per_roc <- do.call(rbind, rows_roc)
write.csv(per_auc, file.path(out, "per_gene_auc_85.csv"), row.names = FALSE)
write.csv(per_roc, file.path(out, "per_gene_roc_85.csv"), row.names = FALSE)
cat("[sig] per-gene AUC:", paste(sprintf("%s=%.3f", per_auc$gene, per_auc$AUC), collapse=" "), "\n")
cat("[sig] DONE\n")
