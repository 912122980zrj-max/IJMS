#!/usr/bin/env Rscript
# 08_consensus_clustering.R —— 机械基因分子分型 C1/C2（真实实现）
# 输入：GSE66494 CKD 样本（53）+ signature 基因；ConsensusClusterPlus
# （pItem=pFeature=0.8, maxK=9, reps=1000, pearson 距离, seed=42），PAC 判 K。

suppressPackageStartupMessages({
  library(optparse)
  library(ConsensusClusterPlus)
  library(ggplot2)
  library(ggpubr)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "08_consensus_clustering.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--expr"), type = "character",
              default = file.path(DATA_PROC, "expr_human_main.rds"), help = "GSE115857 表达矩阵"),
  make_option(c("--meta"), type = "character",
              default = file.path(DATA_PROC, "meta_human_main.csv"), help = "GSE115857 表型"),
  make_option(c("--sig"), type = "character",
              default = file.path(DATA_PROC, "progression_signature.rds"), help = "06 输出 signature"),
  make_option(c("--score"), type = "character",
              default = file.path(DATA_PROC, "ssgsea_mech_scores.csv"), help = "04 机械打分"),
  make_option(c("--outdir"), type = "character", default = DATA_PROC, help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

pac_score <- function(cm) {
  x <- cm[lower.tri(cm)]
  mean(x > 0.1 & x < 0.9)
}

main <- function() {
  ensure_dirs()
  expr <- readRDS(opt$expr)
  meta <- read.csv(opt$meta, stringsAsFactors = FALSE)
  ckd <- meta$status == "IgAN patient"
  sig <- readRDS(opt$sig)
  dat <- as.matrix(expr[intersect(sig, rownames(expr)), ckd, drop = FALSE])
  log_msg("cc", paste0("features = ", nrow(dat), " | samples = ", ncol(dat)))

  set.seed(42)
  res <- ConsensusClusterPlus(dat, maxK = 9, reps = 1000, pItem = 0.8, pFeature = 0.8,
                              clusterAlg = "hc", distance = "pearson", seed = 42)
  pacs <- sapply(res[2:9], function(x) pac_score(x$consensusMatrix))
  log_msg("cc", paste0("PAC (K=2..9) = ", paste(round(pacs, 3), collapse = ", ")))
  cls <- res[[2]]$consensusClass
  safe_write(file.path(opt$outdir, "subtype_assignment.csv"), function(p) {
    write_csv_utf8(data.frame(sample = names(cls), subtype = paste0("C", cls)), p)
  })

  score <- read.csv(opt$score, stringsAsFactors = FALSE, check.names = FALSE)
  score <- score[match(names(cls), score$sample), ]
  df <- data.frame(subtype = factor(paste0("C", cls), levels = c("C1", "C2")),
                   score = score$MRG_up)
  p <- ggboxplot(df, x = "subtype", y = "score", color = "subtype",
                 palette = CB_PALETTE, add = "jitter", xlab = NULL,
                 ylab = "ssGSEA score (MRG up)") +
    stat_compare_means(method = "wilcox.test") + pub_theme()
  save_pub_fig(p, file.path(FIG_DIR, "fig_subtype_mech_score"))
  log_msg("cc", paste0("C1 n=", sum(cls == 1), " | C2 n=", sum(cls == 2),
                       " | wilcox p=", signif(wilcox.test(score ~ subtype, data = df)$p.value, 3)))
  log_msg("cc", "DONE")
}

if (!interactive()) main()
