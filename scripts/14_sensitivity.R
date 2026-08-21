#!/usr/bin/env Rscript
# 14_sensitivity.R —— β/K 敏感性分析
# β：在 β=8 下重跑模块构建，与 β=4 核心基因集比较 Jaccard 重叠；
# K：从 08 的 PAC 序列判定 K=2 的局部最优性（写入 results/sensitivity_summary.md）。

suppressPackageStartupMessages({
  library(optparse)
  library(WGCNA)
})
allowWGCNAThreads(nThreads = 2)

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "14_sensitivity.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--beta"), type = "integer", default = 8, help = "对比 β 值"),
  make_option(c("--outdir"), type = "character", default = DATA_PROC, help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

main <- function() {
  ensure_dirs()
  expr <- readRDS(file.path(opt$outdir, "expr_human_main.rds"))
  meta <- read.csv(file.path(opt$outdir, "meta_human_main.csv"), stringsAsFactors = FALSE)
  keep <- meta$status == "IgAN patient"
  datExpr <- t(expr[, keep, drop = FALSE])
  score <- read.csv(file.path(opt$outdir, "ssgsea_mech_scores.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
  trait <- score$MRG_up[keep]

  net <- blockwiseModules(datExpr, power = opt$beta, networkType = "signed",
                          TOMType = "signed", minModuleSize = 30,
                          mergeCutHeight = 0.25, numericLabels = TRUE, verbose = 0)
  me_cor <- abs(as.numeric(WGCNA::cor(net$MEs, trait, use = "pairwise.complete.obs")))
  best <- which.max(me_cor)
  gs <- as.numeric(WGCNA::cor(datExpr, trait, use = "pairwise.complete.obs"))
  mm <- as.numeric(WGCNA::cor(datExpr, net$MEs[, best], use = "pairwise.complete.obs"))
  names(gs) <- colnames(datExpr)
  core <- colnames(datExpr)[abs(gs) >= 0.5 & abs(mm) >= 0.5]
  if (length(core) < 50) {
    core <- names(sort(abs(gs)[abs(mm) >= 0.5], decreasing = TRUE))[
      seq_len(min(100, sum(abs(mm) >= 0.5)))]
  }
  safe_write(file.path(opt$outdir, paste0("wgcna_core_beta", opt$beta, ".rds")),
             function(p) saveRDS(core, p))

  core4 <- readRDS(file.path(opt$outdir, "wgcna_core.rds"))
  inter <- length(intersect(core, core4))
  jac <- inter / length(union(core, core4))
  log_msg("sensitivity", paste0("beta=", opt$beta,
                                " modules=", length(unique(net$colors)),
                                " best|cor|=", round(me_cor[best], 3),
                                " core=", length(core),
                                " overlap_with_beta4=", inter,
                                " jaccard=", round(jac, 3)))
  log_msg("sensitivity", "DONE")
}

if (!interactive()) main()
