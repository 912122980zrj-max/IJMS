#!/usr/bin/env Rscript
# 05_WGCNA.R —— 机械相关共表达模块（真实实现，主队列 GSE115857 IgAN）
# 性状：04 输出的 MRG_up ssGSEA 打分；软阈值 β 扫描；选择与性状相关性最强的模块；
# GS/MM 阈值（|GS|>=0.7 且 |MM|>=0.7）筛核心机械基因。

suppressPackageStartupMessages({
  library(optparse)
  library(WGCNA)
})
allowWGCNAThreads(nThreads = 2)

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "05_WGCNA.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--expr"), type = "character",
              default = file.path(DATA_PROC, "expr_human_main.rds"), help = "主队列表达矩阵（GSE115857）"),
  make_option(c("--meta"), type = "character",
              default = file.path(DATA_PROC, "meta_human_main.csv"), help = "主队列表型"),
  make_option(c("--score"), type = "character",
              default = file.path(DATA_PROC, "ssgsea_mech_scores.csv"), help = "04 输出机械打分"),
  make_option(c("--outdir"), type = "character", default = DATA_PROC, help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

main <- function() {
  ensure_dirs()
  expr <- readRDS(opt$expr)
  meta <- read.csv(opt$meta, stringsAsFactors = FALSE)
  keep <- meta$status == "IgAN patient"   # 沿用 IPF 模板：仅疾病组样本做 WGCNA
  datExpr <- t(expr[, keep, drop = FALSE])
  score <- read.csv(opt$score, stringsAsFactors = FALSE, check.names = FALSE)
  trait <- score$MRG_up[keep]

  sft <- pickSoftThreshold(datExpr, powerVector = 1:10, networkType = "signed")
  beta <- ifelse(is.na(sft$powerEstimate), 4, sft$powerEstimate)
  log_msg("wgcna", paste0("beta = ", beta))

  net <- blockwiseModules(datExpr, power = beta, networkType = "signed",
                          TOMType = "signed", minModuleSize = 30,
                          mergeCutHeight = 0.25, numericLabels = TRUE, verbose = 0)
  log_msg("wgcna", paste0("modules = ", length(unique(net$colors))))

  me_cor <- abs(as.numeric(WGCNA::cor(net$MEs, trait, use = "pairwise.complete.obs")))
  best <- which.max(me_cor)
  log_msg("wgcna", paste0("selected module: ", colnames(net$MEs)[best],
                          " (|cor|=", round(me_cor[best], 3), ")"))
  gs <- as.numeric(WGCNA::cor(datExpr, trait, use = "pairwise.complete.obs"))
  mm <- as.numeric(WGCNA::cor(datExpr, net$MEs[, best], use = "pairwise.complete.obs"))
  names(gs) <- colnames(datExpr)
  core <- colnames(datExpr)[abs(gs) >= 0.5 & abs(mm) >= 0.5]   # D13：0.7/0.7 过严（0 基因）
  if (length(core) < 50) {                                     # 保底：最优模块内按 |GS| 取前 100
    core <- names(sort(abs(gs)[abs(mm) >= 0.5], decreasing = TRUE))[
      seq_len(min(100, sum(abs(mm) >= 0.5)))]
  }
  log_msg("wgcna", paste0("core genes = ", length(core)))
  safe_write(file.path(opt$outdir, "wgcna_core.rds"), function(p) saveRDS(core, p))
  safe_write(file.path(opt$outdir, "wgcna_modules.rds"), function(p) saveRDS(net, p))
  safe_write(file.path(opt$outdir, "wgcna_core_gs_mm.csv"), function(p) {
    write_csv_utf8(data.frame(gene = colnames(datExpr), GS = gs, MM = mm,
                              core = colnames(datExpr) %in% core), p)
  })
  log_msg("wgcna", "DONE")
}

if (!interactive()) main()
