#!/usr/bin/env Rscript
# 17_immune_consistency.R —— Charoentong 28 vs Bindea 24 反卷积一致性对比
# 指标：① 组均值差异方向一致（IgAN-control 符号）；② 逐样本评分 Spearman 相关。

suppressPackageStartupMessages({
  library(optparse)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "17_immune_consistency.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--outdir"), type = "character", default = DATA_PROC, help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

main <- function() {
  ensure_dirs()
  ch <- read.csv(file.path(opt$outdir, "immune28_charoentong_scores.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
  bn <- read.csv(file.path(opt$outdir, "immune_bindea_scores.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)

  get_score <- function(df, nm) {
    as.numeric(df[[nm]])
  }

  pairs <- list(
    list(c("Activated B cell"), "B_cells", "exact"),
    list(c("Activated CD8 T cell"), "CD8_T_cells", "exact"),
    list(c("Gamma delta T cell"), "Tgd", "exact"),
    list(c("Immature dendritic cell"), "iDC", "exact"),
    list(c("Activated dendritic cell"), "aDC", "exact"),
    list(c("Plasmacytoid dendritic cell"), "pDC", "exact"),
    list(c("Eosinophil"), "Eosinophils", "exact"),
    list(c("Macrophage"), "Macrophages", "exact"),
    list(c("Mast cell"), "Mast_cells", "exact"),
    list(c("Natural killer cell"), "NK_cells", "exact"),
    list(c("CD56dim natural killer cell"), "NK_CD56dim_cells", "exact"),
    list(c("Neutrophil"), "Neutrophils", "exact"),
    list(c("Regulatory T cell"), "Treg", "exact"),
    list(c("T follicular helper cell"), "TFH", "exact"),
    list(c("Type 1 T helper cell"), "Th1_cells", "exact"),
    list(c("Type 2 T helper cell"), "Th2_cells", "exact"),
    list(c("Type 17 T helper cell"), "Th17_cells", "exact"),
    list(c("Activated CD4 T cell"), "T_helper_cells", "approx"),
    list(c("Central memory CD4 T cell", "Central memory CD8 T cell"), "Tcm", "approx_pool"),
    list(c("Effector memeory CD4 T cell", "Effector memeory CD8 T cell"), "Tem", "approx_pool"),
    list(c("Memory B cell"), "B_cells", "approx")
  )

  # 组方向：用分组均值差（样本顺序与 meta 对齐）
  meta <- read.csv(file.path(opt$outdir, "meta_human_main.csv"), stringsAsFactors = FALSE)
  meta$disease <- ifelse(meta$status == "IgAN patient", "IgAN", "control")
  m <- match(ch$sample, meta$sample)
  grp <- meta$disease[m]
  ok <- vapply(pairs, function(pp) {
    all(pp[[1]] %in% colnames(ch)) && pp[[2]] %in% colnames(bn)
  }, logical(1))
  log_msg("consistency", paste0("pairs kept = ", sum(ok), " / skipped = ", sum(!ok)))
  pairs <- pairs[ok]
  res <- do.call(rbind, lapply(pairs, function(pp) {
    cnm <- pp[[1]]; bnm <- pp[[2]]; typ <- pp[[3]]
    x <- if (length(cnm) == 1) get_score(ch, cnm) else rowMeans(sapply(cnm, get_score, df = ch))
    y <- get_score(bn, bnm)
    dc <- sign(mean(x[grp == "IgAN"]) - mean(x[grp == "control"])) ==
          sign(mean(y[grp == "IgAN"]) - mean(y[grp == "control"]))
    data.frame(charoentong = paste(cnm, collapse = "+"),
               bindea = bnm, type = typ, dir_concordant = dc,
               rho = round(cor(x, y, method = "spearman"), 3))
  }))
  safe_write(file.path(opt$outdir, "immune_consistency_comparison.csv"),
             function(p) write_csv_utf8(res, p))
  log_msg("consistency", paste0("pairs = ", nrow(res),
                                " | direction concordant = ",
                                sum(res$dir_concordant), "/", nrow(res),
                                " | median rho = ", round(median(res$rho), 3)))
  log_msg("consistency", "DONE")
}

if (!interactive()) main()
