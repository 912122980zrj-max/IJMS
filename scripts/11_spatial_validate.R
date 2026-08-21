#!/usr/bin/env Rscript
# 11_spatial_validate.R —— 人肾 CosMx 空间验证（轻量实现，GSE282059）
# 检验：机械核心基因评分（MRG）与纤维化 ECM 评分是否在空间上共定位——
#   细胞级 Spearman 相关 + FOV 级相关 + 高 ECM 细胞的高 MRG 富集（wilcoxon）。

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "11_spatial_validate.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--exprmat"), type = "character",
              default = file.path(DATA_RAW, "GSE282059", "suppl", "GSE282059_2B_FFPE_exprMat_file.csv.gz"),
              help = "CosMx 表达矩阵（fov, cell_ID, genes...）"),
  make_option(c("--meta"), type = "character",
              default = file.path(DATA_RAW, "GSE282059", "suppl", "GSE282059_2B_FFPE_metadata_file.csv.gz"),
              help = "CosMx 元数据"),
  make_option(c("--core"), type = "character",
              default = file.path(DATA_PROC, "wgcna_core.rds"), help = "05 核心机械基因"),
  make_option(c("--outdir"), type = "character", default = DATA_PROC, help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

main <- function() {
  ensure_dirs()
  core <- readRDS(opt$core)
  ecm <- c("COL1A1", "COL3A1", "DCN", "FN1", "ACTA2", "COL4A1", "COL6A1", "VCAN")
  first <- fread(opt$exprmat, nrows = 1)
  genes <- setdiff(names(first), c("fov", "cell_ID"))
  use <- intersect(union(core, ecm), genes)
  log_msg("load", paste0("genes used: MRG ", sum(core %in% use),
                         " | ECM ", sum(ecm %in% use)))
  x <- fread(opt$exprmat, select = c("fov", "cell_ID", use))
  log_msg("load", paste0("cells = ", nrow(x)))
  mrg_cols <- use[use %in% core]
  ecm_cols <- use[use %in% ecm]

  z <- function(m) {
    m <- log1p(as.matrix(m))
    s <- apply(m, 2, stats::sd)
    keep <- s > 0
    if (any(!keep)) m <- m[, keep, drop = FALSE]
    rowMeans(scale(m))
  }
  mrg <- z(x[, ..mrg_cols])
  ecmv <- z(x[, ..ecm_cols])
  res <- data.table(fov = x$fov, cell_ID = x$cell_ID, MRG = mrg, ECM = ecmv)

  rho_cell <- cor(res$MRG, res$ECM, method = "spearman")
  fov_agg <- res[, .(MRG = mean(MRG), ECM = mean(ECM), n = .N), by = fov]
  rho_fov <- cor(fov_agg$MRG, fov_agg$ECM, method = "spearman")
  hi <- res$ECM >= quantile(res$ECM, 0.75)
  p_hi <- wilcox.test(res$MRG[hi], res$MRG[!hi])$p.value
  log_msg("spatial", paste0("cell-level spearman rho = ", round(rho_cell, 3),
                            " | FOV-level rho = ", round(rho_fov, 3),
                            " | high-ECM enrichment p = ", signif(p_hi, 3)))

  safe_write(file.path(opt$outdir, "spatial_validation.csv"), function(p) {
    write_csv_utf8(data.frame(metric = c("cell_spearman", "fov_spearman", "high_ecm_wilcox_p"),
                              value = c(rho_cell, rho_fov, p_hi)), p)
  })
  set.seed(42)
  plt <- res[sample(.N, min(.N, 20000))]
  p <- ggplot(plt, aes(ECM, MRG)) +
    geom_point(alpha = 0.15, size = 0.2, color = CB_PALETTE[1]) +
    geom_smooth(method = "lm", color = CB_PALETTE[2], linewidth = 0.5) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.4, size = 2.5,
             label = paste0("rho = ", round(rho_cell, 3))) +
    labs(x = "Fibrotic ECM score", y = "Mechanical core score",
         title = "CosMx: mechanical-fibrotic co-localization") + pub_theme()
  save_pub_fig(p, file.path(FIG_DIR, "fig_spatial_mrg_ecm"))
  log_msg("spatial", "DONE")
}

if (!interactive()) main()
