#!/usr/bin/env Rscript
# 29_timecourse_panel_data.R
#
# 为 Fig 4D 制备逐样本机械核心评分（与 13_UUO_timecourse.R 完全同口径，只读原始数据）：
#   GSE118339 15 个 TPM 文件 -> Ensembl ID 单次 biomaRt 映射到小鼠符号 ->
#   保留人类机械核心基因的 one-to-one 小鼠同源集 -> log2(TPM+1) 按基因 z 化 -> 逐样本均值。
# 结果写入 data/processed/panel_export/timecourse_samples.csv（新文件，不覆盖）。

suppressPackageStartupMessages({
  library(optparse)
  library(biomaRt)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "29_timecourse_panel_data.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--indir"), type = "character",
              default = file.path(DATA_RAW, "GSE118339"), help = "TPM 文件目录"),
  make_option(c("--ortholog"), type = "character",
              default = file.path(DATA_PROC, "ortholog_map.csv"), help = "同源映射表"),
  make_option(c("--core"), type = "character",
              default = file.path(DATA_PROC, "wgcna_core_gs_mm.csv"),
              help = "WGCNA 核心基因表（含 core 标志）"),
  make_option(c("--outdir"), type = "character",
              default = file.path(DATA_PROC, "panel_export"), help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

main <- function() {
  ensure_dirs()
  if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

  files <- list.files(opt$indir, pattern = "tpm.txt.gz$", full.names = TRUE)
  mats <- lapply(files, function(f) {
    d <- read.csv(gzfile(f), header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE)
    setNames(d[[2]], d[[1]])
  })
  genes <- sort(unique(unlist(lapply(mats, names))))
  M <- do.call(cbind, lapply(mats, function(m) m[genes]))
  colnames(M) <- sub("GSM\\d+_(.+)\\.tpm\\.txt\\.gz", "\\1", basename(files))
  log_msg("tc-panel", paste0("samples = ", ncol(M), " | genes = ", nrow(M)))

  ensembl <- NULL
  for (mirror in c("useast", "asia", "www")) {
    ensembl <- tryCatch(
      useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror = mirror),
      error = function(e) NULL)
    if (!is.null(ensembl)) break
  }
  if (is.null(ensembl)) stop("Ensembl 镜像均不可达")
  ann <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
               filters = "ensembl_gene_id", values = genes, mart = ensembl)
  log_msg("tc-panel", paste0("ensembl annotation rows = ", nrow(ann)))

  orth <- read.csv(opt$ortholog, stringsAsFactors = FALSE)
  orth <- orth[orth$mmusculus_homolog_orthology_type == "ortholog_one2one", ]
  core <- read.csv(opt$core, stringsAsFactors = FALSE)
  core_genes <- core$gene[core$core]
  mouse_core <- unique(orth$mmusculus_homolog_associated_gene_name[
    orth$hgnc_symbol %in% core_genes])
  keep_ids <- ann$ensembl_gene_id[ann$external_gene_name %in% mouse_core]
  row_sel <- rownames(M) %in% keep_ids
  log_msg("tc-panel", paste0("core mouse genes present = ", sum(row_sel)))

  v <- log2(M[row_sel, , drop = FALSE] + 1)
  zmat <- t(scale(t(v)))
  score <- colMeans(zmat, na.rm = TRUE)

  day_map <- c("normal_1" = 0, "normal_2" = 0, "normal_3" = 0,
               "day3_4" = 3, "day3_5" = 3, "day3_6" = 3, "day3_7" = 3,
               "day7_8" = 7, "day7_9" = 7, "day7_10" = 7, "day7_11" = 7,
               "day14_12" = 14, "day14_13" = 14, "day14_14" = 14,
               "day14_15" = 14)
  out <- data.frame(sample = colnames(M),
                    day = unname(day_map[colnames(M)]),
                    score = score)
  out_path <- file.path(opt$outdir, "timecourse_samples.csv")
  if (file.exists(out_path)) stop("Refusing to overwrite: ", out_path)
  write.csv(out, out_path, row.names = FALSE, fileEncoding = "UTF-8")
  log_msg("tc-panel", paste0("written -> ", out_path))
  log_msg("tc-panel", "DONE")
}

if (!interactive()) main()
