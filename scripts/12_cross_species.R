#!/usr/bin/env Rscript
# 12_cross_species.R —— 人鼠同源映射 + UUO 方向验证（真实实现）
# 将人类 signature 经 biomaRt one2one 映射到小鼠，在 GSE299417（UUO_Veh vs Sham_Veh）
# 中核对 log2FC 方向一致性（人类 CKD 上调 → 小鼠 UUO 上调）。

suppressPackageStartupMessages({
  library(optparse)
  library(biomaRt)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "12_cross_species.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--sig"), type = "character",
              default = file.path(DATA_PROC, "progression_signature.rds"), help = "06 signature"),
  make_option(c("--core"), type = "character",
              default = file.path(DATA_PROC, "wgcna_core.rds"),
              help = "05 核心基因（跨物种验证主体，功效更高）"),
  make_option(c("--deg-human"), dest = "deg_human", type = "character",
              default = file.path(DATA_PROC, "deg_igag_vs_control.csv"), help = "人类 IgAN DEG 表（取方向）"),
  make_option(c("--deg-mouse"), dest = "deg_mouse", type = "character",
              default = file.path(DATA_PROC, "deg_mouse_uuo.csv"), help = "小鼠 UUO DEG 表"),
  make_option(c("--outdir"), type = "character", default = DATA_PROC, help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

main <- function() {
  ensure_dirs()
  sig <- readRDS(opt$sig)
  core <- readRDS(opt$core)
  genes_use <- union(sig, core)
  dh <- read.csv(opt$deg_human, stringsAsFactors = FALSE, check.names = FALSE)
  dm <- read.csv(opt$deg_mouse, stringsAsFactors = FALSE, check.names = FALSE)
  human_dir <- setNames(sign(dh$logFC), dh$gene)

  log_msg("ortholog", "biomaRt 人→鼠 one2one 映射（联网）")
  ensembl <- NULL
  for (mirror in c("useast", "asia", "www")) {
    ensembl <- tryCatch(useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",
                                   mirror = mirror),
                        error = function(e) NULL)
    if (!is.null(ensembl)) break
  }
  if (is.null(ensembl)) stop("所有 Ensembl 镜像均不可达（网络受限），请设代理后重试")
  orth <- do.call(rbind, lapply(genes_use, function(g) {
    tryCatch(
      getBM(attributes = c("mmusculus_homolog_associated_gene_name",
                           "mmusculus_homolog_orthology_type"),
            filters = "hgnc_symbol", values = g, mart = ensembl) |>
        transform(hgnc_symbol = g),
      error = function(e) NULL)
  }))
  orth <- unique(orth[orth$mmusculus_homolog_orthology_type == "ortholog_one2one", ])
  log_msg("ortholog", paste0("mapped ", length(unique(orth$hgnc_symbol)), "/", length(genes_use)))
  safe_write(file.path(opt$outdir, "ortholog_map.csv"), function(p) write_csv_utf8(orth, p))

  # 4v4 小鼠功效低：主指标为全部映射基因的趋势一致性；另报告达 p<0.05 的个数
  mm <- dm[dm$gene %in% orth$mmusculus_homolog_associated_gene_name,
           c("gene", "logFC", "adj.P.Val")]
  mm$human <- orth$hgnc_symbol[match(mm$gene, orth$mmusculus_homolog_associated_gene_name)]
  mm$human_dir <- human_dir[mm$human]
  mm$consistent <- sign(mm$logFC) == mm$human_dir
  mm <- mm[!is.na(mm$consistent), ]
  safe_write(file.path(opt$outdir, "ortholog_concordance.csv"), function(p) write_csv_utf8(mm, p))
  log_msg("cross_species", paste0("concordant ", sum(mm$consistent), "/", nrow(mm),
                                  " (", round(100 * mean(mm$consistent), 1), "%)"))
  log_msg("cross_species", paste0("significant (padj<0.05): ",
                                  sum(mm$adj.P.Val < 0.05), "/", nrow(mm)))
  log_msg("cross_species", "DONE")
}

if (!interactive()) main()
