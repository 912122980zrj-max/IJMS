#!/usr/bin/env Rscript
# 13_UUO_timecourse.R —— UUO 时间序列验证（GSE118339：normal/D3/D7/D14，TPM，n=3-4/组）
# 检验：人类机械核心基因（经 12 同源映射）在小鼠 UUO 进程中的时间轨迹——
#   样本级核心评分与天数 Spearman 趋势、D14 vs D0 wilcox、核心基因正向斜率占比。

suppressPackageStartupMessages({
  library(optparse)
  library(biomaRt)
  library(ggplot2)
  library(ggpubr)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "13_UUO_timecourse.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--indir"), type = "character",
              default = file.path(DATA_RAW, "GSE118339"), help = "TPM 文件目录"),
  make_option(c("--ortholog"), type = "character",
              default = file.path(DATA_PROC, "ortholog_map.csv"), help = "12 号同源映射表"),
  make_option(c("--outdir"), type = "character", default = DATA_PROC, help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

main <- function() {
  ensure_dirs()
  files <- list.files(opt$indir, pattern = "tpm.txt.gz$", full.names = TRUE)
  mats <- lapply(files, function(f) {
    d <- read.csv(gzfile(f), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    setNames(d[[2]], d[[1]])
  })
  genes <- sort(unique(unlist(lapply(mats, names))))
  M <- do.call(cbind, lapply(mats, function(m) m[genes]))
  colnames(M) <- sub("GSM\\d+_(.+)\\.tpm\\.txt\\.gz", "\\1", basename(files))
  day_map <- c("normal_1" = 0, "normal_2" = 0, "normal_3" = 0,
               "day3_4" = 3, "day3_5" = 3, "day3_6" = 3, "day3_7" = 3,
               "day7_8" = 7, "day7_9" = 7, "day7_10" = 7, "day7_11" = 7,
               "day14_12" = 14, "day14_13" = 14, "day14_14" = 14, "day14_15" = 14)
  day <- unname(day_map[colnames(M)])
  log_msg("tc", paste0("samples = ", ncol(M), " | genes = ", nrow(M)))

  # Ensembl -> 小鼠符号
  ensembl <- NULL
  for (mirror in c("useast", "asia", "www")) {
    ensembl <- tryCatch(useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl",
                                   mirror = mirror), error = function(e) NULL)
    if (!is.null(ensembl)) break
  }
  if (is.null(ensembl)) stop("Ensembl 镜像均不可达")
  ann <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
               filters = "ensembl_gene_id", values = genes, mart = ensembl)

  orth <- read.csv(opt$ortholog, stringsAsFactors = FALSE)
  mouse_names <- unique(orth$mmusculus_homolog_associated_gene_name)
  keep <- ann$external_gene_name %in% mouse_names
  row_sel <- rownames(M) %in% ann$ensembl_gene_id[keep]
  log_msg("tc", paste0("core mouse genes present = ", sum(row_sel)))

  zmat <- t(scale(t(log2(M[row_sel, , drop = FALSE] + 1))))
  score <- colMeans(zmat, na.rm = TRUE)
  df <- data.frame(day = day, score = score)
  rho <- cor(df$day, df$score, method = "spearman")
  p14 <- wilcox.test(score[day == 14], score[day == 0])$p.value
  log_msg("tc", paste0("spearman(score~day) = ", round(rho, 3),
                       " | D14 vs D0 wilcox p = ", signif(p14, 3)))

  slopes <- apply(zmat, 1, function(x) {
    if (all(is.na(x))) return(NA)
    coef(lm(x ~ day))[2]
  })
  log_msg("tc", paste0("core genes positive slope: ",
                       round(100 * mean(slopes > 0, na.rm = TRUE), 1), "% (",
                       sum(slopes > 0, na.rm = TRUE), "/", sum(!is.na(slopes)), ")"))

  safe_write(file.path(opt$outdir, "uuo_timecourse_matrix.rds"), function(p) saveRDS(M, p))
  safe_write(file.path(opt$outdir, "uuo_timecourse_results.csv"), function(p) {
    write_csv_utf8(data.frame(metric = c("spearman_day", "D14_vs_D0_wilcox_p",
                                         "positive_slope_frac"),
                              value = c(rho, p14, mean(slopes > 0, na.rm = TRUE))), p)
  })

  p <- ggplot(df, aes(factor(day), score, color = factor(day))) +
    geom_boxplot(outlier.size = 0.5, width = 0.5) +
    geom_jitter(width = 0.12, size = 0.8) +
    scale_color_manual(values = CB_PALETTE, guide = "none") +
    labs(x = "Days post-UUO", y = "Mechanical core score",
         title = "UUO time course: mechanical core trajectory") +
    annotate("text", x = 1, y = max(df$score) * 0.95, size = 2.5, hjust = 0,
             label = paste0("spearman rho = ", round(rho, 3),
                            " ; D14 vs D0 p = ", signif(p14, 2))) +
    pub_theme()
  save_pub_fig(p, file.path(FIG_DIR, "fig_uuo_timecourse"))
  log_msg("tc", "DONE")
}

if (!interactive()) main()
