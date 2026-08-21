#!/usr/bin/env Rscript
# 16_immune28_charoentong.R —— Charoentong 28 免疫细胞反卷积（ssGSEA）
# 基因集来源：用户提供 Charoentong 2017 Cell Reports Table S6（mmc3.xlsx），
#   已转换为 data/reference/Charoentong_28_immune.gmt（28 细胞类型、782 基因）。
# 输出：IgAN vs 对照 28 细胞评分差异（wilcox+BH）、C1/C2 差异、28 细胞热图。

suppressPackageStartupMessages({
  library(optparse)
  library(GSVA)
  library(ggpubr)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "16_immune28_charoentong.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--expr"), type = "character",
              default = file.path(DATA_PROC, "expr_human_main.rds"), help = "主队列表达矩阵"),
  make_option(c("--meta"), type = "character",
              default = file.path(DATA_PROC, "meta_human_main.csv"), help = "主队列表型"),
  make_option(c("--subtype"), type = "character",
              default = file.path(DATA_PROC, "subtype_assignment.csv"), help = "08 分型"),
  make_option(c("--gmt"), type = "character",
              default = file.path(DATA_REF, "Charoentong_28_immune.gmt"), help = "28 细胞 GMT"),
  make_option(c("--outdir"), type = "character", default = DATA_PROC, help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

read_gmt <- function(path) {
  lines <- strsplit(readLines(path), "\t")
  sets <- lapply(lines, function(x) x[-(1:2)])
  names(sets) <- vapply(lines, `[[`, character(1), 1)
  sets
}

main <- function() {
  ensure_dirs()
  expr <- readRDS(opt$expr)
  meta <- read.csv(opt$meta, stringsAsFactors = FALSE)
  meta$disease <- ifelse(meta$status == "IgAN patient", "IgAN", "control")
  sets <- lapply(read_gmt(opt$gmt), function(g) intersect(g, rownames(expr)))
  sets <- sets[lengths(sets) >= 5]
  log_msg("immune28", paste0("cell types = ", length(sets),
                             " | genes used = ", sum(lengths(sets))))

  score <- gsva(ssgseaParam(exprData = as.matrix(expr), geneSets = sets))
  safe_write(file.path(opt$outdir, "immune28_charoentong_scores.csv"), function(p) {
    out <- data.frame(sample = colnames(score), t(score), check.names = FALSE)
    write_csv_utf8(out, p)
  })

  res <- do.call(rbind, lapply(rownames(score), function(nm) {
    x <- as.numeric(score[nm, ])
    data.frame(celltype = nm,
               IgAN = mean(x[meta$disease == "IgAN"]),
               control = mean(x[meta$disease == "control"]),
               p = wilcox.test(x[meta$disease == "IgAN"], x[meta$disease == "control"])$p.value)
  }))
  res$padj <- p.adjust(res$p, method = "BH")
  res <- res[order(res$p), ]
  safe_write(file.path(opt$outdir, "immune28_charoentong_results.csv"), function(p) {
    write_csv_utf8(res, p)
  })
  log_msg("immune28", paste0("top5: ", paste(head(res$celltype, 5), collapse = ", ")))
  log_msg("immune28", paste0("significant (BH<0.05) n=", sum(res$padj < 0.05)))

  sub <- read.csv(opt$subtype, stringsAsFactors = FALSE)
  common <- intersect(sub$sample, colnames(expr))
  top2 <- head(res$celltype, 2)
  for (ct in top2) {
    df <- data.frame(subtype = factor(sub$subtype[match(common, sub$sample)]),
                     score = as.numeric(score[ct, common]))
    p <- ggboxplot(df, x = "subtype", y = "score", color = "subtype",
                   palette = CB_PALETTE, add = "jitter", xlab = NULL, ylab = ct) +
      stat_compare_means(method = "wilcox.test") + pub_theme()
    save_pub_fig(p, file.path(FIG_DIR, paste0("fig_subtype_", gsub("[ /]", "_", ct))))
  }

  hm <- sapply(rownames(score), function(nm) {
    c(IgAN = mean(as.numeric(score[nm, meta$disease == "IgAN"])),
      control = mean(as.numeric(score[nm, meta$disease == "control"])))
  })
  hp <- file.path(FIG_DIR, "fig_immune28_heatmap.pdf")
  hp_png <- file.path(FIG_DIR, "fig_immune28_heatmap.png")
  for (f in c(hp, hp_png)) if (file.exists(f)) unlink(f)
  pdf(hp, width = 4.5, height = 9)
  pheatmap(hm, color = colorRampPalette(c("#0072B2", "white", "#E69F00"))(50),
           fontsize = 6, cluster_cols = FALSE)
  dev.off()
  png(hp_png, width = 4.5, height = 9, units = "in", res = 300)
  pheatmap(hm, color = colorRampPalette(c("#0072B2", "white", "#E69F00"))(50),
           fontsize = 6, cluster_cols = FALSE)
  dev.off()
  log_msg("immune28", "DONE")
}

if (!interactive()) main()
