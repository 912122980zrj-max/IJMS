#!/usr/bin/env Rscript
# 07_nomogram.R —— CKD 风险 nomogram（真实实现；终点 = CKD vs normal）
# 逐步 logistic（AIC 双向）→ rms::lrm → nomogram（CKD 概率）→ val.prob 校准。

suppressPackageStartupMessages({
  library(optparse)
  library(rms)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "07_nomogram.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--expr"), type = "character",
              default = file.path(DATA_PROC, "expr_human_main.rds"), help = "GSE115857 表达矩阵"),
  make_option(c("--meta"), type = "character",
              default = file.path(DATA_PROC, "meta_human_main.csv"), help = "GSE115857 表型"),
  make_option(c("--sig"), type = "character",
              default = file.path(DATA_PROC, "progression_signature.rds"), help = "06 输出 signature"),
  make_option(c("--outdir"), type = "character", default = DATA_PROC, help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

main <- function() {
  ensure_dirs()
  expr <- readRDS(opt$expr)
  meta <- read.csv(opt$meta, stringsAsFactors = FALSE)
  sig <- readRDS(opt$sig)
  meta <- meta[meta$status %in% c("IgAN patient", "minimal change disease",
                                  "Membranous glomerulonephritis", "Living donor"), ]
  expr <- expr[, meta$sample]
  dat <- data.frame(IgAN = as.integer(meta$status == "IgAN patient"),
                    t(expr[intersect(sig, rownames(expr)), , drop = FALSE]))

  full <- glm(IgAN ~ ., data = dat, family = binomial())
  step <- MASS::stepAIC(full, direction = "both", trace = 0)
  log_msg("glm", paste0("stepwise keeps: ", paste(names(coef(step))[-1], collapse = ", ")))
  safe_write(file.path(opt$outdir, "stepwise_glm.rds"), function(p) saveRDS(step, p))

  assign("dd", datadist(dat), envir = .GlobalEnv)
  options(datadist = "dd")
  lrm_fit <- lrm(IgAN ~ ., data = dat, x = TRUE, y = TRUE)
  nomo <- nomogram(lrm_fit, fun = plogis, funlabel = "IgAN probability")
  pdf_path <- file.path(FIG_DIR, "fig_nomogram_igag.pdf")
  png_path <- file.path(FIG_DIR, "fig_nomogram_igag.png")
  for (p in c(pdf_path, png_path)) if (file.exists(p)) unlink(p)
  pdf(pdf_path, width = 6.5, height = 5)
    plot(nomo)
  dev.off()
  png(png_path, width = 6.5, height = 5, units = "in", res = 300)
  plot(nomo)
  dev.off()

  pred <- predict(lrm_fit, type = "fitted")
  cal <- val.prob(pred, dat$IgAN, pl = FALSE)
  log_msg("calib", paste0("C-index (Dxy) = ", round(cal["Dxy"], 3),
                          " | slope = ", round(cal["Slope"], 3)))
  log_msg("nomogram", "DONE")
}

if (!interactive()) main()
