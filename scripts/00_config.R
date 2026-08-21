# 00_config.R —— 选题三全局配置、日志、安全写入与发表级绘图规范
# 沿用 IPF 模板（template/repro_IPF_MRG/scripts/00_config.R）的工程约定。

suppressPackageStartupMessages({
  library(ggplot2)
})

raw_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", raw_args[grep("^--file=", raw_args)])
if (length(file_arg) == 0L) file_arg <- "00_config.R"
SCRIPT_DIR <- normalizePath(dirname(file_arg))
PROJECT_ROOT <- normalizePath(dirname(SCRIPT_DIR))

DATA_RAW  <- file.path(PROJECT_ROOT, "data", "raw")
DATA_PROC <- file.path(PROJECT_ROOT, "data", "processed")
DATA_REF  <- file.path(PROJECT_ROOT, "data", "reference")
RES_DIR   <- file.path(PROJECT_ROOT, "results")
FIG_DIR   <- file.path(PROJECT_ROOT, "figures")

set.seed(42)  # 全局复现种子（论文未报告，项目约定）

ensure_dirs <- function() {
  for (d in c(DATA_RAW, DATA_PROC, DATA_REF, RES_DIR, FIG_DIR)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  invisible(TRUE)
}

log_msg <- function(stage, msg) {
  message(sprintf("[%s] %s | %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), stage, msg))
}

safe_write <- function(path, writer) {
  if (file.exists(path)) stop("Refusing to overwrite existing file: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writer(path)
  invisible(path)
}

write_csv_utf8 <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

CB_PALETTE <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00")

pub_theme <- function(base_size = 6, single_column = TRUE) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title   = element_text(size = 7, face = "bold", hjust = 0),
      axis.title   = element_text(size = 6),
      axis.text    = element_text(size = 5),
      legend.title = element_text(size = 6),
      legend.text  = element_text(size = 5),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
}

save_pub_fig <- function(p, stem, single_column = TRUE, dpi = 300,
                         width_mm = if (single_column) 89 else 183,
                         height_mm = if (single_column) 80 else 120) {
  # 图为可再生派生物：允许覆盖；数据文件仍由 safe_write 防覆盖
  for (path in c(paste0(stem, ".pdf"), paste0(stem, ".png"))) {
    if (file.exists(path)) unlink(path)
  }
  ggsave(paste0(stem, ".pdf"), plot = p, width = width_mm, height = height_mm, units = "mm",
         device = grDevices::cairo_pdf)
  ggsave(paste0(stem, ".png"), plot = p, width = width_mm, height = height_mm, units = "mm",
         dpi = dpi)
  invisible(stem)
}

log_msg("config", paste("PROJECT_ROOT =", PROJECT_ROOT))
