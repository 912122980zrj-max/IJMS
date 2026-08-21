#!/usr/bin/env Rscript
# 10_scRNA_UUO.R —— 小鼠 UUO 单细胞定位（真实实现，GSE175412：1 sham + 2 UUO）
# QC（mt<20%、nFeature>200）→ 归一/HVG/PCA/聚类/UMAP → 标记基因模块打分注释细胞类型 →
# 机械核心基因（人源经 12 号同源映射）模块打分 → sham vs UUO 比例与机械评分比较。

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0L) script_path <- "10_scRNA_UUO.R"
source(file.path(dirname(normalizePath(script_path)), "00_config.R"))

option_list <- list(
  make_option(c("--indir"), type = "character",
              default = file.path(DATA_RAW, "GSE175412", "RAW"), help = "每样本一个 10X 目录"),
  make_option(c("--ortholog"), type = "character",
              default = file.path(DATA_PROC, "ortholog_map.csv"), help = "12 号输出的人鼠同源表"),
  make_option(c("--outdir"), type = "character", default = DATA_PROC, help = "输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

main <- function() {
  ensure_dirs()
  dirs <- list.dirs(opt$indir, recursive = FALSE)
  objs <- lapply(dirs, function(d) {
    CreateSeuratObject(counts = Read10X(d), project = basename(d),
                       min.cells = 3, min.features = 200)
  })
  obj <- merge(objs[[1]], y = objs[-1])
  obj$group <- ifelse(grepl("^CS", obj$orig.ident), "sham", "UUO")
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj <- subset(obj, subset = percent.mt < 20 & nFeature_RNA > 200)
  log_msg("qc", paste0("cells after QC = ", ncol(obj)))

  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, features = VariableFeatures(obj))
  obj <- FindNeighbors(obj, dims = 1:20)
  obj <- FindClusters(obj, resolution = 0.8)
  obj <- RunUMAP(obj, dims = 1:20, seed.use = 42)

  markers <- list(
    PT    = c("Slc34a1", "Slc22a6", "Slc22a8", "Lrp2"),
    TAL   = c("Slc12a1", "Umod"),
    DCT   = c("Slc12a3", "Pvalb"),
    CD    = c("Aqp2", "Aqp3"),
    Podo  = c("Nphs1", "Nphs2"),
    Endo  = c("Pecam1", "Emcn"),
    Macro = c("Adgre1", "Cd68", "Lyz2"),
    Tcell = c("Cd3d", "Cd3e", "Nkg7"),
    Bcell = c("Ms4a1", "Cd79a"),
    Fibro = c("Pdgfrb", "Col1a1", "Acta2")
  )
  for (nm in names(markers)) {
    obj <- AddModuleScore(obj, features = list(markers[[nm]]), name = paste0("sc_", nm))
  }
  sc_cols <- paste0("sc_", names(markers), "1")
  sc_mat <- FetchData(obj, vars = sc_cols)
  obj$celltype <- apply(sc_mat, 1, function(x) names(markers)[which.max(x)])

  orth <- read.csv(opt$ortholog, stringsAsFactors = FALSE)
  mrg_mouse <- intersect(orth$mmusculus_homolog_associated_gene_name, rownames(obj))
  obj <- AddModuleScore(obj, features = list(mrg_mouse), name = "MRG")
  log_msg("mrg", paste0("mouse MRG genes used = ", length(mrg_mouse)))

  safe_write(file.path(opt$outdir, "uuo_seurat.rds"), function(p) saveRDS(obj, p))
  safe_write(file.path(opt$outdir, "uuo_meta.csv"), function(p) {
    write_csv_utf8(obj@meta.data, p)
  })
  tab <- prop.table(table(obj$group, obj$celltype), margin = 1)
  safe_write(file.path(opt$outdir, "uuo_celltype_proportions.csv"), function(p) {
    write_csv_utf8(as.data.frame.matrix(tab), p)
  })

  p1 <- DimPlot(obj, group.by = "celltype", label = TRUE, repel = TRUE) +
    scale_color_manual(values = rep(CB_PALETTE, 2)) + pub_theme()
  save_pub_fig(p1, file.path(FIG_DIR, "fig_uuo_umap_celltype"))

  df <- obj@meta.data
  p2 <- ggplot(df, aes(group, MRG1, fill = group)) +
    geom_violin(scale = "width") + geom_boxplot(width = 0.15, outlier.size = 0.2) +
    scale_fill_manual(values = CB_PALETTE) +
    labs(x = NULL, y = "Mechanical core module score") +
    stat_compare_means(method = "wilcox.test") + pub_theme()
  save_pub_fig(p2, file.path(FIG_DIR, "fig_uuo_mrg_score"))
  log_msg("scRNA", paste0("MRG score sham vs UUO wilcox p=",
                          signif(wilcox.test(MRG1 ~ group, data = df)$p.value, 3)))
  log_msg("scRNA", "DONE")
}

if (!interactive()) main()
