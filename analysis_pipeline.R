############################################################
## Spatial (Visium HD) workflow for OsiR project
## Steps: RCTD -> tumor-dominant spots -> Harmony -> clustering
##        -> markers -> paired DE -> KEGG enrichment
############################################################

suppressPackageStartupMessages({
  library(arrow)
  library(spacexr)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(Seurat)
  library(harmony)            # Harmony backend (used by Seurat HarmonyIntegration)
  library(org.Hs.eg.db)
  library(clusterProfiler)
})

set.seed(1234)

# -----------------------------
# 0) Config (EDIT HERE ONLY)
# -----------------------------
cfg <- list(
  # inputs
  seurat_rctd_pp   = "./seurat_RCTD_pp.rds",
  seurat_rctd_all  = "./seurat_RCTD_all.rds",
  sc_ref_rds       = "./GSE131907_Lung_Cancer_raw_UMI_matrix.rds",
  tissue_positions = "./tissue_positions.parquet",
  
  # RCTD params
  rctd_mode  = "doublet",
  max_cores  = 10,
  
  # tumor-dominant spot selection
  tumor_ct   = c("tS1", "tS2", "tS3", "Malignant-cells"),
  weight_col = "top1_weight",
  ct_col     = "top1_ct",
  
  # Harmony / clustering
  batch_var  = "sampleid",
  dims_use   = 1:10,
  resolution = 0.3,
  
  # outputs
  out_dir    = "./result",
  rds_rctd   = "RCTD_result.rds",
  rds_harmony= "0.3cluster_after_harmony.rds"
)

dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(
  file.exists(cfg$seurat_rctd_pp),
  file.exists(cfg$seurat_rctd_all),
  file.exists(cfg$sc_ref_rds),
  file.exists(cfg$tissue_positions)
)

# Small helper
write_csv2 <- function(x, filename) {
  write.csv(x, file = file.path(cfg$out_dir, filename), row.names = TRUE)
}

# -----------------------------
# 1) Load spatial Seurat object
# -----------------------------
seurat_pp <- readRDS(cfg$seurat_rctd_pp)

# -----------------------------
# 2) Build reference for RCTD (scRNA-seq)
# -----------------------------
ref <- readRDS(cfg$sc_ref_rds)
sc_counts <- ref@assays$RNA@counts
sc_cell_types <- ref@meta.data$cell_type

reference_se <- SummarizedExperiment(
  assays  = list(counts = sc_counts),
  colData = data.frame(cell_type = sc_cell_types)
)

rm(ref, sc_counts, sc_cell_types)

# -----------------------------
# 3) Build SpatialExperiment for RCTD
# -----------------------------
sp_counts <- GetAssayData(seurat_pp, assay = "Spatial", layer = "counts")

pos <- as.data.frame(read_parquet(cfg$tissue_positions))
pos <- pos[match(colnames(sp_counts), pos$barcode), c("pxl_col_in_fullres", "pxl_row_in_fullres")]
rownames(pos) <- colnames(sp_counts)
colnames(pos) <- c("x", "y")

sp_spe <- SpatialExperiment(
  assays        = list(counts = sp_counts),
  spatialCoords = as.matrix(pos)
)

rm(pos)

# -----------------------------
# 4) Run RCTD
# -----------------------------
rctd_obj <- createRctd(sp_spe, reference_se, pixel_count_min = 0, UMI_min = 0)
rctd_res <- runRctd(rctd_obj, rctd_mode = cfg$rctd_mode, max_cores = cfg$max_cores)

w <- rctd_res@results$weights

top1_ct <- rownames(w)[apply(w, 2, which.max)]
names(top1_ct) <- colnames(seurat_pp)

top1_weight <- apply(w, 2, max)
names(top1_weight) <- colnames(seurat_pp)

seurat_pp <- AddMetaData(seurat_pp, top1_weight, col.name = cfg$weight_col)
seurat_pp <- AddMetaData(seurat_pp, top1_ct,     col.name = cfg$ct_col)

saveRDS(seurat_pp, file = cfg$rds_rctd)

# -----------------------------
# 5) Load merged spatial object + keep tumor-dominant spots
# -----------------------------
seurat_all <- readRDS(cfg$seurat_rctd_all)
DefaultAssay(seurat_all) <- "Spatial"

keep <- seurat_all[[cfg$ct_col]][, 1] %in% cfg$tumor_ct
seurat_all <- seurat_all[, keep]

# Scale spot counts by tumor confidence score (top1_weight)
ca_score <- seurat_all[[cfg$weight_col]][, 1]
counts <- GetAssayData(seurat_all, assay = "Spatial", layer = "counts")
counts <- sweep(counts, MARGIN = 2, STATS = ca_score, FUN = "*")
seurat_all <- SetAssayData(seurat_all, assay = "Spatial", layer = "counts", new.data = counts)

rm(counts, keep, ca_score)

# -----------------------------
# 6) Seurat standard workflow: Normalize -> HVG -> Scale -> PCA
# -----------------------------
seurat_all <- NormalizeData(seurat_all)
seurat_all <- FindVariableFeatures(seurat_all, selection.method = "vst", nfeatures = 3000)
seurat_all <- ScaleData(seurat_all)
seurat_all <- RunPCA(seurat_all, features = VariableFeatures(seurat_all))
seurat_all <- RunUMAP(seurat_all, dims = 1:30)

# -----------------------------
# 7) Harmony batch correction (Seurat v5 compatible)
# -----------------------------
seurat_all[["Spatial"]] <- split(seurat_all[["Spatial"]], f = seurat_all[[cfg$batch_var]][, 1])

seurat_all <- IntegrateLayers(
  object         = seurat_all,
  method         = HarmonyIntegration,
  orig.reduction = "umap",
  new.reduction  = "harmony",
  verbose        = FALSE
)

# Downstream graph/cluster on Harmony reduction
seurat_all <- FindNeighbors(seurat_all, reduction = "harmony", dims = cfg$dims_use)
seurat_all <- FindClusters(seurat_all, resolution = cfg$resolution)

saveRDS(seurat_all, file = cfg$rds_harmony)

# -----------------------------
# 8) Cluster markers (all clusters)
# -----------------------------
all_markers <- FindAllMarkers(seurat_all, group.by = "Spatial_snn_res.0.3", logfc.threshold = -1e10)
write_csv2(all_markers, "all_markers.csv")

# -----------------------------
# 9) "Paired" subset DE
# -----------------------------
seurat_sub <- seurat_all
seurat_sub@images <- list()

# keep selected sample pair
keep_pair <- seurat_sub$sampleid %in% c("S1320641", "S1492809")
seurat_sub <- seurat_sub[, keep_pair]

# keep selected clusters
keep_clu <- seurat_sub$Spatial_snn_res.0.3 %in% c("0", "1", "3", "5", "9")
seurat_sub <- seurat_sub[, keep_clu]

Idents(seurat_sub) <- "Spatial_snn_res.0.3"

markers <- FindMarkers(
  seurat_sub,
  ident.1 = "1", ident.2 = "0",
  min.pct = 0.25,
  logfc.threshold = -1e-8
)
write_csv2(markers, "markers_S1320641_S1492809_1_vs_0.csv")

# -----------------------------
# 10) KEGG enrichment on up markers + add SYMBOL gene list
# -----------------------------
markers_up <- markers[markers$avg_log2FC > 0.25, , drop = FALSE]
cluster_genes <- rownames(markers_up)

gene_entrez <- mapIds(
  org.Hs.eg.db,
  keys     = cluster_genes,
  keytype  = "SYMBOL",
  column   = "ENTREZID"
)
gene_entrez <- unique(na.omit(as.character(gene_entrez)))

kegg <- enrichKEGG(
  gene     = gene_entrez,
  organism = "hsa",
  keyType  = "ncbi-geneid"   # Entrez IDs
)

# Convert KEGG result geneID to SYMBOL
if (!is.null(kegg) && nrow(kegg@result) > 0) {
  to_symbol <- function(geneID_string) {
    ids <- strsplit(geneID_string, "/", fixed = TRUE)[[1]]
    sym <- mapIds(org.Hs.eg.db, keys = ids, keytype = "ENTREZID", column = "SYMBOL")
    sym <- unique(na.omit(as.character(sym)))
    paste(sym, collapse = "/")
  }
  kegg@result$geneID_symbol <- vapply(kegg@result$geneID, to_symbol, FUN.VALUE = character(1))
}

write_csv2(kegg@result, "kegg_S1320641_S1492809_1_vs_0_up_markers.csv")

# -----------------------------
# 11) Session info (good practice for GitHub)
# -----------------------------
writeLines(capture.output(sessionInfo()), con = file.path(cfg$out_dir, "sessionInfo.txt"))
