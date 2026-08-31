#!/usr/bin/env Rscript

# Step 144: cluster one shared-orthobin WCG matrix with Seurat.
# Mapping-derived GM12878/PK15 labels are metadata added only after PCA,
# neighbors, clusters and UMAP have been constructed.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

PIPELINE_SIGNATURE <- "snmar_dna_ortholog_100kb_seurat_20260825"
VALID_METHODS <- c(
  "ratio_paper",
  "beta_posterior",
  "pearson_residual",
  "coverage_control"
)

parse_arguments <- function(values) {
  result <- list(
    project = "/share/org/YZWL/yzwl_wujx/Projects/2607_methylation",
    run_id = NA_character_,
    method = "ratio_paper",
    resolution = 0.2,
    neighbors = 30L,
    min_dist = 0.30,
    seed = 2607L,
    regress_depth = FALSE,
    force = FALSE
  )
  index <- 1L
  while (index <= length(values)) {
    key <- values[[index]]
    if (key %in% c("--force", "--regress-depth")) {
      result[[gsub("-", "_", sub("^--", "", key))]] <- TRUE
      index <- index + 1L
      next
    }
    if (index == length(values)) stop("missing value after ", key)
    value <- values[[index + 1L]]
    name <- gsub("-", "_", sub("^--", "", key))
    if (!name %in% names(result)) stop("unknown argument: ", key)
    result[[name]] <- value
    index <- index + 2L
  }
  result$resolution <- as.numeric(result$resolution)
  result$neighbors <- as.integer(result$neighbors)
  result$min_dist <- as.numeric(result$min_dist)
  result$seed <- as.integer(result$seed)
  if (
    is.na(result$run_id) || !nzchar(result$run_id) ||
    !result$method %in% VALID_METHODS ||
    result$resolution <= 0 || result$neighbors < 2L ||
    result$min_dist < 0 || result$min_dist > 1
  ) stop("invalid Seurat parameter")
  result
}

read_metrics <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop("missing/empty summary: ", path)
  }
  table <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(names(table), c("metric", "value"))) {
    stop("unexpected metric table header: ", path)
  }
  if (anyDuplicated(table$metric)) stop("duplicate metric: ", path)
  setNames(table$value, table$metric)
}

write_tsv_atomic <- function(table, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  write.table(
    table, temporary, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = "NA"
  )
  if (!file.rename(temporary, path)) stop("cannot replace output: ", path)
}

write_metrics <- function(values, path) {
  write_tsv_atomic(
    data.frame(metric = names(values), value = unname(values), check.names = FALSE),
    path
  )
}

read_gzip_tsv <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection))
  read.delim(connection, stringsAsFactors = FALSE, check.names = FALSE)
}

read_gzip_matrix_market <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection))
  as(readMM(connection), "dgCMatrix")
}

make_data_object <- function(matrix, cells) {
  placeholder <- Matrix(
    1, nrow = 1L, ncol = ncol(matrix), sparse = TRUE,
    dimnames = list("placeholder", colnames(matrix))
  )
  object <- CreateSeuratObject(
    counts = placeholder,
    assay = "placeholder",
    min.cells = 0,
    min.features = 0
  )
  object[["WCG"]] <- CreateAssayObject(data = matrix)
  DefaultAssay(object) <- "WCG"
  object <- AddMetaData(object, cells[colnames(object), , drop = FALSE])
  object
}

scatter_panel <- function(xy, values, title, palette = NULL, continuous = FALSE) {
  if (continuous) {
    colors <- hcl.colors(100, "Viridis")
    range_value <- range(values, finite = TRUE)
    if (diff(range_value) == 0) {
      point_colors <- colors[[50]]
    } else {
      bins <- pmax(1L, pmin(100L, as.integer(
        1 + 99 * (values - range_value[[1]]) / diff(range_value)
      )))
      point_colors <- colors[bins]
    }
    plot(xy, pch = 16, cex = 0.45, col = point_colors,
         xlab = "UMAP1", ylab = "UMAP2", main = title)
    return(invisible(NULL))
  }
  levels_value <- sort(unique(as.character(values)))
  if (is.null(palette)) {
    palette <- setNames(hcl.colors(length(levels_value), "Dark 3"), levels_value)
  }
  missing <- setdiff(levels_value, names(palette))
  if (length(missing)) {
    palette <- c(palette, setNames(hcl.colors(length(missing), "Dark 3"), missing))
  }
  plot(xy, pch = 16, cex = 0.45, col = palette[as.character(values)],
       xlab = "UMAP1", ylab = "UMAP2", main = title)
  legend("topright", legend = levels_value, col = palette[levels_value],
         pch = 16, cex = 0.65, bty = "n")
}

save_overview <- function(object, path, title) {
  xy <- Embeddings(object, "umap")
  metadata <- object@meta.data
  png(path, width = 2200, height = 1800, res = 180)
  old <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  on.exit({par(old); dev.off()}, add = TRUE)
  scatter_panel(xy, Idents(object), "Seurat clusters")
  scatter_panel(
    xy,
    metadata$DNA_mapping_derived_cell_line_posthoc,
    "mapping-derived label (post hoc)",
    c(GM12878 = "#2563EB", PK15 = "#DC2626", ambiguous = "#9CA3AF")
  )
  scatter_panel(
    xy,
    metadata$well_design,
    "physical-well design",
    c(GM12878_only_wells = "#2A9D8F", GM12878_PK15_mix_wells = "#BDBDBD")
  )
  scatter_panel(
    xy,
    log10(metadata$selected_shared_orthobin_WCG_covered_calls + 1),
    "log10(selected WCG calls + 1)",
    continuous = TRUE
  )
  mtext(title, side = 3, outer = TRUE, line = -1.2, cex = 1.2)
}

args <- parse_arguments(commandArgs(trailingOnly = TRUE))
project <- normalizePath(args$project, mustWork = TRUE)
input_root <- file.path(
  project, "10_downstream_analysis/11_dna_ortholog_100kb_seurat_input",
  args$run_id
)
input_summary_path <- file.path(
  project, "08_qc_summary/41_dna_ortholog_100kb_seurat_input",
  args$run_id, "dna_ortholog_100kb_seurat_input_summary.tsv"
)
input_summary <- read_metrics(input_summary_path)
if (
  input_summary[["status"]] != "validated" ||
  input_summary[["pipeline_signature"]] !=
    "snmar_dna_ortholog_100kb_seurat_input_20260825"
) stop("shared-orthobin Seurat input is not current/validated")

suffix <- if (args$regress_depth) "_depth_regressed" else ""
analysis_id <- paste0(args$run_id, "_", args$method, suffix)
output_root <- file.path(
  project, "10_downstream_analysis/12_dna_ortholog_100kb_seurat", analysis_id
)
qc_root <- file.path(
  project, "08_qc_summary/42_dna_ortholog_100kb_seurat", analysis_id
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)
summary_path <- file.path(qc_root, "dna_ortholog_100kb_seurat_summary.tsv")
if (file.exists(summary_path) && !args$force) {
  existing <- read_metrics(summary_path)
  if (
    existing[["status"]] == "validated" &&
    existing[["pipeline_signature"]] == PIPELINE_SIGNATURE
  ) {
    message("[SKIP] validated Seurat result exists: ", summary_path)
    quit(save = "no", status = 0)
  }
  stop("existing non-current result; inspect or rerun with --force: ", summary_path)
}

matrix_path <- input_summary[[paste0("matrix_", args$method)]]
features_path <- input_summary[["features"]]
cells_path <- input_summary[["cells"]]
for (path in c(matrix_path, features_path, cells_path)) {
  if (!file.exists(path) || file.info(path)$size == 0) stop("missing/empty input: ", path)
}
matrix <- read_gzip_matrix_market(matrix_path)
features <- read_gzip_tsv(features_path)
cells <- read_gzip_tsv(cells_path)
if (
  nrow(matrix) != nrow(features) || ncol(matrix) != nrow(cells) ||
  anyDuplicated(features$orthobin_id) || anyDuplicated(cells$cell_id)
) stop("matrix/feature/cell dimensions or identifiers are invalid")
rownames(matrix) <- features$orthobin_id
colnames(matrix) <- cells$cell_id
rownames(cells) <- cells$cell_id

set.seed(args$seed)
if (args$method == "ratio_paper") {
  object <- CreateSeuratObject(
    counts = matrix, assay = "WCG", min.cells = 0, min.features = 0,
    meta.data = cells[colnames(matrix), , drop = FALSE]
  )
  object <- NormalizeData(
    object, normalization.method = "LogNormalize", scale.factor = 10000,
    verbose = FALSE
  )
  object <- FindVariableFeatures(
    object, selection.method = "vst",
    nfeatures = min(2000L, nrow(object)), verbose = FALSE
  )
} else {
  object <- make_data_object(matrix, cells)
  VariableFeatures(object) <- rownames(matrix)
}

regressors <- NULL
if (args$regress_depth && args$method != "coverage_control") {
  object$log_selected_WCG_calls <- log1p(
    object$selected_shared_orthobin_WCG_covered_calls
  )
  regressors <- "log_selected_WCG_calls"
}
object <- ScaleData(
  object, features = VariableFeatures(object), vars.to.regress = regressors,
  verbose = FALSE
)
npcs <- min(30L, ncol(object) - 1L, length(VariableFeatures(object)) - 1L)
if (npcs < 2L) stop("too few dimensions for PCA")
object <- RunPCA(
  object, features = VariableFeatures(object), npcs = npcs,
  seed.use = args$seed, verbose = FALSE
)
dimensions <- seq_len(npcs)
object <- FindNeighbors(
  object, dims = dimensions,
  k.param = min(args$neighbors, ncol(object) - 1L), verbose = FALSE
)
object <- FindClusters(
  object, resolution = args$resolution, random.seed = args$seed,
  verbose = FALSE
)
object <- RunUMAP(
  object, dims = dimensions,
  n.neighbors = min(args$neighbors, ncol(object) - 1L),
  min.dist = args$min_dist, seed.use = args$seed, verbose = FALSE
)

coordinates <- as.data.frame(Embeddings(object, "umap"))
coordinates$cell_id <- rownames(coordinates)
coordinates$seurat_cluster <- as.character(Idents(object)[coordinates$cell_id])
metadata_columns <- setdiff(colnames(object@meta.data), "cell_id")
coordinates <- cbind(
  coordinates,
  object@meta.data[coordinates$cell_id, metadata_columns, drop = FALSE]
)
coordinates_path <- file.path(output_root, "dna_ortholog_100kb_seurat_coordinates.tsv")
write_tsv_atomic(coordinates, coordinates_path)

contingency <- as.data.frame.matrix(table(
  cluster = as.character(Idents(object)),
  mapping_label = object$DNA_mapping_derived_cell_line_posthoc
))
contingency$seurat_cluster <- rownames(contingency)
rownames(contingency) <- NULL
contingency_path <- file.path(output_root, "cluster_by_mapping_label.tsv")
write_tsv_atomic(contingency, contingency_path)

object_path <- file.path(output_root, "dna_ortholog_100kb_seurat.rds")
saveRDS(object, object_path, compress = "xz")
plot_path <- file.path(output_root, "dna_ortholog_100kb_seurat_overview.png")
save_overview(
  object, plot_path,
  paste("snMAR shared 100 kb WCG Seurat:", args$method)
)

labels <- object$DNA_mapping_derived_cell_line_posthoc
metrics <- c(
  status = "validated",
  pipeline_signature = PIPELINE_SIGNATURE,
  analysis_id = analysis_id,
  run_id = args$run_id,
  method = args$method,
  mode = "unsupervised_shared_orthobin_WCG_labels_posthoc_only",
  retained_cells = ncol(object),
  input_features = nrow(matrix),
  variable_features = length(VariableFeatures(object)),
  pca_components = npcs,
  seurat_clusters = length(levels(Idents(object))),
  resolution = args$resolution,
  neighbors = args$neighbors,
  min_dist = args$min_dist,
  regress_depth = args$regress_depth,
  GM12878_posthoc_cells = sum(labels == "GM12878"),
  PK15_posthoc_cells = sum(labels == "PK15"),
  ambiguous_posthoc_cells = sum(labels == "ambiguous"),
  posthoc_labels_used_before_embedding = "no",
  input_summary = input_summary_path,
  input_matrix = matrix_path,
  coordinates = coordinates_path,
  cluster_mapping_contingency = contingency_path,
  seurat_object = object_path,
  overview_png = plot_path
)
write_metrics(metrics, summary_path)
message(
  "[OK] Seurat ", analysis_id, ": cells=", ncol(object),
  " features=", nrow(matrix), " clusters=", length(levels(Idents(object)))
)
