.libPaths("~/.guix-profile/site-library")

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(SeuratDisk)
library(SingleCellExperiment)

# Set options for large datasets and numerical output
options(future.globals.maxSize = 3e+09)
options(scipen = 999)

args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
ref <- as.numeric(args[1])
method <- args[2]
n_samples <- as.numeric(args[3])
rep <- as.numeric(args[4])
seed <- as.numeric(args[5])

# Set the seed for reproducibility
set.seed(seed)

# Define your base path
base_path <- "/fast/AG_Haghverdi/Ehsan_Karimiara/facs_sampling/lcmv/benchmark/"

# Initialize variables
base_seed <- 1000

# Set the global seed for reproducibility
set.seed(base_seed)

# Path construction
current_path <- paste0(base_path, ref, "/")
input_csv <- paste0(current_path, "adata_x.csv")
input_csv_meta <- paste0(current_path, "obs.csv")
output_csv_prefix <- paste0(current_path, method, "/", n_samples, "/", rep, "/")

# Read data and metadata
message("Reading counts...")
x <- read.csv(input_csv, header = TRUE, row.names = 1)

message("Reading metadata...")
m <- read.csv(input_csv_meta, header = TRUE, row.names = 1)
colnames(m)[1] <- "sample"

# Create Seurat object
message("Creating Seurat object...")
obj <- CreateSeuratObject(counts = t(x), meta.data = m, project = "seurat", min.cells = 0, min.features = 0)
obj[["RNA"]]$data <- obj[["RNA"]]$counts

# Find variable features
obj <- FindVariableFeatures(obj, verbose = TRUE)

# Add batch metadata for reproducible splitting
obj$batch <- sample(x = 1:50, size = ncol(obj), replace = TRUE)
obj.list <- SplitObject(obj, split.by = "batch")

# Apply SketchData to each batch with a batch-specific seed
sketched_list <- list()

for (i in seq_along(obj.list)) {
  batch_specific_seed <- seed + i
  message(paste("Sketching batch", i, "with seed", batch_specific_seed))
  obj.list[[i]] <- SketchData(object = obj.list[[i]], ncells = ceiling(n_samples / length(obj.list)), method = "LeverageScore", sketched.assay = "sketch", seed = batch_specific_seed)
  cells_data <- obj.list[[i]][["sketch"]]@cells@.Data[, 1]
  cell_numbers <- names(cells_data)
  
  # Save sampled indices for each batch
  output_csv <- paste0(output_csv_prefix, "results_batch_", i, ".csv")
  write.csv(cell_numbers, file = output_csv, row.names = FALSE)
}

message(paste("Finished processing", ref, "with n_samples:", n_samples, "repetition", rep))
