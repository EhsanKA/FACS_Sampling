.libPaths("~/.guix-profile/site-library")

library(Seurat)
library(BPCells)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(SeuratDisk)
library(SingleCellExperiment)

# set this option when analyzing large datasets
options(future.globals.maxSize = 100e+09)
options(scipen = 999)


args <- commandArgs(trailingOnly = TRUE)

ref <- as.numeric(args[1])
method <- args[2]
n_samples <- as.numeric(args[3])
rep <- as.numeric(args[4])
seed <- as.numeric(args[5])

# ref <- 1
# method <- 'atomic'
# n_samples <- 50000
# rep <- 0
# seed <- 25636

# Define your base path
base_path <- "/fast/AG_Ohler/ekarimi/projects/FACS_Sampling/data/mcc/benchmark/"

# Base seed for reproducibility
base_seed <- 1000
runtimes <- list()
set.seed(seed)

current_path <- paste0(base_path, ref, "/")

input_csv <- paste0(current_path, "adata_x.csv")
input_csv_meta <- paste0(current_path, "obs.csv")
output_csv_prefix <- paste0(current_path, method, "/", n_samples, "/", rep, "/")
output_rds <- paste0(current_path, "object.rds")
output_matrix_count <- paste0(current_path, "counts")


if (file.exists(output_rds)) {
  message("Loading existing Seurat object...")
  obj <- readRDS(output_rds)
  obj[["RNA"]]$counts <- open_matrix_dir(dir = output_matrix_count)#, buffer_size = 1048576L)

} else {
  
  message("Reading counts...")
  x <- read.csv(input_csv, header = TRUE)
  rownames(x) <- x[, 1]
  x[, 1] <- NULL
  
  message("Reading metadata...")
  m <- read.csv(input_csv_meta, header = TRUE)
  rownames(m) <- m[, 1]
  colnames(m)[1] <- "sample"
  
  message("Creating Seurat object...")
  obj <- CreateSeuratObject(counts = t(x), meta.data = m, project = "seurat",
                            min.cells = 0, min.features = 0)
  rm(x)
  rm(m)
  gc()

  # # Write the counts layer to a directory
  write_matrix_dir(mat = obj[["RNA"]]$counts, dir = output_matrix_count)#, buffer_size = 1048576L)
  obj[["RNA"]]$counts <- open_matrix_dir(dir = output_matrix_count)#, buffer_size = 1048576L)
  saveRDS(obj, file = output_rds)
}

# obj[["RNA"]]$counts <- as.matrix(obj[["RNA"]]$counts)
obj[["RNA"]]$data <- obj[["RNA"]]$counts

obj <- FindVariableFeatures(obj, verbose = TRUE)

message(paste("Starting SketchData for ref:", ref, "n_samples:", n_samples, "rep:", rep))
start.time <- Sys.time()

obj <- SketchData(object = obj, ncells = n_samples, method = "LeverageScore", sketched.assay = "sketch", seed = seed)

end.time <- Sys.time()
time.taken <- as.numeric(end.time - start.time, units = "secs")
message(paste("time taken is: ", time.taken))
runtimes[["runtime"]] <- time.taken

cells_data <- obj[["sketch"]]@cells@.Data[, 1]
cell_numbers <- names(cells_data)

# Save sampled indices and runtimes for each repetition
output_csv <- paste0(output_csv_prefix, "results.csv")
write.csv(cell_numbers, file = output_csv, row.names = FALSE)

message(paste("Finished processing", ref, "with n_samples:", n_samples, "repetition", rep, "in", time.taken, "seconds"))

# Save the runtimes and results to files
runtime_output <- paste0(output_csv_prefix, "runtimes.csv")
write.csv(runtimes, file = runtime_output, row.names = TRUE)



