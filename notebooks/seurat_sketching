.libPaths("~/.guix-profile/site-library")

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(SeuratDisk)
library(SingleCellExperiment)

# set this option when analyzing large datasets
options(future.globals.maxSize = 3e+09)

# args <- commandArgs(trailingOnly = TRUE)
# 
# # Assigning arguments to variables
# input_rds <- args[1]  # Input RDS file name, e.g., '80_4_sample.rds'
# output_csv <- args[2]  # Output CSV file name
# 
# # Validate arguments
# if(length(args) < 2) {
#   stop("Not enough arguments. Usage: Rscript script.R inputRDS outputCSV", call. = FALSE)
# }

base_path <- "/fast/AG_Haghverdi/Ehsan_Karimiara/facs_sampling/sara_data/batches/"
unique_id <- "80_4"
n_samples <- 40000

base_path <- paste0(base_path, unique_id, "/")

input_csv <- paste0(base_path, "raw_counts.csv")
input_csv_meta <- paste0(base_path, "metadata.csv")
output_rds <- paste0(base_path, "sample.rds")
output_csv <- paste0(base_path, "atomic_indices_", n_samples, ".csv")

message("Reading counts...")

x <- read.csv(input_csv, header=TRUE)
rownames(x) <- x[,1]
x[,1] <- NULL
print(dim(x))
print(x[1:5,1:5])

message("Reading metadata...")
m <- read.csv(input_csv_meta, header=TRUE)
rownames(m) <- m[,1]
colnames(m)[1] <- "sample"
print(dim(m))
print(head(m))


message("Writing seurat object...")
saveRDS(
  CreateSeuratObject(counts=t(x),meta.data=m,project="seurat",min.cells=0,min.features=0),
  output_rds
)


# Read .h5ad file
input_address <- 
obj <- NULL
obj <- readRDS(output_rds)

obj[["RNA"]]$data <- obj[["RNA"]]$counts
obj <- FindVariableFeatures(obj, verbose = T)

start.time <- Sys.time()

obj <- SketchData( object = obj, ncells = n_samples, method = "LeverageScore", 
                   sketched.assay = "sketch")

end.time <- Sys.time()
time.taken <- as.numeric(end.time - start.time, units = "secs")
time.taken

cells_data <- obj[["sketch"]]@cells@.Data[,1]
# Assuming cells_data contains your data
# Extract only the cell numbers
cell_numbers <- names(cells_data)


# Alternatively, save to a CSV file
write.csv(cell_numbers, file = output_csv, row.names = FALSE)

