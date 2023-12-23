.libPaths("~/.guix-profile/site-library")
library(Seurat)

message("Reading counts...")

x <- read.csv("/fast/AG_Haghverdi/Ehsan_Karimiara/facs_sampling/sara_data/adata_ref_sara_2M_raw_counts.csv",header=TRUE)
rownames(x) <- x[,1]
x[,1] <- NULL
print(dim(x))
print(x[1:5,1:5])

message("Reading metadata...")
m <- read.csv("/fast/AG_Haghverdi/Ehsan_Karimiara/facs_sampling/sara_data/adata_ref_sara_2M_metadata.csv",header=TRUE)
rownames(m) <- m[,1]
colnames(m)[1] <- "sample"
print(dim(m))
print(head(m))


message("Writing seurat object...")
saveRDS(
  CreateSeuratObject(counts=t(x),meta.data=m,project="seurat",min.cells=0,min.features=0),
  "/fast/AG_Haghverdi/Ehsan_Karimiara/facs_sampling/sara_data/adata_ref_sara_2M.rds"
)
