# Downloading and preparing the Human HEAD5 file for music and remove redunant genes from the counts files. 
```
# Install necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("SingleCellExperiment", "zellkonverter", "HDF5Array"))

# Load required libraries
library(SingleCellExperiment)
library(zellkonverter)
library(HDF5Array)


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("HDF5Array")

# Now try loading the library
library(HDF5Array)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("zellkonverter")
library(zellkonverter)

print("HDF5Array installation and loading complete.")
# Download the file - note it did not work file too big
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE185224&format=file&file=GSE185224%5Fclustered%5Fannotated%5Fadata%5Fk10%5Flr0%2E92%5Fv1%2E7%2Eh5ad%2Egz"
destfile <- "GSE185224_clustered_annotated_adata_k10_lr0.92_v1.7.h5ad.gz"
download.file(url, destfile)

# Decompress the file
R.utils::gunzip(destfile)

# Read the h5ad file and convert to SingleCellExperiment
sce <- readH5AD("GSE185224_clustered_annotated_adata_k10_lr0.92_v1.7.h5ad")

# Print summary of the SingleCellExperiment object
print(sce)

# Display the first few cell type annotations
print(head(sce$cell_type))

# Display the first few gene names
print(head(rownames(sce)))

print("Conversion to SingleCellExperiment object complete.")

# Load necessary libraries
library(readxl)
library(Biobase)

# Read the bulk RNA-seq data
bulk_data <- read_excel("RNA seq counts all.xlsx")

# Identify duplicate gene names
gene_names <- bulk_data[[1]]
duplicates <- gene_names[duplicated(gene_names)]
cat("Duplicate gene names found:", duplicates, "\n")

# Remove duplicates by keeping the first occurrence
bulk_data_unique <- bulk_data[!duplicated(gene_names), ]

# Convert bulk data to ExpressionSet
gene_names_unique <- bulk_data_unique[[1]]
bulk_counts_unique <- as.matrix(bulk_data_unique[,-1])
rownames(bulk_counts_unique) <- gene_names_unique
colnames(bulk_counts_unique) <- paste0("Sample_", 1:ncol(bulk_counts_unique))

bulk_eset <- ExpressionSet(assayData = bulk_counts_unique)




# Load necessary libraries
library(Biobase)

# Assuming 'data' is already read from the "28.count" file
# Convert the data to a matrix
gene_names <- data[[1]]  # Assuming the first column contains gene names
bulk_counts <- as.matrix(data[,-1])  # Assuming the second column contains counts
rownames(bulk_counts) <- gene_names
colnames(bulk_counts) <- "Sample_1"  # Name the column as a sample

# Create an ExpressionSet object
bulk_eset <- ExpressionSet(assayData = bulk_counts)

# Check the ExpressionSet object
print(bulk_eset)



# Load necessary libraries
library(Biobase)

# Read the data from the "28.count" file
bulk_data <- read.table("28.count", header = TRUE, sep = "\t", quote = "")

# Extract gene names and counts
gene_names <- bulk_data[[1]]  # Assuming the first column contains gene names
bulk_counts <- as.matrix(bulk_data[,-1])  # Assuming the second column contains counts

# Set row names to gene names
rownames(bulk_counts) <- gene_names

# Name the column as a sample
colnames(bulk_counts) <- "Sample_1"  # You can change this to a more descriptive name if needed

# Create an ExpressionSet object
bulk_eset <- ExpressionSet(assayData = bulk_counts)

# Check the ExpressionSet object
print(bulk_eset)

# Read the 28.count file
data <- read.table("28.count", header = TRUE, sep = "\t", quote = "")

# Extract gene names from the first column
gene_names <- data[[1]]

# Convert the gene names to a list
gene_list <- as.list(gene_names)

# Display the first few elements of the list
print(head(gene_list))

# Display the total number of genes
cat
("Total number of genes:", length(gene_list), "\n")

```
