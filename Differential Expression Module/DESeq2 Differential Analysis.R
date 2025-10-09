library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(HDF5Array)
library(SummarizedExperiment)
library(DESeq2)
library(tibble)
library(tidyr)
library(stringr)

# Initialize
rm(list = ls())

# Load the file (Change the path according the actual path)
barcodes_file <- "/Users/a1234/Desktop/Processed_L2_3 IT/barcodes.tsv"
genes_file <- "/Users/a1234/Desktop/Processed_L2_3 IT/genes.tsv"
matrix_file <- "/Users/a1234/Desktop/Processed_L2_3 IT/matrix.mtx"
metadata_file <- "/Users/a1234/Desktop/Processed_L2_3 IT/metadata.csv"

barcodes <- readLines(barcodes_file)
genes <- read.table(genes_file, stringsAsFactors = FALSE)

# Check the uniqueness of the genes
gene_ids_unique <- make.unique(genes[, 1])

# Read the expression matrix in blocks
con <- gzfile(matrix_file, "r")
readLines(con, n = 2)  

# Set the blocks
chunk_size <- 300000
col_index <- integer(0)
row_index <- integer(0)
values <- numeric(0)

# Read the matrix
while (length(chunk <- readLines(con, n = chunk_size)) > 0) {
  chunk_data <- read.table(text = chunk)
  col_index <- c(col_index, chunk_data$V2)
  row_index <- c(row_index, chunk_data$V1)
  values <- c(values, chunk_data$V3)
  rm(chunk_data)
  gc()
}
close(con)

# Create a sparse expression matrix
sparse_matrix <- sparseMatrix(
  i = row_index, 
  j = col_index, 
  x = values, 
  dims = c(nrow(genes), length(barcodes))
)
colnames(sparse_matrix) <- barcodes
rownames(sparse_matrix) <- gene_ids_unique

# Filter out low-expression genes
cat("Number of original genes:", nrow(sparse_matrix), "\n")
min_cells_expressing <- 3
cells_per_gene <- rowSums(sparse_matrix > 0)
genes_to_keep <- cells_per_gene >= min_cells_expressing
sparse_matrix <- sparse_matrix[genes_to_keep, ]
gene_ids_unique <- gene_ids_unique[genes_to_keep]
cat("Number of filtered genes:", nrow(sparse_matrix), "\n")

# Transform the sparse matrix into the format of HDF5 file
h5file <- tempfile(fileext = ".h5")
data <- writeHDF5Array(sparse_matrix, filepath = h5file)
rm(sparse_matrix)
gc()

# Read the metadata
metadata <- read.csv(
  metadata_file, 
  row.names = 1, 
  check.names = TRUE, 
  stringsAsFactors = FALSE
)
cat("The column names of metadata:\n"); print(colnames(metadata))

# Downsampling the cells
set.seed(123)
common_cells <- intersect(colnames(data), rownames(metadata))
cat("Common cell count:", length(common_cells), "\n")

n_sample_cells <- 10000
if (length(common_cells) < n_sample_cells) {
  selected_cells <- common_cells
  warning("If the number of common cells is less than 10,000, all common cells should be used for downsampling.")
} else {
  selected_cells <- sample(common_cells, size = n_sample_cells)
}

sub_data <- data[, selected_cells]
sub_metadata <- metadata[selected_cells, , drop = FALSE]

# Normalize the column names
if(!("Cognitive.Status" %in% colnames(sub_metadata))) {
  if("Cognitive Status" %in% colnames(sub_metadata)) {
    sub_metadata$Cognitive.Status <- sub_metadata[["Cognitive Status"]]
  } else {
    stop("The Cognitive.Status column cannot be found in the metadata. Please check the metadata.csv column name.")
  }
}

# Convert Cognitive States into factors
sub_metadata$Cognitive.Status <- factor(
  sub_metadata$Cognitive.Status, 
  levels = c("No dementia", "Dementia")
)

# Create Seurat object
obj <- CreateSeuratObject(
  counts = sub_data, 
  meta.data = sub_metadata
)
rm(sub_data, sub_metadata)
gc()

# Convert the HDF5 matrix to the dgCMatrix
counts_delayed <- LayerData(obj, assay = "RNA", layer = "counts")
counts_matrix <- as(counts_delayed, "dgCMatrix")
obj <- SetAssayData(
  obj, 
  assay = "RNA", 
  layer = "counts", 
  new.data = counts_matrix
)
rm(counts_delayed, counts_matrix)
gc()


# Obtain the metadata column name of the Seurat object
meta_names <- colnames(obj@meta.data)

# Determine the donor ID column name
donor_col_name <- if("Donor.ID" %in% meta_names) "Donor.ID" else 
  if("Donor ID" %in% meta_names) "Donor ID" else 
    stop("Cannot find the column of Donor ID，Check the column names of Seurat meta.data.")

# Wash the format of Donor ID
obj$Donor.ID.Clean <- obj[[donor_col_name]]

if(!("Cognitive.Status" %in% meta_names)) {
  if("Cognitive Status" %in% meta_names) {
    obj$Cognitive.Status <- obj[["Cognitive Status"]]
  } else {
    stop("Cannot find the column of Cognitive.Status / Cognitive Status in Seurat meta.data.")
  }
}

# Wash the format of Cognitive Status
obj$Cognitive.Status.Clean <- as.character(obj$Cognitive.Status)
obj$Cognitive.Status.Clean[is.na(obj$Cognitive.Status.Clean)] <- "Unknown"
obj$Cognitive.Status.Clean <- gsub(" ", ".", obj$Cognitive.Status.Clean)

obj$Cognitive.Status.Clean <- factor(
  obj$Cognitive.Status.Clean, 
  levels = c("No.dementia", "Dementia", "Unknown")
)

# ------------Create pseudobulk data---------
pseudobulk_counts <- AggregateExpression(
  obj, 
  group.by = c("Cognitive.Status.Clean", "Donor.ID.Clean"), 
  assays = "RNA", 
  slot = "counts", 
  return.seurat = FALSE
)$RNA

# ------------Create colData for DESeq2---------------
col_names <- colnames(pseudobulk_counts)

cat("Check the colnames (top 40)：\n"); print(head(col_names, 40))
cat("Total numbers：", length(col_names), "，the number of columns, including '-'：", sum(grepl("-", col_names)), "\n")

colData <- tibble(orig = col_names) %>%
  extract(
    col = "orig", 
    into = c("Cognitive.Status.Clean", "Donor.ID"), 
    regex = "^([^_]+)_(.*)$", 
    remove = FALSE
  ) %>%
  as.data.frame()
rownames(colData) <- colData$orig

# Transform into factors
colData$Cognitive.Status <- factor(
  colData$Cognitive.Status.Clean, 
  levels = c("No.dementia", "Dementia")
)
colData$Donor.ID <- factor(colData$Donor.ID)

cat("The column names of the original metadata (for aggregating donor-level)：\n"); print(colnames(metadata))

# Match the covariate column names
meta_colnames <- colnames(metadata)
donor_meta_col <- if("Donor.ID" %in% meta_colnames) "Donor.ID" else 
  if("Donor ID" %in% meta_colnames) "Donor ID" else 
    stop("Cannot find the 'Donor ID' in metadata.")

age_col <- if("Age.at.Death" %in% meta_colnames) "Age.at.Death" else 
  if("Age at Death" %in% meta_colnames) "Age at Death" else 
    stop("Cannot find the 'Age at Death' in metadata.")

sex_col <- if("Sex" %in% meta_colnames) "Sex" else 
  stop("Cannot find the 'Sex' in metadata.")

pmi_col <- if("PMI" %in% meta_colnames) "PMI" else 
  stop("Cannot find the 'PMI' in metadata.")

# Wash the column of "PMI"
metadata <- metadata %>%
  mutate(!!pmi_col := ifelse(grepl("^[0-9.]+$", .data[[pmi_col]]), .data[[pmi_col]], NA))

# Define the function for safely converting to numeric types
safe_as_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

# Aggregate covariates at the donor level (calculated by donor grouping)
# Rename the "Donor ID" into "DonorID_for_agg"
metadata2 <- metadata %>%
  mutate(DonorID_for_agg = as.character(.data[[donor_meta_col]]))

# Aggregation by donor group
donor_cov <- metadata2 %>%
  group_by(DonorID_for_agg) %>%
  summarise(
    Donor.ID = {
      id_raw <- dplyr::first(DonorID_for_agg)
      id_clean <- gsub("[^A-Za-z0-9.]", "", id_raw)
      id_final <- make.names(id_clean, unique = FALSE)
      id_final 
    },
    Age_at_Death = {
      # Mean
      v <- safe_as_numeric(.data[[age_col]])
      v <- v[!is.na(v)]
      if (length(v) == 0) NA_real_ else mean(v, na.rm = TRUE)
    },
    Sex = {
      # Mode
      s <- as.character(.data[[sex_col]])
      s <- s[!is.na(s) & s != ""]
      if (length(s) == 0) NA_character_ else names(which.max(table(s)))
    },
    PMI = {
      # Mean
      v <- safe_as_numeric(.data[[pmi_col]])
      v <- v[!is.na(v)]
      if (length(v) == 0) NA_real_ else mean(v, na.rm = TRUE)
    },
    .groups = "drop"
  )

cat("The first few lines of the aggregated donor_cov (with the ID cleaned)：\n"); print(head(donor_cov))

match_count <- sum(colData$Donor.ID %in% donor_cov$Donor.ID)
cat("The number of Donor ids found in colData that can match donor_cov: ", match_count, "/", nrow(colData), "\n")

colData <- as.data.frame(colData) %>%
  left_join(donor_cov, by = "Donor.ID")
rownames(colData) <- colData$orig


# Deal with the missing value
na_counts <- colSums(is.na(colData[, c("Age_at_Death","Sex","PMI")]))
na_prop <- na_counts / nrow(colData)
cat("Missing proportion: \n"); print(na_prop)

covariates <- c("Age_at_Death","Sex","PMI")
for (cov in covariates) {
  if (na_prop[cov] > 0.2) {
    cat("Variable", cov, "has too high missing rate (", round(na_prop[cov]*100,1),"%)，Remove. \n")
    colData[[cov]] <- NULL
  } else {
    if (cov %in% c("Age_at_Death","PMI")) {
      median_val <- median(colData[[cov]], na.rm = TRUE)
      colData[[cov]][is.na(colData[[cov]])] <- median_val
      cat("Numerical variable", cov, "The missing values have been filled with the median", median_val, "\n")
    } else if (cov == "Sex") {
      n_before <- nrow(colData)
      colData <- colData[!is.na(colData$Sex), ]
      n_after <- nrow(colData)
      cat("Categorical variable Sex：has been deleted", n_before - n_after, "missing sample(s)。\n")
    }
  }
}

# Convert variable types
if ("Age_at_Death" %in% colnames(colData)) colData$Age_at_Death <- as.numeric(colData$Age_at_Death)
if ("PMI" %in% colnames(colData)) colData$PMI <- as.numeric(colData$PMI)
if ("Sex" %in% colnames(colData)) colData$Sex <- factor(colData$Sex)

# Remove the samples with a cognitive status of "Unknown"
colData <- subset(colData, Cognitive.Status.Clean != "Unknown")
pseudobulk_counts <- pseudobulk_counts[, rownames(colData), drop = FALSE]
cat("Alignment successful. Ultimate dimension: ", dim(pseudobulk_counts), "\n")
print(table(colData$Cognitive.Status))

# ----------DESeq2 Differential Analysis----------
# Filter out low-expression genes
keep <- rowSums(pseudobulk_counts) >= 10
pseudobulk_counts <- pseudobulk_counts[keep, ]

# Extract the final gene list for DESeq2 analysis as the background set
background_genes <- rownames(pseudobulk_counts)
write.table(
  background_genes,
  file = "/Users/a1234/Desktop/Dementia_Background_Genes.txt",
  sep = "\n",
  row.names = FALSE,
  col.names = FALSE
)

# Make sure the expression is an integer
pseudobulk_counts <- round(pseudobulk_counts)

# Select the analysis mode
use_pair_mode <- FALSE  

if(use_pair_mode) {
  message("Run the pairing mode: design = ~ Donor.ID + Cognitive.Status")
  dds 
  dds <- DESeqDataSetFromMatrix(
    countData = pseudobulk_counts,
    colData = colData,
    design = ~ Donor.ID + Cognitive.Status
  )
} else {
  na_covariates <- colData %>%
    select(Age_at_Death, Sex, PMI) %>%
    summarise_all(~sum(is.na(.)))
  cat("The number of nas per covariate: \n"); print(na_covariates)
  
  if(any(unlist(na_covariates) > 0)) {
    warning("There are samples with covariate containing NA. 
            DESeq2 cannot accept covariates containing NA. 
            Please complete or delete these samples before running.")
    complete_idx <- complete.cases(colData[, c("Age_at_Death", "Sex", "PMI")])
    if(sum(!complete_idx) > 0) {
      warning("The number of samples containing NA to be deleted: ", sum(!complete_idx))
      colData <- colData[complete_idx, , drop = FALSE]
      pseudobulk_counts <- pseudobulk_counts[, rownames(colData), drop = FALSE]
    }
  }
  
  # Create DESeq2 object
  message("Run in standard mode: design = ~ Age_at_Death + Sex + PMI + Cognitive.Status")
  dds <- DESeqDataSetFromMatrix(
    countData = pseudobulk_counts,
    colData = colData,
    design = ~ Age_at_Death + Sex + PMI + Cognitive.Status
  )
}

# Conduct a DESeq2 Difference Analysis
dds <- DESeq(dds)

# Extract the difference results (Dementia vs. No dementia)
res <- results(dds, contrast = c("Cognitive.Status", "Dementia", "No.dementia"))

# Deal with the result
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[complete.cases(res_df), ]

# Calculate the 95% confidence interval(CI)
res_df <- res_df %>%
  mutate(
    CI_lower = log2FoldChange - 1.96 * lfcSE, 
    CI_upper = log2FoldChange + 1.96 * lfcSE  
  )


class.de <- dplyr::rename(
  res_df,
  base_mean = baseMean,
  avg_log2FC = log2FoldChange,
  lfcSE = lfcSE,
  test_stat = stat,
  p_val = pvalue,
  p_val_adj = padj
) %>%
  dplyr::select(gene, base_mean, avg_log2FC, lfcSE, CI_lower, CI_upper, test_stat, p_val, p_val_adj) %>%
  dplyr::arrange(p_val_adj)


# Screen for genes with significant differences (adjusted p value <0.05)
sig_genes <- subset(class.de, p_val_adj < 0.05)
head(class.de, 10)


write.csv(
  sig_genes,
  file = "/Users/a1234/Desktop/Dementia_Significant_DEGs_Full_Report.csv",
  row.names = FALSE
)
cat("The complete report (including logFC and 95% CI) has been exported as a CSV file.\n")


# -------Visualize--------
volcano_data <- class.de %>%
  mutate(
    log10_padj = -log10(p_val_adj + 1e-300),
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.5 ~ "Up in Dementia",
      p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Down in Dementia",
      TRUE ~ "Not significant"
    )
  )

ggplot(volcano_data, aes(x = avg_log2FC, y = log10_padj, color = significance)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(
    x = "Log2 Fold Change (Dementia vs Control)",
    y = "-Log10(Adjusted P-value)"
  ) +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed")

# Extract up-down-regulated genes
de_genes <- sig_genes$gene
up_genes <- sig_genes$gene[sig_genes$avg_log2FC > 0]
down_genes <- sig_genes$gene[sig_genes$avg_log2FC < 0]

# Export the results
write.table(de_genes, file = "/Users/a1234/Desktop/Dementia_DEGs.txt", 
            sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(up_genes, file = "/Users/a1234/Desktop/Dementia_Up_Genes.txt", 
            sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(down_genes, file = "/Users/a1234/Desktop/Dementia_Down_Genes.txt", 
            sep = "\n", row.names = FALSE, col.names = FALSE)

cat("The Pipeline is completed. The significant genes have been exported to the desktop.\n")

