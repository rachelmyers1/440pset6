# Load required libraries
library(edgeR)
library(dplyr)




# Specify the file path
file_path <- "/Users/Rachel/Downloads/CLASSES/20.440/project/440project/R code/RawCounts.csv"

# Read the CSV file into a dataframe
counts_df <- read.csv(file_path, sep = ",", stringsAsFactors = FALSE, header = TRUE)  # Set header = TRUE if the file has column names

# Identify and remove rows with duplicated gene names in counts_df
unique_counts_df <- counts_df[!duplicated(counts_df$X), ]

# Print the dimensions of the original and unique counts dataframes
cat("Original counts dataframe dimensions:", dim(counts_df), "\n")
cat("Unique counts dataframe dimensions:", dim(unique_counts_df), "\n")

# Identify rows with all zero values across columns (starting from the second column)
zero_rows <- apply(counts_df[, -1], 1, function(row) all(row == 0))

# Subset counts_df to keep rows where not all values are zero across columns (starting from the second column)
filtered_counts_df <- counts_df[!zero_rows, ]

# Print the dimensions of the original and filtered counts dataframes
cat("Original counts dataframe dimensions:", dim(counts_df), "\n")
cat("Filtered counts dataframe dimensions:", dim(filtered_counts_df), "\n")

# Filter out genes with low counts across all samples
keep <- rowSums(filtered_counts_df[, -1] > 10) >= 3  # Keep genes with at least 3 samples having counts > 10
filtered_rawcounts_df <- filtered_counts_df[keep, ]

cat("Filtered counts dataframe dimensions:", dim(filtered_rawcounts_df), "\n")

# Remove first column with the ENS gene names
filtered_rawcounts_df_withoutENS <- filtered_rawcounts_df[, -1]  # Exclude first column (Gene names)

# Add 1 to every numerical cell
rawcounts_filtered <- filtered_rawcounts_df_withoutENS %>%
  mutate_if(is.numeric, ~ . + 1)

# Transpose the data frame
#rawcounts_filtered <- t(rawcounts_filtered)





# Load pdata file
pdata <- read.csv("/Users/Rachel/Downloads/CLASSES/20.440/project/440project/R code/pdata.csv")

# Filter pdata to include only patients with DCM and NF
pdata_dcm_nf <- pdata[pdata$etiology %in% c("DCM", "NF"), ]

# Subset counts to include only samples corresponding to pdata_filtered
filtered_rawcounts_bysample <- rawcounts_filtered[, colnames(rawcounts_filtered) %in% pdata_dcm_nf$sample_name]



# Create DGEList object
d <- DGEList(counts = filtered_rawcounts_bysample)

# Filter out lowly expressed genes
#keep <- rowSums(cpm(d) > 1) >= 2
#d <- d[keep,]

# Perform TMM normalization
d <- calcNormFactors(d)

# Estimate dispersion
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)

# Create design matrix with explicit column names
design <- model.matrix(~ etiology, data = pdata_dcm_nf)

# Rename columns to ensure valid R names
colnames(design) <- make.names(colnames(design))

# Find the column index corresponding to "DCM" and "NF" levels
dcm_col_index <- which(colnames(design) == "etiologyDCM")
nf_col_index <- which(colnames(design) == "etiologyNF")

# Create contrast matrix
contrast_matrix <- matrix(0, nrow = ncol(design), ncol = 1)
contrast_matrix[dcm_col_index] <- 1
contrast_matrix[nf_col_index] <- -1

# Perform differential expression analysis
fit <- glmFit(d, design)
contrast_fit <- glmLRT(fit, contrast = contrast_matrix)

# Get differentially expressed genes
DE_genes_dcm_nf <- topTags(contrast_fit, n = Inf)$table
DE_genes_dcm_nf <- DE_genes_dcm_nf[order(DE_genes_dcm_nf$PValue), ]

# Adjust p-values for multiple testing
DE_genes_dcm_nf$AdjustedPValue <- p.adjust(DE_genes_dcm_nf$PValue, method = "fdr")

# Filter significant differentially expressed genes (optional)
significant_genes_dcm_nf <- DE_genes_dcm_nf[DE_genes_dcm_nf$AdjustedPValue < 0.05, ]

# View the results
head(significant_genes_dcm_nf)

# Define the file path
output_folder <- "/Users/Rachel/Downloads/CLASSES/20.440/project/440project/R code outputs/"
output_file <- paste0(output_folder, "rawDEGs_DCMvsNF_RAW")

# Write the ENSgene column to a CSV file
write.csv(DE_genes_dcm_nf, file = output_file, row.names = FALSE)




# Extract the row names (ENS gene names) of the significant_genes data frame
significant_ens_genes_dcm_nf <- rownames(significant_genes_dcm_nf)

# Convert significant_ens_genes from characters to numbers
significant_ens_genes_numeric_dcm_nf <- as.numeric(significant_ens_genes_dcm_nf)

# Filter rows in combined_data based on row indexes in significant_ens_genes_numeric
filtered_rawcounts_dcm_nf <- filtered_rawcounts_df[significant_ens_genes_numeric_dcm_nf, ]

filtered_rawcounts_dcm_nf <- filtered_rawcounts_dcm_nf[!is.na(filtered_rawcounts_dcm_nf$X), ]

# Display the first few rows of the ENSGene column in the filtered combined_data
print(filtered_rawcounts_dcm_nf$X)



# Define the file path
output_folder <- "/Users/Rachel/Downloads/CLASSES/20.440/project/440project/R code outputs/"
output_file <- paste0(output_folder, "SignficantGenes_DCMvsNF")

# Write the ENSgene column to a CSV file
write.csv(filtered_combined_data_dcm_nf$ENSGene, file = output_file, row.names = FALSE)
