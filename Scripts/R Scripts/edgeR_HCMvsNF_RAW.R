# Load required libraries
library(edgeR)

# Load pdata file
pdata <- read.csv("/Users/Rachel/Downloads/CLASSES/20.440/project/440project/R code/pdata.csv")

# Load counts file
counts <- read.csv("/Users/Rachel/Downloads/CLASSES/20.440/project/440project/R code/baySeq/combined_data_without_ENSGene.csv")

# Filter pdata to include only patients with HCM and NF
pdata_hcm_nf <- pdata[pdata$etiology %in% c("HCM", "NF"), ]

# Subset counts to include only samples corresponding to pdata_hcm_nf
counts_filtered <- counts[, colnames(counts) %in% pdata_hcm_nf$sample_name]

# Create DGEList object
d <- DGEList(counts = counts_filtered)

# Filter out lowly expressed genes
keep <- rowSums(cpm(d) > 1) >= 2
d <- d[keep,]

# Perform TMM normalization
d <- calcNormFactors(d)

# Estimate dispersion
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)

# Create design matrix with explicit column names
design <- model.matrix(~ etiology, data = pdata_hcm_nf)

# Rename columns to ensure valid R names
colnames(design) <- make.names(colnames(design))

# Find the column index corresponding to "HCM" and "NF" levels
hcm_col_index <- which(colnames(design) == "etiologyHCM")
nf_col_index <- which(colnames(design) == "etiologyNF")

# Create contrast matrix
contrast_matrix <- matrix(0, nrow = ncol(design), ncol = 1)
contrast_matrix[hcm_col_index] <- 1
contrast_matrix[nf_col_index] <- -1

# Perform differential expression analysis
fit <- glmFit(d, design)
contrast_fit <- glmLRT(fit, contrast = contrast_matrix)

# Get differentially expressed genes
DE_genes_hcm_nf <- topTags(contrast_fit, n = Inf)$table
DE_genes_hcm_nf <- DE_genes_hcm_nf[order(DE_genes_hcm_nf$PValue), ]

# Adjust p-values for multiple testing
DE_genes_hcm_nf$AdjustedPValue <- p.adjust(DE_genes_hcm_nf$PValue, method = "fdr")

# Filter significant differentially expressed genes (optional)
significant_genes_hcm_nf <- DE_genes_hcm_nf[DE_genes_hcm_nf$AdjustedPValue < 0.05, ]

# View the results
head(significant_genes_hcm_nf)

# Define the file path
output_folder <- "/Users/Rachel/Downloads/CLASSES/20.440/project/440project/R code outputs/"
output_file <- paste0(output_folder, "DEGs_HCMvsNF")

# Write the ENSgene column to a CSV file
write.csv(DE_genes_hcm_nf, file = output_file, row.names = FALSE)




# Extract the row names (ENS gene names) of the significant_genes data frame
significant_ens_genes_hcm_nf <- rownames(significant_genes_hcm_nf)

# Convert significant_ens_genes from characters to numbers
significant_ens_genes_numeric_hcm_nf <- as.numeric(significant_ens_genes_hcm_nf)

# Filter rows in combined_data based on row indexes in significant_ens_genes_numeric
filtered_combined_data_hcm_nf <- combined_data[significant_ens_genes_numeric_hcm_nf, ]

# Display the first few rows of the ENSGene column in the filtered combined_data
print(filtered_combined_data_hcm_nf$ENSGene)

# Define the file path
output_folder <- "/Users/Rachel/Downloads/CLASSES/20.440/project/440project/R code outputs/"
output_file <- paste0(output_folder, "SignficantGenes_HCMvsNF")

# Write the ENSgene column to a CSV file
write.csv(filtered_combined_data_hcm_nf$ENSGene, file = output_file, row.names = FALSE)
