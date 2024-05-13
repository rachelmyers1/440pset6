# Load required libraries
library(edgeR)

# Load pdata file
pdata <- read.csv("/Users/Rachel/Downloads/CLASSES/20.440/project/440project/R code/pdata.csv")

# Load counts file
counts <- read.csv("/Users/Rachel/Downloads/CLASSES/20.440/project/440project/R code/baySeq/combined_data_without_ENSGene.csv")

# Filter pdata to include only patients with DCM
pdata_dcm <- pdata[pdata$etiology == "DCM", ]

# Filter pdata_dcm to include only African American and Caucasian patients
pdata_filtered <- pdata_dcm[pdata_dcm$race %in% c("AA", "Caucasian"), ]

# Subset counts to include only samples corresponding to pdata_filtered
counts_filtered <- counts[, colnames(counts) %in% pdata_filtered$sample_name]


# Create DGEList object
d <- DGEList(counts = counts_filtered)

# Perform filtering and normalization
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)




# Create design matrix with explicit column names
design <- model.matrix(~ race, data = pdata_filtered)

# Rename columns to ensure valid R names
colnames(design) <- make.names(colnames(design))

# Find the column index corresponding to "Caucasian" and "AA" levels
caucasian_col_index <- which(colnames(design) == "raceCaucasian")
aa_col_index <- which(colnames(design) == "raceAA")

# Create contrast matrix
contrast_matrix <- matrix(0, nrow = ncol(design), ncol = 1)
contrast_matrix[caucasian_col_index] <- -1
contrast_matrix[aa_col_index] <- 1

# Perform differential expression analysis
fit <- glmFit(d, design)
contrast_fit <- glmLRT(fit, contrast = contrast_matrix)

# Get differentially expressed genes
DE_genes <- topTags(contrast_fit, n = Inf)$table
DE_genes <- DE_genes[order(DE_genes$PValue), ]

# Adjust p-values for multiple testing
DE_genes$AdjustedPValue <- p.adjust(DE_genes$PValue, method = "fdr")

# Filter significant differentially expressed genes (optional)
significant_genes <- DE_genes[DE_genes$AdjustedPValue < 0.05, ]

# View the results
head(significant_genes)






# Read the combined_data.csv file to get ENS gene names
combined_data <- read.csv("combined_data.csv")

# Extract the row names (ENS gene names) of the significant_genes data frame
significant_ens_genes <- rownames(significant_genes)

# Convert significant_ens_genes from characters to numbers
significant_ens_genes_numeric <- as.numeric(significant_ens_genes)

# Filter rows in combined_data based on row indexes in significant_ens_genes_numeric
filtered_combined_data <- combined_data[significant_ens_genes_numeric, ]

# Display the first few rows of the ENSgene column in the filtered combined_data
print(filtered_combined_data$ENSGene)

# Define the file path
output_folder <- "/Users/Rachel/Downloads/CLASSES/20.440/project/440project/R code outputs/"
output_file <- paste0(output_folder, "SignficantGenes_DCM_AAvsCaucasian")

# Write the ENSgene column to a CSV file
write.csv(filtered_combined_data$ENSGene, file = output_file, row.names = FALSE)

# Load required library
library(openxlsx)

# Read the Excel file with common gene names
common_gene_names <- read.xlsx("/Users/Rachel/Downloads/CLASSES/20.440/project/440project/R code outputs/signficantgenes_DCM_AAvsCaucasian.xlsx")







# Create a volcano plot with custom coloring and labeling
volcano_plot <- function(DE_genes, threshold = 0.05, logFC_threshold = 0.3) {
  # Create volcano plot
  plot(DE_genes$logFC, -log10(DE_genes$AdjustedPValue), pch = 20, col = ifelse(DE_genes$AdjustedPValue > threshold, "grey",
                                                                               ifelse(abs(DE_genes$logFC) > logFC_threshold, ifelse(DE_genes$logFC > 0, "#22bcbc", "#4457e3"), "black")),
       xlab = expression(paste("log"[2], " Fold Change")), ylab = expression(paste("-log"[10], "(Adjusted P-value)")),
       main = "DEGs in AA DCM Patients Relative to Caucasian DCM Patients")
  
  # Add horizontal line for significance threshold
  abline(h = -log10(threshold), col = "grey", lty = 2)
  
  # Add vertical lines for logFC threshold
  abline(v = c(-logFC_threshold, logFC_threshold), col = "grey", lty = 2)
  
  # Define significant genes based on threshold
  significant_genes <- DE_genes[DE_genes$AdjustedPValue < threshold, ]
  
  # Add gene labels for significant genes with common names
  text(significant_genes$logFC, -log10(significant_genes$AdjustedPValue), labels = common_gene_names$Common.Gene.Name, cex = 0.7, pos = 4, col = ifelse(abs(significant_genes$logFC) > logFC_threshold, ifelse(significant_genes$logFC > 0, "#22bcbc", "#013746"), "black"))
  
  # Add legend
  legend("top", legend = c("Upregulated", "Downregulated", "Non-significant"), pch = c(20, 20, 20), col = c("#22bcbc", "#4457e3", "grey"), title = "Gene Expression", bty = "n")
}

# Call the function to create the volcano plot
volcano_plot(DE_genes)








# Create a volcano plot with custom coloring and labeling
volcano_plot <- function(DE_genes, threshold = 0.05, logFC_threshold = 0.3) {
  # Create volcano plot
  plot(DE_genes$logFC, -log10(DE_genes$AdjustedPValue), pch = 20, col = ifelse(DE_genes$AdjustedPValue > threshold, "grey",
                                                                               ifelse(abs(DE_genes$logFC) > logFC_threshold, ifelse(DE_genes$logFC > 0, "#22bcbc", "#4457e3"), "black")),
       xlab = expression(paste("log"[2], " Fold Change")), ylab = expression(paste("-log"[10], "(Adjusted P-value)")),
       main = "DEGs in AA DCM Patients Relative to Caucasian DCM Patients")
  
  # Add horizontal line for significance threshold
  abline(h = -log10(threshold), col = "grey", lty = 2)
  
  # Add vertical lines for logFC threshold
  abline(v = c(-logFC_threshold, logFC_threshold), col = "grey", lty = 2)
  
  # Define significant genes based on threshold
  significant_genes <- DE_genes[DE_genes$AdjustedPValue < threshold, ]
  
  # Add gene labels for significant genes with common names
  text(significant_genes$logFC, -log10(significant_genes$AdjustedPValue), 
       labels = common_gene_names$Common.Gene.Name, 
       cex = 0.7, 
       pos = ifelse(significant_genes$logFC > 0, 2, 4), 
       col = ifelse(abs(significant_genes$logFC) > logFC_threshold, ifelse(significant_genes$logFC > 0, "#22bcbc", "#013746"), "black"))
  
  # Add legend
  legend("top", legend = c("Upregulated", "Downregulated", "Non-significant"), pch = c(20, 20, 20), col = c("#22bcbc", "#4457e3", "grey"), title = "Gene Expression", bty = "n")
}

# Call the function to create the volcano plot
volcano_plot(DE_genes)
