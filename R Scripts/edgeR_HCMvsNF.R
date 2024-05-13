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
DE_genes_hcm_nf$AdjustedPValue <- p.adjust(DE_genes_hcm_nf$PValue, method = "BH")

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



# Create a volcano plot with custom coloring and labeling
volcano_plot <- function(DE_genes_hcm_nf, threshold = 0.05, logFC_threshold = 0.3, top_n = 5) {
  # Create volcano plot
  plot(DE_genes_hcm_nf$logFC, -log10(DE_genes_hcm_nf$AdjustedPValue), pch = 20, col = ifelse(DE_genes_hcm_nf$AdjustedPValue > threshold, "grey",
                                                                                             ifelse(abs(DE_genes_hcm_nf$logFC) > logFC_threshold, ifelse(DE_genes_hcm_nf$logFC > 0, "#22bcbc", "#4457e3"), "black")),
       xlab = expression(paste("log"[2], " Fold Change")), ylab = expression(paste("-log"[10], "(Adjusted P-value)")),
       main = "DEGs in HCM vs NF Samples")
  
  # Add horizontal line for significance threshold
  abline(h = -log10(threshold), col = "grey", lty = 2)
  
  # Add vertical lines for logFC threshold
  abline(v = c(-logFC_threshold, logFC_threshold), col = "grey", lty = 2)
  
  # Define significant genes based on threshold
  significant_genes <- DE_genes_hcm_nf[DE_genes_hcm_nf$AdjustedPValue < threshold, ]
  
  # Sort significant genes by log-fold change
  significant_genes <- significant_genes[order(abs(significant_genes$logFC), decreasing = TRUE), ]
  
  # Select top 5 upregulated and top 5 downregulated genes
  top_upregulated <- head(significant_genes[significant_genes$logFC > 0, ], top_n)
  top_downregulated <- head(significant_genes[significant_genes$logFC < 0, ], top_n)
  
  # Add gene labels for top upregulated genes
  text(top_upregulated$logFC, -log10(top_upregulated$AdjustedPValue), 
       labels = top_upregulated$GeneID, 
       cex = 0.7, 
       pos = 2, 
       col = "#22bcbc")
  
  # Add gene labels for top downregulated genes
  text(top_downregulated$logFC, -log10(top_downregulated$AdjustedPValue), 
       labels = top_downregulated$GeneID, 
       cex = 0.7, 
       pos = 4, 
       col = "#4457e3")
  
  # Add legend
  legend("top", legend = c("Upregulated", "Downregulated", "Non-significant"), pch = c(20, 20, 20), col = c("#22bcbc", "#4457e3", "grey"), title = "Gene Expression", bty = "n")
}

# Call the function to create the volcano plot with top 5 upregulated and top 5 downregulated labels
volcano_plot(DE_genes_hcm_nf)
