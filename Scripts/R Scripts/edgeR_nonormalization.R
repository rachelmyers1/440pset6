# Load required libraries
library(edgeR)

# Load pdata file
pdata <- read.csv("/Users/Rachel/Downloads/CLASSES/20.440/project/440project/R code/pdata.csv")

# Filter pdata to include only patients with DCM
pdata_dcm <- pdata[pdata$etiology == "DCM", ]

# Filter pdata_dcm to include only African American and Caucasian patients
pdata_filtered <- pdata_dcm[pdata_dcm$race %in% c("AA", "Caucasian"), ]

# Load counts file
counts <- read.csv("/Users/Rachel/Downloads/CLASSES/20.440/project/440project/R code/baySeq/combined_data_without_ENSGene.csv")





# Subset counts to include only samples corresponding to pdata_filtered
counts_filtered <- counts[, colnames(counts) %in% pdata_filtered$sample_name]

# Create DGEList object
d <- DGEList(counts = counts_filtered)

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





# Set dispersion values manually
num_rows <- nrow(counts_filtered)
d$common.dispersion <- rep(1, num_rows) # Set a common dispersion value
d$tagwise.dispersion <- rep(0, num_rows)  # Set tagwise dispersion values





# Perform differential expression analysis
fit <- glmFit(d, design)

# Estimate contrasts
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

# Create a volcano plot
volcano_plot <- function(DE_genes, threshold = 0.05, logFC_threshold = 1) {
  # Define significant genes based on threshold
  significant_genes <- DE_genes[DE_genes$AdjustedPValue < threshold & abs(DE_genes$logFC) > logFC_threshold, ]
  
  # Create volcano plot
  plot(DE_genes$logFC, -log10(DE_genes$AdjustedPValue), pch = 20, col = ifelse(abs(DE_genes$logFC) > logFC_threshold & DE_genes$AdjustedPValue < threshold, "red", "black"),
       xlab = "log2 Fold Change", ylab = "-log10(Adjusted P-value)",
       main = "Volcano Plot")
  
  # Add horizontal line for significance threshold
  abline(h = -log10(threshold), col = "grey", lty = 2)
  
  # Add vertical lines for logFC threshold
  abline(v = c(-0.3, 0.3), col = "grey", lty = 2)
}

# Call the function to create the volcano plot
volcano_plot(DE_genes)
