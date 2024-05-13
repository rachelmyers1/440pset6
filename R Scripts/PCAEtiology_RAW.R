# Load required library for plotting
library(ggplot2)



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













# Extract numerical data for PCA
data_for_pca <- filtered_rawcounts_df_with_annotation[, -c(1:4)] # Exclude first three columns (Sample_title, Sex, Race)

# Check the data type of data for pca
data_type <- class(data_for_pca)
print(data_type)
# Check if all columns are numeric
all_numeric <- all(sapply(data_for_pca, is.numeric))
#print(all_numeric)

# Display the structure of the dataframe
str(data_for_pca)

# Convert character columns to numeric
data_for_pca_numeric <- data_for_pca
data_for_pca_numeric[] <- lapply(data_for_pca_numeric, as.numeric)

# Check the structure of the new dataframe
str(data_for_pca_numeric)

# Check for NA values in each column
na_columns <- colSums(is.na(data_for_pca_numeric))
# Print columns with NA values
print(names(na_columns[na_columns > 0]))
# Remove rows with NA values only from the column "Sample_title.1"
data_for_pca_numeric <- data_for_pca_numeric[, !colnames(data_for_pca_numeric) %in% "Sample_title.1"]

# Perform PCA
pca_result <- prcomp(data_for_pca_numeric, scale. = TRUE)


# Extract PC scores
pc_scores <- as.data.frame(pca_result$x)  # Convert PC scores to dataframe

# Add Etiology information to PC scores dataframe
pc_scores$Etiology <- combined_data_with_annotation$Etiology

# Plot PCA results
ggplot(pc_scores, aes(x = PC1, y = PC2, color = Etiology)) +
  geom_point() +
  labs(title = "PCA of Samples by Etiology",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +
  theme(legend.position = "right")


# In the context of the code we've developed:
#
# Each sample corresponds to a dot. The variables (genes) define the dimensions
# of the space. The principal components are linear combinations of the original
# variables. So, each dot represents a sample (e.g., an individual), and its
# position in the plot is determined by its scores on the principal components.
# The PCA plot visualizes the relationships and patterns in the data by showing
# how samples cluster together based on their similarities and differences in
# gene expression profiles.