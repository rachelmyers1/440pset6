# Load required library for plotting
library(ggplot2)
library(dplyr)

# Load the data
# combined_data_with_annotation_PCA <- read.csv("/Users/Rachel/Downloads/CLASSES/20.440/project/440project/R code/combined_data_with_annotation.csv")

# Filter rows where 'Etiology' column contains 'Dilated cardiomyopathy (DCM)'
filtered_data <- combined_data_with_annotation %>%
  filter(Etiology == "Dilated cardiomyopathy (DCM)")

# Now filtered_data contains only the rows where Etiology is 'Dilated cardiomyopathy (DCM)'

# Extract numerical data for PCA
data_for_pca <- filtered_data[, -c(1:4)] # Exclude first four columns (Sample_title, Sex, Race, Etiology)

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


# Filter pdata
filtered_pdata <- pdata %>%
  filter(etiology == "DCM")

# Add diabetes information to PC scores dataframe
pc_scores$Diabetes <- filtered_pdata$Diabetes

# Plot PCA results
ggplot(pc_scores, aes(x = PC1, y = PC2, color = Diabetes)) +
  geom_point() +
  labs(title = "PCA of DCM Patients by Diabetes Status",
       x = "PC1",
       y = "PC2") +
  scale_color_manual(values = c("purple", "#45a989")) +
  theme_minimal() +
  theme(legend.position = "top")


# In the context of the code we've developed:
#
# Each sample corresponds to a dot. The variables (genes) define the dimensions
# of the space. The principal components are linear combinations of the original
# variables. So, each dot represents a sample (e.g., an individual), and its
# position in the plot is determined by its scores on the principal components.
# The PCA plot visualizes the relationships and patterns in the data by showing
# how samples cluster together based on their similarities and differences in
# gene expression profiles.