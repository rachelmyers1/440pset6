# Merge annotation with combined data frame
library(stringr)

# Sample titles from combined_data
sample_titles_combined <- colnames(combined_data)[-1]

# # Extract the sample titles without quotation marks
# sample_titles_combined <- str_replace_all(sample_titles_combined, "\"", "")

# Transpose the combined_data dataframe
combined_data_transposed <- t(combined_data)

# Set the first row as the column names
colnames(combined_data_transposed) <- combined_data_transposed[1, ]
# Remove the first row
combined_data_transposed <- combined_data_transposed[-1, ]
# Check to make sure that the original column label "ENS Gene" is removed. Both lines should output "FALSE"
has_ENS_gene_column <- any(grepl("ENSGene", colnames(combined_data_transposed)))
has_ENS_gene_column <- any(grepl("ENS Gene", colnames(combined_data_transposed)))

# Save a copy of transposed data frame
combined_data_transposed_working <- combined_data_transposed

# Assign sample_titles_combined as a new column to combined_data_transposed dataframe
combined_data_transposed_working <- data.frame(combined_data_transposed_working, Sample_title = sample_titles_combined)

# Merge annotation with combined data frame
combined_data_with_annotation <- merge(combined_data_transposed_working, annotation_df, by = "Sample_title")

# Reorder the columns
combined_data_with_annotation <- combined_data_with_annotation[, c("Sample_title", "Sex", "Race", "Etiology", colnames(combined_data_transposed_working))]
