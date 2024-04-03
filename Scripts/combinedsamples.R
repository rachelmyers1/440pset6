library(dplyr)
library(tidyr)
library(purrr)

#NOTE: Raw gene counts were transformed using SVA and the 'LIMMA-Voom' method

# Specify the directory containing CSV files
directory <- "/Users/Rachel/Downloads/CLASSES/20.440/project/440project/unzipped_files/"

# Get a list of all CSV files in the directory
csv_files <- list.files(path = directory, pattern = "\\.csv$", full.names = TRUE)

# Function to read CSV files correctly and add column headers
read_csv_with_header <- function(file_path) {
  # Read the CSV file
  data <- read.csv(file_path, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  # Add a column for 'ENSGene' with appropriate values
  data$ENSGene <- rownames(data)
  # Reset row names
  rownames(data) <- NULL
  return(data)
}

# Read and combine the CSV files into a list of data frames
count_lists <- map(csv_files, read_csv_with_header)

# Join the data frames in the list
combined_data <- reduce(count_lists, full_join, by = 'ENSGene')

# Unite the columns containing gene names into one column
combined_data <- combined_data %>%
  unite(ENSGene, starts_with("X"), sep = ",", remove = TRUE)

# Remove duplicate gene names within each row
combined_data$ENSGene <- sapply(combined_data$ENSGene, function(x) toString(unique(unlist(strsplit(x, ",")))))

# Reorder columns with 'ENSGene' as the first column
combined_data <- combined_data %>%
  select(ENSGene, everything())

# head(combined_data)
