library(stringr)

# Read the text file
file_path <- "/Users/Rachel/Downloads/CLASSES/20.440/project/440project/GSE141910_series_matrix.txt"
txt_data <- readLines(file_path)

# Find the lines containing sample titles
sample_title_lines <- grep("!Sample_title", txt_data, value = TRUE)

# Extract sample titles
sample_titles <- str_extract_all(sample_title_lines, "\\\"(.*?)\\\"")[[1]]
sample_titles <- gsub("\"", "", sample_titles)  # Remove quotation marks

# Find the lines containing sex information
sex_lines <- grep("!Sample_characteristics_ch1.*Sex", txt_data, value = TRUE)

# Extract sex information
sex_info <- str_extract_all(sex_lines, "Sex: (.*?)\"")[[1]]
sex_info <- gsub("Sex: ", "", sex_info)  # Remove "Sex: " part
sex_info <- gsub("\"", "", sex_info)  # Remove quotation marks

# Find the lines containing race information
race_lines <- grep("!Sample_characteristics_ch1.*race", txt_data, value = TRUE)

# Extract race information
race_info <- str_extract_all(race_lines, "race: (.*?)\"")[[1]]
race_info <- gsub("race: ", "", race_info)  # Remove "race: " part
race_info <- gsub("\"", "", race_info)  # Remove quotation marks

# Find the lines containing etiology information
etiology_lines <- grep("!Sample_characteristics_ch1.*etiology", txt_data, value = TRUE)

# Extract etiology information
etiology_info <- str_extract_all(etiology_lines, "etiology: (.*?)\"")[[1]]
etiology_info <- gsub("etiology: ", "", etiology_info)  # Remove "etiology: " part
etiology_info <- gsub("\"", "", etiology_info)  # Remove quotation marks

# Create a data frame for annotation
annotation_df <- data.frame(
  Sample_title = sample_titles,
  Sex = sex_info,
  Race = race_info,
  Etiology = etiology_info,  # Add etiology column
  stringsAsFactors = FALSE
)

head(annotation_df)
