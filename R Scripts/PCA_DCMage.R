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



# Convert 'Age' column to numeric
filtered_pdata$age <- as.numeric(as.character(filtered_pdata$age))

# Check for NAs or non-numeric values
na_or_non_numeric <- filtered_pdata$age[!is.na(filtered_pdata$age) & !is.numeric(filtered_pdata$age)]
if (length(na_or_non_numeric) > 0) {
  print("Warning: There are NA or non-numeric values in the 'age' column.")
  print(na_or_non_numeric)
}






cat("Minimum Age:", min(pdata$age), "\nMaximum Age:", max(pdata$age))


# Create 4 age intervals
age_intervals <- cut(filtered_pdata$age, breaks = c(0, 16, 30, 55, Inf), labels = c("0-16", "17-30", "31-55", "55+"))

# Add age intervals to PC scores dataframe
pc_scores$Age <- age_intervals

# Plot PCA results with age group coloring
ggplot(pc_scores, aes(x = PC1, y = PC2, color = Age)) +
  geom_point() +
  labs(title = "PCA of DCM Patients by Age Group",
       x = "PC1",
       y = "PC2") +
  scale_color_manual(values = c("#45a989", "purple", "blue", "red")) + # Define colors for age groups
  theme_minimal() +
  theme(legend.position = "top")







# Now try creating eight age intervals
age_intervals <- cut(filtered_pdata$age, breaks = c(0, 16, 30, 55, 60, 65, 70, 75, Inf), labels = c("0-16", "17-30", "31-55", "56-60", "61-65", "66-70", "71-75", "75+"))

# Add age intervals to PC scores dataframe
pc_scores$age <- age_intervals

# Plot PCA results with age group coloring
ggplot(pc_scores, aes(x = PC1, y = PC2, color = age)) +
  geom_point() +
  labs(title = "PCA of DCM Patients by Age Group",
       x = "PC1",
       y = "PC2") +
  scale_color_manual(values = c("#45a989", "purple", "blue", "orange", "yellow", "green", "red", "black")) + # Define colors for age groups
  theme_minimal() +
  theme(legend.position = "top")




# Create two age intervals
age_intervals <- cut(filtered_pdata$age, breaks = c(0, 56, Inf), labels = c("0-56", "56+"))

# Add age intervals to PC scores dataframe
pc_scores$age <- age_intervals

# Plot PCA results with age group coloring
ggplot(pc_scores, aes(x = PC1, y = PC2, color = age)) +
  geom_point() +
  labs(title = "PCA of DCM Patients by Age Group",
       x = "PC1",
       y = "PC2") +
  scale_color_manual(values = c("#45a989", "purple", "blue", "orange", "yellow", "green", "red", "black")) + # Define colors for age groups
  theme_minimal() +
  theme(legend.position = "top")







# Convert 'age' column to numeric
filtered_pdata$age <- as.numeric(as.character(filtered_pdata$age))

# Check for NAs or non-numeric values
na_or_non_numeric <- filtered_pdata$age[!is.na(filtered_pdata$age) & !is.numeric(filtered_pdata$age)]
if (length(na_or_non_numeric) > 0) {
  print("Warning: There are NA or non-numeric values in the 'age' column.")
  print(na_or_non_numeric)
}

cat("Minimum Age:", min(filtered_pdata$age), "\nMaximum Age:", max(filtered_pdata$age))

# Define color scale
library(viridis)

# Now try creating four age intervals
age_intervals <- cut(filtered_pdata$age, breaks = c(0, 35, 56, 80, Inf), labels = c("0-35", "35-56", "57-80", "81+"))

# Add age intervals to PC scores dataframe
pc_scores$Age <- filtered_pdata$age

# Plot PCA results with age gradient coloring
ggplot(pc_scores, aes(x = PC1, y = PC2, color = Age)) +
  geom_point() +
  labs(title = "PCA of DCM Patients by Age Group",
       x = "PC1",
       y = "PC2") +
  scale_color_viridis_c(option = "A", direction = -1) + # Use blueish turquoise color scale
  theme_minimal() +
  theme(legend.position = "top")

