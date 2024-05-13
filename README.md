**Overview:**

This repository contains code for analyzing RNA expression data related to cardiomyopathy. It includes scripts to compile the sample data into a singular data frame, assign phenotypic annotations to the data, PCA and UMAPs that labels the samples by disease etiology, race, and other phenotypic factors. Most importantly, it contains the code to perform DESeq2 and edgeR differential gene analysis in python and R, respectively.

**Data:**

The gene expression data were generated using RNA sequencing techniques from samples collected from individuals with various forms of cardiomyopathy. The RNA-seq libraries were prepared with the Illumina TruSEq stranded mrNA kit and Ovation amplification kit. The original data can be found publicly in the NCBI GEO repository under the accession number GSE141910 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141910).

**Folder Structure:**

- **Data:** Contains the RNA-seq data files from each sample.
- **Scripts:** Contains R and Python scripts for preprocessing, organizing annotations, compiling data, dimensionality-reduction plots, and differential gene expression analysis.
- **Results:** Contains outputs generated from the scripts.

**Installation:**

To run the code, you will need R and python installed on your system along with the following R packages:
- limma
- ggplot2
- dplyr
- readr
- tidyr
- edgeR

You can install these packages using the following R commands:

```R
install.packages("limma")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("readr")
install.packages("tidyr")
```
These scripts were generated and ran with R version 4.1.2 (Bird Hippie).

After installing the required packages, you can clone or download this repository and run the R scripts provided in the "Scripts" folder. Make sure to update file paths if necessary to point to the location of your data files. To generate the PCA plot, run the scripts in the "Scripts" folder in the following order:

1. combinedsamples.R
2. justtheannotations.R
3. combinedwithannotations.R
4. PCAattemptEtiology.R

Python Packages Used:
- pandas (version 1.5.3): Used for data manipulation and analysis, providing DataFrame structures to handle structured data efficiently.
- numpy (version 1.24.3): Essential for numerical operations, providing arrays and mathematical functions for efficient computation.
- google.colab (python version 3.10.12)
- condacolab (python version 3.10.12)
- pyDESeq2 (version 0.4.8): Used for differential expression analysis of RNA-seq data, performing statistical tests to identify significant gene expression changes.
- os (built-in)
- pickle (built-in)
- statsmodels (version 0.14.0): Provides statistical modeling and testing functionalities, used here for statistical analysis and hypothesis testing.
- matplotlib (version 3.7.1): Used for data visualization, particularly for creating plots like scatter plots and volcano plots to visualize gene expression and statistical results.




