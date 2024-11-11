## Integrative analysis of aberrant epigenomic landscape in chronic lymphocytic leukemia
The epigenome of chronic lymphocytic leukemia (CLL) is characterized by the accumulation of epigenetic aberrations. Annotating CLL-specific epigenetic signatures and those that reflect normal B-cell development remains challenging due to heterogeneity within B-cell populations and CLL cases. Here, we analyzed 50 reference epigenomes of CLL and B-cell subtypes to map the interconnected layers of epigenetic aberrations in CLL. This repository contains the source code for the integrative analysis of B-cell and CLL epigenomes. 

### 1. Rscripts for figures 
- Find the scripts and plots [here](https://rashedul.github.io/bcell_cll_epigenomes/)

### 2. Generate figures

#### Requirements

- Ubuntu >= 22.04.4 LTS  
- R version >= 4.3.1
- pandoc >= 2.9.2.1
- conda >= 24.7.1
- RAM ~10GB 
- Runtime ~10 minutes 

#### Create R environment and install R packages

In a Linux/Unix terminal, execute the following commands to create an R environment (`renv`) and install the necessary R packages. 

```
# Install R packages from conda
conda create -n renv r-base=4.3 -y
conda activate renv
conda install -c conda-forge -y r-ggplot2 r-dplyr r-pandoc r-tidyverse r-reshape2 r-factoextra r-pheatmap r-UpSetR r-data.table r-survival r-survminer r-readxl r-ggrepel r-patchwork r-matrixStats r-codetools

# Install R packages outside of conda
R -e "install.packages('circlize', repos='http://cran.r-project.org', dependencies = TRUE)"
```

#### Download the project repository and run Rscript to generate plots

In a terminal, run the command to generate plots. This will create .md and .html files that include the scripts and resulting plots.

```
# Clone the repo
git clone https://github.com/Rashedul/bcell_cll_epigenomes.git
cd ./bcell_cll_epigenomes/docs

Rscript -e "rmarkdown::render('index.Rmd', output_format = 'all')"
```

#### Alternatively, use RStudio to knit the R Markdown file on Windows and macOS
 
To generate the output document from the R Markdown file `index.Rmd`, 

Step 1: Install listed packages through RStudio.

Step 2: Knit the document by clicking the Knit button in the toolbar of RStudio. This will process the `index.Rmd` file and generate the output documents (such as HTML, MD).

### 3. Pseudocode for data processing
- Find the pseudocodes for data processing pipelines [here](https://github.com/Rashedul/bcell_cll_epigenomes/blob/main/docs/pseudocode.md)

### 4. diffER pipeline 

diffER identifies differential enrichment of ChIP-seq peaks between two groups.

- Find the diffER pipeline [here](https://github.com/Rashedul/diffER)

### 5. CRIS pipeline 

CRIS uses RNA-seq data to detect perccentage of somatic hypermutations of IGHV genes to classify CLL cases into uCLL and mCLL.

- Find the CRIS pipeline [here](https://github.com/Rashedul/CRIS)

### 6. Citation 
Islam R. et al., Integrative analysis of aberrant epigenomic landscape in chronic lymphocytic leukemia, (2024).