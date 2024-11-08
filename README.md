## Integrative analysis of aberrant epigenomic landscape in chronic lymphocytic leukemia
The epigenome of chronic lymphocytic leukemia (CLL) is characterized by the accumulation of epigenetic aberrations. Annotating CLL-specific epigenetic signatures and those that reflect normal B-cell development remains challenging due to heterogeneity within B-cell populations and CLL cases. Here, we analyzed 50 reference epigenomes of CLL and B-cell subtypes to map the interconnected layers of epigenetic aberrations in CLL. This repository contains the source code for the integrative analysis of B-cell and CLL epigenomes. 

### 1. Rscripts for figures 
- Please find the scripts and plots [here](https://rashedul.github.io/bcell_cll_epigenomes/)

### 2. Generate figures

#### Requirements

- Ubuntu >= 22.04.4 LTS or macOS >= 13.7 
- R version >= 4.1.2
- pandoc >= 2.9.2.1
- conda >= 24.7.1
- RAM ~32GB 
- Runtime ~60 minutes 

#### Create R environment and install R packages

In a Linux/Unix terminal, execute the following commands to create an R environment (`renv`) and install the necessary R packages. Make sure that conda is already installed.

```
# Install packages from conda
conda create -n renv r-base -y
conda activate renv
conda install -c conda-forge -y r-ggplot2 r-dplyr

# Install packages outside of conda

# did not install dependencies
R -e "install.packages('tidyverse', repos='http://cran.r-project.org', dependencies = TRUE)"
```

#### Run Rscript to generate plots

In a Linux/Unix terminal, run the command to generate plots. This will create .md and .html files that include the scripts and resulting plots.

```
# Clone the repo
git clone https://github.com/Rashedul/bcell_cll_epigenomes.git
cd ./bcell_cll_epigenomes/docs

Rscript -e "rmarkdown::render('index.Rmd', output_format = 'all')"
```

### 3. Pseudocode for data processing
- Find the pseudocodes for data processing pipelines [here](https://github.com/Rashedul/bcell_cll_epigenomes/blob/main/docs/pseudocode.md)

### 4. diffER pipeline 

diffER identifies differential enrichment of ChIP-seq peaks between two groups.

- Please find the diffER pipeline [here](https://github.com/Rashedul/diffER)

### 5. CRIS pipeline 

CRIS uses RNA-seq data to detect perccentage of somatic hypermutations of IGHV genes to classify CLL cases into uCLL and mCLL.

- Please find the CRIS pipeline [here](https://github.com/Rashedul/CRIS)

### 6. Citation 
Islam R. et al., Integrative analysis of aberrant epigenomic landscape in chronic lymphocytic leukemia, (2024).