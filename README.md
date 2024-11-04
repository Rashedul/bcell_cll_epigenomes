## B-cell and CLL epigenomes
Source code of integrative B-cell and CLL epigenomes 

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

On Linux/Unix terminal, please run the command to creat a R environment and to install R packages. This will create R environment named `renv` and install all packages within this environment. Note that you need to have conda installed. 

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

On Linux/Unix terminal, please run the command to generate plots. This will generate .md and .html files containing scripts and plots. 

```
cd docs
Rscript -e "rmarkdown::render('bcell_cll_epigenomes_plots.Rmd', output_format = 'all')"
```

### 3. Pseudocode for data processing
- Please find the pseudocodes for data processing pipelines [here](https://github.com/Rashedul/bcell_cll_epigenomes/blob/main/docs/pseudocode.md)

### 4. diffER pipeline 

diffER identifies differential enrichment of ChIP-seq peaks between two groups.

- Please find the diffER pipeline [here](https://github.com/Rashedul/diffER)

### 5. CRIS pipeline 

CRIS uses RNA-seq data to detect perccentage of somatic hypermutations of IGHV genes to classify CLL cases.

- Please find the CRIS pipeline [here](https://github.com/Rashedul/CRIS)

### 6. Citation 