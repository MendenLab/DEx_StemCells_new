# Differential expression to investigate ferropstosis and antioxidant effect on stem cells

## Description

In this repository is presented a script with the necessary steps to replicate the differential expression analysis and visualizations of a time and compound dependent dataset of stem cells. 

## Usage 

The script considers a data format of raw RNA-seq count data with IDs and genes represented in columns and rows, respectively.

Please check the manuscript to access the data used for the study.

This project was conducted in [R software](https://www.r-project.org). 

All the necessary R package dependencies are

* ggplot2
* ggpubr
* readxl
* readr
* dplyr
* tidyverse
* plyr
* ggrepel

These packages and dependencies should be installed a priori using the install.packages() function and complemented by the library() function to be ready to use. 

Furthermore, we also leveraged the DESeq2 package available in the Bioconductor.

Similarly, to employ this, first install the BiocManager package using install.packages("BiocManager"), and later the packages above using BiocManager::install(). 


## Contact

For any inquiries related to this work, please contact me via e-mail ana.galhoz@helmholtz-muenchen.de
