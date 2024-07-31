# Comparative Gene Coexpression Analysis between Sugarcane and Sorghum

This repository contains scripts and data for the comparative gene coexpression analysis between sugarcane and sorghum. The project aims to identify and compare gene coexpression networks in these two important crop species.

## Table of Contents
- [Introduction](#introduction)
- [Data](#data)
- [Methods](#methods)
- [Results](#results)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgements](#acknowledgements)

## Introduction
Gene coexpression analysis is a powerful tool for understanding the functional relationships between genes. This project focuses on comparing gene coexpression networks in sugarcane and sorghum to identify conserved and species-specific coexpression patterns in the sugar-fiber accumulation proccess.

## Data
The data used in this analysis include:
- Gene expression datasets retrieved from SRA (https://www.ncbi.nlm.nih.gov/sra) for Sugarcane and Sorghum. for Sorghum  it includes 115 accessions and six genotypes, three grain and three sweet. For Sugarcane it includes 109 accessions including seven gentypes, three rich in sugar and four rich in fiber, metadata in data/sugarcane_metadata.csv and data/sorghum_metadata.csv

- Sorghum Pantranscriptome (Transcript of longest cds per orthogroup):
- Sugarcane Pantranscriptome (Transcript of longest cds per orthogroup):


| Genotype | Number of Accessions | Group | Species |
|----------|----------------------|-------|---------|
| BTx623   | 20                   |sweet  |Sorghum  |
| RTx430   | 20                   |sweet  |Sorghum  |
| BTx642   | 20                   |sweet  |Sorghum  |
| della    | 23                   |grain  |Sorghum  |
| keller   | 20                   |grain  |Sorghum  |
| rio      | 12                   |grain  |Sorghum  |
| SRA5     | 15                   |sugar  |Sugarcane|
| SRA3     | 15                   |sugar  |Sugarcane|
| Q238     | 15                   |sugar  |Sugarcane|
| KQB09-20432| 19                 |fiber  |Sugarcane|
| Q157     | 15                   |fiber  |Sugarcane|
| Q151     | 15                   |fiber  |Sugarcane|
| Q135     | 15                   |fiber  |Sugarcane|

complete list of accession in accessions directory

## Methods
The analysis was performed using the following steps:
1. Preprocessing of gene expression data.
    - Download sequences using sra-toolkit (sra-tools 2.10.0)
    - bbduk: remove poor quality sequences short sequences and remove rRNA (bbmap 38.84)
    - Read mapping salmon (v1.10.2) in Sorghum and Sugarcane pantranscriptomes references (transcript of longest cds per orthogroup) ( salmon 0.12.0)

There is a snakefile for this steps in snakefile_preprocessing.py, download_sra.sh launches snakefile in SLURM type system.
get_results.sh is useful to join the salmon data
plot_quant.r is useful to take the results of get_results.sh and plot the quantfification data

2. Data normalization
    - Load quantification dsalmon (v1.10.2)ata and plot exploratory PCA of vst transformed data (pca_for_rnaseq.r)
    - Correct for batch effects using RUVr from RUVseq package (RUV.r) (https://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.html)

3. Construction of gene coexpression networks.
    - filter low expressed data
 
4. Module detection and correlation analysis.
5. Comparative analysis of the networks to identify conserved and unique coexpression modules.

## Results

We download a total of 14122223702 raw sequences, 6651867490 from sorghum and 7470356212 from sugarcane, after removing poor quality, rRNA, very short sequences and mapping with salmon against pan-transcriptome references we got 10698352940 sequences (75.6% of total download sequences, 84.5% sorghum and 68.0% for sugarcane), see table and plots in bbduk results. 

Sorghum quantification
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/quant/quant_sorghum.png?raw=true)

Sugarcane quantification
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/quant/quant_sugarcane.png?raw=true)

Raw PCA for sorghum
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/PCA/pca_sorghum_group.png?raw=true)

Raw PCA For sugarcane
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/PCA/pca_sugarcane_group.png?raw=true)

Raw facet PCA for sugarcane
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/PCA/facet_pca_sugarcane_repetition_group.png?raw=true)

Normalized PCA for sorghum - k1
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/RUVr/sorghum/k1_RUVr_groups.png?raw=true)

Normalized PCA for sugarcane - k2
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/RUVr/cane/k2_RUVr_groups.png?raw=true)

Count matrixes: Matrices directory 
## Installation
You can find yml configuration files for conda envs in conda_envs directory
