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
    - Read mapping salmon (v1.10.2) in Sorghum and Sugarcane pantranscriptomes references (transcript of longest cds per orthogroup)

There is a snakefile for this steps in snakefile_preprocessing.py, run_snakemake.sh launches snakefile in SLURM type system.
plot_quant.r is useful to take the results of get_results.sh and plot the quantfification data
get_bbduk_results.sh is useful to join bbduk results
plot_bbduk can plot the results from get_bbduk_results.sh



2. Quantification
    - salmon (v1.10.2) and plot exploratory PCA of vst transformed data (pca_for_rnaseq.r)
    - get_results.sh collect salmon results
    - plot_quant.r plot salmon results

3. Plot PCA
    - pca_for_rnaseq.r plot PCA using DESeq2

4. Correct for batch effects
    - RUV.r implements RUVr (https://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.html) method and plot PCA  	 

5. Make network
    - pearson_cor.r ultra fast pearson correlation method in R using HiClimR (https://cran.r-project.org/web/packages/HiClimR/index.html)
    - histogram_cor.py make plot of the distribution of the correlations
    - very_fancy_python_code_for_melt.py select the correlations greater than 0.6 and makes an edge list from the pearson correlations matrix (multicore, very fast and memory efficient)

6. Find modules in network
    - run_mcl.sh run mcl (i=2.7) edge list generated by very_fancy_python_code_for_melt.py
    - format_modules.sh format mcl output in a convinient way (gene,module)

7. Correlation analysis of modules
    - corr_modules.r calculate the correlation of the eigengene of each module with the sugar-fiber trait

8. Annotate data
    - run_Panzzer2.sh assing GO to proteins
    - run_tf.sh run hmmer to find TF proteind domais (ref) pfam 37
    - assign_family_membership.pl make table with tf annotations from run_tf.sh (ref)

9. Get orthologues
    - run_orthofinder.sh find othologues with orthofinder (2.5.4) on sugcarcane, sorghum, maize and rice proteins.

10. Make heatmaps of expression data
    - make_heatmaps_for_all_modules.r make heatmaps for all modules and for module eigengenes correlated to sugar-fiber trait
    - make heatmap of orthologues genes in modules correlated to sugar-fiber trait

11. Enrichment GO analysis
    - topGO.r makes overrepresentation test on each module
    - compare_modules.r make overrepresentation analyisis of groups genes of modules correlated to sugar-fiber trait

12. Comparatve transcription factor analysis
    - comparative_tf.r detect transcription factors correlated to sugar-fiber trait and see if familias of tf conserve expression pattern through species.

## Results

1. Preprocessing of gene expression data

    - We downloaded a total sequences 13043204134 from sorghum and 7470356212 from sugarcane, after removing poor quality, rRNA, very short sequences and mapping with salmon against pan-transcriptome references we got 90.66% sorghum and 82.9% for sugarcane. See full results in results/summary_cleaning.csv and results/bbduk/ . 

2. Quantification

    We identified 600891 genes for sugarcane and 327598 for sorghum

    - See results/quant/pretty_table_cane.csv and results/quant/pretty_table_sorghum.csv

Sorghum quantification
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/quant/quant_sorghum.png?raw=true)

Sugarcane quantification
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/quant/quant_cane.png?raw=true)

3. Plot PCA for quant data
PCA for sorghum
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/PCA/pca_sorghum_group.png?raw=true)

 PCA For sugarcane
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/PCA/pca_sugarcane_group.png?raw=true)

facet PCA for sugarcane
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/PCA/pca_sugarcane_facet.png?raw=true)

4. Correct for batch effects (RUVr)
 PCA for sorghum - k1
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/RUVr/sorghum_k1_RUVr_groups.png?raw=true)

 PCA for sugarcane - k2
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/RUVr/sugarcane_k2_RUVr_groups.png?raw=true)

after batch effects correction we ended up with 147418 genes in sugarcane and 76309 in sorghum

5. Make network

Correlations distribution in Sugarcane 	
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/networks/histogram_sugarcane_cor_.png?raw=true)

140833 vertices
29233350 edges
0.02105017 clustering coefficient (transitivity)
21 connected components

Correlations distribution in Sorghum
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/networks/histogram_sorghum_cor_.png?raw=true)

76047 vertices
58469192 edges
0.1591881 transitivity
1 connected components

6. Find modules

   - See modules in results/networks/abs_mcl_cane_p0.6_i2.7.mcl.formated.csv for sugarcane and results/networks/abc_mcl_p0.6_i2.7.mcl.formated.csvfor sorghum

7. Correlation analysis for modules
   - Table for correlated modules for sugarcane results/networks/abs_cor_sugarcane_sweet_rho80_padj001.csv and for sorghum results/correlation_analysis/abs_cor_sorghum_sweet_rho80_padj001.csv

8. Annotate data 
   - GO annotations: For sugarcane data/GO_sugarcane/GO_cane_formated.txt, for sorghum results/panzzer/GO_formated.txt
   - Transcription factors annotation: For sugarcane results/tf_annot/tf_class_sugarcane.txt and for sorghum results/tf_annot/tf_class_sorghum.txt

9. Orthofinder results
   - See results/orthofinder/Results_Aug30/Orthogroups/Orthogroups.tsv

10. Make heatmaps of expression data
 
Sorghum correated modules eigengene
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/correlation_analysis/abs/abs_heatmap_sorghum_eigen_cor_mod_sweet.png?raw=true)

Sugarcane correlated modules eigengene
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/correlation_analysis/abs/abs_heatmap_sugarcane_eigen_cor_mod_sweet.png?raw=true)

11. Enrichment analysis
   - Enrichment analysis for sorghum modules: results/GO_enrichment/abs/sugarcane/modules, for sorghum: /home/dmpachon/jorge/comparative_cane/results/GO_enrichment/abs/sorghum/modules

     group 1: low expressed in sugar genotpypes

     group 2: high expressed in sugar genotypes		
  

for gene group 1 in sorghum
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/GO_enrichment/abs/correlation/group1_sorghum_GO.png?raw=true)

for gene group 2 in sorghum
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/GO_enrichment/abs/correlation/group2_sorghum_GO.png?raw=true)

for gene group 1 in sugarcane
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/GO_enrichment/abs/correlation/group1_sugarcane_GO.png?raw=true)

for gene group 2 in sugarcane
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/GO_enrichment/abs/correlation/group2_sugarcane_GO.png?raw=true)

12. Comparative transcription factor analisys: expression patterns

Orthologues TAP correlation sign across species
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/correlation_analysis/stacked_bar_comparative_tf.png?raw=true)


Conserved families
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/correlation_analysis/conserved_tf_families.png?raw=true)

Not conserved families
![alt text](https://github.com/jomare1188/glowing-palm-tree/blob/master/results/correlation_analysis/not_conserved_tf_families.png?raw=true)

## Installation

You can find yml configuration files for conda envs in conda_envs directory
