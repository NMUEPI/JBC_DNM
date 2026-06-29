# Analysis of De Novo Mutations from Trio-based Whole-genome Sequencing

This repository contains R scripts used to analyze de novo mutations (DNMs) identified from trio-based whole-genome sequencing data. The analyses focus on parental-origin DNMs, embryonic post-zygotic mutations (EPZMs), parental age, assisted reproductive technology (ART), and mutation-spectrum differences.

## Overview

The scripts generate summary tables, regression results, and manuscript figures for DNM analyses. Most analyses are implemented in R using Fisher's exact tests, linear regression, restricted cubic spline models, and formatted forest plots.

## Main analysis topics include:

- Mutational spectrum comparison among paternal DNMs, maternal DNMs, and EPZMs
- Trinucleotide-context enrichment analysis
- Functional annotation distribution of DNMs
- Parental-age associations with DNM burden
- Mutation-type-specific parental-age effects
- ART-associated differences in DNM burden
- Forest plots for regression summary tables

## Files

| File | Description |
| --- | --- |
| `DNM_analysis_final.r` | Main analysis script. Reads DNM datasets, generates figures and figure-related results. |
| `DNM_mutclass_barplot.r` | Generates Figure 1B: comparison of mutation-class proportions among paternal DNMs, maternal DNMs, and EPZMs using Fisher's exact tests. |
| `DNM_mutclass_heatmap.r` | Generates Figure 1C: trinucleotide-context comparison and enrichment heatmap for parental-origin DNMs and EPZMs. |
| `age_linear_fig.r` | Generates Figure 2A: linear association between parental age at conception and standardized DNM counts. |
| `age_mutclass_fig.r` | Generates Figure 2B: mutation-type-specific maternal and paternal age effects. |
| `age_RCS_fig.r` | Generates Figure 2C: restricted cubic spline analysis of maternal age and maternal-origin dnSNVs. |
| `ART_boxplot.r` | Generates Figure 3A: violin/box plots comparing standardized dnSNV burden between ART and naturally conceived groups. |
| `forest_fig.r` | Provides helper functions to prepare and draw forest plots from regression summary tables. |

## Usage

Run the main analysis script from the project analysis directory:
```r
source("DNM_analysis_final.r")
```

Individual figure scripts can be run separately after the required input objects and helper functions are loaded:
```r
source("DNM_mutclass_barplot.r")
source("DNM_mutclass_heatmap.r")
source("age_linear_fig.r")
source("age_mutclass_fig.r")
source("age_RCS_fig.r")
source("ART_boxplot.r")
source("forest_fig.r")
```
