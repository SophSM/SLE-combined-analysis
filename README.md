# Differential gene expression in Peripheral Blood Mononuclear Cells from people with Systemic Lupus Erythematosus

This is a meta-analysis on expression data conducted on different SLE patients versus control experiments that were found in the Gene Expression Omnibus (GEO) database. This repository consists of three parts: **Monorail pipeline, Differential expression and analysis & WGCNA co-expression analysis**

### 1- Monorail pipeline

The pipeline used for uniform processing of bulk RNA-seq expression data, using the `monorail external` software, which was used for uniformly process data in the Recount3 project.

The scripts used are located on the `monorail` folder, see for README file for further information about the pipeline.

### 2- Differential expression and analysis

The differential expression analysis for the raw count matrices derived from the previous `monorail external` pipeline, as well as the scripts used to analyze the data.

The scripts used for this part are in the `scripts/` folder and follow this order.

-   `1-getRSE.Rmd`: With this script, each study is retrieved from the Recount3 database and turned into a RangedSummarizedExperiment-class (rse) object.

**Output:** `downloads.xz`

-   `2-curateRSE.Rmd`: With this script, for each study a new column `DISEASE` is added in order to get uniform information about the disease status of each sample.

**Input:** `downloads.xz`

**Output:** `curated_rse.xz`

-   `3-diffExpr.Rmd`: In this script we get the curated RSE object for each experiment, get the counts and combine them into a single data frame, visualize raw counts with a PCA, perform differential expression analysis and obtain a vsd count matrix.

**Input:** `curated_rse.xz`

**Output:** `LRT-dds.RData`, `vsd2.RData`, `all_data.csv`

-   `4-annotateDGElist`: In this script we annotate the gene name for each transcript in the DGE list, and filter out non-coding transcripts. Then we include some non-coding transcripts.

**Input:** `LRT-dds.RData`

**Output:** `namedDGElist.RData`, `DGElist_withNames_noncoding.csv`

Some data that was used or output from these scripts can be found in the `data/` folder.

**The next scripts are used for data analysis, they were roughly used in this order but can be used in any order, their input is some of the past output files**

-   `5-plot_PCA.R`

-   `6-plot_VolcanoPlot.R`

-   `7-plot_barDGEs.R`

-   `8-plot_Heatmap_fig.R`

-   `9-plot_ViolinPlots.R`

-   `10-plot_GOenrichment.R`

-   `11-plot_HeatmapGO.R`

-   `12-plot_AcyclicGraph.R`

-   `plot_figures.R`: This script generates the final figures used in the manuscript, the figures are in the `final_figs/` folder

### WGCNA co-expression analysis

The scripts used for the co-expression analysis can be found in the `WGCNA/` folder.
