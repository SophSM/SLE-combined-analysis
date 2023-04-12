# SLE-combined-analysis

This is a meta-analysis on expression data conducted on different SLE patients versus control experiments that were found in the Gene Expression Omnibus (GEO) database. This meta-analysis contains the pipeline used for uniform processing of bulk RNA-seq expression data, using the `monorail external` software, which was used for uniformly process data in the Recount3 project. (See the `monorail` folder README file for further information about the pipeline).

The differential expression analysis for the raw count matrices derived from the previous `monorail external` pipeline is available in `combined_counts.Rmd`, were counts from all studies were combined into a single count matrix.

Functional enrichment analysis, differential expression visualization and further exploring of the data was performed and scripts are available for the plots generated at the `scripts` folder. All code the code used while exploring the data is also all together available in the markdown `combined-meta-analysis.Rmd`.