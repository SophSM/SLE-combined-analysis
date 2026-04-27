## Weighted gene co-expression analysis with WGCNA

This directory contains the scripts used for the co-expression analysis on the SLE samples.

The order of the scripts goes as follows:

`1-reducedVSD.R`: Makes a VSD count matrix with only SLE samples and desired transcript annotations

- **Input**: `vsd2.RData`, `LRT-dds.RData`
- **Output**: `named-sle.vsd.RData`

`2-softPowersWGCNA.R`: Clean samples and choose soft-power threshold for module detection.

- **Input**: `named-sle.vsd.RData`
- **Output**: `less_samplesDatExpr.RData`

`3-module-detection.R`: Compute adjacency and TOM matrices, perform module detection with previously chosen soft power.

- **Input**: `less_samplesDatExpr.RData`
- **Output**: `named-reduced-sle-TOM-WGCNA.RData`
