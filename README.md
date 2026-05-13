# Cross-TF Generalization Benchmark

**Evaluating cross-gene-category generalization of single-cell perturbation response predictors**

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)

## Overview

This repository contains the analysis code for evaluating how perturbation response predictors generalize across different transcription factor (TF) categories using a Leave-One-Gene-Out (LOGO) evaluation protocol.

### Key Findings (Preliminary)

- **Ridge Regression marginally outperforms SCP** (mean Pearson r = 0.455 vs 0.387), but the difference is **not statistically significant** (bootstrap 95% CI: [-0.01, 0.14], permutation p = 0.09).
- **PCA-SCP fails entirely** (mean r ≈ 0), suggesting linear dimensionality reduction of delta vectors is not helpful with n << p.
- **Perturbation effect size correlates with predictability** (r = 0.51, p = 0.13), though this is an exploratory observation with n=10 TFs.
- **ETS-family TFs show correlated responses** (ETS1-GABPA: r = 0.78), but this requires validation on larger TF panels.

> **Note:** This is a small-scale pilot study (n=10 TFs). Findings should be treated as hypothesis-generating rather than definitive.

## Dataset

- **Dixit et al. 2016 CRISPR Perturb-seq** (K562 cells, 10 TF knockdowns)
- Source: [Zenodo scPerturb](https://zenodo.org/record/13350497)
- Cells: 4,639 (3,643 perturbed + 996 control)
- Features: 2,000 highly variable genes

## Installation

```bash
pip install scanpy scikit-learn scipy pandas numpy matplotlib seaborn
```

## Usage

```bash
# Download the data (121 MB) from Zenodo
# https://zenodo.org/record/13350497

# Run the benchmark
python benchmark_logo.py
```

## Methods Evaluated

| Method | Description |
|--------|-------------|
| SCP | Simple Control Prediction: mean of all training deltas |
| Ridge Regression | L2-regularized linear regression |
| Random Forest | Ensemble of 100 regression trees |
| PCA-SCP | SCP after PCA dimensionality reduction |

## Evaluation Protocol

Leave-One-Gene-Out (LOGO):
1. For each held-out TF H (10 total):
   - Train on remaining 9 TFs' perturbation deltas
   - Predict delta for H
   - Evaluate using Pearson correlation

## Results

See `benchmark_results_extended.csv` for per-TF results across all methods.

```
Method      | Mean r |  SD
------------|--------|-------
SCP         |  0.387 | 0.093
Ridge       |  0.455 | 0.094
Random Forest|  0.402 | 0.092
PCA-SCP     | -0.001 | 0.032
```

## Limitations

1. **Small sample size** (n=10 TFs) limits statistical power
2. **Single cell type** (K562) — may not generalize to other cell types
3. **TF knockdown only** — CRISPRa, CRISPRko, chemical perturbations untested
4. **No deep learning baselines** — scGen, CPA, scFoundation not evaluated
5. **Mean delta only** — cell-level heterogeneity not evaluated

## Citation

If you use this code, please cite:

```
Xue J. et al. (2025). A Pilot Evaluation of Cross-Transcription-Factor
Generalization in Single-Cell Perturbation Response Prediction.
GitHub: https://github.com/JunbiaoXue/tf-perturbation-benchmark
```

## References

- Dixit, A. et al. (2016). Perturb-seq: dissecting molecular circuits with scalable single cell RNA profiling of pooled genetic screens. *Cell* 167, 1853-1866.
- Weinreb, C. et al. (2018). SCAPE: single-cell analysis of progenitor cell identity. *bioRxiv*.
- Lotfollahi, M. et al. (2019). scGen predicts single-cell perturbation responses. *Nature Methods* 16, 715-721.
- Lotfollahi, M. et al. (2023). Predicting cellular responses to complex perturbations. *Molecular Systems Biology* 18, e11517.
- Zheng, G.R. et al. (2023). Learning single-cell perturbation responses using neural optimal transport. *Cell Systems* 14, 978-994.

## Author

**Junbiao Xue**  
GitHub: [@JunbiaoXue](https://github.com/JunbiaoXue)

## License

MIT License
