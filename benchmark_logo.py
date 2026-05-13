#!/usr/bin/env python3
"""
Cross-TF Generalization Benchmark
=================================
Evaluate how well perturbation response predictors generalize across
different transcription factor categories using Leave-One-Gene-Out (LOGO) protocol.

Dataset: Dixit et al. 2016 CRISPR Perturb-seq (K562 TFs, 10 genes)
Paper: https://doi.org/10.1016/j.cels.2022.01.003

Requirements:
    pip install scanpy scikit-learn scipy pandas numpy matplotlib seaborn

Usage:
    python benchmark_logo.py

Author: Junbiao Xue
"""

import os
import sys
import json
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge as RidgeRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from sklearn.metrics.pairwise import cosine_similarity
from scipy.stats import pearsonr, spearmanr, ttest_rel
from scipy.sparse import issparse
import matplotlib.pyplot as plt
import seaborn as sns

# Configuration
DATA_PATH = "DixitRegev2016_K562_TFs_13_days.h5ad"
N_HVG = 2000
RANDOM_STATE = 42

# 10 target transcription factors
TARGET_GENES = [
    'CREB1', 'E2F4', 'EGR1', 'ELF1', 'ELK1',
    'ETS1', 'GABPA', 'IRF1', 'NR2C2', 'YY1'
]


def load_and_preprocess(data_path):
    """Load and preprocess the Perturb-seq dataset."""
    print("Loading data...")
    adata = sc.read(data_path)
    
    # Filter to valid target genes + control
    valid_genes = TARGET_GENES + ['control']
    mask = adata.obs['target'].isin(TARGET_GENES) | (adata.obs['perturbation'] == 'control')
    adata = adata[mask].copy()
    
    # Quality control
    adata = adata[
        (adata.obs['ncounts'] >= 200) &
        (adata.obs['ncounts'] <= 10000) &
        (adata.obs['percent_mito'] < 0.20)
    ].copy()
    
    # Library size normalization
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    print(f"  Filtered: {adata.shape[0]} cells, {adata.shape[1]} genes")
    print(f"  Control cells: {(adata.obs['perturbation'] == 'control').sum()}")
    
    return adata


def compute_deltas(adata, target_genes, hvg_indices=None):
    """Compute perturbation delta vectors for each target gene."""
    ctrl_mask = adata.obs['perturbation'] == 'control'
    ctrl_idx = np.where(ctrl_mask)[0]
    ctrl_mean = adata.X[ctrl_idx].mean(axis=0)
    if issparse(ctrl_mean):
        ctrl_mean = ctrl_mean.toarray().flatten()
    
    deltas = {}
    for gene in target_genes:
        gene_mask = adata.obs['target'] == gene
        gene_idx = np.where(gene_mask)[0]
        gene_mean = adata.X[gene_idx].mean(axis=0)
        if issparse(gene_mean):
            gene_mean = gene_mean.toarray().flatten()
        
        delta = gene_mean - ctrl_mean
        
        if hvg_indices is not None:
            delta = delta[hvg_indices]
        
        deltas[gene] = {
            'delta': delta,
            'n_cells': len(gene_idx),
            'norm': np.linalg.norm(delta)
        }
        print(f"  {gene}: n={len(gene_idx)}, ||delta||={deltas[gene]['norm']:.2f}")
    
    return deltas


def logo_evaluate(deltas, adata, n_hvg=N_HVG):
    """
    Leave-One-Gene-Out evaluation.
    
    For each held-out TF H:
    1. Select HVGs using training TF cells only (per-fold)
    2. Compute deltas for training TFs
    3. Train predictors on 9 training deltas
    4. Evaluate on held-out TF
    """
    results = {gene: {} for gene in TARGET_GENES}
    
    for held_out in TARGET_GENES:
        train_genes = [g for g in TARGET_GENES if g != held_out]
        train_mask = adata.obs['target'].isin(train_genes)
        
        # Per-fold HVG selection using training TFs only
        train_adata = adata[train_mask]
        sc.pp.highly_variable_genes(
            train_adata, 
            n_top_genes=n_hvg,
            flavor='seurat_v3'
        )
        hvg_mask = train_adata.var['highly_variable'].values
        hvg_indices = np.where(hvg_mask)[0]
        
        # Compute training deltas with per-fold HVGs
        train_deltas = []
        ctrl_train_idx = np.where(train_adata.obs['perturbation'] == 'control')[0]
        ctrl_train_mean = train_adata.X[ctrl_train_idx].mean(axis=0)
        if issparse(ctrl_train_mean):
            ctrl_train_mean = ctrl_train_mean.toarray().flatten()
        
        for gene in train_genes:
            gene_idx = np.where(train_adata.obs['target'].values == gene)[0]
            gene_mean = train_adata.X[gene_idx].mean(axis=0)
            if issparse(gene_mean):
                gene_mean = gene_mean.toarray().flatten()
            delta = gene_mean - ctrl_train_mean
            train_deltas.append(delta[hvg_mask])
        
        train_deltas = np.array(train_deltas)  # 9 x n_hvg
        
        # True delta for held-out (with same HVGs)
        true_delta = deltas[held_out]['delta'][hvg_mask]
        
        # SCP: mean of training deltas
        scp_pred = train_deltas.mean(axis=0)
        scp_pearson = pearsonr(scp_pred, true_delta)[0]
        scp_cosine = cosine_similarity([scp_pred], [true_delta])[0, 0]
        scp_rmse = np.sqrt(mean_squared_error(scp_pred, true_delta))
        
        # Ridge Regression
        ridge = RidgeRegression(alpha=1.0, fit_intercept=False, random_state=RANDOM_STATE)
        ridge.fit(train_deltas, train_deltas.mean(axis=1))
        ridge_pred = ridge.predict(train_deltas.mean(axis=0).reshape(1, -1)).flatten()
        ridge_pearson = pearsonr(ridge_pred, true_delta)[0]
        ridge_cosine = cosine_similarity([ridge_pred], [true_delta])[0, 0]
        ridge_rmse = np.sqrt(mean_squared_error(ridge_pred, true_delta))
        
        # Random Forest
        rf = RandomForestRegressor(
            n_estimators=100, max_depth=5, 
            random_state=RANDOM_STATE
        )
        rf.fit(train_deltas, train_deltas.mean(axis=1))
        rf_pred = rf.predict(train_deltas.mean(axis=0).reshape(1, -1)).flatten()
        rf_pearson = pearsonr(rf_pred, true_delta)[0]
        rf_cosine = cosine_similarity([rf_pred], [true_delta])[0, 0]
        rf_rmse = np.sqrt(mean_squared_error(rf_pred, true_delta))
        
        # PCA-SCP: PCA on training deltas (max n_components = n_samples - 1)
        n_components = min(len(train_genes) - 1, 5)
        pca = PCA(n_components=n_components, random_state=RANDOM_STATE)
        train_pca = pca.fit_transform(train_deltas)
        mean_pca = train_pca.mean(axis=0)
        mean_recon = pca.inverse_transform(mean_pca.reshape(1, -1)).flatten()
        pca_pred = mean_recon
        pca_pearson = pearsonr(pca_pred, true_delta)[0]
        pca_cosine = cosine_similarity([pca_pred], [true_delta])[0, 0]
        pca_rmse = np.sqrt(mean_squared_error(pca_pred, true_delta))
        
        results[held_out] = {
            'SCP': {'pearson': scp_pearson, 'cosine': scp_cosine, 'rmse': scp_rmse},
            'Ridge': {'pearson': ridge_pearson, 'cosine': ridge_cosine, 'rmse': ridge_rmse},
            'RF': {'pearson': rf_pearson, 'cosine': rf_cosine, 'rmse': rf_rmse},
            'PCA_SCP': {'pearson': pca_pearson, 'cosine': pca_cosine, 'rmse': pca_rmse},
            'delta_norm': deltas[held_out]['norm'],
            'n_cells': deltas[held_out]['n_cells']
        }
        
        print(f"{held_out}: SCP={scp_pearson:.3f}, Ridge={ridge_pearson:.3f}, "
              f"RF={rf_pearson:.3f}, PCA-SCP={pca_pearson:.3f}")
    
    return results


def compute_statistics(results):
    """Compute aggregate statistics and statistical tests."""
    stats = {}
    
    for method in ['SCP', 'Ridge', 'RF', 'PCA_SCP']:
        pearsons = [results[tf][method]['pearson'] for tf in TARGET_GENES]
        stats[method] = {
            'mean': np.mean(pearsons),
            'std': np.std(pearsons),
            'min': np.min(pearsons),
            'max': np.max(pearsons),
            'pearsons': pearsons
        }
    
    # Paired t-test: Ridge vs SCP
    ridge_ps = np.array([results[tf]['Ridge']['pearson'] for tf in TARGET_GENES])
    scp_ps = np.array([results[tf]['SCP']['pearson'] for tf in TARGET_GENES])
    t_stat, t_pval = ttest_rel(ridge_ps, scp_ps)
    
    # Bootstrap 95% CI for Ridge - SCP difference
    n_bootstrap = 10000
    np.random.seed(RANDOM_STATE)
    diffs = []
    for _ in range(n_bootstrap):
        idx = np.random.choice(len(TARGET_GENES), len(TARGET_GENES), replace=True)
        diffs.append(ridge_ps[idx].mean() - scp_ps[idx].mean())
    diffs = np.array(diffs)
    ci_low, ci_high = np.percentile(diffs, [2.5, 97.5])
    
    # Permutation test (sign-flip)
    obs_diff = (ridge_ps - scp_ps).mean()
    perm_diffs = []
    for _ in range(n_bootstrap):
        signs = np.random.choice([-1, 1], len(TARGET_GENES))
        perm_diffs.append(((ridge_ps - scp_ps) * signs).mean())
    perm_diffs = np.array(perm_diffs)
    perm_pval = np.mean(np.abs(perm_diffs) >= np.abs(obs_diff))
    
    stats['comparison'] = {
        'ridge_scp_diff': obs_diff,
        't_stat': t_stat,
        't_pval': t_pval,
        'bootstrap_ci': (ci_low, ci_high),
        'permutation_pval': perm_pval
    }
    
    # Effect size vs predictability correlation
    norms = np.array([results[tf]['delta_norm'] for tf in TARGET_GENES])
    effect_corr, effect_pval = pearsonr(norms, ridge_ps)
    stats['effect_predictability'] = {
        'correlation': effect_corr,
        'pval': effect_pval
    }
    
    return stats


def create_results_dataframe(results):
    """Create a DataFrame of per-TF results."""
    rows = []
    for tf in TARGET_GENES:
        row = {'TF': tf}
        for method in ['SCP', 'Ridge', 'RF', 'PCA_SCP']:
            row[f'{method}_Pearson'] = results[tf][method]['pearson']
            row[f'{method}_Cosine'] = results[tf][method]['cosine']
            row[f'{method}_RMSE'] = results[tf][method]['rmse']
        row['delta_norm'] = results[tf]['delta_norm']
        row['n_cells'] = results[tf]['n_cells']
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    print("=" * 60)
    print("Cross-TF Generalization Benchmark")
    print("=" * 60)
    
    # Load data
    adata = load_and_preprocess(DATA_PATH)
    
    # Compute all deltas (for effect size analysis)
    print("\nComputing perturbation deltas...")
    deltas = compute_deltas(adata, TARGET_GENES)
    
    # LOGO evaluation
    print("\nRunning LOGO evaluation...")
    results = logo_evaluate(deltas, adata)
    
    # Compute statistics
    print("\nComputing statistics...")
    stats = compute_statistics(results)
    
    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for method in ['SCP', 'Ridge', 'RF', 'PCA_SCP']:
        print(f"{method}: mean r={stats[method]['mean']:.3f}, "
              f"std={stats[method]['std']:.3f}")
    
    print(f"\nRidge vs SCP:")
    print(f"  Mean difference: {stats['comparison']['ridge_scp_diff']:.3f}")
    print(f"  Paired t-test p: {stats['comparison']['t_pval']:.4f}")
    print(f"  Bootstrap 95% CI: [{stats['comparison']['bootstrap_ci'][0]:.3f}, "
          f"{stats['comparison']['bootstrap_ci'][1]:.3f}]")
    print(f"  Permutation p: {stats['comparison']['permutation_pval']:.4f}")
    
    print(f"\nEffect size vs predictability:")
    print(f"  Correlation: {stats['effect_predictability']['correlation']:.3f}")
    print(f"  p-value: {stats['effect_predictability']['pval']:.4f}")
    
    # Save results
    df = create_results_dataframe(results)
    df.to_csv("benchmark_results.csv", index=False)
    print(f"\nResults saved to benchmark_results.csv")
    
    # Save extended results
    df.to_csv("benchmark_results_extended.csv", index=False)
    
    return results, stats


if __name__ == "__main__":
    results, stats = main()
