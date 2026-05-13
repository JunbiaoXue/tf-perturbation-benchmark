# 转录因子跨类别泛化基准测试

**利用 Leave-One-Gene-Out（LOGO）评估协议，评估单细胞扰动响应预测模型在不同转录因子类别间的泛化能力**

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)

## 简介

本仓库包含分析代码，用于评估扰动响应预测模型在不同转录因子（TF）类别间的泛化能力，采用 Leave-One-Gene-Out（LOGO）评估协议。

### 主要发现（初步）

- **Ridge 回归略优于 SCP**（平均 Pearson r = 0.455 vs 0.387），但差异**未达统计显著**（bootstrap 95% CI: [-0.01, 0.14]，置换检验 p = 0.09）。
- **PCA-SCP 完全失效**（平均 r ≈ 0），表明在 n << p 条件下，对 delta 向量进行线性降维没有帮助。
- **扰动效应大小与可预测性相关**（r = 0.51，p = 0.13），但这是探索性观察（n=10 TFs）。
- **ETS 家族 TF 表现出相关性**（ETS1-GABPA: r = 0.78），但需要更大 TF 面板验证。

> **注意：** 这是一项小规模试点研究（n=10 TFs）。结论应视为"产生假设"而非"定论"。

## 数据集

- **Dixit et al. 2016 CRISPR Perturb-seq**（K562 细胞，10 个 TF 敲降）
- 来源：[Zenodo scPerturb](https://zenodo.org/record/13350497)
- 细胞数：4,639（3,643 扰动 + 996 对照）
- 特征：2,000 个高变异基因

## 安装

```bash
pip install scanpy scikit-learn scipy pandas numpy matplotlib seaborn
```

## 使用方法

```bash
# 从 Zenodo 下载数据（121 MB）
# https://zenodo.org/record/13350497

# 运行基准测试
python benchmark_logo.py
```

## 评估的方法

| 方法 | 描述 |
|------|------|
| SCP | 简单对照预测：所有训练 delta 的均值 |
| Ridge 回归 | L2 正则化线性回归 |
| 随机森林 | 100 棵回归树的集成 |
| PCA-SCP | PCA 降维后的 SCP |

## 评估协议

Leave-One-Gene-Out（LOGO）：
1. 对每个待预测 TF H（共 10 个）：
   - 使用其余 9 个 TF 的扰动 delta 进行训练
   - 预测 H 的 delta
   - 使用 Pearson 相关系数评估

## 结果

详见 `benchmark_results_extended.csv`，包含所有方法在每个 TF 上的结果。

```
方法       | 平均 r | 标准差
----------|--------|-------
SCP       |  0.387 | 0.093
Ridge     |  0.455 | 0.094
随机森林   |  0.402 | 0.092
PCA-SCP   | -0.001 | 0.032
```

## 局限性

1. **样本量小**（n=10 TFs），统计力有限
2. **仅一种细胞系**（K562），可能不泛化至其他细胞类型
3. **仅 TF 敲降**，CRISPRa、CRISPRko、化学扰动未测试
4. **未评估深度学习基线**（scGen、CPA、scFoundation 等）
5. **仅预测平均 delta**，未评估细胞水平异质性

## 引用

如使用本代码，请引用：

```
Xue J. et al. (2025). A Pilot Evaluation of Cross-Transcription-Factor
Generalization in Single-Cell Perturbation Response Prediction.
GitHub: https://github.com/JunbiaoXue/tf-perturbation-benchmark
```

## 参考文献

- Dixit, A. et al. (2016). Perturb-seq: dissecting molecular circuits with scalable single cell RNA profiling of pooled genetic screens. *Cell* 167, 1853-1866.
- Weinreb, C. et al. (2018). SCAPE: single-cell analysis of progenitor cell identity. *bioRxiv*.
- Lotfollahi, M. et al. (2019). scGen predicts single-cell perturbation responses. *Nature Methods* 16, 715-721.
- Lotfollahi, M. et al. (2023). Predicting cellular responses to complex perturbations. *Molecular Systems Biology* 18, e11517.
- Zheng, G.R. et al. (2023). Learning single-cell perturbation responses using neural optimal transport. *Cell Systems* 14, 978-994.

## 作者

**Junbiao Xue**  
GitHub: [@JunbiaoXue](https://github.com/JunbiaoXue)

## 许可

MIT License
