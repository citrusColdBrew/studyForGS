# CLAUDE.md

本文件为Claude Code (claude.ai/code) 在此代码库中工作提供指导。

## 项目概述

这是一个专注于全基因组关联分析（GWAS）和基因组选择（GS）的统计基因组学研究代码库。该仓库包含讲座材料、实现代码和统计遗传学方法的分析工具。

## 仓库结构

- **讲座目录** (`Lecture01_phenotype/`, `Lecture02_GWASbyR/` 等)：涵盖不同统计遗传学主题的独立课程模块，每个都包含R脚本和支持材料
- **`function/`**：核心可重用函数和算法，包括：
  - `gapit_functions.R`：GAPIT（基因组关联和预测集成工具）框架，包含多种GWAS方法（GLM、MLM、BLUP等）
  - `G2P.R`：基因型到表型模拟函数
  - `NN-GWAS.R` 和 `NN-GS.R`：用于GWAS和基因组选择的神经网络实现
  - `GWASbyCor.R`：基于相关性的GWAS方法
  - `bernoulli.py`：用于统计建模的Python工具
- **`data/`**：标准格式的基因组数据集（mdp_numeric.txt为基因型，mdp_traits.txt为表型，mdp_SNP_information.txt为标记信息）
- **`studyForGS/`**：重复的学习材料目录

## 关键技术和依赖项

### R依赖项
代码库主要使用R语言和以下关键包：
- **统计遗传学**：genetics, EMMREML, lme4, multtest, snpStats (Bioconductor)
- **机器学习**：keras, tensorflow（用于神经网络实现）
- **数据处理**：data.table, dplyr, caret
- **可视化**：ggplot2, gplots, scatterplot3d
- **工具包**：bigmemory, ape, compiler, grid

### Python依赖项
有限的Python使用：
- 用于统计建模的标准科学计算库

## 开发工作流程

### 运行分析脚本
- 导航到特定讲座目录运行单个课程
- 引用主函数库：`source("../function/gapit_functions.R")`
- 使用相对路径从`data/`目录加载数据文件
- 大多数脚本是自包含的，内置了包安装检查

### 数据格式要求
- **基因型数据**：数值格式（0,1,2编码）的制表符分隔文件
- **表型数据**：包含个体ID和性状值的制表符分隔文件
- **标记信息**：SNP位置和染色体信息

### 关键函数用法
- **GAPIT**：支持11种以上GWAS方法的主要框架（GLM、MLM、BLUP、EMMA等）
- **G2P**：根据指定的遗传力和QTN效应从基因型模拟表型
- **NN-GWAS/NN-GS**：带交叉验证的神经网络方法

## 已实现的分析方法

代码库涵盖了全面的统计遗传学方法：
1. **GWAS方法**：GLM（Q方法）、MLM（Q+K）、gBLUP、PCA、EMMA、CMLM、P3D、FaST-LMM、SUPER
2. **模拟**：可配置遗传力和QTN分布的表型模拟
3. **机器学习**：用于关联分析和基因组预测的神经网络
4. **统计建模**：贝叶斯方法、混合线性模型、BLUP变体

## 使用代码库

- 脚本包含大量中英文注释，解释统计遗传学概念
- 每个讲座建立在之前的概念基础上 - 按数字顺序学习以获得最佳效果
- 函数文件有详细文档，包含参数说明和使用示例
- GAPIT框架版本为4.0（截至2025.03.03）