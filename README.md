# Nutrient-Based Mendelian Randomization and Polygenic Risk Score Analysis Workflow

## Table of Contents
- [Introduction](#introduction)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Directory Structure](#directory-structure)
- [Usage](#usage)
  - [1. Format GWAS Results](#1-format-gwas-results)
  - [2. Mendelian Randomization (MR) Analysis](#2-mendelian-randomization-mr-analysis)
  - [3. Polygenic Risk Score (PRS) Computation and Association with Parkinson's Disease (PD)](#3-polygenic-risk-score-prs-computation-and-association-with-parkinsons-disease-pd)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)
- [License](#license)
- [Contact](#contact)

## Introduction

This repository contains a comprehensive workflow for conducting Nutrient-Based Mendelian Randomization (MR) analyses and Polygenic Risk Score (PRS) computations to investigate the associations between various nutrient exposures and health outcomes, such as Bone Mineral Density (BMD), Ferritin levels, Blood Pressure, HDL cholesterol, Homocysteine levels, C-reactive Protein (CRP), and Parkinson's Disease (PD).

The workflow is divided into three main components:
1. **GWAS Result Formatting**: A Python script to standardize GWAS summary statistics.
2. **Mendelian Randomization Analysis**: An R script to perform MR analyses across multiple exposure-outcome pairs.
3. **Polygenic Risk Score Computation and Association Analysis**: A Bash script to compute PRS and assess their association with PD using logistic regression.

## Prerequisites

Before running the scripts, ensure that the following software and packages are installed on your system:

### General
- **Operating System**: Linux-based (tested on Ubuntu)
- **Python**: Version 3.6 or higher
- **R**: Version 4.0 or higher
- **PLINK**: Version 1.9 or higher
- **Git**: For cloning repositories

### Python Dependencies
- pandas
- argparse

Install Python dependencies using `pip`:
```bash
pip install pandas argparse
