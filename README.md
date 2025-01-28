# Network MR Workflow

This repository contains three main components for genetic analysis:

1. **Python GWAS Formatting Script**  
   A Python script (`format_gwas.py`) to convert various GWAS result files into a standardized format.

2. **R-based Network Mendelian Randomization Script**  
   An R script that runs a network-style Mendelian Randomization pipeline, including data harmonization, Steiger filtering, MR-PRESSO, and MR summary analyses.

3. **Bash-based PRS Workflow**  
   A shell script that performs polygenic risk score (PRS) calculations for nutrient-related GWAS results and tests their association with a given phenotype (e.g., Parkinson’s disease).

Below is a comprehensive overview, including dependencies, usage instructions, and expected outputs.

---

## Requirements and Dependencies

### Python (for `format_gwas.py`)
- Python 3.x
- `pandas` (>= 1.0)
- Command-line usage of `argparse`

### R and R Packages (for Network MR)
- R (>= 4.0)
- `TwoSampleMR`
- `ieugwasr`
- `MRInstruments`
- `MRPRESSO`
- `ggplot2`
- `dplyr`
- `data.table`
- `devtools` or `remotes` (to install packages from GitHub)

### Shell Environment (for PRS workflow)
- Unix-like environment (Linux or macOS)
- `plink` (version 1.9+ recommended)
- R (for the final logistic regression step)

Ensure that `plink` is available on your system path if you plan to run the PRS workflow.

---

### Usage

#### 1. Summary stats foramting

**File:** `format_gwas.py`  
**Description:** Converts an input GWAS file (tab-delimited) into a standardized format with columns:



```bash
python scripts/format_gwas.py \
  --input <input_gwas_file.tsv> \
  --output <formatted_output_file.tsv> \
  --N <Sample Size
```
#### PRS Analysis

In this step, you will run PLINK to:

1. Filter out variants with minor allele frequency (MAF) below 0.01.
2. Exclude samples/variants with more than 5% missingness (`--geno 0.05`).
3. Exclude variants failing the Hardy-Weinberg equilibrium threshold of `1e-6` (`--hwe 1e-6 midp`).
4. Generate a clean PLINK dataset in binary format (`.bed/.bim/.fam`).

**Command:**
```bash
plink --bfile PD_genotype/PD \
      --maf 0.01 \
      --geno 0.05 \
      --hwe 1e-6 midp \
      --make-bed \
      --out PD_sup_results/PD.QC
```

#### MR analisys
```bash
source("network_MR_script.R")
Rscript scripts/network_MR_script.R
```
