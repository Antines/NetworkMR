# Network MR Workflow

This repository contains several components for performing genetic analyses, including formatting GWAS results, conducting Mendelian Randomization (MR) analyses, and running polygenic risk score (PRS) workflows. A key component is an LCV (Latent Causal Variable) pipeline that integrates nutrient GWAS data with Parkinson’s disease (PD) data to assess potential partial genetic causality.

## Repository Structure

```
NetworkMR/
├── format_gwas.py
├── LCV_analysis
│   ├── lcv_pipeline.sh         # Bash script to run the complete LCV pipeline for multiple nutrients
│   ├── parse_results.py        # Python script to parse and summarize LCV outputs
│   └── scripts                 # R scripts for the LCV pipeline
│       ├── call_LCV.R          # Wrapper to call LCV functions (if needed)
│       ├── harmonize.R         # Script to merge and harmonize nutrient and PD summary statistics
│       ├── merge_LD_scores.R   # Script to merge precomputed LD scores with merged summary statistics
│       ├── MomentFunctions.R   # Functions for moment calculations (derived and modified from [LCV](https://github.com/lukejoconnor/LCV))
│       ├── prepare_LCV_input.R # Script to extract LCV input vectors from merged data with LD scores
│       └── RunLCV.R            # Main LCV function (derived and modified from [LCV](https://github.com/lukejoconnor/LCV))
├── MRscript.R
├── README.md
└── RunPRS.sh
```

## Components

### 1. Python GWAS Formatting Script
- **Script:** `format_gwas.py`
- **Description:** This script converts various GWAS result files into a standardized format suitable for downstream analyses.
- **Usage:**
  ```bash
  python format_gwas.py \
    --input <input_gwas_file.tsv> \
    --output <formatted_output_file.tsv> \
    --N <Sample Size>
  ```

### 2. R-based Network Mendelian Randomization (MR) Pipeline
- **Script:** `MRscript.R`
- **Description:** An R script that executes a network-style Mendelian Randomization pipeline. It performs data harmonization, Steiger filtering, MR-PRESSO, and MR summary analyses.
- **Usage:**  
  Run the script in R:
  ```r
  source("MRscript.R")
  ```

### 3. Bash-based PRS Workflow
- **Script:** `RunPRS.sh`
- **Description:** This shell script performs polygenic risk score (PRS) calculations for nutrient-related GWAS results and tests their association with a phenotype (e.g., Parkinson’s disease).
- **Usage:**  
  Execute the script in your shell:
  ```bash
  bash RunPRS.sh
  ```
  *Note:* Ensure that `plink` (v1.9+) is installed and available in your system path.

### 4. LCV Analysis Pipeline
The LCV pipeline integrates nutrient GWAS summary statistics with PD GWAS data to investigate partial genetic causality between traits. It involves several steps implemented through R scripts (located in `LCV_analysis/scripts/`):

#### Step 1: Harmonize Nutrient GWAS with PD GWAS
```bash
Rscript scripts/harmonize.R <nutrient_summary_stats> <PD_summary_stats>
```
- Merges and aligns summary statistics based on SNP rsIDs.
- Outputs a file named `<nutrient>_merged_for_LCV.tsv`.

#### Step 2: Merge with LD Scores
```bash
Rscript scripts/merge_LD_scores.R <nutrient_merged_for_LCV.tsv>
```
- Merges precomputed LD scores with the harmonized data.
- Outputs `<nutrient>_merged_for_LCV_merged_with_ld.tsv`.

#### Step 3: Prepare LCV Input Data
```bash
Rscript scripts/prepare_LCV_input.R <nutrient_merged_for_LCV_merged_with_ld.tsv>
```
- Extracts relevant variables (`ell`, `z_pd`, and `z_nutr`) for the LCV model.
- Saves an RData object, e.g., `<nutrient>_merged_for_LCV_merged_with_ld_LCV_input.RData`.

#### Step 4: Run LCV Analysis
```bash
Rscript scripts/call_LCV.R <LCV_input_file> <PD_sample_size> <Nutrient_sample_size>
```
- Runs the LCV model using the prepared data and specified sample sizes.
- Outputs results to a text file (e.g., `<nutrient>_PD_LCV.txt`).

#### Automated Pipeline Execution
To run the complete LCV analysis for all nutrient datasets, use the provided shell pipeline script:
```bash
bash LCV_analysis/lcv_pipeline.sh
```
This script loops over each nutrient file in the nutrient directory, performs harmonization, LD score merging, LCV input preparation, and finally runs the LCV analysis. Logs for each step are saved in a dedicated `logs/` directory.

---

## Requirements and Dependencies

### Python Dependencies (for `format_gwas.py`)
- Python 3.x
- `pandas`
- `argparse`

### R Dependencies (for MR and LCV Analysis)
- R (>= 4.0)
- R packages: `data.table`, `TwoSampleMR`, `ieugwasr`, `MRInstruments`, `MRPRESSO`, `ggplot2`, `dplyr`
- Install missing packages via CRAN or GitHub (using `remotes::install_github()` if needed)

### External Tools
- **PLINK:** Version 1.9+ (for PRS workflow)
- **LD Score Files:** Precomputed LD score files should be located in `data/eur_w_ld_chr/` (accessible at https://alkesgroup.broadinstitute.org/LDSCORE/).

---

## Troubleshooting and Notes

- **LCV Analysis Quality:**  
  Reliable LCV analysis requires well-powered GWAS datasets. Nutrients with low sample sizes or poor heritability estimates may yield unstable or non-interpretable LCV results.
  
- **File Paths:**  
  Ensure that file paths in the scripts (especially for LD score files and GWAS summary statistics) are correctly specified and accessible from your working environment.

- **Derived Scripts:**  
  The LCV scripts `MomentFunctions.R` and `RunLCV.R` in the `LCV_analysis/scripts/` directory are derived and modified from the [LCV GitHub repository](https://github.com/lukejoconnor/LCV) to integrate seamlessly with this pipeline.

- **Debugging:**  
  If any step fails or produces errors, check the log files in the `logs/` directory (when using the automated pipeline) or the console output. Common issues include file path mismatches, memory limitations, or low-quality GWAS data.

---

## Acknowledgments

- The original LCV methodology and code are derived from the [LCV GitHub repository](https://github.com/lukejoconnor/LCV).

---

