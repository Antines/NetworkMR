#!/usr/bin/env bash

# =============================================================================
# A workflow to compute nutrient-based PRS and test association with PD,
# using final summary stats in 'nutrients' and writing outputs to 'PD_sup_results'.
# =============================================================================

# 1. Paths & Parameters
# -----------------------------------------------------------------------------
GENO_PREFIX="PD_genotype/PD" 
NUTRIENTS_DIR="nutrients_data"
RESULTS_DIR="PD_sup_results"

# Clumping thresholds
CLUMP_P1=5e-8
CLUMP_P2=1e-2
CLUMP_R2=0.1
CLUMP_KB=250

THREADS=4

# Create the results directory if missing
mkdir -p "${RESULTS_DIR}"

# 2. Basic QC on Genotype Data
# -----------------------------------------------------------------------------
echo "Step 2: Genotype Quality Control (QC)"
plink \
  --bfile "${GENO_PREFIX}" \
  --maf 0.01 \
  --geno 0.05 \
  --hwe 1e-6 midp \
  --make-bed \
  --out "${RESULTS_DIR}/PD.QC"

# 3. Loop Over Each Nutrient .fin.tsv File & Build a PRS
# -----------------------------------------------------------------------------
for FILE in "${NUTRIENTS_DIR}"/*.tsv; do

  BASENAME=$(basename "${FILE}")
  TRAIT="${BASENAME%.*}"   # e.g. iron.fin => iron

  echo "Processing trait: ${TRAIT}"

  # 3a. Reformat Summary Stats for PLINK
  #    Adjust these columns if Beta/P are in different columns in your .fin.tsv
  FORMATTED="${RESULTS_DIR}/${TRAIT}_sumstats_for_plink.txt"
  echo "  [3a] Formatting summary stats..."
  
  awk 'BEGIN{OFS="\t"; print "SNP","ALT","BETA","P"} 
       FNR>1 {
         snp=$1; alt=$2; beta=$7; p=$5;
         print snp,alt,beta,p
       }' "${FILE}" > "${FORMATTED}"

  # # 3b. Clump to Select Independent SNPs
  echo "  [3b] Clumping..."
  plink \
    --bfile "${RESULTS_DIR}/PD.QC" \
    --clump "${FORMATTED}" \
    --clump-p1 ${CLUMP_P1} \
    --clump-p2 ${CLUMP_P2} \
    --clump-r2 ${CLUMP_R2} \
    --clump-kb ${CLUMP_KB} \
    --out "${RESULTS_DIR}/${TRAIT}_clumped"

  CLUMP_FILE="${RESULTS_DIR}/${TRAIT}_clumped.clumped"
  if [ ! -f "${CLUMP_FILE}" ]; then
    echo "  [WARNING] No .clumped file for ${TRAIT}; skipping PRS scoring."
    continue
  fi

  # # Extract the clumped SNPs
  CLUMP_SNPS="${RESULTS_DIR}/${TRAIT}_clumped.snplist"
  awk 'NR>1 {print $3}' "${CLUMP_FILE}" > "${CLUMP_SNPS}"

  # # 3c. Create Scoring File (SNP ALT BETA)
  SCOREFILE="${RESULTS_DIR}/${TRAIT}_scorefile.txt"
  echo "  [3c] Creating score file..."
  awk 'NR==FNR {clumped[$1]=1; next}
       (FNR==1){print $1,$2,$3} 
       (FNR>1 && ($1 in clumped)){print $1,$2,$3}' \
       "${CLUMP_SNPS}" "${FORMATTED}" \
       > "${SCOREFILE}"

  # 3d. Compute PRS
  echo "  Computing PRS..."
  plink \
    --bfile "${RESULTS_DIR}/PD.QC" \
    --score "${SCOREFILE}" 1 2 3 header \
    --out "${RESULTS_DIR}/${TRAIT}_PRS"

done

# 4. Merge All Nutrient PRS & Associate with PD in R
# -----------------------------------------------------------------------------

echo "Step 4: Merge PRSs and run logistic regression on PD"
Rscript -e '
library(dplyr)

# 4a. Collect all _PRS.profile files in PD_sup_results/
prs_files <- list.files("PD_sup_results", pattern="_PRS.profile$", full.names=TRUE)
if(length(prs_files) == 0) {
  stop("No PRS.profile files found in PD_sup_results. Did PRS scoring fail?")
}

prs_list <- list()
for (f in prs_files) {
  # e.g., PD_sup_results/iron.fin_PRS.profile => iron.fin
  trait <- sub("_PRS.profile","", basename(f))
  df <- read.table(f, header=TRUE)
  df <- df %>%
    select(FID, IID, SCORESUM) %>%
    rename(!!paste0(trait,"_PRS") := SCORESUM)
  prs_list[[trait]] <- df
}

# Merge all PRS data frames
merged_prs <- Reduce(function(x,y) merge(x,y, by=c("FID","IID"), all=TRUE), prs_list)
write.table(merged_prs, file="PD_sup_results/merged_nutrient_prs.txt", 
            sep="\t", quote=FALSE, row.names=FALSE)

# 4b. Read PD phenotype from PD_sup_results/PD.QC.fam
phen <- read.table("PD_sup_results/PD.QC.fam", header=FALSE)
colnames(phen) <- c("FID","IID","PID","MID","Sex","Pheno")

# Convert PLINK-coded Pheno (1=control, 2=case) to PD_status (0/1)
phen$PD_status <- ifelse(phen$Pheno == 2, 1, 0)

# Merge phen with PRS
data_merged <- merge(phen, merged_prs, by=c("FID","IID"))

# 4c. Simple logistic regression of PD_status ~ multiple PRSs (+ Sex)
prs_vars <- grep("_PRS$", names(data_merged), value=TRUE)
if(length(prs_vars) == 0) {
  stop("No PRS columns found. Check your PRS computation step.")
}

formula_str <- paste("PD_status ~", paste(prs_vars, collapse=" + "), "+ Sex")
model <- glm(formula_str, data=data_merged, family=binomial())

cat("\n=== Logistic Regression Results ===\n")
print(summary(model))
'
