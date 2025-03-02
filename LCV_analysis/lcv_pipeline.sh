#!/bin/bash

NUTRIENT_DIR="nutrients_data"
PD_FILE="pd_data/IPDGC_sum_stats_no_UKB.txt"
LOG_DIR="logs"
mkdir -p ${LOG_DIR}

N_PD=27693

# Loop over nutrient files
for nutrient in ${NUTRIENT_DIR}/*.tsv; do
    nutrient_base=$(basename "${nutrient}" .tsv)
    echo "Processing nutrient: ${nutrient_base}"
    
    N_NUTR=$(awk -F"\t" 'NR==2 {print $6; exit}' "${nutrient}")
    
    # Harmonize PD and nutrient data
    echo "Merging PD and ${nutrient_base} data..."
    Rscript scripts/harmonize.R "${nutrient}" "${PD_FILE}" > "${LOG_DIR}/${nutrient_base}_merge.log" 2>&1
    if [ $? -ne 0 ]; then
        echo "Error during merging for ${nutrient_base}. Check ${LOG_DIR}/${nutrient_base}_merge.log"
        continue
    fi

    MERGED_FILE="${nutrient_base}_merged_for_LCV.tsv"
    
    # Merge with LD scores
    echo "Merging with LD scores for ${nutrient_base}..."
    Rscript scripts/merge_LD_scores.R "${MERGED_FILE}" > "${LOG_DIR}/${nutrient_base}_ld_merge.log" 2>&1
    if [ $? -ne 0 ]; then
        echo "Error during LD score merging for ${nutrient_base}. Check ${LOG_DIR}/${nutrient_base}_ld_merge.log"
        continue
    fi

    MERGED_LD_FILE="${nutrient_base}_merged_for_LCV_merged_with_ld.tsv"
    
    # Prepare LCV input vectors
    echo "Preparing LCV input for ${nutrient_base}..."
    Rscript scripts/prepare_LCV_input.R "${MERGED_LD_FILE}"
    if [ $? -ne 0 ]; then
        echo "Error during LCV input preparation for ${nutrient_base}. Check ${LOG_DIR}/${nutrient_base}_prepare_LCV.log"
        continue
    fi

    LCV_INPUT_FILE="${nutrient_base}_merged_for_LCV_merged_with_ld_LCV_input.RData"
    
    # Run LCV analysis
    echo "Running LCV analysis for ${nutrient_base}..."
    Rscript scripts/call_LCV.R "${LCV_INPUT_FILE}" ${N_PD} ${N_NUTR} > "${LOG_DIR}/${nutrient_base}_PD_LCV.txt" 2>&1
    if [ $? -ne 0 ]; then
        echo "Error during LCV analysis for ${nutrient_base}. Check ${LOG_DIR}/${nutrient_base}_PD_LCV.txt"
        continue
    fi
    
    echo "Completed processing for ${nutrient_base}. LCV results in ${LOG_DIR}/${nutrient_base}_PD_LCV.txt"
done

echo "Pipeline completed for all nutrients."