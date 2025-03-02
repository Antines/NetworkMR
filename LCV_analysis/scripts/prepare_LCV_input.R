#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {
  stop("Usage: Rscript prepare_LCV_input.R <merged_ld_file>")
}
merged_ld_file <- args[1]
nutrient_base <- tools::file_path_sans_ext(basename(merged_ld_file))
cat("Reading merged file with LD scores:", merged_ld_file, "\n")
merged_ld <- fread(merged_ld_file)

# Sort by CHR and BP
setorder(merged_ld, CHR, BP)

# Extract LCV input vectors:
# ell: LD scores (from column L2)
# z_pd: PD z-scores (from column z_pd)
# z_nutr: Nutrient z-scores (from column Zscore)
ell <- merged_ld$L2
z_pd <- merged_ld$z_pd
z_nutr <- merged_ld$Zscore

output_file <- paste0(nutrient_base, "_LCV_input.RData")
save(ell, z_pd, z_nutr, file = output_file)
cat("LCV input vectors saved to", output_file, "\n")
