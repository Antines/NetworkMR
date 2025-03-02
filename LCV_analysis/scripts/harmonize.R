#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) {
  stop("Usage: Rscript harmonize_for_LCV.R <nutrient_file> <pd_file>")
}
nutrient_file <- args[1]
pd_file <- args[2]
nutrient_base <- tools::file_path_sans_ext(basename(nutrient_file))

cat("Reading nutrient file:", nutrient_file, "\n")
cat("Reading PD file:", pd_file, "\n")

nutrient <- fread(nutrient_file)
pd  <- fread(pd_file)

cat("Nutrient data head:\n")
print(head(nutrient))
cat("\nPD data head:\n")
print(head(pd))

# Merge the two datasets by rsID
merged <- merge(pd, nutrient, by = "SNP")

# Extract CHR and BP from PD's MarkerName column (format: "chrX:position")
merged[, c("CHR", "BP") := tstrsplit(MarkerName, ":", fixed = TRUE)]
merged[, CHR := as.numeric(gsub("chr", "", CHR, ignore.case = TRUE))]
merged[, BP  := as.numeric(BP)]


merged[, Allele1 := toupper(Allele1)]
merged[, Allele2 := toupper(Allele2)]
merged[, A1 := toupper(A1)]
merged[, A2 := toupper(A2)]
setorder(merged, CHR, BP)

# Compute PD z-score
merged[, z_pd := Effect / StdErr]

# Check allele alignment
merged[, allele_match := (Allele1 == A1) | (Allele1 == A2)]
if (any(!merged$allele_match)) {
  warning("Some SNP alleles appear misaligned between PD and nutrient data.")
} else {
  cat("Allele alignment check passed: All PD Allele1 values match one of nutrient alleles.\n")
}

output_file <- paste0(nutrient_base, "_merged_for_LCV.tsv")
fwrite(merged, output_file, sep = "\t")
cat("Merged and sorted summary statistics saved to:", output_file, "\n")
