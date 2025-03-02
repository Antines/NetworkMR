#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {
  stop("Usage: Rscript merge_LD_scores.R <merged_file>")
}
merged_file <- args[1]
nutrient_base <- tools::file_path_sans_ext(basename(merged_file))

ld_dir <- path.expand("data/eur_w_ld_chr/")
output_file <- paste0(nutrient_base, "_merged_with_ld.tsv")

cat("Reading merged file:", merged_file, "\n")
merged <- fread(merged_file)
total_snps <- nrow(merged)
cat("Total SNPs in merged file:", total_snps, "\n")

ld_list <- list()
for(chr in 1:22) {
  ld_file <- file.path(ld_dir, paste0(chr, ".l2.ldscore.gz"))
  if(file.exists(ld_file)){
    ld_chr <- fread(ld_file)
    # Keep only SNP and L2 (LD score) columns
    ld_chr <- ld_chr[, .(SNP, L2)]
    ld_list[[chr]] <- ld_chr
    cat("Chromosome", chr, "processed:", nrow(ld_chr), "SNPs\n")
  } else {
    warning(paste("LD score file for chromosome", chr, "not found."))
  }
}
ld_all <- rbindlist(ld_list)

# Merge on SNP ID
merged_ld <- merge(merged, ld_all, by = "SNP", all.x = TRUE)
missing_count <- sum(is.na(merged_ld$L2))
percentage_missing <- round((missing_count / total_snps) * 100, 2)
cat("Number of SNPs with missing LD scores:", missing_count, 
    "out of", total_snps, "(", percentage_missing, "%)\n")

# Filter out SNPs missing LD scores
merged_ld <- merged_ld[!is.na(L2)]
fwrite(merged_ld, output_file, sep = "\t")
cat("Merged file with LD scores saved to", output_file, "\n")
