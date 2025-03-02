#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3) {
  stop("Usage: Rscript run_LCV.R <LCV_input_file> <n_pd> <n_nutr>")
}
input_file <- args[1]
n_pd   <- as.numeric(args[2])
n_nutr <- as.numeric(args[3])

cat("Loading LCV input from", input_file, "\n")
load(input_file)

source("MomentFunctions.R")
source("RunLCV.R")

no_blocks <- 100  # number of jackknife blocks

lcv_result <- RunLCV(
  ell,    # LD score vector
  z_pd,   # PD z-scores (trait 1)
  z_nutr, # Nutrient z-scores (trait 2)
  no.blocks = no_blocks,
  n.1 = n_pd,
  n.2 = n_nutr
)

print(lcv_result)
