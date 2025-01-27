remotes::install_github("MRCIEU/TwoSampleMR")
remotes::install_github("MRCIEU/ieugwasr")  
devtools::install_github("MRCIEU/MRInstruments")
devtools::install_github("rondolab/MR-PRESSO")
library(MRInstruments)
library(TwoSampleMR)
require(ggplot2)	
library(dplyr)
library(MRPRESSO)
library(data.table)

#######################################################################################################################
# EXPOSURE | MEDIATORS Network MR RUN
#######################################################################################################################

exposure_outcome_pairs <- list(
  "calcium_GWAS.csv"       = "BMD.txt",
  "phosphorus_GWAS.csv"    = "BMD.txt",
  "vitamin_d_GWAS.csv"     = "BMD.txt",
  "vitamin_k_GWAS.csv"     = "BMD.txt",
  "iron_GWAS.csv"          = "ferritin.txt",
  "magnesium_GWAS.csv"     = "blood_pressure.txt",
  "potassium_GWAS.csv"     = "blood_pressure.txt",
  "sodium_GWAS.csv"        = "blood_pressure.txt",
  "fat_intake_GWAS.csv"    = "HDL.txt",
  "vitamin_b1_GWAS.csv"    = "homocesteine.txt",
  "vitamin_a_GWAS.csv"     = "CRP.txt",
  "vitamin_c_GWAS.csv"     = "CRP.txt",
  "vitamin_e_GWAS.csv"     = "CRP.txt"
)

for (exposure_file in names(exposure_outcome_pairs)) {
  
  tryCatch({
    # Get the outcome file from the list
    outcome_file <- exposure_outcome_pairs[[exposure_file]]
    
    message("=== Now processing ===")
    message("Exposure: ", exposure_file)
    message("Outcome:  ", outcome_file)
    
    # Create a folder
    exposure_name <- tools::file_path_sans_ext(exposure_file)
    if (!dir.exists(exposure_name)) {
      dir.create(exposure_name)
    }
    
    # Read the exposure data
    exp_data <- read_exposure_data(
      filename = file.path("data", "exposure", "exp", exposure_file), 
      sep = ",", 
      snp_col = "SNP", 
      beta_col = "beta.exposure",
      eaf_col = "eaf.exposure",
      se_col = "se.exposure",
      effect_allele_col = "effect_allele.exposure",
      other_allele_col = "other_allele.exposure",
      pval_col = "pval.exposure",
      samplesize_col = "samplesize.exposure",
      clump = FALSE
    )
    
    # If all key columns are NA, skip this exposure
    if (all(is.na(exp_data$effect_allele.exposure)) && 
        all(is.na(exp_data$other_allele.exposure)) || 
        all(is.na(exp_data$eaf.exposure))) {
      message("Skipping ", exposure_file, ": no valid exposure data.")
      next
    }
    
    # Read the outcome data
    out_data <- read_outcome_data(
      snps                 = exp_data$SNP,
      filename             = file.path("data", "mediator", "out", outcome_file),
      sep                  = "\t",
      snp_col              = "SNP",
      beta_col             = "beta",
      se_col               = "se",
      eaf_col              = "EAF",
      effect_allele_col    = "A1",
      other_allele_col     = "A2",
      pval_col             = "Pvalue",
      samplesize_col       = "N"
    )
    out_data$r.outcome<- get_r_from_pn(out_data$pval.outcome, out_data$samplesize.outcome)
    
    # If no overlapping SNPs, skip this pair
    if (nrow(out_data) == 0) {
      message("No overlapping SNPs in outcome data for ", exposure_file)
      next
    }
    
    # Harmonise the data
    harmonise_ok <- TRUE
    dat <- NULL
    tryCatch({
      dat <- harmonise_data(exposure_dat = exp_data, outcome_dat = out_data, action = 2)
    }, error = function(e) {
      message("harmonise_data() failed for ", exposure_file, " with error: ", e$message)
      harmonise_ok <<- FALSE
    })
    
    # Skip if harmonization fails or returns empty data
    if (!harmonise_ok || is.null(dat) || nrow(dat) == 0) {
      message("Skipping ", exposure_file, ": harmonise_data() returned empty or errored out.")
      next
    }
    
    #dat$units.outcome <- "log odds"
    
    # Subset to valid data
    dat1 <- subset(dat, !is.na(eaf.exposure))
    dat1$r.exposure <- get_r_from_pn(dat1$pval.exposure, dat1$samplesize.exposure)
    
    # Steiger filtering
    steiger_res <- steiger_filtering(dat1)
    sig <- subset(steiger_res, steiger_dir == TRUE)
    if (is.null(sig) || nrow(sig) == 0) {
      message("No valid SNPs remain after Steiger filtering for ", exposure_file)
      next
    }
    
    # Run MR-PRESSO
    tryCatch({
      presso <- mr_presso(
        BetaOutcome     = "beta.outcome",
        BetaExposure    = "beta.exposure",
        SdOutcome       = "se.outcome",
        SdExposure      = "se.exposure",
        OUTLIERtest     = TRUE,
        DISTORTIONtest  = TRUE,
        data            = sig,
        NbDistribution  = 1000,
        SignifThreshold = 0.05
      )
      capture.output(print(presso), file = file.path(exposure_name, "presso.txt"))
    }, error = function(e) {
      message("MR-PRESSO failed for ", exposure_file, " with error: ", e$message)
      write("Not enough IV occurrence", file = file.path(exposure_name, "Fail_to_presso.txt"))
    })
    
    # Run MR
    tryCatch({
      res <- mr(sig)
      capture.output(print(res), file = file.path(exposure_name, "mr_main_results.txt"))
    }, error = function(e) {
      message("MR analysis failed for ", exposure_file, " with error: ", e$message)
    })
    
    # Calculate R2 and F-statistics
    R2 <- mean(sig$r.exposure, na.rm = TRUE)
    n  <- mean(sig$samplesize.exposure, na.rm = TRUE)
    k  <- nrow(subset(sig, ambiguous == FALSE))
    
    # Avoid division by zero in F-stat calculation
    if (R2 < 1) {
      F_stat <- (R2 * (n - 1 - k)) / ((1 - R2) * k)
    } else {
      F_stat <- NA
    }
    
    capture.output(print(R2), file = file.path(exposure_name, "r2.txt"))
    capture.output(print(F_stat), file = file.path(exposure_name, "f.txt"))
    
    # Generate MR report
    tryCatch({
      mr_report(sig, study = exposure_name, output_path = exposure_name)
    }, error = function(e) {
      message("MR report failed for ", exposure_file, " with error: ", e$message)
    })
    
    # Single SNP results and forest plot
    tryCatch({
      res_single <- mr_singlesnp(sig)
      if (nrow(res_single) < 2) {
        message("Cannot create forest plot with <2 SNPs for ", exposure_file)
      } else {
        p5 <- mr_forest_plot(res_single)
        ggsave(
          filename = "forest_plot.jpg",
          plot     = p5[[1]],
          path     = exposure_name,
          width    = 7,
          height   = 12
        )
      }
    }, error = function(e) {
      message("Single SNP analysis or forest plot failed for ", exposure_file, " with error: ", e$message)
    })
    
    message("Completed MR for exposure=", exposure_file, " outcome=", outcome_file)
    
  }, error = function(e) {
    message("Error processing exposure=", exposure_file, 
            " outcome=", outcome_file, ": ", e$message)
  })
}

#######################################################################################################################
# PARSE RESULTS
#######################################################################################################################

setwd("$resultsDir")

# List all immediate subdirectories
folders <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
results_list <- list()

for (folder in folders) {
  r2_file          <- file.path(folder, "r2.txt")
  mr_report_file   <- file.path(folder, "mr_report.md")
  presso_file      <- file.path(folder, "presso.txt")
  fail_presso_file <- file.path(folder, "Fail_to_presso.txt")
  
  if (!(file.exists(r2_file) && file.exists(mr_report_file))) {
    message("Skipping folder '", folder, 
            "' because 'r2.txt' or 'mr_report.md' not found.")
    next
  }
  
  # Parse R^2
  r2_content <- readLines(r2_file, warn = FALSE)
  r2_value <- gsub("^\\[\\d+\\]\\s+", "", r2_content[1])
  r2_value <- as.numeric(r2_value)
  
  # Parse IVW from mr_report.md
  report_lines <- readLines(mr_report_file, warn = FALSE)
  
  ivw_line_idx <- grep("inverse variance weighted", report_lines, ignore.case = TRUE)
  n_snps   <- NA_integer_
  ivw_se   <- NA_real_
  ivw_pval <- NA_real_
  ivw_beta <- NA_real_  # Initialize beta as NA_real_
  
  if (length(ivw_line_idx) > 0) {
    parsed_line <- FALSE
    for (idx in ivw_line_idx) {
      line_clean <- report_lines[idx] |> 
        trimws() |>
        sub("^\\|+", "", x = _) |>
        sub("\\|+$", "", x = _) |>
        trimws()
      
      splitted <- strsplit(line_clean, "\\|")[[1]]
      splitted <- trimws(splitted)
      
      if (length(splitted) >= 5) {
        n_snps   <- as.integer(splitted[2])
        ivw_beta <- as.numeric(splitted[3])  # Extract beta from the 3rd column
        ivw_se   <- as.numeric(splitted[4])
        ivw_pval <- as.numeric(splitted[5])
        parsed_line <- TRUE
        break
      }
    }
    if (!parsed_line) {
      message("IVW method line found but could not be parsed in folder '", folder, "'.")
    }
  } else {
    message("IVW method not found in '", mr_report_file, "' for folder '", folder, "'.")
  }
  
  # Parse MR-PRESSO p-value
  mr_presso_val <- NA_character_  # Force character from the start
  
  if (file.exists(presso_file)) {
    # parse numeric p-value
    presso_lines <- readLines(presso_file, warn = FALSE)
    pval_header_idx <- grep("Pvalue$", presso_lines)
    if (length(pval_header_idx) > 0) {
      target_idx <- pval_header_idx[length(pval_header_idx)] + 1
      if (target_idx <= length(presso_lines)) {
        line_val <- presso_lines[target_idx]            # e.g. "[1] 0.953"
        line_val <- gsub("^\\[\\d+\\]\\s+", "", line_val) # remove "[1] "
        numeric_val <- as.numeric(line_val)
        
        # Convert numeric to character
        if (!is.na(numeric_val)) {
          mr_presso_val <- as.character(numeric_val)  # e.g. "0.953"
        }
      }
    }
  } else if (file.exists(fail_presso_file)) {
    # read entire error message as character
    fail_lines     <- readLines(fail_presso_file, warn = FALSE)
    mr_presso_val  <- paste(fail_lines, collapse = " ")
  }
  
  # Build data frame row
  one_res <- data.frame(
    "Exposure Reported trait" = basename(folder),
    "N_SNPs"                  = n_snps,
    "se"                      = ivw_se,
    "beta"                    = ivw_beta,
    "IVW_P-value"             = ivw_pval,
    "R2"                      = r2_value,
    "MR-PRESSO P-value"       = mr_presso_val,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  results_list[[folder]] <- one_res
}

# Combine all results into a single data frame
final_results <- bind_rows(results_list)
final_results <- final_results %>%
  arrange(`Exposure Reported trait`)

write.csv(final_results, "MR_summary.csv", row.names = FALSE)
print(final_results)
