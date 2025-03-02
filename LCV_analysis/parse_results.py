#!/usr/bin/env python3

import glob
import re
import csv

def parse_merged_log(filename):

    total_snps = None
    missing_ld = None
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            m_total = re.search(r"Total SNPs in merged file:\s+(\d+)", line)
            if m_total:
                total_snps = int(m_total.group(1))
            m_missing = re.search(r"Number of SNPs with missing LD scores:\s+(\d+)\s+out of", line)
            if m_missing:
                missing_ld = int(m_missing.group(1))
    return total_snps, missing_ld

def parse_pd_lcv(filename):

    data = {}
    warning = ""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    for i, line in enumerate(lines):
        if "Warning message:" in line:
            for j in range(i+1, len(lines)):
                if lines[j].strip().startswith("In RunLCV("):
                    warning = lines[j].strip()
                    break
            break

    keys = ["zscore", "pval.gcpzero.2tailed", "gcp.pm", "gcp.pse", 
            "rho.est", "rho.err", "pval.fullycausal", "h2.zscore"]
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('$'):
            key = line[1:]
            if key in keys:
                j = i + 1
                while j < len(lines) and lines[j].strip() == "":
                    j += 1
                if j < len(lines):
                    num_line = lines[j].strip()
                    # Remove an initial "[1]" if present and split the numbers
                    num_line = re.sub(r'^\[1\]\s*', '', num_line)
                    nums = num_line.split()
                    if len(nums) == 1:
                        data[key] = float(nums[0])
                    else:
                        data[key] = tuple(float(x) for x in nums)
                i = j
        i += 1
    data["warning"] = warning
    return data

def main():
    results = {}
    ld_files = glob.glob("logs/*_ld_merge.log")
    for filepath in ld_files:
        base = filepath.split("/")[-1]
        nutrient = base.replace("_ld_merge.log", "")
        total_snps, missing_ld = parse_merged_log(filepath)
        if nutrient not in results:
            results[nutrient] = {}
        results[nutrient]["total_snps"] = total_snps
        results[nutrient]["missing_ld"] = missing_ld

    pd_files = glob.glob("logs/*_PD_LCV.txt")
    for filepath in pd_files:
        base = filepath.split("/")[-1]
        nutrient = base.replace("_PD_LCV.txt", "")
        pd_data = parse_pd_lcv(filepath)
        if nutrient not in results:
            results[nutrient] = {}
        results[nutrient].update(pd_data)

    headers = ["Nutrient", "Total SNPs", "Missing LD", "zscore", "pval.gcpzero.2tailed", 
               "gcp.pm", "gcp.pse", "rho.est", "rho.err", 
               "pval.fullycausal_1", "pval.fullycausal_2", 
               "h2.zscore_1", "h2.zscore_2", "Warning"]

    with open("summary.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        
        for nutrient in sorted(results.keys()):
            row = results[nutrient]
            pval_fc = row.get("pval.fullycausal", ("", ""))
            if isinstance(pval_fc, float):
                pval_fc = (pval_fc, "")
            h2_z = row.get("h2.zscore", ("", ""))
            if isinstance(h2_z, float):
                h2_z = (h2_z, "")
            
            row_values = [
                nutrient,
                row.get("total_snps", ""),
                row.get("missing_ld", ""),
                row.get("zscore", ""),
                row.get("pval.gcpzero.2tailed", ""),
                row.get("gcp.pm", ""),
                row.get("gcp.pse", ""),
                row.get("rho.est", ""),
                row.get("rho.err", ""),
                pval_fc[0],
                pval_fc[1],
                h2_z[0],
                h2_z[1],
                row.get("warning", "")
            ]
            writer.writerow(row_values)

if __name__ == "__main__":
    main()
