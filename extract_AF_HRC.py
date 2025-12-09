#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]

library(data.table)

cat("Reading HRC sites...\n")
dt <- fread(infile)

# HRC format (columns vary slightly but typically include):
# CHR POS ID REF ALT AC AN ...
#
# Compute AF = AC / AN

if(!all(c("ID","REF","ALT","AC","AN","CHR") %in% names(dt))) {
  stop("HRC format missing required columns (CHR, ID, REF, ALT, AC, AN)")
}

dt[, AF_ref := AC / AN]

out <- dt[, .(SNP = ID,
              CHR = CHR,
              REF = REF,
              ALT = ALT,
              AF_ref = AF_ref)]

fwrite(out, outfile, sep="\t")

cat("Done. Written:", outfile, "\n")
