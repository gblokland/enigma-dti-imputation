#!/usr/bin/env Rscript

# ---------- Packages ----------
suppressPackageStartupMessages({
  library(vcfR)
  library(ggplot2)
  library(dplyr)
  library(optparse)
  library(stringr)
})

# ---------- Command-line arguments ----------
option_list <- list(
  make_option(c("--vcf_folder"), type="character", help="Folder containing imputed VCFs (.vcf.gz)"),
  make_option(c("--ref"), type="character", help="Reference AF file (SNP, AF_ref, CHR)"),
  make_option(c("--out_plot"), type="character", default="AF_QQ_plot.png", help="Output QQ plot file"),
  make_option(c("--flip_threshold"), type="numeric", default=0.2, help="AF difference threshold to flag potential flips")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$vcf_folder) | is.null(opt$ref)) {
  print_help(opt_parser)
  stop("You must specify --vcf_folder and --ref", call.=FALSE)
}

# ---------- Step 1: Read reference panel ----------
ref_af <- read.table(opt$ref, header = TRUE, stringsAsFactors = FALSE)

# ---------- Step 2: Function to extract allele frequencies ----------
get_af_from_vcf <- function(vcf_path) {
  vcf <- read.vcfR(vcf_path, verbose = FALSE)
  gt <- extract.gt(vcf, element = "GT")
  
  af <- apply(gt, 1, function(g) {
    alleles <- unlist(strsplit(g, "[/|]"))
    alleles <- alleles[alleles != "."]
    if(length(alleles) == 0) return(NA)
    mean(as.numeric(alleles))
  })
  
  data.frame(
    SNP = vcf@fix[, "ID"],
    AF_imputed = af,
    CHR = vcf@fix[, "CHROM"],
    stringsAsFactors = FALSE
  )
}

# ---------- Step 3: Read all VCFs ----------
vcf_files <- list.files(opt$vcf_folder, pattern = "\\.vcf\\.gz$", full.names = TRUE)
if(length(vcf_files) == 0) stop("No VCF files found in folder")

imputed_af_list <- lapply(vcf_files, get_af_from_vcf)
imputed_af <- bind_rows(imputed_af_list)

# ---------- Step 4: Merge with reference panel ----------
merged_af <- imputed_af %>%
  inner_join(ref_af, by = c("SNP", "CHR")) %>%
  filter(!is.na(AF_imputed) & !is.na(AF_ref))

# ---------- Step 5: Flag potential strand flips ----------
merged_af <- merged_af %>%
  mutate(
    AF_diff = abs(AF_imputed - AF_ref),
    flip_flag = ifelse(
      AF_diff > opt$flip_threshold & abs(AF_imputed + AF_ref - 1) < 0.1,
      TRUE,
      FALSE
    )
  )

# ---------- Step 6: QQ plot ----------
p <- ggplot(merged_af, aes(x = AF_ref, y = AF_imputed, color = flip_flag)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "orange")) +
  labs(
    x = "Reference Panel Allele Frequency",
    y = "Imputed Data Allele Frequency",
    title = "Observed vs Expected Allele Frequencies",
    color = "Potential Strand Flip"
  ) +
  theme_minimal() +
  coord_equal()

ggsave(opt$out_plot, p, width = 6, height = 6, dpi = 300)
print(p)

# ---------- Step 7: Summary table per chromosome ----------
summary_table <- merged_af %>%
  group_by(CHR) %>%
  summarize(
    total_SNPs = n(),
    flagged_flips = sum(flip_flag, na.rm = TRUE),
    pct_flagged = 100 * flagged_flips / total_SNPs
  )

print(summary_table)
write.csv(summary_table, "strand_flip_summary_by_chr.csv", row.names = FALSE)

# Save all flagged SNPs
flagged_snps <- merged_af %>% filter(flip_flag)
write.csv(flagged_snps, "flagged_strand_flips.csv", row.names = FALSE)
