#!/usr/bin/env python3

import sys, gzip, os

POP = sys.argv[1]           # e.g. EUR or AFR or ALL for global AF
vcf_dir = sys.argv[2]       # folder containing downloaded VCFs
outfile = sys.argv[3]

vcf_files = sorted([x for x in os.listdir(vcf_dir) if x.endswith(".vcf.gz")])

af_key = "AF=" if POP == "ALL" else f"{POP}_AF="

with open(outfile, "w") as out:
    out.write("SNP\tCHR\tREF\tALT\tAF_ref\n")

    for vcf in vcf_files:
        print(f"Processing {vcf} ...")
        with gzip.open(os.path.join(vcf_dir, vcf), "rt") as f:
            for line in f:
                if line.startswith("#"): 
                    continue
                
                fields = line.strip().split("\t")
                chrom, pos, snp, ref, alt, qual, filt, info = fields[:8]

                # Skip multiallelic ALT
                if "," in alt:
                    continue

                # Find AF in INFO field
                af = None
                for item in info.split(";"):
                    if item.startswith(af_key):
                        af = item.replace(af_key, "")
                        break
                if af is None:
                    continue

                try:
                    af = float(af)
                except ValueError:
                    continue

                out.write(f"{snp}\t{chrom}\t{ref}\t{alt}\t{af}\n")
