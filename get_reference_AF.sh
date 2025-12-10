#!/bin/bash

# Usage:
#   bash get_reference_AF.sh 1000G EUR out/ref_AF.txt
#   bash get_reference_AF.sh HRC NONE out/ref_AF.txt

PANEL=$1              # "1000G" or "HRC"
POP=$2                # e.g. "EUR", "AFR", "ALL", or "NONE" for HRC
OUTFILE=$3           # Output ref_AF.txt

mkdir -p tmp_ref/
cd tmp_ref

############################################
### OPTION 1: 1000 Genomes Phase 3 AFs   ###
############################################
if [ "$PANEL" == "1000G" ]; then
    echo "Downloading 1000 Genomes Phase 3 VCFs..."
    for CHR in {1..22}; do
        wget -q https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        wget -q https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
    done
    wget -q -O ALL.chr23.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz
    wget -q -O ALL.chr23.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz.tbi

    echo "Extracting allele frequencies (1000G, POP=$POP)..."
    python3 ../extract_AF_1000G.py "$POP" . "$OUTFILE"
    exit 0
fi

############################################
### OPTION 2: HRC (mac5 sites table)     ###
############################################
if [ "$PANEL" == "HRC" ]; then
    echo "Downloading HRC sites table (GRCh37)"
    wget -q ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

    echo "Extracting allele frequencies using R..."
    Rscript ../extract_AF_HRC.R HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz "$OUTFILE"
    exit 0
fi

echo "ERROR: PANEL must be '1000G' or 'HRC'"
exit 1
