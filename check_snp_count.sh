#!/bin/bash

check_snp_count() {
    local snp_count="$1"
    local min_ok=8000000      # absolute minimum acceptable
    local expected=10000000   # expected ideal

    if [[ -z "$snp_count" ]]; then
        echo "ERROR: No SNP count provided to check_snp_count()"
        return 2
    fi

    echo "SNPs after filtering: $snp_count"

    # Hard failure: below minimum acceptable threshold
    if (( snp_count < min_ok )); then
        echo "❌ ERROR: SNP count ($snp_count) is BELOW the acceptable minimum ($min_ok)."
        echo "   This strongly suggests issues with strand alignment, MAF filtering, INFO filtering, or input QC."
        return 1
    fi

    # Warning: below expected range but above acceptable minimum
    if (( snp_count < expected )); then
        echo "⚠️ WARNING: SNP count is lower than expected (~10M), but above the minimum."
        echo "   Investigate potential SNP loss (e.g., over-aggressive filtering, mismatched build/strand)."
        return 0
    fi

    # All good
    echo "✅ SNP count is within the expected range."
    return 0
}
