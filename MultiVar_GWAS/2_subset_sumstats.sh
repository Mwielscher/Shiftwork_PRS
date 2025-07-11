#!/bin/bash


hapref='/gpfs/data/fs71407/epiwielsch/ukb_2/application/exploratory_analysis/genomicSEM/multiVAR_GWAS/lookup_ref'
indir='/gpfs/data/fs71407/epiwielsch/ukb_2/application/exploratory_analysis/input_sumstats/result/'
outdir='/gpfs/data/fs71407/epiwielsch/ukb_2/application/exploratory_analysis/genomicSEM/multiVAR_GWAS/input_sumstats/'


while read -r trait; do
echo "subsetting $trait now"
# Filter based on HapMap list and prepare for LD score
awk 'NR==FNR {rsids[$1]; next} $1 in rsids' "$hapref" "${indir}${trait}_allSNPs.txt" > "${outdir}${trait}_SNPs_for_genomicSEM.txt"


done<traits
