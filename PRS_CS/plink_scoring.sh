#!/bin/bash
#SBATCH -J merge
#SBATCH -N 1
#SBATCH --mem=90GB
#SBATCH -p skylake_0096
#SBATCH --qos=skylake_0096

source activate plink2

dir='/gpfs/data/fs71407/epiwielsch/ukb_2/application/exploratory_analysis/genomicSEM/PRS_CS_from_factors/'
qc='/home/fs71407/epiwielsch/data_ukb/application/ukb48576/genotyped/final_lists/'
out='/gpfs/data/fs71407/epiwielsch/ukb_2/application/exploratory_analysis/genomicSEM/PRS_CS_from_factors/PLINK_SCORES/final_scores/'

## this assumes ID in column 1; Allele in column 2; and score (i.e. beta) in column 3

while read trait; do

plink \
        --bfile ${dir}plink_dataset/per_chrom_bed/gSEM_SCORE_data \
		--keep ${qc}ukb48576-basic-qc-all.fam \
		--score ${dir}/PLINK_SCORES/score_input/aligned_${trait}.txt 1 2 3 header \
		--out ${out}${trait}_plinkSCORE

done < /gpfs/data/fs71407/epiwielsch/ukb_2/application/exploratory_analysis/genomicSEM/PRS_CS_from_factors/PLINK_SCORES/score_TRAITS