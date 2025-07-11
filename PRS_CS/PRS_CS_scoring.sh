#!/bin/bash
#SBATCH -J score5
#SBATCH --ntasks=15
#SBATCH --mem=70GB
#SBATCH -p zen3_0512
#SBATCH --qos=zen3_0512


ldREF='/gpfs/data/fs71407/epiwielsch/ukb_2/application/exploratory_analysis/genomicSEM/PRS_CS_from_factors/refPanels/ldblk_1kg_eur'
wdir='/gpfs/data/fs71407/epiwielsch/ukb_2/application/exploratory_analysis/genomicSEM/PRS_CS_from_factors/'


/home/fs71407/epiwielsch/PRScs/./PRScs.py --ref_dir=${ldREF} \
--bim_prefix=${wdir}plink_dataset/per_chrom_bed/gSEM_SCORE_data \
--sst_file=${wdir}factor_sumstats/PRS_CS_READY_sleep_factor2_multiGWAS_FINAL.txt \
--n_gwas=3172 \
--out_dir=${wdir}/SCORES_FINAL/PRS_CS_sleep_factor2
