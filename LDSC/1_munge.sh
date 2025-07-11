#!/bin/bash

wdir='/gpfs/data/fs71407/epiwielsch/ukb_2/application/exploratory_analysis/LDSCORE_regression/'
ddir='/gpfs/data/fs71407/epiwielsch/ukb_2/application/exploratory_analysis//input_sumstats/result/'

while IFS=$'\t' read -r trait number; do
    echo "------ trait is ${trait}, number is ${number} ----"

cat <<EOF >run_munge_${trait}.sh
#!/bin/bash
#SBATCH -J ${trait}
#SBATCH -N 1
#SBATCH -p skylake_0096
#SBATCH --qos=skylake_0096
#SBATCH --tasks-per-node=20

source ~/anaconda3/etc/profile.d/conda.sh
source activate ldsc


##  -- adujst MAF

awk 'BEGIN {OFS="\t"; print "RSID","CHROM","POS","REF","ALT","BUILD","BETA","SE","P_value","MAF"} {if (NF == 10 && \$10 > 0.5) \$10 = 1 - \$10; print}' ${ddir}${trait}_SNPs_for_LDscore.txt > tmp_${trait}_for_munge.txt



/home/fs71407/epiwielsch/ldsc/./munge_sumstats.py \
--signed-sumstats beta,0 \
--out ${wdir}/munge_out/munge_${trait} \
--merge-alleles ${wdir}/w_hm3.snplist \
--N ${number} \
--a1 REF \
--a2 ALT \
--frq MAF \
--snp RSID \
--sumstats tmp_${trait}_for_munge.txt \
--p P_value

rm tmp_${trait}_for_munge.txt

EOF

chmod +x run_munge_${trait}.sh
sbatch run_munge_${trait}.sh


done < trait_names
