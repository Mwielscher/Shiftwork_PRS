#!/bin/bash

while read -r fac; do
echo "------ this is ${fac} ----"
cat <<EOF >run_GCTA_${fac}.sh
#!/bin/bash

###  ------------------ 
wdir='/gpfs/data/fs71707/mwielsch1/ukb_2/application/exploratory_analysis/genomicSEM/multiVAR_GWAS/sum_result/'
imp_dir='/gpfs/data/fs71707/mwielsch1/ukb_2/application/ukb48576/imputed/plink_for_GCTA/'
chrlist='/gpfs/data/fs71707/mwielsch1/ukb_2/application/JAN_traits/CHR_LIST'

chr=\$(awk -v var="\$SLURM_ARRAY_TASK_ID" 'FNR==var {print \$1}' \${chrlist} )
## ----------------------------------------------------------------------
##---    set threshhold

if [ "${fac}" = "sleep_factor2" ] || [ "${fac}" = "sleep_factor5" ] || [ "${fac}" = "sleep_factor6" ]; then
    thresh=1e-5
else
    thresh=5e-8
fi


##   ---------------------    prep sumStats  --- sumStats are output form script 3
awk '{print \$1, \$5, \$6, \$4, \$7, \$8, \$10, \$17}' \${wdir}${fac}_multiGWAS_sorted_plusN_MAF_fixed.txt > ${fac}_chr\${chr}GCTA_input.ma


/home/fs71707/mwielsch1/gcta-1.94.1/./gcta-1.94.1 --bfile \${imp_dir}EUR/GCTA_EUR_chr\${chr} \
--chr \${chr} \
--cojo-file ${fac}_chr\${chr}GCTA_input.ma --cojo-slct \
--cojo-wind 10000 --cojo-collinear 0.9 \
--cojo-p \${thresh} --out \${wdir}${fac}_chr\${chr}_sentinel_snps


## --------- GCTA clean up
rm ${fac}_chr\${chr}GCTA_input.ma
#rm \${wdir}${fac}_chr\${chr}_sentinel_snps.badsnps 
rm \${wdir}${fac}_chr\${chr}_sentinel_snps.cma.cojo 
#rm \${wdir}${fac}_chr\${chr}_sentinel_snps.badsnps 
#rm \${wdir}${fac}_chr\${chr}_sentinel_snps.freq.badsnps 
rm \${wdir}${fac}_chr\${chr}_sentinel_snps.ldr.cojo
rm \${wdir}${fac}_chr\${chr}_sentinel_snps.log

EOF

chmod +x run_GCTA_${fac}.sh

cat <<EOFSUBMIT >submit_arrayGCTA_${fac}.sh
#!/bin/bash
#SBATCH -J small-chr
#SBATCH --ntasks=15
#SBATCH --mem=50GB
#SBATCH -p skylake_0096
#SBATCH --qos=skylake_0096
#SBATCH --array=13-22

./run_GCTA_${fac}.sh \$SLURM_ARRAY_TASK_ID

EOFSUBMIT

echo "------ conditional analysis for ${fac} submitted  ----"

chmod +x submit_arrayGCTA_${fac}.sh
sbatch submit_arrayGCTA_${fac}.sh

done<FAC_LIST

