#!/bin/bash


wdir='/gpfs/data/fs71407/epiwielsch/ukb_2/application/exploratory_analysis/genomicSEM/multiVAR_GWAS/'
chrlist='/gpfs/data/fs71407/epiwielsch/ukb_2/application/JAN_traits/CHR_LIST'
chr=$(awk -v var="$SLURM_ARRAY_TASK_ID" 'FNR==var {print $1}' ${chrlist} )

## -------------------------------
## prep sumstats

while read -r trait; do
awk -v chrom="$chr" 'BEGIN {print "rsID chrom pos A1 A2 beta SE P INFO"}
{
    if ($2 == chrom && $7 ~ /^[+-]?[0-9]+([.][0-9]+)?$/ ) { print $1, $2, $3, $4, $5, $7, $8, $9, "0.9" }
}' ${wdir}input_sumstats/${trait}_SNPs_for_genomicSEM.txt > ${wdir}/tmp_sumstats/${trait}_sumstat_chr${chr}.txt

done <${wdir}traits


cat <<EOF >gwas_${chr}.R
library(lavaan,lib.loc = "/gpfs/data/fs71407/epiwielsch/ukb_2/application/exploratory_analysis/genomicSEM/multiVAR_GWAS/lavaan")
library(GenomicSEM,lib.loc = "/gpfs/data/fs71407/epiwielsch/ukb_2/application/initial_4_traits/genomic_SEM")
wdir=c("${wdir}")
dat_dir=c("${wdir}tmp_sumstats/")
trait_f=read.table("${wdir}traits")
ref=paste0(wdir,"reference.1000G.maf.0.005.txt")
file1=paste0(dat_dir,trait_f\$V1,"_sumstat_chr${chr}.txt")
files <- c(file1)
## update names to fit LDSC object !! also order must be the same !!
trait_f\$V1[4:10]=paste0(trait_f\$V1[4:10],"_SIN")
trait_f\$V1[21]=paste0(trait_f\$V1[21],"_SIN")


## set parameters:
N=c(446118,272683.71,487268.57,84757,85205,85670,85449,85449
    ,84810,84810,22700,452071,452071,452071,446118,697828,71500,19733,22700,22700,84810,344615.4)
N=as.numeric(N)

se_logit=c(T,T,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T)

OLS=c(F,F,F,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,F)


test_sumstats <- sumstats(files=files,ref=ref,trait.names=c(trait_f\$V1),
    se.logit=se_logit, OLS = OLS ,N=N,betas=NULL, maf.filter=0.01,info.filter=0.6,
    keep.indel=FALSE,parallel=FALSE,cores=1)
save(test_sumstats,file="${wdir}test_sumstats_${chr}.RData")

##  --------------------------------------------------
##### run the GWAS !!!

load("${wdir}LDSCoutput_for_GWAS_20240523.RData")
#specify the model

model <- 'F1 =~ NA*SR_chronotype + SR_morning_person + ACC_L5_TIME_RAW_SIN + ACC_M10_TIME_RAW_SIN + ACC_SLEEP_MIDP_RAW_SIN
             F2 =~ NA*relativeAmplitude + ACC_SLEEP_DUR_SD_RAW_SIN
             F3 =~ NA*SR_SleepDuration + ACC_SLEEP_DUR_RAW_SIN + SR_LongSleep + SR_ShortSleep
             F4 =~ NA*SR_DaytimeNapping + SR_DaytimeSleepiness + SR_DaytimeSleepiness_adjBMI + ACC_DIURNAL_INACT_RAW_SIN
             F5=~ NA*Per90 + SnoringBMIadj + sleep_apnea + average_spo2 + minimum_spo2
	     F6 =~ NA*ACC_N_SLEEP_EPISODES_RAW_SIN + ACC_SLEEP_EFF_RAW_SIN
	     
F1~~F2
F1~~F3
F1~~F4
F1~~F5
F1~~F6
F2~~F3
F2~~F4
F2~~F5
F2~~F6
F3~~F4
F3~~F5
F3~~F6
F4~~F5
F4~~F6
F5~~F6
F1 ~ SNP
F2 ~ SNP
F3 ~ SNP
F4 ~ SNP
F5 ~ SNP
F6 ~ SNP
'

print ("model_building sucessful ")
gwas_result<-userGWAS(covstruc = LDSCoutput, SNPs = test_sumstats, estimation = "DWLS", model = model, printwarn = TRUE, 
    sub=c("F1~SNP","F2~SNP","F3~SNP","F4~SNP","F5~SNP","F6~SNP"), 
    cores = 1, toler = FALSE, SNPSE = FALSE, 
    parallel = FALSE,GC="standard",MPI=FALSE,
    smooth_check=TRUE,fix_measurement=TRUE)


save(gwas_result,file=paste0(wdir,"result/sleep_factorGWAS_chr${chr}_20240523.RData"))


write.table(as.data.frame(gwas_result[[1]]),file=paste0(wdir,"result/sleep_factor1_chr${chr}_20240523.txt"), sep="\t",col.names=T,row.names=F,quote=F)

write.table(as.data.frame(gwas_result[[2]]),file=paste0(wdir,"result/sleep_factor2_chr${chr}_20240523.txt"), sep="\t",col.names=T,row.names=F,quote=F)

write.table(as.data.frame(gwas_result[[3]]),file=paste0(wdir,"result/sleep_factor3_chr${chr}_20240523.txt"),sep="\t",col.names=T,row.names=F,quote=F)

write.table(as.data.frame(gwas_result[[4]]),file=paste0(wdir,"result/sleep_factor4_chr${chr}_2020523.txt"),sep="\t",col.names=T,row.names=F,quote=F)

write.table(as.data.frame(gwas_result[[5]]),file=paste0(wdir,"result/sleep_factor5_chr${chr}_20240523.txt"), sep="\t",col.names=T,row.names=F,quote=F)

write.table(as.data.frame(gwas_result[[6]]),file=paste0(wdir,"result/sleep_factor6_chr${chr}_20240523.txt"), sep="\t",col.names=T,row.names=F,quote=F)


EOF
chmod +x gwas_${chr}.R

singularity exec --bind /gpfs/data/fs71407/epiwielsch/ukb_2/application /gpfs/data/fs71407/epiwielsch/shortcake.sif Rscript gwas_${chr}.R
