conda activate mtag
sumdir='input_sumstats/'

/home/fs71707/mwielsch1/mtag/./mtag.py \
  --sumstats ${sumdir}sleep_factor4_orig_SNPs_for_MTAG.txt,${sumdir}SWB_SNPs_for_MTAG.txt,${sumdir}NEUR_SNPs_for_MTAG.txt,${sumdir}ADHD_SNPs_for_MTAG.txt,${sumdir}MDD_SNPs_for_MTAG.txt,${sumdir}PTST_SNPs_for_MTAG.txt \
  --out result/sleep_factor_4_ORIG_inclPsych \
  --n_min 0.0 \
  --cores 12 \
  --force \
  --stream_stdout &



### -------------------- factor 2

/home/fs71707/mwielsch1/mtag/./mtag.py \
  --sumstats ${sumdir}sleep_factor2_SNPs_for_MTAG.txt,${sumdir}SWB_SNPs_for_MTAG.txt,${sumdir}NEUR_SNPs_for_MTAG.txt,${sumdir}ADHD_SNPs_for_MTAG.txt,${sumdir}MDD_SNPs_for_MTAG.txt,${sumdir}PTST_SNPs_for_MTAG.txt \
  --out result/sleep_factor_2_inclPsych \
  --n_min 0.0 \
  --cores 12 \
  --force \
  --stream_stdout &
