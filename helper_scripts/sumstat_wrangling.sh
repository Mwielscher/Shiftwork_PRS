#!/bin/bash

# Define the input file name
input_file="../input/accel_GWAS_all_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt"

# Define the list of traits in an array
traits=(
  "ACC_DIURNAL_INACT_RAW"
  "ACC_DIURNAL_INACT_RAW_SIN"
  "ACC_L5_TIME_RAW"
  "ACC_L5_TIME_RAW_SIN"
  "ACC_M10_TIME_RAW"
  "ACC_M10_TIME_RAW_SIN"
  "ACC_N_SLEEP_EPISODES_RAW"
  "ACC_N_SLEEP_EPISODES_RAW_SIN"
  "ACC_SLEEP_DUR_RAW"
  "ACC_SLEEP_DUR_RAW_SIN"
  "ACC_SLEEP_DUR_SD_RAW"
  "ACC_SLEEP_DUR_SD_RAW_SIN"
  "ACC_SLEEP_EFF_RAW"
  "ACC_SLEEP_EFF_RAW_SIN"
  "ACC_SLEEP_MIDP_RAW"
  "ACC_SLEEP_MIDP_RAW_SIN"
)

# Process each trait
for trait in "${traits[@]}"; do
  # Use awk to process the file
  awk -v trait="$trait" '
  BEGIN {
    FS=OFS="\t";
    print "RSID", "CHROM", "POS", "REF", "ALT", "BUILD", "BETA", "SE", "P", "A1FREQ";
  }
  NR==1 {
    for (i=1; i<=NF; i++) {
      if ($i == trait "_BETA") beta = i;
      if ($i == trait "_SE") se = i;
      if ($i == trait "_P") p = i;
      if ($i == trait "_A1FREQ") a1freq = i;
    }
  }
  NR>1 {
    print $1, $2, $3, $4, $5, "GRCh37", $(beta), $(se), $(p), $(a1freq);
  }
  ' "$input_file" > "${trait}_step1.txt"
done

