#!/bin/bash


reference_file="../plink_dataset/per_chrom_bed/gSEM_SCORE_data.bim"



# Step 1: Create a key-value pair file from the reference, with rsID as key and allele as value
awk '{print $2 "\t" $5}' $reference_file > ref_key_value.txt


while read trait; do

target_file="score_input/${trait}.txt"
output_file="score_input/aligned_${trait}.txt"

awk -v OFS='\t' '{print $2 , $4, $6 }' ${target_file} > inter_FILE

# Step 2: Process the target file against the key-value pair file
awk 'BEGIN {
    FS=OFS="\t"
    print "ID", "effect_allele", "beta" > "'$output_file'"
}
NR==FNR {
    allele[$1] = $2
    next
}
($1 in allele) {
    if ($2 != allele[$1]) {
        print "Flipping allele for " $1 > "/dev/stderr"
        $3 = -$3  # Negate the beta value
        $2 = allele[$1]  # Update the effect allele to match the reference
    }
    print $1, $2, $3 > "'$output_file'"
}' ref_key_value.txt inter_FILE

rm inter_FILE
done < traits