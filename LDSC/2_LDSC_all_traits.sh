#!/bin/bash
#SBATCH -J LDSC
#SBATCH -N 1
#SBATCH -p skylake_0096
#SBATCH --qos=skylake_0096
#SBATCH --tasks-per-node=20


source ~/anaconda3/etc/profile.d/conda.sh
source activate ldsc

wdir='/gpfs/data/fs71407/epiwielsch/ukb_2/application/exploratory_analysis/LDSCORE_regression/'
indir='/gpfs/data/fs71407/epiwielsch/ukb_2/application/exploratory_analysis/LDSCORE_regression/munge_out/'


# The file containing the list of traits -- acutally filenames after munge script was run
traits_file="${wdir}traits_for_LDscore"

# Read each line from the traits file
while IFS= read -r trait; do
    # Extract the trait name for naming the output file
    trait_name=$(basename "${trait}" ".sumstats.gz")

    # Prepare the --rg argument with the current trait as the first entry
    rg_arg="${indir}${trait}"

    # Append the rest of the traits, excluding the current one
    while IFS= read -r other_trait; do
        if [ "$trait" != "$other_trait" ]; then
            rg_arg+=",${indir}${other_trait}"
        fi
    done < "$traits_file"

    # Define the output file name dynamically based on the current trait
    output="${wdir}LD_score_reg_${trait_name}_vs_others"

    # Run the ldsc.py command with the current setup
    /home/fs71407/epiwielsch/ldsc/./ldsc.py \
    --rg $rg_arg \
    --ref-ld-chr ${wdir}eur_w_ld_chr/ \
    --w-ld-chr ${wdir}eur_w_ld_chr/ \
    --out $output

    # Feedback
    echo "Processed $trait_name, output in $output"

done < "$traits_file"

