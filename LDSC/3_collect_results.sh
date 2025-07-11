#!/bin/bash

# Directory where result files are 
results_dir="/gpfs/data/fs71407/epiwielsch/ukb_2/application/exploratory_analysis/LDSCORE_regression/"

# Output file for the concatenated results
output_file="Exploratory_analysis_genetic_correlations_summary.tsv"

# Temporary file 
temp_file="temp_summary_results.tsv"

# Ensure the output file is empty
> "$output_file"

# First, we find all .log files in the results directory
# Then, for each file, we extract the summary table and append it to the temp_file
find "$results_dir" -name "*.log" | while read -r result_file; do
    awk '/Summary of Genetic Correlation Results/,/Analysis finished/' "$result_file" | \
    sed '/Summary of Genetic Correlation Results/d;/Analysis finished/d' >> "$temp_file"
done

# Now, clean up the data to ensure the header is included only once and append to the final output file
# Print the header
head -n 1 "$temp_file" > "$output_file"
# Print the rest, excluding the header from subsequent sections
awk 'NR > 1 {print}' "$temp_file" | sort | uniq >> "$output_file"

# Remove the temporary file
rm "$temp_file"

echo "All summaries have been concatenated into $output_file"

