#!/bin/bash

# Create or clear the output file
output_file="tumor_names.txt"
> $output_file

# Loop through each BAM file in the directory
for file in *.bam; do
    # Extract the full sample name (remove the .recal.bam extension)
    sample_name=$(basename "$file" .recal.bam)
    
    # Append the sample name to the output file
    echo $sample_name >> $output_file
done



     
    