#!/bin/bash
# bam_counts.sh
#
# This script takes a file list (each line containing the full path to a BAM file)
# as input, counts the number of reads for each BAM file using samtools, and
# outputs the results in CSV format to "output/bam_read_counts.csv" by default.
#
# Usage:
#   ./bam_counts.sh <file_list> [output_csv]
#
# Example:
#   ./bam_counts.sh output/bam_list.txt output/bam_read_counts.csv
#
# Note: The second argument (output CSV file path) is optional.
#       If omitted, the default output is "output/bam_read_counts.csv".

# Check if the number of arguments is valid.
if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Usage: $0 <file_list> [output_csv]"
    exit 1
fi

# Input file list containing full paths to BAM files.
file_list="$1"

# Set output CSV file path (use second argument if provided, else default).
if [ "$#" -eq 2 ]; then
    output_csv="$2"
else
    output_csv="output/bam_read_counts.csv"
fi

# Verify that the file list exists.
if [ ! -f "$file_list" ]; then
    echo "Error: File list '$file_list' does not exist."
    exit 1
fi

# Create the output directory if it doesn't exist.
mkdir -p "$(dirname "$output_csv")"

# Write CSV header.
echo "Sample,ReadCount" > "$output_csv"

# Process each BAM file in the file list.
while IFS= read -r bam; do
    # Check if the BAM file exists.
    if [ ! -f "$bam" ]; then
        echo "Warning: '$bam' not found. Skipping."
        continue
    fi

    # Extract the sample name by removing the .bam extension.
    sample_name=$(basename "$bam" .bam)

    # Count reads using samtools idxstats and awk.
    read_count=$(samtools idxstats "$bam" | awk '{s+=$3} END {print s}')

    # Append the result to the CSV file.
    echo "$sample_name,$read_count" >> "$output_csv"
done < "$file_list"

echo "Processing complete. Output file: $output_csv"
