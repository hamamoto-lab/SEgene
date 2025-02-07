#!/bin/bash
# generate_file_list.sh
#
# This script generates a list of full paths for all BAM files within a specified directory,
# including its subdirectories, and writes the list to "output/bam_list.txt".
#
# Usage:
#   ./generate_file_list.sh /path/to/your/directory
#
# Example:
#   ./generate_file_list.sh /home/user/data/bam_files

# Check if exactly one argument (the target directory) is provided.
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

directory="$1"

# Verify that the provided directory exists.
if [ ! -d "$directory" ]; then
    echo "Error: Directory '$directory' does not exist."
    exit 1
fi

# Define the output file name and ensure the output directory exists.
output_file="output/bam_list.txt"
mkdir -p "$(dirname "$output_file")"

# Find all files ending with .bam in the specified directory (and subdirectories),
# sort the results, and write the full paths to the output file.
find "$directory" -type f -name "*.bam" | sort > "$output_file"

# Count the number of lines (BAM files) in the output file.
file_count=$(wc -l < "$output_file")

echo "The list of BAM file full paths has been created in '${output_file}'."
echo "Total number of BAM files: $file_count"
