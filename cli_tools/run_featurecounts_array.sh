#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -jc your_job_class   # Optional: specify your job class if needed
#$ -o ./output/featurecounts_output/logs/
#$ -e ./output/featurecounts_output/logs/

# Ensure that the log directory exists
log_dir="./output/featurecounts_output/logs"
mkdir -p "$log_dir"

# Retrieve the BAM file corresponding to the current SGE task ID from 'output/bam_list.txt'
bam_file=$(sed -n "${SGE_TASK_ID}p" output/bam_list.txt)

# Create the output directory for featureCounts results if it does not exist
output_dir="./output/featurecounts_output"
mkdir -p "$output_dir"

# Set the output file name based on the input file name (removing the .bam extension)
base_name=$(basename "$bam_file" .bam)
output_file="$output_dir/featurecounts_${base_name}.txt"

# Set the annotation file to the default GTF file in the output folder
annotation_file="output/modified_for_featurecounts.gtf"

# Run featureCounts on the current BAM file using the annotation file
featureCounts -a "$annotation_file" -o "$output_file" -F GTF -T 1 -p "$bam_file"

# Log the processing
echo "Processed $bam_file with featureCounts"
