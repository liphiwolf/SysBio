#!/bin/bash
<<docstring
This script searches for all files in a given path, decompresses them if necessary, sorts them by chromosome and start position. The sorted file is saved, the original is deleted.
docstring

# Record the start time
start_time=$(date +%s)

# Check if the folder path is provided as an argument
if [ $# -eq 0 ]; then
  echo "Please provide the relative folder path of the files to be sorted:"
  read folder_path
else
  folder_path=$1
fi

# decompress all files
find "$folder_path" -name "*.gz" -exec gunzip {} +

# Get the list of BED files in the specified folder, excluding those containing "sorted"
bed_files=$(find "$folder_path" -maxdepth 1 -name "*.bed" ! -name "*sorted*")

# Count the number of BED files
num_files=$(echo "$bed_files" | wc -l)
counter=0

# Loop through each BED file
while IFS= read -r bed_file; do
  ((counter++))

  # Extract the file name without extension
  base_name=$(basename "$bed_file" .bed)

  # Sort the BED file using sort command and save the sorted output
  sort -k1,1V -k2,2n "$bed_file" > "$folder_path/$base_name.sorted.bed"

  # Compress the sorted BED file using gzip and overwrite if it already exists
  #gzip -f "$folder_path/$base_name.sorted.bed"

  # Print the path to the sorted BED file
  echo "Sorted BED file: $folder_path/$base_name.sorted.bed"

  # Display the progress as a fraction
  progress=$((counter * 100 / num_files))
  echo -ne "Progress: $counter/$num_files ($progress%)\r"
done <<< "$bed_files"

# Delete all non-sorted files in one operation
find "$folder_path" -maxdepth 1 -type f ! -name "*sorted*" -delete

# Record the end time
end_time=$(date +%s)

# Calculate and print the total execution time
total_time=$((end_time - start_time))
echo "Total execution time: $total_time seconds"
