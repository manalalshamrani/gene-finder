#!/bin/bash


# Loop through each directory and copy .gbk files
for dir in */; do
  if [[ -d "$dir" ]]; then
    cp "${dir}"*.gbk combined_prokka_results/ 2>/dev/null
  fi
done

echo "All .gbk files have been copied to combined_output."

