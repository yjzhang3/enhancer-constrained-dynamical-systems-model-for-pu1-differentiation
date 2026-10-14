#!/bin/bash
output_file="combined_Y.txt"
> "$output_file"

for file in ODENOTSfit_PS_fsolve_extraP_[0-9][0-9][0-9]PCn22_Y.txt; do
    if [[ -f "$file" ]]; then
        cat "$file" >> "$output_file"
        echo "" >> "$output_file"
    fi
done

echo "All files combined into $output_file"
