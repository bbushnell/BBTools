#!/bin/bash

# Exceptions - these stay as aligner.
exceptions=("alltoall.sh" "alignrandom.sh" "microalign.sh")

cd /mnt/c/releases/bbmap

# Loop through all shell scripts containing aligner references
for file in *aligner*.sh; do
    # Check if this file is in the exceptions list
    skip=false
    for exception in "${exceptions[@]}"; do
        if [[ "$file" == "$exception" ]]; then
            skip=true
            echo "Skipping $file (exception)"
            break
        fi
    done

    # If not an exception, do the replacement
    if [[ "$skip" == false ]]; then
        echo "Updating $file: aligner. -> idaligner."
        sed -i 's/\$CP aligner\./\$CP idaligner\./g' "$file"
    fi
done

echo "Done!"
