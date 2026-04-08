#!/bin/bash

# Fix the SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")" bug in all BBTools shell scripts
# This ensures scripts can find javasetup.sh even when invoked as "sh script.sh"

echo "Fixing script path resolution in BBTools shell scripts..."
echo ""

count=0
for script in *.sh; do
    if grep -q 'SCRIPT="\$0"' "$script"; then
        # Apply the fix: resolve $0 to absolute path before symlink handling
        sed -i 's|SCRIPT="\$0"|SCRIPT="$(cd "$(dirname "\$0")" \&\& pwd)/$(basename "\$0")"|g' "$script"
        echo "Fixed: $script"
        ((count++))
    fi
done

echo ""
echo "Fixed $count shell scripts"
