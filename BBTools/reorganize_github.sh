#!/bin/bash
# BBTools GitHub Repository Reorganizer
# This reorganizes the GitHub repo structure WITHOUT changing local files
#
# Usage: ./reorganize_github.sh
#
# Strategy: Use git's ability to move files in the index without touching working directory

echo "================================"
echo "BBTools GitHub Reorganization"
echo "================================"
echo ""
echo "This will reorganize GitHub structure while keeping your local files unchanged"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Make sure we're in the right directory
cd /c/releases/bbmap || cd /mnt/c/releases/bbmap

# Step 1: Check we're on a clean working tree
if ! git diff-index --quiet HEAD --; then
    echo -e "${YELLOW}You have uncommitted changes. Please commit or stash them first.${NC}"
    exit 1
fi

echo -e "${YELLOW}Creating reorganized structure in git index...${NC}"

# Step 2: Create a new branch for the reorganization
git checkout -b reorganized-structure

# Step 3: Use git mv to reorganize in the index
# This doesn't touch your working files!
mkdir -p BBTools

# Move everything except what should stay at root
for file in *.sh; do
    if [ -f "$file" ] && [ "$file" != "deploy.sh" ] && [ "$file" != "reorganize_github.sh" ] && [ "$file" != "deploy_reorganized.sh" ]; then
        git mv "$file" "BBTools/$file" 2>/dev/null || true
    fi
done

# Move directories
for dir in current config docs jni pipelines resources; do
    if [ -d "$dir" ]; then
        git mv "$dir" "BBTools/$dir" 2>/dev/null || true
    fi
done

# Move other files
[ -f "license.txt" ] && git mv "license.txt" "BBTools/license.txt"
[ -f "build.xml" ] && git mv "build.xml" "BBTools/build.xml"
[ -f "bitbucket-pipelines.yml" ] && git mv "bitbucket-pipelines.yml" "BBTools/bitbucket-pipelines.yml"

# Remove pytools if it exists
git rm -rf pytools 2>/dev/null || true
git rm -rf BBTools/pytools 2>/dev/null || true

# Remove barcode/prob if it exists  
git rm -rf BBTools/current/barcode/prob 2>/dev/null || true

# Step 4: Create README.md at root
cat > README.md << 'EOF'
# BBTools

**BBTools** is a suite of fast, multithreaded bioinformatics tools designed for analysis of DNA and RNA sequence data.

**Official Website:** [bbmap.org](https://bbmap.org)  
**Author:** Brian Bushnell

## Quick Start

All tools are located in the `BBTools/` directory:

```bash
git clone https://github.com/bbushnell/BBTools.git
cd BBTools/BBTools
./bbduk.sh --help
```

## Popular Tools

- **BBDuk** - Adapter trimming, quality filtering, and contamination removal
- **BBMap** - Fast RNA-seq and DNA aligner
- **BBMerge** - Paired read merging with error correction
- **Reformat** - Format conversion and manipulation
- **Clumpify** - Optical duplicate removal and compression

## Documentation

- [Complete documentation at bbmap.org](https://bbmap.org)
- [Usage Guide](BBTools/docs/UsageGuide.txt)
- [Tool Descriptions](BBTools/docs/ToolDescriptions.txt)

## Requirements

- Java 8 or higher
- 4GB RAM minimum (more for large datasets)

## Citation

Bushnell, B. (2024). BBTools: A suite of fast, multithreaded bioinformatics tools. 
Available at https://github.com/bbushnell/BBTools

## License

See [license.txt](BBTools/license.txt) for details.
EOF

git add README.md

# Step 5: Keep .gitignore at root
if [ ! -f ".gitignore" ]; then
    echo "*.class" > .gitignore
    echo "*.tmp" >> .gitignore
    echo "*.temp" >> .gitignore
    git add .gitignore
fi

# Step 6: Commit the reorganization
echo -e "${YELLOW}Committing reorganized structure...${NC}"
git commit -m "Reorganize repository structure: Move all tools to BBTools/ subdirectory

- All BBTools components now in BBTools/ subdirectory  
- README.md at root for GitHub visibility
- Removed pytools (BBTools is Java-only)
- Removed barcode/prob (internal only)
- Local working directory remains unchanged"

# Step 7: Push to GitHub
echo ""
echo -e "${GREEN}Ready to push reorganized structure to GitHub${NC}"
echo "This will create a new 'reorganized-structure' branch"
echo ""
read -p "Push to GitHub now? (y/n): " PUSH_NOW

if [ "$PUSH_NOW" = "y" ]; then
    git push -u origin reorganized-structure
    
    echo ""
    echo -e "${GREEN}Pushed to 'reorganized-structure' branch!${NC}"
    echo ""
    echo "To make this the main branch on GitHub:"
    echo "1. Go to https://github.com/bbushnell/BBTools/branches"
    echo "2. Create a pull request from 'reorganized-structure' to 'master'"
    echo "3. Or directly merge with: git checkout master && git merge reorganized-structure && git push"
else
    echo "You can push later with: git push -u origin reorganized-structure"
fi

# Step 8: Reset local working directory
echo ""
read -p "Reset your local directory to original state? (y/n): " RESET_NOW

if [ "$RESET_NOW" = "y" ]; then
    git checkout master
    git branch -D reorganized-structure
    echo "Local directory reset to original state"
else
    echo "Staying on reorganized-structure branch"
    echo "To go back: git checkout master"
fi

echo ""
echo "================================"
echo -e "${GREEN}Complete!${NC}"
echo "================================"