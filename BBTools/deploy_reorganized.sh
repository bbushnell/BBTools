#!/bin/bash
# BBTools Deployment Script with Repository Reorganization
# This creates a reorganized structure ONLY for GitHub, leaving local files untouched
#
# Usage: ./deploy_reorganized.sh
#
# This script will:
# 1. Create a temporary clone with reorganized structure
# 2. Move everything except README.md into BBTools/ subdirectory
# 3. Push the reorganized structure to GitHub
# 4. Clean up temporary files
#
# Your local filesystem remains unchanged!

echo "============================================"
echo "BBTools GitHub Reorganization Deployment"
echo "============================================"
echo ""
echo "This will reorganize GitHub WITHOUT changing your local files!"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Detect if running in WSL or Windows Git Bash
if [ -d "/mnt/c" ]; then
    # WSL paths
    BBTOOLS_DIR="/mnt/c/releases/bbmap"
    TEMP_DIR="/mnt/c/releases/bbmap_temp_deploy"
else
    # Git Bash paths
    BBTOOLS_DIR="/c/releases/bbmap"
    TEMP_DIR="/c/releases/bbmap_temp_deploy"
fi

# Step 1: Clean up any previous temp directory
if [ -d "$TEMP_DIR" ]; then
    echo -e "${YELLOW}Removing old temporary directory...${NC}"
    rm -rf "$TEMP_DIR"
fi

# Step 2: Create temporary clone
echo -e "${YELLOW}Creating temporary reorganized structure...${NC}"
mkdir -p "$TEMP_DIR"
cd "$TEMP_DIR"

# Initialize git repo
git init
git remote add origin git@github.com:bbushnell/BBTools.git

# Step 3: Create the new structure
echo "Creating BBTools/ subdirectory structure..."
mkdir -p BBTools

# Copy everything from original location EXCEPT certain files
echo "Copying files to reorganized structure..."
cd "$BBTOOLS_DIR"

# Copy all shell scripts
cp *.sh "$TEMP_DIR/BBTools/" 2>/dev/null || true

# Copy all directories
for dir in current config docs jni pipelines resources; do
    if [ -d "$dir" ]; then
        echo "  Copying $dir/..."
        cp -r "$dir" "$TEMP_DIR/BBTools/"
    fi
done

# Copy other important files
cp -r .gitignore "$TEMP_DIR/" 2>/dev/null || true
cp -r license.txt "$TEMP_DIR/BBTools/" 2>/dev/null || true
cp -r *.yml "$TEMP_DIR/BBTools/" 2>/dev/null || true
cp -r *.xml "$TEMP_DIR/BBTools/" 2>/dev/null || true

# Step 4: Create README.md at root
cd "$TEMP_DIR"
cat > README.md << 'EOF'
# BBTools

**BBTools** is a suite of fast, multithreaded bioinformatics tools designed for analysis of DNA and RNA sequence data.

**Official Website:** [bbmap.org](https://bbmap.org)

## Quick Navigation

All BBTools programs and resources are located in the `BBTools/` directory:
- Shell scripts (e.g., `bbduk.sh`, `bbmap.sh`) are in `BBTools/`
- Java source code is in `BBTools/current/`
- Resources and adapters are in `BBTools/resources/`
- Documentation is in `BBTools/docs/`

## Installation

1. **Download**: Clone this repository or download the latest release
2. **Java Requirement**: Ensure Java 8 or higher is installed
3. **Run**: Execute shell scripts from the `BBTools/` directory

Example:
```bash
cd BBTools
./bbduk.sh in=reads.fq out=clean.fq ref=resources/adapters.fa ktrim=r k=23 mink=11 hdist=1
```

## Popular Tools

- **BBDuk**: Adapter trimming, quality filtering, contamination removal
- **BBMap**: Fast splice-aware aligner for RNA-seq and DNA
- **BBMerge**: Paired read merging with error correction
- **Reformat**: Format conversion and basic operations
- **BBNorm**: Kmer-based error correction and normalization

## Documentation

- [Usage Guide](BBTools/docs/UsageGuide.txt)
- [Tool Descriptions](BBTools/docs/ToolDescriptions.txt)
- Full documentation at [bbmap.org](https://bbmap.org)

## Citation

If you use BBTools in your work, please cite:
> Bushnell, B. (2024). BBTools: A suite of fast, multithreaded bioinformatics tools. Available at https://github.com/bbushnell/BBTools

## Author

**Brian Bushnell** - Creator and maintainer of BBTools

## License

See [license.txt](BBTools/license.txt) for details.
EOF

# Step 5: Remove forbidden directories from the reorganized structure
echo -e "${YELLOW}Cleaning up forbidden directories...${NC}"
rm -rf "$TEMP_DIR/BBTools/current/barcode/prob" 2>/dev/null || true
rm -rf "$TEMP_DIR/BBTools/pytools" 2>/dev/null || true

# Step 6: Add compiled class files
echo "Adding compiled class files..."
cd "$TEMP_DIR/BBTools"
# Force-add class files from original location
find "$BBTOOLS_DIR/current" -name "*.class" -exec cp --parents {} . \; 2>/dev/null || true

# Step 7: Commit everything
cd "$TEMP_DIR"
git add -A
git add -f BBTools/current/*/*.class 2>/dev/null || true
git add -f BBTools/current/*/*/*.class 2>/dev/null || true

git commit -m "Reorganized repository structure: all tools now in BBTools/ subdirectory

- Moved all BBTools components to BBTools/ subdirectory
- README.md remains at repository root for GitHub display
- Removed pytools (BBTools is Java-only)
- Removed barcode/prob (internal development only)
- Improved organization for easier navigation"

# Step 8: Push to GitHub
echo ""
echo -e "${YELLOW}Ready to push reorganized structure to GitHub${NC}"
echo "This will reorganize the GitHub repository while leaving your local files unchanged."
echo ""
read -p "Push to GitHub now? (y/n): " PUSH_NOW

if [ "$PUSH_NOW" = "y" ]; then
    echo "Pushing to GitHub..."
    git push -f origin master
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}Successfully pushed reorganized structure to GitHub!${NC}"
        echo ""
        echo "GitHub repository has been reorganized:"
        echo "  - All tools are now in BBTools/ subdirectory"
        echo "  - README.md is at the root for visibility"
        echo "  - Your local files remain unchanged"
    else
        echo -e "${RED}Failed to push to GitHub${NC}"
        echo "You may need to pull first or check your credentials"
    fi
else
    echo "Deployment cancelled. Temporary files preserved at: $TEMP_DIR"
    echo "You can manually push later from that directory."
fi

# Step 9: Optional cleanup
echo ""
read -p "Remove temporary directory? (y/n): " CLEANUP

if [ "$CLEANUP" = "y" ]; then
    rm -rf "$TEMP_DIR"
    echo "Temporary directory removed."
else
    echo "Temporary directory preserved at: $TEMP_DIR"
fi

echo ""
echo "============================================"
echo -e "${GREEN}Reorganization Complete!${NC}"
echo "============================================"
echo ""
echo "Your local filesystem structure is UNCHANGED"
echo "GitHub repository now has organized structure"
echo ""
echo "Users will clone like this:"
echo "  git clone https://github.com/bbushnell/BBTools.git"
echo "  cd BBTools/BBTools"
echo "  ./bbduk.sh ..."
echo ""