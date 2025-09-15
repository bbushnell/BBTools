#!/bin/bash
# BBTools Deployment Script
# Automates deployment to GitHub and website
#
# Usage: ./deploy.sh
#
# This script will:
# 1. Remove sensitive/dev-only packages (barcode/prob)
# 2. Check for uncommitted changes
# 3. Optionally update version numbers
# 4. Push to GitHub
# 5. Optionally update website
#
# Prerequisites:
# - Git configured with GitHub credentials
# - Write access to bbushnell/BBTools repository
# - (Optional) Access to bbmap-website repository

echo "================================"
echo "BBTools Deployment Script v1.0"
echo "================================"
echo ""

# Colors for output (makes errors/success obvious)
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration paths
# Detect if running in WSL or Windows Git Bash
if [ -d "/mnt/c" ]; then
    # WSL paths
    BBTOOLS_DIR="/mnt/c/releases/bbmap"
    WEBSITE_DIR="/mnt/c/playground/Chloe/bbmap_website"
    GH_CLI="/mnt/c/Program Files/GitHub CLI/gh.exe"
else
    # Git Bash paths
    BBTOOLS_DIR="/c/releases/bbmap"
    WEBSITE_DIR="/c/playground/Chloe/bbmap_website"
    GH_CLI="/c/Program Files/GitHub CLI/gh.exe"
fi

# CRITICAL: Remove dev-only packages that should NEVER be deployed
FORBIDDEN_DIR="$BBTOOLS_DIR/current/barcode/prob"
if [ -d "$FORBIDDEN_DIR" ]; then
    echo -e "${RED}WARNING: Found barcode/prob package - removing before deployment${NC}"
    rm -rf "$FORBIDDEN_DIR"
    echo "Removed: $FORBIDDEN_DIR"
fi

# Remove Python tools - BBTools is Java-only!
PYTOOLS_DIR="$BBTOOLS_DIR/pytools"
if [ -d "$PYTOOLS_DIR" ]; then
    echo -e "${RED}WARNING: Found pytools directory - removing before deployment${NC}"
    rm -rf "$PYTOOLS_DIR"
    echo "Removed: $PYTOOLS_DIR (BBTools is Java-only, no Python!)"
fi

# Ensure we're in the BBTools directory
if [ "$PWD" != "$BBTOOLS_DIR" ]; then
    echo -e "${YELLOW}Switching to BBTools directory...${NC}"
    cd "$BBTOOLS_DIR"
fi

# Step 1: Check for uncommitted changes and add class files
# This prevents accidentally deploying incomplete work
echo "Checking for uncommitted changes..."

# First, force-add class files (they're gitignored)
echo "Adding compiled class files..."
git add -f current/*/*.class 2>/dev/null
git add -f current/*/*/*.class 2>/dev/null

# Now check for all uncommitted changes including class files
if ! git diff-index --quiet HEAD --; then
    echo -e "${YELLOW}Uncommitted changes found:${NC}"
    echo ""
    git status --short
    echo ""
    read -p "Do you want to commit them now? (y/n): " COMMIT_NOW
    
    if [ "$COMMIT_NOW" = "y" ]; then
        read -p "Enter commit message: " COMMIT_MSG
        git add .
        git commit -m "$COMMIT_MSG"
    else
        echo -e "${RED}Aborting deployment. Please commit changes first.${NC}"
        exit 1
    fi
fi

# Step 2: Version management
# Shows current version and optionally updates it
echo ""
echo "Current version in README:"
CURRENT_VERSION=$(grep "Current Version:" README.md | sed 's/.*\*\*\(.*\)\*\*.*/\1/')
echo -e "${GREEN}$CURRENT_VERSION${NC}"
echo ""

# Ask about version update
read -p "Enter new version number (or press Enter to keep $CURRENT_VERSION): " NEW_VERSION

if [ ! -z "$NEW_VERSION" ]; then
    echo "Updating version to $NEW_VERSION..."
    
    # Update README.md
    sed -i "s/Current Version: \*\*.*\*\*/Current Version: **$NEW_VERSION**/" README.md
    
    # Update shell scripts if major version
    read -p "Update version in shell scripts? (y/n): " UPDATE_SCRIPTS
    if [ "$UPDATE_SCRIPTS" = "y" ]; then
        echo "This would update 252 shell scripts - not implemented yet"
        # TODO: Implement shell script version updates
    fi
    
    # Commit version change
    git add README.md
    git commit -m "Update version to $NEW_VERSION"
    
    VERSION_TO_USE="$NEW_VERSION"
else
    VERSION_TO_USE="$CURRENT_VERSION"
fi

# Step 3: Add class files and push to GitHub
# This is the main deployment - pushes all changes INCLUDING compiled class files
echo ""
echo -e "${YELLOW}Adding class files and pushing to GitHub...${NC}"

# Add all regular files
git add -A

# Explicitly add class files (in case gitignore excludes them)
echo "Adding compiled class files..."
git add -f current/*/*.class 2>/dev/null || true
git add -f current/*/*/*.class 2>/dev/null || true

# Commit if there are changes
git commit -m "Include compiled class files" 2>/dev/null || true

# Push to GitHub
git push origin master

if [ $? -eq 0 ]; then
    echo -e "${GREEN}Successfully pushed to GitHub${NC}"
else
    echo -e "${RED}Failed to push to GitHub${NC}"
    exit 1
fi

# Step 4: Website deployment (optional)
# Updates bbmap.org via GitHub -> Netlify auto-deploy
echo ""
read -p "Update bbmap.org website? (y/n): " UPDATE_WEBSITE

if [ "$UPDATE_WEBSITE" = "y" ]; then
    echo -e "${YELLOW}Updating website...${NC}"
    
    cd "$WEBSITE_DIR"
    
    # Update version on website if changed
    if [ ! -z "$NEW_VERSION" ]; then
        sed -i "s/<strong>Current Version:<\/strong> .*</<strong>Current Version:<\/strong> $NEW_VERSION</" index.html
        git add index.html
        git commit -m "Update version to $NEW_VERSION on website"
    fi
    
    # Check for other website changes
    if ! git diff-index --quiet HEAD --; then
        echo "Website has uncommitted changes:"
        git status --short
        read -p "Commit and deploy these changes? (y/n): " DEPLOY_WEB
        
        if [ "$DEPLOY_WEB" = "y" ]; then
            git add .
            read -p "Commit message for website changes: " WEB_MSG
            git commit -m "$WEB_MSG"
        fi
    fi
    
    # Push website
    git push origin main
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}Website pushed - Netlify will auto-deploy${NC}"
        echo "  Website will be live at https://bbmap.org in ~1-2 minutes"
    else
        echo -e "${RED}Failed to push website${NC}"
    fi
    
    cd "$BBTOOLS_DIR"
fi

# Summary
echo ""
echo "================================"
echo -e "${GREEN}Deployment Complete!${NC}"
echo "================================"
echo ""
echo "Version: $VERSION_TO_USE"
echo ""
echo "Checklist of things to update manually:"
echo "  [ ] JGI BBTools page"
echo "  [ ] SourceForge project page" 
echo "  [ ] Biostar/SeqAnswers announcement (if major version)"
echo "  [ ] Google Scholar profile (if new publication)"
echo ""
echo "URLs to check:"
echo "  - https://github.com/bbushnell/BBTools"
echo "  - https://bbmap.org"
echo "  - https://sourceforge.net/projects/bbmap/"
echo ""
echo "To run this script again: ./deploy.sh"
echo ""

# Make this script executable for next time
chmod +x "$BBTOOLS_DIR/deploy.sh" 2>/dev/null

echo "Deployment complete."