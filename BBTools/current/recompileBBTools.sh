#!/bin/bash
# BBTools Recompilation Script
# Location: current/recompileBBTools.sh

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo "================================"
echo "BBTools Recompilation Script"
echo "================================"
echo ""
echo -e "${YELLOW}WARNING: Recompilation is generally NOT necessary!${NC}"
echo "BBTools comes with pre-compiled class files that are optimized"
echo "and ready to use. Only recompile if you have modified the Java source."
echo ""
read -p "Are you sure you want to recompile? (yes/no): " CONFIRM

if [ "$CONFIRM" != "yes" ]; then
    echo "Recompilation cancelled."
    exit 0
fi

# Check Java version
echo ""
echo "Checking Java version..."
java_version=$(java -version 2>&1 | head -n 1 | cut -d'"' -f2 | cut -d'.' -f1)
if [ "$java_version" -lt 17 ]; then
    echo -e "${RED}Error: Java 17 or higher is required${NC}"
    echo "Your version: $(java -version 2>&1 | head -n 1)"
    exit 1
fi
echo "Java version OK: $(java -version 2>&1 | head -n 1)"

# Check if we're in the right directory
if [ ! -f "driver/Driver.java" ]; then
    echo -e "${RED}Error: This script must be run from the 'current' directory${NC}"
    echo "Current directory: $(pwd)"
    exit 1
fi

# Clean old class files
echo ""
echo "Removing old class files..."
find . -name "*.class" -delete
echo "Old class files removed."

# Compile
echo ""
echo "Compiling BBTools (this may take 1-2 minutes)..."
echo "Including: */*.java and */*/*.java"
echo ""

# Compile with incubator.vector for SIMD optimizations
javac --add-modules jdk.incubator.vector \
      --enable-preview \
      -source 17 -target 17 \
      -cp ".:../jni/*:../resources/*" \
      */*.java */*/*.java 2>&1 | tee compile.log

# Check if compilation succeeded
if [ ${PIPESTATUS[0]} -eq 0 ]; then
    echo ""
    echo -e "${GREEN}Compilation successful!${NC}"
    
    # Count class files
    class_count=$(find . -name "*.class" | wc -l)
    echo "Created $class_count class files"
    echo ""
    echo -e "${GREEN}BBTools has been recompiled successfully.${NC}"
    echo "The class files are in their respective package directories."
    echo "You can now use BBTools normally with the shell scripts."
    
    # Clean up
    rm -f compile.log
else
    echo ""
    echo -e "${RED}Compilation failed!${NC}"
    echo "Check compile.log for error details."
    echo ""
    echo "Common issues:"
    echo "  - Missing Java 17+ (you have: $java_version)"
    echo "  - Missing jni libraries in ../jni/"
    echo "  - Syntax errors in modified Java files"
    exit 1
fi

echo ""
echo -e "${YELLOW}Note:${NC} The recompiled bytecode may be slower than the"
echo "original Eclipse-compiled version. For best performance,"
echo "use the pre-compiled class files from the official release."