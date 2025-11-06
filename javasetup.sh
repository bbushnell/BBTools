#!/bin/bash

# javasetup.sh v1.2
# Parses Java command-line arguments and sets up paths
# Authors: Brian Bushnell, Doug Jacobsen, Alex Copeland, Bryce Foster, Isla
# Date: November 6, 2025

# Source memory detection script
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/memdetect.sh"

# Initialize global variables
XMX=""
XMS=""
EA="-ea"
EOOM=""
SIMD=""
SIMD_AUTO=""
json=0
silent=0

# Detect if CPU supports AVX2 (or ARM NEON equivalent)
function detectCPUVectorSupport() {
    if [ -f /proc/cpuinfo ]; then
        # x86_64: Check for AVX2 (256-bit minimum)
        if grep -q "avx2" /proc/cpuinfo; then
            return 0
        fi
        # x86_64: Check for AVX-512 (even better)
        if grep -q "avx512" /proc/cpuinfo; then
            return 0
        fi
        # ARM: Check for SVE (scalable, can support 256-bit+)
        # Note: SVE vector length is runtime-configurable, but typically 256-bit+
        if grep -q "sve" /proc/cpuinfo; then
            return 0
        fi
    fi
    
    # macOS x86_64: Use sysctl for AVX2
    if command -v sysctl >/dev/null 2>&1; then
        if sysctl -a 2>/dev/null | grep -q "hw.optional.avx2.*: 1"; then
            return 0
        fi
        # TODO: macOS ARM (Apple Silicon) - need to verify SVE support detection
    fi
    
    return 1
}

# Detect Java version (need 17+)
function detectJavaVersion() {
    if ! command -v java >/dev/null 2>&1; then
        return 1  # Java not found
    fi
    
    # Get version string
    local version_output=$(java --version 2>&1 | head -n 1)
    
    # Extract major version number
    # Handles formats like:
    #   "openjdk 25 2025-09-16"
    #   "openjdk 17.0.1 2021-10-19"
    #   "java version "1.8.0_291""
    local major_version=0
    
    if echo "$version_output" | grep -q "openjdk [0-9]"; then
        # Modern format: "openjdk 25" or "openjdk 17.0.1"
        major_version=$(echo "$version_output" | sed -n 's/.*openjdk \([0-9]\+\).*/\1/p')
    elif echo "$version_output" | grep -q "java version"; then
        # Old format: "java version "1.8.0_291""
        major_version=$(echo "$version_output" | sed -n 's/.*"\(1\.\)\?\([0-9]\+\).*/\2/p')
    fi
    
    # Check if version >= 17
    if [ "$major_version" -ge 17 ] 2>/dev/null; then
        return 0
    fi
    
    return 1
}

# Auto-detect SIMD support
function autoDetectSIMD() {
    if detectCPUVectorSupport && detectJavaVersion; then
        SIMD="--add-modules jdk.incubator.vector"
        SIMD_AUTO="simd"
        return 0
    fi
    return 1
}

# Normalize memory specification with intelligent defaults
# Args: 
#   $1 - memory string (e.g., "2", "4g", "512m")
#   $2 - prefix ("-Xmx" or "-Xms")
# Returns: normalized Java memory flag
normalizeMemory() {
    local mem="$1"
    local prefix="$2"
    
    # Strip excessive leading dashes (handle ---, ----, etc.)
    while [[ "$mem" == --* ]]; do
        mem="${mem#-}"
    done
    
    # Strip any existing prefix variations
    mem="${mem#-Xmx}"
    mem="${mem#-Xms}"
    mem="${mem#Xmx}"
    mem="${mem#Xms}"
    mem="${mem#-xmx}"
    mem="${mem#-xms}"
    mem="${mem#xmx}"
    mem="${mem#xms}"
    
    # Check if already has suffix (case-insensitive)
    if [[ "$mem" =~ ^[0-9]+[gGmMkK]$ ]]; then
        # Already has suffix, normalize to lowercase
        mem=$(echo "$mem" | tr '[:upper:]' '[:lower:]')
        echo "${prefix}${mem}"
        return
    fi
    
    # Pure number - apply heuristic
    if [[ "$mem" =~ ^[0-9]+$ ]]; then
        # Get TOTAL installed physical memory in GB (not available)
        local physicalMemGB=0
        if [ -e /proc/meminfo ]; then
            local memTotal=$(grep '^MemTotal:' /proc/meminfo | awk '{print $2}')
            physicalMemGB=$((memTotal / 1024 / 1024))
        elif command -v sysctl >/dev/null 2>&1; then
            local totalMemBytes=$(sysctl -n hw.memsize 2>/dev/null)
            physicalMemGB=$((totalMemBytes / 1024 / 1024 / 1024))
        fi
        
        # Calculate 90% of physical memory (round up)
        local ninetyPercent=$(((physicalMemGB * 9 + 9) / 10))
        
        # Decision logic:
        # - Values < 20: always GB (JVM can't bootstrap with <20MB)
        # - Values > 1000: always MB (terabyte users should know better)
        # - Values 20-1000: GB if <= 90% of physical memory, else MB
        if [ "$mem" -lt 20 ]; then
            echo "${prefix}${mem}g"
        elif [ "$mem" -gt 1000 ]; then
            echo "${prefix}${mem}m"
        elif [ "$physicalMemGB" -eq 0 ] || [ "$mem" -le "$ninetyPercent" ]; then
            # If we can't detect memory, assume GB (safer)
            echo "${prefix}${mem}g"
        else
            echo "${prefix}${mem}m"
        fi
        return
    fi
    
    # Fallback: use as-is with prefix
    echo "${prefix}${mem}"
}

# Parse Java memory and other flags
# Arguments:
#   All command-line arguments
function parseJavaArgs() {
    local setxmx=0
    local setxms=0
    local defaultMem="4g"
    local memPercent=84
    local memMode="auto"
    local simd_specified=0
    
    # Process all arguments
    for arg in "$@"; do
        # Tool-specific memory settings
        if [ "${arg%%=*}" = "--mem" ]; then
            defaultMem="$(echo "$arg" | cut -d= -f2)"
        elif [ "${arg%%=*}" = "--percent" ]; then
            memPercent="$(echo "$arg" | cut -d= -f2)"
        elif [ "${arg%%=*}" = "--mode" ]; then
            memMode="$(echo "$arg" | cut -d= -f2)"
            
        # Fix broken Xmx flags
        elif [[ "$arg" =~ ^-*[xX][mM][xX]=?([0-9].*)$ ]]; then
            XMX=$(normalizeMemory "${BASH_REMATCH[1]}" "-Xmx")
            setxmx=1
        elif [[ "$arg" =~ ^-*[xX][mM][sS]=?([0-9].*)$ ]]; then
            XMS=$(normalizeMemory "${BASH_REMATCH[1]}" "-Xms")
            setxms=1
        
        # Assertion settings
        elif [ "$arg" = "-da" ] || [ "$arg" = "-ea" ]; then
            EA="$arg"
        elif [ "$arg" = "da" ] || [ "$arg" = "ea" ]; then
            EA="-$arg"
        
        # Out of memory handling
        elif [ "$arg" = "ExitOnOutOfMemoryError" ] || [ "$arg" = "exitonoutofmemoryerror" ] || [ "$arg" = "eoom" ]; then
            EOOM="-XX:+ExitOnOutOfMemoryError"
        elif [ "$arg" = "-ExitOnOutOfMemoryError" ] || [ "$arg" = "-exitonoutofmemoryerror" ] || [ "$arg" = "-eoom" ]; then
            EOOM="-XX:+ExitOnOutOfMemoryError"
        
        # SIMD instructions
        elif [ "$arg" = "simd" ] || [ "$arg" = "SIMD" ] || [ "$arg" = "simd=t" ] || [ "$arg" = "simd=true" ]; then
            SIMD="--add-modules jdk.incubator.vector"
            SIMD_AUTO="simd"
            simd_specified=1
        elif [ "$arg" = "simd=f" ] || [ "$arg" = "simd=false" ] || [ "$arg" = "nosimd" ]; then
            SIMD=""
            SIMD_AUTO=""
            simd_specified=1
        
        # Output format
        elif [ "$arg" = "json" ] || [ "$arg" = "json=t" ] || [ "$arg" = "json=true" ] || [ "$arg" = "format=json" ]; then
            json=1
        
        # Silence output
        elif [ "$arg" = "silent" ] || [ "$arg" = "silent=t" ] || [ "$arg" = "silent=true" ]; then
            silent=1
        fi
    done
    
    # Auto-detect SIMD if not explicitly specified
    if [ "$simd_specified" = "0" ]; then
        autoDetectSIMD
    fi
    
    # If Xmx was set but Xms wasn't, make Xms = Xmx
    if [ "$setxmx" = "1" ] && [ "$setxms" = "0" ]; then
        local substring=$(echo $XMX | cut -d'x' -f 2)
        XMS="-Xms$substring"
    # If Xms was set but Xmx wasn't, make Xmx = Xms
    elif [ "$setxmx" = "0" ] && [ "$setxms" = "1" ]; then
        local substring=$(echo $XMS | cut -d's' -f 2)
        XMX="-Xmx$substring"
    fi
    
    # If no Xmx was specified, determine it automatically
    if [ "$setxmx" = "0" ]; then
        detectMemory "$defaultMem" "$memPercent" "$memMode"
        XMX="-Xmx${RAM}m"
        XMS="-Xms${RAM}m"
    fi
    z="$XMX"
    z2="$XMS"
}

# Setup environment paths based on the execution environment
function setEnvironment() {
    if [ "$SHIFTER_RUNTIME" = "1" ]; then
        shifter=1
    elif [ -n "$EC2_HOME" ]; then
        PATH=/test1/binaries/bgzip:$PATH
        PATH=/test1/binaries/lbzip2/bin:$PATH
        PATH=/test1/binaries/sambamba:$PATH
        PATH=/test1/binaries/pigz2/pigz-2.4:$PATH
    elif [ -n "$NERSC_HOST" ]; then
        PATH=/global/cfs/cdirs/bbtools/bgzip:$PATH
        PATH=/global/cfs/cdirs/bbtools/lbzip2/bin:$PATH
        PATH=/global/cfs/cdirs/bbtools/samtools116/samtools-1.16.1:$PATH
        PATH=/global/cfs/cdirs/bbtools/java/jdk-17/bin:$PATH
        PATH=/global/cfs/cdirs/bbtools/pigz2/pigz-2.4:$PATH
    fi
}

# Get Java command with all the appropriate flags
# Arguments:
#   $@ - All command-line arguments
# Returns:
#   Echoes the complete Java command
function getJavaCommand() {
    parseJavaArgs "$@"
    setEnvironment
    
    local JAVA_CMD="java $EA $EOOM $XMX $XMS $SIMD"
    echo "$JAVA_CMD"
}

# Check if this script is being sourced or run directly
if [ "$0" != "$BASH_SOURCE" ] && [ "$BASH_SOURCE" != "" ]; then
    :
else
    getJavaCommand "$@"
fi
