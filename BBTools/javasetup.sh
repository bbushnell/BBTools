#!/bin/bash

# javasetup.sh v1.0.3
# Parses Java command-line arguments and sets up paths
# Authors: Brian Bushnell, Doug Jacobsen, Alex Copeland, Bryce Foster, Claude AI
# Date: April 20, 2025

# Source memory detection script - try to be more compatible with different shells
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/memdetect.sh"

# Initialize global variables
XMX=""
XMS=""
EA="-ea"
EOOM=""
SIMD=""
json=0
silent=0

# Parse Java memory and other flags
# Arguments:
#   All command-line arguments
function parseJavaArgs() {
    local setxmx=0
    local setxms=0
    local defaultMem="4g"  # Default memory size
    local memPercent=84    # Default percentage
    local memMode="auto"   # Default mode (auto, partial, fixed)
    
    # Process all arguments
    for arg in "$@"; do
        # Tool-specific memory settings (process first)
        if [ "${arg%%=*}" = "--mem" ]; then
            defaultMem="$(echo "$arg" | cut -d= -f2)"
        elif [ "${arg%%=*}" = "--percent" ]; then
            memPercent="$(echo "$arg" | cut -d= -f2)"
        elif [ "${arg%%=*}" = "--mode" ]; then
            memMode="$(echo "$arg" | cut -d= -f2)"
            
        # Memory settings
        elif [ "${arg%%=*}" = "Xmx" ] || [ "${arg%%=*}" = "xmx" ]; then
            XMX="-Xmx$(echo "$arg" | cut -d= -f2)"
            setxmx=1
        elif [ "${arg%%=*}" = "-Xmx" ] || [ "${arg%%=*}" = "-xmx" ]; then
            XMX="-Xmx$(echo "$arg" | cut -d= -f2)"
            setxmx=1
        elif echo "$arg" | grep -q "^-Xmx"; then
            XMX="$arg"
            setxmx=1
        elif echo "$arg" | grep -q "^Xmx"; then
            XMX="-$arg"
            setxmx=1
        elif echo "$arg" | grep -q "^-Xms"; then
            XMS="$arg"
            setxms=1
        elif echo "$arg" | grep -q "^Xms"; then
            XMS="-$arg"
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
        
        # Output format
        elif [ "$arg" = "json" ] || [ "$arg" = "json=t" ] || [ "$arg" = "json=true" ] || [ "$arg" = "format=json" ]; then
            json=1
        
        # Silence output
        elif [ "$arg" = "silent" ] || [ "$arg" = "silent=t" ] || [ "$arg" = "silent=true" ]; then
            silent=1
        fi
    done
    
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
}

# Setup environment paths based on the execution environment
function setEnvironment() {
    # Check for specific execution environments
    
    if [ "$SHIFTER_RUNTIME" = "1" ]; then
        # Running in Shifter container
        shifter=1
    elif [ -n "$EC2_HOME" ]; then
        # Running on AWS
        PATH=/test1/binaries/bgzip:$PATH
        PATH=/test1/binaries/lbzip2/bin:$PATH
        PATH=/test1/binaries/sambamba:$PATH
        PATH=/test1/binaries/pigz2/pigz-2.4:$PATH
    elif [ -n "$NERSC_HOST" ]; then
        # Running on NERSC
        PATH=/global/cfs/cdirs/bbtools/bgzip:$PATH
        PATH=/global/cfs/cdirs/bbtools/lbzip2/bin:$PATH
        PATH=/global/cfs/cdirs/bbtools/samtools116/samtools-1.16.1:$PATH
        PATH=/global/cfs/cdirs/bbtools/java/jdk-17/bin:$PATH
        PATH=/global/cfs/cdirs/bbtools/pigz2/pigz-2.4:$PATH
    fi
    
    # Add other environment-specific setups here
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
# This is a more portable way than using BASH_SOURCE
if [ "$0" != "$BASH_SOURCE" ] && [ "$BASH_SOURCE" != "" ]; then
    # Being sourced
    :
else
    # Being run directly
    getJavaCommand "$@"
fi
