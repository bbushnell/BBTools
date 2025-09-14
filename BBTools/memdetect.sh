#!/bin/bash

# memdetect.sh v1.0.4
# Detects available memory for Java applications across various environments
# Authors: Brian Bushnell, Doug Jacobsen, Alex Copeland, Bryce Foster, Isla (Claude Sonnet)
# Date: May 24, 2025

# Constants
DEFAULT_MEM_MB=3200
DEFAULT_MEM_PERCENT_SHARED=45
DEFAULT_MEM_PERCENT_EXCLUSIVE=85
RESERVED_MEM_KB=500000  # 500MB reserved for OS and other processes

# External flags imported from javasetup.sh
# silent=0 # This will be imported

# Detect available memory in kilobytes
# Arguments:
#   $1 - Default memory allocation (optional: e.g. "3200m", "4g")
#   $2 - Percentage of memory to use (optional: e.g. 75)
#   $3 - Memory mode: "auto", "partial", or "fixed"
# Returns:
#   Sets RAM variable with memory in megabytes
function detectMemory() {
    RAM=0
    local defaultMem=$DEFAULT_MEM_MB
    local defaultMemKB=$(($defaultMem * 1024))
    local memPercent=$DEFAULT_MEM_PERCENT_SHARED
    local memMode="auto"  # Default: "auto", "partial", or "fixed"
    local isExclusiveNode=0
    local userSpecifiedMemWithUnit=""
    
    # Process custom memory specification if provided
    if [ $# -gt 0 ]; then
        userSpecifiedMemWithUnit=$1  # Save original with unit
        defaultMem=$1
        case $defaultMem in
            *g)
            defaultMem=$(echo $defaultMem | cut -d'g' -f 1)
            defaultMemKB=$(( $defaultMem * 1024 * 1024 ))
            ;;
            *m)
            defaultMem=$(echo $defaultMem | cut -d'm' -f 1)
            defaultMemKB=$(( $defaultMem * 1024 ))
            ;;
            *k)
            defaultMemKB=$(echo $defaultMem | cut -d'k' -f 1)
            ;;
        esac
    fi

    # Process memory percentage if provided
    if [ $# -gt 1 ]; then
        memPercent=$2
    fi
    
    # Process memory mode if provided
    if [ $# -gt 2 ]; then
        memMode=$3
    fi
    
    # For fixed mode, just use the default and return
    if [ "$memMode" = "fixed" ]; then
        # Handle the original memory specification with unit
        case $userSpecifiedMemWithUnit in
            *g)
                # If specified in gigabytes, convert to MB (1g = 1024m)
                local numericPart=$(echo "$userSpecifiedMemWithUnit" | sed 's/g$//')
                RAM=$(( numericPart * 1024 ))
                ;;
            *m)
                # If specified in megabytes, use directly
                RAM=$(echo "$userSpecifiedMemWithUnit" | sed 's/m$//')
                ;;
            *k)
                # If specified in kilobytes, convert to MB
                local numericPart=$(echo "$userSpecifiedMemWithUnit" | sed 's/k$//')
                RAM=$(( numericPart / 1024 ))
                ;;
            *)
                # If no unit, assume MB
                RAM=$defaultMem
                ;;
        esac
        
        #if [ "$silent" != "1" ]; then
            #echo "Using fixed memory allocation: ${RAM}MB" >&2
        #fi
        return 0
    fi
    
    # Check for user-specified memory
    if [ -n "$RQCMEM" ] && [ $RQCMEM -gt 0 ]; then
        #if [ "$silent" != "1" ]; then
            #echo "Using manually specified memory: ${RQCMEM}MB" >&2
        #fi
        RAM=$RQCMEM
        return 0
    fi
    
    # Check for ulimit constraints
    local ulimit=$(ulimit -v)
    ulimit="${ulimit:-0}"
    if [ "$ulimit" = "unlimited" ]; then ulimit=0; fi
    
    # Check for scheduler-specific memory constraints
    local slurmMem=0
    if [ -n "$SLURM_MEM_PER_NODE" ]; then
        slurmMem=$(( SLURM_MEM_PER_NODE * 1024 ))  # Convert MB to KB
    fi
    
    # Platform-specific memory detection
    if echo "$OSTYPE" | grep -q "darwin"; then
        detectMacMemory
    elif [ -e /proc/meminfo ]; then
        detectLinuxMemory
    else
        # Fall back to defaults
        RAM=$(($defaultMem))
        if [ "$silent" != "1" ]; then
            echo "WARNING: Cannot detect system memory. Using default $RAM MB." >&2
        fi
        return 0
    fi
    
    # Apply scheduler constraints if applicable
    if [ $slurmMem -gt 0 ]; then
        if [ $availableMemKB -gt $slurmMem ] || [ $availableMemKB -eq 0 ]; then
            availableMemKB=$slurmMem
        fi
    fi
    
    # Apply ulimit if applicable
    if [ $ulimit -gt 0 ]; then
        if [ $availableMemKB -gt $ulimit ] || [ $availableMemKB -eq 0 ]; then
            availableMemKB=$ulimit
        fi
    fi
    
    # Calculate final memory allocation
    if [ $availableMemKB -gt 0 ]; then
        # For "partial" mode, use a lower percentage
        if [ "$memMode" = "partial" ]; then
            # (available RAM - reserved) * percentage / 100 / 1024 = MB
            RAM=$(( ((availableMemKB - RESERVED_MEM_KB) * memPercent / 100) / 1024 ))
        else
            # For "auto" mode, use the full percentage
            RAM=$(( ((availableMemKB - RESERVED_MEM_KB) * memPercent / 100) / 1024 ))
        fi
    else
        # Fall back to defaults
        RAM=$(($defaultMem))
        if [ "$silent" != "1" ]; then
            echo "WARNING: Memory detection failed. Using default $RAM MB." >&2
        fi
    fi
    
    # Final sanity check - don't return negative or tiny values
    if [ $RAM -lt 100 ]; then
        RAM=$(($defaultMem))
        if [ "$silent" != "1" ]; then
            echo "WARNING: Detected memory too low. Using default $RAM MB." >&2
        fi
    fi
    
    return 0
}

# Detect available memory on Linux systems
function detectLinuxMemory() {
    # Get memory statistics from /proc/meminfo
    # This accounts for both free memory and reclaimable buffer/cache memory
    availableMemKB=0
    
    # Try using "MemAvailable" field first (most accurate, available in newer kernels)
    local memAvailable=$(grep '^MemAvailable:' /proc/meminfo | awk '{print $2}')
    
    if [ -n "$memAvailable" ] && [ $memAvailable -gt 0 ]; then
        # Modern kernels have MemAvailable which is the most accurate metric
        availableMemKB=$memAvailable
        #if [ "$silent" != "1" ]; then
            #echo "Detected ${availableMemKB}KB available memory (MemAvailable)" >&2
        #fi
    else
        # Older kernels - calculate available memory as free + buffers + cached
        local memFree=$(grep '^MemFree:' /proc/meminfo | awk '{print $2}')
        local memBuffers=$(grep '^Buffers:' /proc/meminfo | awk '{print $2}')
        local memCached=$(grep '^Cached:' /proc/meminfo | awk '{print $2}')
        
        if [ -n "$memFree" ] && [ -n "$memBuffers" ] && [ -n "$memCached" ]; then
            availableMemKB=$(( memFree + memBuffers + memCached ))
            #if [ "$silent" != "1" ]; then
                #echo "Detected ${availableMemKB}KB available memory (Free+Buffers+Cached)" >&2
            #fi
        fi
    fi
    
    # Also check virtual memory as a safety measure
    local commitLimit=$(grep '^CommitLimit:' /proc/meminfo | awk '{print $2}')
    local committedAS=$(grep '^Committed_AS:' /proc/meminfo | awk '{print $2}')
    
    if [ -n "$commitLimit" ] && [ -n "$committedAS" ] && [ $(( commitLimit - committedAS )) -gt 0 ]; then
        local availableVirtKB=$(( commitLimit - committedAS ))
        #if [ "$silent" != "1" ]; then
            #echo "Detected ${availableVirtKB}KB available virtual memory" >&2
        #fi
        
        # Use the more conservative of physical vs virtual available
        if [ $availableVirtKB -lt $availableMemKB ] || [ $availableMemKB -eq 0 ]; then
            availableMemKB=$availableVirtKB
        fi
    fi
}

# Detect available memory on macOS systems
function detectMacMemory() {
    availableMemKB=0
    
    # Get total physical memory
    if command -v sysctl >/dev/null 2>&1; then
        local totalMemBytes=$(sysctl -n hw.memsize 2>/dev/null)
        local totalMemKB=$(( totalMemBytes / 1024 ))
        
        # macOS doesn't easily report free memory, so estimate conservatively
        if [ $totalMemKB -gt 0 ]; then
            # Use 65% of total memory by default on macOS (conservative)
            availableMemKB=$(( totalMemKB * 65 / 100 ))
            if [ "$silent" != "1" ]; then
                echo "Detected ${totalMemKB}KB total memory on macOS, estimating ${availableMemKB}KB available" >&2
            fi
        fi
    fi
    
    # Try to get more accurate free memory using vm_stat if available
    if command -v vm_stat >/dev/null 2>&1; then
        # Extract page size and free pages
        local pageSize=$(vm_stat | grep "page size" | awk '{print $8}')
        if [ -z "$pageSize" ]; then
            # Try alternate format
            pageSize=$(vm_stat | grep "page size" | awk '{print $4}' | tr -d '.')
        fi
        
        local freePages=$(vm_stat | grep "Pages free" | awk '{print $3}' | tr -d '.')
        local inactivePages=$(vm_stat | grep "Pages inactive" | awk '{print $3}' | tr -d '.')
        
        if [ -n "$pageSize" ] && [ -n "$freePages" ] && [ -n "$inactivePages" ]; then
            # Convert pages to KB
            local freeMemKB=$(( (freePages + inactivePages) * pageSize / 1024 ))
            if [ "$silent" != "1" ]; then
                echo "Detected ${freeMemKB}KB free memory on macOS (vm_stat)" >&2
            fi
            
            # Use whichever is more conservative
            if [ $freeMemKB -lt $availableMemKB ] || [ $availableMemKB -eq 0 ]; then
                availableMemKB=$freeMemKB
            fi
        fi
    fi
}

# Check if this script is being sourced or run directly
# This is a more portable way than using BASH_SOURCE
if [ "$0" != "$BASH_SOURCE" ] && [ "$BASH_SOURCE" != "" ]; then
    # Being sourced
    :
else
    # Being run directly
    detectMemory "$@"
    echo "Detected memory: ${RAM}MB"
fi
