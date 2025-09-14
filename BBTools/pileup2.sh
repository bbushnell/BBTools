#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified April 17, 2024

Description:  This is a version of Pileup designed to process multiple files
concurrently.  If you only have one input file just use regular Pileup.
Other than those mentioned here, the flags are the same as in pileup.sh.

Usage:        pileup2.sh in=<file,file,file> out=<file>
Alternate:    pileup2.sh *.sam out=<file>

Parameters:
in=<file,file>     The input sam/bam files.  Omit the 'in=' if a wildcard
                   is used.
streams=-1         If positive, use at most this many concurrent streams.
                   Default is half the number of logical processors.  Note 
                   that each stream uses multiple threads.
atomic=false       Use atomic arrays instead of locks.
prealloc=false     Preallocate coverage arrays instead of creating them 
                   as needed.

Java Parameters:
-Xmx               This will set Java's memory usage, overriding 
                   autodetection. -Xmx20g will 
                   specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  
                   The max is typically 85% of physical memory.
-eoom              This flag will cause the process to exit if an out-of-memory
                   exception occurs.  Requires Java 8u92+.
-da                Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
"
}

#This block allows symlinked shellscripts to correctly set classpath.
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx1g"
z2="-Xms1g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

pileup() {
	local CMD="java $EA $EOOM $z -cp $CP jgi.CoveragePileupMT $@"
	echo $CMD >&2
	eval $CMD
}

pileup "$@"
