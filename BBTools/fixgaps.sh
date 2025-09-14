#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified September 11, 2019

Description:  Uses paired read insert sizes to estimate the correct
length of scaffold gaps, and resizes incorrectly-sized gaps.

Usage:  fixgaps.sh in=mapped.sam ref=scaffolds.fa out=fixed.fa

Standard parameters:
in=<file>       Reads mapped to the reference; should be sam or bam.
ref=<file>      Reference; may be fasta or fastq.
out=<file>      Modified reference; may be fasta or fastq.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Processing parameters:
gap=10          Consider any consecutive streak of Ns at least this long to
                be a scaffold break.  Gaps will not be resized to less than
                this.
border=0.4      Ignore the outermost (border*readlen) of an insert (read pair)
                when incrementing coverage.  A higher value is more accurate 
                but requires more coverage and/or longer inserts.  Range: 0-1.
mindepth=10     Minimum spanning read pairs to correct a gap.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

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

z="-Xmx4g"
z2="-Xms4g"
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
	freeRam 4000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

fixgaps() {
	local CMD="java $EA $EOOM $z -cp $CP consensus.FixScaffoldGaps $@"
	echo $CMD >&2
	eval $CMD
}

fixgaps "$@"
