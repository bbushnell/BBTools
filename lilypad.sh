#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified September 13, 2019

Description:  Uses mapped paired reads to generate scaffolds from contigs.
Designed for use with ordinary paired-end Illumina libraries.

Usage:  lilypad.sh in=mapped.sam ref=contigs.fa out=scaffolds.fa

Standard parameters:
in=<file>       Reads mapped to the reference; should be sam or bam.
ref=<file>      Reference; may be fasta or fastq.
out=<file>      Modified reference; should be fasta.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Processing parameters:
gap=10          Pad gaps with a minimum of this many Ns.
mindepth=4      Minimum spanning read pairs to join contigs.
maxinsert=3000  Maximum allowed insert size for proper pairs.
mincontig=200   Ignore contigs under this length if there is a 
                longer alternative.
minwr=0.8       (minWeightRatio) Minimum fraction of outgoing edges 
                pointing to the same contig.  Lower values will increase
                continuity at a risk of misassemblies.
minsr=0.8       (minStrandRatio) Minimum fraction of outgoing edges 
                indicating the same orientation.  Lower values will increase
                continuity at a possible risk of inversions.
passes=8        More passes may increase continuity.
samestrandpairs=f  Read pairs map to the same strand.  Currently untested.

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

lilypad() {
	local CMD="java $EA $EOOM $z -cp $CP consensus.Lilypad $@"
	echo $CMD >&2
	eval $CMD
}

lilypad "$@"
