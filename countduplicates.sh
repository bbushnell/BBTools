#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 14, 2022

Description:  Counts duplicate sequences probabilistically,
using around 20 bytes per unique read.  Read pairs are treated
as a single read.  Reads are converted to a hashcode and only
the hashcode is stored when tracking duplicates, so (rare) hash
collisions will result in false positive duplicate detection.
Optionally outputs the deduplicated and/or duplicate reads.

Usage:  countduplicates.sh in=<input file>

Input may be fasta, fastq, or sam, compressed or uncompressed.
in2, out2, and outd2 are accepted for paired files.

Standard parameters:
in=<file>       Primary input, or read 1 input.
out=<file>      Optional output for deduplicated reads.
outd=<file>     Optional output for duplicate reads.  An extension like .fq
                will output reads; .txt will output headers only.
stats=stdout    May be replaced by a filename to write stats to a file.
showspeed=t     (ss) Set to 'f' to suppress display of processing speed.

Processing parameters (these are NOT mutually exclusive):
bases=t         Include bases when generating hashcodes. 
names=f         Include names (headers) when generating hashcodes.
qualities=f     Include qualities when generating hashcodes.
maxfraction=-1.0  Set to a positive number 0-1 to FAIL input
                  that exceeds this fraction of reads with duplicates.
maxrate=-1.0    Set to a positive number >=1 to FAIL input that exceeds this
                average duplication rate (the number of copies per read).
failcode=0      Set to some other number like 1 to produce a
                non-zero exit code for failed input.
samplerate=1.0  Fraction of reads to subsample, to conserve memory.  Sampling
                is deterministic - if a read is sampled, copies will be too.
                Unsampled reads are not sent to any output stream or counted 
                in statistics.

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

countDupes() {
	local CMD="java $EA $EOOM $z -cp $CP jgi.CountDuplicates $@"
	echo $CMD >&2
	eval $CMD
}

countDupes "$@"
