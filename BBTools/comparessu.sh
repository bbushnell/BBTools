#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified December 4, 2019

Description:  Aligns SSUs to each other and reports identity.
This requires sequences annotated with a taxID in their header.

Usage:  comparessu.sh in=<input file> out=<output file>

Input may be fasta or fastq, compressed or uncompressed.

Standard parameters:
in=<file>       Input sequences.
out=<file>      Output data.
t=              Set the number of threads; default is logical processors.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
showspeed=t     (ss) Set to 'f' to suppress display of processing speed.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression
                level; lower compression is faster.
reads=-1        If positive, quit after this many sequences.

Processing parameters:
ata=f           Do an all-to-all comparison.  Otherwise, each sequence will
                only be compared to one other randomly-selected sequence
                per taxonomic level.
minlen=0        Ignore sequences shorter than this.
maxlen=BIG      Ignore sequences longer than this.
maxns=-1        If positive, ignore sequences with more than this many Ns.


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

comparessu() {
	local CMD="java $EA $EOOM $z -cp $CP sketch.CompareSSU $@"
	echo $CMD >&2
	eval $CMD
}

comparessu "$@"
