#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified September 6, 2023

Description:  Splits a sequence file evenly into multiple files.

Usage:  partition.sh in=<file> in2=<file2> out=<outfile> out2=<outfile2> ways=<number>

in2 and out2 are for paired reads and are optional.
If input is paired and out2 is not specified, data will be written interleaved.
Output filenames MUST contain a '%' symbol.  This will be replaced by a number.

Parameters and their defaults:

in=<file>       Input file.
out=<file>      Output file pattern (containing a % symbol, like 'part%.fa').
ways=-1         The number of output files to create; must be positive.
pacbio=f        Set to true to keep PacBio subreads together.
bp=f            Optimize for an even split by base pairs instead of sequences.
                Not compatible with PacBio mode.

ow=f            (overwrite) Overwrites files that already exist.
app=f           (append) Append to files that already exist.
zl=4            (ziplevel) Set compression level, 1 (low) to 9 (max).
int=f           (interleaved) Determines whether INPUT file is considered interleaved.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
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

z="-Xmx1g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
}
calcXmx "$@"

function partition() {
	local CMD="java $EA $EOOM $z -cp $CP jgi.PartitionReads $@"
	echo $CMD >&2
	eval $CMD
}

partition "$@"
