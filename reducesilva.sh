#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified July 31, 2015

Description:  Reduces Silva entries down to one entry per taxa (the first).
This is accomplished by splitting the semicolon-delimited name on semicolons,
and assuming everything is in the form of:
kingdom;phylum;class;order;family;genus;species
...so it's not very reliable.

Usage:  reducesilva.sh in=<file> out=<file> column=<1>

Parameters:
column          The taxonomic level.  0=species, 1=genus, etc.
ow=f            (overwrite) Overwrites files that already exist.
zl=4            (ziplevel) Set compression level, 1 (low) to 9 (max).
fastawrap=70    Length of lines in fasta output.

Sampling parameters:
reads=-1        Set to a positive number to only process this many INPUT sequences, then quit.

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

function reducesilva() {
	local CMD="java $EA $EOOM $z -cp $CP driver.ReduceSilva $@"
	echo $CMD >&2
	eval $CMD
}

reducesilva "$@"
