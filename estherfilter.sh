#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 21, 2015

Description:  BLASTs queries against reference, and filters out hits with
              scores less than 'cutoff'.  The score is taken from column 12
              of the BLAST output.  The specific BLAST command is:
              blastall -p blastn -i QUERY -d REFERENCE -e 0.00001 -m 8

Usage:  estherfilter.sh <query> <reference> <cutoff>

For example:

estherfilter.sh reads.fasta genes.fasta 1000 > results.txt

'fasta' can be used as a fourth argument to get output in Fasta format.  Requires more memory.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
"
}

if [ -z "$1" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	usage
	exit
fi

resolveSymlinks(){
	SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
	CP="$DIR/current/"
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=3200m" "--xms=3200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP driver.EstherFilter $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
