#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 2, 2017

Description:  Remove Smart Bell adapters from PacBio reads.

Usage:        removesmartbell in=<input> out=<output> split=t

Input may be fasta or fastq, compressed or uncompressed (not H5 files).

Parameters:
in=file         Specify the input file, or stdin.
out=file        Specify the output file, or stdout.
adapter=        Specify the adapter sequence (default is normal SmrtBell).
split=t            t: Splits reads at adapters.
                   f: Masks adapters with X symbols.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
"
}

if [ -z "$1" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	usage
	exit
fi

resolveSymlinks(){
	SCRIPT="$0"
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

	parseJavaArgs "--xmx=400m" "--xms=400m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP pacbio.RemoveAdapters2 $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
