#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 5, 2022

Description:  Parses verbose output from bloomfilter.sh for a specific paper.
Irrelevant for most people, but useful for reproducing published results.
You use it to parse output from bloomfilter.sh and tabulate it.

Usage:  bloomfilterparser.sh in=<input file> out=<output file>

...where the input file is whatever bloomfilter.sh prints to the screen.  E.G.
in=slurm-3249652.out out=summary.txt

You get details of calls to increment() if you add the verbose flag.

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

	parseJavaArgs "--xmx=300m" "--xms=300m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP bloom.ParseBloomFilter $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
