#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 13, 2025

Description:  Calculates EST (expressed sequence tags) capture by an assembly from a sam file.
Designed to use BBMap output generated with these flags:
k=13 maxindel=100000 customtag ordered nodisk

Usage:          bbest.sh in=<sam file> out=<stats file>

Parameters:
in=<file>       Specify a sam file (or stdin) containing mapped ests.
                If a fastq file is specified it will be mapped to a temporary
                sam file using BBMap, then deleted.
out=<file>      Specify the output stats file (default is stdout).
ref=<file>      Specify the reference file (optional).
est=<file>      Specify the est fasta file (optional).
fraction=0.98   Min fraction of bases mapped to ref to be 
                considered 'all mapped'.

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

	parseJavaArgs "--xmx=120m" "--xms=120m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.SamToEst $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
