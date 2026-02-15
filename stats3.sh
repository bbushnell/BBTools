#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 21, 2025

Description:  In progress.
Generates some assembly stats for multiple files.

Usage:        stats3.sh in=file
Or:           stats3.sh in=file,file
Or:           stats3.sh file file file

Parameters:
in=file         Specify the input fasta file(s), or stdin.
                Multiple files can be listed without a 'in=' flag.
out=stdout      Destination of primary output; may be directed to a file.

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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.AssemblyStats3 $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
