#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified April 25, 2019

Description:  Compresses or decompresses files based on extensions.
This only exists because the syntax and default behavior of many
compression utilities is unintuitive; it is just a wrapper, and
relies on existing executables in the command line (pigz, lbzip, etc.)
Does not delete the input file.
Does not untar files.

Usage:  unzip.sh in=<file> out=<file>

Parameters:
in=<file>       Input file.
out=<file>      Output file for good reads.
zl=             Set the compression level; 0-9 or 11.

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

	parseJavaArgs "--xmx=80m" "--xms=80m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.Unzip $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
