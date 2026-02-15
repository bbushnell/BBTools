#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 18, 2016

Description:  Generates multiple assemblies with Tadpole
to estimate the optimal kmer length.

Usage:
tadwrapper.sh in=reads.fq out=contigs%.fa k=31,62,93

Parameters:
out=<file>      Output file name.  Must contain a % symbol.
outfinal=<file> Optional.  If set, the best assembly file
                will be renamed to this.
k=31            Comma-delimited list of kmer lengths.
delete=f        Delete assemblies before terminating.
quitearly=f     Quit once metrics stop improving with longer kmers.
bisect=f        Recursively assemble with the kmer midway between
                the two best kmers until improvement halts.
expand=f        Recursively assemble with kmers shorter or longer
                than the current best until improvement halts.

All other parameters are passed to Tadpole.
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

	parseJavaArgs "--xmx=15000m" "--xms=15000m" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP assemble.TadpoleWrapper $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
