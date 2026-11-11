#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 8, 2026

Description:  Generates multiple assemblies with Tadpole
to estimate the optimal kmer length.

Usage:
tadwrapper.sh in=reads.fq out=contigs%.fa k=31,62,93

Parameters:
out=<file>      Output file name.  Must contain a % symbol.
outfinal=<file> Optional destination for the selected assembly.
                Never replaces an input.  Existing outputs require ow=t.
                A failed publication retains the selected assembly.
k=31            Comma-delimited list of kmer lengths.
delete=f        Delete nonselected intermediates after success,
                excluding the final destination.
quitearly=f     Stop the initial K scan on a worse trial; ties continue.
bisect=f        Test K values between the current winner and its
                adjacent tested K values, refining where they improve.
expand=f        Recursively assemble with kmers shorter or longer
                than the current best until improvement halts.

Keep showstats=t: selection uses each trial's reported assembly statistics.
All other parameters are passed to Tadpole.  The search is heuristic.
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
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
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
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" assemble.TadpoleWrapper "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
