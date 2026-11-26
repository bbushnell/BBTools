#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 22, 2024

Description:  Filters a sam file to retain reads passing clipping bounds and
at least one of three alternate retention routes: substitutions, counted indels,
or no substitutions and no counted indels.  The substitution route uses all
substitutions for minsubs and requires at least one substitution in the inclusive
minq/maxq quality interval unless minsubs is below 1.  Optional counted-indel
and perfect-read routes are controlled by their flags and are not gated by
minsubs, minq, or maxq.  Used for manually examining specific reads that may
have incorrectly calibrated quality scores.

Usage:  filtersubs.sh in=<file> out=<file> minq=<number> maxq=<number>

Parameters:
in=<file>       Input sam or bam file.
out=<file>      Output file.
minq=0          Substitution route: inclusive minimum for one qualifying sub.
maxq=99         Substitution route: inclusive maximum for one qualifying sub.
countindels=t   Optional route: keep reads with counted indels; not quality-gated.
minsubs=1       Substitution route: require this many total substitutions.
minclips=0      Require at least this many clipped bases.
maxclips=-1     If nonnegative, discard reads with more clipped bases.
keepperfect=f   Optional route: keep reads with no substitutions or counted indels.

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
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=120m" "--xms=120m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.FilterReadsWithSubs $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" jgi.FilterReadsWithSubs "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
