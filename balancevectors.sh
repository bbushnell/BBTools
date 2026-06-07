#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and UMP45
Last modified June 7, 2026

Description:  Balances a labeled training-vector TSV (from VcfToTrainingVectors)
into train/validation splits for neural-network training.  Keeps all positives and
performs stratified negative sampling so rare-but-hard false positives are not
drowned out: a small enriched fraction is drawn evenly across category axes
(variant type x depth/score/event-length, plus artifact axes - homopolymer,
allele fraction, strand ratio, mapq - where false positives concentrate), and the
remainder is a representative random sample.  Deterministic given a seed.

Usage:  balancevectors.sh in=<vectors.tsv> outtrain=<file> outval=<file>

I/O parameters:
in=<file>       Input labeled vector TSV (last column is the 0-1 label).
outtrain=<file> Output training split.
outval=<file>   Output validation split.

Balancing parameters:
posfraction=0.3 Target fraction of positives in the output (0.3 = 30% pos / 70% neg).
valfraction=0.1 Fraction of the balanced set held out for validation.
enrich=t        Use stratified category sampling for negatives (f = purely random).
newcats=t       Include the artifact category axes (homopolymer/allele-fraction/
                strand-ratio/mapq) in addition to the type/depth/score/length axes.
newweight=0.3   Sample the artifact axes at this fraction of the base per-category
                quota (a smaller fraction, not an equal split).
quota=          Override the per-existing-category sample size (default: derived).
noscore=t       Zero the composite-score column (29) on output, so the network must
                learn its own scoring instead of copying the existing score.
seed=42         Random seed (sampling and shuffle are deterministic given this).

Java Parameters:
-Xmx            Set Java memory, e.g. -Xmx64g.  Large inputs (tens of millions of
                vectors) need enough heap to hold all lines; budget ~1.5x the file size.
-eoom           Exit on out-of-memory.
-da             Disable assertions.

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

	parseJavaArgs "--xmx=4g" "--xms=4g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP var2.BalanceVectors $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
