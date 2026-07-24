#!/bin/bash

usage(){
echo "
Written by Brian Bushnell, Noire
Last modified July 24, 2026

Description:  Manipulates neural-network training-vector files (the tab-delimited
'#dims'-header format read by train.sh / ml.DataLoader).  Combines any number of
inputs and writes one or more outputs, applying (in order): deduplicate, subsample,
shuffle, partition by fraction, and positive/negative balance.  The '#dims' header
is preserved as line 1 of each output (unlike a plain 'shuf', which corrupts it).

Usage:  vectorutils.sh in=<files> out=<file[:frac],...> [flags]

Example:
  vectorutils.sh a.tsv b.tsv c.tsv out=training.tsv:0.9,validation.tsv:0.1 shuffle samplerate=0.5 balance=0.3

Parameters:
in=<a,b,c>      Input vector files (also accepted as bare positional filenames).
out=<f1:0.9,f2:0.1>  Output files, each with a partition fraction.  Fractions are
                normalized to sum to 1; a single out= with no fraction gets the whole set.
shuffle         Randomly shuffle rows (needed for a meaningful random partition).
samplerate=1    Keep a random fraction of rows (subsample).
balance=0       Upsample the MINORITY class (label = last column >= 0.5) with random
                duplicates until it is this fraction of each output's total (e.g. 0.3).
deduplicate     (dedupe) Remove exact-duplicate rows (sort-based).
seed=1          RNG seed (>=0 deterministic; -1 random).
overwrite=t     (ow) Overwrite existing output files.

Java Parameters:
-Xmx            Set max memory; vectors are held in RAM, so size to the data
                (e.g. -Xmx200g).  Autodetects if unset.
-eoom           Exit on out-of-memory.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
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

	parseJavaArgs "--xmx=8g" "--xms=8g" "--percent=60" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ml.VectorUtils $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
