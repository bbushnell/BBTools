#!/bin/bash

usage(){
echo "
Written by Eru
Last modified August 10, 2026

Description:  Computes per-organism dinucleotide-composition features (HH and
CAGA, from tracker.KmerTracker) from a contig FASTA whose headers embed the
taxon id (..._tid_<int>), then min-max scales each feature so 1 = the maximum
and 0 = the minimum observed across the TRAINING organisms only (validation
organisms are scaled by the same train-derived range, clamped to [0,1]).  The
train/val split reproduces MagQCVectorMaker's split exactly (same seed and
valfrac => same held-out organisms).  Output is consumed by
magqcvectormaker.sh subnet=ncrna snhhcaga=t kmerfile=<output>.

This is a DEVELOPMENT tool (training-feature generation); it is not part of
the end-user MAG-QC path.

Usage:  genomedimerfeatures.sh in=<renamed.fa> cache=<percontig_cache.tsv> \\
          out=<kmerfeat.tsv> [seed=1] [valfrac=0.10]

Required parameters:
in=<file>       Contig FASTA with _tid_<int> header suffixes.
cache=<file>    Per-contig precompute cache TSV (defines the usable organisms).
out=<file>      Output tid<tab>HH<tab>CAGA TSV (train-scaled).

Optional parameters (and their defaults):
seed=1          RNG seed; MUST match the MagQCVectorMaker run it feeds.
valfrac=0.10    Held-out organism fraction; MUST match likewise.

Java Parameters:
-Xmx            Set Java heap (overrides autodetection).
-eoom           Exit if an out-of-memory exception occurs.
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
	DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=8g" "--xms=1g" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prot.GenomeDimerFeatures $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
