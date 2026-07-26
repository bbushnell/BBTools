#!/bin/bash

usage(){
echo "
Written by Noire

Description:  Fits probability-calibration constants [K,a,b,c] for the map
                  p = K * sigmoid(a*logit(x) + b)^c
              converting a raw score (e.g. a neural-net's continuous output)
              into a calibrated probability that the label is 1.  Fit per
              dataset, independently, by minimizing log-loss (Nelder-Mead).
              Generic: works on any net + vectors, or on raw (score,label)
              pairs, or either.  Prints log-loss and ECE before/after so the
              improvement is visible, and a '#cal K a b c' line ready to paste
              into a .bbnets file (clade.SerialNNLoader / CladeConfidence).

Usage:  calibrate.sh net=<net.bbnet> in=<vectors.tsv> [out=<cal.txt>]
        calibrate.sh in=<score_label_pairs.tsv> [out=<cal.txt>]

Vectors mode (net= given): 'in' is a '#dims IN OUT'-headed tab file; each row's
    first IN columns are run through the net, and column IN is the target label.
Pairs mode (no net): 'in' has a raw score and a label in two columns.

Parameters:
in=<file>       Tab-delimited vectors (with net=) or (score,label) pairs.
                Should be a HELD-OUT set, not the training data.
net=<file>      Optional single-output .bbnet; run it on the vectors in 'in'.
out=<file>      Write the '#cal K a b c' line here (also printed to stderr).
scorecol=0      Pairs mode: column holding the raw score.
labelcol=1      Pairs mode: column holding the label (>=0.5 -> 1).
bins=20         Bins for the ECE (expected calibration error) diagnostic.
restarts=5      Deterministic multi-start fits; best (lowest log-loss) is kept.
iters=2000      Max Nelder-Mead iterations per start.

Java Parameters:
-Xmx            Set Java memory, e.g. -Xmx4g.  Overrides autodetection.
-eoom           Exit on out-of-memory exception.
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

	parseJavaArgs "--xmx=4g" "--xms=4g" "--percent=42" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ml.Calibrate $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
