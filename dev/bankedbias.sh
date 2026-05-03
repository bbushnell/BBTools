#!/bin/bash

usage(){
echo "
Written by Chloe
Last modified April 2026

Description:  Multi-bucket per-state bias simulator for BDLL5 (banked 3-bit
exponent + 2-bit history, 6 buckets per 32-bit word).

Unlike mantissacompare.sh (single-bucket, no banking), this runs a full BDLL5
tracker (default 1536 buckets = 256 words = 1 KB, the standard BDLL5 size).
Bank promotion is a word-level operation that correlates 6 buckets, and
changes per-state collision dynamics vs plain UDLL6.

Outputs:
  stderr: per-tier verbose report (matches mantissacompare.sh format).
  stdout: HSB table TSV ready for StateTable.loadHsbTable().

Usage:  bankedbias.sh buckets=1536 inner=4m outer=1024 t=32

Parameters:
buckets=1536    Bucket count. Real modBuckets = roundToWords(b)*6.
                Standard BDLL5 size is 1536 (256 words = 1 KB).
inner=4000000   Elements per trial.
outer=1024      Number of independent trials.
maxtier=11      Highest NLZ tier to record.
sstier=11       Tier to use for steady-state row (default 11 matches
                MantissaCompare2 convention).
pf=0.004        PROMOTE_FRAC for global tier advancement.
t=1             Threads.
avg=geo         Averaging mode for HSB table output (lin, geo, harm).
                geo matched the UDLL6 empirical winner.

Workflow:
  1. Run bankedbias.sh to generate HSB TSV on stdout.
  2. Pass the TSV to ddlcalibrate.sh via hsbtable2=file.tsv.
  3. Run ddlcalibrate.sh cf=f to generate a BDLL5 v5 CF table with new biases.
  4. Run ddlcalibrate.sh cf=t to verify Mean+H beats Mean.

Java Parameters:
-Xmx            Heap size, e.g. -Xmx8g.
-eoom           Exit if out of memory.
-da             Disable assertions.

Contact: bbushnell@lbl.gov
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

	parseJavaArgs "--xmx=2g" "--xms=200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP cardinality.BankedBiasSimulator $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
