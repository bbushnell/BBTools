#!/bin/bash

usage(){
echo "
Written by UMP45
Last modified June 13, 2026

Description:  Grades a VCF against a truth set and produces concordance metrics
and cumulative histograms.  Replaces a multi-step filtervcf+comparevcf pipeline
with a single fast tool.

Usage:  gradevcf.sh in=<file> truth=<file>

I/O parameters:
in=<file>       Input candidate VCF (from CallVariants or filtervcf clearfilters).
truth=<file>    Truth VCF (e.g. GIAB, normalized).
out=<file>      Optional: write passing variants to this file.
ref=<file>      Reference fasta (optional; used for NN vector or normalize).
overwrite=f     Set to false to force the program to abort rather than
                overwrite an existing file.

Region parameters:
bed=<file>      Restrict evaluation to variants inside these BED intervals
                (e.g. a high-confidence benchmark region set).

Quality cutoff parameters:
minscore=0.0    Reject variants with QUAL below this threshold.

Neural network parameters:
net=<file>      Neural network file (.bbnet).
netmode=        Feature vector mode (ump45, elba, lawrence, donovan).
netcutoff=      NN output threshold; variants below this are rejected.
                Use 'auto' to read cutoff from the .bbnet file.
cutoff=         Alias for netcutoff=.
includescore=f  Include composite score in NN input vector.

Histogram output:
histnn=<file>   Write the NN score histogram (raw per-bin counts, then cumulative
                threshold/TP/FP/FN) at 0.002 resolution over 0-2.0.
histqual=<file> Write the QUAL histogram (raw per-bin counts, then cumulative
                threshold/TP/FP/FN) at 0.1 resolution over 0-60.

Operating-point parameters:
crossover=0     If >0, also report the cutoff at the FN=crossover*FP operating
                point (e.g. crossover=4 -> the FN=4*FP point) for BOTH the QUAL and
                NN score axes, plus FP, FN, and FP+FN there.  Emitted as
                CROSSOVER_QUAL / CROSSOVER_NN lines.  0 (default) disables.

Variant-quality pre-filtering parameters:
minreads=0              Ignore variants seen in fewer reads.
minqualitymax=0         Ignore variants with lower max base quality.
minedistmax=0           Ignore variants with lower max distance from read ends.
minmapqmax=0            Ignore variants with lower max mapq.
minidmax=0              Ignore variants with lower max read identity.
minscore=0.0            Ignore variants with lower Phred-scaled score.
clearfilters            Reset all variant filters to zero.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP var2.GradeVCF $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
