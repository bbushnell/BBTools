#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified December 19, 2018

Description:  Calls peaks from a whitespace-delimited count histogram.
              One-column rows contain counts at successive indices; indexed rows
              contain an x-coordinate and count. An optional final GC-percent
              column must contain a decimal point (for example, 50.0).

Usage:        callpeaks.sh in=<histogram file> out=<output file>

Peak-calling parameters:
in=<file>           'in=stdin.fq' will pipe from standard in.
out=<file>          Write the peaks to this file.  Default is stdout.
minHeight=2         (h) Ignore peaks shorter than this.
minVolume=5         (v) Ignore peaks with less area than this.
minWidth=3          (w) Minimum width in histogram positions, not bases.
minPeak=2           (minp) Ignore peaks with an X-value below this.
                    Useful when low-count kmers are filtered).
maxPeak=BIG         (maxp) Ignore peaks with an X-value above this.
maxPeakCount=10     (maxpc) Target count when condensing peaks by height and volume.
                    Separate cutoffs and ties can retain more peaks.
countColumn=1       (col) For multi-column input, this column, zero-based,
                    contains the counts.
ploidy=-1           Specify ploidy; otherwise it will be autodetected.
logscale=f          Transform to log-scale prior to peak-calling.  Useful
                    for kmer-frequency histograms.

Smoothing parameters:
smoothradius=0      Integer radius of triangle filter.  Set above zero to
                    smooth data prior to peak-calling.  Higher values are
                    smoother.
smoothprogressive=t Set to true to widen the filter as the x-coordinate
                    increases.  Useful for kmer-frequency histograms.
maxradius=10        Stop further radius growth after crossing this value;
                    the initial smoothradius is not clamped.
progressivemult=2   Spacing factor for positions where the radius increases by one.

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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.CallPeaks $@"
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" jgi.CallPeaks "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
