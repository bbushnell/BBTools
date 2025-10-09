#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 9, 2025

Description:  Visualizes 3D compositional metrics (GC, HH, CAGA) as 2D scatter plots.
Supports both TSV interval data and FASTA input (via ScalarIntervals).
Generates PNG images with configurable scaling and point sizes.

Usage:  cloudplot.sh in=<input file> out=<output file>
e.g.
cloudplot.sh in=data.tsv out=plot.png
or
cloudplot.sh in=ecoli.fasta out=plot.png shred=5k

Standard parameters:
in=<file>       Primary input; TSV (GC/HH/CAGA columns) or FASTA/FASTQ.
out=<file>      Output PNG image file.

Rendering parameters:
scale=1         Image scale multiplier (1=800x600, 2=1600x1200, etc).
pointsize=2     Radius of plotted points in pixels.
xmin=-1         X-axis minimum (GC). Negative = autoscale.
xmax=-1         X-axis maximum (GC). Negative = autoscale.
ymin=-1         Y-axis minimum (HH). Negative = autoscale.
ymax=-1         Y-axis maximum (HH). Negative = autoscale.
zmin=-1         Z-axis/color minimum (CAGA). Negative = autoscale.
zmax=-1         Z-axis/color maximum (CAGA). Negative = autoscale.

FASTA processing parameters (only used with FASTA input):
window=0        If nonzero, calculate metrics over sliding windows.
                Otherwise calculate per contig.
interval=5000   Generate a data point every this many bp.
shred=-1        If positive, set window and interval to the same size.
break=t         Reset metrics at contig boundaries.
minlen=500      Minimum interval length to generate a point.
maxreads=-1     Maximum number of reads/contigs to process.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

#This block allows symlinked shellscripts to correctly set classpath.
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

calcXmx () {
    # Source the new scripts
    source "$DIR""/memdetect.sh"
    source "$DIR""/javasetup.sh"

    parseJavaArgs "--mem=2g" "--mode=auto" "$@"

    # Set environment paths
    setEnvironment
}
calcXmx "$@"

cloudplot() {
	if [[ $# -eq 0 ]]; then
		usage
		return
	fi
	local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP scalar.CloudPlot $@"
	#echo $CMD >&2
	eval $CMD
}

cloudplot "$@"
