#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 12, 2025

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
order=caga,hh,gc  Plotting order of dimensions as x,y,z.
scale=1         Image scale multiplier (1=1024x768).
pointsize=3.5   Width of plotted points in pixels.
autoscale=t     Autoscale dimensions with negative values based on data.
                If false, they will be scaled to 0-1.
xmin=-1         X-axis minimum.
xmax=-1         X-axis maximum.
ymin=-1         Y-axis minimum.
ymax=-1         Y-axis maximum.
zmin=-1         Z-axis (rotation/color) minimum.
zmax=-1         Z-axis (rotation/color) maximum.
xpct=0.998      Percentile of x-axis values to use for autoscaling.
ypct=0.998      Percentile of y-axis values to use for autoscaling.
zpct=0.99       Percentile of z-axis values to use for autoscaling.

Taxonomy/Coloring parameters:
colorbytax=f    Color by taxonomy.  Default coloring is by the 
colorbyname=f   Color by contig name, so points on the same contig have
                the same, random color.
level=          Raise taxonomy to this level before assigning color.
                Requires a taxonomic tree.  e.g. 'level=genus'
                See https://sourceforge.net/projects/bbmap/files/Resources/
parsetid=f      Parse TaxIDs from file and sequence headers.
sketch=f        Use BBSketch (SendSketch) to assign taxonomy per contig.
clade=f         Use QuickClade to assign taxonomy per contig.

Decorrelation parameters:
decorrelate=t   Modify plotted data to reduce inter-dimension correlation.
GChh=-0.5       Correlation between GC and HH.
GChhs=0.2       (GChhStrength) Modify HH by -GChhs*GC*GChh.
hhGCs=1.4       (hhGCStrength) Modify GC by -hhGCs*hh*GChh.
GCcaga=0.1      Correlation between GC and CAGA.
GCcagas=0.5     (GCcagaStrength) Modify CAGA by -GCcagas*GC*GCcaga.
cagaGCs=0.0     (cagaGCStrength) Modify GC by -cagaGCs*caga*GCcaga.

Sequence processing parameters (not used with TSV input):
window=50000    If nonzero, calculate metrics over sliding windows.
                Otherwise calculate per contig.
interval=10000  Generate a data point every this many bp.
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
