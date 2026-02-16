#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 15, 2026

Description:  Visualizes up to 5D compositional metrics as 2D scatter plots.
X, Y, Z (rotation), Size, and Color channels can display GC, HH, CAGA, Depth, or Length.
Supports both TSV interval data and FASTA input (via ScalarIntervals).
Generates PNG images with configurable scaling, sizes, and colors.

Usage:  cloudplot.sh in=<input file> out=<output file> [options]

Examples:
# Basic 3D plot (legacy format)
cloudplot.sh in=contigs.fa out=plot.png order=gc,hh,caga

# 5D plot: GC vs HH, rotated by CAGA, sized by depth, colored by taxonomy
cloudplot.sh in=contigs.fa out=plot.png order=gc,hh,caga,depth colorby=tax cov=coverage.txt

# Alternative: depth on X-axis, length as size
cloudplot.sh in=contigs.fa out=plot.png order=depth,hh,caga,length colorby=tax cov=coverage.txt

# Fixed size with taxonomic coloring
cloudplot.sh in=contigs.fa out=plot.png order=gc,hh,caga colorby=taxonomy

Standard parameters:
in=<file>       Primary input; TSV (GC/HH/CAGA columns) or FASTA/FASTQ.
out=<file>      Output PNG image file.

Dimension assignment:
order=gc,hh,caga        Assign metrics to X, Y, Z(rotation) dimensions (3D mode).
order=gc,hh,caga,depth  Assign metrics to X, Y, Z(rotation), Size dimensions (5D mode).
                        Available metrics: gc, hh, caga, depth, length, taxonomy, none
                        Note: Taxonomy can ONLY be used for Z (rotation) or colorby.

colorby=<metric>        Metric for point color (default: caga gradient).
                        Options: gc, hh, caga, depth, length, taxonomy
                        Use 'taxonomy' or 'tax' for taxonomic coloring.

Depth/coverage sources:
cov=<file>      Coverage file from pileup.sh (format: #ID, Avg_fold) or
                covmaker.sh (format: #Contigs, AvgFold).
depth=<file>    SAM/BAM file for depth calculation (not yet implemented -
                use pileup.sh to generate coverage file instead).

Rendering parameters:
scale=1         Image scale multiplier (1=1024x768).
pointsize=3.5   Base point size in pixels.
                When size dimension is set, this is the reference size.
                When size dimension is NOT set, this is the fixed size.

Size scaling:
minsize=-1      Minimum point size for variable sizing (pixels).
                Default: 0.8 * pointsize (e.g., 2.8 pixels if pointsize=3.5).
                Negative value triggers autoscaling.
maxsize=-1      Maximum point size for variable sizing (pixels).
                Default: 3.0 * pointsize (e.g., 10.5 pixels if pointsize=3.5).
                Negative value triggers autoscaling.
spct=0.998      Percentile of size values to use for autoscaling.
                Note: Depth and Length use logarithmic scaling for size.

Axis scaling:
autoscale=t     Autoscale dimensions with negative min/max based on data percentiles.
                If false, dimensions are scaled to 0-1 range.
xmin=-1         X-axis minimum (negative = autoscale from data).
xmax=-1         X-axis maximum (negative = autoscale from data).
ymin=-1         Y-axis minimum (negative = autoscale from data).
ymax=-1         Y-axis maximum (negative = autoscale from data).
zmin=-1         Z-axis (rotation) minimum (negative = autoscale from data).
zmax=-1         Z-axis (rotation) maximum (negative = autoscale from data).
smin=-1         Size minimum (negative = autoscale from data).
smax=-1         Size maximum (negative = autoscale from data).
xpct=0.998      Percentile of x-axis values to use for autoscaling.
ypct=0.998      Percentile of y-axis values to use for autoscaling.
zpct=0.99       Percentile of z-axis values to use for autoscaling.

Taxonomy/Coloring parameters:
colorbytax=f    (Legacy) Color by taxonomy. Use colorby=tax instead.
colorbyname=f   (Legacy) Color by contig name. Not compatible with colorby parameter.
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

Dimension Usage Guidelines:
- X, Y axes: Best for GC, HH, CAGA (0-1 range, easy to interpret)
- Z (rotation): Works for any metric, but CAGA or Taxonomy recommended
- Size: Works well for Depth or Length (high dynamic range with log scaling)
- Color: Taxonomy (categorical) or any continuous metric (gradient)

Validation Rules:
- Taxonomy can ONLY be assigned to Z (rotation) or colorby
- Attempting to assign Taxonomy to X, Y, or Size will produce an error
- All other metrics can be assigned to any dimension

Notes:
- When size dimension is enabled, point length no longer varies with Y position
- Depth and Length metrics use logarithmic scaling for size dimension
- GC, HH, CAGA use linear scaling for size dimension
- Color gradient uses the cagaToColor6 palette (Red → Purple → Blue → Cyan → Green → Yellow)

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
	CP="$DIR/current/"
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=2g" "--xms=2g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	# Detect if X11 is available for antialiased rendering
	if command -v xset >/dev/null 2>&1 && xset q &>/dev/null; then
		HEADLESS=""  # X11 available - use normal rendering with antialiasing
	else
		HEADLESS="-Djava.awt.headless=true"  # No X11 - use headless mode
	fi

	CMD="java $HEADLESS $EA $EOOM $SIMD $XMX $XMS -cp $CP scalar.CloudPlot $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
