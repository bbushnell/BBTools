#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 19, 2026

Description:  Calculates compositional scalar metrics from nucleotide sequences.
Computes GC content, HH (homo/hetero dimer ratio), CAGA (transition preference),
depth, and length for sequence intervals. Outputs data in TSV format for 
visualization with CloudPlot or external analysis tools.

Output TSV Format:
#Name   Length  GC      HH      CAGA    Depth   Start   TaxID   TaxID2
- Name: Sequence/contig name
- Length: Interval length in bases
- GC: GC content (0-1)
- HH: Homopolymer-heteropolymer ratio (0-1, GC-independent)
- CAGA: Compositional asymmetry (0-1, GC-independent)
- Depth: Read coverage depth (from cov/depth file or header)
- Start: Start position within contig (0 for whole contigs)
- TaxID: Primary taxonomy ID (from clade, sketch, or header)
- TaxID2: Secondary taxonomy ID (for concordance checking)

Usage:  scalarintervals.sh in=<input file> out=<output file> [options]

Examples:
# Basic interval generation
scalarintervals.sh in=assembly.fa out=data.tsv shred=20k header=t

# With coverage from BAM file
scalarintervals.sh in=contigs.fa out=data.tsv shred=20k depth=mapped.bam header=t

# With coverage file and taxonomy
scalarintervals.sh in=assembly.fa out=data.tsv shred=20k cov=coverage.txt clade=t header=t

# Multiple input files
scalarintervals.sh *.fa.gz out=combined.tsv shred=10k header=t printname=t

Standard parameters:
in=<file>       Primary input; FASTA or FASTQ.
                Can be a directory or comma-delimited list.
                Filenames can also be used without in=
out=stdout      Output TSV file. Mean and stdev printed to stderr.

Depth/Coverage parameters:
cov=<file>      Coverage file from pileup.sh (format: #ID, Avg_fold) or
                covmaker.sh (format: #Contigs, AvgFold).
depth=<file>    SAM/BAM file for depth calculation.
                Calculates depth from aligned bases in the file.

Processing parameters:
header=f        Print TSV header line.
window=50000    If nonzero, calculate metrics over sliding windows.
                Otherwise calculate per contig. Larger has lower variance.
interval=10000  Generate a data point every this many bp.
shred=-1        If positive, set window and interval to the same size.
                Example: shred=20k sets both window and interval to 20000.
break=t         Reset metrics at contig boundaries.
minlen=500      Minimum interval length to generate a point.
maxreads=-1     Maximum number of reads/contigs to process.
printname=f     Print contig names in output.
printpos=f      Print start position in output (same as Start column).
printtime=t     Print timing information to stderr.

Taxonomy parameters:
parsetid=f      Parse TaxIDs from file and sequence headers.
sketch=f        Use BBSketch (SendSketch) to assign taxonomy per contig.
                Assigns TaxID2 field.
clade=f         Use QuickClade to assign taxonomy per contig.
                Assigns TaxID field.

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

	parseJavaArgs "--xmx=800m" "--xms=800m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP scalar.ScalarIntervals $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
