#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 13, 2025

Description:  Calculates some scalars from nucleotide sequence data.
Writes them periodically as a tsv.

Usage:  scalarintervals.sh in=<input file> out=<output file>
e.g.
scalarintervals.sh in=ecoli.fasta out=data.tsv shred=5k
or
scalarintervals.sh *.fa.gz out=data.tsv shred=5k

Standard parameters:
in=<file>       Primary input; fasta or fastq.
                This can also be a directory or comma-delimited list.
		Filenames can also be used without in=
out=stdout      Set to a file to redirect tsv output.  The mean and stdev
                will be printed to stderr.

Processing parameters:
header=f        Print a header line.
window=50000    If nonzero, calculate metrics over sliding windows.
                Otherwise calculate per contig.  Larger has lower variance.
interval=10000  Generate a data point every this many bp.
shred=-1        If positive, set window and interval to the same size.
break=t         Reset metrics at contig boundaries.
minlen=500      Minimum interval length to generate a point.
maxreads=-1     Maximum number of reads/contigs to process.
printname=f     Print contig names in output.
printpos=f      Print contig position in output.
printtime=t     Print timing information to screen.
parsetid=f      Parse TaxIDs from file and sequence headers.
sketch=f        Use BBSketch (SendSketch) to assign taxonomy per contig.
clade=f         Use QuickClade to assign taxonomy per contig.

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

    parseJavaArgs "--mem=800m" "--mode=fixed" "$@"

    # Set environment paths
    setEnvironment
}
calcXmx "$@"

intervals() {
	if [[ $# -eq 0 ]]; then
		usage
		return
	fi
	local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP scalar.ScalarIntervals $@"
	#echo $CMD >&2
	eval $CMD
}

intervals "$@"
