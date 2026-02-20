#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 1, 2017

Description:  Runs stats.sh on multiple assemblies to produce one output line per file.

Usage:  statswrapper.sh in=<input file>

Parameters:
in=<file>       Specify the input fasta file, or stdin.  For multiple files a, b, and c: 'statswrapper.sh in=a,b,c'.
                'in=' may be omitted if this is the first arg, and asterisks may be used; e.g. statswrapper.sh *.fa
gc=<file>       Writes ACGTN content per scaffold to a file.
gchist=<file>   Filename to output scaffold gc content histogram.
gcbins=<200>    Number of bins for gc histogram.
n=<10>          Number of contiguous Ns to signify a break between contigs.
k=<13>          Estimate memory usage of BBMap with this kmer length.
minscaf=<0>     Ignore scaffolds shorter than this.
n_=<t>          This flag will prefix the terms 'contigs' and 'scaffolds' with 'n_' in formats 3-6.
addname=<t>     Adds a column for input file name, for formats 3-6.

format=<1 through 6>   Format of the stats information.  Default is format=3.
   format=1 uses variable units like MB and KB, and is designed for compatibility with existing tools.
   format=2 uses only whole numbers of bases, with no commas in numbers, and is designed for machine parsing.
   format=3 outputs stats in 2 rows of tab-delimited columns: a header row and a data row.
   format=4 is like 3 but with scaffold data only.
   format=5 is like 3 but with contig data only.
   format=6 is like 3 but the header starts with a #.

gcformat=<1 or 2>   Select GC output format.
   gcformat=1:   name   start   stop   A   C   G   T   N   GC
   gcformat=2:   name   GC
   Note that in gcformat 1, A+C+G+T=1 even when N is nonzero.
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

	parseJavaArgs "--xmx=200m" "--xms=200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.AssemblyStatsWrapper format=3 $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
