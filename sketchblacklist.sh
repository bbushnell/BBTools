#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 7, 2019

Description:  Creates a blacklist sketch from common kmers, 
which occur in at least X different sequences or taxa.
Please read bbmap/docs/guides/BBSketchGuide.txt for more information.

Usage:  sketchblacklist.sh in=<fasta file> out=<sketch file>

Standard parameters:
in=<file>           A fasta file containing one or more sequences.
out=<file>          Output filename.
mintaxcount=100     Sketch kmers occuring in at least this many taxa.
k=31                Kmer length, 1-32.  To maximize sensitivity and 
                    specificity, dual kmer lengths may be used:  k=31,24
mode=sequence       Possible modes:
                       sequence: Count kmers once per sequence.
                       taxa: Count kmers once per taxonomic unit.
name=               Set the blacklist sketch name.
delta=t             Delta-compress sketches.
a48=t               Encode sketches as ASCII-48 rather than hex.
amino=f             Amino-acid mode.
entropy=0.66        Ignore sequence with entropy below this value.
keyfraction=0.16    Smaller values reduce blacklist size by ignoring a
                    a fraction of the key space.  Range: 0.0001-0.5.

Taxonomy-specific parameters:
tree=               Specify a taxtree file.  On Genepool, use 'auto'.
gi=                 Specify a gitable file.  On Genepool, use 'auto'.
accession=          Specify one or more comma-delimited NCBI accession to
                    taxid files.  On Genepool, use 'auto'.
taxlevel=subspecies Taxa hits below this rank will be promoted and merged
                    with others.
prefilter=t         Use a bloom filter to ignore low-count kmers.
prepasses=2         Number of prefilter passes.
prehashes=2         Number of prefilter hashes.
prebits=-1          Manually override number of prefilter cell bits.
tossjunk=t          For taxa mode, discard taxonomically uninformative
                    sequences.  This includes sequences with no taxid,
                    with a tax level NO_RANK, of parent taxid of LIFE.
silva=f             Parse headers using Silva or semicolon-delimited syntax.

Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.

For more detailed information, please read /bbmap/docs/guides/BBSketchGuide.txt.
Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
"
}

if [ -z "$1" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	usage
	exit
fi

resolveSymlinks(){
	SCRIPT="$0"
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

	parseJavaArgs "--xmx=4g" "--xms=4g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP sketch.BlacklistMaker $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
