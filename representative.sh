#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified September 4, 2019

Description:  Makes a representative set of taxa from all-to-all identity
comparison.  Input should be in 3+ column TSV format (first 3 are required):
(query, ref, ANI, qsize, rsize, qbases, rbases)
...as produced by CompareSketch with format=3 and usetaxidname.
Additional columns are allowed and will be ignored.

Usage:  representative.sh in=<input file> out=<output file>

Parameters:
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
threshold=0     Ignore edges under threshold value.  This also affects the
                choice of centroids; a high threshold gives more weight to 
                higher-value edges.
minratio=0      Ignores edges with a ratio below this value.
invertratio=f   Invert the ratio when greater than 1.
printheader=t   Print a header line in the output.
printsize=t     Print the size of retained nodes.
printclusters=t Print the nodes subsumed by each retained node.
minsize=0       Ignore nodes under this size (in unique kmers).
maxsize=0       If positive, ignore nodes over this size (unique kmers).
minbases=0      Ignore nodes under this size (in total bases).
maxbases=0      If positive, ignore nodes over this size (total bases).

Taxonomy parameters:
level=          Taxonomic level, such as phylum.  Filtering will operate on
                sequences within the same taxonomic level as specified ids.
                If not set, only matches to a node or its descendants will 
                be considered.
ids=            Comma-delimited list of NCBI numeric IDs.  Can also be a
                file with one taxID per line.
names=          Alternately, a list of names (such as 'Homo sapiens').
                Note that spaces need special handling.
include=f       'f' will discard filtered sequences, 't' will keep them.
tree=<file>     Specify a TaxTree file like tree.taxtree.gz.  
                On Genepool, use 'auto'.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will
                specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                The max is typically around 85% of physical memory.
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
	CP="$DIR/current/"
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=4g" "--xms=4g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.RepresentativeSet $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
