#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 12, 2019

Description:  Creates a blacklist sketch from common kmers, 
which occur in at least X different sketches or taxa.
BlacklistMaker2 makes blacklists from sketches rather than sequences.
It is advisable to make the input sketches larger than normal,
e.g. sizemult=2, because new kmers will be introduced in the final
sketches to replace the blacklisted kmers.

Usage:  sketchblacklist2.sh ref=<sketch files> out=<sketch file>
or      sketchblacklist2.sh *.sketch out=<sketch file>
or      sketchblacklist2.sh ref=taxa#.sketch out=<sketch file>

Standard parameters:
ref=<file>          Sketch files.
out=<file>          Output filename.
mintaxcount=20      Retain keys occuring in at least this many taxa.
length=300000       Retain at most this many keys (prioritizing high count).
k=32,24             Kmer lengths, 1-32.
mode=taxa           Possible modes:
                       sequence: Count kmers once per sketch.
                       taxa: Count kmers once per taxonomic unit.
name=               Set the blacklist sketch name.
delta=t             Delta-compress sketches.
a48=t               Encode sketches as ASCII-48 rather than hex.
amino=f             Amino-acid mode.

Taxonomy-specific parameters:
tree=               Specify a taxtree file.  On Genepool, use 'auto'.
gi=                 Specify a gitable file.  On Genepool, use 'auto'.
accession=          Specify one or more comma-delimited NCBI accession to
                    taxid files.  On Genepool, use 'auto'.
taxlevel=subspecies Taxa hits below this rank will be promoted and merged
                    with others.
tossjunk=t          For taxa mode, discard taxonomically uninformative
                    sequences.  This includes sequences with no taxid,
                    with a tax level NO_RANK, of parent taxid of LIFE.

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

z="-Xmx31g"
z2="-Xms31g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 4000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

sketchblacklist() {
	local CMD="java $EA $SIMD $EOOM $z $z2 -cp $CP sketch.BlacklistMaker2 $@"
	echo $CMD >&2
	eval $CMD
}

sketchblacklist "$@"
