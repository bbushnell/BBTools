#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified June 18, 2018

Description:   Filters sequences according to their taxonomy,
as determined by the sequence name.  Sequences should
be labeled with a gi number, NCBI taxID, or species name.

Usage:  filterbytaxa.sh in=<input file> out=<output file> tree=<tree file> table=<table file> ids=<numbers> level=<name or number>

I/O parameters:
in=<file>       Primary input, or read 1 input.
out=<file>      Primary output, or read 1 output.
results=<file>  Optional; prints a list indicating which taxa were retained.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
showspeed=t     (ss) Set to 'f' to suppress display of processing speed.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression
                level; lower compression is faster.

Processing parameters:
level=          Taxonomic level, such as phylum.  Filtering will operate on
                sequences within the same taxonomic level as specified ids.
                If not set, only matches to a node or its descendants will 
                be considered.
reqlevel=       Require nodes to have ancestors at these levels.  For example,
                reqlevel=species,genus would ban nodes that are not defined
                at both the species and genus levels.
ids=            Comma-delimited list of NCBI numeric IDs.  Can also be a
                file with one taxID per line.  Names (like bacteria) are also
                acceptable.
include=f       'f' will discard filtered sequences, 't' will keep them.
besteffort=f    Intended for include mode.  Iteratively increases level
                while the input file has no hits to the tax list.
tree=<file>     Specify a TaxTree file like tree.taxtree.gz.  
                On Genepool, use 'auto'.
gi=<file>       Specify a gitable file like gitable.int1d.gz. Only needed
                if gi numbers will be used.  On Genepool, use 'auto'.
accession=      Specify one or more comma-delimited NCBI accession to taxid
                files.  Only needed if accessions will be used; requires ~45GB
                of memory.  On Genepool, use 'auto'.
printnodes=t    Print the names of nodes added to the filter.
requirepresent=t   Crash with an error message if a header cannot be resolved
                   to a taxid.

String-matching parameters:
regex=          Filter names matching this Java regular expression.
contains=       Filter names containing this substring (case-insensitive).

* Note *
Tree and table files are in /global/projectb/sandbox/gaag/bbtools/tax
For non-Genepool users, or to make new ones, use taxtree.sh and gitable.sh

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

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

z="-Xmx4g"
z2="-Xms4g"
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
	freeRam 1000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

filterbytaxa() {
	local CMD="java $EA $SIMD $EOOM $z -cp $CP tax.FilterByTaxa $@"
	echo $CMD >&2
	eval $CMD
}

filterbytaxa "$@"
