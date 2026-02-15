#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified Jan 7, 2020

Description:   Prints the full taxonomy of a string.
String may be a gi number, NCBI taxID, or Latin name.
An NCBI identifier should just be a number or ncbi|number.
A gi number should be gi|number.
Please read bbmap/docs/guides/TaxonomyGuide.txt for more information.
Not: It is more convenient to use taxonomy.jgi-psf.org.

Usage:  taxonomy.sh tree=<tree file> <identifier>
Alternate usage: taxonomy.sh tree=<tree file> in=<file>

Usage examples:
taxonomy.sh tree=tree.taxtree.gz homo_sapiens canis_lupus 9606
taxonomy.sh tree=tree.taxtree.gz gi=gitable.int1.d.gz in=refseq.fasta

Processing parameters:
in=<file>       A file containing named sequences, or just the names.
out=<file>      Output file.  If blank, use stdout.
tree=<file>     Specify a TaxTree file like tree.taxtree.gz.
                On Genepool, use 'auto'.
gi=<file>       Specify a gitable file like gitable.int1d.gz. Only needed
                if gi numbers will be used.  On Genepool, use 'auto'.
accession=      Specify one or more comma-delimited NCBI accession to taxid
                files.  Only needed if accessions will be used; requires ~45GB
                of memory.  On Genepool, use 'auto'.
level=null      Set to a taxonomic level like phylum to just print that level.
minlevel=-1     For multi-level printing, do not print levels below this.
maxlevel=life   For multi-level printing, do not print levels above this.
silva=f         Parse headers using Silva or semicolon-delimited syntax.
taxpath=auto    Set the path to taxonomy files; auto only works at NERSC.

Parameters without an '=' symbol will be considered organism identifiers.

* Note *
Tree and table files are in /global/projectb/sandbox/gaag/bbtools/tax
For non-Genepool users, or to make new ones, use taxtree.sh and gitable.sh

Java Parameters:
-Xmx            This will set Java's memory usage,
                overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify
                200 megs.  The max is typically 85% of physical memory.
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

	parseJavaArgs "--xmx=2g" "--xms=2g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP tax.PrintTaxonomy $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
