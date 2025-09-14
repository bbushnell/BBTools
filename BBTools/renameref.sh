#!/bin/bash

usage(){
echo "
Written by Isla and Brian Bushnell
Last modified July 15, 2025

Description:  Converts reference sequence names in genomics files,
supporting SAM, BAM, FASTA, VCF, and GFF.  Updates reference names in headers
and data records according to a mapping file.  Useful for converting between
reference naming conventions (e.g. HG19 <-> GRCh37).
Sequence names not in the mapping file are kept as-is.  Name mapping will
first be attempted using the full header, and secondly using the prefix
of the original name up to the first whitespace.

Usage:
renameref.sh in=<input file> out=<output file> mapping=<ref_mapping.tsv>

Examples:
renameref.sh in=aligned.sam out=converted.sam mapping=hg19_to_grch37.tsv
renameref.sh in=data.sam out=renamed.sam mapping=refs.tsv strict=true

Parameters:
in=<file>       Input file to process
out=<file>      Output file with converted reference names  
map=<file>      Tab-delimited file with old_name<tab>new_name mappings
invert=<bool>   Reverse the order of names in the map file.
strict=<bool>   Crash on unknown references (default: false)
verbose=<bool>  Print detailed progress information (default: false)

Mapping file format:
chr1	1
chr2	2
chrX	X
chrM	MT

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

calcXmx () {
    # Source the new scripts
    source "$DIR""/memdetect.sh"
    source "$DIR""/javasetup.sh"
    
    parseJavaArgs "--mem=1g" "--mode=fixed" "$@"
    
    # Set environment paths
    setEnvironment
}
calcXmx "$@"

renameref() {
	local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.RefRenamer $@"
	#echo $CMD >&2
	eval $CMD
}

renameref "$@"
