#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 6, 2022

Description:  Sorts reads by name or other keys such as length,
quality, mapping position, flowcell coordinates, or taxonomy.
Writes temp files if memory is exceeded.

Usage:   sortbyname.sh in=<file> out=<file>

Input may be fasta, fastq, or sam, compressed or uncompressed.
Temp files will use the same format as the output.  Pairs are
kept together if reads are paired, and in2/out2 may be used for that.

Example 1 - sort by name:
sortbyname.sh in=raw.fq out=sorted.fq
Example 2 - sort by sequence:
sortbyname.sh in=raw.fq out=sorted.fq sequence
Example 3 - sort by mapping position:
sortbyname.sh in=mapped.sam out=sorted.sam position

Parameters:

in=<file>       Input file.
out=<file>      Output file.
name=t          Sort reads by name.
length=f        Sort reads by length.
quality=f       Sort reads by quality.
position=f      Sort reads by position (for mapped reads).
taxa=f          Sort reads by taxonomy (for NCBI naming convention).
sequence=f      Sort reads by sequence, alphabetically.
clump=f         Sort reads by shared kmers, like Clumpify.
flowcell=f      Sort reads by flowcell coordinates.
shuffle=f       Shuffle reads randomly (untested).
list=<file>     Sort reads according to this list of names.
ascending=t     Sort ascending.

Memory parameters (you might reduce these if you experience a crash)
memmult=0.30    Write a temp file when used memory exceeds this fraction
                of available memory.
memlimit=0.65   Wait for temp files to finish writing until used memory
                drops below this fraction of available memory.
delete=t        Delete temporary files.
allowtemp=t     Allow writing temporary files.

Taxonomy-sorting parameters (for taxa mode only):
tree=           Specify a taxtree file.  On Genepool, use 'auto'.
gi=             Specify a gitable file.  On Genepool, use 'auto'.
accession=      Specify one or more comma-delimited NCBI accession to
                taxid files.  On Genepool, use 'auto'.

Note: name, length, and quality are mutually exclusive.
Sorting by quality actually sorts by average expected error rate,
so ascending will place the highest-quality reads first.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding
                autodetection.  -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an 
                out-of-memory exception occurs.  Requires Java 8u92+.
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

z="-Xmx2g"
z2="-Xms2g"
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
	freeRam 2000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

sortbyname() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP sort.SortByName $@"
	echo $CMD >&2
	eval $CMD
}

sortbyname "$@"
