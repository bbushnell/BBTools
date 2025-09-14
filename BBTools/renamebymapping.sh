#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 21, 2025

Description:  Renames contigs based on mapping information.
Appends coverage and optionally taxID from parsing sam line headers.
For taxID renaming, read headers should contain a term like 'tid_1234';
output will be named as 'original tid_1234 cov_45.67' with potentially
multiple coverage entries (if there are multiple sam files) but
only one tid entry based on the highest-coverage sam file.
Designed for metagenome binning evaluation and synthetic read generation.

Usage:  renamebymapping.sh in=contigs.fa out=renamed.fa *.sam

Parameters:
in=<file>        Assembly to rename.
out=<file>       Renamed assembly.
sam=<file>       This can be a file, directory, or comma-delimited list.
                 Unrecognized arguments that are existing files will also
                 be treated as sam files.  Bam is acceptable too.
delimiter=space  Delimiter between appended fields.
wipe=f           Replace the original header with contig_#.
depth=t          Add a depth field.
tid=t            Add a tid field (if not already present).

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
	freeRam 4000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

rename() {
	local CMD="java $EA $EOOM $z -cp $CP bin.ContigRenamer $@"
	echo $CMD >&2
	eval $CMD
}

rename "$@"
