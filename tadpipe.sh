#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified September 18, 2019

Description:  Runs TadpoleWrapper after some preprocessing,
to allow optimal assemblies using long kmers.
Only paired reads are supported.

Usage:
tadpipe.sh in=reads.fq out=contigs.fa


Parameters:
in=<file>           Input reads.
in2=<file>          Optional read 2, if reads are in two files.
out=contigs.fa      Output file name.
temp=$TMPDIR        Path to a directory for temp files.
delete=t            Delete intermediate files.
gz=f                Gzip intermediate files.

Other parameters can be passed to individual phases like this:

assemble_k=200,250  Set kmer lengths for assembly phase.
merge_strict        Set the strict flag in merge phase.
extend_el=120       Set the left-extension distance in the extension phase.

Valid prefixes:

filter_             PhiX and contaminant filtering.
trim_               Adapter trimmming.
merge_              Paired-read merging.
correct_            Error correction.
extend_             Read extension.
assemble_           Final assembly.
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

z="-Xmx14g"
z2="-Xms14g"
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
	freeRam 15000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

tadpipe() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP assemble.TadPipe $@"
	echo $CMD >&2
	eval $CMD
}

tadpipe "$@"
