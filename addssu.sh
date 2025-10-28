#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified Jan 29, 2020

Description:  Adds, removes, or replaces SSU sequence of existing sketches.
Sketches and SSU fasta files must be annotated with TaxIDs.

Usage:           addssu.sh in=a.sketch out=b.sketch 16S=16S.fa 18S=18S.fa

Standard parameters:
in=<file>       Input sketch file.
out=<file>      Output sketch file.

Additional file parameters (optional):
16S=<file>      A fasta file of 16S sequences.  These should be renamed
                so that they start with tid|# where # is the taxID.
                Should not contain organelle rRNA.
18S=<file>      A fasta file of 18S sequences.  These should be renamed
                so that they start with tid|# where # is the taxID.
                Should not contain organelle rRNA.
tree=auto       Path to TaxTree, if performing prok/euk-specific operations.

Processing parameters:
preferSSUMap=f
preferSSUMapEuks=f
preferSSUMapProks=f
SSUMapOnly=f
SSUMapOnlyEuks=f
SSUMapOnlyProks=f
clear16S=f
clear18S=f
clear16SEuks=f
clear18SEuks=f
clear16SProks=f
clear18SProks=f


Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-da             Disable assertions.

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
	#freeRam 3200m 84
	#z="-Xmx${RAM}m"
	#z2="-Xms${RAM}m"
}
calcXmx "$@"

sendsketch() {
	local CMD="java $EA $SIMD $EOOM $z -cp $CP sketch.AddSSU $@"
	echo $CMD >&2
	eval $CMD
}

sendsketch "$@"
