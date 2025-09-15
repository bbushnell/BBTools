#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified December 22, 2021
This script requires at least 52GB RAM.
It is designed for NERSC and uses hard-coded paths.

Description:  Removes all reads that map to the cat, dog, mouse, or human genome with at least 95% identity after quality trimming.
Removes approximately 98.6% of human 2x150bp reads, with zero false-positives to non-animals.
NOTE!  This program uses hard-coded paths and will only run on Nersc systems.

Usage:  removecatdogmousehuman.sh in=<input file> outu=<clean output file>

Input may be fasta or fastq, compressed or uncompressed.

Parameters:
threads=auto        (t) Set number of threads to use; default is number of logical processors.
overwrite=t         (ow) Set to false to force the program to abort rather than overwrite an existing file.
interleaved=auto    (int) If true, forces fastq input to be paired and interleaved.
trim=t              Trim read ends to remove bases with quality below minq.
                    Values: t (trim both ends), f (neither end), r (right end only), l (left end only).
untrim=t            Undo the trimming after mapping.
minq=4              Trim quality threshold.
ziplevel=2          (zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster.
outm=<file>         File to output the reads that mapped to human.

***** All BBMap parameters can be used; run bbmap.sh for more details. *****

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
JNI="-Djava.library.path=""$DIR""jni/"
JNI=""

z="-Xmx50g"
z2="-Xms50g"
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
}
calcXmx "$@"

function removecatdogmousehuman() {
	local CMD="java $EA $EOOM $z $z2 $JNI -cp $CP align2.BBMap minratio=0.9 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path=/global/cfs/cdirs/bbtools/mousecatdoghuman/ pigz unpigz zl=6 qtrim=r trimq=10 untrim idtag usemodulo printunmappedcount ztd=2 kfilter=25 maxsites=1 k=14 bloomfilter $@"
	echo $CMD >&2
	eval $CMD
}

removecatdogmousehuman "$@"
