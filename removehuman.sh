#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified December 22, 2021
This script requires at least 16GB RAM.
It is designed for NERSC and uses hard-coded paths.

Description:  Removes all reads that map to the human genome with at least 95% identity after quality trimming.
Removes approximately 98.6% of human 2x150bp reads, with zero false-positives to non-animals.
NOTE!  This program uses hard-coded paths and will only run on Nersc systems unless you change the path.

Usage:  removehuman.sh in=<input file> outu=<clean output file>

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
path=               Set the path to an indexed human genome.

***** All BBMap parameters can be used; run bbmap.sh for more details. *****

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

	parseJavaArgs "--xmx=15000m" "--xms=15000m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP align2.BBMap minratio=0.9 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path=/global/cfs/cdirs/bbtools/hg19 pigz unpigz zl=6 qtrim=r trimq=10 untrim idtag usemodulo printunmappedcount ztd=2 kfilter=25 maxsites=1 k=14 bloomfilter $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
