#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified June 20, 2025

Description:  Grades metagenome bins for completeness and contamination.
The contigs can be labeled with their taxID; in which case the header should
contain 'tid_X' somewhere where X is a number unique to their proper genome.
Alternately, CheckM2 and/or EukCC output can be fed to it.
Do not include a 'chaff' file (for unbinned contigs) when grading.
Completeness Score is (sum of completeness*size)/(total size) for all bins.
Contamination Score is (sum of contam*size)/(total size) for all bins.
Total Score is (sum of (completeness-5*contam)^2) for all bins.
Bin Definitions:
UHQ: >=99% complete and <=1% contam (subset of VHQ)
VHQ: >=95% complete and <=2% contam (subset of HQ)
HQ:  >=90% complete and <=5% contam
MQ:  >=50% complete and <=10% contam, but not HQ
LQ:  <50% complete or >10% contam
VLQ: <20% complete or >5% contam    (subset of LQ)

Usage:  gradebins.sh ref=assembly bin*.fa
or
gradebins.sh ref=assembly.fa in=bin_directory
or
gradebins.sh taxin=tax.txt in=bins

Input parameters:
ref=<file>      The original assembly that was binned.
in=<directory>  Location of bin fastas.
checkm=<file>   Optional CheckM2 quality_report.tsv file or directory.
eukcc=<file>    Optional EukCC eukcc.csv file or directory.
cami=<file>     Optional binning file from CAMI which indicates contig TaxIDs.
taxin=<file>    Optional file with taxIDs and sizes (instead of loading ref).
                Does not need to include taxIDs.  The tax file loads faster.
gtdb=<file>     Optional gtdbtk file.
gff=<file>      Optional gff file.
imgmap=<file>   Optional IMG map file, for renamed IMG gff input.
spectra=<file>  Optional path to QuickClade index.
cov=<file>      Optional path to QuickBin coverage file.
loadmt=t        Load bins multithreaded.

Output parameters:
report=<file>   Report on bin size, quality, and taxonomy.
taxout=<file>   Generate a tax file from the reference (for use with taxin).
hist=<file>     Cumulative bin size and contamination histogram.
ccplot=<file>   Per-bin completeness/contam data.
contamhist=<file> Histogram plotting #bins or bases vs %contam.

Processing parameters:
userna=f        Require rRNAs and tRNAs for HQ genomes.  This needs either
                a gff file or the callgenes flag.  Specifically, HQ and
                subtypes require at least 1 16S, 23S, and 5S, plus 18 tRNAs.
callgenes=f     Call rRNAs and tRNAs.  Suboptimal for some RNA types.
aligner=ssa2    Do not change this.
quickclade=f    Assign taxonomy using QuickClade.


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
}
calcXmx "$@"

gradeBins() {
	local CMD="java $EA $EOOM $SIMD $z -cp $CP bin.GradeBins $@"
	echo $CMD >&2
	eval $CMD
}

gradeBins "$@"
