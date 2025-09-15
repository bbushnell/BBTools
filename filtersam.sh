#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified September 20, 2022

Description:  Filters a sam file to remove reads with variations unsupported
by other reads (bad vars, aka bad subs).  For particularly bad data,
it may be advisable to iteratively re-call variants and re-run FilterSam.
Calling variants may be performed like this:

callvariants.sh in=mapped.sam ref=ref.fa out=vars.vcf clearfilters minreads=2

Usage:  filtersam.sh in=<file> out=<file> vcf=<file>

Parameters:
in=<file>       Input sam or bam file.
ref=<file>      Optional fasta reference file.
out=<file>      Output file for good reads.
outb=<file>     Output file for bad reads.
vcf=<file>      VCF file of variants called from these reads.
vars=<file>     Alternatively, variants can be provided in CallVariants'
                native output format.
mbv=2           (maxbadvars) Discarded reads with more bad vars than this.
mbad=2          (maxbadalleledepth) A var is bad if the allele depth is at
                most this much.
mbaf=0.01       (maxbadalleledepth) A var is bad if the allele fraction is at
                most this much.  The more stringent of mbad or mbaf is used,
                so in low depth regions mbad is dominant while in high depth 
                regions mbaf is more important.  Vars are considered bad if
                they fail either threshold (meaning ad<=mbad or af<=mbaf).
mbrd=2          (minbadreaddepth) Substitutions may only be considered
                bad if the total read depth spanning the variant is
                at least this much.
border=5        (minenddist) Ignore vars within this distance of read ends.
sub=t           Consider bad substitutions.
ins=f           Consider bad insertions.
del=f           Consider bad deletions.

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

z="-Xmx8g"
z2="-Xms8g"
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
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

filtersam() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP var2.FilterSam $@"
#	echo $CMD >&2
	eval $CMD
}

filtersam "$@"
