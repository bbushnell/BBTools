#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 10, 2020

Description:  Generates the consensus sequence of a reference
using aligned sequences.  This can be used for polishing assemblies,
making representative ribosomal sub units, correcting PacBio reads, etc.

If unaligned sequences are used as input, they should be in fasta or fastq
format, and they will be aligned to the first reference sequence.

Usage:  consensus.sh in=mapped.sam ref=ref.fa out=consensus.fa

Recommended settings for assembly polishing via Illumina reads:  mafsub=0.5


Standard parameters:
in=<file>       Reads mapped to the reference; should be sam or bam.
ref=<file>      Reference; may be fasta or fastq.
out=<file>      Modified reference; may be fasta or fastq.
outm=<file>     Optional output for binary model file.
                Preferred extension is .alm.
inm=<file>      Optional input model file for statistics.
hist=<file>     Optional score histogram output.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Processing parameters:
mindepth=2      Do not change to alleles present at depth below this.
mafsub=0.25     Do not incorporate substitutions below this allele fraction.
mafdel=0.50     Do not incorporate deletions below this allele fraction.
mafins=0.50     Do not incorporate insertions below this allele fraction.
mafn=0.40       Do not change Ns (noref) to calls below this allele fraction.
usemapq=f       Include mapq as a positive factor in edge weight.
nonly=f         Only change Ns to different bases.
noindels=f      Don't allow indels.
ceiling=        If set, alignments will be weighted by their inverse identity.
                For example, at ceiling=105, a read with 96% identity will get
                bonus weight of 105-96=9 while a read with 70% identity will
                get 105-70=35.  This favors low-identity reads.
name=           Set the output sequence name (for a single output sequence).

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

#Note that this needs -Xss flag to prevent serialization stack overflow
consensus() {
	local CMD="java $EA $EOOM $z -Xss8m -cp $CP consensus.ConsensusMaker $@"
	echo $CMD >&2
	eval $CMD
}

consensus "$@"
