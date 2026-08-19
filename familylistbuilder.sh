#!/bin/bash
#familylistbuilder cluster=<tsv> reps=<fasta> out=<familylist.tsv> repsout=<top_reps.fasta>

usage(){
echo "
Written by Brian Bushnell and Eru
Last modified August 17, 2026

Description:  Builds a per-phylum-selected family list from mmseqs cluster output.
For each phylum with enough genomes, ranks its candidate families by in-phylum
prevalence (tiebreak: total prevalence, then rep id), takes the top N per phylum,
unions across phyla, and pads with the globally most prevalent remaining families
up to a floor count. Excludes tids on excluded= from all prevalence counting.
Output familylist.tsv contains ONLY the selected set (not every cluster).

Usage:  familylistbuilder.sh cluster=<cluster.tsv> reps=<rep_seq.fasta> taxpgm=<taxpgm.tsv> excluded=<EXCLUDED_TIDS.tsv> out=<familylist.tsv> repsout=<top_reps.fasta> [perphylumtop=100] [floor=6000] [minphylumgenomes=3] [tidsout=<tids.txt>]
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
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=16g" "--xms=16g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prot.FamilyListBuilder $@"
	echo $CMD >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
