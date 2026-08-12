#!/bin/bash

usage(){
echo "
Written by Eru
Last modified August 10, 2026

Description:  Builds one consensus sequence per pre-computed protein family (the
hybrid foundation step of MAG-QC).  Keeps an existing clustering's membership
(e.g. mmseqs), but replaces each family's representative with the CONSENSUS of
its members, so the search representative fits divergent members instead of
being one arbitrary longest member.  The consensus reps are then re-searched
(by mmseqs) to rebuild the feature table.

This is a DEVELOPMENT tool (feature-database construction); it is not part of
the end-user MAG-QC path.

Usage:  consensusrepbuilder.sh seqs=<members.faa> clusters=<cluster.tsv> \\
          familylist=<familylist.tsv> out=<consensus_reps.fasta> [options]

Required parameters:
seqs=<file>     Member protein sequences (faa= and sequences= also accepted).
clusters=<file> Cluster membership TSV (rep<tab>member, mmseqs cluster format).
familylist=<f>  Family list TSV (rank<tab>cluster_rep<tab>...); defines which
                families get a consensus (families= also accepted).
out=<file>      Output consensus representative FASTA.

Optional parameters (and their defaults):
maxmembers=300  Max members aligned into one family's consensus (subsampled
                above this).
minmembers=1    Skip families with fewer members than this.
maxmemberlen=6000  Skip members longer than this (aligner cost guard).
passes=1        Consensus refinement passes.
pad=20          Graph end padding.
ani=f           Weight members by identity to the seed when voting.
ceiling=40      Identity ceiling (percent) for weighting saturation.
mafsub=0.25     Minor-allele fraction to call a substitution column.
mafdel=0.5      Minor-allele fraction to call a deletion column.
mafins=0.5      Minor-allele fraction to call an insertion column.
trimdepth=0.1   Trim consensus ends whose member depth fraction is below this.
t=auto          Worker threads.
ow=t            Overwrite existing output.

Java Parameters:
-Xmx            Set Java heap (overrides autodetection).
-eoom           Exit if an out-of-memory exception occurs.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
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
	DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=8g" "--xms=1g" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prot.ConsensusRepBuilder $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
