#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Neptune
Last modified August 6, 2026

Description:  Builds per-anticodon tRNA consensus sequences and HBM
(Hidden Brian Model) profiles from a collection of tRNA sequences.
Uses consensus-centroid clustering with reassignment for optimal
cluster representatives.  Output can be used with callgenes.sh
via the trnalib= and trnamodel= flags.

Usage:  trnaconsensus.sh in=trnas.fa out=consensus.fa outmodel=models.hbm

File parameters:
in=<file>       Input fasta of tRNA sequences (with anticodon in header).
out=<file>      Output consensus fasta.
outmodel=<file> Output HBM model file (.hbm for text, .bgm for binary).

Clustering parameters:
cluster=t       Enable clustering (default true).
clusterid=0.75  Clustering identity threshold.
mincluster=3    Minimum cluster size to keep.
mingroup=10     Minimum sequences per anticodon group.

Recruitment parameters:
recruit=t       Recruit orphans into clusters via consensus alignment.
recruitid=0.70  Minimum identity to consensus for recruitment.

Consensus parameters:
passes=2        Number of consensus refinement passes.
minid=0.3       Minimum alignment identity for consensus building.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM.  The max is typically 85% of
                physical memory.  The human genome requires around 24g minimum.

Please contact Brian at bbushnell@lbl.gov if you encounter any problems.
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

resolveSymlinks(){
	if [ -d "$DIR""current/" ]; then
		CP="$DIR""current/"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=8g" "--xms=8g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prok.TrnaConsensusBuilder $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
