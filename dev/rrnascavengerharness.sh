#!/bin/bash

usage(){
echo "
Development-only direct rRNA scavenger runner.

Usage: rrnascavengerharness.sh in=<fasta> out=<gff> consensus=<fa> kmers=<fa>
       [family=5S k=15 minlen=90 windowpad=100 indextopn=12 adaptive=t fixedminhits=12 quantumthresh=120 normalizecase=f collapsefrac=.85 idpass=.60 idborderline=.60 hbm=<models.hbm> hbmpass=.75 workloadout=<tsv>]

This command deliberately does not register a resource with CallGenes.  It runs
NcrnaScavenger directly, with no HBM or boundary network, for reproducible
seed/consensus A/B canaries.  Output is generic RNA GFF with model:<family>.
Input records must be complete contigs, not pre-shredded windows.  workloadout=
emits raw zero-based inclusive per-strand coordinates; strand 1 is in the
reverse-complement frame, not contig-absolute/GFF coordinates.
When hbm= is supplied, its model records must be index-aligned with consensus=
and borderline candidates may be accepted at hbmpass= (default .75).
normalizecase=t clones and uppercases each input contig for an explicitly
opt-in soft-mask compatibility canary; it does not alter CallGenes.
"
}

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
	elif [ -d "$DIR""../current/" ]; then
		CP="$DIR""../current/"
	fi
}

setEnv(){
	DIR="$DIR""../"
	. "$DIR""javasetup.sh"
	parseJavaArgs "--xmx=1g" "--xms=1g" "--mode=fixed" "$@"
	setEnvironment
}

launch(){
	if [ "$#" -eq 0 ]; then usage; return 0; fi
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prok.RrnaScavengerHarness $@"
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" prok.RrnaScavengerHarness "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
