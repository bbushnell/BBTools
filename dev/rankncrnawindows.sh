#!/bin/bash
# Member-voted DDL diagnostic on one ordinary native caller execution per contig.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/rankncrnawindows.sh caller=callgenes.config reference=reference.tsv familymap=families.tsv subtypes=subtypes.tsv witnesses=witnesses.tsv events=events.tsv members=members.fa ledger=train.members.tsv assignments=assignments.tsv census=models.tsv membercounts=membership.summary.tsv out=fresh_prefix t=1 -Xmx8g -eoom
Reads the exact configured genomes once through CallGenes. Requires taxonomy=f.
RF00174 models keep their native order; full construction headers bind every vote.
Defaults: expectedmodels=508 ks=5,7,9,11,13,15,19 widths=32,64,128,256,512,1024,2048
Each model is fully verified once per captured window. All three ranking methods
reuse that result: native SSU float score, equal sketch buckets, all raw window keys.
Seed, native-index, oracle, ranking and report CPU are recorded separately.
Per-k/width profile construction, query/scoring and rank-query CPU are retained
in .geometry_cpu.tsv; all three score methods share each query/scoring measurement.
Ordinary final paths are written to .calls.gff; whole-command time includes analysis.
Reports retain window bytes, complete native model counts and oracle outcomes,
every accepted span/parent link, all actual attempts,
per-event ranks and per-target minima, empty inputs, ambiguous spans and missing ranks.
The .summary.tsv gives both all-ranked and common-native-rankable distributions.
All scored cobalamin parents compete for unique attribution. Masks never consume it.
Minima across events are diagnostic and do not describe one actual caller walk.
Existing outputs, native CPU profiling, sampling and incomplete bindings reject.
Read limits or sampling inside the caller configuration also reject before output.
This is an in-corpus ranking experiment, not a resource rebuild or caller adoption.
Compile on Dori. Bulk inputs and traces belong on Dori, never on the laptop.
HELP
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ]; do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
if [ -f "${DIR}bbtools.jar" ]; then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=8g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.RankNcrnaWindows "$@"
