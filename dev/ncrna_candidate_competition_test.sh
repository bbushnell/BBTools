#!/bin/bash
# Permanent full-caller TPP/CDS competition and candidate-reuse regression.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: ncrna_candidate_competition_test.sh in=original_genome.fa.gz reference=registered.tsv familymap=TPP_rules.tsv out=fresh_directory
Default fixture requires unchanged Microvirgula aerodenitrificans DSM15089 GCF_000620105.1.
A separately identified SAM-I fixture supplies target=sam1 gradefamily=SAM_I,
contig= strand= rnastart= rnastop= cdsstart= cdsstop=, mincdslen= baselinecap=128 floor=4.
Coordinates must be established independently from ordinary-caller controls.
Runs fresh ordinary callers first, then reuse off/on and native grading.
Uses a fixed artificial overlap/score/CDS-length configuration, never a production setting.
Parent heap defaults to 512 MiB; sequential caller children use 4096 MiB.
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
parseJavaArgs --xmx=512m --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.NcrnaCandidateCompetitionTest "$@"
