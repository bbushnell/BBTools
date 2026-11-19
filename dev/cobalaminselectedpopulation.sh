#!/bin/bash
# Selected assembly/source-bound native-CM population and reference comparison.
set -eo pipefail
if [[ $# == 0 || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/cobalaminselectedpopulation.sh in=MANIFEST.tsv provenance=original_MANIFEST.tsv cohorts=bindings.tsv out=fresh_directory expectedgenomes=N t=16 -Xmx4g -eoom
Recounts RF00174/RF01482/RF01689 from saved native-GA events for every selected
assembly. No DNA reads or model fitting. Native truncations remain in recall;
only loci without an untruncated event are excluded from construction counts.
Reference parent counts are separate: fresh scored CM, scored protected prior,
and don't-care. Saved witnesses compare each corresponding native event to its
scored GFF parent. Nonidentical links block construction pending interpretation.
cohorts= has columns cohort, reference, events, witnesses (tab-separated).
provenance= binds selected input IDs and exact physical source FASTA/GFF identities
to their original scan cohorts. It is an identity map, never another input roster.
Outputs include event/locus ledgers, per-input/summary counts, reference categories,
all raw/reference links, bounded mismatch examples, and result.tsv written last.
HELP
 exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [[ -h $SCRIPT ]];do
 cd "$(dirname "$SCRIPT")"
 nextScript="$(readlink "$(basename "$SCRIPT")")"
 SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=4g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "${DIR}current/" prok.CobalaminSelectedPopulation "$@"
