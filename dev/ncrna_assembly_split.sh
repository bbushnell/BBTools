#!/bin/bash
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_assembly_split.sh in=original_corpus.genomes.tsv
 groups=resolved.inputs.tsv.gz taxonomy=bac120_taxonomy.tsv.gz,ar53_taxonomy.tsv.gz
 out=fresh_prefix seed=20261113 expectedgenomes=39979 t=1 -Xmx512m -eoom
 Optional coassign=reviewed_pair.tsv has pool_group, partner_group, reason columns.

Assigns resolved assembly groups once for all families, without using source
locus counts or historical panel membership. Strata use GTDB genus, falling
back to family; fewer than five groups all train as coverage-critical.
Other strata hold out floor(group_count/5), using sorted strata and group IDs
then Java Random shuffle with the recorded seed. Missing/conflicting taxonomy
or unresolved grouping fails. One explicitly reviewed unresolved pair may share
a side, retaining both pool groups and labelling any outside-pool aliases. No
identity is claimed. If two pool partners land on opposite sides or strata,
the tool stops for a new disposition rather than silently changing quotas.
This is not member admission or sequence-seen grading.
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
CP="${DIR}current/"
if [ -f "${DIR}bbtools.jar" ]; then CP="${DIR}bbtools.jar"; fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=512m --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.NcrnaAssemblySplit "$@"
