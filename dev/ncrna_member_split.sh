#!/bin/bash
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_member_split.sh in=saved_member_prefix split=common_split_prefix
 groups=original_resolved_inputs.tsv.gz out=fresh_prefix expectedinputs=41082
 expectedpool=39979 expectedold=26896 expectedevents=27666 expectedcores=27631
 t=1 bgzip=t zl=9 -Xmx1g -eoom

Preserves every historical admission/evidence field and appends common-group
methodology assignments and eligibility. Outside paths in pool groups follow
their group's side. Truly outside-pool old members retain eligible roles;
native-only outside cores remain recorded without expanding the training pool.
No admission replay, sequence reads, deduplication, consensus or voting.
Training candidate row counts are not deduplicated sequence weights.
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
parseJavaArgs --xmx=1g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.NcrnaMemberSplit "$@"
