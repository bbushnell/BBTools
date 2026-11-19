#!/bin/bash
# Complete recorded-frame baseline preparation; metadata only.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]];then
	cat <<'HELP'
Usage: dev/prepare_ncrna_roster_baseline.sh in=source_frames.tsv.gz originals=original_calls.tsv.gz out=fresh_prefix expectedgenes=N expectedframes=N t=1 -Xmx1g -eoom
Retains all recorded canonical source frames, including older scanned assemblies,
zero-gene assemblies and unavailable inputs. No assembly-group deduplication.
Outputs .frames.tsv, .original_calls.tsv.gz, .reference.tsv, .inputs.config and
.summary.tsv. Only readable recorded frames enter the caller input list; missing
frames and their genes retain an explicit unresolved state in both ledgers.
The capture manifest's GFF column is unused and names placeholders, not truth.
The .inputs.config contains only in=; combine it with a reviewed caller recipe.
Current file availability/size checks are metadata, not proof of DNA identity.
No sequence read, resource generation, caller, grading or GC fit is performed.
HELP
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ];do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=1g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "${DIR}current/" prok.PrepareNcrnaRosterBaseline "$@"
