#!/bin/bash
# Two-length member-voted ranking/statistics on saved windows only.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/rank_dual_ncrna_windows.sh panels=table selection=table baseline=corrected.loci.tsv.gz members=fa ledger=tsv assignments=tsv census=tsv summary=tsv ref=fa out=fresh_prefix t=16 -Xmx8g -eoom
Writes fixed native-model identities and quartiles before new ranks, the seeded
assembly-group split, sixteen two-length geometries, saved partition half-counts,
equal-locus signal/background moments, fixed/cross-fitted weights, complete target
rows and named last-reachable-input losses. No genomes, oracle or caller rerun.
Fitting and judging rows stay separate; missing identity/weights remain explicit.
The corrected ambiguity replay supplies all original native/masked10/hybrid
controls. Results describe saved ordering, not caller recall or adoption.
Optional identities=accepted_prefix reads .identity.tsv, .identity_bins.tsv and
.groups.tsv instead of measuring identities again. Original parent/model/event/
gene witnesses, float bits, bins, cutoff population and group assignment must
match. This path runs no identity alignments; omit it only for a first measurement.
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
parseJavaArgs --xmx=8g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$CP" prok.RankDualNcrnaWindows "$@"
