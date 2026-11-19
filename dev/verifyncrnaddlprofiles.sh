#!/bin/bash
# Check the actual compact caller loader/scorer against saved parent ranks.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/verifyncrnaddlprofiles.sh in=profiles512.tsv in2=profiles1024.tsv \
  panels=panels.tsv selection=selection.tsv ref=consensus.fa \
  ranks=accepted_partition.loci.tsv.gz out=fresh_parity.tsv \
  expectedmodels=508 expectedgroups=3153 t=1 -Xmx2g -eoom

Loads both real resources against the complete original model roster and uses
the caller scorer on existing captured windows. Compares minimum acceptable-model
rank for EVERY original parent with the selected plain7/masked10/shared-consensus
512/1024 rows in the saved table, including unrankable/unknown parents.
Original panel/input/parent identity and order must match; zero-target inputs
remain in the reported population. Writes two small count rows only after exact
rank equality. Identity-bin labels are retained in the accepted report, not refit.
No member voting, resource rebuild, identities, alignments or genome caller.
The temporary indexparity flag has been retired after the accepted index
comparison. Supplying it now fails clearly. Ordinary saved-rank verification
uses the same shared-index scorer as the resource consumer.
Use the accepted input/output bindings in the reviewed execution recipe.
This is resource/rank parity, not biological caller acceptance or adoption.
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
parseJavaArgs --xmx=2g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.VerifyNcrnaDdlProfiles "$@"
