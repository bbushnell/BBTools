#!/usr/bin/env bash
# Source-relative CDS metrics and changed-start evidence from paired saved GFFs.
set -eo pipefail
if [[ $# == 0 || ${1:-} == --help ]];then
 cat <<'HELP'
Usage: dev/gradecdsstartchanges.sh in=before.calls.tsv in2=after.calls.tsv corpus=formal.MANIFEST.tsv expectedinputs=N out=fresh_prefix t=8 -Xmx4g
Reads original CDS-containing annotations; no DNA or caller replay.
Optional rnacallchanges=t writes .rna_call_changes.tsv in the same saved-GFF pass.
It requires familymap= with the same explicit model-to-scoring-family projection
as RNA grading. Included target families are excluded from other-RNA losses;
SRP/6S subtype switches retain their shared scoring-family identity.
Exact other-RNA identity includes input, contig, type, family, strand and both
coordinates. Duplicates retain multiplicity; model-only changes are retained.
An endpoint shift is loss plus gain,
not proven biological gene loss. Default is false and preserves all historical
output tables and schemas.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=4g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.GradeCdsStartChanges "$@"
