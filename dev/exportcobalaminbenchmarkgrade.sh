#!/usr/bin/env bash
# Export completed union grades; no caller or matching operation runs.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/exportcobalaminbenchmarkgrade.sh in=completed_grade out=fresh_prefix
For historical CobalaminPanelGrade: legacy=t panel=original_panel.loci.tsv

Writes .members.tsv, .genomes.tsv, .witnesses.tsv, .matches.tsv,
.queries.tsv, .references.tsv and .metrics.tsv. The Java metrics row gives total
calls and scored/masked reference counts separately from masked query calls.
Optional timing=native.time.tsv callerlog=complete.stderr copies actual
user/system/wall seconds and the three cobalamin alignment counters into metrics.
Query/model/coordinates remain available for
overlap-only loss classification and matched versus unmatched length diagnostics.
Emitted members partition UNION outcomes; these are not biological subtype grades.
The nine-column matches contain clan and member views for ncrnaendpointsummary.sh.
Offsets retain gene-sense signs and actual matched reference coordinates.
Zero-call genomes survive. No calling, rematching, DNA read or truth change.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.ExportCobalaminBenchmarkGrade "$@"
