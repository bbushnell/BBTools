#!/usr/bin/env bash
# Saved-only held-out endpoint comparison; no caller, grader or fit.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Render accepted saved boundary/detection tables as one HTML page and source TSV.
Usage: reportncrnaboundarycomparison.sh in=manifest.tsv out=fresh_directory \
 family=cobalamin mainmember=RF00174
Manifest columns: configuration, panel, state, prefix, comparison (tab-separated).
state=MEASURED uses normalized prefix files and comparison.loci/summary.tsv.
state=NOT_RUN uses literal . for both paths and renders explicit missing rows.
The prefix's per-family endpoint summary/radii tables must already exist.
Rates use the current matched emitted-member cohort. Lost/absent references use
baseline matches, with physical overlap and grading losses shown separately.
No input filtering, overwrite, calling, matching or training. One writer/thread.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.ReportNcrnaBoundaryComparison "$@"
