#!/bin/bash
# Select an identical explicit saved-caller scope from two hypothetical references.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden
Last modified November 8, 2026

Usage: dev/selectncrnascenarios.sh in=global_ga.tsv.gz in2=global_margin.tsv.gz \
       calls=calls.tsv bindings=bindings.tsv models=models.tsv scope=scope.tsv \
       callerconfig=original.config out=new_directory

Selects the same complete input population and enabled/supported families on both
sides. Included rows are copied unchanged. This is diagnostic scope selection,
not biological adjudication, threshold selection or an accuracy estimate.
No FASTA payload, caller or covariance-model scan is run.

in=/in2= Paired global dev/draftncrnascenarios.sh tables, in identical decision order.
    Keys are (input_id, decision_id); duplicate keys fail even among exclusions.
    Keys/families use canonical evidence TSV escaping (%25, %09, %0D, %0A).
    Other row validation remains the downstream comparator's responsibility.
calls= Complete three-column saved-call manifest: input_id, fasta, gff.
bindings= The dev/bindncrnacalls.sh bindings.tsv for that exact complete manifest.
    Original configuration inputs, ordered bindings and canonical paths must agree.
    Zero-call inputs remain in the denominator; no shortened manifest is accepted.
models= Tab-separated family, model_prefix, display_label with that literal header.
    Every generic RNA call must match exactly one prefix. A present family label
    must agree exactly. These are historical caller names, not current defaults.
scope= Tab-separated family, required_flags, supported with that literal header.
    required_flags is a comma-separated list of canonical lowercase gate names.
    All gates for a family must be true for it to be enabled. supported is t/f.
    Declared but unsupported families remain explicit counted exclusions.
callerconfig= Original flat caller configuration, with one absolute in= per input
    and every scope gate explicit exactly once. euk=f is mandatory. Nested configs,
    paired inputs, unsupported aliases and *only modes fail. ncrnafamily controls
    diagnostics/overrides and does not itself disable the other families.
    This is a restricted explicit-gate profile, not an arbitrary CallGenes replay.
out= Fresh directory, never overwritten. A failed run may leave partial reports;
    require exit0 and summary.tsv before using them.

Outputs:
ga.tsv.gz / margin.tsv.gz  Same selected paired keys, original row bytes.
familymap.tsv   Explicit four-column include/exclude map for dev/comparencrnascenarios.sh.
selection_counts.tsv Family, exclusion reason and paired before/after dispositions.
inputs.tsv      Complete original input population, query counts and selected rows.
scope.tsv       Explicit enabled/supported status and query counts for each family.
sources.tsv     Original global/manifest/config/mapping companion paths.
summary.tsv     Conserved global/selected/excluded counts, written last.

Queries from disabled families fail. Unknown global context families remain counted
exclusions; scored/masked global families require an explicit scope declaration.
Keep original global tables and their source-link/witness companions. Downstream
comparison uses this family's map, both selected tables, and the same original
calls manifest and adapter contigs.tsv.gz. Empty queries can still contribute FNs.
One analysis worker; -Xmx1g default; retained identity keys scale with global rows.
HELP
}
if [ $# -eq 0 ] || [ "$1" = -h ] || [ "$1" = --help ]; then
	usage
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ]; do
	cd "$(dirname "$SCRIPT")"
	RESOLVED_SCRIPT="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$RESOLVED_SCRIPT"
done
DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)/"
popd > /dev/null
if [ -f "${DIR}bbtools.jar" ]; then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=1g --xms=256m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.SelectNcrnaScenarios "$@"
