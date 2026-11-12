#!/bin/bash
# Paired accounting for completed native Infernal search pilots.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden
Last modified November 9, 2026

Usage: dev/comparecmsearchcosts.sh in=completed_pilot census=full_summary.tsv \
       population=sample_populations.tsv stage=64 out=fresh_prefix

Compares saved native cmsearch costs and GA-included hit multisets. Runs no search.
Requires the pilot bindings.tsv and successful root/point completion records.
Control commands must use --cut_ga; candidate commands must use -T and --incT
bound to the original model GA. All other search settings/input paths must match.
Simple unescaped recorded argv is required; arbitrary shell commands are refused.

in=          Completed pilot directory with bindings.tsv and search_STAGE/model/.
census=      Five-column full-genome population summary.
population=  Four-column stage/genomes/nucleotide_bases/fasta_records sample table.
stage=64     Selects the sample and native point directory.
out=         Fresh output prefix. Overwrite and append are unsupported.

Outputs: .model_costs.tsv, .hit_differences.tsv, .rounding.tsv, .cost_summary.tsv.
Hit identity uses literal target, mode, model/target bounds, strand and truncation.
Repeated target names in concatenated genomes cannot distinguish individual genomes.
Inclusion resolves only the native GA cutoff, not full length or biological truth.
Displayed scores are rounded; conflicts with the bound native flag are reported.
Lost/gained witnesses preserve duplicate multiplicity; scores/E-values are not keys.

CPU projections require at least ten CPU seconds in both paired timings. The total
new-pass cost is distinct from its growth over the historical control. Linear
base scaling is not measured full-scan cost or a worst-case bound. The eight-worker
wall scenario assumes perfect balance and excludes queue/startup/I/O/contention.
Report bytes include sample headers; a separate generous planning allowance assigns
each genome the whole sample report. File counts exclude per-input metadata and
scratch. Historical controls are nonconcurrent; no sensitivity/accuracy claim follows.
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
parseJavaArgs --xmx=512m --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.CompareCmSearchCosts "$@"
