#!/usr/bin/env bash
# Exports all native-score-passing domain evidence for unnamed BBGAP families.
# Shared JVM setup permits unset optional variables; match other BBTools launchers.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
	cat <<'HELP'
Usage: bbgapcoverageaudit.sh decisions=decisions.tsv targets=targets.compact.tsv \
    domains=domains.compact.tsv out=NEW_PREFIX

Inventories families with decision_class=unnamed and
rule_or_reason=insufficient_native_ncbi_evidence. Reads the original compact
target/domain tables, not the one-domain-per-model candidate table.
Targets must contain HMMER's reported (rep) count, not its detected (dom) count.
Recompact older evidence with the corrected BbgapLabelEvidence producer.
Missing or inconsistent reported-domain rows fail before any output is written.

All native-score-passing domains are retained, including non-naming models and
alternative domains. No model ranking or new naming rule is applied.
out=PREFIX       Writes .domains.tsv, .families.tsv and .summary.tsv.
mincov=0.8       Classifies each domain as none/model_only/query_only/both failing.
ow=f            Existing outputs are protected by default.

Family output retains original decision fields, including available Swiss
evidence. Domain-count categories may overlap at family level. Model type and
for_naming remain explicit; no best-raw-score row stands in for alternatives.
HELP
	exit 0
fi
script=$(readlink -f -- "${BASH_SOURCE[0]}")
repo=$(cd "$(dirname "$script")/.." && pwd)
DIR=$repo
source "$repo/javasetup.sh"
source "$repo/memdetect.sh"
args=(); jvm_args=()
for arg in "$@"; do
	case ${arg,,} in
		--xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
		-ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
		simd|simd=t|simd=true|simd=f|simd=false|nosimd) jvm_args+=("$arg");;
		*) args+=("$arg");;
	esac
done
parseJavaArgs --xmx=1g --xms=128m --mode=fixed "${jvm_args[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$repo/current" prot.BbgapCoverageAudit "${args[@]}"
