#!/bin/bash
# Reconcile existing ncRNA calls and member coordinates; never read genome payloads.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]];then
	cat <<'HELP'
Usage: dev/reconcilencrnadata.sh in=cohorts.tsv mapping=members.tsv groups=groups.inputs.tsv.gz roots=root1,root2,root3,root4 genomes=corpus.genomes.tsv assignments=assignments.tsv ledger=train.members.tsv proxy=corpus.loci.tsv.gz recovered=recovered.ledger.tsv overlays=overlay1,overlay2 out=fresh_directory t=1 -Xmx2g -eoom
Reads existing calls, alias tables and corrected member-coordinate metadata.
Lists immediate .fna.gz source-root members, without opening genome payloads.
Preserves raw RF00174 events, partial/ambiguous evidence and protected priors.
Exact-coordinate groups and same-frame overlapping inclusive components are
reported separately. Components never merge adjacency or unproved aliases.
Native RF00174 events with valid coordinates and native_included=true qualify,
including partial and prior-covered hits. Saved annotations and protected-prior
support are reported separately through the existing subtype reader. Member exact bounds and
reciprocal-80-percent overlap are reported separately; overlap is not identity.
Whole-assembly equality or accession-looking names cannot map coordinates across
different source files. Those comparisons stay unresolved without a contig map.
Writes assembly coverage, unscanned paths, annotations, events, witnesses, exact
loci, merged_loci.tsv.gz, component links, saved-census source-frame cost inputs,
member-coordinate matches and summary tables. No model/DDL build or caller run.
The cohort manifest has seven tab-separated columns:
cohort,census,reference,events,witnesses,expected_inputs,expected_scored_RF00174.
Set expectedcohorts=3 for the three-cohort extension (default2). Name the original
and audited prior-import cohorts original/import; new cohort names are distinct.
All census paths bind before live coverage is computed. Historical population
checks default to 40,723 alias-audit groups, 25,916 original
scored loci, 57,723 confident imported events and 85,726 kept library members.
Outputs must be fresh. Compile only on Dori.
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
if [ -f "${DIR}bbtools.jar" ];then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=2g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.NcrnaDataReconciliation "$@"
