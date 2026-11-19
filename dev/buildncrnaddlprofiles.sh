#!/bin/bash
# Build one fixed model-bound RF00174 profile resource before caller experiments.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/buildncrnaddlprofiles.sh in=members.fa ledger=routing.tsv
  assignments=assignments.tsv census=membership.models.tsv
  summary=membership.summary.tsv ref=consensus.fa out=fresh_profiles.tsv
  expectedmodels=508 expectedmembers=85726 dual=f partitionbins=0 t=1 -Xmx2g -eoom

Joins the original frozen construction members, validates retained voter counts,
and votes one forward middle-masked k10/2048 profile per current model.
dual=t instead builds the fixed partitioned plain9+masked10/4096 resource,
scored as long + 0.5*short in sketch-only caller order.
Without selected=t, partitionbins=512 or 1024 builds the historical plain-7 member quarter,
masked-10 member half and shared consensus quarter in one sketch-only index.
It cannot be combined with dual=t. Other geometries and weights are not enabled
by this builder. Undefined sequence bases break words.
Seed12345/exponent4 and raw-query-key scoring are fixed by the versioned format.
The resource binds ordered model IDs and uppercase consensus bytes; changing
the consensus library requires rebuilding it. Voting occurs once, not per
worker or genome. No genome, verification, alignment or caller invocation.
This experimental resource does not enable a caller or change its defaults.

New all-gene mode:
  selected=t in=construction_members.fa assignments=builder.assignments.tsv
  ref=new_consensus.fa expectedinput=N expectedmembers=KEPT expectedmodels=M
  partitionbins=1024 paddingblacklist=none out=fresh_profiles.tsv
Omit the historical ledger/census/summary inputs. Full FASTA headers bind each
actual new builder assignment, including nonvoting small/orphan records.
One native maximum key per member/bin votes, with100extra votes for every
distinct eligible consensus bin/key; ties use the greatest key. The version4
resource fixes forward8/masked10 quarter-short routing and native4+12 encoding.
Queries count matching bin maxima, preserving a caller's global DDL exponent.
Use explicit paddingblacklist=none initially, or a complete256short/2048long
table plus paddingfit=scope later. The resource retains those words so query
and model filtering agree. Both512 and1024 total bins are supported explicitly.
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
if [ -f "${DIR}bbtools.jar" ]; then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=2g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.BuildNcrnaDdlProfiles "$@"
