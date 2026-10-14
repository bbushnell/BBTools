#!/bin/bash
# Development-only: collect raw eligible unitig connections using one live table.
# task=collect in=reads.fq.gz k=64 out=candidates.tsv rawcontigs=unitigs.fa t=8
# Collect fixes mode=contig, pop=f, sweeplen=0, processcontigs=f,
# retainshortcontigs=t, trimends=0, trimcircular=f, hashkmers=f. Native seed,
# minextension, coverage, washing and mce settings are otherwise retained.
# This is not a production assembly recipe; no candidates are joined.
# task=label in=candidates.tsv ref=source.fa out=labels.tsv circular=f
# Labeling runs separately; circular=t applies to ALL reference records.
# Local circular Kmer truth wraps origins. Initial bounded pair-context labels
# do not resolve origin-spanning pairs and never call those absent pairs false.
# context=200 sets retained bases outside the K-1 overlap. Outputs must be new.
# Optional collect flags: features=vectors.tsv controls=controls.tsv
# controlcount=-1 (default) samples floor(candidates/9) eligible internal cuts;
# controlcount=N overrides that count; controlseed=1 fixes the sample seed.
# controls= requires features=. Both remain diagnostic and do not alter contigs.
# Feature rows have kind=candidate|control, event, then 27 shared normalized inputs.
# Comment metadata records live-table moments and pre-label sampling population.
# Label the distinct control file with task=label in=controls.tsv ref=source.fa
# out=control_labels.tsv. Control labels are never accepted as merger labels.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=4g --xms=64m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.SameKJoinDiagnostic "$@"
