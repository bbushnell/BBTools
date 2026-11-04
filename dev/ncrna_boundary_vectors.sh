#!/usr/bin/env bash
# Native shared-boundary vectors; offsets, tables and family policy are explicit.
usage(){
cat <<'EOF'
NcrnaBoundaryVectorGen: labeled endpoint vectors from annotated, flanked RNA.

Usage: bash dev/ncrna_boundary_vectors.sh family=rnasep fasta=contexts.fa \
         tablestart=start.tsv tablestop=stop.tsv out=vectors.tsv.gz

Required: family=, fasta=, tablestart=, tablestop=, out=.
Headers need lflank=, rflank= and contig_gc=; GC is the source contig's fraction.
Table endpoint types and word geometry must match the intended caller.

boundaryfeatures=v1  Feature layout: v1, v2 or v3. Extra inputs for v2/v3 are
                     described in docs/guides/NcrnaBoundaryVectorGenGuide.txt.
meanlen=0            Nonpositive selects the upper median consensus length.
leaveoneout=t        Subtract this locus's true-position count before smoothing.
stoponly=f           Emit both ends by default; true emits only stop candidates.
countonly=f          Count candidates without writing vectors; needs statsout=.
statsout=            Optional plain TSV of counts and skipped-candidate reasons.
locusout=            Optional locus/assembly/model/start_rows/stop_rows TSV.
                     Requires both ends and actual vectors; incompatible with
                     stoponly=t and countonly=t. Use a fresh destination.
singlemodel=f        With true, require exactly one consensus and assign model 0
                     without shortlisting. Alignment features are still computed.
                     Mutually exclusive with fixedmodel=.
r58boundarystartrange=  For family=r58 only: inclusive minimum,maximum offsets.
r58boundarystoprange=   Each range includes zero, stays within -256..256, and
                       contains at most 257 positions. Omission keeps defaults.
zl= / pigz= / bgzip= / unpigz=  Native compression controls. Compressed suffixes
                     apply to vector/locus outputs; statsout remains plain TSV.

This low-level launcher currently fixes a 4 GB heap and one active processor
and requires the JVM vector module. It does not accept config= or JVM flags.
It generates vectors; it does not train, select or install boundary networks.
EOF
}
if [[ $# -eq 0 || "$1" = -h || "$1" = --help ]]; then
  usage
  exit 0
fi
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea -Xmx4g -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector -cp "$root/current" prok.NcrnaBoundaryVectorGen "$@"
