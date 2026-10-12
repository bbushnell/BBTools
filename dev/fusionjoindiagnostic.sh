#!/bin/bash
# Development-only fusion depth census; in=reads trace=log ref=fa assembly=fa out=tsv k=N.
# Optional outvectors=tsv writes fusion_join_v1 (107 numeric features after 3 ID fields).
# Optional outhist=tsv preserves all depth bins, including explicit empty regions.
# IDs are local to one trace: keep genome/trace provenance when pooling these files.
# Labels and final/reference hits exist only in out=, never in numeric features.
# Use query_k==phase_k for the same-K view; selected pairs are not all candidates.
# plan in=genomes.tsv out=samples.tsv expected=1000 seed=350194 randomizes recipes.
# label trace=joins.trace.tsv ref=reference.fa out=labels.tsv labels collected pairs.
# Label mode runs after assembly, without any read-count table; verify=t checks
# indexed matches against exhaustive matching for small validation fixtures.
# entropy in=samples.tsv root=corpus_root out=new_directory appends tip entropy
# as fusion_join_v2 (109 inputs) and writes group-disjoint native trainer tables.
# entropytest runs the native sequence-complexity fixtures without a read table.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=8g --xms=256m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.FusionJoinDiagnostic "$@"
