#!/bin/bash
# Deduplicates a v2 junction trace and labels it against a known source reference.
set -eo pipefail
if [[ $# -eq 0 ]]; then
    echo 'Usage: livebranchcorpus.sh in=trace.tsv ref=source.fa out=NEW_directory circular=f'
    echo 'Writes nodes.tsv (all labels), vectors.tsv (known grey/control rows), rows.tsv, census.tsv and complete.txt.'
    echo 'Use circular=t only when every supplied reference record is explicitly circular.'
    echo 'controlfraction=1 keeps all controls (default); 0.1 caps controls at 10% of emitted rows.'
    echo 'controlseed=1 ranks known controls by canonical sequence only. Use a frozen per-genome seed.'
    echo 'All known grey nodes and the full audit/census are retained. Do not sample validation/test exports.'
    exit 0
fi
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=2g --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.LiveBranchCorpus "$@"
