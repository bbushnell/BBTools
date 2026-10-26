#!/bin/bash
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
if [[ $# == 0 || "$1" == --help || "$1" == -h ]]; then
  echo "Audits RandomReadsMG Q20 evidence or grades right-only extensions of 200bp Shred pieces."
  echo "Audit: hollowextensiontruth.sh mode=audit ref=original.fa mutant=mutant.fa shreds=shreds.fa in=reads.fq.gz out=audit.tsv"
  echo "Grade: hollowextensiontruth.sh mode=grade ref=original.fa shreds=shreds.fa in=extended.fa maxextension=20 label=solid20 out=summary.tsv details=errors.tsv"
  echo "Reports only added-base errors; verifies original prefixes and complete unique shred IDs. Outputs must be new."
  exit
fi
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.assemble.HollowExtensionTruth "$@"
