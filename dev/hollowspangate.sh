#!/bin/bash
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
if [[ $# == 0 || "$1" == --help || "$1" == -h ]]; then
  echo "Experimental post-hoc start veto; does not extend reads or correct bases."
  echo "in=saved_extensions.fa out=guarded.fa details=flags.tsv minus=k127.tsv nominal=k128.tsv plus=k129.tsv"
  echo "k=128 ratio=4 minsupport=3 (defaults). Profile rows must match extension order, IDs and seed anchors."
  echo "Veto when max(alternate depths)>=minsupport and >=ratio*nominal; truncate only the added suffix."
  echo "Requires fresh outputs; test=t runs boundary/overflow checks. Audit against truth separately."
  exit
fi
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.assemble.HollowSpanGate "$@"
