#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'Characterize the accepted mitochondrial construction audit without searches or membership changes.'
  echo 'in=AUDIT.tsv.gz source=ORIGINAL_CANDIDATES.fa.gz verified=SHA80_PATHS outprefix=NEW_PREFIX expected=37716 expectedcm=718 expectednohit=7050 expected16=4831 expected18=6586 expectedssu=11420'
  echo 'lengthbins=200,500,1000,1500,2000,3000,5000 t=1 -Xmx2g'
  echo 'Outputs: .nohit.tsv, .margins.tsv, exact-copy .nohit.fa.gz and unsplit .ssu.fa.gz including preference ties. Annotation-copy weighting, signed16S-minus18S margins, nearest-rank quantiles. Below-GA scores are NA.'
  exit 0
fi
DIR=$root
set +u
source "$DIR/javasetup.sh"
if ! parseJavaArgs --xmx=2g --mode=fixed "$@";then exit 1;fi
setEnvironment
set -u
read -r -a java_options <<< "$EA $EOOM $SIMD $XMX $XMS $XSS"
command=(java "${java_options[@]}" -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=4 -cp "$root/current" prok.NcrnaMitoCharacterization "$@")
printf '%q ' "${command[@]}" >&2
printf '\n' >&2
exec "${command[@]}"
