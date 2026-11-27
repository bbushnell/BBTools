#!/bin/bash
# Learn training-only RF00174 padding through the existing real caller/oracle.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]];then
	cat <<'HELP'
Usage: dev/buildncrnaddlpadding.sh caller=ordinary.config sample=training.tsv
  binding=worker_verified_execution.tsv
  first=0 count=128 expectedmodels=354 makelist=f fit=TRAINING
  out=fresh_pilot t=16 -Xmx8g -eoom

Calls the explicit training-sample slice once using the ordinary main-subfamily
caller. The existing exhaustive oracle verifies every current main model on
each real window. Count only the two flanks outside the envelope of ALL accepted
spans, separately: forward plain8 and middle-masked10 words, strict ambiguity
rejection before masking. No annotation-context shortcut or old-model envelopes.
No ordinary GFF or per-model oracle dump is written. Other caller controls stay
at the supplied flat config; inputs are replaced by the declared sample slice,
taxonomy is off, and no contig sampling or selected DDL/GC is allowed.

Outputs: .frequencies.tsv.gz (all positive raw-word counts), .inputs.tsv,
.sample.tsv (complete ordered sample), .binding.tsv (worker execution binding),
.models.fa (exact model binding), .summary.tsv, and with makelist=t,
.padding.tsv (256short/2048long occurrence-frequency prefixes; raw-word tie).
For the remainder, use first=128 count=3872 reuse=fresh_pilot makelist=t.
The pilot's frequencies and input audit are retained exactly once; overlapping
inputs, another sample, or different ordered consensus identities/bytes fail.
The reviewed experiment worker owns binding= creation and verification: exact
source/build, effective caller options and every loaded resource. The tool
compares that retained binding before importing counts; a final list requires
complete sample coverage. Supplying an unverified label is not that check.
This is resource preparation, not caller accuracy or the nine-cell comparison.
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
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=8g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError \
	-cp "${DIR}current/" prok.BuildNcrnaDdlPadding "$@"
