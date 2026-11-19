#!/usr/bin/env bash
# Fixed-donor generic comparison plus ordinary file-load/worker-copy checks.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Compare every native network output with saved GPU predictions or an exact net.
Usage: subnetgpuroundtriptest.sh net=export.bbnet vectors=raw.tsv \
 predictions=reference.tsv rows=1000 tolerance=0.00002 ordinary=t argmax=t
Use reference=other.bbnet instead of predictions= for exact native output bits.

vectors= starts #dims TABinputs TABtargets TAB0; targets match outputs, or two
stored value labels serve the existing four/six-head subnet case. Only raw input
columns enter the net. Predictions contain row_id followed by every output,
without a required header; blank/comment rows are skipped. Row IDs start at0.
No normalization is recomputed. Serialized preprocessing applies through CellNet.
forceJavaParseDouble=t supports exponent-formatted numeric fixtures.

outbits=checked.tsv optionally saves the actual first-pass Java outputs without
another inference pass. The plain TSV starts #prediction_bits TAB1 TABoutputs;
each data row contains its zero-based row ID and signed decimal IEEE754 float-bit
integers, decoded with Float.intBitsToFloat. This avoids formatting changing ties.
The final path is created only after every requested check and output close passes.
Existing output is never replaced. A failed run retains checked.tsv.partial as
unaccepted evidence; neither an old partial nor a final output is overwritten.

ordinary=t also checks the ordinary loader and the caller's copy(false) path.
With predictions=, copied forward parameters must be exact, every accelerated
path must satisfy the exported-reference tolerance, and a separate scalar
original/copy control must match exactly. Native path errors and first-maximum
differences are reported: vector lane sums need not repeat identical bits.
With reference=, native output bits remain exact. argmax=t reports first-maximum
changes against the reference; a small numeric error can change the selected
class. This is not a confidence-gated caller.
selftest=t runs the finite six-output/schema/load/copy/argmax oracle instead.
HELP
 exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}") || exit 1
DIR=$(cd "$(dirname "$SCRIPT")" && pwd) || exit 1
if [[ ! -f "$DIR/javasetup.sh" ]]; then DIR=$(dirname "$DIR"); fi
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=64m --mode=fixed "$@"
setEnvironment
ARGS=()
CLASS=test.ml.SubnetGpuRoundTripTest
for arg in "$@";do
 case "$arg" in
  selftest=t) CLASS=test.ml.SubnetGpuRoundTripFixtureTest ;;
  --xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-[xX][mM][xX]*|[xX][mM][xX]=*|-[xX][mM][sS]*|[xX][mM][sS]=*|-da|-ea|da|ea|ExitOnOutOfMemoryError|exitonoutofmemoryerror|eoom|-ExitOnOutOfMemoryError|-exitonoutofmemoryerror|-eoom) continue ;;
  *) ARGS+=("$arg") ;;
 esac
done
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" "$CLASS" "${ARGS[@]}"
