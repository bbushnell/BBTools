#!/usr/bin/env bash
# Paired native NPY value-output evaluation; explicit preprocessing mode.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
  echo 'Usage: cellnetweightbitseval.sh reference=32.bbnet candidate=24.bbnet x=X.npy y=Y.npy rows=400000 originalmse=0.00267 out=fresh.tsv t=16 store=bfloat16'
  echo 'Self-test: cellnetweightbitseval.sh selftest=t -Xmx256m'
  echo 'Candidate precision comes from its critical #weightbits 18 or24 header; the reference remains the matched32-bit net.'
  echo 'Optional limitmse= sets an absolute nonnegative MSE bar (D249 composite:0.00254). Otherwise the bar is1.005*originalmse.'
  echo 'originalmse= always records the actual original measurement. Explicit limits label the verdict within_absolute_bar; partial runs remain CANARY_ONLY.'
  exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
if [[ -f $DIR/bbtools.jar ]]; then CP=$DIR/bbtools.jar; else CP=$DIR/current; fi
JAVA_BIN=$(command -v java)
source "$DIR/javasetup.sh"
CLASS=ml.CellNetWeightBitsEval
ARGS=()
JVM_ARGS=()
for arg in "$@"; do
  case ${arg,,} in
    selftest=t|selftest=true) CLASS=ml.CellNetWeightBitsEvalTest;;
    -xmx*|-xms*|-ea|-da|ea|da|eoom|-eoom|simd=t|simd=f) JVM_ARGS+=("$arg");;
    *) ARGS+=("$arg");;
  esac
done
parseJavaArgs --xmx=4g --xms=64m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec "$JAVA_BIN" $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" "$CLASS" "${ARGS[@]}"
