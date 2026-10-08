#!/usr/bin/env bash
# Preserve declared two-layer input normalization before weight quantization.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
  echo 'Usage: cellnetaffineinput.sh in=explicit_affine.bbnet out=fresh_input_headers.bbnet x=validation.npy rows=1000 -Xmx4g'
  exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
if [[ -f $DIR/bbtools.jar ]]; then CP=$DIR/bbtools.jar; else CP=$DIR/current; fi
JAVA_BIN=$(command -v java)
source "$DIR/javasetup.sh"
ARGS=()
JVM_ARGS=()
for arg in "$@"; do
  case ${arg,,} in
    -xmx*|-xms*|-ea|-da|ea|da|eoom|-eoom|simd=t|simd=f) JVM_ARGS+=("$arg");;
    *) ARGS+=("$arg");;
  esac
done
parseJavaArgs --xmx=4g --xms=64m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec "$JAVA_BIN" $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" ml.CellNetAffineInput "${ARGS[@]}"
