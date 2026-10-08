#!/usr/bin/env bash
# Convert only neural edge precision; preserve bias, normalization and topology.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
  echo 'Usage: cellnetweightbits.sh in=model.bbnet out=fresh.bbnet bits=18|24|32 -Xmx2g'
  echo '       cellnetweightbits.sh selftest=t -Xmx256m'
  echo '18/24-bit A48 truncates the low fourteen/eight edge-weight bits. Biases and normalization remain float32.'
  exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
if [[ -f $DIR/bbtools.jar ]]; then CP=$DIR/bbtools.jar; else CP=$DIR/current; fi
JAVA_BIN=$(command -v java)
source "$DIR/javasetup.sh"
CLASS=ml.CellNetWeightBits
ARGS=(); JVM_ARGS=(); SELFTEST=''
for arg in "$@"; do
  case ${arg,,} in
    selftest=*)
      [[ -z $SELFTEST ]] || exit 2
      SELFTEST=${arg#*=}
      case ${SELFTEST,,} in t|true) CLASS=ml.CellNetWeightBitsTest;; f|false) ;; *) exit 2;; esac;;
    -xmx*|-xms*|-ea|-da|ea|da|eoom|-eoom|simd=t|simd=f) JVM_ARGS+=("$arg");;
    *) ARGS+=("$arg");;
  esac
done
parseJavaArgs --xmx=2g --xms=64m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec "$JAVA_BIN" $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" "$CLASS" "${ARGS[@]}"
