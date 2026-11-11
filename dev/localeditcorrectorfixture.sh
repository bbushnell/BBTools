#!/bin/bash
# Generate the existing fixed CLI panel; do not forward JVM options to its positional main.
set -eo pipefail
if [[ $# != 1 || $1 = -* ]];then
 echo 'Usage: localeditcorrectorfixture.sh NEW_directory' >&2
 exit 2
fi
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed -ea -eoom
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.assemble.LocalEditCorrectorTest "$1"
