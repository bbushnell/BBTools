#!/bin/bash
# Finite source/strand/identity conservation controls; all files stay in out=.
set -eo pipefail
if [[ $# == 0 || $1 == --help ]];then
 echo 'Usage: dev/cobalaminconstructioninputstest.sh out=fresh_fixture_directory -Xmx2g -eoom'
 exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [[ -h $SCRIPT ]];do
 cd "$(dirname "$SCRIPT")"
 nextScript="$(readlink "$(basename "$SCRIPT")")"
 SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=2g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "${DIR}current/" test.prok.CobalaminConstructionInputsTest "$@"
