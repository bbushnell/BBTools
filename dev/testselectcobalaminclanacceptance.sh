#!/usr/bin/env bash
set -eo pipefail
if [[ $# == 0 || $1 == --help ]];then
 echo 'Usage: dev/testselectcobalaminclanacceptance.sh out=fresh_directory -Xmx1g -ea -eoom'
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
parseJavaArgs --xmx=1g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError \
 -cp "${DIR}current/" test.prok.SelectCobalaminClanAcceptanceTest "$@"
