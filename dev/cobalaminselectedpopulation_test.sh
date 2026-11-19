#!/bin/bash
# Metadata population conservation and source-coordinate collision fixtures.
set -eo pipefail
[[ $# -gt 0 ]] || { echo 'Supply out=fresh_fixture_directory' >&2;exit 1; }
pushd . > /dev/null
cd "$(dirname "${BASH_SOURCE[0]}")/.."
DIR="$(pwd)/"
popd > /dev/null
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=1g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "${DIR}current/" test.prok.CobalaminSelectedPopulationTest "$@"
