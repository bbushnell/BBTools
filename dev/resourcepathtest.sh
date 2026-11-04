#!/usr/bin/env bash
# Developer regression; requires exclusive access to the actual compiled runtime.
# Creates audit-prefixed fixtures; follow DataResourcePathTest's absence preflight.
set -eo pipefail
if [[ $# == 0 ]]; then
 echo 'Usage: resourcepathtest.sh root=/absolute/exclusive/runtime [baseline=t] [jar=t]'
 exit 0
fi
script=$(readlink -f "${BASH_SOURCE[0]}")
dir=$(cd "$(dirname "$script")/.." && pwd)
DIR=$dir
source "$dir/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS -cp "$dir/current${RESOURCE_TEST_CLASSPATH:+:$RESOURCE_TEST_CLASSPATH}" test.dna.DataResourcePathTest "$@"
