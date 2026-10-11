#!/usr/bin/env bash
# Isolated developer regression; root must hold a disposable compiled runtime.
set -eo pipefail
if [[ $# == 0 ]]; then
 echo 'Usage: resourcepathtest.sh root=/absolute/disposable/runtime [baseline=t] [jar=t]'
 exit 0
fi
script=$(readlink -f "${BASH_SOURCE[0]}")
dir=$(cd "$(dirname "$script")/.." && pwd)
DIR=$dir
source "$dir/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS -cp "$dir/current${RESOURCE_TEST_CLASSPATH:+:$RESOURCE_TEST_CLASSPATH}" dna.DataResourcePathTest "$@"
