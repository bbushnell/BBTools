#!/bin/bash
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "${TADPOLE_POLICY_TEST_CLASSPATH:-$DIR/current/}" assemble.TadpoleAssemblyPolicyTest "$@"
