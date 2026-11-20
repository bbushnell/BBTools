#!/usr/bin/env bash
# Run the finite static-report checks with assertions and a supplied classpath.
set -euo pipefail
test "$#" = 4
exec "$1" -ea --add-modules jdk.incubator.vector -Xmx128m -XX:ActiveProcessorCount=1 \
    -cp "$2" test.jgi.CallPeaksKNormalizationTest "$3" "$4"
