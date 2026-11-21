#!/usr/bin/env bash
# Run the finite public smoothing checks with the supplied compiled classpath.
set -euo pipefail
test "$#" = 3
exec "$1" -ea --add-modules jdk.incubator.vector -Xmx128m -XX:ActiveProcessorCount=1 \
    -cp "$2" test.jgi.CallPeaksSmoothingTest "$3"
