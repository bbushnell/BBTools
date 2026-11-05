#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
exec java -ea --add-modules jdk.incubator.vector -Xmx1g -cp "$root/current" test.jgi.Dedupe2RegressionTest "$@"
