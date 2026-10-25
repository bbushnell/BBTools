#!/usr/bin/env bash
# A new diagnostic class on the original accepted runtime/resources; no caller changes.
set -euo pipefail
overlay=${1:?Diagnostic_classes_only_directory}
original=${2:?Original_candidate_runtime}
shift 2
cd "$original"
exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx2g -cp "$overlay:$original/current" prok.MitoSsuShortlistControl "$@"
