#!/usr/bin/env bash
# Selected retained correctness build FIRST; shared timing helper SECOND; dependencies LAST.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
selected_build=; timing_build=
for arg in "$@"; do
 case "$arg" in classes=*) selected_build=${arg#classes=};; timingclasses=*) timing_build=${arg#timingclasses=};; esac
done
[[ $selected_build = /* && $timing_build = /* && -d $selected_build && -d $timing_build ]]
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed -ea simd=f "$@"
setEnvironment
exec java -Dfile.encoding=UTF-8 -XX:ActiveProcessorCount=1 $EA $EOOM $SIMD $XMX $XMS -cp "$selected_build:$timing_build:$DIR/current/" test.jgi.BBMergeErrorFilterTiming "$@"
