#!/usr/bin/env bash
# Build a plain frozen public-base or16S-release snapshot, never Brian's source tree.
set -euo pipefail
revision=${1:?exact_committed_revision}
mode=${2:?baseline_or_release}
[[ "$revision" =~ ^[0-9a-f]{40}$ ]]
[[ "$mode" = baseline || "$mode" = release ]]
code=/mnt/c/codex-lbl2/Raiden/workspace/BBTools_canonical_20260924
out=/home/bbushnell/.cache/raiden/prok16s_release_20261012/$revision
sha=/home/bbushnell/.local/bin/sha80
test ! -e "$out"
mkdir -p "$out/runtime" "$out/logs"
trap 'printf "%s\n" "$?" > "$out/exit_code.txt"' EXIT
git -C "$code" archive "$revision" > "$out/source.tar"
tar -xf "$out/source.tar" -C "$out/runtime"
find "$out/runtime/current" -name '*.class' -delete
javac -J-Xmx2g -J-XX:ActiveProcessorCount=2 -Xprefer:source --add-modules jdk.incubator.vector -cp "$out/runtime/current" "$out/runtime/current/prok/"*.java "$out/runtime/current/idaligner/"*.java "$out/runtime/current/ml/"*.java "$out/runtime/current/jgi/"*.java "$out/runtime/current/consensus/"*.java "$out/runtime/current/dna/"*.java "$out/runtime/current/shared/"*.java > "$out/logs/compile.stdout" 2> "$out/logs/compile.stderr"
cat "$out/logs/compile.stderr"
if [[ "$mode" = release ]];then
  bash "$out/runtime/dev/prok16s_release_test.sh" > "$out/logs/tests.stdout" 2> "$out/logs/tests.stderr"
  cat "$out/logs/tests.stdout" "$out/logs/tests.stderr"
fi
tar -cf "$out/runtime.tar" -C "$out" runtime
"$sha" "$out/runtime.tar" > "$out/RUNTIME.sha80"
printf '%s\n' "$revision" > "$out/SOURCE.txt"
printf 'PROK16S_RELEASE_BUILD_PASS\t%s\n' "$mode" > "$out/COMPLETE.txt"
cat "$out/COMPLETE.txt" "$out/RUNTIME.sha80"
