#!/usr/bin/env bash
# After full prok/test.prok compilation on the sole Dori build:
# unit FRESH_DIR
# run FRESH_PARENT INPUT_PATH_LIST 1|16 off|on  (four independent invocations)
# verify PARENT INPUT_PATH_LIST             (after all four finish)
# Caller resources must be installed. INPUT_PATH_LIST has >=16 frozen genome paths.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
mode=${1:?unit|run|verify}
if [[ "$mode" = unit ]];then
  exec java -ea --add-modules jdk.incubator.vector -XX:ActiveProcessorCount=1 -XX:+ExitOnOutOfMemoryError -Xmx2g -cp "$root/current" test.prok.CallGenesFileCostTest "outdir=${2:?fresh_dir}"
fi
parent=${2:?parent_dir}
list=${3:?input_path_list}
mapfile -t inputs < "$list"
test "${#inputs[@]}" -ge 16
for path in "${inputs[@]}";do test -f "$path";done
if [[ "$mode" = run ]];then
  threads=${4:?1|16}
  diagnostic=${5:?off|on}
  case "$threads" in 1|16) ;; *) exit 2;; esac
  case "$diagnostic" in off|on) ;; *) exit 2;; esac
  out=$parent/t${threads}_${diagnostic}
  test ! -e "$out";mkdir -p "$out"
  trap 'printf "%s\n" "$?" > "$out/exit_code.txt"' EXIT
  # JVM option tokenization is deliberately explicit in this development fixture.
  [[ "$out" != *[[:space:]\"]* ]]
  unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS
  if [[ "$diagnostic" = on ]];then
    export JDK_JAVA_OPTIONS="-Dbbtools.diagnostic.callGenesFileCostOut=$out/cost -Dbbtools.diagnostic.ncrnaCpuOut=$out/cpu.tsv"
  fi
  printf '%s\n' "${JDK_JAVA_OPTIONS:-}" > "$out/jvm_options.txt"
  printf -v input_list '%s,' "${inputs[@]}"
  printf 'taxonomy=f\ncallcds=t\neuk=f\nncrna=t\ntmrna=t\nsixs=t\ncobalamin=t\ncm=f\now=f\nordered=t\nextendedstats=t\nncrnafamily=cobalamin\nncrnaindexk=7\nncrnaindextopn=100\nncrnaindexminhits=1\nncrnaadaptiveminhits=f\nncrnaidpass=0.65\nncrnaidborderline=0.65\nt=%s\n' "$threads" > "$out/common.config"
  command=(bash "$root/callgenes.sh" "config=$out/common.config" "in=${input_list%,}" "out=$out/calls.gff" -Xmx8g)
  printf '%q ' "${command[@]}" > "$out/COMMAND.txt";printf '\n' >> "$out/COMMAND.txt"
  "${command[@]}" > "$out/caller.stdout" 2> "$out/caller.stderr"
  exit 0
fi
test "$mode" = verify
for threads in 1 16;do
  for diagnostic in off on;do test "$(cat "$parent/t${threads}_${diagnostic}/exit_code.txt")" = 0;done
  cmp "$parent/t${threads}_off/calls.gff" "$parent/t${threads}_on/calls.gff"
  cmp <(grep '^ncRNA Alignments (' "$parent/t${threads}_off/caller.stderr") <(grep '^ncRNA Alignments (' "$parent/t${threads}_on/caller.stderr")
  test ! -e "$parent/t${threads}_off/cost.inputs.tsv"
  test "$(wc -l < "$parent/t${threads}_on/cost.inputs.tsv")" = "$((${#inputs[@]}+2))"
  grep -q $'^B12_RESOURCES\tCOMPLETE\t' "$parent/t${threads}_on/cost.startup.tsv"
  grep -q $'^SHARED_SEED_INDEX\tCOMPLETE\t' "$parent/t${threads}_on/cost.startup.tsv"
done
cmp "$parent/t1_on/calls.gff" "$parent/t16_on/calls.gff"
cmp <(cut -f1-9 "$parent/t1_on/cost.families.tsv") <(cut -f1-9 "$parent/t16_on/cost.families.tsv")
printf 'CALLGENES_FILE_COST_CALLER_PARITY_PASS\tinputs=%s\tt1_and_t16\tdiagnostics_off_on\n' "${#inputs[@]}"
