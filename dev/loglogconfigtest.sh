#!/bin/bash
# Compare config expansion, CLI precedence, and redirected LogLog diagnostics.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
work=$(mktemp -d /tmp/loglog-config-test.XXXXXX)
printf 'LOGLOG_CONFIG_TEST_DIR %s\n' "$work"
common=(synth=t reads=200 t=1 seed=36 seed2=57 buckets=2048 histcolumns=3 overwrite=f)
for k in 31 62; do
  bash "$root/loglog.sh" -Xmx128m "${common[@]}" k="$k" \
    out="$work/direct$k.card" khist="$work/direct$k.tsv" > "$work/direct$k.stdout" 2> "$work/direct$k.stderr"
  bash "$root/loglog.sh" -Xmx128m config="$root/testdata/loglog_config_test.config" k="$k" \
    out="$work/config$k.card" khist="$work/config$k.tsv" outstream="$work/config$k.diagnostics" \
    > "$work/config$k.stdout" 2> "$work/config$k.stderr"
  cmp "$work/direct$k.card" "$work/config$k.card"
  cmp "$work/direct$k.tsv" "$work/config$k.tsv"
  grep -q '^Cardinality:' "$work/config$k.diagnostics"
  if grep -q '^Cardinality:' "$work/config$k.stderr"; then exit 1; fi
done
# A config after a conflicting CLI option follows normal in-order precedence.
bash "$root/loglog.sh" -Xmx128m k=31 config="$root/testdata/loglog_config_test.config" \
  out="$work/configlast.card" khist="$work/configlast.tsv" > "$work/configlast.stdout" 2> "$work/configlast.stderr"
cmp "$work/direct62.card" "$work/configlast.card"
cmp "$work/direct62.tsv" "$work/configlast.tsv"
if bash "$root/loglog.sh" -Xmx128m config="$work/missing.config" > "$work/missing.stdout" 2> "$work/missing.stderr"; then
  printf 'Missing config unexpectedly succeeded.\n' >&2
  exit 1
fi
printf 'LOGLOG_CONFIG_TEST_PASS\n'
