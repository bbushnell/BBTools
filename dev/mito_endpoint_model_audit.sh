#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'mode=counts|support|test in=ASSIGNMENTS enrollment=ENROLLMENT consensus=FASTA outprefix=FRESH_PREFIX [corpus=WINDOWS model=5 end=3 expected=5575 edits=DECLARED_EDITS.tsv t=4]'
  echo 'Support tests one unedited model. Test mode requires a declared source-supported edit roster (model_index,model,end,trim,extend_bases); only its explicit alternatives are aligned, with unchanged-model controls. No packaged model is modified.'
  exit 0
fi
exec java -ea -Xmx8g -XX:ActiveProcessorCount=4 -XX:+ExitOnOutOfMemoryError --add-modules jdk.incubator.vector -cp "$root/current" test.prok.MitoEndpointModelAudit "$@"
