#!/usr/bin/env bash
# Internal owned-session launcher for NcrnaChildProcess; never submits a job.
set -euo pipefail
if [[ "${1-}" != --session ]];then
  exec setsid --wait bash "$0" --session "$@"
fi
shift
pidfile=${1:?fresh_process_identity_file}
nonce=${2:?parent_nonce}
shift 2
[[ $# -gt 0 ]]
raw_pgid=$(ps -o pgid= -p "$$")
pgid=${raw_pgid// /}
[[ "$pgid" = "$$" ]]
set -C
printf '%s\t%s\n' "$nonce" "$$" > "$pidfile.pending"
test ! -e "$pidfile"
mv -- "$pidfile.pending" "$pidfile"
# No workload starts until Java has recorded cleanup ownership. A cancelled or
# abandoned startup must not leave a waiting launcher indefinitely.
for attempt in {1..1000};do
  [[ ! -e "$pidfile.cancelled" ]] || exit 143
  if [[ -f "$pidfile.ready" ]];then
    read -r ready_nonce < "$pidfile.ready"
    [[ "$ready_nonce" = "$nonce" ]]
    [[ ! -e "$pidfile.cancelled" ]] || exit 143
    exec "$@"
  fi
  sleep 0.01
done
printf 'Parent did not confirm child ownership\n' >&2
exit 70
