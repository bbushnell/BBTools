#!/usr/bin/env bash
# Finite native process-lifecycle fixture; not a training entry point.
set -euo pipefail
mode=
payload=
childfile=
control=
denydir=
for arg in "$@";do
  case "$arg" in
    mode=*) mode=${arg#mode=};;
    payload=*) payload=${arg#payload=};;
    childfile=*) childfile=${arg#childfile=};;
    control=*) control=${arg#control=};;
    denydir=*) denydir=${arg#denydir=};;
    t=*|-Xmx*|-Xms*|-eoom) ;;
    *) exit 64;;
  esac
done
case "$mode" in
  echo)
    printf '%s\n' "$payload"
    cat
    printf 'stderr captured\n' >&2
    ;;
  cwd) pwd -P;;
  sentinel) [[ -n "$childfile" ]];printf 'workload executed\n' > "$childfile";;
  fail) printf 'expected child failure\n' >&2;exit 7;;
  failcleanup)
    [[ -d "$denydir" ]]
    chmod a-w -- "$denydir"
    exit 7
    ;;
  readwait)
    [[ -n "$control" ]]
    for attempt in {1..3000};do
      [[ ! -e "$control.read" ]] || break
      sleep 0.01
    done
    [[ -e "$control.read" ]]
    cat > /dev/null
    touch -- "$control.consumed"
    for attempt in {1..3000};do
      [[ ! -e "$control.exit" ]] || exit 0
      sleep 0.01
    done
    exit 70
    ;;
  ignoreterm|holdstdin)
    [[ -n "$childfile" ]]
    trap '' TERM
    sleep 300 <&0 &
    child=$!
    printf '%s\n' "$child" > "$childfile"
    if [[ "$mode" = ignoreterm ]];then wait "$child";fi
    ;;
  *) exit 64;;
esac
