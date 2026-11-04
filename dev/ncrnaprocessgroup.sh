#!/usr/bin/env bash
# Internal NcrnaProcess startup protocol. Invoke through setsid, never directly.
# Stdout/stderr stay inherited. Parent removes the ID file before releasing work.
set -euo pipefail
identity=${1:?temporary_group_identity}
shift
[[ $# -gt 0 ]]
printf '%s\n' "$$" > "$identity"
IFS= read -r release
[[ "$release" = start ]]
exec "$@"
