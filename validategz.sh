#!/bin/bash
# validategz.sh - Report gzip files that are empty (0 bytes) or fail gzip -t.
# Written by Amber, August 16, 2026.
#
# Usage:  validategz.sh <dir-or-file> [more dirs/files ...] [-jN]
#   Scans each argument (recursively for directories) for *.gz files and, for
#   each, reports it if it is 0 bytes (EMPTY) or fails a gzip integrity test
#   (CORRUPT).  Prints one line per bad file plus a final summary; exit status
#   is nonzero if any bad file was found (so it can gate a pipeline).
#     -jN   run N gzip tests in parallel (default 4).
#
# This is the compression-level companion to validatepairs.sh, which checks the
# decompressed CONTENT (partial GFF records, fna/gff header mismatch).

usage(){
	echo "Usage: validategz.sh <dir-or-file> [more ...] [-jN]" >&2
	echo "  Reports *.gz files that are 0 bytes (EMPTY) or fail gzip -t (CORRUPT)." >&2
}

if [ -z "$1" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	usage
	exit 1
fi

JOBS=4
TARGETS=()
for arg in "$@"; do
	case "$arg" in
		-j*) JOBS="${arg#-j}";;
		*)   TARGETS+=("$arg");;
	esac
done
if [ "${#TARGETS[@]}" -eq 0 ]; then usage; exit 1; fi

# One-file test used by the parallel workers: prints "REASON<TAB>path" for bad files.
checkone(){
	f="$1"
	if [ ! -s "$f" ]; then
		echo -e "EMPTY\t$f"
	elif ! gzip -t "$f" 2>/dev/null; then
		echo -e "CORRUPT\t$f"
	fi
}
export -f checkone

# Collect the file list (recursive for dirs, literal for files), then test in parallel.
# NOTE: workers APPEND (>>) to $tmp, never a shared truncating redirect (>).  Each worker
# prints at most one short line per bad file; a single write under PIPE_BUF to an O_APPEND
# fd is atomic, so parallel workers cannot interleave or DROP lines.  A shared '>' redirect
# would race and silently lose failures (a real gzip-sweep bug this tool exists to avoid).
tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT
: > "$tmp"
for tgt in "${TARGETS[@]}"; do
	if [ -d "$tgt" ]; then
		find "$tgt" -type f -name '*.gz' -print0
	elif [ -f "$tgt" ]; then
		printf '%s\0' "$tgt"
	else
		echo "Warning: not found: $tgt" >&2
	fi
done | TMPOUT="$tmp" nice xargs -0 -P "$JOBS" -I{} bash -c 'checkone "$@" >> "$TMPOUT"' _ {}

total_bad=$(wc -l < "$tmp")
empty=$(grep -c '^EMPTY' "$tmp" 2>/dev/null || echo 0)
corrupt=$(grep -c '^CORRUPT' "$tmp" 2>/dev/null || echo 0)

if [ "$total_bad" -gt 0 ]; then
	cat "$tmp"
fi
echo "----" >&2
echo "Bad gzip files: $total_bad  (EMPTY: $empty, CORRUPT: $corrupt)" >&2

[ "$total_bad" -eq 0 ]
