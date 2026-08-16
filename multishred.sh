#!/bin/bash
#multishred in=<directory> out=<file> length=<number>

usage(){
	echo "
Written by Brian Bushnell and Eru
Last modified August 12, 2026

Description:  Shreds all sequences in a directory into one output file,
using a single JVM invocation.  Equivalent to running shred.sh on every
file individually but orders of magnitude faster for large directories.

Usage:  multishred.sh indir=<directory> out=<file> length=<number>

File Parameters:
indir=<dir>     Directory containing input sequence files.
out=<file>      Destination of output shreds.
ext=fna.gz      Extension filter for input files (default fna.gz).

All other parameters are passed through to the shredding engine.
Run shred.sh with no arguments for the full parameter list.
"
	exit 0
}

pushd . > /dev/null
DIR="$(cd "$(dirname "$0")" && pwd)"
popd > /dev/null

INDIR=""
EXT="fna.gz"
ARGS=()
for a in "$@"; do
	case "$a" in
		indir=*) INDIR="${a#indir=}" ;;
		ext=*)   EXT="${a#ext=}" ;;
		in=*)    echo "Use indir= instead of in= with multishred.sh" >&2; exit 1 ;;
		*)       ARGS+=("$a") ;;
	esac
done

if [ -z "$INDIR" ] || [ ! -d "$INDIR" ]; then
	usage
fi

shopt -s nullglob
FILES=("$INDIR"/*."$EXT")
N=${#FILES[@]}
if [ "$N" -eq 0 ]; then
	echo "No *.$EXT files found in $INDIR" >&2
	exit 1
fi
echo "Shredding $N files from $INDIR (*.$EXT)" >&2

# Bash array expansion above avoids invoking `ls`/exec on the full file list,
# which silently truncates ("Argument list too long") once a directory holds
# enough files to exceed ARG_MAX (observed: fine at 9,734 files, broken at
# 35,395 -- caught 2026-08-16 when the count check exited 1 with no `set -e`
# in the caller, so the bacteria shred step silently no-op'd).
(
	for f in "${FILES[@]}"; do
		[ -f "$f" ] || continue
		case "$f" in
			*.gz) unpigz -c "$f" 2>/dev/null || zcat "$f" ;;
			*)    cat "$f" ;;
		esac
	done
) | "$DIR/shred.sh" in=stdin.fa "${ARGS[@]}"
