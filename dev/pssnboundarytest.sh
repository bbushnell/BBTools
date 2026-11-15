#!/usr/bin/env bash
# Bounded ProcessStackedSitesNormalized baseline/fixed boundary fixture.
set -eo pipefail
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
outdir=
args=()
for arg in "$@"; do
	case "$arg" in
	outdir=*) outdir=${arg#outdir=}; args+=("$arg");;
	phase=baseline|phase=fixed) args+=("$arg");;
	child=*) args+=("$arg");;
	*) echo 'Use outdir=directory [phase=baseline|fixed] [child=name]' >&2; exit 1;;
	esac
done
[[ -n $outdir && $outdir = /* && ! -e $outdir && -d ${outdir%/*} ]]
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed -ea -eoom simd=f
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS -XX:ActiveProcessorCount=1 \
	-Dpssn.wrapper="$SCRIPT" -cp "$DIR/current" \
	test.pacbio.ProcessStackedSitesNormalizedBoundaryTest "${args[@]}"
