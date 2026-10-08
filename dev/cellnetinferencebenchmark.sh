#!/bin/bash
# Compare native inference representations on one subnet and its real input row.
set -o pipefail
if [[ $# -eq 0 || "$1" = -h || "$1" = --help ]]; then
	echo 'Usage: cellnetinferencebenchmark.sh bundle=<bbnets> subnet=<id> in=<subnet vectors> seconds=10 -Xmx4g
Reports source/dense outputs, real SIMD flags, rows/s and MACs/s. Each timing lasts >=10s.
selftest=t checks dense-copy arithmetic, isolation and dispatch without running a benchmark.'
	exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}") || exit 1
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd) || exit 1
if [[ -f "$DIR/bbtools.jar" ]]; then CP=$DIR/bbtools.jar; else CP=$DIR/current/; fi
source "$DIR/javasetup.sh" || exit 1
ARGS=(); JVM_ARGS=()
for arg in "$@"; do
	case "${arg,,}" in
	--xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx[0-9]*|xmx[0-9]*|-xmx=*|xmx=*|\
	-xms[0-9]*|xms[0-9]*|-xms=*|xms=*|\
	-ea|-da|ea|da|-eoom|eoom|-exitonoutofmemoryerror|exitonoutofmemoryerror|\
	simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("${arg,,}");;
	*) ARGS+=("$arg");;
	esac
done
parseJavaArgs --xmx=4g --xms=64m --mode=fixed "${JVM_ARGS[@]}" || exit 1
setEnvironment || exit 1
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" ml.CellNetInferenceBenchmark "${ARGS[@]}"
