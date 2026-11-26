#!/usr/bin/env bash
# Direct positive-window CoveragePileup fixture and timing carrier.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"

app_args=()
jvm_args=()
for arg in "$@"; do
	case "$arg" in
		--xmx=*|--xms=*|--xss=*|--mem=*|--percent=*|--mode=*|\
		*[xX][mM][xX]*|*[xX][mM][sS]*|-[xX][sS][sS]*|[xX][sS][sS]*|\
		-da|-ea|da|ea|\
		ExitOnOutOfMemoryError|exitonoutofmemoryerror|eoom|-ExitOnOutOfMemoryError|-exitonoutofmemoryerror|-eoom|\
		simd|SIMD|simd=t|simd=true|simd=f|simd=false|nosimd)
			jvm_args+=("$arg")
			;;
		*)
			app_args+=("$arg")
			;;
	esac
done

parseJavaArgs --xmx=128m --xms=32m --mode=fixed -ea simd=f "${jvm_args[@]}"
setEnvironment

exec java -Dfile.encoding=UTF-8 -XX:ActiveProcessorCount=1 $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" test.jgi.CoveragePileupWindowTest "${app_args[@]}"
