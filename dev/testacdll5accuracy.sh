#!/usr/bin/env bash
# Internal fixed-workload cardinality diagnostic; no biological inputs or outputs.
set -euo pipefail
usage(){
cat <<'HELP'
Usage: dev/testacdll5accuracy.sh run=t [-Xmx1g] [-Xms64m] [t=1]

Compare compressed and arithmetic-compressed five-bit cardinality trackers.
Ten million insertions, scanning 2048 buckets before and after each insertion.
This is an expensive full diagnostic, not a short smoke test.
Results go to stdout; divergences are reported, not converted into failure exits.
The existing diagnostic examines 2048 of the trackers' 2052 actual registers.

run=t            Explicitly run the fixed workload. No arguments, run=f or --help
                 prints this help without starting Java. There are no data inputs,
                 output-file flags, iteration limits or random-seed overrides.
t=1, threads=1   One diagnostic worker; other values are rejected. JVM service
                 threads are separate from the serial diagnostic loop.
-Xmx1g           Default maximum heap. Standard BBTools memory overrides apply.
-Xms64m          Default initial heap; setting only -Xmx also sets -Xms to that value.
eoom             Enable the JVM exit-on-out-of-memory flag.
simd=f           Default; simd=t is optional on a JVM with the Vector API.
                 The diagnostic's estimator updates are scalar.
library=DIR      Development override for the BBTools directory containing
                 javasetup.sh and memdetect.sh; default is this script's parent.
classpath=PATH   Development override: one class directory or JAR, with no implicit
                 fallback. Default is library/current. Use the matching built code.

Assertions are always enabled. The workload and checkpoint values are fixed.
HELP
}
if [[ $# == 0 ]]; then usage; exit 0; fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
CP=
RUN=false
JVM_ARGS=()
for arg in "$@"; do
 case "${arg,,}" in
  -h|--help|help|help=t|help=true) usage; exit 0;;
  run=t|run=true) RUN=true;;
  run=f|run=false) RUN=false;;
  t=1|threads=1) ;;
  t=*|threads=*) printf '%s\n' 'This diagnostic uses one worker; require t=1.' >&2; exit 2;;
  library=*) DIR=${arg#*=};;
  classpath=*) CP=${arg#*=};;
  -xmx*|-xms*) JVM_ARGS+=("$arg");;
  -ea|ea|-eoom|eoom|simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("${arg,,}");;
  -da|da) printf '%s\n' 'Assertions must remain enabled.' >&2; exit 2;;
  *) printf 'Unknown option: %s\n' "$arg" >&2; exit 2;;
 esac
done
if [[ $RUN != true ]]; then usage; exit 0; fi
[[ -f $DIR/javasetup.sh && -f $DIR/memdetect.sh ]] || { printf '%s\n' 'Missing BBTools Java setup scripts.' >&2; exit 2; }
if [[ -z $CP ]]; then CP=$DIR/current; fi
[[ -d $CP || -f $CP ]] || { printf 'Missing classpath: %s\n' "$CP" >&2; exit 2; }
: "${SHIFTER_RUNTIME:=}" "${EC2_HOME:=}" "${NERSC_HOST:=}"
source "$DIR/javasetup.sh"
# An explicit scalar default also avoids auto-detection failure under Java 8/set -e.
parseJavaArgs --xmx=1g --xms=64m --mode=fixed simd=f "${JVM_ARGS[@]}"
setEnvironment
read -r -a JAVA_FLAGS <<< "$EA $EOOM $SIMD $XMX $XMS $PROXY"
exec java "${JAVA_FLAGS[@]}" -ea -cp "$CP" test.cardinality.TestACDLL5Accuracy
