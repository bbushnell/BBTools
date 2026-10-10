#!/usr/bin/env bash
# Experimental full-top50 positional assignment; production remains unchanged.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
 cat <<'HELP'
Usage: hbmcompetitiveassign.sh config=FILE in=QUERY_FASTA insha80=PIN
  expected=RECORD_COUNT out=NEW_DIRECTORY
Required model inputs, each with its KEYsha80=PIN argument:
  consensus, identity, hbm, provenance, cores, background, sidecar, kmers
Uses the full F4 top50, positional logodds beta.01 clipped[-4,11], gap4,
then identity>=40.612846% and paired-core coverage>=.8 (float32).
Chooses the highest score passing both gates, with ASCII family-ID ties.
Every protein gets an assignment or an explicit NO_ELIGIBLE_FAMILY row.
Only malformed edge-stop markers are repaired; internal or unsupported residues
fail the run. Raw inputs remain intact, and repairs are recorded.
Runs one worker per query shard. selftest=t runs the input/top50 fixtures.
familylimit=0 uses all families. A positive limit restricts candidates to that
many leading roster entries before selecting top50, retaining the full sidecar
and its F4 scores. This is an experimental old-roster-only control.
HELP
 exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
CP=$DIR/current
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
ARGS=(); JVM_ARGS=(); main=prot.HbmCompetitiveAssign
for arg in "$@"; do
 case ${arg,,} in
  selftest=t) main=prot.HbmCompetitiveAssignTest;;
  --xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
  -ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
  simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("$arg");;
  *) ARGS+=("$arg");;
 esac
done
parseJavaArgs --xmx=8g --xms=256m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -XX:ActiveProcessorCount=1 -cp "$CP" "$main" "${ARGS[@]}"
