#!/usr/bin/env bash
# Construct reviewable comparison candidates; no model filter is changed.
set -euo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Usage: proteinkmercandidates.sh real=/accepted/real_count_root \
       shadow=/accepted/shadow_count_root legacy=/old/full.faa out=/new/candidates

Create candidate FASTAs from independently accepted K5/K6 counts.
Each root contains real_k5/real_k6 or shadow_k5/shadow_k6 phases with owner
acceptance, successful native/outer/watcher receipts, and frequency outputs.
Real-only, shadow-only and alternating-rank mixed lists each contain2048 words.
Full unions are separate different-size controls. Legacy K5 is re-emitted.
k6sizes=2048          Default: nine accepted-top-list candidates.
k6sizes=2049         Add real/shadow/mixed lists of exactly2049 words.
k6sizes=2048,8192,40960  Any unique positive sizes, bounded by available K6 words.
                       Retain the nine default controls; add three lists per
                       other size. Read full counts once per population when
                       needed, retain only the largest requested rank set, and
                       verify exact top2048 parity. Sizes are emitted ascending.
                       Outputs above200KiB raw size use BGZF9.
                       Requires external bgzip with -l 9 and -@ 1 support;
                       fails rather than falling back to native level5 compression.
No model is modified. Input verification is structural; no hashes are computed.
Alternatively replace legacy/k6sizes with explicit normalized treatments:
arms=raw:shadow:16384,raw:weighted:16384,reduced:shadow:8192,reduced:shadow:16384,reduced:weighted:16384
alphabet=all:A/C/DE/FWY/G/IV/LM/NH/P/Q/RK/ST
negativeweight=4 positiveweight=1
Each population frequency is normalized by its total valid K6 windows.
Reduced words aggregate ALL raw counts before ranking. Rank ties use packed
word order. Only K6 count phases are consumed; output includes exact selected
counts/scores and an explicit alphabet/space/weight manifest. The weighting
and size choices are parameters, not a separate script per comparison.
Output must be a new absolute directory. t=1; default heap1g.

rankedexport=/accepted/candidate_run out=/new/ranked_lists
Export existing raw shadow/weighted K6 rankings as BGZF9 FASTAs in score order.
Reads accepted manifests, score tables and FASTAs only; no recount or reselection.
Supports one list per selection, each with up to16384 words. Validates ranks,
scores, canonical words and exact FASTA membership. Do not combine with selection
parameters. Preserve first N whole FASTA records to obtain a smaller raw list.

library=DIR          Optional BBTools runtime; default parent of dev directory.
modulejar=JAR        Development module before runtime classes; testjar alias.
dependencyjar=JARS   Colon-separated dependency jars between module and runtime.
test=t testdir=DIR   Focused fixture in a fresh absolute directory (Dori only).
HELP
 exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
LIBRARY=$(cd "$(dirname "$SCRIPT")/.." && pwd)
MODULE= DEPENDENCIES= TEST=false TESTDIR=
ARGS=();JVM_ARGS=()
for arg in "$@";do
 case ${arg,,} in
  library=*) LIBRARY=${arg#*=};;
  modulejar=*|testjar=*) [[ -z $MODULE ]];MODULE=${arg#*=};;
  dependencyjar=*) [[ -z $DEPENDENCIES ]];DEPENDENCIES=${arg#*=};;
  test=t|test=true) TEST=true;;
  test=f|test=false) TEST=false;;
  testdir=*) [[ -z $TESTDIR ]];TESTDIR=${arg#*=};;
  -da|da|-da:*) printf '%s\n' 'Assertions must remain enabled' >&2;exit 2;;
  -xmx*|-xms*|-ea|ea|-eoom|eoom|simd|simd=t|simd=f|nosimd) JVM_ARGS+=("$arg");;
  *) ARGS+=("$arg");;
 esac
done
CLASS=ddl.ProteinKmerCandidates
if [[ $TEST == true ]];then
 [[ ${#ARGS[@]} == 0 && $TESTDIR == /* && ! -e $TESTDIR && ! -L $TESTDIR ]]
 CLASS=test.ddl.ProteinKmerCandidatesTest;ARGS=("$TESTDIR")
else
 [[ -z $TESTDIR ]]
fi
DIR=$LIBRARY
if [[ -f $DIR/bbtools.jar ]];then CP=$DIR/bbtools.jar;else [[ -d $DIR/current ]];CP=$DIR/current;fi
if [[ -n $DEPENDENCIES ]];then
 IFS=: read -r -a jars <<< "$DEPENDENCIES"
 for jar in "${jars[@]}";do [[ -f $jar ]];done
 CP=$DEPENDENCIES:$CP
fi
if [[ -n $MODULE ]];then [[ -f $MODULE ]];CP=$MODULE:$CP;fi
: "${SHIFTER_RUNTIME:=}" "${EC2_HOME:=}" "${NERSC_HOST:=}"
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -ea -XX:ActiveProcessorCount=2 -cp "$CP" "$CLASS" "${ARGS[@]}"
