#!/usr/bin/env bash
set -euo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
 printf '%s\n' \
  'Usage: dev/familymemberlengthcensus.sh in=/absolute/preparation.config out=NEW_DIR' \
  '       dev/familymemberlengthcensus.sh histogram=lengths.tsv members=summary.tsv expected=N out=NEW_DIR' \
  '       dev/familymemberlengthcensus.sh reference=/absolute/references.tsv.gz marginscale=M out=NEW_DIR' \
  '' \
  'Summarize retained member lengths, join length coverage, or report saved identity gates.' \
  'No classification, clustering, alignment, consensus update or model fitting occurs.' \
  '' \
  'Member mode:' \
  'in=FILE               Native family-preparation configuration. Structural mode requires' \
  '                      an absolute, nonsymlink path to a regular configuration file.' \
  'verification=structural  Default: check structure and metadata stability without hashing.' \
  'insha80=unmeasured    Optional identity field in structural mode; an inherited 20-hex' \
  '                      identity is also accepted but does not certify current bytes.' \
  'verification=sha80   Explicit exact-match check; requires an existing matching' \
  '                      insha80 and hashes the configuration before and after reading.' \
  '                      Explicit minmemberaa/maxmemberaa allow 40..200 residues;' \
  '                      absent bounds retain the legacy 100..200 member range.' \
  '                      Writes summary.tsv, family_lengths.tsv, family_mean_lengths.tsv,' \
  '                      quantiles.tsv and PASS. Lengths describe retained peptides.' \
  '' \
  'Coverage mode:' \
  'histogram=FILE        Native Reformat length/count histogram, at most 200 residues.' \
  'members=FILE          Member-mode summary.tsv containing the four legacy length bins.' \
  'expected=N            Positive total population count; all four bins must be populated.' \
  '                      Writes coverage.tsv and PASS; conserves population/member totals.' \
  '' \
  'Threshold mode:' \
  'reference=FILE        Fitted five-line family references (legacy or metric-bound).' \
  '                      Structural mode requires an absolute,' \
  '                      nonsymlink regular file and checks its metadata before/after.' \
  'verification=structural  Default; insha80 is optional and is not remeasured.' \
  'verification=sha80   Explicit reference byte check against a supplied insha80 identity.' \
  'marginscale=M         Explicit finite, nonnegative scale; no value is chosen by the tool.' \
  '                      Writes thresholds.tsv and PASS from saved gates without refitting.' \
  '                      Quantiles use measured live families only; none measured gives NA.' \
  '                      Total, measured, unavailable and empty counts are separate columns;' \
  '                      total = measured + unavailable + empty. Empty ranks are not' \
  '                      unavailable alignment families. Nonzero unavailable counts warn' \
  '                      and carry UNEXPECTED_UNAVAILABLE; zero is reported explicitly.' \
  '                      Threshold min/median/p90/max use nearest rank. The table records' \
  '                      the bound alignment geometry, span/pair minima and matrix.' \
  '' \
  'Shared options:' \
  'out=NEW_DIR           Fresh output directory under an existing parent; existing output' \
  '                      is rejected. Summary tables are also printed to stdout.' \
  '-Xmx4g -Xms64m        Default heap limits; standard BBTools memory flags override them.' \
  '                      Supplying only -Xmx also sets -Xms to that value.' \
  'Threads               JVM processor count is fixed at 2; member reading sets 2 threads.' \
  '                      There is no t= or threads= option in this diagnostic.' \
  'eoom, simd=t|f        Standard Java setup flags. Assertions are always enabled.' \
  'library=DIR           Development override for the BBTools library directory.' \
  'modulejar=JAR         Development module placed before the library class directory.' \
  'dependencyjar=PATHS   Development dependency JARs, separated by colons.' \
  'test=t                Development-module-only boundary/aggregation fixture instead of' \
  '                      the census. Requires a separately supplied test module on the' \
  '                      classpath containing test.benchmark.FamilyMemberLengthCensusTest;' \
  '                      this test class is not included in the shipped selection.' \
  '                      Do not combine test=t with census input or output arguments.'
 exit
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
LIBRARY=$(cd "$(dirname "$SCRIPT")/.." && pwd)
MODULE= DEPENDENCY= CLASS=benchmark.FamilyMemberLengthCensus
ARGS=(); JVM_ARGS=()
for arg in "$@"; do
 case "${arg,,}" in
  library=*) LIBRARY=${arg#*=};;
  modulejar=*) MODULE=${arg#*=};;
  dependencyjar=*) DEPENDENCY=${arg#*=};;
  test=t|test=true) CLASS=test.benchmark.FamilyMemberLengthCensusTest;;
  -xmx*|-xms*|-ea|ea|-eoom|eoom|simd|simd=t|simd=f|nosimd) JVM_ARGS+=("$arg");;
  -da|da) printf '%s\n' 'Assertions must remain enabled' >&2; exit 2;;
  *) ARGS+=("$arg");;
 esac
done
DIR=$LIBRARY
[[ -d $DIR/current ]]
: "${SHIFTER_RUNTIME:=}" "${EC2_HOME:=}" "${NERSC_HOST:=}"
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
parseJavaArgs --xmx=4g --xms=64m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
CP=$DIR/current
if [[ -n $DEPENDENCY ]]; then
 IFS=: read -r -a dependencies <<< "$DEPENDENCY"
 for path in "${dependencies[@]}"; do [[ -f $path ]]; done
 CP=$DEPENDENCY:$CP
fi
if [[ -n $MODULE ]]; then [[ -f $MODULE ]]; CP=$MODULE:$CP; fi
exec java -XX:ActiveProcessorCount=2 $EA $EOOM $SIMD $XMX $XMS $PROXY -ea -cp "$CP" "$CLASS" "${ARGS[@]}"
