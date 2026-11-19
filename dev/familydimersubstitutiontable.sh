#!/usr/bin/env bash
# Empirical overlapping dimer substitutions; explicit structural verification.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 printf '%s\n' \
  'familydimersubstitutiontable.sh families=TSV familiessha80=PIN chunks=TSV chunkssha80=PIN consensus=FAA consensussha80=PIN expectedfamilies=N expectedmembers=N batches=1 buckets=256 maxreferenceaa=1000 maxbucketmembers=N maxbucketbytes=N verification=structural monomercontrol=TSV out=NEW_DIR t=64 [-Xmx8g]' \
  'familydimersubstitutiontable.sh test [accepted_symmetric_tables.tsv accepted_merge_groups.tsv] runs independent trace and merge-oracle fixtures.' \
  'familydimersubstitutiontable.sh cuts in=symmetric_tables.tsv.gz ladder=merge_ladder.tsv groups=merge_groups.tsv maps=maps.tsv cuts=140,120,80 out=/absolute/NEW_DIR t=1 reuses original odds without recounting.' \
  'Cut export requires byte-identical accepted 300-merge prefix and 300/200/100 maps; these controls are also exported.' \
  'familydimersubstitutiontable.sh cutstest /absolute/NEW_DIR runs cut-export fixtures on Dori.' \
  'Development-module-only: test and cutstest require a separately supplied test module on the classpath; test.prot.FamilyDimerSubstitutionTableTest is not included in the shipped selection.' \
  'The Java test entry also accepts fixture=NEW_DIR launcher=SCRIPT library=ROOT together for full driver/config and rejected-control cases.' \
  'identitystudy=t counts four disjoint bands and exports cumulative all/<80/<60/<40 alphabets: monomer12/16 and ordered-dimer120.' \
  'Optional bands=0to40,25to40,25to45,0to35 emits overlapping lower-inclusive, upper-exclusive identity intervals and writes monomer_counts.tsv without requiring monomercontrol.' \
  'Counts overlapping ordered amino-acid dimers in consecutive standard aligned columns; gaps and unknown residues break adjacency.' \
  'Input byte identity is inherited, not recomputed. Exact accepted monomer-table parity is required before emitting 300/200/100 dimer maps.' \
  'Optional library=ROOT, modulejar=JAR, dependencyjar=JAR[:JAR...] select the reviewed runtime.' \
  'Alternatively classpathfile=PATH supplies one exact classpath line; every absolute provider and its order are preserved.' \
  'classpathfile is mutually exclusive with modulejar/dependencyjar. Tests and builds run on Dori only.'
 exit
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
CLASS=prot.FamilyDimerSubstitutionTable
if [[ $1 == test ]];then CLASS=test.prot.FamilyDimerSubstitutionTableTest;shift;fi
cutstest=false
if [[ ${1:-} == cutstest ]];then CLASS=test.prot.FamilyDimerSubstitutionTableTest;cutstest=true;shift;fi
MODULE=
DEPENDENCIES=
CLASSPATH_FILE=;CLASSPATH_SEEN=false;MODULE_SEEN=false;DEPENDENCIES_SEEN=false
ARGS=();JVM_ARGS=();THREADS=1
if $cutstest;then ARGS+=(cutstest);fi
for arg in "$@";do
 case ${arg,,} in
  library=*) DIR=${arg#*=};;
  modulejar=*) MODULE_SEEN=true;MODULE=${arg#*=};;
  dependencyjar=*) DEPENDENCIES_SEEN=true;DEPENDENCIES=${arg#*=};[[ -n $DEPENDENCIES ]] || { printf 'dependencyjar requires an explicit jar list\n' >&2;exit 2; };;
  classpathfile=*)
   [[ $CLASSPATH_SEEN == false ]] || { printf 'classpathfile may be supplied only once\n' >&2;exit 2; }
   CLASSPATH_SEEN=true;CLASSPATH_FILE=${arg#*=};;
  t=*|threads=*) THREADS=${arg#*=};[[ $THREADS =~ ^[1-9][0-9]*$ && $THREADS -le 64 ]] || { printf 'threads must be 1..64: %s\n' "$arg" >&2;exit 2; };if [[ $CLASS == prot.* ]];then ARGS+=("t=$THREADS");fi;;
  -xmx*|-xms*|-ea|ea|-eoom|eoom|simd|simd=t|simd=f|nosimd) JVM_ARGS+=("$arg");;
  -da) printf 'Assertions must remain enabled\n' >&2;exit 2;;
  *) ARGS+=("$arg");;
 esac
done
if $CLASSPATH_SEEN;then
 [[ $MODULE_SEEN == false && $DEPENDENCIES_SEEN == false ]] || {
  printf 'classpathfile cannot be combined with modulejar or dependencyjar\n' >&2;exit 2;
 }
 [[ -f $CLASSPATH_FILE && -s $CLASSPATH_FILE ]] || { printf 'classpathfile requires an existing nonempty file\n' >&2;exit 2; }
 mapfile -t CLASSPATH_LINES < "$CLASSPATH_FILE"
 [[ ${#CLASSPATH_LINES[@]} == 1 ]] || { printf 'classpathfile requires exactly one line\n' >&2;exit 2; }
 EXPLICIT_CP=${CLASSPATH_LINES[0]}
 [[ -n $EXPLICIT_CP && $EXPLICIT_CP != :* && $EXPLICIT_CP != *: && $EXPLICIT_CP != *::* && $EXPLICIT_CP != *$'\r'* ]] || {
  printf 'classpathfile has an empty provider or invalid line ending\n' >&2;exit 2;
 }
 IFS=: read -r -a CLASSPATH_PROVIDERS <<< "$EXPLICIT_CP"
 for provider in "${CLASSPATH_PROVIDERS[@]}";do
  [[ $provider == /* && ( -f $provider || -d $provider ) ]] || {
   printf 'Classpath provider must be an existing absolute file or directory: %s\n' "$provider" >&2;exit 2;
  }
 done
fi
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
parseJavaArgs --xmx=8g --xms=64m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
CP=$DIR/current
if $CLASSPATH_SEEN;then
 CP=$EXPLICIT_CP
elif [[ -n $DEPENDENCIES ]];then
 [[ $DEPENDENCIES != :* && $DEPENDENCIES != *: && $DEPENDENCIES != *::* ]] || { printf 'Empty dependency classpath component\n' >&2;exit 2; }
 IFS=: read -r -a DEPENDENCY_PATHS <<< "$DEPENDENCIES"
 for jar in "${DEPENDENCY_PATHS[@]}";do [[ $jar == /* && -f $jar ]] || { printf 'Missing absolute dependency jar: %s\n' "$jar" >&2;exit 2; };done
 CP=$DEPENDENCIES:$CP
fi
if [[ -n $MODULE ]];then [[ $MODULE == /* && -f $MODULE ]];CP=$MODULE:$CP;fi
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -ea -XX:ActiveProcessorCount="$THREADS" -cp "$CP" "$CLASS" "${ARGS[@]}"
