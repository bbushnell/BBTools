#!/usr/bin/env bash
# Pool per-view empirical monomer and dimer substitution counts with integer weights.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 printf '%s\n' \
  'familyweightedsubstitutiontables.sh monomers=SHORT.tsv,LONG.tsv dimers=SHORT.tsv,LONG.tsv weights=SHORT_N,LONG_N out=NEW_DIR [-Xmx1g]' \
  'familyweightedsubstitutiontables.sh test runs a small weighted pooling fixture.' \
  'Weights are positive integers applied to every monomer and dimer observation before odds or matrix normalization.'
 exit
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
CLASS=prot.FamilyWeightedSubstitutionTables
if [[ $1 == test ]];then ARGS=(test); shift; else ARGS=(); fi
MODULE=
DEPENDENCIES=
CLASSPATH_FILE=;CLASSPATH_SEEN=false;MODULE_SEEN=false;DEPENDENCIES_SEEN=false
JVM_ARGS=()
for arg in "$@";do
 case ${arg,,} in
  library=*) DIR=${arg#*=};;
  modulejar=*) MODULE_SEEN=true;MODULE=${arg#*=};;
  dependencyjar=*) DEPENDENCIES_SEEN=true;DEPENDENCIES=${arg#*=};[[ -n $DEPENDENCIES ]] || { printf 'dependencyjar requires an explicit jar list\n' >&2;exit 2; };;
  classpathfile=*)
   [[ $CLASSPATH_SEEN == false ]] || { printf 'classpathfile may be supplied only once\n' >&2;exit 2; }
   CLASSPATH_SEEN=true;CLASSPATH_FILE=${arg#*=};;
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
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "${JVM_ARGS[@]}"
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
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -ea -XX:ActiveProcessorCount=1 -cp "$CP" "$CLASS" "${ARGS[@]}"
