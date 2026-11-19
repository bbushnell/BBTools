#!/usr/bin/env bash
# Shell entrypoint for rarity-aware BBGAP reduced-alphabet sweeps.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 printf '%s\n' 'rarityawarealphabetsweep.sh in=DIMER_SYMMETRIC_TABLE out=DIR [insha80=PIN] [verification=sha80|structural] [band=<40] [lambdas=0,0.1,0.25,0.5,1] [targets=120] [-Xmx4g]' \
  'Reads retained symmetric dimer odds, keeps complete-link original-odds scoring, and adds a normalized group-frequency uniformity gain to each greedy merge.' \
  'Optional library=ROOT, classpathfile=FILE, modulejar=JAR and dependencyjar=JAR[:JAR...] select the installed runtime and reviewed additive modules.' \
  'verification=structural retains input metadata and all table/partition checks, emits no digests, and rejects insha80. Default sha80 behavior is unchanged.' \
  'format=member_pair_counts reads the pooled five-column ordered-count export directly, using the original symmetric odds arithmetic. Default format=symmetric reads nine-column saved odds.'
 exit
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
CLASS=prot.RarityAwareAlphabetSweep
MODULE=
DEPENDENCIES=
CP_FILE=
ARGS=();JVM_ARGS=();THREADS=1
for arg in "$@";do
 case ${arg,,} in
  library=*) DIR=${arg#*=};;
  classpathfile=*) [[ -z $CP_FILE ]] || exit 2; CP_FILE=${arg#*=};;
  modulejar=*) MODULE=${arg#*=};;
  dependencyjar=*) DEPENDENCIES=${arg#*=};[[ -n $DEPENDENCIES ]] || { printf 'dependencyjar requires an explicit jar list\n' >&2;exit 2; };;
  t=*|threads=*) THREADS=${arg#*=};[[ $THREADS =~ ^[1-9][0-9]*$ ]] || { printf 'threads must be positive: %s\n' "$arg" >&2;exit 2; };;
  -xmx*|-xms*|-ea|ea|-eoom|eoom|simd|simd=t|simd=f|nosimd) JVM_ARGS+=("$arg");;
  -da) printf 'Assertions must remain enabled\n' >&2;exit 2;;
  *) ARGS+=("$arg");;
 esac
done
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
parseJavaArgs --xmx=4g --xms=64m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
CP=$DIR/current
if [[ -n $CP_FILE ]];then
 [[ -z $MODULE && -z $DEPENDENCIES && -f $CP_FILE && -s $CP_FILE ]] || exit 2
 mapfile -t lines < "$CP_FILE"
 [[ ${#lines[@]} == 1 ]] || exit 2
 CP=${lines[0]}
 [[ -n $CP && $CP != :* && $CP != *: && $CP != *::* && $CP != *$'\r'* ]] || exit 2
 IFS=: read -r -a entries <<< "$CP"
 for entry in "${entries[@]}";do [[ $entry == /* && ( -f $entry || -d $entry ) ]] || exit 2;done
elif [[ -n $DEPENDENCIES ]];then
 [[ $DEPENDENCIES != :* && $DEPENDENCIES != *: && $DEPENDENCIES != *::* ]] || { printf 'Empty dependency classpath component\n' >&2;exit 2; }
 IFS=: read -r -a DEPENDENCY_PATHS <<< "$DEPENDENCIES"
 for jar in "${DEPENDENCY_PATHS[@]}";do [[ $jar == /* && -f $jar ]] || { printf 'Missing absolute dependency jar: %s\n' "$jar" >&2;exit 2; };done
 CP=$DEPENDENCIES:$CP
fi
if [[ -n $MODULE ]];then [[ $MODULE == /* && -f $MODULE ]];CP=$MODULE:$CP;fi
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -ea -XX:ActiveProcessorCount="$THREADS" -cp "$CP" "$CLASS" "${ARGS[@]}"
