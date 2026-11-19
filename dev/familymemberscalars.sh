#!/usr/bin/env bash
# Shell entrypoint for the installed all-member scalar-moment producer.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 printf '%s\n' 'familymemberscalars.sh families=TSV familiessha80=IDENTITY chunks=TSV chunkssha80=IDENTITY expectedfamilies=N expectedmembers=N batches=1 buckets=256 shards=256 shard=N out=TSV ow=f t=1 [-Xmx2g]' \
  'familymemberscalars.sh collect families=TSV familiessha80=IDENTITY parts=TSV partssha80=IDENTITY expectedfamilies=N expectedmembers=N shards=256 out=TSV ow=f [-Xmx2g]' \
  'Reads all assigned members in the selected rank shard; emits family N and five scalar mean/SD/min/max sets. No consensus-only substitute.' \
  'verification=structural is the default: validates structure, counts and stable file metadata without hashing. IDENTITY must be unmeasured or a lowercase 20-hex label; manifest identities follow the same rule.' \
  'verification=sha80 checks exact stored bytes against supplied lowercase 20-hex identities for the roster, manifest and every selected chunk or declared part; unmeasured is rejected.' \
  'Producer member bounds default to 100..200 amino acids. Supply minmemberaa=N maxmemberaa=N together for another declared interval within 40..200; out-of-range members fail, not filter.' \
  'scalars=hydrophobicity emits only family_rank, family_id, N, hydrophobicity_mean, equally weighting retained-member KD means; empty means are NA.' \
  'The default scalars=all retains the five-scalar schema. Use shards=1 shard=0 for one complete KD table; collect supports only the default schema.' \
  'collect validates and joins shard outputs from a manifest with the literal header shard/file/sha80 (tab-separated; sha80 holds the identity label), without recomputing moments; it requires a fresh absolute out= path and ow=f.' \
  'Optional library=ROOT, modulejar=JAR and dependencyjar=JAR[:JAR...] select the installed runtime and reviewed additive modules.' \
  'Development-module-only: familymemberscalars.sh test root=NEW_FIXTURE_DIR runs the finite scalar checks.' \
  'Requires a separately supplied test module on the classpath; test.prot.FamilyMemberScalarsTest is not included in the shipped selection.'
 exit
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
CLASS=prot.FamilyMemberScalars
if [[ ${1,,} == collect ]];then CLASS=prot.FamilyMemberScalarsCollect;shift;fi
if [[ ${1,,} == test ]];then CLASS=test.prot.FamilyMemberScalarsTest;shift;fi
MODULE=
DEPENDENCIES=
ARGS=();JVM_ARGS=();THREADS=1
for arg in "$@";do
 case ${arg,,} in
  library=*) DIR=${arg#*=};;
  modulejar=*) MODULE=${arg#*=};;
  dependencyjar=*) DEPENDENCIES=${arg#*=};[[ -n $DEPENDENCIES ]] || { printf 'dependencyjar requires an explicit jar list\n' >&2;exit 2; };;
  t=*|threads=*) THREADS=${arg#*=};[[ $THREADS =~ ^[1-9][0-9]*$ ]];ARGS+=("$arg");;
  -xmx*|-xms*|-ea|ea|-eoom|eoom|simd|simd=t|simd=f|nosimd) JVM_ARGS+=("$arg");;
  -da) printf 'Assertions must remain enabled\n' >&2;exit 2;;
  *) ARGS+=("$arg");;
 esac
done
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
parseJavaArgs --xmx=2g --xms=64m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
CP=$DIR/current
if [[ -n $DEPENDENCIES ]];then
 [[ $DEPENDENCIES != :* && $DEPENDENCIES != *: && $DEPENDENCIES != *::* ]] || { printf 'Empty dependency classpath component\n' >&2;exit 2; }
 IFS=: read -r -a DEPENDENCY_PATHS <<< "$DEPENDENCIES"
 for jar in "${DEPENDENCY_PATHS[@]}";do [[ $jar == /* && -f $jar ]] || { printf 'Missing absolute dependency jar: %s\n' "$jar" >&2;exit 2; };done
 CP=$DEPENDENCIES:$CP
fi
if [[ -n $MODULE ]];then [[ $MODULE == /* && -f $MODULE ]];CP=$MODULE:$CP;fi
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -ea -XX:ActiveProcessorCount="$THREADS" -cp "$CP" "$CLASS" "${ARGS[@]}"
