#!/usr/bin/env bash
# Converts final DDLClust members to FamilyCalibration reference buckets.
set -euo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
 printf '%s\n' \
  'Usage: ddlclustcalibrationadapter.sh members=final.members.tsv input=sampled.faa codec=final.pddl out=NEW_DIR [buckets=256]' \
  'Emits families.tsv, family_map.tsv, chunks.tsv, buckets.tsv, models.tsv, consensus.faa, and bucket_###.tsv.gz files.' \
  'Default storage is bounded resident memory; keepwork=t is rejected.' \
  'Optional maxresidentrows=ROWS maxresidentbytes=BYTES must be positive and paired; omitted caps are automatic.' \
  'compressionthreads=1..96 controls bucket_###.tsv.gz compression; values above 1 require threaded bgzip.' \
  'Use consensus=consensus.faa for FamilyCalibration and models=models.tsv consensus=consensus.faa for FamilyNegativeUpperCalibration.' \
  'Development-module-only fixture: test=t fixture=NEW_DIR runs the finite adapter regression.' \
  'Requires a separately supplied test module on the classpath; test.prot.DdlClustCalibrationAdapterTest is not included in the shipped selection.'
 exit 0
fi
script=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$script")" && pwd)
library=$DIR; modulejar=; clusterjar=; testmode=f; fixture=; args=(); jvm=()
for arg in "$@"; do
 case $arg in
  library=*) library=${arg#*=};;
  modulejar=*) modulejar=${arg#*=};;
  clusterjar=*) clusterjar=${arg#*=};;
  test=*) testmode=${arg#*=};;
  fixture=*) fixture=${arg#*=};;
  -Xmx*|-Xms*|-ea|-da|eoom|-eoom) jvm+=("$arg");;
  simd|SIMD|simd=t|simd=true|simd=f|simd=false|nosimd) jvm+=("$arg");;
  *) args+=("$arg");;
 esac
done
[[ $library == /* && -d $library/current && ( $testmode == t || $testmode == f ) ]]
[[ -z $modulejar || ( $modulejar == /* && -f $modulejar ) ]]
[[ -z $clusterjar || ( $clusterjar == /* && -f $clusterjar ) ]]
DIR=$library
: "${SHIFTER_RUNTIME:=}" "${EC2_HOME:=}" "${NERSC_HOST:=}"
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
parseJavaArgs --xmx=4g --xms=32m --mode=fixed "${jvm[@]}"
setEnvironment
class=prot.DdlClustCalibrationAdapter
if [[ $testmode == t ]]; then
 [[ -n $fixture && ${#args[@]} == 0 ]]
 class=test.prot.DdlClustCalibrationAdapterTest
 args=("$fixture")
else
 [[ -z $fixture ]]
fi
cp="$library/current"
if [[ -n $modulejar ]]; then cp="$modulejar:$cp"; fi
if [[ -n $clusterjar ]]; then cp="${modulejar:+$modulejar:}$clusterjar:$library/current"; fi
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -ea -XX:ActiveProcessorCount=1 -cp "$cp" "$class" "${args[@]}"
