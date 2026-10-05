#!/usr/bin/env bash
# GradeBins' report-only ProkCC import regressions; no model loading.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
  echo 'Usage: gradebinsprokcctest.sh mode=prepare dir=new_dir real=prokcc_report.tsv'
  echo 'After the saved CLI comparisons: mode=verify dir=fixture_dir'
  echo 'RNA import comparisons: mode=rna dir=fixture_dir; standalone formatter: reportselftest=t'
  exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
if [[ -f $DIR/bbtools.jar ]]; then CP=$DIR/bbtools.jar; else CP=$DIR/current; fi
JAVA=$(command -v java)
source "$DIR/javasetup.sh"
ARGS=(); JVM_ARGS=()
CLASS=bin.GradeBinsProkCCTest
for arg in "$@"; do
  case ${arg,,} in
    reportselftest=t) CLASS=prot.MagQCAssemblyReportTest;;
    -ea|-da|-eoom|eoom|-xmx*|-xms*|simd|simd=t|simd=f) JVM_ARGS+=("$arg");;
    *) ARGS+=("$arg");;
  esac
done
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec "$JAVA" $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" "$CLASS" "${ARGS[@]}"
