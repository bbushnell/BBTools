#!/usr/bin/env bash
# Convert or test nucleotide HBM1 / lossless compact HBMC models.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
 echo 'Usage: nucleotidehbm.sh in=source.hbm out=fresh.hbmc[.gz]'
 echo 'Reverse: in=source.hbmc[.gz] out=fresh.hbm[.gz]'
 echo 'Tests: nucleotidehbm.sh selftest=t dir=fresh_directory'
 echo 'Parity: nucleotidehbm.sh verify=t in=source.hbm other=candidate.hbmc[.gz]'
 echo 'Counts and metadata are retained. No rare-state filtering or precision loss.'
 exit 0
fi
SCRIPT=$(readlink -f "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
source "$DIR/javasetup.sh"
ARGS=(); JVM_ARGS=(); CLASS=prok.NucleotideHbmCodec
for arg in "$@"; do
 case ${arg,,} in
  selftest=t|verify=t) CLASS=test.prok.NucleotideHbmCodecTest;;
  -xmx*|-xms*|-ea|-da|-eoom|eoom|simd=*|nosimd) JVM_ARGS+=("$arg");;
  *) ARGS+=("$arg");;
 esac
done
parseJavaArgs --xmx=2g --xms=64m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current" "$CLASS" "${ARGS[@]}"
