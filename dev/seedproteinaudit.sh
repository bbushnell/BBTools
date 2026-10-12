#!/usr/bin/env bash
# Audits raw seed residues; explicit trimedges mode writes a derived repaired copy.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
 cat <<'HELP'
Usage: seedproteinaudit.sh mode=audit manifest=family_manifest.tsv manifestsha80=PIN source=DIR
  out=NEW_DIR
Reads pinned seed FASTAs and reports edge stops, internal stops and unsupported
residues. This does not modify or filter proteins. Required: manifest=,
manifestsha80=, source= and out=NEW_DIR.
mode=trimedges writes members/ and a new manifest, repairing only records with
edge-only stops. Every record is retained; internal stops/unsupported residues
fail. Already valid records, including one terminal stop marker, are unchanged.
HELP
 exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
if [[ -d $DIR/current ]]; then CP=$DIR/current; else CP=$DIR/bbtools.jar; fi
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
ARGS=(); JVM_ARGS=()
for arg in "$@"; do
 case ${arg,,} in
  --xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
  -ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
  simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("$arg");;
  *) ARGS+=("$arg");;
 esac
done
parseJavaArgs --xmx=4g --xms=256m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -XX:ActiveProcessorCount=1 -cp "$CP" prot.SeedProteinAudit "${ARGS[@]}"
