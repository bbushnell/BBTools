#!/usr/bin/env bash
# Input preparation for an experimental extension of the protein HBM library.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
 cat <<'HELP'
Usage: hbmfamilyexpand.sh mode=select in=clusters.tsv original=familylist.tsv
         roster=roster.tsv.gz out=prefix count=4000
       hbmfamilyexpand.sh mode=extract in=clusters.tsv selected=prefix.selected.tsv
         seqs=families_all_seqs.fasta out=NEW_DIRECTORY
Selection ranks unrepresented original clusters by membership count, then ID.
Extraction validates every selected ID against MMseqs segmented FASTA records.
The original empty cluster-boundary headers are required in seqs=.
Default heap: 4g. Existing outputs are refused. One input-processing thread.
HELP
 exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
if [[ -d $DIR/current ]]; then CP=$DIR/current; else CP=$DIR/bbtools.jar; fi
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
ARGS=(); JVM_ARGS=(); TASK=prot.AdditionalFamilySelector
for arg in "$@"; do
 case ${arg,,} in
  mode=select) TASK=prot.AdditionalFamilySelector;;
  mode=extract) TASK=prot.ClusterMemberExtractor;;
  --xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
  -ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
  simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("$arg");;
  *) ARGS+=("$arg");;
 esac
done
parseJavaArgs --xmx=4g --xms=256m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -XX:ActiveProcessorCount=1 -cp "$CP" "$TASK" "${ARGS[@]}"
