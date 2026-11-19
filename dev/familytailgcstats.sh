#!/usr/bin/env bash
# Measure family tail GC from retained memberships, extraction metadata and DNA.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Usage: bash dev/familytailgcstats.sh inputs=ordered_inputs.tsv \
 members=members.tsv.gz familymap=family_map.tsv genomes=genomes.tsv \
 readers=readers.tsv records=N assigned=N families=N minaa=40 maxaa=99 \
 out=/absolute/new_directory [-Xmx8g]

Compute equal-member GC mean and sample standard deviation for each family.
Join the exact ordered FASTA inputs and memberships with retained extraction
metadata and original DNA. Verify each measured tail translates to its peptide.
The measured interval contains sense codons and excludes the stop codon.

inputs=       Two-column ordinal/path manifest in original clustering order.
members=      Four-column input_index/input_id/family_id/status membership table.
familymap=    Five-column stable_family_id/roster_rank/family_id/assigned/residues.
genomes=      Original ten-column genome manifest, including source DNA paths.
readers=      Six-column extraction-reader manifest, including metadata paths.
records=      Exact total input records, including unassigned records.
assigned=     Exact assigned-member count.
families=     Exact positive family count in the map.
minaa=40 maxaa=99    Whole short-protein view.
minaa=100 maxaa=100  Terminal 100-amino-acid view of retained long proteins.
out=          Fresh absolute output directory; writes family_gc.tsv and PASS.

library=DIR        BBTools installation root; default parent of this dev directory.
classpathfile=FILE Optional exact one-line classpath of existing absolute providers.
                   Otherwise use library/bbtools.jar if present, else library/current.
JVM options include -Xmx/-Xms, eoom and simd=t/f. Assertions remain enabled.
This wrapper calls the existing GC producer only; it does not cluster, select
inputs, submit jobs or run a GC-threshold sweep. No thread-count flag is exposed.
HELP
 exit 0
fi

SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
CLASSPATH_FILE=;CLASSPATH_SEEN=false
ARGS=();JVM_ARGS=()
for arg in "$@";do
 case ${arg,,} in
  library=*) DIR=${arg#*=};;
  classpathfile=*)
   [[ $CLASSPATH_SEEN == false ]] || { printf 'classpathfile may be supplied only once\n' >&2;exit 2; }
   CLASSPATH_SEEN=true;CLASSPATH_FILE=${arg#*=};;
  -xmx*|-xms*|-ea|ea|-eoom|eoom|simd|simd=t|simd=f|nosimd) JVM_ARGS+=("$arg");;
  -da|da|-da:*) printf 'Assertions must remain enabled\n' >&2;exit 2;;
  *) ARGS+=("$arg");;
 esac
done

if $CLASSPATH_SEEN;then
 [[ -f $CLASSPATH_FILE && -s $CLASSPATH_FILE ]] || { printf 'classpathfile requires an existing nonempty file\n' >&2;exit 2; }
 mapfile -t lines < "$CLASSPATH_FILE"
 [[ ${#lines[@]} == 1 ]] || { printf 'classpathfile requires exactly one line\n' >&2;exit 2; }
 CP=${lines[0]}
 [[ -n $CP && $CP != :* && $CP != *: && $CP != *::* && $CP != *$'\r'* ]] || {
  printf 'classpathfile has an empty provider or invalid line ending\n' >&2;exit 2;
 }
 IFS=: read -r -a providers <<< "$CP"
 for provider in "${providers[@]}";do
  [[ $provider == /* && ( -f $provider || -d $provider ) ]] || {
   printf 'Classpath provider must be an existing absolute file or directory: %s\n' "$provider" >&2;exit 2;
  }
 done
elif [[ -f $DIR/bbtools.jar ]];then CP=$DIR/bbtools.jar
else CP=$DIR/current
fi

source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
parseJavaArgs --xmx=8g --xms=64m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -ea -XX:ActiveProcessorCount=4 -cp "$CP" prot.cluster.FamilyTailGcStats "${ARGS[@]}"
