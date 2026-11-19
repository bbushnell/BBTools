#!/usr/bin/env bash
usage(){
cat <<'HELP'
Extract paired-caller scorer metadata from four saved tables. No sequence reads or caller runs.
Usage: dev/callgenesmetadata.sh accepted=accepted_inputs.tsv raw=by_genome.tsv \
  retained=retained_gffs.tsv summary=scorecard.tsv out=/fresh/metadata policy=cached expectedgenomes=983
Outputs: grades.tsv, inputs.tsv, expected.tsv, PASS. Existing output roots are rejected.
The cohort size and safe policy prefix are explicit; native/formula aliases and saved order are retained.

library= overrides the detected installation root.
classpathfile= supplies one literal classpath line of existing absolute providers.
-Xmx/-Xms, assertions and ordinary BBTools JVM options are supported.
Installation root is detected automatically from top-level or dev/ placement.
HELP
}
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 usage
 exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")" && pwd)
if [[ ! -f "$DIR/javasetup.sh" ]]; then DIR=$(dirname "$DIR"); fi
LIBRARY=$DIR
CLASSPATH_FILE=
EXACT_CP=false
ARGS=()
JVM_ARGS=()
for arg in "$@";do
 case ${arg,,} in
  library=*) LIBRARY=${arg#*=};;
  classpathfile=*)
   [[ $EXACT_CP == false ]] || { printf 'Duplicate classpathfile\n' >&2;exit 2; }
   CLASSPATH_FILE=${arg#*=}
   EXACT_CP=true;;
  -xmx*|-xms*|-ea|ea|-eoom|eoom|simd|simd=t|simd=f|nosimd) JVM_ARGS+=("$arg");;
  -da) printf 'Assertions must remain enabled\n' >&2;exit 2;;
  *) ARGS+=("$arg");;
 esac
done
DIR=$LIBRARY
if [[ $EXACT_CP == true ]];then
 [[ -f $CLASSPATH_FILE && -s $CLASSPATH_FILE ]] || { printf 'Nonempty classpathfile required\n' >&2;exit 2; }
 mapfile -t cp_lines < "$CLASSPATH_FILE" || exit 2
 [[ ${#cp_lines[@]} == 1 ]] || { printf 'classpathfile requires one line\n' >&2;exit 2; }
 CP=${cp_lines[0]}
 [[ -n $CP && $CP != :* && $CP != *: && $CP != *::* && $CP != *$'\r'* ]] || { printf 'Invalid classpath line\n' >&2;exit 2; }
 IFS=: read -r -a cp_entries <<< "$CP"
 for path in "${cp_entries[@]}";do
  [[ $path == /* && ( -f $path || -d $path ) ]] || { printf 'Existing absolute classpath provider required: %s\n' "$path" >&2;exit 2; }
 done
elif [[ -f $DIR/bbtools.jar ]];then
 CP=$DIR/bbtools.jar
else
 CP=$DIR/current
fi
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=32m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -ea -cp "$CP" benchmark.CallGenesSavedMetadata "${ARGS[@]}"
