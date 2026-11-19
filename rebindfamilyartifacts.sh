#!/usr/bin/env bash
usage(){
cat <<'HELP'
Rebind explicitly listed legacy artifact paths, preserving other bytes and resolved dependencies.
Usage: rebindfamilyartifacts.sh in=bindings.tsv out=/fresh/package
Input magic: #bbgap_family_rebind<TAB>1
Columns: artifact, format, source, target, key, expected, replacement (tab-separated).
Formats: classifier or cascade. Keys: blacklist; cascade also allows short40manifest.
Sources/replacements are absolute; targets are safe basenames in the fresh output root.
SHORT edits require their composite root in the same batch. Diagnostic-matrix, fitted-floor,
scalar-mean and cyclic graphs are rejected before structural readers. No payload hashing or refitting.
Outputs: rebound artifacts, bindings.tsv, PASS. Failure has no PASS; existing roots are never overwritten.

library= selects an installed BBTools directory (default: this script's directory).
classpathfile= supplies one literal classpath line of existing absolute providers.
-Xmx/-Xms, assertions and ordinary BBTools JVM options are supported.
HELP
}
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 usage
 exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")" && pwd)
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
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -ea -cp "$CP" prot.FamilyArtifactRebind "${ARGS[@]}"
