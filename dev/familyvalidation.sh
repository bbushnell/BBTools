#!/usr/bin/env bash
# One portable entry point for the public family fixture programs.
set -eo pipefail
usage(){
 cat <<'HELP'
Usage: familyvalidation.sh mode=paired|inputs|admission|selector [root=NEW_ABSOLUTE_DIR]

paired    Build a tiny two-family parent and test paired evaluation, ordinary
          behavior, malformed policies and one/four-worker parity.
inputs    Generate the paired and hit-floor inputs, manifests and malformed
          variants used by default-hash wrapper validation. Does not run those
          wrapper comparisons; generation alone is not their acceptance.
admission Check public family hit-floor admission, binding, observations,
          callback handling and unchanged fixed-floor controls.
selector  Check the public covering selector against an independent full-sort
          oracle over 41,878 finite cases. Creates no files; omit root=.

paired, inputs and admission require a fresh absolute root with an existing parent.
No real population, historical private exporter or private benchmark is used.
The paired and inputs modes share one public parent-fixture implementation.
Defaults: mode=paired; admission uses -Xmx4g, other modes use -Xmx1g.
Explicit heap options override these defaults. SIMD and exit-on-OOM options are accepted;
assertions always remain enabled for validation.
HELP
}
die(){ printf '%s\n' "$*" >&2; exit 2; }
if [[ $# == 0 || $1 == -h || $1 == --help || $1 == help ]];then usage; exit 0;fi
fixture_mode=; fixture_root=
jvm_args=()
for arg in "$@";do
 case ${arg,,} in
  mode=*) [[ -z $fixture_mode ]] || die 'Duplicate mode'; fixture_mode=${arg#*=}; [[ -n $fixture_mode ]] || die 'Empty mode';;
  root=*) [[ -z $fixture_root ]] || die 'Duplicate root'; fixture_root=${arg#*=}; [[ -n $fixture_root ]] || die 'Empty root';;
  --xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
  -ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
  simd|simd=t|simd=true|simd=f|simd=false|nosimd) jvm_args+=("$arg");;
  *) die "Unknown argument: $arg";;
 esac
done
fixture_mode=${fixture_mode:-paired}
program_args=()
case ${fixture_mode,,} in
 paired) fixture_class=test.prot.FamilyPairedEvalTest;;
 inputs) fixture_class=test.prot.FamilyDefaultHashInputFixture;;
 admission) fixture_class=test.prot.FamilyHitFloorIntegrationTest;;
 selector) fixture_class=test.prot.HbmTailCoveringSelectionTest;;
 *) die "Unknown mode: $fixture_mode";;
esac
if [[ ${fixture_mode,,} == selector ]];then
 [[ -z $fixture_root ]] || die 'selector does not take root='
else
 [[ $fixture_root == /* ]] || die 'paired, inputs and admission require an absolute root='
 [[ ! -e $fixture_root && ! -L $fixture_root ]] || die 'root must be fresh'
 program_args+=("$fixture_root")
fi
script=$(readlink -f -- "${BASH_SOURCE[0]}")
bbtools_root=$(cd "$(dirname "$script")/.." && pwd)
if [[ -f $bbtools_root/bbtools.jar ]];then cp_path=$bbtools_root/bbtools.jar; else cp_path=$bbtools_root/current;fi
DIR="$bbtools_root"
source "$bbtools_root/javasetup.sh"
default_heap=1g
if [[ ${fixture_mode,,} == admission ]];then default_heap=4g;fi
parseJavaArgs "--xmx=$default_heap" --xms=64m --mode=fixed "${jvm_args[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -ea -cp "$cp_path" "$fixture_class" "${program_args[@]}"
