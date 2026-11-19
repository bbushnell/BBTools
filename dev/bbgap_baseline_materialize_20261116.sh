#!/usr/bin/env bash
# Fixed full-LONG OLDSTRICT config/command materializer. Never runs science.
# Required: mode=collection|lengths|hydro|prepare|fit out=NEW_CONTROL run=NEW_RUN
# java=ABS_JAVA (or jdk=ABS_JDK_HOME) classpathfile=ONE_LINE_CPP custom40=ABS_FALLBACK
# Optional toolpath=ABS_DIR_LIST; defaults to /usr/bin:/bin after the selected Java bin.
# Config stages: parent=ACCEPTED_NATIVE_CONFIG. Metadata/prep: collection=RESULT.
# Prepare also needs lengthmedians=FILE familyhydro=FILE resourceroot=DIR.
# Fit needs prepared=CLASSIFIER resourceroot=DIR. Optional prepare-only blacklist
# and blacklistprovenance come from the accepted header/resource handoff.
# Fit requires negk5=FILE negshadow=FILE negreversed=FILE negshuffled=FILE.
# These relocate the four accepted cohorts; names, identity labels and ceilings stay fixed.
# GNU realpath checks destination identity; generated paths retain the supplied spelling.
set -euo pipefail
export LC_ALL=C TZ=UTC

usage(){
	cat <<'EOF'
Usage: bbgap_baseline_materialize_20261116.sh key=value ...

Write fixed full-LONG baseline configs and commands; do not run them.
All file and directory arguments must be absolute paths.

Required for every mode:
  mode=collection|lengths|hydro|prepare|fit
  out=NEW_CONTROL_DIR run=NEW_RUN_DIR
  java=EXECUTABLE or jdk=JDK_HOME
  classpathfile=FILE custom40=FALLBACK_MATRIX

Additional required arguments by mode:
  collection  parent=ACCEPTED_NATIVE_CONFIG
  lengths     parent=ACCEPTED_NATIVE_CONFIG collection=RESULT_DIR
  hydro       collection=RESULT_DIR
  prepare     parent=ACCEPTED_NATIVE_CONFIG collection=RESULT_DIR
              lengthmedians=FILE familyhydro=FILE resourceroot=DIR
  fit         prepared=CLASSIFIER resourceroot=DIR
              negk5=FILE negshadow=FILE negreversed=FILE negshuffled=FILE

Optional:
  toolpath=DIRS         Colon-separated absolute directories; /usr/bin:/bin
                       by default, after the selected Java bin directory.
  blacklist=FILE        Prepare-only accepted blacklist replacement.
  blacklistprovenance=TEXT  Prepare-only accepted resource provenance.

The classpath file contains exactly one line. Control and run directories must
be new and disjoint. Writes bindings.tsv, classpath.txt, native.argv and
native.command.sh under out, plus native.config for collection, lengths and
prepare. The generated command requires a Slurm allocation. Fit preserves the
four accepted negative cohorts and their fixed names, labels and ceilings.
No arguments, -h or --help prints this help and exits without creating files.
EOF
}
if [[ $# == 0 || ${1:-} == -h || ${1:-} == --help ]];then
	usage
	exit 0
fi

die(){ printf '%s\n' "$*" >&2;exit 2; }
declare -A option=()
for arg in "$@";do
	[[ $arg == *=* ]] || die 'Use key=value arguments'
	key=${arg%%=*};key=${key,,};value=${arg#*=}
	case $key in mode|out|run|java|jdk|toolpath|classpathfile|custom40|parent|collection|lengthmedians|familyhydro|resourceroot|prepared|blacklist|blacklistprovenance|negk5|negshadow|negreversed|negshuffled);;*)die "Unknown option: $key";;esac
	[[ -n $value && ! ${option[$key]+present} && $value != *[$'\t\r\n\\']* ]] || die "Empty, duplicate or escaped option: $key"
	option[$key]=$value
done
required(){ [[ ${option[$1]+present} ]] || die "Required option: $1"; }
absolute(){ required "$1";[[ ${option[$1]} == /* ]] || die "Absolute path required: $1"; }
if [[ ${option[jdk]+present} ]];then
	absolute jdk
	jdk_java=${option[jdk]%/}/bin/java
	if [[ ${option[java]+present} ]];then
		[[ ${option[java]} == "$jdk_java" ]] || die 'java and jdk must name the same Java executable'
	else
		option[java]=$jdk_java
	fi
fi
option[toolpath]=${option[toolpath]-/usr/bin:/bin}
[[ ${option[toolpath]} != :* && ${option[toolpath]} != *: && ${option[toolpath]} != *::* ]] || die 'Malformed toolpath'
IFS=: read -r -a tool_directories <<< "${option[toolpath]}"
for directory in "${tool_directories[@]}";do [[ $directory == /* ]] || die 'toolpath entries must be absolute';done
for key in out run java classpathfile custom40;do absolute "$key";done
java_bin=${option[java]%/*}
[[ -n $java_bin ]] || java_bin=/
runtime_path=$java_bin:${option[toolpath]}
required mode
mode=${option[mode]};out=${option[out]};run=${option[run]}
[[ ! -e $out && ! -L $out && ! -e $run && ! -L $run ]] || die 'Control and run must both be fresh paths'
physical_out=$(realpath -m -- "$out") || die 'Cannot resolve control path'
physical_run=$(realpath -m -- "$run") || die 'Cannot resolve run path'
[[ ! -e $physical_out && ! -L $physical_out && ! -e $physical_run && ! -L $physical_run ]] || die 'Control and run must both resolve to fresh paths'
[[ $physical_out != "$physical_run" && $physical_out != "$physical_run/"* && $physical_run != "$physical_out/"* ]] || die 'Control and run paths must resolve to disjoint destinations'
[[ -x ${option[java]} && -s ${option[classpathfile]} && -s ${option[custom40]} ]] || die 'Missing Java, classpath or fallback binding'
mapfile -t cpp_lines < "${option[classpathfile]}"
[[ ${#cpp_lines[@]} == 1 ]] || die 'Classpath file must contain exactly one line'
cpp=${cpp_lines[0]}
[[ -n $cpp && $cpp != :* && $cpp != *: && $cpp != *::* && $cpp != *[$'\r\n']* ]] || die 'Malformed classpath'
case $mode in
	collection)absolute parent;heap=64g;threads=1;;
	lengths)absolute parent;absolute collection;heap=4g;threads=2;;
	hydro)absolute collection;heap=4g;threads=4;;
	prepare)absolute parent;absolute collection;absolute lengthmedians;absolute familyhydro;absolute resourceroot;heap=128g;threads=64;;
	fit)absolute prepared;absolute resourceroot;heap=64g;threads=8;;
	*)die 'Unknown mode';;
esac
for key in negk5 negshadow negreversed negshuffled;do
	if [[ $mode == fit ]];then
		absolute "$key"
		[[ ${option[$key]} != *,* ]] || die "Negative input path cannot contain a comma: $key"
	else
		[[ ! ${option[$key]+present} ]] || die "$key is fit-only"
	fi
done
if [[ ${option[parent]+present} ]];then [[ -s ${option[parent]} ]] || die 'Missing parent native config';fi
if [[ ${option[blacklist]+present} ]];then absolute blacklist;[[ $mode == prepare ]] || die 'blacklist is prepare-only';fi
if [[ ${option[blacklistprovenance]+present} ]];then [[ $mode == prepare ]] || die 'blacklistprovenance is prepare-only';fi
mkdir -p -- "${out%/*}"
mkdir -- "$out"
cp -- "${option[classpathfile]}" "$out/classpath.txt"
printf 'format\tbbgap_baseline_materialization_v1\n' > "$out/bindings.tsv"
for key in mode out run java jdk toolpath classpathfile custom40 parent collection lengthmedians familyhydro resourceroot prepared blacklist blacklistprovenance negk5 negshadow negreversed negshuffled;do
	[[ ! ${option[$key]+present} ]] || printf '%s\t%s\n' "$key" "${option[$key]}" >> "$out/bindings.tsv"
done
if [[ $mode == collection || $mode == lengths || $mode == prepare ]];then
	awk -F= -v mode="$mode" -v run="$run" -v collection="${option[collection]-}" \
	 -v medians="${option[lengthmedians]-}" -v hydro="${option[familyhydro]-}" \
	 -v blacklist="${option[blacklist]-}" -v provenance="${option[blacklistprovenance]-}" '
	 BEGIN{
	  if(mode=="collection"){
	   change["out"]=run "/result";change["maxresidentrows"]="117316888";change["maxresidentbytes"]="45097156608"
	  }else{
	   change["families"]=collection "/families.tsv";change["chunks"]=collection "/chunks.tsv"
	   change["consensus"]=collection "/consensus.faa";change["expectedmembers"]="117316888"
	   if(mode=="prepare"){
	    change["lengthmedians"]=medians;change["familyhydro"]=hydro
	    if(length(blacklist))change["blacklist"]=blacklist
	    if(length(provenance))change["blacklistprovenance"]=provenance
	   }
	  }
	 }
	 /^#/||!length{print;next}
	 {
	  k=$1;if(index($0,"=")<2||seen[k]++){bad=1;exit 2}
	  if(k=="expectedmembers"&&substr($0,index($0,"=")+1)!="117316888"){bad=1;exit 2}
	  if(k in change){print k "=" change[k];replaced[k]++}else print
	 }
	 END{for(k in change)if(replaced[k]!=1)bad=1;exit bad}
	' "${option[parent]}" > "$out/native.config"
fi
command=("${option[java]}" -ea -XX:+ExitOnOutOfMemoryError -Xms64m "-Xmx$heap" "-XX:ActiveProcessorCount=$threads"
 "-Dbbgap.custom40.resource=${option[custom40]}" --add-modules jdk.incubator.vector -cp "$cpp")
case $mode in
	collection)command+=(prot.DdlClustCalibrationAdapter "config=$out/native.config");;
	lengths)command+=(benchmark.FamilyMemberLengthCensus "in=$out/native.config" insha80=unmeasured "out=$run/result" verification=structural);;
	hydro)
		command+=(prot.FamilyMemberScalars "families=${option[collection]}/families.tsv" familiessha80=unmeasured
		 "chunks=${option[collection]}/chunks.tsv" chunkssha80=unmeasured expectedfamilies=120000 expectedmembers=117316888
		 buckets=256 batches=1 shards=1 shard=0 minmemberaa=100 maxmemberaa=100 scalars=hydrophobicity t=4 ow=f
		 "out=$run/result/member_hydrophobicity.tsv" verification=structural);;
	prepare)
		command+=(prot.cluster.ClusteringLauncher target=calibrate "resourceroot=${option[resourceroot]}"
		 "prepare=$out/native.config" verification=structural targetfpr=0.0001 perclass=t rounds=4 t=64 threadsin=1 marginscale=1 "out=$run/result");;
	fit)
		command+=(prot.cluster.ClusteringLauncher "resourceroot=${option[resourceroot]}" "in=${option[prepared]}"
		 verification=structural marginscale=1 t=8 threadsin=1 target=calibrate
		 "neg=${option[negk5]},${option[negshadow]},${option[negreversed]},${option[negshuffled]}"
		 negnames=k5,shadow,reversed,shuffled negsha80=b54bd7c3f60261d9bbfa,4f8f156aa3a388d512a5,cdf5e374051d832abfcd,a3991495ce9695607a8c
		 expectedpanelqueries=1000000,1000000,1000000,1000000 targetacceptsbyclass=0,0,0,0 perclass=t rounds=4 "out=$run/result");;
esac
printf '%q ' "${command[@]}" > "$out/native.argv"
printf '\n' >> "$out/native.argv"
{
	printf '#!/usr/bin/env bash\nset -euo pipefail\nexport LC_ALL=C TZ=UTC\n'
	printf 'unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS CLASSPATH\n'
	printf 'export PATH=%q\n' "$runtime_path"
	printf ': "${SLURM_JOB_ID:?Run on an allocated Dori node}"\n'
	printf 'run=%q\n' "$run"
	printf '[[ ! -e $run && ! -L $run ]]\nmkdir -- "$run"\n'
	printf 'trap '\''rc=$?;printf "%%s\\n" "$rc" > "$run/exit"'\'' EXIT\n'
	[[ $mode != hydro ]] || printf 'mkdir -- "$run/result"\n'
	printf 'rc=0\n/usr/bin/time -f '\''%%e\t%%U\t%%S\t%%x'\'' -o "$run/native.time.tsv" '
	printf '%q ' "${command[@]}"
	printf '> "$run/native.stdout" 2> "$run/native.stderr" || rc=$?\n'
	printf 'printf "%%s\\n" "$rc" > "$run/native.exit"\nexit "$rc"\n'
} > "$out/native.command.sh"
chmod a-w -- "$out/classpath.txt" "$out/bindings.tsv" "$out/native.argv" "$out/native.command.sh"
[[ ! -e $out/native.config ]] || chmod a-w -- "$out/native.config"
for file in "$out/classpath.txt" "$out/bindings.tsv" "$out/native.argv" "$out/native.command.sh" "$out/native.config";do
	[[ -e $file ]] || continue
	bits=$(stat -c %a -- "$file")
	(( (8#$bits & 8#222)==0 )) || die "Output is writable: $file"
done
printf 'MATERIALIZED mode=%s control=%s run=%s science_executed=false outputs_validated=false\n' "$mode" "$out" "$run"
