#!/usr/bin/env bash
# Materialize accepted Step03 all-tails construction configs and commands.
# Required: component=shared|full view=short|long phase=partition|seedN|shardN|final
# recipe=recipe_v1.tsv resources=resources_v1.config inputs=VIEW.inputs.tsv
# template=PHASE.config partition=ABS_DIR work=ABS_DIR out=NEW_CONTROL run=NEW_RUN
# java=ABS_JAVA or jdk=ABS_JDK_HOME, classpathfile=ONE_LINE_CPP,
# resourceroot=ABS_DIR, custom40=ABS_RESOURCE_MATRIX, matrixfile=ABS_RESOURCE_MATRIX.
# Optional: addmodules=COMMA_SEPARATED_JAVA_MODULES.
set -euo pipefail
export LC_ALL=C TZ=UTC

usage(){
	cat <<'EOF'
Usage: bbgap_baseline_construct_materialize_20261117.sh key=value ...

Write accepted Step03 construction configs and commands; do not run them.
All file and directory arguments must be absolute paths.

Required:
  component=shared|full view=short|long phase=partition|seedN|shardN|final
  recipe=FILE resources=FILE inputs=FILE template=FILE
  partition=DIR work=DIR out=NEW_CONTROL_DIR run=NEW_RUN_DIR
  java=EXECUTABLE or jdk=JDK_HOME
  classpathfile=FILE resourceroot=DIR custom40=FILE matrixfile=FILE

N is a zero-based partition ordinal from the recipe. The partition phase is
only for full/long. Use the matching VIEW_partition.config, VIEW_seedN.config,
VIEW_shard_cakeN.config or VIEW_final_cake.config template. The inputs file has
ordinal/path columns. The classpath file contains one line of absolute entries.
Both matrix arguments must match the matrix in the accepted resources config.
The control and run directories must be new and disjoint.

Optional:
  addmodules=MODULES    Comma-separated Java modules.
  toolpath=DIRS         Colon-separated absolute directories; /usr/bin:/bin
                       by default, after the selected Java bin directory.

Writes bindings.tsv, classpath.txt, native.config, native.argv and
native.command.sh under out. The generated command requires a Slurm allocation.
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
	case $key in
		component|view|recipe|resources|inputs|template|phase|partition|work|out|run|java|jdk|classpathfile|resourceroot|custom40|matrixfile|addmodules|toolpath);;
		*)die "Unknown option: $key";;
	esac
	[[ -n $value && ! ${option[$key]+present} && $value != *[$'\t\r\n\\']* ]] || die "Empty, duplicate or escaped option: $key"
	option[$key]=$value
done

required(){ [[ ${option[$1]+present} ]] || die "Required option: $1"; }
absolute(){ required "$1";[[ ${option[$1]} == /* ]] || die "Absolute path required: $1"; }
scalar(){
	local file=$1 key=$2 value
	value=$(awk -F '\t' -v key="$key" 'NR==1{if($0!="#final_all_tails_recipe\t1")bad=1;next}
		NR==2{if($0!="key\tvalue\tsource")bad=1;next}
		$1==key{if(NF!=3)bad=1;print $2;n++}
		END{exit bad||n!=1}' "$file") || die "Missing recipe key: $key"
	printf '%s' "$value"
}
config_value(){
	local file=$1 key=$2 value
	value=$(awk -F= -v key="$key" '$1==key{print substr($0,index($0,"=")+1);n++}END{exit n!=1}' "$file") ||
		die "Missing resource key: $key"
	printf '%s' "$value"
}
checked_uint(){
	local value=$1 min=$2 max=$3 label=$4 normalized
	[[ $value =~ ^0*([0-9]+)$ ]] || die "Bad $label"
	normalized=${BASH_REMATCH[1]}
	if [[ ${#normalized} -gt ${#max} || ( ${#normalized} -eq ${#max} && $normalized > $max ) ]];then
		die "Bad $label"
	fi
	normalized=$((10#$normalized))
	[[ $normalized -ge $min ]] || die "Bad $label"
	printf '%s' "$normalized"
}
phase=${option[phase]-}
case $phase in
	partition)phase_kind=partition;ordinal=0;;
	final)phase_kind=final;ordinal=0;;
	*)
		if [[ $phase =~ ^seed([0-9]+)$ ]];then phase_kind=seed;ordinal=${BASH_REMATCH[1]}
		elif [[ $phase =~ ^shard([0-9]+)$ ]];then phase_kind=shard;ordinal=${BASH_REMATCH[1]}
		else die 'Unknown phase'
		fi;;
esac

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
for key in recipe resources inputs template partition work out run java classpathfile resourceroot custom40 matrixfile;do absolute "$key";done
for key in recipe resources inputs template classpathfile custom40 matrixfile;do [[ -s ${option[$key]} ]] || die "Missing input file: $key";done
[[ -x ${option[java]} && -d ${option[resourceroot]} ]] || die 'Missing Java executable or resource root'
java_bin=${option[java]%/*}
[[ -n $java_bin ]] || java_bin=/
runtime_path=$java_bin:${option[toolpath]}

component=${option[component]-};view=${option[view]-}
case $component in shared|full);;*)die 'component must be shared or full';;esac
case $view in short|long);;*)die 'view must be short or long';;esac
[[ $(scalar "${option[recipe]}" experiment_arm) == "$component" ]] || die 'Recipe component mismatch'
case $view in short)[[ $(scalar "${option[recipe]}" views) == short ]] || die 'Recipe view mismatch';;long)[[ $(scalar "${option[recipe]}" views) == long ]] || die 'Recipe view mismatch';;esac
partition_ways=$(scalar "${option[recipe]}" partition_ways)
cluster_seed=$(scalar "${option[recipe]}" clustering_seeds)
sketch_seed=$(scalar "${option[recipe]}" sketchseed)
blacklist_words=$(scalar "${option[recipe]}" blacklistwords)
blacklist_k=$(scalar "${option[recipe]}" blacklistk)
blacklist_mode=$(scalar "${option[recipe]}" blacklistmode)
keytiebreak=$(scalar "${option[recipe]}" keytiebreak)
matrix_band=$(scalar "${option[recipe]}" matrix_band_policy)
resource_policy=$(scalar "${option[recipe]}" resource_policy)
alignmentmode=$(scalar "${option[recipe]}" alignmentmode)
pinright=$(scalar "${option[recipe]}" pinright)
trimleft=$(scalar "${option[recipe]}" trimleft)
denominator=$(scalar "${option[recipe]}" denominator)
consensustrimleft=$(scalar "${option[recipe]}" consensustrimleft)
trimdepth=$(scalar "${option[recipe]}" trimdepth)
token_re='^[A-Za-z0-9_.-]+$'
decimal_re='^((0|[1-9][0-9]*)([.][0-9]+)?|[.][0-9]+)$'
partition_ways=$(checked_uint "$partition_ways" 1 1024 'partition count')
if [[ $phase_kind == seed || $phase_kind == shard ]];then
	ordinal=$(checked_uint "$ordinal" 0 1023 'phase ordinal')
fi
[[ $ordinal -lt $partition_ways ]] || die 'Phase ordinal exceeds partition count'
cluster_seed=$(checked_uint "$cluster_seed" 0 2147483647 'cluster seed')
sketch_seed=$(checked_uint "$sketch_seed" 0 2147483647 'sketch seed')
[[ $blacklist_words =~ ^[1-9][0-9]*$ ]] || die 'Bad blacklist word count'
blacklist_k=$(checked_uint "$blacklist_k" 1 31 'blacklist K')
[[ $blacklist_mode =~ $token_re ]] || die 'Bad blacklist mode'
[[ $keytiebreak =~ $token_re ]] || die 'Bad key tiebreak'
[[ $resource_policy =~ $token_re ]] || die 'Bad resource policy'
[[ $matrix_band =~ $token_re ]] || die 'Bad matrix band'
[[ $alignmentmode =~ $token_re ]] || die 'Bad alignment mode'
[[ $pinright =~ ^[tf]$ && $trimleft =~ ^[tf]$ && $denominator =~ $token_re &&
	$consensustrimleft =~ ^[tf]$ && $trimdepth =~ ^[0-9]+$ ]] || die 'Bad trim gate'
[[ $(scalar "${option[recipe]}" construction_ratio) =~ $decimal_re ]] || die 'Bad construction ratio'
[[ $(config_value "${option[resources]}" schema) == alltails_resources_v1 ]] || die 'Unsupported resource schema'
[[ $(config_value "${option[resources]}" resource_acceptance) == ACCEPTED ]] || die 'Resources are not accepted'
[[ $(config_value "${option[resources]}" matrix_name) =~ $token_re ]] || die 'Bad matrix name'
[[ $(config_value "${option[resources]}" matrix_band) == "$matrix_band" ]] || die 'Resource matrix band mismatch'

case $phase_kind in
	partition)
		[[ $component == full && $view == long ]] || die 'Only full/long has a public partition phase'
		expected_template=${view}_partition.config;cpus=$(scalar "${option[recipe]}" cake_threads);heap=$(scalar "${option[recipe]}" cake_heap);;
	seed)
		expected_template=${view}_seed${ordinal}.config;cpus=$(scalar "${option[recipe]}" greedy_threads);heap=$(scalar "${option[recipe]}" greedy_heap);;
	shard)
		expected_template=${view}_shard_cake${ordinal}.config;cpus=$(scalar "${option[recipe]}" cake_threads);heap=$(scalar "${option[recipe]}" cake_heap);;
	final)
		expected_template=${view}_final_cake.config;cpus=$(scalar "${option[recipe]}" cake_threads);heap=$(scalar "${option[recipe]}" cake_heap);;
esac
[[ ${option[template]##*/} == "$expected_template" ]] || die "Template must be $expected_template"
cpus=$(checked_uint "$cpus" 1 64 'thread count')
[[ $heap =~ ^[0-9]+[gGmMkK]?$ ]] || die 'Bad heap'

out=${option[out]};run=${option[run]}
[[ ! -e $out && ! -L $out && ! -e $run && ! -L $run ]] || die 'Control and run must both be fresh paths'
physical_out=$(realpath -m -- "$out") || die 'Cannot resolve control path'
physical_run=$(realpath -m -- "$run") || die 'Cannot resolve run path'
[[ ! -e $physical_out && ! -L $physical_out && ! -e $physical_run && ! -L $physical_run ]] || die 'Control and run must both resolve to fresh paths'
[[ $physical_out != "$physical_run" && $physical_out != "$physical_run/"* && $physical_run != "$physical_out/"* ]] ||
	die 'Control and run paths must resolve to disjoint destinations'

mapfile -t cpp_lines < "${option[classpathfile]}"
[[ ${#cpp_lines[@]} == 1 ]] || die 'Classpath file must contain exactly one line'
cpp=${cpp_lines[0]}
[[ -n $cpp && $cpp != :* && $cpp != *: && $cpp != *::* && $cpp != *[$'\r\n']* ]] || die 'Malformed classpath'
IFS=: read -r -a cpp_entries <<< "$cpp"
for entry in "${cpp_entries[@]}";do [[ $entry == /* && -e $entry ]] || die "Missing absolute classpath entry: $entry";done

minpaired=$(scalar "${option[recipe]}" construction_${view}_minpaired)
ratio=$(scalar "${option[recipe]}" construction_ratio)
tailcrop=$(scalar "${option[recipe]}" ${view}_tailcrop)
blacklist=$(config_value "${option[resources]}" blacklist)
alphabet=$(config_value "${option[resources]}" alphabet)
matrix=$(config_value "${option[resources]}" matrix)
matrix_name=$(config_value "${option[resources]}" matrix_name)
[[ ${option[custom40]} == "$matrix" && ${option[matrixfile]} == "$matrix" ]] || die 'Matrix arguments must match accepted resources'
[[ -s $blacklist && -s $alphabet && -s $matrix ]] || die 'Missing alphabet, blacklist or matrix resource'
if [[ ${option[addmodules]+present} ]];then
	[[ ${option[addmodules]} =~ ^[A-Za-z0-9_.]+(,[A-Za-z0-9_.]+)*$ ]] || die 'Bad Java module list'
fi

mkdir -p -- "${out%/*}"
mkdir -- "$out"
cp -- "${option[classpathfile]}" "$out/classpath.txt"
printf 'format\tbbgap_step03_materialization_v1\n' > "$out/bindings.tsv"
for key in component view recipe resources inputs template phase partition work out run java jdk classpathfile resourceroot custom40 matrixfile addmodules toolpath;do
	[[ ! ${option[$key]+present} ]] || printf '%s\t%s\n' "$key" "${option[$key]}" >> "$out/bindings.tsv"
done

awk -F= -v OFS='=' -v phase="$phase" -v kind="$phase_kind" -v shard="$ordinal" -v run="$run" \
 -v work="${option[work]}" -v part="${option[partition]}" -v blacklist="$blacklist" -v cpus="$cpus" \
 -v ways="$partition_ways" -v seed="$cluster_seed" -v sketch="$sketch_seed" \
 -v words="$blacklist_words" -v bk="$blacklist_k" -v bmode="$blacklist_mode" \
 -v kt="$keytiebreak" -v align="$alignmentmode" -v pin="$pinright" -v tl="$trimleft" -v denom="$denominator" \
 -v ctl="$consensustrimleft" -v td="$trimdepth" \
 -v provenance="policy=$resource_policy;component=$component;binding=${option[resources]}" \
 -v matrix_name="$matrix_name" -v ratio="$ratio" -v paired="$minpaired" -v crop="$tailcrop" '
 FILENAME==ARGV[1]{
  if(FNR==1){if($0!="ordinal\tpath")bad=1;next}
  n=split($0,p,"\t");if(n!=2||p[1]!=count||p[2]!~/^\//||p[2]~/[,\r\n]/)bad=1
  vector=vector (count?",":"") p[2];count++;next
 }
 FILENAME==ARGV[2]{
  if($1!="aabits"||++maps!=1)bad=1;map=substr($0,index($0,"=")+1);next
 }
 /^#/ || /^$/ {next}
 {
  key=$1
  if(key=="checkpointevery"){
   ckpt=substr($0,index($0,"=")+1);if(ckpt!~/^[0-9]+$/||ckpt<1)bad=1;ckptcount++;next
  }
  if(key=="ways"){
   value=substr($0,index($0,"=")+1);if(kind!="partition"||value!=ways)bad=1;wayscount++;next
  }
  if(key=="blacklistsha80"||key=="outdir"||key=="out"||key=="map"||key=="ready"||
     key=="ddlout"||key=="members"||key=="stats"||key=="consensus"||key=="refinement"||
     key=="checkpoints"||key=="hist"||key=="pairhist"||key=="codecs"||key=="keytiebreak"||
     key=="aabits"||key=="blacklist"||key=="blacklistwords"||key=="verification"||
     key=="blacklistk"||key=="blacklistmode"||key=="blacklistprovenance"||key=="matrix"||
     key=="matrixfile"||key=="matrixfilesha80"||key=="tailcrop"||key=="alphabet"||key=="minblosumratio"||
     key=="t"||key=="greedythreads"||key=="frozenindexthreads"||key=="checkpointevery"||
     key=="alignmentmode"||key=="pinright"||key=="minpaired"||key=="trimleft"||key=="ways"||
     key=="denominator"||key=="consensustrimleft"||key=="trimdepth"){next}
  if(key=="in"){
   if(kind=="seed"||kind=="shard"){print key,part"/part"shard".faa.gz"}else{print key,vector}
   seenin++;next
  }
  if(key=="seed"){if(substr($0,index($0,"=")+1)!=seed)bad=1;seedcount++}
  if(key=="sketchseed"){if(substr($0,index($0,"=")+1)!=sketch)bad=1;sketchcount++}
  print
 }
 END{
  if(bad||count<1||(kind!="partition"&&maps!=1)||seenin!=1||seedcount!=1||(kind!="partition"&&sketchcount!=1)||
     (kind=="seed"&&ckptcount!=1)||(kind!="seed"&&ckptcount>0)||
     (kind=="partition"&&wayscount!=1)||(kind!="partition"&&wayscount>0))exit 2
  if(kind=="partition"){
   print "ways",ways
   print "out",run"/part#.faa.gz";print "map",run"/partition.tsv.gz";print "ready",run"/part#.ready"
   if(crop>0)print "tailcrop",crop
  }else{
   print "minpaired",paired;print "minblosumratio",ratio;print "matrix",matrix_name;print "keytiebreak",kt
   print "alignmentmode",align
   print "pinright",pin;print "trimleft",tl;print "denominator",denom;print "consensustrimleft",ctl;print "trimdepth",td
   print "aabits",map;print "blacklist",blacklist;print "blacklistsha80","unmeasured"
   print "blacklistwords",words;print "blacklistk",bk;print "blacklistmode",bmode;print "blacklistprovenance",provenance
   print "verification","structural";print "t",cpus
   if(crop>0)print "tailcrop",crop
   if(kind=="seed"){
    print "greedythreads",cpus;print "checkpointevery",ckpt;print "out","null"
    print "ddlout",run"/result/families.pddl.gz";print "members",run"/result/members.tsv.gz";print "stats",run"/result/stats.tsv"
    print "consensus",run"/result/consensus.faa.gz";print "refinement",run"/result/refinement.tsv.gz"
    print "checkpoints",run"/result/checkpoints.tsv.gz";print "hist",run"/result/hist.tsv.gz";print "pairhist",run"/result/pairs"
   }else{
    print "frozenindexthreads",cpus;print "outdir",run"/result"
    if(kind=="final"){
     for(i=0;i<ways;i++){c=c (i?",":"") work"/shard"i"/result/final/families.pddl.gz"}
    }else{c=work"/seed"shard"/result/families.pddl.gz"}
    print "codecs",c
   }
  }
 }' "${option[inputs]}" "$alphabet" "${option[template]}" > "$out/native.config"

command=("${option[java]}" -ea -XX:+ExitOnOutOfMemoryError -Xms64m "-Xmx$heap" "-XX:ActiveProcessorCount=$cpus"
 "-Dbbgap.custom40.resource=${option[custom40]}")
[[ ! ${option[addmodules]+present} ]] || command+=(--add-modules "${option[addmodules]}")
command+=(-cp "$cpp")
if [[ $phase_kind == partition ]];then
	command+=(prot.cluster.DdlClustPartition "config=$out/native.config")
else
	command+=(prot.cluster.ClusteringLauncher target=cluster "resourceroot=${option[resourceroot]}"
		"config=$out/native.config" "matrixfile=${option[matrixfile]}")
	[[ $phase_kind == seed ]] || command+=(mode=cake)
fi
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
	[[ $phase_kind == partition || $phase_kind == shard || $phase_kind == final ]] || printf 'mkdir -- "$run/result"\n'
	printf 'rc=0\n/usr/bin/time -f '\''%%e\t%%U\t%%S\t%%x'\'' -o "$run/native.time.tsv" '
	printf '%q ' "${command[@]}"
	printf '> "$run/native.stdout" 2> "$run/native.stderr" || rc=$?\n'
	printf 'printf "%%s\\n" "$rc" > "$run/native.exit"\nexit "$rc"\n'
} > "$out/native.command.sh"
chmod a-w -- "$out/classpath.txt" "$out/bindings.tsv" "$out/native.config" "$out/native.argv" "$out/native.command.sh"
for file in "$out/classpath.txt" "$out/bindings.tsv" "$out/native.config" "$out/native.argv" "$out/native.command.sh";do
	bits=$(stat -c %a -- "$file")
	(( (8#$bits & 8#222)==0 )) || die "Output is writable: $file"
done
printf 'MATERIALIZED_STEP03 component=%s view=%s phase=%s control=%s run=%s science_executed=false\n' "$component" "$view" "$phase" "$out" "$run"
