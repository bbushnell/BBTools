#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 12, 2026

Description:  Assembles preprocessed paired Illumina reads with deduplication,
overlap/Clumpify/Tadpole correction, staged merging, extension, and ordered
multi-K Tadpole assembly. Inputs must already be adapter-trimmed, filtered,
and quality-score calibrated. Each stage runs in a separate JVM.
The starting recipe targets roughly 150bp bacterial isolate libraries;
shorter reads or unusual coverage require adjusted stage parameters.

Usage:
tadpipe.sh in=reads_R1.fq.gz in2=reads_R2.fq.gz out=contigs.fa -Xmx16g t=8
tadpipe.sh in=interleaved.fq.gz out=contigs.fa -Xmx16g t=8
tadpipe.sh in=interleaved.fq.gz out=contigs.fa k=96,124,300,64,32 assemblek=96


Parameters:
in=<file>           Paired input reads; interleaved when in2 is omitted.
in2=<file>          Optional read 2, if reads are in two files.
out=contigs.fa      Output file name.
temp=<directory>    Parent for a unique work directory (default Java tmpdir).
delete=t            Delete pipeline FASTQs only after successful assembly.
                    Commands and per-stage logs are always retained.
                    All intermediates are retained on failure.
gz=t                Compress intermediates through native BBTools I/O.
overwrite=f         Allow replacement of an existing final output.
dryrun=f            Validate paths and print/save commands without running them.
t=<auto>            Threads per child. Stages run serially to reuse the heap
                    budget and avoid keeping two kmer tables in memory.
-Xmx16g             Maximum heap PER CHILD. Coordinator uses only 256 MB.
childheap=<auto>    Alternative explicit child heap (e.g. childheap=16g).
zl=4                Intermediate compression level.

Stage switches (all default true):
dedupe=t            Sequence-based paired deduplication, NOT optical-only.
                    Both mates must match under Clumpify's duplicate criteria.
ecco=t              BBMerge overlap correction.
clump=t             Clumpify correction of both mates, with pairing restored.
ecc=t               Tadpole kmer correction at K62.
merge=t             Overlap merging, then optional REM passes.
rem=t               Enable the remk schedule after ordinary overlap merging.
qtrim=t             Quality-trim/filter the residual unmerged reads.
extend=t            Extend merged and unmerged reads with both sets as evidence.
extendm=t           Extend merged reads; requires extend=t and merge=t.
extendu=t           Extend unmerged reads; requires extend=t.
nn=t                Use the bundled fusion NN at cutoff 0.667098.

REM schedule controls (do not change ECC, extension, or final assembly K):
remk=124,145,93      Ordered REM K list; arbitrary positive K values and repeats.
                    remk=none or rem=f retains overlap merging only.
remextend=auto       Extension limit per REM pass. One integer applies to all,
                    or give one comma-separated integer per remk entry.
                    Zero is allowed. Auto uses 120 for K124, 140 for K145,
                    and 100 for other K values, preserving the old recipe.
                    Each pass consumes residual pairs and earlier merged reads.
                    Repeated K passes have distinct outputs and logs.
The scheduled K93 pass starts with strict=t; other REM passes use native
stringency unless overridden. merge_= all merging, merge<K>_= scheduled K,
rem<N>_= one-based REM pass (excluding overlap). More specific overrides win:
remextend defaults, then merge_, then merge<K>_, then rem<N>_.

For example, remk=64,96,64 remextend=40,80,40 rem3_strict=t runs three REM
passes, making only the last pass strict. These are controls, not tuned presets.

Final assembly controls (never applied to correction, merging, or extension):
k=124,300,96,64,32  Ordered final assembly/traversal schedule; repeats preserved.
assemblek=124       Initial contigging K for the final assembly.
graphk=<last K>     Final graph K; defaults to the last entry in k.
bridgek=<auto>      Override the eligible bridging K list.
fusek=<auto>        Override the eligible fusion K list.
These are aliases for assemble_k, assemble_assemblek, assemble_graphk,
assemble_bridgek, and assemble_fusek. Last occurrence wins across aliases.
Changing k alone keeps assemblek=124; set both to change the initial K.

Coverage guards are explicitly set to 1.75; they can reject genuine joins
across uneven coverage. Use assemble_fusecoverageratio=0 and/or
assemble_graphmergecoverageratio=0 to disable them independently.
For reads lacking long merged/extended support, use k=124,96,64,32.

Other parameters can be passed to individual phases like this:

assemble_k=124,96,64,32   Set the ORDERED assembly/traversal schedule.
assemble_assemblek=124   Set initial contigging K (also accepts assemblek=124).
assemble_mincontig=500   Override final assembly minimum contig length.
assemble_mcs=2           Override final assembly minimum seed kmer count.
assemble_mce=1           Override final assembly minimum extension kmer count.
assemble_prefilter=t     Enable the final assembly prefilter.
dedupe_subs=0            Require exact duplicates (native default permits 2).
clump_passes=6           Set Clumpify correction passes (pipeline default 4).
ecco_k=31               Set overlap-correction kmer filtering K (default 31).
ecc_k=62                Set Tadpole correction K (default 62).
ecc_wash=t              Clean the correction graph (pipeline default true).
ecc_ldf=0.15            Set the low-depth-read discard fraction (default 0.15).
ecc_tossdepth=1         Set the low-depth kmer threshold (pipeline default 1).
                        Use -1 to disable this low-depth-read discard test.
ecc_rbb=t               Require both mates to be bad before pair rejection
                        (native default false).
ecc_tu=f                Keep uncorrectable reads (pipeline default true).
ecc_tossjunk=t          Enable native junk-read rejection (default false).
merge_strict=t           Apply to all merging passes, not overlap correction.
merge145_extend2=100     Override only the K145 REM pass.
rem2_strict=f            Override stringency of only the second REM pass.
qtrim_trimq=15          Residual-read trimming quality (pipeline default 15).
qtrim_maq=14            Residual-read minimum mean quality (default 14).
qtrim_minlen=90         Residual-read minimum length (default 90).
extend_mce=4            Minimum extension kmer count for both sets (default 4).
extend_el=20             Set left extension for both read sets.
extendu_k=96             Override only unmerged-read extension.
Merged extension defaults: extendm_k=145 extendm_el=10 extendm_er=10.
Unmerged defaults: extendu_k=124 extendu_el=10 extendu_er=0.
Each extension sees the unextended opposite read set, even if that set's
extension is disabled. extend=f disables both sets regardless of their switches.

Valid prefixes:

dedupe_, ecco_, clump_ (clumpify_), ecc_ (correct_),
merge_, merge0_ (overlap), merge<K>_, rem<N>_, qtrim_,
extend_ (extend1_), extendm_, extendu_, assemble_.

Stage-specific overrides follow group-wide overrides. TadPipe owns filenames
and pairing; phase-prefixed in/out/extra/config/overwrite options are rejected.
Use assemble_<Tadpole option>=<value> for other final assembly controls;
unknown options are rejected by the receiving tool, not silently ignored.
Use ecc_k, extendm_k, or extendu_k to adjust those earlier phases explicitly.
Custom assemble_fusenet requires an explicit assemble_fusencutoff.
Adapter trimming, contaminant filtering, and the old unreachable extend2_
stage are not part of this pipeline. Named files and Bash are required.

Workflow experiments can be stored as ordinary config=files, one flag=value
per line; e.g. config=low_depth.config. Keep JVM heap flags on the command line.
Inspect dryrun=t first. No validated high/mid/low/metagenome presets or automatic
profile selection are supplied yet; current no-override defaults are unchanged.
"
}

if [ -z "$1" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	usage
	exit
fi

resolveSymlinks(){
	SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=14g" "--xms=14g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD -Xmx256m -Xms32m -cp $CP assemble.TadPipe childheap=${XMX#-Xmx} $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD -Xmx256m -Xms32m -cp "$CP" assemble.TadPipe "bbtools=$DIR" "childheap=${XMX#-Xmx}" "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
