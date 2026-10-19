#!/bin/bash
# Conservative standalone contig-tip polishing; not enabled by TadPipe.
usage(){
echo "
Written by Brian Bushnell and Fischl
Last modified October 18, 2026

Usage: tippolish.sh assembly=contigs.fa in=reads.fq.gz out=polished.fa
       tippolish.sh assembly=contigs.fa in=R1.fq.gz in2=R2.fq.gz out=polished.fa
       tippolish.sh assembly=contigs.fa in=paired.fq.gz interleaved=t out=polished.fa

Repairs supported substitutions near contig ends using unique assembly-end
overlaps as temporary mapping context. Does not join, extend, trim or remove
contigs, or edit indels. Experimental: validated on synthetic bacterial isolates,
not established safe for natural mixed-strain libraries or metagenomes.
Use adapter-trimmed, filtered, calibrated, preferably deduplicated reads BEFORE
error correction, merging or extension; inferred bases are not independent support.
No reference genome is needed. The input assembly remains the unpolished copy.

assembly=       Input FASTA (ref= is an alias); unique whitespace-free names.
in=             Input reads, single-end unless interleaved=t or in2 is given.
in2=            Optional second mate file; incompatible with interleaved=t.
interleaved=f   Set true for paired reads in one file.
out=            Output FASTA; .gz is supported.
overwrite=f     (ow) Allow replacement of final output/report; never input files.
report=         New edit ledger; default <out>.edits.tsv.
temp=           Existing work parent; default Java temporary directory.
t=auto          Threads per child. Stages run sequentially.
-Xmx8g          Heap per child, not the coordinator (fixed 256 MB).
childheap=8g    Alternative per-child heap setting; parsed config values supported.
minoverlap=33   Minimum candidate overlap; must be at least 33.
maxoverlap=200  Maximum candidate overlap.
indexed=t      Exact 16mer prefix index accelerates candidate discovery.
config=        Standard BBTools configuration file.

Fixed safeguards: 16 exact bases on each side of the candidate mismatch;
context support >=2, original mapped support <2, and no overlapping original
32mer supported >=2 times. Missing evidence abstains. Targets containing ambiguous
bases are retained unchanged. One-read true alleles can still resemble errors.

Whole-assembly mapping runs once; each context maps separately. This first driver
still counts read kmers per eligible target, so many contexts can be expensive.
All evidence, commands, per-stage timing and logs remain in a unique work directory,
on success or failure. No automatic cleanup. No stdin/stdout or comma-list inputs.
See docs/guides/TipContextGuide.txt for the experimental evidence and limitations.
";
}
if [[ $# == 0 || "$1" == -h || "$1" == --help || "$1" == help ]]; then usage; exit 0; fi
resolveSymlinks(){
	SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then CP="$DIR/bbtools.jar"
	else CP="$DIR/current/"; fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"
	parseJavaArgs --xmx=8g --xms=32m --mode=fixed "$@"
	setEnvironment
}

launch(){
	CMD="java $EA $EOOM $SIMD -Xmx256m -Xms32m -cp $CP consensus.TipPolish childheap=${XMX#-Xmx} $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD -Xmx256m -Xms32m -cp "$CP" consensus.TipPolish "bbtools=$DIR" "childheap=${XMX#-Xmx}" "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
