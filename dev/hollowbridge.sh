#!/bin/bash

usage(){
echo "Experimental read-supported connections between unique contig-tip anchors.
Candidate report by default; joining requires explicit join=t. No Tadpole defaults change.

Usage: hollowbridge.sh contigs=assembly.fa in=reads.fq.gz out=candidates.tsv

contigs=        Input FASTA assembly; no reference is used for discovery.
in=             Primary read file. Each read supplies separate evidence.
in2=            Optional mate file. Mates never form a cross-read kmer.
extra=          Comma-delimited additional read files, each loaded once.
out=            Candidate TSV; indices and anchor starts are zero-based.
join=f          Enable a second read pass and conservative linear-path joining.
fillgaps=f      Accepted links become scaffolds with exactly10 Ns between full
                contigs, even for predicted overlaps. Ns mark unknown sequence,
                not a10bp gap estimate. Set t for experimental sequence filling.
outassembly=    Required FASTA destination when join=t; never overwrites inputs.
outjoins=       Required decision TSV when join=t, including rejected candidates.
gapratio=4      Winning whole-gap sequence needs this ratio over the runner-up.
flankcheck=f    Check observed read flanks OUTSIDE the selected span during joining.
flankbases=128  Maximum bases checked per outer flank; at most4096.
flankmin=16     Minimum usable bases per side; shorter or N-containing sides abstain.
k=224           Hollow span, greater than64; each arm is exactly32bp.
tipbases=192    Consider anchors wholly inside this many bases of either tip.
tipmin=0        Minimum seed distance from the near edge of its32bp arm to tip.
tipmax=-1       Maximum seed distance; -1 means tipbases-32. Both ends qualify.
slide=0         Audit the selected pair at every offset -slide through +slide.
                Requires full known-contig margin on both sides, supported
                expected pairs, and no supported competing partner either way.
minsupport=3    Minimum observations per arm pair; not independent molecules.
reads=-1        Maximum fragments PER input stream; -1 means all.
t=auto          Native I/O/compression threads; table counting is serial.
ow=f            Allow overwriting the report only when explicitly enabled.

Only assembly-wide-unique nonpalindromic anchors qualify. All supported
partners count toward ambiguity, including those absent from the assembly.
Both endpoints must have a sole supported partner. Competing end/gap
geometries are reported as conflicted, not selected. Support is the maximum
single pair count, never a sum of correlated overlapping anchor windows.
Ns invalidate a span, including its masked interior. Assembly uniqueness is
not proof of genomic uniqueness.

For example, tipmin=10 tipmax=160 slide=10 with k=224 leaves sliding room.
The inland band restricts SEEDS, not negative evidence: terminal, nonunique,
and palindromic probe words still retain read partners when slide>0.
A failed sliding check vetoes the candidate; it never retries a different seed.
The candidate report distinguishes branches, low support, Ns and short margins.
slide=0 retains the original behavior; all options remain experimental.

With join=t, the selected pair must be observed in at least minsupport read
records whose known flanks exactly match both contigs. The winning whole-gap
sequence must satisfy gapratio; no per-base mosaic consensus is manufactured.
More than16 gap variants vetoes the candidate. Repeated pair occurrences in
one read invalidate that record's vote. Records need not represent independent
molecules. Negative gaps require exact contig overlaps. Only unambiguous
linear paths are joined; cycles and rejected contigs remain unchanged. No
reference sequence enters this experimental joining decision.
With fillgaps=f, gap votes select links but never supply output bases or trim
overlaps. Preserve the original assembly and reports for later gap resolution.

With flankcheck=t, an anchored comparison allows one substitution, insertion
or deletion per outer flank. Conflicting records do not vote for gap sequence.
At least minsupport contradictory records veto the candidate unless confirmed
outer-flank records outnumber them by gapratio. Reads too short to test either
side do not count as confirmed. A repeated selected pair revokes both kinds of
evidence. No extra read pass or solid table is built. The join report appends
outer_confirmed and outer_contradictory counts; this safeguard is experimental.
"
}

if [[ $# == 0 || "$1" == -h || "$1" == --help ]]; then usage; exit 0; fi
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=2g --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.HollowBridge "$@"
