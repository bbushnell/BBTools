#!/usr/bin/env bash
# Serial greedy protein clustering with parallel frozen passes.
usage(){
cat <<'HELP'
Written by Keqing, Yelan, and Fischl.

Clusters proteins using modal DDL sketches, monomer cosine ranking, and selected
substitution-matrix alignments (POSSUM55 by default). Greedy admission is serial;
frozen reassignment and per-family rebuilding can use multiple workers.
This is not a calibrated detector.
It loads proteins into memory and ranks live families exactly through an index
or the optional exhaustive reference scan.

Usage: ddlclust.sh in=tails.faa out=cluster_#.faa members=members.tsv.gz
                  consensus=consensus.faa.gz ddlout=families.pddl.gz

Saved-greedy structural arm (explicit opt-in):
  ddlclust.sh mode=cake in=tails.faa codecs=shard0.pddl.gz,shard1.pddl.gz
    outdir=new_arm familybudget=80000 initialfraction=1.5
    fractions=.25,.4375,.625,.8125,1 sizelambdas=.25,.1875,.125,.0625,0
    pricemode=floor belowmode=priorscaled passes=3 t=64

mode=cake          Select saved greedy families globally, assign all tails once
                   without rebuilding, run the specified layers, then select
                   up to familybudget positive families and refine. No negatives
                   or calibration occur inside this command. Normal mode is unchanged.
codecs=            Comma-separated greedy codec files with matching geometry.
outdir=            New directory; existing/partial runs are never overwritten.
familybudget=80000 Final maximum family count.
initialfraction=1.5 Initial count is floor(familybudget*initialfraction).
setupslotcost=auto  Constant/linear/sqrt fallback B0: auto uses the initial
                   post-assignment census; N sets a nonnegative integer.
                   Prior modes use the latest above-edge B instead.
fractions=         Required layer slot-cost fractions; none gives a no-cake control.
lambda=0.25        Constant size pressure; sizelambda is an alias.
sizelambdas=       One constant or exactly one nonnegative finite value per layer.
                   Do not combine with lambda/sizelambda. Zero is supported.
pricemode=floor    floor or fraction; belowmode=priorscaled, prior, constant,
                   linear or sqrt.
edgefraction=1.01  Floor-price population edge. Priorscaled updates the prior B
                   above the edge and scales it by active families / budget below.
pricedfraction=.99 Rank fraction when pricemode=fraction.
fissiontopfraction=1 Largest layer-start fraction eligible for fission. Below1
                   runs one frontier; 1 preserves recursive child fission.
fissiontopfractions= Cake-only: one constant or one value per layer, finite in [0,1].
                   Do not combine with fissiontopfraction. Zero skips fission in
                   that layer while preserving fusion and pricing. For example,
                   1,0,0,0,0 enables fission only in the first of five layers.
                   Each layer/run.config records its effective arguments.
fusionbottomfraction=1 Smallest post-fission fraction eligible for fusion pairs.
                   Both endpoints must qualify; ties use stable family IDs.
stopchangefraction=.001 Final refinement stops below this changed/input fraction
                   from pass2. passes=3 sets the maximum, with no interlayer refinement.
quality=auto       In cake mode, final/quality is the default quality prefix;
                   another prefix overrides it, null disables it. Initial assignment
                   never runs quality scoring. Final standard outputs and a layers.tsv
                   census retain initial/final coverage and each layer's settings.

in=                 Protein FASTA. Short proteins remain in the membership table.
out=                Optional per-family FASTA pattern with one # for stable ID.
members=            Optional TSV of every input ordinal and its final disposition.
consensus=          Optional consensus FASTA, ordered by final family size.
ddlout=             Optional uncalibrated sketch/consensus/statistics file.
seeds=              Optional existing ddlout model file for all-input reassignment.
                   Skips greedy and its initial rebuild. Requires passes>=1 and
                   identical sketch geometry/blacklist; imported member counts are
                   source metadata, never invented assignments to the new input.
                   No concurrentgreedy, scavenge, ramp, pairhist or decision audits.
                   Ordinary reassignment, pruning, rebuilding and quality reports
                   operate on actual new memberships. No new families are seeded.
stats=ddlclust.stats.tsv          Final counters and elapsed seconds.
refinement=ddlclust.refinement.tsv Consensus proposal outcomes and rejection reasons.
checkpoints=        Optional progress TSV with stage-local costs and sampled heap.
hist=               Optional shortlist-size and accepted-rank histograms.
assignmetrics=      Optional final-member identity/BLOSUM mean and median TSV.
                   Uses existing final-pass winner decisions, without realignment.
                   Requires passes>=1; includes only members retained after pruning
                   and topn. Values describe the last assignment reference, before
                   its final consensus rebuild. Empty results are explicitly NA.
quality=            Optional prefix for final-consensus cohesion and neighbor reports.
                   Reuses existing final rebuild alignments for per-family mean
                   identity and BLOSUM ratio. Searches ten shortlisted other final
                   consensuses for descriptive separation; this is approximate,
                   not a guaranteed global nearest neighbor. Reports a-b and
                   (a-b)/max(a,b), cohort quantiles and explicit missing values.
                   This output adds neighbor alignments, without changing clusters.
pairhist=           Optional prefix for decision-time candidate diagnostics:
                   .bins.tsv, .quantiles.tsv, .counts.tsv, .wrong_side.tsv.
                   Chosen assignments count tails; gate rejections count pairs.
                   Records DDL hits, cosine, scalar differences and winning rank
                   before updates. Quantiles are bin intervals, not exact scores.
                   No assignment cutoff is changed by these diagnostics.
decisionaudit=      Optional final-pass decision trace prefix (requires t=1,
                   passes>=1 and auditqueries=). Writes raw .queries.tsv,
                   .candidates.tsv and .references.tsv for verified packaging.
auditqueries=       Validated lost-query/proxy map from the membership reporter.
greedyaudit=        Optional serial-greedy decision trace prefix. Requires
                   auditqueries= and auditoldmembers=; width1024 and greedy hit
                   floors0-4, flex or ramp; shortlist must be at most500.
                   Rejects concurrentgreedy=t. Writes
                   .greedy_queries.tsv, .greedy_candidates.tsv and aggregate
                   .greedy_screens.tsv. Includes target hits, hit requirement,
                   occupied bins and live-family count even outside the shortlist.
                   Targets meeting the hit requirement get a rank capped at501,
                   where501 means greater than500. No extra alignment is performed.
                   Does not alter clustering decisions.
auditoldmembers=    Complete unfiltered membership table for live shared-member
                   correspondence; numeric family IDs alone are not equivalent.
                   Trace observations do not change ranking or gate decisions.
checkpointevery=10000 Progress reporting interval in eligible inputs.
progress=t         Emit and flush progress diagnostics at each checkpoint.
                   Index primitive bytes exclude object and temporary-build overhead.
                   With parallel frozen indexing, they describe row storage only;
                   the separate immutable query index is not included.

minlen=100          Minimum input protein length after removing a terminal '*'.
minspan=70          Minimum retained span on each side of an alignment.
mincluster=5        Prune smaller families after greedy and each reassignment.
prunefinal=f       With seeds=, retain all reference models between passes, even
                   those with zero assigned members. Empty models keep their
                   existing sequence and sketch until they regain members.
                   Apply mincluster only after the final reassignment, before
                   its final rebuild and topn selection. Requires passes>=1 and
                   merge=f. Default retains ordinary per-pass pruning.
topn=80000          Keep this many largest families at the end; 0 keeps all.
shortlist=50        Rank this many candidates; screen them before alignment.
indexed=t          Exact indexed ranking; f retains the exhaustive reference scan.
frozenindexthreads=1 Index-building workers for fixed-seed reassignment, in [1,t].
                   Values above one build the immutable index in parallel each
                   pass and skip mutable index construction. Requires seeds=,
                   indexed=t and merge=f. Default keeps the existing serial build.
concurrentgreedy=f Experimental asynchronous greedy admission using greedythreads and
                   one additional background index builder. Requires indexed=t,
                   pinright=t, and no decisionaudit. Results can vary by scheduling.
greedythreads=     Concurrent greedy workers; defaults to min(t,16). Explicit
                   values must be in [1,t] and require concurrentgreedy=t.
                   Frozen reassignment and family rebuilding still use t workers.
scavenge=f         After greedy, remove initial singleton families and replay their
                   tails against families of size at least two, using greedy gates.
                   Successful matches join existing families; failures stay
                   unassigned and never seed. Later reassignment still considers
                   every input. Runs before mincluster pruning, serially in seeded
                   input order. Reports rescued counts and separate stage costs.
                   Unsupported with pairhist or greedyaudit; defaults remain off.
snapshotfraction=0.02 Rebuild the immutable prefix index when recent families exceed
                   this fraction of its size; only used with concurrentgreedy=t.
merge=f            After each ordinary reassignment/prune/rebuild pass, cluster all
                   live consensuses and unite the original members of merged groups.
                   Each union gets one inclusive consensus rebuild. A rejected
                   proposal retains the seed and every member until reassignment;
                   it is not a claim that all members pass the retained seed.
                   Smallest source family ID survives. Pass zero is unchanged.
                   Defaults to the same alignment gate. mergeminid and
                   mergeminratio may override only consensus-vs-consensus
                   clustering; later assignment and union validation keep minid
                   and minblosumratio. Consensus inputs use seeded shuffle and
                   mincluster=1 with the same sketch geometry and hit/hydro screens.
                   Statistics and merge checkpoints report groups, counts and costs.
                   Nested model and union alignments count as construction work,
                   separate from original-input assignment attempts.
                   Requires passes>=1, pinright=t; incompatible with pairhist,
                   decision audits, assignmetrics and greedy ramp. Default off.
mergeminid=        Optional consensus-vs-consensus identity floor for merge=t.
                   Defaults to minid.
mergeminratio=     Optional consensus-vs-consensus BLOSUM ratio floor for merge=t.
                   Defaults to minblosumratio.
mergeroster=      Stage-B collect roster binding each consensus to its original
                   tail count. Requires positive mergemaxsize, passes=0,
                   mincluster=1, topn=0, all inputs, and no seeds/interpass merge,
                   scavenging, ramp or decision/pair diagnostics.
mergemaxsize=0    Restrict greedy consensus joins: at least one side must have
                   fewer than this many original tails. The target count grows
                   with each join and is rechecked atomically under concurrency.
                   For example, 20 allows 19+20 but forbids 20+20. Uses the ordinary
                   minid/minblosumratio gate; the subsequent original-tail union
                   has its own member gate. Zero disables this mode (no roster).
                   Shortlist counts/ranks in this mode describe size-eligible
                   candidates before alignment; merge_size_filtered and
                   merge_seed_size_filtered separately count size exclusions.
                   Lost size reservations try already-scored alternatives without
                   rerunning the query; merge_reservation_failures records them.
minhits=4          Minimum DDL hits before reassignment alignment.
passhits=4         Alias for minhits=N, or flex to apply the occupied-bin formula
                   during reassignment too. Default remains fixed4. If both
                   aliases are supplied, the last one wins. Summary minhits
                   records the effective number or flex. Flexible passes require
                   width1024 and do not support the fixed-floor decision audit.
greedyhits=flex    Initial greedy hit floor: floor(2+(filled reference bins+200)/512).
                   Requires width=1024; gives 2/3/4 hits at 0-311/312-823/824-1024 bins.
                   A nonnegative integer supplies a fixed floor; 0 disables it.
                   ramp starts one hit lower (1/2/3), then switches to flex once
                   rampfamilies qualifying families exist. Default remains flex.
rampfamilies=0     Required positive family target when greedyhits=ramp.
rampminmembers=1   A family qualifies once its published membership reaches this
                   size. One ramp state is captured per query, including seed
                   competition; concurrent queries may observe the switch later.
order=shuffle      Seeded random input order (unchanged default). hydro sorts by
                   ascending native hydrophobicity, retaining shuffled ties;
                   hydrobin groups bins of width0.05, shuffling within each bin.
                   Undefined hydrophobicity sorts last. Concurrent workers claim
                   this order but their admission completion order may differ.
maxhydrodelta=0.55 Maximum absolute query/consensus hydrophobicity difference in reassignment.
greedymaxhydrodelta=0.53 Same screen for initial greedy. Missing values leave this
                   screen unavailable; hit and alignment gates still apply.
                   Both caps accept finite values in [0,9]; 9 disables screening.
                   To recover unfiltered behavior set both hit floors=0 and both caps=9.
minid=0.35          Minimum identities / alignment columns, including gaps.
minblosumratio=0.14 Minimum selected-matrix score / query self-score.
matrix=POSSUM55     Alignment matrix for greedy, consensus rebuilding, structural
                   layers and reassignment. Without matrixfile, POSSUM55 loads
                   the installed resources/bbgap/POSSUM55.tsv resource.
                   Use matrix=blosum62 for the built-in BLOSUM62 scores.
matrixfile=         Optional versioned replacement TSV for POSSUM55; required for
                   external matrices without an installed default; forbidden with
                   blosum62. The matrix name stored in the TSV must match matrix=.
                   Its scores, scale and linear gap stay together.
                   This selects clustering scores, not downstream calibration.
denominator=aligned Use retained query suffix self-score; full uses all query AAs.
pinright=t         Both sequence ends must be consumed; terminal gaps are penalized.
                   With f, align the whole query to an unpinned reference interval.
trimleft=t         Permit left clipping during pinned alignment.
maxalignlen=2147483647 Explicit pinned alignment window cap; default is no truncation.
consensustrimleft=t Trim left consensus X/low-depth columns independently of alignment.
trimdepth=0        Consensus depth threshold as a fraction of maximum column depth.
                   Pinned mode never applies right-edge X/depth trimming.
passes=3           Number of frozen-reference, all-input reassignment rounds.
consensuspasses=2  Maximum internal consensus rounds after each retained-family pass.
                   Online updates propose one consensus at sizes 5,9,17,... .
seed=1             Random input-order seed; independent of sketchseed.
k=5 width=1024 sketchseed=0  DDL sketch configuration.
positional=        Defaults to t with pinright, f otherwise.
partitions=4 overlap=2       Positional K5/1024 stop-relative geometry.
blacklist=          Optional pentamer FASTA excluded from both family/query sketches.
blacklistsha80=     Optional expected final 20 hex digest; otherwise measured on load.

Sketches and member statistics refresh at sizes 1,2,4,8,... and after consensus
updates. Ranking is matching nonzero registers plus normalized monomer cosine.
Alignment ties use normalized score, identity, raw score, then stable family ID.
Hit and hydrophobicity screens do not change ranking or shortlist membership.
Assignment alignment counts include only actual alignments, not screened pairs.
pairhist also reports chosen matches to reference families of size >=5 and >=10:
pre-admission size in greedy, previous-pass frozen size during reassignment.
All assigned members contribute to a consensus proposal. A proposal replaces the
old consensus only if every member still passes the same gate; otherwise it is
rejected and recorded. Membership changes in greedy joins, optional singleton
scavenging, pruning and global reassignment.
All inputs, including members of pruned families, are reconsidered each pass.
Unmatched reassignment inputs remain unassigned; they do not seed new families.
Family IDs are stable and separate from final size ranks.

Concurrent greedy uses an immutable indexed prefix plus recent families and
may rank against older sketches while aligning to a current model. Regeneration
constructs from one fixed membership snapshot using positive right-pinned traces;
the proposed consensus must pass the original gate for every snapshot member.
Later joins remain members and never veto that proposal. Subsequent global
reassignment repairs membership. Diagnostic late_member_failures counts failed
checks, not unique failed genes; seed and version collision counters expose races.
Initial-shortlist histograms exclude competitive seed checks; separate seed_*
counters report those checks. All actual alignments still enter total work.

Original member FASTA headers, residue case, ambiguity symbols and terminal '*'
are preserved on output (line wrapping may change). Alignment folds B/J/Z to X
and strips a terminal '*'; internal stops or unsupported residues fail loudly.
Duplicate input names are distinguished by input ordinal.

config=file        Standard BBTools key=value config files are supported.
t=                 Workers for frozen reassignment and per-family rebuilding.
                   Sets the upper bound for greedythreads. The default
                   greedy path stays serial; the experimental one-worker path uses
                   the same asynchronous machinery as larger worker counts.
ow=t               Overwrite outputs. append mode and paired inputs are rejected.
-Xmx               Override Java memory allocation. -eoom exits on out-of-memory.
HELP
}
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then usage; exit 0; fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")" && pwd)
if [[ -f $DIR/bbtools.jar ]]; then CP=$DIR/bbtools.jar; else CP=$DIR/current; fi
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
parseJavaArgs --xmx=4g --xms=32m --mode=auto "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -ea -cp "$CP" prot.cluster.DdlClust "$@"
