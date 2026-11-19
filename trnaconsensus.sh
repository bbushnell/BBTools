#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Neptune
Last modified August 11, 2026

Description:  Builds per-anticodon tRNA consensus sequences and HBM
(Hidden Brian Model) profiles from a collection of tRNA sequences.
Uses consensus-centroid clustering with reassignment for optimal
cluster representatives.  Output can be used with callgenes.sh
via the trnalib= and trnamodel= flags.

Usage:  trnaconsensus.sh in=trnas.fa out=consensus.fa outmodel=models.hbmc

File parameters:
in=<file>       Input fasta of tRNA sequences (with anticodon in header).
out=<file>      Output consensus fasta.
outmodel=<file> Output model: .hbmc for compact text, .hbm for legacy text,
                or .bgm for serialized graphs. Text formats also support .gz.

Clustering parameters:
cluster=t       Enable clustering (default true).
clusterid=0.75  Clustering identity threshold.
mincluster=3    Minimum cluster size to keep.
mingroup=10     Minimum sequences per anticodon group.
reassignrounds=2  Reassignment iterations; consensus is rebuilt between
                rounds, with early exit when no sequences move.
lentilt=0.03    Slight reassignment incentive toward the cluster whose median
                member length is closest to the sequence's own; 0 disables.
clusterpolicy=legacy  Optional nucleotide clustering: top50, patience or all.
                top50 chooses the best passing model in the first50 DDL ranks;
                only if none passes does it try every remaining model.
                patience walks DDL ranks; all exhausts every model. These three
                modes use Quantum for assignment, consensus and HBM traces.
clusterpatience=100  After the first pass, stop after this many alignments
                without a better passing identity;0 returns the first pass.
ddlbins=        Required512 or1024 for DDL policies; no selected default.
ddlseed=12345   Shared query/member/model DDL hash seed.
paddingblacklist=  Existing exported padding-word TSV; requires the complete
               256 SHORT8 and2048 MASKED10 prefixes for the selected fit_side.
                Use paddingblacklist=none explicitly for unfiltered initial
                DDL clustering; omit paddingfit in that mode. A missing file
                argument does not silently select an unfiltered sketch.
paddingfit=     Explicit literal fit_side label from that file. A historical
                fit-fold list is not an all-gene refit. No list learning occurs.
                The all policy needs no DDL/blacklist inputs and performs no
                DDL construction or indexing. DDL profiles include one vote per
                member/bin plus100 once per distinct consensus bin/key.
                New policies use lentilt=0 and stochastic=f; the legacy mode
                retains its existing settings and behavior.
outbuildstats=  Optional TSV for a new nucleotide policy with prune=0. Reports
                executed assignment, clustering-consensus, final-consensus and
                optional census alignments separately and in total. Also writes
                final model-size and member-identity histograms (identity rounded
                to nearest0.001) from existing final placements, without aligning
                again for reporting. Includes input/modeled/unmodeled members.
                With these policies, outassignments preserves absent source
                taxonomy as -1; the original full record ID remains authoritative.
                Native cmfrom/cmto/cmlen header fields additionally give each
                model's member-end ranges and joint left-missing44/right-missing20
                classes. Missing CM metadata remains explicitly counted.

Recruitment parameters:
recruit=t       Recruit orphans into clusters via consensus alignment.
recruitid=0.70  Minimum identity to consensus for recruitment.

Consensus parameters:
passes=2        Number of consensus refinement passes.
minid=0.3       Minimum alignment identity for consensus building.
endtrim=0.3     Truncate consensus termini at the modal member start/stop
                boundaries; 0 disables.
census=f        Print a per-model 5'-flushness census line (aligns members
                back to the final consensus; diagnostic).

Processing parameters:
t=              Worker threads for parallel cluster building (autodetected).

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM.  The max is typically 85% of
                physical memory.  The human genome requires around 24g minimum.

Please contact Brian at bbushnell@lbl.gov if you encounter any problems.
"
}

#This block allows symlinked shellscripts to correctly set classpath.
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  nextDir="$(readlink "$(basename "$DIR")")"
  DIR="$nextDir"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

resolveSymlinks(){
	if [ -d "$DIR""current/" ]; then
		CP="$DIR""current/"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=8g" "--xms=8g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prok.TrnaConsensusBuilder $@"
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" prok.TrnaConsensusBuilder "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
