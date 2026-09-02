#!/bin/bash

set -e

resolveSymlinks(){
	SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

runTest(){
	local class="$1"
	java -ea --add-modules jdk.incubator.vector -cp "$CP" "$class"
}

runCrossKLeftBridgeTest(){
	local temp expected actual rc
	temp="$(mktemp -d)"
	trap 'rm -rf "$temp"' RETURN
	expected="CATACTCGCCTCGGCATATGTAGCCCGGAATCATCACACTCCTCAGGGAGCTCCTACTAGCTTGGAGCTCATCGCCGTTCACTCGGACTAACTCTGCAGCAAATGACGCCCGACCACGGTCCCAACTGACATCAGCCGCTTTGATGCGAAGCCAATTTCATTCCTCCTCTATCCTAAGCATCGATCTCTCCAACATCTCGACTCTTGGAG"
	rc="$(printf '%s' "$expected" | tr ACGT TGCA | rev)"
	for spec in "crossk_left_bridge_reads.fa:60,40:Tadpole2" \
		"crossk_left_bridge_reads_k31.fa:31,20:Tadpole1"; do
		local input="${spec%%:*}"
		local rest="${spec#*:}"
		local kmers="${rest%%:*}"
		local implementation="${rest##*:}"
		"$DIR/tadpolemulti.sh" in="$DIR/testdata/$input" out="$temp/out.fa.gz" \
			k="$kmers" mcs=1 mce=1 mincontig=1 sweeplen=0 prefilter=f pop=f ow=t t=1 showstats=f \
			1>"$temp/stdout" 2>"$temp/stderr"
		actual="$(zcat "$temp/out.fa.gz" | awk '!/^>/{printf "%s", $0} END{print ""}')"
		if [ "$actual" != "$expected" ] && [ "$actual" != "$rc" ]; then
			echo "FAIL: $implementation reciprocal cross-k left bridge did not reconstruct the expected sequence" >&2
			cat "$temp/stderr" >&2
			exit 1
		fi
		echo "PASS: reciprocalCrossKLeftBridge$implementation"
	done
}

runLongKBridgeEligibilityTest(){
	local temp
	temp="$(mktemp -d)"
	trap 'rm -rf "$temp"' RETURN
	"$DIR/tadpole.sh" in="$DIR/testdata/crossk_long_bridge_short_contig.fa" out="$temp/out.fa.gz" \
		assemblek=31 fusek=none bridgek=60,20 mcs=1 mce=1 mincontig=1 sweeplen=0 prefilter=f pop=f \
		ow=t t=1 showstats=f verbose=t 1>"$temp/stdout" 2>"$temp/stderr"
	if ! grep -Fq "Cross-k endpoints: eligible=2" "$temp/stderr"; then
		echo "FAIL: long bridge K did not exclude the short contig safely" >&2
		cat "$temp/stderr" >&2
		exit 1
	fi
	if ! grep -Fq "Cross-k endpoints: eligible=4" "$temp/stderr"; then
		echo "FAIL: later short bridge K lost the previously ineligible endpoints" >&2
		cat "$temp/stderr" >&2
		exit 1
	fi
	echo "PASS: longKBridgeSkipsShortContigsWithoutLosingEndpoints"
}

runLongKBridgeResolutionTest(){
	local temp before after
	temp="$(mktemp -d)"
	trap 'rm -rf "$temp"' RETURN
	"$DIR/tadpole.sh" in="$DIR/testdata/crossk_long_bridge_x.fa" out="$temp/baseline.fa.gz" \
		k=31 mcs=1 mce=1 mincontig=1 sweeplen=0 prefilter=f pop=f ow=t t=1 showstats=f \
		1>"$temp/baseline.stdout" 2>"$temp/baseline.stderr"
	"$DIR/tadpole.sh" in="$DIR/testdata/crossk_long_bridge_x.fa" out="$temp/bridged.fa.gz" \
		assemblek=31 fusek=none bridgek=60 mcs=1 mce=1 mincontig=1 sweeplen=0 prefilter=f pop=f \
		ow=t t=1 showstats=f 1>"$temp/bridged.stdout" 2>"$temp/bridged.stderr"
	before="$(zgrep -c '^>' "$temp/baseline.fa.gz")"
	after="$(zgrep -c '^>' "$temp/bridged.fa.gz")"
	if [ "$before" != 5 ] || [ "$after" != 3 ]; then
		echo "FAIL: long-K branch resolution produced $before -> $after contigs instead of 5 -> 3" >&2
		cat "$temp/bridged.stderr" >&2
		exit 1
	fi
	echo "PASS: longKBridgeResolvesShortKBranch"
}

runLowDepthFateTrackingTest(){
	local temp expected
	temp="$(mktemp -d)"
	trap 'rm -rf "$temp"' RETURN
	"$DIR/tadpole.sh" in="$DIR/testdata/crossk_left_bridge_reads.fa" out="$temp/out.fa.gz" \
		k=60,40 mcs=1 mce=1 mincontig=1 sweeplen=0 prefilter=f pop=f lowdepthcontigdiag=both \
		ldcmaxcov=100 ldcfrac=1 ldcmaxlen=1000 ow=t t=1 showstats=f \
		1>"$temp/stdout" 2>"$temp/stderr"
	expected="Early low-depth candidate fates: initial=2/200 contigs/bases, present=1, retiredOrAbsorbed=1, grown=1/+110 bases, stillCandidates=1, reclassified=0."
	if ! grep -Fqx "$expected" "$temp/stderr"; then
		echo "FAIL: early low-depth candidate fate accounting was incorrect" >&2
		cat "$temp/stderr" >&2
		exit 1
	fi
	echo "PASS: earlyLowDepthCandidateFatesTrackCrossKMerge"
}

runShortContigRetentionTest(){
	local temp common
	temp="$(mktemp -d)"
	trap 'rm -rf "$temp"' RETURN
	common="in=$DIR/testdata/crossk_left_bridge_reads.fa k=60 mincontig=500 sweeplen=0 mcs=1 mce=1 pop=t lowdepthcontigdiag=early ldcmaxcov=100 ldcfrac=1 ldcmaxlen=1000 showstats=f ow=t t=1"
	"$DIR/tadpole.sh" $common out="$temp/keep.fa.gz" retainshortcontigs=t \
		1>"$temp/keep.stdout" 2>"$temp/keep.stderr"
	"$DIR/tadpole.sh" $common out="$temp/drop.fa.gz" retainshortcontigs=f \
		1>"$temp/drop.stdout" 2>"$temp/drop.stderr"
	"$DIR/tadpole.sh" $common out="$temp/evict.fa.gz" evictlowdepthcontigs=t \
		1>"$temp/evict.stdout" 2>"$temp/evict.stderr"
	"$DIR/tadpole.sh" $common out="$temp/classdrop.fa.gz" graphclassify=t retainshortcontigs=f \
		1>"$temp/classdrop.stdout" 2>"$temp/classdrop.stderr"
	"$DIR/tadpole.sh" in="$DIR/testdata/crossk_left_bridge_reads.fa" out="$temp/multidrop.fa.gz" \
		assemblek=60 fusek=none bridgek=none graphk=60 mincontig=500 sweeplen=0 mcs=1 mce=1 pop=t \
		graphclassify=t retainshortcontigs=f lowdepthcontigdiag=early ldcmaxcov=100 ldcfrac=1 \
		ldcmaxlen=1000 showstats=f ow=t t=1 1>"$temp/multidrop.stdout" 2>"$temp/multidrop.stderr"
	if ! grep -Fq 'Candidate filters (assemble-k): live=2,' "$temp/keep.stderr" \
			|| ! grep -Fq 'skipped; no positive contig coverage' "$temp/drop.stderr" \
			|| ! grep -Fq 'Low-depth contigs evicted (pre-graph): 2 contigs, 200 bases; retained=0.' "$temp/evict.stderr" \
			|| ! grep -Fq 'retainshortcontigs=false' "$temp/classdrop.stderr" \
			|| ! grep -Fq 'skipped; no positive contig coverage' "$temp/classdrop.stderr" \
			|| ! grep -Fq 'retainshortcontigs=false' "$temp/multidrop.stderr" \
			|| ! grep -Fq 'skipped; no positive contig coverage' "$temp/multidrop.stderr"; then
		echo "FAIL: short-contig retain/drop/selective-eviction controls diverged" >&2
		cat "$temp/keep.stderr" "$temp/drop.stderr" "$temp/evict.stderr" \
			"$temp/classdrop.stderr" "$temp/multidrop.stderr" >&2
		exit 1
	fi
	echo "PASS: shortContigRetentionAndEvictionControls"
}

runHashedBridgeParityTest(){
	local temp
	temp="$(mktemp -d)"
	trap 'rm -rf "$temp"' RETURN
	for mode in t f; do
		"$DIR/tadpole.sh" in="$DIR/testdata/crossk_long_bridge_x.fa" out="$temp/$mode.fa" \
			assemblek=31 fusek=none bridgek=95 bridgehash=$mode maskcore=t mcs=1 mce=1 mincontig=1 sweeplen=0 \
			prefilter=f pop=f ow=t t=1 showstats=f 1>"$temp/$mode.stdout" 2>"$temp/$mode.stderr"
	done
	if ! cmp -s "$temp/t.fa" "$temp/f.fa"; then
		echo "FAIL: fixed-width and explicit k=95 bridge tables emitted different contigs" >&2
		diff -u "$temp/f.fa" "$temp/t.fa" >&2 || true
		exit 1
	fi
	if ! grep -Eq '^hashkmers[[:space:]]+pair$' "$temp/t.stderr"; then
		echo "FAIL: compact bridge-table execution plan was not reported" >&2
		cat "$temp/t.stderr" >&2
		exit 1
	fi
	echo "PASS: hashedBridgeMatchesExplicitK95"
}

runHashedSequenceModeParityTest(){
	local temp
	temp="$(mktemp -d)"
	trap 'rm -rf "$temp"' RETURN
	"$DIR/tadpole.sh" -Xmx400m in="$DIR/testdata/crossk_long_bridge_x.fa" out="$temp/explicit.fa" \
		mode=extend k=95 er=10 el=10 ecc=t hashkmers=f prealloc=f mcs=1 mce=1 mincr=2 \
		prefilter=f pop=f ow=t t=1 showstats=f 1>"$temp/explicit.stdout" 2>"$temp/explicit.stderr"
	"$DIR/tadpole.sh" -Xmx400m in="$DIR/testdata/crossk_long_bridge_x.fa" out="$temp/pair.fa" \
		mode=extend k=95 er=10 el=10 ecc=t hashkmers=pair prealloc=f mcs=1 mce=1 mincr=2 \
		prefilter=f pop=f ow=t t=1 showstats=f 1>"$temp/pair.stdout" 2>"$temp/pair.stderr"
	"$DIR/tadpole.sh" -Xmx400m in="$DIR/testdata/crossk_long_bridge_x.fa" out="$temp/fixed.fa" \
		mode=extend k=95 er=10 el=10 ecc=t hashkmers=fixed prealloc=t mcs=1 mce=1 mincr=2 \
		prefilter=f pop=f ow=t t=1 showstats=f 1>"$temp/fixed.stdout" 2>"$temp/fixed.stderr"
	if ! cmp -s "$temp/explicit.fa" "$temp/pair.fa" || ! cmp -s "$temp/explicit.fa" "$temp/fixed.fa"; then
		echo "FAIL: explicit, hash-pair, and fixed-fingerprint extension/ECC outputs differed" >&2
		diff -u "$temp/explicit.fa" "$temp/pair.fa" >&2 || true
		diff -u "$temp/explicit.fa" "$temp/fixed.fa" >&2 || true
		exit 1
	fi
	echo "PASS: hashedExtensionAndCorrectionMatchExplicitK95"
	"$DIR/tadpole.sh" -Xmx400m in="$DIR/testdata/crossk_long_bridge_x.fa" out="$temp/k64_explicit.fa" \
		mode=correct k=64 hashkmers=f prealloc=f mcs=1 mce=1 prefilter=f pop=f ow=t t=1 showstats=f \
		1>"$temp/k64_explicit.stdout" 2>"$temp/k64_explicit.stderr"
	"$DIR/tadpole.sh" -Xmx400m in="$DIR/testdata/crossk_long_bridge_x.fa" out="$temp/k64_fixed.fa" \
		mode=correct k=64 hashkmers=auto prealloc=t mcs=1 mce=1 prefilter=f pop=f ow=t t=1 showstats=f \
		1>"$temp/k64_fixed.stdout" 2>"$temp/k64_fixed.stderr"
	if ! cmp -s "$temp/k64_explicit.fa" "$temp/k64_fixed.fa" \
			|| ! grep -Eq '^hashkmers[[:space:]]+fixed$' "$temp/k64_fixed.stderr"; then
		echo "FAIL: automatic preallocated k=64 fingerprint storage differed or was not selected" >&2
		exit 1
	fi
	echo "PASS: preallocatedK64UsesMatchingFixedFingerprint"
}

runHashedGraphParityTest(){
	local temp mode opts
	temp="$(mktemp -d)"
	trap 'rm -rf "$temp"' RETURN
	for mode in explicit pair fixed; do
		if [ "$mode" = explicit ]; then opts="hashkmers=f prealloc=f"
		elif [ "$mode" = pair ]; then opts="hashkmers=pair prealloc=f"
		else opts="hashkmers=fixed prealloc=t"
		fi
		"$DIR/tadpole.sh" -Xmx400m in="$DIR/testdata/crossk_long_bridge_x.fa" out="$temp/$mode.fa" \
			gfa="$temp/$mode.gfa" assemblek=31 fusek=none bridgek=none graphk=95 $opts \
			mcs=1 mce=1 mincr=2 mincontig=1 sweeplen=0 prefilter=f pop=f ow=t t=1 showstats=f \
			1>"$temp/$mode.stdout" 2>"$temp/$mode.stderr"
	done
	if ! cmp -s "$temp/explicit.fa" "$temp/pair.fa" || ! cmp -s "$temp/explicit.fa" "$temp/fixed.fa" \
			|| ! cmp -s "$temp/explicit.gfa" "$temp/pair.gfa" || ! cmp -s "$temp/explicit.gfa" "$temp/fixed.gfa"; then
		echo "FAIL: explicit, hash-pair, and fixed-fingerprint graph outputs differed" >&2
		exit 1
	fi
	echo "PASS: hashedPrunedGraphExplorationMatchesExplicitK95"
}

runAssemblyHashGuardTest(){
	local temp
	temp="$(mktemp -d)"
	trap 'rm -rf "$temp"' RETURN
	"$DIR/tadpole.sh" -Xmx400m in="$DIR/testdata/crossk_long_bridge_x.fa" out="$temp/explicit.fa" \
		k=95 hashkmers=f prealloc=t mcs=1 mce=1 mincontig=1 sweeplen=0 buildthreads=1 \
		prefilter=f pop=f ow=t t=1 showstats=f 1>"$temp/explicit.stdout" 2>"$temp/explicit.stderr"
	"$DIR/tadpole.sh" -Xmx400m in="$DIR/testdata/crossk_long_bridge_x.fa" out="$temp/auto.fa" \
		k=95 hashkmers=auto prealloc=t mcs=1 mce=1 mincontig=1 sweeplen=0 buildthreads=1 \
		prefilter=f pop=f ow=t t=1 showstats=f 1>"$temp/auto.stdout" 2>"$temp/auto.stderr"
	if ! cmp -s "$temp/explicit.fa" "$temp/auto.fa" || grep -Eq '^hashkmers[[:space:]]' "$temp/auto.stderr"; then
		echo "FAIL: initial contig assembly did not retain explicit kmer storage" >&2
		exit 1
	fi
	if "$DIR/tadpole.sh" -Xmx400m in="$DIR/testdata/crossk_long_bridge_x.fa" out="$temp/invalid.fa" \
			k=95 hashkmers=pair mcs=1 mce=1 mincontig=1 sweeplen=0 buildthreads=1 prefilter=f pop=f ow=t t=1 showstats=f \
			1>"$temp/invalid.stdout" 2>"$temp/invalid.stderr"; then
		echo "FAIL: explicit compact storage was accepted for table-seeded contig construction" >&2
		exit 1
	fi
	if ! grep -Fq 'requires explicit kmer enumeration' "$temp/invalid.stderr"; then
		echo "FAIL: incompatible compact assembly did not explain the enumeration requirement" >&2
		exit 1
	fi
	echo "PASS: initialAssemblyRetainsExplicitKmers"
}

runMultiKGfaTest(){
	local temp headers writes
	temp="$(mktemp -d)"
	trap 'rm -rf "$temp"' RETURN
	"$DIR/tadpole.sh" in="$DIR/testdata/crossk_long_bridge_x.fa" out="$temp/out.fa.gz" \
		gfa="$temp/out.gfa" assemblek=31 fusek=none bridgek=60 graphk=60 \
		mcs=1 mce=1 mincontig=1 sweeplen=0 prefilter=f pop=f ow=t t=1 showstats=f \
		1>"$temp/stdout" 2>"$temp/stderr"
	headers="$(grep -c '^H' "$temp/out.gfa")"
	writes="$(grep -c 'Writing GFA contig graph' "$temp/stderr")"
	if [ "$headers" != 1 ] || [ "$writes" != 1 ]; then
		echo "FAIL: multi-K GFA was emitted $writes times with $headers headers" >&2
		cat "$temp/stderr" >&2
		exit 1
	fi
	echo "PASS: multiKGfaWritesFinalGraphOnce"
}

runExecutionPlanTest(){
	local temp actual expected
	temp="$(mktemp -d)"
	trap 'rm -rf "$temp"' RETURN
	"$DIR/tadpole.sh" in="$DIR/testdata/crossk_long_bridge_x.fa" out="$temp/single.fa.gz" \
		k=31 mcs=1 mce=1 mincontig=1 sweeplen=0 prefilter=f pop=f ow=t t=1 showstats=f \
		1>"$temp/single.stdout" 2>"$temp/single.stderr"
	actual="$(awk '/^mode +/{p=1} p{if($0==""){exit} print}' "$temp/single.stderr")"
	expected="$(printf 'mode       assemble\nassemblek  31')"
	if [ "$actual" != "$expected" ]; then
		echo "FAIL: default startup plan was not minimal" >&2
		printf 'Expected:\n%s\nObserved:\n%s\n' "$expected" "$actual" >&2
		exit 1
	fi
	"$DIR/tadpole.sh" in="$DIR/testdata/crossk_long_bridge_x.fa" out="$temp/multi.fa.gz" \
		k=96,64,32 wash=t simpleomnitigs=t graphk=96 mcs=1 mce=1 mincontig=1 sweeplen=0 \
		prefilter=f pop=f ow=t t=1 showstats=f 1>"$temp/multi.stdout" 2>"$temp/multi.stderr"
	actual="$(awk '/^mode +/{p=1} p{if($0==""){exit} print}' "$temp/multi.stderr")"
	expected="$(printf 'mode       assemble\nextra      shave rinse simpleomnitigs\nassemblek  96\nfusek      64,32\nbridgek    64,32\nhashkmers  pair\ngraphk     96')"
	if [ "$actual" != "$expected" ]; then
		echo "FAIL: multi-K startup plan did not report the resolved phases" >&2
		printf 'Expected:\n%s\nObserved:\n%s\n' "$expected" "$actual" >&2
		exit 1
	fi
	echo "PASS: startupPlanReportsOnlyRelevantPhases"
}

runUnifiedMultiGraphTest(){
	local temp expected rc actual
	temp="$(mktemp -d)"
	trap 'rm -rf "$temp"' RETURN
	expected="CATACTCGCCTCGGCATATGTAGCCCGGAATCATCACACTCCTCAGGGAGCTCCTACTAGCTTGGAGCTCATCGCCGTTCACTCGGACTAACTCTGCAGCAAATGACGCCCGACCACGGTCCCAACTGACATCAGCCGCTTTGATGCGAAGCCAATTTCATTCCTCCTCTATCCTAAGCATCGATCTCTCCAACATCTCGACTCTTGGAG"
	rc="$(printf '%s' "$expected" | tr ACGT TGCA | rev)"
	for spec in "default::Loading graph-k table at k=31" \
		"intermediate:graphk=25:Loading graph-k table at k=25" \
		"washed:wash=t:Removing dead ends and error bubbles"; do
		local name="${spec%%:*}"
		local rest="${spec#*:}"
		local grapharg="${rest%%:*}"
		local logtext="${rest#*:}"
		"$DIR/tadpole.sh" in="$DIR/testdata/crossk_left_bridge_reads_k31.fa" out="$temp/$name.fa.gz" \
			k=20,31 $grapharg simpleomnitigs=t mcs=1 mce=1 mincontig=1 sweeplen=0 prefilter=f pop=f \
			ow=t t=1 showstats=f verbose=t 1>"$temp/$name.stdout" 2>"$temp/$name.stderr"
		actual="$(zcat "$temp/$name.fa.gz" | awk '!/^>/{printf "%s", $0} END{print ""}')"
		if [ "$actual" != "$expected" ] && [ "$actual" != "$rc" ]; then
			echo "FAIL: unified multi-k $name graph extraction changed the expected sequence" >&2
			cat "$temp/$name.stderr" >&2
			exit 1
		fi
		if ! grep -Fq "$logtext" "$temp/$name.stderr"; then
			echo "FAIL: unified multi-k $name did not use the expected graph-k lifecycle" >&2
			cat "$temp/$name.stderr" >&2
			exit 1
		fi
		if ! grep -Fq "Graph-k endpoints refreshed:" "$temp/$name.stderr"; then
			echo "FAIL: unified multi-k $name did not refresh final graph-k endpoints" >&2
			cat "$temp/$name.stderr" >&2
			exit 1
		fi
		for diagnostic in "Graph-k endpoint topology:" "Graph-k unique topology:" \
			"Graph-k ambiguous topology:" "Graph-k missing topology:" \
			"Graph-k endpoint seeds:" "Graph-k traversal exits:"; do
			if ! grep -Fq "$diagnostic" "$temp/$name.stderr"; then
				echo "FAIL: unified multi-k $name did not report '$diagnostic'" >&2
				cat "$temp/$name.stderr" >&2
				exit 1
			fi
		done
		if [ "$name" = "intermediate" ]; then
			for diagnostic in "Graph-k tip overlaps:" "Graph-k overlaps"; do
				if ! grep -Fq "$diagnostic" "$temp/$name.stderr"; then
					echo "FAIL: unified multi-k $name did not report '$diagnostic'" >&2
					cat "$temp/$name.stderr" >&2
					exit 1
				fi
			done
		fi
		echo "PASS: unifiedMultiKGraph${name^}"
	done
}

runQuietMultiGraphTest(){
	local temp diagnostic
	temp="$(mktemp -d)"
	trap 'rm -rf "$temp"' RETURN
	"$DIR/tadpole.sh" in="$DIR/testdata/crossk_left_bridge_reads_k31.fa" out="$temp/out.fa.gz" \
		k=20,31 graphk=20 simpleomnitigs=t resolverepeats=t mcs=1 mce=1 mincontig=1 sweeplen=0 \
		prefilter=f pop=f ow=t t=1 showstats=f 1>"$temp/stdout" 2>"$temp/stderr"
	if ! grep -Fq "X repeats resolved:" "$temp/stderr"; then
		echo "FAIL: unified multi-k repeat resolution did not run on the longest-k assembly" >&2
		cat "$temp/stderr" >&2
		exit 1
	fi
	for diagnostic in "Cross-k merge:" "Cross-k tip overlaps:" "Graph-k tip overlaps:" \
		"Cross-k endpoints:" "Graph-k endpoints refreshed:" "Graph-k endpoint topology:" \
		"Graph-k unique topology:" "Graph-k ambiguous topology:" "Graph-k missing topology:" \
		"Graph-k endpoint seeds:" "Graph-k traversal exits:" "Cross-k merge evaluations:" \
		"Popping bubbles; contigs="; do
		if grep -Fq "$diagnostic" "$temp/stderr"; then
			echo "FAIL: default multi-k output contained internal diagnostic '$diagnostic'" >&2
			cat "$temp/stderr" >&2
			exit 1
		fi
	done
	echo "PASS: quietMultiKRepeatResolution"
}

resolveSymlinks
runTest assemble.BubblePopperUnitTest
runTest assemble.PathPreservingBubbleSimplifierSpec
runTest assemble.ReadThreadedXResolverUnitTest
runTest assemble.CrossKTipOverlapperUnitTest
runTest assemble.SimpleOmnitigExtractorUnitTest
runTest assemble.TadpoleMultiUnitTest
runTest assemble.ContigGraphClassifierUnitTest
runCrossKLeftBridgeTest
runLongKBridgeEligibilityTest
runLongKBridgeResolutionTest
runLowDepthFateTrackingTest
runShortContigRetentionTest
runHashedBridgeParityTest
runHashedSequenceModeParityTest
runHashedGraphParityTest
runAssemblyHashGuardTest
runMultiKGfaTest
runExecutionPlanTest
runUnifiedMultiGraphTest
runQuietMultiGraphTest
