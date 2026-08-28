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
			k="$kmers" mcs=1 mce=1 mincontig=1 prefilter=f pop=f ow=t t=1 showstats=f \
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

runUnifiedMultiGraphTest(){
	local temp expected rc actual
	temp="$(mktemp -d)"
	trap 'rm -rf "$temp"' RETURN
	expected="CATACTCGCCTCGGCATATGTAGCCCGGAATCATCACACTCCTCAGGGAGCTCCTACTAGCTTGGAGCTCATCGCCGTTCACTCGGACTAACTCTGCAGCAAATGACGCCCGACCACGGTCCCAACTGACATCAGCCGCTTTGATGCGAAGCCAATTTCATTCCTCCTCTATCCTAAGCATCGATCTCTCCAACATCTCGACTCTTGGAG"
	rc="$(printf '%s' "$expected" | tr ACGT TGCA | rev)"
	for spec in "default::Reusing shortest-k table for final graph at k=20" \
		"intermediate:graphk=25:Loading graph-k table at k=25" \
		"washed:wash=t:Removing dead ends and error bubbles"; do
		local name="${spec%%:*}"
		local rest="${spec#*:}"
		local grapharg="${rest%%:*}"
		local logtext="${rest#*:}"
		"$DIR/tadpole.sh" in="$DIR/testdata/crossk_left_bridge_reads_k31.fa" out="$temp/$name.fa.gz" \
			k=20,31 $grapharg simpleomnitigs=t mcs=1 mce=1 mincontig=1 prefilter=f pop=f \
			ow=t t=1 showstats=f 1>"$temp/$name.stdout" 2>"$temp/$name.stderr"
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
			"Graph-k endpoint seeds:" "Graph-k traversal exits:" "Graph-k tip overlaps:" \
			"Graph-k overlaps"; do
			if ! grep -Fq "$diagnostic" "$temp/$name.stderr"; then
				echo "FAIL: unified multi-k $name did not report '$diagnostic'" >&2
				cat "$temp/$name.stderr" >&2
				exit 1
			fi
		done
		echo "PASS: unifiedMultiKGraph${name^}"
	done
}

resolveSymlinks
runTest assemble.BubblePopperUnitTest
runTest assemble.PathPreservingBubbleSimplifierSpec
runTest assemble.ReadThreadedXResolverUnitTest
runTest assemble.CrossKTipOverlapperUnitTest
runTest assemble.SimpleOmnitigExtractorUnitTest
runTest assemble.TadpoleMultiUnitTest
runCrossKLeftBridgeTest
runUnifiedMultiGraphTest
