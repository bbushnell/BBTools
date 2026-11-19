package test.benchmark;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import benchmark.CallGenesPartialScorecard;
import benchmark.CallGenesPartialScorecard.ContigLengths;
import benchmark.CallGenesPartialScorecard.References;
import benchmark.CallGenesPartialScorecard.Result;
import benchmark.CallGenesBenchmarkSupport;
import gff.GffLine;
import shared.Shared;

/** Literal partial/collision oracle and path-only/legacy saved-input CLI cases. */
public final class CallGenesPartialScorecardTest{
	/** Checks literal partial/edge grades and saved-input CLI contracts under one fresh fixture directory. */
	public static void main(String[] args)throws Exception{
		check(args.length==1, "One fresh fixture directory is required");
		final Path root=Paths.get(args[0]);
		Files.createDirectory(root);
		Shared.setThreads(1);
		CallGenesBenchmarkSupport.cdsOnly();
		final String refs=g(10, 99, 0, "partial=10")+g(13, 99, 0, "partial=00")
			+g(200, 299, 0, "partial=false")+g(203, 299, 0, "partial=true")
			+g(400, 489, 1, "partial=10")+g(500, 589, 1, "partial=01")+g(600, 689, 0, "partial=11")
			+g(700, 789, 0, "partial=true;start_range=.,700")+g(800, 889, 0, "partial=false;start_range=.,800")
			+g(900, 989, 0, "partial=00")+g(1000, 1089, 0, "partial=true;end_range=1089,.")
			+g(1200, 1289, 1, "partial=true")+g(1300, 1389, 0, "ID=complete");
		final String calls=g(10, 99, 0, ".")+g(10, 99, 0, ".")+g(200, 299, 0, ".")+g(400, 480, 1, ".")
			+g(500, 589, 1, ".")+g(500, 589, 1, ".")+g(600, 689, 0, ".")+g(700, 789, 0, ".")+g(800, 889, 0, ".")
			+g(1089, 1100, 1, ".")+g(1500, 1589, 0, ".")+g(1000, 1089, 0, ".").replace("\tCDS\t", "\tmisc_feature\t");
		final References reference=new References(lines(refs));
		final Result result=analyze(reference, lines(calls));
		check(Arrays.equals(reference.raw, new long[]{13, 5, 3, 2, 1, 2}), "Raw partial categories include both shadowed rows");
		check(Arrays.equals(reference.effective, new long[]{11, 4, 2, 2, 1, 2}), "Effective census preserves original winners before partition");
		check(reference.all.duplicateKeys==2 && reference.mixedPartialKeys==2, "Mixed flag collisions are explicit");
		check(reference.complete.uniqueStops()==4 && reference.excluded.uniqueStops()==7, "Unknown-sided true remains excluded");
		check(Arrays.equals(result.all, new long[]{5, 9, 7, 4, 5, 9, 7, 3}), "Literal original eight-field grade oracle");
		check(Arrays.equals(result.complete, new long[]{1, 3, 3, 2, 1, 3, 5, 3}), "Wrong-start partial matches are dropped, wrong strand/type and unmatched calls remain extras");
		check(Arrays.equals(result.dropped, new long[]{3, 1, 1, 1}), "Six dropped occurrences retain duplicate multiplicity and biological sides");
		check(result.allKeys-result.all[3]==7 && result.all[1]==9, "Recall counts seven unique stops while precision counts nine matching occurrences");
		check(result.completeKeys-result.complete[3]==2 && result.complete[1]==3, "Filtered recall and precision retain distinct denominators");
		check(CallGenesPartialScorecard.ratio(0, 0).equals("NA"), "Empty denominators are unavailable");
		check(CallGenesPartialScorecard.ratio(2, 4).equals("0.5"), "Literal rate");
		final References terminalRefs=new References(lines(refs), true);
		final Result terminal=analyze(terminalRefs, lines(calls));
		check(terminalRefs.complete.uniqueStops()==8 && terminalRefs.excluded.uniqueStops()==3,
			"Only three-prime and both-end winners are excluded; unknown and five-only remain eligible");
		check(Arrays.equals(terminal.all, result.all), "Selecting terminal-stop eligibility preserves all-gene grades");
		check(Arrays.equals(terminal.complete, new long[]{4, 7, 5, 3, 4, 7, 6, 3}),
			"Literal terminal-stop grade retains unknown/five-only wrong starts and duplicate calls");
		check(Arrays.equals(terminal.dropped, new long[]{0, 1, 1, 0}), "Only the three-prime and both-end matching calls are dropped");
		collisionTests();
		sideTests();
		normalizationTests();
		final Result empty=analyze(new References(Collections.<GffLine>emptyList()), lines(g(10, 99, 0, ".")));
		check(empty.completeKeys==0 && empty.complete[7]==1, "Empty reference retains unmatched query extras");
		final Result noCalls=analyze(reference, Collections.<GffLine>emptyList());
		check(noCalls.all[3]==11 && noCalls.complete[3]==4 && noCalls.complete[7]==0, "Empty calls miss every retained key");
		write(root.resolve("reference.fa"), ">c description with spaces\n"+repeat('A', 1600)+"\n");
		write(root.resolve("reference.gff"), refs);
		write(root.resolve("calls.gff"), calls);
		final String gradeHeader="genome\tarm\t"+String.join("\t", CallGenesBenchmarkSupport.GRADE_NAMES)+"\n";
		write(root.resolve("grades.tsv"), gradeHeader+"g0001\tbonus400\t5\t9\t7\t4\t5\t9\t7\t3\n");
		final Path manifest=root.resolve("inputs.tsv"), expected=root.resolve("expected.tsv");
		manifest(root, manifest);
		write(expected, CallGenesPartialScorecard.EXPECTED_HEADER+"\nfixture\t1\t11\t9\t4\t3\t7\t14\n");
		new CallGenesPartialScorecard().run(manifest, expected, root.resolve("result"));
		final List<String> summary=Files.readAllLines(root.resolve("result/scorecard.tsv"), StandardCharsets.US_ASCII);
		check(summary.size()==3, "Both all-gene and partial-excluded scorecards emitted");
		check(summary.get(2).equals("fixture\tALL\texclude_partial\t4\t3\t2\t3\t5\t8\t2\t4\t0.5\t3\t6\t0.5\t1\t3\t2\t3\t5\t8\t2\t4\t0.5\t3\t6\t0.5\t1\t1\t3\t3\t2\t1\t3\t5\t3"), "Literal complete scorecard including rates, duplicate count and all eight grader fields");
		check(Files.readAllLines(root.resolve("result/dropped_calls.tsv"), StandardCharsets.US_ASCII).contains("fixture\tALL\t6\t3\t1\t1\t1\t4\t2"), "Both-end count contributes to inclusive side totals only");
		check(Files.readAllLines(root.resolve("result/normalized_dropped_calls.tsv"), StandardCharsets.US_ASCII).contains("fixture\tALL\t6\t3\t1\t1\t1\t4\t2"),
			"Adjusted dropped-call census is emitted separately");
		check(Files.readAllLines(root.resolve("result/start_only_partial_starts.tsv"), StandardCharsets.US_ASCII).get(0)
			.startsWith("policy\tgenome\tfive_prime_references"), "Start-only adjusted-start cohort is emitted");
		check(Files.readAllLines(root.resolve("result/partial_counts.tsv"), StandardCharsets.US_ASCII).contains("ALL\teffective_keys\t11\t4\t7\t2\t2\t1\t2\t3\t3\t2\t2"), "Effective annotation census retains exclusive and inclusive sides");
		CallGenesPartialScorecard.main(new String[]{"in="+manifest, "expected="+expected,
			"out="+root.resolve("terminal"), "exclude=threeprime"});
		final List<String> terminalSummary=Files.readAllLines(root.resolve("terminal/scorecard.tsv"), StandardCharsets.US_ASCII);
		check(terminalSummary.size()==3 && terminalSummary.get(1).equals(summary.get(1)), "CLI mode retains the exact all-gene row");
		check(terminalSummary.get(2).equals("fixture\tALL\tterminal_stop_eligible_by_annotation\t8\t7\t3\t3\t6\t11\t5\t8\t0.625\t7\t10\t0.7\t2\t7\t3\t3\t6\t11\t5\t8\t0.625\t7\t10\t0.7\t2\t4\t7\t5\t3\t4\t7\t6\t3"),
			"Terminal-stop CLI emits the literal counts and explicit occurrence/key rates");
		check(Files.readAllLines(root.resolve("terminal/dropped_calls.tsv"), StandardCharsets.US_ASCII).contains("fixture\tALL\t2\t0\t1\t1\t0\t1\t2"),
			"Unknown-side and five-prime-only calls are explicitly retained in the dropped-call census");
		check(Arrays.equals(Files.readAllBytes(root.resolve("result/partial_counts.tsv")), Files.readAllBytes(root.resolve("terminal/partial_counts.tsv"))),
			"The annotation census, including unknown cases, is independent of exclusion mode");
		boolean invalidMode=false;
		try{CallGenesPartialScorecard.main(new String[]{"exclude=unknown"});}catch(IllegalArgumentException e){invalidMode=true;}
		check(invalidMode, "Unknown exclusion policy fails rather than changing the reference population silently");
		final Path legacy=root.resolve("legacy_inputs.tsv");
		final StringBuilder old=new StringBuilder(CallGenesPartialScorecard.LEGACY_INPUT_HEADER+"\nfixture\tg0001");
		for(String file:new String[]{"reference.fa", "reference.gff", "calls.gff", "grades.tsv"}){
			old.append('\t').append(root.resolve(file)).append("\tignored");
		}
		write(legacy, old.append("\tbonus400\n").toString());
		CallGenesPartialScorecard.main(new String[]{"in="+legacy, "insha80=ignored", "expected="+expected,
			"expectedsha80=ignored", "out="+root.resolve("legacy")});
		check(Arrays.equals(Files.readAllBytes(root.resolve("result/scorecard.tsv")),
			Files.readAllBytes(root.resolve("legacy/scorecard.tsv"))), "Legacy pin fields are ignored without changing science");
		new CallGenesPartialScorecard().run(manifest, "ignored", expected, "ignored", root.resolve("legacy_api"));
		check(Arrays.equals(Files.readAllBytes(root.resolve("result/scorecard.tsv")),
			Files.readAllBytes(root.resolve("legacy_api/scorecard.tsv"))), "Legacy API delegates to path-only scoring");
		final String one=Files.readAllLines(manifest, StandardCharsets.US_ASCII).get(1);
		final String expectedText=new String(Files.readAllBytes(expected), StandardCharsets.US_ASCII);
		final Path copy=root.resolve("same_reference_other_path.gff");
		Files.copy(root.resolve("reference.gff"), copy);
		write(expected, expectedText+"other\t1\t11\t9\t4\t3\t7\t14\n");
		write(manifest, CallGenesPartialScorecard.INPUT_HEADER+"\n"+one+"\n"+
			one.replace("fixture\t", "other\t").replace(root.resolve("reference.gff").toString(), copy.toString())+"\n");
		final String conflict=failure(manifest, expected, root.resolve("path_conflict"));
		check(conflict!=null && conflict.contains("same literal reference path"), "Different literal reference paths fail even with identical bytes");
		write(expected, expectedText);
		write(manifest, CallGenesPartialScorecard.INPUT_HEADER+"\n"+one+"\n"+one+"\n");
		check(fails(manifest, expected, root.resolve("duplicate")), "Duplicate policy/genome rows fail");
		manifest(root, manifest);
		write(expected, CallGenesPartialScorecard.EXPECTED_HEADER+"\nfixture\t2\t22\t18\t8\t6\t14\t28\n");
		check(fails(manifest, expected, root.resolve("missing")), "Missing genome cannot pass expected roster size");
		write(expected, CallGenesPartialScorecard.EXPECTED_HEADER+"\nfixture\t1\t11\t9\t4\t3\t7\t14\n");
		write(root.resolve("grades.tsv"), gradeHeader+"g0001\tbonus400\t6\t9\t7\t4\t5\t9\t7\t3\n");
		manifest(root, manifest);
		check(fails(manifest, expected, root.resolve("badgrade")), "One wrong all-gene grade field fails before filtered publication");
		write(root.resolve("grades.tsv"), gradeHeader+"g0001\tbonus400\t5\t9\t7\t4\t5\t9\t7\t3\n");
		write(root.resolve("reference.fa"), ">c one\nA\n>c two\nA\n");
		manifest(root, manifest);
		final String duplicateFasta=failure(manifest, expected, root.resolve("duplicateFasta"));
		check(duplicateFasta!=null && duplicateFasta.contains("Duplicate FASTA contig first token: c"),
			"Duplicate FASTA first-token contig ids fail with their own diagnostic before scoring");
		System.out.println("CallGenesPartialScorecardTest PASS checks="+checks);
	}
	/*--------------------------------------------------------------*/
	/*----------------         Edge fixtures        ----------------*/
	/*--------------------------------------------------------------*/

	/** Checks original last-wins semantics in both mixed-partial orders and across strands. */
	private static void collisionTests(){
		final String partial=g(10, 99, 0, "partial=10"), complete=g(13, 99, 0, "partial=00"), call=g(13, 99, 0, ".");
		Result r=analyze(new References(lines(complete+partial)), lines(call));
		check(r.completeKeys==0 && r.complete[7]==0 && r.dropped[0]==1, "Partial winner hides earlier complete annotation; no resurrection");
		r=analyze(new References(lines(partial+complete)), lines(call));
		check(r.completeKeys==1 && r.complete[0]==1 && r.dropped[0]==0, "Complete winner keeps call despite shadowed partial");
		final String otherStrand=g(99, 180, 1, "partial=true");
		r=analyze(new References(lines(complete+otherStrand)), lines(call));
		check(r.all[7]==1 && r.complete[7]==1 && r.completeKeys==0, "Wrong-strand last winner remains authoritative even after exclusion");
		final String three=g(10, 99, 0, "partial=01");
		r=analyze(new References(lines(complete+three), true), lines(call));
		check(r.completeKeys==0 && r.complete[7]==0 && r.dropped[1]==1, "A three-prime winner cannot resurrect a shadowed complete reference");
		r=analyze(new References(lines(three+complete), true), lines(call));
		check(r.completeKeys==1 && r.complete[1]==1 && r.dropped[1]==0, "A complete last winner survives a shadowed three-prime row");
		r=analyze(new References(lines(complete+partial), true), lines(call));
		check(r.completeKeys==1 && r.complete[1]==1, "A five-prime-only last winner remains authoritative and eligible");
	}
	/** Checks explicit flags, biological sides and malformed attributes. */
	private static void sideTests(){
		for(int strand=0; strand<2; strand++){
			for(String flag:new String[]{"10", "01", "11", "true", "false", "00"}){
				final References refs=new References(lines(g(10, 99, strand, "partial="+flag)));
				check(refs.excluded.uniqueStops()==(flag.equals("false") || flag.equals("00")?0:1), "Exclusion follows flag on either strand: "+flag);
				final References terminal=new References(lines(g(10, 99, strand, "partial="+flag)), true);
				final boolean remove=flag.equals("11") || flag.equals(strand==0 ? "01" : "10");
				check(terminal.excluded.uniqueStops()==(remove?1:0), "Terminal-stop exclusion uses biological strand and retains unknown flags: "+flag);
			}
		}
		final References reverse=new References(lines(g(10, 99, 1, "partial=true;end_range=99,.")));
		check(reverse.effective[2]==1, "High coordinate range means five-prime on minus strand");
		check(new References(lines(g(10, 99, 1, "partial=true;end_range=99,.")), true).complete.uniqueStops()==1,
			"Minus-strand high-coordinate partiality is five-prime and retained");
		check(new References(lines(g(10, 99, 1, "partial=true;start_range=.,10")), true).excluded.uniqueStops()==1,
			"Minus-strand low-coordinate partiality is three-prime and excluded");
		check(new References(lines(g(10, 99, 0, "partial=10;end_range=99,.")), true).complete.uniqueStops()==1,
			"Contradictory side evidence remains unknown and is retained");
		final References absent=new References(lines(g(10, 99, 0, "start_range=.,10")));
		check(absent.complete.uniqueStops()==1, "Range alone does not override absence of explicit partial flag");
		boolean bad=false;
		try{new References(lines(g(10, 99, 0, "partial=oops")));}
		catch(IllegalArgumentException e){bad=true;}
		check(bad, "Malformed flag fails loudly");
		bad=false;
		try{new References(lines(g(10, 99, 0, "partial=true;partial=false")));}
		catch(IllegalArgumentException e){bad=true;}
		check(bad, "Duplicate flag fails loudly");
	}
	/** Checks Brian's symmetric narrow physical-edge 1-2 base normalization. */
	private static void normalizationTests(){
		final String plusStart=gPhase("c", 1, 92, 0, 2, "partial=10;start_range=.,1");
		final String plusStop=gPhase("c", 10, 100, 0, 0, "partial=01;end_range=100,.");
		final String minusStart=gPhase("d", 10, 101, 1, 2, "partial=01;end_range=101,.");
		final String minusStop=gPhase("c", 1, 91, 1, 0, "partial=10;start_range=.,1");
		final References refs=new References(lines(plusStart+plusStop+minusStart+minusStop));
		final String calls=gPhase("c", 3, 92, 0, 0, ".")+gPhase("c", 10, 99, 0, 0, ".")+gPhase("d", 10, 99, 1, 0, ".")+gPhase("c", 2, 91, 1, 0, ".")+
			gPhase(10, 97, 0, 0, "partial=01");
		final Result result=CallGenesPartialScorecard.analyze(refs, lines(calls), lengths(100, 101));
		check(Arrays.equals(result.all, new long[]{0, 2, 4, 2, 0, 2, 5, 3}), "Exact scorecard still charges 1-2 base edge endpoints");
		check(Arrays.equals(result.allNormalized, new long[]{4, 4, 0, 0, 4, 4, 1, 1}), "Adjusted scorecard credits start and stop leftovers on both strands");
		check(Arrays.equals(result.startOnly, new long[]{2, 2, 0, 2, 2, 2, 0, 2}),
			"Five-prime-only cohort counts matched-stop start repair separately from native all-key start errors");
		final ContigLengths lengths=lengths(100);
		final GffLine ref=lines(plusStart).get(0), near=lines(gPhase(3, 92, 0, 0, ".")).get(0), far=lines(gPhase(4, 92, 0, 0, ".")).get(0);
		check(CallGenesPartialScorecard.sameStartOrPhysicalEdge(ref, near, lengths), "A two-base plus-strand fake start is adjusted");
		check(!CallGenesPartialScorecard.sameStartOrPhysicalEdge(ref, far, lengths), "A three-base plus-strand start gap is outside the narrow adjustment");
		check(!CallGenesPartialScorecard.sameStartOrPhysicalEdge(lines(gPhase(1, 92, 0, 2, "partial=10")).get(0), near, lengths),
			"Physical-edge normalization requires the literal GFF open-boundary range");
		check(!CallGenesPartialScorecard.sameStartOrPhysicalEdge(ref, lines(gPhase(2, 92, 0, 1, ".")).get(0), lengths),
			"Plus-strand starts require the actual first full codon, not a compensating same-frame phase");
		check(!CallGenesPartialScorecard.sameStartOrPhysicalEdge(ref, lines(gPhase(3, 92, 0, 1, ".")).get(0), lengths),
			"Plus-strand starts require frame agreement at the first full codon");
		check(!CallGenesPartialScorecard.sameStartOrPhysicalEdge(lines(gPhase("d", 10, 101, 1, 2, "partial=01;end_range=101,.")).get(0),
			lines(gPhase("d", 10, 100, 1, 1, ".")).get(0), lengths(100, 101)),
			"Minus-strand starts require the actual first full codon, not a compensating same-frame phase");
		check(!CallGenesPartialScorecard.sameStartOrPhysicalEdge(lines(gPhase("d", 10, 101, 1, 2, "partial=01;end_range=101,.")).get(0),
			lines(gPhase("d", 10, 99, 1, 1, ".")).get(0), lengths(100, 101)),
			"Minus-strand starts require frame agreement at the first full codon");
		check(!CallGenesPartialScorecard.sameStartOrPhysicalEdge(lines(gPhase(5, 90, 0, 2, "partial=10;start_range=.,5")).get(0),
			lines(gPhase(7, 90, 0, 0, ".")).get(0), lengths), "Internal start_range gaps are not physical edges");
		check(!CallGenesPartialScorecard.sameStopOrPhysicalEdge(lines(gPhase(1, 101, 0, 0, "partial=01;end_range=101,.")).get(0),
			lines(gPhase(1, 100, 0, 0, ".")).get(0), lengths(101)), "A one-base plus-strand stop leftover with 100 coding bases is not credited");
		check(!CallGenesPartialScorecard.sameStopOrPhysicalEdge(lines(gPhase(3, 101, 0, 0, "partial=01;end_range=101,.")).get(0),
			lines(gPhase(3, 100, 0, 0, ".")).get(0), lengths(101)), "A one-base plus-strand stop leftover with 98 coding bases is not credited");
		check(!CallGenesPartialScorecard.sameStopOrPhysicalEdge(lines(gPhase(1, 92, 1, 0, "partial=10;start_range=.,1")).get(0),
			lines(gPhase(2, 92, 1, 0, ".")).get(0), lengths(92)), "A one-base minus-strand stop leftover with 91 coding bases is not credited");
		check(!CallGenesPartialScorecard.sameStopOrPhysicalEdge(lines(gPhase(1, 93, 1, 0, "partial=10;start_range=.,1")).get(0),
			lines(gPhase(2, 93, 1, 0, ".")).get(0), lengths(93)), "A one-base minus-strand stop leftover with 92 coding bases is not credited");
		check(!CallGenesPartialScorecard.sameStopOrPhysicalEdge(lines(gPhase(9, 99, 0, 0, "partial=01;end_range=99,.")).get(0),
			lines(gPhase(9, 98, 0, 0, ".")).get(0), lengths), "Right open ranges inside a contig are gap-associated, not physical edges");
		check(!CallGenesPartialScorecard.sameStopOrPhysicalEdge(lines(gPhase(1, 101, 0, 0, "partial=01;end_range=101,.")).get(0),
			lines(gPhase(2, 99, 0, 0, ".")).get(0), lengths(101)), "Plus-strand stops require frame agreement after endpoint normalization");
		check(!CallGenesPartialScorecard.sameStopOrPhysicalEdge(lines(gPhase(1, 101, 1, 0, "partial=10;start_range=.,1")).get(0),
			lines(gPhase(3, 100, 1, 0, ".")).get(0), lengths(101)), "Minus-strand stops require frame agreement after endpoint normalization");
		final String included=gPhase(11, 100, 0, 0, "partial=00"), excluded=gPhase(11, 101, 0, 0, "partial=01;end_range=101,.");
		final Result exactIncluded=CallGenesPartialScorecard.analyze(new References(lines(included+excluded)),
			lines(gPhase(11, 100, 0, 0, ".")), lengths(101));
		check(Arrays.equals(exactIncluded.completeNormalized, new long[]{1, 1, 0, 0, 1, 1, 0, 0}) && exactIncluded.droppedNormalized[1]==0,
			"Global exact retained matches take precedence over nearby excluded physical-edge matches");
		final GffLine unknownExact=lines(gPhase(1, 90, ".", 2, "partial=10;start_range=.,1")).get(0);
		final GffLine qmarkExact=lines(gPhase(1, 90, "?", 2, "partial=10;start_range=.,1")).get(0);
		check(unknownExact.strand==GffLine.DOT && qmarkExact.strand==GffLine.QMARK,
			"Fixture must parse both unknown strand spellings");
		check(CallGenesPartialScorecard.sameStartOrPhysicalEdge(unknownExact, lines(gPhase(1, 90, ".", 0, ".")).get(0), lengths),
			"Unknown-strand exact start equality remains an exact match");
		check(CallGenesPartialScorecard.sameStopOrPhysicalEdge(qmarkExact, lines(gPhase(1, 90, "?", 0, ".")).get(0), lengths),
			"Unknown-strand exact stop equality remains an exact match");
		final ContigLengths longLength=lengths(101);
		final GffLine minusStartRef=lines(gPhase(1, 101, 1, 2, "partial=01;end_range=101,.")).get(0);
		final GffLine dotStart=lines(gPhase(1, 101, ".", 2, "partial=01;end_range=101,.")).get(0);
		final GffLine startQuery=lines(gPhase(1, 99, 1, 0, ".")).get(0);
		final GffLine dotStartQuery=lines(gPhase(1, 99, ".", 0, ".")).get(0);
		check(CallGenesPartialScorecard.sameStartOrPhysicalEdge(minusStartRef, startQuery, longLength),
			"Minus-strand physical-edge start fixture is otherwise valid");
		check(!CallGenesPartialScorecard.sameStartOrPhysicalEdge(dotStart, dotStartQuery, longLength),
			"Start physical-edge normalization requires an explicit plus or minus strand");
		final GffLine minusStopRef=lines(gPhase(1, 101, 1, 0, "partial=10;start_range=.,1")).get(0);
		final GffLine qmarkStop=lines(gPhase(1, 101, "?", 0, "partial=10;start_range=.,1")).get(0);
		final GffLine stopQuery=lines(gPhase(3, 101, 1, 0, ".")).get(0);
		final GffLine qmarkStopQuery=lines(gPhase(3, 101, "?", 0, ".")).get(0);
		check(CallGenesPartialScorecard.sameStopOrPhysicalEdge(minusStopRef, stopQuery, longLength),
			"Minus-strand physical-edge stop fixture is otherwise valid");
		check(!CallGenesPartialScorecard.sameStopOrPhysicalEdge(qmarkStop, qmarkStopQuery, longLength),
			"Stop physical-edge normalization requires an explicit plus or minus strand");
	}
	/*--------------------------------------------------------------*/
	/*----------------        Fixture helpers       ----------------*/
	/*--------------------------------------------------------------*/

	/** Requires rejection without a success marker for a malformed input fixture. */
	private static boolean fails(Path in, Path expected, Path out)throws Exception{
		return failure(in, expected, out)!=null;
	}
	/** Returns the failure diagnostic only when a malformed input fixture leaves no success marker. */
	private static String failure(Path in, Path expected, Path out)throws Exception{
		try{new CallGenesPartialScorecard().run(in, expected, out);}
		catch(IllegalArgumentException e){return Files.exists(out.resolve("PASS")) ? null : e.getMessage();}
		return null;
	}
	/** Grades a literal reference/call pair against the fixture contig. */
	private static Result analyze(References refs, List<GffLine> calls){
		return CallGenesPartialScorecard.analyze(refs, calls, lengths(100));
	}
	/** Creates a one-contig FASTA length map for physical edge tests. */
	private static ContigLengths lengths(int c){
		final java.util.HashMap<String, Integer> map=new java.util.HashMap<String, Integer>();
		map.put("c", c);
		return new ContigLengths(map);
	}
	/** Creates a two-contig FASTA length map for asymmetric physical edge tests. */
	private static ContigLengths lengths(int c, int d){
		final java.util.HashMap<String, Integer> map=new java.util.HashMap<String, Integer>();
		map.put("c", c);
		map.put("d", d);
		return new ContigLengths(map);
	}
	/** Binds fixture FASTA, reference, calls and grade literal paths after each edit. */
	private static void manifest(Path root, Path out)throws Exception{
		final StringBuilder b=new StringBuilder(CallGenesPartialScorecard.INPUT_HEADER+"\nfixture\tg0001");
		for(String file:new String[]{"reference.fa", "reference.gff", "calls.gff", "grades.tsv"}){
			final Path p=root.resolve(file);
			b.append('\t').append(p);
		}
		write(out, b.append("\tbonus400\n").toString());
	}
	/** Parses literal GFF rows through the production GFF parser. */
	private static ArrayList<GffLine> lines(String text){
		final ArrayList<GffLine> lines=new ArrayList<GffLine>();
		for(String line:text.split("\n")){
			if(!line.isEmpty()){lines.add(new GffLine(line.getBytes(StandardCharsets.US_ASCII)));}
		}
		return lines;
	}
	/** Creates one one-based CDS fixture row with an explicit strand. */
	private static String g(int start, int end, int strand, String attributes){return gPhase(start, end, strand, 0, attributes);}
	/** Creates a CDS row on contig c with an explicit numeric phase. */
	private static String gPhase(int start, int end, int strand, int phase, String attributes){return "c\tfixture\tCDS\t"+start+"\t"+end+"\t.\t"+(strand==0?"+":"-")+"\t"+phase+"\t"+attributes+"\n";}
	/** Creates a CDS row on the supplied contig with an explicit numeric phase. */
	private static String gPhase(String seqid, int start, int end, int strand, int phase, String attributes){return seqid+"\tfixture\tCDS\t"+start+"\t"+end+"\t.\t"+(strand==0?"+":"-")+"\t"+phase+"\t"+attributes+"\n";}
	/** Preserves literal strand text, including unknown spellings, in a CDS fixture row. */
	private static String gPhase(int start, int end, String strand, int phase, String attributes){return "c\tfixture\tCDS\t"+start+"\t"+end+"\t.\t"+strand+"\t"+phase+"\t"+attributes+"\n";}
	/** Builds a literal ASCII homopolymer FASTA contig. */
	private static String repeat(char c, int times){
		final char[] a=new char[times];
		Arrays.fill(a, c);
		return new String(a);
	}
	/** Writes exact fixture bytes. */
	private static void write(Path file, String text)throws Exception{Files.write(file, text.getBytes(StandardCharsets.US_ASCII));}
	/** Fails on a violated expectation before incrementing the successful-check count. */
	private static void check(boolean ok, String why){
		if(!ok){throw new AssertionError(why);}
		checks++;
	}
	/** Successful checks in this process; main does not reset it. */
	private static int checks;
}
