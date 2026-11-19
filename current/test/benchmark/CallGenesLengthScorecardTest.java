package test.benchmark;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import benchmark.CallGenesLengthScorecard;
import benchmark.CallGenesBenchmarkSupport;
import benchmark.CallGenesStopAudit;
import gff.GffLine;
import prot.DigestSuffix;

/** Literal independent grade and routing oracles. Runs only in Dori validation. */
public final class CallGenesLengthScorecardTest{
	/** Builds literal saved-call fixtures and checks grades, routing, input rejection and named query caps. */
	public static void main(String[] args)throws Exception{
		check(args.length==1 || args.length==2 && args[1].equals("verification=structural"), "Fresh fixture directory and optional verification=structural required");
		structural=args.length==2;
		final Path root=Paths.get(args[0]);
		Files.createDirectory(root);
		boundaries();
		final String[] ref={g(10, 132, 0, "."), g(400, 597, 1, "."), g(700, 1002, 0, "."), g(1100, 1552, 0, "."),
			g(1700, 2302, 1, "."), g(2400, 2492, 0, "."), g(2600, 2749, 0, "partial=10"), g(2800, 2994, 0, ".")};
		// The shadowed start has the same stop but must never become a winner again.
		write(root.resolve("reference.gff"), g(9, 132, 0, "partial=10")+String.join("", ref));
		final String extras=g(3200, 3397, 0, ".")+g(132, 254, 1, ".")+ref[3].replace("\tCDS\t", "\tmisc_feature\t");
		write(root.resolve("native.gff"), ref[0]+ref[0]+g(400, 594, 1, ".")+ref[2]+ref[4]+ref[6]+extras);
		final String assisted=ref[0]+ref[0]+ref[1]+ref[3]+ref[4]+ref[6]+g(2803, 2994, 0, ".")+extras;
		write(root.resolve("cap150.gff"), assisted);
		write(root.resolve("k5.gff"), assisted);
		final String gh="genome\tarm\t"+String.join("\t", CallGenesBenchmarkSupport.GRADE_NAMES)+"\n";
		final String nativeGrade=gh+"g0001\tnative\t5\t6\t4\t3\t5\t6\t4\t3\n";
		final String assistedGrade=gh+"g0001\tbonus400\t6\t7\t3\t2\t6\t7\t4\t3\n";
		write(root.resolve("native.grades.tsv"), nativeGrade);
		write(root.resolve("cap150.grades.tsv"), assistedGrade);
		write(root.resolve("k5.grades.tsv"), assistedGrade);
		final String expected=CallGenesLengthScorecard.EXPECTED_HEADER+"\nnative\t1\t8\t6\t3\t3\t6\t8\ncap150\t1\t8\t7\t2\t3\t5\t7\nk5\t1\t8\t7\t2\t3\t5\t7\n";
		write(root.resolve("expected.tsv"), expected);
		loci(root, ref);
		diagnostics(root, "cap150", ref);
		diagnostics(root, "k5", ref);
		manifests(root);
		run(root, "result", false);
		final List<String> scores=read(root.resolve("result/scorecard.tsv"));
		check(scores.size()==28, "Three arms by nine explicit bands plus header");
		// Each row is literal: TPstart,TPstop,FNstart,FNstop,TPstartQ,TPstopQ,FPstart,FPstop.
		final long[][] nativeBands={{0, 0, 1, 1, 0, 0, 0, 0}, {2, 2, 1, 1, 2, 2, 1, 1}, {0, 1, 1, 0, 0, 1, 2, 1},
			{1, 1, 0, 0, 1, 1, 0, 0}, {0, 0, 1, 1, 0, 0, 0, 0}, {1, 1, 0, 0, 1, 1, 0, 0}, {1, 1, 0, 0, 1, 1, 1, 1}};
		final long[][] assistedBands={{0, 0, 1, 1, 0, 0, 0, 0}, {2, 3, 1, 0, 2, 3, 2, 1}, {1, 1, 0, 0, 1, 1, 1, 1},
			{0, 0, 1, 1, 0, 0, 0, 0}, {1, 1, 0, 0, 1, 1, 0, 0}, {1, 1, 0, 0, 1, 1, 0, 0}, {1, 1, 0, 0, 1, 1, 1, 1}};
		for(String policy:POLICIES){
			final long[][] want=policy.equals("native")?nativeBands:assistedBands;
			for(int b=0; b<7; b++){
				final String[] row=row(scores, policy, BANDS[b]);
				check(Long.parseLong(row[3])==(b==1?2:1), "Frozen reference census: "+policy+"/"+BANDS[b]);
				for(int j=0; j<8; j++){check(Long.parseLong(row[11+j])==want[b][j], "Literal grade field "+j+": "+policy+"/"+BANDS[b]);}
			}
		}
		check(row(scores, "native", "ALL")[10].equals("1"), "Duplicate matching occurrences remain distinct from unique stops");
		check(row(scores, "cap150", "40_99")[4].equals("4"), "Combined short band includes annotated matches rather than query lengths");
		final List<String> routes=read(root.resolve("result/routing.tsv"));
		check(Arrays.equals(numbers(row(routes, "native", "ALL"), 3), new long[]{8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 5, 0, 0}), "Native has no family searches; wrong type at D does not count as a CDS call");
		for(String policy:new String[]{"cap150", "k5"}){
			check(Arrays.equals(numbers(row(routes, policy, "ALL"), 3), new long[]{8, 0, 2, 6, 6, 4, 4, 3, 2, 1, 4, 2, 2, 2, 5, 6, 2, 1}), "Unique routing, both accepted long precedence, two recoveries and one loss: "+policy);
			check(Arrays.equals(numbers(row(routes, policy, "40_64"), 3), new long[]{2, 0, 0, 2, 2, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 2, 1, 0}), "Rejected short query can gain a call without asserting causal acceptance");
			check(row(routes, policy, "65_99")[8].equals("1"), "Annotated65-aa locus really searched the long view too");
			check(row(routes, policy, "unknown")[5].equals("1"), "Unknown annotation is unscanned, not searched and rejected");
		}
		check(read(root.resolve("result/deltas.tsv")).contains("k5\tnative\tALL\t1\t-1\t0\t-1\t-1"), "Literal full-population delta");
		check(Files.exists(root.resolve("result/PASS")), "Accepted output has success marker");
		savedBands(root);
		nativeComparatorOnly(root);
		check(fails(root, "badpin", true), "Malformed input binding fails before output");
		mutateFail(root, "loci.tsv", 1, 4, "9", "shadow_start");
		mutateFail(root, "loci.tsv", 1, 7, "1", "phase");
		mutateFail(root, "cap150.sets.tsv", 1, 10, "long", "primary");
		mutateFail(root, "cap150.sets.tsv", 1, 3, "130", "stop_conversion");
		mutateFail(root, "cap150.annotations.tsv", 1, 17, "0", "call_flag");
		mutateFail(root, "cap150.annotations.tsv", 1, 11, "0", "scan_flag");
		mutateFail(root, "native.grades.tsv", 1, 2, "6", "grade");
		mutateFail(root, "expected.tsv", 1, 3, "7", "aggregate");
		final Path setFile=root.resolve("cap150.sets.tsv");
		final String sets=contents(setFile);
		write(setFile, dropRow(sets, 1));
		manifests(root);
		check(fails(root, "missing_set", false), "Missing applicable search fails");
		write(setFile, sets);
		final Path annotationFile=root.resolve("cap150.annotations.tsv");
		final String annotations=contents(annotationFile);
		write(annotationFile, dropRow(annotations, 1));
		manifests(root);
		check(fails(root, "missing_annotation", false), "Incomplete reference routing fails");
		write(annotationFile, annotations);
		final Path stopFile=root.resolve("cap150.stops.tsv");
		final String stops=contents(stopFile);
		write(stopFile, stops+stops.split("\n")[1]+"\n");
		manifests(root);
		check(fails(root, "duplicate_stop", false), "Duplicate bonus stop fails");
		write(stopFile, stops);
		write(setFile, sets+sets.split("\n")[1]+"\n");
		manifests(root);
		check(fails(root, "duplicate_set", false), "Duplicate family-set row fails");
		write(setFile, sets);
		manifests(root);
		namedCaps(root, ref);
		System.out.println("CallGenesLengthScorecardTest PASS checks="+checks);
	}

	/** The lean facade must reproduce every existing count, not just an objective. */
	private static void savedBands(Path root)throws Exception{
		final java.util.Set<String> cohort=java.util.Collections.singleton("g0001");
		final CallGenesLengthScorecard.SavedBands saved=CallGenesLengthScorecard.readSavedBands(root.resolve("loci.tsv"), cohort);
		final CallGenesStopAudit.ReferenceSet refs=CallGenesStopAudit.loadReference(root.resolve("reference.gff").toString());
		final CallGenesLengthScorecard.ReferenceBands bound=saved.bind("g0001", refs);
		final StringBuilder actual=new StringBuilder(CallGenesLengthScorecard.scoreHeader());
		final ArrayList<CallGenesLengthScorecard.ScoreCounts> nativeScores=new ArrayList<>();
		for(String policy:POLICIES){
			final ArrayList<GffLine> cds=new ArrayList<>();
			final List<GffLine> mixedCalls=GffLine.loadGffFile(root.resolve(policy+".gff").toString(), null, false);
			int excluded=0;
			for(GffLine call:mixedCalls){
				if(call.type.equals("CDS")){cds.add(call);}else{
					check(call.type.equals("misc_feature") && call.start==1100 && call.stop==1552,
						"Only the deliberate wrong-type fixture record is excluded: "+policy);
					excluded++;
				}
			}
			check(excluded==1, "Exactly one deliberate misc_feature per policy: "+policy);
			boolean mixedRejected=false;
			try{bound.grade(mixedCalls);}catch(IllegalArgumentException e){mixedRejected=true;}
			check(mixedRejected, "Lean grading rejects the deliberate mixed-feature input: "+policy);
			final CallGenesLengthScorecard.ScoreCounts counts=bound.grade(cds);
			check(Arrays.equals(counts.grade(), CallGenesStopAudit.grade(refs, cds)), "Adoption fixture checks original eight-field parity: "+policy);
			actual.append(counts.rows(policy, "g0001"));
			if(policy.equals("native")){
				nativeScores.add(counts);
				nativeScores.add(counts);
				final long[] original=counts.grade(), copy=counts.grade();
				copy[0]=-1;
				check(Arrays.equals(counts.grade(), original), "Returned grade arrays cannot mutate stored counts");
				cds.add(line(g(3300, 3398, 0, ".").replace("\tCDS\t", "\trRNA\t")));
				boolean bad=false;
				try{bound.grade(cds);}catch(IllegalArgumentException e){bad=true;}
				check(bad, "Lean grading rejects an unfiltered RNA call list without changing flags");
			}
		}
		final String strictExpected=strictCdsExpected(read(root.resolve("result/by_genome.tsv")));
		if(!actual.toString().equals(strictExpected)){
			System.err.println("STRICT_CDS_EXPECTED_BEGIN\n"+strictExpected+"STRICT_CDS_EXPECTED_END");
			System.err.println("STRICT_CDS_ACTUAL_BEGIN\n"+actual+"STRICT_CDS_ACTUAL_END");
		}
		check(actual.toString().equals(strictExpected), "Lean facade preserves all16counts/everyband/duplicate winners and occurrences for the same strict-CDS input");
		final List<String> doubled=Arrays.asList(CallGenesLengthScorecard.combineScores(nativeScores).rows("native", "ALL").toString().split("\n"));
		final List<String> original=Arrays.asList(strictCdsExpected(read(root.resolve("result/scorecard.tsv"))).split("\n"));
		for(String band:new String[]{"lt40", "40_64", "65_99", "100_149", "150_199", "ge200", "unknown", "ALL", "40_99"}){
			final long[] a=numbers(row(doubled, "native", band), 3), b=numbers(row(original, "native", band), 3);
			for(int i=0; i<16; i++){check(a[i]==2*b[i], "Combined score count "+i+" must conserve band "+band);}
		}
		boolean bad=false;
		try{CallGenesLengthScorecard.readSavedBands(root.resolve("loci.tsv"), java.util.Collections.singleton("g0983"));}catch(IllegalArgumentException e){bad=true;}
		check(bad, "Wrong saved locus cohort fails before calling/grading");
		final ArrayList<GffLine> changed=new ArrayList<>(refs.sourceRecords());
		changed.add(line(g(11, 132, 0, ".")));
		bad=false;
		try{saved.bind("g0001", CallGenesStopAudit.fromLines(changed));}catch(IllegalArgumentException e){bad=true;}
		check(bad, "Changed last-wins reference coordinates fail the existing locus verifier");
	}

	/** The mixed-input oracle deliberately counts one misc_feature as an unmatched
	 * unknown-band call. Strict CDS removes only its FPstart/FPstop contribution;
	 * preserve every other cell, and keep the mixed-input regression unchanged. */
	private static String strictCdsExpected(List<String> mixed){
		check(mixed.size()==28, "Strict oracle requires three policies by nine bands plus header");
		final StringBuilder strict=new StringBuilder(mixed.get(0)).append('\n');
		final String[] bands={"lt40", "40_64", "65_99", "100_149", "150_199", "ge200", "unknown", "ALL", "40_99"};
		int adjusted=0;
		for(int i=1; i<mixed.size(); i++){
			final String[] cells=mixed.get(i).split("\t", -1);
			check(cells.length==19 && cells[0].equals(POLICIES[(i-1)/9]) && cells[2].equals(bands[(i-1)%9]),
				"Strict oracle retains the exact mixed scorecard schema and row order");
			if(cells[2].equals("unknown") || cells[2].equals("ALL")){
				// extras, stop_errors, start_errors, fp_start, fp_stop: one fewer each.
				for(int column:new int[]{6, 7, 8, 17, 18}){
					final long count=Long.parseLong(cells[column]);
					check(count>=1, "Mixed oracle must contain the deliberate false positive");
					cells[column]=Long.toString(count-1);
				}
				adjusted++;
			}
			strict.append(String.join("\t", cells)).append('\n');
		}
		check(adjusted==6, "Only unknown and ALL change for each of three policies");
		return strict.toString();
	}

	/** Exact boundaries include a stop codon; partial/nontriplet/phase/type stay unknown. */
	private static void boundaries(){
		final int[] aa={0, 39, 40, 64, 65, 99, 100, 149, 150, 199, 200}, band={0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5};
		for(int strand=0; strand<2; strand++){
			for(int i=0; i<aa.length; i++){check(CallGenesLengthScorecard.calledBand(line(g(1, 3*(aa[i]+1), strand, "partial=00")))==band[i], "Called boundary "+aa[i]+" on strand "+strand);}
			for(String partial:new String[]{"10", "01", "11", "true"}){check(CallGenesLengthScorecard.calledBand(line(g(1, 198, strand, "partial="+partial)))==6, "Incomplete span cannot establish full length");}
			check(CallGenesLengthScorecard.calledBand(line(g(1, 198, strand, "partial=false")))==2, "Explicit complete legacy flag");
			check(CallGenesLengthScorecard.calledBand(line(g(1, 197, strand, ".")))==6, "Nontriplet span unknown");
			check(CallGenesLengthScorecard.calledBand(line(g(1, 198, strand, ".").replace("\t0\t", "\t1\t")))==6, "Nonzero phase unknown");
			check(CallGenesLengthScorecard.calledBand(line(g(1, 198, strand, ".").replace("\tCDS\t", "\tmisc_feature\t")))==6, "Wrong feature type unknown");
		}
		check(CallGenesLengthScorecard.trueStop(129, 0)==132, "Plus scanner coordinate adds three");
		check(CallGenesLengthScorecard.trueStop(399, 1)==400, "Minus scanner coordinate adds one");
		boolean bad=false;
		try{CallGenesLengthScorecard.trueStop(1, 2);}catch(IllegalArgumentException e){bad=true;}
		check(bad, "Unknown strand fails");
	}
	/** Writes the literal winner roster in the saved 42-column locus schema. */
	private static void loci(Path root, String[] references)throws Exception{
		final String[] header=new String[42];
		Arrays.fill(header, "unused");
		final String[] names={"genome", "contig", "true_stop_1based", "strand", "reference_start_1based", "reference_left_1based", "reference_right_1based", "phase", "annotation_band"};
		System.arraycopy(names, 0, header, 0, names.length);
		for(int i=names.length; i<42; i++){header[i]="unused"+i;}
		final StringBuilder b=new StringBuilder(String.join("\t", header)).append('\n');
		final String[] bands={"40_64", "65_99", "100_149", "150_199", "ge200", "lt40", "unknown", "40_64"};
		for(int i=0; i<references.length; i++){
			final GffLine g=line(references[i]);
			final String[] cells=new String[42];
			Arrays.fill(cells, "NA");
			final String[] known={"g0001", "c", ""+g.trueStop(), ""+g.strand, ""+g.trueStart(), ""+g.start, ""+g.stop, "0", bands[i]};
			System.arraycopy(known, 0, cells, 0, known.length);
			b.append(String.join("\t", cells)).append('\n');
		}
		write(root.resolve("loci.tsv"), b.toString());
	}
	/** Selects the legacy cap150 or cap100 fixture routing policy. */
	private static void diagnostics(Path root, String policy, String[] references)throws Exception{
		final int cap=policy.equals("cap150")?150:100;
		diagnostics(root, policy, references, cap);
	}
	/** Writes explicit stop, family-set and annotation rows for the requested query cap. */
	private static void diagnostics(Path root, String policy, String[] references, int cap)throws Exception{
		final int[] queries={40, 100, cap, Math.min(120, cap), cap, 0, 0, 64};
		final boolean[] shortAccept={true, true, false, true, false, false, false, false}, longAccept={false, false, false, true, true, false, false, false};
		final int[] nativeCall={1, 1, 1, 0, 1, 0, 1, 0}, called={1, 1, 0, 1, 1, 0, 1, 1};
		final StringBuilder stops=new StringBuilder(CallGenesLengthScorecard.STOPS_HEADER).append('\n'),
			sets=new StringBuilder(CallGenesLengthScorecard.SETS_HEADER).append('\n'), annotations=new StringBuilder(CallGenesLengthScorecard.ANNOTATIONS_HEADER).append('\n');
		for(int i=0; i<references.length; i++){
			final GffLine g=line(references[i]);
			final int stopLeft=g.trueStop()-(g.strand==0?3:1), q=queries[i];
			final boolean accept=shortAccept[i] || longAccept[i];
			final String primary=longAccept[i]?"long":shortAccept[i]?"short":"none";
			if(q>0){
				stops.append(tsv("g0001", "c", "bonus400", stopLeft, g.strand, 220, bit(accept), 1, "NA", 0, 0, nativeCall[i], called[i], q, "stop"));
				sets.append(tsv("g0001", "c description", "bonus400", stopLeft, g.strand, "short", Math.min(q, 99), 220, bit(shortAccept[i]), shortAccept[i]?"familyS":"NA", primary, 1, nativeCall[i], called[i], 4, 0, 0));
				if(q>=100){sets.append(tsv("g0001", "c description", "bonus400", stopLeft, g.strand, "long", q, 220, bit(longAccept[i]), longAccept[i]?"familyL":"NA", primary, 1, nativeCall[i], called[i], 4, 0, 0));}
			}
			annotations.append(tsv("g0001", "c", "bonus400", g.trueStop(), g.strand, g.start, g.stop, g.length()/3, i==6?"10":"00", 0.5, q>0?"scanned":"below_minimum", bit(q>0), bit(accept), "NA", 0, 0, nativeCall[i], called[i]));
		}
		// Identical zero-effect replay rows are a separate diagnostic arm, never additional observations.
		for(StringBuilder table:new StringBuilder[]{stops, sets, annotations}){
			final String data=table.substring(table.indexOf("\n")+1);
			table.append(data.replace("\tbonus400\t", "\tzero_effect\t"));
		}
		write(root.resolve(policy+".stops.tsv"), stops.toString());
		write(root.resolve(policy+".sets.tsv"), sets.toString());
		write(root.resolve(policy+".annotations.tsv"), annotations.toString());
	}

	/** Real selected setting names use explicit150 bindings alongside legacy100/150 arms.
	 * Two genomes expose inconsistent per-policy bindings; the literal data includes150-aa queries. */
	private static void namedCaps(Path original, String[] references)throws Exception{
		final Path root=original.resolve("named_caps");
		Files.createDirectory(root);
		final String[] policies={"native", "cap150", "k5", "K5_512_f8", "K5_1024_f4", "K4b2_1024_f3"};
		final int[] caps={0, 150, 100, 150, 150, 150};
		final StringBuilder input=new StringBuilder(CallGenesLengthScorecard.INPUT_HEADER+"\n"), routing=new StringBuilder(CallGenesLengthScorecard.ROUTE_INPUT_HEADER+"\n");
		final StringBuilder expected=new StringBuilder(CallGenesLengthScorecard.EXPECTED_HEADER+"\n");
		for(String policy:policies){expected.append(policy).append(policy.equals("native")?"\t2\t16\t12\t6\t6\t12\t16\n":"\t2\t16\t14\t4\t6\t10\t14\n");}
		write(root.resolve("expected.tsv"), expected.toString());
		final String oldLoci=contents(original.resolve("loci.tsv"));
		write(root.resolve("loci.tsv"), oldLoci+oldLoci.substring(oldLoci.indexOf('\n')+1).replace("g0001\t", "g0002\t"));
		for(String genome:new String[]{"g0001", "g0002"}){
			final Path dir=root.resolve(genome);
			Files.createDirectory(dir);
			for(int i=0; i<policies.length; i++){
				final String policy=policies[i], source=policy.equals("native")?"native":"cap150";
				final Path ref=original.resolve("reference.gff"), calls=original.resolve(source+".gff"), grades=dir.resolve(policy+".grades.tsv");
				write(grades, contents(original.resolve(source+".grades.tsv")).replace("g0001\t", genome+"\t"));
				input.append(tsv(policy, genome, ref, pin(ref), calls, pin(calls), grades, pin(grades), source.equals("native")?"native":"bonus400"));
				if(policy.equals("native")){continue;}
				diagnostics(dir, policy, references, caps[i]);
				final Path stops=dir.resolve(policy+".stops.tsv"), annotations=dir.resolve(policy+".annotations.tsv"), sets=dir.resolve(policy+".sets.tsv");
				for(Path p:Arrays.asList(stops, annotations, sets)){write(p, contents(p).replace("g0001\t", genome+"\t"));}
				routing.append(tsv(policy, genome, stops, pin(stops), annotations, pin(annotations), sets, pin(sets), caps[i]));
			}
		}
		write(root.resolve("inputs.tsv"), input.toString());
		write(root.resolve("routing.tsv"), routing.toString());
		run(root, "result", false);
		final List<String> scores=read(root.resolve("result/scorecard.tsv")), routes=read(root.resolve("result/routing.tsv"));
		check(scores.size()==55 && routes.size()==55, "Six policies by nine bands, two genomes");
		final long[] want={16, 0, 4, 12, 12, 8, 8, 6, 4, 2, 8, 4, 4, 4, 10, 12, 4, 2};
		for(String policy:new String[]{"K5_512_f8", "K5_1024_f4", "K4b2_1024_f3"}){
			check(Arrays.equals(numbers(row(routes, policy, "ALL"), 3), want), "Literal recovered/lost and routing with named150 policy: "+policy);
			final String[] score=row(scores, policy, "ALL");
			check(Arrays.equals(Arrays.copyOfRange(numbers(score, 3), 0, 8), new long[]{16, 14, 4, 6, 10, 14, 12, 2}), "Literal whole-GFF grades preserve duplicate occurrences: "+policy);
		}
		final String[] rows=routing.toString().split("\n");
		int target=-1;
		for(int i=1; i<rows.length; i++){if(rows[i].startsWith("K5_512_f8\tg0002\t")){target=i;}}
		check(target>0, "Second-genome binding exists");
		String[] changed=rows.clone();
		changed[target]=replaceColumn(changed[target], 8, "100");
		capFailure(root, String.join("\n", changed)+"\n", "inconsistent_cap", "Inconsistent query cap for policy");
		for(int bad:new int[]{0, 99, 149, 151}){
			changed=rows.clone();
			changed[target]=replaceColumn(changed[target], 8, Integer.toString(bad));
			capFailure(root, String.join("\n", changed)+"\n", "invalid_cap_"+bad, "Unsupported query cap binding");
		}
		changed=rows.clone();
		for(int i=1; i<changed.length; i++){if(changed[i].startsWith("K5_512_f8\t")){changed[i]=replaceColumn(changed[i], 8, "100");}}
		capFailure(root, String.join("\n", changed)+"\n", "query_exceeds_cap", "Saved query outside bound view");
		changed=rows.clone();
		for(int i=1; i<changed.length; i++){if(changed[i].startsWith("cap150\t")){changed[i]=replaceColumn(changed[i], 8, "100");}}
		capFailure(root, String.join("\n", changed)+"\n", "baseline_cap", "cap150 comparator requires query cap150");
		capFailure(root, dropRow(routing.toString(), target), "missing_binding", "Missing routing binding");
		capFailure(root, routing+rows[target]+"\n", "duplicate_binding", "Duplicate routing input");
		changed=rows.clone();
		changed[target]=replaceColumn(changed[target], 0, "unknown_policy");
		capFailure(root, String.join("\n", changed)+"\n", "unknown_binding", "Unexpected routing identity");
		final List<String> reversed=new ArrayList<String>(Arrays.asList(rows));
		java.util.Collections.reverse(reversed.subList(1, reversed.size()));
		write(root.resolve("routing.tsv"), String.join("\n", reversed)+"\n");
		run(root, "permuted", false);
		for(String file:new String[]{"scorecard.tsv", "by_genome.tsv", "routing.tsv", "routing_by_genome.tsv", "deltas.tsv"}){
			check(Arrays.equals(Files.readAllBytes(root.resolve("result/"+file)), Files.readAllBytes(root.resolve("permuted/"+file))), "Binding row order preserves output bytes: "+file);
		}
	}
	/** Replaces one tab-separated field while retaining empty trailing fields. */
	private static String replaceColumn(String row, int column, String value){
		final String[] cells=row.split("\t", -1);
		cells[column]=value;
		return String.join("\t", cells);
	}
	/** Checks a specific routing rejection and restores the original table after the expected failure. */
	private static void capFailure(Path root, String changed, String output, String reason)throws Exception{
		final Path routing=root.resolve("routing.tsv");
		final String before=contents(routing);
		write(routing, changed);
		boolean rejected=false;
		try{run(root, output, false);}catch(IllegalArgumentException e){
			rejected=true;
			check(e.getMessage().startsWith(reason), "Specific cap-binding rejection: "+output);
		}
		check(rejected && !Files.exists(root.resolve(output)), "Bad cap binding fails before output: "+output);
		write(routing, before);
	}
	/** Rebinds every policy's current fixture inputs and routing paths after a mutation or restoration. */
	private static void manifests(Path root)throws Exception{
		final StringBuilder in=new StringBuilder(CallGenesLengthScorecard.INPUT_HEADER).append('\n'), routes=new StringBuilder(CallGenesLengthScorecard.ROUTE_INPUT_HEADER).append('\n');
		for(String p:POLICIES){
			final Path ref=root.resolve("reference.gff"), call=root.resolve(p+".gff"), grade=root.resolve(p+".grades.tsv");
			in.append(tsv(p, "g0001", ref, pin(ref), call, pin(call), grade, pin(grade), p.equals("native")?"native":"bonus400"));
			if(!p.equals("native")){
				final Path stop=root.resolve(p+".stops.tsv"), annotation=root.resolve(p+".annotations.tsv"), set=root.resolve(p+".sets.tsv");
				routes.append(tsv(p, "g0001", stop, pin(stop), annotation, pin(annotation), set, pin(set), p.equals("cap150")?150:100));
			}
		}
		write(root.resolve("inputs.tsv"), in.toString());
		write(root.resolve("routing.tsv"), routes.toString());
	}
	/** No invented cap150 row when a new RNA mode has a different native path. */
	private static void nativeComparatorOnly(Path root)throws Exception{
		run(root, "explicit_legacy", false, "comparators=native,cap150");
		for(String file:new String[]{"scorecard.tsv", "by_genome.tsv", "routing.tsv", "routing_by_genome.tsv", "deltas.tsv", "PASS"}){
			check(Arrays.equals(Files.readAllBytes(root.resolve("result/"+file)), Files.readAllBytes(root.resolve("explicit_legacy/"+file))),
				"Explicit legacy comparator preserves bytes: "+file);
		}
		final String[] names={"inputs.tsv", "expected.tsv", "routing.tsv"};
		final String[] original=new String[names.length];
		try{
			for(int i=0; i<names.length; i++){
				original[i]=contents(root.resolve(names[i]));
				final StringBuilder filtered=new StringBuilder();
				for(String line:original[i].split("\n")){if(!line.startsWith("cap150\t")){filtered.append(line).append('\n');}}
				write(root.resolve(names[i]), filtered.toString());
			}
			check(fails(root, "missing_default_cap150", false), "Default still requires historical comparator");
			run(root, "native_only", false, "comparators=native");
			for(String file:new String[]{"scorecard.tsv", "by_genome.tsv", "routing.tsv", "routing_by_genome.tsv"}){
				final ArrayList<String> expected=new ArrayList<String>();
				for(String line:read(root.resolve("result/"+file))){if(!line.startsWith("cap150\t")){expected.add(line);}}
				check(read(root.resolve("native_only/"+file)).equals(expected), "Native-only retains exact grade/routing rows: "+file);
			}
			final List<String> deltas=read(root.resolve("native_only/deltas.tsv"));
			check(deltas.size()==19 && deltas.contains("k5\tnative\tALL\t1\t-1\t0\t-1\t-1"), "Only native deltas emitted");
			for(String bad:new String[]{"cap150", "native,unknown", "native,native"}){
				boolean failed=false;
				try{run(root, "invalid_comparators", false, "comparators="+bad);}catch(IllegalArgumentException expected){failed=true;}
				check(failed && !Files.exists(root.resolve("invalid_comparators")), "Invalid comparator selection fails before output");
			}
		}finally{
			for(int i=0; i<names.length; i++){if(original[i]!=null){write(root.resolve(names[i]), original[i]);}}
		}
	}
	/** Runs the scorecard with current fixture identities and optional explicit comparator arguments. */
	private static void run(Path root, String name, boolean badPin, String... options)throws Exception{
		final ArrayList<String> args=new ArrayList<String>();
		for(String key:new String[]{"in", "expected", "loci", "routing"}){
			final Path p=root.resolve(key.equals("in")?"inputs.tsv":key+".tsv");
			args.add(key+"="+p);
			args.add((key.equals("loci")?"locisha80":key+"sha80")+"="+(badPin && key.equals("in")?(structural?"invalid":"00000000000000000000"):pin(p)));
		}
		if(structural){args.add("verification=structural");}
		args.addAll(Arrays.asList(options));
		args.add("out="+root.resolve(name));
		CallGenesLengthScorecard.main(args.toArray(new String[0]));
	}
	/** Returns true only for an argument failure that leaves no partial output directory. */
	private static boolean fails(Path root, String name, boolean badPin)throws Exception{
		try{run(root, name, badPin);}catch(IllegalArgumentException e){
			check(!Files.exists(root.resolve(name)), "Failed fixture emits no partial result: "+name);
			return true;
		}
		return false;
	}
	/** Changes one fixture cell, requires rejection, then restores the file and its manifest bindings. */
	private static void mutateFail(Path root, String file, int row, int column, String value, String name)throws Exception{
		final Path p=root.resolve(file);
		final String before=contents(p);
		final String[] lines=before.split("\n"), cells=lines[row].split("\t", -1);
		cells[column]=value;
		lines[row]=String.join("\t", cells);
		write(p, String.join("\n", lines)+"\n");
		manifests(root);
		check(fails(root, name, false), "Corrupt "+name+" fails closed");
		write(p, before);
		manifests(root);
	}
	/** Drops one zero-based line, including the header in line numbering, and retains a final newline. */
	private static String dropRow(String text, int row){
		final ArrayList<String> rows=new ArrayList<String>(Arrays.asList(text.split("\n")));
		rows.remove(row);
		return String.join("\n", rows)+"\n";
	}
	/** Finds the aggregate row for a policy and length band or fails with the missing identity. */
	private static String[] row(List<String> rows, String policy, String band){
		for(String s:rows){
			final String[] c=s.split("\t");
			if(c[0].equals(policy) && c[1].equals("ALL") && c[2].equals(band)){return c;}
		}
		throw new AssertionError("Missing row "+policy+"/"+band);
	}
	/** Parses all numeric cells starting at the specified column. */
	private static long[] numbers(String[] cells, int from){
		final long[] n=new long[cells.length-from];
		for(int i=0; i<n.length; i++){n[i]=Long.parseLong(cells[i+from]);}
		return n;
	}
	/** Encodes a fixture Boolean as its literal zero/one table value. */
	private static int bit(boolean b){return b?1:0;}
	/** Formats one tab-separated fixture row with a final newline. */
	private static String tsv(Object... cells){
		final StringBuilder b=new StringBuilder();
		for(int i=0; i<cells.length; i++){
			if(i>0){b.append('\t');}
			b.append(cells[i]);
		}
		return b.append('\n').toString();
	}
	/** Creates a one-based, phase-zero CDS row on fixture contig c. */
	private static String g(int start, int stop, int strand, String attrs){return "c\tfixture\tCDS\t"+start+"\t"+stop+"\t.\t"+(strand==0?"+":"-")+"\t0\t"+attrs+"\n";}
	/** Parses one literal GFF row after trimming its trailing newline. */
	private static GffLine line(String text){return new GffLine(text.trim().getBytes(StandardCharsets.US_ASCII));}
	/** Writes or replaces exact ASCII fixture text. */
	private static void write(Path p, String text)throws Exception{Files.write(p, text.getBytes(StandardCharsets.US_ASCII));}
	/** Reads complete ASCII fixture text for mutation and restoration. */
	private static String contents(Path p)throws Exception{return new String(Files.readAllBytes(p), StandardCharsets.US_ASCII);}
	/** Reads an ASCII result table into ordered lines. */
	private static List<String> read(Path p)throws Exception{return Files.readAllLines(p, StandardCharsets.US_ASCII);}
	/** Returns an explicit unmeasured marker in structural mode, otherwise the exact fixture identity. */
	private static String pin(Path p)throws Exception{return structural?"unmeasured":DigestSuffix.file(p.toString());}
	/** Counts an attempted expectation and fails immediately if it is false. */
	private static void check(boolean ok, String why){
		checks++;
		if(!ok){throw new AssertionError(why);}
	}
	/** Ordered fixture policies and their seven elementary length bands. */
	private static final String[] POLICIES={"native", "cap150", "k5"}, BANDS={"lt40", "40_64", "65_99", "100_149", "150_199", "ge200", "unknown"};
	/** Verification mode selected once by main for all subsequent fixture invocations. */
	private static boolean structural;
	/** Attempted checks in this process; main does not reset the counter. */
	private static int checks;
}
