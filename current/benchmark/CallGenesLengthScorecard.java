package benchmark;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;

import benchmark.CallGenesStopAudit.ReferenceSet;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import gff.GffLine;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import benchmark.CallGenesBenchmarkSupport.InputVerification;
import prok.ProkObject;
import shared.Shared;
import structures.ByteBuilder;
import structures.StringNum;

/** Saved-call grades by frozen annotation length, with called-length extras and
 * observed per-set routing. No caller, classifier or reference filtering runs.
 * @author Sayu
 */
public final class CallGenesLengthScorecard{

	public static void main(String[] args)throws Exception{
		final PreParser pp=new PreParser(args, CallGenesLengthScorecard.class, false);
		try{
			final HashMap<String, String> o=new HashMap<String, String>();
			for(String arg:pp.args){
				final String[] kv=Parse.splitOnFirst(arg, '=');
				require(kv.length==2 && !kv[1].isEmpty()
					&& o.put(kv[0].toLowerCase(Locale.ROOT), kv[1])==null,
					"Unique key=value options required");
			}
			final String[] keys={"in", "insha80", "expected", "expectedsha80", "loci", "locisha80",
				"routing", "routingsha80", "out"};
			final HashSet<String> allowed=new HashSet<String>(Arrays.asList(keys));
			allowed.add("verification");
			allowed.add("comparators");
			require(allowed.containsAll(o.keySet()) && o.keySet().containsAll(Arrays.asList(keys)),
				"Explicit saved input/control/locus/routing pins and fresh output required");
			final String comparatorChoice=o.getOrDefault("comparators", "native,cap150");
			require(comparatorChoice.equals("native") || comparatorChoice.equals("native,cap150"),
				"comparators must be native or native,cap150");
			Shared.setThreads(1);
			CallGenesBenchmarkSupport.cdsOnly();
			new CallGenesLengthScorecard().run(o);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Read the existing locus schema once for a caller benchmark's exact cohort.
	 * The caller binds the same FASTA/reference/code/scanner contract that produced
	 * these bands; PGM or formula policy alone does not define annotation length. */
	public static SavedBands readSavedBands(Path file, java.util.Set<String> genomes){
		require(file!=null && genomes!=null && !genomes.isEmpty(),
			"Saved bands require an explicit nonempty cohort");
		final HashMap<String, HashMap<StringNum, Locus>> loci=loadLoci(file);
		require(loci.keySet().equals(genomes),
			"Saved bands and caller cohort differ;100-genome loci cannot grade983 genomes");
		return new SavedBands(loci);
	}

	/** Read-only after construction. No locus or mutable map escapes this facade. */
	public static final class SavedBands{
		private SavedBands(HashMap<String, HashMap<StringNum, Locus>> loci_){loci=loci_;}
		/** Verify the exact saved winners once; their feature objects remain read-only. */
		public ReferenceBands bind(String genome, ReferenceSet refs){
			final HashMap<StringNum, Locus> bands=loci.get(genome);
			require(refs!=null && bands!=null,
				"Reference winners must belong to the bound saved cohort: "+genome);
			verifyLoci(refs, bands);
			return new ReferenceBands(refs, bands);
		}
		private final HashMap<String, HashMap<StringNum, Locus>> loci;
	}

	/** One genome's immutable reference binding; callers must not mutate features.
	 * Grade after the calling barrier while CDS is enabled; no global flags change. */
	public static final class ReferenceBands{
		private ReferenceBands(ReferenceSet refs_, HashMap<StringNum, Locus> bands_){
			refs=refs_;
			bands=bands_;
		}
		/** Frozen annotation band of an original effective winner, including retained
		 * subsets. Never derive annotation length from a call or adjusted endpoint. */
		public int annotationBand(GffLine ref){
			require(ref!=null && refs.lookup(ref.seqid, ref.trueStop())==ref,
				"Annotation band requires an original bound reference winner");
			return bands.get(key(ref)).band;
		}
		/** Grade strict CDS calls after caller workers join. Whole-run conservation
		 * belongs to the harness; this path does not run a second parity grader. */
		public ScoreCounts grade(List<GffLine> calls){
			require(calls!=null && ProkObject.callCDS,
				"Lean grading requires CDS enabled and an explicit strict-CDS call list");
			for(GffLine call:calls){
				require(call!=null && "CDS".equals(call.type),
					"Filter RNA and other features before lean CDS grading");
			}
			final long[][] counts=gradeBands(refs, calls, bands);
			return new ScoreCounts(counts);
		}
		private final ReferenceSet refs;
		private final HashMap<StringNum, Locus> bands;
	}

	/** Private counters retain the unchanged16-count formatting and band algebra. */
	public static final class ScoreCounts{
		private ScoreCounts(long[][] counts_){counts=counts_;}
		/** Format all nine bands with the existing16-count schema; the caller owns the returned buffer. */
		public ByteBuilder rows(String policy, String genome){
			require(policy!=null && policy.matches("[A-Za-z][A-Za-z0-9_]*") && genome!=null
				&& genome.matches("[A-Za-z0-9_-]+"), "Score labels must be safe table keys");
			final ByteBuilder rows=new ByteBuilder();
			score(rows, policy, genome, counts);
			return rows;
		}
		/** Return a defensive copy of the aggregate eight-field grade. */
		public long[] grade(){return Arrays.copyOf(sum(counts), 8);}
		private final long[][] counts;
	}
	/** Shared header for both guarded and lean16-count tables. */
	public static String scoreHeader(){return SCORE_HEADER;}
	/** Sum completed genome counters without mutating any component result. */
	public static ScoreCounts combineScores(List<ScoreCounts> scores){
		require(scores!=null && !scores.isEmpty(), "Aggregate counts require completed genome results");
		final long[][] total=new long[7][9];
		for(ScoreCounts score:scores){
			require(score!=null, "A missing genome result must fail whole-run conservation");
			add(total, score.counts);
		}
		return new ScoreCounts(total);
	}

	/** One genome's reference winners and six call sets are resident at a time. */
	private void run(Map<String, String> o)throws Exception{
		final String[] comparators=o.getOrDefault("comparators", "native,cap150").split(",");
		verification=new InputVerification(o.get("verification"));
		final Path input=Paths.get(o.get("in")), expectedFile=Paths.get(o.get("expected")),
			lociFile=Paths.get(o.get("loci")), routingFile=Paths.get(o.get("routing")),
			out=Paths.get(o.get("out"));
		require(!Files.exists(out) && Files.isDirectory(out.toAbsolutePath().getParent()),
			"Fresh output under an existing directory required");
		pin(input, o.get("insha80"));
		pin(expectedFile, o.get("expectedsha80"));
		pin(lociFile, o.get("locisha80"));
		pin(routingFile, o.get("routingsha80"));
		final LinkedHashMap<String, long[]> expected=expected(expectedFile);
		for(String comparator:comparators){
			require(expected.containsKey(comparator), "Missing required comparator: "+comparator);
		}
		final TreeMap<String, LinkedHashMap<String, Input>> roster=inputs(input, expected);
		final HashMap<String, HashMap<StringNum, Locus>> loci=loadLoci(lociFile);
		require(loci.keySet().equals(roster.keySet()), "Frozen loci and input genome roster differ");
		final HashMap<String, RouteInput> routeInputs=routingInputs(routingFile, roster);
		final LinkedHashMap<String, long[][]> totals=new LinkedHashMap<String, long[][]>(),
			routeTotals=new LinkedHashMap<String, long[][]>();
		for(String policy:expected.keySet()){
			totals.put(policy, new long[7][9]);
			routeTotals.put(policy, new long[7][18]);
		}
		final ByteBuilder byGenome=new ByteBuilder(SCORE_HEADER), routingRows=new ByteBuilder(ROUTING_HEADER);
		long referenceCount=0;
		for(Map.Entry<String, LinkedHashMap<String, Input>> entry:roster.entrySet()){
			final String genome=entry.getKey();
			final LinkedHashMap<String, Input> arms=entry.getValue();
			final Input nativeInput=arms.get("native");
			pin(nativeInput.reference, nativeInput.referencePin);
			final ReferenceSet refs=CallGenesStopAudit.loadReference(nativeInput.reference.toString());
			final HashMap<StringNum, Locus> bands=loci.get(genome);
			verifyLoci(refs, bands);
			referenceCount+=refs.uniqueStops();
			final LinkedHashMap<String, List<GffLine>> calls=new LinkedHashMap<String, List<GffLine>>();
			for(String policy:expected.keySet()){
				final Input in=arms.get(policy);
				pin(in.calls, in.callsPin);
				pin(in.grades, in.gradesPin);
				calls.put(policy, GffLine.loadGffFile(in.calls.toString(), null, false));
			}
			final HashSet<StringNum> nativeStops=calledStops(calls.get("native"));
			for(String policy:expected.keySet()){
				final Input in=arms.get(policy);
				final List<GffLine> query=calls.get(policy);
				final long[] grade=CallGenesStopAudit.grade(refs, query),
					saved=grades(in.grades).get(genome+"/"+in.arm);
				require(saved!=null && Arrays.equals(grade, saved),
					"Original eight-field grade changed: "+policy+"/"+genome);
				final long[][] split=gradeBands(refs, query, bands);
				require(Arrays.equals(grade, Arrays.copyOf(sum(split), 8)),
					"Band counts do not conserve native grade");
				add(totals.get(policy), split);
				score(byGenome, policy, genome, split);
				final long[][] routes;
				if(policy.equals("native")){
					routes=new long[7][18];
					for(GffLine ref:refs.records()){
						final long[] n=routes[bands.get(key(ref)).band];
						n[0]++;
						n[1]++;
						if(nativeStops.contains(routeKey(ref.seqid, ref.trueStop(), ref.strand))){
							n[14]++;
							n[15]++;
						}
					}
				}else{
					routes=routes(genome, policy, in.arm, refs, bands, routeInputs.get(policy+"/"+genome),
						nativeStops, calledStops(query));
				}
				add(routeTotals.get(policy), routes);
				routing(routingRows, policy, genome, routes);
			}
			System.out.println("LENGTH_SCORECARD_GENOME_PASS genome="+genome+" references="+
				refs.uniqueStops()+" policies="+expected.size());
		}
		final ByteBuilder summary=new ByteBuilder(SCORE_HEADER),
			routeSummary=new ByteBuilder(ROUTING_HEADER), deltas=new ByteBuilder(DELTA_HEADER);
		for(String policy:expected.keySet()){
			final long[] g=sum(totals.get(policy)), want=expected.get(policy);
			final long[] measured={roster.size(), g[8], g[1], g[3], g[7], g[3]+g[7], g[2]+g[6]};
			require(Arrays.equals(want, measured) && g[8]==referenceCount,
				"Independent aggregate control failed: "+policy);
			score(summary, policy, "ALL", totals.get(policy));
			routing(routeSummary, policy, "ALL", routeTotals.get(policy));
			for(String comparator:comparators){
				for(int b=0;b<9;b++){
					final long[] a=view(totals.get(policy), b), c=view(totals.get(comparator), b);
					deltas.append(policy).tab().append(comparator).tab().append(label(b))
						.tab().append(a[1]-c[1]).tab().append(a[3]-c[3])
						.tab().append(a[7]-c[7]).tab().append(a[3]+a[7]-c[3]-c[7])
						.tab().append(a[2]+a[6]-c[2]-c[6]).nl();
				}
			}
		}
		verification.verifyComplete();
		Files.createDirectory(out);
		write(out.resolve("scorecard.tsv"), summary);
		write(out.resolve("by_genome.tsv"), byGenome);
		write(out.resolve("routing.tsv"), routeSummary);
		write(out.resolve("routing_by_genome.tsv"), routingRows);
		write(out.resolve("deltas.tsv"), deltas);
		final String pass="LENGTH_SCORECARD_PASS policies="+expected.size()+" genomes="+roster.size()+" reference_keys="+referenceCount+" grade_fields=8 loci_reads=1 caller_runs=0 classifier_runs=0 verification="+verification.mode()+"\n";
		write(out.resolve("PASS"), new ByteBuilder(pass));
		System.out.print(pass);
	}

	/** Preserve every original winner; the saved annotation band is never recomputed by arm. */
	private static void verifyLoci(ReferenceSet refs, HashMap<StringNum, Locus> loci){
		require(loci.size()==refs.uniqueStops(), "Exactly one saved locus per effective reference required");
		for(GffLine ref:refs.records()){
			final Locus l=loci.get(key(ref));
			require(l!=null && ref.type.equals("CDS") && l.strand==ref.strand && l.left==ref.start && l.right==ref.stop
				&& l.start==ref.trueStart() && l.phase==ref.phase,
				"Saved locus differs from original winning reference");
		}
	}

	/** TP values count occurrences; FN values count missing unique reference keys. */
	private static long[][] gradeBands(ReferenceSet refs, List<GffLine> calls,
		HashMap<StringNum, Locus> bands){
		final long[][] counts=new long[7][9];
		final HashSet<StringNum> starts=new HashSet<StringNum>(), stops=new HashSet<StringNum>();
		for(GffLine q:calls){
			if(!ProkObject.processType(q.prokType())){continue;}
			if(!refs.matches(q)){
				final long[] n=counts[calledBand(q)];
				n[6]++;
				n[7]++;
				continue;
			}
			final StringNum k=key(q);
			final GffLine ref=refs.lookup(q.seqid, q.trueStop());
			final long[] n=counts[bands.get(k).band];
			stops.add(k);
			n[1]++;
			n[5]++;
			if(q.trueStart()==ref.trueStart()){
				starts.add(k);
				n[0]++;
				n[4]++;
			}else{n[6]++;}
		}
		for(GffLine ref:refs.records()){
			final StringNum k=key(ref);
			final long[] n=counts[bands.get(k).band];
			n[8]++;
			if(!starts.contains(k)){n[2]++;}
			if(!stops.contains(k)){n[3]++;}
		}
		return counts;
	}

	/** Called CDS span excludes its terminating codon. Incomplete or nontriplet
	 * calls cannot establish that full protein length and remain explicitly unknown. */
	public static int calledBand(GffLine call){
		final int bases=call.length();
		if(!call.type.equals("CDS") || call.partial() || call.phase!=0 || bases<3 || bases%3!=0){return 6;}
		final int aa=bases/3-1;
		return aa<40 ? 0 : aa<65 ? 1 : aa<100 ? 2 : aa<150 ? 3 : aa<200 ? 4 : 5;
	}

	/** Reads all saved routes but counts only unique original winning reference stops. */
	private long[][] routes(String genome, String policy, String gradeArm, ReferenceSet refs,
		HashMap<StringNum, Locus> bands, RouteInput files, HashSet<StringNum> nativeStops,
		HashSet<StringNum> calls)throws Exception{
		require(files!=null, "Missing saved routing binding");
		for(int i=0;i<3;i++){pin(files.paths[i], files.hashes[i]);}
		final HashMap<StringNum, Route> routes=new HashMap<StringNum, Route>();
		try(Table t=new Table(files.paths[0], STOPS_HEADER)){
			while(t.next()){
				if(!bonus(t, genome, gradeArm)){continue;}
				final Route r=new Route();
				r.contig=t.s(1);
				r.strand=t.i(4);
				r.stop=trueStop(t.i(3), r.strand);
				r.upstream=t.i(5);
				r.accepted=t.bit(6);
				r.annotated=t.bit(7);
				r.nativeCall=t.bit(11);
				r.call=t.bit(12);
				r.query=t.i(13);
				require(r.query>=40 && r.query<=files.cap && r.query<=r.upstream,
					"Saved query outside bound view");
				final StringNum k=routeKey(r.contig, r.stop, r.strand);
				require(routes.put(k, r)==null, "Duplicate selected-arm stop row");
				final GffLine ref=refs.lookup(r.contig, r.stop);
				require(r.annotated==(ref!=null && ref.strand==r.strand && ref.type.equals("CDS")),
					"Saved annotation membership differs from native reference");
				require(r.nativeCall==nativeStops.contains(k) && r.call==calls.contains(k),
					"Saved call flags differ from actual GFF stops");
			}
		}
		try(Table t=new Table(files.paths[2], SETS_HEADER)){
			while(t.next()){
				if(!bonus(t, genome, gradeArm)){continue;}
				final int strand=t.i(4);
				final String name=firstToken(t.s(1));
				final Route r=routes.get(routeKey(name, trueStop(t.i(3), strand), strand));
				require(r!=null, "Per-set row has no matching saved stop");
				final boolean shortSet=t.s(5).equals("short");
				require(shortSet || t.s(5).equals("long"), "Unknown family set");
				final int expected=shortSet ? Math.min(99, r.query) : r.query;
				require(t.i(6)==expected && (shortSet || r.query>=100) && t.i(7)==r.upstream,
					"Per-set query geometry mismatch");
				require(t.bit(11)==r.annotated && t.bit(12)==r.nativeCall && t.bit(13)==r.call,
					"Per-set membership/call flags changed");
				final boolean accepted=t.bit(8);
				require(accepted ? !t.s(9).equals("NA") : t.s(9).equals("NA"),
					"Family identity and acceptance disagree");
				final String primary=t.s(10);
				require(primary.equals("none") || primary.equals("short") || primary.equals("long"),
					"Unknown primary source");
				require(r.primary==null || r.primary.equals(primary), "Primary differs across per-set rows");
				r.primary=primary;
				if(shortSet){
					require(!r.shortSeen, "Duplicate short-set row");
					r.shortSeen=true;
					r.shortAccept=accepted;
				}else{
					require(!r.longSeen, "Duplicate long-set row");
					r.longSeen=true;
					r.longAccept=accepted;
				}
			}
		}
		for(Route r:routes.values()){
			require(r.shortSeen && r.longSeen==(r.query>=100), "Missing applicable family-set search");
			require(r.accepted==(r.shortAccept || r.longAccept), "Union acceptance mismatch");
			final String primary=r.longAccept ? "long" : r.shortAccept ? "short" : "none";
			require(primary.equals(r.primary), "Recorded primary violates long-first acceptance precedence");
		}
		final long[][] counts=new long[7][18];
		final HashSet<StringNum> seen=new HashSet<StringNum>();
		try(Table t=new Table(files.paths[1], ANNOTATIONS_HEADER)){
			while(t.next()){
				if(!bonus(t, genome, gradeArm)){continue;}
				final String name=t.s(1);
				final int stop=t.i(3), strand=t.i(4);
				final GffLine ref=refs.lookup(name, stop);
				final StringNum k=new StringNum(name, stop);
				require(ref!=null && ref.strand==strand && ref.start==t.i(5) && ref.stop==t.i(6)
					&& seen.add(k), "Annotation must join exactly one original winner");
				final StringNum rk=routeKey(name, stop, strand);
				final Route r=routes.get(rk);
				final boolean scanned=t.bit(11), accepted=t.bit(12);
				require(scanned==(r!=null) && scanned==t.s(10).equals("scanned")
					&& accepted==(r!=null && r.accepted),
					"Annotation scan/acceptance disagrees with saved stop");
				final boolean nativeCall=t.bit(16), called=t.bit(17);
				require(nativeCall==nativeStops.contains(rk) && called==calls.contains(rk),
					"Annotation called flags differ from saved GFF");
				final long[] n=counts[bands.get(k).band];
				n[0]++;
				if(r==null){
					n[2]++;
				}else{
					n[3]++;
					n[4]+=bit(r.shortSeen);
					n[5]+=bit(r.longSeen);
					n[6]+=bit(r.shortSeen && r.longSeen);
					n[7]+=bit(r.shortAccept);
					n[8]+=bit(r.longAccept);
					n[9]+=bit(r.shortAccept && r.longAccept);
					n[10]+=bit(r.accepted);
					n[11]+=bit(!r.accepted);
					n[12]+=bit(r.primary.equals("short"));
					n[13]+=bit(r.primary.equals("long"));
				}
				n[14]+=bit(nativeCall);
				n[15]+=bit(called);
				n[16]+=bit(!nativeCall && called);
				n[17]+=bit(nativeCall && !called);
			}
		}
		require(seen.size()==refs.uniqueStops(), "Saved annotation routing omitted reference winners");
		return counts;
	}

	/** Zero-effect replay rows are validated for identity then skipped, never counted twice. */
	private static boolean bonus(Table t, String genome, String gradeArm){
		require(t.s(0).equals(genome) && (t.s(2).equals("bonus400") || t.s(2).equals("formula")
			|| t.s(2).equals("zero_effect")), "Unexpected diagnostic identity/arm");
		return t.s(2).equals(gradeArm);
	}
	public static int trueStop(int stopLeft, int strand){
		require(stopLeft>=0 && (strand==0 || strand==1), "Invalid scanner stop coordinate/strand");
		return Math.addExact(stopLeft, strand==0 ? 3 : 1);
	}
	private static String firstToken(String s){
		int i=0;
		while(i<s.length() && !Character.isWhitespace(s.charAt(i))){i++;}
		require(i>0, "Empty contig");
		return s.substring(0, i);
	}
	private static StringNum key(GffLine g){return new StringNum(g.seqid, g.trueStop());}
	private static StringNum routeKey(String contig, int stop, int strand){
		require(stop>0 && (strand==0 || strand==1), "Invalid genomic stop");
		return new StringNum(contig, ((long)stop<<1)|strand);
	}
	private static HashSet<StringNum> calledStops(List<GffLine> calls){
		final HashSet<StringNum> result=new HashSet<StringNum>();
		for(GffLine q:calls){
			if(q.type.equals("CDS")){result.add(routeKey(q.seqid, q.trueStop(), q.strand));}
		}
		return result;
	}

	private TreeMap<String, LinkedHashMap<String, Input>> inputs(Path file,
		Map<String, long[]> expected){
		final TreeMap<String, LinkedHashMap<String, Input>> result=
			new TreeMap<String, LinkedHashMap<String, Input>>();
		final HashMap<String, String> policyArms=new HashMap<String, String>();
		try(Table t=new Table(file, INPUT_HEADER)){
			while(t.next()){
				final String policy=t.s(0), genome=t.s(1);
				require(expected.containsKey(policy) && genome.matches("g[0-9]{4}"),
					"Unexpected policy/genome");
				final Input in=new Input(t);
				require(policy.equals("native") ? in.arm.equals("native") :
					in.arm.equals("bonus400") || !policy.equals("cap150") && in.arm.equals("formula"),
					"Wrong saved grade arm");
				final String priorArm=policyArms.putIfAbsent(policy, in.arm);
				require(priorArm==null || priorArm.equals(in.arm),
					"Saved grade arm changed within one policy");
				LinkedHashMap<String, Input> map=result.get(genome);
				if(map==null){
					map=new LinkedHashMap<String, Input>();
					result.put(genome, map);
				}
				require(map.put(policy, in)==null, "Duplicate policy/genome input");
			}
		}
		for(LinkedHashMap<String, Input> map:result.values()){
			require(map.keySet().equals(expected.keySet()), "Incomplete policy roster");
			final Input n=map.get("native");
			for(Input in:map.values()){
				require(in.reference.equals(n.reference) && in.referencePin.equals(n.referencePin),
					"References differ across arms");
			}
		}
		for(long[] e:expected.values()){require(e[0]==result.size(), "Expected genome count differs");}
		return result;
	}
	private static LinkedHashMap<String, long[]> expected(Path file){
		final LinkedHashMap<String, long[]> result=new LinkedHashMap<String, long[]>();
		try(Table t=new Table(file, EXPECTED_HEADER)){
			while(t.next()){
				final String name=t.s(0);
				require(name.matches("[A-Za-z][A-Za-z0-9_]*"), "Unsafe policy name");
				final long[] n=new long[7];
				for(int i=0;i<7;i++){
					n[i]=t.n(i+1);
					require(n[i]>=0, "Negative expected count");
				}
				require(n[0]>0 && result.put(name, n)==null, "Duplicate policy or empty population");
			}
		}
		return result;
	}
	private static HashMap<String, HashMap<StringNum, Locus>> loadLoci(Path file){
		final HashMap<String, HashMap<StringNum, Locus>> result=
			new HashMap<String, HashMap<StringNum, Locus>>();
		try(Table t=new Table(file, null)){
			require(t.width==42, "Accepted locus schema must have42 fields");
			while(t.next()){
				final String genome=t.named("genome");
				require(genome.matches("g[0-9]{4}"), "Unsafe locus genome");
				final Locus l=new Locus(t);
				HashMap<StringNum, Locus> map=result.get(genome);
				if(map==null){
					map=new HashMap<StringNum, Locus>();
					result.put(genome, map);
				}
				require(map.put(new StringNum(l.contig, l.stop), l)==null,
					"Duplicate frozen reference winner");
			}
		}
		return result;
	}
	private static HashMap<String, RouteInput> routingInputs(Path file,
		Map<String, LinkedHashMap<String, Input>> roster){
		final HashMap<String, RouteInput> result=new HashMap<String, RouteInput>();
		final HashMap<String, Integer> policyCaps=new HashMap<String, Integer>();
		try(Table t=new Table(file, ROUTE_INPUT_HEADER)){
			while(t.next()){
				final String policy=t.s(0), genome=t.s(1);
				require(!policy.equals("native") && roster.containsKey(genome)
					&& roster.get(genome).containsKey(policy), "Unexpected routing identity");
				// The pinned manifest binds geometry; policy IDs identify experiments, not cap values.
				final RouteInput in=new RouteInput(t);
				require(in.cap==100 || in.cap==150, "Unsupported query cap binding");
				require(!policy.equals("cap150") || in.cap==150, "cap150 comparator requires query cap150");
				final Integer previousCap=policyCaps.putIfAbsent(policy, in.cap);
				require(previousCap==null || previousCap.intValue()==in.cap,
					"Inconsistent query cap for policy: "+policy);
				require(result.put(policy+"/"+genome, in)==null, "Duplicate routing input");
			}
		}
		for(Map.Entry<String, LinkedHashMap<String, Input>> e:roster.entrySet()){
			for(String policy:e.getValue().keySet()){
				require(policy.equals("native") || result.containsKey(policy+"/"+e.getKey()),
					"Missing routing binding");
			}
		}
		return result;
	}
	private HashMap<String, long[]> grades(Path file){
		HashMap<String, long[]> result=gradeCache.get(file);
		if(result!=null){return result;}
		result=new HashMap<String, long[]>();
		try(Table t=new Table(file, GRADE_HEADER)){
			while(t.next()){
				final long[] n=new long[8];
				for(int i=0;i<8;i++){
					n[i]=t.n(i+2);
					require(n[i]>=0, "Negative saved grade");
				}
				require(result.put(t.s(0)+"/"+t.s(1), n)==null, "Duplicate saved grade");
			}
		}
		gradeCache.put(file, result);
		return result;
	}
	private void pin(Path path, String hash)throws Exception{verification.pin(path, hash);}

	private static void score(ByteBuilder b, String policy, String genome, long[][] counts){
		for(int band=0;band<9;band++){
			final long[] n=view(counts, band);
			final long matched=n[8]-n[3];
			require(n[0]==n[4] && n[1]==n[5] && n[2]<=n[8] && matched>=0 && matched<=n[1],
				"Invalid band grade conservation");
			b.append(policy).tab().append(genome).tab().append(label(band)).tab().append(n[8])
				.tab().append(n[1]).tab().append(n[3]).tab().append(n[7]).tab().append(n[3]+n[7])
				.tab().append(n[2]+n[6]).tab().append(matched).tab().append(n[1]-matched);
			for(int i=0;i<8;i++){b.tab().append(n[i]);}
			b.nl();
		}
	}
	private static void routing(ByteBuilder b, String policy, String genome, long[][] counts){
		for(int band=0;band<9;band++){
			final long[] n=view(counts, band);
			require(n[0]==n[1]+n[2]+n[3] && n[3]==n[10]+n[11] && n[7]<=n[4] && n[8]<=n[5]
				&& n[9]<=n[6] && n[10]==n[7]+n[8]-n[9] && n[10]==n[12]+n[13],
				"Routing set/key conservation failed");
			b.append(policy).tab().append(genome).tab().append(label(band));
			for(long x:n){b.tab().append(x);}
			b.nl();
		}
	}
	private static String label(int band){return band<7 ? BANDS[band] : band==7 ? "ALL" : "40_99";}
	private static long[] view(long[][] counts, int band){
		if(band<7){return counts[band];}
		if(band==7){return sum(counts);}
		final long[] n=counts[1].clone();
		add(n, counts[2]);
		return n;
	}
	private static long[] sum(long[][] rows){
		final long[] n=new long[rows[0].length];
		for(long[] row:rows){add(n, row);}
		return n;
	}
	private static void add(long[][] a, long[][] b){
		require(a.length==b.length, "Counter bands differ");
		for(int i=0;i<a.length;i++){add(a[i], b[i]);}
	}
	private static void add(long[] a, long[] b){
		require(a.length==b.length, "Counter shapes differ");
		for(int i=0;i<a.length;i++){a[i]=Math.addExact(a[i], b[i]);}
	}
	private static int bit(boolean b){return b ? 1 : 0;}
	private static void write(Path path, ByteBuilder data){
		final ByteStreamWriter w=new ByteStreamWriter(path.toString(), false, false, false);
		w.start();
		w.print(data);
		require(!w.poisonAndWait(), "Output write failed");
	}
	private static void require(boolean b, String why){if(!b){throw new IllegalArgumentException(why);}}

	private static final class Input{
		Input(Table t){
			reference=Paths.get(t.s(2));
			referencePin=t.s(3);
			calls=Paths.get(t.s(4));
			callsPin=t.s(5);
			grades=Paths.get(t.s(6));
			gradesPin=t.s(7);
			arm=t.s(8);
		}
		final Path reference, calls, grades;
		final String referencePin, callsPin, gradesPin, arm;
	}
	private static final class RouteInput{
		RouteInput(Table t){
			for(int i=0;i<3;i++){
				paths[i]=Paths.get(t.s(2+2*i));
				hashes[i]=t.s(3+2*i);
			}
			cap=t.i(8);
		}
		final Path[] paths=new Path[3];
		final String[] hashes=new String[3];
		final int cap;
	}
	private static final class Locus{
		Locus(Table t){
			contig=t.named("contig");
			stop=t.integer("true_stop_1based");
			strand=t.integer("strand");
			start=t.integer("reference_start_1based");
			left=t.integer("reference_left_1based");
			right=t.integer("reference_right_1based");
			phase=t.integer("phase");
			final String value=t.named("annotation_band");
			int b=-1;
			for(int i=0;i<BANDS.length;i++){
				if(BANDS[i].equals(value)){b=i;}
			}
			require(b>=0, "Unknown frozen annotation band: "+value);
			band=b;
		}
		final String contig;
		final int stop, strand, start, left, right, phase, band;
	}
	private static final class Route{
		String contig, primary;
		int stop, strand, query, upstream;
		boolean accepted, annotated, nativeCall, call, shortSeen, longSeen, shortAccept, longAccept;
	}
	private static final class Table implements AutoCloseable{
		Table(Path path, String expected){
			in=ByteFile.makeByteFile1(path.toString(), false);
			try{
				final byte[] h=in.nextLine();
				require(h!=null, "Missing header: "+path);
				final String header=new String(h, StandardCharsets.US_ASCII);
				require(expected==null || expected.equals(header), "Wrong table schema: "+path);
				final String[] names=header.split("\t", -1);
				width=names.length;
				for(int i=0;i<width;i++){require(columns.put(names[i], i)==null, "Duplicate header field");}
			}catch(RuntimeException|Error e){
				in.close();
				throw e;
			}
		}
		boolean next(){
			final byte[] row=in.nextLine();
			if(row==null){return false;}
			parser.set(row);
			require(parser.terms()==width, "Malformed table width");
			return true;
		}
		String s(int i){return parser.parseString(i);}
		int i(int n){return parser.parseInt(n);}
		long n(int i){return parser.parseLong(i);}
		boolean bit(int i){
			final String v=s(i);
			require(v.equals("0") || v.equals("1"), "Expected boolean0/1");
			return v.equals("1");
		}
		String named(String key){
			require(columns.containsKey(key), "Missing column "+key);
			return s(columns.get(key));
		}
		int integer(String key){
			require(columns.containsKey(key), "Missing column "+key);
			return i(columns.get(key));
		}
		@Override public void close(){require(!in.close(), "Input read failed");}
		final ByteFile in;
		final int width;
		final LineParser1 parser=new LineParser1('\t');
		final HashMap<String, Integer> columns=new HashMap<String, Integer>();
	}
	private InputVerification verification;
	private final HashMap<Path, HashMap<String, long[]>> gradeCache=new HashMap<Path, HashMap<String, long[]>>();
	private static final String[] BANDS={"lt40", "40_64", "65_99", "100_149", "150_199", "ge200", "unknown"};
	public static final String INPUT_HEADER="policy\tgenome\treference\treference_sha80\tcalls\tcalls_sha80\tgrades\tgrades_sha80\tgrade_arm";
	public static final String EXPECTED_HEADER="policy\tgenomes\treference_keys\tcorrect\tmissed\textras\tstop_errors\tstart_errors";
	public static final String ROUTE_INPUT_HEADER="policy\tgenome\tstops\tstops_sha80\tannotations\tannotations_sha80\tfamily_sets\tfamily_sets_sha80\tquery_cap";
	public static final String STOPS_HEADER="genome\tcontig\tarm\tstop_left_0based\tstrand\tupstream_codons\taccepted\tannotated\tparent\tadmitted_candidates\tretained_candidates\tcalled_native\tcalled_arm\tquery_length\tinterval_boundary";
	public static final String SETS_HEADER="genome\tcontig\tarm\tstop_left_0based\tstrand\tfamily_set\tquery_length\tupstream_codons\taccepted\tfamily_id\tprimary_source\tannotated\tcalled_native\tcalled_arm\traw_hit_floor\tadmitted_families\tdp_attempts";
	public static final String ANNOTATIONS_HEADER="genome\tcontig\tarm\ttrue_stop_1based\tstrand\tstart_1based\tend_1based\tcds_codons_including_stop\tpartial\tcontig_gc\tscanner_status\tscanned\taccepted\tparent\tadmitted_candidates\tretained_candidates\tcalled_native\tcalled_arm";
	private static final String GRADE_HEADER="genome\tarm\t"+String.join("\t", CallGenesBenchmarkSupport.GRADE_NAMES);
	private static final String SCORE_HEADER="policy\tgenome\tannotation_or_called_length_band\treference_keys\tcorrect\tmissed\textras\tstop_errors\tstart_errors\tunique_matched_keys\tduplicate_matching_occurrences\t"+String.join("\t", CallGenesBenchmarkSupport.GRADE_NAMES)+"\n";
	private static final String ROUTING_HEADER="policy\tgenome\tannotation_band\treference_stops\tnative_no_family\tunscanned\tscanned\tshort_searched\tlong_searched\tboth_searched\tshort_accepted\tlong_accepted\tboth_accepted\tany_accepted\tsearched_rejected\tprimary_short\tprimary_long\tnative_called\tcalled\trecovered\tlost\n";
	private static final String DELTA_HEADER="policy\tcomparator\tannotation_or_called_length_band\tcorrect_delta\tmissed_delta\textras_delta\tstop_errors_delta\tstart_errors_delta\n";
}
