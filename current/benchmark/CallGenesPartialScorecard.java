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

import benchmark.CallGenesStopAudit.ReferenceSet;
import benchmark.CallGenesLengthScorecard.ReferenceBands;
import benchmark.CallGenesLengthScorecard.SavedBands;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import gff.GffLine;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import prok.ProkObject;
import shared.Shared;
import stream.FastaReadInputStream;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;
import structures.StringNum;

/** Regrades saved calls against all effective references and a chosen partial-excluded subset.
 * Original last-wins winners are fixed before filtering; no caller/classifier runs.
 * Exact and physical-edge-adjusted views retain query multiplicity, while recall
 * counts unique reference stop keys. Optional saved loci supply fixed length bands.
 * Instances retain path-keyed grade and contig-length caches and are not thread-safe;
 * callers must keep cached input files unchanged and configure the desired feature types.
 * The command-line entry point selects CDS-only grading and sets the shared thread limit to one.
 * @author Qiqi
 */
public final class CallGenesPartialScorecard{

	/** Parses saved-input paths and exclusion policy, selects CDS-only globals, and closes the parser stream. */
	public static void main(String[] args)throws Exception{
		final PreParser pp=new PreParser(args, CallGenesPartialScorecard.class, false);
		try{
			final HashMap<String, String> options=new HashMap<String, String>();
			for(String arg:pp.args){
				final String[] kv=Parse.splitOnFirst(arg, '=');
				require(kv.length==2 && !kv[1].isEmpty() && options.put(kv[0].toLowerCase(Locale.ROOT), kv[1])==null, "Require unique key=value arguments");
			}
			final String exclusion=options.remove("exclude");
			final String loci=options.remove("loci");
			// Historical pin arguments are accepted as metadata only; no input is hashed.
			options.remove("insha80");
			options.remove("expectedsha80");
			require(exclusion==null || exclusion.equals("any") || exclusion.equals("threeprime"),
				"exclude must be any (default) or threeprime; unknown partial sides remain eligible in threeprime mode");
			require(options.size()==3 && options.keySet().containsAll(Arrays.asList("in", "expected", "out")),
				"Require in/expected/out; optional loci binds saved length bands; legacy pins are ignored");
			Shared.setThreads(1);
			CallGenesBenchmarkSupport.cdsOnly();
			new CallGenesPartialScorecard().run(Paths.get(options.get("in")),
				Paths.get(options.get("expected")), Paths.get(options.get("out")), "threeprime".equals(exclusion), loci==null?null:Paths.get(loci));
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------    Reference partitions    ----------------*/
	/*--------------------------------------------------------------*/

	/** Count raw rows separately from effective native grader keys. The default
	 * removes every flagged partial; three-prime mode retains five-only and unknown. */
	public static final class References{
		/** Builds the default partition excluding every explicitly flagged partial. */
		public References(List<GffLine> rows){
			this(rows, false);
		}
		/** Retains original GFF objects, fixes last-wins keys, then partitions effective winners. */
		public References(List<GffLine> rows, boolean threePrimeOnly_){
			require(rows!=null, "Original reference rows are required to preserve last-wins collisions");
			threePrimeOnly=threePrimeOnly_;
			all=CallGenesStopAudit.fromLines(rows);
			final ArrayList<GffLine> kept=new ArrayList<GffLine>(), removed=new ArrayList<GffLine>();
			final HashSet<StringNum> mixed=new HashSet<StringNum>();
			for(GffLine ref:rows){
				countPartial(raw, ref);
				final GffLine winner=all.lookup(ref.seqid, ref.trueStop());
				if(ref.partial()!=winner.partial()){mixed.add(new StringNum(ref.seqid, ref.trueStop()));}
			}
			for(GffLine ref:all.records()){
				countPartial(effective, ref);
				final boolean remove=ref.partial() && (!threePrimeOnly || side(ref)==1 || side(ref)==2);
				if(remove){removed.add(ref);}else{kept.add(ref);}
			}
			complete=CallGenesStopAudit.fromLines(kept);
			excluded=CallGenesStopAudit.fromLines(removed);
			mixedPartialKeys=mixed.size();
			require(complete.uniqueStops()+excluded.uniqueStops()==all.uniqueStops() && complete.duplicateKeys==0 && excluded.duplicateKeys==0,
				"Partition only original effective winners; filtering raw rows could resurrect shadowed references");
		}
		/** complete is the retained population, not a biological completeness claim. */
		public final ReferenceSet all, complete, excluded;
		/** True retains five-prime-only and unknown-side partials. */
		public final boolean threePrimeOnly;
		/** All, complete, five-only, three-only, both and unknown counts, before/after key collapse. */
		public final long[] raw=new long[6], effective=new long[6];
		/** Raw stop keys containing a partial flag that differs from their effective winner. */
		public final int mixedPartialKeys;
	}

	/** Results preserve query occurrences for native grade/precision; recall uses
	 * unique observed reference keys, derived from the native FN-stop count. */
	public static final class Result{
		/** Retains the supplied counter arrays and verifies partition and occurrence conservation. */
		private Result(long[] all_, long[] complete_, long[] allNormalized_, long[] completeNormalized_, long[] dropped_,
				long[] droppedNormalized_, long[] startOnly_, References refs, BandCounts bands_){
			all=all_;
			complete=complete_;
			allNormalized=allNormalized_;
			completeNormalized=completeNormalized_;
			dropped=dropped_;
			droppedNormalized=droppedNormalized_;
			startOnly=startOnly_;
			bands=bands_;
			allKeys=refs.all.uniqueStops();
			completeKeys=refs.complete.uniqueStops();
			verify(all, allKeys);
			verify(complete, completeKeys);
			verify(allNormalized, allKeys);
			verify(completeNormalized, completeKeys);
			require(startOnly.length==8, "Start-only physical-edge counters must preserve exact and adjusted starts");
			final long droppedTotal=sum(dropped), droppedNormalizedTotal=sum(droppedNormalized);
			require(all[1]-complete[1]==droppedTotal && all[7]==complete[7],
				"Dropping only matched partial calls must conserve extras and remove exactly their TP occurrences");
			require(allNormalized[1]-completeNormalized[1]==droppedNormalizedTotal && allNormalized[7]==completeNormalized[7],
				"Dropping adjusted partial calls must conserve extras and remove exactly their adjusted TP occurrences");
			require(allKeys-all[3]-(completeKeys-complete[3])<=droppedTotal,
				"Removed matched reference keys cannot outnumber dropped matching query occurrences");
			require(allKeys-allNormalized[3]-(completeKeys-completeNormalized[3])<=droppedNormalizedTotal,
				"Removed adjusted reference keys cannot outnumber dropped adjusted query occurrences");
			require(!refs.threePrimeOnly || dropped[0]==0 && dropped[3]==0,
				"Terminal-stop exclusion must retain every five-prime-only and unknown-side matching call");
			require(!refs.threePrimeOnly || droppedNormalized[0]==0 && droppedNormalized[3]==0,
				"Adjusted terminal-stop exclusion must retain every five-prime-only and unknown-side matching call");
			if(bands!=null){bands.verify(this);}
		}
		/** Four eight-counter grade views followed by the two four-side dropped-occurrence arrays. */
		public final long[] all, complete, allNormalized, completeNormalized, dropped, droppedNormalized;
		/** Five-prime-only cohort size, exact/adjusted stop and start counts, errors and correction count. */
		public final long[] startOnly;
		/** Unique effective reference-stop denominators for all and retained scopes. */
		public final long allKeys, completeKeys;
		private final BandCounts bands;
	}
	/** Grades both scopes after dropping only calls matching excluded effective winners. */
	public static Result analyze(References refs, List<GffLine> calls, ContigLengths lengths){
		return analyze(refs, calls, lengths, null);
	}
	/** Grades borrowed features without mutating them; exact matches precede adjusted matches when excluding calls. */
	private static Result analyze(References refs, List<GffLine> calls, ContigLengths lengths, ReferenceBands binding){
		require(refs!=null && calls!=null && lengths!=null, "Saved reference/query records and FASTA contig lengths are required");
		final ArrayList<GffLine> kept=new ArrayList<GffLine>(calls.size()), keptNormalized=new ArrayList<GffLine>(calls.size());
		final long[] dropped=new long[4], droppedNormalized=new long[4];
		for(GffLine query:calls){
			if(!ProkObject.processType(query.prokType())){
				kept.add(query);
				keptNormalized.add(query);
				continue;
			}
			if(refs.excluded.matches(query)){
				require(!refs.complete.matches(query), "An original effective winner cannot be both included and excluded");
				final GffLine ref=refs.excluded.lookup(query.seqid, query.trueStop());
				dropped[side(ref)]++;
			}else{kept.add(query);}
			final GffLine completeExact=stopMatch(refs.complete, query, lengths, false);
			if(completeExact!=null){
				keptNormalized.add(query);
				continue;
			}
			final GffLine excludedExact=stopMatch(refs.excluded, query, lengths, false);
			if(excludedExact!=null){
				droppedNormalized[side(excludedExact)]++;
				continue;
			}
			final GffLine completeNormalized=stopMatch(refs.complete, query, lengths, true);
			if(completeNormalized!=null){
				keptNormalized.add(query);
				continue;
			}
			final GffLine excludedNormalized=stopMatch(refs.excluded, query, lengths, true);
			if(excludedNormalized==null){keptNormalized.add(query);}else{droppedNormalized[side(excludedNormalized)]++;}
		}
		final BandCounts bands=binding==null?null:new BandCounts();
		return new Result(binding==null?CallGenesStopAudit.grade(refs.all, calls):grade(refs.all, calls, lengths, false, binding, bands.counts[0]),
			binding==null?CallGenesStopAudit.grade(refs.complete, kept):grade(refs.complete, kept, lengths, false, binding, bands.counts[1]),
			grade(refs.all, calls, lengths, true, binding, bands==null?null:bands.counts[2]),
			grade(refs.complete, keptNormalized, lengths, true, binding, bands==null?null:bands.counts[3]),
			dropped, droppedNormalized, startOnly(refs, calls, lengths), refs, bands);
	}
	/** Same eight counters as the native grader, with only physical contig-edge 1-2 nt endpoints normalized. */
	public static long[] gradeNormalized(ReferenceSet refs, List<GffLine> queries, ContigLengths lengths){
		return grade(refs, queries, lengths, true, null, null);
	}
	/** Match against the full scope before assigning a frozen annotation or called
	 * band. The optional counters add no matching policy and never filter inputs. */
	private static long[] grade(ReferenceSet refs, List<GffLine> queries, ContigLengths lengths, boolean normalize,
			ReferenceBands binding, long[][] bands){
		require(refs!=null && queries!=null && lengths!=null, "Reference/query features and FASTA contig lengths are required");
		final long[] n=new long[8];
		final HashSet<StringNum> starts=new HashSet<StringNum>(), stops=new HashSet<StringNum>();
		for(GffLine q:queries){
			if(!ProkObject.processType(q.prokType())){continue;}
			final GffLine ref=stopMatch(refs, q, lengths, normalize);
			final long[] band=bands==null?null:bands[ref==null?CallGenesLengthScorecard.calledBand(q):binding.annotationBand(ref)];
			if(ref==null){
				n[6]++;
				n[7]++;
				if(band!=null){
					band[6]++;
					band[7]++;
				}
				continue;
			}
			final StringNum key=new StringNum(ref.seqid, ref.trueStop());
			stops.add(key);
			n[1]++;
			n[5]++;
			if(band!=null){
				band[1]++;
				band[5]++;
			}
			if(startMatch(ref, q, lengths, normalize)){
				starts.add(key);
				n[0]++;
				n[4]++;
				if(band!=null){
					band[0]++;
					band[4]++;
				}
			}else{
				n[6]++;
				if(band!=null){band[6]++;}
			}
		}
		if(bands!=null){
			for(GffLine ref:refs.records()){
				final long[] band=bands[binding.annotationBand(ref)];
				final StringNum key=new StringNum(ref.seqid, ref.trueStop());
				band[8]++;
				if(!starts.contains(key)){band[2]++;}
				if(!stops.contains(key)){band[3]++;}
			}
		}
		n[2]=refs.uniqueStops()-starts.size();
		n[3]=refs.uniqueStops()-stops.size();
		return n;
	}
	/** Five-prime-only partial cohort requested as the adjusted-start denominator.
	 * Start-error fields are restricted to references whose stop key matched. */
	private static long[] startOnly(References refs, List<GffLine> queries, ContigLengths lengths){
		final long[] n=new long[8];
		final HashSet<StringNum> exactStarts=new HashSet<StringNum>(), exactStops=new HashSet<StringNum>();
		final HashSet<StringNum> normalizedStarts=new HashSet<StringNum>(), normalizedStops=new HashSet<StringNum>();
		for(GffLine ref:refs.all.records()){
			if(startOnly(ref)){n[0]++;}
		}
		for(GffLine query:queries){
			if(!ProkObject.processType(query.prokType())){continue;}
			final GffLine exact=stopMatch(refs.all, query, lengths, false);
			if(startOnly(exact)){
				final StringNum key=new StringNum(exact.seqid, exact.trueStop());
				exactStops.add(key);
				if(startMatch(exact, query, lengths, false)){exactStarts.add(key);}
			}
			final GffLine normalized=stopMatch(refs.all, query, lengths, true);
			if(startOnly(normalized)){
				final StringNum key=new StringNum(normalized.seqid, normalized.trueStop());
				normalizedStops.add(key);
				if(startMatch(normalized, query, lengths, true)){normalizedStarts.add(key);}
			}
		}
		n[1]=exactStops.size();
		n[2]=exactStarts.size();
		n[3]=n[1]-n[2];
		n[4]=normalizedStops.size();
		n[5]=normalizedStarts.size();
		n[6]=n[4]-n[5];
		n[7]=n[5]-n[2];
		require(n[2]<=n[1] && n[5]<=n[4] && n[2]<=n[5], "Adjusted start-only counts cannot lose exact matches");
		return n;
	}
	/** Checks native occurrence counts against unseen-reference-key counts. */
	private static void verify(long[] grade, long keys){
		require(grade.length==8 && keys>=0, "Native grade must contain eight counters and a nonnegative key denominator");
		for(long n:grade){require(n>=0, "Native grade counters cannot be negative");}
		require(grade[0]==grade[4] && grade[1]==grade[5] && grade[2]<=keys && grade[3]<=keys && grade[1]>=keys-grade[3],
			"Native occurrence TPs, unseen-key FNs and unique matched keys must remain consistent");
	}
	/** Validated contig lengths used by edge matching; the caller binds them to the saved GFF inputs. */
	public static final class ContigLengths{
		/** Copies lengths under first-space/tab-token IDs; does not retain reads or their bases. */
		public ContigLengths(List<Read> reads){
			require(!FastaReadInputStream.SPLIT_READS && FastaReadInputStream.MIN_READ_LEN<=1,
				"FASTA length loading requires whole unsplit contigs");
			require(reads!=null && !reads.isEmpty(), "Physical-edge scorecards require nonempty FASTA contig lengths");
			for(Read r:reads){
				require(r!=null && r.id!=null && r.bases!=null && r.bases.length>0, "Every FASTA contig must have an id and positive length");
				final String id=firstToken(r.id);
				require(!id.isEmpty(), "Every FASTA contig must have a first-whitespace token id");
				require(lengths.put(id, r.bases.length)==null, "Duplicate FASTA contig first token: "+id);
			}
		}
		/** Copies positive lengths under the supplied IDs without truncating those IDs. */
		public ContigLengths(HashMap<String, Integer> lengths_){
			require(lengths_!=null && !lengths_.isEmpty(), "Physical-edge scorecards require nonempty FASTA contig lengths");
			for(Map.Entry<String, Integer> e:lengths_.entrySet()){
				require(e.getKey()!=null && e.getValue()!=null && e.getValue()>0, "Every FASTA contig length must be positive");
				require(lengths.put(e.getKey(), e.getValue())==null, "Duplicate FASTA contig id: "+e.getKey());
			}
		}
		/** Returns a known contig length in nucleotides; an absent ID is an input error. */
		int length(String seqid){
			final Integer length=lengths.get(seqid);
			require(length!=null, "Reference GFF contig is absent from the paired FASTA: "+seqid);
			return length.intValue();
		}
		private final HashMap<String, Integer> lengths=new HashMap<String, Integer>();
	}
	/** Maps the accepted biological-side classifier into four exclusive partial bins. */
	private static int side(GffLine ref){
		require(ref.partial(), "Side counts apply only to explicitly flagged partial records");
		final String side=CallGenesBenchmarkSupport.partial(ref);
		if(side.equals("five_prime")){return 0;}
		if(side.equals("three_prime")){return 1;}
		if(side.equals("both")){return 2;}
		require(side.equals("unknown"), "A flagged partial must be sided or explicitly unknown, never silently complete");
		return 3;
	}
	/** Tests whether a nonnull annotation is explicitly five-prime-only partial. */
	private static boolean startOnly(GffLine ref){return ref!=null && ref.partial() && side(ref)==0;}
	/** Counts complete flags and mutually exclusive partial sides for one annotation. */
	private static void countPartial(long[] n, GffLine ref){
		n[0]++;
		if(ref.partial()){n[side(ref)+2]++;}else{n[1]++;}
	}
	/** Returns an explicit unavailable value for a zero denominator. */
	public static String ratio(long numerator, long denominator){
		require(numerator>=0 && numerator<=denominator, "A rate numerator must be within its stated denominator");
		return denominator==0?"NA":Double.toString((double)numerator/denominator);
	}
	/** Stop equality plus the approved 0-2 nt physical three-prime edge adjustment. */
	public static boolean sameStopOrPhysicalEdge(GffLine ref, GffLine query, ContigLengths lengths){
		return startOrStopMatch(ref, query, lengths, false, true);
	}
	/** Start equality plus the approved 0-2 nt physical five-prime edge adjustment. */
	public static boolean sameStartOrPhysicalEdge(GffLine ref, GffLine query, ContigLengths lengths){
		return startOrStopMatch(ref, query, lengths, true, true);
	}
	/** Returns an exact same-type/strand winner first, otherwise a unique permitted physical-edge match. */
	private static GffLine stopMatch(ReferenceSet refs, GffLine query, ContigLengths lengths, boolean normalize){
		require(refs!=null && query!=null, "Reference set and query feature are required");
		final GffLine exact=refs.lookup(query.seqid, query.trueStop());
		if(exact!=null && exact.strand==query.strand && exact.type.equals(query.type)){return exact;}
		if(!normalize || !nearPhysicalStop(query, lengths)){return null;}
		GffLine match=null;
		for(GffLine ref:refs.records(query.seqid)){
			if(ref.strand==query.strand && ref.type.equals(query.type) && normalizedBoundary(ref, query, lengths, false)){
				require(match==null, "Ambiguous adjusted stop match");
				match=ref;
			}
		}
		return match;
	}
	/** Compares biological start coordinates after the caller has established feature identity. */
	private static boolean startMatch(GffLine ref, GffLine query, ContigLengths lengths, boolean normalize){
		require(ref!=null && query!=null, "Reference and query feature are required");
		return ref.trueStart()==query.trueStart() || normalize && normalizedBoundary(ref, query, lengths, true);
	}
	/** Requires matching contig, strand and feature type before comparing the selected endpoint. */
	private static boolean startOrStopMatch(GffLine ref, GffLine query, ContigLengths lengths, boolean start, boolean normalize){
		require(ref!=null && query!=null && lengths!=null, "Reference/query features and FASTA contig lengths are required");
		return ref.seqid.equals(query.seqid) && ref.strand==query.strand && ref.type.equals(query.type) &&
			(start ? startMatch(ref, query, lengths, normalize) : ref.trueStop()==query.trueStop() || normalize && normalizedBoundary(ref, query, lengths, false));
	}
	/** Checks phase/frame and an inward zero-to-two-base shift from an explicitly open contig boundary. */
	private static boolean normalizedBoundary(GffLine ref, GffLine query, ContigLengths lengths, boolean start){
		if(ref.seqid==null || !ref.seqid.equals(query.seqid) || ref.strand!=query.strand || !knownStrand(ref.strand) ||
			ref.phase<0 || ref.phase>2 || query.phase<0 || query.phase>2){return false;}
		final int delta=start ? startDelta(ref, query) : stopDelta(ref, query);
		return delta>=0 && delta<=2 && sameFrame(ref, query) && (start ? firstFullStart(ref)==query.trueStart() : fullCodonStop(ref, query)) &&
			(start ? startBoundary(ref, lengths) : stopBoundary(ref, lengths));
	}
	/** Accepts only the two oriented strands required by physical-boundary arithmetic. */
	private static boolean knownStrand(int strand){return strand==GffLine.PLUS || strand==GffLine.MINUS;}
	/** Returns the query start displacement inward from the reference, in nucleotides. */
	private static int startDelta(GffLine ref, GffLine query){return ref.strand==GffLine.PLUS ? query.trueStart()-ref.trueStart() : ref.trueStart()-query.trueStart();}
	/** Returns the query stop displacement inward from the reference, in nucleotides. */
	private static int stopDelta(GffLine ref, GffLine query){return ref.strand==GffLine.PLUS ? ref.trueStop()-query.trueStop() : query.trueStop()-ref.trueStop();}
	/** Tests whether the query biological stop is within two bases of the physical contig end. */
	private static boolean nearPhysicalStop(GffLine query, ContigLengths lengths){
		if(query.seqid==null || !knownStrand(query.strand)){return false;}
		final int delta=query.strand==GffLine.PLUS ? lengths.length(query.seqid)-query.trueStop() : query.trueStop()-1;
		return delta>=0 && delta<=2;
	}
	/** Requires an explicitly partial annotation open at its physical five-prime boundary. */
	private static boolean startBoundary(GffLine ref, ContigLengths lengths){return ref.partial() && (ref.strand==GffLine.PLUS ? openLeft(ref) : openRight(ref, lengths));}
	/** Requires an explicitly partial annotation open at its physical three-prime boundary. */
	private static boolean stopBoundary(GffLine ref, ContigLengths lengths){return ref.partial() && (ref.strand==GffLine.PLUS ? openRight(ref, lengths) : openLeft(ref));}
	/** Checks coordinate one and the exact open-left start_range attribute. */
	private static boolean openLeft(GffLine ref){return ref.start==1 && (".,"+ref.start).equals(CallGenesStopAudit.attribute(ref, "start_range"));}
	/** Checks the final contig coordinate and the exact open-right end_range attribute. */
	private static boolean openRight(GffLine ref, ContigLengths lengths){return ref.stop==lengths.length(ref.seqid) && (ref.stop+",.").equals(CallGenesStopAudit.attribute(ref, "end_range"));}
	/** Returns the first complete-codon coordinate after applying the annotation phase. */
	private static int firstFullStart(GffLine ref){return ref.strand==GffLine.PLUS ? ref.start+ref.phase : ref.stop-ref.phase;}
	/** Compares phase-adjusted start coordinates modulo three on the oriented strand. */
	private static boolean sameFrame(GffLine a, GffLine b){return Math.floorMod((a.strand==GffLine.PLUS?a.start+a.phase:a.stop-a.phase)-1, 3)==Math.floorMod((b.strand==GffLine.PLUS?b.start+b.phase:b.stop-b.phase)-1, 3);}
	/** Checks that the query stop ends a positive whole number of codons from the phased reference start. */
	private static boolean fullCodonStop(GffLine ref, GffLine query){
		final int codingBases=ref.strand==GffLine.PLUS ? query.trueStop()-ref.start-ref.phase+1 : ref.stop-ref.phase-query.trueStop()+1;
		return codingBases>0 && codingBases%3==0;
	}
	/** Returns the prefix before the first literal space or tab, preserving other characters. */
	private static String firstToken(String s){
		final int space=s.indexOf(' '), tab=s.indexOf('\t');
		final int cut=space<0 ? tab : (tab<0 ? space : Math.min(space, tab));
		return cut<0 ? s : s.substring(0, cut);
	}
	/** Sums the fixed partial-side counters used by conservation checks. */
	private static long sum(long[] values){
		long total=0;
		for(long n:values){total+=n;}
		return total;
	}

	/*--------------------------------------------------------------*/
	/*----------------    Saved input reduction    ----------------*/
	/*--------------------------------------------------------------*/

	/** Legacy API: pin strings are ignored; rows still bind the original per-genome grade oracle. */
	public void run(Path manifest, String manifestPin, Path expectedFile, String expectedPin, Path out)throws Exception{
		run(manifest, expectedFile, out, false);
	}
	/** Legacy API compatibility; pin strings are ignored metadata, not verified identities. */
	public void run(Path manifest, String manifestPin, Path expectedFile, String expectedPin, Path out, boolean threePrimeOnly)throws Exception{
		run(manifest, expectedFile, out, threePrimeOnly);
	}
	/** Writes both grading views with all flagged partials excluded from the retained subset. */
	public void run(Path manifest, Path expectedFile, Path out)throws Exception{run(manifest, expectedFile, out, false);}
	/** Three-prime mode measures terminal-stop eligibility by annotation, without
	 * checking literal genomic codons or treating unknown partial sides as absent stops. */
	public void run(Path manifest, Path expectedFile, Path out, boolean threePrimeOnly)throws Exception{
		run(manifest, expectedFile, out, threePrimeOnly, null);
	}
	/** Optional saved loci add length-band sidecars without changing existing tables. */
	public void run(Path manifest, Path expectedFile, Path out, boolean threePrimeOnly, Path loci)throws Exception{
		require(!Files.exists(out) && Files.isDirectory(out.toAbsolutePath().getParent()), "Output must be fresh under an existing parent");
		final String retainedScope=threePrimeOnly ? "terminal_stop_eligible_by_annotation" : "exclude_partial";
		requireInput(manifest);
		requireInput(expectedFile);
		final SavedBands savedBands;
		if(loci==null){savedBands=null;}else{
			requireInput(loci);
			final HashSet<String> genomes=new HashSet<String>();
			try(Table input=new Table(manifest, INPUT_HEADER, LEGACY_INPUT_HEADER)){
				while(input.next()){
					final String genome=input.s(1);
					require(genome.matches("g[0-9]{4}"), "Safe saved-band genome key required");
					genomes.add(genome);
				}
			}
			savedBands=CallGenesLengthScorecard.readSavedBands(loci, genomes);
		}
		final LinkedHashMap<String, long[]> expected=expected(expectedFile);
		final LinkedHashMap<String, Totals> totals=new LinkedHashMap<String, Totals>();
		for(String policy:expected.keySet()){totals.put(policy, new Totals());}
		final HashSet<String> identities=new HashSet<String>();
		final HashMap<String, String> referencePaths=new HashMap<String, String>(), fastaPaths=new HashMap<String, String>();
		final HashMap<String, long[]> censusByGenome=new HashMap<String, long[]>();
		final long[] rawTotal=new long[6], effectiveTotal=new long[6];
		long duplicateKeys=0, mixedKeys=0;
		final ByteBuilder rows=new ByteBuilder(SCORE_HEADER), drops=new ByteBuilder(DROP_HEADER), normalizedDrops=new ByteBuilder(DROP_HEADER);
		final ByteBuilder startRows=new ByteBuilder(START_ONLY_HEADER), census=new ByteBuilder(CENSUS_HEADER);
		final ByteBuilder bandRows=savedBands==null?null:new ByteBuilder(BAND_SCORE_HEADER);
		try(Table input=new Table(manifest, INPUT_HEADER, LEGACY_INPUT_HEADER)){
			final boolean legacy=input.width==11;
			while(input.next()){
				final String policy=input.s(0), genome=input.s(1), arm=input.s(legacy ? 10 : 6);
				require(expected.containsKey(policy) && genome.matches("g[0-9]{4}") && identities.add(policy+"/"+genome), "Unique explicitly expected policy/genome rows required");
				require(arm.equals("native") || arm.equals("bonus400"), "Only retained native/bonus400 grades are allowed");
				final String fastaName=input.s(2), referenceName=input.s(legacy ? 4 : 3);
				final Path fasta=Paths.get(fastaName), reference=Paths.get(referenceName),
					calls=Paths.get(input.s(legacy ? 6 : 4)), grades=Paths.get(input.s(legacy ? 8 : 5));
				requireInput(fasta);
				requireInput(reference);
				requireInput(calls);
				requireInput(grades);
				final String previousFasta=fastaPaths.putIfAbsent(genome, fastaName);
				require(previousFasta==null || previousFasta.equals(fastaName), "All policies must use the same literal FASTA path per genome");
				final String previous=referencePaths.putIfAbsent(genome, referenceName);
				require(previous==null || previous.equals(referenceName), "All policies must use the same literal reference path per genome");
				final References refs=new References(CallGenesStopAudit.loadReference(reference.toString()).sourceRecords(), threePrimeOnly);
				final ReferenceBands binding=savedBands==null?null:savedBands.bind(genome, refs.all);
				final Result result=analyze(refs, GffLine.loadGffFile(calls.toString(), null, false), contigLengths(fasta), binding);
				final long[] saved=grades(grades).get(genome+"/"+arm);
				require(saved!=null && Arrays.equals(saved, result.all), "All eight all-gene counts must reproduce saved grades: "+policy+"/"+genome);
				final long[] signature=new long[14];
				System.arraycopy(refs.raw, 0, signature, 0, 6);
				System.arraycopy(refs.effective, 0, signature, 6, 6);
				signature[12]=refs.all.duplicateKeys;
				signature[13]=refs.mixedPartialKeys;
				final long[] prior=censusByGenome.putIfAbsent(genome, signature);
				require(prior==null || Arrays.equals(prior, signature), "Reference partial/collision census must be identical across policies");
				if(prior==null){
					census(census, genome, "source_rows", refs.raw, refs.all.duplicateKeys, refs.mixedPartialKeys);
					census(census, genome, "effective_keys", refs.effective, refs.all.duplicateKeys, refs.mixedPartialKeys);
					add(rawTotal, refs.raw);
					add(effectiveTotal, refs.effective);
					duplicateKeys+=refs.all.duplicateKeys;
					mixedKeys+=refs.mixedPartialKeys;
				}
				score(rows, policy, genome, "all_genes", result.allKeys, result.all, result.allNormalized);
				score(rows, policy, genome, retainedScope, result.completeKeys, result.complete, result.completeNormalized);
				if(result.bands!=null){scoreBands(bandRows, policy, genome, retainedScope, result.bands);}
				dropped(drops, policy, genome, result.dropped);
				dropped(normalizedDrops, policy, genome, result.droppedNormalized);
				startOnly(startRows, policy, genome, result.startOnly);
				final Totals total=totals.get(policy);
				total.genomes++;
				total.allKeys+=result.allKeys;
				total.completeKeys+=result.completeKeys;
				add(total.all, result.all);
				add(total.complete, result.complete);
				add(total.allNormalized, result.allNormalized);
				add(total.completeNormalized, result.completeNormalized);
				add(total.dropped, result.dropped);
				add(total.droppedNormalized, result.droppedNormalized);
				add(total.startOnly, result.startOnly);
				if(result.bands!=null){
					if(total.bands==null){total.bands=new BandCounts();}
					total.bands.add(result.bands);
				}
			}
		}
		final ByteBuilder summary=new ByteBuilder(SCORE_HEADER);
		final ByteBuilder bandSummary=savedBands==null?null:new ByteBuilder(BAND_SCORE_HEADER);
		for(Map.Entry<String, long[]> e:expected.entrySet()){
			final String policy=e.getKey();
			final long[] want=e.getValue();
			final Totals t=totals.get(policy);
			require(t.genomes==want[0] && t.genomes==referencePaths.size(), "Every policy must cover the same complete reference roster");
			for(String genome:referencePaths.keySet()){require(identities.contains(policy+"/"+genome), "Missing policy/genome input");}
			final long[] measured={t.genomes, t.allKeys, t.all[1], t.all[3], t.all[7], t.all[3]+t.all[7], t.all[2]+t.all[6]};
			require(Arrays.equals(want, measured), "All-gene aggregate controls disagree for "+policy+": "+Arrays.toString(measured));
			score(summary, policy, "ALL", "all_genes", t.allKeys, t.all, t.allNormalized);
			score(summary, policy, "ALL", retainedScope, t.completeKeys, t.complete, t.completeNormalized);
			if(t.bands!=null){
				t.bands.verify(t.all, t.complete, t.allNormalized, t.completeNormalized, t.allKeys, t.completeKeys);
				scoreBands(bandSummary, policy, "ALL", retainedScope, t.bands);
			}
			dropped(drops, policy, "ALL", t.dropped);
			dropped(normalizedDrops, policy, "ALL", t.droppedNormalized);
			startOnly(startRows, policy, "ALL", t.startOnly);
		}
		census(census, "ALL", "source_rows", rawTotal, duplicateKeys, mixedKeys);
		census(census, "ALL", "effective_keys", effectiveTotal, duplicateKeys, mixedKeys);
		Files.createDirectory(out);
		write(out.resolve("scorecard.tsv"), summary);
		write(out.resolve("by_genome.tsv"), rows);
		write(out.resolve("partial_counts.tsv"), census);
		write(out.resolve("dropped_calls.tsv"), drops);
		write(out.resolve("normalized_dropped_calls.tsv"), normalizedDrops);
		write(out.resolve("start_only_partial_starts.tsv"), startRows);
		if(savedBands!=null){
			write(out.resolve("band_scorecard.tsv"), bandSummary);
			write(out.resolve("band_by_genome.tsv"), bandRows);
		}
		final String pass="PARTIAL_SCORECARD_PASS policies="+expected.size()+" genomes_per_policy="+referencePaths.size()+" grade_fields=8 caller_runs=0 classifier_runs=0"+
			(threePrimeOnly ? " exclusion=threeprime" : "")+"\n";
		write(out.resolve("PASS"), new ByteBuilder(pass));
		System.out.print(pass);
	}
	/*--------------------------------------------------------------*/
	/*----------------    Output tables    ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes native counts with unique-key recall and query-occurrence precision. */
	private static void score(ByteBuilder b, String policy, String genome, String scope, long keys, long[] g){
		score(b, policy, genome, scope, keys, g, g);
	}
	/** Writes native exact counts followed by the symmetric physical-edge adjusted view. */
	private static void score(ByteBuilder b, String policy, String genome, String scope, long keys, long[] g, long[] adjusted){
		score(b, policy, genome, scope, null, keys, g, adjusted);
	}
	/** Appends one checked score row, deriving rates from counts and emitting the eight raw grade fields. */
	private static void score(ByteBuilder b, String policy, String genome, String scope, String band, long keys, long[] g, long[] adjusted){
		verify(g, keys);
		verify(adjusted, keys);
		final long matched=keys-g[3], queries=g[5]+g[7], matchedNormalized=keys-adjusted[3], queriesNormalized=adjusted[5]+adjusted[7];
		b.append(policy).tab().append(genome).tab().append(scope);
		if(band!=null){b.tab().append(band);}
		b.tab().append(keys).tab().append(g[1]).tab().append(g[3]).tab().append(g[7])
			.tab().append(g[3]+g[7]).tab().append(g[2]+g[6]).tab().append(matched).tab().append(keys).tab().append(ratio(matched, keys))
			.tab().append(g[5]).tab().append(queries).tab().append(ratio(g[5], queries)).tab().append(g[1]-matched)
			.tab().append(adjusted[1]).tab().append(adjusted[3]).tab().append(adjusted[7]).tab().append(adjusted[3]+adjusted[7]).tab().append(adjusted[2]+adjusted[6])
			.tab().append(matchedNormalized).tab().append(keys).tab().append(ratio(matchedNormalized, keys))
			.tab().append(adjusted[5]).tab().append(queriesNormalized).tab().append(ratio(adjusted[5], queriesNormalized)).tab().append(adjusted[1]-matchedNormalized);
		for(long n:g){b.tab().append(n);}
		b.nl();
	}
	/** ALL is the sum of seven disjoint bins; 40_99 is exactly bins 1 and 2. */
	private static void scoreBands(ByteBuilder out, String policy, String genome, String retainedScope, BandCounts bands){
		for(int scope=0; scope<2; scope++){
			for(int b=0; b<9; b++){
				final long[] raw=bands.view(scope, b), normalized=bands.view(scope+2, b);
				require(raw[8]==normalized[8], "Normalization cannot change frozen annotation-band denominators");
				score(out, policy, genome, scope==0?"all_genes":retainedScope, BAND_NAMES[b], raw[8],
					Arrays.copyOf(raw, 8), Arrays.copyOf(normalized, 8));
			}
		}
	}
	/** Writes exclusive partial bins and overlapping inclusive end counts. */
	private static void census(ByteBuilder b, String genome, String scope, long[] n, long duplicate, long mixed){
		final long any=n[2]+n[3]+n[4]+n[5];
		require(n[0]==n[1]+any, "Exclusive partial categories must conserve reference population");
		b.append(genome).tab().append(scope).tab().append(n[0]).tab().append(n[1]).tab().append(any);
		for(int i=2; i<6; i++){b.tab().append(n[i]);}
		b.tab().append(n[2]+n[4]).tab().append(n[3]+n[4]).tab().append(duplicate).tab().append(mixed).nl();
	}
	/** Writes dropped query occurrences, retaining duplicate calls. */
	private static void dropped(ByteBuilder b, String policy, String genome, long[] n){
		require(n.length==4, "Drop categories must preserve five-only,three-only,both,unknown");
		b.append(policy).tab().append(genome).tab().append(n[0]+n[1]+n[2]+n[3]);
		for(long x:n){b.tab().append(x);}
		b.tab().append(n[0]+n[2]).tab().append(n[1]+n[2]).nl();
	}
	/** Writes exact and adjusted starts for five-prime-only partials with a matched stop. */
	private static void startOnly(ByteBuilder b, String policy, String genome, long[] n){
		require(n.length==8 && n[2]<=n[1] && n[5]<=n[4] && n[7]==n[5]-n[2],
			"Start-only physical-edge counters must preserve their exact/adjusted shape");
		b.append(policy).tab().append(genome);
		for(long x:n){b.tab().append(x);}
		b.nl();
	}
	/** Adds identically shaped counters with overflow detection. */
	private static void add(long[] sum, long[] values){
		require(sum.length==values.length, "Counter shapes must agree");
		for(int i=0; i<sum.length; i++){sum[i]=Math.addExact(sum[i], values[i]);}
	}
	/** Reads independent all-gene aggregate controls for each named policy. */
	private static LinkedHashMap<String, long[]> expected(Path path){
		final LinkedHashMap<String, long[]> result=new LinkedHashMap<String, long[]>();
		try(Table input=new Table(path, EXPECTED_HEADER)){
			while(input.next()){
				final String policy=input.s(0);
				require(policy.matches("[A-Za-z][A-Za-z0-9_]*"), "Safe named policy required");
				final long[] values=new long[7];
				for(int i=0; i<7; i++){
					values[i]=input.n(i+1);
					require(values[i]>=0, "Expected controls cannot be negative");
				}
				require(values[0]>0 && result.put(policy, values)==null, "Unique policy and positive genome count required");
			}
		}
		require(!result.isEmpty(), "Expected policies cannot be empty");
		return result;
	}
	/** Caches saved per-genome eight-field grade oracles without changing multiplicity. */
	private HashMap<String, long[]> grades(Path path){
		HashMap<String, long[]> result=gradeCache.get(path);
		if(result!=null){return result;}
		result=new HashMap<String, long[]>();
		try(Table input=new Table(path, GRADE_HEADER)){
			while(input.next()){
				final String genome=input.s(0), arm=input.s(1);
				require(genome.matches("g[0-9]{4}") && (arm.equals("native") || arm.equals("bonus400") || arm.equals("zero_effect")),
					"Unexpected saved grade identity");
				final long[] values=new long[8];
				for(int i=0; i<8; i++){
					values[i]=input.n(i+2);
					require(values[i]>=0, "Negative saved grade");
				}
				require(result.put(genome+"/"+arm, values)==null, "Duplicate saved grade identity");
			}
		}
		gradeCache.put(path, result);
		return result;
	}
	/** Loads listed FASTA lengths once per normalized path for physical-boundary matching. */
	private ContigLengths contigLengths(Path path){
		final Path key=path.toAbsolutePath().normalize();
		ContigLengths lengths=lengthCache.get(key);
		if(lengths==null){
			lengths=new ContigLengths(ReadInputStream.toReads(key.toString(), FileFormat.FASTA, -1));
			lengthCache.put(key, lengths);
		}
		return lengths;
	}
	/** Path metadata only; the committed manifest is the input identity contract. */
	private static void requireInput(Path path){
		require(path.isAbsolute() && Files.isRegularFile(path) && !Files.isSymbolicLink(path), "Explicit regular input paths required: "+path);
	}
	/** Writes one complete table and propagates asynchronous writer failure. */
	private static void write(Path path, ByteBuilder data){
		final ByteStreamWriter out=new ByteStreamWriter(path.toString(), false, false, false);
		out.start();
		out.print(data);
		require(!out.poisonAndWait(), "Scorecard output write failed: "+path);
	}
	/*--------------------------------------------------------------*/
	/*----------------       Table and totals       ----------------*/
	/*--------------------------------------------------------------*/

	/** Strict native byte-table reader with exact header and field-count checks. */
	private static final class Table implements AutoCloseable{
		/** Opens the table and closes it if header validation fails. */
		Table(Path path, String... headers){
			in=ByteFile.makeByteFile1(path.toString(), false);
			try{
				final byte[] actual=in.nextLine();
				String matched=null;
				for(String header:headers){
					if(Arrays.equals(actual, header.getBytes(StandardCharsets.US_ASCII))){
						matched=header;
						break;
					}
				}
				require(matched!=null, "Unexpected table schema: "+path);
				width=matched.split("\t", -1).length;
			}catch(RuntimeException | Error e){
				in.close();
				throw e;
			}
		}
		/** Advances to one complete row; malformed widths cannot be silently parsed. */
		boolean next(){
			final byte[] line=in.nextLine();
			if(line==null){return false;}
			parser.set(line);
			require(parser.terms()==width, "Table row width changed");
			return true;
		}
		/** Materializes one text field from the current borrowed input line. */
		String s(int i){return parser.parseString(i);}
		/** Parses one decimal counter from the current input line. */
		long n(int i){return parser.parseLong(i);}
		/** Closes the byte reader and propagates its error flag. */
		@Override
		public void close(){require(!in.close(), "Table read failed");}
		final ByteFile in;
		final int width;
		final LineParser1 parser=new LineParser1('\t');
	}
	/** Per-policy native counts and dropped-call categories. */
	private static final class Totals{
		long genomes, allKeys, completeKeys;
		final long[] all=new long[8], complete=new long[8], allNormalized=new long[8], completeNormalized=new long[8];
		final long[] dropped=new long[4], droppedNormalized=new long[4], startOnly=new long[8];
		BandCounts bands;
	}
	/** Four existing grading views, seven disjoint bands, eight counters plus keys. */
	private static final class BandCounts{
		final long[][][] counts=new long[4][7][9];
		/** Adds all four views without changing the seven disjoint band assignments. */
		void add(BandCounts other){
			for(int scope=0; scope<4; scope++){
				for(int band=0; band<7; band++){CallGenesPartialScorecard.add(counts[scope][band], other.counts[scope][band]);}
			}
		}
		/** Borrows a disjoint row for bands 0–6; returns a new aggregate for ALL or 40_99. */
		long[] view(int scope, int band){
			if(band<7){return counts[scope][band];}
			final long[] out=new long[9];
			for(int b=band==7?0:1; b<(band==7?7:3); b++){CallGenesPartialScorecard.add(out, counts[scope][b]);}
			return out;
		}
		/** Checks band conservation against one genome result. */
		void verify(Result r){verify(r.all, r.complete, r.allNormalized, r.completeNormalized, r.allKeys, r.completeKeys);}
		/** Requires every view to sum to its eight aggregate counters and reference-key denominator. */
		void verify(long[] all, long[] complete, long[] normalized, long[] retainedNormalized, long allKeys, long retainedKeys){
			final long[][] expected={all, complete, normalized, retainedNormalized};
			for(int scope=0; scope<4; scope++){
				final long[] actual=view(scope, 7);
				require(Arrays.equals(Arrays.copyOf(actual, 8), expected[scope]) && actual[8]==(scope%2==0?allKeys:retainedKeys),
					"Disjoint annotation/called bands must conserve every existing ALL counter and reference key");
			}
		}
	}
	/** Rejects malformed inputs or violated scorecard invariants even when Java assertions are disabled. */
	private static void require(boolean valid, String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}
	/*--------------------------------------------------------------*/
	/*----------------    Fields    ----------------*/
	/*--------------------------------------------------------------*/

	/** Saved grade arrays cached by the supplied path; inputs must remain unchanged. */
	private final HashMap<Path, HashMap<String, long[]>> gradeCache=new HashMap<Path, HashMap<String, long[]>>();
	/** Contig lengths cached by absolute normalized FASTA path. */
	private final HashMap<Path, ContigLengths> lengthCache=new HashMap<Path, ContigLengths>();
	/*--------------------------------------------------------------*/
	/*----------------    Constants    ----------------*/
	/*--------------------------------------------------------------*/

	public static final String INPUT_HEADER="policy\tgenome\tfasta\treference\tcalls\tgrades\tgrade_arm";
	/** Compatibility schema only: legacy SHA fields are ignored and never verified. */
	public static final String LEGACY_INPUT_HEADER="policy\tgenome\tfasta\tfasta_sha80\treference\treference_sha80\tcalls\tcalls_sha80\tgrades\tgrades_sha80\tgrade_arm";
	public static final String EXPECTED_HEADER="policy\tgenomes\treference_keys\tcorrect\tmissed\textras\tstop_errors\tstart_errors";
	private static final String GRADE_HEADER="genome\tarm\t"+String.join("\t", CallGenesBenchmarkSupport.GRADE_NAMES);
	private static final String SCORE_HEADER="policy\tgenome\tscope\treference_keys\tcorrect\tmissed\textras\tstop_errors\tstart_errors\trecall_numerator\trecall_denominator\trecall\tprecision_numerator\tprecision_denominator\tprecision\tduplicate_matching_occurrences\tnormalized_correct\tnormalized_missed\tnormalized_extras\tnormalized_stop_errors\tnormalized_start_errors\tnormalized_recall_numerator\tnormalized_recall_denominator\tnormalized_recall\tnormalized_precision_numerator\tnormalized_precision_denominator\tnormalized_precision\tnormalized_duplicate_matching_occurrences\t"+String.join("\t", CallGenesBenchmarkSupport.GRADE_NAMES)+"\n";
	private static final String BAND_SCORE_HEADER=SCORE_HEADER.replace("scope\t", "scope\tannotation_or_called_length_band\t");
	private static final String[] BAND_NAMES={"lt40", "40_64", "65_99", "100_149", "150_199", "ge200", "unknown", "ALL", "40_99"};
	private static final String DROP_HEADER="policy\tgenome\tdropped_occurrences\tfive_prime_only\tthree_prime_only\tboth\tunknown_side\tfive_prime_inclusive\tthree_prime_inclusive\n";
	private static final String START_ONLY_HEADER="policy\tgenome\tfive_prime_references\texact_stop_keys\texact_start_keys\texact_matched_stop_start_errors\tnormalized_stop_keys\tnormalized_start_keys\tnormalized_matched_stop_start_errors\tnormalized_start_corrections\n";
	private static final String CENSUS_HEADER="genome\tscope\tall\tcomplete_flag\tany_partial\tfive_prime_only\tthree_prime_only\tboth\tunknown_side\tfive_prime_inclusive\tthree_prime_inclusive\tduplicate_raw_keys\tmixed_partial_raw_keys\n";
}
