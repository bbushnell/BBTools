package test.processor;

import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashSet;

import parse.Parser;
import processor.ReformatProcessor;
import stream.Read;
import tracker.ReadStats;

/** Literal name, output-bit and clone ownership controls. @author Sucrose */
public final class ReformatProcessorBoundaryTest{

	/** Optional scopes isolate non-name, Unicode conversion or force-trim statistics controls. */
	public static void main(final String[] args){
		if(args.length<1 || args.length>2 || (!args[0].equals("phase=baseline") && !args[0].equals("phase=fixed")) ||
			(args.length==2 && !args[1].equals("scope=nonname") && !args[1].equals("scope=all") &&
				!args[1].equals("scope=unicodeheader") && !args[1].equals("scope=forcetrimstats"))){
			throw new IllegalArgumentException("Expected phase=baseline|fixed [scope=all|nonname|unicodeheader|forcetrimstats]");
		}
		final boolean baseline=args[0].equals("phase=baseline");
		final boolean nonname=args.length==2 && args[1].equals("scope=nonname");
		final boolean unicode=args.length==2 && args[1].equals("scope=unicodeheader");
		final boolean forcetrim=args.length==2 && args[1].equals("scope=forcetrimstats");
		if(forcetrim){
			forceTrimStats(baseline);
		}else if(unicode){
			unicodeHeaders(baseline);
		}else{
			if(!nonname){
				singleNames(baseline);
				pairedNames(baseline);
				sharedNames(baseline);
				collisionNames(baseline);
			}
			outputBits();
			cloneState();
			emptyStats(baseline);
			minimumSubs(baseline);
		}
		System.out.println("REFORMAT_PROCESSOR_BOUNDARY_PASS phase="+(baseline ? "baseline" : "fixed")+" checks="+checks);
	}

	/** First names remain unchanged; only repeats acquire the existing second-occurrence suffix. */
	private static void singleNames(final boolean baseline){
		final ReformatProcessor proc=processor(true);
		final String[] input={"x", "x", "y"};
		final String[] expected=baseline ? new String[]{"x_0", "x_1", "y_0"} : new String[]{"x", "x_2", "y"};
		for(int i=0; i<input.length; i++){
			final Read r=read(input[i], "ACGT", i);
			check(proc.processReadPair(r, null)==1, "A valid single read must retain only bit one.");
			check(expected[i].equals(r.id), "Single name "+i+": expected "+expected[i]+", observed "+r.id);
		}
		check(proc.readsProcessedT==3 && proc.basesProcessedT==12 && proc.singlesOut==3 && proc.singleBasesOut==12,
			"Three four-base single reads must be counted once each.");
		System.out.println("unique_single\t"+expected[0]+","+expected[1]+","+expected[2]);
	}

	/** Mate histories remain independent; the fixed policy inserts suffixes before terminal mate markers. */
	private static void pairedNames(final boolean baseline){
		final ReformatProcessor proc=processor(true);
		final String[] input={"p", "p", "q"};
		final String[] suffix=baseline ? new String[]{"_0", "_1", "_0"} : new String[]{"", "_2", ""};
		for(int i=0; i<input.length; i++){
			final Read left=read(input[i]+"/1", "ACGT", i), right=read(input[i]+"/2", "TGCA", i);
			pair(left, right);
			check(proc.processReadPair(left, right)==3, "A valid pair must retain both bits.");
			final String expectedLeft=baseline ? input[i]+"/1"+suffix[i] : input[i]+suffix[i]+"/1";
			final String expectedRight=baseline ? input[i]+"/2"+suffix[i] : input[i]+suffix[i]+"/2";
			check(expectedLeft.equals(left.id), "Unexpected first mate name: "+left.id);
			check(expectedRight.equals(right.id), "Unexpected second mate name: "+right.id);
			check(left.mate==right && right.mate==left, "The processor must leave mate routing to its driver.");
		}
		check(proc.readsProcessedT==6 && proc.basesProcessedT==24 && proc.pairsOut==3 && proc.pairBasesOut==24,
			"Three pairs contain six four-base reads and twenty-four bases.");
	}

	/** Cloning must not split the serial name history into competing namespaces. */
	private static void sharedNames(final boolean baseline){
		final ReformatProcessor parent=processor(true), clone=parent.clone();
		final Read first=read("shared", "ACGT", 0), second=read("shared", "ACGT", 1);
		check(clone.processReadPair(first, null)==1 && parent.processReadPair(second, null)==1,
			"Serial parent/clone calls must retain both valid reads.");
		check((baseline ? "shared_0" : "shared").equals(first.id), "First shared-map name: "+first.id);
		check((baseline ? "shared_1" : "shared_2").equals(second.id), "Second shared-map name: "+second.id);
		check(parent.recommendedWorkers()==1 && clone.recommendedWorkers()==1,
			"Shared name maps require the existing single-worker recommendation.");
	}

	/** Original names and earlier generated suffixes must not collide in either mate namespace. */
	private static void collisionNames(final boolean baseline){
		final String[][] input={
			{"x", "x", "x_2"},
			{"x_2", "x", "x"},
			{"x_3", "x", "x", "x", "x_2"},
			{"x", "x", "x_02", "x_1", "x_2_2", "x_2"}
		};
		final String[][] before={
			{"x_0", "x_1", "x_2_0"},
			{"x_2_0", "x_0", "x_1"},
			{"x_3_0", "x_0", "x_1", "x_2", "x_2_0"},
			{"x_0", "x_1", "x_02_0", "x_1_0", "x_2_2_0", "x_2_0"}
		};
		final String[][] after={
			{"x", "x_2", "x_2_2"},
			{"x_2", "x", "x_3"},
			{"x_3", "x", "x_2", "x_4", "x_2_2"},
			{"x", "x_2", "x_02", "x_1", "x_2_2", "x_2_3"}
		};
		final String[][] expected=baseline ? before : after;
		for(int row=0; row<input.length; row++){
			final ReformatProcessor proc=processor(true);
			for(int i=0; i<input[row].length; i++){
				final Read r=read(input[row][i], "ACGT", i);
				check(proc.processReadPair(r, null)==1, "Name collisions must not discard valid reads.");
				check(expected[row][i].equals(r.id), "Collision row "+row+", item "+i+": expected "+expected[row][i]+", observed "+r.id);
			}
		}
		final ReformatProcessor paired=processor(true);
		for(int i=0; i<3; i++){
			final Read left=read(input[0][i], "ACGT", i), right=read(input[1][i], "TGCA", i);
			pair(left, right);
			check(paired.processReadPair(left, right)==3, "Independent mate name collisions must retain both reads.");
			check(expected[0][i].equals(left.id), "First-mate collision namespace: "+left.id);
			check(expected[1][i].equals(right.id), "Second-mate collision namespace: "+right.id);
		}
	}

	/** The processor returns read-selection bits; the driver owns list/singleton routing. */
	private static void outputBits(){
		for(int mask=0; mask<4; mask++){
			final ReformatProcessor proc=processor(false);
			proc.minReadLength=3;
			proc.requireBothBad=true;
			final Read left=read("left", (mask&1)==0 ? "AC" : "ACGT", mask);
			final Read right=read("right", (mask&2)==0 ? "TG" : "TGCA", mask);
			pair(left, right);
			check(proc.processReadPair(left, right)==mask, "Length controls require literal retained mask "+mask);
		}
		final ReformatProcessor proc=processor(false);
		proc.minReadLength=3;
		proc.requireBothBad=false;
		final Read left=read("left", "ACGT", 0), right=read("right", "TG", 0);
		pair(left, right);
		check(proc.processReadPair(left, right)==0, "Either-bad policy must remove a mixed passing/failing pair.");
	}

	/** Fresh worker scratch coexists with copied scalar totals and borrowed frozen configuration. */
	private static void cloneState(){
		final ReformatProcessor parent=processor(false);
		parent.readsProcessedT=7;
		parent.basesProcessedT=11;
		parent.remap1=new byte[]{0, 1};
		parent.remap2=parent.remap1;
		parent.barcodes=new HashSet<String>();
		parent.barcodes.add("ACGT");
		parent.readstats=new ReadStats(false);
		final ArrayList<ReadStats> previous=ReadStats.objectList;
		try{
			ReadStats.objectList=new ArrayList<ReadStats>();
			final ReformatProcessor clone=parent.clone();
			check(clone!=parent && clone.rhp!=parent.rhp, "A worker needs a distinct processor and header parser.");
			check(clone.readstats!=parent.readstats && ReadStats.objectList.size()==1 && ReadStats.objectList.get(0)==clone.readstats,
				"A cloned stats object must be independent and registered for ReadStats aggregation.");
			check(clone.remap1==parent.remap1 && clone.remap2==parent.remap2 && clone.barcodes==parent.barcodes,
				"Configured remaps/barcodes retain their established borrowed-reference semantics.");
			check(clone.readsProcessedT==7 && clone.basesProcessedT==11,
				"clone copies existing scalar counters; it is not a counter-reset operation.");
			parent.add(clone);
			check(parent.readsProcessedT==14 && parent.basesProcessedT==22,
				"add sums copied scalar totals; drivers must clone before processing and merge each worker once.");
		}finally{ReadStats.objectList=previous;}
	}

	/** Empty denominators have explicit baseline defects and literal finite fixed percentages. */
	private static void emptyStats(final boolean baseline){
		final String percent=baseline ? "NaN" : "0.00";
		final ReformatProcessor empty=processor(false);
		empty.trimBadSequence=true;
		final String expected="QTrimmed:               \t0 reads ("+percent+"%) \t0 bases ("+percent+"%)\n";
		check(expected.equals(empty.toStats().toString()), "Empty trimming statistics must match literal percentages: "+empty.toStats());
		final ReformatProcessor shortRead=processor(false);
		shortRead.k=5;
		check(shortRead.processReadPair(read("short", "ACGT", 0), null)==1,
			"A read shorter than k remains a valid output read.");
		check(shortRead.kmersProcessed==0 && shortRead.correctKmers==0,
			"A four-base read contains no five-mers.");
		final String kStats="5-mers processed:        \t0\nCorrect 5-mers:          \t"+percent+"%\n";
		check(kStats.equals(shortRead.toStats().toString()), "Zero-kmer statistics must match literal percentages: "+shortRead.toStats());
	}

	/** A perfect alignment has zero substitutions and must obey a positive minimum too. */
	private static void minimumSubs(final boolean baseline){
		final ReformatProcessor proc=processor(false);
		final Read perfect=read("perfect", "ACGT", 0);
		perfect.match="mmmm".getBytes(StandardCharsets.US_ASCII);
		perfect.setMapped(true);
		perfect.setPerfect(true);
		proc.minsubs=1;
		check(proc.passesEditFilter(perfect, false)==baseline,
			"A stored-perfect read has zero substitutions; baseline accepts it, the repair must reject minsubs=1.");
		proc.minsubs=0;
		check(proc.passesEditFilter(perfect, false), "Zero minimum substitutions must still accept a perfect alignment.");
		proc.minsubs=1;
		perfect.setPerfect(false);
		check(!proc.passesEditFilter(perfect, false), "A zero-substitution match without the perfect bit must fail the same minimum.");
		perfect.match="mSmm".getBytes(StandardCharsets.US_ASCII);
		check(proc.passesEditFilter(perfect, false), "One literal substitution meets minsubs=1.");
	}

	/** JGI-021: UTF-16 offsets must not index UTF-8 bytes; ASCII remains a nearby control. */
	private static void unicodeHeaders(final boolean baseline){
		final Charset charset=Charset.defaultCharset();
		System.out.println("header_encoding\t"+charset.name()+"\tuniquenames=false");
		if(!charset.equals(StandardCharsets.UTF_8)){
			throw new IllegalStateException("Unicode baseline requires UTF-8 default encoding, observed "+charset);
		}
		final String[] input={"\u00e9_a_b", "x_a_b"};
		final String[] option={"underscoretospace", "lastunderscoretospace"};
		final String[] expected={baseline ? "\u00e9_ _b" : "\u00e9_a b", "x_a b"};
		for(int i=0; i<input.length; i++){
			final ReformatProcessor proc=new ReformatProcessor();
			proc.setFromParser(new Parser());
			proc.minReadLength=0;
			proc.maxReadLength=0;
			if(!proc.parse("uniquenames=f", "uniquenames", "f") || !proc.parse(option[i]+"=t", option[i], "t")){
				throw new IllegalStateException("The real processor parser must recognize header options.");
			}
			proc.postParse();
			if(proc.uniqueNames || !proc.lastUnderscoreToSpace){throw new IllegalStateException("Enable only underscore conversion.");}
			final Read read=new Read("ACGT".getBytes(StandardCharsets.US_ASCII), new byte[]{40, 40, 40, 40}, input[i], i);
			check(proc.processReadPair(read, null)==1, "Header conversion must retain the valid four-base single read.");
			check(expected[i].equals(read.id), "Header conversion expected "+expected[i]+", observed "+read.id);
			System.out.println("header_result\t"+input[i]+"\t"+read.id);
		}
	}

	/** JGI-022: force-trim failures and ordinary short reads must share the short-discard report.
	 * Exact percentages retain input-base denominators and count remaining discarded bases.
	 * Only short-after-trim cases differ by phase; output bits and sequences must not change.
	 */
	private static void forceTrimStats(final boolean baseline){
		final String forcedSix="FTrimmed:               \t1 reads (100.00%) \t6 bases (75.00%)\n";
		final String shortTwo="Short Read Discards:    \t1 reads (100.00%) \t2 bases (25.00%)\n";
		forceTrimStatsCase("short after trim", "ACGTACGT", null, 6, 0, "GT", null,
			forcedSix+(baseline ? "" : shortTwo));
		forceTrimStatsCase("untrimmed short", "GT", null, 0, 0, "GT", null,
			"Short Read Discards:    \t1 reads (100.00%) \t2 bases (100.00%)\n");
		forceTrimStatsCase("passing trim", "ACGTACGT", null, 4, 1, "ACGT", null,
			"FTrimmed:               \t1 reads (100.00%) \t4 bases (50.00%)\n");
		final String pairTrim="FTrimmed:               \t2 reads (100.00%) \t12 bases (66.67%)\n";
		final String pairShort="Short Read Discards:    \t1 reads (50.00%) \t2 bases (11.11%)\n";
		forceTrimStatsCase("second mate short after trim", "ACGTACGTAC", "ACGTACGT", 6, 1,
			"GTAC", "GT", pairTrim+(baseline ? "" : pairShort));
	}

	/** Creates an independent processor per case; the driver remains responsible for mate routing. */
	private static void forceTrimStatsCase(final String label, final String bases1, final String bases2,
		final int left, final int expectedKeep, final String expected1, final String expected2,
		final String expectedStats){
		final ReformatProcessor proc=processor(false);
		proc.minReadLength=3;
		proc.forceTrimLeft=left;
		proc.forceTrimRight=-1;
		proc.forceTrimRight2=0;
		proc.forceTrimModulo=0;
		proc.requireBothBad=true;
		final Read first=read(label+"/1", bases1, 0);
		final Read second=bases2==null ? null : read(label+"/2", bases2, 0);
		if(second!=null){pair(first, second);}
		final int keep=proc.processReadPair(first, second);
		check(keep==expectedKeep, label+": expected keep bits "+expectedKeep+", observed "+keep);
		final String actual1=new String(first.bases, StandardCharsets.US_ASCII);
		final String actual2=second==null ? null : new String(second.bases, StandardCharsets.US_ASCII);
		check(expected1.equals(actual1) && (expected2==null ? actual2==null : expected2.equals(actual2)),
			label+": unexpected remaining bases "+actual1+" / "+actual2);
		final String stats=proc.toStats().toString();
		check(expectedStats.equals(stats), label+": expected statistics\n"+expectedStats+"observed\n"+stats);
		System.out.println("force_trim_case\t"+label+"\tkeep="+keep+"\tleft="+actual1+"\tright="+actual2);
		System.out.print(stats);
	}

	/** Uses the actual parser transfer/post-parse lifecycle with explicit unrestricted lengths. */
	private static ReformatProcessor processor(final boolean unique){
		final ReformatProcessor proc=new ReformatProcessor();
		proc.setFromParser(new Parser());
		proc.minReadLength=0;
		proc.maxReadLength=0;
		proc.uniqueNames=unique;
		proc.postParse();
		return proc;
	}

	/** Makes a fresh finite read so each production call owns its input mutations. */
	private static Read read(final String name, final String bases, final long id){
		return new Read(bases.getBytes(StandardCharsets.US_ASCII), null, name, id);
	}

	/** Installs the same bidirectional mate/ordinal relationship used by read streamers. */
	private static void pair(final Read left, final Read right){
		left.mate=right;
		right.mate=left;
		left.setPairnum(0);
		right.setPairnum(1);
	}

	/** Independent literal controls must fail loudly even if assertions were accidentally disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	private static int checks;
}
