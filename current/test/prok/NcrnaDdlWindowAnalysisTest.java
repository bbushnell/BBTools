package test.prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import cardinality.DynamicDemiLog;
import dna.AminoAcid;
import map.LongHashSet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlOracle;
import prok.NcrnaDdlProfile;
import prok.NcrnaDdlWindowAnalysis;
import prok.NcrnaFamily;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaModelAttemptInstrumentSink;
import prok.NcrnaScavenger;
import prok.Orf;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;

/** Checks synchronous oracle/ranking observation against an ordinary control.
 * Synthetic genes have independent literal spans; full GFF and actual initial
 * attempt rows must agree with observation disabled. Two orientations exercise
 * the same Read's in-place reversal, and an empty contig keeps zero events.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlWindowAnalysisTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaDdlWindowAnalysisTest.class, false);
		final int previous=DynamicDemiLog.exponentBits();
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown window-analysis fixture parameter: "+arg);
			}
			check(parser.in1==null && parser.in2==null && parser.out1==null && parser.out2==null
				&& !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"The finite window-analysis fixture owns synthetic inputs without file output or sampling");
			NcrnaDdlWindowAnalysis.checkEnvironment();
			DynamicDemiLog.setExponent(4);
			new NcrnaDdlWindowAnalysisTest().run();
			pp.outstream.println("NCRNA_DDL_WINDOW_ANALYSIS_TEST_PASS");
		}finally{
			DynamicDemiLog.setExponent(previous);
			Shared.closeStream(pp.outstream);
		}
	}

	private NcrnaDdlWindowAnalysisTest(){
		final LongHashSet seeds=new LongHashSet();
		for(final byte[] ref:refs){
			seeds.add(pack(ref, 70));
		}
		family=new NcrnaFamily("cobalamin_rf00174", refs, null, NAMES,
			seeds, 17, 100, 100, 9, 2, false, 0f, 0f, 0f, 1,
			0f, 1f, .64f, .64f);
		family.scavengePass2=false;
		family.strictIndexCutoff=true;
		profiles=new NcrnaDdlProfile[KS.length][refs.length];
		for(int geometry=0; geometry<KS.length; geometry++){
			final DynamicDemiLog sketch=DynamicDemiLog.create(WIDTHS[geometry], KS[geometry], 12345, 0, false, false);
			for(int model=0; model<refs.length; model++){
				sketch.clear();
				sketch.hash(refs[model], null);
				final NcrnaDdlProfile.Votes vote=new NcrnaDdlProfile.Votes(WIDTHS[geometry]);
				vote.add(sketch.maxArray());
				profiles[geometry][model]=vote.finish();
			}
		}
	}

	private void run(){
		checkOrdinaryBinding();
		final NcrnaScavenger control=NcrnaDdlWindowAnalysis.oracleWorker(family);
		final NcrnaScavenger observed=NcrnaDdlWindowAnalysis.oracleWorker(family);
		final Attempts baselineAttempts=new Attempts();
		control.setModelAttemptSink(baselineAttempts);
		final Observations sink=new Observations();
		final NcrnaDdlWindowAnalysis analysis=new NcrnaDdlWindowAnalysis("I", family, observed, profiles, KS, sink);
		observed.setWindowCaptureSink((name, strand, pass, start, stop, source, bases, window)->{
			final byte[] before=window.clone();
			analysis.capture(name, strand, pass, start, stop, source, bases, window);
			check(Arrays.equals(before, window), "The synchronous oracle and ranking must preserve every borrowed window byte");
		});
		final byte[] plus=new byte[500];
		Arrays.fill(plus, (byte)'N');
		System.arraycopy(refs[1], 0, plus, 120, refs[1].length);
		final byte[] minus=plus.clone();
		AminoAcid.reverseComplementBasesInPlace(minus);
		final byte[] both=new byte[1200], reverse=refs[0].clone();
		Arrays.fill(both, (byte)'N');
		System.arraycopy(refs[1], 0, both, 120, 160);
		AminoAcid.reverseComplementBasesInPlace(reverse);
		System.arraycopy(reverse, 0, both, 700, 160);
		final ByteBuilder expected=new ByteBuilder(), actual=new ByteBuilder();
		for(int i=0; i<4; i++){
			final byte[] bases=i==0 ? plus : i==1 ? minus : i==2 ? new byte[]{'N', 'N', 'N'} : both;
			final String title=i==0 ? "plus description" : i==1 ? "minus\tdescription" : i==2 ? "empty" : "both description";
			final Read baseline=new Read(bases.clone(), null, title, i);
			final Read measured=new Read(bases.clone(), null, title, i);
			call(control, baseline, expected, null);
			call(observed, measured, actual, analysis);
			check(Arrays.equals(bases, baseline.bases) && Arrays.equals(bases, measured.bases),
				"Both ordinary calls must restore each original Read's exact bytes and orientation");
		}
		check(Arrays.equals(expected.toBytes(), actual.toBytes()),
			"Synchronous exhaustive verification/ranking must preserve complete ordinary GFF bytes");
		check(Arrays.equals(baselineAttempts.rows.toBytes(), sink.attempts.rows.toBytes()),
			"Observation must preserve every actual native attempt, rank, shared count, identity and initial-pass flag");
		check(analysis.events()==4 && sink.windows==4 && sink.ended==4 && sink.nonemptyEnds==4
			&& sink.bothStrands==3 && analysis.attempts()>0 && analysis.seedScans()==4,
			"Both orientations of one Read need distinct scans; the empty contig creates no window, attempt or seed scan");
		check(analysis.oracleCpu()>=0 && analysis.rankingCpu()>=0 && analysis.nativeCpu()>=0,
			"Analysis CPU counters must remain separate nonnegative current-worker measurements");
		boolean duplicate=false, stale=false;
		try{analysis.begin(new Read(plus.clone(), null, "plus another description", 4));}
		catch(IllegalArgumentException expectedFailure){duplicate=true;}
		try{analysis.oracle();}catch(IllegalArgumentException expectedFailure){stale=true;}
		check(duplicate && stale, "Duplicate native seqid tokens and out-of-callback oracle access must reject");
		// These two manual captures test only event identity, not ordinary calling.
		final Read repeatRead=new Read(plus.clone(), null, "repeat description", 5);
		analysis.begin(repeatRead);
		final byte[] window=Arrays.copyOfRange(plus, 20, 310);
		analysis.capture("repeat description", 0, 1, 20, 309, 0, repeatRead.bases, window);
		analysis.capture("repeat description", 0, 1, 20, 309, 0, repeatRead.bases, window);
		check(Arrays.equals(window, Arrays.copyOfRange(plus, 20, 310)), "Analysis must preserve borrowed window bytes across replacement");
		analysis.end();
		check(analysis.events()==6 && sink.ended==6 && sink.nonemptyEnds==4 && sink.lastEvent==5,
			"Identical coordinate tuples retain separate event ordinals, including windows with zero ordinary attempts");
		analysis.finish();
		check(control.alignmentCount()==observed.alignmentCount(),
			"Exhaustive oracle work must not increment the ordinary target worker's alignment counter");
	}

	/** An ordinary worker with no accepted window must still expose its wrong
	 * model population before any observation hook is installed. */
	private void checkOrdinaryBinding(){
		for(final String[] names:new String[][]{{NAMES[1], NAMES[0]}, {NAMES[0]}}){
			final byte[][] library=names.length==1 ? new byte[][]{refs[0]} : refs;
			final NcrnaFamily wrong=new NcrnaFamily("cobalamin_rf00174", library, null, names,
				family.kmerSet, 17, 100, 100, 9, 2, false, 0f, 0f, 0f, 1, 0f, 1f, .64f, .64f);
			final NcrnaScavenger ordinary=NcrnaDdlWindowAnalysis.oracleWorker(wrong);
			check(ordinary.scavenge("seedless", new byte[]{'N', 'N', 'N'}, 0, new ArrayList<int[]>()).isEmpty()
				&& ordinary.alignmentCount()==0, "The mismatched ordinary control must genuinely make zero model attempts");
			boolean rejected=false;
			try{new NcrnaDdlWindowAnalysis("I", family, ordinary, profiles, KS, new Observations());}
			catch(IllegalArgumentException expected){rejected=true;}
			check(rejected, "Missing or reordered ordinary model names must reject even on a seedless, zero-attempt input");
		}
	}

	/** Calls each strand once and restores the original Read before ending it. */
	private static void call(final NcrnaScavenger caller, final Read read, final ByteBuilder out,
			final NcrnaDdlWindowAnalysis analysis){
		if(analysis!=null){analysis.begin(read);}
		for(int strand=0; strand<2; strand++){
			if(strand==1){AminoAcid.reverseComplementBasesInPlace(read.bases);}
			final ArrayList<Orf> found=caller.scavenge(read.id, read.bases, strand, new ArrayList<int[]>());
			for(final Orf orf:found){orf.appendGff(out).nl();}
		}
		AminoAcid.reverseComplementBasesInPlace(read.bases);
		if(analysis!=null){analysis.end();}
	}

	private static final class Attempts implements NcrnaModelAttemptInstrumentSink{
		@Override
		public void modelAttempt(final String name, final int strand, final int pass, final int start, final int stop,
				final int rank, final int model, final int shared, final float identity, final int length, final boolean initialPass){
			rows.append(name).tab().append(strand).tab().append(pass).tab().append(start).tab().append(stop).tab()
				.append(rank).tab().append(model).tab().append(shared).tab().append(Float.floatToRawIntBits(identity)).tab().append(length).tab().append(initialPass).nl();
		}
		final ByteBuilder rows=new ByteBuilder();
	}

	private static final class Observations implements NcrnaDdlWindowAnalysis.Sink{
		@Override
		public void window(final NcrnaDdlWindowAnalysis view, final byte[] window){
			check(view.event()==windows && view.event()>lastEvent, "Event ordinals must increase independently of coordinate tuples");
			lastEvent=view.event();
			windows++;
			final boolean both=view.header().startsWith("both"), minus=both ? view.strand()==1 : view.header().startsWith("minus");
			final int model=both && minus ? 0 : 1;
			final int seed=both && minus ? 418 : 198;
			check(Arrays.equals(view.seedPositions(), new int[]{seed}), "One Read's strands must use their own literal oriented seed center198 or418");
			check(view.seedValues().length==1, "The retained native seed key aligns with its one position");
			if(both){bothStrands|=1<<view.strand();}
			final int first=both && minus ? 701 : minus ? 221 : 121;
			final Locus parent=new Locus(view.seqid(), "cobalamin", first, first+159, minus ? '-' : '+');
			final NcrnaDdlOracle oracle=view.oracle();
			check(oracle.compatible("I", parent, compatible)>0 && compatible[model],
				"The present model must fully verify the literal physical parent on the original strand");
			view.nativeBest(compatible, rank);
			check(rank[0]==1 && rank[3]==model, "The ordinary index must rank the exact synthetic model first");
			for(int geometry=0; geometry<KS.length; geometry++){
				for(int method=0; method<3; method++){
					view.best(geometry, method, compatible, rank);
					check(rank[0]==1 && rank[3]==model, "All profile methods must rank the exact model first on this unambiguous window");
				}
			}
			check(window.length==view.stop()-view.start()+1, "The borrowed window must retain its inclusive native length");
		}
		@Override
		public void attempt(final NcrnaDdlWindowAnalysis view, final int rank, final int model, final int shared,
				final float identity, final int length, final boolean initialPass){
			check(view.nativeShortlistRank(model)==rank, "Actual attempts must use the same-window native shortlist");
			attempts.modelAttempt(view.header(), view.strand(), view.pass(), view.start(), view.stop(),
				rank, model, shared, identity, length, initialPass);
		}
		@Override
		public void endWindow(final NcrnaDdlWindowAnalysis view){
			ended++;
			if(view.eventAttempts()>0){nonemptyEnds++;}
		}
		final Attempts attempts=new Attempts();
		final boolean[] compatible=new boolean[2];
		final int[] rank=new int[4];
		int windows, ended, nonemptyEnds, bothStrands;
		long lastEvent=-1;
	}

	private static byte[] dna(final int length, long value){
		final byte[] bases=new byte[length];
		for(int i=0; i<length; i++){
			value^=value<<13;
			value^=value>>>7;
			value^=value<<17;
			bases[i]=AminoAcid.numberToBase[(int)(value&3)];
		}
		return bases;
	}
	private static long pack(final byte[] bases, final int from){
		long value=0;
		for(int i=from; i<from+17; i++){value=(value<<2)|AminoAcid.baseToNumber[bases[i]];}
		return value;
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private final byte[][] refs={dna(160, 12345), dna(160, 98765)};
	private final NcrnaFamily family;
	private final NcrnaDdlProfile[][] profiles;
	private static final int[] KS={5, 9}, WIDTHS={32, 64};
	private static final String[] NAMES={"cobalamin_RF00174_fixture0", "cobalamin_RF00174_fixture1"};
}
