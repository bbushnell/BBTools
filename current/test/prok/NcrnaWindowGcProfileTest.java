package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Locale;

import map.LongHashSet;
import cardinality.DynamicDemiLog;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlInterleave;
import prok.NcrnaDdlOrf;
import prok.BuildNcrnaReference;
import prok.NcrnaDdlProfileResource;
import prok.NcrnaDdlShortlist;
import prok.NcrnaTextIO.Output;
import prok.NcrnaWindowGcProfile;
import prok.GeneCaller;
import prok.GeneModel;
import prok.NcrnaFamily;
import prok.NcrnaScavenger;
import prok.Orf;
import shared.Shared;
import structures.ByteBuilder;

/** Finite GC resource/scorer/ordering fixture, adapted from
 * template.A_SampleBasic2's argument and output lifecycle. Model rows are
 * deliberately reordered. Tests physical zero/unknown states, exact fractional
 * reload, scientific notation, hard cuts, negative/tied ordering and scratch
 * reuse through the real selected profile scorer and shortlist.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaWindowGcProfileTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaWindowGcProfileTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown GC fixture argument: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append,
				"Supply out=fresh_fixture_directory without input, append or overwrite");
			new NcrnaWindowGcProfileTest(new File(parser.out1).getAbsoluteFile()).run();
			pp.outstream.println("NCRNA_WINDOW_GC_PROFILE_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private NcrnaWindowGcProfileTest(final File root_){
		root=root_;
		check(root.mkdir(), "The GC fixture requires a fresh output directory");
	}

	private void run(){
		final ByteBuilder escaped=BuildNcrnaReference.field(new ByteBuilder(), E);
		final String encodedName=new String(escaped.array, 0, escaped.length(), StandardCharsets.UTF_8);
		final String table=NcrnaWindowGcProfile.HEADER+"\n"
			+encodedName+"\t2\t0\t0.5\t1.0E-7\t0.4999999\t0.5000001\tFIT\n"
			+C+"\t0\t1\tNA\tNA\tNA\tNA\tNO_KNOWN_GC_MEMBERS\n"
			+A+"\t2\t0\t0.5\t0.1\t0.4\t0.6\tFIT\n"
			+D+"\t1\t0\t0.3333333333333333\t0.0\t0.3333333333333333\t0.3333333333333333\tZERO_SD\n"
			+B+"\t3\t0\t0.5\t0.0\t0.5\t0.5\tZERO_SD\n";
		write("gc.tsv", table);
		final String[] names={A+" description", B, C, D, E};
		final NcrnaWindowGcProfile sd=NcrnaWindowGcProfile.load(path("gc.tsv"), names, 2, 1, false);
		check(sd.models()==5 && sd.knownModels()==4, "Reordered rows must bind all actual models and retain the explicit unfitted row");
		final float[] scores=new float[5];
		final boolean[] eligible=new boolean[5];
		check(sd.apply(bytes("GGGAAN"), scores, eligible)==3, "GC0.6 should retain broad model A and unknown model C only");
		check(eligible[0] && !eligible[1] && eligible[2] && !eligible[3] && !eligible[4]
			&& Math.abs(scores[0]+2)<1e-6 && scores[2]==0,
			"The physical SD cut must preserve a negative adjusted score and leave unknown-model scores unchanged");
		final NcrnaDdlInterleave order=new NcrnaDdlInterleave(5);
		check(order.buildFilteredSketchCaller(scores, eligible, 5)==2 && order.modelAt(0)==2 && order.modelAt(1)==0,
			"A valid negative-score proposal must follow the unchanged zero-score unknown model");
		Arrays.fill(scores, 7);
		check(sd.apply(bytes("acgtNN"), scores, eligible)==1 && eligible[0] && eligible[1] && eligible[2] && !eligible[3] && eligible[4],
			"Known-base denominator and exact zero-SD matching must agree for upper/lowercase ACGT");
		check(scores[0]==7 && scores[1]==7 && scores[4]==7, "Exact means incur no GC penalty, including scientific-notation SD");
		Arrays.fill(scores, 3);
		check(sd.apply(bytes("gAA"), scores, eligible)==3 && eligible[3], "The one-third profile must admit its original nonterminating fraction exactly");
		Arrays.fill(scores, 9);
		Arrays.fill(eligible, false);
		check(sd.apply(bytes("NNRYU"), scores, eligible)==0, "A window without ACGT has unknown GC, not zero GC");
		for(int i=0; i<5; i++){check(eligible[i] && scores[i]==9, "Unknown-window fallback must replace stale eligibility and preserve all raw scores");}
		final NcrnaWindowGcProfile range=NcrnaWindowGcProfile.load(path("gc.tsv"), names, .5, Double.NaN, true);
		Arrays.fill(scores, 0);
		range.apply(bytes("GGGAAN"), scores, eligible);
		check(eligible[0] && Math.abs(scores[0]+.5)<1e-6, "The inclusive observed maximum must survive while its penalty is retained");
		final byte[] outside=new byte[100];
		Arrays.fill(outside, (byte)'A');
		Arrays.fill(outside, 0, 61, (byte)'G');
		Arrays.fill(scores, 0);
		range.apply(outside, scores, eligible);
		check(!eligible[0], "GC0.61 lies outside the saved0.4..0.6range");
		checkSignedOrder();
		checkResourceScorer();
		checkRealVerifier();
		reject("duplicate", table+A+"\t2\t0\t0.5\t0.1\t0.4\t0.6\tFIT\n", names);
		reject("missing", table.substring(0, table.lastIndexOf(B)), names);
		reject("foreign", table.replace(encodedName+"\t", "cobalamin_RF00174_foreign\t"), names);
		reject("negative", table.replace("1.0E-7", "-0.1"), names);
		reject("nan", table.replace("1.0E-7", "NaN"), names);
		reject("false_zero", table.replace("0.0\t0.5\t0.5\tZERO_SD", "0.0\t0.4\t0.6\tZERO_SD"), names);
	}

	/** The real GeneCaller factory passes the shared optional resource to its
	 * actual Quantum verifier. Both exact models pass without GC; excluding the
	 * first must remove exactly one executed alignment and change the winner. */
	private void checkRealVerifier(){
		final byte[] query=Prok16sPartialRescueTest.random(320, 20261116);
		final byte[][] sequences={query.clone(), query.clone()};
		final String[] names={A, C};
		NcrnaDdlProfileResource.writeSelected(path("verifier_profiles.tsv"), new char[2][512],
			new int[]{1,1}, names, sequences, 512, new LongHashSet(4));
		final NcrnaDdlProfileResource raw=NcrnaDdlProfileResource.load(path("verifier_profiles.tsv"), names, sequences);
		int gc=0;
		for(final byte base:query){if(base=='G' || base=='C'){gc++;}}
		final double wrongMean=gc/(double)query.length==.5 ? 0 : .5;
		final double fraction=gc/(double)query.length;
		write("verifier_gc.tsv", NcrnaWindowGcProfile.HEADER+"\n"+A+"\t1\t0\t"+wrongMean
			+"\t0.0\t"+wrongMean+"\t"+wrongMean+"\tZERO_SD\n"
			+C+"\t2\t0\t"+(fraction-.05)+"\t0.05\t"+(fraction-.1)+"\t"+fraction+"\tFIT\n");
		final NcrnaDdlProfileResource filtered=raw.withGc(NcrnaWindowGcProfile.load(path("verifier_gc.tsv"), names, 1, 3, false));
		final ArrayList<NcrnaFamily> saved=new ArrayList<NcrnaFamily>(GeneCaller.ncrnaFamilies);
		final int exponent=DynamicDemiLog.exponentBits();
		try{
			DynamicDemiLog.setExponent(5);
			final LongHashSet seeds=new LongHashSet(4);
			long seed=0;
			for(int base=0; base<17; base++){seed=(seed<<2)|dna.AminoAcid.baseToNumber[query[base]];}
			seeds.add(seed);
			final NcrnaFamily family=new NcrnaFamily("cobalamin_rf00174", sequences, null, names, seeds,
				17, 109, query.length, 9, Integer.MAX_VALUE, true, 2f, 0f, 0f, 0, 0f, 20f, .99f, .90f);
			family.modelPatience=Integer.MAX_VALUE;
			GeneCaller.ncrnaFamilies.clear();
			GeneCaller.ncrnaFamilies.add(family);
			final NcrnaScavenger before=worker(raw), after=worker(filtered);
			final ArrayList<Orf> unadjusted=before.scavenge("without_gc", query, 0, new ArrayList<int[]>());
			final ArrayList<Orf> adjusted=after.scavenge("with_gc", query, 0, new ArrayList<int[]>());
			check(unadjusted.size()==1 && adjusted.size()==1 && A.equals(unadjusted.get(0).trnaModel)
				&& C.equals(adjusted.get(0).trnaModel) && before.alignmentCount()==2 && after.alignmentCount()==1,
				"The real optional GC route must skip the excluded exact model before alignment, execute the retained model and preserve one acceptance");
			final NcrnaDdlOrf winner=(NcrnaDdlOrf)adjusted.get(0);
			check(winner.proposalScore()<0 && Math.abs(winner.proposalScore()+1)<1e-5
				&& winner.toGff().contains(";ddl_gc_adjusted=true"),
				"An actually accepted negative-score GC winner must retain finite signed evidence through GFF output");
			check(DynamicDemiLog.exponentBits()==5, "GC scoring must not alter the caller's process-global key encoding");
		}finally{
			GeneCaller.ncrnaFamilies.clear();
			GeneCaller.ncrnaFamilies.addAll(saved);
			DynamicDemiLog.setExponent(exponent);
		}
	}

	private static NcrnaScavenger worker(final NcrnaDdlProfileResource profiles){
		final GeneCaller caller=new GeneCaller(60, 10, 10, -1f, -1f, -1f, -1f, -1f, new GeneModel(false));
		caller.setCobalaminDdl(profiles, Integer.MAX_VALUE);
		final NcrnaScavenger worker=caller.ncrnaScavenger(0);
		worker.trimAlignmentExtent=false;
		worker.scavengePass2=false;
		return worker;
	}

	private void checkSignedOrder(){
		final NcrnaDdlInterleave order=new NcrnaDdlInterleave(8);
		final float[] scores={-1, 0, -0.0f, 2, -Float.MIN_VALUE, -Float.MAX_VALUE, Float.MIN_VALUE, -2};
		final boolean[] eligible=new boolean[8];
		Arrays.fill(eligible, true);
		final int[] expected={3,6,1,2,4,0,7,5};
		check(order.buildFilteredSketchCaller(scores, eligible, 8)==8, "Every finite signed score must reach its proper order");
		for(int i=0; i<8; i++){check(order.modelAt(i)==expected[i], "Signed float order, zero ties and subnormal ordering must preserve numerical rank");}
		check(order.buildFilteredSketchCaller(scores, eligible, 3)==3 && order.modelAt(2)==1,
			"The finite cap counts eligible model positions after GC filtering");
		Arrays.fill(eligible, false);
		check(order.buildFilteredSketchCaller(scores, eligible, 8)==0, "A completely filtered query has an empty live prefix");
		for(int i=0; i<8; i++){check(order.rank(i)==0, "An empty query must clear every earlier emitted rank");}
		boolean rejected=false;
		try{order.buildSketchCaller(scores, 8);}catch(IllegalArgumentException expectedFailure){rejected=true;}
		check(rejected, "The original nonnegative sketch-only API must retain its validation");
	}

	/** Real resource/scorer/shortlist plumbing: zero-GC profiles reject one actual
	 * selected model while a named unknown model stays; omitted GC stays unchanged. */
	private void checkResourceScorer(){
		final String[] names={A, C};
		final byte[][] bases={bytes("ACGTACGTACGTACGTACGT"), bytes("TGCATGCATGCATGCATGCA")};
		final char[][] keys=new char[2][1024];
		Arrays.fill(keys[0], (char)17);
		Arrays.fill(keys[1], (char)31);
		NcrnaDdlProfileResource.writeSelected(path("profiles.tsv"), keys, new int[]{1,1}, names, bases, 1024, new LongHashSet(4));
		final NcrnaDdlProfileResource raw=NcrnaDdlProfileResource.load(path("profiles.tsv"), names, bases);
		write("two_gc.tsv", NcrnaWindowGcProfile.HEADER+"\n"+C+"\t0\t0\tNA\tNA\tNA\tNA\tNO_KNOWN_GC_MEMBERS\n"
			+A+"\t1\t0\t0.5\t0.0\t0.5\t0.5\tZERO_SD\n");
		final NcrnaDdlProfileResource filtered=raw.withGc(NcrnaWindowGcProfile.load(path("two_gc.tsv"), names, 1, 3, false));
		final NcrnaDdlShortlist before=new NcrnaDdlShortlist(raw, Integer.MAX_VALUE);
		final NcrnaDdlShortlist after=new NcrnaDdlShortlist(filtered, Integer.MAX_VALUE);
		before.build(bytes("GGGGGGGGGGGGGGGGGGGG"), null);
		final int[] emitted=after.build(bytes("GGGGGGGGGGGGGGGGGGGG"), null);
		check(before.size()==2 && after.size()==1 && emitted[0]==1 && after.rank(0)==0 && after.rank(1)==1,
			"The optional immutable GC resource must reach the actual selected scorer and emitted live prefix");
		final int[] reused=after.build(bytes("NNNNNNNNNNNNNNNNNNNN"), null);
		check(reused==emitted && after.size()==2, "Unknown-window scoring must reuse scratch and clear the previous GC exclusion");
		check(raw.key(0, 0)==17 && filtered.key(0, 0)==17, "Adding GC must share immutable profile content without changing the original resource");
	}

	private void reject(final String name, final String table, final String[] names){
		write(name+".tsv", table);
		boolean rejected=false;
		try{NcrnaWindowGcProfile.load(path(name+".tsv"), names, 1, 3, false);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Malformed/incomplete GC resource must fail loudly: "+name);
	}
	private void write(final String name, final String value){
		// ByteBuilder(String) and toString() are ASCII APIs; metadata fixture names
		// must cross the same explicit UTF-8 boundary as the production writer.
		try(Output out=new Output(path(name))){out.writer.print(new ByteBuilder(value.getBytes(StandardCharsets.UTF_8)));}
	}
	private String path(final String name){return new File(root, name).getPath();}
	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private final File root;
	private static final String A="cobalamin_RF00174_A", B="cobalamin_RF00174_B", C="cobalamin_RF00174_C",
		D="cobalamin_RF00174_D", E="cobalamin_RF00174_E%25Ω";
}
