package test.prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import consensus.BaseGraph;
import consensus.BaseGraphHelper;
import dna.AminoAcid;
import map.LongHashSet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.ConservedRnaSeedIndex;
import prok.NcrnaDdlOracle;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaScavenger;
import prok.NcrnaWorkloadInstrumentSink;
import prok.Orf;
import shared.Shared;

/** Independent fixed-model verification must retain ordinary accepted spans,
 * including physical-edge clipping, without inheriting mapping restrictions or
 * leaking its single-model choice into later ordinary calls.
 * Synthetic references and literal coordinates avoid a learned-resource oracle.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaModelOracleTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaModelOracleTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
			}
			check(parser.in1==null && parser.in2==null && parser.out1==null && parser.out2==null
				&& !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"The finite oracle fixture owns synthetic inputs and reports only its test result");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Model-oracle validation requires -ea for the requested-model invariant");
			new NcrnaModelOracleTest().run();
			pp.outstream.println("NCRNA_MODEL_ORACLE_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private NcrnaModelOracleTest(){
		refs=new byte[][]{dna(160, 12345), dna(160, 98765)};
		seeds.add(pack(refs[0], 0));
		seeds.add(pack(refs[0], 70));
		seeds.add(pack(refs[1], 0));
		seeds.add(pack(refs[1], 70));
		index=new ConservedRnaSeedIndex(new LongHashSet[]{seeds});
	}

	private void run(){
		for(int model=0; model<refs.length; model++){
			final byte[] interior=new byte[500], edge=new byte[400];
			Arrays.fill(interior, (byte)'N');
			Arrays.fill(edge, (byte)'N');
			System.arraycopy(refs[model], 0, interior, 120, 160);
			System.arraycopy(refs[model], 60, edge, 0, 100);
			final ConservedRnaSeedIndex.ScanResult full=index.scan(interior), clipped=index.scan(edge);
			final int[] interiorHits=full.hits(0), edgeHits=clipped.hits(0);
			check(Arrays.equals(interiorHits, new int[]{128, 198}), "Synthetic full-model seed centers are128 and198");
			check(Arrays.equals(edgeHits, new int[]{18}), "The retained seed at model70 maps to contig center18");
			for(int strand=0; strand<2; strand++){
				checkCase("interior_"+model+"_"+strand, strand, model, interior, interiorHits, full.keys(0), 98, 315, 120, 279, 0);
				checkCase("edge_"+model+"_"+strand, strand, model, edge, edgeHits, clipped.keys(0), 0, 135, 0, 99, 1);
			}
		}
		checkHbmRescue();
		checkModelOrder();
		checkWindowResult();
		checkTandemAttribution();
	}

	/** A mutant fails identity-only verification but passes an initialized HBM. */
	private void checkHbmRescue(){
		final BaseGraph[] graphs=new BaseGraph[refs.length];
		for(int i=0; i<graphs.length; i++){
			graphs[i]=new BaseGraph(NAMES[i], refs[i], null, i, 0);
			BaseGraphHelper.initForScoring(graphs[i]);
		}
		for(int selected=0; selected<refs.length; selected++){
			final byte[] mutant=refs[selected].clone(), bases=new byte[500];
			for(final int position:new int[]{30, 50, 90, 110, 130, 150}){
				mutant[position]=(byte)(mutant[position]=='A' ? 'C' : 'A');
			}
			Arrays.fill(bases, (byte)'N');
			System.arraycopy(mutant, 0, bases, 120, mutant.length);
			final ConservedRnaSeedIndex.ScanResult scan=index.scan(bases);
			check(Arrays.equals(scan.hits(0), new int[]{128, 198}), "HBM mutations must preserve both independent seed anchors");
			final int[] model={selected};
			final NcrnaScavenger identityOnly=caller(1000, null, 160, .99f, .5f);
			check(identityOnly.verifyOracleModel("hbm", bases, 0, 1, 0, 305, 0, scan.hits(0), scan.keys(0), model)==null,
				"Neither primary nor Scrabble identity may explain the mutant's HBM rescue");
			final NcrnaScavenger ordinary=caller(1, graphs, 160, .99f, .5f);
			final ArrayList<Orf> expected=ordinary.scavenge("hbm", bases, 0, new ArrayList<int[]>());
			check(expected.size()==1 && ordinary.hbmScoreCalls()>0,
				"The separate ordinary control must genuinely score and accept the HBM-rescued mutant");
			final NcrnaScavenger oracle=caller(1000, graphs, 160, .99f, .5f);
			final Orf found=oracle.verifyOracleModel("hbm", bases, 0, 1, 0, 305, 0, scan.hits(0), scan.keys(0), model);
			check(found!=null && found.start==120 && found.stop==279 && found.trnaModel.equals(NAMES[selected])
				&& found.start==expected.get(0).start && found.stop==expected.get(0).stop && oracle.hbmScoreCalls()==1,
				"Fixed-model HBM rescue must preserve its requested original model and literal accepted span");
			final NcrnaDdlOracle batch=new NcrnaDdlOracle("hbm_input", caller(1000, graphs, 160, .99f, .5f), NAMES);
			check(batch.verify("hbm", bases, 0, 1, 0, 305, 0, scan.hits(0), scan.keys(0))>=1
				&& batch.start(selected)==121 && batch.stop(selected)==280
				&& Float.compare(batch.score(selected), found.orfScore)==0,
				"Cached exhaustive results must retain full HBM acceptance and one-based forward coordinates");
		}
	}

	/** Seedless windows must not conceal a missing or reordered model identity. */
	private void checkModelOrder(){
		final byte[] bases=new byte[500];
		Arrays.fill(bases, (byte)'N');
		final ConservedRnaSeedIndex.ScanResult scan=index.scan(bases);
		final NcrnaScavenger worker=caller(1000);
		final NcrnaDdlOracle complete=new NcrnaDdlOracle("order_input", worker, NAMES);
		check(complete.verify("seedless", bases, 0, 1, 98, 315, 0, scan.hits(0), scan.keys(0))==0,
			"The complete library control must reject this window before any accepted-model assertion can fire");
		for(final String[] names:new String[][]{{NAMES[0]}, {NAMES[1], NAMES[0]},
				{NAMES[0], NAMES[1], "extra"}, {NAMES[0], null}}){
			boolean rejected=false;
			try{new NcrnaDdlOracle("order_input", worker, names);}
			catch(IllegalArgumentException expected){rejected=true;}
			check(rejected, "Missing, reordered, extra or unnamed models must fail at construction, independently of window acceptance");
		}
		check(worker.alignmentCount()==0, "The seedless model-order fixture must not rely on any positive alignment");
	}

	/** Two model alternatives, two strand orientations and repeated replacements
	 * use independent literal spans. Tandem links remain links, not recoveries. */
	private void checkWindowResult(){
		final String name="chr description\twith fields";
		check(NcrnaDdlOracle.seqid(name).equals("chr") && NcrnaDdlOracle.seqid("chr\tdescription").equals("chr")
			&& NcrnaDdlOracle.seqid("plain").equals("plain"),
			"Reference joins must use native first-space/tab tokens while the caller retains its full title");
		final NcrnaDdlOracle batch=new NcrnaDdlOracle("exact_input", caller(1000), NAMES);
		boolean rejected=false;
		try{batch.start(0);}catch(IllegalStateException expected){rejected=true;}
		check(rejected, "An unverified result cannot look like a zero-coordinate acceptance");
		final boolean[] compatible=new boolean[2];
		for(int strand=0; strand<2; strand++){
			for(int selected=0; selected<2; selected++){
				final byte[] bases=new byte[500];
				Arrays.fill(bases, (byte)'N');
				System.arraycopy(refs[selected], 0, bases, 120, 160);
				final byte[] saved=bases.clone();
				final ConservedRnaSeedIndex.ScanResult scan=index.scan(bases);
				check(batch.verify(name, bases, strand, selected+1, 98, 315, 0, scan.hits(0), scan.keys(0))==1,
					"Exactly the matching original model passes each exhaustive window");
				final int start=strand==0 ? 121 : 221, stop=strand==0 ? 280 : 380;
				final char direction=strand==0 ? '+' : '-';
				check(batch.models()==2 && batch.start(selected)==start && batch.stop(selected)==stop
					&& batch.partial(selected)==0 && batch.start(1-selected)==-1 && Float.isNaN(batch.score(1-selected)),
					"Replacement must clear the preceding model and convert minus coordinates exactly once");
				rejected=false;
				try{batch.parentCount(selected);}catch(IllegalStateException expected){rejected=true;}
				check(rejected, "A new window must not inherit the preceding window's parent attribution");
				check(batch.attribute("exact_input", new Locus[]{new Locus("chr", "cobalamin", start, stop, direction)})==0
					&& batch.parentCount(selected)==1 && batch.uniqueParent(selected)==0,
					"One exact scored parent provides one uniquely attributable compatibility link");
				check(batch.compatible("exact_input", new Locus("chr", "cobalamin", start, stop, direction), compatible)==1
					&& compatible[selected] && !compatible[1-selected], "Qualified exact parent must bind one actual accepted model");
				check(batch.compatible("exact_input", new Locus("other", "cobalamin", start, stop, direction), compatible)==0
					&& !compatible[0] && !compatible[1], "A different contig must clear all preceding compatibility flags");
				check(batch.compatible("exact_input", new Locus("chr", "cobalamin", start, start+79, direction), compatible)==1
					&& batch.compatible("exact_input", new Locus("chr", "cobalamin", start+80, stop, direction), compatible)==1,
					"The exact half-overlap tandem ambiguity must retain both compatibility links without claiming two recoveries");
				rejected=false;
				try{batch.compatible("another_input", new Locus("chr", "cobalamin", start, stop, direction), compatible);}
				catch(IllegalArgumentException expected){rejected=true;}
				check(rejected && Arrays.equals(bases, saved), "Exact genome identity is mandatory and original contig bytes remain read-only");
			}
		}
		final byte[] empty=new byte[500];
		Arrays.fill(empty, (byte)'N');
		final ConservedRnaSeedIndex.ScanResult noSeeds=index.scan(empty);
		check(batch.verify(name, empty, 0, 1, 98, 315, 0, noSeeds.hits(0), noSeeds.keys(0))==0
			&& batch.start(0)==-1 && batch.start(1)==-1 && Float.isNaN(batch.score(0)) && Float.isNaN(batch.score(1)),
			"A complete seedless replacement must clear every preceding accepted model, not retain a previous span");
		rejected=false;
		try{batch.verify(null, new byte[10], 0, 1, 0, 9, 0, new int[0], null);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Malformed replacement must fail before it can publish partial oracle state");
		rejected=false;
		try{batch.start(0);}catch(IllegalStateException expected){rejected=true;}
		check(rejected, "A failed replacement must invalidate the preceding complete window result");
	}

	/** One broad compatible span cannot uniquely recover both tandem parents.
	 * A second parent-specific span can recover its own parent, in either order. */
	private static void checkTandemAttribution(){
		final Locus left=new Locus("chr", "cobalamin", 101, 150, '+');
		final Locus right=new Locus("chr", "cobalamin", 151, 200, '+');
		final int[] starts={101, -1}, stops={200, -1}, counts=new int[2], unique=new int[2];
		for(final Locus[] parents:new Locus[][]{{left, right}, {right, left}}){
			check(NcrnaDdlOracle.attribute("chr", '+', starts, stops, parents, counts, unique)==1
				&& counts[0]==2 && unique[0]==-1 && counts[1]==0 && unique[1]==-1,
				"A shared span keeps both compatibility links and neither unique recovery, independent of reference order");
			starts[1]=151;
			stops[1]=200;
			check(NcrnaDdlOracle.attribute("chr", '+', starts, stops, parents, counts, unique)==1
				&& counts[0]==2 && unique[0]==-1 && counts[1]==1 && parents[unique[1]]==right,
				"A parent-specific alternative can uniquely recover the right parent without crediting the broad span twice");
			starts[1]=stops[1]=-1;
		}
		check(NcrnaDdlOracle.attribute("other", '+', starts, stops, new Locus[]{left, right}, counts, unique)==0
			&& Arrays.equals(counts, new int[2]) && Arrays.equals(unique, new int[]{-1, -1}),
			"Attribution replacement must clear links from another contig");
	}

	/** Ordinary calling supplies a separate control; exact coordinates are literal. */
	private void checkCase(final String name, final int strand, final int selected, final byte[] bases,
			final int[] hits, final long[] keys,
			final int start, final int stop, final int expectedStart, final int expectedStop, final int partial){
		final NcrnaScavenger ordinary=caller(1);
		final SeedCapture observed=new SeedCapture(start, stop);
		ordinary.setWorkloadSink(observed);
		final ArrayList<Orf> expected=ordinary.scavenge(name, bases, strand, new ArrayList<int[]>());
		check(Arrays.equals(observed.positions, hits) && observed.scheduled,
			"Recomputed shared-index seed positions must equal the ordinary caller's independent scan for this window");
		check(keys.length==0,
			"GeneCaller retains no seed keys for RF00174 without vote/join routes; the oracle must preserve that empty stream");
		check(expected.size()==1, "The ordinary control must produce exactly one synthetic locus: "+name);
		check(expected.get(0).start==expectedStart && expected.get(0).stop==expectedStop
			&& expected.get(0).ncrnaPartial==partial, "Native control must match independently specified bounds and clipping: "+name);
		final NcrnaScavenger oracle=caller(1000);
		check(oracle.scavenge(name, bases, strand, new ArrayList<int[]>(), hits).isEmpty(),
			"The intentionally rejecting index must block ordinary calling before oracle verification");
		final int[] model={0};
		for(int round=0; round<3; round++){
			model[0]=selected;
			final Orf found=oracle.verifyOracleModel(name, bases, strand, 1, start, stop, 0, hits, keys, model);
			check(found!=null && found.start==expectedStart && found.stop==expectedStop && found.ncrnaPartial==partial,
				"Full fixed-model verification must preserve native interior/edge semantics: "+name);
			check(found.trnaModel.equals(NAMES[selected]) && Float.compare(found.orfScore, expected.get(0).orfScore)==0,
				"Fixed-model acceptance must retain its original model identity and ordinary score");
			model[0]=1-selected;
			check(oracle.verifyOracleModel(name, bases, strand, 1, start, stop, 0, hits, keys, model)==null,
				"An unrelated model must not inherit the preceding model's acceptance or position scratch");
		}
		check(oracle.scavenge(name, bases, strand, new ArrayList<int[]>(), hits).isEmpty(),
			"The oracle must not disable index filtering in later ordinary calls");
	}

	private NcrnaScavenger caller(final int minHits){
		return caller(minHits, null, 100, .95f, .95f);
	}

	private NcrnaScavenger caller(final int minHits, final BaseGraph[] models, final int pad,
			final float pass, final float borderline){
		final NcrnaScavenger caller=new NcrnaScavenger(refs, models, NAMES, seeds, 17, 20, pad,
			9, 2, false, 0f, 0f, 0f, minHits, 0f, 1f, pass, borderline);
		caller.family="cobalamin_rf00174";
		caller.minKmerHits=1;
		caller.strictIndexCutoff=true;
		caller.physicalEdgeMinLen=20;
		caller.hbmPass=models==null ? 1.01f : .5f;
		return caller;
	}

	/** Uses the seed index's A/C/G/T encoding, not cardinality's alternate encoding. */
	private static long pack(final byte[] bases, final int start){
		check(start>=0 && start+17<=bases.length, "A fixture seed must be a complete17-mer");
		long key=0;
		for(int i=start; i<start+17; i++){key=(key<<2)|AminoAcid.baseToNumber[bases[i]];}
		return key;
	}

	private static byte[] dna(final int length, final long seed){
		check(length>=17, "Synthetic references must contain native17-mer seeds");
		final byte[] bases=new byte[length];
		long state=seed;
		for(int i=0; i<length; i++){
			state^=state<<13;
			state^=state>>>7;
			state^=state<<17;
			bases[i]=AminoAcid.numberToBase[(int)(state&3)];
		}
		return bases;
	}

	private static void check(final boolean valid, final String why){
		if(!valid){throw new AssertionError(why);}
	}

	/** Observes the ordinary per-family scan, not the recomputed shared-index input. */
	private static final class SeedCapture implements NcrnaWorkloadInstrumentSink{
		SeedCapture(final int start_, final int stop_){
			start=start_;
			stop=stop_;
		}
		@Override
		public void seedHits(final String name, final int strand, final int[] hits){positions=hits;}
		@Override
		public void scheduledWindow(final String name, final int strand, final int pass, final int a, final int b){
			if(pass==1 && a==start && b==stop){scheduled=true;}
		}
		final int start, stop;
		int[] positions;
		boolean scheduled;
	}
	private final byte[][] refs;
	private final LongHashSet seeds=new LongHashSet();
	private final ConservedRnaSeedIndex index;
	private static final String[] NAMES={"cobalamin_RF00174_fixture_0", "cobalamin_RF00174_fixture_1"};
}
