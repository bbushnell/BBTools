package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import consensus.BaseGraph;
import consensus.BaseGraphHelper;
import idaligner.QuantumAligner;
import idaligner.ScrabbleAligner;
import ml.CellNet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.CallGenes;
import prok.ConservedRnaSeedIndex;
import prok.FamilyResourceResolver;
import prok.GeneCaller;
import prok.NcrnaFamily;
import prok.NcrnaFamilyLoader;
import prok.NcrnaFamilyRegistry;
import prok.NcrnaScavenger;
import prok.NcrnaReportedOrf;
import prok.NcrnaTextIO.Output;
import prok.NcrnaWorkloadInstrumentSink;
import prok.NucleotideHbmCodec;
import prok.Orf;
import prok.TrnaBoundaryFeatures;
import shared.Shared;
import structures.ByteBuilder;

/** Finite loading controls for the default-off glycine caller.
 * Writes three literal paired models, then checks loading, overrides and failures.
 * Run in a fresh JVM. This fixture does not measure scientific caller accuracy.
 * @author Brian Bushnell, Raiden
 */
public final class CallGenesGlycineTest{

	/** Uses the native parser for out= and closes the preparser-owned log. */
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CallGenesGlycineTest.class, false);
		try{
			final Parser parser=new Parser();
			for(String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown test option: "+arg);}
			}
			check(parser.out1!=null && parser.in1==null, "Supply out=fresh_directory for the literal fixture");
			final File root=new File(parser.out1);
			check(root.mkdir(), "Fixture output must be a new directory: "+root);
			run(root);
			pp.outstream.println("GLYCINE_LOADING_CONTROLS_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Checks the direct and caller factories against the same independent fixture. */
	private static void run(final File root){
		check(!CallGenes.GLYCINE_ENABLED && GeneCaller.ncrnaFamilies.isEmpty(), "Fresh startup must leave glycine disabled");
		final String absent=new File(root, "absent").getPath();
		CallGenes.parseGlycineFlag("glycineroot", absent);
		CallGenes.loadGlycineResources();
		check(GeneCaller.ncrnaFamilies.isEmpty(), "Disabled glycine must not open an explicitly missing root");
		check("glycine".equals(CallGenes.parseNcrnaFamily("GLYCINE")), "The selector must bind the registered family");
		check("Glycine".equals(NcrnaFamilyRegistry.displayName("glycine")), "GFF attribution must identify glycine aptamers");
		check("glycine".equals(NcrnaFamilyRegistry.modelOverrideMember("glycine").id), "Paired-model tools must recognize glycine");
		check("glycine".equals(NcrnaFamilyRegistry.sharedBoundaryMember("glycine").id), "Later boundary tools must share the registration");

		final File family=new File(root, "family");
		writeFixture(family, "glycine_");
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
		options.resources.root=family.getPath();
		final NcrnaFamily direct=NcrnaFamilyLoader.load("glycine", options, new FamilyResourceResolver());
		check(direct.library.length==3 && direct.indexTopN==58 && direct.windowPad==222,
			"A smaller explicit library must retain the selected production cap and padding");
		check(direct.adaptive && direct.adaptFloor==9f && direct.adaptTopFrac==.390625f && direct.adaptQFrac==.030f,
			"The selected adaptive mapping policy must reach direct loading, including replacement libraries");
		check(direct.kLong==ConservedRnaSeedIndex.K && direct.seedMinHits==2 && direct.seedDistinct,
			"Direct loading must retain the shared seed scan and distinct two-key gate");
		check(direct.minLen==17 && direct.boundary5NetTemplate==null && direct.cdsStopConstraintFactory==null,
			"The experimental seed-length floor must not acquire endpoint nets or a CDS-stop constraint");
		CallGenes.parseGlycineFlag("glycine", "t");
		CallGenes.NCRNA_FAMILIES_ENABLED=true;
		expectMissing();
		CallGenes.parseGlycineFlag("glycineroot", family.getPath());
		CallGenes.loadGlycineResources();
		check(GeneCaller.ncrnaFamilies.size()==1, "The factory must publish exactly one glycine family");
		final NcrnaFamily called=GeneCaller.ncrnaFamilies.get(0);
		check("glycine".equals(called.name) && called.windowPad==direct.windowPad && called.indexTopN==direct.indexTopN
			&& called.adaptive==direct.adaptive && called.adaptFloor==direct.adaptFloor
			&& called.adaptTopFrac==direct.adaptTopFrac && called.adaptQFrac==direct.adaptQFrac
			&& called.seedMinHits==2 && called.seedDistinct
			&& Arrays.deepEquals(called.library, direct.library) && Arrays.equals(called.modelNames, direct.modelNames),
			"Caller and direct loading must consume the same paired library and controls");
		verifyQuantum(called, false);
		verifyQuantum(called, true);
		verifyBoundary(called, false, 0);
		verifyBoundary(called, true, 0);
		verifyBoundary(called, false, 1);
		verifyBoundary(called, true, 1);
		verifyReportedBounds();
		GeneCaller.ncrnaFamilies.clear();

		options.windowPad=123;
		options.topN=1;
		options.seedMinHits=1;
		final NcrnaFamily overridden=NcrnaFamilyLoader.load("glycine", options, new FamilyResourceResolver());
		check(overridden.windowPad==123 && overridden.indexTopN==1
			&& overridden.seedMinHits==1 && !overridden.seedDistinct,
			"Explicit experimental controls must reach the direct loader");
		final File wrong=new File(root, "wrong_prefix");
		writeFixture(wrong, "other_");
		CallGenes.parseGlycineFlag("glycineroot", wrong.getPath());
		expectMissing();
	}

	/** Runs the caller's scavenger and control handoff, not a standalone aligner replay. */
	private static void verifyQuantum(final NcrnaFamily family, final boolean padded){
		final NcrnaScavenger caller=new NcrnaScavenger(family.library, family.models, family.modelNames,
			family.kmerSet, family.kLong, family.minLen, family.windowPad,
			family.indexK, family.indexTopN, family.adaptive, family.adaptFloor,
			family.adaptTopFrac, family.adaptQFrac, family.fixedMinHits,
			family.scoreA, family.scoreB, family.idPass, family.idBorderline);
		GeneCaller.applyNcrnaFamilyControls(caller, family);
		final WindowWidths widths=new WindowWidths();
		caller.setWorkloadSink(widths);
		final byte[] bases=new byte[padded ? 224 : 64];
		Arrays.fill(bases, (byte)'N');
		final int start=padded ? 80 : 0;
		System.arraycopy(family.library[0], 0, bases, start, 64);
		final QuantumAligner quantum=new QuantumAligner();
		final ScrabbleAligner scrabble=new ScrabbleAligner();
		final long quantumBefore=quantum.loops(), scrabbleBefore=scrabble.loops();
		final ArrayList<Orf> calls=caller.scavenge(padded ? "long" : "short", bases, 0, new ArrayList<int[]>());
		check(calls.size()==1 && calls.get(0).start==start && calls.get(0).stop==start+63,
			"The literal64-base gene must survive seed/window/verification with its exact coordinates");
		check(widths.primary>0 && (padded ? widths.maximum>120 : widths.maximum<=120),
			"Both sides of the former120-base aligner boundary must actually be scheduled");
		final long quantumCells=quantum.loops()-quantumBefore;
		final long scrabbleCells=scrabble.loops()-scrabbleBefore;
		check(quantumCells>0 && scrabbleCells==0,
			"Glycine verification and endpoints require Quantum only: padded="+padded
			+" quantumCells="+quantumCells+" scrabbleCells="+scrabbleCells);
		System.err.println("GLYCINE_QUANTUM_PATH_PASS padded="+padded+" maximumWindow="+widths.maximum
			+" quantumCells="+quantumCells+" scrabbleCells="+scrabbleCells);
	}

	/** A length-preferring scorer must extend both ends by exactly one base.
	 * Uses the real loaded-family controls and seeded Quantum caller. The null-net
	 * control keeps the same raw endpoints, so eligibility alone cannot enable a net.
	 */
	private static void verifyBoundary(final NcrnaFamily family, final boolean enabled, final int strand){
		final LengthNet net=new LengthNet();
		final CellNet selected=enabled ? net : null;
		final TrnaBoundaryFeatures.NinemerTable table=(bases, position)->1;
		final int[] offsets={-1, 0, 1};
		final NcrnaScavenger caller=new NcrnaScavenger(family.library, family.models, family.modelNames,
			family.kmerSet, family.kLong, family.minLen, family.windowPad,
			family.indexK, family.indexTopN, family.adaptive, family.adaptFloor,
			family.adaptTopFrac, family.adaptQFrac, family.fixedMinHits,
			family.scoreA, family.scoreB, family.idPass, family.idBorderline,
			selected, selected, table, table, 7, 0, 7, 0, 64f, offsets, offsets);
		GeneCaller.applyNcrnaFamilyControls(caller, family);
		final byte[] bases=new byte[224];
		Arrays.fill(bases, (byte)'N');
		System.arraycopy(family.library[0], 0, bases, 80, 64);
		final ArrayList<int[]> claims=new ArrayList<int[]>();
		final ArrayList<Orf> calls=caller.scavenge("boundary", bases, strand, claims);
		check(calls.size()==1 && calls.get(0).start==80 && calls.get(0).stop==143,
			"Both NN states must return the raw accepted span to candidate selection");
		check(claims.size()==1 && claims.get(0)[0]==80 && claims.get(0)[1]==143,
			"Neighbor discovery must consume the same raw claim with or without endpoint movement");
		final float rawScore=calls.get(0).orfScore;
		if(strand==1){calls.get(0).flip();}
		NcrnaReportedOrf.applyTo(calls);
		final int expectedStart=enabled ? 79 : 80;
		final int expectedStop=enabled ? 144 : 143;
		check(calls.size()==1 && calls.get(0).start==expectedStart && calls.get(0).stop==expectedStop,
			"The loaded glycine policy must reach reported endpoint movement: enabled="+enabled
			+" expected="+expectedStart+"-"+expectedStop+" calls="+calls);
		NcrnaReportedOrf.applyTo(calls);
		check(calls.get(0).start==expectedStart && calls.get(0).stop==expectedStop
			&& calls.get(0).orfScore==rawScore,
			"Repeated reporting handoff must preserve coordinates and raw selection score");
		check(enabled ? net.calls>0 : net.calls==0,
			"Only an installed scorer may run on the retained Quantum endpoints: enabled="+enabled+" scores="+net.calls);
		System.err.println("GLYCINE_BOUNDARY_INVOCATION_PASS enabled="+enabled
			+" strand="+strand+" scores="+net.calls);
	}

	/** Asymmetric literals expose flip errors and demonstrate raw overlap decisions. */
	private static void verifyReportedBounds(){
		final byte[] bases=new byte[100];
		Arrays.fill(bases, (byte)'A');
		for(int strand=0; strand<2; strand++){
			final NcrnaReportedOrf first=new NcrnaReportedOrf("pair", 20, 40, strand, bases);
			first.start=19;
			first.stop=42;
			first.captureRefinement(20, 40, 17, 80);
			final Orf next=new Orf("pair", 41, 60, strand, 0, bases, false, Orf.RNA);
			check(first.start==20 && first.stop==40 && next.isValidPrev(first, 0),
				"A refined overlap must not enter the final-path predecessor decision");
			if(strand==1){first.flip();}
			first.applyReportedBounds();
			check(first.start==(strand==0 ? 19 : 57) && first.stop==(strand==0 ? 42 : 80),
				"Reported bounds must use current flip parity, not biological strand alone");
			first.flip();
			first.applyReportedBounds();
			check(first.start==(strand==0 ? 57 : 19) && first.stop==(strand==0 ? 80 : 42),
				"Reporting remains correct after a later coordinate flip");
		}
		final NcrnaReportedOrf shortProposal=new NcrnaReportedOrf("short", 20, 40, 0, bases);
		shortProposal.start=30;
		shortProposal.stop=40;
		shortProposal.captureRefinement(20, 40, 17, 80);
		shortProposal.applyReportedBounds();
		check(shortProposal.start==20 && shortProposal.stop==40,
			"An undersized refinement must retain the raw call rather than remove detection");
		final NcrnaReportedOrf contraction=new NcrnaReportedOrf("contract", 20, 40, 0, bases);
		contraction.start=22;
		contraction.stop=39;
		contraction.captureRefinement(20, 40, 17, 80);
		check(contraction.start==20 && contraction.stop==40,
			"Contraction must not release seed positions into later discovery");
		contraction.applyReportedBounds();
		check(contraction.start==22 && contraction.stop==39,
			"A valid contraction remains available for the final report");
	}

	/** Single-thread fixture: copies share the counter so actual scoring is observable. */
	private static final class LengthNet extends CellNet{
		@Override
		public CellNet copy(final boolean copyWeights){return this;}
		@Override
		public CellNet applyInput(final float[] input){
			CallGenesGlycineTest.check(input.length==10 && Float.isFinite(input[8]),
				"The shared v1 scorer supplies its finite length-ratio feature at index8");
			lengthRatio=input[8];
			return this;
		}
		@Override
		public float feedForward(){
			calls++;
			return lengthRatio;
		}
		private float lengthRatio;
		private int calls;
	}

	/** Records actual scheduled primary windows without changing calling decisions. */
	private static final class WindowWidths implements NcrnaWorkloadInstrumentSink{
		@Override
		public void seedHits(String name, int strand, int[] positions){}
		@Override
		public void scheduledWindow(String name, int strand, int pass, int start, int stop){
			check(stop>=start, "A scheduled inclusive window must have valid bounds");
			if(pass==1){
				primary++;
				maximum=Math.max(maximum, stop-start+1);
			}
		}
		int primary, maximum;
	}

	/** A selected missing or misattributed family must fail without publishing partial state. */
	private static void expectMissing(){
		boolean failed=false;
		try{CallGenes.loadGlycineResources();}
		catch(IllegalArgumentException expected){failed=true;}
		check(failed && GeneCaller.ncrnaFamilies.isEmpty(), "Invalid selected resources must fail before family publication");
	}

	/** Produces exact64/68/72-base paired records using the native model codec. */
	private static void writeFixture(final File root, final String prefix){
		check(root.mkdir(), "Each fixture family needs a fresh directory: "+root);
		final String base="GATCCGTACGAGTCCATGCTAGACCTGATCGAGCTTAGGCATACGTCATGGTACCGATCAGTGC";
		final String[] sequences={base, base+"AACC", base+"AACCGGTT"};
		check(base.length()==64, "The literal fixture must preserve its declared 64-base gene");
		final ArrayList<BaseGraph> models=new ArrayList<BaseGraph>();
		final ByteBuilder bb=new ByteBuilder();
		try(Output fasta=new Output(new File(root, "consensus.fa").getPath())){
			for(int i=0; i<sequences.length; i++){
				final String name=prefix+i;
				final byte[] bases=sequences[i].getBytes(StandardCharsets.US_ASCII);
				bb.clear();
				fasta.writer.print(bb.append('>').append(name).nl().append(bases).nl());
				final BaseGraph graph=new BaseGraph(name, bases, null, 0, 0);
				BaseGraphHelper.initForScoring(graph);
				models.add(graph);
			}
		}
		NucleotideHbmCodec.writeModels(models, new File(root, "models.hbmc").getPath(), false);
		try(Output seeds=new Output(new File(root, "seeds.fa").getPath());
			Output manifest=new Output(new File(root, "manifest.tsv").getPath())){
			seeds.writer.print(">first\n"+base.substring(0, 17)+"\n>second\n"+base.substring(17, 34)+"\n");
			manifest.writer.print("#ncrna_family\t1\tglycine\nrole\tmodel\tend\tpath\tformat\tbytes\tsha80\n"
				+"consensus\t.\t.\tconsensus.fa\tfasta\t-\t-\nmodels\t.\t.\tmodels.hbmc\thbm\t-\t-\n"
				+"seeds\t.\t.\tseeds.fa\tfasta\t-\t-\n");
		}
	}

	private static void check(final boolean condition, final String why){
		if(!condition){throw new AssertionError(why);}
	}
}
