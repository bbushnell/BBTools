package test.prok;

import java.io.File;
import java.util.Arrays;
import java.util.Locale;

import ml.CellNet;
import ml.CellNetParser;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.FamilyResourceResolver;
import prok.NcrnaEndpointModelOffsets;
import prok.NcrnaFamily;
import prok.NcrnaFamilyLoader;
import prok.NcrnaFixedEndpointTrials;
import prok.NcrnaPositionalBoundaryFeatures;
import prok.NcrnaPositionalBoundaryRefiner;
import prok.NcrnaPositionalBoundaryRefiner.GatePolicy;
import prok.NcrnaReportedOrf;
import prok.NcrnaTextIO.Output;
import prok.NcrnaWeightedEndpointTable;
import prok.RrnaPositionalKmerTable.Table;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Existing real packages on bounded synthetic spans: no fitting or genome scan.
 * Independent35-input arithmetic checks the caller against its exported nets;
 * literal nets distinguish whole-move rejection from a residual-only gate.
 * @author Brian, Ganyu
 */
public final class NcrnaCenteredGateServingTest{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaCenteredGateServingTest.class, false);
		try{
			new NcrnaCenteredGateServingTest(pp.args).run();
			pp.outstream.println("NCRNA_CENTERED_GATE_SERVING_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private NcrnaCenteredGateServingTest(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("familyroot")){root=b;}
			else if(a.equals("family")){name=b;}
			else if(a.equals("cutoff5")){setCutoff(0, b);}
			else if(a.equals("cutoff3")){setCutoff(1, b);}
			else{check(parser.parse(arg, a, b), "Unknown centered serving fixture option: "+arg);}
		}
		check(root!=null && ("cobalamin_rf00174".equals(name) || "cobalamin_rf01482".equals(name))
			&& cutoffSet[0] && cutoffSet[1] && parser.out1!=null
			&& parser.in1==null && parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append
			&& parser.maxReads<0 && parser.samplerate==1, "Bind a real centered package, its exact cutoffs and out=fresh_fixture_directory");
		out=new File(parser.out1);
		check(out.mkdir(), "Finite serving outputs must be fresh");
	}

	private void run(){
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
		options.resources.root=root;
		options.endpoints=false;
		family=NcrnaFamilyLoader.load(name, options, new FamilyResourceResolver());
		check(family.boundaryOnRawEndpoints && family.deferBoundaryReport, "The loaded experimental package must reach both actual caller controls");
		final NcrnaPositionalBoundaryRefiner.Resources resources=(NcrnaPositionalBoundaryRefiner.Resources)family.boundaryRefinerFactory;
		check(resources.centeredGate() && resources.fullErrors() && resources.radius()==10, "Only the explicit centered35 serving route is under test");
		final String member=name.substring("cobalamin_".length()).toUpperCase(Locale.ROOT);
		final File payload=new File(root, member);
		for(int end=0; end<2; end++){
			for(int probe=0; probe<3; probe++){
				weighted[end][probe]=NcrnaWeightedEndpointTable.loadScores(new File(payload, (end==0 ? "five" : "three")+probe+".tsv").getPath(), member);
				geometry[end][probe]=weighted[end][probe].geometry;
			}
			nets[end]=CellNetParser.load(new File(payload, end==0 ? "five.bbnet" : "three.bbnet").getPath());
			gates[end]=CellNetParser.load(new File(payload, end==0 ? "gate5.bbnet" : "gate3.bbnet").getPath());
			check(nets[end].numInputs()==69 && nets[end].numOutputs()==22 && gates[end].numInputs()==35 && gates[end].numOutputs()==1,
				"Real first-stage and whole-move exports must retain their distinct widths");
		}
		features=new NcrnaPositionalBoundaryFeatures(geometry[0], geometry[1], weighted[0], weighted[1], false);
		final NcrnaPositionalBoundaryRefiner worker=(NcrnaPositionalBoundaryRefiner)resources.create(family);
		for(final int model:new int[]{0, family.library.length-1}){
			final byte[] bases=bases(model);
			for(final int shift:new int[]{-3, 0, 4}){
				for(int fixed=0; fixed<4; fixed++){
					for(int strand=0; strand<2; strand++){
						parity(worker, resources.modelOffsets(), bases, model, PAD+shift, bases.length-PAD-1-shift, fixed, strand);
					}
				}
			}
		}
		check(observedGate[0]>0 && observedGate[1]>0, "Both real gate exports must actually be exercised; successful loading alone is insufficient");
		literal(member);
		try(Output output=new Output(new File(out, "summary.tsv").getPath())){
			output.writer.println("family\tscalar_cases\tsimd_cases\tstart_gate_checks\tstop_gate_checks\tqualification");
			output.writer.print(new ByteBuilder().append(name).tab().append(scalarCases).tab().append(simdCases)
				.tab().append(observedGate[0]).tab().append(observedGate[1]).append("\tFINITE_SYNTHETIC_SERVING_PARITY_NOT_ACCURACY\n"));
		}
	}

	/** The oracle builds35inputs literally, without calling the production gate feature helper. */
	private void parity(final NcrnaPositionalBoundaryRefiner worker, final NcrnaEndpointModelOffsets offsets,
			final byte[] bases, final int model, final int start, final int stop, final int fixed, final int strand){
		final boolean saved=Shared.SIMD;
		int scalarStart=-1;
		int scalarStop=-1;
		try{
			for(int mode=0; mode<2; mode++){
				if(mode==1 && (!saved || !simd.Vector.simd128)){continue;}
				Shared.SIMD=mode==1;
				final boolean vector=mode==1;
				final long left=(long)start-((fixed&1)==0 ? offsets.median(0, model) : 0);
				final long right=(long)stop-((fixed&2)==0 ? offsets.median(1, model) : 0);
				final boolean valid=left>=0 && right<bases.length && right-left+1>=family.minLen && right-left+1<=family.maxLen;
				final int centerStart=valid ? (int)left : start;
				final int centerStop=valid ? (int)right : stop;
				final int usable=valid ? features.fill(bases, centerStart, centerStop, family.library[model], family.models[model], plateau(model), .37f, fixed) : 0;
				final float[][] expectedInputs=new float[2][];
				final float[][] expectedScores=new float[2][];
				final long[] proposals={start, stop};
				final long[] chosen={start, stop};
				int expectedMask=0;
				for(int end=0; end<2; end++){
					if((usable&(1<<end))==0){continue;}
					nets[end].applyInput(features.values[end]);
					nets[end].feedForward();
					final float[] scores=nets[end].getOutput();
					expectedScores[end]=scores;
					final int adjustment=adjustment(scores);
					final int raw=end==0 ? start : stop;
					final int center=end==0 ? centerStart : centerStop;
					proposals[end]=(long)center+adjustment;
					if(proposals[end]==raw){continue;}
					expectedMask|=1<<end;
					final float[] input=new float[35];
					System.arraycopy(scores, 0, input, 0, 22);
					input[22]=features.values[end][63];
					input[23]=((long)centerStop-centerStart+1+(end==0 ? -adjustment : adjustment))/(float)family.library[model].length;
					input[24]=features.values[end][64];
					input[25]=gc(bases, centerStart, centerStop);
					input[26]=gc(family.library[model], 0, family.library[model].length-1);
					input[27]=input[25]-input[26];
					input[28]=2*features.values[end][65]-1;
					System.arraycopy(features.values[end], 66, input, 29, 3);
					input[32]=Math.abs((long)center-raw);
					input[33]=(long)center-raw;
					input[34]=adjustment==0 ? 0 : 1;
					expectedInputs[end]=input;
				}
				final int[] observed={0, 0};
				worker.setScoreObserver((end, input, scores, selected)->{
					check(Arrays.equals(input, features.values[end]), "Both first-stage inputs must describe the same guarded center");
					for(int i=0; i<22; i++){check(same(scores[i], expectedScores[end][i], vector), "First-stage output mismatch");}
					check(selected==adjustment(scores)+10, "Centered first-stage observer reports its neural class before the whole-move gate");
					observed[0]|=1<<end;
					Arrays.fill(input, Float.NaN);
					Arrays.fill(scores, Float.NaN);
				});
				worker.setWholeMoveObserver((end, input, score, proposal, accepted)->{
					check(proposal==proposals[end] && input.length==35, "Whole proposal must retain its original raw anchor");
					for(int i=0; i<35; i++){check(same(input[i], expectedInputs[end][i], vector && i<22), "Independent35-input oracle differs at field "+i);}
					gates[end].applyInput(input);
					gates[end].feedForward();
					final float direct=gates[end].getOutput(0);
					final boolean take=policies[end]==GatePolicy.ALWAYS || policies[end]==GatePolicy.THRESHOLD && direct>=cutoffs[end];
					check(same(score, direct, vector) && accepted==take, "Same-input exported gate score and explicit policy must agree");
					chosen[end]=accepted ? proposal : end==0 ? start : stop;
					check((observed[1]&(1<<end))==0, "Each distinct proposal invokes its gate once");
					observed[1]|=1<<end;
					observedGate[end]++;
					Arrays.fill(input, Float.NaN);
				});
				final NcrnaReportedOrf call=call(bases, model, start, stop, strand);
				worker.refine(call, bases, model, .37f, .8f, true, fixed, family.minLen, family.maxLen);
				check(observed[0]==usable && observed[1]==expectedMask, "Missing/fixed/invalid-center and raw-equal ends must not invent gate work");
				if(chosen[0]<0 || chosen[1]>=bases.length || chosen[1]-chosen[0]+1<family.minLen || chosen[1]-chosen[0]+1>family.maxLen){
					chosen[0]=start;
					chosen[1]=stop;
				}
				check(call.start==chosen[0] && call.stop==chosen[1] && call.orfScore==42, "Joint invalid or rejected movement must restore the path anchor");
				if(vector){check(call.start==scalarStart && call.stop==scalarStop, "Scalar/SIMD whole-move decisions must agree");simdCases++;}
				else{scalarStart=call.start;scalarStop=call.stop;scalarCases++;}
				call.captureRefinement(start, stop, family.minLen, family.maxLen);
				check(call.start==start && call.stop==stop, "Deferred capture must leave ordinary selection coordinates intact");
				if(strand==1){call.flip();}
				call.applyReportedBounds();
				check(call.start==(strand==0 ? chosen[0] : bases.length-1-chosen[1])
					&& call.stop==(strand==0 ? chosen[1] : bases.length-1-chosen[0]), "Deferred whole-move report must honor genomic flip parity");
			}
		}finally{Shared.SIMD=saved;worker.setScoreObserver(null);worker.setWholeMoveObserver(null);}
	}

	/** Literal +5/-3 training errors give applied median movements -5/+3. */
	private void literal(final String member){
		final File file=new File(out, "literal_offsets.tsv");
		try(Output output=new Output(file.getPath())){
			output.writer.println(NcrnaFixedEndpointTrials.OFFSET_HEADER);
			final ByteBuilder b=new ByteBuilder();
			for(int i=0; i<family.modelNames.length; i++){
				output.writer.print(b.clear().append(member).tab().append(Tools.trimToWhitespace(family.modelNames[i]))
					.append(i==0 ? "\t10\t10\t5\t-3\n" : "\t0\t10\t0\t0\n"));
			}
		}
		final NcrnaEndpointModelOffsets offsets=NcrnaEndpointModelOffsets.load(file.getPath(), family, member, 10);
		fails(()->new NcrnaPositionalBoundaryRefiner.Resources(family, geometry[0], geometry[1], weighted[0], weighted[1],
			new FixedNet(10, false), new FixedNet(10, false), 0, 0, false,
			new FixedGate(new int[1], 32, .1f), new FixedGate(new int[1], 32, .1f), .1f, .1f, offsets));
		final byte[] bases=bases(0);
		final int start=PAD;
		final int stop=bases.length-PAD-1;
		for(int policyIndex=0; policyIndex<6; policyIndex++){
			final float threshold=policyIndex==0 ? Math.nextDown(.1f) : policyIndex==1 ? .1f : policyIndex==2 ? Math.nextUp(.1f) : 0;
			final GatePolicy policy=policyIndex==3 ? GatePolicy.ALWAYS : policyIndex==4 ? GatePolicy.NEVER : GatePolicy.THRESHOLD;
			for(final boolean tied:new boolean[]{false, true}){
				final int[] count={0};
				final NcrnaPositionalBoundaryRefiner worker=literalWorker(offsets, tied ? 12 : 10, tied ? 8 : 10, tied, threshold, policy, count);
				final NcrnaReportedOrf call=call(bases, 0, start, stop, 0);
				worker.refine(call, bases, 0, .37f, .8f, true, 0, 16, bases.length);
				final boolean accepted=policy==GatePolicy.ALWAYS || policy==GatePolicy.THRESHOLD && .1f>=threshold;
				check(count[0]==2 && call.start==(accepted ? start-5 : start) && call.stop==(accepted ? stop+3 : stop),
					"Current-best/tied neural scores still gate center-only moves; rejection/never must restore raw");
			}
		}
		final int[] count={0};
		final NcrnaPositionalBoundaryRefiner inward=literalWorker(offsets, 20, 0, false, .1f, count);
		final NcrnaReportedOrf joint=call(bases, 0, start, stop, 0);
		inward.refine(joint, bases, 0, .37f, .8f, true, 0, stop-start+1-10, bases.length);
		check(count[0]==2 && joint.start==start && joint.stop==stop, "Two separately legal moves that jointly shorten too far restore BOTH raw ends");
		count[0]=0;
		final NcrnaPositionalBoundaryRefiner centered=literalWorker(offsets, 10, 10, false, .1f, count);
		final NcrnaReportedOrf invalid=call(bases, 0, 2, stop, 0);
		centered.refine(invalid, bases, 0, .37f, .8f, true, 0, 16, bases.length);
		check(count[0]==0 && invalid.start==2 && invalid.stop==stop, "Invalid center performs no gate inference and keeps raw");
		final NcrnaReportedOrf fixed=call(bases, 0, 0, stop, 0);
		centered.refine(fixed, bases, 0, .37f, .8f, true, 1, 16, bases.length);
		check(count[0]==1 && fixed.start==0 && fixed.stop==stop+3, "Fixed start bypasses its gate while the free stop may move");
		count[0]=0;
		bases[start-5]='N';
		final NcrnaReportedOrf missing=call(bases, 0, start, stop, 0);
		centered.refine(missing, bases, 0, .37f, .8f, true, 0, 16, bases.length);
		check(count[0]==1 && missing.start==start && missing.stop==stop+3, "Unavailable centered start features restore RAW start, not center");
		final byte[] clean=bases(0);
		final NcrnaPositionalBoundaryRefiner nonfinite=(NcrnaPositionalBoundaryRefiner)new NcrnaPositionalBoundaryRefiner.Resources(
			family, geometry[0], geometry[1], weighted[0], weighted[1], new FixedNet(10, false), new FixedNet(10, false),
			0, 0, false, new FixedGate(new int[1], 35, Float.NaN), new FixedGate(new int[1]),
			0, 0, offsets, GatePolicy.ALWAYS, GatePolicy.ALWAYS).create(family);
		fails(()->nonfinite.refine(call(clean, 0, start, stop, 0), clean, 0, .37f, .8f, true, 0, 16, clean.length));
		if(family.library.length>1){
			count[0]=0;
			final int model=family.library.length-1;
			final byte[] unsupported=bases(model);
			final NcrnaReportedOrf unchanged=call(unsupported, model, PAD, unsupported.length-PAD-1, 0);
			centered.refine(unchanged, unsupported, model, .37f, .8f, true, 0, 16, unsupported.length);
			check(count[0]==0 && unchanged.start==PAD && unchanged.stop==unsupported.length-PAD-1,
				"Unsupported offsets plus a zero neural adjustment are raw-equal holds without gate work");
		}
	}

	private NcrnaPositionalBoundaryRefiner literalWorker(final NcrnaEndpointModelOffsets offsets, final int five,
			final int three, final boolean tied, final float cutoff, final int[] count){
		return literalWorker(offsets, five, three, tied, cutoff, GatePolicy.THRESHOLD, count);
	}
	private NcrnaPositionalBoundaryRefiner literalWorker(final NcrnaEndpointModelOffsets offsets, final int five,
			final int three, final boolean tied, final float cutoff, final GatePolicy policy, final int[] count){
		return (NcrnaPositionalBoundaryRefiner)new NcrnaPositionalBoundaryRefiner.Resources(family, geometry[0], geometry[1],
			weighted[0], weighted[1], new FixedNet(five, tied), new FixedNet(three, tied), 0, 0, false,
			new FixedGate(count), new FixedGate(count), cutoff, cutoff, offsets, policy, policy).create(family);
	}
	private byte[] bases(final int model){
		final byte[] result=new byte[family.library[model].length+2*PAD];
		final byte[] alphabet={'A','C','G','T'};
		for(int i=0; i<result.length; i++){result[i]=alphabet[(i*17+i/7)&3];}
		for(int i=0; i<family.library[model].length; i++){
			final byte b=family.library[model][i];
			result[PAD+i]=b=='A' || b=='C' || b=='G' || b=='T' ? b : (byte)'A';
		}
		return result;
	}
	private float plateau(final int model){
		int maximum=0;
		for(int i=0; i<family.models[model].ref.length; i++){
			maximum=Math.max(maximum, family.models[model].ref[i].countSum+family.models[model].del[i].countSum);
		}
		check(maximum>0, "Real paired HBM must have positive coverage");
		return maximum;
	}
	private NcrnaReportedOrf call(final byte[] bases, final int model, final int start, final int stop, final int strand){
		final NcrnaReportedOrf call=new NcrnaReportedOrf("fixture", start, stop, strand, bases);
		call.ncrnaFamily=family.name;
		call.trnaModel=family.modelNames[model];
		call.orfScore=42;
		return call;
	}
	private void setCutoff(final int end, final String value){
		check(!cutoffSet[end] && value!=null, "Each fixture policy must be supplied once");
		policies[end]="never".equals(value) ? GatePolicy.NEVER : "always".equals(value) ? GatePolicy.ALWAYS : GatePolicy.THRESHOLD;
		cutoffs[end]=policies[end]==GatePolicy.THRESHOLD ? Float.parseFloat(value) : 0;
		check(Float.isFinite(cutoffs[end]), "Fixture numeric cutoffs must be finite");
		cutoffSet[end]=true;
	}
	private static int adjustment(final float[] scores){
		int best=0;
		for(int i=1; i<21; i++){if(scores[i]>scores[best]){best=i;}}
		for(int i=0; i<21; i++){if(i!=best && scores[i]==scores[best]){return 0;}}
		return best-10;
	}
	private static float gc(final byte[] bases, final int start, final int stop){
		int recognized=0;
		int gc=0;
		for(int i=start; i<=stop; i++){
			final byte b=bases[i];
			if(b=='A' || b=='C' || b=='G' || b=='T'){recognized++;}
			if(b=='C' || b=='G'){gc++;}
		}
		check(recognized>0, "Finite GC oracle needs recognized bases");
		return gc/(float)recognized;
	}
	private static boolean same(final float a, final float b, final boolean vector){return vector ? Math.abs((double)a-b)<=2.4e-7 : Float.floatToRawIntBits(a)==Float.floatToRawIntBits(b);}
	private static void check(final boolean value, final String why){if(!value){throw new AssertionError(why);}}
	private static void fails(final Runnable action){
		try{action.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Invalid gate dimensions or nonfinite scores must fail loudly");
	}
	private static final class FixedNet extends CellNet{
		FixedNet(final int peak_, final boolean tied_){
			peak=peak_;
			tied=tied_;
			scores[peak]=1;
			scores[21]=2;
			if(tied){scores[7]=1;}
		}
		@Override public int numInputs(){return 69;}
		@Override public int numOutputs(){return 22;}
		@Override public CellNet copy(final boolean deltas){return new FixedNet(peak, tied);}
		@Override public CellNet applyInput(final float[] input){
			NcrnaCenteredGateServingTest.check(input.length==69, "First-stage fixture width");
			return this;
		}
		@Override public float feedForward(){return 1;}
		@Override public float getOutput(final int index){return scores[index];}
		private final int peak;
		private final boolean tied;
		private final float[] scores=new float[22];
	}
	private static final class FixedGate extends CellNet{
		FixedGate(final int[] count_){this(count_, 35, .1f);}
		FixedGate(final int[] count_, final int inputs_, final float value_){
			count=count_;
			inputs=inputs_;
			value=value_;
		}
		@Override public int numInputs(){return inputs;}
		@Override public int numOutputs(){return 1;}
		@Override public CellNet copy(final boolean deltas){return new FixedGate(count, inputs, value);}
		@Override public CellNet applyInput(final float[] input){
			NcrnaCenteredGateServingTest.check(input.length==inputs, "Whole-move fixture width");
			return this;
		}
		@Override public float feedForward(){
			count[0]++;
			return value;
		}
		@Override public float getOutput(final int index){
			NcrnaCenteredGateServingTest.check(index==0, "One gate score");
			return value;
		}
		private final int[] count;
		private final int inputs;
		private final float value;
	}
	private String root, name;
	private File out;
	private NcrnaFamily family;
	private NcrnaPositionalBoundaryFeatures features;
	private final Table[][] geometry=new Table[2][3];
	private final NcrnaWeightedEndpointTable[][] weighted=new NcrnaWeightedEndpointTable[2][3];
	private final CellNet[] nets=new CellNet[2], gates=new CellNet[2];
	private final float[] cutoffs=new float[2];
	private final boolean[] cutoffSet=new boolean[2];
	private final GatePolicy[] policies=new GatePolicy[2];
	private final long[] observedGate=new long[2];
	private long scalarCases,simdCases;
	private static final int PAD=80;
}
