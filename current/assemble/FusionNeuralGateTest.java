package assemble;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import dna.AminoAcid;
import fileIO.ByteFile;
import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;
import structures.ByteBuilder;
import ukmer.Kmer;

/** Deterministic v2 feature, inference, veto-order and parser fixtures. @author Fischl */
final class FusionNeuralGateTest {

	/** Tests the actual candidate model without treating a classifier score as assembly truth. */
	static void run(final String model){
		final boolean dense=CellNet.DENSE;
		final FusionNeuralGate gate=new FusionNeuralGate(model, .5f);
		check(CellNet.DENSE==dense, "Loading the fusion model changed global DENSE mode.");
		final CellNet independent=CellNetParser.loadInferenceFromLines(ByteFile.toLines(model));
		independent.simdFF=true;
		final TadpoleGraph.Counts counts=new TadpoleGraph.Counts(){
			@Override
			public int count(final Kmer key){return (int)((key.xor()&Long.MAX_VALUE)%71);}
		};
		int cases=0;
		float maxScoreDifference=0;
		for(int k : new int[]{3, 24, 35, 64}){
			gate.beginPhase(k, counts);
			final FusionJoinDiagnostic offline=new FusionJoinDiagnostic(k, counts);
			for(int orientation=0; orientation<4; orientation++){
				for(int trim : new int[]{0, 7}){
					for(int flank : new int[]{5, 31, 200, 350}){
						final int overlap=k+5;
						final byte[] sequence=randomBases(2*flank+overlap, 196+cases);
						final byte[] a=Arrays.copyOf(sequence, flank+overlap+trim);
						final byte[] b=new byte[trim+overlap+flank];
						Arrays.fill(b, 0, trim, (byte)'N');
						System.arraycopy(sequence, flank, b, trim, overlap+flank);
						Arrays.fill(a, flank+overlap, a.length, (byte)'N');
						final boolean sr=(orientation&1)!=0, dr=(orientation&2)!=0;
						final Contig source=contig(a, sr, 0), dest=contig(b, dr, 1);
						final byte[] savedA=Arrays.copyOfRange(a, Math.max(0, a.length-trim-overlap-200), a.length);
						final byte[] savedB=Arrays.copyOf(b, Math.min(b.length, trim+overlap+200));
						final FusionJoinDiagnostic.Join join=new FusionJoinDiagnostic.Join(1, k, 0, 1,
								overlap, trim, trim, savedA, savedB);
						final float[] core=offline.measure(join.product, join.product.length,
								join.left.length-overlap, overlap, trim, trim, new ByteBuilder());
						final ByteBuilder serialized=new ByteBuilder();
						for(float value : core){
							if(serialized.length()>0){serialized.tab();}
							serialized.append(value, 9);
						}
						new FusionTipEntropy().append(serialized, savedA, savedB);
						final LineParser1 parser=new LineParser1('\t');
						parser.set(serialized.toBytes());
						check(parser.terms()==109, "The offline v2 row changed width.");
						final float[] expected=new float[109];
						final float[] actual=gate.features(source, sr, trim, dest, dr, trim, overlap);
						for(int i=0; i<expected.length; i++){
							expected[i]=parser.parseFloat(i);
							check(Float.floatToIntBits(expected[i])==Float.floatToIntBits(actual[i]),
									"Training/live feature mismatch at input "+i+", case "+cases);
						}
						final float score=independent.applyInput(expected).feedForward();
						final float liveScore=gate.score(actual);
						final float difference=Math.abs(score-liveScore);
						maxScoreDifference=Math.max(maxScoreDifference, difference);
						// Vector reductions may round differently across interpreter/JIT transitions
						// (simd.Vector.dotRows documents this). Inputs above remain bit-exact;
						// testVeto separately verifies exact cutoff equality with fixed scores.
						check(Float.isFinite(score) && difference<=0.000001f,
								"Training/live probability mismatch: case="+cases+", k="+k+
								", offline="+score+", live="+liveScore+", difference="+difference);
						check(Arrays.equals(source.bases, contig(a, sr, 0).bases), "Gate mutated the source.");
						check(Arrays.equals(dest.bases, contig(b, dr, 1).bases), "Gate mutated the destination.");
						cases++;
					}
				}
			}
			gate.endPhase();
		}
		testVeto(counts);
		testParser();
		System.err.println("FUSION_NEURAL_TEST_PASS feature_score_cases="+cases+
				" maxScoreDifference="+maxScoreDifference+" veto_cycle_parser");
	}

	/** Equality accepts, a stricter gate rejects, and cycles never reach the model. */
	private static void testVeto(final TadpoleGraph.Counts counts){
		final FixedNet net=new FixedNet(.5f);
		final FusionNeuralGate accept=new FusionNeuralGate(net, .5f);
		final FusionNeuralGate reject=new FusionNeuralGate(net, .6f);
		for(FusionNeuralGate gate : new FusionNeuralGate[]{accept, reject}){
			gate.beginPhase(3, counts);
			final ArrayList<Contig> pair=contigs("AGGTTAAC", "AACGGCCT");
			final CrossKTipOverlapper overlapper=new CrossKTipOverlapper(pair, 3, 3);
			overlapper.neural=gate;
			check(overlapper.addEdges()==(gate==accept ? 1 : 0), "Inclusive cutoff or rejection changed.");
			check(gate.evaluations==1, "A frozen pair was scored more than once or not scored.");
			final ArrayList<Contig> cycle=contigs("GGATTAAC", "AACGGCCT", "CCTCCGGA");
			final CrossKTipOverlapper cyclic=new CrossKTipOverlapper(cycle, 3, 3);
			cyclic.neural=gate;
			check(cyclic.addEdges()==0 && gate.evaluations==1, "Neural filtering modified cycle rejection.");
			gate.endPhase();
		}
		net.value=Float.NaN;
		boolean failed=false;
		try{accept.score(new float[109]);}catch(IllegalStateException expected){failed=true;}
		check(failed, "Nonfinite inference must fail loudly.");
	}

	/** NN flags route to multi-K only and never leak into underlying table arguments. */
	private static void testParser(){
		final TadpoleMulti.Config config=new TadpoleMulti.Config(new String[]{"k=95,35", "out=out.fa",
				"fusenet=model.bbnet", "fusencutoff=0.667098"});
		check(TadpoleMulti.hasMultipleK(new String[]{"fusenet=model.bbnet"}), "Network flag did not route to multi-K.");
		check(config.fuseCutoff==.667098f && config.fuseNet.equals("model.bbnet"), "Explicit gate options lost.");
		final String[] args=new TadpoleMulti(config).makeFusionSupportArgs(35);
		check(!Arrays.asList(args).contains("minprob=0"), "Neural gate must not change a shared path guard's quality filtering.");
		for(String arg : args){check(!arg.startsWith("fusen"), "Neural policy leaked into table parsing.");}
		for(String[] bad : new String[][]{{"fusenet=model"}, {"fusencutoff=.5"},
				{"fusenet=model", "fusencutoff=NaN"}, {"fusenet=model", "fusencutoff=1.1"}}){
			final ArrayList<String> all=new ArrayList<String>(Arrays.asList("k=95,35", "out=out.fa"));
			all.addAll(Arrays.asList(bad));
			boolean failed=false;
			try{new TadpoleMulti.Config(all.toArray(new String[0]));}catch(IllegalArgumentException expected){failed=true;}
			check(failed, "Invalid neural policy was silently accepted.");
		}
	}

	/** Creates original-oriented inputs while preserving the separately stored trace bases. */
	private static Contig contig(final byte[] bases, final boolean reverse, final int id){
		final byte[] copy=bases.clone();
		if(reverse){AminoAcid.reverseComplementBasesInPlace(copy);}
		return new Contig(copy, id);
	}

	/** Sets only endpoint metadata needed by the overlapper's reciprocal fixtures. */
	private static ArrayList<Contig> contigs(final String... sequences){
		final ArrayList<Contig> answer=new ArrayList<Contig>();
		for(String sequence : sequences){
			final Contig c=contig(sequence.getBytes(StandardCharsets.US_ASCII), false, answer.size());
			c.leftBridgeEndpoint=c.rightBridgeEndpoint=true;
			c.leftCode=c.rightCode=Tadpole.DEAD_END;
			c.coverage=20;
			answer.add(c);
		}
		return answer;
	}

	/** Reproducible unique-looking flanks cover both short and capped context cases. */
	private static byte[] randomBases(final int length, final long seed){
		assert(length>0) : "A feature fixture must contain a retained sequence.";
		final byte[] bases=new byte[length];
		final Random random=new Random(seed);
		for(int i=0; i<length; i++){bases[i]=(byte)"ACGT".charAt(random.nextInt(4));}
		return bases;
	}

	/** Fixture failures must remain visible even if assertions are disabled accidentally. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** An explicit score stub separates decision-policy fixtures from learned weights. */
	private static final class FixedNet extends CellNet {
		/** No matrices or global network state are needed for a fixed-output fixture. */
		FixedNet(final float value_){value=value_;}
		@Override
		public int numInputs(){return 109;}
		@Override
		public int numOutputs(){return 1;}
		@Override
		public CellNet applyInput(final float[] values){return this;}
		@Override
		public float feedForward(){return value;}
		float value;
	}
}
