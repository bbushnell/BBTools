package assemble;

import java.util.ArrayList;
import java.util.Arrays;

import ml.CellNet;
import parse.LineParser1;
import structures.ByteBuilder;
import structures.IntList;
import ukmer.Kmer;

/** Exact decile, scaling, serialization and lazy inference parity checks.
 * @author Fischl */
public final class LocalEditReadDepthTest {

	private LocalEditReadDepthTest(){}

	/** Run small deterministic fixtures without a production model or cluster. */
	public static void main(final String[] args){
		testDeciles();
		testRolling();
		testGate();
		testDualGate();
		System.out.println("LOCAL_EDIT_READ_DEPTH_TEST_OK deciles=11 features=39 dual_boundary=30/31");
	}

	/** Verify nearest-rank endpoints, ties, zero clamping and model-input parity. */
	private static void testDeciles(){
		final IntList depths=new IntList();
		for(int i=99; i>=0; i--){depths.add(i);}
		final int[] percentiles=new int[11];
		LocalEditNeuralFeatures.selectDepthPercentiles(depths, percentiles);
		for(int i=0; i<=10; i++){
			check(percentiles[i]==(i==0 ? 0 : i*10-1), "Incorrect nearest-rank decile "+i);
		}
		depths.clear();
		depths.add(256);
		LocalEditNeuralFeatures.selectDepthPercentiles(depths, percentiles);
		final double[] features=new double[39];
		check(LocalEditNeuralFeatures.appendReadDepths(features, percentiles), "Singleton depth rejected");
		for(int i=28; i<39; i++){check(Math.abs(features[i]-1)<1e-12, "Depth256 must scale to1");}
		for(int i=0; i<=10; i++){percentiles[i]=i==0 ? 0 : 1<<(i-1);}
		check(LocalEditNeuralFeatures.appendReadDepths(features, percentiles), "Ordered depths rejected");
		for(int i=0; i<=10; i++){
			check(Math.abs(features[28+i]-(i==0 ? 0 : (i-1)*0.125))<1e-12, "Incorrect log2 scaling");
		}
		final float[] input=new float[39];
		check(LocalEditNeuralFeatures.toModelInput(features, input), "Extended vector rejected");
		final ByteBuilder bb=new ByteBuilder();
		LocalEditNeuralFeatures.append(bb, features, 1);
		final LineParser1 parser=new LineParser1('\t');
		parser.set(bb.toBytes());
		check(parser.terms()==40, "Extended vector must contain39 features plus label");
		for(int i=0; i<39; i++){
			check(Float.floatToIntBits(input[i])==Float.floatToIntBits(parser.parseFloat(i)), "Serialization drift at "+i);
		}
		percentiles[5]=-1;
		check(!LocalEditNeuralFeatures.appendReadDepths(features, percentiles), "Negative depth accepted");
		check(!LocalEditNeuralFeatures.toModelInput(features, new float[28]), "Mixed feature dimensions accepted");
	}

	/** Confirm that N breaks the rolling key and absent valid kmers retain zero. */
	private static void testRolling(){
		final byte[] bases="AAAAANAAAAA".getBytes(java.nio.charset.StandardCharsets.US_ASCII);
		final IntList scratch=new IntList();
		final int[] percentiles=new int[11];
		final HomopolymerIndelProposal.CountLookup absent=new HomopolymerIndelProposal.CountLookup(){
			@Override public int count(final Kmer key){return -1;}
		};
		check(LocalEditNeuralFeatures.readDepthPercentiles(bases, new Kmer(5), absent, scratch, percentiles)==2,
				"Only two valid5-mers surround the N");
		for(final int depth:percentiles){check(depth==0, "Absent kmers must have depth0");}
		check(LocalEditNeuralFeatures.readDepthPercentiles(new byte[]{'N'}, new Kmer(5), absent, scratch, percentiles)==0,
				"Empty profile must not fabricate a lookup");
		check(scratch.size==0, "Scratch leaked prior-read depths");
	}

	/** Match a direct39-input score and prove whole-read counting is cached/reset. */
	private static void testGate(){
		ml.Function.normalizeTypeRates();
		final CellNet model=new CellNet(new int[]{39, 8, 1}, 7, 1, 1, 1, new ArrayList<String>());
		model.randomize();
		final HomopolymerIndelProposal.CountLookup lookup=new HomopolymerIndelProposal.CountLookup(){
			@Override public int count(final Kmer key){return 5;}
		};
		final LocalEditNeuralGate gate=new LocalEditNeuralGate(62, lookup, model.copy(false), 0.5f);
		final byte[] bases=new byte[140], alphabet={'A', 'C', 'G', 'T'};
		for(int i=0; i<bases.length; i++){bases[i]=alphabet[i&3];}
		gate.beginRead(bases);
		final float observed=gate.score(bases, 70, LocalSingleBaseEdit.Operation.SUBSTITUTION, (byte)'A');
		check(Float.isFinite(observed), "Complete depth-aware evidence must score");
		check(gate.windowQueries==79+3, "Whole-read profile must query each of79 valid windows once");
		final double[] features=new double[39];
		final float[] input=new float[39];
		final int[] depths=new int[8], percentiles=new int[11];
		Arrays.fill(depths, 5);
		Arrays.fill(percentiles, 5);
		check(LocalEditNeuralFeatures.fill(features, 5, depths, 0, 2, 5, 5, 1, false, 140, 70), "Legacy prefix rejected");
		check(LocalEditNeuralFeatures.appendReadDepths(features, percentiles), "Depth suffix rejected");
		check(LocalEditNeuralFeatures.toModelInput(features, input), "Direct input rejected");
		model.simdFF=true;
		final float expected=model.applyInput(input).feedForwardDense();
		check(Float.floatToIntBits(observed)==Float.floatToIntBits(expected), "Training/inference39-input score differs");
		gate.score(bases, 70, LocalSingleBaseEdit.Operation.SUBSTITUTION, (byte)'A');
		check(gate.windowQueries==79+6, "A second candidate rescanned the read");
		gate.beginRead(bases);
		gate.score(bases, 70, LocalSingleBaseEdit.Operation.SUBSTITUTION, (byte)'A');
		check(gate.windowQueries==158+9, "A new transaction reused stale read depths");
	}

	/** Verify original-read P50 routing, low/high cutoffs and per-transaction reset. */
	private static void testDualGate(){
		final CellNet low=new CellNet(new int[]{39, 8, 1}, 17, 1, 1, 1, new ArrayList<String>());
		final CellNet high=new CellNet(new int[]{39, 8, 1}, 19, 1, 1, 1, new ArrayList<String>());
		low.randomize(); high.randomize();
		final int[] depth=new int[]{30};
		final HomopolymerIndelProposal.CountLookup lookup=new HomopolymerIndelProposal.CountLookup(){
			@Override public int count(final Kmer key){return depth[0];}
		};
		final LocalEditNeuralGate gate=new LocalEditNeuralGate(62, lookup, low.copy(false), 0f, high.copy(false), 1f);
		final byte[] bases=new byte[140], alphabet={'A', 'C', 'G', 'T'};
		for(int i=0; i<bases.length; i++){bases[i]=alphabet[i&3];}
		gate.beginRead(bases);
		check(gate.accept(bases, 70, LocalSingleBaseEdit.Operation.SUBSTITUTION, (byte)'A'),
				"P50=30 must route to permissive low-depth cutoff");
		check(!gate.lastHighDepth && Float.floatToIntBits(gate.lastCutoff)==Float.floatToIntBits(0f),
				"P50=30 selected the high-depth model");
		check(gate.windowQueries==79+3, "Dual low-depth route did not cache whole-read percentiles");
		depth[0]=31;
		check(gate.accept(bases, 70, LocalSingleBaseEdit.Operation.SUBSTITUTION, (byte)'A'),
				"Route changed before beginRead reset");
		gate.beginRead(bases);
		check(!gate.accept(bases, 70, LocalSingleBaseEdit.Operation.SUBSTITUTION, (byte)'A'),
				"P50=31 must route to rejecting high-depth cutoff");
		check(gate.lastHighDepth && Float.floatToIntBits(gate.lastCutoff)==Float.floatToIntBits(1f),
				"P50=31 selected the low-depth model");
	}

	/** Keep fixture failures active under both assertion modes. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
