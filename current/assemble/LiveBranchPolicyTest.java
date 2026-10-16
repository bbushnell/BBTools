package assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/** Deterministic invariance and guard fixtures for the new live-node schema. @author Fischl */
public final class LiveBranchPolicyTest {

	/*--------------------------------------------------------------*/
	/*----------------            Entry             ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs local unit tests only; does not load a model or change an assembly. */
	public static void main(final String[] args){
		if(args.length!=0){throw new IllegalArgumentException("LiveBranchPolicyTest accepts no program arguments.");}
		boolean assertions=false;
		assert(assertions=true);
		if(!assertions){throw new IllegalStateException("Run the test with assertions enabled through its launcher.");}
		guards(); features(); symmetry(); invalid();
		System.err.println("LIVE_BRANCH_POLICY_PASS checks="+checks);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Hard stops cannot be softened by swapping direction or by a score. */
	private static void guards(){
		final int[] clear={10, 0, 0, 0}, grey={10, 2, 0, 0};
		check(level(clear, grey, false, false)==LiveBranchPolicy.NATIVE_ACCEPT, "Native accepts bypass the new model.");
		check(level(clear, grey, false, true)==LiveBranchPolicy.GREY_ZONE, "10/2 is an eligible experimental grey zone.");
		check(level(grey, grey, true, true)==LiveBranchPolicy.GREY_ZONE, "Both sides use one node-level decision.");
		check(level(clear, new int[]{10, 5, 0, 0}, false, true)==LiveBranchPolicy.HARD_STOP, "Half-depth equality remains a hard stop.");
		check(level(clear, new int[]{10, 4, 0, 0}, false, true)==LiveBranchPolicy.GREY_ZONE, "Below-screen competitor is eligible.");
		check(level(clear, new int[]{10, 10, 0, 0}, false, true)==LiveBranchPolicy.HARD_STOP, "Equal best routes are ambiguous.");
		check(level(new int[]{0, 0, 0, 0}, grey, false, true)==LiveBranchPolicy.HARD_STOP, "No predecessor cannot be rescued.");
		check(level(new int[]{1, 0, 0, 0}, grey, false, true)==LiveBranchPolicy.HARD_STOP, "Preferred predecessor must meet mce.");
		check(LiveBranchPolicy.level(1, clear, grey, false, true, 2, 0.5)==LiveBranchPolicy.HARD_STOP, "Node depth must meet mce.");
		check(LiveBranchPolicy.level(20, clear, grey, false, true, 2, 0.2)==LiveBranchPolicy.HARD_STOP, "Configured equality is a stop.");
		check(LiveBranchPolicy.level(20, clear, grey, false, true, 2, 1)==LiveBranchPolicy.GREY_ZONE, "Screen is configurable without permitting ties.");
		check(LiveBranchPolicy.preferredBase(new int[]{0, 0, 0, 10})==3, "Raw preferred base survives independent sorting.");
		check(LiveBranchPolicy.preferredBase(new int[]{0, 0, 0, 0})==-1, "All-zero side has no preferred edge.");
		check(LiveBranchPolicy.preferredBase(new int[]{10, 10, 0, 0})==-1, "Tied raw edges have no preferred base.");
		check(level(new int[]{Integer.MAX_VALUE, 2, 0, 0}, grey, true, true)==LiveBranchPolicy.GREY_ZONE, "High counts do not overflow the ratio test.");
	}

	/** Fixed slots and transforms are checked independently of permutation testing. */
	private static void features(){
		final int[] left={0, 1, 3, 7}, right={15, 31, 63, 127};
		final float[] out=new float[LiveBranchPolicy.INPUTS];
		LiveBranchPolicy.fill(out, 64, 255, left, right, 511, 1023);
		final float[] expected={1, 0.375f, 0.25f, 0.125f, 0, 0.875f, 0.75f, 0.625f, 0.5f, 0.75f, 1.125f, 1.25f};
		for(int i=0; i<out.length; i++){check(Math.abs(out[i]-expected[i])<0.000001, "Feature slot/normalization mismatch at "+i);}
		check(out[3]>out[4], "Singleton and absence must differ.");
		check(Arrays.equals(left, new int[]{0, 1, 3, 7}) && Arrays.equals(right, new int[]{15, 31, 63, 127}), "Features must not mutate caller counts.");
		LiveBranchPolicy.fill(out, 127, Integer.MAX_VALUE, right, left, Integer.MAX_VALUE, Integer.MAX_VALUE);
		for(float value : out){check(Float.isFinite(value), "Extreme finite depths must remain finite.");}
	}

	/** Tests all 24x24 nucleotide permutations plus side swap and reverse complement. */
	private static void symmetry(){
		final ArrayList<int[]> permutations=new ArrayList<int[]>();
		permutations(new int[4], 0, 0, permutations);
		check(permutations.size()==24, "Permutation fixture must cover all nucleotide assignments.");
		final Random random=new Random(381401);
		final float[] expected=new float[12], actual=new float[12];
		for(int sample=0; sample<64; sample++){
			final int k=24+random.nextInt(104), count=2+random.nextInt(500);
			final int[] left=new int[4], right=new int[4];
			for(int i=0; i<4; i++){left[i]=random.nextInt(200); right[i]=random.nextInt(200);}
			if(sample%4==0){left[0]=Integer.MAX_VALUE; right[0]=Integer.MAX_VALUE-1;}
			final boolean lb=(sample&1)!=0, rb=(sample&2)!=0;
			final int decision=LiveBranchPolicy.level(count, left, right, lb, rb, 2, 0.5);
			LiveBranchPolicy.fill(expected, k, count, left, right, 20.5, 80.25);
			for(int[] lp : permutations){for(int[] rp : permutations){
				final int[] a=permute(left, lp), b=permute(right, rp);
				LiveBranchPolicy.fill(actual, k, count, a, b, 20.5, 80.25);
				check(Arrays.equals(expected, actual), "Base assignment changed the feature vector.");
				check(LiveBranchPolicy.level(count, a, b, lb, rb, 2, 0.5)==decision, "Base assignment changed the screen.");
				LiveBranchPolicy.fill(actual, k, count, b, a, 20.5, 80.25);
				check(Arrays.equals(expected, actual), "Traversal reversal changed the feature vector.");
				check(LiveBranchPolicy.level(count, b, a, rb, lb, 2, 0.5)==decision, "Traversal reversal changed the screen.");
			}}
			final int[] rcLeft={right[3], right[2], right[1], right[0]};
			final int[] rcRight={left[3], left[2], left[1], left[0]};
			LiveBranchPolicy.fill(actual, k, count, rcLeft, rcRight, 20.5, 80.25);
			check(Arrays.equals(expected, actual), "Reverse complement changed the feature vector.");
			final int preferred=LiveBranchPolicy.preferredBase(right);
			check(LiveBranchPolicy.preferredBase(rcLeft)==(preferred<0 ? -1 : 3-preferred), "Preferred raw edge did not complement correctly.");
		}
	}

	/** Supplies bad evidence deliberately; configuration/data errors must be loud. */
	private static void invalid(){
		final int[] counts={10, 2, 0, 0};
		expectInvalid(()->LiveBranchPolicy.preferredBase(new int[]{1, 2}));
		expectInvalid(()->LiveBranchPolicy.preferredBase(new int[]{1, -1, 0, 0}));
		expectInvalid(()->LiveBranchPolicy.level(10, counts, counts, true, true, 0, 0.5));
		expectInvalid(()->LiveBranchPolicy.level(10, counts, counts, true, true, 2, Double.NaN));
		expectInvalid(()->LiveBranchPolicy.level(10, counts, counts, true, true, 2, 0));
		expectInvalid(()->LiveBranchPolicy.level(10, counts, counts, true, true, 2, 1.1));
		expectInvalid(()->LiveBranchPolicy.fill(new float[27], 32, 10, counts, counts, 10, 20));
		expectInvalid(()->LiveBranchPolicy.fill(new float[12], 1, 10, counts, counts, 10, 20));
		expectInvalid(()->LiveBranchPolicy.fill(new float[12], 32, 10, counts, counts, -1, 20));
		expectInvalid(()->LiveBranchPolicy.fill(new float[12], 32, 10, counts, counts, 10, Double.POSITIVE_INFINITY));
	}

	/** Small fixture wrapper fixes node depth/mce without hiding native flags. */
	private static int level(final int[] left, final int[] right, final boolean lb, final boolean rb){
		assert(left!=null && right!=null) : "Guard fixtures require explicit independent side evidence.";
		return LiveBranchPolicy.level(20, left, right, lb, rb, 2, 0.5);
	}

	/** Enumerates permutations once; allocations belong to test fixtures, not inference. */
	private static void permutations(final int[] value, final int position, final int used,
			final ArrayList<int[]> out){
		assert(position>=0 && position<=4) : "Exactly four nucleotide indices form one permutation.";
		if(position==4){out.add(value.clone()); return;}
		for(int i=0; i<4; i++){
			if((used&(1<<i))==0){value[position]=i; permutations(value, position+1, used|(1<<i), out);}
		}
	}

	/** Allocates only a temporary test case; raw counts remain unchanged. */
	private static int[] permute(final int[] input, final int[] order){
		assert(input.length==4 && order.length==4) : "Fixture permutations must cover the whole nucleotide alphabet.";
		return new int[]{input[order[0]], input[order[1]], input[order[2]], input[order[3]]};
	}

	/** Requires the intended validation exception, not any arbitrary failure. */
	private static void expectInvalid(final Runnable test){
		boolean failed=false;
		try{test.run();}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Malformed evidence/configuration was accepted.");
	}

	/** Counts checks even with Java assertions disabled, though main also requires -ea. */
	private static void check(final boolean condition, final String reason){
		checks++;
		if(!condition){throw new AssertionError(reason);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private static long checks;
}
