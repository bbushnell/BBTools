package test.benchmark;

import java.util.Arrays;
import benchmark.FamilyMemberLengthCensus;
import prot.MemberLengthRange;

/** Boundary, histogram and family-aggregation oracles for the member census.
 * Uses literal synthetic lengths and independent expected bins, medians and sums;
 * it reads no corpus and writes no evidence files. A failed check aborts immediately.
 * @author Fischl
 */
public final class FamilyMemberLengthCensusTest{

	/** Runs legacy and explicit-range checks, then prints one success marker.
	 * @param args Empty argument array; this fixture takes no external inputs */
	public static void main(final String[] args){
		if(args.length!=0){throw new IllegalArgumentException("No test arguments expected");}
		final int[] lengths={100, 124, 125, 149, 150, 199, 200};
		final int[] expected={0, 0, 1, 1, 2, 2, 3};
		for(int i=0; i<lengths.length; i++){
			check(FamilyMemberLengthCensus.lengthBin(lengths[i])==expected[i], "Length boundary "+lengths[i]);
		}
		for(int length:new int[]{99, 201}){
			boolean failed=false;
			try{FamilyMemberLengthCensus.lengthBin(length);}catch(IllegalArgumentException e){failed=true;}
			check(failed, "Out-of-contract length must fail: "+length);
		}
		final long[] bins=new long[4];
		check(FamilyMemberLengthCensus.addFamily(new byte[][]{new byte[100], new byte[124], new byte[125],
			new byte[150], new byte[199], new byte[200]}, bins)==5,
			"Five short members qualify regardless of the200aa member");
		check(FamilyMemberLengthCensus.addFamily(new byte[][]{new byte[149], new byte[200]}, bins)==1,
			"Count each family independently");
		check(FamilyMemberLengthCensus.addFamily(new byte[0][], bins)==0, "An empty family contributes no members");
		check(Arrays.equals(bins, new long[]{2, 2, 2, 2}), "All eight members must enter exactly one bin");
		final int[] histogram=new int[201];
		FamilyMemberLengthCensus.addFamily(new byte[][]{new byte[100], new byte[200]}, new long[4], histogram);
		check(FamilyMemberLengthCensus.quantile(histogram, 2, 10)==100
			&& FamilyMemberLengthCensus.quantile(histogram, 2, 50)==150
			&& FamilyMemberLengthCensus.quantile(histogram, 2, 90)==200,
			"Even-family median is the mean of the two central lengths, with nearest-rank deciles");
		histogram[150]=1;
		check(FamilyMemberLengthCensus.quantile(histogram, 3, 50)==150, "Odd-family median is its central member length");
		Arrays.fill(histogram, 0);
		histogram[124]=10;
		histogram[199]=80;
		histogram[200]=10;
		check(FamilyMemberLengthCensus.quantile(histogram, 100, 10)==124
			&& FamilyMemberLengthCensus.quantile(histogram, 100, 50)==199
			&& FamilyMemberLengthCensus.quantile(histogram, 100, 90)==199,
			"Exact decile endpoints must use complete empirical family counts");
		boolean wrongCount=false;
		try{FamilyMemberLengthCensus.quantile(histogram, 99, 50);}catch(IllegalArgumentException expectedFailure){wrongCount=true;}
		check(wrongCount, "A different denominator must fail rather than produce plausible but wrong percentiles");
		shortRanges();
		memberMeans();
		System.out.println("FAMILY_MEMBER_LENGTH_CENSUS_TEST_PASS");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Extended ranges        ----------------*/
	/*--------------------------------------------------------------*/

	/** Literal populations distinguish arithmetic means from medians and preserve exact boundary numerators. */
	private static void memberMeans(){
		final MemberLengthRange range=new MemberLengthRange(40, 200);
		final int[] histogram=new int[201];
		check(FamilyMemberLengthCensus.lengthSum(histogram, 0, range)==0, "Empty family has zero residues, not a measured mean");
		histogram[40]=1;
		histogram[64]=1;
		histogram[99]=1;
		check(FamilyMemberLengthCensus.lengthSum(histogram, 3, range)==203
			&& FamilyMemberLengthCensus.quantile(histogram, 3, 50, range)==64,
			"Lengths40,64,99 have mean203/3 above65 but median64 below65");
		Arrays.fill(histogram, 0);
		histogram[64]=1;
		histogram[65]=2;
		check(FamilyMemberLengthCensus.lengthSum(histogram, 3, range)==194, "Mean194/3 stays below65 without rounding");
		histogram[64]=0;
		histogram[65]=3;
		check(FamilyMemberLengthCensus.lengthSum(histogram, 3, range)==195, "Mean195/3 is exactly65");
		Arrays.fill(histogram, 0);
		histogram[99]=1;
		histogram[100]=1;
		check(FamilyMemberLengthCensus.lengthSum(histogram, 2, range)==199, "Mean199/2 stays below100");
		histogram[99]=0;
		histogram[100]=2;
		check(FamilyMemberLengthCensus.lengthSum(histogram, 2, range)==200, "Mean200/2 is exactly100");
		reject(()->FamilyMemberLengthCensus.lengthSum(histogram, 3, range),
			"A wrong denominator must not change the family band silently");
		reject(()->FamilyMemberLengthCensus.lengthSum(histogram, 2, new MemberLengthRange(40, 99)),
			"Out-of-range lengths cannot enter family means");
		histogram[100]=-1;
		reject(()->FamilyMemberLengthCensus.lengthSum(histogram, 0, range),
			"Negative histogram counts cannot cancel valid members");
		Arrays.fill(histogram, 0);
		histogram[200]=Integer.MAX_VALUE;
		check(FamilyMemberLengthCensus.lengthSum(histogram, Integer.MAX_VALUE, range)==429496729400L,
			"Residue sums must use long arithmetic even though each histogram count is an int");
	}

	/** Literal boundary and order-statistic oracles cover the newly admitted range and legacy rejection. */
	private static void shortRanges(){
		final MemberLengthRange range=new MemberLengthRange(40, 99);
		final int[] lengths={40, 59, 60, 79, 80, 99}, expected={0, 0, 1, 1, 2, 2};
		for(int i=0; i<lengths.length; i++){
			check(FamilyMemberLengthCensus.lengthBin(lengths[i], range)==expected[i], "Short bin boundary "+lengths[i]);
		}
		for(int length:new int[]{39, 100}){
			reject(()->FamilyMemberLengthCensus.lengthBin(length, range), "Declared40..99 must reject "+length);
		}
		final long[] bins=new long[7];
		final int[] histogram=new int[201];
		final byte[][] members=new byte[lengths.length][];
		for(int i=0; i<lengths.length; i++){members[i]=new byte[lengths[i]];}
		check(FamilyMemberLengthCensus.addFamily(members, bins, histogram, range)==6, "Every short-range member is below200");
		check(Arrays.equals(bins, new long[]{2, 2, 2, 0, 0, 0, 0}), "Six boundary members partition into three short bins");
		check(FamilyMemberLengthCensus.quantile(histogram, 6, 10, range)==40
			&& FamilyMemberLengthCensus.quantile(histogram, 6, 50, range)==69.5
			&& FamilyMemberLengthCensus.quantile(histogram, 6, 90, range)==99,
			"Sorted40,59,60,79,80,99 has median69.5 and nearest-rank endpoints40/99");
		reject(()->FamilyMemberLengthCensus.quantile(histogram, 6, 50), "Legacy quantiles still reject short observations");
		reject(()->FamilyMemberLengthCensus.addFamily(members, new long[4], null, range),
			"Wrong counter layout cannot discard short bins");
		histogram[60]++;
		check(FamilyMemberLengthCensus.quantile(histogram, 7, 50, range)==60, "Odd short median uses the fourth of seven members");
		histogram[100]=1;
		reject(()->FamilyMemberLengthCensus.quantile(histogram, 8, 50, range),
			"Histogram cannot include a length above its declared maximum");
		final long[] mixed=new long[7];
		check(FamilyMemberLengthCensus.addFamily(new byte[][]{new byte[40], new byte[100], new byte[199], new byte[200]},
			mixed, null, new MemberLengthRange(40, 200))==3, "Mixed-range below200 count is independent of bin numbering");
		check(Arrays.equals(mixed, new long[]{1, 0, 0, 1, 0, 1, 1}), "Mixed-range histogram retains legacy boundaries");
	}

	/** Failure is required rather than an apparently valid statistic from out-of-contract input. */
	private static void reject(final Runnable action, final String why){
		try{action.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError(why);
	}

	/** Requires the literal expectation independently of whether JVM assertions are enabled. */
	private static void check(final boolean valid, final String why){
		if(!valid){throw new AssertionError(why);}
	}
}
