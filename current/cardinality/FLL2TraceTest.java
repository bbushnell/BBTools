package cardinality;

import rand.FastRandomXoshiro;

/**
 * Traces a single FutureLogLog2 word (6 buckets): prints register state,
 * CF-table lookup, and compares with the real estimator output at each step.
 *
 * @author Brian Bushnell
 */
public class FLL2TraceTest {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		FutureLogLog2.loadCFTable();

		final long seed=42;
		final FutureLogLog2 fll=new FutureLogLog2(6, 31, seed, 0);
		final FastRandomXoshiro rng=new FastRandomXoshiro(seed+1000000);

		System.err.printf("%-8s %-4s %-4s %-6s %-6s %-12s %-12s %-12s %-12s %-10s%n",
			"TrueC", "GE", "LE", "Hist", "Idx84", "TierAvg", "Mult", "CF_Est", "FLL2_Est", "Error%");

		for(long trueCard=1; trueCard<=500; trueCard++){
			fll.add(rng.nextLong());

			final int w=fll.getWord(0);
			final int ge=fll.getGlobalExp();
			final int le=(w>>>12)&0xF;
			final int hist=w&0xFFF;
			final int absNLZ=ge+le;
			final int tier=Math.min(absNLZ, FutureLogLog2.CF_TABLE_TIERS-1);
			final int eq=FutureLogLog2.idx84(hist);

			final double ta=(absNLZ<FutureLogLog2.tierAvg.length) ?
				FutureLogLog2.tierAvg[absNLZ] :
				FutureLogLog2.tierAvg[FutureLogLog2.tierAvg.length-1]*
				Math.pow(2.0, absNLZ-FutureLogLog2.tierAvg.length+1);
			final double mult=FutureLogLog2.cfTable[tier][eq];
			final double cfEst=(w==0) ? 0 : ta*mult;
			final long fllEst=fll.cardinality();
			final double err=(fllEst-trueCard)/(double)trueCard;

			final boolean print=trueCard<=40 || trueCard%50==0;
			if(print){
				System.err.printf("%-8d %-4d %-4d 0x%03X  %-6d %-12.2f %-12.6f %-12.2f %-12d %-10.4f%n",
					trueCard, ge, le, hist, eq, ta, mult, cfEst, fllEst, err);
			}
		}
	}
}
