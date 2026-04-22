package cardinality;

import rand.FastRandomXoshiro;
import shared.Tools;

/**
 * Measures DLL4 branch rates at various cardinalities.
 * Inserts elements via Xoshiro RNG and reports minZeros and both branch rates
 * at each target cardinality.
 *
 * @author Brian Bushnell
 */
public class BranchRateTest {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		final int buckets=2048;
		final long[] cards={100, 1000, 10000, 100000, 1000000, 10000000};
		for(long target : cards){
			final DynamicLogLog4 dll=new DynamicLogLog4(buckets, 31, 12345L, 0);
			final FastRandomXoshiro rng=new FastRandomXoshiro(99);
			for(long i=0; i<target; i++){
				dll.add(rng.nextLong());
			}
			System.out.println("card="+target
				+" minZeros="+dll.getMinZeros()
				+" branch1Rate="+String.format("%.6f%%", dll.branch1Rate()*100)
				+" branch2Rate="+String.format("%.6f%%", dll.branch2Rate()*100));
		}
	}

}
