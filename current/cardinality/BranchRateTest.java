package cardinality;

import rand.FastRandomXoshiro;
import shared.Tools;

/** Quick test to measure DLL4 branch rates at various cardinalities. */
public class BranchRateTest {
	public static void main(String[] args){
		int buckets=2048;
		long[] cards={100, 1000, 10000, 100000, 1000000, 10000000};
		for(long target : cards){
			DynamicLogLog4 dll=new DynamicLogLog4(buckets, 31, 12345L, 0);
			FastRandomXoshiro rng=new FastRandomXoshiro(99);
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
