package cardinality;

/**
 * Compare UDLL6 vs ULLc converted register distributions at card=16384, 512 buckets.
 */
public class UDLL6Test {
	public static void main(String[] args){
		final int buckets=512;
		final int k=31;
		final long seed=12345L;

		CardinalityTracker.clampToAdded=false;
		//TRACE removed in UDLL6 rewrite
		UltraDynamicLogLog6 udll=new UltraDynamicLogLog6(buckets, k, seed, 0);
		UltraDynamicLogLog6 ullc=new UltraDynamicLogLog6(buckets, k, seed, 0);

		for(long i=1; i<=16384; i++){
			udll.add(i);
			ullc.add(i);
		}
		udll.lastCardinality=-1;
		ullc.lastCardinality=-1;

		System.out.println("UDLL6 mz="+udll.getMinZeros()+" ULLc mz="+ullc.getMinZeros());
		System.out.println("UDLL6 fgra="+udll.fgraEstimate()+" ULLc fgra="+ullc.fgraEstimate());

		// Compare every register (converted to absolute Ertl format)
		int p=Integer.numberOfTrailingZeros(buckets); // log2(buckets)
		int regOffset=4*(udll.getMinZeros()+p-1-2);
		int diffs=0;
		int maxDiff=0;
		for(int b=0; b<buckets; b++){
			int ur=udll.getRegister(b)&0x3F;
			int urAbs=(ur==0) ? 0 : Math.min(ur+regOffset, 255);
			int cr=ullc.getRegister(b)&0xFF;
			if(urAbs!=cr){
				diffs++;
				int d=Math.abs(urAbs-cr);
				if(d>maxDiff){maxDiff=d;}
				if(diffs<=10){
					System.out.printf("  bucket=%d  udll6_rel=%d  udll6_abs=%d  ullc_abs=%d  diff=%d%n",
						b, ur, urAbs, cr, urAbs-cr);
				}
			}
		}
		System.out.println("Total differing buckets: "+diffs+"/"+buckets+" maxDiff="+maxDiff);
	}
}
