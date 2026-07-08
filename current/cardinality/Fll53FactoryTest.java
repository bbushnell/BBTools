package cardinality;

/** One-off: verify makeTracker("fll53") end-to-end. */
public class Fll53FactoryTest {
	public static void main(String[] args){
		parse.Parser.loglogType="fll53";
		final CardinalityTracker t=CardinalityTracker.makeTracker("fll53");
		System.out.println("factory returned: "+t.getClass().getSimpleName());
		final java.util.SplittableRandom r=new java.util.SplittableRandom(777);
		final long N=1_000_000;
		for(long i=0; i<N; i++){t.add(r.nextLong());}
		final long est=t.cardinality();
		System.out.printf("unique 1M: est=%d err=%+.2f%%%n", est, 100.0*(est-N)/N);
		// dup-heavy: 4x reuse
		final CardinalityTracker t2=CardinalityTracker.makeTracker("fll53");
		final long[] pool=new long[250_000];
		for(int i=0; i<pool.length; i++){pool[i]=r.nextLong();}
		for(long i=0; i<N; i++){t2.add(pool[r.nextInt(pool.length)]);}
		final long est2=t2.cardinality();
		System.out.printf("dup4 1M adds ~250k distinct: est=%d err=%+.2f%%%n",
			est2, 100.0*(est2-250_000)/250_000);
		try{t.add(t2); System.out.println("MERGE DID NOT THROW — BUG");}
		catch(UnsupportedOperationException e){System.out.println("merge throws as designed");}
	}
}
