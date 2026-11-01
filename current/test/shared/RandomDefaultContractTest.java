package test.shared;

/** Finite checks that exercise the shared.Random defaults rather than rand overrides.
 * A scripted long source makes the wide-interval collapse reproducible and caps
 * consumption so a rejection-loop regression fails instead of hanging.
 * @author Sucrose
 */
public final class RandomDefaultContractTest{

	/** Runs zero-bound, wide-interval and valid bounded controls; arguments are unused. */
	public static void main(String[] args){
		int failures=0;
		try{zeroBound();System.out.println("zero-bound PASS");}
		catch(RuntimeException | AssertionError e){failures++;System.out.println("zero-bound FAIL: "+e);}
		try{wideInterval();System.out.println("wide-interval PASS");}
		catch(RuntimeException | AssertionError e){failures++;System.out.println("wide-interval FAIL: "+e);}
		try{boundedControls();System.out.println("bounded-controls PASS");}
		catch(RuntimeException | AssertionError e){failures++;System.out.println("bounded-controls FAIL: "+e);}
		System.out.println("RandomDefaultContractTest cases=3 failures="+failures);
		if(failures>0){throw new AssertionError("Default random contract failures="+failures);}
	}

	/** An empty interval must be rejected rather than returning an arbitrary integer. */
	private static void zeroBound(){
		final ScriptedRandom random=new ScriptedRandom();
		boolean rejected=false;
		try{random.nextInt(0);}
		catch(IllegalArgumentException | AssertionError expected){rejected=true;}
		require(rejected, "nextInt(0) returned a value for an empty interval");
	}

	/** A valid nearly full signed-int interval must not reject or collapse every draw.
	 * This checks nondegeneracy, not an exact seeded sequence or a statistical claim.
	 * The fixed scripted source varies both signs of the low word; the broken width
	 * calculation passes -1 to the bounded default, yielding the same origin each time.
	 */
	private static void wideInterval(){
		final ScriptedRandom random=new ScriptedRandom();
		final int first=random.nextInt(Integer.MIN_VALUE, Integer.MAX_VALUE);
		boolean distinct=false;
		require(first<Integer.MAX_VALUE, "exclusive upper bound returned");
		for(int i=1; i<128; i++){
			final int value=random.nextInt(Integer.MIN_VALUE, Integer.MAX_VALUE);
			require(value<Integer.MAX_VALUE, "exclusive upper bound returned");
			distinct|=value!=first;
		}
		require(distinct, "128 varying source words collapsed to one value: "+first);
	}

	/** Power-of-two, rejection-sampled and shifted narrow intervals remain valid. */
	private static void boundedControls(){
		final ScriptedRandom random=new ScriptedRandom();
		for(int i=0; i<32; i++){
			final int power=random.nextInt(16), odd=random.nextInt(7), shifted=random.nextInt(-3, 4);
			require(power>=0 && power<16, "power-of-two bounds");
			require(odd>=0 && odd<7, "odd bounds");
			require(shifted>=-3 && shifted<4, "shifted bounds");
		}
	}

	/** Fails independently of JVM assertion settings. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Implements only required primitives so all int calls exercise interface defaults. */
	private static final class ScriptedRandom implements shared.Random{
		/** Finite repeating source; high bits keep the old -1-bound path out of its retry loop. */
		@Override
		public long nextLong(){
			if(calls>=4096){throw new IllegalStateException("Random source consumption exceeded finite test budget");}
			return WORDS[calls++&3];
		}
		/** Resets test consumption; seed values are not interpreted by this scripted source. */
		@Override
		public void setSeed(long seed){calls=0;}
		/** Supplies a valid bounded-long primitive if an int implementation delegates to it. */
		@Override
		public long nextLong(long bound){
			if(bound<=0){throw new IllegalArgumentException("Positive long bound required");}
			long bits, value;
			do{bits=nextLong()>>>1;value=bits%bound;}while(bits-value+(bound-1)<0);
			return value;
		}
		/** Unused primitive. */
		@Override
		public double nextGaussian(){throw new UnsupportedOperationException("Not used by int defaults");}
		/** Number of long words consumed. */
		private int calls;
		/** Distinct positive/negative low words; no word equals the excluded upper int bound. */
		private static final long[] WORDS={0x4000000000000001L, 0x40000000FFFFFFFEL, 0x4000000040000000L, 0x4000000080000000L};
	}
}
