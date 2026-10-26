package test.map;

import map.AbstractBitSet;
import map.AtomicBitSet;

/** Checks backing-array reuse at the first capacity requiring long arithmetic.
 * The largest case uses two 512 MiB arrays on the unfixed implementation.
 * @author Shinobu
 */
public final class BitSetCapacityTest{

	public static void main(String[] args){
		if(args.length!=2){throw new IllegalArgumentException("Expected type=raw|multi|atomic expect=baseline|fixed");}
		final boolean fixed=args[1].equals("expect=fixed");
		check(fixed || args[1].equals("expect=baseline"), "Unknown expectation");
		final AbstractBitSet set;
		if(args[0].equals("type=raw")){set=AbstractBitSet.make(1, 1);}
		else if(args[0].equals("type=multi")){set=AbstractBitSet.make(1, 2);}
		else if(args[0].equals("type=atomic")){set=new AtomicBitSet(1);}
		else{throw new IllegalArgumentException(args[0]);}
		final long capacity=1L<<31;
		set.setCapacity(capacity, 0);
		set.increment(Integer.MAX_VALUE);
		check(set.getCount(Integer.MAX_VALUE)==1, "Highest int-indexed bit must initially be set");
		set.setCapacity(capacity, 0);
		check(set.getCount(Integer.MAX_VALUE)==(fixed ? 1 : 0), "Same-capacity reuse must preserve the backing array");
		check(set.capacity()==capacity, "Logical capacity must remain unchanged");
		System.out.println("PASS\t"+args[0]+"\tfixed="+fixed+"\tretained="+set.getCount(Integer.MAX_VALUE));
	}

	private static void check(boolean ok, String message){
		if(!ok){throw new AssertionError(message);}
	}
}
