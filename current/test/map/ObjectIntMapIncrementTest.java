package test.map;

import map.ObjectIntMap;
import shared.Tools;

/** Finite regressions for absent counts and increment returns across rehashing.
 * @author Shinobu
 */
public final class ObjectIntMapIncrementTest{

	/** Runs explicit oracles under either assertion setting. */
	public static void main(String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final boolean fixed="expect=fixed".equals(args[0]);
		check(fixed || "expect=baseline".equals(args[0]), "Unknown expectation");
		final ObjectIntMap<String> fresh=new ObjectIntMap<String>(64, String.class);
		check(fresh.get("absent")==-1, "Absent get sentinel must remain -1");
		check(fresh.put("sentinel", 7)==-1, "Fresh absent put sentinel must remain -1");
		check(fresh.increment("count")==(fixed ? 1 : 0), "Fresh count origin");
		check(fresh.increment("amount", 5)==(fixed ? 5 : 4), "Absent arbitrary increment");
		check(fresh.increment("negative", -3)==(fixed ? -3 : -4), "Absent negative increment");
		fresh.put("existing", -1);
		check(fresh.increment("existing")==0, "Present -1 must be added, not treated as absence");
		fresh.put("maximum", Integer.MAX_VALUE-1);
		check(fresh.increment("maximum", 5)==Integer.MAX_VALUE, "Positive overflow clamp");
		fresh.clear();
		check(fresh.increment("count")==1, "Cleared slot count origin");
		check(fresh.put("cleared-slot", 9)==0, "Retain cleared-slot put behavior");
		check(fresh.remove("count"), "Remove existing count");
		check(fresh.increment("count")==1, "Removed slot count origin");
		checkResizeReturn(fixed);
		checkMerge(fixed);
		if(fixed){checkCollisions();}
		System.out.println("PASS\tfixed="+fixed+"\tget/put sentinels preserved");
	}

	/** Forces the first insertion to move from cell zero to cell one during growth. */
	private static void checkResizeReturn(final boolean fixed){
		int id=0;
		while((Tools.hash32plus(id)&1)==0){id++;}
		final Key key=new Key(id, id);
		final ObjectIntMap<Key> map=new ObjectIntMap<Key>(1, Key.class);
		final int returned=map.increment(key, 7);
		final int stored=map.get(key);
		check(stored==(fixed ? 7 : 6), "Growth must retain computed count");
		check(returned==(fixed ? 7 : -1), "Return must not read an obsolete cell after growth");
		check(map.size()==1, "Growth must preserve size");
		System.out.println("resize\treturned="+returned+"\tstored="+stored);
	}

	/** Merging count maps must not subtract an absent-slot sentinel. */
	private static void checkMerge(final boolean fixed){
		final ObjectIntMap<String> source=new ObjectIntMap<String>(64, String.class);
		final ObjectIntMap<String> dest=new ObjectIntMap<String>(64, String.class);
		source.increment("one");
		source.increment("many", 4);
		dest.incrementAll(source);
		check(dest.get("one")==(fixed ? 1 : -1), "Merged singleton count");
		check(dest.get("many")==(fixed ? 4 : 2), "Merged larger count");
		check(dest.size()==2, "Merged key count");
	}

	/** Exercises collision chains, repeated growth, deletion and signed existing values. */
	private static void checkCollisions(){
		final ObjectIntMap<Key> map=new ObjectIntMap<Key>(2, Key.class);
		final Key[] keys=new Key[80];
		for(int i=0; i<keys.length; i++){
			keys[i]=new Key(i, 11);
			check(map.increment(keys[i], i)==i, "Insertion return at "+i);
			check(map.get(keys[i])==i, "Inserted count at "+i);
		}
		for(int i=0; i<keys.length; i+=2){check(map.remove(keys[i]), "Removal at "+i);}
		for(int i=0; i<keys.length; i++){
			final int expected=(i%2==0 ? -2 : i-2);
			check(map.increment(keys[i], -2)==expected, "Collision increment at "+i);
			check(map.get(keys[i])==expected, "Collision stored value at "+i);
		}
		check(map.size()==keys.length, "Collision repair must preserve all keys");
	}

	/** Uses unconditional test checks so -da validates the same contract. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Supplies deterministic equality and independently chosen hash collisions. */
	private static final class Key{
		Key(final int id_, final int hash_){id=id_; hash=hash_;}
		@Override
		public int hashCode(){return hash;}
		@Override
		public boolean equals(final Object other){return other instanceof Key && ((Key)other).id==id;}
		final int id, hash;
	}
}
