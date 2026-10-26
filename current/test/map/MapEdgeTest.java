package test.map;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import map.AbstractBitSet;
import map.ObjectDoubleMap;
import map.ObjectIntMap;
import map.SeqMap;
import map.SeqPos;
import map.SeqPosM;

/** Small independent edge cases discovered during the map package review.
 * @author Shinobu
 */
public final class MapEdgeTest{

	public static void main(String[] args){
		check(args.length==2, "Expected case=name expect=baseline|fixed");
		final String name=args[0].substring("case=".length());
		final boolean fixed=args[1].equals("expect=fixed");
		check(args[0].startsWith("case=") && (fixed || args[1].equals("expect=baseline")), "Invalid fixture arguments");
		if(name.equals("double-count")){doubleCount(fixed);}
		else if(name.equals("render")){render(fixed);}
		else if(name.equals("bit-overflow")){bitOverflow(fixed);}
		else if(name.equals("nan-order")){nanOrder(fixed);}
		else if(name.equals("double-fetch")){doubleFetch(fixed);}
		else if(name.equals("gc-copy")){gcCopy(fixed);}
		else if(name.equals("double-fetch-gc")){doubleFetchGc(fixed);}
		else{throw new IllegalArgumentException(name);}
		System.out.println("PASS\t"+name+"\tfixed="+fixed);
	}

	/** An absent count starts at zero regardless of backing-slot initialization. */
	private static void doubleCount(final boolean fixed){
		final ObjectDoubleMap<String> map=new ObjectDoubleMap<String>(64, String.class);
		check(map.increment("one")== (fixed ? 1 : 0), "Fresh singleton count");
		check(map.increment("fraction", 2.5)==(fixed ? 2.5 : 1.5), "Fresh fractional count");
		map.put("signed", -1);
		check(map.increment("signed", .5)==-.5, "Stored -1 is not absence");
		map.clear();
		check(map.increment("one")==1, "Clear/reuse control");
		final ObjectDoubleMap<String> growth=new ObjectDoubleMap<String>(1, String.class);
		String key=null;
		for(int i=0; i<1024 && key==null; i++){
			final String candidate="key"+i;
			if((shared.Tools.hash32plus(candidate.hashCode())&1)==1){key=candidate;}
		}
		check(key!=null, "Fixture must select an odd mixed hash");
		final double returned=growth.increment(key, 7);
		check(growth.get(key)==(fixed ? 7 : 6), "Growth stored count");
		check(returned==(fixed ? 7 : -1), "Growth return must follow relocated key");
		System.out.println("double-growth\treturned="+returned+"\tstored="+growth.get(key));
	}

	/** The maps' ordinary debug rendering must close the brace it opens. */
	private static void render(final boolean fixed){
		final ObjectIntMap<String> ints=new ObjectIntMap<String>(String.class);
		final ObjectDoubleMap<String> doubles=new ObjectDoubleMap<String>(String.class);
		check(ints.toString().equals(fixed ? "{}" : "{"), "Empty integer rendering");
		check(doubles.toString().equals(fixed ? "{}" : "{"), "Empty double rendering");
		ints.put("a", 3); doubles.put("a", 3);
		check(ints.toString().equals(fixed ? "{a=3}" : "{a=3"), "Singleton integer rendering");
		check(doubles.toString().endsWith("}")==fixed, "Singleton double rendering");
	}

	/** Large positive increments must saturate one two-bit element, never set its neighbor. */
	private static void bitOverflow(final boolean fixed){
		final AbstractBitSet set=AbstractBitSet.make(16, 2);
		set.increment(0);
		set.increment(0, Integer.MAX_VALUE);
		check(set.getCount(0)==(fixed ? 3 : 0), "Target must saturate at three");
		check(set.getCount(15)==(fixed ? 0 : 2), "Overflow must not change another packed element");
		for(int i=1; i<15; i++){check(set.getCount(i)==0, "Untouched neighbor "+i);}
		System.out.println("bit-overflow\tfirst="+set.getCount(0)+"\tlast="+set.getCount(15));
	}

	/** Sorting must preserve comparator antisymmetry even with a NaN score. */
	private static void nanOrder(final boolean fixed){
		final SeqMap map=new SeqMap(1, 0, 0);
		map.add(bytes("A"), 1, Float.NaN);
		map.add(bytes("A"), 1, 1);
		final ArrayList<SeqPos> entries=map.get(0);
		check(entries!=null && entries.size()==2, "Both scored records must share the A bucket");
		final int ab=Integer.signum(entries.get(0).compareTo(entries.get(1)));
		final int ba=Integer.signum(entries.get(1).compareTo(entries.get(0)));
		check((ab==-ba)==fixed, "NaN ordering must be antisymmetric: "+ab+", "+ba);
		System.out.println("nan-order\tab="+ab+"\tba="+ba);
	}

	/** The second region may contain a valid all-T kmer when the first ends in N. */
	private static void doubleFetch(final boolean fixed){
		final SeqMap map=new SeqMap(3, 0, 0);
		map.add(bytes("TTT"), 1, 1);
		final byte[] query=bytes("NNTTTT");
		final ArrayList<SeqPosM> found=map.doubleFetch(query, 0, 2, 3, 5, 3, 0, 0, 0, 1, false);
		final int count=found==null ? 0 : found.size();
		check(count==(fixed ? 1 : 0), "Valid second-region kmer must be queried: "+count);
		if(found!=null){check(found.get(0).pos()==0, "Second-region candidate remains first-region relative");}
		System.out.println("double-fetch\tcandidates="+count);
	}

	/** Copying a position must replace its GC metadata along with its sequence. */
	private static void gcCopy(final boolean fixed){
		final SeqMap map=new SeqMap(1, 0, 0);
		map.add(bytes("C"), 1, 1);
		final SeqPos source=map.get(1).get(0);
		final SeqPosM copy=new SeqPosM(bytes("A"), 0);
		check(source.gc==1 && copy.gc==0, "Fixture must start with distinct GC fractions");
		copy.setFrom(source);
		check(copy.gc==(fixed ? 1 : 0), "GC must follow the copied sequence");
		check(copy.seq==source.seq && copy.count==source.count && copy.hashcode==source.hashcode, "Other copy metadata");
	}

	/** A shared kmer still needs the second region's distinct GC filter. */
	private static void doubleFetchGc(final boolean fixed){
		final SeqMap map=new SeqMap(3, 0, 0);
		map.add(bytes("NNNNNNNNNNNNNNNNNNNNCCC"), 1, 1);
		final byte[] query=bytes("AAAAAAAAAAAAAAAAAAAACCCGGGGGGGGGGGGGGGGGGGGCCC");
		final ArrayList<SeqPosM> found=map.doubleFetch(query, 0, 22, 23, 45, 23, 0, 0, 0, 1, false);
		final int count=found==null ? 0 : found.size();
		check(count==(fixed ? 1 : 0), "Second-region GC filter must be considered for a shared kmer: "+count);
		if(found!=null){check(found.get(0).pos()==0, "Shared-key candidate stays first-region relative");}
	}

	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static void check(final boolean ok, final String message){
		if(!ok){throw new AssertionError(message);}
	}
}
