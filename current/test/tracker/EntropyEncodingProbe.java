package test.tracker;

import java.lang.reflect.Array;
import java.lang.reflect.Field;

import tracker.EntropyTracker;

/** JT011 arithmetic evidence and explicitly bounded real amino constructors.
 * Arithmetic rows do not allocate k-mer arrays or claim constructor execution.
 * @author Brian Bushnell, Jean */
public final class EntropyEncodingProbe {

	public static void main(String[] args) throws ReflectiveOperationException{
		assert(args.length==0) : "JT011 encoding fixtures are fixed; no arbitrary constructor sizes are accepted.";
		System.out.println("kind\tcase\texpected\tobserved\tagrees");
		arithmetic("amino_k2", 5, 2, "10,10,1024,1023,true");
		arithmetic("amino_k6", 5, 6, "30,30,1073741824,1073741823,true");
		arithmetic("amino_k7", 5, 7, "35,3,8,7,false");
		arithmetic("amino_k13", 5, 13, "65,1,2,1,false");
		arithmetic("dna_k15", 2, 15, "30,30,1073741824,1073741823,true");
		arithmetic("dna_k16", 2, 16, "32,0,1,0,false");
		construct("default_amino_control", 2, true, false, "32,1024,26,26,25");
		construct("amino_k1_control", 1, false, false, "32,32,27,27,25");
		construct("amino_k2_control", 2, false, false, "32,1024,26,26,25");
		construct("amino_k7_reject", 7, false, true, "32,8,21,21,25");
		construct("amino_k13_reject", 13, false, true, "32,2,15,15,25");
	}

	/** Integer-only reproduction of the source expressions; expected tuples were counted separately. */
	private static void arithmetic(String name, int bitsPerBase, int k, String expected){
		final int bits=bitsPerBase*k;
		final String observed=bits+","+(bits&31)+","+(1<<bits)+","+(~((-1)<<bits))+","+(bits<31);
		show("source_arithmetic", name, expected, observed, expected.equals(observed));
	}

	private static void construct(String name, int k, boolean defaults, boolean reject, String expectedSizes)
			throws ReflectiveOperationException{
		//Whitelist before calling production: never allocate the enormous amino-k6 table.
		if(k!=1 && k!=2 && k!=7 && k!=13){throw new IllegalArgumentException("Unapproved constructor k="+k);}
		EntropyTracker tracker;
		try{
			tracker=defaults ? new EntropyTracker(true, -1, true) : new EntropyTracker(k, 25, true);
		}catch(AssertionError e){
			System.err.println("JT011 encoding case="+name);
			e.printStackTrace(System.err);
			show("constructor", name, reject ? "REJECT" : "ACCEPT", "REJECT:java.lang.AssertionError", reject);
			return;
		}
		final String sizes=length(tracker, "baseCounts")+","+length(tracker, "counts")+","+
			length(tracker, "countCounts")+","+length(tracker, "entropy")+","+length(tracker, "baseRingBuffer");
		if(!sizes.equals(expectedSizes)){
			throw new AssertionError("Actual constructor allocations differ from reviewed counts: "+name+" "+sizes);
		}
		String observed="ACCEPT:k="+tracker.k()+",window="+tracker.windowBases()+",arrays="+sizes;
		boolean agrees=false;
		if(!reject){
			for(int i=0; i<26; i++){tracker.add((byte)'A');}
			final float entropy=tracker.calcEntropy(), monomer=tracker.calcMaxMonomerFraction();
			observed+=",entropy="+entropy+",monomer="+monomer;
			agrees=tracker.k()==k && tracker.windowBases()==25 && Math.abs(entropy)<1e-6 && Math.abs(monomer-1)<1e-6;
		}
		show("constructor", name, reject ? "REJECT" : "ACCEPT:entropy0,monomer1", observed, agrees);
	}

	/** Test-only observation of the five actual instance arrays, without any production helper. */
	private static int length(EntropyTracker tracker, String name) throws ReflectiveOperationException{
		final Field field=EntropyTracker.class.getDeclaredField(name);
		field.setAccessible(true);
		return Array.getLength(field.get(tracker));
	}

	private static void show(String kind, String name, String expected, String observed, boolean agrees){
		System.out.println(kind+"\t"+name+"\t"+expected+"\t"+observed+"\t"+agrees);
	}
}
