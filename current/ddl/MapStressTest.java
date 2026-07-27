package ddl;

import map.IntHashSet;
import map.LongObjectMap;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * Stress test for LongObjectMap + IntHashSet at the scales DDLBlacklistMaker uses.
 *
 * Motivation: building the 32k blacklist, the SAME kmer reports a different family count and
 * some high-genus-count kmers vanish entirely, depending only on how many total kmers are in the
 * map (2^29 capacity clean, 2^30 corrupt, 2^31 clean).  That should be impossible: per-kmer genus
 * sets are independent of what else is in the map.  This isolates the collections from the
 * sketches, the taxonomy, and the prefilter.
 *
 * For each target size N: insert N distinct keys, each mapped to an IntHashSet holding a
 * deterministic set of values (sized like real genus sets: mostly small, some in the thousands).
 * Then read every key back and verify the set is exactly what was stored.
 *
 * Usage: java ddl.MapStressTest <sizeInMillions> [setSizeCap]
 *
 * @author Noire
 * @date July 14, 2026
 */
public class MapStressTest {

	public static void main(String[] args){
		final long n=(args.length>0 ? Long.parseLong(args[0]) : 100)*1000000L;
		final int cap=(args.length>1 ? Integer.parseInt(args[1]) : 6000);
		Timer t=new Timer();

		System.err.println("Inserting "+n+" keys (set sizes up to "+cap+")...");
		final LongObjectMap<IntHashSet> map=new LongObjectMap<IntHashSet>(IntHashSet.class);
		for(long i=0; i<n; i++){
			final long key=key(i);
			IntHashSet set=map.get(key);
			if(set==null){set=new IntHashSet(4); map.put(key, set);}
			final int size=setSize(i, cap);
			for(int j=0; j<size; j++){set.add(value(i, j));}
		}
		System.err.println("Inserted.  map.size()="+map.size()+"  expected="+n);
		if(map.size()!=n){System.err.println("*** SIZE MISMATCH: lost "+(n-map.size())+" keys ***");}

		System.err.println("Verifying...");
		long missing=0, wrongSize=0, wrongContents=0, ok=0;
		for(long i=0; i<n; i++){
			final long key=key(i);
			final IntHashSet set=map.get(key);
			if(set==null){missing++; continue;}
			final int expected=setSize(i, cap);
			if(set.size()!=expected){
				if(wrongSize<5){
					System.err.println("  key "+Long.toHexString(key)+" (i="+i+"): set.size()="
						+set.size()+" expected "+expected);
				}
				wrongSize++;
				continue;
			}
			boolean good=true;
			for(int j=0; j<expected && good; j++){
				if(!set.contains(value(i, j))){good=false;}
			}
			if(good){ok++;}else{wrongContents++;}
		}

		t.stop();
		System.err.println("\n=== RESULT for N="+n+" ===");
		System.err.println("  correct:         "+ok);
		System.err.println("  missing keys:    "+missing);
		System.err.println("  wrong set size:  "+wrongSize);
		System.err.println("  wrong contents:  "+wrongContents);
		System.err.println(missing+wrongSize+wrongContents==0 ? "  VERDICT: CLEAN" : "  VERDICT: *** CORRUPTION ***");
		System.err.println("Time: \t"+t);
	}

	/** Distinct, well-spread 64-bit keys. */
	private static long key(long i){return Tools.hash64shift(i*0x9E3779B97F4A7C15L)|1L;}

	/** Genus-set-like size distribution: most tiny, a heavy tail into the thousands. */
	private static int setSize(long i, int cap){
		final int m=(int)(i%1000);
		if(m<900){return 1+(int)(i%5);}       //90%: 1-5, like most kmers
		if(m<990){return 5+(int)(i%50);}      //9%:  5-54
		return 5+(int)(i%cap);                 //1%:  up to cap, the hot kmers
	}

	/** Deterministic positive member values, DISTINCT for a given i by construction.
	 * (Hashing each (i,j) into a small space instead would collide by the birthday bound for large
	 * sets, so IntHashSet would rightly dedupe and size() would fall below the add() count --
	 * an artifact of the generator, not the collection.) */
	private static int value(long i, int j){
		final int base=(int)(Tools.hash64shift(i)&0x3FFFFFF);//0..67M, leaves headroom for +j
		return 1+base+j;
	}
}
