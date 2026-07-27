package ddl;

import map.IntHashSet;
import shared.Timer;

/**
 * Perfect-recall test for IntHashSet (Brian's method): feed it X values from a seeded RNG, then
 * RESET the seed, regenerate the same X values, and assert every one is still retrievable.
 *
 * Two regimes, because the real workload does both:
 *   A) FRESH sets of growing size (like the per-kmer genus sets: 1 .. thousands of members).
 *   B) ONE set, reused across many rounds via clear() (like `promoted` in DDLBlacklistMaker, which
 *      is a single IntHashSet cleared once per kmer per rank -- tens of millions of clear() calls).
 * Regime B is the one that has never been validated, and IntHashSet's `invalid` sentinel is a
 * RANDOM NEGATIVE int per instance with a resetInvalid() rehash path, so reuse is worth proving.
 *
 * Duplicates are expected from an RNG, so the invariant is not size==X; it is:
 *   every generated value is contained, AND size == number of DISTINCT generated values.
 *
 * Usage: java ddl.SetRecallTest [rounds] [maxSetSize]
 *
 * @author Noire
 * @date July 14, 2026
 */
public class SetRecallTest {

	public static void main(String[] args){
		final int rounds=(args.length>0 ? Integer.parseInt(args[0]) : 200000);
		final int maxSize=(args.length>1 ? Integer.parseInt(args[1]) : 6000);
		Timer t=new Timer();

		long badFresh=0, badReused=0, badSizeFresh=0, badSizeReused=0;
		final IntHashSet reused=new IntHashSet(256);//Regime B: the SAME set, cleared every round

		for(int r=0; r<rounds; r++){
			final int n=1+(int)(mix(r)%maxSize);
			final long seed=0x5DEECE66DL^r;

			//Regime A: a fresh set, like a per-kmer genus set
			final IntHashSet fresh=new IntHashSet(4);
			long rng=seed;
			for(int i=0; i<n; i++){rng=next(rng); fresh.add(val(rng));}

			//Regime B: the reused set, like `promoted`
			reused.clear();
			rng=seed;
			for(int i=0; i<n; i++){rng=next(rng); reused.add(val(rng));}

			//Reset the RNG and demand perfect recall of every value we fed in
			rng=seed;
			for(int i=0; i<n; i++){
				rng=next(rng);
				final int v=val(rng);
				if(!fresh.contains(v)){badFresh++;}
				if(!reused.contains(v)){badReused++;}
			}

			//Distinct-count invariant: size must equal the number of distinct values generated
			rng=seed;
			final IntHashSet truth=new IntHashSet(4);
			for(int i=0; i<n; i++){rng=next(rng); truth.add(val(rng));}
			if(fresh.size()!=truth.size()){badSizeFresh++;}
			if(reused.size()!=truth.size()){badSizeReused++;}

			if(r>0 && r%50000==0){System.err.println("  round "+r+"...");}
		}

		t.stop();
		System.err.println("\n=== IntHashSet recall over "+rounds+" rounds (sets up to "+maxSize+") ===");
		System.err.println("  fresh  set: missing values="+badFresh+"  wrong size="+badSizeFresh);
		System.err.println("  reused set: missing values="+badReused+"  wrong size="+badSizeReused);
		System.err.println((badFresh+badReused+badSizeFresh+badSizeReused)==0
			? "  VERDICT: CLEAN" : "  VERDICT: *** CORRUPTION ***");
		System.err.println("Time: \t"+t);
	}

	/** Positive values, wide range, so duplicates are possible but not dominant. */
	private static int val(long rng){return 1+(int)((rng>>>17)&0x7FFFFFF);}
	private static long next(long rng){return rng*0x5DEECE66DL+0xBL;}
	private static long mix(long x){x^=x>>>33; x*=0xff51afd7ed558ccdL; x^=x>>>33; return x&0x7FFFFFFF;}
}
