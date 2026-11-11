package test.assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import assemble.HomopolymerIndelProposal;
import assemble.LocalEditKmerProbe;
import ukmer.Kmer;

/**
 * Checks cached probe reuse against fresh instances of the same probe algorithm.
 * Compares ordered query keys, query counts, depths and unchanged input bytes;
 * this is not an independent edited-sequence oracle or correction accuracy test.
 * Synthetic lookup counts are deterministic functions of encoded keys. Both PACKED
 * layouts are exercised; MASK_CORE is inherited rather than varied by this fixture.
 * Run in an isolated process because the main method temporarily changes global
 * Kmer layout state. Explicit checks remain active without JVM assertions.
 * @author Fischl */
public final class LocalEditProbeReuseTest{

	/** Runs forward substitution/indel scans and reverse direct-indel scans.
	 * Each layout uses effective K derived from requested K31/62/63/94 and requested
	 * maxima of one, three or five windows per edit family. A second pass adds undefined
	 * bases and begins a new
	 * immutable-array scope. Restores PACKED on exit; ignores arguments and retains
	 * the cumulative check counter if invoked repeatedly within one JVM.
	 */
	public static void main(String[] args){
		final boolean packed=Kmer.PACKED;
		try{
			for(boolean layout:new boolean[]{false, true}){
				Kmer.PACKED=layout;
				for(int requested:new int[]{31, 62, 63, 94}){
					for(int n:new int[]{1, 3, 5}){
						final int k=new Kmer(requested).kbig;
						final byte[] bases=new byte[3*k+5];
						final Random random=new Random(551+k);
						for(int i=0; i<bases.length; i++){bases[i]=(byte)"ACGT".charAt(random.nextInt(4));}
						final Trace trace=new Trace();
						final LocalEditKmerProbe cached=new LocalEditKmerProbe(k, trace, n);
						for(int pass=0; pass<2; pass++){
							if(pass==1){bases[k+2]='N'; bases[bases.length-1]='N';}
							cached.beginRead(bases);
							for(int p=0; p<bases.length; p++){
								if(bases[p]!='N'){
									final int w=Math.max(0, Math.min(p-k/2, bases.length-k));
									compare(cached, trace, k, n, bases, w, p, false);
								}
							}
							for(int p=bases.length-1; p>=0; p--){
								if(bases[p]!='N'){
									final int w=Math.max(0, Math.min(p-k/2, bases.length-k));
									compare(cached, trace, k, n, bases, w, p, true);
								}
							}
						}
						budget(k, n);
					}
				}
			}
			System.out.println("LOCAL_EDIT_PROBE_REUSE_TEST_OK checks="+checks);
		}finally{Kmer.PACKED=packed;}
	}

	/** Compares one call sequence against a newly constructed, unscoped probe.
	 * Clears only the cached lookup's trace, preserving its probe cache and failure setting.
	 * Availability flags are not compared separately; array sentinels and depths are.
	 * @param cached Reused probe with any preceding immutable-array scope
	 * @param trace Lookup belonging to cached, with failure injection disabled for this comparison
	 * @param k Effective K shared by both probes
	 * @param n Requested odd maximum number of adjacent candidate windows per edit family
	 * @param bases Borrowed read, required to survive both call sequences unchanged
	 * @param w In-bounds original window start
	 * @param p Called base position inside that window
	 * @param direct True for indels alone, false for substitutions followed by indels
	 */
	private static void compare(final LocalEditKmerProbe cached, final Trace trace, final int k, final int n, final byte[] bases, final int w, final int p, final boolean direct){
		final Trace expectedTrace=new Trace();
		final LocalEditKmerProbe fresh=new LocalEditKmerProbe(k, expectedTrace, n);
		trace.keys.clear(); final byte[] before=bases.clone();
		if(direct){
			fresh.indelsAt(bases, w, p); cached.indelsAt(bases, w, p);
		}else{
			fresh.substitutions(bases, w, p); fresh.indels();
			cached.substitutions(bases, w, p); cached.indels();
		}
		check(cached.queries==fresh.queries && trace.keys.size()==expectedTrace.keys.size(), "Cached query budget differs from fresh probe.");
		for(int i=0; i<trace.keys.size(); i++){check(Arrays.equals(trace.keys.get(i), expectedTrace.keys.get(i)), "Cached key/query order differs from fresh probe.");}
		check(Arrays.equals(cached.substitutionDepth, fresh.substitutionDepth) && Arrays.equals(cached.insertionDepth, fresh.insertionDepth) && cached.deletionDepth==fresh.deletionDepth, "Cached candidate depths differ.");
		check(Arrays.equals(before, bases), "Window reuse must never mutate the input sequence.");
	}

	/** Checks full-load reuse, rolling, same-array invalidation and recovery after lookup errors.
	 * Uses deterministic counters rather than timings. The single-window right-edge case
	 * cannot delete, ensuring that its injected exception occurs during insertion instead.
	 */
	private static void budget(final int k, final int n){
		final byte[] bases=new byte[k+10]; Arrays.fill(bases, (byte)'A');
		final Trace trace=new Trace();
		final LocalEditKmerProbe probe=new LocalEditKmerProbe(k, trace, n);
		probe.beginRead(bases);
		for(int p=k/2; p<k/2+5; p++){probe.substitutions(bases, 3, p); probe.indelsAt(bases, 3, p);}
		check(probe.windowBuilds()==1, "One immutable original window group should need one full build across positions/families.");
		probe.substitutions(bases, 4, k/2);
		check(probe.windowBuilds()==1 && probe.windowRolls()>0, "Adjacent groups must roll, not rebuild.");
		bases[k/2]='C'; probe.beginRead(bases);
		compare(probe, trace, k, n, bases, 3, k/2, false);
		check(probe.windowBuilds()==2, "Beginning a new read scope must invalidate even same-array cached content.");
		trace.failAt=0; boolean failed=false;
		try{probe.substitutions(bases, 3, k/2);}catch(IllegalStateException expected){failed=true;}
		check(failed, "Throwing lookup fixture must execute."); trace.failAt=-1;
		compare(probe, trace, k, n, bases, 3, k/2, false);
		trace.failAt=0; failed=false;
		try{probe.indelsAt(bases, 3, k/2);}catch(IllegalStateException expected){failed=true;}
		check(failed, "Deletion lookup failure fixture must execute."); trace.failAt=-1;
		compare(probe, trace, k, n, bases, 3, k/2, false);
		// The direct single-window last window cannot delete; fail during insertion.
		if(n==1){
			trace.failAt=0; failed=false;
			try{probe.indelsAt(bases, 10, k/2);}catch(IllegalStateException expected){failed=true;}
			check(failed, "Insertion lookup failure fixture must execute."); trace.failAt=-1;
			compare(probe, trace, k, n, bases, 10, k/2, true);
		}
	}

	/** Synthetic count provider with ordered key snapshots and an optional failure countdown. */
	private static final class Trace implements HomopolymerIndelProposal.CountLookup{
		/** Clones each successful query key and returns its hash masked to 0..1023.
		 * An injected failure occurs before recording a key and stays armed until reset.
		 */
		@Override public int count(final Kmer key){
			if(failAt==0){throw new IllegalStateException("Intentional lookup failure.");}
			if(failAt>0){failAt--;}
			final long[] value=key.key().clone(); keys.add(value);
			return Arrays.hashCode(value)&1023;
		}
		/** Owned snapshots in query order; comparisons clear this list between loci. */
		final ArrayList<long[]> keys=new ArrayList<long[]>();
		/** Negative disables failure, zero throws, positive counts successful queries before throwing. */
		int failAt=-1;
	}

	/** Counts every check attempt and throws independently of JVM assertion settings. */
	private static void check(final boolean ok, final String why){checks++; if(!ok){throw new AssertionError(why);}}
	/** Cumulative check attempts across main invocations in this JVM. */
	private static long checks;
}
