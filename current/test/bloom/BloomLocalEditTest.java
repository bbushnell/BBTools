package test.bloom;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import assemble.LocalEditEngine;
import bloom.BloomFilter;
import bloom.BloomFilterCorrectorWrapper;
import bloom.KCountArray;
import bloom.KCountArray7MTA;
import dna.AminoAcid;
import stream.Read;
import ukmer.Kmer;

/** Small real-count-min-table correction/geometry regression and load diagnostic.
 * Test-only exact string keys deliberately do not approximate the Bloom oracle.
 * The sweep reports adverse results; it does not assert noisy tables are safe.
 * @author Fischl */
public final class BloomLocalEditTest{
	/** Runs six low-load configurations by default, or a larger load diagnostic.
	 * Restores the cell override; the hash seed remains changed. Counters accumulate
	 * across calls in one JVM, so the usual invocation uses a fresh process.
	 * @param args Empty for regression checks, or positional K and cell count */
	public static void main(final String[] args){
		final long previousCells=BloomFilter.OVERRIDE_CELLS;
		try{
			if(args.length==0){
				for(int k:new int[]{31, 62, 63}){
					for(boolean dual:new boolean[]{false, true}){run(k, 2000003, dual, 6, true);}
				}
				System.out.println("BLOOM_LOCAL_EDIT_TEST_OK checks="+checks);
			}else{
				if(args.length!=2){throw new IllegalArgumentException("Load diagnostic: K cells, or no arguments for unit gates.");}
				run(Integer.parseInt(args[0]), Long.parseLong(args[1]), true, 96, false);
			}
		}finally{BloomFilter.OVERRIDE_CELLS=previousCells;}
	}
	/** Loads deterministic references/queries, then compares Bloom and exact correction.
	 * Checks original-array ownership, output quality/header shape and strand symmetry
	 * for every read; strict mode also requires exact repairs and unchanged controls.
	 * @param k Exact short/long k-mer length used by both tables and query engines
	 * @param cells Requested Bloom counter capacity
	 * @param dual Whether Bloom lookup uses the second hash
	 * @param samples Number of generated reads, cycling through six error/control types
	 * @param strict Require the low-load repair oracle; otherwise add noise and report outcomes */
	private static void run(final int k, final long cells, final boolean dual, final int samples, final boolean strict){
		KCountArray7MTA.setSeed(20260911); BloomFilter.OVERRIDE_CELLS=cells;
		final BloomFilter filter=new BloomFilter(null, null, null, k, k, 4, 3, 1, true, false, false, 0.01f, dual);
		final Exact counts=new Exact(); final ArrayList<String> truth=new ArrayList<String>(), queries=new ArrayList<String>();
		final Random random=new Random(2026091102L);
		for(int sample=0; sample<samples; sample++){
			final StringBuilder s=new StringBuilder(); for(int i=0; i<4*k; i++){s.append("ACGT".charAt(random.nextInt(4)));}
			final int p=2*k;
			if((sample&1)==0){for(int i=p-3; i<p+3; i++){s.setCharAt(i, 'A');}}
			final String ref=s.toString(); add(filter, counts, ref, 12, k);
			// Three error classes and three untouched controls, including supported
			// real +/-3 alleles in a higher-depth reference background.
			final int type=sample%6;
			String q=ref, expected=ref;
			if(type==0){
				q=ref.substring(0, p)+(ref.charAt(p)=='A' ? 'C' : 'A')+ref.substring(p+1);
			}else if(type==1){
				q=ref.substring(0, p)+ref.substring(p+1);
			}else if(type==2){
				q=ref.substring(0, p)+ref.charAt(p)+ref.substring(p);
			}else if(type==4){
				q=ref.substring(0, p)+"ACG"+ref.substring(p); expected=q; add(filter, counts, q, 6, k);
			}else if(type==5){
				q=ref.substring(0, p)+ref.substring(p+3); expected=q; add(filter, counts, q, 6, k);
			}
			add(filter, counts, q, 1, k); queries.add(q); truth.add(expected);
		}
		// Same deterministic key population at every capacity. Never insert a
		// second copy merely because the query is evaluated at multiple N values.
		if(!strict){for(int i=0; i<10000; i++){
			final StringBuilder s=new StringBuilder(); for(int j=0; j<k; j++){s.append("ACGT".charAt(random.nextInt(4)));}
			add(filter, counts, s.toString(), i%5==0 ? 8 : 1, k);
		}}
		System.out.println("TABLE k="+k+" cells="+filter.filter.cells+" dual="+dual+" unique_exact="+counts.map.size()+
				" used1="+filter.filter.usedFraction(1)+" used2="+filter.filter.usedFraction(2)+" used4="+filter.filter.usedFraction(4));
		for(int n:new int[]{1, 3, 5}){
			final LocalEditEngine bloom=BloomFilterCorrectorWrapper.makeLocalEditEngine(filter, n);
			final LocalEditEngine exact=new LocalEditEngine(k, counts, n);
			int repaired=0, wrong=0, abstained=0, controlsChanged=0, exactRepaired=0, disagree=0;
			for(int i=0; i<queries.size(); i++){
				final String q=queries.get(i), expected=truth.get(i); final boolean control=q.equals(expected);
				final Read read=read(q), reference=read(q); final byte[] old=read.bases, oldQ=read.quality, before=old.clone(), beforeQ=oldQ.clone();
				final int edits=bloom.correct(read, 1); exact.correct(reference, 1);
				final boolean correct=Arrays.equals(read.bases, expected.getBytes(StandardCharsets.US_ASCII));
				if(control){if(edits>0){controlsChanged++;}}else if(edits==0){abstained++;}else if(correct){repaired++;}else{wrong++;}
				if(!control && Arrays.equals(reference.bases, expected.getBytes(StandardCharsets.US_ASCII))){exactRepaired++;}
				if(!Arrays.equals(read.bases, reference.bases)){disagree++;}
				check(Arrays.equals(old, before) && Arrays.equals(oldQ, beforeQ), "Correction must not mutate the original arrays.");
				check(read.quality.length==read.bases.length && read.id.equals("fixture"), "Output quality/header contract broken.");
				final Read reverse=read(new String(AminoAcid.reverseComplementBases(before), StandardCharsets.US_ASCII));
				for(int j=0; j<beforeQ.length; j++){reverse.quality[j]=beforeQ[beforeQ.length-1-j];}
				final int reverseEdits=bloom.correct(reverse, 1);
				check(reverseEdits==edits && Arrays.equals(AminoAcid.reverseComplementBases(read.bases), reverse.bases), "Bloom correction must commute with reverse complementation.");
				for(int j=0; j<read.quality.length; j++){check(read.quality[j]==reverse.quality[read.quality.length-1-j], "Retained/edit quality placement must commute with reverse complementation.");}
			}
			System.out.println("RESULT k="+k+" n="+n+" repaired="+repaired+" wrong="+wrong+" abstained="+abstained+
					" controls_changed="+controlsChanged+" exact_repaired="+exactRepaired+" disagreed="+disagree+
					" profile_queries="+bloom.profileQueries+" probe_queries="+bloom.probeQueries+" verification_queries="+bloom.verificationQueries);
			if(strict){check(repaired==samples/2 && wrong==0 && controlsChanged==0 && disagree==0, "Low-load Bloom oracle must repair all three error classes and preserve all three controls.");}
		}
	}
	/** Creates an independent read with deterministic position-dependent qualities.
	 * @param s DNA bases copied into a new byte array
	 * @return Unpaired fixture read with qualities cycling from zero through forty */
	private static Read read(final String s){
		final byte[] q=new byte[s.length()]; for(int i=0; i<q.length; i++){q[i]=(byte)(i%41);}
		return new Read(s.getBytes(StandardCharsets.US_ASCII), q, "fixture", 0, false);
	}
	/** Adds each complete window to the real Bloom table and independent exact counts.
	 * Uses primitive canonical keys for short k and the long-key hash pair otherwise.
	 * @param filter Mutable Bloom table constructed with the supplied k
	 * @param exact Mutable exact-key multiplicity oracle
	 * @param s Sequence whose windows are inserted
	 * @param copies Positive multiplicity added for every window
	 * @param k Same exact k used to construct filter and the query engines */
	private static void add(final BloomFilter filter, final Exact exact, final String s, final int copies, final int k){
		assert(copies>0) : "The fixture loader adds known positive multiplicities to both independent count stores.";
		final Kmer key=new Kmer(k);
		for(int i=0; i<s.length(); i++){
			key.addRight(s.charAt(i)); if(key.len()<k){continue;}
			final String identity=Arrays.toString(key.key()); final Integer previous=exact.map.get(identity);
			exact.map.put(identity, (previous==null ? 0 : previous)+copies);
			final long first=k<=31 ? filter.toKey(key.array1()[0], key.array2()[0]) : key.xor();
			if(filter.dualHash){filter.filter.increment(first, k<=31 ? KCountArray.secondaryHash(first) : key.xor2(), copies);}else{for(int j=0; j<copies; j++){filter.filter.increment(first);}}
		}
	}
	/** Exact canonical-key lookup for a small test oracle, independent of Bloom hashing. */
	private static final class Exact implements LocalEditEngine.CountLookup{
		/** Returns the loaded multiplicity of the complete canonical key, or zero if absent.
		 * @param key Complete k-mer to query without mutation
		 * @return Exact accumulated count */
		@Override
		public int count(final Kmer key){final Integer v=map.get(Arrays.toString(key.key())); return v==null ? 0 : v;}
		/** Canonical array text to unsaturated multiplicity; deliberately test-only allocation. */
		final HashMap<String, Integer> map=new HashMap<String, Integer>();
	}
	/** Counts every attempted check and throws on failure in either JVM assertion mode.
	 * @param ok Independent expected condition
	 * @param why Diagnostic describing the violated contract */
	private static void check(final boolean ok, final String why){checks++; if(!ok){throw new AssertionError(why);}}
	/** Number of attempted checks across this process; not reset between configurations. */
	private static long checks;
}
