package test.prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import idaligner.QuantumAligner;
import idaligner.ScrabbleAligner;
import map.LongHashSet;
import prok.NcrnaAlignmentCache;
import prok.NcrnaScavenger;
import prok.Orf;
import shared.Shared;
import structures.ByteBuilder;

/** Exact-kernel and real-scavenger comparisons for bounded alignment reuse.
 * These finite cases are not whole-caller or biological accuracy validation.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaAlignmentCacheTest{

	public static void main(final String[] args){
		check(args.length==0, "This finite suite takes no arguments");
		final boolean saved=Shared.SIMD;
		try{
			Shared.SIMD=false;
			run();
			if(simd.Vector.simd256){
				Shared.SIMD=true;
				run();
			}
		}finally{Shared.SIMD=saved;}
		System.out.println("NCRNA_ALIGNMENT_CACHE_TEST_PASS");
	}

	private static void run(){
		final byte[][] library={sequence(96, 31), sequence(83, 73)};
		final byte[] window=join(sequence(17, 5), library[0], sequence(29, 101));
		final NcrnaAlignmentCache cache=new NcrnaAlignmentCache(library, 2, 512);
		compareQuantum(cache, library, window, 0, 0);
		compareQuantum(cache, library, window.clone(), 0, 0);
		check(cache.hits()==1 && cache.fills()==1, "Equal bytes in different arrays must reuse one result");
		compareQuantum(cache, library, window, 0, 23);
		check(cache.fills()==2, "A changed incoming Quantum offset must execute the kernel");
		compareQuantum(cache, library, window, 0, 23);
		compareQuantum(cache, library, window, 1, 0);
		compareScrabble(cache, library, window, 0, false);
		compareScrabble(cache, library, window, 0, true);
		compareScrabble(cache, library, window, 0, false);
		compareScrabble(cache, library, window, 0, true);

		final long beforeMutation=cache.fills();
		window[27]=window[27]=='A' ? (byte)'C' : (byte)'A';
		compareQuantum(cache, library, window, 0, 23);
		check(cache.fills()==beforeMutation+1, "Mutating the caller's input must not mutate a retained key");
		final long beforeClear=cache.fills();
		cache.clear();
		compareQuantum(cache, library, window, 0, 23);
		check(cache.fills()==beforeClear+1, "Clearing an input must invalidate retained results");

		final NcrnaAlignmentCache tiny=new NcrnaAlignmentCache(library, 1, 512);
		compareQuantum(tiny, library, window, 0, 0);
		compareQuantum(tiny, library, library[1], 1, 0);
		compareQuantum(tiny, library, window, 0, 0);
		check(tiny.hits()==0 && tiny.evictions()==2, "Evicted windows must recompute rather than reuse stale model slots");
		final NcrnaAlignmentCache disabled=new NcrnaAlignmentCache(library, 0, 512);
		compareQuantum(disabled, library, window, 0, 0);
		compareQuantum(disabled, library, window, 0, 0);
		check(disabled.hits()==0 && disabled.fills()==2, "Zero retention is the real-kernel control");
		final NcrnaAlignmentCache shortOnly=new NcrnaAlignmentCache(library, 1, 50);
		compareQuantum(shortOnly, library, window, 0, 0);
		compareQuantum(shortOnly, library, window, 0, 0);
		check(shortOnly.hits()==0, "An oversized window must not evade the retention bound");
		collision(library);
		scavenger(library);
	}

	private static void compareQuantum(final NcrnaAlignmentCache cache, final byte[][] library,
			final byte[] sequence, final int model, final int offset){
		final int[] expected={offset, 0, 0, 0};
		final int[] actual=expected.clone();
		final float identity=QuantumAligner.alignStatic(library[model], sequence, expected);
		check(Float.floatToIntBits(identity)==Float.floatToIntBits(cache.quantum(model, sequence, actual))
			&& Arrays.equals(expected, actual), "Cached Quantum must preserve identity and all four positions");
		Arrays.fill(actual, Integer.MIN_VALUE);
	}

	private static void compareScrabble(final NcrnaAlignmentCache cache, final byte[][] library,
			final byte[] sequence, final int model, final boolean positioned){
		final int[] expected=positioned ? new int[2] : null;
		final int[] actual=positioned ? new int[2] : null;
		final float identity=positioned ? ScrabbleAligner.alignStatic(library[model], sequence, expected)
			: ScrabbleAligner.alignStatic(sequence, library[model], null);
		check(Float.floatToIntBits(identity)==Float.floatToIntBits(cache.scrabble(model, sequence, actual))
			&& Arrays.equals(expected, actual), "Scrabble direction and positioned/unpositioned modes must remain separate");
	}

	/** Aa and BB have equal polynomial hashes but different bytes; surrounding
	 * sequence preserves that collision. This is a storage fixture, not a gene.
	 */
	private static void collision(final byte[][] library){
		final byte[] a=join(library[0], "Aa".getBytes(StandardCharsets.US_ASCII), library[1]);
		final byte[] b=a.clone();
		b[library[0].length]='B';
		b[library[0].length+1]='B';
		check(Arrays.hashCode(a)==Arrays.hashCode(b) && !Arrays.equals(a, b), "Fixture must contain a real hash collision");
		final NcrnaAlignmentCache cache=new NcrnaAlignmentCache(library, 2, 512);
		compareQuantum(cache, library, a, 0, 0);
		compareQuantum(cache, library, b, 0, 0);
		compareQuantum(cache, library, a, 0, 0);
		check(cache.fills()==2 && cache.hits()==1, "Hash equality must never replace exact sequence equality");
	}

	private static void scavenger(final byte[][] library){
		final byte[] bases=join(sequence(40, 11), library[0], sequence(180, 29), library[0], sequence(40, 13));
		final int[] hits={48, 324};
		long seed=0;
		for(int i=0; i<17; i++){seed=(seed<<2)|dna.AminoAcid.baseToNumber[library[0][i]];}
		final LongHashSet seeds=new LongHashSet();
		seeds.add(seed);
		final NcrnaAlignmentCache cache=new NcrnaAlignmentCache(library, 32, 512);
		for(int cap=2; cap>=1; cap--){
			final NcrnaScavenger ordinary=caller(library, seeds, cap);
			final NcrnaScavenger reused=caller(library, seeds, cap);
			reused.setAlignmentCache(cache);
			final ArrayList<Orf> expected=ordinary.scavenge("finite", bases, 0, new ArrayList<int[]>(), hits);
			final ArrayList<Orf> actual=reused.scavenge("finite", bases, 0, new ArrayList<int[]>(), hits);
			check(expected.size()==2, "Finite separated exact loci must both reach the scavenger output");
			check(Arrays.equals(gff(expected), gff(actual)), "Reused kernels must preserve actual emitted scavenger records");
			check(ordinary.alignmentCount()==reused.alignmentCount(), "Reuse must preserve logical alignment attempts");
		}
		check(cache.hits()>0, "The second policy must exercise result reuse rather than only misses");
	}

	private static NcrnaScavenger caller(final byte[][] library, final LongHashSet seeds, final int cap){
		final NcrnaScavenger caller=new NcrnaScavenger(library, null, null, seeds, 17, 20, 80,
			7, cap, true, 0f, 0f, 0f, 1);
		caller.family="tpp";
		return caller;
	}

	private static byte[] gff(final ArrayList<Orf> orfs){
		final ByteBuilder out=new ByteBuilder();
		for(Orf orf:orfs){orf.appendGff(out).nl();}
		return out.toBytes();
	}

	private static byte[] sequence(final int length, final int seed){
		final byte[] out=new byte[length];
		int value=seed;
		for(int i=0; i<length; i++){
			value=value*1664525+1013904223;
			out[i]=(byte)"ACGT".charAt((value>>>24)&3);
		}
		return out;
	}

	private static byte[] join(final byte[] a, final byte[] b, final byte[] c){
		final byte[] out=Arrays.copyOf(a, a.length+b.length+c.length);
		System.arraycopy(b, 0, out, a.length, b.length);
		System.arraycopy(c, 0, out, a.length+b.length, c.length);
		return out;
	}

	private static byte[] join(final byte[] a, final byte[] b, final byte[] c, final byte[] d, final byte[] e){
		return join(join(a, b, c), d, e);
	}
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
