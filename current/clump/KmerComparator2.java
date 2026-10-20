package clump;

import java.util.Comparator;

import shared.Tools;
import stream.Read;

/**
 * Base comparator for optical duplicate passes on reads that already have
 * {@link ReadKey} objects attached.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public abstract class KmerComparator2 implements Comparator<Read>{

	/**
	 * Compares attached ReadKeys, then read sequence, quality, then read ID.
	 *
	 * @param a first read to compare
	 * @param b second read to compare
	 * @return negative if a &lt; b, positive if a &gt; b, zero if equal
	 */
	@Override
	public final int compare(Read a, Read b){
		// a.obj/b.obj are ReadKeys populated during hashing
		// (KmerComparator.hash). The cast is unguarded because the X/Y
		// subclass comparators only run after keys are set.
		final ReadKey keyA=(ReadKey)a.obj;
		final ReadKey keyB=(ReadKey)b.obj;

		// The X/Y comparator runs first; sequence, quality, and ID break ties.
		int x=compare(keyA, keyB);
		if(x==0){
			x=compareSequence(a, b, 0);
		}
		return x==0 ? a.id.compareTo(b.id) : x;
	}
	
	/** Compares two attached read keys on the X or Y optical axis. */
	public abstract int compare(ReadKey a, ReadKey b);
	
	/**
	 * Breaks key ties by bases, then mate bases, then quality-sum.
	 *
	 * @param a first read to compare
	 * @param b second read to compare
	 * @param depth unused recursion vestige
	 * @return sequence comparison, or quality comparison if bases tie
	 */
	// The mate comparison is guarded only by a.mate!=null. Callers must pass
	// reads with a consistent mate convention; mixed paired/unpaired reads
	// can still dereference b.mate here.
	public static final int compareSequence(Read a, Read b, int depth){
		int x=compareSequence(a.bases, b.bases);
		if(x!=0){return x;}
		if(a.mate!=null){x=compareSequence(a.mate.bases, b.mate.bases);}
		if(x!=0){return x;}
		return compareQuality(a, b);
	}
	
	//CLEVER [verified]: total order = longer-first (L<below> b.length-a.length),
	//then byte-ascending (a[i]-b[i], byte->int promotion bounded -255..255,
	//no overflow). The length compare short-circuits the O(n) byte loop for the
	//common unequal-length case. Deterministic -> reproducible within-clump
	//ordering. null sorts last.
	/** Orders null arrays last, then longer arrays first, then signed byte values
	 * in ascending order. Equal lengths and bytes compare equal. */
	public static final int compareSequence(final byte[] a, final byte[] b){
		if(a==null || b==null){
			if(a==null && b!=null){return 1;}
			if(a!=null && b==null){return -1;}
			return 0;
		}
		if(a.length!=b.length){
			return b.length-a.length;
		}
		for(int i=0, lim=a.length; i<lim; i++){
			int x=a[i]-b[i];
			if(x!=0){return x;}
		}
		return 0;
	}
	
	//Final ordering fallback among byte-identical reads: higher first-read
	//total quality sorts EARLIER (qb-qa; qa>qb -> negative -> a first), before
	//the final ID tie-break. Duplicate retention is decided later by Clump's
	//expected-error, containment, and scan policy. Asymmetric null handling
	//BOTH ways: if only a is null it returns 0 (non-antisymmetric vs the
	//reverse call); if only b is null it NPEs at sumInt(b.quality).
	/**
	 * Breaks byte-identical sequence ties by total quality.
	 *
	 * @param a first read to compare
	 * @param b second read to compare
	 * @return negative when a has higher quality, positive when b has higher
	 *         quality, or zero when a has no quality scores
	 */
	// Callers must pass reads with a consistent quality convention; mixed
	// quality/no-quality reads can give an asymmetric zero or dereference
	// b.quality here.
	public static final int compareQuality(Read a, Read b){
		if(a.quality==null){return 0;}
		int qa=Tools.sumInt(a.quality);
		int qb=Tools.sumInt(b.quality);
		return qb-qa;
	}
	
}
