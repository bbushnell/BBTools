package prok;

import java.util.ArrayList;

import structures.IntHashMap;
import structures.IntList;

/** Exact sense-oriented DNA membership with case folding and U treated as T.
 * The primitive lookup key only selects a collision chain. Every match compares
 * the complete normalized sequence, including length; keys never prove identity.
 * Retains borrowed read-only member arrays without copying or String conversion.
 * Duplicate sequences collapse presence, not the original construction votes.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaExactSequenceSet{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaExactSequenceSet(){}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Adds one immutable borrowed sequence; returns false for an existing sequence.
	 * Undefined bases fail loudly rather than claiming an ambiguous exact identity.
	 */
	public boolean add(final byte[] bases){
		final int key=key(bases), old=heads.get(key);
		for(int i=old; i>=0; i=next.get(i)){
			if(same(sequences.get(i), bases)){return false;}
		}
		final int index=sequences.size();
		sequences.add(bases);
		next.add(old);
		heads.put(key, index);
		assert(sequences.size()==next.size) : "Every retained exact sequence needs its collision-chain link";
		return true;
	}

	/** Queries the complete sense-oriented sequence, never its reverse complement.
	 * Callers with partial or undefined query sequence must retain UNKNOWN instead.
	 */
	public boolean contains(final byte[] bases){
		final int key=key(bases);
		for(int i=heads.get(key); i>=0; i=next.get(i)){
			if(same(sequences.get(i), bases)){return true;}
		}
		return false;
	}

	public int size(){return sequences.size();}

	/*--------------------------------------------------------------*/
	/*----------------         Static Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Polynomial lookup key; intentional integer wrapping is not sequence proof. */
	private static int key(final byte[] bases){
		if(bases==null || bases.length==0){throw new IllegalArgumentException("Exact gene membership requires a nonempty complete sequence");}
		int key=bases.length;
		for(final byte value:bases){
			final int base=RrnaResourceIO.base(value);
			if(base<0){throw new IllegalArgumentException("Undefined bases cannot establish exact gene membership");}
			key=31*key+base;
		}
		return key;
	}

	/** Full comparison disambiguates any lookup-key collision. */
	private static boolean same(final byte[] a, final byte[] b){
		assert(a!=null && b!=null) : "Exact-sequence collision chains retain validated complete arrays";
		if(a.length!=b.length){return false;}
		for(int i=0; i<a.length; i++){
			if(RrnaResourceIO.base(a[i])!=RrnaResourceIO.base(b[i])){return false;}
		}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	private final IntHashMap heads=new IntHashMap();
	private final IntList next=new IntList();
	private final ArrayList<byte[]> sequences=new ArrayList<byte[]>();
}
