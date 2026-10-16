package structures;

import java.util.Arrays;

import dna.AminoAcid;
import shared.Tools;

/**
 * Canonicalized sequence key with a fixed public base-array reference and
 * cached hash.
 * <p>
 * Constructors canonize the stored sequence to the smaller strand under
 * Tools.canonical's extended base-number ordering. One constructor copies a
 * half-open slice; the whole-array constructor
 * retains and mutates its input. Equality first checks the cached hash, then
 * compares the current base bytes. Instances are safe as hash keys only while
 * the public base array is left unchanged. Ordering uses count, length, score
 * and bases; equality and hashing remain sequence-only.
 *
 * @author Brian Bushnell
 */
public class SeqCount implements Comparable<SeqCount>, Cloneable{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Copies a sequence slice and stores it in canonical orientation.
	 * @param s Source bases.
	 * @param start First included source index.
	 * @param stop First excluded source index.
	 */
	public SeqCount(byte[] s, int start, int stop){
		synchronized(this){
			bases=Arrays.copyOfRange(s, start, stop);
			synchronized(bases){
				Tools.canonize(bases);
				hashcode=Tools.hash(bases, 22);
			}
		}
	}

	/**
	 * Stores the supplied array in canonical orientation.
	 * The array is retained and canonized in place.
	 * @param bases_ Sequence bytes to retain and canonize.
	 */
	public SeqCount(byte[] bases_){
		synchronized(this){
			bases=bases_;
			synchronized(bases){
				Tools.canonize(bases);
				hashcode=Tools.hash(bases, 22);
			}
		}
	}

	/** Returns a shallow clone that shares the final base array. */
	@Override
	public SeqCount clone(){
		synchronized(this){
			try{
				SeqCount clone=(SeqCount) super.clone();
//				assert(clone.equals(this));
				return clone;
			}catch(CloneNotSupportedException e){
				e.printStackTrace();
				return null;
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

//	public void add(SeqCount s) {
//		throw new RuntimeException("This class is immutable.");
//	}
//
//	public void increment(int x) {
//		throw new RuntimeException("This class is immutable.");
//	}

	/** Returns the implicit occurrence count for this sequence key. */
	public int count(){return 1;}

	/** Returns the implicit comparison score for this sequence key. */
	public float score(){return -1;}

	/*--------------------------------------------------------------*/
	/*----------------       Inherited Methods      ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final boolean equals(Object s){
		if(!(s instanceof SeqCount)){return false;}
		return equals((SeqCount)s);
	}

	/** Tests whether another SeqCount stores the same canonicalized bases. */
	public final boolean equals(SeqCount s){
		if(s==null){return false;}
//		synchronized(this) {
//			synchronized(s) {
//				synchronized(bases) {
//					synchronized(s.bases) {
		if(hashcode!=s.hashcode){
//			assert(!Tools.equals(bases, s.bases)) : new String(bases)+", "+new String(s.bases)+", "+hashcode+", "+s.hashcode;
			return false;
		}
		return Tools.equals(bases, s.bases);
//					}
//				}
//			}
//		}
	}

	@Override
	public int compareTo(SeqCount s){
		int x=Integer.compare(count(), s.count());
		if(x!=0){return x;}
		x=Integer.compare(bases.length, s.bases.length);
		if(x!=0){return x;}
		x=compareScores(score(), s.score());
		if(x!=0){return x;}
		return Tools.compare(bases, s.bases);
	}

	/** Compares scores with finite values below NaN and +0/-0 as ties. */
	private static int compareScores(final float a, final float b){
		if(a==b){return 0;}
		if(a>b){return 1;}
		if(a<b){return -1;}
		return Float.isNaN(a) ? (Float.isNaN(b) ? 0 : 1) : -1;
	}

	/** Returns the precomputed hash code for the canonicalized bases. */
	@Override
	public final int hashCode(){
		return hashcode;
	}

	/** Returns the runtime class, identity hash and current count. */
	@Override
	public String toString(){
		return getClass()+"@"+super.hashCode()+"="+Integer.toString(count());
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Canonicalized sequence bytes retained by this key. */
	public final byte[] bases;

	/** Precomputed hash code for the canonicalized bases. */
	public final int hashcode;

	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/

	/** Base-to-number lookup. */
	public static final byte[] symbolToNumber=AminoAcid.baseToNumber;

	/** Base-to-complement-number lookup. */
	public static final byte[] symbolToComplementNumber=AminoAcid.baseToComplementNumber;

	/** Base-to-number lookup mapping A/C/G/T/U to 0-3 and other symbols to 0. */
	public static final byte[] symbolToNumber0=AminoAcid.baseToNumber0;

}
