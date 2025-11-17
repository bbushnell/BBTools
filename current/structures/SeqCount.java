package structures;

import java.util.Arrays;

import dna.AminoAcid;
import shared.Tools;

/**
 * Represents a sequence with occurrence counting functionality.
 * Designed for counting identical sequences in genomic datasets with canonical
 * form normalization and efficient hash-based equality testing.
 * Immutable after construction with thread-safe initialization.
 *
 * @author Brian Bushnell
 */
public class SeqCount implements Comparable<SeqCount>, Cloneable {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructs SeqCount from a subsequence of bases.
	 * Creates a copy of the specified range, canonizes it, and computes hash.
	 *
	 * @param s Source sequence array
	 * @param start Starting position (inclusive)
	 * @param stop Ending position (exclusive)
	 */
	public SeqCount(byte[] s, int start, int stop) {
		synchronized(this) {
			bases=Arrays.copyOfRange(s, start, stop);
			synchronized(bases) {
				Tools.canonize(bases);
				hashcode=Tools.hash(bases, 22);
			}
		}
	}
	
	/**
	 * Constructs SeqCount from a complete sequence.
	 * Stores reference to input array, canonizes it, and computes hash.
	 * @param bases_ Sequence bases array
	 */
	public SeqCount(byte[] bases_) {
		synchronized(this) {
			bases=bases_;
			synchronized(bases) {
				Tools.canonize(bases);
				hashcode=Tools.hash(bases, 22);
			}
		}
	}
	
	@Override
	public SeqCount clone() {
		synchronized(this) {
			try {
				SeqCount clone=(SeqCount) super.clone();
//				assert(clone.equals(this));
				return clone;
			} catch (CloneNotSupportedException e) {
				// TODO Auto-generated catch block
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
	
	/**
	 * Returns the occurrence count for this sequence.
	 * Always returns 1 for base SeqCount instances.
	 * @return Count value (always 1)
	 */
	public int count() {return 1;}
	
	/*--------------------------------------------------------------*/
	/*----------------       Inherited Methods      ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final boolean equals(Object s) {
		return equals((SeqCount)s);
	}
	
	/**
	 * Tests equality with another SeqCount.
	 * Uses hash codes for fast inequality testing, then compares bases arrays.
	 * @param s SeqCount to compare against
	 * @return true if sequences are identical, false otherwise
	 */
	public final boolean equals(SeqCount s) {
		if(s==null) {return false;}
//		synchronized(this) {
//			synchronized(s) {
//				synchronized(bases) {
//					synchronized(s.bases) {
		if(hashcode!=s.hashcode) {
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
	public int compareTo(SeqCount s) {
		if(count()!=s.count()) {return count()-s.count();}
		if(bases.length!=s.bases.length) {return bases.length-s.bases.length;}
		return Tools.compare(bases, s.bases);
	}
	
	@Override
	public final int hashCode() {
		return hashcode;
	}
	
	@Override
	public String toString() {
		return getClass()+"@"+super.hashCode()+"="+Integer.toString(count());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Canonized sequence bases array */
	public final byte[] bases;
	/** Pre-computed hash code for the canonized sequence */
	public final int hashcode;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/

	/** Base-to-number conversion table from AminoAcid class */
	public static final byte[] symbolToNumber=AminoAcid.baseToNumber;
	/** Base-to-complement-number conversion table from AminoAcid class */
	public static final byte[] symbolToComplementNumber=AminoAcid.baseToComplementNumber;
	/**
	 * Base-to-number conversion table with 0-based indexing from AminoAcid class
	 */
	public static final byte[] symbolToNumber0=AminoAcid.baseToNumber0;
	
}