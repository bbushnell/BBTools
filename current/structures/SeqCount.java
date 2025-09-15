package structures;

import java.util.Arrays;

import dna.AminoAcid;
import shared.Tools;

public class SeqCount implements Comparable<SeqCount>, Cloneable {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	public SeqCount(byte[] s, int start, int stop) {
		synchronized(this) {
			bases=Arrays.copyOfRange(s, start, stop);
			synchronized(bases) {
				Tools.canonize(bases);
				hashcode=Tools.hash(bases, 22);
			}
		}
	}
	
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
	
	public int count() {return 1;}
	
	/*--------------------------------------------------------------*/
	/*----------------       Inherited Methods      ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final boolean equals(Object s) {
		return equals((SeqCount)s);
	}
	
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
	
	public final byte[] bases;
	public final int hashcode;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/

	public static final byte[] symbolToNumber=AminoAcid.baseToNumber;
	public static final byte[] symbolToComplementNumber=AminoAcid.baseToComplementNumber;
	public static final byte[] symbolToNumber0=AminoAcid.baseToNumber0;
	
}