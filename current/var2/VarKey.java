package var2;

/**
 * Simplified key representation of genetic variants for efficient hashing and comparison.
 * Provides a lightweight alternative to full Var objects when only basic variant 
 * identification is needed, without storing complete allele sequences.
 * 
 * Note: This class is not currently thought to be used in the active codebase.
 * 
 * @author Brian Bushnell
 * @contributor Isla Winglet
 */
public class VarKey implements Comparable<VarKey> {
	
	/**
	 * Creates a VarKey from a full Var object.
	 * Extracts essential identifying information while discarding detailed statistics.
	 * 
	 * @param v Source Var object to convert
	 * @return VarKey representing the same variant
	 */
	public static VarKey toVarKey(Var v){
		if(v.type==Var.INS){
			return new VarKey(v.scafnum, v.start, v.allele.length, v.type, v.allele[0]);
		}else if(v.type==Var.DEL){
			return new VarKey(v.scafnum, v.start, v.reflen(), v.type, 0);
		}
		return new VarKey(v.scafnum, v.start, v.reflen(), v.type, v.allele[0]);
	}
	
	/**
	 * Constructs a VarKey with the specified variant parameters.
	 * 
	 * @param scafNum_ Scaffold number (chromosome identifier)
	 * @param start_ Start position on the scaffold
	 * @param length_ Length of the variant (reference length for DEL, allele length for INS)
	 * @param type_ Variant type (SUB, INS, DEL, etc.)
	 * @param allele_ First nucleotide of allele (0 for deletions)
	 */
	public VarKey(int scafNum_, int start_, int length_, int type_, int allele_){
		scafNum=scafNum_;
		start=start_;
		length=length_;
		type=type_;
		allele=allele_;
	}
	
	/**
	 * Computes hash code for efficient storage in hash-based collections.
	 * Uses bit rotation to distribute hash values across all fields.
	 * 
	 * @return Hash code for this VarKey
	 */
	@Override
	public int hashCode(){
		return scafNum^Integer.rotateLeft(start, 4)^Integer.rotateRight(start, 18)^Integer.rotateLeft(type, 8)^Integer.rotateLeft(allele, 12);
	}
	
	/**
	 * Tests equality with another object.
	 * 
	 * @param b Object to compare with
	 * @return true if objects represent the same variant
	 */
	@Override
	public boolean equals(Object b){
		return equals((VarKey)b);
	}
	
	/**
	 * Tests equality with another VarKey.
	 * Two VarKeys are equal if all identifying fields match.
	 * 
	 * @param b VarKey to compare with
	 * @return true if both VarKeys represent the same variant
	 */
	public boolean equals(VarKey b){
		if(b==null){return false;}
		return scafNum==b.scafNum && start==b.start && length==b.length && type==b.type && allele==b.allele;
	}
	
	/**
	 * Compares this VarKey with another for sorting purposes.
	 * Orders by scaffold, then position, then length, then type, then allele.
	 * 
	 * @param b VarKey to compare with
	 * @return Negative, zero, or positive integer for less-than, equal, or greater-than
	 */
	@Override
	public int compareTo(VarKey b){
		if(b==null){return -1;}
		if(scafNum!=b.scafNum){return scafNum-b.scafNum;}
		if(start!=b.start){return start-b.start;}
		if(length!=b.length){return length-b.length;}
		if(type!=b.type){return type-b.type;}
		if(allele!=b.allele){return allele-b.allele;}
		return 0;
	}
	
	/** Scaffold number (chromosome identifier) */
	int scafNum;
	
	/** Start position of the variant on the scaffold */
	int start;
	
	/** Length of the variant (reference length for DEL, allele length for INS) */
	int length;
	
	/** Type of variant (SUB, INS, DEL, etc.) */
	int type;
	
	/** First nucleotide of allele as integer (0 for deletions) */
	int allele;
}