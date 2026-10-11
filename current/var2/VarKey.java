package var2;

/**
 * Legacy five-field variant key that retains only one allele byte.
 * Distinct multi-base alleles can therefore produce equal keys. Fields are mutable
 * within the package; callers using keys in collections must keep their values stable.
 * No active caller was found in the current Java sources during the package review.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @contributor Sucrose
 */
public class VarKey implements Comparable<VarKey>{

	/**
	 * Copies scaffold, start and type, with type-dependent length and allele fields.
	 * Insertions use allele length and first allele byte; deletions use reference
	 * length and allele zero; other types use reference length and first allele byte.
	 * Does not preserve the complete allele identity or any observation statistics.
	 * @param v Nonnull variant; all non-deletion paths require a nonempty allele
	 * @return New reduced key
	 */
	public static VarKey toVarKey(Var v){
		if(v.type==Var.INS){
			return new VarKey(v.scafnum, v.start, v.allele.length, v.type, v.allele[0]);
		}else if(v.type==Var.DEL){
			return new VarKey(v.scafnum, v.start, v.reflen(), v.type, 0);
		}
		return new VarKey(v.scafnum, v.start, v.reflen(), v.type, v.allele[0]);
	}

	/** Stores the five values without validation or normalization.
	 * @param scafNum_ Scaffold number
	 * @param start_ Variant start coordinate
	 * @param length_ Allele length for insertions, reference length otherwise
	 * @param type_ Variant type code
	 * @param allele_ Promoted first allele byte, or zero for deletions */
	public VarKey(int scafNum_, int start_, int length_, int type_, int allele_){
		scafNum=scafNum_;
		start=start_;
		length=length_;
		type=type_;
		allele=allele_;
	}

	/**
	 * Combines scaffold, start, type and allele using XOR and rotations.
	 * Length participates in equality but is not included in this hash.
	 * @return Hash of the current field values
	 */
	@Override
	public int hashCode(){
		return scafNum^Integer.rotateLeft(start, 4)^Integer.rotateRight(start, 18)^Integer.rotateLeft(type, 8)^Integer.rotateLeft(allele, 12);
	}

	/** Delegates through an unchecked VarKey cast; null compares unequal.
	 * @param b Object to compare
	 * @return Equality of all five fields for a VarKey, false for null
	 * @throws ClassCastException If b is a nonnull object of another type */
	@Override
	public boolean equals(Object b){
		//TODO: Probable bug - unrelated object types throw instead of comparing unequal (VAR2-004).
		return equals((VarKey)b);
	}

	/**
	 * Tests equality with another VarKey.
	 * Two VarKeys are equal if all identifying fields match.
	 * @param b VarKey to compare with
	 * @return True if all stored fields match, not necessarily the complete allele sequence
	 */
	public boolean equals(VarKey b){
		if(b==null){return false;}
		return scafNum==b.scafNum && start==b.start && length==b.length && type==b.type && allele==b.allele;
	}

	/**
	 * Compares this VarKey with another for sorting purposes.
	 * Orders by scaffold, then position, then length, then type, then allele.
	 * Uses int subtraction, so ordering assumes those differences do not overflow.
	 * @param b Key to compare, or null
	 * @return First nonzero field difference, zero for equal fields, or -1 for null
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

	/** Scaffold number. */
	int scafNum;

	/** Variant start coordinate. */
	int start;

	/** Insertion allele length or other variant's reference length. */
	int length;

	/** Variant type code. */
	int type;

	/** First allele byte promoted to int, or zero for deletion keys. */
	int allele;
}
