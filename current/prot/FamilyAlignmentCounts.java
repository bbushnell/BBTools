package prot;

import java.util.Arrays;

/** Reusable development-only counts from one existing alignment traceback.
 * Reference counts include canonical residues opposite mismatches and deletions;
 * match counts include only exact canonical residue pairs. Not a scoring policy.
 * @author Sayu
 */
public final class FamilyAlignmentCounts{

	public void clear(){Arrays.fill(reference,0);Arrays.fill(matches,0);measured=false;}
	public void paired(byte query,byte ref){
		assert(!measured) : "A completed trace must be cleared before counting another candidate";
		if(Blosum62.isStandard(ref)){reference[ref]++;if(query==ref){matches[ref]++;}}
	}
	public void deleted(byte ref){
		assert(!measured) : "A completed trace must be cleared before counting another candidate";
		if(Blosum62.isStandard(ref)){reference[ref]++;}
	}
	public void finish(){
		assert(!measured) : "One traceback has exactly one completion event";
		measured=true;
	}
	public boolean measured(){return measured;}
	/** Copy a completed borrowed trace into independently owned worker scratch. */
	public void copyFrom(FamilyAlignmentCounts other){
		if(other==null || !other.measured){throw new IllegalArgumentException("Only measured traces can be retained for primary-view selection");}
		System.arraycopy(other.reference,0,reference,0,20);
		System.arraycopy(other.matches,0,matches,0,20);measured=true;
	}
	public int referenceCount(int code){check(code);return reference[code];}
	public int matchCount(int code){check(code);return matches[code];}
	private void check(int code){
		if(!measured || code<0 || code>=20){
			throw new IllegalStateException("Residue counts require a measured trace and a canonical code: "+code);
		}
	}
	private final int[] reference=new int[20],matches=new int[20];
	private boolean measured;
}
