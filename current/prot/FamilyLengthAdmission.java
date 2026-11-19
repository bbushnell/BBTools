package prot;

/** Immutable raw-hit floors from exact retained-member counts, never query or consensus length.
 * The first released form requires4 hits below65aa,5 below100aa, and6 otherwise.
 * @author Keqing
 */
public final class FamilyLengthAdmission{

	FamilyLengthAdmission(final FamilyReference[] families_){
		if(families_==null || families_.length==0){throw new IllegalArgumentException("Family-length admission needs a nonempty live roster");}
		families=families_; floors=new byte[families.length];
		int min=6, max=4;
		for(int i=0;i<families.length;i++){
			final FamilyReference family=families[i];
			if(family==null){throw new IllegalArgumentException("Missing live family for length admission at ordinal "+i);}
			final int floor=rejectFloor(family.members,family.residues);
			floors[i]=(byte)floor; min=Math.min(min,floor); max=Math.max(max,floor);
		}
		minimum=min; maximum=max;
	}

	/** Compares exact integers at the band boundaries; no rounded mean enters admission.
	 * Supported retained peptides are40..200aa. Reject absent/overflowing counts loudly. */
	public static int rejectFloor(final long members,final long residues){
		if(members<=0 || members>Long.MAX_VALUE/200 || residues<40*members || residues>200*members){
			throw new IllegalArgumentException("Family admission requires positive retained membership and exact residue totals within40..200aa per member: "+members+"/"+residues);
		}
		return residues<65*members ? 4 : residues<100*members ? 5 : 6;
	}

	int floor(final int ordinal){
		assert(ordinal>=0 && ordinal<floors.length) : "Admission ordinals must match the immutable live reference roster";
		return floors[ordinal];
	}
	boolean boundTo(final FamilyReference[] candidate){return candidate==families;}

	final int minimum,maximum;
	private final byte[] floors;
	private final FamilyReference[] families;
}
