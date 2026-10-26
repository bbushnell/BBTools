package assemble;

import dna.AminoAcid;
import ukmer.HollowKmerLookup;
import ukmer.KmerTableSetU;

/** Worker-local proposals for a substitution entering a Hollow left anchor.
 * Counts are observations, not proof of sequence truth. This class never edits
 * a sequence; the extension caller owns the transaction and rollback.
 * @author Brian Bushnell, Fischl */
public final class HollowAnchorRescue {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrows immutable Hollow counts without creating another evidence table. */
	public HollowAnchorRescue(final KmerTableSetU source){
		this(source, false);
	}

	/** Optionally rejects supported competing ingress alleles in the same count matrix. */
	public HollowAnchorRescue(final KmerTableSetU source, final boolean rejectCompetitor_){
		lookup=new HollowKmerLookup(source);
		k=source.kbig;
		rejectCompetitor=rejectCompetitor_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns a replacement ingress base, or -1; records whether a trial must roll back.
	 * Call only at an extension stopping point with a defined terminal K window.
	 * A supported ordinary successor must never be replaced by this rescue.
	 * The rare stop-path packs 64 bases, then queries four ordinary successors;
	 * only a strong collapse pays for the full16-choice matrix. */
	public int propose(final byte[] bases, final int length, final int minExtend){
		if(bases==null || length<k || length>bases.length || minExtend<1){
			throw new IllegalArgumentException("Hollow rescue needs a complete terminal window and positive extension threshold.");
		}
		collapsed=false;
		final long left=arm(bases, length-k), right=arm(bases, length-32);
		final int ingress=AminoAcid.baseToNumber[bases[length-k+32]];
		assert(ingress>=0) : "rightmostKmer must validate the full terminal window before a rescue can edit its interior.";
		final int context=lookup.count(left, right);
		if(context<FLOOR){
			// An initial weak context cannot propose a correction. After a trial,
			// however, the same weakness must reject it rather than silently count
			// as confirmation. A deletion-seed failure retained a false correction
			// after support fell14->2; the old early return bypassed rollback.
			collapsed=true;
			return -1;
		}
		final long nextLeft=(left<<2)|ingress, nextRight=right<<2;
		int ordinary=0;
		for(int base=0; base<4; base++){
			ordinary=Math.max(ordinary, lookup.count(nextLeft, nextRight|base));
		}
		if(ordinary>=minExtend || ordinary>1 || context<RATIO*(long)Math.max(1, ordinary)){return -1;}
		collapsed=true;
		lookup.fillRightCounts16(left, right, choices);
		int best=-1, high=0, second=0;
		for(int i=0; i<choices.length; i++){
			final int count=choices[i];
			if(count>high){second=high; high=count; best=i;}
			else if(count>second){second=count;}
		}
		// The winner must explain the observed depth, not merely win a weak or
		// repeat-pooled competition. A tie is rejected regardless of ACGT order.
		//TODO: Repeat copies can still satisfy these ratios and replace a correct
		// seed base; this occurred with mixed-error reads from Dactylococcopsis
		// at K192. Strong terminal support did not distinguish the repeat copy.
		// Keep recovery opt-in until original-seed support is independently guarded.
		if(best<0 || (best>>>2)==ingress || high<minExtend || high<FLOOR
				|| high<RATIO*(long)Math.max(1, second) || 2L*high<context || high>2L*context){return -1;}
		if(rejectCompetitor){
			// A repeat-majority allele can pass the ratio while the true allele
			// still has independent support. Exclude the entire winning ingress
			// row: alternative NEXT bases do not imply a competing seed edit.
			for(int i=0; i<choices.length; i++){
				if((i>>>2)!=(best>>>2) && choices[i]>=FLOOR){return -1;}
			}
		}
		return best>>>2;
	}

	/** Encodes only the two counted arms; never reconstructs the masked interior. */
	private static long arm(final byte[] bases, final int start){
		assert(start>=0 && start+32<=bases.length) : "A count anchor must contain exactly32 existing sequence bases.";
		long value=0;
		for(int i=start; i<start+32; i++){
			final int base=AminoAcid.baseToNumber[bases[i]];
			if(base<0){throw new IllegalArgumentException("Undefined base in Hollow rescue anchor at "+i);}
			value=(value<<2)|base;
		}
		return value;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** True for inadequate terminal support or a strong ordinary-support collapse.
	 * Initial proposals still return -1; after a tentative correction this flag
	 * requires restoring its original byte and extension checkpoint. */
	public boolean collapsed;
	/** Committed correction in the current oriented extension call; -1 after rollback. */
	public int correctedPosition=-1;
	private final HollowKmerLookup lookup;
	private final int k;
	private final boolean rejectCompetitor;
	private final int[] choices=new int[16];
	private static final int FLOOR=4, RATIO=4;
}
