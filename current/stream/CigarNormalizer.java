package stream;

import dna.AminoAcid;
import shared.Tools;

/**
 * Deterministic left-alignment (canonicalization) of a SAM record's CIGAR against a reference.
 *
 * <p>This is <b>alignment</b> normalization, not <b>call</b> (variant) normalization: indels are
 * slid to their leftmost reference position while <i>conserving identity</i> (a base that matched
 * still matches; a mismatch stays a mismatch; a '=' is never turned into an 'X'), but they are kept
 * <b>strictly interior</b> to the read. An indel is never rolled into the outermost {@link #MARGIN}
 * bases of the alignment, and never across a soft/hard clip or an intron (N) — a CIGAR that begins
 * with a deletion, or an indel jammed against a read end, is nonsense.
 *
 * <p>Work is done on the <b>long match string</b> (one symbol per alignment column, always
 * reference-forward because it is derived from the CIGAR, which the SAM spec mandates be
 * reference-forward). "Left" on that string is therefore reference-left for both strands. The only
 * strand-sensitive input is the read sequence used to test insertion rolls: under {@code FLIP_ON_LOAD}
 * a minus-strand read's stored {@code seq} is in native orientation, so it is reverse-complemented to
 * a reference-forward view here. Deletions need only reference bases and are strand-safe.
 *
 * <p>Coalescing of SW-fragmented gaps (e.g. {@code 1M3D1M3D1M5D1M1D1M3D -> 1M15D4M}) is an emergent
 * consequence of left-aligning each indel through equal bases until adjacent gaps meet; it is not a
 * separate merge, and it can never convert a match to a mismatch. Adjacent same-type ops are merged
 * by {@link SamLine#toCigar14} when the CIGAR is rebuilt.
 *
 * <p>Requires a loaded reference (the scaffold's bases), because rolling can need reference context
 * beyond the read's own footprint.
 *
 * @author Furina
 */
public class CigarNormalizer {

	/** Outermost number of alignment columns at each end that indels may not be pushed into. */
	public static final int MARGIN=3;

	/**
	 * Left-aligns the indels of {@code sl}'s CIGAR against {@code refBases}, rebuilding the CIGAR in
	 * =/X form if anything moved. The CIGAR must already be in =/X form (M resolved) — i.e. run eqx first.
	 *
	 * @param sl A mapped SAM record with a =/X CIGAR.
	 * @param refBases The full reference scaffold bases for {@code sl}'s rname (reference-forward).
	 * @return true if the CIGAR changed.
	 */
	public static boolean normalize(SamLine sl, byte[] refBases){
		if(sl==null || refBases==null || sl.cigar==null || !sl.mapped()){return false;}

		final byte[] shortMatch=sl.toShortMatch(false); //M already resolved to m/S via ref/MD
		if(shortMatch==null){return false;}
		final byte[] lm=Read.toLongMatchString(shortMatch); //one symbol per column, reference-forward
		if(lm==null || lm.length==0){return false;}

		final byte[] readFwd=refForwardRead(sl); //read bases in reference-forward orientation
		final int refStart=sl.pos-1;             //0-based reference coordinate of the first aligned column

		final boolean changed=leftAlign(lm, readFwd, refBases, refStart);
		if(!changed){return false;}

		final int start=sl.pos-1;
		final int stop=start+Read.calcMatchLength(lm)-1;
		sl.setCigar(SamLine.toCigar14(lm, start, stop, Integer.MAX_VALUE, sl.seq));
		return true;
	}

	/** Returns the read bases in reference-forward orientation (reverse-complemented for a flipped minus read). */
	private static byte[] refForwardRead(SamLine sl){
		final byte[] seq=sl.seq;
		if(seq==null){return null;}
		if(SamLine.FLIP_ON_LOAD && sl.mapped() && sl.strand()==shared.Shared.MINUS){
			return AminoAcid.reverseComplementBases(seq); //returns a new array; sl.seq is untouched
		}
		return seq;
	}

	/**
	 * Repeatedly rolls every indel run as far left as it will go (to the fixpoint), conserving identity
	 * and respecting the interior/clip/intron boundaries. Operates in place on {@code lm}.
	 * @return true if any column moved.
	 */
	private static boolean leftAlign(byte[] lm, byte[] readFwd, byte[] ref, int refStart){
		boolean anyChange=false;
		boolean pass=true;
		while(pass){
			pass=false;
			for(int c=0; c<lm.length; c++){
				final byte sym=lm[c];
				if(sym=='D'){
					if(rollDeletionLeft(lm, c, ref, refStart, readFwd)){pass=true; anyChange=true;}
				}else if(isInsertion(sym)){
					if(rollInsertionLeft(lm, c, readFwd, ref, refStart)){pass=true; anyChange=true;}
				}
			}
		}
		return anyChange;
	}

	/**
	 * Rolls the deletion run beginning at column {@code k} one position left, if the column to its left
	 * is a true match, the move stays outside the margin, and the reference base entering on the
	 * left equals the reference base leaving on the right (so the displaced base's match status is preserved).
	 */
	private static boolean rollDeletionLeft(byte[] lm, int k, byte[] ref, int refStart, byte[] readFwd){
		if(k==0 || !isTrueMatch(lm[k-1])){return false;} //roll only PAST a true match (never through a mismatch)
		int k2=k; while(k2<lm.length && lm[k2]=='D'){k2++;} //run is [k, k2)

		final int refBefore=refConsumedBefore(lm, k-1);     //ref index consumed by the left match column
		if(refBefore<MARGIN){return false;}                 //outermost MARGIN columns are sacrosanct
		final int refLeftCoord=refStart+refBefore;          //reference position of the left match column
		final int refLastDelCoord=refStart+refConsumedBefore(lm, k2-1); //ref position of the last D column

		if(refLeftCoord<0 || refLastDelCoord>=ref.length){return false;}
		final byte rL=upper(ref[refLeftCoord]), rR=upper(ref[refLastDelCoord]);
		if(!isACGT(rL) || rL!=rR){return false;} //roll only through equal, defined bases

		//Roll: [M, D, ..., D] -> [D, ..., D, M]. The match column moves to k2-1; re-derive its m/S vs the new ref base.
		final int readIdx=readConsumedBefore(lm, k-1);
		final byte readBase=(readFwd==null ? 0 : upper(readFwd[readIdx]));
		lm[k-1]='D';
		lm[k2-1]=(readFwd!=null && readBase==rR) ? (byte)'m' : (byte)'S';
		return true;
	}

	/**
	 * Rolls the insertion run beginning at column {@code k} one position left, if the column to its left
	 * is a true match, the move stays outside the margin, and the read base entering on the left
	 * equals the read base leaving on the right (so the displaced base's match status is preserved).
	 */
	private static boolean rollInsertionLeft(byte[] lm, int k, byte[] readFwd, byte[] ref, int refStart){
		if(k==0 || !isTrueMatch(lm[k-1]) || readFwd==null){return false;} //roll only PAST a true match
		int k2=k; while(k2<lm.length && isInsertion(lm[k2])){k2++;} //run is [k, k2)

		final int readBefore=readConsumedBefore(lm, k-1);
		final int refBefore=refConsumedBefore(lm, k-1);
		if(readBefore<MARGIN || refBefore<MARGIN){return false;} //margin on both axes

		final int readLeftIdx=readBefore;                     //read base of the left match column
		final int readLastInsIdx=readConsumedBefore(lm, k2-1);//read base of the last inserted column
		if(readLeftIdx>=readFwd.length || readLastInsIdx>=readFwd.length){return false;}
		final byte qL=upper(readFwd[readLeftIdx]), qR=upper(readFwd[readLastInsIdx]);
		if(!isACGT(qL) || qL!=qR){return false;}

		//Roll: [M, I, ..., I] -> [I, ..., I, M]. The match column moves to k2-1; its ref position is unchanged.
		final int refCoord=refStart+refBefore;
		final byte refBase=(refCoord>=0 && refCoord<ref.length) ? upper(ref[refCoord]) : 0;
		lm[k-1]='I';
		lm[k2-1]=(qR==refBase) ? (byte)'m' : (byte)'S';
		return true;
	}

	/** Number of reference-consuming columns strictly before column {@code c}. */
	private static int refConsumedBefore(byte[] lm, int c){
		int n=0;
		for(int i=0; i<c; i++){if(consumesRef(lm[i])){n++;}}
		return n;
	}

	/** Number of read-consuming columns strictly before column {@code c}. */
	private static int readConsumedBefore(byte[] lm, int c){
		int n=0;
		for(int i=0; i<c; i++){if(consumesRead(lm[i])){n++;}}
		return n;
	}

	/** True MATCH column (read base equals reference) — the only kind an indel may be rolled past. */
	private static boolean isTrueMatch(byte s){return s=='m' || s=='s';}
	/** Aligned column: match OR substitution (consumes both read and reference). */
	private static boolean isAligned(byte s){return s=='m' || s=='s' || s=='S' || s=='V';}
	/** Insertion column (consumes read only). */
	private static boolean isInsertion(byte s){return s=='I' || s=='X' || s=='Y';}
	/** True if the column consumes a reference base (match/sub/deletion/intron). */
	private static boolean consumesRef(byte s){return isAligned(s) || s=='D' || s=='N' || s=='B';}
	/** True if the column consumes a read base (match/sub/insertion/soft-clip). */
	private static boolean consumesRead(byte s){return isAligned(s) || isInsertion(s) || s=='C';}

	private static byte upper(byte b){return (b>='a' && b<='z') ? (byte)(b-32) : b;}
	private static boolean isACGT(byte b){return b=='A' || b=='C' || b=='G' || b=='T';}

}
