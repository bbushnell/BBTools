package stream;

import dna.AminoAcid;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Repeated reference-left shifting of indels in a SAM record's alignment.
 * This changes alignment placement within the CIGAR, not variant-call representation
 * or reference assembly. Intended inputs have resolved =/X operations, read bases
 * consistent with those operations and the supplied reference, and a storage orientation
 * consistent with SamLine.FLIP_ON_LOAD.
 *
 * <p>Each roll crosses an immediately preceding true-match symbol in the converted
 * long match string, requiring equal defined bases at the affected boundary. The
 * margin is a consumed-prefix bound: reference bases for deletions, both read and
 * reference bases for insertions. There is no separate right-end margin check or
 * repair of already terminal indels. Repeated allowed rolls can bring gaps together;
 * SamLine.toCigar14 groups adjacent operations when rebuilding the CIGAR.
 *
 * <p>The match-string conversion does not retain all original CIGAR distinctions:
 * SamLine.toShortMatch maps N to D and omits H. N and D are treated equivalently
 * for alignment normalization; retaining the original N/D annotation is optional.
 * Original terminal H tokens are retained around the rebuilt CIGAR, outside any soft clips.
 * Rebuilding also applies SamLine.INTRON_LIMIT when choosing D versus N for a gap.
 * Coordinates and optional tags are not refreshed by this class.
 * @author Furina, Shinobu
 */
public class CigarNormalizer{

	/** Minimum consumed reference prefix before a roll; insertion rolls also require this read prefix. */
	public static final int MARGIN=3;

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Rolls eligible indels left and assigns a rebuilt CIGAR if any roll occurred.
	 * Resolve M operations first. Supply read bases consistent with the resolved alignment,
	 * and keep FLIP_ON_LOAD consistent with the record's stored sequence orientation.
	 * The conversion/rebuild limitations described on this class also apply here.
	 * @param sl Mapped SAM record with a resolved =/X CIGAR and available read sequence
	 * @param refBases Full reference scaffold bases for the record, in reference-forward orientation
	 * @return true after at least one roll and CIGAR reassignment; false for null record/reference,
	 * missing CIGAR, unmapped record, unavailable/empty match string or no eligible roll
	 */
	public static boolean normalize(SamLine sl, byte[] refBases){
		if(sl==null || refBases==null || sl.cigar==null || !sl.mapped()){return false;}

		//N/D annotation may change during normalization; preserving that label is not required.
		//FIXED [stream/CigarNormalizer#001]: the codec omits H; rebuilding 2H5=1I5= used
		//to yield 3=1I7= and lose the terminal clip. Retain original H text around the rebuilt
		//alignment without changing the shared codec or the accepted N/D normalization policy.
		final String originalCigar=sl.cigar;
		final byte[] shortMatch=sl.toShortMatch(false);//M already resolved to m/S via ref/MD
		if(shortMatch==null){return false;}
		final byte[] lm=Read.toLongMatchString(shortMatch);//one symbol per column, reference-forward
		if(lm==null || lm.length==0){return false;}

		final byte[] readFwd=refForwardRead(sl);//read bases in reference-forward orientation
		final int refStart=sl.pos-1;//0-based reference coordinate of the first aligned column

		final boolean changed=leftAlign(lm, readFwd, refBases, refStart);
		if(!changed){return false;}

		final int start=sl.pos-1;
		final int stop=start+Read.calcMatchLength(lm)-1;
		final String rebuilt=SamLine.toCigar14(lm, start, stop, Integer.MAX_VALUE, sl.seq);
		sl.setCigar(retainHardClips(originalCigar, rebuilt));
		return true;
	}

	/**
	 * Restores original terminal hard-clip tokens around a rebuilt alignment.
	 * Copies token text without converting clip counts; hard clips stay outside soft clips.
	 * @param original Nonnull original valid CIGAR, captured before match-string conversion
	 * @param rebuilt Rebuilt CIGAR, or null if conversion produced no CIGAR
	 * @return Combined CIGAR, or rebuilt unchanged when null or no terminal H is present
	 */
	private static String retainHardClips(final String original, final String rebuilt){
		assert(original!=null) : "normalize checks for a CIGAR before retaining its terminal clips";
		if(rebuilt==null){return null;}
		final int len=original.length();
		int left=0;
		for(int end=0; end<len;){
			while(end<len && Tools.isDigit(original.charAt(end))){end++;}
			if(end<len && original.charAt(end)=='H'){left=++end;}
			else{break;}
		}
		int right=len;
		while(right>left && original.charAt(right-1)=='H'){
			int start=right-1;
			while(start>left && Tools.isDigit(original.charAt(start-1))){start--;}
			right=start;
		}
		if(left==0 && right==len){return rebuilt;}
		final ByteBuilder bb=new ByteBuilder(left+rebuilt.length()+len-right);
		//ByteBuilder's substring append writes nullBytes for an empty span; skip absent clips.
		if(left>0){bb.append(original, 0, left);}
		bb.append(rebuilt);
		if(right<len){bb.append(original, right, len);}
		return bb.toString();
	}

	/** Returns a new reverse-complement array for a flipped mapped minus record; otherwise borrows seq.
	 * Returns null when the record has no sequence; does not modify the stored sequence. */
	private static byte[] refForwardRead(SamLine sl){
		final byte[] seq=sl.seq;
		if(seq==null){return null;}
		if(SamLine.FLIP_ON_LOAD && sl.mapped() && sl.strand()==shared.Shared.MINUS){
			return AminoAcid.reverseComplementBases(seq);//returns a new array; sl.seq is untouched
		}
		return seq;
	}

	/**
	 * Repeats permitted one-position rolls until a pass makes no change, modifying lm in place.
	 * Tests operate on converted symbols and consumed prefixes, not original CIGAR boundaries.
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
	 * is a true-match symbol, the consumed reference prefix meets MARGIN and the reference
	 * boundary bases are equal A/C/G/T. Reclassifies the displaced column using readFwd.
	 * @return true when one roll was applied, false when its checks prevented movement
	 */
	private static boolean rollDeletionLeft(byte[] lm, int k, byte[] ref, int refStart, byte[] readFwd){
		if(k==0 || !isTrueMatch(lm[k-1])){return false;} //roll only PAST a true match (never through a mismatch)
		int k2=k; while(k2<lm.length && lm[k2]=='D'){k2++;} //run is [k, k2)

		final int refBefore=refConsumedBefore(lm, k-1);//ref index consumed by the left match column
		if(refBefore<MARGIN){return false;}                 //left consumed-reference prefix, not a right-end or column-index bound
		final int refLeftCoord=refStart+refBefore;//reference position of the left match column
		final int refLastDelCoord=refStart+refConsumedBefore(lm, k2-1);//ref position of the last D column

		if(refLeftCoord<0 || refLastDelCoord>=ref.length){return false;}
		final byte rL=upper(ref[refLeftCoord]), rR=upper(ref[refLastDelCoord]);
		if(!isACGT(rL) || rL!=rR){return false;} //roll only through equal, defined bases

		//Roll: [M, D, ..., D] -> [D, ..., D, M]. The match column moves to k2-1; re-derive its m/S vs the new ref base.
		final int readIdx=readConsumedBefore(lm, k-1);
		//TODO: Probable bug [stream/CigarNormalizer#002]: without read bases, the assignment
		//below emits S after moving a true-match column. Missing-sequence handling is unverified;
		//do not treat this path as an identity-preservation guarantee.
		final byte readBase=(readFwd==null ? 0 : upper(readFwd[readIdx]));
		lm[k-1]='D';
		lm[k2-1]=(readFwd!=null && readBase==rR) ? (byte)'m' : (byte)'S';
		return true;
	}

	/**
	 * Rolls the insertion run beginning at column {@code k} one position left, if the column to its left
	 * is a true-match symbol, both consumed prefixes meet MARGIN and the read boundary bases
	 * are equal A/C/G/T. Reclassifies the displaced column against the reference base.
	 * @return true when one roll was applied, false when its checks prevented movement
	 */
	private static boolean rollInsertionLeft(byte[] lm, int k, byte[] readFwd, byte[] ref, int refStart){
		if(k==0 || !isTrueMatch(lm[k-1]) || readFwd==null){return false;} //roll only PAST a true match
		int k2=k; while(k2<lm.length && isInsertion(lm[k2])){k2++;} //run is [k, k2)

		final int readBefore=readConsumedBefore(lm, k-1);
		final int refBefore=refConsumedBefore(lm, k-1);
		if(readBefore<MARGIN || refBefore<MARGIN){return false;} //margin on both axes

		final int readLeftIdx=readBefore;//read base of the left match column
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

	/** Recognizes internal true-match labels m/s; rolls do not cross any other preceding symbol. */
	private static boolean isTrueMatch(byte s){return s=='m' || s=='s';}
	/** Aligned column: match OR substitution (consumes both read and reference). */
	private static boolean isAligned(byte s){return s=='m' || s=='s' || s=='S' || s=='V';}
	/** Insertion column (consumes read only). */
	private static boolean isInsertion(byte s){return s=='I' || s=='X' || s=='Y';}
	/** Recognizes reference-consuming internal symbols: aligned bases, D and N/B placeholders. */
	private static boolean consumesRef(byte s){return isAligned(s) || s=='D' || s=='N' || s=='B';}
	/** True if the column consumes a read base (match/sub/insertion/soft-clip). */
	private static boolean consumesRead(byte s){return isAligned(s) || isInsertion(s) || s=='C';}

	/** Uppercases ASCII letters without changing other bytes. */
	private static byte upper(byte b){return (b>='a' && b<='z') ? (byte)(b-32) : b;}
	/** Tests membership in the uppercase unambiguous DNA alphabet. */
	private static boolean isACGT(byte b){return b=='A' || b=='C' || b=='G' || b=='T';}

}
