package stream;

import dna.AminoAcid;
import shared.Tools;

/**
 * Safer MD walker that maps MD:Z tag symbols to a long match string.
 *
 * Differences from MDWalker:
 * - Removes brittle assertions around encountering 'D' in longmatch when CIGAR lacks 'N'.
 *   Deletions are valid regardless of intron encoding and must be handled.
 * - Adds bounds checks around longmatch access to avoid ArrayIndexOutOfBounds.
 * - Keeps original semantics: digits advance matches; '^' enters deletion; letters
 *   mark substitutions; 'I' in longmatch are skipped when counting matches.
 */
public class MDWalker2 {

	MDWalker2(String tag, String cigar_, byte[] longmatch_, SamLine sl_){
		mdTag=tag;
		cigar=cigar_;
		longmatch=longmatch_;
		sl=sl_;
		mdPos=(mdTag.startsWith("MD:Z:") ? 5 : mdTag.startsWith("Z:") ? 2 : 0);

		matchPos=0;
		bpos=0;
		rpos=0;
		sym=0;
		current=0;
		mode=0;

		// Skip leading clipping in longmatch
		while(matchPos<longmatch.length && longmatch[matchPos]=='C'){
			matchPos++;
			bpos++;
		}
	}

	/**
	 * Apply MD tag edits to longmatch in-place. Marks 'S' (or 'N') for subs.
	 */
	void fixMatch(byte[] bases){
		sym=0;
		while(mdPos<mdTag.length()){
			char c=mdTag.charAt(mdPos++);
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
				mode=NORMAL;
				continue;
			}

			int target=matchPos;
			if(current>0){
				target=matchPos+current;
				assert(mode==NORMAL) : mode+", "+current;
				current=0;
			}

			// If prior token was a deletion and the next longmatch symbol is not D,
			// treat upcoming letter as a substitution context.
			if(mode==DEL && matchPos<longmatch.length && longmatch[matchPos]!='D'){
				mode=SUB;
			}

			// Advance across matches until reaching target, skipping inserts ('I')
			// and consuming deletions ('D') by advancing reference only.
			while(matchPos<target || (matchPos<longmatch.length && longmatch[matchPos]=='I')){
				if(matchPos>=longmatch.length){break;}
				byte lm=longmatch[matchPos];
				if(lm=='I'){
					target++; // keep match count aligned ignoring insertions
					bpos++;
					matchPos++;
				}else if(lm=='D'){
					// deletion: advance reference only
					rpos++;
					matchPos++;
				}else{
					rpos++;
					bpos++;
					matchPos++;
				}
			}

			if(c=='^'){
				mode=DEL; // entering deletion block; upcoming letters enumerate ref bases
				continue;
			}

			if(mode==DEL){
				// Within deletion: consume ref and longmatch position; do not touch bpos
				if(matchPos<longmatch.length){
					rpos++;
					matchPos++;
				}else{
					rpos++;
				}
				sym=c;
				continue;
			}

			// Substitution (or continuation thereof)
			if(mode==NORMAL || mode==SUB){
				// If longmatch currently points at a run of deletions, consume them first
				while(matchPos<longmatch.length && longmatch[matchPos]=='D'){
					rpos++;
					matchPos++;
				}
				if(matchPos>=longmatch.length){
					// Nothing left to mark; out-of-bounds MD tag or mismatch vs. longmatch
					break;
				}
				longmatch[matchPos]=(byte)'S';
				if((bases!=null && !AminoAcid.isFullyDefined(bases[bpos])) || !AminoAcid.isFullyDefined(c)){
					longmatch[matchPos]='N';
				}
				mode=SUB;
				bpos++;
				rpos++;
				matchPos++;
				sym=c;
			}else{
				assert(false) : "Unexpected mode "+mode;
			}
		}
	}

	boolean nextSub(){
		sym=0;
		while(mdPos<mdTag.length()){
			char c=mdTag.charAt(mdPos++);
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
				mode=NORMAL;
				continue;
			}
			if(current>0){
				bpos+=current;
				rpos+=current;
				matchPos+=current;
				assert(mode==NORMAL) : mode+", "+current;
				current=0;
			}
			if(c=='^'){
				mode=DEL;
			}else if(mode==DEL){
				rpos++;
				matchPos++;
				sym=c;
			}else if(matchPos<longmatch.length && longmatch[matchPos]=='I'){
				mode=INS;
				bpos++;
				matchPos++;
				sym=c;
			}else if(mode==NORMAL || mode==SUB || mode==INS){
				mode=SUB;
				bpos++;
				rpos++;
				matchPos++;
				sym=c;
				return true;
			}
		}
		return false;
	}

	public int matchPosition(){return matchPos-1;}
	public int basePosition(){return bpos-1;}
	public int refPosition(){return rpos-1;}
	public char symbol(){assert(sym!=0); return sym;}

	/** Position in match string (excluding clipping and insertions) */
	private int matchPos;
	/** Position in read bases (excluding clipping and insertions) */
	private int bpos;
	/** Position in reference bases (excluding clipping) */
	private int rpos;
	private char sym;

	private final String mdTag;
	private final String cigar; //Optional; for debugging
	private final byte[] longmatch;
	private int mdPos;
	private int current;
	private int mode;

	private final SamLine sl; // For debugging

	private static final int NORMAL=0, SUB=1, DEL=2, INS=3;
}

