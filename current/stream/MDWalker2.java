package stream;

import dna.AminoAcid;
import shared.Tools;

/**
 * Retained reference implementation for applying MD-tag events to an expanded match array.
 * SamLine uses MDWalker; the current Java tree constructs this class only in TestMDWalker.
 * Mutable traversal state and the borrowed match array belong to the caller. Use one
 * traversal strategy per instance; there is no reset, copying or synchronization.
 * Particular longmatch accesses are guarded, but this is not a complete MD/query validator
 * and may leave partial edits. Consistent aligned inputs remain a caller precondition.
 *
 * Historical status note (2026-06-25): this alternative removed restrictions on D entries
 * when CIGAR lacks N and added longmatch bounds checks. That note reported matching output
 * on tested cases and a successful live MDWalker run over 3,000,000 Roche bwa long reads
 * (M+MD, about 9% error), favoring MDWalker's crash-loud checks over this class's partial
 * correction behavior. The report is retained as historical rationale, not a guarantee of
 * equivalence or safety for all inputs. TestMDWalker prints five cases for inspection;
 * it does not prove universal equivalence. This class remains unwired reference code.
 */
public class MDWalker2{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Retains inputs and initializes traversal after an optional MD:Z: or Z: prefix.
	 * Leading C entries advance the match and query cursors, not the reference cursor.
	 * @param tag Nonnull MD text, with either recognized prefix or none
	 * @param cigar_ Optional CIGAR retained but not interpreted here
	 * @param longmatch_ Nonnull caller-owned expanded match array, modified by fixMatch
	 * @param sl_ Optional SAM record retained but not consulted here
	 */
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

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Consumes remaining MD events and changes substitution entries to S or N in place.
	 * Numeric runs advance over match entries while accounting for I/D; a terminal
	 * numeric run remains pending and does not establish final cursor positions.
	 * Undefined reference or supplied query bases select N. CIGAR is not traversed.
	 * Shares state with nextSub; do not interleave the two algorithms or expect rollback.
	 * @param bases Optional aligned query bases, consistent with the expanded match array
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

	/** Advances the legacy substitution iterator without editing longmatch.
	 * Numeric runs ignore intervening I entries (#001); this is not equivalent to fixMatch.
	 * A terminal numeric run remains pending. No iterator callers were found in the Java
	 * tree; the comparison driver uses fixMatch. Both algorithms share all traversal state.
	 * @return true after a reported substitution, false when MD text is exhausted
	 */
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
				//TODO: Probable bug #001 (also MDWalker#001) - direct MD-run advances ignore I entries.
				//MD2A0 over mImS reports match/query position2 instead of3, with reference2.
				//No iterator callers found; retained source concern, not runtime-verified.
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

	/** Returns match cursor minus one, including leading clips and processed events. */
	public int matchPosition(){return matchPos-1;}
	/** Returns query cursor minus one; leading clips and processed insertions count. */
	public int basePosition(){return bpos-1;}
	/** Returns reference cursor minus one, relative to alignment start, not SamLine.pos. */
	public int refPosition(){return rpos-1;}
	/** Returns the most recently stored MD event; asserts if no nonzero symbol is stored. */
	public char symbol(){assert(sym!=0); return sym;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Next expanded-match entry to consume, including leading clips and insertions. */
	private int matchPos;
	/** Next query offset, including leading clips and insertions. */
	private int bpos;
	/** Next reference offset relative to alignment start. */
	private int rpos;
	/** Last stored event character, reset to zero on entry to either walker method. */
	private char sym;

	/** Retained MD text, including any recognized prefix. */
	private final String mdTag;
	/** Optional CIGAR retained for reference; not otherwise used here. */
	private final String cigar;
	/** Caller-owned expanded match array; fixMatch edits substitution entries in place. */
	private final byte[] longmatch;
	/** Index of the next MD character to consume. */
	private int mdPos;
	/** Accumulated numeric run; trailing digits are not applied to cursors. */
	private int current;
	/** Current interpretation of symbolic MD characters. */
	private int mode;

	/** Optional SAM record retained for reference; not otherwise used here. */
	private final SamLine sl;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Parser modes; INS is used only by the legacy iterator. */
	private static final int NORMAL=0, SUB=1, DEL=2, INS=3;
}
