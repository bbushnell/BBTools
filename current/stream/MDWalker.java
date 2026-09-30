package stream;

import dna.AminoAcid;
import shared.KillSwitch;
import shared.Tools;

/**
 * Consumes SAM MD text against an expanded BBTools match array.
 * The live correction path marks substitutions in the caller-owned array. Match, query
 * and reference cursors are relative to the supplied alignment, not absolute coordinates.
 * CIGAR text is retained but its operations are not parsed here.
 * Each walker owns mutable traversal state with no reset; use a fresh instance per walk.
 * fixMatch and the legacy nextSub iterator share that state and have different gap handling.
 * Trailing numeric runs are accumulated without advancing cursors, so getters do not
 * report final alignment lengths.
 *
 * @author Brian Bushnell
 * @date May 5, 2016
 */
public class MDWalker{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Retains the supplied references and positions cursors after leading clip entries.
	 * Accepts MD text with no prefix, {@code MD:Z:}, or {@code Z:}. Leading C entries
	 * advance match/query cursors only; the reference cursor starts at zero.
	 * @param tag Nonnull MD text
	 * @param cigar_ Optional CIGAR text; its N-presence scan is unused by correction
	 * @param longmatch_ Nonnull expanded match array, retained and modified by fixMatch
	 * @param sl_ Optional SAM record used only in diagnostics
	 */
	MDWalker(String tag, String cigar_, byte[] longmatch_, SamLine sl_){//SamLine is just for debugging
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

		while(matchPos<longmatch.length && longmatch[matchPos]=='C'){
			matchPos++;
			bpos++;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Consumes remaining MD text and marks described substitutions S or N in place.
	 * N is used when the MD symbol or corresponding supplied read base is undefined.
	 * Insertions advance query/match cursors; deletions advance reference/match cursors.
	 * A final numeric run remains pending, rather than advancing to the alignment end.
	 * Inputs must describe a consistent alignment; this is not a complete MD validator.
	 * Detected length inconsistencies use the existing terminating diagnostics under -ea.
	 * @param bases Optional query bases for distinguishing defined substitutions from no-calls
	 */
	void fixMatch(byte[] bases){
		final boolean cigarContainsN=(cigar!=null && cigar.indexOf('N')>=0);
		sym=0;
		while(mdPos<mdTag.length()){
			char c=mdTag.charAt(mdPos);
			mdPos++;

			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
				mode=NORMAL;
			}else{
				int matchPos2=matchPos;
				if(current>0){
					matchPos2=matchPos+current;
					assert(mode==NORMAL) : mode+", "+current;
					current=0;
				}

				//Fixes subs after dels getting ignored
				if(mode==DEL && (matchPos<longmatch.length && longmatch[matchPos]!='D')){
					mode=SUB;
				}

				while(matchPos<matchPos2 || (matchPos<longmatch.length && longmatch[matchPos]=='I')){
					//Crash-loud under -ea on an inconsistent MD/longmatch (matchPos overran a match run); -da degrades (break).
					if(matchPos>=longmatch.length){
						assert(false) : KillSwitch.assertDie("MD/longmatch inconsistent: matchPos overran the match string in a match run.\n"+(sl==null ? "" : sl.toString()+"\n")+new String(longmatch));
						break;
					}
					if(longmatch[matchPos]=='I'){
						matchPos2++;
						bpos++;
						matchPos++;
					}else if(longmatch[matchPos]=='D'){
						// Advance reference for deletions regardless of CIGAR N presence
						while(matchPos<longmatch.length && longmatch[matchPos]=='D'){
							rpos++;
							matchPos++;
						}
					}else{
						rpos++;
						bpos++;
						matchPos++;
					}
				}

				if(c=='^'){
					mode=DEL;
				}else if(mode==DEL){
					//Crash-loud under -ea if the deletion runs past the match string; -da advances reference only (best-effort).
					if(matchPos<longmatch.length){
						rpos++;
						matchPos++;
					}else{
						assert(false) : KillSwitch.assertDie("MD/longmatch inconsistent: deletion ran past the match string.\n"+(sl==null ? "" : sl.toString()+"\n")+new String(longmatch));
						rpos++;
					}
					sym=c;
				}else if(mode==NORMAL || mode==SUB){
					// Consume any pending deletions at current position
					while(matchPos<longmatch.length && longmatch[matchPos]=='D'){
						rpos++;
						matchPos++;
					}
					//Crash-loud under -ea if the MD names a sub past the end of the match string; -da degrades (break).
					if(matchPos>=longmatch.length){
						assert(false) : KillSwitch.assertDie("MD names a substitution past the end of the match string (length-inconsistent MD vs CIGAR).\n"+(sl==null ? "" : sl.toString()+"\n")+new String(longmatch));
						break;
					}
					longmatch[matchPos]=(byte)'S';
					//Crash-loud under -ea if the MD names a sub past the end of read bases; -da treats it as base-defined (best-effort).
					final boolean basesOver=(bases!=null && bpos>=bases.length);
					assert(!basesOver) : KillSwitch.assertDie("MD names a substitution past the end of read bases (length-inconsistent MD vs SEQ): bpos="+bpos+", bases.length="+(bases==null ? -1 : bases.length));
					if((bases!=null && !basesOver && !AminoAcid.isFullyDefined(bases[bpos])) || !AminoAcid.isFullyDefined(c)){longmatch[matchPos]='N';}
					mode=SUB;
					bpos++;
					rpos++;
					matchPos++;
					sym=c;
				}else{
					assert(false);
				}

			}

		}
	}

	/**
	 * Advances the legacy substitution iterator without modifying the match array.
	 * No callers were found in the current Java tree; SamLine uses fixMatch instead.
	 * Numeric runs advance cursors directly and do not account for intervening insertions
	 * (#001), so this is not equivalent to fixMatch. A terminal numeric run stays pending.
	 * Shares all cursor/mode state with correction; do not interleave the two algorithms.
	 * @return true after a reported substitution, false when the MD text is exhausted
	 */
	boolean nextSub(){
		sym=0;
		while(mdPos<mdTag.length()){
			char c=mdTag.charAt(mdPos);
			mdPos++;

			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
				mode=NORMAL;
			}else{
				if(current>0){
					//TODO: Probable bug #001 - direct MD-run advances ignore intervening I entries in longmatch.
					//For MD2A0 over mImS, nextSub reports match/query position2 instead of3 (reference2).
					//No iterator callers found; retained as dormant source behavior, not runtime-verified.
					bpos+=current;
					rpos+=current;
					matchPos+=current;
					assert(mode==NORMAL) : mode+", "+current;
					current=0;
				}
				if(c=='^'){mode=DEL;}else if(mode==DEL){
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

		}
		return false;
	}

	/** Returns match cursor minus one, including leading clips and processed events.
	 * @return Relative zero-based predecessor index, or -1 before any advance */
	public int matchPosition(){return matchPos-1;}

	/** Returns query cursor minus one; leading clips count toward this cursor.
	 * @return Relative zero-based predecessor index, or -1 before any advance */
	public int basePosition(){return bpos-1;}

	/** Returns reference cursor minus one, relative to alignment start rather than SamLine.pos.
	 * @return Relative zero-based predecessor offset, or -1 before any reference advance */
	public int refPosition(){return rpos-1;}

	/** Returns the most recently stored MD event character.
	 * @return Stored MD character from the most recent symbolic event
	 * @throws AssertionError With assertions enabled when no nonzero symbol is stored */
	public char symbol(){
		assert(sym!=0);
		return sym;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Next expanded-match entry to consume. */
	private int matchPos;
	/** Next query base offset, including leading clips and insertions. */
	private int bpos;
	/** Next reference offset relative to alignment start. */
	private int rpos;
	/** Last stored event character; reset to zero on entry to either walker method. */
	private char sym;

	/** Retained MD text, including any recognized prefix. */
	private String mdTag;
	/** Optional retained CIGAR; fixMatch scans for N but does not use the result. */
	private String cigar;
	/** Caller-owned expanded match array; fixMatch changes substitution entries in place. */
	private byte[] longmatch;
	/** Index of the next MD character to consume. */
	private int mdPos;
	/** Accumulated numeric run; terminal digits are not applied to cursors. */
	private int current;
	/** Current interpretation of symbolic MD characters. */
	private int mode;

	/** Optional SAM record retained for length-inconsistency diagnostics. */
	private SamLine sl;

	/** Parser modes; INS is used only by the legacy iterator. */
	private static final int NORMAL=0, SUB=1, DEL=2, INS=3;

}
