package stream;

import dna.AminoAcid;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date May 5, 2016
 * 
 *
 */
public class MDWalker {

	/**
	 * Constructs an MDWalker for parsing alignment information.
	 * Initializes position tracking variables and advances past initial clipping.
	 *
	 * @param tag MD tag string containing mismatch information
	 * @param cigar_ CIGAR string for debugging (optional, may be null)
	 * @param longmatch_ Match array to be corrected with MD tag information
	 * @param sl_ SamLine object for debugging purposes
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

	/**
	 * Corrects the match array using MD tag information and read bases.
	 * Processes the entire MD tag, updating the longmatch array to reflect
	 * actual substitutions, deletions, and insertions in the alignment.
	 * @param bases Read sequence bases for ambiguous base checking
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
					//					System.err.println(mpos+", "+current+", "+mpos2);
					assert(mode==NORMAL) : mode+", "+current;
					current=0;
				}

				//Fixes subs after dels getting ignored
				if(mode==DEL && (matchPos<longmatch.length && longmatch[matchPos]!='D')){
					mode=SUB;
				}

				while(matchPos<matchPos2 || (matchPos<longmatch.length && (longmatch[matchPos]=='I' || false))){
					assert(matchPos<longmatch.length) : longmatch.length+"\n"+sl.toString()+"\n"+new String(longmatch);
					if(longmatch[matchPos]=='I'){
						//						System.err.println("I: mpos="+mpos+", bpos="+bpos);
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
						//						System.err.println("M: mpos="+mpos+", bpos="+bpos);
						rpos++;
						bpos++;
						matchPos++;
					}
				}

				//				while(mpos<longmatch.length && longmatch[mpos]=='I'){
				//					System.err.println("I2: mpos="+mpos+", bpos="+bpos);
				//					mpos++;
				//					bpos++;
				//				}

				if(c=='^'){
					mode=DEL;
					//					System.err.println("c="+((char)c)+", mpos="+mpos+", rpos="+rpos+", bpos="+bpos+", mode="+mode+(mode==NORMAL ? "" : ", match="+(char)longmatch[mpos-1])+"\n"+new String(longmatch));
				}else if(mode==DEL){
					//					System.err.println("c="+((char)c)+", mpos="+mpos+", rpos="+rpos+", bpos="+bpos+", mode="+mode+(mode==NORMAL ? "" : ", match="+(char)longmatch[mpos-1])+"\n"+new String(longmatch));
					rpos++;
					matchPos++;
					sym=c;
				}
				//				else if(longmatch[mpos]=='I'){
				//					mode=INS;
				//					bpos++;
				//					mpos++;
				//					sym=c;
				//				}
				else if(mode==NORMAL || mode==SUB){
					// Consume any pending deletions at current position
					while(matchPos<longmatch.length && longmatch[matchPos]=='D'){
						rpos++;
						matchPos++;
					}
					if(matchPos>=longmatch.length){break;}
					longmatch[matchPos]=(byte)'S';
					if((bases!=null && !AminoAcid.isFullyDefined(bases[bpos])) || !AminoAcid.isFullyDefined(c)){longmatch[matchPos]='N';}
					mode=SUB;
					//					System.err.println("c="+((char)c)+", mpos="+mpos+", rpos="+rpos+", bpos="+bpos+", mode="+mode+(mode==NORMAL ? "" : ", match="+(char)longmatch[mpos-1])+"\n"+new String(longmatch));
					bpos++;
					rpos++;
					matchPos++;
					sym=c;
				}else{
					assert(false);
				}

			}

		}
		//		System.err.println();
		//		assert((bases==null || Read.calcMatchLength(longmatch)==bases.length)) :
		//			bases.length+", "+Read.calcMatchLength(longmatch)+"\n"+new String(longmatch)+"\n"
		//					+ new String(Read.toShortMatchString(longmatch))+"\n"+mdTag;
	}

	/**
	 * Advances to the next substitution in the MD tag.
	 * Updates position counters and returns true if a substitution is found.
	 * @return true if another substitution exists, false if end of MD tag reached
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
					bpos+=current;
					rpos+=current;
					matchPos+=current;
					assert(mode==NORMAL) : mode+", "+current;
					current=0;
				}
				if(c=='^'){mode=DEL;}
				else if(mode==DEL){
					rpos++;
					matchPos++;
					sym=c;
				}else if(longmatch[matchPos]=='I'){
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
					//					System.err.println("c="+((char)c)+", mpos="+mpos+", rpos="+rpos+", bpos="+bpos+", mode="+mode+(mode==NORMAL ? "" : ", match="+(char)longmatch[mpos-1])+"\n"+new String(longmatch));
					return true;
				}
			}

			//			System.err.println("c="+((char)c)+", mpos="+mpos+", rpos="+rpos+", bpos="+bpos+", mode="+mode+(mode==NORMAL ? "" : ", match="+(char)longmatch[mpos-1])+"\n"+new String(longmatch));
		}
		return false;
	}

	/** Gets the current position in the match string.
	 * @return Zero-based position in match array (last processed position) */
	public int matchPosition(){
		return matchPos-1;
	}

	/** Gets the current position in read bases.
	 * @return Zero-based position in read sequence (last processed position) */
	public int basePosition(){
		return bpos-1;
	}

	/** Gets the current position in reference sequence.
	 * @return Zero-based position in reference (last processed position) */
	public int refPosition(){
		return rpos-1;
	}

	/** Gets the current substitution symbol from MD tag.
	 * @return Reference base character at current substitution position */
	public char symbol(){
		assert(sym!=0);
		return sym;
	}

	/** Position in match string (excluding clipping and insertions) */
	private int matchPos;
	/** Position in read bases (excluding clipping and insertions) */
	private int bpos;
	/** Position in reference bases (excluding clipping) */
	private int rpos;
	/** Current substitution symbol from MD tag */
	private char sym;

	/** MD tag string containing mismatch and deletion information */
	private String mdTag;
	/** CIGAR string for debugging purposes (optional) */
	private String cigar; //Optional; for debugging
	/** Match array to be corrected with MD tag information */
	private byte[] longmatch;
	/** Current parsing position within the MD tag string */
	private int mdPos;
	/** Current numeric value being accumulated from MD tag digits */
	private int current;
	/** Current parsing mode (NORMAL, SUB, DEL, or INS) */
	private int mode;

	/** SamLine object for debugging purposes */
	private SamLine sl;

	//	private int dels=0, subs=0, normals=0;
	/** Mode constant for normal matching bases */
	private static final int NORMAL=0, SUB=1, DEL=2, INS=3;

}
