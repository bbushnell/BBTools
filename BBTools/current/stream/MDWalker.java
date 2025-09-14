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
	
	MDWalker(String tag, String cigar_, byte[] longmatch_, SamLine sl_){//SamLine is just for debugging
		mdTag=tag;
		cigar=cigar_;
		longmatch=longmatch_;
		sl=sl_;
		mdPos=(mdTag.startsWith("MD:Z:") ? 5 : 0);
		
		matchPos=0;
		bpos=0;
		rpos=0;
		sym=0;
		current=0;
		mode=0;
		
		while(longmatch[matchPos]=='C'){
			matchPos++;
			bpos++;
		}
	}
	
	void fixMatch(byte[] bases){
		final boolean cigarContainsN=(cigar==null || cigar.indexOf('N')>=0);
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
					}else if(longmatch[matchPos]=='D'){//This should only happen if the cigar string contains an 'N'
						assert(cigarContainsN) : matchPos+"\n"+new String(longmatch)+"\n"+mdTag+"\n"+cigar;
//						System.err.println("M: mpos="+mpos+", bpos="+bpos);
						while(longmatch[matchPos]=='D' && matchPos<longmatch.length){
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
					if(longmatch[matchPos]=='D'){
						assert(cigarContainsN) : matchPos+"\n"+new String(longmatch)+"\n"+mdTag+"\n"+cigar;
						while(longmatch[matchPos]=='D' && matchPos<longmatch.length){
							rpos++;
							matchPos++;
						}
					}
					assert(longmatch[matchPos]!='D') : matchPos+"\n"+new String(longmatch)+"\n"+mdTag+"\n"+cigar;
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
	
	public int matchPosition(){
		return matchPos-1;
	}
	
	public int basePosition(){
		return bpos-1;
	}
	
	public int refPosition(){
		return rpos-1;
	}
	
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
	private char sym;
	
	private String mdTag;
	private String cigar; //Optional; for debugging
	private byte[] longmatch;
	private int mdPos;
	private int current;
	private int mode;
	
	private SamLine sl;
	
//	private int dels=0, subs=0, normals=0;
	private static final int NORMAL=0, SUB=1, DEL=2, INS=3;
	
}
