package stream;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import align2.GapTools;
import align2.QualityTools;
import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import shared.TrimRead;
import structures.ByteBuilder;
import structures.FloatList;
import ukmer.Kmer;

public final class Read implements Comparable<Read>, Cloneable, Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -1026645233407290096L;

	public static void main(String[] args){
		byte[] a=args[0].getBytes();
		System.out.println(new String(a));
		byte[] b=toShortMatchString(a);
		System.out.println(new String(b));
		byte[] c=toLongMatchString(b);
		System.out.println(new String(c));
		byte[] d=toLongMatchString(c);
		System.out.println(new String(d));
//		byte[] e=toShortMatchString(b);
//		System.out.println(new String(e));
		
	}
	
	public Read(byte[] bases_, byte[] quals_, long id_){
		this(bases_, quals_, Long.toString(id_), id_);
	}
	
	public Read(byte[] bases_, byte[] quals_, String name_, long id_){
		this(bases_, quals_, name_, id_, 0, -1, -1, -1);
	}
	
	public Read(byte[] bases_, byte[] quals_, String name_, long id_, int flag_){
		this(bases_, quals_, name_, id_, flag_, -1, -1, -1);
	}
	
	public Read(byte[] s_, byte[] quals_, long id_, int chrom_, int start_, int stop_, byte strand_){
		this(s_, quals_, Long.toString(id_), id_, (int)strand_, chrom_, start_, stop_);
	}
	
//	public Read(byte[] bases_, byte[] quals_, String id_, long numericID_, byte strand_, int chrom_, int start_, int stop_){
//		this(bases_, quals_, id_, numericID_, (int)strand_, chrom_, start_, stop_);
//		assert(strand_==0 || strand_==1);
//		assert(start_<=stop_) : chrom_+", "+start_+", "+stop_+", "+numericID_;
//	}
	
	/** Note that strand can be used as flag */
	public Read(byte[] bases_, byte[] quals_, String id_, long numericID_, int flags_, int chrom_, int start_, int stop_){
		flags=flags_&~VALIDATEDMASK;
		bases=bases_;
		quality=quals_;

		chrom=chrom_;
		start=start_;
		stop=stop_;
		
		id=id_;
		numericID=numericID_;

//		assert(amino()) : Shared.AMINO_IN;//123
		if(VALIDATE_IN_CONSTRUCTOR){validate(true);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Validation          ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean validate(final boolean processAssertions){
		assert(!validated());
		
//		if(false){//This causes problems with error-corrected PacBio reads.
//			boolean x=(quality==null || quality.length<1 || quality[0]<=80 || !FASTQ.DETECT_QUALITY || FASTQ.IGNORE_BAD_QUALITY);
//			if(!x){
//				if(processAssertions){
//					KillSwitch.kill("Quality value ("+quality[0]+") appears too high.\n"+Arrays.toString(quality)+
//							"\n"+Arrays.toString(bases)+"\n"+numericID+"\n"+id+"\n"+FASTQ.ASCII_OFFSET);
//				}
//				return false;
//			}
//		}
		
		if(bases==null){
			quality=null; //I could require this to be true
			if(FIX_HEADER){fixHeader(processAssertions);}
			setValidated(true);
			return true;
		}
		
		validateQualityLength(processAssertions);
		
		final boolean passesJunk;
		
//		assert(false) : SKIP_SLOW_VALIDATION+","+VALIDATE_BRANCHLESS+","+JUNK_MODE;
		if(SKIP_SLOW_VALIDATION){
			passesJunk=true;
		}else if(!aminoacid()){
			if(U_TO_T){uToT();}
			if(VALIDATE_BRANCHLESS){
				passesJunk=validateCommonCase_branchless(processAssertions);
			}else{
				passesJunk=validateCommonCase(processAssertions);
			}
		}else{
//			if(U_TO_T){uToT();}  This is amino, so...
			fixCase();
			passesJunk=validateJunk(processAssertions);
			if(CHANGE_QUALITY){fixQuality();}
		}
		
		if(FIX_HEADER){fixHeader(processAssertions);}
		
		setValidated(true);
		
		return true;
	}
	

	
	public boolean checkQuality(){
		if(quality==null){return true;}
		for(int i=0; quality!=null && i<bases.length; i++){
			if((bases[i]=='N')!=(quality[i]<1)){
				assert(false) : (char)bases[i]+", "+quality[i]+", "+i+", "+CHANGE_QUALITY+", "+MIN_CALLED_QUALITY+"\n"+toFastq();
				return false;
			}
		}
		return true;
	}
	
	private void uToT(){
		if(bases==null){return;}
		for(int i=0; i<bases.length; i++){
			bases[i]=AminoAcid.uToT[bases[i]];
		}
	}
	
	private void tToU(){
		if(bases==null){return;}
		for(int i=0; i<bases.length; i++){
			bases[i]=AminoAcid.tToU[bases[i]];
		}
	}
	
	private boolean validateJunk(boolean processAssertions){
		assert(bases!=null);
		if(JUNK_MODE==IGNORE_JUNK){return true;}
		final byte nocall;
		final byte[] toNum;
		final boolean aa=aminoacid();
		if(aa){
			nocall='X';
			toNum=AminoAcid.acidToNumberExtended;
		}else{
			nocall='N';
			toNum=AminoAcid.baseToNumberExtended;
		}
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int num=toNum[b];
//			System.err.println(Character.toString(b)+" -> "+num);
			if(num<0){
				if(JUNK_MODE==FIX_JUNK){
					bases[i]=nocall;
				}else if(JUNK_MODE==FLAG_JUNK){
					setJunk(true);
					return false;
				}else{
					if(processAssertions){
						KillSwitch.kill("\nAn input file appears to be misformatted:\n"
							+ "The character with ASCII code "+b+" appeared where "+(aa ? "an amino acid" : "a base")+" was expected"
									+ (b>31 && b<127 ? ": '"+Character.toString((char)b)+"'\n" : ".\n")
							+ "Sequence #"+numericID+"\n"
							+ "Sequence ID '"+id+"'\n"
							+ "Sequence: "+Tools.toStringSafe(bases)+"\n"
							+ "Flags: "+Long.toBinaryString(flags)+"\n\n"
							+ "This can be bypassed with the flag 'tossjunk', 'fixjunk', or 'ignorejunk'");
					}
					setJunk(true);
					return false;
				}
			}
		}
		return true;
	}
	
	private void validateQualityLength(boolean processAssertions){
		if(quality==null || quality.length==bases.length){return;}
		if(TOSS_BROKEN_QUALITY){
			quality=null;
			setDiscarded(true);
			setJunk(true);
		}else if(NULLIFY_BROKEN_QUALITY){
			quality=null;
			setJunk(true);
		}else if(FLAG_BROKEN_QUALITY){
			setJunk(true);
		}else{
			boolean x=false;
			assert(x=processAssertions);
			if(x){
				KillSwitch.kill("\nMismatch between length of bases and qualities for read "+numericID+" (id="+id+").\n"+
						"# qualities="+quality.length+", # bases="+bases.length+"\n\n"+
						FASTQ.qualToString(quality)+"\n"+new String(bases)+"\n\n"
						+ "This can be bypassed with the flag 'tossbrokenreads' or 'nullifybrokenquality'");
			}
		}
	}
	
	private void fixQuality(){
		if(quality==null || !CHANGE_QUALITY){return;}

		final byte[] toNumber=aminoacid() ? AminoAcid.acidToNumber : AminoAcid.baseToNumber;
		
		for(int i=0; i<quality.length; i++){
			byte b=bases[i];
			byte q=quality[i];
			quality[i]=capQuality(q, b);
//			if(toNumber[b]>=0){
//				if(q<MIN_CALLED_QUALITY){
//					quality[i]=MIN_CALLED_QUALITY;
//				}else if(q>MAX_CALLED_QUALITY){
//					quality[i]=MAX_CALLED_QUALITY;
//				}
//			}else{
//				quality[i]=0;
//			}
		}
	}
	
	private void fixCase(){
		if(bases==null || (!DOT_DASH_X_TO_N && !TO_UPPER_CASE && !LOWER_CASE_TO_N)){return;}
		final boolean aa=aminoacid();
		
		final byte[] caseMap, ddxMap;
		if(!aa){
			caseMap=TO_UPPER_CASE ? AminoAcid.toUpperCase : LOWER_CASE_TO_N ? AminoAcid.lowerCaseToNocall : null;
			ddxMap=DOT_DASH_X_TO_N ? AminoAcid.dotDashXToNocall : null;
		}else{
			caseMap=TO_UPPER_CASE ? AminoAcid.toUpperCase : LOWER_CASE_TO_N ? AminoAcid.lowerCaseToNocallAA : null;
			ddxMap=DOT_DASH_X_TO_N ? AminoAcid.dotDashXToNocallAA : null;
		}
		
//		assert(false) : (AminoAcid.toUpperCase==caseMap)+", "+ddxMap;
		
		if(DOT_DASH_X_TO_N){
			if(TO_UPPER_CASE || LOWER_CASE_TO_N){
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					bases[i]=caseMap[ddxMap[b]];
				}
			}else{
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					bases[i]=ddxMap[b];
				}
			}
		}else{
			if(TO_UPPER_CASE || LOWER_CASE_TO_N){
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					bases[i]=caseMap[b];
				}
			}else{
				assert(false);
			}
		}
	}
	
	private boolean validateCommonCase_branchless(boolean processAssertions){
		
		assert(!aminoacid());
		assert(bases!=null);
		
		if(TO_UPPER_CASE || LOWER_CASE_TO_N){fixCase();}
		
		final byte nocall='N';
		final byte[] toNum=AminoAcid.baseToNumber;
		final byte[] map=(DOT_DASH_X_TO_N && IUPAC_TO_N ? AminoAcid.baseToACGTN : 
			DOT_DASH_X_TO_N ? AminoAcid.dotDashXToNocall : IUPAC_TO_N ? AminoAcid.iupacToNocall : null);
		
		if(JUNK_MODE==IGNORE_JUNK){
			if(quality!=null && CHANGE_QUALITY){
				if(DOT_DASH_X_TO_N || IUPAC_TO_N){
					for(int i=0; i<bases.length; i++){
						final byte b=map[bases[i]];
						final byte q=quality[i];
						final int num=toNum[b];
						bases[i]=b;
						quality[i]=(num>=0 ? qMap[q] : 0);
					}
				}else{
					for(int i=0; i<bases.length; i++){
						final byte b=bases[i];
						final byte q=quality[i];
						final int num=toNum[b];
						quality[i]=(num>=0 ? qMap[q] : 0);
					}
				}
			}else if(DOT_DASH_X_TO_N || IUPAC_TO_N){
				for(int i=0; i<bases.length; i++){
					byte b=map[bases[i]];
					bases[i]=b;
				}
			}
			return true;
		}
		
		int junkOr=0;
//		int iupacOr=0;
		final byte[] toNumE=AminoAcid.baseToNumberExtended;
		
//		asdf
		
		if(DOT_DASH_X_TO_N || IUPAC_TO_N){
			if(quality!=null && CHANGE_QUALITY){
				for(int i=0; i<bases.length; i++){
					final byte b=map[bases[i]];
					final byte q=quality[i];
					final int numE=toNumE[b];
					final int num=toNum[b];
					junkOr|=numE;
//					iupacOr|=num;
					bases[i]=b;
					quality[i]=(num>=0 ? qMap[q] : 0);
				}
			}else{
				for(int i=0; i<bases.length; i++){
					final byte b=map[bases[i]];
					final int numE=toNumE[b];
//					final int num=toNum[b];
					junkOr|=numE;
//					iupacOr|=num;
					bases[i]=b;
				}
			}
		}else{
			if(quality!=null && CHANGE_QUALITY){
				for(int i=0; i<bases.length; i++){
					final byte b=bases[i];
					final byte q=quality[i];
					final int numE=toNumE[b];
					final int num=toNum[b];
					junkOr|=numE;
//					iupacOr|=num;
					quality[i]=(num>=0 ? qMap[q] : 0);
				}
			}else{
				for(int i=0; i<bases.length; i++){
					final byte b=bases[i];
					final int numE=toNumE[b];
//					final int num=toNum[b];
					junkOr|=numE;
//					iupacOr|=num;
				}
			}
		}
		
//		System.err.println(junkOr+", "+JUNK_MODE);
		
		//Common case
		if(junkOr>=0){return true;}
//		if(junkOr>=0 && (JUNK_MODE!=FIX_JUNK_AND_IUPAC || iupacOr>=0)){return true;}
//		
//		assert(junkOr<0 || (JUNK_MODE==FIX_JUNK_AND_IUPAC && iupacOr<0));
		
		//TODO: I could disable VALIDATE_BRANCHLESS here, if it's not final
		//VALIDATE_BRANCHLESS=false;
		if(JUNK_MODE==FIX_JUNK){
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				final int numE=toNumE[b];

				if(numE<0){
					bases[i]=nocall;
					if(quality!=null){quality[i]=0;}
				}
			}
			return true;
		}else if(JUNK_MODE==FLAG_JUNK){
			setJunk(true);
			return false;
		}else if(JUNK_MODE==FIX_JUNK_AND_IUPAC){
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				final byte c=AminoAcid.baseToACGTN[b];
				bases[i]=c;
			}
			return true;
		}else{
			if(processAssertions){
				int i=0;
				for(; i<bases.length; i++){
					if(toNumE[bases[i]]<0){break;}
				}
				byte b=bases[i];
				KillSwitch.kill("\nAn input file appears to be misformatted:\n"
						+ "The character with ASCII code "+b+" appeared where a base was expected"
							+ (b>31 && b<127 ? ": '"+Character.toString((char)b)+"'\n" : ".\n")
						+ "Sequence #"+numericID+"\n"
						+ "Sequence ID: '"+id+"'\n"
						+ "Sequence: '"+Tools.toStringSafe(bases)+"'\n\n"
						+ "This can be bypassed with the flag 'tossjunk', 'fixjunk', or 'ignorejunk'");
			}
			setJunk(true);
			return false;
		}
	}
	
	private boolean validateCommonCase(boolean processAssertions){
		
		assert(!aminoacid());
		assert(bases!=null);
		
		if(TO_UPPER_CASE || LOWER_CASE_TO_N){fixCase();}
		
		final byte nocall='N';
		final byte[] toNum=AminoAcid.baseToNumber;
		final byte[] ddxMap=AminoAcid.dotDashXToNocall;
		
		if(JUNK_MODE==IGNORE_JUNK){
			if(quality!=null && CHANGE_QUALITY){
				if(DOT_DASH_X_TO_N){
					for(int i=0; i<bases.length; i++){
						final byte b=ddxMap[bases[i]];
						final byte q=quality[i];
						final int num=toNum[b];
						bases[i]=b;
						quality[i]=(num>=0 ? qMap[q] : 0);
					}
				}else{
					for(int i=0; i<bases.length; i++){
						final byte b=bases[i];
						final byte q=quality[i];
						final int num=toNum[b];
						quality[i]=(num>=0 ? qMap[q] : 0);
					}
				}
			}else if(DOT_DASH_X_TO_N){
				for(int i=0; i<bases.length; i++){
					byte b=ddxMap[bases[i]];
					bases[i]=b;
				}
			}
		}else if(DOT_DASH_X_TO_N){
			final byte[] toNumE=AminoAcid.baseToNumberExtended;
			if(quality!=null && CHANGE_QUALITY){
				for(int i=0; i<bases.length; i++){
					byte b=ddxMap[bases[i]];
					final byte q=quality[i];
					final int numE=toNumE[b];

					if(numE<0){
						if(JUNK_MODE==FIX_JUNK){
							b=nocall;
						}else if(JUNK_MODE==FLAG_JUNK){
							setJunk(true);
							return false;
						}else{
							if(processAssertions){
								KillSwitch.kill("\nAn input file appears to be misformatted:\n"
										+ "The character with ASCII code "+bases[1]+" appeared where a base was expected.\n"
										+ "Sequence #"+numericID+"\n"
										+ "Sequence ID '"+id+"'\n"
										+ "Sequence: "+Tools.toStringSafe(bases)+"\n\n"
										+ "This can be bypassed with the flag 'tossjunk', 'fixjunk', or 'ignorejunk'");
							}
							setJunk(true);
							return false;
						}
					}

					final int num=toNum[b];
					bases[i]=b;
					quality[i]=(num>=0 ? qMap[q] : 0);
				}
			}else{
				for(int i=0; i<bases.length; i++){
					byte b=ddxMap[bases[i]];
					final int numE=toNumE[b];

					if(numE<0){
						if(JUNK_MODE==FIX_JUNK){
							b=nocall;
						}else if(JUNK_MODE==FLAG_JUNK){
							setJunk(true);
							return false;
						}else{
							if(processAssertions){
								KillSwitch.kill("\nAn input file appears to be misformatted:\n"
										+ "The character with ASCII code "+bases[1]+" appeared where a base was expected.\n"
										+ "Sequence #"+numericID+"\n"
										+ "Sequence ID '"+id+"'\n"
										+ "Sequence: "+Tools.toStringSafe(bases)+"\n\n"
										+ "This can be bypassed with the flag 'tossjunk', 'fixjunk', or 'ignorejunk'");
							}
							setJunk(true);
							return false;
						}
					}

					bases[i]=b;
				}
			}
		}else{
			final byte[] toNumE=AminoAcid.baseToNumberExtended;
			if(quality!=null && CHANGE_QUALITY){
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					final byte q=quality[i];
					final int numE=toNumE[b];

					if(numE<0){
						if(JUNK_MODE==FIX_JUNK){
							bases[i]=b=nocall;
						}else if(JUNK_MODE==FLAG_JUNK){
							setJunk(true);
							return false;
						}else{
							if(processAssertions){
								KillSwitch.kill("\nAn input file appears to be misformatted:\n"
										+ "The character with ASCII code "+bases[1]+" appeared where a base was expected.\n"
										+ "Sequence #"+numericID+"\n"
										+ "Sequence ID '"+id+"'\n"
										+ "Sequence: "+Tools.toStringSafe(bases)+"\n\n"
										+ "This can be bypassed with the flag 'tossjunk', 'fixjunk', or 'ignorejunk'");
							}
							setJunk(true);
							return false;
						}
					}

					final int num=toNum[b];
					bases[i]=b;
					quality[i]=(num>=0 ? qMap[q] : 0);
				}
			}else{
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					final int numE=toNumE[b];

					if(numE<0){
						if(JUNK_MODE==FIX_JUNK){
							bases[i]=b=nocall;
						}else if(JUNK_MODE==FLAG_JUNK){
							setJunk(true);
							return false;
						}else{
							if(processAssertions){
								KillSwitch.kill("\nAn input file appears to be misformatted:\n"
										+ "The character with ASCII code "+bases[1]+" appeared where a base was expected.\n"
										+ "Sequence #"+numericID+"\n"
										+ "Sequence ID '"+id+"'\n"
										+ "Sequence: "+Tools.toStringSafe(bases)+"\n\n"
										+ "This can be bypassed with the flag 'tossjunk', 'fixjunk', or 'ignorejunk'");
							}
							setJunk(true);
							return false;
						}
					}
				}
			}
		}
		return true;
	}
	
	private final void fixHeader(boolean processAssertions){
		id=Tools.fixHeader(id, ALLOW_NULL_HEADER, processAssertions);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Various            ----------------*/
	/*--------------------------------------------------------------*/
	
	
	private static final int absdif(int a, int b){
		return a>b ? a-b : b-a;
	}
	
	/** Returns true if these reads are identical, allowing at most n no-calls and m mismatches of max quality q*/
	public boolean isDuplicateByBases(Read r, int nmax, int mmax, byte qmax, boolean banSameQualityMismatch){
		return isDuplicateByBases(r, nmax, mmax, qmax, banSameQualityMismatch, false);
	}
	
	
	
	/** Returns true if these reads are identical, allowing at most n no-calls and m mismatches of max quality q*/
	public boolean isDuplicateByBases(Read r, int nmax, int mmax, byte qmax, boolean banSameQualityMismatch, boolean allowDifferentLength){
		int n=0, m=0;
		assert(r.length()==bases.length) : "Merging different-length reads is supported but seems to be not useful.";
		if(!allowDifferentLength && r.length()!=bases.length){return false;}
		int minLen=Tools.min(bases.length, r.length());
		for(int i=0; i<minLen; i++){
			byte b1=bases[i];
			byte b2=r.bases[i];
			if(b1=='N' || b2=='N'){
				n++;
				if(n>nmax){return false;}
			}else if(b1!=b2){
				m++;
				if(m>mmax){return false;}
				if(quality[i]>qmax && r.quality[i]>qmax){return false;}
				if(banSameQualityMismatch && quality[i]==r.quality[i]){return false;}
			}
		}
		return true;
	}
	
	public boolean isDuplicateByMapping(Read r, boolean bothEnds, boolean checkAlignment){
		if(bases.length!=r.length()){
			return isDuplicateByMappingDifferentLength(r, bothEnds, checkAlignment);
		}
		assert(this!=r && mate!=r);
		assert(!bothEnds || bases.length==r.length());
		if(!mapped() || !r.mapped()){return false;}
//		if(chrom==-1 && start==-1){return false;}
		if(chrom<1 && start<1){return false;}
		
//		if(chrom!=r.chrom || strand()!=r.strand() || start!=r.start){return false;}
////		if(mate==null && stop!=r.stop){return false;} //For unpaired reads, require both ends match
//		if(stop!=r.stop){return false;} //For unpaired reads, require both ends match
//		return true;
		
		if(chrom!=r.chrom || strand()!=r.strand()){return false;}
		if(bothEnds){
			if(start!=r.start || stop!=r.stop){return false;}
		}else{
			if(strand()==Shared.PLUS){
				if(start!=r.start){return false;}
			}else{
				if(stop!=r.stop){return false;}
			}
		}
		if(checkAlignment){
			if(perfect() && r.perfect()){return true;}
			if(match!=null && r.match!=null){
				if(match.length!=r.match.length){return false;}
				for(int i=0; i<match.length; i++){
					byte a=match[i];
					byte b=r.match[i];
					if(a!=b){
						if((a=='D') != (b=='D')){return false;}
						if((a=='I' || a=='X' || a=='Y') != (b=='I' || b=='X' || b=='Y')){return false;}
					}
				}
			}
		}
		return true;
	}
	
	public boolean isDuplicateByMappingDifferentLength(Read r, boolean bothEnds, boolean checkAlignment){
		assert(this!=r && mate!=r);
		assert(bases.length!=r.length());
		if(bothEnds){return false;}
//		assert(!bothEnds || bases.length==r.length());
		if(!mapped() || !r.mapped()){return false;}
//		if(chrom==-1 && start==-1){return false;}
		if(chrom<1 && start<1){return false;}
		
//		if(chrom!=r.chrom || strand()!=r.strand() || start!=r.start){return false;}
////		if(mate==null && stop!=r.stop){return false;} //For unpaired reads, require both ends match
//		if(stop!=r.stop){return false;} //For unpaired reads, require both ends match
//		return true;
		
		if(chrom!=r.chrom || strand()!=r.strand()){return false;}

		if(strand()==Shared.PLUS){
			if(start!=r.start){return false;}
		}else{
			if(stop!=r.stop){return false;}
		}
		
		if(checkAlignment){
			if(perfect() && r.perfect()){return true;}
			if(match!=null && r.match!=null){
				int minLen=Tools.min(match.length, r.match.length);
				for(int i=0; i<minLen; i++){
					byte a=match[i];
					byte b=r.match[i];
					if(a!=b){
						if((a=='D') != (b=='D')){return false;}
						if((a=='I' || a=='X' || a=='Y') != (b=='I' || b=='X' || b=='Y')){return false;}
					}
				}
			}
		}
		return true;
	}
	
	public void merge(Read r, boolean mergeVectors, boolean mergeN){mergePrivate(r, mergeVectors, mergeN, true);}
	
	private void mergePrivate(Read r, boolean mergeVectors, boolean mergeN, boolean mergeMate){
		assert(r!=this);
		assert(r!=this.mate);
		assert(r!=r.mate);
		assert(this!=this.mate);
		assert(r.mate==null || r.mate.mate==r);
		assert(this.mate==null || this.mate.mate==this);
		assert(r.mate==null || r.numericID==r.mate.numericID);
		assert(mate==null || numericID==mate.numericID);
		mergeN=(mergeN||mergeVectors);
		
		assert(r.length()==bases.length) : "Merging different-length reads is supported but seems to be not useful.";
		
		if((mergeN || mergeVectors) && bases.length<r.length()){
			int oldLenB=bases.length;
			start=Tools.min(start, r.start);
			stop=Tools.max(stop, r.stop);
			mapScore=Tools.max(mapScore, r.mapScore);
			
			bases=KillSwitch.copyOfRange(bases, 0, r.length());
			quality=KillSwitch.copyOfRange(quality, 0, r.quality.length);
			for(int i=oldLenB; i<bases.length; i++){
				bases[i]='N';
				quality[i]=0;
			}
			match=null;
			r.match=null;
		}
		
		copies+=r.copies;
		
		
//		if(numericID==11063941 || r.numericID==11063941 || numericID==8715632){
//			System.err.println("***************");
//			System.err.println(this.toText()+"\n");
//			System.err.println(r.toText()+"\n");
//			System.err.println(mergeVectors+", "+mergeN+", "+mergeMate+"\n");
//		}
		
		boolean pflag1=perfect();
		boolean pflag2=r.perfect();

		final int minLenB=Tools.min(bases.length, r.length());
		
		if(mergeN){
			if(quality==null){
				for(int i=0; i<minLenB; i++){
					byte b=r.bases[i];
					if(bases[i]=='N' && b!='N'){bases[i]=b;}
				}
			}else{
				for(int i=0; i<minLenB; i++){
					final byte b1=bases[i];
					final byte b2=r.bases[i];
					final byte q1=Tools.max((byte)0, quality[i]);
					final byte q2=Tools.max((byte)0, r.quality[i]);
					if(b1==b2){
						if(b1=='N'){
							//do nothing
						}else if(mergeVectors){
							//merge qualities
							//						quality[i]=(byte) Tools.min(40, q1+q2);
							if(q1>=q2){
								quality[i]=(byte) Tools.min(48, q1+1+q2/4);
							}else{
								quality[i]=(byte) Tools.min(48, q2+1+q1/4);
							}
						}
					}else if(b1=='N'){
						bases[i]=b2;
						quality[i]=q2;
					}else if(b2=='N'){
						//do nothing
					}else if(mergeVectors){
						if(q1<1 && q2<1){
							//Special case - e.g. Illumina calls bases at 0 quality.
							//Possibly best to keep the matching allele if one matches the ref.
							//But for now, do nothing.
							//This was causing problems changing perfect match strings into imperfect matches.
						}else if(q1==q2){
							assert(b1!=b2);
							bases[i]='N';
							quality[i]=0;
						}else if(q1>q2){
							bases[i]=b1;
							quality[i]=(byte)(q1-q2/2);
						}else{
							bases[i]=b2;
							quality[i]=(byte)(q2-q1/2);
						}
						assert(quality[i]>=0 && quality[i]<=48);
					}
				}
			}
		}
		
		//TODO:
		//Note that the read may need to be realigned after merging, so the match string may be rendered incorrect.
		
		if(mergeN && match!=null){
			if(r.match==null){match=null;}
			else{
				if(match.length!=r.match.length){match=null;}
				else{
					boolean ok=true;
					for(int i=0; i<match.length && ok; i++){
						byte a=match[i], b=r.match[i];
						if(a!=b){
							if((a=='m' || a=='S') && b=='N'){
								//do nothing;
							}else if(a=='N' && (b=='m' || b=='S')){
								match[i]=b;
							}else{
								ok=false;
							}
						}
					}
					if(!ok){match=null;}
				}
			}
		}
		
		if(mergeMate && mate!=null){
			mate.mergePrivate(r.mate, mergeVectors, mergeN, false);
			assert(copies==mate.copies);
		}
		assert(copies>1);
		
		assert(r!=this);
		assert(r!=this.mate);
		assert(r!=r.mate);
		assert(this!=this.mate);
		assert(r.mate==null || r.mate.mate==r);
		assert(this.mate==null || this.mate.mate==this);
		assert(r.mate==null || r.numericID==r.mate.numericID);
		assert(mate==null || numericID==mate.numericID);
	}
	
	@Override
	public String toString(){return toText(false).toString();}
	
	public ByteBuilder toSites(){return toSites((ByteBuilder)null);}
	
	public ByteBuilder toSites(ByteBuilder sb){
		if(numSites()==0){
			if(sb==null){sb=new ByteBuilder(2);}
			sb.append('.');
		}else{
			if(sb==null){sb=new ByteBuilder(sites.size()*20);}
			int appended=0;
			for(SiteScore ss : sites){
				if(appended>0){sb.append('\t');}
				if(ss!=null){
					ss.toBytes(sb);
					appended++;
				}
			}
			if(appended==0){sb.append('.');}
		}
		return sb;
	}
	
	public ByteBuilder toInfo(){
		if(obj==null){return new ByteBuilder();}
		if(obj.getClass()==ByteBuilder.class){return (ByteBuilder)obj;}
		return new ByteBuilder(obj.toString());
	}
	
	public ByteBuilder toInfo(ByteBuilder bb){
		if(obj==null){return bb;}
		if(obj.getClass()==ByteBuilder.class){return bb.append((ByteBuilder)obj);}
		return bb.append(obj.toString());
	}
	
	public ByteBuilder toFastq(){
		return FASTQ.toFASTQ(this, (ByteBuilder)null);
	}
	
	public ByteBuilder toFastq(ByteBuilder bb){
		return FASTQ.toFASTQ(this, bb);
	}
	
	public ByteBuilder toFasta(){return toFasta(Shared.FASTA_WRAP);}
	public ByteBuilder toFasta(ByteBuilder bb){return toFasta(Shared.FASTA_WRAP, bb);}
	
	public ByteBuilder toFasta(int wrap){
		return toFasta(wrap, (ByteBuilder)null);
	}
	
	public ByteBuilder toFasta(int wrap, ByteBuilder bb){
		if(wrap<1){wrap=Integer.MAX_VALUE;}
		int len=(id==null ? Tools.stringLength(numericID) : id.length())+(bases==null ? 0 : bases.length+bases.length/wrap)+5;
		if(bb==null){bb=new ByteBuilder(len+1);}
		bb.append('>');
		if(id==null){bb.append(numericID);}
		else{bb.append(id);}
		if(bases!=null){
			int pos=0;
			while(pos<bases.length-wrap){
				bb.append('\n');
				bb.append(bases, pos, wrap);
				pos+=wrap;
			}
			if(pos<bases.length){
				bb.append('\n');
				bb.append(bases, pos, bases.length-pos);
			}
		}
		return bb;
	}
	
	public ByteBuilder toSam(){
		return toSam((ByteBuilder)null);
	}
	
	public ByteBuilder toSam(ByteBuilder bb){
		SamLine sl=new SamLine(this, pairnum());
//		System.err.println("Called toSam on read "+id+"; num="+numericID+", pairnum="+pairnum()+"; result="+sl.toString());
		return sl.toBytes(bb);
	}
	
	public static CharSequence header(){

		StringBuilder sb=new StringBuilder();
		sb.append("id");
		sb.append('\t');
		sb.append("numericID");
		sb.append('\t');
		sb.append("chrom");
		sb.append('\t');
		sb.append("strand");
		sb.append('\t');
		sb.append("start");
		sb.append('\t');
		sb.append("stop");
		sb.append('\t');

		sb.append("flags");
		sb.append('\t');
		
		sb.append("copies");
		sb.append('\t');
		
		sb.append("errors,fixed");
		sb.append('\t');
		sb.append("mapScore");
		sb.append('\t');
		sb.append("length");
		sb.append('\t');
		
		sb.append("bases");
		sb.append('\t');
		sb.append("quality");
		sb.append('\t');
		
		sb.append("insert");
		sb.append('\t');
		{
			//These are not really necessary...
			sb.append("avgQual");
			sb.append('\t');
		}
		
		sb.append("match");
		sb.append('\t');
		sb.append("SiteScores: "+SiteScore.header());
		return sb;
	}
	
	public ByteBuilder toText(boolean okToCompressMatch){
		return toText(okToCompressMatch, (ByteBuilder)null);
	}
	
	public ByteBuilder toText(boolean okToCompressMatch, ByteBuilder bb){
		
		final byte[] oldmatch=match;
		final boolean oldshortmatch=this.shortmatch();
		if(COMPRESS_MATCH_BEFORE_WRITING && !shortmatch() && okToCompressMatch){
			match=toShortMatchString(match);
			setShortMatch(true);
		}
		
		if(bb==null){bb=new ByteBuilder();}
		bb.append(id);
		bb.tab();
		bb.append(numericID);
		bb.tab();
		bb.append(chrom);
		bb.tab();
		bb.append(Shared.strandCodes2[strand()]);
		bb.tab();
		bb.append(start);
		bb.tab();
		bb.append(stop);
		bb.tab();
		
		for(int i=maskArray.length-1; i>=0; i--){
			bb.append(flagToNumber(maskArray[i]));
		}
		bb.tab();
		
		bb.append(copies);
		bb.tab();

		bb.append(errors);
		bb.tab();
		bb.append(mapScore);
		bb.tab();
		
		if(bases==null){bb.append('.');}
		else{bb.append(bases);}
		bb.tab();
		
//		int qualSum=0;
//		int qualMin=99999;
		
		if(quality==null){
			bb.append('.');
		}else{
			bb.ensureExtra(quality.length);
			for(int i=0, j=bb.length; i<quality.length; i++, j++){
				byte q=quality[i];
				bb.array[j]=(byte)(q+ASCII_OFFSET);
//				qualSum+=q;
//				qualMin=Tools.min(q, qualMin);
			}
			bb.length+=quality.length;
		}
		bb.tab();
		
		if(insert<1){bb.append('.');}else{bb.append(insert);};
		bb.tab();
		
		if(true || quality==null){
			bb.append('.');
			bb.tab();
		}else{
//			//These are not really necessary...
//			sb.append(qualSum/quality.length);
//			sb.append('\t');
		}
		
		if(match==null){bb.append('.');}
		else{bb.append(match);}
		bb.tab();
		
		if(gaps==null){
			bb.append('.');
		}else{
			for(int i=0; i<gaps.length; i++){
				if(i>0){bb.append('~');}
				bb.append(gaps[i]);
			}
		}
		
		if(sites!=null && sites.size()>0){
			
			assert(absdif(start, stop)<3000 || (gaps==null) == (sites.get(0).gaps==null)) :
				"\n"+this.numericID+"\n"+Arrays.toString(gaps)+"\n"+sites.toString()+"\n";
			
			for(SiteScore ss : sites){
				bb.tab();
				if(ss==null){
					bb.append((byte[])null);
				}else{
					ss.toBytes(bb);
				}
				bb.append(ss==null ? "null" : ss.toText());
			}
		}
		
		if(originalSite!=null){
			bb.tab();
			bb.append('*');
			originalSite.toBytes(bb);
		}
		
		match=oldmatch;
		setShortMatch(oldshortmatch);
		
		return bb;
	}
	
	public static Read fromText(String line){
		if(line.length()==1 && line.charAt(0)=='.'){return null;}
		
		String[] split=line.split("\t");
		
//		if(split.length<17){
//			throw new RuntimeException("Error parsing read from text.\n\n" +
//					"This may be caused be attempting to parse the wrong format.\n" +
//					"Please ensure that the file extension is correct:\n" +
//					"\tFASTQ should end in .fastq or .fq\n" +
//					"\tFASTA should end in .fasta or .fa, .fas, .fna, .ffn, .frn, .seq, .fsa\n" +
//					"\tSAM should end in .sam\n" +
//					"\tNative format should end in .txt or .bread\n" +
//					"If a file is compressed, there must be a compression extension after the format extension:\n" +
//					"\tgzipped files should end in .gz or .gzip\n" +
//					"\tzipped files should end in .zip and have only 1 file per archive\n" +
//					"\tbz2 files should end in .bz2\n");
//		}
		
		final String id=new String(split[0]);
		long numericID=Long.parseLong(split[1]);
		int chrom=Byte.parseByte(split[2]);
//		byte strand=Byte.parseByte(split[3]);
		int start=Integer.parseInt(split[4]);
		int stop=Integer.parseInt(split[5]);
		
		int flags=Integer.parseInt(split[6], 2);
		
		int copies=Integer.parseInt(split[7]);

		int errors;
		int errorsCorrected;
		if(split[8].indexOf(',')>=0){
			String[] estring=split[8].split(",");
			errors=Integer.parseInt(estring[0]);
			errorsCorrected=Integer.parseInt(estring[1]);
		}else{
			errors=Integer.parseInt(split[8]);
			errorsCorrected=0;
		}
		
		int mapScore=Integer.parseInt(split[9]);
		
		byte[] basesOriginal=split[10].getBytes();
		byte[] qualityOriginal=(split[11].equals(".") ? null : split[11].getBytes());
		
		if(qualityOriginal!=null){
			for(int i=0; i<qualityOriginal.length; i++){
				byte b=qualityOriginal[i];
				b=(byte) (b-ASCII_OFFSET);
				assert(b>=-1) : b;
				qualityOriginal[i]=b;
			}
		}
		
		int insert=-1;
		if(!split[12].equals(".")){insert=Integer.parseInt(split[12]);}
		
		byte[] match=null;
		if(!split[14].equals(".")){match=split[14].getBytes();}
		int[] gaps=null;
		if(!split[15].equals(".")){
			
			String[] gstring=split[16].split("~");
			gaps=new int[gstring.length];
			for(int i=0; i<gstring.length; i++){
				gaps[i]=Integer.parseInt(gstring[i]);
			}
		}
		
//		assert(false) : split[16];
		
		Read r=new Read(basesOriginal, qualityOriginal, id, numericID, flags, chrom, start, stop);
		r.match=match;
		r.errors=errors;
		r.mapScore=mapScore;
		r.copies=copies;
		r.gaps=gaps;
		r.insert=insert;
		
		int firstScore=(ADD_BEST_SITE_TO_LIST_FROM_TEXT) ? 17 : 18;
		
		int scores=split.length-firstScore;
		
		int mSites=0;
		for(int i=firstScore; i<split.length; i++){
			if(split[i].charAt(0)!='*'){mSites++;}
		}
		
		//This can be disabled to handle very old text format.
		if(mSites>0){r.sites=new ArrayList<SiteScore>(mSites);}
		for(int i=firstScore; i<split.length; i++){
			SiteScore ss=SiteScore.fromText(split[i]);
			if(split[i].charAt(0)=='*'){r.originalSite=ss;}
			else{r.sites.add(ss);}
		}
		
		if(DECOMPRESS_MATCH_ON_LOAD && r.shortmatch()){
			r.toLongMatchString(true);
		}

		assert(r.numSites()==0 || absdif(r.start, r.stop)<3000 || (r.gaps==null) == (r.topSite().gaps==null)) :
			"\n"+r.numericID+", "+r.chrom+", "+r.strand()+", "+r.start+", "+r.stop+", "+Arrays.toString(r.gaps)+"\n"+r.sites+"\n"+line+"\n";
		
		return r;
	}

	/** Inflates gaps between contigs in a scaffold. */
	public void inflateGaps(int minGapIn, int minGapOut) {
		assert(minGapIn>0);
		if(!containsNocalls()){return;}
		final ByteBuilder bbb=new ByteBuilder();
		final ByteBuilder bbq=(quality==null ? null : new ByteBuilder());
		
		int gap=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			byte q=(quality==null ? 0 : quality[i]);
			if(b=='N'){
				gap++;
			}else{
				while(gap>=minGapIn && gap<minGapOut){
					gap++;
					bbb.append('N');
					if(bbq!=null){bbq.append(0);}
				}
				gap=0;
			}
			bbb.append(b);
			if(bbq!=null){bbq.append(q);}
		}
		
		while(gap>=minGapIn && gap<minGapOut){//Handle trailing bases
			gap++;
			bbb.append('N');
			if(bbq!=null){bbq.append(0);}
		}
		
		assert(bbb.length()>=bases.length);
		if(bbb.length()>bases.length){
			bases=bbb.toBytes();
			if(bbq!=null){quality=bbq.toBytes();}
		}
	}
	
	public ArrayList<Read> breakAtGaps(final boolean agp, final int minContig){
		ArrayList<Read> list=new ArrayList<Read>();
		byte prev='N';
		int lastN=-1, lastBase=-1;
		int contignum=1;
		long feature=1;
		ByteBuilder bb=(agp ? new ByteBuilder() : null);
		assert(obj==null);
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			if(b=='N'){
				if(prev!='N'){
					final int start=lastN+1, stop=i;
					byte[] b2=KillSwitch.copyOfRange(bases, start, stop);
					byte[] q2=(quality==null ? null : KillSwitch.copyOfRange(quality, start, stop));
					Read r=new Read(b2, q2, id+"_c"+contignum, numericID);
					if(r.length()>=minContig){list.add(r);}
					contignum++;
					
					if(bb!=null){
						bb.append(id).append('\t');
						bb.append(start+1).append('\t');
						bb.append(stop).append('\t');
						bb.append(feature).append('\t');
						feature++;
						bb.append('W').append('\t');
						bb.append(r.id).append('\t');
						bb.append(1).append('\t');
						bb.append(r.length()).append('\t');
						bb.append('+').append('\n');
					}
				}
				lastN=i;
			}else{
				if(bb!=null && prev=='N' && lastBase>=0){
					bb.append(id).append('\t');
					bb.append(lastBase+2).append('\t');
					bb.append(i).append('\t');
					bb.append(feature).append('\t');
					feature++;
					bb.append('N').append('\t');
					bb.append((i-lastBase-1)).append('\t');
					bb.append("scaffold").append('\t');
					bb.append("yes").append('\t');
					bb.append("paired-ends").append('\n');
				}
				lastBase=i;
			}
			prev=b;
		}
		if(prev!='N'){
			final int start=lastN+1, stop=bases.length;
			byte[] b2=KillSwitch.copyOfRange(bases, start, stop);
			byte[] q2=(quality==null ? null : KillSwitch.copyOfRange(quality, start, stop));
			Read r=new Read(b2, q2, id+"_c"+contignum, numericID);
			if(r.length()>=minContig){list.add(r);}
			contignum++;
			
			if(bb!=null){
				bb.append(id).append('\t');
				bb.append(start+1).append('\t');
				bb.append(stop).append('\t');
				bb.append(feature).append('\t');
				feature++;
				bb.append('W').append('\t');
				bb.append(r.id).append('\t');
				bb.append(1).append('\t');
				bb.append(r.length()).append('\t');
				bb.append('+').append('\n');
			}
		}else{
			if(bb!=null && prev=='N' && lastBase>=0){
				bb.append(id).append('\t');
				bb.append(lastBase+2).append('\t');
				bb.append(bases.length).append('\t');
				bb.append(feature).append('\t');
				feature++;
				bb.append('N').append('\t');
				bb.append((bases.length-lastBase-1)).append('\t');
				bb.append("scaffold").append('\t');
				bb.append("yes").append('\t');
				bb.append("paired-ends").append('\n');
			}
			lastBase=bases.length;
		}
		if(bb!=null){obj=bb.toBytes();}
		return list;
	}

	/** Reverse-complements the read. */
	public Read reverseComplement() {
		AminoAcid.reverseComplementBasesInPlace(bases);
		Tools.reverseInPlace(quality);
		setStrand(strand()^1);
		return this;
	}

	/** Complements the read. */
	public void complement() {
		AminoAcid.reverseComplementBasesInPlace(bases);
	}
	
	@Override
	public int compareTo(Read o) {
		if(chrom!=o.chrom){return chrom-o.chrom;}
		if(start!=o.start){return start-o.start;}
		if(stop!=o.stop){return stop-o.stop;}
		if(strand()!=o.strand()){return strand()-o.strand();}
		return 0;
	}
	
	public SiteScore toSite(){
		assert(start<=stop) : this.toText(false);
		SiteScore ss=new SiteScore(chrom, strand(), start, stop, 0, 0, rescued(), perfect());
		if(paired()){
			ss.setSlowPairedScore(mapScore-1, mapScore);
		}else{
			ss.setSlowPairedScore(mapScore, 0);
		}
		ss.setScore(mapScore);
		ss.gaps=gaps;
		ss.match=match;
		originalSite=ss;
		return ss;
	}
	
	public SiteScore topSite(){
		final SiteScore ss=(sites==null || sites.isEmpty()) ? null : sites.get(0);
		assert(sites==null || sites.isEmpty() || ss!=null) : "Top site is null for read "+this;
		return ss;
	}
	
	public int numSites(){
		return (sites==null ? 0 : sites.size());
	}
	
	public SiteScore makeOriginalSite(){
		originalSite=toSite();
		return originalSite;
	}
	
	public void setFromSite(SiteScore ss){
		assert(ss!=null);
		chrom=ss.chrom;
		setStrand(ss.strand);
		start=ss.start;
		stop=ss.stop;
		mapScore=ss.slowScore;
		setRescued(ss.rescued);
		gaps=ss.gaps;
		setPerfect(ss.perfect);
		
		match=ss.match;
		
		if(gaps!=null){
			gaps=ss.gaps=GapTools.fixGaps(start, stop, gaps, Shared.MINGAP);
//			gaps[0]=Tools.min(gaps[0], start);
//			gaps[gaps.length-1]=Tools.max(gaps[gaps.length-1], stop);
		}
	}
	
//	public static int[] fixGaps(int a, int b, int[] gaps, int minGap){
////		System.err.println("fixGaps input: "+a+", "+b+", "+Arrays.toString(gaps)+", "+minGap);
//		int[] r=GapTools.fixGaps(a, b, gaps, minGap);
////		System.err.println("fixGaps output: "+Arrays.toString(r));
//		return r;
//	}

	public void setFromOriginalSite(){
		setFromSite(originalSite);
	}
	public void setFromTopSite(){
		final SiteScore ss=topSite();
		if(ss==null){
			clearSite();
			setMapped(false);
			return;
		}
		setMapped(true);
		setFromSite(ss);
	}
	
	public void setFromTopSite(boolean randomIfAmbiguous, boolean primary, int maxPairDist){
		final SiteScore ss0=topSite();
		if(ss0==null){
			clearSite();
			setMapped(false);
			return;
		}
		setMapped(true);
		
		if(sites.size()==1 || !randomIfAmbiguous || !ambiguous()){
			setFromSite(ss0);
			return;
		}
		
		if(primary || mate==null || !mate.mapped() || !mate.paired()){
			int count=1;
			for(int i=1; i<sites.size(); i++){
				SiteScore ss=sites.get(i);
				if(ss.score<ss0.score || (ss0.perfect && !ss.perfect) || (ss0.semiperfect && !ss.semiperfect)){break;}
				count++;
			}

			int x=(int)(numericID%count);
			if(x>0){
				SiteScore ss=sites.get(x);
				sites.set(0, ss);
				sites.set(x, ss0);
			}
			setFromSite(sites.get(0));
			return;
		}
		
//		assert(false) : "TODO: Proper strand orientation, and more.";
		//TODO: Also, this code appears to sometimes duplicate sitescores(?)
//		for(int i=0; i<list.size(); i++){
//			SiteScore ss=list.get(i);
//			if(ss.chrom==mate.chrom && Tools.min(Tools.absdifUnsigned(ss.start, mate.stop), Tools.absdifUnsigned(ss.stop, mate.start))<=maxPairDist){
//				list.set(0, ss);
//				list.set(i, ss0);
//				setFromSite(ss);
//				return;
//			}
//		}
		
		//If unsuccessful, recur unpaired.
		
		this.setPaired(false);
		mate.setPaired(false);
		setFromTopSite(randomIfAmbiguous, true, maxPairDist);
	}
	
	public void clearPairMapping(){
		clearMapping();
		if(mate!=null){mate.clearMapping();}
	}
	
	public void clearMapping(){
		clearSite();
		match=null;
		sites=null;
		setMapped(false);
		setPaired(false);
		if(mate!=null){mate.setPaired(false);}
	}
	
	public void clearSite(){
		chrom=-1;
		setStrand(0);
		start=-1;
		stop=-1;
//		errors=0;
		mapScore=0;
		gaps=null;
	}


	public void clearAnswers(boolean clearMate) {
//		assert(mate==null || (pairnum()==0 && mate.pairnum()==1)) : pairnum()+", "+mate.pairnum();
		clearSite();
		match=null;
		sites=null;
		flags=(flags&(SYNTHMASK|PAIRNUMMASK|SWAPMASK));
		if(clearMate && mate!=null){
			mate.clearSite();
			mate.match=null;
			mate.sites=null;
			mate.flags=(mate.flags&(SYNTHMASK|PAIRNUMMASK|SWAPMASK));
		}
//		assert(mate==null || (pairnum()==0 && mate.pairnum()==1)) : pairnum()+", "+mate.pairnum();
	}
	
	
	public boolean isBadPair(boolean requireCorrectStrands, boolean sameStrandPairs, int maxdist){
		if(mate==null || paired()){return false;}
		if(!mapped() || !mate.mapped()){return false;}
		if(chrom!=mate.chrom){return true;}
		
		{
			int inner;
			if(start<=mate.start){inner=mate.start-stop;}
			else{inner=start-mate.stop;}
			if(inner>maxdist){return true;}
		}
//		if(absdif(start, mate.start)>maxdist){return true;}
		if(requireCorrectStrands){
			if((strand()==mate.strand())!=sameStrandPairs){return true;}
		}
		if(!sameStrandPairs){
			if(strand()==Shared.PLUS && mate.strand()==Shared.MINUS){
				if(start>=mate.stop){return true;}
			}else if(strand()==Shared.MINUS && mate.strand()==Shared.PLUS){
				if(mate.start>=stop){return true;}
			}
		}
		return false;
	}
	
	public int countMismatches(){
		assert(match!=null);
		int x=0;
		for(byte b : match){
			if(b=='S'){x++;}
		}
		return x;
	}

	/**
	 * @param k
	 * @return Number of valid kmers
	 */
	public int numValidPairKmers(int k) {
		return numValidKmers(k)+(mate==null ? 0 : mate.numValidKmers(k));
	}

	/**
	 * @param k
	 * @return Number of valid kmers
	 */
	public int numValidKmers(int k) {
		if(bases==null){return 0;}
		int len=0, counted=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			if(x<0){len=0;}else{len++;}
			if(len>=k){counted++;}
		}
		return counted;
	}
	
	/**
	 * @param match string
	 * @return Total number of match, sub, del, ins, or clip symbols
	 */
	public static final int[] matchToMsdicn(byte[] match) {
		if(match==null || match.length<1){return null;}
		int[] msdicn=KillSwitch.allocInt1D(6);
		
		byte mode='0', c='0';
		int current=0;
		for(int i=0; i<match.length; i++){
			c=match[i];
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(mode==c){
					current=Tools.max(current+1, 2);
				}else{
					current=Tools.max(current, 1);

					if(mode=='m'){
						msdicn[0]+=current;
					}else if(mode=='S'){
						msdicn[1]+=current;
					}else if(mode=='D'){
						msdicn[2]+=current;
					}else if(mode=='I'){
						msdicn[3]+=current;
					}else if(mode=='C' || mode=='X' || mode=='Y'){
						msdicn[4]+=current;
					}else if(mode=='N' || mode=='R'){
						msdicn[5]+=current;
					}
					mode=c;
					current=0;
				}
			}
		}
		if(current>0 || !Tools.isDigit(c)){
			current=Tools.max(current, 1);
			if(mode=='m'){
				msdicn[0]+=current;
			}else if(mode=='S'){
				msdicn[1]+=current;
			}else if(mode=='D'){
				msdicn[2]+=current;
			}else if(mode=='I'){
				msdicn[3]+=current;
			}else if(mode=='C' || mode=='X' || mode=='Y'){
				msdicn[4]+=current;
			}else if(mode=='N' || mode=='R'){
				msdicn[5]+=current;
			}
		}
		return msdicn;
	}

	
	/**
	 * @param match string
	 * @return Ref length of match string
	 */
	public static final int calcMatchLength(byte[] match) {
		if(match==null || match.length<1){return 0;}
		
		byte mode='0', c='0';
		int current=0;
		int len=0;
		for(int i=0; i<match.length; i++){
			c=match[i];
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(mode==c){
					current=Tools.max(current+1, 2);
				}else{
					current=Tools.max(current, 1);

					if(mode=='m'){
						len+=current;
					}else if(mode=='S'){
						len+=current;
					}else if(mode=='D'){
						len+=current;
					}else if(mode=='I'){ //Do nothing
						//len+=current;
					}else if(mode=='C'){
						len+=current;
					}else if(mode=='X'){ //Not sure about this case, but adding seems fine
						len+=current;
//						assert(false) : new String(match);
					}else if(mode=='Y'){ //Do nothing
						//len+=current;
//						assert(false) : new String(match);
					}else if(mode=='N' || mode=='R'){
						len+=current;
					}
					mode=c;
					current=0;
				}
			}
		}
		if(current>0 || !Tools.isDigit(c)){
			current=Tools.max(current, 1);
			if(mode=='m'){
				len+=current;
			}else if(mode=='S'){
				len+=current;
			}else if(mode=='D'){
				len+=current;
			}else if(mode=='I'){ //Do nothing
				//len+=current;
			}else if(mode=='C'){
				len+=current;
			}else if(mode=='X'){ //Not sure about this case, but adding seems fine
				len+=current;
				assert(false) : new String(match);
			}else if(mode=='Y'){ //Do nothing
				//len+=current;
//				assert(false) : new String(match);
			}else if(mode=='N' || mode=='R'){
				len+=current;
			}
		}
		return len;
	}

	public final float identity() {return identity(match);}
	
	public static final float identity(byte[] match) {
		if(FLAT_IDENTITY){
			return identityFlat(match, true);
		}else{
			return identitySkewed(match, true, true, false, false);
		}
	}
	
	public final boolean hasLongInsertion(int maxlen){
		return hasLongInsertion(match, maxlen);
	}
	
	public final boolean hasLongDeletion(int maxlen){
		return hasLongDeletion(match, maxlen);
	}
	
	public static final boolean hasLongInsertion(byte[] match, int maxlen){
		if(match==null || match.length<maxlen){return false;}
		byte prev='0';
		int len=0;
		for(byte b : match){
			if(b=='I' || b=='X' || b=='Y'){
				if(b==prev){len++;}
				else{len=1;}
				if(len>maxlen){return true;}
			}else{
				len=0;
			}
			prev=b;
		}
		return false;
	}
	
	public static final boolean hasLongDeletion(byte[] match, int maxlen){
		if(match==null || match.length<maxlen){return false;}
		byte prev='0';
		int len=0;
		for(byte b : match){
			if(b=='D'){
				if(b==prev){len++;}
				else{len=1;}
				if(len>maxlen){return true;}
			}else{
				len=0;
			}
			prev=b;
		}
		return false;
	}
	
	public int mappedNonClippedBases() {
		if(!mapped() || match==null || bases==null){return 0;}
		
		int len=0;
		byte mode='0', c='0';
		int current=0;
		for(int i=0; i<match.length; i++){
			c=match[i];
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(mode==c){
					current=Tools.max(current+1, 2);
				}else{
					current=Tools.max(current, 1);

					if(mode=='D' || mode=='C' || mode=='X' || mode=='Y' || mode=='d'){
						
					}else{
						len+=current;
					}
					mode=c;
					current=0;
				}
			}
		}
		if(current>0 || !Tools.isDigit(c)){
			current=Tools.max(current, 1);
			if(mode=='D' || mode=='C' || mode=='X' || mode=='Y' || mode=='d'){
				
			}else{
				len+=current;
			}
			mode=c;
			current=0;
		}
		return len;
	}
	
	/**
	 * Handles short or long mode.
	 * @param match string
	 * @return Identity based on number of match, sub, del, ins, or N symbols
	 */
	public static final float identityFlat(byte[] match, boolean penalizeN) {
//		assert(false) : new String(match);
		if(match==null || match.length<1){return 0;}
		
		int good=0, bad=0, n=0;
		
		byte mode='0', c='0';
		int current=0;
		for(int i=0; i<match.length; i++){
			c=match[i];
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(mode==c){
					current=Tools.max(current+1, 2);
				}else{
					current=Tools.max(current, 1);

					if(mode=='m'){
						good+=current;
//						System.out.println("G: mode="+(char)mode+", c="+(char)c+", current="+current+", good="+good+", bad="+bad);
					}else if(mode=='R' || mode=='N'){
						n+=current;
					}else if(mode=='C' || mode=='V'){
						//Do nothing
						//I assume this is clipped because it went off the end of a scaffold, and thus is irrelevant to identity
					}else if(mode!='0'){
						assert(mode=='S' || mode=='D' || mode=='I' || mode=='X' || mode=='Y' || mode=='i' || mode=='d') : (char)mode;
						if(mode!='D' || current<SamLine.INTRON_LIMIT){
							bad+=current;
						}
//						System.out.println("B: mode="+(char)mode+", c="+(char)c+", current="+current+", good="+good+", bad="+bad);
					}
					mode=c;
					current=0;
				}
			}
		}
		if(current>0 || !Tools.isDigit(c)){
			current=Tools.max(current, 1);
			if(mode=='m'){
				good+=current;
			}else if(mode=='R' || mode=='N'){
				n+=current;
			}else if(mode=='C' || mode=='V'){
				//Do nothing
				//I assume this is clipped because it went off the end of a scaffold, and thus is irrelevant to identity
			}else if(mode!='0'){
				assert(mode=='S' || mode=='I' || mode=='X' || mode=='Y') : (char)mode;
				if(mode!='D' || current<SamLine.INTRON_LIMIT){
					bad+=current;
				}
//				System.out.println("B: mode="+(char)mode+", c="+(char)c+", current="+current+", good="+good+", bad="+bad);
			}
		}
		

		float good2=good+n*(penalizeN ? 0.25f : 0);
		float bad2=bad+n*(penalizeN ? 0.75f : 0);
		float r=good2/Tools.max(good2+bad2, 1);
//		assert(false) : new String(match)+"\nmode='"+(char)mode+"', current="+current+", good="+good+", bad="+bad;

//		System.out.println("match="+new String(match)+"\ngood="+good+", bad="+bad+", r="+r);
//		System.out.println(Arrays.toString(matchToMsdicn(match)));
		
		return r;
	}
	
	/**
	 * Handles short or long mode.
	 * @param match string
	 * @return Identity based on number of match, sub, del, ins, or N symbols
	 */
	public static final float identitySkewed(byte[] match, boolean penalizeN, boolean sqrt, boolean log, boolean single) {
//		assert(false) : new String(match);
		if(match==null || match.length<1){return 0;}
		
		int good=0, bad=0, n=0;
		
		byte mode='0', c='0';
		int current=0;
		for(int i=0; i<match.length; i++){
			c=match[i];
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(mode==c){
					current=Tools.max(current+1, 2);
				}else{
					current=Tools.max(current, 1);

					if(mode=='m'){
						good+=current;
//						System.out.println("G: mode="+(char)mode+", c="+(char)c+", current="+current+", good="+good+", bad="+bad);
					}else if(mode=='D'){
						if(current<SamLine.INTRON_LIMIT){
							int x;
							
							if(sqrt){x=(int)Math.ceil(Math.sqrt(current));}
							else if(log){x=(int)Math.ceil(Tools.log2(current));}
							else{x=1;}
							
							bad+=(Tools.min(x, current));
						}
						
//						System.out.println("D: mode="+(char)mode+", c="+(char)c+", current="+current+", good="+good+", bad="+bad+", x="+x);
					}else if(mode=='R' || mode=='N'){
						n+=current;
					}else if(mode=='C' || mode=='V'){
						//Do nothing
						//I assume this is clipped because it went off the end of a scaffold, and thus is irrelevant to identity
					}else if(mode!='0'){
						assert(mode=='S' || mode=='I' || mode=='X' || mode=='Y' || mode=='i' || mode=='d') : (char)mode;
						bad+=current;
//						System.out.println("B: mode="+(char)mode+", c="+(char)c+", current="+current+", good="+good+", bad="+bad);
					}
					mode=c;
					current=0;
				}
			}
		}
		if(current>0 || !Tools.isDigit(c)){
			current=Tools.max(current, 1);
			if(mode=='m'){
				good+=current;
			}else if(mode=='R' || mode=='N'){
				n+=current;
			}else if(mode=='C' || mode=='V'){
				//Do nothing
				//I assume this is clipped because it went off the end of a scaffold, and thus is irrelevant to identity
			}else if(mode!='0'){
				assert(mode=='S' || mode=='I' || mode=='X' || mode=='Y') : (char)mode;
				if(mode!='D' || current<SamLine.INTRON_LIMIT){
					bad+=current;
				}
//				System.out.println("B: mode="+(char)mode+", c="+(char)c+", current="+current+", good="+good+", bad="+bad);
			}
		}
		
		
		float good2=good+n*(penalizeN ? 0.25f : 0);
		float bad2=bad+n*(penalizeN ? 0.75f : 0);
		float r=good2/Tools.max(good2+bad2, 1);
//		assert(false) : new String(match)+"\nmode='"+(char)mode+"', current="+current+", good="+good+", bad="+bad;

//		System.out.println("match="+new String(match)+"\ngood="+good+", bad="+bad+", r="+r);
//		System.out.println(Arrays.toString(matchToMsdicn(match)));
		
		return r;
	}
	
	public boolean failsChastity(){
		return failsChastity(true);
	}
	
	public boolean failsChastity(boolean processAssertions){
		if(id==null){return false;}
		int space=id.indexOf(' ');
		if(space<0 || space+5>id.length()){return false;}
		char a=id.charAt(space+1);
		char b=id.charAt(space+2);
		char c=id.charAt(space+3);
		char d=id.charAt(space+4);
		
		if(a=='/'){
			if(b<'1' || b>'4' || c!=':'){
				if(!processAssertions){return false;}
				KillSwitch.kill("Strangely formatted read.  Please disable chastityfilter with the flag chastityfilter=f.  id:"+id);
			}
			return d=='Y';
		}else{
			if(processAssertions){
				assert(a=='1' || a=='2' || a=='3' || a=='4') : id;
				assert(b==':') : id;
				assert(d==':');
			}
			if(a<'1' || a>'4' || b!=':' || d!=':'){
				if(!processAssertions){return false;}
				KillSwitch.kill("Strangely formatted read.  Please disable chastityfilter with the flag chastityfilter=f.  id:"+id);
			}
			return c=='Y';
		}
	}
	
	public boolean failsBarcode(HashSet<String> set, boolean failIfNoBarcode){
		if(id==null){return false;}
		
		final int loc=id.lastIndexOf(':');
		final int loc2=Tools.max(id.indexOf(' '), id.indexOf('/'));
		if(loc<0 || loc<=loc2 || loc>=id.length()-1){
			return failIfNoBarcode;
		}
		
		if(set==null){
			for(int i=loc+1; i<id.length(); i++){
				char c=id.charAt(i);
				boolean ok=(c=='+' || AminoAcid.isFullyDefined(c));
				if(!ok){return true;}
			}
			return false;
		}else{
			String code=id.substring(loc+1);
			return !set.contains(code);
		}
	}

	public String barcode(boolean failIfNoBarcode){
		return headerToBarcode(id, failIfNoBarcode);
	}
	
	/** 
	 * Parse the barcode from an Illumina header.
	 * @param failIfNoBarcode Terminate the JVM if no barcode is present. 
	 * @return The barcode
	 */
	public static String headerToBarcode(String id, boolean failIfNoBarcode){
		
		if(id==null){
			if(failIfNoBarcode && Shared.EA()){
				KillSwitch.kill("Encountered a read header without a barcode:\n"+id+"\n");
			}
			return null;
		}
		
		final int loc=id.lastIndexOf(':');
		final int loc2=Tools.max(id.indexOf(' '), id.indexOf('/'));
		if(loc<0 || loc<=loc2 || loc>=id.length()-1){
			if(failIfNoBarcode && Shared.EA()){
				KillSwitch.kill("Encountered a read header without a barcode:\n"+id+"\n");
			}
			return null;
		}
		
		//This section allows for comments after the barcode
		final int bcStart=loc+1;
		int bcStop=id.indexOf(' ', bcStart);
		bcStop=(bcStop>=0 ? bcStop : id.indexOf('\t', bcStart));
		String code;
		if(bcStop<0) {
			code=id.substring(bcStart);
		}else {
			code=id.substring(bcStart, bcStop);
		}
		
		return code;
	}
	
	/**
	 * @return The rname of this Read's SamLine, if present and mapped.
	 */
	public String rnameS() {
		if(samline==null || !samline.mapped()){return null;}
		return samline.rnameS();
	}

	/** Average based on summing quality scores */
	public double avgQuality(boolean countUndefined, int maxBases){
		return AVERAGE_QUALITY_BY_PROBABILITY ? avgQualityByProbabilityDouble(countUndefined, maxBases) : avgQualityByScoreDouble(maxBases);
	}

	/** Average based on summing quality scores */
	public int avgQualityInt(boolean countUndefined, int maxBases){
		return AVERAGE_QUALITY_BY_PROBABILITY ? avgQualityByProbabilityInt(countUndefined, maxBases) : avgQualityByScoreInt(maxBases);
	}
	
	/** Average based on summing error probabilities */
	public int avgQualityByProbabilityInt(boolean countUndefined, int maxBases){
		if(bases==null || bases.length==0){return 0;}
		return avgQualityByProbabilityInt(bases, quality, countUndefined, maxBases);
	}
	
	/** Average based on summing error probabilities */
	public double avgQualityByProbabilityDouble(boolean countUndefined, int maxBases){
		if(bases==null || bases.length==0){return 0;}
		return avgQualityByProbabilityDouble(bases, quality, countUndefined, maxBases);
	}
	
	/** Average based on summing error probabilities */
	public double probabilityErrorFree(boolean countUndefined, int maxBases){
		if(bases==null || bases.length==0){return 0;}
		return probabilityErrorFree(bases, quality, countUndefined, maxBases);
	}
	
	/** Average based on summing error probabilities */
	public static int avgQualityByProbabilityInt(byte[] bases, byte[] quality, boolean countUndefined, int maxBases){
		if(quality==null){return 40;}
		if(quality.length==0){return 0;}
		float e=expectedErrors(bases, quality, countUndefined, maxBases);
		final int div=(maxBases<1 ? quality.length : Tools.min(maxBases, quality.length));
		float p=e/div;
		return QualityTools.probErrorToPhred(p);
	}
	
	/** Average based on summing error probabilities */
	public static double avgQualityByProbabilityDouble(byte[] bases, byte[] quality, boolean countUndefined, int maxBases){
		if(quality==null){return 40;}
		if(quality.length==0){return 0;}
		float e=expectedErrors(bases, quality, countUndefined, maxBases);
		final int div=(maxBases<1 ? quality.length : Tools.min(maxBases, quality.length));
		float p=e/div;
		return QualityTools.probErrorToPhredDouble(p);
	}

	/** Average based on summing quality scores */
	public int avgQualityByScoreInt(int maxBases){
		if(bases==null || bases.length==0){return 0;}
		if(quality==null){return 40;}
		int x=0, limit=(maxBases<1 ? quality.length : Tools.min(maxBases, quality.length));
		for(int i=0; i<limit; i++){
			byte b=quality[i];
			x+=(b<0 ? 0 : b);
		}
		return x/limit;
	}

	/** Average based on summing quality scores */
	public double avgQualityByScoreDouble(int maxBases){
		if(bases==null || bases.length==0){return 0;}
		if(quality==null){return 40;}
		int x=0, limit=(maxBases<1 ? quality.length : Tools.min(maxBases, quality.length));
		for(int i=0; i<limit; i++){
			byte b=quality[i];
			x+=(b<0 ? 0 : b);
		}
		return x/(double)limit;
	}
	
	/** Used by BBMap tipsearch. */
	public int avgQualityFirstNBases(int n){
		if(bases==null || bases.length==0){return 0;}
		if(quality==null || n<1){return 40;}
		assert(quality!=null);
		int x=0;
		if(n>quality.length){return 0;}
		for(int i=0; i<n; i++){
			byte b=quality[i];
			x+=(b<0 ? 0 : b);
		}
		return x/n;
	}
	
	/** Used by BBMap tipsearch. */
	public int avgQualityLastNBases(int n){
		if(bases==null || bases.length==0){return 0;}
		if(quality==null || n<1){return 40;}
		assert(quality!=null);
		int x=0;
		if(n>quality.length){return 0;}
		for(int i=bases.length-n; i<bases.length; i++){
			byte b=quality[i];
			x+=(b<0 ? 0 : b);
		}
		return x/n;
	}
	
	/** Used by BBDuk. */
	public int minQuality(){
		byte min=41;
		if(bases!=null && quality!=null){
			for(byte q : quality){
				min=Tools.min(min, q);
			}
		}
		return min;
	}
	
	/** Used by BBMap tipsearch. */
	public byte minQualityFirstNBases(int n){
		if(bases==null || bases.length==0){return 0;}
		if(quality==null || n<1){return 41;}
		assert(quality!=null && n>0);
		if(n>quality.length){return 0;}
		byte x=quality[0];
		for(int i=1; i<n; i++){
			byte b=quality[i];
			if(b<x){x=b;}
		}
		return x;
	}
	
	/** Used by BBMap tipsearch. */
	public byte minQualityLastNBases(int n){
		if(bases==null || bases.length==0){return 0;}
		if(quality==null || n<1){return 41;}
		assert(quality!=null && n>0);
		if(n>quality.length){return 0;}
		byte x=quality[bases.length-n];
		for(int i=bases.length-n; i<bases.length; i++){
			byte b=quality[i];
			if(b<x){x=b;}
		}
		return x;
	}
	
	public boolean containsNonM(){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b>'9' && b!='m'){return true;}
		}
		return false;
	}
	
	public boolean containsNonNM(){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b>'9' && b!='m' && b!='N'){return true;}
		}
		return false;
	}
	
	public boolean containsVariants(){
		assert(match!=null && valid()) : (match==null)+", "+(valid())+"\n"+samline+"\n";
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b>'9' && b!='m' && b!='N' && b!='C'){return true;}
		}
		return false;
	}
	
	public boolean containsClipping(){
		assert(match!=null && valid()) : (match==null)+", "+(valid())+"\n"+samline+"\n";
		if(match.length<1){return false;}
		if(match[0]=='C'){return true;}
		for(int i=match.length-1; i>0; i--){
			if(match[i]=='C'){return true;}
			if(match[i]>'9'){break;}
		}
		return false;
	}
	
	public int countAlignedBases() {
		return countAlignedBases(match);
	}
	
	public static int countAlignedBases(byte[] match) {
		//Short version; correct but slow
//		int[] sym=countMatchSymbols(match);
//		return sym[0]+sym[1]+sym[4];
		
		//Long version; faster
		int msi=0;
		int current=0;
		byte last='?';
		for(byte b : match){
			if(Tools.isDigit(b)){
				current=current*10+b-'0';
			}else{
				current=Tools.max(current, 1);
				if(last=='m' || last=='S' || last=='I'){
					msi+=current;
				}
				else{
					//Ignore
					assert(last=='?' || last=='C' || last=='N' || last=='D') : "Unhandled symbol "+(char)last+"\n"+new String(match);
				}
				current=0;
				last=b;
			}
		}
		current=Tools.max(current, 1);
		if(last=='m' || last=='S' || last=='I'){
			msi+=current;
		}
		else{
			//Ignore
			assert(last=='?' || last=='C' || last=='N' || last=='D') : "Unhandled symbol "+(char)last+"\n"+new String(match);
		}
		current=0;
		return msi;
	}
	
	//This ignores N and clip, and counts each deletion event as a single error
	public int countErrors(){
		return countErrors(match);
	}
	
	//This ignores N and clip, and counts each deletion event as a single error
	public static int countErrors(byte[] match){
		int m=0, S=0, C=0, N=0, I=0, D=0;
		int current=0;
		byte last='?';
		for(byte b : match){
			if(Tools.isDigit(b)){
				current=current*10+b-'0';
			}else{
				current=Tools.max(current, 1);
				if(last=='m'){
					m+=current;
				}else if(last=='S'){
					S+=current;
				}else if(last=='C'){
					C+=current;
				}else if(last=='N'){
					N+=current;
				}else if(last=='I'){
					I+=current;
				}else if(last=='D'){
					D++;
				}else{
					assert(last=='?') : "Unhandled symbol "+(char)last+"\n"+new String(match);
				}
				current=0;
				last=b;
			}
		}
		current=Tools.max(current, 1);
		if(last=='m'){
			m+=current;
		}else if(last=='S'){
			S+=current;
		}else if(last=='C'){
			C+=current;
		}else if(last=='N'){
			N+=current;
		}else if(last=='I'){
			I+=current;
		}else if(last=='D'){
			D++;
		}else{
			assert(last=='?') : "Unhandled symbol "+(char)last+"\n"+new String(match);
		}
		current=0;
		int errors=S+I+D;
//		if(errors>30) {System.err.println(errors+": "+new String(match));}
		return errors;
	}
	
	/**
	 * @return {m,S,C,N,I,D};
	 */
	public int[] countMatchSymbols(){
		return countMatchSymbols(match);
	}
	
	/**
	 * @return {m,S,C,N,I,D};
	 */
	public static int[] countMatchSymbols(byte[] match){
		int m=0, S=0, C=0, N=0, I=0, D=0;
		int current=0;
		byte last='?';
		for(byte b : match){
			if(Tools.isDigit(b)){
				current=current*10+b-'0';
			}else{
				current=Tools.max(current, 1);
				if(last=='m'){
					m+=current;
				}else if(last=='S'){
					S+=current;
				}else if(last=='C'){
					C+=current;
				}else if(last=='N'){
					N+=current;
				}else if(last=='I'){
					I+=current;
				}else if(last=='D'){
					D+=current;
				}else{
					assert(last=='?') : "Unhandled symbol "+(char)last+"\n"+new String(match);
				}
				current=0;
				last=b;
			}
		}
		current=Tools.max(current, 1);
		if(last=='m'){
			m+=current;
		}else if(last=='S'){
			S+=current;
		}else if(last=='C'){
			C+=current;
		}else if(last=='N'){
			N+=current;
		}else if(last=='I'){
			I+=current;
		}else if(last=='D'){
			D+=current;
		}else{
			assert(last=='?') : "Unhandled symbol "+(char)last+"\n"+new String(match);
		}
		current=0;
		return new int[] {m,S,C,N,I,D};
	}
	
	/**
	 * Here, consecutive symbols are collapsed, so mmmDDDmmmm would yield 2 m and 1 D.
	 * @return {m,S,C,N,I,D};
	 */
	public static int[] countMatchEvents(byte[] match){
		int m=0, S=0, C=0, N=0, I=0, D=0;
		int current=0;
		byte last='?';
		for(byte b : match){
			if(Tools.isDigit(b)){
				current=current*10+b-'0';
			}else{
				current=Tools.max(current, 1);
				if(last=='m'){
					m++;
				}else if(last=='S'){
					S++;
				}else if(last=='C'){
					C++;
				}else if(last=='N'){
					N++;
				}else if(last=='I'){
					I++;
				}else if(last=='D'){
					D++;
				}else{
					assert(last=='?') : "Unhandled symbol "+(char)last+"\n"+new String(match);
				}
				current=0;
				last=b;
			}
		}
		current=Tools.max(current, 1);
		if(last=='m'){
			m++;
		}else if(last=='S'){
			S++;
		}else if(last=='C'){
			C++;
		}else if(last=='N'){
			N++;
		}else if(last=='I'){
			I++;
		}else if(last=='D'){
			D++;
		}else{
			assert(last=='?') : "Unhandled symbol "+(char)last+"\n"+new String(match);
		}
		current=0;
		return new int[] {m,S,C,N,I,D};
	}
	
	public boolean containsNonNMXY(){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b>'9' && b!='m' && b!='N' && b!='X' && b!='Y'){return true;}
		}
		return false;
	}
	
	public boolean containsSDI(){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b=='S' || b=='s' || b=='D' || b=='I'){return true;}
		}
		return false;
	}
	
	public boolean containsNonNMS(){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b>'9' && b!='m' && b!='s' && b!='N' && b!='S'){return true;}
		}
		return false;
	}
	
	public boolean containsConsecutiveS(int num){
		assert(match!=null && valid() && !shortmatch());
		int cnt=0;
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b=='S'){
				cnt++;
				if(cnt>=num){return true;}
			}else{
				cnt=0;
			}
		}
		return false;
	}
	
	public boolean containsIndels(){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='I' || b=='D' || b=='X' || b=='Y'){return true;}
		}
		return false;
	}
	
	public int countSubs(){
		assert(match!=null && valid()) : (match!=null)+", "+valid()+", "+shortmatch();
		return countSubs(match);
	}
	
	public boolean containsInMatch(char c){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b==c){return true;}
		}
		return false;
	}
	
	public boolean containsNocalls(){
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			if(b=='N'){return true;}
		}
		return false;
	}
	
	public int countNocalls(){
		return countNocalls(bases);
	}
	
	public int longestHomopolymer(){
		return longestHomopolymer(bases);
	}
	
	public static int countSubs(byte[] match){
		int S=0;
		int current=0;
		byte last='?';
		for(byte b : match){
			if(Tools.isDigit(b)){
				current=current*10+b-'0';
			}else{
				if(last=='S'){S+=Tools.max(1, current);}
				current=0;
				last=b;
			}
		}
		if(last=='S'){S+=Tools.max(1, current);}
//		assert(S==0) : S+"\t"+new String(match);
		return S;
//		int x=0;
//		assert(match!=null);
//		for(int i=0; i<match.length; i++){
//			byte b=match[i];
//			if(b=='S'){x++;}
//			assert(!Tools.isDigit(b));
//		}
//		return x;
	}
	
	public static int countVars(byte[] match){
		return countVars(match, true, true, true);
	}
	
	public static int countVars(byte[] match, boolean sub, boolean ins, boolean del){
		int S=0, I=0, D=0;
		int current=0;
		byte last='?';
		for(byte b : match){
			if(Tools.isDigit(b)){
				current=current*10+b-'0';
			}else{
				if(last=='S'){S+=Tools.max(1, current);}
				else if(last=='I'){I++;}
				else if(last=='D'){D++;}
				current=0;
				last=b;
			}
		}
		if(last=='S'){S+=Tools.max(1, current);}
		else if(last=='I'){I++;}
		else if(last=='D'){D++;}
		return (sub ? S : 0)+(ins ? I : 0)+(del ? D : 0);
	}
	
	public static boolean containsSubs(byte[] match){
		int x=0;
		assert(match!=null);
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='S'){return true;}
		}
		return false;
	}
	
	public static boolean containsVars(byte[] match){
		int x=0;
		assert(match!=null);
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='S' || b=='I' || b=='D'){return true;}
		}
		return false;
	}
	
	public static int countNocalls(byte[] match){
		int n=0;
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='N'){n++;}
		}
		return n;
	}
	
	public static int countUndefined(byte[] bases){
		return countUndefined(bases, 0, bases.length-1);
	}
	
	public static int countUndefined(byte[] bases, int from, int to){
		int n=0;
		if(bases==null) {return 0;}
		for(int i=from; i<=to; i++){
			n+=(AminoAcid.isFullyDefined(bases[i]) ? 0 : 1);
		}
		return n;
	}
	
	public static int longestHomopolymer(byte[] bases){
		return longestHomopolymer(bases, 0, bases.length-1);
	}
	
	public static int longestHomopolymer(byte[] bases, int from, int to){
		int max=0, streak=0;
		if(bases==null) {return 0;}
		for(int i=from, prev=-1; i<=to; i++){
			final byte b=bases[i];
			if(b==prev && AminoAcid.isFullyDefined(b)) {
				streak++;
			}else {
				max=Tools.max(max, streak);
				streak=1;
				prev=b;
			}
		}
		return Tools.max(max, streak);
	}
	
	public static int countInsertions(byte[] match){
		int n=0;
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='I'){n++;}
		}
		return n;
	}
	
	public static int countDeletions(byte[] match){
		int n=0;
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='D'){n++;}
		}
		return n;
	}
	
	public static int countInsertionEvents(byte[] match){
		int n=0;
		byte prev='N';
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='I' && prev!=b){n++;}
			prev=b;
		}
		return n;
	}
	
	public static int countDeletionEvents(byte[] match){
		int n=0;
		byte prev='N';
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='D' && prev!=b){n++;}
			prev=b;
		}
		return n;
	}
	
	public boolean containsNonACGTN(){
		if(bases==null){return false;}
		for(byte b : bases){
			if(AminoAcid.baseToNumberACGTN[b]<0){return true;}
		}
		return false;
	}
	
	public boolean containsUndefined(){
		if(bases==null){return false;}
		final byte[] symbolToNumber=AminoAcid.symbolToNumber(amino());
		for(byte b : bases){
			if(symbolToNumber[b]<0){return true;}
		}
		return false;
	}
	
	public boolean containsLowercase(){
		if(bases==null){return false;}
		for(byte b : bases){
			if(Tools.isLowerCase(b)){return true;}
		}
		return false;
	}
	
	public int countUndefined(){
		if(bases==null){return 0;}
		final byte[] symbolToNumber=AminoAcid.symbolToNumber(amino());
		int n=0;
		for(byte b : bases){
			n+=(symbolToNumber[b]>=0 ? 0 : 1);
//			if(symbolToNumber[b]<0){n++;}
		}
		return n;
	}
	
	public boolean hasMinConsecutiveBases(final int min){
		if(bases==null){return min<=0;}
		final byte[] symbolToNumber=AminoAcid.symbolToNumber(amino());
		int len=0;
		for(byte b : bases){
			if(symbolToNumber[b]<0){len=0;}
			else{
				len++;
				if(len>=min){return true;}
			}
		}
		return false;
	}
	
	
	/**
	 * @return The number of occurrences of the rarest base.
	 */
	public int minBaseCount(){
		if(bases==null){return 0;}
		int a=0, c=0, g=0, t=0;
		for(byte b : bases){
			if(b=='A'){a++;}
			else if(b=='C'){c++;}
			else if(b=='G'){g++;}
			else if(b=='T'){t++;}
		}
		return Tools.min(a, c, g, t);
	}
	
	public boolean containsXY(){
		assert(match!=null && valid());
		return containsXY(match);
	}
	
	public static boolean containsXY(byte[] match){
		if(match==null){return false;}
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='X' || b=='Y'){return true;}
		}
		return false;
	}
	
	public boolean containsXY2(){
		if(match==null || match.length<1){return false;}
		boolean b=(match[0]=='X' || match[match.length-1]=='Y');
		assert(!valid() || b==containsXY());
		return b;
	}
	
	public boolean containsXYC(){
		if(match==null || match.length<1){return false;}
		boolean b=(match[0]=='X' || match[match.length-1]=='Y');
		assert(!valid() || b==containsXY());
		return b || match[0]=='C' || match[match.length-1]=='C';
	}
	
	/** Replaces 'B' in match string with 'S', 'm', or 'N' */
	public boolean fixMatchB(){
		assert(match!=null);
		final ChromosomeArray ca;
		if(Data.GENOME_BUILD>=0){
			ca=Data.getChromosome(chrom);
		}else{
			ca=null;
		}
		boolean originallyShort=shortmatch();
		if(originallyShort){match=toLongMatchString(match);}
		int mloc=0, cloc=0, rloc=start;
		for(; mloc<match.length; mloc++){
			byte m=match[mloc];
			
			if(m=='B'){
				byte r=(ca==null ? (byte)'?' : ca.get(rloc));
				byte c=bases[cloc];
				if(r=='N' || c=='N'){
					match[mloc]='N';
				}else if(r==c || Tools.toUpperCase(r)==Tools.toUpperCase(c)){
					match[mloc]='m';
				}else{
					if(ca==null){
						if(originallyShort){
							match=toShortMatchString(match);
						}
						for(int i=0; i<match.length; i++){
							if(match[i]=='B'){match[i]='N';}
						}
						return false;
					}
					match[mloc]='S';
				}
				cloc++;
				rloc++;
			}else if(m=='m' || m=='S' || m=='N' || m=='s' || m=='C'){
				cloc++;
				rloc++;
			}else if(m=='D'){
				rloc++;
			}else if(m=='I' || m=='X' || m=='Y'){
				cloc++;
			}
		}
		if(originallyShort){match=toShortMatchString(match);}
		return true;
	}
	
	public float expectedTipErrors(boolean countUndefined, int maxBases){
		return expectedTipErrors(bases, quality, countUndefined, maxBases);
	}
	
	public float expectedErrorsIncludingMate(boolean countUndefined){
		float a=expectedErrors(countUndefined, length());
		float b=(mate==null ? 0 : mate.expectedErrors(countUndefined, mate.length()));
		return a+b;
	}
	
	public float expectedErrors(boolean countUndefined, int maxBases){
		return expectedErrors(bases, quality, countUndefined, maxBases);
	}
	
	public static float probabilityErrorFree(byte[] bases, byte[] quality, boolean countUndefined, int maxBases){
		if(quality==null){return 0;}
		final int limit=(maxBases<1 ? quality.length : Tools.min(maxBases, quality.length));
		final float[] array=QualityTools.PROB_CORRECT;
		assert(array[0]>0 && array[0]<1);
		float product=1;
		for(int i=0; i<limit; i++){
			byte b=bases[i];
			byte q=quality[i];
			if(AminoAcid.isFullyDefined(b)){
				product*=array[q];
			}else if(countUndefined){
				return 0;
			}
		}
		return product;
	}
	
	public static float expectedErrors(byte[] bases, byte[] quality, boolean countUndefined, int maxBases){
		if(quality==null){return 0;}
		final int limit=(maxBases<1 ? quality.length : Tools.min(maxBases, quality.length));
		final float[] array=QualityTools.PROB_ERROR;
		assert(array[0]>0 && array[0]<1);
		float sum=0;
		for(int i=0; i<limit; i++){
			byte b=bases[i];
			boolean d=AminoAcid.isFullyDefined(b);
			//assert((quality[i]==0)==d) : "Q="+quality[i]+" for base "+(char)b;
			if(d || countUndefined){
				byte q=(d ? quality[i] : 0);
				sum+=array[q];
			}
		}
		return sum;
	}
	
	/** Runs backwards instead of forwards */
	public static float expectedTipErrors(byte[] bases, byte[] quality, boolean countUndefined, int maxBases){
		if(quality==null){return 0;}
		final int limit;
		{
			final int limit0=(maxBases<1 ? quality.length : Tools.min(maxBases, quality.length));
			limit=quality.length-limit0;
		}
		final float[] array=QualityTools.PROB_ERROR;
		assert(array[0]>0 && array[0]<1);
		float sum=0;
		for(int i=quality.length-1; i>=limit; i--){
			byte b=bases[i];
			byte q=quality[i];
			if(AminoAcid.isFullyDefined(b)){
				sum+=array[q];
			}else{
				assert(q==0);
				if(countUndefined){sum+=0.75f;}
			}
		}
		return sum;
	}

	public int estimateErrors() {
		if(quality==null){return 0;}
		assert(match!=null) : this.toText(false);
		
		int count=0;
		for(int ci=0, mi=0; ci<bases.length && mi<match.length; mi++){
			
//			byte b=bases[ci];
			byte q=quality[ci];
			byte m=match[mi];
			if(m=='m' || m=='s' || m=='N'){
				ci++;
			}else if(m=='X' || m=='Y'){
				ci++;
				count++;
			}else if(m=='I'){
				ci++;
			}else if(m=='D'){
				
			}else if(m=='S'){
				ci++;
				if(q<19){
					count++;
				}
			}
			
		}
		return count;
	}
	
	/** {M, S, D, I, N, splice} */
	public int[] countErrors(int minSplice) {
		assert(match!=null) : this.toText(false);
		int m=0;
		int s=0;
		int d=0;
		int i=0;
		int n=0;
		int splice=0;
		
		byte prev=' ';
		int streak=0;
		minSplice=Tools.max(minSplice, 1);
		
		for(int pos=0; pos<match.length; pos++){
			final byte b=match[pos];
			
			if(b==prev){streak++;}else{streak=1;}
			
			if(b=='m'){
				m++;
			}else if(b=='N' || b=='C'){
				n++;
			}else if(b=='X' || b=='Y'){
				i++;
			}else if(b=='I'){
				i++;
			}else if(b=='D'){
				d++;
				if(streak==minSplice){splice++;}
			}else if(b=='S'){
				s++;
			}else{
				if(Tools.isDigit(b) && shortmatch()){
					System.err.println("Warning! Found read in shortmatch form during countErrors():\n"+this); //Usually caused by verbose output.
					if(mate!=null){System.err.println("mate:\n"+mate.id+"\t"+new String(mate.bases));}
					System.err.println("Stack trace: ");
					new Exception().printStackTrace();
					match=toLongMatchString(match);
					setShortMatch(false);
					return countErrors(minSplice);
				}else{
					throw new RuntimeException("\nUnknown symbol "+(char)b+":\n"+new String(match)+"\n"+this+"\nshortmatch="+this.shortmatch());
				}
			}
			
			prev=b;
		}
		
//		assert(i==0) : i+"\n"+this+"\n"+new String(match)+"\n"+Arrays.toString(new int[] {m, s, d, i, n, splice});
		
		return new int[] {m, s, d, i, n, splice};
	}
	
	public static boolean isShortMatchString(byte[] match){
		byte last=' ';
		int streak=0;
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(Tools.isDigit(b)){return true;}
			if(b==last){
				streak++;
				if(streak>3){return true;}
			}else{
				streak=0;
				last=b;
			}
		}
		return false;
	}
	
	public void toShortMatchString(boolean doAssertion){
		if(shortmatch()){
			assert(!doAssertion);
			return;
		}
		match=toShortMatchString(match);
		setShortMatch(true);
	}
	
	public static byte[] toShortMatchString(byte[] match){
		if(match==null){return null;}
		assert(match.length>0);
		ByteBuilder sb=new ByteBuilder(10);
		
		byte prev=match[0];
		int count=1;
		for(int i=1; i<match.length; i++){
			byte m=match[i];
			assert(Tools.isLetter(m) || m==0) : new String(match);
			if(m==0){System.err.println("Warning! Converting empty match string to short form.");}
			if(m==prev){count++;}
			else{
				sb.append(prev);
				if(count>1){sb.append(count);}
//				else if(count==2){sb.append(prev);}
				prev=m;
				count=1;
			}
		}
		sb.append(prev);
		if(count>1){sb.append(count);}
//		else if(count==2){sb.append(prev);}
		
		byte[] r=sb.toBytes();
		return r;
	}
	
	public void toLongMatchString(boolean doAssertion){
		if(!shortmatch()){
			assert(!doAssertion);
			return;
		}
		match=toLongMatchString(match);
		setShortMatch(false);
	}
	
	public static byte[] toLongMatchString(byte[] shortmatch){
		if(shortmatch==null){return null;}
		assert(shortmatch.length>0);
		
		int count=0;
		int current=0;
		for(int i=0; i<shortmatch.length; i++){
			byte m=shortmatch[i];
			if(Tools.isLetter(m)){
				count++;
				count+=(current>0 ? current-1 : 0);
				current=0;
			}else{
				assert(Tools.isDigit(m));
				current=(current*10)+(m-48); //48 == '0'
			}
		}
		count+=(current>0 ? current-1 : 0);
		
		
		byte[] r=new byte[count];
		current=0;
		byte lastLetter='?';
		int j=0;
		for(int i=0; i<shortmatch.length; i++){
			byte m=shortmatch[i];
			if(Tools.isLetter(m)){
				while(current>1){
					r[j]=lastLetter;
					current--;
					j++;
				}
				current=0;
				
				r[j]=m;
				j++;
				lastLetter=m;
			}else{
				assert(Tools.isDigit(m));
				current=(current*10)+(m-48); //48 == '0'
			}
		}
		while(current>1){
			r[j]=lastLetter;
			current--;
			j++;
		}
		
		assert(r[r.length-1]>0);
		return r;
	}
	
	public String parseCustomRname(){
		assert(id.startsWith("SYN")) : "Can't parse header "+id;
		return new CustomHeader(id, pairnum()).rname;
	}
	
	public FloatList fetchVector() {
		FloatList x=(FloatList)obj;
		return x;
	}
	
	public ByteBuilder fetchBB() {
		ByteBuilder x=(ByteBuilder)obj;
		return x;
	}
	
	public Object obj() {return obj;}
	
	public Object setObj(Object x) {
		assert(x!=null);
		Object old=obj;
		assert(old==null) : obj.getClass()+", "+x.getClass(); //Just a warning.  This should probably not happen.
		obj=x;
		return obj;
	}
	
	public Object nullifyObj() {
		Object old=obj;
		obj=null;
		return old;
	}
	
	public Object nullifyObject() {
		Object old=obj;
		obj=null;
		return old;
	}
	
	public Object setObjNull() {
		return nullifyObject();
	}
	
	public TrimRead fetchTrimRead() {
		return (TrimRead)obj;
	}
	
	public Class getObjectClass() {
		return obj==null ? null : obj.getClass();
	}
	
	public Object fetchObject() {
		return obj;
	}

	/** Bases of the read. */
	public byte[] bases;
	
	/** Quality of the read. */
	public byte[] quality;
	
	/** Alignment string.  E.G. mmmmDDDmmm would have 4 matching bases, then a 3-base deletion, then 3 matching bases. */
	public byte[] match;
	
	public int[] gaps;
	
	public String name() {return id;}
	public String id;
	public long numericID;
	public int chrom;
	public int start;
	public int stop;
	
	public int copies=1;

	/** Errors detected */
	public int errors=0;
	
	/** Alignment score from BBMap.  Assumed to max at approx 100*bases.length */
	public int mapScore=0;
	
	public ArrayList<SiteScore> sites;
	public SiteScore originalSite; //Origin site for synthetic reads
	public Object obj=null;//Don't set this to a SamLine; it's for other things.
	public SamLine samline=null;
	public Read mate;
	
	public int flags;
	
	/** -1 if invalid.  TODO: Currently not retained through most processes. */
	private int insert=-1;
	
	/** A random number for deterministic usage.
	 * May decrease speed in multithreaded applications.
	 */
	public double rand=-1;

	public long time(){
		assert(obj!=null && obj.getClass()==Long.class) : obj;
		return ((Long)obj).longValue();
	}
	/** Number of bases in this pair, including the mate if present. */
	public int pairLength(){return length()+mateLength();}
	/** Number of bases in this pair, including the mate if present. */
	public int numPairKmers(int k){return numKmers(k)+numMateKmers(k);}
	/** Number of reads in this pair.  Returns 1 if the read has no mate, and 2 if it does. */
	public int pairCount(){return 1+mateCount();}
	public int pairMappedCount(){return (mapped() ? 1 : 0)+(mate==null || !mate.mapped() ? 0 : 1);}
	public int length(){return bases==null ? 0 : bases.length;}
	public int numKmers(int k){return bases==null ? 0 : Tools.max(0, bases.length-k+1);}
	public int qlength(){return quality==null ? 0 : quality.length;}
	public int mateLength(){return mate==null ? 0 : mate.length();}
	public int numMateKmers(int k){return mate==null ? 0 : mate.numKmers(k);}
	public String mateId(){return mate==null ? null : mate.id;}
	public int mateCount(){return mate==null ? 0 : 1;}
	public boolean mateMapped(){return mate==null ? false : mate.mapped();}
	public boolean eitherMapped(){return mapped() || mateMapped();}
	public long countMateBytes(){return mate==null ? 0 : mate.countBytes();}
	public long countMateFastqBytes(){return mate==null ? 0 : mate.countFastqBytes();}
	/** Number of bytes this read pair uses in memory, approximately */
	public long countPairBytes(){return countBytes()+(mate==null ? 0 : mate.countBytes());}
	
	/** Number of bytes this read uses in memory, approximately */
	public long countBytes(){
		long sum=144; //Approximate per-read overhead
		sum+=(bases==null ? 0 : bases.length+16);
		sum+=(quality==null ? 0 : quality.length+16);
		sum+=(id==null ? 0 : id.length()*2+16);
		sum+=(match==null ? 0 : match.length+16);
		sum+=(samline==null ? 0 : samline.countBytes());
		sum+=(obj==null ? 0 : 32);
		return sum;
	}
	
	/** Number of bytes this read uses in on disk in Fastq format */
	public long countFastqBytes(){
		long sum=6;//4 newlines, +, @
		sum+=(bases==null ? 0 : bases.length);
		sum+=(quality==null ? 0 : quality.length);
		sum+=(id==null ? 0 : id.length());
		return sum;
	}

	public int countLeading(final char base){return countLeft((byte)base);}
	public int countTrailing(final char base){return countRight((byte)base);}
	public int countLeft(final char base){return countLeft((byte)base);}
	public int countRight(final char base){return countRight((byte)base);}
	
	public int countLeft(final byte base){
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			if(b!=base){return i;}
		}
		return bases.length;
	}
	
	public int countRight(final byte base){
		for(int i=bases.length-1; i>=0; i--){
			final byte b=bases[i];
			if(b!=base){return bases.length-i-1;}
		}
		return bases.length;
	}
	
	public boolean untrim(){
		if(obj==null || obj.getClass()!=TrimRead.class){return false;}
		((TrimRead)obj).untrim();
		obj=null;
		return true;
	}
	
	public int trailingLowerCase(){
		for(int i=bases.length-1; i>=0;){
			if(Tools.isLowerCase(bases[i])){
				i--;
			}else{
				return bases.length-i-1;
			}
		}
		return bases.length;
	}
	public int leadingLowerCase(){
		for(int i=0; i<bases.length; i++){
			if(!Tools.isLowerCase(bases[i])){return i;}
		}
		return bases.length;
	}

	public char strandChar(){return Shared.strandCodes2[strand()];}
	public byte strand(){return (byte)(flags&1);}
	public boolean mapped(){return (flags&MAPPEDMASK)==MAPPEDMASK;}
	public boolean paired(){return (flags&PAIREDMASK)==PAIREDMASK;}
	public boolean synthetic(){return (flags&SYNTHMASK)==SYNTHMASK;}
	public boolean ambiguous(){return (flags&AMBIMASK)==AMBIMASK;}
	public boolean perfect(){return (flags&PERFECTMASK)==PERFECTMASK;}
//	public boolean semiperfect(){return perfect() ? true : list!=null && list.size()>0 ? list.get(0).semiperfect : false;} //TODO: This is a hack.  Add a semiperfect flag.
	public boolean rescued(){return (flags&RESCUEDMASK)==RESCUEDMASK;}
	public boolean discarded(){return (flags&DISCARDMASK)==DISCARDMASK;}
	public boolean invalid(){return (flags&INVALIDMASK)==INVALIDMASK;}
	public boolean swapped(){return (flags&SWAPMASK)==SWAPMASK;}
	public boolean shortmatch(){return (flags&SHORTMATCHMASK)==SHORTMATCHMASK;}
	public boolean insertvalid(){return (flags&INSERTMASK)==INSERTMASK;}
	public boolean hasAdapter(){return (flags&ADAPTERMASK)==ADAPTERMASK;}
	public boolean secondary(){return (flags&SECONDARYMASK)==SECONDARYMASK;}
	public boolean aminoacid(){return (flags&AAMASK)==AAMASK;}
	public boolean amino(){return (flags&AAMASK)==AAMASK;}
	public boolean junk(){return (flags&JUNKMASK)==JUNKMASK;}
	public boolean validated(){return (flags&VALIDATEDMASK)==VALIDATEDMASK;}
	public boolean tested(){return (flags&TESTEDMASK)==TESTEDMASK;}
	public boolean invertedRepeat(){return (flags&IRMASK)==IRMASK;}
	public boolean trimmed(){return (flags&TRIMMEDMASK)==TRIMMEDMASK;}
	
	/** For paired ends: 0 for read1, 1 for read2 */
	public int pairnum(){return (flags&PAIRNUMMASK)>>PAIRNUMSHIFT;}
	public boolean valid(){return !invalid();}

	public boolean getFlag(int mask){return (flags&mask)==mask;}
	public int flagToNumber(int mask){return (flags&mask)==mask ? 1 : 0;}
	
	public void setFlag(int mask, boolean b){
		flags=(flags&~mask);
		if(b){flags|=mask;}
	}
	
	public void setStrand(int b){
		assert(b==1 || b==0);
		flags=(flags&(~1))|b;
	}
	
	/** For paired ends: 0 for read1, 1 for read2 */
	public void setPairnum(int b){
//		System.err.println("Setting pairnum to "+b+" for "+id);
//		assert(!id.equals("2_chr1_0_1853883_1853982_1845883_ecoli_K12") || b==1);
		assert(b==1 || b==0);
		flags=(flags&(~PAIRNUMMASK))|(b<<PAIRNUMSHIFT);
//		assert(pairnum()==b);
	}
	
	public void setPaired(boolean b){
		flags=(flags&~PAIREDMASK);
		if(b){flags|=PAIREDMASK;}
	}
	
	public void setSynthetic(boolean b){
		flags=(flags&~SYNTHMASK);
		if(b){flags|=SYNTHMASK;}
	}
	
	public void setAmbiguous(boolean b){
		flags=(flags&~AMBIMASK);
		if(b){flags|=AMBIMASK;}
	}
	
	public boolean setPerfectFlag(int maxScore){
		final SiteScore ss=topSite();
		if(ss==null){
			setPerfect(false);
		}else{
			assert(ss.slowScore<=maxScore) : maxScore+", "+ss.slowScore+", "+ss.toText();
			
			if(ss.slowScore==maxScore || ss.perfect){
				assert(testMatchPerfection(true)) : "\n"+ss+"\n"+maxScore+"\n"+this+"\n"+mate+"\n";
				setPerfect(true);
			}else{
				boolean flag=testMatchPerfection(false);
				setPerfect(flag);
				assert(flag || !ss.perfect) : "flag="+flag+", ss.perfect="+ss.perfect+"\nmatch="+new String(match)+"\n"+this.toText(false);
				assert(!flag || ss.slowScore>=maxScore) : "\n"+ss+"\n"+maxScore+"\n"+this+"\n"+mate+"\n";
			}
		}
		return perfect();
	}
	
	private boolean testMatchPerfection(boolean returnIfNoMatch){
		if(match==null){return returnIfNoMatch;}
		boolean flag=(match.length==bases.length);
		if(shortmatch()){
			flag=(match.length==0 || match[0]=='m');
			for(int i=0; i<match.length && flag; i++){flag=(match[i]=='m' || Tools.isDigit(match[i]));}
		}else{
			for(int i=0; i<match.length && flag; i++){flag=(match[i]=='m');}
		}
		for(int i=0; i<bases.length && flag; i++){flag=(bases[i]!='N');}
		return flag;
	}

	/**
	 * @return GC fraction
	 */
	public float gc() {
		if(bases==null || bases.length<1){return 0;}
		int at=0, gc=0;
		for(byte b : bases){
			int x=AminoAcid.baseToNumber[b];
			if(x>-1){
				if(x==0 || x==3){at++;}
				else{gc++;}
			}
		}
		if(gc<1){return 0;}
		return gc*1f/(at+gc);
	}

	/**
	 * @param swapFrom
	 * @param swapTo
	 * @return number of swaps
	 */
	public int swapBase(byte swapFrom, byte swapTo) {
		if(bases==null){return 0;}
		int swaps=0;
		for(int i=0; i<bases.length; i++){
			if(bases[i]==swapFrom){
				bases[i]=swapTo;
				swaps++;
			}
		}
		return swaps;
	}

	/**
	 * @param remap Table of new values
	 */
	public void remap(byte[] remap) {
		if(bases==null){return;}
		for(int i=0; i<bases.length; i++){
			bases[i]=remap[bases[i]];
		}
	}

	/**
	 * @param remap Table of new values
	 */
	public int remapAndCount(byte[] remap) {
		if(bases==null){return 0;}
		int swaps=0;
		for(int i=0; i<bases.length; i++){
			byte a=bases[i];
			byte b=remap[a];
			if(a!=b){
				bases[i]=b;
				swaps++;
			}
		}
		return swaps;
	}
	
	public int convertUndefinedTo(byte b){
		if(bases==null){return 0;}
		int changed=0;
		for(int i=0; i<bases.length; i++){
			if(b<0 || AminoAcid.baseToNumberACGTN[bases[i]]<0){
				changed++;
				bases[i]=b;
				if(quality!=null){quality[i]=0;}
			}
		}
		return changed;
	}
	
	public void swapBasesWithMate(){
		if(mate==null){
			assert(false);
			return;
		}
		byte[] temp=bases;
		bases=mate.bases;
		mate.bases=temp;
		temp=quality;
		quality=mate.quality;
		mate.quality=temp;
	}
	
	public int insert(){
		return insertvalid() ? insert : -1;
	}
	
	public int insertSizeMapped(boolean ignoreStrand){
		return insertSizeMapped(this, mate, ignoreStrand);
	}
	
	public static int insertSizeMapped(Read r1, Read r2, boolean ignoreStrand){
//		assert(false) : ignoreStrand+", "+(r2==null)+", "+(r1.mapped())+", "+(r2.mapped())+", "+(r1.strand()==r2.strand())+", "+r1.strand()+", "+r2.strand();
		if(ignoreStrand || r2==null || !r1.mapped() || !r2.mapped() || r1.strand()==r2.strand()){return insertSizeMapped_Unstranded(r1, r2);}
		return insertSizeMapped_PlusLeft(r1, r2);
	}
	
	/** TODO: This is not correct when the insert is shorter than a read's bases with same-strand reads */
	public static int insertSizeMapped_PlusLeft(Read r1, Read r2){
		if(r1.strand()>r2.strand()){return insertSizeMapped_PlusLeft(r2, r1);}
		if(r1.strand()==r2.strand() || r1.start>r2.stop){return insertSizeMapped_Unstranded(r2, r1);} //So r1 is always on the left.
//		if(!mapped() || !mate.mapped()){return 0;}
		if(r1.chrom!=r2.chrom){return 0;}
		if(r1.start==r1.stop || r2.start==r2.stop){return 0;} //???
		
		int a=r1.length();
		int b=r2.length();
		int mid=r2.start-r1.stop-1;
		if(-mid>=a+b){return insertSizeMapped_Unstranded(r1, r2);} //Not properly oriented; plus read is to the right of minus read
		return mid+a+b;
	}
	
	public static int insertSizeMapped_Unstranded(Read r1, Read r2){
		if(r2==null){return r1.start==r1.stop ? 0 : r1.stop-r1.start+1;}
		
		if(r1.start>r2.start){return insertSizeMapped_Unstranded(r2, r1);} //So r1 is always on the left side.
		
//		if(!mapped() || !mate.mapped()){return 0;}
		if(r1.start==r1.stop || r2.start==r2.stop){return 0;} //???
		
		if(r1.chrom!=r2.chrom){return 0;}
		int a=r1.length();
		int b=r2.length();
		if(false && Tools.overlap(r1.start, r1.stop, r2.start, r2.stop)){
			//This does not handle very short inserts
			return Tools.max(r1.stop, r2.stop)-Tools.min(r1.start, r2.start)+1;
			
		}else{
			if(r1.start<r2.start){
				int mid=r2.start-r1.stop-1;
//				assert(false) : mid+", "+a+", "+b;
//				if(-mid>a && -mid>b){return Tools.min(a, b);} //Strange situation, no way to guess insert size
				if(-mid>=a+b){return 0;} //Strange situation, no way to guess insert size
				return mid+a+b;
			}else{
				assert(r1.start==r2.start);
				return Tools.min(a, b);
			}
		}
	}
	
	public int insertSizeOriginalSite(){
		if(mate==null){
//			System.err.println("A: "+(originalSite==null ? "null" : (originalSite.stop-originalSite.start+1)));
			return (originalSite==null ? -1 : originalSite.stop-originalSite.start+1);
		}
		
		final SiteScore ssa=originalSite, ssb=mate.originalSite;
		final int x;
		if(ssa==null || ssb==null){
//			System.err.println("B: 0");
			x=0;
		}else{
			x=insertSize(ssa, ssb, bases.length, mate.length());
		}
		
		assert(pairnum()>=mate.pairnum() || x==mate.insertSizeOriginalSite());
		return x;
	}
	
	public static int insertSize(SiteScore ssa, SiteScore ssb, int lena, int lenb){
		return insertSize(ssa.chrom, ssb.chrom, ssa.start, ssb.start, ssa.stop, ssb.stop, lena, lenb);
	}
	
	public static int insertSize(int chroma, int chromb, int starta, int startb, int stopa, int stopb, int lena, int lenb){
		
		final int x;

		//		if(mate==null || ){return bases==null ? 0 : bases.length;}
		if(chroma!=chromb){x=0;}
		else{

			if(Tools.overlap(starta, stopa, startb, stopb)){
				x=Tools.max(stopa, stopb)-Tools.min(starta, startb)+1;
//				System.err.println("C: "+x);
			}else{
				if(starta<=startb){
					int mid=startb-stopa-1;
					//				assert(false) : mid+", "+a+", "+b;
					x=mid+lena+lenb;
//					System.err.println("D: "+x);
				}else{
					int mid=starta-stopb-1;
					//				assert(false) : mid+", "+a+", "+b;
					x=mid+lena+lenb;
//					System.err.println("E: "+x);
				}
			}
		}
		return x;
	}
	
	public Read subRead(int from, int to){
		Read r=this.clone();
		r.bases=KillSwitch.copyOfRange(bases, from, to);
		r.quality=(quality==null ? null : KillSwitch.copyOfRange(quality, from, to));
		r.mate=null;
//		assert(Tools.indexOf(r.bases, (byte)'J')<0);
		return r;
	}
	
	public Read joinRead(){
		if(insert<1 || mate==null || !insertvalid()){return this;}
		assert(insert>9 || bases.length<20) : "Perhaps old read format is being used?  This appears to be a quality value, not an insert.\n"+this+"\n\n"+mate+"\n";
		return joinRead(this, mate, insert);
	}
	
	public Read joinRead(int x){
		if(x<1 || mate==null){return this;}
		//assert(x>9 || bases.length<20) : "Perhaps old read format is being used?  This appears to be a quality value, not an insert.\n"+this+"\n\n"+mate+"\n";
		return joinRead(this, mate, x);
	}
	
	public static Read joinRead(Read a, Read b, int insert){
		assert(a!=null && b!=null && insert>0);
		final int lengthSum=a.length()+b.length();
		final int overlap=Tools.min(insert, lengthSum-insert);
		
//		System.err.println(insert);
		final byte[] bases=new byte[insert], abases=a.bases, bbases=b.bases;
		final byte[] aquals=a.quality, bquals=b.quality;
		final byte[] quals=(aquals==null || bquals==null ? null : new byte[insert]);
		assert(aquals==null || (aquals.length==abases.length && bquals.length==bbases.length));
		
		int mismatches=0;
		
		int start, stop;
		
		if(overlap<=0){//Simple join in which there is no overlap
			int lim=insert-b.length();
			if(quals==null){
				for(int i=0; i<a.length(); i++){
					bases[i]=abases[i];
				}
				for(int i=a.length(); i<lim; i++){
					bases[i]='N';
				}
				for(int i=0; i<b.length(); i++){
					bases[i+lim]=bbases[i];
				}
			}else{
				for(int i=0; i<a.length(); i++){
					bases[i]=abases[i];
					quals[i]=aquals[i];
				}
				for(int i=a.length(); i<lim; i++){
					bases[i]='N';
					quals[i]=0;
				}
				for(int i=0; i<b.length(); i++){
					bases[i+lim]=bbases[i];
					quals[i+lim]=bquals[i];
				}
			}
			
			start=Tools.min(a.start, b.start);
//			stop=start+insert-1;
			stop=Tools.max(a.stop, b.stop);
			
//		}else if(insert>=a.length() && insert>=b.length()){ //Overlapped join, proper orientation
//			final int lim1=a.length()-overlap;
//			final int lim2=a.length();
//			for(int i=0; i<lim1; i++){
//				bases[i]=abases[i];
//				quals[i]=aquals[i];
//			}
//			for(int i=lim1, j=0; i<lim2; i++, j++){
//				assert(false) : "TODO";
//				bases[i]='N';
//				quals[i]=0;
//			}
//			for(int i=lim2, j=overlap; i<bases.length; i++, j++){
//				bases[i]=bbases[j];
//				quals[i]=bquals[j];
//			}
		}else{ //reads go off ends of molecule.
			if(quals==null){
				for(int i=0; i<a.length() && i<bases.length; i++){
					bases[i]=abases[i];
				}
				for(int i=bases.length-1, j=b.length()-1; i>=0 && j>=0; i--, j--){
					byte ca=bases[i], cb=bbases[j];
					if(ca==0 || ca=='N'){
						bases[i]=cb;
					}else if(ca==cb){
					}else{
						bases[i]=(ca>=cb ? ca : cb);
						if(ca!='N' && cb!='N'){mismatches++;}
					}
				}
			}else{
				for(int i=0; i<a.length() && i<bases.length; i++){
					bases[i]=abases[i];
					quals[i]=aquals[i];
				}
				for(int i=bases.length-1, j=b.length()-1; i>=0 && j>=0; i--, j--){
					byte ca=bases[i], cb=bbases[j];
					byte qa=quals[i], qb=bquals[j];
					if(ca==0 || ca=='N'){
						bases[i]=cb;
						quals[i]=qb;
					}else if(cb==0 || cb=='N'){
						//do nothing
					}else if(ca==cb){
						quals[i]=(byte)Tools.min((Tools.max(qa, qb)+Tools.min(qa, qb)/4), MAX_MERGE_QUALITY);
					}else{
						bases[i]=(qa>qb ? ca : qa<qb ? cb : (byte)'N');
						quals[i]=(byte)(Tools.max(qa, qb)-Tools.min(qa, qb));
						if(ca!='N' && cb!='N'){mismatches++;}
					}
				}
			}
			
			if(a.strand()==0){
				start=a.start;
//				stop=start+insert-1;
				stop=b.stop;
			}else{
				stop=a.stop;
//				start=stop-insert+1;
				start=b.start;
			}
			if(start>stop){
				start=Tools.min(a.start, b.start);
				stop=Tools.max(a.stop, b.stop);
			}
		}
//		assert(mismatches>=countMismatches(a, b, insert, 999));
//		System.err.println(mismatches);
		if(a.chrom==0 || start==stop || (!a.mapped() && !a.synthetic())){start=stop=a.chrom=0;}
		
//		System.err.println(bases.length+", "+start+", "+stop);
		
		Read r=new Read(bases, null, a.id, a.numericID, a.flags, a.chrom, start, stop);
		r.quality=quals; //This prevents quality from getting capped.
		if(a.chrom==0 || start==stop || (!a.mapped() && !a.synthetic())){r.setMapped(true);}
		r.setInsert(insert);
		r.setPaired(false);
		r.copies=a.copies;
		r.mapScore=a.mapScore+b.mapScore;
		if(overlap<=0){
			r.mapScore=a.mapScore+b.mapScore;
			r.errors=a.errors+b.errors;
			//TODO r.gaps=?
		}else{//Hard to calculate
			r.mapScore=(int)((a.mapScore*(long)a.length()+b.mapScore*(long)b.length())/insert);
			r.errors=a.errors;
		}
		
		
		assert(r.insertvalid()) : "\n\n"+a.toText(false)+"\n\n"+b.toText(false)+"\n\n"+r.toText(false)+"\n\n";
		assert(r.insert()==r.length()) : r.insert()+"\n\n"+a.toText(false)+"\n\n"+b.toText(false)+"\n\n"+r.toText(false)+"\n\n";
//		assert(false) : "\n\n"+a.toText(false)+"\n\n"+b.toText(false)+"\n\n"+r.toText(false)+"\n\n";
		
		//TODO: Triggered by BBMerge in useratio mode for some reason.
//		assert(Shared.anomaly || (a.insertSizeMapped(false)>0 == r.insertSizeMapped(false)>0)) :
//			"\n"+r.length()+"\n"+r.insert()+"\n"+r.insertSizeMapped(false)+"\n"+a.insert()+"\n"+a.insertSizeMapped(false)+
//			"\n\n"+a.toText(false)+"\n\n"+b.toText(false)+"\n\n"+r.toText(false)+"\n\n";
		
		return r;
	}

	/**
	 * @param minlen
	 * @param maxlen
	 * @return A list of read fragments
	 */
	public ArrayList<Read> split(int minlen, int maxlen) {
		int len=bases==null ? 0 : bases.length;
		if(len<minlen){return null;}
		int parts=(len+maxlen-1)/maxlen;
		ArrayList<Read> subreads=new ArrayList<Read>(parts);
		if(len<=maxlen){
			subreads.add(this);
		}else{
			float ideal=Tools.max(minlen, len/(float)parts);
			int actual=(int)ideal;
			assert(false) : "TODO"; //Some assertion goes here, I forget what
			for(int i=0; i<parts; i++){
				int a=i*actual;
				int b=(i+1)*actual;
				if(b>bases.length){b=bases.length;}
//				if(b-a<)
				byte[] subbases=KillSwitch.copyOfRange(bases, a, b);
				byte[] subquals=(quality==null ? null : KillSwitch.copyOfRange(quality, a, b+1));
				Read r=new Read(subbases, subquals, id+"_"+i, numericID, flags);
				subreads.add(r);
			}
		}
		return subreads;
	}
	
	/** Generate and return an array of canonical kmers for this read */
	public long[] toKmers(final int k, final int gap, long[] kmers, boolean makeCanonical, Kmer longkmer) {
		if(gap>0){throw new RuntimeException("Gapped reads: TODO");}
		if(k>31){return toLongKmers(k, kmers, makeCanonical, longkmer);}
		if(bases==null || bases.length<k+gap){return null;}
		
		final int arraylen=bases.length-k+1;
		if(kmers==null || kmers.length!=arraylen){kmers=new long[arraylen];}
		Arrays.fill(kmers, -1);
		
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, rkmer=0;
		int len=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(x<0){len=0; rkmer=0;}else{len++;}
			if(len>=k){
				kmers[i-k+1]=makeCanonical ? Tools.max(kmer, rkmer) : kmer;
			}
		}
		return kmers;
	}
	
//	/** Generate and return an array of canonical kmers for this read */
//	public long[] toKmers(final int k, final int gap, long[] kmers, boolean makeCanonical, Kmer longkmer) {
//		if(gap>0){throw new RuntimeException("Gapped reads: TODO");}
//		if(k>31){return toLongKmers(k, kmers, makeCanonical, longkmer);}
//		if(bases==null || bases.length<k+gap){return null;}
//
//		final int kbits=2*k;
//		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
//
//		int len=0;
//		long kmer=0;
//		final int arraylen=bases.length-k+1;
//		if(kmers==null || kmers.length!=arraylen){kmers=new long[arraylen];}
//		Arrays.fill(kmers, -1);
//
//		for(int i=0; i<bases.length; i++){
//			byte b=bases[i];
//			int x=AminoAcid.baseToNumber[b];
//			if(x<0){
//				len=0;
//				kmer=0;
//			}else{
//				kmer=((kmer<<2)|x)&mask;
//				len++;
//
//				if(len>=k){
//					kmers[i-k+1]=kmer;
//				}
//			}
//		}
//
////		System.out.println(new String(bases));
////		System.out.println(Arrays.toString(kmers));
//
//		if(makeCanonical){
//			this.reverseComplement();
//			len=0;
//			kmer=0;
//			for(int i=0, j=bases.length-1; i<bases.length; i++, j--){
//				byte b=bases[i];
//				int x=AminoAcid.baseToNumber[b];
//				if(x<0){
//					len=0;
//					kmer=0;
//				}else{
//					kmer=((kmer<<2)|x)&mask;
//					len++;
//
//					if(len>=k){
//						assert(kmer==AminoAcid.reverseComplementBinaryFast(kmers[j], k));
//						kmers[j]=Tools.max(kmers[j], kmer);
//					}
//				}
//			}
//			this.reverseComplement();
//
////			System.out.println(Arrays.toString(kmers));
//		}
//
//
//		return kmers;
//	}
	
	/** Generate and return an array of canonical kmers for this read */
	public long[] toLongKmers(final int k, long[] kmers, boolean makeCanonical, Kmer kmer) {
		assert(k>31) : k;
		assert(makeCanonical);
		if(bases==null || bases.length<k){return null;}
		kmer.clear();
		
		final int arraylen=bases.length-k+1;
		if(kmers==null || kmers.length!=arraylen){kmers=new long[arraylen];}
		Arrays.fill(kmers, -1);
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			kmer.addRight(b);
			if(!AminoAcid.isFullyDefined(b)){kmer.clear();}
			if(kmer.len>=k){
				kmers[i-k+1]=kmer.xor();
			}
		}
		
		return kmers;
	}
	
//	/** Generate and return an array of canonical kmers for this read */
//	public long[] toLongKmers(final int k, long[] kmers, boolean makeCanonical, Kmer longkmer) {
//		assert(k>31) : k;
//		if(bases==null || bases.length<k){return null;}
//
//		final int kbits=2*k;
//		final long mask=Long.MAX_VALUE;
//
//		int len=0;
//		long kmer=0;
//		final int arraylen=bases.length-k+1;
//		if(kmers==null || kmers.length!=arraylen){kmers=new long[arraylen];}
//		Arrays.fill(kmers, -1);
//
//
//		final int tailshift=k%32;
//		final int tailshiftbits=tailshift*2;
//
//		for(int i=0; i<bases.length; i++){
//			byte b=bases[i];
//			int x=AminoAcid.baseToNumber[b];
//			if(x<0){
//				len=0;
//				kmer=0;
//			}else{
//				kmer=Long.rotateLeft(kmer, 2);
//				kmer=kmer^x;
//				len++;
//
//				if(len>=k){
//					long x2=AminoAcid.baseToNumber[bases[i-k]];
//					kmer=kmer^(x2<<tailshiftbits);
//					kmers[i-k+1]=kmer;
//				}
//			}
//		}
//		if(makeCanonical){
//			this.reverseComplement();
//			len=0;
//			kmer=0;
//			for(int i=0, j=bases.length-1; i<bases.length; i++, j--){
//				byte b=bases[i];
//				int x=AminoAcid.baseToNumber[b];
//				if(x<0){
//					len=0;
//					kmer=0;
//				}else{
//					kmer=Long.rotateLeft(kmer, 2);
//					kmer=kmer^x;
//					len++;
//
//					if(len>=k){
//						long x2=AminoAcid.baseToNumber[bases[i-k]];
//						kmer=kmer^(x2<<tailshiftbits);
//						kmers[j]=mask&(Tools.max(kmers[j], kmer));
//					}
//				}
//			}
//			this.reverseComplement();
//		}else{
//			assert(false) : "Long kmers should be made canonical here because they cannot be canonicized later.";
//		}
//
//		return kmers;
//	}
	
	public static final boolean CHECKSITES(Read r, byte[] basesM){
		return CHECKSITES(r.sites, r.bases, basesM, r.numericID, true);
	}
	
	public static final boolean CHECKSITES(Read r, byte[] basesM, boolean verifySorted){
		return CHECKSITES(r.sites, r.bases, basesM, r.numericID, verifySorted);
	}
	
	public static final boolean CHECKSITES(ArrayList<SiteScore> list, byte[] basesP, byte[] basesM, long id){
		return CHECKSITES(list, basesP, basesM, id, true);
	}
	
	public static final boolean CHECKSITES(ArrayList<SiteScore> list, byte[] basesP, byte[] basesM, long id, boolean verifySorted){
		return true; //Temporarily disabled
//		if(list==null || list.isEmpty()){return true;}
//		SiteScore prev=null;
//		for(int i=0; i<list.size(); i++){
//			SiteScore ss=list.get(i);
//			if(ss.strand==Gene.MINUS && basesM==null && basesP!=null){basesM=AminoAcid.reverseComplementBases(basesP);}
//			byte[] bases=(ss.strand==Gene.PLUS ? basesP : basesM);
//			if(verbose){System.err.println("Checking site "+i+": "+ss);}
//			boolean b=CHECKSITE(ss, bases, id);
//			assert(b) : id+"\n"+new String(basesP)+"\n"+ss+"\n";
//			if(verbose){System.err.println("Checked site "+i+" = "+ss+"\nss.p="+ss.perfect+", ss.sp="+ss.semiperfect);}
//			if(!b){
////				System.err.println("Error at SiteScore "+i+": ss.p="+ss.perfect+", ss.sp="+ss.semiperfect);
//				return false;
//			}
//			if(verifySorted && prev!=null && ss.score>prev.score){
//				if(verbose){System.err.println("verifySorted failed.");}
//				return false;
//			}
//			prev=ss;
//		}
//		return true;
	}
	
	/** Makes sure 'bases' is for correct strand. */
	public static final boolean CHECKORDER(ArrayList<SiteScore> list){
		if(list==null || list.size()<2){return true;}
		SiteScore prev=list.get(0);
		for(int i=0; i<list.size(); i++){
			SiteScore ss=list.get(i);
			if(ss.score>prev.score){return false;}
			prev=ss;
		}
		return true;
	}
	
	/** Makes sure 'bases' is for correct strand. */
	public static final boolean CHECKSITE(SiteScore ss, byte[] basesP, byte[] basesM, long id){
		return CHECKSITE(ss, ss.plus() ? basesP : basesM, id);
	}
	
	/** Make sure 'bases' is for correct strand! */
	public static final boolean CHECKSITE(SiteScore ss, byte[] bases, long id){
		return true; //Temporarily disabled
//		if(ss==null){return true;}
////		System.err.println("Checking site "+ss+"\nss.p="+ss.perfect+", ss.sp="+ss.semiperfect+", bases="+new String(bases));
//		if(ss.perfect){assert(ss.semiperfect) : ss+"\n"+new String(bases);}
//		if(ss.gaps!=null){
//			if(ss.gaps[0]!=ss.start || ss.gaps[ss.gaps.length-1]!=ss.stop){return false;}
////			assert(ss.gaps[0]==ss.start && ss.gaps[ss.gaps.length-1]==ss.stop);
//		}
//
//		if(!(ss.pairedScore<1 || (ss.slowScore<=0 && ss.pairedScore>ss.quickScore ) || ss.pairedScore>ss.slowScore)){
//			System.err.println("Site paired score violation: "+ss.quickScore+", "+ss.slowScore+", "+ss.pairedScore);
//			return false;
//		}
//
//		final boolean xy=ss.matchContainsXY();
//		if(bases!=null){
//
//			final boolean p0=ss.perfect;
//			final boolean sp0=ss.semiperfect;
//			final boolean p1=ss.isPerfect(bases);
//			final boolean sp1=(p1 ? true : ss.isSemiPerfect(bases));
//
//			assert(p0==p1 || (xy && p1)) : p0+"->"+p1+", "+sp0+"->"+sp1+", "+ss.isSemiPerfect(bases)+
//				"\nnumericID="+id+"\n"+new String(bases)+"\n\n"+Data.getChromosome(ss.chrom).getString(ss.start, ss.stop)+"\n\n"+ss+"\n\n";
//			assert(sp0==sp1 || (xy && sp1)) : p0+"->"+p1+", "+sp0+"->"+sp1+", "+ss.isSemiPerfect(bases)+
//				"\nnumericID="+id+"\n"+new String(bases)+"\n\n"+Data.getChromosome(ss.chrom).getString(ss.start, ss.stop)+"\n\n"+ss+"\n\n";
//
////			ss.setPerfect(bases, false);
//
//			assert(p0==ss.perfect) :
//				p0+"->"+ss.perfect+", "+sp0+"->"+ss.semiperfect+", "+ss.isSemiPerfect(bases)+"\nnumericID="+id+"\n\n"+new String(bases)+"\n\n"+
//				Data.getChromosome(ss.chrom).getString(ss.start, ss.stop)+"\n"+ss+"\n\n";
//			assert(sp0==ss.semiperfect) :
//				p0+"->"+ss.perfect+", "+sp0+"->"+ss.semiperfect+", "+ss.isSemiPerfect(bases)+"\nnumericID="+id+"\n\n"+new String(bases)+"\n\n"+
//				Data.getChromosome(ss.chrom).getString(ss.start, ss.stop)+"\n"+ss+"\n\n";
//			if(ss.perfect){assert(ss.semiperfect);}
//		}
//		if(ss.match!=null && ss.matchLength()!=ss.mappedLength()){
//			if(verbose){System.err.println("Returning false because matchLength!=mappedLength:\n"+ss.matchLength()+", "+ss.mappedLength()+"\n"+ss);}
//			return false;
//		}
//		return true;
	}
	
	public void setPerfect(boolean b){
		flags=(flags&~PERFECTMASK);
		if(b){flags|=PERFECTMASK;}
	}
	
	public void setRescued(boolean b){
		flags=(flags&~RESCUEDMASK);
		if(b){flags|=RESCUEDMASK;}
	}
	
	public void setMapped(boolean b){
		flags=(flags&~MAPPEDMASK);
		if(b){flags|=MAPPEDMASK;}
	}
	
	public void setPairDiscarded(boolean b){
		setDiscarded(b);
		if(mate!=null) {mate.setDiscarded(b);}
	}
	public void setDiscarded(boolean b){
		flags=(flags&~DISCARDMASK);
		if(b){flags|=DISCARDMASK;}
	}
	
	public void setInvalid(boolean b){
		flags=(flags&~INVALIDMASK);
		if(b){flags|=INVALIDMASK;}
	}
	
	public void setSwapped(boolean b){
		flags=(flags&~SWAPMASK);
		if(b){flags|=SWAPMASK;}
	}
	
	public void setShortMatch(boolean b){
		flags=(flags&~SHORTMATCHMASK);
		if(b){flags|=SHORTMATCHMASK;}
	}
	
	public void setInsertValid(boolean b){
		flags=(flags&~INSERTMASK);
		if(b){flags|=INSERTMASK;}
	}
	
	public void setHasAdapter(boolean b){
		flags=(flags&~ADAPTERMASK);
		if(b){flags|=ADAPTERMASK;}
	}
	
	public void setSecondary(boolean b){
		flags=(flags&~SECONDARYMASK);
		if(b){flags|=SECONDARYMASK;}
	}
	
	public void setAminoAcid(boolean b){
		flags=(flags&~AAMASK);
		if(b){flags|=AAMASK;}
	}
	
	public void setJunk(boolean b){
		flags=(flags&~JUNKMASK);
		if(b){flags|=JUNKMASK;}
	}
	
	public void setValidated(boolean b){
		flags=(flags&~VALIDATEDMASK);
		if(b){flags|=VALIDATEDMASK;}
	}
	
	public void setTested(boolean b){
		flags=(flags&~TESTEDMASK);
		if(b){flags|=TESTEDMASK;}
	}
	
	public void setInvertedRepeat(boolean b){
		flags=(flags&~IRMASK);
		if(b){flags|=IRMASK;}
	}
	
	public void setTrimmed(boolean b){
		flags=(flags&~TRIMMEDMASK);
		if(b){flags|=TRIMMEDMASK;}
	}
	
	public void setInsert(int x){
		if(x<1){x=-1;}
//		assert(x==-1 || x>9 || length()<20) : x+", "+length(); //Invalid assertion for synthetic reads.
		insert=x;
		setInsertValid(x>0);
		if(mate!=null){
			mate.insert=x;
			mate.setInsertValid(x>0);
		}
	}

	private static int[] makeMaskArray(int max) {
		int[] r=new int[max+1];
		for(int i=0; i<r.length; i++){r[i]=(1<<i);}
		return r;
	}
	


	public static byte[] getFakeQuality(int len){
		if(len>=QUALCACHE.length){
			byte[] r=KillSwitch.allocByte1D(len);
			Arrays.fill(r, (byte)30);
			return r;
		}
		if(QUALCACHE[len]==null){
			synchronized(QUALCACHE){
				if(QUALCACHE[len]==null){
					QUALCACHE[len]=KillSwitch.allocByte1D(len);
					Arrays.fill(QUALCACHE[len], (byte)30);
				}
			}
		}
		return QUALCACHE[len];
	}
	
	public byte[] getScaffoldName(boolean requireSingleScaffold){
		byte[] name=null;
		if(mapped()){
			if(!requireSingleScaffold || Data.isSingleScaffold(chrom, start, stop)){
				int idx=Data.scaffoldIndex(chrom, (start+stop)/2);
				name=Data.scaffoldNames[chrom][idx];
//				int scaflen=Data.scaffoldLengths[chrom][idx];
//				a1=Data.scaffoldRelativeLoc(chrom, start, idx);
//				b1=a1-start1+stop1;
			}
		}
		return name;
	}
	
	public void bisulfite(boolean AtoG, boolean CtoT, boolean GtoA, boolean TtoC){
		for(int i=0; i<bases.length; i++){
			final int x=AminoAcid.baseToNumber[bases[i]];
			if(x==0 && AtoG){bases[i]='G';}
			else if(x==1 && CtoT){bases[i]='T';}
			else if(x==2 && GtoA){bases[i]='A';}
			else if(x==3 && TtoC){bases[i]='C';}
		}
	}
	
	public Read copy(){
		Read r=clone();
		r.bases=(r.bases==null ? null : r.bases.clone());
		r.quality=(r.quality==null ? null : r.quality.clone());
		r.match=(r.match==null ? null : r.match.clone());
		r.gaps=(r.gaps==null ? null : r.gaps.clone());
		r.originalSite=(r.originalSite==null ? null : r.originalSite.clone());
		r.sites=(ArrayList<SiteScore>) (r.sites==null ? null : r.sites.clone());
		r.mate=null;
		
		if(r.sites!=null){
			for(int i=0; i<r.sites.size(); i++){
				r.sites.set(i, r.sites.get(i).clone());
			}
		}
		return r;
	}
	
	@Override
	public Read clone(){
		try {
			return (Read) super.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		throw new RuntimeException();
	}
	
	/**
	 * @return This protein in canonical nucleotide space.
	 */
	public Read aminoToNucleic() {
		assert(aminoacid()) : "This read is not flagged as an amino acid sequence.";
		Read r=this.clone();
		r.setAminoAcid(false);
		r.bases=AminoAcid.toNTs(r.bases);
		if(quality!=null){
			byte[] ntquals=new byte[r.quality.length*3];
			for(int i=0; i<quality.length; i++){
				byte q=quality[i];
				byte q2=(byte)Tools.min(q+5, MAX_CALLED_QUALITY);
				ntquals[i]=ntquals[i+1]=ntquals[i+2]=q2;
			}
			r.quality=ntquals;
		}
		return r;
	}
	
	private static final byte[][] QUALCACHE=new byte[1000][];
	

	public static final int STRANDMASK=1;
	public static final int MAPPEDMASK=(1<<1);
	public static final int PAIREDMASK=(1<<2);
	public static final int PERFECTMASK=(1<<3);
	public static final int AMBIMASK=(1<<4);
	public static final int RESCUEDMASK=(1<<5);
//	public static final int COLORMASK=(1<<6); //TODO:  Change to semiperfectmask?
	public static final int SYNTHMASK=(1<<7);
	public static final int DISCARDMASK=(1<<8);
	public static final int INVALIDMASK=(1<<9);
	public static final int SWAPMASK=(1<<10);
	public static final int SHORTMATCHMASK=(1<<11);
	
	public static final int PAIRNUMSHIFT=12;
	public static final int PAIRNUMMASK=(1<<PAIRNUMSHIFT);

	public static final int INSERTMASK=(1<<13);
	public static final int ADAPTERMASK=(1<<14);
	public static final int SECONDARYMASK=(1<<15);
	public static final int AAMASK=(1<<16);
	public static final int JUNKMASK=(1<<17);
	public static final int VALIDATEDMASK=(1<<18);
	public static final int TESTEDMASK=(1<<19);
	public static final int IRMASK=(1<<20);
	public static final int TRIMMEDMASK=(1<<21);
	
	private static final int[] maskArray=makeMaskArray(22); //Be sure this is big enough for all flags!

	public static boolean TO_UPPER_CASE=false;
	public static boolean LOWER_CASE_TO_N=false;
	public static boolean DOT_DASH_X_TO_N=false;
	public static boolean AVERAGE_QUALITY_BY_PROBABILITY=true;
	public static boolean FIX_HEADER=false;
	public static boolean ALLOW_NULL_HEADER=false;
	public static boolean SKIP_SLOW_VALIDATION=false;
	public static final boolean VALIDATE_BRANCHLESS=true;

	public static final int IGNORE_JUNK=0;
	public static final int FLAG_JUNK=1;
	public static final int FIX_JUNK=2;
	public static final int CRASH_JUNK=3;
	public static final int FIX_JUNK_AND_IUPAC=4;
	public static int JUNK_MODE=CRASH_JUNK;
	public static boolean IUPAC_TO_N=false;
	
	public static boolean U_TO_T=false;
	public static boolean COMPRESS_MATCH_BEFORE_WRITING=true;
	public static boolean DECOMPRESS_MATCH_ON_LOAD=true; //Set to false for some applications, like sorting, perhaps
	
	public static boolean ADD_BEST_SITE_TO_LIST_FROM_TEXT=true;
	public static boolean NULLIFY_BROKEN_QUALITY=false;
	public static boolean TOSS_BROKEN_QUALITY=false;
	public static boolean FLAG_BROKEN_QUALITY=false;
	public static boolean FLAT_IDENTITY=true;
	public static boolean VALIDATE_IN_CONSTRUCTOR=true;
	
	public static boolean verbose=false;

	/*--------------------------------------------------------------*/
	
	/** Offset for quality scores, normally 33, or 64 for some old Illumina data */
	private static final byte ASCII_OFFSET=33;
	
	/** Allow quality scores to be changed, typically for corrupt Illumina data */
	public static boolean CHANGE_QUALITY=true; //Cap all quality values between MIN_CALLED_QUALITY and MAX_CALLED_QUALITY
	
	/** Minimum allowed quality score for called (ACGT) bases */
	private static byte MIN_CALLED_QUALITY=2;
	
	/** Maximum allowed quality score for called (ACGT) bases */
	private static byte MAX_CALLED_QUALITY=50; //TODO: Find, replace, and test all instances of 41 (old value).
	
	/** Maximum allowed quality score for merged reads, which otherwise would normally be very high */
	public static byte MAX_MERGE_QUALITY=50;
	public static byte[] qMap=makeQmap(MIN_CALLED_QUALITY, MAX_CALLED_QUALITY);

	public static byte MIN_CALLED_QUALITY(){return MIN_CALLED_QUALITY;}
	public static byte MAX_CALLED_QUALITY(){return MAX_CALLED_QUALITY;}
	
	public static void setMaxCalledQuality(int x){
		x=Tools.mid(1, x, 93);
		if(x!=MAX_CALLED_QUALITY){
			MAX_CALLED_QUALITY=(byte)x;
			qMap=makeQmap(MIN_CALLED_QUALITY, MAX_CALLED_QUALITY);
		}
	}
	
	public static void setMinCalledQuality(int x){
		x=Tools.mid(0, x, 93);
		if(x!=MIN_CALLED_QUALITY){
			MIN_CALLED_QUALITY=(byte)x;
			qMap=makeQmap(MIN_CALLED_QUALITY, MAX_CALLED_QUALITY);
		}
	}

	public static byte capQuality(long q){
		return (byte)Tools.mid(MIN_CALLED_QUALITY, q, MAX_CALLED_QUALITY);
	}

	public static byte capQuality(byte q){
		return qMap[q];
	}
	
	public static byte capQuality(byte q, byte b){
		return AminoAcid.isFullyDefined(b) ? qMap[q] : 0;
	}
	
	private static byte[] makeQmap(byte min, byte max){
		byte[] r=(qMap==null ? new byte[128] : qMap);
		for(int i=0; i<r.length; i++){
			r[i]=(byte) Tools.mid(min, i, max);
		}
		return r;
	}
}
