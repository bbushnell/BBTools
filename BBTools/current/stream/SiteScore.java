package stream;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

import align2.GapTools;
import align2.MSA;
import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;



public final class SiteScore implements Comparable<SiteScore>, Cloneable, Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -8096245242590075081L;

	public SiteScore(int chrom_, byte strand_, int start_, int stop_, int hits_, int quickScore_){
		start=start_;
		stop=stop_;
		hits=hits_;
		quickScore=quickScore_;
		score=quickScore_;
		chrom=chrom_;
		strand=strand_;
//		assert(chrom_>=0) : this.toText()+"\nchrom_="+chrom_+", strand_="+strand_+", start_="+start_+", stop_="+stop_+", hits_="+hits_+", quickScore_="+quickScore_;
		assert(start_<=stop_) : this.toText()+"\nchrom_="+chrom_+", strand_="+strand_+", start_="+start_+", stop_="+stop_+", hits_="+hits_+", quickScore_="+quickScore_;
	}
	
	public SiteScore(int chrom_, byte strand_, int start_, int stop_, int hits_, int quickScore_, boolean rescued_, boolean perfect_){
		start=start_;
		stop=stop_;
		hits=hits_;
		quickScore=quickScore_;
		score=quickScore_;
		chrom=chrom_;
		strand=strand_;
		rescued=rescued_;
		perfect=perfect_;
		semiperfect=perfect;
		assert(start_<=stop_) : this.toText();
	}
	
	@Override
	public int compareTo(SiteScore other) {
		int x=other.score-score;
		if(x!=0){return x;}
		
		x=other.slowScore-slowScore;
		if(x!=0){return x;}
		
		x=other.pairedScore-pairedScore;
		if(x!=0){return x;}
		
		x=other.quickScore-quickScore;
		if(x!=0){return x;}
		
		x=chrom-other.chrom;
		if(x!=0){return x;}
		
		x=start-other.start;
		return x;
	}
	
	@Override
	public boolean equals(Object other){
		return compareTo((SiteScore)other)==0;
	}
	
	@Override
	public int hashCode() {
		assert(false) : "This class should not be hashed.";
		return super.hashCode();
	}
	
	@Override
	public String toString(){
		return toText().toString();
	}
	
//	9+2+1+9+9+1+1+4+4+4+4+gaps
	public CharSequence toText(){
		StringBuilder sb=new StringBuilder(53+(gaps==null ? 0 : gaps.length*10));
		sb.append(chrom);
		sb.append(',');
		sb.append(strand);
		sb.append(',');
		sb.append(start);
		sb.append(',');
		sb.append(stop);
		sb.append(',');
		sb.append((rescued ? 1 : 0));
		sb.append(',');
		sb.append((semiperfect ? 1 : 0));
		sb.append((perfect ? 1 : 0));
		sb.append(',');
		sb.append(hits);
		sb.append(',');
		sb.append(quickScore);
		sb.append(',');
		sb.append(slowScore);
		sb.append(',');
		sb.append(pairedScore);
		sb.append(',');
		sb.append(score);
		
		if(gaps!=null){
			sb.append(',');
			for(int i=0; i<gaps.length; i++){
				if(i>0){sb.append('~');}
				sb.append(gaps[i]);
			}
		}
		
		if(match!=null){
			if(gaps==null){sb.append(',');}
			sb.append(',');
			final char[] buffer=Shared.getTLCB(match.length);
			for(int i=0; i<match.length; i++){buffer[i]=(char)match[i];}
			sb.append(buffer, 0, match.length);
		}
		
		return sb;
//		chrom+","+strand+","+start+","+stop+","+(rescued ? 1 : 0)+","+
//		(perfect ? 1 : 0)+","+quickScore+","+slowScore+","+pairedScore+","+score;
	}
	
//	9+2+1+9+9+1+1+4+4+4+4+gaps
	public ByteBuilder toBytes(ByteBuilder sb){
		if(sb==null){sb=new ByteBuilder(53+(gaps==null ? 0 : gaps.length*10));}
		sb.append(chrom);
		sb.append(',');
		sb.append((int)strand);
		sb.append(',');
		sb.append(start);
		sb.append(',');
		sb.append(stop);
		sb.append(',');
		sb.append((rescued ? 1 : 0));
		sb.append(',');
		sb.append((semiperfect ? 1 : 0));
		sb.append((perfect ? 1 : 0));
		sb.append(',');
		sb.append(hits);
		sb.append(',');
		sb.append(quickScore);
		sb.append(',');
		sb.append(slowScore);
		sb.append(',');
		sb.append(pairedScore);
		sb.append(',');
		sb.append(score);
		
		if(gaps!=null){
			sb.append(',');
			for(int i=0; i<gaps.length; i++){
				if(i>0){sb.append('~');}
				sb.append(gaps[i]);
			}
		}
		
		if(match!=null){
			if(gaps==null){sb.append(',');}
			sb.append(',');
			sb.append(match);
		}
		
		return sb;
//		chrom+","+strand+","+start+","+stop+","+(rescued ? 1 : 0)+","+
//		(perfect ? 1 : 0)+","+quickScore+","+slowScore+","+pairedScore+","+score;
	}
	
	public boolean isSemiPerfect(byte[] bases){
		if(bases.length!=stop-start+1){return false;}
		byte[] ref=Data.getChromosome(chrom).array;

		//This block handles cases where the read runs outside the reference
		//Of course, padding the reference with 'N' would be better, but...
		int readStart=0;
		int readStop=bases.length;
		final int refStop=start+bases.length;
		int maxNoref=bases.length/2;

		if(start<0){
			readStart=0-start;
		}
		if(refStop>ref.length){
			int dif=(refStop-ref.length);
			readStop-=dif;
		}

		for(int i=readStart; i<readStop; i++){
			byte c=bases[i];
			byte r=ref[start+i];
			
//			assert(Tools.isUpperCase(c) && Tools.isUpperCase(r));
			if(c=='N'){return false;}
			if(c!=r){
				maxNoref--;
				if(maxNoref<0 || r!='N'){return false;}
			}
		}
		return true;
	}
	
	public boolean isPerfect(byte[] bases){
		if(bases.length!=stop-start+1 || start<0){return false;}
		byte[] ref=Data.getChromosome(chrom).array;
		if(stop>=ref.length){return false;}
		
		for(int i=0; i<bases.length; i++){
			byte c=bases[i];
			byte r=ref[start+i];
			assert(Tools.isUpperCase(c) && Tools.isUpperCase(r)) : "Lowercase letters detected: ref="+(char)r+", read="+(char)c+"\n"+new String(bases)+"\n"+
					"Please re-run with the 'tuc=t' flag (touppercase=true).";

			if((c!=r /* && (Tools.toUpperCase(c)!=Tools.toUpperCase(r))*/) || c=='N'){
				return false;
			}
		}
		return true;
	}
	
	
	public boolean setPerfectFlag(int maxScore, byte[] bases){
		if(maxScore==slowScore){
			assert(isPerfect(bases));
			return perfect=semiperfect=true;
		}
		return setPerfect(bases, false);
	}
	
	/** Sets "perfect" and "semiperfect" flags */
	public boolean setPerfect(byte[] bases){return setPerfect(bases, false);}
	
	/** Sets "perfect" and "semiperfect" flags, optionally assuming "perfect" flag is correct. */
	public boolean setPerfect(byte[] bases, boolean assumePerfectCorrect){
		if(bases.length!=stop-start+1){
			assert(!perfect || !assumePerfectCorrect) : perfect+", "+toString()+", "+
					new String(Data.getChromosome(chrom).array, Tools.max(0, start), (Tools.min(Data.chromLengths[chrom], stop)-start));
			perfect=false;
			semiperfect=false;
			assert(Read.CHECKSITE(this, bases, 0)) : new String(bases)+"\n"+this+"\n"; //123
			return perfect;
		}
		byte[] ref=Data.getChromosome(chrom).array;
		
		perfect=semiperfect=true;
		int refloc=start, readloc=0, N=0, max=Tools.min(stop, ref.length-1), nlimit=bases.length/2;
		if(start<0){
			N-=start;
			readloc-=start;
			refloc-=start;
			assert(!perfect || !assumePerfectCorrect);
			perfect=false;
		}
		if(stop>=ref.length){
			N+=(stop-ref.length+1);
			assert(!perfect || !assumePerfectCorrect);
			perfect=false;
		}
		if(N>nlimit){
			perfect=semiperfect=false;
			assert(Read.CHECKSITE(this, bases, 0)); //123
			return perfect;
		}
		
		final byte bn=(byte)'N';
		for(; refloc<=max; refloc++, readloc++){
			final byte c=bases[readloc];
			final byte r=ref[refloc];
			assert(Tools.isUpperCase(r) && Tools.isUpperCase(c)) :
				"\nAn input read appears to contain a non-upper-case base.  Please rerun with the 'touppercase' flag.\n"+
				"ref base = "+r+", read base = "+c+", TO_UPPER_CASE = "+Read.TO_UPPER_CASE+"\n"+(bases.length<=500 ? new String(bases) : "")+"\n";
			if(c!=r || c==bn){
				perfect=false;
				if(c==bn){semiperfect=false;}
				if(r!=bn || (N=N+1)>nlimit){
					semiperfect=false;
					assert(Read.CHECKSITE(this, bases, 0)); //123
					return semiperfect;
				}
			}
		}
		
		semiperfect=(semiperfect && (N<=nlimit));
		perfect=(perfect && semiperfect && (N==0));
		assert(Read.CHECKSITE(this, bases, 0)); //123
		return perfect;
	}
	
	public final boolean overlaps(SiteScore ss){
		return chrom==ss.chrom && strand==ss.strand && overlap(start, stop, ss.start, ss.stop);
	}
	public final boolean overlaps(SiteScore ss, boolean ignoreStrand){
		return chrom==ss.chrom && (ignoreStrand || strand==ss.strand) && overlap(start, stop, ss.start, ss.stop);
	}
	private static boolean overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1;
	}
	
	public static String header() {
		return "chrom,strand,start,stop,rescued,semiperfect+perfect,hits,quickScore,slowScore,pairedScore,score,match";
	}
	
	public static SiteScore fromText(String s){
//		System.err.println("Trying to make a SS from "+s);
		String line[]=s.split(",");
		
		SiteScore ss;

		assert(line.length>=11 && line.length<=13) : "\n"+line.length+"\n"+s+"\n"+Arrays.toString(line)+"\n";
		int chrom=Byte.parseByte(line[0].charAt(0)=='*' ? line[0].substring(1) : line[0]);
		byte strand=Byte.parseByte(line[1]);
		int start=Integer.parseInt(line[2]);
		int stop=Integer.parseInt(line[3]);
		boolean rescued=Integer.parseInt(line[4])==1;
//		[1, 1, 9397398, 9398220, 0, 00, 20, 8701, 9084, 0, 9084, 9397398~9397471~9398145~9398220]
		int p=Integer.parseInt(line[5], 2);
//		assert(false) : line[5]+"->"+p;
		boolean perfect=(p&1)==1;
		boolean semiperfect=(p&2)==2;
		int hits=Integer.parseInt(line[6]);
		int quickScore=Integer.parseInt(line[7]);
		int swscore=Integer.parseInt(line[8]);
		int pairedScore=Integer.parseInt(line[9]);
		int score=Integer.parseInt(line[10]);
		ss=new SiteScore(chrom, strand, start, stop, hits, quickScore, rescued, perfect);
		ss.setScore(score);
		ss.setSlowPairedScore(swscore, pairedScore);
		ss.semiperfect=semiperfect;
		
		if(line.length>11){
			if(line[11]!=null && line[11].length()>0){
				String[] gstring=line[11].split("~");
				ss.gaps=new int[gstring.length];
				for(int i=0; i<gstring.length; i++){
					ss.gaps[i]=Integer.parseInt(gstring[i]);
				}
			}
		}
		
		if(line.length>12){
			ss.match=line[12].getBytes();
		}
		
		return ss;
	}
	
	public boolean positionalMatch(SiteScore b, boolean testGaps){
//		return chrom==b.chrom && strand==b.strand && start==b.start && stop==b.stop;
		if(chrom!=b.chrom || strand!=b.strand || start!=b.start || stop!=b.stop){
			return false;
		}
		if(!testGaps || (gaps==null && b.gaps==null)){return true;}
		if((gaps==null) != (b.gaps==null)){return false;}
		if(gaps.length!=b.gaps.length){return false;}
		for(int i=0; i<gaps.length; i++){
			if(gaps[i]!=b.gaps[i]){return false;}
		}
		return true;
	}
	
	public byte[] getScaffoldName(boolean requireSingleScaffold){
		byte[] name=null;
		if(!requireSingleScaffold || Data.isSingleScaffold(chrom, start, stop)){
			int idx=Data.scaffoldIndex(chrom, (start+stop)/2);
			name=Data.scaffoldNames[chrom][idx];
			//				int scaflen=Data.scaffoldLengths[chrom][idx];
			//				a1=Data.scaffoldRelativeLoc(chrom, start, idx);
			//				b1=a1-start1+stop1;
		}
		return name;
	}
	
	public static class PositionComparator implements Comparator<SiteScore>{
		
		private PositionComparator(){}
		
		@Override
		public int compare(SiteScore a, SiteScore b) {
			if(a.chrom!=b.chrom){return a.chrom-b.chrom;}
			if(a.start!=b.start){return a.start-b.start;}
			if(a.stop!=b.stop){return a.stop-b.stop;}
			if(a.strand!=b.strand){return a.strand-b.strand;}
			if(a.score!=b.score){return b.score-a.score;}
			if(a.slowScore!=b.slowScore){return b.slowScore-a.slowScore;}
			if(a.quickScore!=b.quickScore){return b.quickScore-a.quickScore;}
			if(a.perfect!=b.perfect){return a.perfect ? -1 : 1;}
			if(a.rescued!=b.rescued){return a.rescued ? 1 : -1;}
			return 0;
		}
		
		public void sort(ArrayList<SiteScore> list){
			if(list==null || list.size()<2){return;}
			Shared.sort(list, this);
		}
		
		public void sort(SiteScore[] list){
			if(list==null || list.length<2){return;}
			Arrays.sort(list, this);
		}
		
	}
	
	public SiteScore copy(){
		SiteScore ss2=this.clone();
		if(gaps!=null){ss2.gaps=ss2.gaps.clone();}
		return ss2;
	}
	
	@Override
	public SiteScore clone(){
		try {
			return (SiteScore)super.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		throw new RuntimeException();
	}
	
	public boolean isInBounds(){
		ChromosomeArray cha=Data.getChromosome(chrom);
		return (start>=0 && stop<=cha.maxIndex);
	}
	
	public boolean matchContainsXY(){
		if(match==null || match.length<1){return false;}
		final byte a=match[0], b=match[match.length-1];
		return (a=='X' ||a=='Y' || b=='X' || b=='Y');
	}
	
	public boolean matchContainsAB(){
		if(match==null || match.length<1){return false;}
		final byte a=match[0], b=match[match.length-1];
		return (a=='A' ||a=='B' || b=='A' || b=='B');
	}
	
	public boolean matchContainsC(){
		if(match==null || match.length<1){return false;}
		final byte a=match[0], b=match[match.length-1];
		return (a=='C' || b=='C');
	}
	
	public boolean isCorrect(int chrom_, byte strand_, int start_, int stop_, int thresh){
		if(chrom_!=chrom || strand_!=strand){return false;}
		if(thresh<=0){return start_==start && stop_==stop;}
		return Tools.absdif(start_, start)<=thresh || Tools.absdif(stop_, stop)<=thresh;
	}
	
	public int leftPaddingNeeded(int tiplen, int maxIndel){
		if(match==null || match.length<1){return 0;}
		
		int neutral=0, insertion=0, deletion=0, xy=0;
		{
			int mloc=0;
			for(; mloc<match.length; mloc++){
				byte c=match[mloc];
				if(c=='I'){insertion++;}
				else if(c=='X' || c=='Y'){xy++;}
				else if(c=='D'){return insertion+xy;}
				else{
					neutral++;
					if(mloc>=tiplen){break;}
				}
			}
		}
		
		if(insertion>maxIndel || xy>0 || match[0]=='I'){return insertion+xy;}
		return 0;
	}
	
	public int rightPaddingNeeded(int tiplen, int maxIndel){
		if(match==null || match.length<1){return 0;}
		final int lastIndex=match.length-1;
		
		int neutral=0, insertion=0, deletion=0, xy=0;
		{
			int mloc=lastIndex;
			for(int min=lastIndex-tiplen; mloc>=0; mloc--){
				byte c=match[mloc];
				if(c=='I'){insertion++;}
				else if(c=='X' || c=='Y'){xy++;}
				else if(c=='D'){return insertion+xy;}
				else{
					neutral++;
					if(mloc>=tiplen){break;}
				}
			}
		}
		
		if(insertion>maxIndel || xy>0 || match[lastIndex]=='I'){return insertion+xy;}
		return 0;
	}
	
	/** Assumes bases are rcomped as needed */
	public int unclip(byte[] bases, byte[] rbases){
		int unclipped=0;
		if(match==null || match.length<1){return unclipped;}
		assert(lengthsAgree());
		
		for(int mpos=0, bpos=0, rpos=start; mpos<match.length && rpos<rbases.length; mpos++, bpos++, rpos++){
			if(match[mpos]!='C'){break;}
			if(rpos>=0){
				byte b=bases[bpos];
				byte r=rbases[rpos];
				
				if(!AminoAcid.isFullyDefined(b) || !AminoAcid.isFullyDefined(r)){match[mpos]='N';}
				else if(b==r){match[mpos]='m';}
				else{match[mpos]='S';}
				unclipped++;
			}
		}
		assert(lengthsAgree());
		
		for(int mpos=match.length-1, bpos=bases.length-1, rpos=stop; mpos>=0 && rpos>=0; mpos--, bpos--, rpos--){
			if(match[mpos]!='C'){break;}
			if(rpos<rbases.length){
				byte b=bases[bpos];
				byte r=rbases[rpos];
				if(b=='N' || r=='N'){match[mpos]='N';}
				else if(b==r){match[mpos]='m';}
				else{match[mpos]='S';}
				unclipped++;
			}
		}
		assert(lengthsAgree());
		
		return unclipped;
	}
	
	/** Simply replaces terminal 'I', 'X', and 'Y' with 'C' and adjusts length
	 * TODO: Also clip out-of-bounds. */
	public int clipTipIndels(int rlen){
		int clipped=0;
		if(match==null || match.length<1){return clipped;}
		assert(lengthsAgree());
		
		for(int mpos=0, rpos=start; mpos<match.length; mpos++){
			final byte m=match[mpos];
			if(m=='C'){
				//Do nothing
				rpos++;
			}else if(m=='m' || m=='S'){
				break;
			}else{
				if(m=='I'){
					start--;
				}else if(m=='X'){
					rpos++;
				}else if(m=='Y'){//Should not happen
					start--;
				}else if(m=='D'){
					//Do nothing
				}else if(rpos>=0){
					break;
				}else{
					assert(m=='N') : (char)m+"\n"+mpos+", "+rpos+"\n"+this;
					rpos++;
				}
				clipped++;
				match[mpos]='C';
			}
		}
		assert(clipped==0 || lengthsAgree());
		
		for(int mpos=match.length-1, rpos=stop; mpos>=0; mpos--){
			final byte m=match[mpos];
			if(m=='C'){
				//Do nothing
				rpos--;
			}else if(m=='m' || m=='S'){
				break;
			}else{
				if(m=='I'){
					stop++;
				}else if(m=='X'){//Should not happen
					assert(false);
					rpos--;
				}else if(m=='Y'){
					stop++;
				}else if(m=='D'){
					//Do nothing
				}else if(rpos<rlen){
					break;
				}else{
					assert(m=='N') : (char)m+"\n"+mpos+", "+rpos+"\n"+this;
					rpos--;
				}
				clipped++;
				match[mpos]='C';
			}
		}
		assert(clipped==0 || lengthsAgree());
		
		return clipped;
	}
	
//	/** Simply replaces terminal 'I', 'X', and 'Y' with 'C' and adjusts length
//	 * TODO: Also clip out-of-bounds. */
//	public int clipTipIndels(int rlen){
//		int clipped=0;
//		if(match==null || match.length<1){return clipped;}
//		assert(lengthsAgree());
//
//		for(int mpos=0, rpos=start; mpos<match.length; mpos++){
//			final byte m=match[mpos];
//			if(m=='I'){
//				start++;
//			}else if(m=='X' || m=='Y' || m=='D'){
//				//Do nothing
//			}else{
//				break;
//			}
//			clipped++;
//			match[mpos]='C';
//		}
//		assert(clipped==0 || lengthsAgree());
//
//		for(int mpos=match.length-1, rpos=stop; mpos>=0; mpos--){
//			final byte m=match[mpos];
//			if(m=='I'){
//				start++;
//			}else if(m=='X' || m=='Y' || m=='D'){
//				//Do nothing
//			}else{
//				break;
//			}
//			clipped++;
//			match[mpos]='C';
//		}
//		assert(clipped==0 || lengthsAgree());
//
//		return clipped;
//	}
	
	public boolean clipTipIndels(byte[] bases, byte[] basesM, int tiplen, int maxIndel, MSA msa){
		return this.plus() ? clipTipIndels(bases, tiplen, maxIndel, msa) : clipTipIndels(basesM, tiplen, maxIndel, msa);
	}
	
	public boolean clipTipIndels(byte[] bases, int tiplen, int maxIndel, MSA msa){
		if(match==null || match.length<maxIndel){return false;}
		if(verbose){
			System.err.println("Calling clipTipIndels:\n"+new String(match));
			System.err.println("slowScore="+slowScore+", pairedScore="+pairedScore);
		}
		assert(lengthsAgree());
		boolean left=clipLeftTipIndel(bases, tiplen, maxIndel);
		assert(lengthsAgree());
		boolean right=clipRightTipIndel(bases, tiplen, maxIndel);
		assert(lengthsAgree());

		if(verbose){System.err.println("left="+left+", right="+right+", match="+new String(match));}
		if(left || right){
			unclip(bases);
			assert(lengthsAgree());
			int oldScore=slowScore;
			if(verbose){System.err.println("oldScore="+oldScore+", slowScore="+slowScore+", pairedScore="+pairedScore+", newPairedScore="+(pairedScore+(slowScore-oldScore)));}
			setSlowScore(msa.score(match));
			setScore(score+(slowScore-oldScore));
			this.setPerfect(bases);
		}
		if(verbose){System.err.println("After clipTipIndels:\n"+new String(match));}
		return left | right;
	}
	
	public boolean clipLeftTipIndel(byte[] bases, int tiplen, int maxIndel){
		if(match==null || match.length<maxIndel){return false;}
		if(match[0]=='C' || match[0]=='Y' || match[0]=='X'){return false;}
		
		int neutral=0, insertion=0, deletion=0;
		{
			int mloc=0;
			for(; mloc<match.length; mloc++){
				byte c=match[mloc];
				if(c=='I'){insertion++;}
				else if(c=='D'){deletion++;}
				else{
					neutral++;
					if(mloc>=tiplen){break;}
				}
			}
			while(mloc>=0 && match[mloc]=='m'){mloc--; neutral--;}
		}
		if(insertion<=maxIndel && deletion<=4*maxIndel){return false;}
		assert(mappedLength()==matchLength() || matchContainsXY()) :
			"start="+start+", stop="+stop+", maplen="+mappedLength()+", matchlen="+matchLength()+"\n"+new String(match)+"\n"+new String(bases)+"\n\n"+this;
		
		int sum=neutral+insertion+deletion;
		if(deletion>0){
			byte[] temp=new byte[match.length-deletion];
			int i=0, j=0;
			for(; i<sum; i++){
				byte c=match[i];
				if(c=='D'){
					//Do nothing
				}else{
					temp[j]=match[i];
					j++;
				}
			}
			for(; i<match.length; i++, j++){
				temp[j]=match[i];
			}
			assert(i==match.length && j==temp.length) : i+", "+j+", "+match.length+", "+temp.length+"\n"+new String(match)+"\n"+new String(bases)+"\n"+this+"\n";
			match=temp; //Be sure to percolate this to the read!
		}
		
		sum=neutral+insertion;
		for(int i=0; i<sum; i++){
			match[i]='C';
		}
		
		final int dif=(insertion-deletion);
		incrementStart(-dif);
		assert(mappedLength()==matchLength() || matchContainsXY());
		
		return true;
	}
	
	public boolean clipRightTipIndel(final byte[] bases, final int tiplen, final int maxIndel){
		if(match==null || match.length<maxIndel){return false;}
		final int lastIndex=match.length-1;
		if(match[lastIndex]=='C' || match[lastIndex]=='Y' || match[lastIndex]=='X'){return false;}
		
		if(verbose){System.err.println("mappedLength="+mappedLength()+", matchLength()="+matchLength());}
		
		int neutral=0, insertion=0, deletion=0;
		{
			int mloc=lastIndex;
			for(int min=lastIndex-tiplen; mloc>=0; mloc--){
				byte c=match[mloc];
				if(c=='I'){insertion++;}
				else if(c=='D'){deletion++;}
				else{
					neutral++;
					if(mloc<=min){break;}
				}
			}
			while(mloc<match.length && match[mloc]=='m'){mloc++; neutral--;}
		}
		if(insertion<=maxIndel && deletion<=4*maxIndel){return false;}
		assert(mappedLength()==matchLength() || matchContainsXY()) : mappedLength()+", "+matchLength()+"\n"+new String(match)+"\n"+new String(bases);
		
		int sum=neutral+insertion+deletion;
		final int limit=match.length-sum;
		if(deletion>0){
			byte[] temp=new byte[match.length-deletion];
			int i=0, j=0;
			for(; i<limit; i++, j++){
				temp[j]=match[i];
			}
			for(; i<match.length; i++){
				byte c=match[i];
				if(c=='D'){
					//Do nothing
				}else{
					temp[j]=match[i];
					j++;
				}
			}
			assert(i==match.length && j==temp.length) : i+", "+j+", "+match.length+", "+temp.length+
				"\n"+new String(match)+"\n"+new String(temp)+"\n"+new String(bases)+"\n"+this+"\n";
			match=temp; //Be sure to percolate this to the read!
		}
		
		sum=neutral+insertion;
		for(int i=limit; i<match.length; i++){
			match[i]='C';
		}


		if(verbose){System.err.println("Final: "+new String(match));}
		final int dif=(insertion-deletion);
		if(verbose){System.err.println("mappedLength="+mappedLength()+", matchLength()="+matchLength()+", dif="+dif);}
		incrementStop(dif);
		assert(mappedLength()==matchLength() || matchContainsXY()) : mappedLength()+", "+matchLength()+", "+neutral+", "+insertion+", "+deletion+", "+dif;
		
		return true;
	}
	
	public boolean unclip(final byte[] bases){
		if(match==null || match.length<1){return false;}
		if(verbose){System.err.println("Calling unclip on "+new String(match));}
		{
			byte first=match[0], last=match[match.length-1];
//			if(first=='X' || last=='Y'){return false;}
			if(first!='C' && last!='C'){
				if(verbose){System.err.println("No unclipping needed.");}
				return false;
			}
		}
		assert(lengthsAgree()) : new String(bases)+"\n"+this;
		final ChromosomeArray ca=Data.getChromosome(chrom);
		for(int rloc=start, cloc=0, mloc=0; mloc<match.length; mloc++){
			final byte m=match[mloc];
			
			if(m=='C'){
				final byte c=bases[cloc];
				final byte r=ca.get(rloc);
				if(!AminoAcid.isFullyDefined(c) || !AminoAcid.isFullyDefined(r)){
					match[mloc]='N';
				}else{
					match[mloc]=(byte)(c==r ? 'm' : 'S');
				}
				rloc++;
				cloc++;
			}else if(m=='m' || m=='N' || m=='S'){
				rloc++;
				cloc++;
			}else if(m=='X' || m=='Y'){
				rloc++;
				cloc++;
			}else if(m=='I'){
				cloc++;
			}else if(m=='D'){
				rloc++;
			}else{
				throw new RuntimeException("Unsupported symbol: ASCII "+(int)m);
			}
		}
		if(verbose){System.err.println("After unclip: "+new String(match));}
		assert(lengthsAgree()) : new String(bases)+"\n"+this;
		return true;
	}
	
	/** TODO: Test
	 * Attempt to extend match/N symbols where there are X and Y symbols
	 * */
	public boolean fixXY(byte[] bases, boolean nullifyOnFailure, MSA msa){
		if(verbose && match!=null){System.err.println("Calling fixXY:\n"+new String(match));}
		if(!matchContainsXY()){return true;}
		if(verbose){System.err.println("lengthsAgree: "+this.lengthsAgree());}
		
		boolean disable=false;
		if(disable){
			if(nullifyOnFailure){
				match=null;
			}
			return false;
		}
		
//		if(match==null || match.length<1){return false;} //Already covered
		final ChromosomeArray ca=Data.getChromosome(chrom);
//		final int tip=3;
		boolean success=true;
		final float maxSubRate=0.4f;
		final int maxSubs=5;
		
		{//Process left side
			if(verbose){System.err.println("Processing left side.  Success="+success+", start="+start+", stop="+stop+", match=\n"+new String(match));}
			int mloc=0;
			while(mloc<match.length && (match[mloc]=='X' || match[mloc]=='Y')){mloc++;}
			if(mloc>=match.length || mloc>=bases.length){success=false;}
			else if(mloc>0){
				mloc--;//Location of last X or Y on left side
				final int numX=mloc+1;
				int rloc=start+mloc, cloc=mloc;
				int subs=0, firstSub=-1;
				while(mloc>=0){
					byte m=match[mloc];
					byte c=bases[cloc];
					byte r=ca.get(rloc);
					assert(m=='X' || m=='Y') : (char)m+", "+mloc+", "+(char)c+", "+(char)r+"\n"+new String(bases)+"\n"+this.toString();
					if(r=='N' || c=='N'){match[mloc]='N';}
					else if(c==r){match[mloc]='m';}
//					else if(mloc<=tip){match[mloc]='S';}
					else{
						match[mloc]='S';
						subs++;
						if(subs==1){firstSub=mloc;}
					}
//					else{
//						if(verbose){System.err.println("A: Set success to false.");}
//						success=false;
//						break;
//					}
					mloc--;
					rloc--;
					cloc--;
				}
				if(success && mappedLength()!=matchLength()){incrementStart(-numX);}
				if(subs>maxSubs && subs>numX*maxSubRate){
					if(verbose){System.err.println("Failed to correct alignment; clipping left side of read.");}
					for(int i=0; i<=firstSub; i++){
						match[i]='C';
					}
					if(!lengthsAgree()){
						assert(false);
						incrementStart(firstSub+1);
						assert(lengthsAgree());
					}
				}
			}
			if(verbose){System.err.println("Finished left side.  Success="+success+", start="+start+", stop="+stop+", match=\n"+new String(match));}
			if(verbose){System.err.println("lengthsAgree: "+this.lengthsAgree());}
		}
		
		if(success){//Process right side
			if(verbose){System.err.println("Processing right side.  Success="+success+", start="+start+", stop="+stop+", match=\n"+new String(match));}
			int mloc=match.length-1;
			while(mloc>=0 && (match[mloc]=='X' || match[mloc]=='Y')){mloc--;}
			int dif=match.length-1-mloc;
			if(mloc<0){
				if(verbose){System.err.println("B: Set success to false.");}
				success=false;
			}
			else if(dif>0){
				mloc++;//Location of first X or Y on right side
				final int numX=match.length-mloc;
				int rloc=stop-dif+1, cloc=bases.length-dif;
				int subs=0, firstSub=-1;
				if(cloc<0){
					if(verbose){System.err.println("C: Set success to false.");}
					success=false;
				}else{
//					final int tip2=match.length-tip;
					while(mloc<match.length){
						byte m=match[mloc];
						byte c=bases[cloc];
						byte r=ca.get(rloc);
						assert(m=='X' || m=='Y') : (char)m+", "+mloc+", "+(char)c+", "+(char)r+"\n"+new String(bases)+"\n"+this.toString();
						if(r=='N' || c=='N'){match[mloc]='N';}
						else if(c==r){match[mloc]='m';}
//						else if(mloc>=tip2){match[mloc]='S';}
						else{
							match[mloc]='S';
							subs++;
							if(subs==1){firstSub=mloc;}
						}
//						else{
//							if(verbose){System.err.println("D: Set success to false.");}
//							success=false;
//							break;
//						}
						mloc++;
						rloc++;
						cloc++;
					}
				}
				if(success){
					if(verbose){System.err.println("A: Start="+start+", stop="+stop+", numX="+numX+", lengthsAgree()="+lengthsAgree());}
					if(mappedLength()!=matchLength()){incrementStop(numX);}
					if(verbose){System.err.println("B: Start="+start+", stop="+stop+", numX="+numX+", lengthsAgree()="+lengthsAgree());}
					if(subs>maxSubs && subs>numX*maxSubRate){
						if(verbose){System.err.println("Failed to correct alignment; clipping right side of read.");}
						for(int i=firstSub; i<match.length; i++){
							match[i]='C';
						}
						int clipped=match.length-firstSub+1;
						if(!lengthsAgree()){
							assert(false);
							incrementStop(-clipped);
							assert(lengthsAgree());
						}
						if(verbose){System.err.println("C: Start="+start+", stop="+stop+", numX="+numX+", lengthsAgree()="+lengthsAgree());}
					}
					if(verbose){System.err.println("mappedLength()="+mappedLength()+", matchLength()="+matchLength()+", numX="+numX);}
				}
			}
			if(verbose){System.err.println("Finished right side.  Success="+success+", start="+start+", stop="+stop+", match=\n"+new String(match));}
			if(verbose){System.err.println("lengthsAgree: "+this.lengthsAgree());}
		}
		
		success=success && !matchContainsXY()/* && mappedLength()==matchLength()*/;
		if(!success){
			if(verbose){System.err.println("E: Set success to false.");}
			if(nullifyOnFailure){match=null;}
		}else{
			if(verbose){System.err.println("E: Success!");}
		}
		
		if(match!=null){
			int oldScore=slowScore;
			setSlowScore(msa.score(match));
			setScore(score+(slowScore-oldScore));
		}
		if(verbose){System.err.println("lengthsAgree: "+this.lengthsAgree());}
		
		setPerfect(bases); //Fixes a rare bug
		return success;
	}
	
	public boolean lengthsAgree(){
		return match==null ? true : matchLength()==mappedLength();
	}
	
	public int mappedLength(){
		return stop-start+1;
	}
	
	public int matchLength(){
		assert(match!=null);
		return Read.calcMatchLength(match);
	}

//	public boolean plus(){return strand()==Gene.PLUS;}
//	public boolean minus(){return strand()==Gene.MINUS;}
//
//	public final byte strand(){return (byte)(flags&strandMask);}
//	public boolean rescued(){return (flags&rescuedMask)!=0;}
//	public boolean perfect(){return (flags&perfectMask)!=0;}
//	public boolean semiperfect(){return (flags&semiperfectMask)!=0;}
//
//	public final int setStrand(int x){
//		assert(x==0 || x==1);
//		if(x==0){flags=(flags&~strandMask);}
//		else{flags=(flags|strandMask);}
//		assert(strand()==x);
//		return x;
//	}
//	public boolean setRescued(boolean b){
//		if(b){flags=(flags|rescuedMask);}
//		else{flags=(flags&~rescuedMask);}
//		assert(rescued()==b);
//		return b;
//	}
//	public boolean setPerfect(boolean b){
//		if(b){flags=(flags|semiperfectMask);}
//		else{flags=(flags&~semiperfectMask);}
//		assert(perfect()==b);
//		return b;
//	}
//	public boolean setSemiperfect(boolean b){
//		if(b){flags=(flags|semiperfectMask);}
//		else{flags=(flags&~semiperfectMask);}
//		assert(semiperfect()==b);
//		return b;
//	}

	public boolean plus(){return strand==Shared.PLUS;}
	public boolean minus(){return strand==Shared.MINUS;}
	public boolean perfect(){return perfect;}
	public boolean semiperfect(){return semiperfect;}
	public boolean rescued(){return rescued;}
	public byte strand(){return strand;}
	
	public final byte strand;
	public boolean rescued=false;
	public boolean perfect=false;
	public boolean semiperfect=false;

	public void setPerfect(){
		perfect=semiperfect=true;
		gaps=null;
	}
	public void incrementStart(int x){setStart(start+x);}
	public void incrementStop(int x){setStop(stop+x);}
	public void setLimits(int a, int b){
		start=a;
		stop=b;
		if(gaps!=null){
			gaps[0]=a;
			gaps[gaps.length-1]=b;
			if(!CHECKGAPS()){gaps=GapTools.fixGaps(this);}
			assert(CHECKGAPS()) : Arrays.toString(gaps);
		}
	}
	
	//Seems to be no longer needed after a change to calculating Y symbols.
	@Deprecated
	public void fixLimitsXY(){
//		if(match==null || match.length<1){return;}
//		int x=0, y=0;
//		for(int i=0; i<match.length; i++){
//			if(match[i]=='X'){x++;}else{break;}
//		}
//		for(int i=match.length-1; i>=0; i--){
//			if(match[i]=='Y'){y++;}else{break;}
//		}
////		if((x!=0 || y!=0) && !lengthsAgree()){
////			setLimits(start-x, stop+y);
////		}
//		if((y!=0)){
////			setLimits(start, stop+y);
//		}
	}
	
	public void setStart(int a){
		start=a;
		if(gaps!=null){
			gaps[0]=a;
			if(gaps[0]>gaps[1]){
				gaps=GapTools.fixGaps(this);
			}
			assert(CHECKGAPS()) : Arrays.toString(gaps);
		}
	}
	public void setStop(int b){
		stop=b;
		if(gaps!=null){
			gaps[gaps.length-1]=b;
			if(gaps.length-1>gaps.length-2){gaps=GapTools.fixGaps(this);}
			assert(CHECKGAPS()) : Arrays.toString(gaps);
		}
	}
	public boolean CHECKGAPS(){
		if(gaps==null){return true;}
		if(gaps.length==0 || ((gaps.length&1)==1)){return false;}
		for(int i=1; i<gaps.length; i++){
			if(gaps[i-1]>gaps[i]){return false;}
		}
		return gaps[0]==start && gaps[gaps.length-1]==stop;
	}
	
	public int start(){return start;}
	public int stop(){return stop;}
	public void setSlowScore(int x){
//		assert(x!=-1);
		if(verbose){System.err.println("Before setSlowScore: x="+x+", quick="+quickScore+", slow="+slowScore+", paired="+pairedScore);}
//		assert(slowScore<=0 || pairedScore<=0 || pairedScore>=slowScore) : "x="+x+", quick="+quickScore+", slow="+slowScore+", paired="+pairedScore+"\n"+this;  //Correct, but temporarily disabled for stability
		if(x<=0){
			pairedScore=slowScore=x;
		}else if(pairedScore<=0){
			slowScore=x;
		}else{
			assert(pairedScore>=slowScore) : this;
			if(pairedScore>0){
				if(slowScore>0){
					pairedScore=x+(pairedScore-slowScore);
				}else{
					pairedScore=x+1;
				}
			}
		}
		slowScore=x;
		if(verbose){System.err.println("After setSlowScore: "+this);}
//		assert(pairedScore<=0 || pairedScore>=slowScore) : "quick="+quickScore+", slow="+slowScore+", paired="+pairedScore+"\n"+this;  //Correct, but temporarily disabled for stability
	}
	public void setPairedScore(int x){
//		assert(x==0 || slowScore<=0 || x>=slowScore) : "x="+x+", quick="+quickScore+", slow="+slowScore+", paired="+pairedScore+"\n"+this;  //Correct, but temporarily disabled for stability
		pairedScore=x;
//		assert(slowScore<=0 || pairedScore<=0 || pairedScore>=slowScore) : "x="+x+", quick="+quickScore+", slow="+slowScore+", paired="+pairedScore+"\n"+this;  //Correct, but temporarily disabled for stability
	}
	public void setSlowPairedScore(int x, int y){
//		assert(slowScore<=0 || pairedScore<=0 || pairedScore>=slowScore) : "x="+x+", quick="+quickScore+", slow="+slowScore+", paired="+pairedScore+"\n"+this;  //Correct, but temporarily disabled for stability
		slowScore=x;
		pairedScore=y;
//		assert(pairedScore<=0 || pairedScore>=slowScore) : this;  //Correct, but temporarily disabled for stability
	}
	public void setScore(int x){
		score=x;
	}

	public int start;
	public int stop;
	public int quickScore;
	public int score;
	public int slowScore;
	public int pairedScore;
	public int hits;
	public final int chrom;
	
	public long flags; //TODO Use this instead of fields
	
	public int[] gaps; //Limits of large gaps
	public byte[] match;
	

	public static final PositionComparator PCOMP=new PositionComparator();
	public static final long strandMask=(1L<<0);
	public static final long rescuedMask=(1L<<1);
	public static final long perfectMask=(1L<<2);
	public static final long semiperfectMask=(1L<<3);
	public static boolean verbose=false;
	
}
