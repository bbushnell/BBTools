package pacbio;

import shared.Shared;
import stream.SiteScore;
import stream.SiteScoreR;

/**
 * @author Brian Bushnell
 * @date Jul 24, 2012
 *
 */
public class SiteR {
	
	public SiteR(SiteScoreR ssr){
		this(ssr.start, ssr.stop, ssr.chrom, ssr.strand, ssr.numericID, ssr.pairnum);
	}
	
	public SiteR(int start_, int stop_, int chrom, byte strand, long numericID, int pairnum){
		start=start_;
		stop=stop_;
		if((pairnum&1)==0){
			idPairnum=numericID;
		}else{
			idPairnum=-numericID;
		}
		if(strand==Shared.PLUS){
			chromStrand=chrom;
		}else{
			chromStrand=-chrom;
		}
		assert(chrom==chrom());
		assert(strand==strand());
		assert(numericID==numericID());
		assert(pairnum==pairNum());
	}
	
	public boolean equals(SiteScore other){
		if(other.start!=start){return false;}
		if(other.stop!=stop){return false;}
		if(other.chrom!=chrom()){return false;}
		if(other.strand!=strand()){return false;}
		return true;
	}
	
	public boolean equals(SiteScoreR other){
		if(other.start!=start){return false;}
		if(other.stop!=stop){return false;}
		if(other.chrom!=chrom()){return false;}
		if(other.strand!=strand()){return false;}
		return true;
	}
	
	public StringBuilder toTextRecursive(StringBuilder sb){
		if(sb==null){sb=new StringBuilder();}else{sb.append(" ");}
		sb.append("("+toText()+")");
		if(next!=null){next.toTextRecursive(sb);}
		return sb;
	}
	
	public StringBuilder toText(){
		StringBuilder sb=new StringBuilder();
		sb.append(start).append(',');
		sb.append(stop).append(',');
		sb.append(chrom()).append(',');
		sb.append(strand()).append(',');
		sb.append(numericID()).append(',');
		sb.append(pairNum());
		return sb;
	}
	
	@Override
	public String toString(){
		return toText().toString();
	}
	
	public final int start;
	public final int stop;
	public final int chromStrand;
	public final long idPairnum;
	public SiteR next;

	public long numericID(){return idPairnum>=0 ? idPairnum : -idPairnum;}
	public int pairNum(){return idPairnum>=0 ? 0 : 1;}
	public int chrom(){return chromStrand>=0 ? chromStrand : -chromStrand;}
	public byte strand(){return chromStrand>=0 ? (byte)0 : (byte)1;};
	public int listLength(){
		int i=1;
		SiteR sr=this;
		while(sr.next!=null){
			sr=sr.next;
			i++;
		}
		return i;
	}
	
}
