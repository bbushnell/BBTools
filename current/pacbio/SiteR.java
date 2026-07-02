package pacbio;

import shared.Shared;
import stream.SiteScore;
import stream.SiteScoreR;

/**
 * Compact representation of a genomic alignment site with packed coordinate storage.
 * Stores genomic coordinates, strand, chromosome, and read metadata in a memory-efficient
 * format using bit packing. Forms linked lists for multiple alignments per read.
 *
 * @author Brian Bushnell
 * @date Jul 24, 2012
 */
public class SiteR {
	
	public SiteR(SiteScoreR ssr){
		this(ssr.start, ssr.stop, ssr.chrom, ssr.strand, ssr.numericID, ssr.pairnum);
	}
	
	public SiteR(int start_, int stop_, int chrom, byte strand, long numericID, int pairnum){
		start=start_;
		stop=stop_;
		//TODO: Possible bug [pacbio/SiteR#001] LOW→MEDIUM (shared class): SIGN-MAGNITUDE ZERO COLLISION. pairnum is encoded as
		//the sign of numericID (odd→negate), but -0==0, so for numericID==0 with pairnum==1 (odd) idPairnum=0 → pairNum()
		//(L86, idPairnum>=0?0:1) returns 0, NOT 1. numericID is 0-based, so this is exactly READ 0's MATE — reachable on the
		//very first read pair. The assert(pairnum==pairNum()) at L37 CATCHES this under -ea (AssertionError on that read);
		//under -da it silently returns the wrong pairnum to every consumer (StackSites/SortSites/etc). The chrom/strand packing
		//below has the same -0 flaw for chrom==0, but chrom is 1-based so that leg is unreachable. Fix needs a real spare bit
		//for pairnum rather than the sign of a value whose domain includes 0.
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
