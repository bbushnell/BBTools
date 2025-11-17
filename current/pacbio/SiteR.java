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
	
	/**
	 * Creates a SiteR from a SiteScoreR object.
	 * Extracts genomic coordinates and metadata from the source object.
	 * @param ssr Source SiteScoreR containing alignment information
	 */
	public SiteR(SiteScoreR ssr){
		this(ssr.start, ssr.stop, ssr.chrom, ssr.strand, ssr.numericID, ssr.pairnum);
	}
	
	/**
	 * Creates a SiteR with genomic coordinates and read metadata.
	 * Packs strand information into chromosome field and pair number into numeric ID field
	 * for memory efficiency. Negative values indicate reverse strand and second-in-pair.
	 *
	 * @param start_ Start coordinate of the alignment
	 * @param stop_ Stop coordinate of the alignment
	 * @param chrom Chromosome identifier
	 * @param strand Strand orientation (0=forward, 1=reverse)
	 * @param numericID Unique identifier for the read
	 * @param pairnum Pair number (0=first, 1=second in paired reads)
	 */
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
	
	/**
	 * Compares this SiteR with a SiteScore for equality.
	 * Checks start, stop, chromosome, and strand coordinates.
	 * @param other SiteScore object to compare against
	 * @return true if genomic coordinates match, false otherwise
	 */
	public boolean equals(SiteScore other){
		if(other.start!=start){return false;}
		if(other.stop!=stop){return false;}
		if(other.chrom!=chrom()){return false;}
		if(other.strand!=strand()){return false;}
		return true;
	}
	
	/**
	 * Compares this SiteR with a SiteScoreR for equality.
	 * Checks start, stop, chromosome, and strand coordinates.
	 * @param other SiteScoreR object to compare against
	 * @return true if genomic coordinates match, false otherwise
	 */
	public boolean equals(SiteScoreR other){
		if(other.start!=start){return false;}
		if(other.stop!=stop){return false;}
		if(other.chrom!=chrom()){return false;}
		if(other.strand!=strand()){return false;}
		return true;
	}
	
	/**
	 * Recursively converts this SiteR and all linked sites to text representation.
	 * Traverses the linked list and appends each site's text to the StringBuilder.
	 * @param sb StringBuilder to append to (creates new one if null)
	 * @return StringBuilder containing text representation of all linked sites
	 */
	public StringBuilder toTextRecursive(StringBuilder sb){
		if(sb==null){sb=new StringBuilder();}else{sb.append(" ");}
		sb.append("("+toText()+")");
		if(next!=null){next.toTextRecursive(sb);}
		return sb;
	}
	
	/**
	 * Converts genomic coordinates to comma-separated text format.
	 * Format: start,stop,chromosome,strand,numericID,pairNumber
	 * @return StringBuilder containing comma-separated coordinate values
	 */
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
	
	/** Start coordinate of the genomic alignment */
	public final int start;
	/** Stop coordinate of the genomic alignment */
	public final int stop;
	/** Packed chromosome and strand information (negative for reverse strand) */
	public final int chromStrand;
	/** Packed numeric ID and pair number (negative for second-in-pair) */
	public final long idPairnum;
	/** Next SiteR in linked list for multiple alignment sites */
	public SiteR next;

	/** Extracts the numeric ID from the packed idPairnum field */
	public long numericID(){return idPairnum>=0 ? idPairnum : -idPairnum;}
	/**
	 * Extracts the pair number from the packed idPairnum field (0=first, 1=second)
	 */
	public int pairNum(){return idPairnum>=0 ? 0 : 1;}
	/** Extracts the chromosome identifier from the packed chromStrand field */
	public int chrom(){return chromStrand>=0 ? chromStrand : -chromStrand;}
	/**
	 * Extracts the strand orientation from the packed chromStrand field (0=forward, 1=reverse)
	 */
	public byte strand(){return chromStrand>=0 ? (byte)0 : (byte)1;};
	/**
	 * Calculates the length of the linked list starting from this SiteR.
	 * Traverses all next pointers to count total sites.
	 * @return Number of SiteR objects in the linked list
	 */
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
