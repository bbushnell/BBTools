package icecream;

import java.util.ArrayList;

import shared.Tools;
import stream.Read;
import stream.SamLine;
import structures.IntList;

/**
 * Container for the list of reads from a single
 * PacBio ZMW.
 * @author Brian Bushnell
 * @date June 5, 2020
 */
public class ZMW extends ArrayList<Read> {
	
	/**
	 * For serialization.
	 */
	private static final long serialVersionUID = -2580124131008824113L;

	/** Creates an empty ZMW container */
	public ZMW(){super();}
	
	/** Creates a ZMW container with specified initial capacity
	 * @param initialSize Expected number of reads from this ZMW */
	public ZMW(int initialSize){super(initialSize);}

	/** Counts total bases across all reads in this ZMW
	 * @return Total number of bases in all reads */
	public long countBases(){
		long x=0;
		for(Read r : this){
			x+=r.length();
		}
		return x;
	}
	
	/**
	 * Calculates median read length excluding first and last reads.
	 * First and last reads are typically adapter/primer sequences.
	 * @param includeDiscarded Whether to include discarded reads in calculation
	 * @return Median length of interior reads, or -1 if fewer than 3 reads
	 */
	public int medianLength(boolean includeDiscarded){
		if(size()<3){return -1;}
		IntList lengths=new IntList(size()-2);
		
		for(int i=1; i<size()-1; i++){
			Read r=get(i);
			if(includeDiscarded || !r.discarded()){
				lengths.add(get(i).length());
			}
		}
		lengths.sort();
		int median=lengths.get(lengths.size/2);
		return median;
	}
	
	/**
	 * Finds the length of the longest read in this ZMW
	 * @param includeDiscarded Whether to consider discarded reads
	 * @return Length of longest read, or 0 if no valid reads
	 */
	public int longestLength(boolean includeDiscarded){
		int max=0;
		for(Read r : this){
			if(includeDiscarded || !r.discarded()){
				max=Tools.max(max, r.length());
			}
		}
		return max;
	}
	
	/**
	 * Returns a read with median length from interior reads.
	 * Falls back to longest read if median calculation fails.
	 * @param includeDiscarded Whether to consider discarded reads
	 * @return Read with median length, or null if no suitable read found
	 */
	public Read medianRead(boolean includeDiscarded){
		int len=medianLength(includeDiscarded);
		if(len<0){return longestRead(includeDiscarded);}
		for(int i=1; i<size()-1; i++){
			Read r=get(i);
			if((includeDiscarded || !r.discarded()) && r.length()==len){
				return r;
			}
		}
		return null;
	}
	
	/**
	 * Returns the longest read in this ZMW
	 * @param includeDiscarded Whether to consider discarded reads
	 * @return Longest read, or null if no valid reads
	 */
	public Read longestRead(boolean includeDiscarded){
		Read max=null;
		for(Read r : this){
			if((includeDiscarded || !r.discarded()) && (max==null || r.length()>max.length())){max=r;}
		}
		return max;
	}
	
	/** Gets the PacBio-assigned ZMW ID, parsing from first read if needed
	 * @return ZMW identifier from PacBio headers */
	public int zid(){
		if(zid==-1){parseZID();}
		return zid;
	}
	
	/** Extracts ZMW ID from the first read's header
	 * @return Parsed ZMW ID, or -1 if no reads or parsing fails */
	private int parseZID(){
		return (size()<1 ? -1 : PBHeader.parseZMW(get(0).id));
	}
	
	/**
	 * Updates read header coordinates after trimming operations.
	 * Modifies PacBio header format to reflect new start/stop positions
	 * after left and right trimming. Also updates SAM optional fields.
	 * @param r Read to modify
	 * @param leftTrim Number of bases trimmed from left end
	 * @param rightTrim Number of bases trimmed from right end
	 */
	public static void fixReadHeader(Read r, int leftTrim, int rightTrim){
		leftTrim=Tools.max(0, leftTrim);
		rightTrim=Tools.max(0, rightTrim);
		if(leftTrim<1 && rightTrim<1){return;}
		final int idx=r.id.lastIndexOf('/');
		if(idx>0 && idx<r.id.length()-3){
			String prefix=r.id.substring(0, idx+1);
			String suffix=r.id.substring(idx+1);
			if(suffix.indexOf('_')>0){
				String coords=suffix, comment="";
				int tab=suffix.indexOf('\t');
				if(tab<0){tab=suffix.indexOf(' ');}
				if(tab>0){
					coords=coords.substring(0, tab);
					comment=coords.substring(tab);
				}
				String[] split=Tools.underscorePattern.split(coords);
				int left=Integer.parseInt(split[0]);
				int right=Integer.parseInt(split[1]);
				left+=leftTrim;
				right-=rightTrim;
				if(left>right){left=right;}
				
				if(right-left!=r.length()){right=left+r.length();}
//				System.err.println(r.length()+", "+(right-left));
				
				r.id=prefix+left+"_"+right+comment;
				final SamLine sl=r.samline;
				if(sl!=null){
					sl.qname=r.id;
					if(sl.optional!=null){
						for(int i=0; i<sl.optional.size(); i++){
							String s=sl.optional.get(i);
							if(s.startsWith("qe:i:")){
								s="qe:i:"+right;
								sl.optional.set(i, s);
							}else if(s.startsWith("qs:i:")){
								s="qs:i:"+left;
								sl.optional.set(i, s);
							}
						}
					}
				}
			}
		}
	}
	
	/** Marks all reads in this ZMW as discarded or not discarded
	 * @param b True to discard all reads, false to keep all reads */
	public void setDiscarded(boolean b){
		for(Read r : this){
			r.setDiscarded(b);
		}
	}

	/**
	 * Returns array of read lengths in order
	 * @return Array where each element is length of corresponding read,
	 * or -1 for null reads
	 */
	public int[] lengths() {
		final int size=size();
		int[] array=new int[size];
		for(int i=0; i<size; i++){
			Read r=get(i);
			array[i]=r==null ? -1 : r.length();
		}
		return array;
	}
	
	/**
	 * Estimates total number of sequencing passes through the DNA molecule.
	 * Uses median read length as reference and calculates fractional passes
	 * for first/last reads based on their relative lengths.
	 * @return Estimated number of passes (interior reads + fractional ends)
	 */
	public float estimatePasses(){
		final int size=size();
		if(size<1){return 0;}
		else if(size==1){return 0.25f;}
		else if(size==2){return 0.5f;}
		
		int median=medianLength(true);
		int first=first().length();
		int last=last().length();

		return size-2+estimatePasses(first, median)+estimatePasses(last, median);
	}
	
	/**
	 * Estimates fractional passes for a read based on its length ratio
	 * to median. Uses asymptotic formula to prevent overestimation.
	 * @param len Length of read to estimate
	 * @param median Reference median length
	 * @return Fractional passes (0.0 to 0.99)
	 */
	private float estimatePasses(int len, int median){
		float ratio=len/(float)median;
		//TODO: I want this to be more asymptotic
		return Tools.min(0.99f, ratio/(1+0.05f*ratio));
	}

	/** Checks if all reads in this ZMW are discarded
	 * @return True if all reads are discarded, false if any read is kept */
	public boolean discarded() {
		for(Read r : this){
			if(!r.discarded()){return false;}
		}
		return true;
	}
	
	/** 
	 * Identifier assigned by streamer, not by PacBio.
	 * First identifier is 0, then 1, etc.
	 */
	public long id;
	
	/** 
	 * ZMW ID assigned by PacBio.
	 */
	private int zid=-1;

	/** Returns first read in this ZMW (typically adapter/primer) */
	public Read first(){return get(0);}
	/** Returns last read in this ZMW (typically adapter/primer) */
	public Read last(){return get(size()-1);}
	
}
