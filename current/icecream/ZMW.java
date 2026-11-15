package icecream;

import java.util.ArrayList;

import shared.Tools;
import stream.Read;
import stream.SamLine;
import structures.IntList;

/**
 * Mutable, ordered list of subread references for one PacBio well.
 * Retains the supplied Read objects; discard and header helpers mutate those
 * objects rather than copying them. Most helpers require non-null reads, and
 * header parsing requires compatible IDs. Inherited list mutations and these
 * helpers are unsynchronized; a consumer owns its delivered well.
 * Internal-read statistics exclude the first/last reads as potentially partial
 * passes, while length arrays and longest-read selection can include them.
 * @author Brian Bushnell, Sucrose
 * @date June 5, 2020
 */
public class ZMW extends ArrayList<Read>{
	
	/** For serialization compatibility */
	private static final long serialVersionUID = -2580124131008824113L;

	/** Creates an empty well without parsing or assigning a header identifier. */
	public ZMW(){super();}
	
	/** Creates an empty well with the requested list capacity.
	 * @param initialSize Initial ArrayList capacity
	 */
	public ZMW(int initialSize){super(initialSize);}

	/** @return Total current sequence length, including discarded reads; zero if empty. */
	public long countBases(){
		long x=0;
		for(Read r : this){
			x+=r.length();
		}
		return x;
	}
	
	/** Selects the upper median among eligible internal-read lengths.
	 * @param includeDiscarded Whether discarded internal reads contribute
	 * @return Median length, -1 for fewer than three reads, or zero for an empty eligible middle
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
		if(lengths.size==0){return 0;}
		lengths.sort();
		int median=lengths.get(lengths.size/2);
		return median;
	}
	
	/** @param includeDiscarded Whether discarded reads contribute
	 * @return Maximum eligible length over the whole well, or zero when none qualify
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
	
	/** Selects the first eligible internal read with the median length.
	 * Wells with fewer than three reads or no eligible interior use longestRead
	 * instead. Returned reads are references, without cloning or reordering.
	 * @param includeDiscarded Whether discarded reads are eligible
	 * @return Selected read, or null only when no eligible read exists anywhere in the well
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
		//PIC-103: A discarded interior must not hide surviving endpoints from the CCS caller.
		return longestRead(includeDiscarded);
	}
	
	/** @param includeDiscarded Whether discarded reads are eligible
	 * @return First longest eligible read, or null when none qualify
	 */
	public Read longestRead(boolean includeDiscarded){
		Read max=null;
		for(Read r : this){
			if((includeDiscarded || !r.discarded()) && (max==null || r.length()>max.length())){max=r;}
		}
		return max;
	}
	
	/** Parses the current first header each time, including after inherited list mutations.
	 * @return Current first read's numeric well ID, or -1 when empty
	 */
	public int zid(){
		//PIC-106: Do not cache across arbitrary list/header mutations inherited from ArrayList.
		return parseZID();
	}
	
	/** Parses without populating the reserved cache field.
	 * @return First read's numeric header well ID, or -1 if empty
	 */
	private int parseZID(){
		return (size()<1 ? -1 : PBHeader.parseZMW(get(0).id));
	}
	
	/** Adjusts a coordinate-like suffix after sequence trimming has already occurred.
	 * Clamps negative requested trims to zero, uses the last slash, preserves a
	 * tab/space description, and reconciles the stored interval length to r.length().
	 * Updates SAM qname and existing qs/qe tags; sequence/quality and other alignment
	 * metadata are the caller's responsibility. Noncoordinate suffixes are unchanged.
	 * @param r Non-null mutable read with a non-null header
	 * @param leftTrim Requested bases removed from the left
	 * @param rightTrim Requested bases removed from the right
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
					comment=suffix.substring(tab);
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
	
	/** Changes every contained read's discard flag without removing list entries.
	 * @param b New flag value
	 */
	public void setDiscarded(boolean b){
		for(Read r : this){
			r.setDiscarded(b);
		}
	}

	/** @return New ordered length array; null entries are represented by -1. */
	public int[] lengths(){
		final int size=size();
		int[] array=new int[size];
		for(int i=0; i<size; i++){
			Read r=get(i);
			array[i]=r==null ? -1 : r.length();
		}
		return array;
	}
	
	/** Estimates passes from well size and endpoint lengths, ignoring discard flags.
	 * Returns0,0.25,0.5 for zero, one, or two reads. Larger wells count internal
	 * reads as full passes and soften each endpoint's length/median ratio to at most
	 * 0.99. A nonpositive internal median is invalid for endpoint length ratios.
	 * @return Heuristic pass estimate
	 * @throws IllegalStateException If a well with at least three reads has a nonpositive internal median
	 */
	public float estimatePasses(){
		final int size=size();
		if(size<1){return 0;}
		else if(size==1){return 0.25f;}
		else if(size==2){return 0.5f;}
		
		int median=medianLength(true);
		if(median<=0){
			throw new IllegalStateException("Cannot estimate passes with internal median length "+median+" for "+size+" subreads");
		}
		int first=first().length();
		int last=last().length();

		return size-2+estimatePasses(first, median)+estimatePasses(last, median);
	}
	
	/** @param len Endpoint sequence length
	 * @param median Positive internal median length
	 * @return Softened endpoint contribution, capped at0.99
	 */
	private float estimatePasses(int len, int median){
		float ratio=len/(float)median;
		//TODO: I want this to be more asymptotic
		return Tools.min(0.99f, ratio/(1+0.05f*ratio));
	}

	/** @return Whether every read is discarded; also true for an empty well. */
	public boolean discarded(){
		for(Read r : this){
			if(!r.discarded()){return false;}
		}
		return true;
	}
	
	/** Stream-order ordinal when assigned by ZMWStreamer, distinct from the header well ID. */
	public long id;
	
	/** Legacy unused cache field retained in the serialized representation. */
	private int zid=-1;

	/** @return First read; requires a nonempty list. */
	public Read first(){return get(0);}
	/** @return Last read; requires a nonempty list. */
	public Read last(){return get(size()-1);}
	
}
