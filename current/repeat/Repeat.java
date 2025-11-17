package repeat;

import java.util.Arrays;
import java.util.Comparator;

import dna.AminoAcid;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import structures.CRange;
import tracker.EntropyTracker;

/**
 * Represents a repetitive sequence region found in genomic data.
 * Tracks location, depth, gaps, and statistics for repeats detected by sliding window analysis.
 * Used for identifying tandem repeats, low-complexity regions, and other repetitive elements.
 *
 * @author Brian Bushnell
 * @date 2014
 */
public class Repeat implements Cloneable, Comparable<Repeat> {
	
	/**
	 * Constructs a new repeat region with specified parameters.
	 *
	 * @param contig_ The contig containing this repeat
	 * @param start_ Starting position of the repeat
	 * @param depth_ Repeat depth or coverage level
	 * @param k_ K-mer size or window size used for detection
	 * @param maxGap_ Maximum gap size allowed within the repeat
	 * @param minRepeat_ Minimum repeat length threshold
	 * @param type_ Type of repeat (E for entropy, R for repeat)
	 */
	public Repeat(Read contig_, int start_, int depth_, int k_, int maxGap_, int minRepeat_, char type_) {
		contig=contig_;
		contigNum=(contig==null ? -1 : contig.numericID);
		contigName=(contig==null ? null : contig.id);
		start=start_;
		depth=depth_;
		k=k_;
		maxGap=maxGap_;
		minRepeat=minRepeat_;
		type=type_;
	}
	
	/**
	 * Creates a deep copy of this repeat region.
	 * @return Cloned repeat with identical properties
	 * @throws RuntimeException if cloning fails
	 */
	public Repeat clone() {
		assert(start>=0) : this;
		try {
			return (Repeat)super.clone();
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException(e);
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/** Resets all fields to initial empty state.
	 * Clears position, statistics, and reference data. */
	public void clear() {
		contigNum=start=stop=-1;
		gapLen=gapBP=gapCount=0;
		
		depthSum=tipdist=0;
		minDepth=Integer.MAX_VALUE;
		maxDepth=-1;
		
		contig=null;
		contigName=null;
		gc=entropy=0;
	}
	
	/**
	 * Attempts to extend this repeat to include a new position.
	 * If the gap is too large or contig changes, finalizes current repeat and starts new one.
	 *
	 * @param currentContig The contig being processed
	 * @param pos Current position to potentially include
	 * @param currentDepth Depth at the current position
	 * @return Finalized repeat if gap exceeded, null if successfully extended
	 */
	public Repeat increment(Read currentContig, int pos, int currentDepth) {
		assert(currentDepth>=depth);//TODO: This can be disabled, but then the currentDepth<depth case needs to be handled; currently the function is not called inside gaps
		final int gap=pos-stop-1;
		if(contigNum==currentContig.numericID && gap<=maxGap) {//advance
//			System.err.println("A:"+this);
			stop=pos;
			gapLen+=gap;
			gapCount+=(gap>0 ? 1 : 0);
			gapBP+=(gap>=k ? gap-k+1 : 0);
			depthSum+=currentDepth;
			maxDepth=Tools.max(currentDepth, maxDepth);
			return null;
		}
		
		Repeat r=null;
		if(contigNum>=0 && length()>=minRepeat) {r=this.clone();}
//		System.err.println("B:"+r+", "+this);
		clear();
		contigNum=currentContig.numericID;
		start=pos-k+1;
		assert(start>=0) : start+", "+pos+", "+k;
		stop=pos;
		assert(stop<currentContig.bases.length) : stop+", "+currentContig.length();
		minDepth=maxDepth=depth;
		
		//These are not *strictly* needed and can use a lot of memory.
		//They could be cleared after calculating entropy and gc or removed entirely.
		contig=currentContig;
		contigName=currentContig.name();
		
		return r;
	}
	
	/** Calculates the length of this repeat region.
	 * @return Length in base pairs, or 0 if not initialized */
	public int length() {
		return start<0 ? 0 : stop-start+1;
	}
	
	/**
	 * Tests if this repeat overlaps with another repeat.
	 * @param r The repeat to test for overlap
	 * @return true if repeats overlap on same contig
	 */
	public final boolean overlaps(Repeat r) {
		return (contigNum==r.contigNum && start<=r.stop && stop>=r.start);
	}
	
	/**
	 * Tests if this repeat completely spans another repeat.
	 * @param r The repeat to test for spanning
	 * @return true if this repeat completely contains the other
	 */
	public final boolean spans(Repeat r) {
		return (contigNum==r.contigNum && start<=r.start && stop>=r.stop);
	}
	
	/**
	 * Tests if this repeat subsumes another repeat based on position and quality.
	 * A repeat subsumes another if it spans it, has equal or greater depth, and similar gap length.
	 *
	 * @param r The repeat to test for subsumption
	 * @param weak If true, allows subsumption even with more gaps
	 * @return true if this repeat subsumes the other
	 */
	public final boolean subsumes(Repeat r, boolean weak) {
		return spans(r) && depth>=r.depth && (gapLen<=r.gapLen || weak);
	}
	
	/**
	 * Calculates GC content, entropy, and tip distance statistics for this repeat.
	 * Uses the entropy tracker to compute sequence complexity.
	 * @param et Entropy tracker for complexity calculations
	 */
	void calcStats(EntropyTracker et) {
		et.clear();
		int[] acgtn=new int[5];
		byte[] bases=contig.bases;
		assert(stop<bases.length && start>=0) : start+", "+stop+", "+bases.length+"\n"+this;
		for(int i=start; i<=stop; i++) {
			byte b=bases[i];
			int num=AminoAcid.baseToNumber4[b];
			acgtn[num]++;
		}
		int atCount=acgtn[0]+acgtn[3];
		int gcCount=acgtn[1]+acgtn[2];
		gc=(gcCount)/(float)Tools.max(1, atCount+gcCount);
		entropy=et.averageEntropy(bases, false, start, stop);
		tipdist=Tools.min(start, bases.length-stop-1);
	}
	
	/** Sets the sequence preview for this repeat using the provided builder.
	 * @param bb ByteBuilder to use for sequence construction */
	void setSeq(ByteBuilder bb) {
		bb.clear();
		appendPreview(bb);
		seq=bb.toBytes();
	}
	
	/*--------------------------------------------------------------*/
	
	/** Converts this repeat to a coordinate range object.
	 * @return CRange representing the genomic coordinates of this repeat */
	public CRange toRange() {
		return new CRange(contigNum, start, stop, contig);
	}
	
	@Override
	public int compareTo(Repeat r) {
		if(depth!=r.depth) {return depth-r.depth;}
		int lenDif=length()-r.length();
		if(lenDif!=0) {return lenDif;}
		if(gapLen!=r.gapLen) {return r.gapLen-gapLen;}
		if(entropy!=r.entropy) {return entropy>r.entropy ? 1 : -1;}
		if(gc!=r.gc) {return gc<r.gc ? 1 : -1;}
		if(contigNum!=r.contigNum) {return r.contigNum>contigNum ? 1 : -1;}
		return r.start-start;
	}
	
	@Override
	public String toString() {
		return toBytes().toString();
	}
	
	/** Converts repeat to ByteBuilder format.
	 * @return ByteBuilder containing formatted repeat data */
	public ByteBuilder toBytes() {
		ByteBuilder bb=new ByteBuilder();
		return appendTo(bb);
	}
	
	/**
	 * Appends formatted repeat data to existing ByteBuilder.
	 * Includes depth, length, gap statistics, position, and sequence preview.
	 * @param bb ByteBuilder to append data to
	 * @return The same ByteBuilder for method chaining
	 */
	public ByteBuilder appendTo(ByteBuilder bb) {
		bb.append(depth).tab().append(length()).tab();
		bb.append(gapLen).tab().append(gapBP).tab().append(gapCount).tab();
		bb.append(contigNum).tab().append(start).tab().append(stop).tab().append(maxDepth).tab();
		bb.append(tipdist).tab().append(gc, 3).tab().append(entropy, 3).tab().append(contigName);
		
		if(SEQ_AFFIX_LEN>0) {
			bb.tab();
			appendPreview(bb);
		}
		return bb;
	}
	
	/**
	 * Appends sequence preview to ByteBuilder.
	 * Shows full sequence if short, or prefix/suffix with ellipsis if long.
	 * @param bb ByteBuilder to append sequence to
	 * @return The same ByteBuilder for method chaining
	 */
	ByteBuilder appendPreview(ByteBuilder bb) {
		if(seq!=null) {
			bb.append(seq);
		}else {
			byte[] bases=contig.bases;
			assert(stop<bases.length && start>=0) : start+", "+stop+", "+bases.length;
			int lim=SEQ_AFFIX_LEN;
			if(length()<=2*lim+3) {
				for(int i=start; i<=stop; i++) {bb.append(bases[i]);}
			}else {
				for(int i=0; i<lim; i++) {bb.append(bases[start+i]);}
				bb.append('.').append('.').append('.');
				for(int i=stop-lim+1; i<=stop; i++) {bb.append(bases[i]);}
			}
		}
		return bb;
	}
	
	/** Converts this repeat region to a Read object.
	 * @return Read containing the repeat sequence and metadata */
	public Read toRead(){
		Read r=new Read(fullSequence(), null, readHeader(), 0L);
		return r;
	}
	
	/** Extracts the complete sequence of this repeat region.
	 * @return Byte array containing the full repeat sequence */
	public byte[] fullSequence(){
		return Arrays.copyOfRange(contig.bases, start, stop+1);
	}
	
	/**
	 * Generates a descriptive header for this repeat when converted to a read.
	 * Includes position, depth, gap, and composition statistics.
	 * @return Formatted header string with repeat metadata
	 */
	public String readHeader() {
		ByteBuilder bb=new ByteBuilder();
		bb.append("seq").append(contigNum).under().append(start).append('-').append(stop);
		bb.under().append("depth").append(depth).under().append("max").append(maxDepth);
		bb.under().append("gaplen").append(gapLen);
		bb.under().append("gc").append(gc, 3).under().append("entropy").append(entropy, 3);
		bb.tab().append(contigName);
		return bb.toString();
	}
	
	/** Generates TSV header line for repeat output files.
	 * @return Tab-separated header describing repeat data columns */
	public static String tsvHeader() {
		ByteBuilder bb=new ByteBuilder();
		bb.append("#");
		bb.append("depth").tab().append("length").tab();
		bb.append("gapLen").tab().append("gapBP").tab().append("gaps").tab();
		bb.append("contig").tab().append("start").tab().append("stop").tab().append("maxDp").tab();
		bb.append("tipDist").tab().append("gc").tab().append("entropy").tab().append("cName");
		bb.tab().append("seq");
		return bb.toString();
	}
	
	/*--------------------------------------------------------------*/
	
	/** For depth subsumption */
	public static class PosComparator implements Comparator<Repeat>{

		@Override
		public int compare(Repeat a, Repeat b) {
			//Lower contigs first
			if(a.contigNum!=b.contigNum) {return a.contigNum>b.contigNum ? 1 : -1;}
			//Lower starts first
			if(a.start!=b.start) {return a.start-b.start;}
			//Lower stops first
			if(a.stop!=b.stop) {return a.stop-b.stop;}//Not really necessary since depth already does that 
			//Higher depth first
			return b.depth-a.depth;
		}
		
		/** Singleton instance of PosComparator for depth subsumption sorting */
		public static final PosComparator comparator=new PosComparator();
	}
	
	/** For positional subsumption */
	public static class PosComparator2 implements Comparator<Repeat>{

		@Override
		public int compare(Repeat a, Repeat b) {
			//Lower contigs first
			if(a.contigNum!=b.contigNum) {return a.contigNum>b.contigNum ? 1 : -1;}
			//Lower starts first
			if(a.start!=b.start) {return a.start-b.start;}
			//Higher stops first
			if(a.stop!=b.stop) {return b.stop-a.stop;}
			//Higher depth first
			return b.depth-a.depth;
		}
		
		/** Singleton instance of PosComparator2 for positional subsumption sorting */
		public static final PosComparator2 comparator=new PosComparator2();
	}
	
	/*--------------------------------------------------------------*/
	
	/** K-mer size or window size used for repeat detection */
	public final int k;//or window
	/** Maximum gap size allowed within a repeat region */
	public final int maxGap;//max gap allowed
	/** Minimum repeat length required for reporting */
	public final int minRepeat;//min repeat allowed
	/** Type of repeat: 'E' for entropy-based, 'R' for repeat-based */
	public final char type;//E for entropy, R for repeat
	
	/** Numeric identifier of the contig containing this repeat */
	public long contigNum=-1;
	/** Starting position of the repeat within the contig */
	public int start=-1;
	/** Ending position of the repeat within the contig */
	public int stop=-1;
	
	/** Number of gaps found within this repeat region */
	public int gapCount=0;
	/** Total length of all gaps within this repeat region */
	public int gapLen=0;
	/** Number of gap base pairs that contribute k-mers to the repeat */
	public int gapBP=0;
	/** Sum of all depth values across the repeat region */
	public long depthSum=0;
	/** Minimum depth observed within this repeat region */
	public int minDepth=Integer.MAX_VALUE;
	/** Maximum depth observed within this repeat region */
	public int maxDepth=-1;

	/** Reference to the contig containing this repeat */
	public Read contig;
	/** Name of the contig containing this repeat */
	public String contigName;
	/** Cached sequence data for this repeat region */
	public byte[] seq;
	/** GC content ratio of this repeat region */
	public float gc;
	/** Sequence entropy measure for this repeat region */
	public float entropy;
	/** Distance from nearest contig tip (start or end) */
	int tipdist;
	
	//TODO
//	public int minDepthInCurrentGap=Integer.MAX_VALUE;
//	public long depthSumInCurrentGap=0;
	
	/*--------------------------------------------------------------*/

	/** Depth threshold used to define this repeat region */
	public final int depth;
	
	/*--------------------------------------------------------------*/
	
	/** Number of bases to show at start and end of sequence previews */
	public static int SEQ_AFFIX_LEN=12;
	
}
