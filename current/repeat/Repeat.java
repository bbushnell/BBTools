package repeat;

import java.util.Arrays;
import java.util.Comparator;

import dna.AminoAcid;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import structures.CRange;
import tracker.EntropyTracker;

public class Repeat implements Cloneable, Comparable<Repeat> {
	
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
	
	public Repeat clone() {
		assert(start>=0) : this;
		try {
			return (Repeat)super.clone();
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException(e);
		}
	}
	
	/*--------------------------------------------------------------*/
	
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
	
	public int length() {
		return start<0 ? 0 : stop-start+1;
	}
	
	public final boolean overlaps(Repeat r) {
		return (contigNum==r.contigNum && start<=r.stop && stop>=r.start);
	}
	
	public final boolean spans(Repeat r) {
		return (contigNum==r.contigNum && start<=r.start && stop>=r.stop);
	}
	
	public final boolean subsumes(Repeat r, boolean weak) {
		return spans(r) && depth>=r.depth && (gapLen<=r.gapLen || weak);
	}
	
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
	
	void setSeq(ByteBuilder bb) {
		bb.clear();
		appendPreview(bb);
		seq=bb.toBytes();
	}
	
	/*--------------------------------------------------------------*/
	
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
	
	public ByteBuilder toBytes() {
		ByteBuilder bb=new ByteBuilder();
		return appendTo(bb);
	}
	
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
	
	public Read toRead(){
		Read r=new Read(fullSequence(), null, readHeader(), 0L);
		return r;
	}
	
	public byte[] fullSequence(){
		return Arrays.copyOfRange(contig.bases, start, stop+1);
	}
	
	public String readHeader() {
		ByteBuilder bb=new ByteBuilder();
		bb.append("seq").append(contigNum).under().append(start).append('-').append(stop);
		bb.under().append("depth").append(depth).under().append("max").append(maxDepth);
		bb.under().append("gaplen").append(gapLen);
		bb.under().append("gc").append(gc, 3).under().append("entropy").append(entropy, 3);
		bb.tab().append(contigName);
		return bb.toString();
	}
	
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
		
		public static final PosComparator2 comparator=new PosComparator2();
	}
	
	/*--------------------------------------------------------------*/
	
	public final int k;//or window
	public final int maxGap;//max gap allowed
	public final int minRepeat;//min repeat allowed
	public final char type;//E for entropy, R for repeat
	
	public long contigNum=-1;
	public int start=-1;
	public int stop=-1;
	
	public int gapCount=0;
	public int gapLen=0;
	public int gapBP=0;
	public long depthSum=0;
	public int minDepth=Integer.MAX_VALUE;
	public int maxDepth=-1;

	public Read contig;
	public String contigName;
	public byte[] seq;
	public float gc;
	public float entropy;
	int tipdist;
	
	//TODO
//	public int minDepthInCurrentGap=Integer.MAX_VALUE;
//	public long depthSumInCurrentGap=0;
	
	/*--------------------------------------------------------------*/

	public final int depth;
	
	/*--------------------------------------------------------------*/
	
	public static int SEQ_AFFIX_LEN=12;
	
}
