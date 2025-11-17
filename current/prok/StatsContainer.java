package prok;

import java.util.Arrays;

import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import structures.LongHashSet;

/**
 * Container for collecting and managing statistics for prokaryotic gene elements.
 * Maintains frame statistics for inner regions, start sites, and stop sites
 * along with length distribution data for statistical analysis.
 *
 * @author Brian Bushnell
 * @date 2025
 */
class StatsContainer {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Constructs a StatsContainer with specified parameters for all frame statistics.
	 * Initializes inner, start, and stop statistics with their respective k-mer lengths,
	 * frame counts, and offsets.
	 *
	 * @param type_ Element type identifier from ProkObject
	 * @param kInner K-mer length for inner region statistics
	 * @param framesInner Number of reading frames for inner regions
	 * @param kStart K-mer length for start site statistics
	 * @param framesStart Number of reading frames for start sites
	 * @param offsetStart Position offset for start site analysis
	 * @param kStop K-mer length for stop site statistics
	 * @param framesStop Number of reading frames for stop sites
	 * @param offsetStop Position offset for stop site analysis
	 */
	StatsContainer(int type_, int kInner, int framesInner, int kStart, int framesStart, int offsetStart, int kStop, int framesStop, int offsetStop){
		type=type_;
		name=ProkObject.typeStrings[type];
		setInner(kInner, framesInner);
		setStart(kStart, framesStart, offsetStart);
		setStop(kStop, framesStop, offsetStop);
	}

	/**
	 * Constructs a basic StatsContainer with only type specification.
	 * Frame statistics must be set separately using setter methods.
	 * @param type_ Element type identifier from ProkObject
	 */
	StatsContainer(int type_){
		type=type_;
		name=ProkObject.typeStrings[type];
	}
	
	/**
	 * Initializes inner region frame statistics with specified parameters.
	 * @param kInner K-mer length for inner region analysis
	 * @param framesInner Number of reading frames to analyze
	 */
	void setInner(int kInner, int framesInner){
		assert(inner==null);
		statsArray[0]=inner=new FrameStats(name+" inner", kInner, framesInner, 0);
	}
	
	/**
	 * Initializes start site frame statistics with specified parameters.
	 * @param kStart K-mer length for start site analysis
	 * @param framesStart Number of reading frames to analyze
	 * @param offsetStart Position offset for start site detection
	 */
	void setStart(int kStart, int framesStart, int offsetStart){
		assert(start==null);
		statsArray[1]=start=new FrameStats(name+" start", kStart, framesStart, offsetStart);
	}
	
	/**
	 * Initializes stop site frame statistics with specified parameters.
	 * @param kStop K-mer length for stop site analysis
	 * @param framesStop Number of reading frames to analyze
	 * @param offsetStop Position offset for stop site detection
	 */
	void setStop(int kStop, int framesStop, int offsetStop){
		assert(stop==null);
		statsArray[2]=stop=new FrameStats(name+" stop", kStop, framesStop, offsetStop);
	}
	
	/** Sets the inner region statistics using an existing FrameStats object.
	 * @param fs Pre-configured FrameStats for inner regions */
	void setInner(FrameStats fs){
		assert(inner==null);
		assert(fs!=null);
		statsArray[0]=inner=fs;
	}
	
	/** Sets the start site statistics using an existing FrameStats object.
	 * @param fs Pre-configured FrameStats for start sites */
	void setStart(FrameStats fs){
		assert(start==null);
		statsArray[1]=start=fs;
	}
	
	/** Sets the stop site statistics using an existing FrameStats object.
	 * @param fs Pre-configured FrameStats for stop sites */
	void setStop(FrameStats fs){
		assert(stop==null);
		statsArray[2]=stop=fs;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public String toString(){
		return appendTo(new ByteBuilder()).toString();
	}
	
	/**
	 * Appends formatted statistics data to a ByteBuilder.
	 * Includes header information, counts, and detailed frame statistics
	 * in tab-delimited format suitable for output files.
	 *
	 * @param bb ByteBuilder to append data to
	 * @return The same ByteBuilder with statistics appended
	 */
	public ByteBuilder appendTo(ByteBuilder bb) {
		bb.append("#name\t").append(name).nl();
		bb.append("#type\t").append(type).nl();
		bb.append("#count\t").append(lengthCount).nl();
		bb.append("#lengthSum\t").append(lengthSum).nl();
		//lengths
		bb.append("#contains\t").append(3).nl();
		if(ProkObject.processType(type)) {
			for(FrameStats fs : statsArray){
				fs.appendTo(bb);
			}
		}else{
			for(FrameStats fs : statsArray){
				fs.append0(bb);
			}
		}
		return bb;
	}
	
	/**
	 * Resets all statistics to initial empty state.
	 * Clears frame statistics, length arrays, and counters,
	 * then recalculates derived values.
	 */
	public void clear(){
		for(int i=0; i<statsArray.length; i++){
			if(statsArray[i]!=null){
				statsArray[i].clear();
			}
		}

		assert(inner==statsArray[0]);
		assert(start==statsArray[1]);
		assert(stop==statsArray[2]);
		
		Arrays.fill(lengths, 0);
		lengthSum=0;
		lengthCount=0;
		calculate();
	}
	
	/**
	 * Copies all statistics from another StatsContainer.
	 * Replaces current statistics with those from the source container,
	 * maintaining compatibility for same element types.
	 * @param sc Source StatsContainer to copy from
	 */
	public void setFrom(StatsContainer sc){
		assert(sc.name.equals(name));
		for(int i=0; i<statsArray.length; i++){
			FrameStats fs=sc.statsArray[i];
			if(statsArray[i]==null){
				statsArray[i]=new FrameStats(fs.name, fs.k, fs.frames, fs.leftOffset);
				statsArray[i].add(fs);
			}else{
				statsArray[i].setFrom(fs);
			}
		}
		inner=statsArray[0];
		start=statsArray[1];
		stop=statsArray[2];
		
		for(int i=0; i<lengths.length; i++){lengths[i]=sc.lengths[i];}
		lengthSum=sc.lengthSum;
		lengthCount=sc.lengthCount;
		calculate();
	}
	
	/**
	 * Merges statistics from another StatsContainer into this one.
	 * Adds frame statistics, length counts, and recalculates derived values
	 * to combine data from multiple sources.
	 * @param sc StatsContainer to merge with this one
	 */
	public void add(StatsContainer sc){
		assert(sc.name.equals(name));
		for(int i=0; i<statsArray.length; i++){
			FrameStats fs=sc.statsArray[i];
			if(statsArray[i]==null){
				statsArray[i]=new FrameStats(fs.name, fs.k, fs.frames, fs.leftOffset);
			}
			statsArray[i].add(fs);
		}

		inner=statsArray[0];
		start=statsArray[1];
		stop=statsArray[2];
		
		Tools.add(lengths, sc.lengths);
		lengthSum+=sc.lengthSum;
		lengthCount+=sc.lengthCount;
		calculate();
	}
	
	/**
	 * Scales all statistics by a multiplication factor.
	 * Applies the multiplier to frame statistics, length arrays, and counts
	 * for proportional scaling operations.
	 * @param mult Multiplication factor to apply to all statistics
	 */
	public void multiplyBy(double mult) {
		for(int i=0; i<statsArray.length; i++){
			FrameStats fs=statsArray[i];
			fs.multiplyBy(mult);
		}
		
		Tools.multiplyBy(lengths, mult);
		lengthSum=Math.round(lengthSum*mult);
		lengthCount=Math.round(lengthCount*mult);
		calculate();
	}
	
	/**
	 * Recalculates derived statistics from current data.
	 * Updates frame statistics calculations and computes average length
	 * and inverse average length from current counts.
	 */
	public void calculate(){
		for(int i=0; i<statsArray.length; i++){
			statsArray[i].calculate();
		}
		lengthAvg=(int)(lengthSum/Tools.max(1.0, lengthCount));
		invLengthAvg=1f/Tools.max(1, lengthAvg);
	}
	
	/**
	 * Records a length observation in the statistics.
	 * Updates length sum, count, and distribution histogram
	 * for length-based analysis.
	 * @param x Length value to record
	 */
	public void addLength(int x){
		lengthSum+=x;
		lengthCount++;
		lengths[Tools.min(x, lengths.length-1)]++;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Frame statistics for inner gene regions */
	FrameStats inner;
	/** Frame statistics for start sites */
	FrameStats start;
	/** Frame statistics for stop sites */
	FrameStats stop;
	/** Array containing inner, start, and stop frame statistics */
	final FrameStats[] statsArray=new FrameStats[3];
	
	/** Gets the long k-mer length for this element type */
	int kLongLen(){return ProkObject.kLongLen(type);}
	/** Gets the k-mer set for this element type */
	LongHashSet kmerSet(){return ProkObject.kmerSet(type);}
	/** Gets the consensus sequence reads for this element type */
	Read[] consensusSequence(){return ProkObject.consensusReads(type);}
	/** Gets the minimum identity threshold for this element type */
	float minIdentity(){return ProkObject.minID(type);}

	/** Gets the start position tolerance for this element type */
	public int startSlop(){return ProkObject.startSlop(type);}
	/** Gets the stop position tolerance for this element type */
	public int stopSlop(){return ProkObject.stopSlop(type);}
	
	/** Human-readable name for this element type */
	final String name;
	/** Total sum of all recorded lengths */
	long lengthSum=0;
	/** Number of length observations recorded */
	long lengthCount=0;
	/** Average length calculated from sum and count */
	int lengthAvg=-1;
	/** Inverse of average length for efficient calculations */
	float invLengthAvg;
	
	/** Histogram of length frequencies up to maximum array size */
	int[] lengths=new int[5000];
	/** Element type identifier corresponding to ProkObject types */
	public final int type;
	
}
