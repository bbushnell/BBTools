package prok;

import java.util.Arrays;

import dna.AminoAcid;
import shared.KillSwitch;
import shared.Parse;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Stores frame-relative kmer counts for a type of genomic feature, such as a coding start site.
 * @author Brian Bushnell
 * @date Sep 24, 2018
 */
public class FrameStats {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructs FrameStats for analyzing k-mer patterns around genomic features.
	 * Initializes count matrices for tracking k-mer frequencies in multiple frames.
	 *
	 * @param name_ Identifier for this statistics object
	 * @param k_ Length of k-mers to analyze
	 * @param frames_ Number of frames to analyze relative to reference points
	 * @param leftOffset_ Offset to the left of the reference point where analysis begins
	 */
	public FrameStats(String name_, int k_, int frames_, int leftOffset_){
		name=name_;
		k=k_;
		mask=~((-1)<<(2*k));
		frames=frames_;
		kMax=1<<(2*k);
		invFrames=1.0f/frames;
		leftOffset=leftOffset_;
		
		probs=new float[frames][kMax];
		countsTrue=new long[frames][kMax];
		countsFalse=new long[frames][kMax];
		counts=new long[][][] {countsFalse, countsTrue};
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Increments the count for a specific k-mer in a given frame.
	 * @param kmer The k-mer value (encoded as integer)
	 * @param frame Frame index relative to reference point
	 * @param valid 1 for valid features, 0 for invalid/background
	 */
	public void add(int kmer, int frame, int valid){
		counts[valid][frame][kmer]++;
		validSums[valid]++;
	}
	
	/**
	 * Checks if this FrameStats is compatible with another for merging operations.
	 * Compatibility requires matching name, leftOffset, k-mer length, and frame count.
	 * @param fs FrameStats to compare with
	 * @return true if compatible for merging, false otherwise
	 */
	public boolean compatibleWith(FrameStats fs) {
		return fs.name.equals(name) && fs.leftOffset==leftOffset && fs.k==k && fs.frames==frames;
	}
	
	/** Resets all count matrices and sums to zero.
	 * Used to reinitialize the statistics object for new data. */
	public void clear() {
		Tools.fill(counts, 0);
		Arrays.fill(validSums, 0);
	}
	
	/**
	 * Clears current data and copies all counts from another FrameStats.
	 * Both objects must be compatible as verified by compatibleWith().
	 * @param fs Source FrameStats to copy data from
	 */
	public void setFrom(FrameStats fs) {
		assert(compatibleWith(fs)) : name+", "+frames+", "+fs.name+", "+fs.frames;
		clear();
		add(fs);
	}
	
	/**
	 * Merges counts from another FrameStats into this one.
	 * Adds corresponding count values across all frames and validity states.
	 * Both objects must have matching parameters.
	 * @param fs FrameStats to merge with this one
	 */
	public void add(FrameStats fs){
		assert(fs.name.equals(name));
		assert(fs.leftOffset==leftOffset);
		assert(fs.k==k);
		assert(fs.frames==frames) : name+", "+frames+", "+fs.name+", "+fs.frames;
//		for(int x=0; x<counts.length; x++) {
//			for(int y=0; y<counts[x].length; y++) {
//				for(int z=0; z<counts[x][y].length; z++) {
//					assert(fs.counts[x][y][z]>=0) : counts[x][y][z]+", "+fs.counts[x][y][z]+", "+fs.name;
//					assert(counts[x][y][z]>=0) : counts[x][y][z]+", "+fs.counts[x][y][z];
//				}
//			}
//		}

		Tools.add(counts, fs.counts);
		Tools.add(validSums, fs.validSums);
//		for(int x=0; x<counts.length; x++) {
//			for(int y=0; y<counts[x].length; y++) {
//				for(int z=0; z<counts[x][y].length; z++) {
//					assert(fs.counts[x][y][z]>=0) : counts[x][y][z]+", "+fs.counts[x][y][z];
//					assert(counts[x][y][z]>=0) : counts[x][y][z]+", "+fs.counts[x][y][z];
//				}
//			}
//		}
//		calculate();
	}
	
	/**
	 * Scales all counts and sums by a multiplication factor.
	 * Used for normalizing or weighting statistics.
	 * @param mult Multiplication factor to apply to all counts
	 */
	public void multiplyBy(double mult) {
		Tools.multiplyBy(counts, mult);
		Tools.multiplyBy(validSums, mult);
	}
	
	/**
	 * Calculates probability scores from accumulated k-mer counts.
	 * Computes the overall average probability and frame-specific k-mer probabilities
	 * normalized by the inverse average for scoring purposes.
	 */
	void calculate(){
		average=(float)((validSums[1]+1.0)/(validSums[0]+validSums[1]+1.0));
		invAvg=1.0f/average;
		
		for(int a=0; a<frames; a++){
			for(int b=0; b<kMax; b++){
				long t=countsTrue[a][b];
				long f=countsFalse[a][b];
				probs[a][b]=(float)(t/(t+f+1.0))*invAvg;
			}
		}
	}
	
	/**
	 * Calculates a statistical score for a genomic position based on k-mer patterns.
	 * Analyzes k-mers in multiple frames around the specified point and computes
	 * a cumulative score based on learned probability patterns.
	 *
	 * @param point Reference position in the sequence to score
	 * @param bases DNA sequence as byte array
	 * @return Statistical score indicating likelihood of feature at this position
	 */
	public float scorePoint(int point, byte[] bases){
		final int mask=~((-1)<<(2*k));
		
		int kmer=0;
		int len=0;
		float score=0;

//		outstream.println("k="+k);
//		outstream.println("mask="+mask);
		
		int start=point-leftOffset;
		for(int i=start, frame=0-k+1; i<bases.length && frame<frames; i++, frame++){
			byte b=(i>=0 ? bases[i] : (byte)'A');
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			
//			outstream.println("b="+(char)b+", kmer="+kmer+", len="+(len+1)+", frame="+frame);
			
			if(x>=0){
				len++;
				if(len>=k){
					float prob=probs[frame][kmer];
					float dif=prob-0.99f;
					score+=dif;
					
//					if(name.equals("startStats")){
//						System.err.println("frame="+frame+" kmer="+AminoAcid.kmerToString(kmer, k)+
//								Tools.format(" prob=%.4f\tdif=%.4f\tscore=%.4f", prob, dif, score)+
//								"\tvalid="+counts[1][frame][kmer]+"\tinvalid="+counts[0][frame][kmer]);
//					}
				}
			}else{len=0;}
		}
		
		return score*invFrames;
	}
	
	/**
	 * Processes a sequence to accumulate k-mer counts for CDS frame analysis.
	 * Uses a validFrames array to indicate which frames contain valid coding sequences.
	 * Only processes when ProkObject.callCDS is enabled.
	 *
	 * @param bases DNA sequence to process
	 * @param validFrames Bit array indicating valid coding frames at each position
	 */
	void processCDSFrames(byte[] bases, byte[] validFrames){
		if(!ProkObject.callCDS){return;}
		int kmer=0;
		int len=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x>=0){
				len++;
				if(len>=k){
					int vf=validFrames[i];
					for(int frame=0; frame<frames; frame++){
						int valid=vf&1;
						add(kmer, frame, valid);
						//For CDS start (0-based) of 189, i=192, j=189, vf=1, frame=0 - all as expected.
//						assert(valid==0) : "vf="+vf+", frame="+frame+", len="+len+", kmer="+AminoAcid.kmerToString(kmer, k)+", i="+i+", j="+j;
						vf=(vf>>1);
					}
				}
			}else{len=0;}
		}
	}
	
	/**
	 * Processes k-mers around a specific genomic feature point.
	 * Accumulates counts across multiple frames relative to the specified point.
	 * Skips positions too close to sequence ends to avoid truncated data.
	 *
	 * @param bases DNA sequence containing the point
	 * @param point Position of the genomic feature (0-based)
	 * @param valid 1 if this is a true feature, 0 for background/negative example
	 */
	void processPoint(byte[] bases, int point, int valid){
		
		//Degenerate cases where the point is at the end, possibly from a truncated gene.
		if(point<3){return;}
		if(point>=bases.length-3){return;}
		
		int kmer=0;
		int len=0;

//		outstream.println("k="+k);
//		outstream.println("mask="+mask);
		
		int start=point-leftOffset;
		
		int i=start, frame=0-k+1;
		while(i<0){i++; frame++;}
		for(; i<bases.length && frame<frames; i++, frame++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			
//			outstream.println("b="+(char)b+", kmer="+kmer+", len="+(len+1)+", frame="+frame);
			
			if(x>=0){
				len++;
				if(len>=k){
					add(kmer, frame, valid);
				}
			}else{len=0;}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Text Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses a tab-delimited line of count data to populate the count matrices.
	 * Expected format: valid\tframe\tkmer_count_0\tkmer_count_1\t...
	 * Updates both the count matrix and validity sums.
	 * @param line Tab-delimited byte array containing count data
	 */
	public void parseData(byte[] line) {
		int a=0, b=0;
		final int valid, frame;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		valid=Parse.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(line);
		frame=Parse.parseInt(line, a, b);
		b++;
		a=b;
		
		assert(valid==0 || valid==1);
		assert(frame>=0 && frame<frames);
		try {
			final long[] row=counts[valid][frame];
			long sum=0;
			for(int kmer=0; kmer<row.length; kmer++){
				while(b<line.length && line[b]!='\t'){b++;}
				assert(b>a) : "Missing field 1: "+new String(line);
				long count=Parse.parseLong(line, a, b);
				b++;
				a=b;
				row[kmer]=count;
				sum+=count;
			}
			validSums[valid]+=sum;
		} catch (Exception e) {
			System.err.println(new String(line)+"\n"+name);
			assert(false) : e;
		}
	}
	
	@Override
	public String toString(){
		return appendTo(new ByteBuilder()).toString();
	}
	
	/**
	 * Appends formatted count data to a ByteBuilder for output.
	 * Includes header information and full count matrix in tab-delimited format.
	 * @param bb ByteBuilder to append data to
	 * @return The same ByteBuilder for method chaining
	 */
	public ByteBuilder appendTo(ByteBuilder bb){
		bb.append("#name\t").append(name).nl();
		bb.append("#k\t").append(k).nl();
		bb.append("#frames\t").append(frames).nl();
		bb.append("#offset\t").append(leftOffset).nl();
		bb.append("#valid\tframe");
		for(int i=0; i<kMax; i++){bb.tab().append(AminoAcid.kmerToString(i, k));}
		bb.nl();
		for(int a=0; a<2; a++){//valid
			for(int b=0; b<frames; b++){
				bb.append(a);
				bb.tab().append(b);
				for(int c=0; c<kMax; c++){
					bb.tab().append(counts[a][b][c]);
				}
				bb.nl();
			}
		}
		return bb;
	}
	
	/**
	 * Appends formatted template with zero counts to a ByteBuilder.
	 * Creates the same structure as appendTo() but with all counts set to zero.
	 * @param bb ByteBuilder to append template to
	 * @return The same ByteBuilder for method chaining
	 */
	public ByteBuilder append0(ByteBuilder bb){
		bb.append("#name\t").append(name).nl();
		bb.append("#k\t").append(k).nl();
		bb.append("#frames\t").append(frames).nl();
		bb.append("#offset\t").append(leftOffset).nl();
		bb.append("#valid\tframe");
		for(int i=0; i<kMax; i++){bb.tab().append(AminoAcid.kmerToString(i, k));}
		bb.nl();
		for(int a=0; a<2; a++){//valid
			for(int b=0; b<frames; b++){
				bb.append(a);
				bb.tab().append(b);
				for(int c=0; c<kMax; c++){
					bb.tab().append(0);
				}
				bb.nl();
			}
		}
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Identifier name for this FrameStats object */
	public final String name;
	/** Length of k-mers analyzed by this object */
	public final int k;
	/** Bit mask for k-mer encoding, computed as ~((-1)<<(2*k)) */
	public final int mask;
	/** Total number of frames analyzed relative to reference points */
	public final int frames;
	/** Maximum possible k-mer value, computed as 1<<(2*k) */
	public final int kMax;
	/** Inverse of frame count (1.0/frames) for score normalization */
	public final float invFrames;
	/** Number of positions to the left of reference point where analysis begins */
	public final int leftOffset;
	/**
	 * Calculates the right offset from the reference point based on total frames and left offset
	 */
	public int rightOffset() {return frames-leftOffset-1;}
	
	/** Probability matrix [frames][kmers] computed from count ratios */
	public final float[][] probs;
	/** Count matrix [frames][kmers] for valid/true feature examples */
	public final long[][] countsTrue;
	/** Count matrix [frames][kmers] for invalid/false background examples */
	public final long[][] countsFalse;
	/**
	 * Combined count matrix [valid][frames][kmers] containing both true and false counts
	 */
	public final long[][][] counts;

	/** Sum of all counts for each validity state [false_sum, true_sum] */
	public final long[] validSums=KillSwitch.allocLong1D(2);
	/** Overall average probability of valid features, computed from validSums */
	private float average=-1;
	/** Inverse of average probability (1.0/average) for score normalization */
	private float invAvg=-1;
	
}
