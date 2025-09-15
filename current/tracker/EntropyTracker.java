package tracker;

import java.util.Arrays;

import dna.AminoAcid;
import shared.Parse;
import shared.Tools;

/**
 * Tracks entropy over a sliding window.
 * @author Brian Bushnell
 * @date Oct 6, 2017
 *
 */
public class EntropyTracker {
	
	public static void main(String[] args){
		final int k=args.length>0 ? Integer.parseInt(args[0]) : 2;
		final int window=args.length>1 ? Integer.parseInt(args[1]) : 3;
		final float cutoff=args.length>2 ? Float.parseFloat(args[2]) : 0.7f;
		final boolean highPass=args.length>3 ? Parse.parseBoolean(args[3]) : true;
		
		EntropyTracker et=new EntropyTracker(k, window, false, cutoff, highPass);
		System.err.println(et);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Normal constructor.
	 * @param k_ Kmer length.
	 * @param window_ Window size in bases.
	 */
	public EntropyTracker(int k_, int window_, boolean amino_){
		this(k_, window_, amino_, -1, true);
	}
	
	/**
	 * Allows the use of passes() based on entropy.
	 * @param k_ Kmer length.
	 * @param window_ Window size in bases.
	 * @param cutoff_ Entropy cutoff, 0 (no entropy) to 1 (max entropy).
	 * @param highPass_ True passes entropy of at least cutoff; false fails.
	 */
	public EntropyTracker(boolean amino_, float cutoff_, boolean highPass_){
		this((setDefaultK ? defaultK : amino_ ? 2 : defaultK), 
				(setDefaultWindow ? defaultWindowBases : amino_ ? 25 : 50),
				amino_, cutoff_, highPass_);
	}
	
	/**
	 * Allows the use of passes() based on entropy.
	 * @param k_ Kmer length.
	 * @param window_ Window size in bases.
	 * @param cutoff_ Entropy cutoff, 0 (no entropy) to 1 (max entropy).
	 * @param highPass_ True passes entropy of at least cutoff; false fails.
	 */
	public EntropyTracker(int k_, int window_, boolean amino_, float cutoff_, boolean highPass_){
		
		k=k_;
		windowBases=window_;
		windowKmers=windowBases-k+1;
		amino=amino_;
		entropyCutoff=cutoff_;
		highPass=highPass_;
		
		
		assert(k>0 && k<=15 && k<windowBases) : k+", "+windowBases;
		assert(windowKmers>0 && (entropyCutoff>=0 || entropyCutoff==-1) && entropyCutoff<=1) : k+", "+windowBases+", "+windowKmers+", "+entropyCutoff;
		
		bitsPerBase=(amino ? 5 : 2);
		mask=(k>15 ? -1 : ~((-1)<<(bitsPerBase*k)));
		kmerSpace=(1<<(bitsPerBase*k));//Note: This should be different for amino, but decoding would be a pain
		symbolToNumber=AminoAcid.symbolToNumber(amino);
		symbolToNumber0=AminoAcid.symbolToNumber0(amino);
		
		entropy=makeEntropyArray(windowKmers);
		entropyMult=-1/Math.log(windowKmers);
		baseCountMult=1f/windowBases;
		baseCounts=new short[4];
		counts=new short[kmerSpace];
		countCounts=new short[windowKmers+2];
		countCounts[0]=(short)windowKmers;
		baseRingBuffer=new byte[windowBases];
		
//		assert(false) : "\namino="+amino+", windowBases="+windowBases+", windowKmers="+windowKmers+", cutoff="+entropyCutoff+", bpb="+bitsPerBase+", mask="+
//			Integer.toBinaryString(mask)+", space="+kmerSpace+", mult="+entropyMult+"\n"+"entropy="+Arrays.toString(entropy);
		
//		entropyDeltaPlus=makeEntropyDeltaPlus(entropy, entropyMult);
//		entropyDeltaMinus=makeEntropyDeltaMinus(entropy, entropyMult);
	}
	
	/**
	 * @param maxCount Highest possible count; equal to window size in kmers.
	 * @return Entropy array
	 */
	private static final double[] makeEntropyArray(int maxCount){
		final double[] array=new double[maxCount+2]; //+2 to handle temporary condition of an extra kmer
		final double mult=1d/maxCount;
		for(int i=1; i<array.length; i++){//First element (kmer count of 0) contains zero entropy, and yields NaN, so is skipped.
			double pk=i*mult;
			array[i]=pk*Math.log(pk);
		}
		return array;
	}
	
	/**
	 * @param entropy Entropy array
	 * @param entropyMult Multiplier applied to each array element
	 * @return entropyDeltaPlus array
	 */
	@SuppressWarnings("unused")
	private static final double[] makeEntropyDeltaPlus(double[] entropy, double entropyMult){
		final double[] array=new double[entropy.length];
		for(int i=0; i<entropy.length-1; i++){//Last element is never used.
			array[i]=(entropy[i+1]-entropy[i])*entropyMult;
		}
		return array;
	}
	
	/**
	 * @param entropy Entropy array
	 * @param entropyMult Multiplier applied to each array element
	 * @return entropyDeltaMinus array
	 */
	@SuppressWarnings("unused")
	private static final double[] makeEntropyDeltaMinus(double[] entropy, double entropyMult){
		final double[] array=new double[entropy.length];
		for(int i=1; i<entropy.length; i++){//First element is never used.
			array[i]=(entropy[i-1]-entropy[i])*entropyMult;
		}
		return array;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** @return Number of unique kmers in current window. */
	public int unique(){return unique;}
	
	/** @return Number of undefined bases in current window. */
	public int ns(){return ns;}
	
	/** @return Sequence position of rightmost base in window. */
	public int rightPos(){return len-1;}
	
	/** @return Sequence position of leftmost base in window. */
	public int leftPos(){return len-windowBases;}

	/** @return Window size in bases. */
	public int windowBases() {return windowBases;}

	public int k() {return k;}

	/** @return Entropy cutoff. */
	public float cutoff() {return entropyCutoff;}
	
	/*--------------------------------------------------------------*/
	/*----------------      Entropy Calculation     ----------------*/
	/*--------------------------------------------------------------*/
	
	public float calcMaxMonomerFraction(){
		//int total=sum(baseCounts); should be window length
		int max=Tools.max(baseCounts);
		return baseCountMult*max;
	}
	
	/**
	 * Calculate entropy in current window.
	 * @return Entropy in current window.
	 */
	public float calcEntropy(){
		if(speed==FAST){return calcEntropyFast();}
		else if(speed==MEDIUM){return calcEntropyMedium();}
		else if(speed==SLOW){return calcEntropySlow();}
		else{return calcEntropySuperSlow();}
	}
	
	/** Sub-method of calcEntropy() */
	private float calcEntropyFast(){
		final float f=(float)(currentEsum*entropyMult);
//		final float f=(float)currentEsum; //For delta arrays
		
		//Avoid potential negative numbers due to underflow
		return f>0 ? f : 0;
	}
	
	/** Sub-method of calcEntropy()
	 * Calculates entropy from countCounts using precalculated entropy array and early exit. */
	private float calcEntropyMedium(){
		//Sum of entropy contributions from each kmer
		double eSum=0;
		
		//Sum of kmer counts, used for loop early exit
		int cSum=countCounts[0];
		
		//Simpler loop with no early exit.  Slower for complex sequences, but faster for homopolymonomers.
//		for(int i=1; i<countCounts.length; i++){
		
		//Loop over all nonzero kmer counts
		for(int i=1; cSum<windowKmers; i++){
			//Number of unique kmers with count i
			final int cc=countCounts[i];
			cSum+=cc;
			
			//Entropy contribution for each unique kmer with count i
			double pklogpk=entropy[i];
			
			//Add the entropy contribution for all unique kmers with count i
			eSum+=(cc*pklogpk);
		}
		//eSum now holds negative entropy in bits.
		
		//Adjust entropy to 0-1 scale based on window size
		float e=(float)(eSum*entropyMult);
		assert(e>=0 && e<=1) : e+", "+eSum+", "+entropyMult+"\n"+Arrays.toString(entropy)+"\n"+this;
		
		//Get rid of negative zero
		if(e<=0){e=0;}
		
		if(verbose){
			System.err.println(Tools.format("%.3f", eSum)+"\t"+Tools.format("%.3f", entropyMult)+"\t"
					+Tools.format("%.3f", e)+"\t"+len+"\t"+unique+"\t"+basesToString());
		}
		
		return e;
	}
	
	/** Sub-method of calcEntropy()
	 * Calculates entropy from counts using precalculated entropy array. */
	private float calcEntropySlow(){
		//Sum of entropy contributions from each kmer
		double eSum=0;
		
		//Loop over all nonzero kmer counts
		for(int count : counts){
			
			//Entropy contribution for this kmer's count
			double pklogpk=entropy[count];
			eSum+=pklogpk;
		}
		//eSum now holds negative entropy in bits.
		
		//Adjust entropy to 0-1 scale based on window size
		float e=(float)(eSum*entropyMult);
		assert(e>=0 && e<=1) : e+", "+eSum+", "+entropyMult+"\n"+Arrays.toString(entropy)+"\n"+this;
		
		//Get rid of negative zero
		if(e<=0){e=0;}
		
		if(verbose){
			System.err.println(Tools.format("%.3f", eSum)+"\t"+Tools.format("%.3f", entropyMult)+"\t"
					+Tools.format("%.3f", e)+"\t"+len+"\t"+unique+"\t"+basesToString());
		}
		
		return e;
	}
	
	/**
	 * Sub-method of calcEntropy()
	 * Demonstrates explicit entropy calculation.
	 * 
	 * Definition:
	 * Entropy, or information content, can be calculated using kmer counts of a sequence.
	 * 
	 * Probability of an event (unique kmer), pk, is that kmer's count divided by the number of kmers.
	 * Entropy contribution from that kmer is -pk*log(pk).
	 * Total entropy is sum of (-pk*log(pk)) for all kmer counts.
	 * entropyMult is simply a multiplier to convert the entropy measure (in bits) to a convenient 0-1 scale,
	 * corresponding to the inverse of the maximum possible entropy.
	 */
	private float calcEntropySuperSlow(){
		//Sum of entropy contributions from each kmer
		double eSum=0;
		
		//Loop over all nonzero kmer counts
		for(int count : counts){
			//Prevent NaN and INF
//			System.err.println("count="+count);
			if(count>0){
				//Fraction of total represented by this kmer
				double pk=count/(double)windowKmers;
//				System.err.println("pk="+pk);

				//Entropy contribution for this kmer's count
				double npklogpk=-pk*Math.log(pk);
//				System.err.println("npklogpk="+npklogpk);
				eSum+=npklogpk;
//				System.err.println("eSum="+eSum);
			}
		}
		//eSum now holds entropy in bits.
		
		//Multiplier to convert entropy to 0-1 scale.
		double multiplier=1/Math.log(windowKmers);
		
		//Adjust entropy to 0-1 scale based on window size
		float e=(float)(eSum*multiplier);
//		System.err.println("e="+e);
		assert(e>=0 && e<=1) : e+", "+eSum+", "+entropyMult+"\n"+Arrays.toString(entropy)+"\n"+this;
		
		//Get rid of negative zero
		if(e<=0){e=0;}
		
		if(verbose){
			System.err.println(Tools.format("%.3f", eSum)+"\t"+Tools.format("%.3f", entropyMult)+"\t"
					+Tools.format("%.3f", e)+"\t"+len+"\t"+unique+"\t"+basesToString());
		}
		
		return e;
	}
	
	/*--------------------------------------------------------------*/
	
	public static float strandedness(byte[] bases, int[] counts, int k) {
		if(counts==null) {counts=new int[1<<(2*k)];}
		countKmers(bases, counts, k);
		return strandedness(counts, k);
	}
	
	public static float strandedness(int[] counts, int k) {
		final int mask=~((-1)<<(2*k));
		assert(mask==counts.length-1);
		long lower=0, upper=0;
		for(int kmer=0, limit=counts.length/2; kmer<limit; kmer++) {
			int a=counts[kmer];
			int b=counts[mask&(~kmer)];
			lower+=Math.min(a, b);
			upper+=Math.max(a, b);
		}
		return lower/(float)(Long.max(1, upper));
	}
	
	public static float strandedness(long[] counts, int k) {
		final int mask=~((-1)<<(2*k));
		assert(mask==counts.length-1);
		long lower=0, upper=0;
		for(int kmer=0, limit=counts.length/2; kmer<limit; kmer++) {
			long a=counts[kmer];
			long b=counts[mask&(~kmer)];
			lower+=Math.min(a, b);
			upper+=Math.max(a, b);
		}
		return lower/(float)(Long.max(1, upper));
	}
	
	public static float strandednessK2(int[] counts) {
		final int mask=15;
		assert(counts.length==16);
		long lower=0, upper=0;
		for(int kmer=0; kmer<8; kmer++) {
			int a=counts[kmer];
			int b=counts[mask&(~kmer)];
			lower+=Math.min(a, b);
			upper+=a+b;
		}
		return lower/(float)(Long.max(1, upper-lower));
	}
	
	public static float strandednessWindowed(byte[] bases, int[] counts, int k, int window) {
		if(k==2) {return strandednessWindowedK2(bases, counts, window);}
		assert(k>2);
		if(counts==null) {counts=new int[1<<(2*k)];}
		
		final int shift=2*k;
		final int mask=~((-1)<<shift);
		assert(mask==counts.length-1);

		int valid=0;
		double sum=0;
		int sums=0;
		for(int i=0, j=-window, ikmer=0, jkmer=0, ilen=0, jlen=0; i<bases.length; i++, j++){
			{
				byte b=bases[i];
				int x=AminoAcid.baseToNumber[b];
				ikmer=((ikmer<<2)|x)&mask;
				if(x>=0){
					ilen++;
					if(ilen>=k) {
						valid++;
						counts[ikmer]++;
					}
				}else{ilen=ikmer=0;}
			}
			
			if(j>=0){
				byte b=bases[i];
				int y=AminoAcid.baseToNumber[b];
				jkmer=((jkmer<<2)|y)&mask;
				if(y>=0){
					jlen++;
					if(jlen>=k) {
						valid--;
						counts[jkmer]--;
					}
				}else{jlen=jkmer=0;}
			}
			
			if(i>=window-1) {
				sums++;
				sum+=strandedness(counts, k);
			}
		}
		if(sums<1) {
			return strandedness(counts, k);
		}
		return (float)(sum/sums);
	}
	
	public static float strandednessWindowedK2(byte[] bases, int[] counts, int window) {
		if(counts==null) {counts=new int[16];}
		
		final int mask=15;
		assert(mask==counts.length-1);

		int valid=0;
		double sum=0;
		int sums=0;
		for(int i=0, j=-window, ikmer=0, jkmer=0, ilen=0, jlen=0; i<bases.length; i++, j++){
			{
				byte b=bases[i];
				int x=AminoAcid.baseToNumber[b];
				ikmer=((ikmer<<2)|x)&mask;
				if(x>=0){
					ilen++;
					if(ilen>=2) {
						valid++;
						counts[ikmer]++;
					}
				}else{ilen=ikmer=0;}
			}
			
			if(j>=0){
				byte b=bases[i];
				int y=AminoAcid.baseToNumber[b];
				jkmer=((jkmer<<2)|y)&mask;
				if(y>=0){
					jlen++;
					if(jlen>=2) {
						valid--;
						counts[jkmer]--;
					}
				}else{jlen=jkmer=0;}
			}
			
			if(i>=window-1) {
				sums++;
				sum+=strandednessK2(counts);
			}
		}
		if(sums<1) {
			return strandednessK2(counts);
		}
		return (float)(sum/sums);
	}
	
	/*--------------------------------------------------------------*/
	
	public static float calcEntropy(byte[] bases, int[] counts, int k){
		assert(k<=10) : k;//This is for small kmers
		if(counts==null) {counts=new int[1<<(2*k)];}
		countKmers(bases, counts, k);
		return calcEntropyFromCounts(counts);
	}
	
	/**
	 * Static, non-windowed version.
	 * Demonstrates explicit entropy calculation.
	 * 
	 * Definition:
	 * Entropy, or information content, can be calculated using kmer counts of a sequence.
	 * 
	 * Probability of an event (unique kmer), pk, is that kmer's count divided by the number of kmers.
	 * Entropy contribution from that kmer is -pk*log(pk).
	 * Total entropy is sum of (-pk*log(pk)) for all kmer counts.
	 * entropyMult is simply a multiplier to convert the entropy measure (in bits) to a convenient 0-1 scale,
	 * corresponding to the inverse of the maximum possible entropy.
	 */
	public static float calcEntropyFromCounts(int[] counts){
		//Sum of entropy contributions from each kmer
		double eSum=0;
		
		long windowKmers=Tools.sum(counts);
		double invKmers=1.0/windowKmers;
		//Loop over all nonzero kmer counts
		for(int count : counts){
			//Prevent NaN and INF
//			System.err.println("count="+count);
			if(count>0){
				//Fraction of total represented by this kmer
				double pk=count*invKmers;
//				System.err.println("pk="+pk);

				//Entropy contribution for this kmer's count
				double npklogpk=-pk*Math.log(pk);
//				System.err.println("npklogpk="+npklogpk);
				eSum+=npklogpk;
//				System.err.println("eSum="+eSum);
			}
		}
		//eSum now holds entropy in bits.
		
		//Multiplier to convert entropy to 0-1 scale.
		double multiplier=1/Math.log(windowKmers);
		
		//Adjust entropy to 0-1 scale based on window size
		float e=(float)(eSum*multiplier);
		
		//Get rid of negative zero
		if(e<=0){e=0;}
		
		return e;
	}
	
	public static int countKmers(final byte[] bases, final int[] counts, int k){
		Arrays.fill(counts, 0);
		if(bases==null || bases.length<k){return 0;}
		
		final int shift=2*k;
		final int mask=~((-1)<<shift);
		
		int kmer=0;
		int len=0;
		int valid=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x>=0){
				len++;
				if(len>=k) {
					valid++;
					counts[kmer]++;
				}
			}else{len=kmer=0;}
		}
		return valid;
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Calculate the average entropy of a sequence.
	 * @param bases Sequence as bytes
	 * @param allowNs True if windows containing undefined bases should be included
	 * @return Average entropy
	 */
	public float averageEntropy(byte[] bases, boolean allowNs){
		return averageEntropy(bases, allowNs, 0, bases.length-1);
	}
	
	/**
	 * Calculate the average entropy of a sequence.
	 * @param bases Sequence as bytes
	 * @param allowNs True if windows containing undefined bases should be included
	 * @return Average entropy
	 */
	public float averageEntropy(final byte[] bases, final boolean allowNs, final int from, final int to){
		assert(from>=0 && to<bases.length && from<=to) : from+", "+to+", "+bases.length;
		final int len=to-from+1;
		
		//Reset the initial state
		clear();
		
		//Position in sequence
		int i=from;
		
		//Accumulated entropy
		double sum=0;
		
		//Number of entropy measurements
		int divisor=0;
		
//		System.err.println("\n"+new String(bases, from, len));
//		System.err.println("from="+from+", to="+to+", lim=min("+(bases.length)+","+windowBases+","+(to+1)+")");
		
		//Prefill the first window
		for(final int lim=Tools.min(bases.length, windowBases+from, to+1); i<lim; i++){
			add(bases[i]);
		}
		
		//Calculate entropy for the first window.
		//This allows entropy to pass if it is high enough even though the sequence is shorter than window length.
		if(allowNs || ns==0){
			sum+=calcEntropy();
			divisor++;
		}
		
		//Calculate entropy for remaining windows
		for(; i<=to; i++){
			add(bases[i]);
			if(allowNs || ns==0){
				sum+=calcEntropy();
				divisor++;
			}
		}
		
		if(divisor<0){return -1;}//No valid windows.
		
		//Calculate the average
		double avg=(sum/(Tools.max(1, divisor)));
		return (float)avg;
	}
	
	/**
	 * Reports the longest block of consecutive bases in which all windows
	 * are below the entropy cutoff and at least the (optional) monomer fraction.
	 * @param bases
	 * @param allowNs
	 * @param maxMonomer
	 * @return
	 */
	public int longestLowEntropyBlock(byte[] bases, boolean allowNs, float maxMonomerFraction){
		//Reset the initial state
		clear();
		
		//Position in sequence
		int i=0;
		
		//Number of entropy measurements
		int totalWindows=0;
		
		double sum=0;
		int totalLowWindows=0;
		int currentLowWindows=0;
		int maxLowWindows=0;
		
		//Prefill the first window
		for(final int lim=Tools.min(bases.length, windowBases); i<lim; i++){add(bases[i]);}
		
		//Calculate entropy for the first window.
		//This allows entropy to pass if it is high enough even though the sequence is shorter than window length.
		if(allowNs || ns==0){
			totalWindows++;
			double e=calcEntropy();
			float mmf=calcMaxMonomerFraction();
			sum+=e;
			if(e<entropyCutoff && mmf>=maxMonomerFraction){
				totalLowWindows++;
				currentLowWindows++;
				maxLowWindows=Tools.max(maxLowWindows, currentLowWindows);
			}else{
				currentLowWindows=0;
			}
		}
		
		//Calculate entropy for remaining windows
		for(; i<bases.length; i++){
			add(bases[i]);
			if(allowNs || ns==0){
				totalWindows++;
				double e=calcEntropy();
				float mmf=calcMaxMonomerFraction();
				sum+=e;
				if(e<entropyCutoff && mmf>maxMonomerFraction){
					totalLowWindows++;
					currentLowWindows++;
					maxLowWindows=Tools.max(maxLowWindows, currentLowWindows);
				}else{
					currentLowWindows=0;
				}
			}
		}
		
		//Calculate the average; not needed
		double avg=(sum/(Tools.max(1, totalWindows)));
		
		int maxLowBlock=maxLowWindows<1 ? 0 : Tools.min(bases.length, maxLowWindows+windowBases-1);
		return maxLowBlock;
	}
	
	/**
	 * Calculate entropy in the window and compare to the cutoff.
	 * If Ns are important they should be handled externally with ns().
	 * @return True if the entropy passes the cutoff.
	 */
	public boolean passes(){
		//This function should only be used if entropyCutoff is set.
		assert(entropyCutoff>=0);
		float e=calcEntropy();
		
		//XOR: highPass inverts truth of comparison.
		return highPass ^ (e<entropyCutoff);
	}
	
	/**
	 * Calculate average entropy of the sequence and compare to the cutoff.
	 * @param sequence Sequence to measure.
	 * @param allowNs True if entropy should be calculated in windows containing Ns.
	 * @return True if the average entropy passes the cutoff.
	 */
	public boolean passes(byte[] sequence, boolean allowNs){
		//This function should only be used if entropyCutoff is set.
		assert(entropyCutoff>=0);
		float e=averageEntropy(sequence, allowNs);
		
		//XOR: highPass inverts truth of comparison.
		return highPass ^ (e<entropyCutoff);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Mutating Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Slide the window by adding a new base.
	 * @param b Base to add.
	 */
	public void add(final byte b){
		
		//Test initial state
		assert(!verify || verify()) : this;
		
		final byte oldBase=baseRingBuffer[pos]; //Leftmost base, about to be overwritten
		
		if(verbose){System.err.println("\nAdding "+Character.toString((char)b)+
				"; evicting "+Character.toString((char)oldBase)+"; counts="+Arrays.toString(counts)+"; countcounts="+Arrays.toString(countCounts)+", pos="+pos+", pos2="+pos2);}
		
		//Increment length
		len++;
		
		{//Add a new rightmost base
			baseRingBuffer[pos]=b;
			final int n=symbolToNumber0[b];
			baseCounts[n]++;
			kmer=((kmer<<bitsPerBase)|n)&mask; //Generate new rightmost kmer using the new base
			
			//Update number of Ns in current window
			if(!isFullyDefined(b)){
				ns++;
				assert(ns<=windowBases+1) : "There are more Ns than bases in the window:\n"+this;
			}
			
			if(len>=k){//Add a kmer
				//System.err.println("Adding "+kmer);
				
				final short oldCount=counts[kmer];
				
				/* Update unique kmer count */
				if(oldCount<1){
					assert(oldCount==0) : "An incoming array has negative counts: \n"+this;
					unique++;
				}
				
				/* Decrement the old countCount */
				countCounts[oldCount]--;
				
				/* countCounts[0] could be temporarily -1 at this point; for all others, min is 0. */
				assert(countCounts[oldCount]>=-1) : this;
				
				/* Increment the kmer count */
				final short newCount=counts[kmer]=(short)(oldCount+1);
				
				/* The count could at most be 1 more than the total window kmers here temporarily */
				assert(newCount<=windowKmers+1) : this;
				
				/* Increment the new countCount */
				countCounts[newCount]++;
				
				/* Update entropy */
				currentEsum=currentEsum+entropy[newCount]-entropy[oldCount];
//				currentEsum+=entropyDeltaPlus[oldCount];

				assert(!verify || unique==Tools.cardinality(counts)) : this;
				assert(!verify || (shared.Vector.sum(countCounts)>0 && (shared.Vector.sum(countCounts)<=windowKmers+1))) : this;
			}
		}
		if(verbose){System.err.println("B: counts="+Arrays.toString(counts)+"; countcounts="+Arrays.toString(countCounts));}
		
		//At this point the state is inconsistent as it may have one too many kmers.
		
		if(pos2>=0){//Remove the leftmost base
			final byte b2=(k>1 ? baseRingBuffer[pos2] : oldBase);//This is not the leftmost base, but the base to the right of the leftmost kmer
			final int n2=symbolToNumber0[b2];
			//Generate old leftmost kmer using an internal base near the left end
			kmer2=((kmer2<<bitsPerBase)|n2)&mask;

			if(verbose){System.err.println("B2: pos="+pos+", pos2="+pos2+"; b2="+Character.toString((char)b2)+"; kmer2="+kmer2);}

			if(len>windowBases){//Remove a kmer, only if a base is leaving the window
				baseCounts[n2]--;
				//System.err.println("Removing "+kmer2);
				
				//Update number of Ns in current window
				if(!isFullyDefined(oldBase)){
					ns--;
					assert(ns>=0) : "There are fewer than 0 Ns in the window:\n"+this;
				}
				
				assert(kmer2>=0) : "A negative kmer was observed: "+kmer2+"\n"+this;
				
				final short oldCount=counts[kmer2];
				assert(oldCount>0) : "Attempting to decrement a nonpositive count: \n"+oldCount+"\n"+this;
				
				//Decrement the old countCount
				countCounts[oldCount]--;
				
				assert(countCounts[oldCount]>=0) : "A countCount became negative: \n"+countCounts[oldCount]+"\n"+this;
				
				//Decrement the kmer count
				final short newCount=counts[kmer2]=(short)(oldCount-1);
				
				/* Increment the new countCount */
				countCounts[newCount]++;
				
				/* Update unique kmer count */
				if(newCount<1){
					assert(newCount==0): "An outgoing array has negative counts: \n"+this;
					unique--;
				}
				
				/* Update entropy */
				currentEsum=currentEsum+entropy[newCount]-entropy[oldCount];
//				currentEsum+=entropyDeltaMinus[oldCount];
				
				assert(!verify || unique==Tools.cardinality(counts)) : this;
				assert(!verify || (shared.Vector.sum(countCounts)>=0 && (shared.Vector.sum(countCounts)<=windowKmers))) : this;
			}
		}
		if(verbose){System.err.println("C: counts="+Arrays.toString(counts)+"; countcounts="+Arrays.toString(countCounts));}
		
		//Update position pointers
		//Can use modulo, but this is faster because the branch is normally skipped.
		//Ternary conditionals are also slower.
		pos++;
		pos2++;
		if(pos>=windowBases){pos=0;}
		if(pos2>=windowBases){pos2=0;}
		assert(k==1 || pos!=pos2) : "pos="+pos+", pos2="+pos2;
		
		//Test final state.
		assert(!verify || verify()) : this;
	}
	
	/**
	 * Reset fields to prepare for a new sequence.
	 */
	public void clear(){
		//Reset scalars
		kmer=kmer2=len=pos=unique=ns=0;
		pos=0;
		pos2=0-windowBases+k-1;
		currentEsum=0;
		
		//Clear mutable arrays.  Bases does not need to be cleared.
		Arrays.fill(baseCounts, (short)0); //4 operations
		Arrays.fill(counts, (short)0); //Time proportional to kmer space
		
		//Note - countCounts are only needed for medium speed mode.
		Arrays.fill(countCounts, (short)0); //Time proportional to window size
		
		//Sets the number of kmers with a count of zero to maximum.
		countCounts[0]=(short)windowKmers;
		
		//Verify the state was cleared
		assert(!verify || verifyClear()) : this;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Validation          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Verify that mutable fields were properly cleared.
	 * Throws an assertion error upon failure.
	 * @return True.
	 */
	public boolean verifyClear(){
		for(int c : baseCounts){assert(c==0) : this;}
		for(int c : counts){assert(c==0) : this;}
		for(int i=1; i<countCounts.length; i++){assert(countCounts[i]==0) : this;}
		assert(kmer==0 && kmer2==0) : this;
		assert(pos==0) : this;
		assert(len==pos) : this;
		assert(ns==0) : this;
		assert(unique==0) : this;
		assert(pos2<0) : this;
		assert(currentEsum==0) : this;
		return true;
	}
	
	/**
	 * Verify that internal state is consistent.
	 * Throws an assertion error upon failure.
	 * @return True.
	 */
	public boolean verify(){
		
		//Number of unique kmers in the window
		int existSum=0;
		
		//Total number of kmers in the window
		int countSum=0;
		
		//Number of undefined symbols in the window
		int nSum=0;
		
		//Check the kmer counts
		for(int c : counts){
			assert(c>=0 && c<=windowKmers) : "A kmer count exceeds the possible bounds.\n"+this;
			if(c>0){
				existSum++;
				countSum+=c;
			}
		}
		
		//Check the countCounts
		for(int cc : countCounts){
			assert(cc>=0 && cc<=windowKmers) : "A countCount exceeds the possible bounds.\n"+this;
		}
		
		//Count undefined symbols
		for(byte b : baseRingBuffer){
			if(!isFullyDefined(b)){nSum++;}
		}
		
		//Number of kmers with count 0
		final int cc0=countCounts[0];
		
		//Sum of countCounts
		final int ccSum=(int)shared.Vector.sum(countCounts);
		
		//Sum of nonzero countCounts; should equal the number of unique kmers
		final int ccSum1=ccSum-cc0;
		
		//Do assertions
		assert(len<windowBases || ns==nSum) : this;
		assert(existSum==unique) : this;
		assert(existSum==ccSum1) : this;
		assert(existSum>=0 && existSum<=windowKmers) : this;
		assert(len<windowBases || countSum==windowKmers) : this;
		assert(ccSum==windowKmers);
		assert(pos==len%windowBases) : this;
		
		//Ensure different entropy calculation methods are consistent
		if(len>=windowKmers){
			float a=calcEntropyFast();
			float b=calcEntropyMedium();
			float c=calcEntropySlow();
			float d=calcEntropySuperSlow(); //"d" should be identical to "c".
			assert(Tools.absdif(a, b)<0.0001) : "Fast and Medium entropy differ:\n"+a+"\n"+b+"\n"+c+"\n"+this;
			assert(Tools.absdif(b, c)<0.000001) : "Medium and Slow entropy differ:\n"+a+"\n"+b+"\n"+c+"\n"+this;
			assert(c==d) : "Slow and SuperSlow entropy differ:\n"+a+"\n"+b+"\n"+c+"\n"+d+"\n"+this;
		}
		return true;
	}
	
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append('\n');
		sb.append("kmer\t"+kmer).append('\n');
		sb.append("kmer2\t"+kmer2).append('\n');
		sb.append("pos\t"+pos).append('\n');
		sb.append("pos2\t"+pos2).append('\n');
		sb.append("len\t"+len).append('\n');
		sb.append("unique\t"+unique).append('\n');
		sb.append("ns\t"+ns).append('\n');
		sb.append('\n');
		sb.append("k\t"+k).append('\n');
		sb.append("windowBases\t"+windowBases).append('\n');
		sb.append("windowKmers\t"+windowKmers).append('\n');
		sb.append("mask\t"+mask).append('\n');
		sb.append('\n');
		sb.append("cardinality\t"+Tools.cardinality(counts)).append('\n');
		sb.append("counts\t"+Arrays.toString(counts)).append('\n');
		sb.append("ccounts\t"+Arrays.toString(countCounts)).append('\n');
		sb.append("bases\t"+basesToString()).append('\n');
		sb.append("entropy\t"+Arrays.toString(entropy)).append('\n');
		sb.append("entropyMult\t"+entropyMult).append('\n');
		sb.append("entropySum\t"+currentEsum).append('\n');
		return sb.toString();
	}
	
	/** Returns the ring buffer as a String in its correct order. */
	public String basesToString(){
		StringBuilder sb=new StringBuilder(baseRingBuffer.length);
		for(int i=0; i<baseRingBuffer.length; i++){
			byte b=(baseRingBuffer[(i+pos)%baseRingBuffer.length]);
			if(b==0){b='N';}
			sb.append((char)b);
		}
		return sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Mutable Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Current leading kmer; rightmost k bases of window */
	int kmer=0;
	/** Current trailing kmer; leftmost k-1 bases of window plus the removed base */
	int kmer2=0;
	/** Position in ring buffer to place incoming bases */
	int pos=0;
	/** Position in ring buffer to read next base for kmer2 */
	int pos2=0;
	/** Current number of processed bases.  Equal to pos without being reset at buffer wrap. */
	int len=0;
	/** Number of unique kmers in the current window */
	int unique=0;
	/** Number of undefined bases in the current window */
	int ns;

	/** Current sum of entropy from kmers in the current window */
	double currentEsum=0;
	
	/*--------------------------------------------------------------*/
	/*----------------        Mutable Arrays        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Number of times each base occurs in current window.
	 * Equivalent to counts if k=1. */
	private final short[] baseCounts;
	
	/** Number of times each kmer occurs in current window.
	 * Indexed by the kmer's numeric value.
	 * counts[0] stores the count of the kmer AAA, if k=3. */
	private final short[] counts;
	
	/** Number of instances of each number in counts.
	 * countCounts[0] stores the number of kmers with count 0.
	 * This is only needed in medium speed mode (or verify mode). */
	private final short[] countCounts;
	
	/** Ring buffer of bases in current window.
	 * Not strictly necessary, but convenient. */
	private final byte[] baseRingBuffer;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Kmer length for entropy calculation */
	private final int k;
	/** Window length for entropy calculation */
	private final int windowBases;
	/** Number of kmers in the window */
	private final int windowKmers;
	/** Amino acid mode */
	private final boolean amino;

	/** Bits per symbol */
	private final int bitsPerBase;
	/** Mask for sliding kmers */
	private final int mask;
	
	/** Minimum entropy to be considered "complex", on a scale of 0-1; optional */
	private final float entropyCutoff;
	/** Pass entropy values above the cutoff */
	private final boolean highPass;
	
	/** Number of possible unique kmers */
	private final int kmerSpace;
	/** A precalculated constant */
	private final double entropyMult;
	/** Array of precalculated constants */
	private final double[] entropy;
	
	/** Precalculated constant equal to 1f/windowBases */
	private final float baseCountMult;

	/** Translation table yielding 0 if undefined */
	private final byte[] symbolToNumber0;
	/** Translation table yielding -1 if undefined */
	private final byte[] symbolToNumber;
	
	final boolean isFullyDefined(byte symbol){
		return symbol>=0 && symbolToNumber[symbol]>=0;
	}
	
	//Note:  These incur fewer operations, but in testing, were not faster.
//	/** For calculating entropy running average quickly when adding a kmer.
//	 * entropyDeltaPlus[i] = (entropy[i+1]-entropy[i])*entropyMult */
//	private final double[] entropyDeltaPlus;
//	/** For calculating entropy running average quickly when removing a kmer.
//	 * entropyDeltaMinus[i] = (entropy[i-1]-entropy[i])*entropyMult */
//	private final double[] entropyDeltaMinus;
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Entropy calculation speed constants.
	 * FAST is less precise for long sequences.
	 * MEDIUM is probably most precise. */
	public static final int FAST=0, MEDIUM=1, SLOW=2, SUPERSLOW=3;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Kmer length for entropy calculation */
	public static int defaultK=5;
	public static boolean setDefaultK=false;
	/** Window length for entropy calculation */
	public static int defaultWindowBases=50;
	public static boolean setDefaultWindow=false;
//	/** Minimum entropy to be considered "complex", on a scale of 0-1 */
//	public static float defaultCutoff=-1;
	
	/** Entropy calculation mode */
	public static int speed=FAST;
	
	/** Verify consistency of related data structures (slow) */
	public static boolean verify=false;
	/** Verbose output */
	public static final boolean verbose=false;
	
}
