package hiseq;

import java.util.Arrays;

import align2.QualityTools;
import dna.AminoAcid;
import shared.Tools;
import stream.Read;
import stream.SamLine;
import structures.ByteBuilder;

/**
 * Represents a spatial tile region on HiSeq flow cells for quality assessment.
 * Accumulates quality metrics, error rates, and base composition statistics
 * for reads originating from a defined rectangular coordinate region.
 * Used for spatial quality analysis and tile-based filtering.
 *
 * @author Brian Bushnell
 */
public class MicroTile implements Comparable<MicroTile>{

	/** Creates a MicroTile with default zero coordinates */
	public MicroTile(){this(0,0,0,0,0,0);}

	/**
	 * Creates a MicroTile with specified spatial coordinates.
	 *
	 * @param lane_ HiSeq lane number
	 * @param tile_ Tile identifier within the lane
	 * @param x1_ Left boundary x-coordinate
	 * @param x2_ Right boundary x-coordinate
	 * @param y1_ Bottom boundary y-coordinate
	 * @param y2_ Top boundary y-coordinate
	 */
	public MicroTile(int lane_, int tile_, int x1_, int x2_, int y1_, int y2_){
		lane=lane_;
		tile=tile_;
		x1=x1_;
		x2=x2_;
		y1=y1_;
		y2=y2_;
	}
	
	/** Processes the cycle tracker if cycle tracking is enabled */
	void process(){
		if(tracker!=null){tracker.process();}
	}
	
	/**
	 * Tests if coordinates fall within this tile's boundaries.
	 * @param x X-coordinate to test
	 * @param y Y-coordinate to test
	 * @return true if coordinates are within tile boundaries
	 */
	public boolean contains(int x, int y){
		return x>=x1 && x<=x2 && y>=y1 && y<=y2;
	}
	
	@Override
	public String toString(){
		return lane+", "+tile+", "+x1+", "+x2+", "+y1+", "+y2;
	}
	
	/** Calculates average read quality using probability-weighted scoring.
	 * @return Average quality score by probability, or 0 if no reads processed */
	public double averageReadQualityByProb(){
		return readCount==0 ? 0 : readQualityByProbSum/readCount;
	}
	
	/** Calculates expected base error rate from quality scores.
	 * @return Average expected error probability per base */
	public double averageExpectedBaseErrorRate(){
		return baseCount==0 ? 0 : baseErrorProbSum/baseCount;
	}
	
	/** Converts expected base error rate to Phred scale.
	 * @return Expected base error rate as Phred quality score */
	public double averageExpectedBaseErrorRatePhred(){
		return QualityTools.probErrorToPhredDouble(averageExpectedBaseErrorRate());
	}
	
	/** Calculates percentage of reads expected to be error-free.
	 * @return Percentage of error-free reads based on quality scores */
	public double percentErrorFree(){
		return readCount==0 ? 0 : probErrorFreeSum/readCount;
	}
	
	/** Calculates fraction of k-mers with high quality scores.
	 * @return Ratio of good k-mers to total k-mers */
	public double goodKmerFraction(){
		double kmers=goodKmerSum+validKmerSum;
		return goodKmerSum/Tools.max(kmers, 1);
	}
	
	/** Calculates fraction of reads that aligned successfully.
	 * @return Alignment rate as ratio of aligned to total reads */
	public double alignmentRate(){
		return readCount==0 ? 0 : alignedReadCount/(double)readCount;
	}
	
	//Small sample sizes will drift toward 23. 
	/**
	 * Calculates empirical quality score from alignment error counts.
	 * Applies Bayesian smoothing with pseudocounts for small samples.
	 * Small samples drift toward Phred 23 due to prior assumptions.
	 * @return True quality score in Phred scale based on alignment errors
	 */
	public double trueQuality(){
		double e=baseErrorCount+1;
		double b=alignedBaseCount+200;
		double prob=e/b;
		double phred=QualityTools.probErrorToPhredDouble(prob);
//		System.err.println(baseErrorCount+", "+alignedBaseCount+", "+prob+", "+phred);
		return phred;
	}
	
	//Small sample sizes will drift toward 0.2. 
	/**
	 * Calculates per-read error rate from alignment data.
	 * Applies Bayesian smoothing with 0.2 pseudocount for stability.
	 * Small samples drift toward 0.2 error rate.
	 * @return Error rate per read
	 */
	public double readErrorRate(){
		double e=readErrorCount+0.2;
		double b=alignedReadCount+1;
		double rate=e/b;
		return rate;
	}
	
	//Small sample sizes will drift toward 0.002. 
	/**
	 * Calculates per-base error rate from alignment data.
	 * Applies Bayesian smoothing with 0.002 pseudocount for stability.
	 * Small samples drift toward 0.002 error rate.
	 * @return Error rate per base
	 */
	public double baseErrorRate(){
		double e=baseErrorCount+0.002;
		double b=alignedBaseCount+1;
		double rate=e/b;
		return rate;
	}
	
	/** Calculates insertion rate per read from alignment data.
	 * @return Fraction of reads containing insertions */
	public double readInsRate(){
//		System.err.println(alignedReadCount+", "+readInsCount+", "+
//		(readInsCount/(double)alignedReadCount));
		if(alignedReadCount==0) {return 0;}
//		double e=readInsCount+0.01;
//		double b=alignedReadCount+10;
		double e=readInsCount;
		double b=alignedReadCount;
		double rate=e/b;
		return rate;
	}
	 
	/** Calculates deletion rate per read from alignment data.
	 * @return Fraction of reads containing deletions */
	public double readDelRate(){
		if(alignedReadCount==0) {return 0;}
//		double e=readDelCount+0.01;
//		double b=alignedReadCount+10;
		double e=readDelCount;
		double b=alignedReadCount;
		double rate=e/b;
		return rate;
	}
	 
	/** Calculates k-mer-based error rate per read.
	 * @return Rate of reads with k-mer-detected errors */
	public double kmerErrorRateR(){
		if(readCount==0) {return 0;}
		double e=kmerReadErrorCount;
		double b=readCount;
		double rate=e/b;
		return rate;
	}
	 
	/** Calculates k-mer-based error rate per base.
	 * @return Rate of bases with k-mer-detected errors per read */
	public double kmerErrorRateB(){
		if(readCount==0) {return 0;}
		double e=kmerBaseErrorCount;
		double b=readCount;
		double rate=e/b;
		return rate;
	}
	
	/** Calculates percentage of k-mers found in reference database.
	 * @return Percentage of k-mer hits */
	public double hitPercent(){
		long count=kmerCount();
		return count==0 ? 0 : hits*100.0/count;
	}
	
	/** Calculates percentage of unique k-mers not found in database.
	 * @return Percentage of unique k-mers */
	public double uniquePercent(){
		long count=kmerCount();
		return count==0 ? 0 : misses*100.0/count;
	}
	
	/** Calculates percentage of reads with poly-G contamination.
	 * @return Percentage of reads exceeding minimum poly-G threshold */
	public double polyGPercent(){
		long count=readCount;
		return count==0 ? 0 : homoPolyGCount*100.0/count;
	}
	
	/** Calculates average k-mer depth from database lookups.
	 * @return Average depth per k-mer */
	public double depth(){
		long count=kmerCount();
		return depthSum*1.0/count;
	}
	
	/** Gets average G content from cycle tracker.
	 * @return Average G percentage, or 0 if tracking disabled */
	public double avgG(){
		return tracker==null ? 0 : tracker.avg('G');
	}
	
	/** Gets maximum G content from cycle tracker.
	 * @return Maximum G percentage, or 0 if tracking disabled */
	public double maxG(){
		return tracker==null ? 0 : tracker.max('G');
	}
	
	/**
	 * Calculates implied error rate from base error rate function.
	 * Uses unique k-mer percentage to estimate per-base error rates.
	 * @param berf Base error rate function coefficients [intercept, slope]
	 * @return Estimated per-base error rate
	 */
	public double impliedErrorRate(double[] berf) {
		if(berf==null) {return 0;}
		double rootBer=Tools.mid(0.000001, 0.75, berf[0]+berf[1]*uniquePercent());
		return rootBer*rootBer;
	}

	/** Gets total k-mer count (hits plus misses) */
	public long kmerCount(){return hits+misses;}
	
	/** Resets all statistics counters to zero.
	 * Clears k-mer counts, error counts, base composition, and quality metrics. */
	public void clear() {
		hits=0;
		misses=0;
		depthSum=0;
		readCount=0;
		baseCount=0;
		alignedReadCount=0;
		alignedBaseCount=0;
		readErrorCount=0;
		baseErrorCount=0;
		kmerReadErrorCount=0;
		kmerBaseErrorCount=0;
		readInsCount=0;
		readDelCount=0;
		readQualityByProbSum=0;
		probErrorFreeSum=0;
		baseErrorProbSum=0;
		goodKmerSum=0;
		validKmerSum=0;
		Arrays.fill(acgtn, 0);
		
		homoPolyGCount=0;
		homoPolyGSum=0;
		if(TRACK_CYCLES){
			tracker.clear();
		}
		barcodes=0;
		barcodeHDistSum=0;
		barcodePolymers=0;
		
		mergedReads=0;
		insertSum=0;
		overlapSum=0;
		mergeErrorSum=0;
	}
	
	/**
	 * Adds statistics from another MicroTile to this one.
	 * Merges all counters including reads, errors, k-mers, and composition.
	 * @param mt MicroTile to merge statistics from
	 */
	public void add(MicroTile mt) {
		hits+=mt.hits;
		misses+=mt.misses;
		depthSum+=mt.depthSum;
		readCount+=mt.readCount;
		baseCount+=mt.baseCount;
		alignedReadCount+=mt.alignedReadCount;
		alignedBaseCount+=mt.alignedBaseCount;
		readErrorCount+=mt.readErrorCount;
		baseErrorCount+=mt.baseErrorCount;
		kmerReadErrorCount+=mt.kmerReadErrorCount;
		kmerBaseErrorCount+=mt.kmerBaseErrorCount;
		readInsCount+=mt.readInsCount;
		readDelCount+=mt.readDelCount;
		readQualityByProbSum+=mt.readQualityByProbSum;
		probErrorFreeSum+=mt.probErrorFreeSum;
		baseErrorProbSum+=mt.baseErrorProbSum;
		goodKmerSum+=mt.goodKmerSum;
		validKmerSum+=mt.validKmerSum;

		for(int i=0; i<acgtn.length; i++){
			acgtn[i]+=mt.acgtn[i];
		}
		homoPolyGCount+=mt.homoPolyGCount;
		homoPolyGSum+=mt.homoPolyGSum;
		if(TRACK_CYCLES){
			tracker.add(mt.tracker);
		}
		barcodes+=mt.barcodes;
		barcodeHDistSum+=mt.barcodeHDistSum;
		barcodePolymers+=mt.barcodePolymers;
		
		mergedReads+=mt.mergedReads;
		insertSum+=mt.insertSum;
		overlapSum+=mt.overlapSum;
		mergeErrorSum+=mt.mergeErrorSum;
	}
	
	/**
	 * Multiplies all statistics by a scaling factor.
	 * Used for normalizing tile statistics or creating projections.
	 * Does not handle cycle tracker multiplication.
	 * @param mult Scaling factor to apply to all counters
	 */
	public void multiplyBy(double mult) {
		hits=(long)(mult*hits);
		misses=(long)(mult*misses);
		depthSum=(long)(mult*depthSum);
		readCount=(long)(mult*readCount);
		baseCount=(long)(mult*baseCount);
		alignedReadCount=(long)(mult*alignedReadCount);
		alignedBaseCount=(long)(mult*alignedBaseCount);
		readErrorCount=(long)(mult*readErrorCount);
		baseErrorCount=(long)(mult*baseErrorCount);
		kmerReadErrorCount=(long)(mult*kmerReadErrorCount);
		kmerBaseErrorCount=(long)(mult*kmerBaseErrorCount);
		readInsCount=(long)(mult*readInsCount);
		readDelCount=(long)(mult*readDelCount);
		readQualityByProbSum=(double)(mult*readQualityByProbSum);
		probErrorFreeSum=(double)(mult*probErrorFreeSum);
		baseErrorProbSum=(double)(mult*baseErrorProbSum);
		goodKmerSum=(double)(mult*goodKmerSum);
		validKmerSum=(long)(mult*validKmerSum);

		for(int i=0; i<acgtn.length; i++){
			acgtn[i]=(long)(mult*acgtn[i]);
		}
		homoPolyGCount=(long)(mult*homoPolyGCount);
		homoPolyGSum=(long)(mult*homoPolyGSum);
		if(TRACK_CYCLES){
//			tracker.add(mt.tracker));
			assert(false) : "TRACK_CYCLES: Not handled with multiplyBy.";
		}
		barcodes=(long)(mult*barcodes);
		barcodeHDistSum=(long)(mult*barcodeHDistSum);
		barcodePolymers=(long)(mult*barcodePolymers);
		
		mergedReads=(long)(mult*mergedReads);
		insertSum=(long)(mult*insertSum);
		overlapSum=(long)(mult*overlapSum);
		mergeErrorSum=(long)(mult*mergeErrorSum);
	}
	
	/**
	 * Counts and scores k-mers based on quality values.
	 * Uses sliding window to calculate probability of k-mer correctness.
	 * Updates goodKmerSum and validKmerSum statistics.
	 *
	 * @param quals Quality score array
	 * @param k K-mer length for scoring window
	 * @return Total quality-weighted k-mer score
	 */
	private float countGoodKmers(byte[] quals, int k) {
		int valid=0;
		float good=0;
		float product=1;
		for(int i=0, len=0; i<quals.length; i++) {
			byte q=quals[i];
			if(q>0) {
				len++;
				float pc=QualityTools.PROB_CORRECT[q];
				product=product*pc;
				if(len>=k) {
					if(len>=k) {
						product=product*QualityTools.PROB_CORRECT_INVERSE[quals[i-k]];
					}
					valid++;
					good+=product;
				}
			}else {
				len=0;
				product=1;
			}
		}
		goodKmerSum+=good;
		validKmerSum+=valid;
		return good;
	}
	
	/**
	 * Processes a read and updates all relevant statistics.
	 * Analyzes quality scores, alignment information, base composition,
	 * poly-G content, and k-mer quality metrics.
	 * @param r Read to process and extract statistics from
	 */
	public void add(Read r){
		if(r==null){return;}
		final int len=r.length();
		if(len<1){return;}
		final SamLine sl=r.samline;
		final byte[] match=r.match;
		
		readCount++;
		baseCount+=r.length();
		readQualityByProbSum+=r.avgQualityByProbabilityDouble(true, len);
		probErrorFreeSum+=100*r.probabilityErrorFree(true, len);
		baseErrorProbSum+=r.expectedErrors(true, len);
		
		countGoodKmers(r.quality, 62);
		
//		if(r.mapped() || (sl!=null && sl.mapped()){alignedRead
		
		if(match!=null) {
			int bc=r.countAlignedBases();
			if(bc>0) {
				alignedReadCount++;
				alignedBaseCount+=bc;
				int errors=r.countErrors();
				readErrorCount+=(errors>0 ? 1 : 0);
				baseErrorCount+=errors;
				int[] mSCNID=Read.countMatchEvents(match);
				readInsCount+=(mSCNID[4]>0 ? 1 : 0);
				readDelCount+=(mSCNID[5]>0 ? 1 : 0);
			}
		}else if(sl!=null && sl.mapped()) {
			alignedReadCount++;
		}
		
		final byte[] bases=r.bases;
		int maxPolyG=0, currentPolyG=0;
		for(int i=0; i<len; i++){
			byte b=bases[i];
			byte x=AminoAcid.baseToNumberACGTN[b];
			acgtn[x]++;
			if(b=='G'){
				currentPolyG++;
				maxPolyG=Tools.max(currentPolyG, maxPolyG);
			}else{
				currentPolyG=0;
			}
		}
		final boolean polyg=(maxPolyG>=MIN_POLY_G);
		homoPolyGCount+=(polyg ? 1 : 0);
		homoPolyGSum+=(polyg ? maxPolyG : 0);
		r.setDiscarded(polyg);
		if(TRACK_CYCLES){
			tracker.add(r);
		}
	}
	
	/** 
	 * Has some slow functions offloaded to increase
	 * concurrency, but did not result in a speedup.
	 */
	public void addQuick(Read r){
		if(r==null){return;}
		final int len=r.length();
		if(len<1){return;}
		
		readCount++;
		baseCount+=r.length();
		
		countGoodKmers(r.quality, 62);
		
		final byte[] bases=r.bases;
		int maxPolyG=0, currentPolyG=0;
		for(int i=0; i<len; i++){
			byte b=bases[i];
			byte x=AminoAcid.baseToNumberACGTN[b];
			acgtn[x]++;
			if(b=='G'){
				currentPolyG++;
				maxPolyG=Tools.max(currentPolyG, maxPolyG);
			}else{
				currentPolyG=0;
			}
		}
		final boolean polyg=(maxPolyG>=MIN_POLY_G);
		homoPolyGCount+=(polyg ? 1 : 0);
		homoPolyGSum+=(polyg ? maxPolyG : 0);
		r.setDiscarded(polyg);
		if(TRACK_CYCLES){
			tracker.add(r);
		}
	}
	
	/**
	 * Generates column headers for tabular output.
	 * Returns either short or long format headers based on shortHeader flag.
	 * @return Tab-delimited header string for statistics output
	 */
	public static String header() {
		if(shortHeader) {
			return "lane\ttile\tx1\tx2\ty1\ty2"
					+ "\treads\tbases\talnRead\talnBase\terrAlnR\terrAlnB"
					+ "\terrKR\terrKB"
					+ "\tinsCnt\tdelCnt\tARERate\tABERate"
					+ "\tunique\tavQScor\tprobEF\tavBEPrb\tdepth"
					+ "\tIERate1\tIERate2\tIERate3\tIQScore"
					+ "\talnRate\ttruQual"
					+ "\teKRRate\teKBRate\tinsRate\tdelRate"
					+ "\tdiscard"
					+ "\tA\tC\tG\tT\tN\tpolyG\tplGLen\tplGRate"
					+ "\tBCCount\tBCHDist\tBCPoly\tBCHDAv\tBCPlyAv"
					+ "\tValKmrs\tGdKmrs\tGdKRate"
					+ "\tMerged\tInsert\tOverlap\tMrgErr"
					+ "\tAvInsrt\tMrgRate\tMrgBER";
		}
		
		return "lane\ttile\tx1\tx2\ty1\ty2"
		+ "\treads\tbases\talignedRead\talignedBase\treadsAlignedWithErrors\talignedErrorCount"
		+ "\treadsWithKmerErrors\tkmerErrorCount"
		+ "\tinsertionCount\tdeletionCount\tAligedReadErrorRate\tAlignedBaseErrorRate"
		+ "\tuniqueKmerRate\tavgQScoreByProb\tprobErrorFree\taverageBaseErrorProb\tdepth"
		+ "\tinferredErrorRate1\tinferredErrorRate2\tinferredErrorRate3\timpliedQualityScore"
		+ "\talignmentRate\ttrueQuality"
		+ "\treadsWithKmerErrorsRate\tkmerErrorsPerRead\tinsertionRate\tdeletionRate"
		+ "\tdiscard"
		+ "\tA\tC\tG\tT\tN\tpolyG_Count\tpolyG_Length\tpolyG_Rate"
		+ "\tBarcodeCount\tBarcodeHDistSum\tBarcodePolymers\tBarcodeHDistAvg\tBarcodePolymerRate"
		+ "\tValidKmers\tGoodKmers\tGoodKmerRate"
		+ "\tMergedReads\tInsertSum\tOverlapSum\tMergeErrors"
		+ "\tAvgInsertSize\tMergeRate\tMergeBaseErrorRate";
	}
	
	/**
	 * Formats all tile statistics as tab-delimited text output.
	 * Includes coordinate information, quality metrics, error rates,
	 * base composition, and derived statistics.
	 *
	 * @param bb ByteBuilder to append formatted output to
	 * @param k K-mer length for calculations
	 * @param HG High-depth genomic k-mer fraction
	 * @param rerf Read error rate function coefficients
	 * @param berf Base error rate function coefficients
	 * @return ByteBuilder with appended statistics line
	 */
	public ByteBuilder toText(ByteBuilder bb, int k, double HG, double[] rerf, double[] berf){
		bb.append(lane).tab();
		bb.append(tile).tab();
		bb.append(x1).tab();
		bb.append(x2).tab();
		bb.append(y1).tab();
		bb.append(y2).tab();

		bb.append(readCount).tab();
		bb.append(baseCount).tab();
		bb.append(alignedReadCount).tab();
		bb.append(alignedBaseCount).tab();
		bb.append(readErrorCount).tab();
		bb.append(baseErrorCount).tab();
		bb.append(kmerReadErrorCount).tab();
		bb.append(kmerBaseErrorCount).tab();
		bb.append(readInsCount).tab();
		bb.append(readDelCount).tab();
		bb.append(readErrorRate(), 5).tab();
		bb.append(baseErrorRate(), 5).tab();
		
		final double uniquePercent=uniquePercent();
		bb.append(uniquePercent, 4).tab();
		bb.append(averageReadQualityByProb(), 4).tab();
		bb.append(percentErrorFree(), 4).tab();
		bb.append(averageExpectedBaseErrorRate(), 5).tab();
		final double depth=depth();
		if(depth>10000) {
			bb.append((long)Math.round(depth)).tab();
		}else {
			bb.append(depth, depth>=100 ? 2 : 4).tab();
		}
		final double E1=(HG>0 && hits+misses>0 ? calcErrorRateFromUniqueness(HG, k) : 0);
		bb.append(E1, 5).tab();
		final double avgBases=baseCount/(Tools.max(1.0, readCount));
		assert(rerf!=null || alignedReadCount<1 || depthSum<1) : readCount+", "+alignedReadCount;
		final double rer=(rerf==null ? 0 : Tools.mid(0.000001, 0.999999, rerf[0]+rerf[1]*uniquePercent));
//		assert(rer>0.001) : rerf[0];
		final double E2=1-Math.pow(1-rer, 1/avgBases);
		bb.append(E2, 5).tab();
		final double rootber=(berf==null ? 0 : Tools.mid(0.000001, 0.75, berf[0]+berf[1]*uniquePercent));
		final double ber=rootber*rootber;
		bb.append(ber, 5).tab();
		bb.append(QualityTools.probErrorToPhredDouble(ber), 4).tab();
		
		bb.append(alignmentRate(), 5).tab();
		bb.append(trueQuality(), 4).tab();
		bb.append(kmerErrorRateR(), 5).tab();
		bb.append(kmerErrorRateB(), 5).tab();
		bb.append(readInsRate(), 5).tab();
		bb.append(readDelRate(), 5).tab();
		bb.append(discard);
		
		for(int i=0; i<5; i++){
			bb.tab().append(acgtn[i]);
		}
		bb.tab().append(homoPolyGCount);
		bb.tab().append(homoPolyGSum);
		bb.tab().append(homoPolyGCount/(double)readCount, 5);
//		assert(false) : homoPolyGSum+"\n"+bb;
		
		bb.tab().append(barcodes);
		bb.tab().append(barcodeHDistSum);
		bb.tab().append(barcodePolymers);
		double invBarcodes=1.0/(Tools.max(1.0, barcodes));
		bb.tab().append(barcodeHDistSum*invBarcodes, 5);//avg hdist
		bb.tab().append(barcodePolymers*invBarcodes, 5);//polymer rate

		bb.tab().append(validKmerSum);
		bb.tab().append(goodKmerSum, (goodKmerSum>=999 ? 0 : 2), true);
		bb.tab().append(goodKmerSum/(Tools.max(1.0, validKmerSum)), 5);
		
		if(mergedReads>0) {
			bb.tab().append(mergedReads);
			bb.tab().append(insertSum);
			bb.tab().append(overlapSum);
			bb.tab().append(mergeErrorSum);
			bb.tab().append(insertSum*2/Tools.max(1.0, mergedReads), 3);//Avg insert
			bb.tab().append(mergedReads/Tools.max(1.0, readCount), 5);//Merge rate
			bb.tab().append(mergeErrorSum/Tools.max(1.0, 2*overlapSum), 5);//Merge error rate
		}
		
		return bb.nl();
	}
	
	/**
	 * Although this seems like it should work, it did not work
	 * very well at all in a test of real data.
	 * @param HG_lane High depth genomic kmer fraction of lane.
	 * @return E_tile Per-base error rate of this tile.
	 */
	public float calcErrorRateFromUniqueness(double HG_lane, int k) {
		//Unique kmer rate
		double U_tile=uniquePercent()*0.01;
		
		//Non-unique kmer rate
		double NU_tile=1-U_tile;
		
		//Kmer correctness probablity
		double P_tile=NU_tile/HG_lane;
//		assert(P_tile<=1) : HG_lane+", "+U_tile+", "+NU_tile+", "+P_tile;
		P_tile=Tools.mid(0.000001, P_tile, 0.99999);

		//Per-base correctness probability
		double NE_tile=Math.pow(P_tile, 1.0/k);
		//Per-base error rate
		double E_tile=1-NE_tile;
		
		return (float)E_tile;
	}
	
	/** Sorts better tiles first */
	@Override
	public int compareTo(MicroTile mt) {
		double ua=uniquePercent();
		double ub=mt.uniquePercent();
		if(ua!=ub) {return ua>ub ? 1 : -1;}
		double qa=averageReadQualityByProb();
		double qb=mt.averageReadQualityByProb();
		if(ua!=ub) {return qa<qb ? 1 : -1;}
		if(readCount!=mt.readCount) {return readCount>mt.readCount ? -1 : 1;}
		return 0;
	}
	
	/** Number of k-mers found in reference database */
	public long hits;
	/** Number of k-mers not found in reference database */
	public long misses;
	/** Sum of k-mer depths from database lookups */
	public long depthSum;
	/** Total number of reads processed */
	public long readCount;
	/** Total number of bases in all reads */
	public long baseCount;
	/** Number of reads that aligned successfully */
	public long alignedReadCount;
	/** Number of bases in aligned reads */
	public long alignedBaseCount;
	/** Number of aligned reads containing errors */
	public long readErrorCount;//Reads aligned with errors
	/** Number of bases aligned with errors */
	public long baseErrorCount;//Bases aligned with errors
	/** Number of reads with k-mer-detected errors */
	public long kmerReadErrorCount;//Reads with errors detected
	/** Number of bases detected as errors by k-mer analysis */
	public long kmerBaseErrorCount;//Bases detected as errors
	/** Number of reads containing insertions */
	public long readInsCount;//Number of reads containing insertions
	/** Number of reads containing deletions */
	public long readDelCount;//Number of reads containing deletions
	/** Sum of read qualities weighted by error probability */
	public double readQualityByProbSum;
	/** Sum of error-free probabilities across all reads */
	public double probErrorFreeSum;
	/** Sum of expected error probabilities for all bases */
	public double baseErrorProbSum;
	
	/** Sum of quality-weighted k-mer scores */
	public double goodKmerSum;
	/** Total number of valid k-mers processed */
	public long validKmerSum;
	
	/** Base composition counts for A, C, G, T, N */
	public long[] acgtn=new long[5];
	/** Number of reads with poly-G runs exceeding threshold */
	public long homoPolyGCount;
	/** Total length of poly-G runs across all reads */
	public long homoPolyGSum;
	
	/** Counter for discarded reads */
	public int discard=0;
	
	/** HiSeq lane number for this tile */
	public final int lane;
	/** Tile identifier within the lane */
	public final int tile;
	public final int x1, x2;
	public final int y1, y2;
	
//	long[] barcodeHDist=new long[4];
	/** Number of barcodes processed */
	long barcodes;
	/** Sum of barcode Hamming distances */
	long barcodeHDistSum=0;
	/** Number of barcodes containing homopolymer runs */
	long barcodePolymers=0;
	
	/** Number of successfully merged paired reads */
	long mergedReads=0;
	/** Sum of insert sizes for merged reads */
	long insertSum=0;
	/** Sum of overlap lengths for merged reads */
	long overlapSum=0;
	/** Sum of errors detected during read merging */
	long mergeErrorSum=0;
	
	/** Cycle-by-cycle quality tracking, null if tracking disabled */
	public final CycleTracker tracker=TRACK_CYCLES ? new CycleTracker() : null;

	/** Minimum poly-G length threshold for contamination detection */
	public static int MIN_POLY_G=15;
	/** Whether to enable cycle-by-cycle quality tracking */
	public static boolean TRACK_CYCLES=false;
	/** Whether to use abbreviated column headers in output */
	public static boolean shortHeader=true;
	
}
