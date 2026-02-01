package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import shared.Random;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import kmer.AbstractKmerTable;
import kmer.HashArray1D;
import kmer.HashForest;
import kmer.KmerTable;
import kmer.ScheduleMaker;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * @author Brian Bushnell
 * @date Mar 24, 2014
 *
 */
public class CalcUniqueness {
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	
	/**
	 * Program entry point.
	 * Creates CalcUniqueness instance and processes input with timing.
	 * @param args Command-line arguments for input files and parameters
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		CalcUniqueness x=new CalcUniqueness(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor that parses arguments and initializes data structures.
	 * Sets up hash tables, counters, and file formats for processing.
	 * Validates parameters and configures k-mer tracking modes.
	 *
	 * @param args Command-line arguments containing options and file paths
	 * @throws RuntimeException if invalid parameters or missing required inputs
	 */
	public CalcUniqueness(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean setInterleaved=false; //Whether it was explicitly set.
		
		Shared.capBuffers(4);
		ReadWrite.USE_UNPIGZ=true;
		
		int k_=25;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("printlastbin") || a.equals("plb")){
				printLastBin=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("cumulative")){
				cumulative=Parse.parseBoolean(b);
			}else if(a.equals("offset")){
				singleOffset=Integer.parseInt(b);
			}else if(a.equals("percent") || a.equals("percents")){
				showPercents=Parse.parseBoolean(b);
			}else if(a.equals("count") || a.equals("counts")){
				showCounts=Parse.parseBoolean(b);
			}else if(a.equals("minprob") || a.equals("percents")){
				minprob=Float.parseFloat(b);
			}else if(a.equals("k")){
				k_=Integer.parseInt(b);
			}else if(a.equals("fixpeaks") || a.equals("fixspikes") || a.equals("fs")){
				fixSpikes=Parse.parseBoolean(b);
			}else if(a.equals("bin") || a.equals("interval")){
				interval=Parse.parseKMG(b);
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else{
//				System.err.println("Unknown parameter "+args[i]);
//				assert(false) : "Unknown parameter "+args[i];
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			minAverageQuality=parser.minAvgQuality;
			minAverageQualityBases=parser.minAvgQualityBases;
			
			maxReads=parser.maxReads;
			samplerate=parser.samplerate;
			sampleseed=parser.sampleseed;

			overwrite=parser.overwrite;
			append=parser.append;
			testsize=parser.testsize;
			
			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;

			out=parser.out1;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		setSampleSeed(-1);
		
		k=k_;
		k2=k-1;
		assert(k>0 && k<32) : "k="+k+"; valid range is 1-31";
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){System.err.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out==null){
			out="stdout.txt";
		}
		
		if(!setInterleaved){
			assert(in1!=null && out!=null) : "\nin1="+in1+"\nin2="+in2+"\nout="+out+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}
		}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
		}
		
		ffout=FileFormat.testOutput(out, FileFormat.TEXT, extout, false, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);

		keySets=new AbstractKmerTable[WAYS];

		//Initialize tables
		ScheduleMaker scheduleMaker=new ScheduleMaker(WAYS, 12, false, 1);
		int[] schedule=scheduleMaker.makeSchedule();
		for(int j=0; j<WAYS; j++){
			if(useForest){
				keySets[j]=new HashForest(initialSize, true, false);
			}else if(useTable){
				keySets[j]=new KmerTable(initialSize, true);
			}else if(useArray){
//				keySets[j]=new HashArray1D(initialSize, -1, -1L, true); //TODO: Set maxSize
				keySets[j]=new HashArray1D(schedule, -1L);
			}else{
				throw new RuntimeException("Must use forest, table, or array data structure.");
			}
		}
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Class          ----------------*/
	/*--------------------------------------------------------------*/
	
	private class Counter{
		
		/** Creates a counter with the specified bit mask for k-mer tracking.
		 * @param mask_ Bit mask used to identify this counter's k-mer category */
		Counter(int mask_){
			mask=mask_;
		}
		
		/**
		 * Updates quality statistics for a read.
		 * Calculates average quality and error-free probability.
		 * @param r The read to analyze for quality metrics
		 */
		void incrementQuality(Read r){
			qualCounts++;
			double q=r.avgQualityByProbabilityDouble(true, r.length());
			quality+=q;
			double p=r.probabilityErrorFree(true, r.length());
			perfectProb+=p;
		}
		
		/**
		 * Tracks k-mer occurrence in the appropriate hash table.
		 * Increments hit count if k-mer exists, miss count if new.
		 * Uses bit masking to identify which counter categories have seen this k-mer.
		 * @param kmer The k-mer value to track (negative values are ignored)
		 */
		void increment(final long kmer){
			if(kmer<0){return;}
			AbstractKmerTable table=keySets[(int)(kmer%WAYS)];
			int count=table.getValue(kmer);
			if(count<1){
				table.set(kmer, mask);
				misses++;
				cmisses++;
			}else if((count&mask)==0){
				table.set(kmer, count|mask);
				misses++;
				cmisses++;
			}else{
				hits++;
				chits++;
			}
		}
		
		/** Resets interval counters while preserving cumulative totals.
		 * Saves current values as previous for spike detection. */
		void reset(){
			prevPercent=percent();
			prevHits=hits;
			prevMisses=misses;
			
			hits=misses=0;
			quality=0;
			perfectProb=0;
			qualCounts=0;
		}
		
		/** Calculates average quality score across processed reads.
		 * @return Average quality score, or 0 if no reads processed */
		public double averageQuality() {
			return qualCounts<1 ? 0 : quality/qualCounts;
		}
		
		/** Calculates average probability of error-free reads.
		 * @return Average error-free probability as percentage, or 0 if no reads */
		public double averagePerfectProb() {
			return qualCounts<1 ? 0 : 100*perfectProb/qualCounts;
		}
		
		/**
		 * Calculates uniqueness percentage for this counter.
		 * Applies spike-fixing algorithm if enabled to smooth anomalous values.
		 * @return Percentage of unique k-mers (misses / total)
		 */
		double percent(){
			final long sum=hits()+misses(), prevSum=prevHits+prevMisses;
			if(sum==0){return 0;}
			double percent=misses()*100.0/sum;
			if(cumulative || !fixSpikes || prevSum==0){return percent;}
			assert(!cumulative && fixSpikes);
			return Tools.min(percent, prevPercent+0.2);
		}
		
		/** Formats uniqueness percentage as string with 3 decimal places.
		 * @return Formatted percentage string */
		String percentS(){
			return Tools.format("%.3f",percent());
		}

		/**
		 * Returns hit count (k-mers seen before).
		 * Uses cumulative or interval count based on mode.
		 * @return Number of k-mer hits
		 */
		long hits(){return cumulative ? chits : hits;}
		/**
		 * Returns miss count (unique k-mers).
		 * Uses cumulative or interval count based on mode.
		 * @return Number of k-mer misses
		 */
		long misses(){return cumulative ? cmisses : misses;}
		
		/** Bit mask used to identify which counter categories have seen a k-mer */
		final int mask;
		
		/** Cumulative probability of error-free reads */
		double perfectProb;
		/** Cumulative quality score across all processed reads */
		double quality;
		/** Number of reads used for quality calculations */
		long qualCounts;
		
		/** Per-interval hash hits */
		long hits=0;
		/** Per-interval hash misses */
		long misses=0;
		
		/** Cumulative hash hits */
		long chits=0;
		/** Cumulative hash misses */
		long cmisses=0;
		
		/** Previous interval hit count for spike detection */
		long prevHits=0;
		/** Previous interval miss count for spike detection */
		long prevMisses=0;
		/** Previous interval uniqueness percentage for spike detection */
		double prevPercent=0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Primary Method        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Main processing method that analyzes sequence uniqueness.
	 * Reads input files, extracts k-mers, and tracks occurrence in hash tables.
	 * Outputs statistics at regular intervals and handles paired-end data.
	 *
	 * @param t Timer for tracking execution performance
	 * @throws RuntimeException if processing encounters unrecoverable errors
	 */
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			cris.setSampleRate(samplerate, sampleseed);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		final boolean paired=cris.paired();
		if(verbose){System.err.println("Input is "+(paired ? "paired" : "unpaired"));}

		TextStreamWriter tsw=null;
		if(out!=null){
			tsw=new TextStreamWriter(ffout);
			tsw.start();
			tsw.print("#count");
			if(showPercents){
				tsw.print("\tfirst\trand");
				if(paired){tsw.print("\tr1_first\tr1_rand\tr2_first\tr2_rand\tpair");}
			}
			if(showCounts){
				tsw.print("\tfirst_cnt\trand_cnt");
				if(paired){tsw.print("\tr1_first_cnt\tr1_rand_cnt\tr2_first_cnt\tr2_rand_cnt\tpair_cnt");}
			}
			if(showQuality){
				tsw.print("\tavg_quality\tperfect_prob");
			}
			tsw.print("\n");
		}
		
		//Counters for overall data statistics
		long pairsProcessed=0;
		long readsProcessed=0;
		long basesProcessed=0;
		
		//Counter for display intervals
		long remaining=interval;
		
		final StringBuilder sb=new StringBuilder(1024);
		
		{
			//Fetch initial list
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			/* Process 1 list of reads per loop iteration */
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				/* Process 1 read per loop iteration */
				for(Read r1 : reads){
					if(minAverageQuality<1 || r1.avgQualityFirstNBases(minAverageQualityBases)>=minAverageQuality){
						final Read r2=r1.mate;
						final byte[] bases1=(r1==null ? null : r1.bases);
						final byte[] bases2=(r2==null ? null : r2.bases);
						final byte[] quals1=(r1==null ? null : r1.quality);
						final byte[] quals2=(r2==null ? null : r2.quality);
						final int length1=(bases1==null ? 0 : bases1.length);
						final int length2=(bases2==null ? 0 : bases2.length);

						pairsProcessed++;

						/* Process read 1 */
						if(r1!=null){
							
							bothCounterFirst.incrementQuality(r1);

							readsProcessed++;
							basesProcessed+=length1;

							if(length1>=k){
								if(length1>=k+singleOffset){//Fixed kmer
									final long kmer=toKmer(bases1, quals1, singleOffset, k);
									r1CounterFirst.increment(kmer);
									bothCounterFirst.increment(kmer);
								}
								{//Random kmer
									final long kmer=toKmer(bases1, quals1, randy.nextInt(length1-k2), k);
									r1CounterRand.increment(kmer);
									bothCounterRand.increment(kmer);
								}
							}
						}

						/* Process read 2 */
						if(r2!=null){
							
							bothCounterFirst.incrementQuality(r2);

							readsProcessed++;
							basesProcessed+=length2;

							if(length2>=k){
								if(length2>=k+singleOffset){//Fixed kmer
									final long kmer=toKmer(bases2, quals2, singleOffset, k);
									r2CounterFirst.increment(kmer);
									bothCounterFirst.increment(kmer);
								}
								{//Random kmer
									final long kmer=toKmer(bases2, quals2, randy.nextInt(length2-k2), k);
									r2CounterRand.increment(kmer);
									bothCounterRand.increment(kmer);
								}
							}
						}

						/* Process pair */
						if(r1!=null && r2!=null){

							if(length1>k+PAIR_OFFSET && length2>k+PAIR_OFFSET){
								final long kmer1=toKmer(bases1, quals1, PAIR_OFFSET, k);
								final long kmer2=toKmer(bases2, quals2, PAIR_OFFSET, k);
								if(kmer1!=-1 && kmer2!=-1){
									final long kmer=(~((-1L)>>2))|((kmer1<<(2*(31-k)))^(kmer2));
									assert(kmer>=0) : k+", "+kmer1+", "+kmer2+", "+kmer;
									{//Pair kmer
										pairCounter.increment(kmer);
									}
								}
							}
						}

						remaining--;
						if(remaining<=0){

							printCountsToBuffer(sb, pairsProcessed, paired);

							if(tsw!=null){tsw.print(sb.toString());}

							//Reset things
							sb.setLength(0);
							remaining=interval;
						}
					}
				}
				
				//Fetch a new list
				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){//Return final list
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		if(remaining<interval && printLastBin){
			
			printCountsToBuffer(sb, pairsProcessed, paired);
			
			if(tsw!=null){tsw.print(sb.toString());}
			
			//Reset things
			sb.setLength(0);
			remaining=interval;
		}
		
		
		//Close things
		errorState|=ReadWrite.closeStream(cris);
		if(tsw!=null){
			tsw.poisonAndWait();
			errorState|=tsw.errorState;
		}
		
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		if(testsize){
			long bytesProcessed=(new File(in1).length()+(in2==null ? 0 : new File(in2).length()));
			double xpnano=bytesProcessed/(double)(t.elapsed);
			String xpstring=(bytesProcessed<100000 ? ""+bytesProcessed : bytesProcessed<100000000 ? (bytesProcessed/1000)+"k" : (bytesProcessed/1000000)+"m");
			while(xpstring.length()<8){xpstring=" "+xpstring;}
			outstream.println("Bytes Processed:    "+xpstring+" \t"+Tools.format("%.2fm bytes/sec", xpnano*1000));
		}
		
		if(errorState){
			throw new RuntimeException("CalcUniqueness terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	/**
	 * Formats uniqueness statistics into output buffer.
	 * Includes percentages, counts, and quality metrics based on configuration.
	 * Resets interval counters after printing.
	 *
	 * @param sb StringBuilder to append formatted output
	 * @param pairsProcessed Number of read pairs processed so far
	 * @param paired Whether input data is paired-end
	 */
	private void printCountsToBuffer(StringBuilder sb, long pairsProcessed, boolean paired){
		
		//Display data for the last interval
		sb.append(pairsProcessed);
		
		if(showPercents){
			sb.append('\t');
			sb.append(bothCounterFirst.percentS());
			sb.append('\t');
			sb.append(bothCounterRand.percentS());
			if(paired){
				sb.append('\t');
				sb.append(r1CounterFirst.percentS());
				sb.append('\t');
				sb.append(r1CounterRand.percentS());
				sb.append('\t');
				sb.append(r2CounterFirst.percentS());
				sb.append('\t');
				sb.append(r2CounterRand.percentS());
				sb.append('\t');
				sb.append(pairCounter.percentS());
			}
		}
		
		if(showCounts){
			sb.append('\t');
			sb.append(bothCounterFirst.misses());
			sb.append('\t');
			sb.append(bothCounterRand.misses());
			if(paired){
				sb.append('\t');
				sb.append(r1CounterFirst.misses());
				sb.append('\t');
				sb.append(r1CounterRand.misses());
				sb.append('\t');
				sb.append(r2CounterFirst.misses());
				sb.append('\t');
				sb.append(r2CounterRand.misses());
				sb.append('\t');
				sb.append(pairCounter.misses());
			}
		}
		
		if(showQuality){
			sb.append('\t').append(Tools.format("%.2f", bothCounterFirst.averageQuality()));
			sb.append('\t').append(Tools.format("%.2f", bothCounterFirst.averagePerfectProb()));
		}
		
		sb.append('\n');

		bothCounterFirst.reset();
		bothCounterRand.reset();
		r1CounterFirst.reset();
		r1CounterRand.reset();
		r2CounterFirst.reset();
		r2CounterRand.reset();
		pairCounter.reset();
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Generate a kmer from specified start location
	 * @param bases
	 * @param start
	 * @param klen kmer length
	 * @return kmer
	 */
	private final long toKmer(final byte[] bases, byte[] quals, final int start, final int klen){
		if(minprob>0 && quals!=null){
			float prob=toProb(quals, start, klen);
			if(prob<minprob){return -1;}
		}
		final int stop=start+klen;
		assert(stop<=bases.length);
		long kmer=0;
		
		for(int i=start; i<stop; i++){
			final byte b=bases[i];
			final long x=Dedupe.baseToNumber[b];
			kmer=((kmer<<2)|x);
		}
		return kmer;
	}
	
	/**
	 * Generate the probability a kmer is error-free, from specified start location
	 * @param quals
	 * @param start
	 * @param klen kmer length
	 * @return kmer
	 */
	private final static float toProb(final byte[] quals, final int start, final int klen){
		final int stop=start+klen;
		assert(stop<=quals.length);
		float prob=1f;
		
		for(int i=start; i<stop; i++){
			final byte q=quals[i];
			float pq=probCorrect[q];
			prob*=pq;
		}
		return prob;
	}
	
	/*--------------------------------------------------------------*/
	
	/** Sets random seed for reproducible k-mer sampling.
	 * @param seed Random seed value (-1 for default thread-local random) */
	public void setSampleSeed(long seed){
		randy=Shared.threadLocalRandom(seed);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path for paired-end data */
	private String in2=null;
	
	/** Output file path for uniqueness statistics */
	private String out=null;
	
	/** Input file extension override */
	private String extin=null;
	/** Output file extension override */
	private String extout=null;
	
	/*--------------------------------------------------------------*/
	
	//Counters for hashtable hits and misses of different kmers
	/** Counter for fixed k-mers from read 1 */
	private final Counter r1CounterFirst=new Counter(1);
	/** Counter for random k-mers from read 1 */
	private final Counter r1CounterRand=new Counter(2);
	/** Counter for fixed k-mers from read 2 */
	private final Counter r2CounterFirst=new Counter(4);
	/** Counter for random k-mers from read 2 */
	private final Counter r2CounterRand=new Counter(8);
	/** Counter for k-mers derived from read pairs */
	private final Counter pairCounter=new Counter(16);

	/** Counter for fixed k-mers from all reads */
	private final Counter bothCounterFirst=new Counter(32);
	/** Counter for random k-mers from all reads */
	private final Counter bothCounterRand=new Counter(64);
	
	/*--------------------------------------------------------------*/
	
	
	/** Maximum number of reads to process (-1 for unlimited) */
	private long maxReads=-1;
	/** Fraction of reads to sample (1.0 for all reads) */
	private float samplerate=1f;
	/** Random seed for read sampling */
	private long sampleseed=-1;

	/** Number of reads between uniqueness calculations */
	private long interval=25000;
	/** Minimum probability threshold for k-mer acceptance */
	private float minprob=0;
	/** Minimum average quality score for read acceptance */
	private float minAverageQuality=0;
	/** Number of bases to consider for average quality calculation */
	private int minAverageQualityBases=20;
	/** Fixed position offset for k-mer extraction from reads */
	private int singleOffset=0;
	/** Whether to use cumulative counts instead of interval counts */
	private boolean cumulative=false;
	/** Whether to include percentage columns in output */
	private boolean showPercents=true;
	/** Whether to include count columns in output */
	private boolean showCounts=false;
	/** Whether to print statistics for final incomplete interval */
	private boolean printLastBin=false;
	/** Whether to include quality metrics in output */
	private boolean showQuality=true;
	/** Whether to apply spike-smoothing algorithm to percentages */
	private boolean fixSpikes=false;

	/** K-mer length minus 1 for array indexing */
	/** K-mer length for uniqueness analysis */
	private final int k, k2;
	/** Number of hash tables for k-mer distribution */
	private static final int WAYS=31;
	/** Fixed offset for pair k-mer extraction */
	private static final int PAIR_OFFSET=10;
	
	/** Initial size of data structures */
	private int initialSize=512000;
	
	/** Hold kmers.  A kmer X such that X%WAYS=Y will be stored in keySets[Y] */
	private final AbstractKmerTable[] keySets;
	
	/*--------------------------------------------------------------*/
	
	/** File format for primary input */
	private final FileFormat ffin1;
	/** File format for secondary input */
	private final FileFormat ffin2;

	/** File format for output */
	private final FileFormat ffout;
	
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for messages and statistics */
	private PrintStream outstream=System.err;
	/** Whether to print verbose debugging information */
	public static boolean verbose=false;
	/** Whether processing encountered unrecoverable errors */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files */
	private boolean append=false;
	/** Whether to report file size processing statistics */
	private boolean testsize=false;
	
	/** Whether to use HashArray1D data structure */
	/** Whether to use KmerTable data structure */
	/** Whether to use HashForest data structure */
	private static final boolean useForest=false, useTable=false, useArray=true;
	
	/** Random number generator for k-mer position sampling */
	private Random randy;
	
	/** Lookup table converting quality scores to error probabilities */
	private static final float[] probCorrect=
		{0.0000f, 0.2501f, 0.3690f, 0.4988f, 0.6019f, 0.6838f, 0.7488f, 0.8005f, 0.8415f, 0.8741f, 0.9000f, 0.9206f, 0.9369f, 0.9499f,
		 0.9602f, 0.9684f, 0.9749f, 0.9800f, 0.9842f, 0.9874f, 0.9900f, 0.9921f, 0.9937f, 0.9950f, 0.9960f, 0.9968f, 0.9975f, 0.9980f,
		 0.9984f, 0.9987f, 0.9990f, 0.9992f, 0.9994f, 0.9995f, 0.9996f, 0.9997f, 0.9997f, 0.9998f, 0.9998f, 0.9999f, 0.9999f, 0.9999f,
		 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f,
		 0.9999f, 0.9999f, 0.9999f, 0.9999f, 1f, 1f, 1f, 1f, 1f, 1f, 1f, 1f, 1f, 1f, 1f, 1f, 1f, 1f, 1f, 1f, 1f, 1f, 1f, 1f, 1f, 1f};
	
	
}
