package cardinality;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import shared.Random;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.CallPeaks;
import parse.Parse;
import parse.Parser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Streamer;
import stream.StreamerFactory;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;
import structures.LongList;
import tracker.ReadStats;

/**
 * Wrapper class for LogLog cardinality estimation testing and benchmarking.
 * Provides functionality to test LogLog cardinality tracking algorithms with both
 * synthetic and real sequence data, supporting multiple trials and statistical analysis.
 *
 * @author Brian Bushnell
 * @date 2014
 */
class LogLogWrapper {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Runs multiple cardinality estimation trials and reports statistical summary
	 * including harmonic mean, median, standard deviation, and percentiles.
	 * @param args Command-line arguments for configuration
	 */
	public static final void main(String[] args){
		final LogLogWrapper llw=new LogLogWrapper(args);

		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;

		final Timer t=new Timer();
		final LongList results=new LongList();
		for(int i=0; i<llw.trials; i++){
			final long card=llw.process();
			results.add(card);
		}

		if(llw.trials>1){
			results.sort();

			final long kmers=llw.kmersProcessed;
			final double mean=results.mean();
			final double hmean=results.harmonicMean();
			final double gmean=results.geometricMean();
			final double div=hmean;
//			long expected=(150-31+1)*llw.maxReads;
//			if(!llw.synth){expected=(long)mean;}

			final long median=results.median();

			final long min=results.min();
			final long max=results.max();
			final long p05=results.percentile(0.05);
			final long p95=results.percentile(0.95);
			double stdev=results.stdev();
			double avgDif=results.avgDif(mean);
//			double rmsDif=results.rmsDif(mean);

			final double range=(p95-p05)/div;

//			System.err.println(stdev);
//			System.err.println(rmsDif);

			stdev=stdev/div;
			avgDif=avgDif/div;
//			rmsDif=rmsDif/mean;

			/* Testing indicates that more buckets and fewer trials is more accurate. */
			/* For example, 8 buckets 256 trials < 32 buckets 64 trials < 256 buckets 8 trials, */
			/* from the standpoint of minimizing the standard deviation of the hmean of the estimates over 45 reps. */

			t.stopAndPrint();
			System.err.println("#Trials\tkmers\thmean\tmedian\tmin\tmax\t5%ile\t95%ile\trange\tstdev\tavgDif");
			System.err.println(llw.trials+"\t"+kmers+"\t"+String.format("%.2f", hmean)+"\t"+median+"\t"+min+"\t"+max
					+"\t"+p05+"\t"+p95+"\t"+String.format("%.5f", range)+"\t"+String.format("%.5f", stdev)+"\t"+String.format("%.5f", avgDif));

			if(CardinalityTracker.trackCounts){
				System.err.println("Avg Count:\t"+String.format("%.4f", llw.countSum/(llw.trials*(double)llw.buckets)));
			}
//			System.err.println("#Trials\tkmers\tmean\thmean\tgmean\tmedian\tmin\tmax\t5%ile\t95%ile\trange\tstdev\tavgDif");
//			System.err.println(llw.trials+"\t"+kmers+"\t"+String.format("%.2f", mean)+"\t"+String.format("%.2f", hmean)+"\t"+String.format("%.2f", gmean)+"\t"+median+"\t"+min+"\t"+max
//					+"\t"+p05+"\t"+p95+"\t"+String.format("%.5f", range)+"\t"+String.format("%.5f", stdev)+"\t"+String.format("%.5f", avgDif));

		}

		Read.VALIDATE_IN_CONSTRUCTOR=vic;
	}
	
	/**
	 * Parses command-line arguments to configure buckets, k-mer length, seeds,
	 * minimum probability, file inputs/outputs, and trial parameters.
	 * @param args Command-line arguments
	 */
	public LogLogWrapper(String[] args){

		Shared.capBufferLen(200);
		Shared.capBuffers(8);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		final Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			final String a=split[0].toLowerCase();
			final String b=split.length>1 ? split[1] : null;

			if(a.equals("buckets") || a.equals("loglogbuckets")){
				buckets=Parse.parseIntKMG(b);
			}else if(a.equals("k") || a.equals("loglogk")){
				k=Integer.parseInt(b);
			}else if(a.equals("seed") || a.equals("loglogseed")){
				parser.parseCardinality(arg, "loglogseed", b);
			}else if(a.equals("seed2")){
				seed2=Long.parseLong(b);
			}else if(a.equals("minprob") || a.equals("loglogminprob")){
				minProb=Float.parseFloat(b);
			}else if(a.equals("synth")){
				synth=Parse.parseBoolean(b);
			}else if(a.equals("trials")){
				trials=Parse.parseIntKMG(b);
			}else if(a.equals("printcounts")){
				printCounts=Parse.parseBoolean(b);
			}else if(a.equals("printhist")){
				printHist=Parse.parseBoolean(b);
			}else if(a.equals("khist") || a.equals("hist")){
				khistFile=b;
			}else if(a.equals("peaks") || a.equals("peaksout")){
				peaksFile=b;
			}else if(a.equals("gchist")){
				gcHist=Parse.parseBoolean(b);
			}else if(a.equals("histmax") || a.equals("histlen")){
				histMax=Parse.parseIntKMG(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("loglogcounts") || a.equals("loglogcount") ||
					a.equals("count") || a.equals("counts") || a.equals("trackcounts")){
				CardinalityTracker.trackCounts=Parse.parseBoolean(b);
			}else if(a.equals("clamptoadded") || a.equals("clamp")){
				CardinalityTracker.clampToAdded=Parse.parseBoolean(b);
			}else if(a.equals("cf") || a.equals("loglogcf")){
				CorrectionFactor.USE_CORRECTION=Parse.parseBoolean(b);
			}else if(a.equals("promotethreshold") || a.equals("pt")){
				final int pt=Integer.parseInt(b);
				DynamicLogLog3.PROMOTE_THRESHOLD=pt;
				DynamicLogLog4.PROMOTE_THRESHOLD=pt;
			}else if(a.equals("frozenhistory") || a.equals("frozen")){
//				UltraLogLog8.FROZEN_HISTORY=Parse.parseBoolean(b);
			}else if(a.equals("atomic")){
				assert(false) : "Atomic flag disabled.";
//				CardinalityTracker.atomic=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		{//Process parser fields
			Parser.processQuality();

			maxReads=parser.maxReads;

			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			in1=(parser.in1==null ? null : parser.in1.split(","));
			in2=(parser.in2==null ? null : parser.in2.split(","));
			out=parser.out1;
			seed=parser.loglogseed;
		}

		if(khistFile!=null || peaksFile!=null){CardinalityTracker.trackCounts=true;}

		assert(synth || (in1!=null && in1.length>0)) : "No primary input file specified.";
		if(synth){
			ffin1=ffin2=null;
		}else{
			ffin1=new FileFormat[in1.length];
			ffin2=new FileFormat[in1.length];

			for(int i=0; i<in1.length; i++){
				String a=in1[i];
				String b=(in2!=null && in2.length>i ? in2[i] : null);
				assert(a!=null) : "Null input filename.";
				if(b==null && a.indexOf('#')>-1 && !new File(a).exists()){
					b=a.replace("#", "2");
					a=a.replace("#", "1");
				}

				ffin1[i]=FileFormat.testInput(a, FileFormat.FASTQ, null, true, true);
				ffin2[i]=FileFormat.testInput(b, FileFormat.FASTQ, null, true, true);
			}
		}

		threads=Shared.threads();
		assert(FastaReadInputStream.settingsOK());
	}
	

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Executes a single cardinality estimation trial.
	 * Creates tracker instances, processes input sequences or synthetic data,
	 * and returns the estimated cardinality.
	 * @return Estimated cardinality from the LogLog algorithm
	 */
	long process(){
		final Timer t=new Timer();
		readsProcessed=basesProcessed=kmersProcessed=0;

		final CardinalityTracker log=CardinalityTracker.makeTracker(buckets, k, seed, minProb);

		for(int ffnum=0, max=(synth ? 1 : ffin1.length); ffnum<max; ffnum++){
			Streamer cris=null;
			if(!synth){
				cris=StreamerFactory.getReadInputStream(maxReads, false, ffin1[ffnum], ffin2[ffnum], -1);
				cris.start();
			}

			final LogLogThread[] threadArray=new LogLogThread[threads];
			for(int tid=0; tid<threadArray.length; tid++){
				threadArray[tid]=new LogLogThread(
					CardinalityTracker.makeTracker(buckets, k, seed, minProb), cris, tid);
			}
			for(LogLogThread llt : threadArray){llt.start();}
			for(LogLogThread llt : threadArray){
				while(llt.getState()!=Thread.State.TERMINATED){
					try{
						llt.join();
					}catch(InterruptedException e){
						e.printStackTrace();
					}
				}
				readsProcessed+=llt.readsProcessedT;
				basesProcessed+=llt.basesProcessedT;
				kmersProcessed+=llt.kmersProcessedT;
				log.add(llt.log);
			}

			if(cris!=null){errorState|=ReadWrite.closeStreams(cris);}
		}

//		final int[] max=new int[buckets];
//		if(CardinalityTracker.atomic){
//			for(int i=0; i<log.maxArray.length(); i++){
//				//System.err.println(log.maxArray.get(i));
//				max[i]=log.maxArray.get(i);
//			}
//		}

		t.stop();

		final long cardinality=log.cardinality();
		countSum+=(CardinalityTracker.trackCounts ? log.countSum() : 0);

		if(out!=null){
			ReadWrite.writeString(cardinality+"\n", out);
		}

		if(khistFile!=null || peaksFile!=null){
			makeKhist(log, cardinality);
		}

		if(!Parser.silent){
//			Arrays.sort(copy);
//			System.err.println("Median:        "+copy[Tools.median(copy)]);
//			System.err.println("Mean:          "+Tools.mean(copy));
//			System.err.println("Harmonic Mean: "+Tools.harmonicMean(copy));
			if(log.getClass()==DynamicDemiLog.class){
				final DynamicDemiLog ddl=(DynamicDemiLog)log;
				if(ddl.branch1>0){
					System.err.println("Branch1 Rate:  "+ddl.branch1Rate());
					System.err.println("Branch2 Rate:  "+ddl.branch2Rate());
				}
			}else if(log.getClass()==DynamicDemiLog8.class){
				final DynamicDemiLog8 ddl=(DynamicDemiLog8)log;
				if(ddl.branch1>0){
					System.err.println("Branch1 Rate:  "+ddl.branch1Rate());
					System.err.println("Branch2 Rate:  "+ddl.branch2Rate());
				}
			}else if(log.getClass()==DynamicLogLog4.class){
				final DynamicLogLog4 ddl=(DynamicLogLog4)log;
				if(ddl.branch1>0){
					System.err.println("Branch1 Rate:  "+ddl.branch1Rate());
					System.err.println("Branch2 Rate:  "+ddl.branch2Rate());
				}
			}//assert(false) : log.getClass();
			System.err.println("Cardinality:   "+cardinality);
//			System.err.println("CardinalityH:  "+log.cardinalityH());
//			for(long i : log.counts){System.err.println(i);}

			System.err.println("Time: \t"+t);
			if(printCounts){
				final char[] counts=log.counts16();
				Arrays.sort(counts);
				for(int i=0; i<counts.length; i++){
					System.err.print(((int)counts[i])+" ");
				}
			}
			if(printHist){
				final IntList hist=new IntList();
				final ByteBuilder bb=new ByteBuilder();
				final char[] counts=log.counts16();
				for(char c : counts){hist.increment(c);}
				for(int i=0; i<hist.size; i++){
					final int count=hist.get(i);
					if(count>0){bb.append(i).tab().append(count).nl();}
				}
				System.err.println(bb);
			}
		}

		return cardinality;
	}
	
	private void makeKhist(CardinalityTracker log, long cardinality){
		final char[] counts=log.counts16();
		if(counts==null){return;}
		final byte[] gc=log.gcArray();
		final double scale=cardinality/(double)log.buckets;

		final long[] hist=new long[histMax+1];
		final long[] gcSum=(gc!=null ? new long[histMax+1] : null);
		for(int i=0; i<counts.length; i++){
			final int c=counts[i];
			if(c>0 && c<=histMax){
				hist[c]++;
				if(gcSum!=null){gcSum[c]+=(gc[i]&0xFF);}
			}else if(c>histMax){
				hist[histMax]++;
				if(gcSum!=null){gcSum[histMax]+=(gc[i]&0xFF);}
			}
		}

		final long[] scaled=new long[hist.length];
		for(int i=0; i<hist.length; i++){
			scaled[i]=Math.round(hist[i]*scale);
		}

		if(khistFile!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(khistFile, overwrite, false, true);
			bsw.start();
			bsw.print(gcSum!=null ? "#Depth\tCount\tScaled\tGC%\n" : "#Depth\tCount\tScaled\n");
			for(int i=1; i<hist.length; i++){
				if(hist[i]>0){
					bsw.print(i);
					bsw.print('\t');
					bsw.print(hist[i]);
					bsw.print('\t');
					bsw.print(scaled[i]);
					if(gcSum!=null){
						bsw.print('\t');
						bsw.print(String.format("%.2f", 100.0*gcSum[i]/(hist[i]*k)));
					}
					bsw.print('\n');
				}
			}
			bsw.poisonAndWait();
		}

		if(peaksFile!=null){
			CallPeaks.printClass=false;
			CallPeaks.printPeaks(scaled, null, peaksFile, overwrite,
					2, 2, 2, 2, histMax, 12, k, 0, false, 0, null);
		}
	}

	/**
	 * Generates a synthetic DNA read of specified length with random bases.
	 * Reuses the provided Read object to minimize allocation.
	 * @param len Length of the synthetic read
	 * @param randy Random number generator
	 * @param r Existing Read to reuse, or null to create new
	 * @return Read containing synthetic DNA sequence
	 */
	private static Read makeRead(final int len, final Random randy, Read r){
		if(r==null || r.bases==null || r.bases.length!=len){
			r=new Read(null, null, 0);
			r.bases=new byte[len];
		}
		final byte[] bases=r.bases;

		int pos=0;
		final int basesPerRand=4;//Fewer calls to rand should be faster
		for(int max=bases.length-(bases.length%basesPerRand); pos<max;){
			int x=randy.nextInt()%prime;
			for(int i=0; i<basesPerRand; i++){
				final int num=x&3;
				final byte b=AminoAcid.numberToBase[num];
				bases[pos]=b;
				pos++;
				x>>=2;
			}
		}
		for(; pos<bases.length; pos++){
			final int x=randy.nextInt()%prime;
			final int num=x&3;
			final byte b=AminoAcid.numberToBase[num];
			bases[pos]=b;
		}
		return r;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Worker thread for parallel cardinality estimation.
	 * Handles both real sequence data and synthetic generation. */
	private class LogLogThread extends Thread{

		LogLogThread(CardinalityTracker log_, Streamer cris_, int tid_){
			log=log_;
			cris=cris_;
			tid=tid_;
		}

		@Override
		public void run(){
			if(cris!=null){runCris();}
			else{runSynth();}
		}

		/** Processes real sequence data from the input stream. */
		public void runCris(){
			final int kt=k;
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

				for(Read r : reads){
//					if(!r.validated()){r.validate(true);}
//					if(r.mate!=null && !r.mate.validated()){r.mate.validate(true);}
					log.hash(r);
					readsProcessedT+=r.pairCount();
					basesProcessedT+=r.pairLength();
					kmersProcessedT+=r.numPairKmers(kt);
				}

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln);
		}

		/** Generates and processes synthetic 150bp reads. */
		public void runSynth(){
			final int kt=k;
			assert(maxReads>0 && maxReads<Long.MAX_VALUE);
			long readsLeft=maxReads/threads;
			readsLeft+=(maxReads%threads>tid ? 1 : 0);
			final Random randy=Shared.threadLocalRandom(seed2<0 ? seed2 : seed2+999999L);

			Read r=null;
			while(readsLeft>0){
				r=makeRead(150, randy, r);
				log.hash(r);
				readsProcessedT+=r.pairCount();
				basesProcessedT+=r.pairLength();
				kmersProcessedT+=r.numPairKmers(kt);
				readsLeft--;
			}
		}

		private final CardinalityTracker log;
		private final Streamer cris;
		private final int tid;

		protected long readsProcessedT=0;
		protected long basesProcessedT=0;
		protected long kmersProcessedT=0;

	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Number of buckets for the LogLog cardinality tracker */
	private int buckets=2048;
	/** K-mer length for hashing sequences */
	private int k=31;
	/** Random seed for hash function initialization */
	private long seed=0;
	/** Secondary random seed for synthetic data generation */
	private long seed2=-1;
	/** Minimum probability threshold for k-mer inclusion */
	private float minProb=0;
	/** Cumulative sum of bucket counts across all trials */
	private long countSum=0;
	
	private String[] in1=null;
	private String[] in2=null;
	private String out=null;
	
	/*--------------------------------------------------------------*/
	
	protected long readsProcessed=0;
	protected long basesProcessed=0;
	protected long kmersProcessed=0;
	
	private long maxReads=-1;
	
	boolean overwrite=true;
	boolean append=false;
	boolean errorState=false;

	/** Number of cardinality estimation trials to run */
	private int trials=1;
	/** Whether to use synthetic data generation instead of real input files */
	private boolean synth=false;
	private boolean printCounts=false;
	private boolean printHist=false;
	private String khistFile=null;
	private String peaksFile=null;
	private boolean gcHist=false;
	private int histMax=100000;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat[] ffin1;
	private final FileFormat[] ffin2;
	
	/** Number of worker threads for parallel processing */
	final int threads;
	/** Prime number used for modular arithmetic in synthetic data generation */
	private static final int prime=32452843; //A prime number
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private static PrintStream outstream=System.err;
	public static boolean verbose=false;
}
