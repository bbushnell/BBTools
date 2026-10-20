package clump;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import bloom.KCountArray;
import bloom.ReadCounter;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Counts one selected pivot per read, including linked mates, into a shared
 * approximate count array. The minimum count controls cell width and later pivot
 * selection; it does not filter reads during this counting pass. Short or
 * ambiguous reads use KmerComparator's synthesized fallback pivot.
 * Configuration mutates shared parser/comparator settings and must be serialized.
 *
 * @author Brian Bushnell
 * @date August 19, 2016
 */
public class PivotSet{

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Program entry point that creates a k-mer count array from input reads.
	 * @param args Command line arguments for input files and parameters */
	public static void main(String[] args){
		makeSet(args);
	}
	
	/** Builds and finalizes a table from command-line input settings. Restores only
	 * the saved pigz/unpigz flags on normal return; returns the shared table object. */
	public static KCountArray makeSet(String[] args){
		final boolean pigz=ReadWrite.USE_PIGZ, unpigz=ReadWrite.USE_UNPIGZ;
		Timer t=new Timer();
		PivotSet x=new PivotSet(args);
		KCountArray kca=x.process(t, false);
		ReadWrite.USE_PIGZ=pigz;
		ReadWrite.USE_UNPIGZ=unpigz;
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
		
		return kca;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parses counting/input settings and resolves file descriptors without starting readers. */
	public PivotSet(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=KmerComparator.verbose=Parse.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
				assert(k>0 && k<32);
			}else if(a.equals("ecco")){
				ecco=Parse.parseBoolean(b);
			}else if(a.equals("rename") || a.equals("addname")){
				//do nothing
			}else if(a.equals("rcomp") || a.equals("reversecomplement")){
				//do nothing
			}else if(a.equals("condense") || a.equals("consensus")){
				//do nothing
			//Fixed [clump/PivotSet#001]: accept mincr like the other clump tools. Consensus is handled above;
			//the previous duplicate consensus branch left mincr unrecognized. ClumpTools uses minCount.
			}else if(a.equals("mincount") || a.equals("mincr")){
				minCount=Integer.parseInt(b);
			}else if(a.equals("correct") || a.equals("ecc")){
				//do nothing
			}else if(a.equals("groups") || a.equals("g") || a.equals("sets") || a.equals("ways")){
				//do nothing
			}else if(a.equals("seed")){
				KmerComparator.defaultSeed=Long.parseLong(b);
			}else if(a.equals("hashes")){
				KmerComparator.setHashes(Integer.parseInt(b));
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			in1=parser.in1;
			in2=parser.in2;
			
			extin=parser.extin;
		}
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Estimates the number of cbits-wide cells fitting a fraction of usable heap. */
	private static long getCells(double fraction, int cbits){
		final long memory=Runtime.getRuntime().maxMemory();
		final long usable=(long)Tools.max(((memory-96000000)*.73), memory*0.45);
		final double filterMem=usable*fraction;
		return (long)((filterMem*8)/cbits);
	}
	
	/**
	 * Creates read streams and processes all data to generate k-mer count array.
	 * Configures memory allocation, initializes concurrent processing threads,
	 * reports counts, closes input, and returns the finalized count array.
	 *
	 * @param t Timer for tracking processing duration
	 * @param amino Passed to the count-array builder; pivot selection still uses
	 * KmerComparator's nucleotide encoding
	 * @return Approximate counts of selected pivots, including linked mates
	 */
	public KCountArray process(Timer t, boolean amino){
		int cbits=2;
		while((1L<<cbits)<=minCount){cbits*=2;}
		int filterHashes=2;
		float fraction=0.1f;
		long cells=getCells(fraction, cbits);
		ReadCounter rc=new ReadCounter(k, true, ecco, false, amino);
		KCountArray kca=rc.makeKca(null, null, null, cbits, cells, filterHashes, 0, maxReads, 1, 1, 1, null, 0);
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			cris.start();
			if(verbose){outstream.println("Started cris");}
		}
		
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the read stream
		processInner(cris, kca);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		errorState|=ReadStats.writeAll();
		errorState|=ReadWrite.closeStreams(cris);
		
		t.stop();

		outstream.println("Made filter:     \t"+kca.toShortString(filterHashes));
		outstream.println("Estimated pivots:      \t"+(long)kca.estimateUniqueKmers(filterHashes));
		outstream.println("Estimated pivots >1x:  \t"+(long)kca.estimateUniqueKmers(filterHashes, minCount));
		
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		if(errorState){
			Clumpify.sharedErrorState=true;
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
		return kca;
	}
	
	/** Counts pivots from an already-started stream and finalizes the returned array.
	 * Does not close cris or report its error state; the caller owns that lifecycle.
	 * This entry point does not perform overlap correction or aggregate statistics.
	 * @param cris Borrowed stream, consumed through its terminal list
	 * @param k Nucleotide pivot length
	 * @param minCount Threshold used to select counting-cell width
	 * @param amino Count-array builder mode; does not change pivot encoding
	 * @return Finalized approximate pivot-count table
	 */
	public static KCountArray makeKcaStatic(final ConcurrentReadInputStream cris, int k, int minCount, boolean amino){

		KmerComparator kc=new KmerComparator(k, false, false);
		int cbits=2;
		while((1L<<cbits)<=minCount){cbits*=2;}
		int filterHashes=2;
		float fraction=0.1f;
		long cells=getCells(fraction, cbits);
		ReadCounter rc=new ReadCounter(k, true, false, false, amino);
		KCountArray kca=rc.makeKca(null, null, null, cbits, cells, filterHashes, 0, -1, 1, 1, 1, null, 0);
		
		if(verbose){System.err.println("Making hash threads.");}
		final int threads=Shared.threads();
		ArrayList<HashThread> alht=new ArrayList<HashThread>(threads);
		for(int i=0; i<threads; i++){alht.add(new HashThread(cris, kc, kca, false));}
		
		if(verbose){System.err.println("Starting threads.");}
		for(HashThread ht : alht){ht.start();}
		
		if(verbose){System.err.println("Waiting for threads.");}
		/* Wait for threads to die */
		for(HashThread ht : alht){
			
			/* Wait for a thread to die */
			while(ht.getState()!=Thread.State.TERMINATED){
				try{
					ht.join();
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
		}
		kca.shutdown();
		return kca;
	}
	
	/** Recounts expanded groups with the current pivot policy, keeping legacy input.
	 * Owns each CRIS and joins every started counting worker before closing it.
	 * All groups accumulate into one table, finalized once after successful counting.
	 * @param first Nonempty primary group descriptors
	 * @param second Aligned optional mate descriptors, or null
	 * @param k Nucleotide pivot length
	 * @param minCount Threshold determining counting-cell width
	 * @param amino Count-array builder mode, as in the legacy single-stream helper
	 * @return One completed approximate table for all current groups */
	public static KCountArray makeKcaStatic(FileFormat[] first, FileFormat[] second,
			int k, int minCount, boolean amino){
		if(first==null || first.length==0 || (second!=null && second.length!=first.length)){
			throw new IllegalArgumentException("Grouped pivot counting requires aligned input arrays");
		}
		for(FileFormat input : first){
			if(input==null){throw new IllegalArgumentException("A pivot-count group cannot be null");}
		}
		assert(k>0 && k<32) : "KmerComparator packs nucleotide pivots into 62 bits; k="+k;
		assert(minCount>0) : "A pivot-count threshold must be positive; minCount="+minCount;
		final KmerComparator comparator=new KmerComparator(k, false, false);
		int bits=2;
		while((1L<<bits)<=minCount){bits*=2;}
		final long cells=getCells(0.1f, bits);
		final ReadCounter counter=new ReadCounter(k, true, false, false, amino);
		final KCountArray counts=counter.makeKca(null, null, null, bits, cells, 2,
				0, -1, 1, 1, 1, null, 0);
		for(int group=0; group<first.length; group++){
			ConcurrentReadInputStream input=null;
			Throwable failure=null;
			boolean interrupted=false;
			try{
				input=ConcurrentReadInputStream.getReadInputStream(-1, true, first[group],
						second==null ? null : second[group], null, null);
				input.start();
				final GroupHashThread[] workers=new GroupHashThread[Tools.max(1, Shared.threads())];
				for(int i=0; i<workers.length; i++){
					workers[i]=new GroupHashThread(input, comparator, counts);
				}
				int started=0;
				try{
					for(; started<workers.length; started++){workers[started].start();}
				}catch(RuntimeException | Error e){failure=e;}
				for(int i=0; i<started; i++){
					while(workers[i].isAlive()){
						try{workers[i].join();}catch(InterruptedException e){
							interrupted=true;
							failure=appendFailure(failure, e);
						}
					}
					failure=appendFailure(failure, workers[i].failure);
				}
			}catch(RuntimeException | Error e){failure=appendFailure(failure, e);}
			finally{
				try{
					if(input!=null && ReadWrite.closeStream(input)){
						failure=appendFailure(failure,
								new IllegalStateException("Error recounting pivot group "+first[group].name()));
					}
				}catch(RuntimeException | Error e){failure=appendFailure(failure, e);}
				if(interrupted){Thread.currentThread().interrupt();}
			}
			if(failure instanceof Error){throw (Error)failure;}
			if(failure!=null){throw new RuntimeException("Grouped pivot recount failed", failure);}
		}
		counts.shutdown();
		return counts;
	}

	/** Keeps the first failure while preserving later independent cleanup errors. */
	private static Throwable appendFailure(Throwable first, Throwable next){
		if(first==null){return next;}
		if(next!=null && next!=first){first.addSuppressed(next);}
		return first;
	}

	/**
	 * Manages thread creation and coordination for k-mer processing.
	 * Creates hash threads, joins them, aggregates mate-inclusive statistics and
	 * finalizes kca. The caller retains responsibility for closing cris.
	 * @param cris Already-started input stream providing reads for processing
	 * @param kca K-mer count array to populate with frequency data
	 */
	public void processInner(final ConcurrentReadInputStream cris, KCountArray kca){
		if(verbose){outstream.println("Making comparator.");}
		KmerComparator kc=new KmerComparator(k, false, false);
		
		if(verbose){outstream.println("Making hash threads.");}
		final int threads=Shared.threads();
		ArrayList<HashThread> alht=new ArrayList<HashThread>(threads);
		for(int i=0; i<threads; i++){alht.add(new HashThread(cris, kc, kca, ecco));}
		
		if(verbose){outstream.println("Starting threads.");}
		for(HashThread ht : alht){ht.start();}
		
		if(verbose){outstream.println("Waiting for threads.");}
		/* Wait for threads to die */
		for(HashThread ht : alht){
			
			/* Wait for a thread to die */
			while(ht.getState()!=Thread.State.TERMINATED){
				try{
					ht.join();
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
			readsProcessed+=ht.readsProcessedT;
			basesProcessed+=ht.basesProcessedT;
		}
		kca.shutdown();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Captures failures for the grouped caller without changing legacy worker entry points. */
	private static final class GroupHashThread extends HashThread{
		GroupHashThread(ConcurrentReadInputStream input, KmerComparator comparator, KCountArray counts){
			super(input, comparator, counts, false);
		}
		@Override
		public void run(){
			try{super.run();}catch(Throwable error){failure=error;}
		}
		/** Published to the grouped caller through Thread.join. */
		Throwable failure;
	}

	/** Counts both mates while owning each borrowed input batch. */
	private static class HashThread extends Thread{

		HashThread(ConcurrentReadInputStream cris_, KmerComparator kc_, KCountArray kca_, boolean ecco_){
			cris=cris_;
			kc=kc_;
			kca=kca_;
			ecco=ecco_;
		}

		/**
		 * Main thread execution method that processes reads and updates k-mer counts.
		 * Retrieves read lists, processes each read pair, optionally applies error correction,
		 * hashes sequences, and increments corresponding k-mer counters.
		 */
		@Override
		public void run(){
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				for(Read r1 : reads){
					Read r2=r1.mate;
					readsProcessedT+=r1.pairCount();
					basesProcessedT+=r1.pairLength();
					if(ecco && r2!=null){
						if(r2!=null){BBMerge.findOverlapStrict(r1, r2, true);}//redundant inner null-check (already guarded above); harmless
					}
					{
						//Count the selected pivot, including fillShort's synthesized pivot for short or all-N reads.
						//The nonnegative guard is retained; short reads are not skipped by this counting pass.
						final long kmer=kc.hash(r1, null, 0, false);
						if(kmer>=0){
							kca.increment(kmer);
						}
					}
					if(r2!=null){
						final long kmer=kc.hash(r2, null, 0, false);
						if(kmer>=0){
							kca.increment(kmer);
						}
					}
				}
				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		final ConcurrentReadInputStream cris;
		final KmerComparator kc;
		final KCountArray kca;
		final boolean ecco;
		
		protected long readsProcessedT=0;
		protected long basesProcessedT=0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Nucleotide pivot length, limited to 1 through 31 by the packed representation. */
	private int k=31;
	/** Threshold used to choose counting-cell width and report estimated abundant pivots. */
	private int minCount=2;
	
	/*--------------------------------------------------------------*/
	/*----------------          I/O Fields          ----------------*/
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String in2=null;
	
	private String extin=null;
	
	/*--------------------------------------------------------------*/
	
	/** Input reads counted, including linked mates. */
	protected long readsProcessed=0;
	/** Input bases counted, including linked mates. */
	protected long basesProcessed=0;
	
	/** Reader fragment limit; negative means no configured limit. */
	private long maxReads=-1;
	/** Whether the file-processing path corrects overlapping pairs before counting. */
	private boolean ecco=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	
}
