package kmer;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.concurrent.atomic.AtomicLong;

import bloom.KCountArray;
import bloom.KmerCountAbstract;
import bloom.ReadCounter;
import fileIO.ByteStreamWriter;
import jgi.CallPeaks;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import structures.IntList;
import ukmer.Kmer;

/**
 * Loads and holds kmers for Tadpole
 * @author Brian Bushnell
 * @date Jun 22, 2015
 *
 */
public abstract class AbstractKmerTableSet {
	
	public static final boolean isValidArgument(String a){
			if(a.equals("in") || a.equals("in1")){
			}else if(a.equals("in2")){
			}else if(a.equals("append") || a.equals("app")){
			}else if(a.equals("overwrite") || a.equals("ow")){
			}else if(a.equals("initialsize")){
			}else if(a.equals("showstats") || a.equals("stats")){
			}else if(a.equals("ways")){
			}else if(a.equals("buflen") || a.equals("bufflen") || a.equals("bufferlength")){
			}else if(a.equals("k")){
			}else if(a.equals("threads") || a.equals("t")){
			}else if(a.equals("showspeed") || a.equals("ss")){
			}else if(a.equals("ecco")){
			}else if(a.equals("merge")){
			}else if(a.equals("verbose")){
			}else if(a.equals("verbose2")){
			}else if(a.equals("minprob")){
			}else if(a.equals("reads") || a.startsWith("maxreads")){
			}else if(a.equals("prealloc") || a.equals("preallocate")){
			}else if(a.equals("prefilter")){
			}else if(a.equals("prefiltersize") || a.equals("prefilterfraction") || a.equals("pff")){
			}else if(a.equals("minprobprefilter") || a.equals("mpp")){
			}else if(a.equals("minprobmain") || a.equals("mpm")){
			}else if(a.equals("prefilterpasses") || a.equals("prepasses")){
			}else if(a.equals("prehashes") || a.equals("hashes")){
			}else if(a.equals("onepass")){
			}else if(a.equals("passes")){
			}else if(a.equals("rcomp")){
			}else if(a.equals("maskmiddle")){
			}else if(a.equals("filtermemory") || a.equals("prefiltermemory") || a.equals("filtermem")){
			}else{
				return false;
			}
			return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public final void process(Timer t){
		
		/* Count kmers */
		long added=processInput();
		
		/* Stop timer and calculate speed statistics */
		t.stop();
		
		showStats(t, added);
		
		/* Throw an exception if errors were detected */
		if(errorState){
			throw new RuntimeException(getClass().getSimpleName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	
	public abstract void clear();
	
	
	public final long processInput(){
		
		/* Start phase timer */
		Timer t=new Timer();

//		if(DISPLAY_PROGRESS){
//			outstream.println("Before loading:");
//			Shared.printMemory();
//			outstream.println();
//		}
		
		prefilterArray=makePrefilter(new KCountArray[1], null);
		if(prefilterArray!=null){
			prefilterArray.purgeFilter();
			filterMax2=Tools.min(filterMax, prefilterArray.maxValue-1);
			
			/* This is already getting printed in makePrefilter */
//			if(DISPLAY_PROGRESS){
//				outstream.println("After prefilter:");
//				Shared.printMemory();
//				outstream.println();
//			}
		}
//		assert(false) : prefilterArray.cellBits+", "+prefilterArray.maxValue+", "+filterMax+", "+filterMax2;
		
		if(DISPLAY_STATS){System.err.println("Estimated kmer capacity: \t"+estimatedKmerCapacity());}
		
		assert(!allocated);
		allocateTables();
		allocated=true;
		
		if(DISPLAY_PROGRESS){
			outstream.println("After table allocation:");
			Shared.printMemory();
			outstream.println();
		}
		
		/* Fill tables with kmers */
		long added=loadKmers();
		
		/* Clear prefilter; no longer needed */
		prefilterArray=null;
		
//		long removed=0;
//		if(prefilter && filterMax>0){
//			removed=removeKmersAtMost(filterMax);
//			System.err.println("Removed "+removed+" low-depth kmers.");
//		}
		
		return added;
	}
	
	
	public final KCountArray makePrefilter(final KCountArray[] filter, Timer ht){
//		assert(false) : lastFilter+", "+prefilter+", "+filterMax()+", "+currentPass+", "+filterMemory(currentPass);
		if(!prefilter){return null;}
		assert(!MASK_MIDDLE) : "TODO: MaskMiddle not yet supported by prefilter.";
		
		if(filter[0]!=null){
			filter[0].purgeFilter();
			assert(filter[0].prefilter()==null);
		}
		
		KmerCountAbstract.CANONICAL=true;

		long precells=-1;
		int cbits=1;
		if(onePass){
			while(filterMax>=(1<<cbits)){cbits*=2;}
		}else{
			while(filterMax+1>=(1<<cbits)){cbits*=2;}
		}
		if(prepasses>2 && currentPass==prepasses-1){cbits=1;}
		
		byte minq=0;
		if(precells<1){
			long prebits=(filterMemory(currentPass)-10)*8;
			
//			System.err.println("prebits="+prebits+", currentPass="+currentPass);
			
			precells=prebits/cbits;
			if(precells<100000){ //Not enough memory - no point.
				prefilter=false;
				return null;
			}
		}
		if(prehashes<1){prehashes=2;}

		ReadCounter rc=new ReadCounter(kbig(), true, ecco(), false, Shared.AMINO_IN);
		if(onePass){
			assert(filter==null || filter.length==1) : "Multiple filtering passes are not allowed in onepass mode.\n"+filter.length+","+prepasses+", "+onePass+", "+prefilter;
			filter[0]=rc.makeKca(null, null, null, cbits, precells, prehashes, minq,
					maxReads, 1, 1, 1, null, 0);
		}else{
			if(ht==null){ht=new Timer();}
			ht.start();
			filter[0]=rc.makeKca_als(in1, in2, extra, cbits, precells, prehashes, minq,
					maxReads, 1, 1, 1, filter[0], filterMax);
			assert(filterMax<filter[0].maxValue || (currentPass>0 && currentPass==prepasses-1));
			outstream.println("Made prefilter:   \t"+filter[0].toShortString(prehashes));
			double uf=filter[0].usedFraction();
//			System.err.println("cellsUsed: "+filter[0].cellsUsed(1)+" //123"); //123
			if(uf>0.5){
				outstream.println("Warning:  This table is "+(uf>0.995 ? "totally" : uf>0.99 ? "crazy" : uf>0.95 ? "incredibly" : uf>0.9 ? "extremely" : uf>0.8 ? "very" :
					uf>0.7 ? "rather" : uf>0.6 ? "fairly" : "somewhat")+" full.  Ideal load is under 50% used." +
						"\nFor better accuracy, run on a node with more memory; quality-trim or error-correct reads; or increase prefiltersize.");
			}
			ht.stop();
			currentPass++;
			
			final double kmers=filter[0].estimateUniqueKmers(prehashes, Tools.min(filterMax+1, filter[0].maxValue));
			outstream.println("Estimated valid kmers: \t\t"+(long)kmers);
			
//			outstream.println("Estimated valid kmers 1+: "+(long)filter[0].estimateUniqueKmers(prehashes, 1));
//			outstream.println("Estimated valid kmers 2+: "+(long)filter[0].estimateUniqueKmers(prehashes, 2));
//			outstream.println("Estimated valid kmers 3+: "+(long)filter[0].estimateUniqueKmers(prehashes, 3));
//			outstream.println("Estimated valid kmers 4+: "+(long)filter[0].estimateUniqueKmers(prehashes, 4));
			
			if(prepasses<0){//auto
				if((currentPass&1)==0){
					return makePrefilter(filter, ht);
				}else if(currentPass<5){
					if(kmers>estimatedKmerCapacity()){
						return makePrefilter(filter, ht);
					}
				}
			}else if(currentPass<prepasses){
				return makePrefilter(filter, ht);
			}
			
			if(DISPLAY_PROGRESS){
				outstream.println("Prefilter time:\t"+ht);
				outstream.println("After prefilter:");
				Shared.printMemory();
				outstream.println();
			}
		}
		
		return filter[0];
	}
	
	
	public final void showStats(Timer t, long added){
		
		if(!DISPLAY_STATS){return;}
		
		if(DISPLAY_PROGRESS){
			outstream.println("After loading:");
			Shared.printMemory();
			outstream.println();
		}
		
		t.stop();
		outstream.println("Input:                      \t"+readsIn+" reads \t\t"+basesIn+" bases.");
		outstream.println("Unique Kmers:               \t"+added);
		outstream.println("Load Time:                  \t"+t);
		
		if(showSpeed){
			outstream.println();
			outstream.println(Tools.readsBasesProcessed(t.elapsed, readsIn, basesIn, 8));
		}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private final long loadKmers(){
		//allocateTables();
		assert(allocated);
		kmersLoaded=0;
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=false;
		for(int i=0; i<in1.size(); i++){
			String a=in1.get(i);
			String b=in2.size()>i ? in2.get(i) : null;
			int idx=a.indexOf('#');
			if(idx>=0 && b==null && !new File(a).exists()){
				b=a.replaceFirst("#", "2");
				a=a.replaceFirst("#", "1");
			}
			kmersLoaded+=loadKmers(a, b);
		}
		for(int i=0; i<extra.size(); i++){
			String a=extra.get(i);
			String b=null;
			int idx=a.indexOf('#');
			if(idx>=0 && b==null && !new File(a).exists()){
				b=a.replaceFirst("#", "2");
				a=a.replaceFirst("#", "1");
			}
			kmersLoaded+=loadKmers(a, b);
		}
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		return kmersLoaded;
	}
	
	/**
	 * Load reads into tables, using multiple LoadThread.
	 */
	public abstract long loadKmers(String fname1, String fname2);
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public abstract long regenerate(int limit);
	
	public abstract Object getTable(int tnum);
	
	public abstract long[] fillHistogram(int histMax);

	public abstract void countGC(long[] gcCounts, int max);
	
	public final long[] fillGcCounts(int histMax){
		long[] gcCounts=new long[histMax+1];
		countGC(gcCounts, histMax);
		return gcCounts;
	}
	
	public final float[] makeGcHistogram(long[] counts, long[] gcCounts){
		float[] gcHist=new float[counts.length];
		final long k=kbig();
		for(int i=0; i<counts.length; i++){
			long gc=gcCounts[i];
			double bases=Tools.max(counts[i], 1)*k;
			gcHist[i]=(float)(gc/bases);
		}
		return gcHist;
	}
	
	public abstract void initializeOwnership();
	
	public abstract void clearOwnership();
	
	public abstract int ways();
	
	public final int fillCounts(byte[] bases, IntList counts, Kmer kmer){
		return fillSpecificCounts(bases, counts, null, kmer);
	}
	
	public abstract int fillSpecificCounts(byte[] bases, IntList counts, BitSet positions, Kmer kmer);
	
	public abstract int regenerateCounts(byte[] bases, IntList counts, Kmer kmer, BitSet changed);
	
	/*--------------------------------------------------------------*/
	/*----------------       Printing Methods       ----------------*/
	/*--------------------------------------------------------------*/

	public abstract boolean dumpKmersAsBytes(String fname, int mincount, int maxcount, boolean printTime, AtomicLong remaining);
	public abstract boolean dumpKmersAsBytes_MT(String fname, int mincount, int maxcount, boolean printTime, AtomicLong remaining);
	
	public final long[][] makeKhist(String fname, int cols, int max, boolean printHeader, boolean printZeros, boolean printTime, 
			boolean smooth, boolean calcGC, boolean doLogScale, double logWidth, int logPasses, int smoothRadius){
		Timer t=new Timer();
		
		long[] ca=fillHistogram(max);
		float[] gcHist=null;
		if(calcGC){
//			assert(false) : max+", "+ca.length;
			long[] gc=(calcGC ? fillGcCounts(max) : null);
//			assert(false) : max+", "+ca.length+", "+gc.length;
			gcHist=makeGcHistogram(ca, gc);
		}
		
		long[] logScale=null;
		
		if(smooth){
			ca=CallPeaks.smoothProgressive(ca, smoothRadius);
		}
		if(doLogScale){
			logScale=CallPeaks.logScale(ca, logWidth, 1, logPasses);
		}
		
		long[][] ret=new long[2][];
		ret[0]=ca;
		if(gcHist!=null){
			final int k=kbig();
			ret[1]=new long[ca.length];
			for(int i=1; i<ca.length; i++){
				ret[1][i]=Math.round(ca[i]*gcHist[i]*k);
			}
		}
		
		if(fname==null){return ret;}
		
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, false, true);
		bsw.start();
		if(printHeader){
			bsw.print("#Depth\t"+(cols==3 ? "RawCount\t" : "")+"Count"+(doLogScale ? "\tlogScale" : "")+(calcGC ? "\tGC%\n" : "\n"));
		}
		
		for(int i=1; i<ca.length; i++){
			long count=ca[i];
			if(printZeros || count>0){
				bsw.print(i);
				bsw.print('\t');
				if(cols==3){
					bsw.print(i*count);
					bsw.print('\t');
				}
				bsw.print(count);
				if(doLogScale){
					bsw.print('\t').print(logScale[i]);
				}
				if(gcHist!=null){
					bsw.print(Tools.format("\t%.2f", 100f*gcHist[i]));
				}
				bsw.print('\n');
			}
		}
		bsw.poisonAndWait();
		t.stop();
		if(printTime){outstream.println("Histogram Write Time:       \t"+t);}
		return ret;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean showStats=true;
	
//	public boolean silent=false;
	
	/** Has this class encountered errors while processing? */
	public boolean errorState=false;
	
	/** Use a count-min prefilter for low-depth kmers */
	public boolean prefilter=false;
	/** Fill the prefilter at the same time as the main table */
	public boolean onePass=false;
	
	public boolean amino=false;
	/** Number of hashes used by prefilter */
	public int prehashes=2;
	/** Fraction of memory used by prefilter */
	public double prefilterFraction=0.2;
	
	/** Initial size of data structures */
	public int initialSize=-1;
	/** Fraction of available memory preallocated to arrays */
	public double preallocFraction=1.0;
	
	public KCountArray prefilterArray=null;
	
	public boolean minProbPrefilter=true;
	public boolean minProbMain=true;

	/** Input reads for kmers */
	public ArrayList<String> in1=new ArrayList<String>(), in2=new ArrayList<String>();
	
	/** Extra files for use as kmers */
	public ArrayList<String> extra=new ArrayList<String>();
	
	/** Maximum input reads (or pairs) to process.  Does not apply to references.  -1 means unlimited. */
	public long maxReads=-1;
	
	public int buflen=1000;
	
	/** Filter kmers up to this level; don't store them in primary data structure */
	protected int filterMax=0;
	protected int filterMax2=0;
	
	public long readsIn=0;
	public long basesIn=0;
	public long lowqReads=0;
	public long lowqBases=0;
	public long readsTrimmed=0;
	public long basesTrimmed=0;
	
	public long kmersIn=0;
	public long kmersLoaded=0;
	
	private int currentPass=0;
	protected int prepasses=1;
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/
	
	public abstract int kbig();
	public abstract long filterMemory(int pass);
	public abstract long tableMemory();
	public abstract long estimatedKmerCapacity();
	public abstract boolean ecco();
	public abstract boolean qtrimLeft();
	public abstract boolean qtrimRight();
	public abstract float minAvgQuality();
	public final int filterMax(){return filterMax;}
	public abstract boolean rcomp();
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	protected abstract void allocateTables();
	protected boolean allocated=false;

	/** Print messages to this stream */
	public static PrintStream outstream=System.err;
	/** Permission to overwrite existing files */
	public static boolean overwrite=true;
	/** Permission to append to existing files */
	public static boolean append=false;
	/** Print speed statistics upon completion */
	public static boolean showSpeed=true;
	/** Display progress messages such as memory usage */
	public static boolean DISPLAY_PROGRESS=true;
	/** Display kmer loading information */
	public static boolean DISPLAY_STATS=true;
	/** Verbose messages */
	public static boolean verbose=false;
	/** Debugging verbose messages */
	public static boolean verbose2=false;
	/** Number of ProcessThreads */
	public static int THREADS=Shared.threads();
	
	public static int maxNs=Integer.MAX_VALUE;
	public static int minLen=0;
	
	/** Increment owner by this much to indicate claim is final. */
	public static final int CLAIM_OFFSET=100000;
	
	/** Default initial table size */
	public static final int initialSizeDefault=128000;
	
	public static final float[] PROB_CORRECT=Arrays.copyOf(align2.QualityTools.PROB_CORRECT, 127);
	public static final float[] PROB_CORRECT_INVERSE=Arrays.copyOf(align2.QualityTools.PROB_CORRECT_INVERSE, 127);
	
	public static boolean IGNORE_UNKNOWN_ARGS=true;

	public static final int NOT_PRESENT=AbstractKmerTable.NOT_PRESENT, HASH_COLLISION=AbstractKmerTable.HASH_COLLISION;
	public static final int NO_OWNER=AbstractKmerTable.NO_OWNER;
	
	public static double defaultMinprob=0;
	
	/** IIRC this has to do with allowing 32-mers */
	public static boolean MASK_CORE=false;
	public static boolean MASK_MIDDLE=false;
	public static boolean FAST_FILL=true;
	
}
