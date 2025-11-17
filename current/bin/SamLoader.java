package bin;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicLongArray;
import java.util.concurrent.locks.ReadWriteLock;

import fileIO.FileFormat;
import shared.Shared;
import shared.Tools;
import stream.SamLine;
import stream.Streamer;
import stream.StreamerFactory;
import structures.IntHashMap;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.EntropyTracker;

/**
 * Loads SAM/BAM files to calculate contig coverage depth and build connectivity graphs.
 * Processes multiple SAM files concurrently using thread-based parallelization.
 * Accumulates alignment statistics and creates contig-to-contig linkage information
 * from paired-end reads mapping to different contigs.
 *
 * @author Brian Bushnell
 * @date April 15, 2015
 */
public class SamLoader implements Accumulator<SamLoader.LoadThread> {
	
	/** Constructs a SamLoader with specified output stream.
	 * @param outstream_ Stream for progress and diagnostic output */
	public SamLoader(PrintStream outstream_) {
		outstream=outstream_;
	}
	
	/**
	 * Loads SAM files using deprecated interface.
	 * Converts contig map to sorted list and delegates to main load method.
	 * @param fnames List of SAM file names to process
	 * @param contigMap Map of contig names to Contig objects
	 * @param graph Array of maps for storing contig connectivity
	 * @deprecated Use load(ArrayList, HashMap, ArrayList, IntHashMap[]) instead
	 */
	@Deprecated
	public void load(ArrayList<String> fnames, HashMap<String, Contig> contigMap, IntHashMap[] graph) {
		//Contig list should already be sorted and numbered.
		ArrayList<Contig> list=new ArrayList<Contig>(contigMap.values());
		Collections.sort(list);
		for(int i=0; i<list.size(); i++) {list.get(i).setID(i);}
		load(fnames, contigMap, list, graph);
	}
	
	/** Spawn process threads */
	public void load(ArrayList<String> fnames, HashMap<String, Contig> contigMap, 
			ArrayList<Contig> contigs, IntHashMap[] graph){
		final int files=fnames.size();
		SamLine.RNAME_AS_BYTES=false;
		//Do anything necessary prior to processing
		
		FileFormat ff0=FileFormat.testInput(fnames.get(0), FileFormat.SAM, null, false, false);
		final boolean compressed=ff0.compressed();
		
		//Determine how many threads may be used
		final int availableThreads=Tools.max(1, Shared.threads());
		final int maxThreadsPerFile=(availableThreads+files-1)/files;
		final float[] ideal;
		{//FileRead, Decompress, SamLine, Coverage
			float x=MAX_SAM_LOADER_THREADS_PER_FILE;
			if(ff0.bam()) {//bam
				ideal=new float[] {1f, 4f, 6f, x};
			}else if(ff0.compressed()) {//sam.bgz or bz2
				ideal=new float[] {1f, 4f, 6f, x};
			}else if(false) {//sam.gz, non-bgzip.  Detectable from magic word...
				ideal=new float[] {0.5f, 1f, 2f, 0.5f*x};
			}else {//sam
				assert(ff0.sam());
				ideal=new float[] {1f, 0f, 6f, x};
			}
		}
		final int[] allocation=allocateThreads(ideal, availableThreads);
//		final int availableLoaderThreads=Tools.mid(1, MAX_SAM_LOADER_THREADS, availableThreads);
//		final int maxLoaderThreadsPerFile=(availableLoaderThreads+files-1)/files;
//		
//		int zipTheadsPF=8;
//		int streamerThreadsPF=Tools.min(maxThreadsPerFile, Streamer.DEFAULT_THREADS);
//		int loaderThreadsPF=availableLoaderThreads;
//		final int loaderThreads=loaderThreadsPF*files;

		final int readThreadsPF=allocation[0];
		final int zipThreadsPF=allocation[1];
		final int streamerThreadsPF=allocation[2];
		final int covThreadsPF=allocation[3];
		
		System.err.println("Using "+zipThreadsPF+":"+streamerThreadsPF+":"+covThreadsPF+
			" zip:stream:cov threads for "+files+" files and "+Tools.plural("thread", Shared.threads())+".");
		
		ArrayList<Streamer> sslist=new ArrayList<Streamer>(files);
		AtomicLongArray[] covlist=new AtomicLongArray[files];
		for(int i=0; i<fnames.size(); i++){
			FileFormat ff=FileFormat.testInput(fnames.get(i), FileFormat.SAM, null, true, false);
			Streamer ss=StreamerFactory.makeSamOrBamStreamer(ff, streamerThreadsPF, false, false, -1, false);
			sslist.add(ss);
			ss.start();
			covlist[i]=new AtomicLongArray(contigs.size());
		}
		
		//Fill a list with LoadThreads
		ArrayList<LoadThread> alpt=new ArrayList<LoadThread>(covThreadsPF);
		for(int i=0; i<covThreadsPF; i++){
			final int sample=i%files;
//			System.err.println("Started a thread for sample "+sample+", pid="+i+
//				", threads="+covThreadsPF+", streamerThreads="+streamerThreadsPF);
			final LoadThread lt=new LoadThread(sslist.get(sample), 
				sample, contigMap, contigs, graph, covlist[sample], i);
			alpt.add(lt);
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState|=!success;
		
		for(int i=0; i<files; i++) {postprocess(covlist[i], contigs, i);}
		
		//Do anything necessary after processing
		
	}
	
	private void postprocess(AtomicLongArray depthArray, ArrayList<Contig> contigs, int sample) {
		for(int cnum=0; cnum<depthArray.length(); cnum++) {
			Contig c=contigs.get(cnum);
			float depth=depthArray.get(cnum)*1f/Tools.max(1, c.size());
			synchronized(c) {c.setDepth(depth, sample);}
		}
	}
	
	/** 
	 * Allocates threads optimally given ideal ratios and budget.
	 * Scales ratios proportionally and applies ceiling to ensure no component is starved.
	 * Each component is capped at its ideal maximum (ceiling of ideal value).
	 * @param ideal Array of ideal thread ratios for each component (also serves as max)
	 * @param budget Total threads available for allocation
	 * @return Array of allocated threads, one per component in ideal
	 */
	private static int[] allocateThreads(float[] ideal, int budget){
	    final int terms=ideal.length;
	    double total=Tools.sum(ideal);
	    double scale=budget/total;
	    
	    int[] allocated=new int[terms];
	    for(int i=0; i<terms; i++){
	        allocated[i]=Tools.min((int)Math.ceil(ideal[i]), (int)Math.ceil(ideal[i]*scale));
	    }
	    for(int i=0; i<ideal.length; i++) {
	   	 assert(Math.ceil(ideal[i])>=allocated[i]) : Arrays.toString(allocated)+", "+budget;
	    }
	    assert(Tools.min(allocated)>0 || Tools.min(allocated)<=0) : Arrays.toString(ideal)+", "+budget;
	    assert(Tools.sum(allocated)>=Tools.sum(ideal) || Tools.sum(ideal)>budget) : 
	        Arrays.toString(allocated)+", "+budget;
//	    System.err.println("Budget: "+Arrays.toString(ideal)+", "+budget+" -> "+Arrays.toString(allocated));
	    return allocated;
	}
	
	@Override
	public synchronized void accumulate(LoadThread t) {
		synchronized(t) {
			readsIn+=t.readsInT;
			readsUsed+=t.readsUsedT;
			basesIn+=t.basesInT;
			bytesIn+=t.bytesInT;
			errorState|=(t.success);
		}
	}

	@Override
	public ReadWriteLock rwlock() {return null;}

	@Override
	public synchronized boolean success() {
		return errorState;
	}
	
	/**
	 * Thread for processing a single SAM file.
	 * Calculates contig coverage depths and builds connectivity graph
	 * from paired-end reads mapping to different contigs.
	 */
	class LoadThread extends Thread {
		
		LoadThread(final Streamer ss_, final int sample_, HashMap<String, Contig> contigMap_, 
				ArrayList<Contig> contigs_, IntHashMap[] graph_, AtomicLongArray depth_, int tid_) {
			ss=ss_;
			sample=sample_;
			contigMap=contigMap_;
			contigs=contigs_;
			graph=graph_;
			depthArray=depth_;
			tid=tid_;
		}
		
		@Override
		public void run() {
			synchronized(this) {runInner();}
		}
		
		/**
		 * Main thread processing logic.
		 * Creates SAM streamer, processes alignments, and calculates final depths.
		 * Updates contig depth values after processing all reads.
		 */
		private void runInner() {
			if(tid<=sample) {outstream.println("Loading "+ss.fname());}
			//else {outstream.println("tid "+tid+">sample "+sample);}

//			System.err.println("SamLoader.LoadThread "+tid+" started processSam_Thread.");
			processSam_Thread(ss, depthArray);//They never leave here aside from one
			success=true;
//			System.err.println("SamLoader.LoadThread "+tid+" terminated successfully.");
		}
		
		void processSam_Thread(Streamer ss, AtomicLongArray depthArray) {
			ListNum<SamLine> ln=ss.nextLines();
			ArrayList<SamLine> reads=(ln==null ? null : ln.list);

			while(ln!=null && reads!=null && reads.size()>0){

				for(int idx=0; idx<reads.size(); idx++){
					SamLine sl=reads.get(idx);
					if(sl.mapped()) {
						boolean used=addSamLine(sl, depthArray);
						readsUsedT+=(used ? 1 : 0);
						readsInT++;
						basesInT+=(sl.seq==null ? 0 : sl.length());
						bytesInT+=(sl.countBytes());
					}
				}
				ln=ss.nextLines();
				reads=(ln==null ? null : ln.list);
			}
//			assert(false) : (ln==null ? "null" : "poison="+ln.poison()+", last="+ln.last());
		}
		
		/**
		 * Calculates aligned bases excluding contig tips for small contigs.
		 * For contigs shorter than 1.5*tipLimit, returns full alignment length.
		 * For longer contigs, excludes alignment portions near contig ends.
		 * @param sl SAM line with alignment information
		 * @param contigLen Length of the reference contig
		 * @return Number of aligned bases excluding tip regions
		 */
		private int calcAlignedBases(SamLine sl, int contigLen) {
			int aligned=sl.mappedNonClippedBases();
			if(contigLen<1.5f*tipLimit) {return aligned;}
			int limit=Tools.min(tipLimit, contigLen/4);
			final int lineStart=sl.start(false, false);
			final int lineStop=sl.stop(lineStart, false, false);
			final int contigStart=limit;
			final int contigStop=contigLen-limit;
			if(lineStart>=contigStart && lineStop<=contigStop) {return aligned;}
			return Tools.overlapLength(lineStart, lineStop, contigStart, contigStop);
		}
		
		private boolean addSamLine(SamLine sl, AtomicLongArray depthArray) {
			if(!sl.mapped()) {return false;}
			if(maxSubs<999 && sl.countSubs()>maxSubs) {return false;}
			if(minID>0 && sl.calcIdentity()<minID) {return false;}
			final String rname=ContigRenamer.toShortName(sl.rnameS());
			final Contig c1=contigMap.get(rname);
			if(c1==null) {return false;}//Contig not found; possibly too short
			assert(c1!=null) : "Can't find contig for rname "+rname;
			final int cid=c1.id();
			final int aligned=calcAlignedBases(sl, (int)c1.size());
			depthArray.addAndGet(cid, aligned);
			
			if(graph==null || sl.ambiguous() || !sl.hasMate() || !sl.nextMapped() 
					|| sl.pairedOnSameChrom() || sl.mapq<minMapq || aligned<minAlignedBases) {return true;}
			if(minMateq>0) {
				int mateq=sl.mateq();
				if(mateq>=0 && mateq<minMateq) {
//					System.err.println("mateq too low: "+mateq);
					return true;
				}
			}
			if(minMateID>0){
				float mateid=sl.mateID();
				if(mateid>0 && mateid<100*minMateID) {
//					System.err.println("mateid too low: "+mateid);
					return true;
				}
			}
			final String rnext=ContigRenamer.toShortName(sl.rnext());
			assert(rnext!=null && !"*".equals(rnext) && !"=".equals(rnext));
			
			final Contig c2=contigMap.get(rnext);
			if(c2==null) {
				//System.err.println("Can't find "+rnext);//Happens when using mincontig
				return true;
			}//Contig not found
//			System.err.println("Adding edge: "+rname+" - "+rnext);
			if(minEntropy>0 && sl.seq!=null && !et.passes(sl.seq, true)) {return true;}
//			if(minEntropy>0 && sl.seq!=null && EntropyTracker.calcEntropy(sl.seq, kmerCounts, 4)<minEntropy) {return true;}
			assert(c2!=null) : "Can't find contig for rnext "+rnext;
			
			//TODO:  Try commenting this out to see if it is the source of the nondeterminism.
			final IntHashMap destMap;
			synchronized(graph) {
				if(graph[cid]==null) {graph[cid]=new IntHashMap(5);}
				destMap=graph[cid];
			}
			synchronized(destMap) {
				destMap.increment(c2.id());
			}
			return true;
		}
		
		final Streamer ss;
		/** Sample number for depth assignment */
		final int sample;
		final int tid;
		/** Map from contig names to Contig objects */
		final HashMap<String, Contig> contigMap;
		/** Ordered list of contigs for depth calculation */
		final ArrayList<Contig> contigs;
		/** Array of connectivity maps between contigs */
		final IntHashMap[] graph;
		/** Entropy tracker for sequence complexity filtering */
		final EntropyTracker et=new EntropyTracker(5, 80, false, minEntropy, true);
		final AtomicLongArray depthArray;
		/** Thread-local count of input reads processed */
		long readsInT=0;
		/** Thread-local count of reads passing filters */
		long readsUsedT=0;
		/** Thread-local count of input bases processed */
		long basesInT=0;
		/** Thread-local count of input bytes processed */
		long bytesInT=0;
		/** Thread completion status flag */
		boolean success=false;
	}
	
	/** Output stream for progress messages */
	public PrintStream outstream=System.err;
	/** Total count of input reads across all threads */
	public long readsIn=0;
	/** Total count of reads passing all filters */
	public long readsUsed=0;
	/** Total count of input bases processed */
	public long basesIn=0;
	/** Total count of input bytes processed */
	public long bytesIn=0;
	/** Minimum mapping quality for connectivity graph edges */
	public int minMapq=4;
	/** Minimum mate mapping quality for graph edges */
	public int minMateq=4;
	/** Minimum alignment identity for read acceptance */
	public float minID=0f;
	/** Minimum mate alignment identity for graph edges */
	public float minMateID=0f;
	/** Maximum substitutions allowed for read acceptance */
	public int maxSubs=999;
	/** Distance from contig ends to exclude in depth calculations */
	public int tipLimit=100;
	/** Minimum sequence entropy for complexity filtering */
	public float minEntropy=0;
	/** Minimum aligned bases required for graph edge creation */
	public int minAlignedBases=0;
	
	/** Overall error state flag for processing */
	public boolean errorState=false;
	public static int MAX_SAM_LOADER_THREADS=1024;
	public static int MAX_SAM_LOADER_THREADS_PER_FILE=2;
	public static final boolean verbose=true;
	
}
