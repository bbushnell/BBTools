package bin;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLongArray;
import java.util.concurrent.locks.ReadWriteLock;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.IntHashMap2;
import shared.Shared;
import shared.Tools;
import stream.SamLine;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.EntropyTracker;

/**
 * Loads SAM/BAM files to compute contig coverage and connectivity graphs.
 * Features:
 * 1. Queue-based loading to limit open file handles (MAX_CONCURRENT_FILES).
 * 2. Modulo-based sample merging to limit memory usage (MAX_SAMPLES).
 * e.g., if MAX_SAMPLES=8 and inputs=96, file 0 and file 8 map to sample 0.
 * 
 * @author Brian Bushnell
 * @contributor Collei
 * @date December 17, 2025
 */
public class SamLoader3 implements Accumulator<SamLoader3.Worker> {
	
	/** Builds a loader that reports status to the provided PrintStream. */
	public SamLoader3(PrintStream outstream_) {
		outstream=outstream_;
	}
	
	/**
	 * Deprecated entry point that converts the contig map to a sorted list before loading.
	 */
	@Deprecated
	public void load(ArrayList<String> fnames, HashMap<String, Contig> contigMap, IntHashMap2[] graph) {
		ArrayList<Contig> list=new ArrayList<Contig>(contigMap.values());
		Collections.sort(list);
		for(int i=0; i<list.size(); i++) {list.get(i).setID(i);}
		load(fnames, contigMap, list, graph);
	}
	
	/**
	 * Main loader entry point.
	 * Branches between parallel loading (few files) and queue-based loading (many files).
	 */
	public void load(ArrayList<String> fnames, HashMap<String, Contig> contigMap, 
			ArrayList<Contig> contigs, IntHashMap2[] graph){
		
		SamLine.RNAME_AS_BYTES=false;
		
		//Determine effective number of samples to store in memory
		final int files=fnames.size();
		final int effectiveSamples=Tools.min(files, MAX_SAMPLES);
		
		if(files > effectiveSamples) {
			outstream.println("Input contains "+files+" files. Merging into "+
				effectiveSamples+" samples (Modulo Mode) to conserve memory.");
		}
		
		if(files > MAX_CONCURRENT_FILES) {
			loadQueue(fnames, contigMap, contigs, graph, effectiveSamples);
		} else {
			loadParallel(fnames, contigMap, contigs, graph, effectiveSamples);
		}
	}

	/**
	 * Queue-based strategy: Limits concurrent open files to MAX_CONCURRENT_FILES.
	 * Threads pick the next file from the list, open it, process it, and close it.
	 */
	private void loadQueue(ArrayList<String> fnames, HashMap<String, Contig> contigMap, 
			ArrayList<Contig> contigs, IntHashMap2[] graph, int effectiveSamples) {
		
		final int files=fnames.size();
		final int concurrentFiles=Tools.min(files, MAX_CONCURRENT_FILES);
		
		//Initialize coverage arrays for effective samples only
		AtomicLongArray[] covlist=new AtomicLongArray[effectiveSamples];
		for(int i=0; i<effectiveSamples; i++) {
			covlist[i]=new AtomicLongArray(contigs.size());
		}
		
		//Determine thread allocation per ACTIVE file
		final int availableThreads=Tools.max(1, Shared.threads());
		final int[] threadAllocation = calculateThreadAllocation(fnames.get(0), concurrentFiles, availableThreads);
		final int streamerThreadsPF=threadAllocation[2]; //Index 2 is streamer threads
		
		outstream.println("Processing "+files+" files in queue mode (concurrency="+concurrentFiles+").");
		
		//Shared atomic index for threads to grab the next file
		final AtomicInteger nextFileIndex = new AtomicInteger(0);
		
		ArrayList<Worker> alpt=new ArrayList<Worker>(concurrentFiles);
		for(int i=0; i<concurrentFiles; i++){
			final QueueWorker lt=new QueueWorker(nextFileIndex, fnames, streamerThreadsPF,
				contigMap, contigs, graph, covlist, i);
			alpt.add(lt);
		}
		
		//Start and wait
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState|=!success;
		
		//Post-process the consolidated sample arrays
		for(int i=0; i<effectiveSamples; i++) {postprocess(covlist[i], contigs, i);}
	}

	/**
	 * Original strategy: Opens ALL files at once. Best for small numbers of files.
	 */
	private void loadParallel(ArrayList<String> fnames, HashMap<String, Contig> contigMap, 
			ArrayList<Contig> contigs, IntHashMap2[] graph, int effectiveSamples) {
		
		final int files=fnames.size();
		
		//Determine thread allocation
		final int availableThreads=Tools.max(1, Shared.threads());
		final int[] threadAllocation = calculateThreadAllocation(fnames.get(0), files, availableThreads);
		
		final int streamerThreadsPF=threadAllocation[2];
		final int covThreadsPF=threadAllocation[3];
		final int zipThreadsPF=threadAllocation[1];
		
		System.err.println("Using "+zipThreadsPF+":"+streamerThreadsPF+":"+covThreadsPF+
			" zip:stream:cov threads for "+files+" files and "+Tools.plural("thread", Shared.threads())+".");
		
		//Pre-open all Streamers
		ArrayList<Streamer> sslist=new ArrayList<Streamer>(files);
		
		//Allocated condensed coverage arrays
		AtomicLongArray[] covlist=new AtomicLongArray[effectiveSamples];
		for(int i=0; i<effectiveSamples; i++) {
			covlist[i]=new AtomicLongArray(contigs.size());
		}
		
		for(int i=0; i<fnames.size(); i++){
			FileFormat ff=FileFormat.testInput(fnames.get(i), FileFormat.SAM, null, true, false);
			Streamer ss=StreamerFactory.makeSamOrBamStreamer(ff, streamerThreadsPF, false, false, -1, false);
			sslist.add(ss);
			ss.start();
		}
		
		//Fill a list with ParallelWorkers
		final int threads=covThreadsPF*files;
		assert(threads>=files);
		ArrayList<Worker> alpt=new ArrayList<Worker>(threads);
		for(int i=0; i<threads; i++){
			final int fileIndex=i%files;
			final int sampleIndex=fileIndex%effectiveSamples; //Map file to sample bucket
			final ParallelWorker lt=new ParallelWorker(sslist.get(fileIndex), 
				sampleIndex, contigMap, contigs, graph, covlist[sampleIndex], i);
			alpt.add(lt);
		}
		
		//Start and wait
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState|=!success;
		
		//Close streams
		for(Streamer st : sslist) {
			ReadWrite.closeStreams(st, (Writer[])null);
		}
		
		for(int i=0; i<effectiveSamples; i++) {postprocess(covlist[i], contigs, i);}
	}
	
	/**
	 * Calculates thread distribution (bgzip vs streamer vs processing) based on available threads and active files.
	 */
	private int[] calculateThreadAllocation(String fname, int activeFiles, int availableThreads) {
		FileFormat ff0=FileFormat.testInput(fname, FileFormat.SAM, null, false, false);
		
		final int maxThreadsPerFile=(availableThreads+activeFiles-1)/activeFiles;
		final float[] ideal;
		{//FileRead, Decompress, SamLine, Coverage
			float x=(activeFiles>MAX_CONCURRENT_FILES ? 1 : MAX_SAM_LOADER_THREADS_PER_FILE);
			if(ff0.bam()) {//bam
				ideal=new float[] {1f, 4f, 6f, x};
			}else if(ff0.bgzip() || ff0.bz2()) {//sam.bgz or bz2
				ideal=new float[] {1f, 4f, 6f, x};
			}else if(ff0.gzip()) {//sam.gz, non-bgzip. Detectable from magic word...
				ideal=new float[] {0.5f, 1f, 2f, 0.5f*x};
			}else {//sam
				assert(ff0.sam());
				ideal=new float[] {1f, 0f, 6f, x};
			}
		}
		return allocateThreads(ideal, availableThreads);
	}

	/**
	 * Converts per-contig depth totals into normalized depth values and stores them on the contigs.
	 */
	private void postprocess(AtomicLongArray depthArray, ArrayList<Contig> contigs, int sample) {
		for(int cnum=0; cnum<depthArray.length(); cnum++) {
			Contig c=contigs.get(cnum);
			float depth=depthArray.get(cnum)*1f/Tools.max(1, c.size());
			synchronized(c) {c.setDepth(depth, sample);}
		}
	}
	
	private static int[] allocateThreads(float[] ideal, int budget){
	    final int terms=ideal.length;
	    double total=Tools.sum(ideal);
	    double scale=budget/total;
	    
	    int[] allocated=new int[terms];
	    for(int i=0; i<terms; i++){
	        allocated[i]=Tools.min((int)Math.ceil(ideal[i]), (int)Math.ceil(ideal[i]*scale));
	    }
	    assert(Tools.min(allocated)>0 || Tools.min(allocated)<=0);
	    assert(Tools.sum(allocated)>=Tools.sum(ideal) || Tools.sum(ideal)>budget);
	    return allocated;
	}
	
	@Override
	public synchronized void accumulate(Worker t) {
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
	public synchronized boolean success() {return errorState;}
	
	/**
	 * Base class for load workers. 
	 */
	abstract class Worker extends Thread {
		
		Worker(HashMap<String, Contig> contigMap_, ArrayList<Contig> contigs_, IntHashMap2[] graph_) {
			contigMap=contigMap_;
			contigs=contigs_;
			graph=graph_;
			et=new EntropyTracker(5, 80, false, minEntropy, true);
		}

		void processSam(Streamer ss, AtomicLongArray depthArray) {
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
		}
		
//		private int calcAlignedBases(SamLine sl, int contigLen) {
//			int aligned=sl.mappedNonClippedBases();
//			if(contigLen<1.5f*tipLimit) {return aligned;}
//			int limit=Tools.min(tipLimit, contigLen/4);
//			final int lineStart=sl.start(false, false);
//			final int lineStop=sl.stop(lineStart, false, false);
//			final int contigStart=limit;
//			final int contigStop=contigLen-limit;
//			if(lineStart>=contigStart && lineStop<=contigStop) {return aligned;}
//			return Tools.overlapLength(lineStart, lineStop, contigStart, contigStop);
//		}
		
		private int calcAlignedBases(SamLine sl, int contigLen) {
			int aligned=sl.mappedNonClippedBases();
			if(contigLen<1.5f*tipLimit) {return aligned;}
			
			// Calculate the trim limit just like before
			int limit=Tools.min(tipLimit, contigLen/4);
			
			// Define the "valid" center region
			final int contigStart=limit;
			final int contigStop=contigLen-limit;
			
			// Calculate bases aligned to the center region
			final int lineStart=sl.start(false, false);
			final int lineStop=sl.stop(lineStart, false, false);
			
			// Scale: (Mapped Bases in Center) * (Total Length / Center Length)
			final float scale=contigLen/(float)(contigStop-contigStart+1);
			final int overlap=Tools.overlapLength(lineStart, lineStop, contigStart, contigStop);
			return Math.round(overlap*scale);
		}
		
		private boolean addSamLine(SamLine sl, AtomicLongArray depthArray) {
			if(!sl.mapped()) {return false;}
			if(maxSubs<999 && sl.countSubs()>maxSubs) {return false;}
			if(minID>0 && sl.calcIdentity()<minID) {return false;}
			final String rname=ContigRenamer.toShortName(sl.rnameS());
			
			final Contig c1=contigMap.get(rname);
			if(c1==null) {return false;}
			assert(c1!=null) : "Can't find contig for rname "+rname;
			final int cid=c1.id();
			final int aligned=calcAlignedBases(sl, (int)c1.size());
			
			//Atomic accumulation handles multiple threads hitting the same sample bin
			depthArray.addAndGet(cid, aligned);
			
			if(graph==null || sl.ambiguous() || !sl.hasMate() || !sl.nextMapped() 
					|| sl.pairedOnSameChrom() || sl.mapq<minMapq || aligned<minAlignedBases) {return true;}
			if(minMateq>0) {
				int mateq=sl.mateq();
				if(mateq>=0 && mateq<minMateq) {return true;}
			}
			if(minMateID>0){
				float mateid=sl.mateID();
				if(mateid>0 && mateid<100*minMateID) {return true;}
			}
			final String rnext=ContigRenamer.toShortName(sl.rnext());
			assert(rnext!=null && !"*".equals(rnext) && !"=".equals(rnext));
			
			final Contig c2=contigMap.get(rnext);
			if(c2==null) {return true;}
			
			if(minEntropy>0 && sl.seq!=null && !et.passes(sl.seq, true)) {return true;}
			assert(c2!=null) : "Can't find contig for rnext "+rnext;
			
			final IntHashMap2 destMap;
			synchronized(graph) {
				if(graph[cid]==null) {graph[cid]=new IntHashMap2(5);}
				destMap=graph[cid];
			}
			synchronized(destMap) {
				destMap.increment(c2.id());
			}
			return true;
		}
		
		final HashMap<String, Contig> contigMap;
		final ArrayList<Contig> contigs;
		final IntHashMap2[] graph;
		final EntropyTracker et;
		
		long readsInT=0;
		long readsUsedT=0;
		long basesInT=0;
		long bytesInT=0;
		boolean success=false;
	}
	
	/**
	 * New Style: Worker that grabs files from a queue.
	 */
	class QueueWorker extends Worker {
		
		QueueWorker(AtomicInteger nextFileIndex_, ArrayList<String> fnames_, int streamerThreadsPF_,
				HashMap<String, Contig> contigMap_, ArrayList<Contig> contigs_, 
				IntHashMap2[] graph_, AtomicLongArray[] covlist_, int tid_) {
			super(contigMap_, contigs_, graph_);
			nextFileIndex=nextFileIndex_;
			fnames=fnames_;
			streamerThreadsPF=streamerThreadsPF_;
			covlist=covlist_;
			tid=tid_;
		}

		@Override
		public void run() {
			synchronized(this) {runInner();}
		}
		
		private void runInner() {
			int idx=nextFileIndex.getAndIncrement();
			while(idx < fnames.size()) {
				String fname=fnames.get(idx);
				outstream.println("Thread "+tid+" loading "+fname);
				
				//Create and start streamer
				FileFormat ff=FileFormat.testInput(fname, FileFormat.SAM, null, true, false);
				Streamer ss=StreamerFactory.makeSamOrBamStreamer(ff, streamerThreadsPF, false, false, -1, false);
				ss.start();
				
				//Process - Map file index to sample bucket
				int sampleBucket = idx % covlist.length;
				processSam(ss, covlist[sampleBucket]);
				
				//Close immediately
				ReadWrite.closeStreams(ss, (Writer[])null);
				
				//Next
				idx=nextFileIndex.getAndIncrement();
			}
			success=true;
		}
		
		final AtomicInteger nextFileIndex;
		final ArrayList<String> fnames;
		final int streamerThreadsPF;
		final AtomicLongArray[] covlist;
		final int tid;
	}
	
	/**
	 * Old Style: Worker dedicated to a single pre-opened Streamer.
	 */
	class ParallelWorker extends Worker {
		
		ParallelWorker(final Streamer ss_, final int sample_, HashMap<String, Contig> contigMap_, 
				ArrayList<Contig> contigs_, IntHashMap2[] graph_, AtomicLongArray depth_, int tid_) {
			super(contigMap_, contigs_, graph_);
			ss=ss_;
			sample=sample_;
			depthArray=depth_;
			tid=tid_;
		}
		
		@Override
		public void run() {
			synchronized(this) {runInner();}
		}
		
		private void runInner() {
			if(tid<=sample) {outstream.println("Loading "+ss.fname());}
			processSam(ss, depthArray);
			success=true;
		}
		
		final Streamer ss;
		final int sample;
		final int tid;
		final AtomicLongArray depthArray;
	}
	
	public PrintStream outstream=System.err;
	public long readsIn=0;
	public long readsUsed=0;
	public long basesIn=0;
	public long bytesIn=0;
	public int minMapq=4;
	public int minMateq=4;
	public float minID=0f;
	public float minMateID=0f;
	public int maxSubs=999;
	public int tipLimit=100;
	public float minEntropy=0;
	public int minAlignedBases=0;
	
	public boolean errorState=false;
	public static int MAX_SAM_LOADER_THREADS=1024;
	public static int MAX_SAM_LOADER_THREADS_PER_FILE=2;
	
	/** Max files to open concurrently. */
	public static int MAX_CONCURRENT_FILES=4;
	
	/** Max samples to store in memory. Files beyond this limit are merged modulo N. */
	public static int MAX_SAMPLES=8;
	
	public static final boolean verbose=true;
}