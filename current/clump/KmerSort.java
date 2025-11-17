package clump;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import bloom.KCountArray;
import fileIO.ReadWrite;
import jgi.BBMerge;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sort.ReadComparatorID;
import sort.ReadComparatorName;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.Read;
import structures.ListNum;
import structures.Quantizer;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date June 20, 2014
 *
 */
public abstract class KmerSort {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Count kmers */
	final void preprocess(){
		if(minCount>1){
			if(groups>1){
				table=ClumpTools.table();
				assert(table!=null);
			}else{
				Timer ctimer=new Timer();
				if(verbose){ctimer.start("Counting pivots.");}
				table=ClumpTools.getTable(in1, in2, k, minCount);
				if(verbose){ctimer.stop("Count time: ");}
			}
		}
	}

	/** Create read streams and process all data */
	abstract void process(Timer t);
	
	/**
	 * Prints comprehensive statistics about processing results.
	 * Includes reads processed, clumps formed, errors corrected, duplicates found,
	 * and output counts. Updates last processed counters and validates error state.
	 * @param t Timer containing elapsed processing time
	 */
	final void printStats(Timer t){
		table=null;
		ClumpTools.clearTable();
		
		errorState|=ReadStats.writeAll();
		
		t.stop();
		
		String rpstring2=readsProcessed+"";
		
		String cpstring=""+(groups==1 ? clumpsProcessedThisPass : clumpsProcessedTotal);
		String epstring=""+correctionsTotal;
		String efstring=""+(entryFiltered);
		String dpstring=""+(duplicatesTotal + entryFiltered);

		String rostring=""+readsOut;
		String bostring=""+basesOut;

		lastReadsIn=readsProcessed;
		lastBasesIn=basesProcessed;
		lastReadsOut=readsOut;
		lastBasesOut=basesOut;
		
		while(rpstring2.length()<12){rpstring2=" "+rpstring2;}
		while(cpstring.length()<12){cpstring=" "+cpstring;}
		while(epstring.length()<12){epstring=" "+epstring;}
		while(efstring.length()<12){efstring=" "+efstring;}
		while(dpstring.length()<12){dpstring=" "+dpstring;}
		dpstring+=String.format("\t%.3f%%", ((duplicatesTotal + entryFiltered)*100.0)/readsProcessed);

		while(rostring.length()<12){rostring=" "+rostring;}
		while(bostring.length()<12){bostring=" "+bostring;}
		
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 10));
		outstream.println();

		outstream.println("Reads In:         "+rpstring2);
		outstream.println("Clumps Formed:    "+cpstring);
		if(correct){
			outstream.println("Errors Corrected: "+epstring);
		}
		if(dedupe || entryfilter){
			outstream.println("Duplicates Found: "+dpstring);
			if(entryfilter && verbose && false){
				outstream.println(" -Entry Filtered: "+efstring);
			}
			outstream.println("Reads Out:        "+rostring);
			outstream.println("Bases Out:        "+bostring);
		}
		
		if(errorState){
			Clumpify.sharedErrorState=true;
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * Executes one complete pass of kmer-based read processing.
	 * Counts pivots if needed, hashes reads, sorts by kmer, creates clumps,
	 * then processes for correction or deduplication.
	 * @param reads Input reads to process
	 * @param kc KmerComparator for hashing and sorting operations
	 * @return Processed reads after correction or deduplication
	 */
	final ArrayList<Read> runOnePass(ArrayList<Read> reads, KmerComparator kc){
		Timer t=new Timer();
		
		table=null;
		if(minCount>1){
			if(verbose){t.start("Counting pivots.");}
			table=ClumpTools.getTable(reads, k, minCount);
			if(verbose){t.stop("Count time: ");}
		}
		
		if(verbose){t.start("Hashing.");}
		kc.hashThreaded(reads, table, minCount);
		if(verbose){t.stop("Hash time: ");}
		
		if(verbose){t.start("Sorting.");}
		Shared.sort(reads, kc);
		if(verbose){t.stop("Sort time: ");}
		
		if(verbose){t.start("Making clumps.");}
		readsProcessedThisPass=reads.size();
		ClumpList cl=new ClumpList(reads, k, false);
		reads.clear();
		clumpsProcessedThisPass=cl.size();
		clumpsProcessedTotal+=clumpsProcessedThisPass;
		if(verbose){t.stop("Clump time: ");}
		
		if(correct){
			if(verbose){t.start("Correcting.");}
			reads=processClumps(cl, ClumpList.CORRECT);
			if(verbose){t.stop("Correct time: ");}
		}else{
			assert(dedupe);
			if(verbose){t.start("Deduplicating.");}
			reads=processClumps(cl, ClumpList.DEDUPE);
			if(verbose){t.stop("Dedupe time: ");}
		}
		
		return reads;
	}
	
	/**
	 * Sorts reads by name and optionally pairs mates with matching names.
	 * For pairing, matches reads with identical names or valid FASTQ pair names,
	 * setting mate relationships and pair numbers.
	 * @param list Reads to sort
	 * @param pair Whether to pair mates after sorting
	 * @return Sorted and optionally paired reads
	 */
	static final ArrayList<Read> nameSort(ArrayList<Read> list, boolean pair){
		Shared.sort(list, ReadComparatorName.comparator);
		if(!pair){return list;}
		
		ArrayList<Read> list2=new ArrayList<Read>(1+list.size()/2);
		
		Read prev=null;
		for(Read r : list){
			if(prev==null){
				prev=r;
				assert(prev.mate==null);
			}else{
				if(prev.id.equals(r.id) || FASTQ.testPairNames(prev.id, r.id, true)){
					prev.mate=r;
					r.mate=prev;
					prev.setPairnum(0);
					r.setPairnum(1);
					list2.add(prev);
					prev=null;
				}else{
					list2.add(prev);
					prev=r;
				}
			}
		}
		return list2;
	}
	
	/**
	 * Sorts reads by numeric ID and optionally pairs mates with matching IDs.
	 * For pairing, matches reads with identical numeric IDs, assuming pairnum 0 and 1.
	 * @param list Reads to sort
	 * @param pair Whether to pair mates after sorting
	 * @return Sorted and optionally paired reads
	 */
	static final ArrayList<Read> idSort(ArrayList<Read> list, boolean pair){
		Shared.sort(list, ReadComparatorID.comparator);
		if(!pair){return list;}
		
		ArrayList<Read> list2=new ArrayList<Read>(1+list.size()/2);
		
		Read prev=null;
		for(Read r : list){
			if(prev==null){
				prev=r;
				assert(prev.mate==null);
			}else{
				if(prev.numericID==r.numericID){
					assert(prev.pairnum()==0 && r.pairnum()==1) : prev.id+"\n"+r.id;
					prev.mate=r;
					r.mate=prev;
					prev.setPairnum(0);
					r.setPairnum(1);
					list2.add(prev);
					prev=null;
				}else{
					list2.add(prev);
					prev=r;
				}
			}
		}
		return list2;
	}
	
	/**
	 * Filters paired reads to retain only read1 (pairnum 0) from each pair.
	 * @param list Paired reads to filter
	 * @return List containing only first reads from each pair
	 */
	static final ArrayList<Read> read1Only(ArrayList<Read> list){
		ArrayList<Read> list2=new ArrayList<Read>(1+list.size()/2);
		for(Read r : list){
			assert(r.mate!=null) : r+"\n"+r.mate;
			if(r.pairnum()==0){list2.add(r);}
		}
		return list2;
	}
	
//	@Deprecated
//	//No longer needed
//	public int countClumps(ArrayList<Read> list){
//		int count=0;
//		long currentKmer=-1;
//		for(final Read r : list){
//			final ReadKey key=(ReadKey)r.obj;
//			if(key.kmer!=currentKmer){
//				currentKmer=key.kmer;
//				count++;
//			}
//		}
//		return count;
//	}
	
	/**
	 * Processes clumps for error correction or duplicate removal.
	 * Updates correction and duplicate statistics from processing results.
	 * @param cl ClumpList containing reads grouped by kmers
	 * @param mode Processing mode (CORRECT or DEDUPE)
	 * @return Processed reads after clump operations
	 */
	public final ArrayList<Read> processClumps(ClumpList cl, int mode){
		long[] rvector=KillSwitch.allocLong1D(2);
		ArrayList<Read> out=cl.process(Shared.threads(), mode, rvector);
		correctionsThisPass=rvector[0];
		correctionsTotal+=correctionsThisPass;
		duplicatesThisPass=rvector[1];
		duplicatesTotal+=duplicatesThisPass;
		cl.clear();
		return out;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Hashes reads and splits them into groups using multi-threaded processing.
	 * Each thread processes a subset of reads, hashing and distributing to groups
	 * based on kmer hash values.
	 * @param list Input reads to hash and split
	 * @param kc KmerComparator for hashing operations
	 * @param array Output arrays for each group
	 */
	public final void hashAndSplit(ArrayList<Read> list, KmerComparator kc, ArrayList<Read>[] array){
		int threads=Shared.threads();
		ArrayList<HashSplitThread> alt=new ArrayList<HashSplitThread>(threads);
		for(int i=0; i<threads; i++){alt.add(new HashSplitThread(i, threads, list, kc));}
		for(HashSplitThread ht : alt){ht.start();}
		
		/* Wait for threads to die */
		for(HashSplitThread ht : alt){
			
			/* Wait for a thread to die */
			while(ht.getState()!=Thread.State.TERMINATED){
				try {
					ht.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			for(int i=0; i<groups; i++){
				array[i].addAll(ht.array[i]);
				ht.array[i]=null;
			}
		}
	}
	
	ArrayList<Read> fetchReads1(final ConcurrentReadInputStream cris, final KmerComparator kc){
		Timer t=new Timer();
		if(verbose){t.start("Making fetch threads.");}
		final int threads=Shared.threads();
		ArrayList<FetchThread1> alft=new ArrayList<FetchThread1>(threads);
		for(int i=0; i<threads; i++){alft.add(new FetchThread1(i, cris, kc, unpair));}
		
		readsThisPass=memThisPass=entryFilteredThisPass=0;
		
		if(verbose){outstream.println("Starting threads.");}
		for(FetchThread1 ht : alft){ht.start();}
		
		
		if(verbose){outstream.println("Waiting for threads.");}
		/* Wait for threads to die */
		for(FetchThread1 ht : alft){
			
			/* Wait for a thread to die */
			while(ht.getState()!=Thread.State.TERMINATED){
				try {
					ht.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			entryFilteredThisPass+=ht.entryFilteredT;
			readsThisPass+=ht.readsProcessedT;
			basesProcessed+=ht.basesProcessedT;
			diskProcessed+=ht.diskProcessedT;
			memThisPass+=ht.memProcessedT;
		}
		readsProcessed+=readsThisPass;
		memProcessed+=memThisPass;
		entryFiltered+=entryFilteredThisPass;

		if(verbose){t.stop("Fetch time: ");}
		if(verbose){System.err.println("Closing input stream.");}
		errorState=ReadWrite.closeStream(cris)|errorState;
		
		if(verbose){t.start("Combining thread output.");}
		long readsLeft=readsThisPass-entryFilteredThisPass;
		long slotsLeft=cris.paired() && !unpair ? readsLeft/2 : readsLeft;
		assert(slotsLeft<=Shared.MAX_ARRAY_LEN) :
			"\nThe number of reads is greater than 2 billion, which is the limit for a single group. "
			+ "\nPlease rerun and manually specify 'groups=7' or similar, "
			+ "\nsuch that the number of reads per group is less than 2 billion.";
		ArrayList<Read> list=new ArrayList<Read>((int)(slotsLeft));
		for(int i=0; i<threads; i++){
			FetchThread1 ft=alft.set(i, null);
			list.addAll(ft.storage);
		}
		if(verbose){t.stop("Combine time: ");}
		
		assert(list.size()==slotsLeft) : list.size()+", "+readsThisPass+", "+readsLeft+", "+slotsLeft+", "+cris.paired();
		//assert(list.size()==readsLeft || (cris.paired() && list.size()*2==readsLeft)) : list.size()+", "+readsThisPass+", "+readsLeft+", "+cris.paired();
		ecco=false;
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class FetchThread1 extends Thread{
		
		/**
		 * Constructs thread for fetching reads from concurrent input stream.
		 * @param id_ Thread identifier
		 * @param cris_ Input stream to read from
		 * @param kc_ KmerComparator for hashing operations
		 * @param unpair_ Whether to unpair reads during processing
		 */
		FetchThread1(int id_, ConcurrentReadInputStream cris_, KmerComparator kc_, boolean unpair_){
			id=id_;
			cris=cris_;
			kc=kc_;
			storage=new ArrayList<Read>();
			unpairT=unpair_;
			entryFilterTable=(entryfilter ? new HashMap<Long, Read>() : null);
		}
		
		@Override
		public void run(){
			ListNum<Read> ln=cris.nextList();
			final boolean paired=cris.paired();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				
				for(Read r : reads){
					if(!r.validated()){
						r.validate(true);
						if(r.mate!=null){r.mate.validate(true);}
					}
					readsProcessedT+=1+r.mateCount();
					basesProcessedT+=r.length()+r.mateLength();
//					diskProcessedT+=r.countFastqBytes()+r.countMateFastqBytes();
//					memProcessedT+=r.countBytes()+r.countMateBytes();
					if(shrinkName){
						Clumpify.shrinkName(r);
						Clumpify.shrinkName(r.mate);
					}else if(shortName){
						Clumpify.shortName(r);
						Clumpify.shortName(r.mate);
					}
					
					if(quantizeQuality){
						Quantizer.quantize(r, r.mate);
					}
				}
				
				if(ecco){
					for(Read r : reads){
						Read r2=r.mate;
						assert(r.obj==null) : "TODO: Pivot should not have been generated yet, though it may be OK.";
						assert(r2!=null) : "ecco requires paired reads.";
						if(r2!=null){
							int x=BBMerge.findOverlapStrict(r, r2, true);
							if(x>=0){
								r.obj=null;
								r2.obj=null;
							}
						}
					}
				}
				
				if(entryFilterTable!=null){
					int removed=0;
					for(int i=0; i<reads.size(); i++){
						Read r=reads.get(i);
						final long key=Hasher.hashPair(r);
						final Long key2=Long.valueOf(key);
						final Read old=entryFilterTable.get(key2);
						if(old==null){
							entryFilterTable.put(key2, r);
						}else{
							boolean same=Hasher.equalsPaired(r, old);
							if(same){
								removed++;
								entryFilteredT+=r.pairCount();
								reads.set(i, null);
							}
						}
					}
					if(removed>0){Tools.condenseStrict(reads);}
				}
				
				ArrayList<Read> hashList=reads;
				if(paired && unpairT){
					hashList=new ArrayList<Read>(reads.size()*2);
					for(Read r1 : reads){
						Read r2=r1.mate;
						assert(r2!=null);
						hashList.add(r1);
						hashList.add(r2);
						if(groups>1 || !repair || namesort){
							r1.mate=null;
							r2.mate=null;
						}
					}
				}
				
				kc.hash(hashList, table, minCount, true);
				storage.addAll(hashList);
				cris.returnList(ln.id, false);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
			
			//Optimization for TimSort
			if(parallelSort){
				storage.sort(kc);
//				Shared.sort(storage, kc); //Already threaded; this is not needed.
			}else{
				Collections.sort(storage, kc);
			}
		}

		/** Thread identifier for the fetch thread. */
		final int id;
		/** Concurrent input stream for reading sequences. */
		final ConcurrentReadInputStream cris;
		/** KmerComparator for hashing and comparing reads. */
		final KmerComparator kc;
		/** Storage for processed reads from this thread. */
		final ArrayList<Read> storage;
		/** Whether this thread should unpair reads during processing. */
		final boolean unpairT;
		/** Hash table for entry filtering to detect duplicate reads. */
		final HashMap<Long, Read> entryFilterTable;
		/** Number of entries filtered by this thread. */
		public long entryFilteredT=0;
		
		/** Number of reads processed by this thread. */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread. */
		protected long basesProcessedT=0;
		/** Disk bytes processed by this thread. */
		protected long diskProcessedT=0;
		/** Memory bytes processed by this thread. */
		protected long memProcessedT=0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private final class HashSplitThread extends Thread{
		
		/**
		 * Constructs thread for hashing reads and splitting into groups.
		 * @param id_ Thread identifier
		 * @param threads_ Total number of threads
		 * @param list_ Input reads to process
		 * @param kc_ KmerComparator for hashing operations
		 */
		@SuppressWarnings("unchecked")
		HashSplitThread(int id_, int threads_, ArrayList<Read> list_, KmerComparator kc_){
			id=id_;
			threads=threads_;
			list=list_;
			kc=kc_;
			array=new ArrayList[groups];
			for(int i=0; i<groups; i++){
				array[i]=new ArrayList<Read>();
			}
		}
		
		@Override
		public void run(){
			for(int i=id; i<list.size(); i+=threads){
				Read r=list.get(i);
				kc.hash(r, null, 0, true);
				ReadKey key=(ReadKey)r.obj;
				array[(int)(kc.hash(key.kmer)%groups)].add(r);
			}
		}
		
		/** Thread identifier for the hash-split thread. */
		final int id;
		/** Total number of threads in the thread pool. */
		final int threads;
		/** Input list of reads to be hashed and split. */
		final ArrayList<Read> list;
		/** KmerComparator for hashing operations. */
		final KmerComparator kc;
		/** Output arrays for distributing reads across groups. */
		final ArrayList<Read>[] array;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Kmer length for sequence hashing and comparison operations. */
	int k=31;
	/** Minimum occurrence count for kmers to be considered valid pivots. */
	int minCount=0;
	
	/** Number of groups to split processing into for memory management. */
	int groups=1;
	
	/** Kmer count table for tracking occurrence frequencies during processing. */
	KCountArray table=null;
	
	/*--------------------------------------------------------------*/
	/*----------------          I/O Fields          ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path. */
	String in1=null;
	/** Secondary input file path for paired reads. */
	String in2=null;

	/** Primary output file path. */
	String out1=null;
	/** Secondary output file path for paired reads. */
	String out2=null;
	
	/** Input file extension override. */
	String extin=null;
	/** Output file extension override. */
	String extout=null;
	
	/*--------------------------------------------------------------*/
	
	/** Total number of reads processed across all operations. */
	protected long readsProcessed=0;
	/** Total number of bases processed across all operations. */
	protected long basesProcessed=0;
	/** Total bytes read from disk during processing. */
	protected long diskProcessed=0;
	/** Total bytes processed in memory during operations. */
	protected long memProcessed=0;
	/** Total reads filtered out due to entry filtering (duplicate detection). */
	protected static long entryFiltered=0;

	/** Total number of reads written to output files. */
	protected long readsOut=0;
	/** Total number of bases written to output files. */
	protected long basesOut=0;

	/** Reads filtered during current processing pass. */
	protected long entryFilteredThisPass=0;
	/** Reads processed during current pass. */
	protected long readsThisPass=0;
	/** Memory bytes processed during current pass. */
	protected long memThisPass=0;
	
	/** Number of reads processed in the current sorting pass. */
	protected long readsProcessedThisPass=0;
	/** Number of clumps processed during current pass. */
	protected long clumpsProcessedThisPass=0;
	/** Number of error corrections made during current pass. */
	protected long correctionsThisPass=0;
	
	/** Number of duplicates found during current pass. */
	protected long duplicatesThisPass=0;
	/** Total number of duplicates found across all passes. */
	protected static long duplicatesTotal=0;
	
	/** Total number of clumps processed across all passes. */
	protected long clumpsProcessedTotal=0;
	/** Total number of error corrections made across all passes. */
	protected static long correctionsTotal=0;
	
	/** Number of processing passes to perform. */
	protected int passes=1;
	
	/** Maximum number of reads to process, -1 for unlimited. */
	long maxReads=-1;
	/** Whether to add information to read names during processing. */
	protected boolean addName=false;
	/** Whether to use shortened read names. */
	boolean shortName=false;
	/** Whether to shrink read names to save memory. */
	boolean shrinkName=false;
	/** Whether to consider reverse complement kmers during processing. */
	boolean rcomp=false;
	/** Whether to condense output by removing null entries. */
	boolean condense=false;
	/** Whether to perform error correction on reads. */
	boolean correct=false;
	/** Whether to perform duplicate removal. */
	boolean dedupe=false;
	/** Whether to split input processing across multiple groups. */
	boolean splitInput=false;
	/** Whether to use ECCO (Error Correction and Consensus from Overlaps) mode. */
	boolean ecco=false;
	/** Whether to unpair reads during processing. */
	boolean unpair=false;
	/** Whether to repair read pairing after processing. */
	boolean repair=false;
	/** Whether to sort reads by name during processing. */
	boolean namesort=false;
	/** Whether to filter duplicate entries during input processing. */
	boolean entryfilter=false;
	/** Whether to use parallel sorting algorithms. */
	final boolean parallelSort=Shared.parallelSort;
	/** Whether memory usage warnings have been issued. */
	boolean memWarned=false;
	
	/** Whether to use shared headers for output files. */
	boolean useSharedHeader=false;
	/** Mode for reordering reads after processing. */
	int reorderMode=REORDER_FALSE;
	
	/*--------------------------------------------------------------*/

	/** Number of reads from the last processing run. */
	public static long lastReadsIn=-1;
	/** Number of bases from the last processing run. */
	public static long lastBasesIn=-1;
	/** Number of reads output from the last processing run. */
	public static long lastReadsOut=-1;
	/** Number of bases output from the last processing run. */
	public static long lastBasesOut=-1;
	
	/** Whether to quantize quality scores to reduce memory usage. */
	static boolean quantizeQuality=false;
	static final int REORDER_FALSE=0, REORDER_CONSENSUS=1, REORDER_PAIRED=2, REORDER_AUTO=3;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and statistics. */
	PrintStream outstream=System.err;
	/** Whether to print verbose status and timing information. */
	public static boolean verbose=true;
	/** Whether to use hash-and-split processing for large datasets. */
	public static boolean doHashAndSplit=true;
	/** Whether an error state has been encountered during processing. */
	public boolean errorState=false;
	/** Whether to overwrite existing output files. */
	boolean overwrite=true;
	/** Whether to append to existing output files instead of overwriting. */
	boolean append=false;
	
}
