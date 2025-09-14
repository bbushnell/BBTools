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

		final int id;
		final ConcurrentReadInputStream cris;
		final KmerComparator kc;
		final ArrayList<Read> storage;
		final boolean unpairT;
		final HashMap<Long, Read> entryFilterTable;
		public long entryFilteredT=0;
		
		protected long readsProcessedT=0;
		protected long basesProcessedT=0;
		protected long diskProcessedT=0;
		protected long memProcessedT=0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private final class HashSplitThread extends Thread{
		
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
		
		final int id;
		final int threads;
		final ArrayList<Read> list;
		final KmerComparator kc;
		final ArrayList<Read>[] array;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	int k=31;
	int minCount=0;
	
	int groups=1;
	
	KCountArray table=null;
	
	/*--------------------------------------------------------------*/
	/*----------------          I/O Fields          ----------------*/
	/*--------------------------------------------------------------*/

	String in1=null;
	String in2=null;

	String out1=null;
	String out2=null;
	
	String extin=null;
	String extout=null;
	
	/*--------------------------------------------------------------*/
	
	protected long readsProcessed=0;
	protected long basesProcessed=0;
	protected long diskProcessed=0;
	protected long memProcessed=0;
	protected static long entryFiltered=0;

	protected long readsOut=0;
	protected long basesOut=0;

	protected long entryFilteredThisPass=0;
	protected long readsThisPass=0;
	protected long memThisPass=0;
	
	protected long readsProcessedThisPass=0;
	protected long clumpsProcessedThisPass=0;
	protected long correctionsThisPass=0;
	
	protected long duplicatesThisPass=0;
	protected static long duplicatesTotal=0;
	
	protected long clumpsProcessedTotal=0;
	protected static long correctionsTotal=0;
	
	protected int passes=1;
	
	long maxReads=-1;
	protected boolean addName=false;
	boolean shortName=false;
	boolean shrinkName=false;
	boolean rcomp=false;
	boolean condense=false;
	boolean correct=false;
	boolean dedupe=false;
	boolean splitInput=false;
	boolean ecco=false;
	boolean unpair=false;
	boolean repair=false;
	boolean namesort=false;
	boolean entryfilter=false;
	final boolean parallelSort=Shared.parallelSort;
	boolean memWarned=false;
	
	boolean useSharedHeader=false;
	int reorderMode=REORDER_FALSE;
	
	/*--------------------------------------------------------------*/

	public static long lastReadsIn=-1;
	public static long lastBasesIn=-1;
	public static long lastReadsOut=-1;
	public static long lastBasesOut=-1;
	
	static boolean quantizeQuality=false;
	static final int REORDER_FALSE=0, REORDER_CONSENSUS=1, REORDER_PAIRED=2, REORDER_AUTO=3;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	PrintStream outstream=System.err;
	public static boolean verbose=true;
	public static boolean doHashAndSplit=true;
	public boolean errorState=false;
	boolean overwrite=true;
	boolean append=false;
	
}
