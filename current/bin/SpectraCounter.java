package bin;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.locks.ReadWriteLock;

import prok.CallGenes;
import prok.GeneCaller;
import prok.Orf;
import shared.KillSwitch;
import shared.LineParserS1;
import shared.LineParserS4;
import shared.Shared;
import shared.Tools;
import stream.Streamer;
import stream.Read;
import structures.IntLongHashMap;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.EntropyTracker;
import tracker.KmerTracker;

/**
 * Concurrent processing and analysis of genomic contigs with multi-threaded feature extraction.
 * Manages parallel processing of contigs, extracting genomic metrics like entropy, depth, and potential gene markers.
 * Supports optional entropy and strandedness calculation, 16S/18S ribosomal gene detection, and taxonomic ID parsing.
 *
 * @author Brian Bushnell
 * @date December 2013
 */
public class SpectraCounter extends BinObject implements Accumulator<SpectraCounter.LoadThread> {
	
	/**
	 * Constructs a SpectraCounter with configuration for genomic contig processing.
	 * Initializes entropy adjustment parameters if entropy calculation is enabled.
	 *
	 * @param outstream_ Output stream for progress and error messages
	 * @param parseDepth_ Whether to parse depth information from contig headers
	 * @param parseTID_ Whether to parse taxonomy IDs from contig headers
	 * @param sizeMap_ Map for tracking contig sizes by taxonomy ID
	 */
	public SpectraCounter(PrintStream outstream_, boolean parseDepth_, 
			boolean parseTID_, IntLongHashMap sizeMap_) {
		outstream=outstream_;
		parseDepth=parseDepth_;
		parseTID=parseTID_;
		sizeMap=sizeMap_;
		if(calcEntropy) {
			if(AdjustEntropy.kLoaded!=4 || AdjustEntropy.wLoaded!=150) {
				AdjustEntropy.load(4, 150);
			}
			assert(AdjustEntropy.kLoaded==4 && AdjustEntropy.wLoaded==150) : 
				AdjustEntropy.kLoaded+", "+calcEntropy;
		}
	}
	
	/** Spawn process threads */
	public void makeSpectra(ArrayList<Contig> contigs, Streamer cris, int minlen){
		
		//Do anything necessary prior to processing
//		sizeMap=(parseTax ? new IntLongHashMap(1021) : null);
		
		GeneTools.setMode(call16S, call18S, false, false, false, false);
		//Determine how many threads may be used
		int threads=Tools.mid(1, cris==null ? contigs.size()/4 : 128, Shared.threads());
		if(loadThreadsOverride>0) {threads=loadThreadsOverride;}
		//Fill a list with LoadThreads
		ArrayList<LoadThread> alpt=new ArrayList<LoadThread>(threads);
		for(int i=0; i<threads; i++){
			LoadThread lt=new LoadThread(contigs, cris, minlen, i, threads);
			alpt.add(lt);
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		
	}
	
	@Override
	public synchronized void accumulate(LoadThread t) {
		synchronized(t) {
			errorState|=(t.success);
			contigsLoaded+=t.contigsLoadedT;
			basesLoaded+=t.basesLoadedT;
			contigsRetained+=t.contigsRetainedT;
			basesRetained+=t.basesRetainedT;
		}
	}

	@Override
	public ReadWriteLock rwlock() {return null;}

	@Override
	public synchronized boolean success() {return errorState;}
	
	/**
	 * Worker thread for concurrent contig processing and analysis.
	 * Handles loading contigs from input streams, processing genomic features,
	 * and tracking statistics for the parent SpectraCounter.
	 */
	class LoadThread extends Thread {
		
		/**
		 * Constructs a LoadThread for processing genomic contigs.
		 *
		 * @param contigs_ Shared list of contigs to process or populate
		 * @param cris_ Concurrent read input stream for loading contigs (may be null)
		 * @param minlen_ Minimum contig length threshold
		 * @param tid_ Thread ID for work distribution
		 * @param threads_ Total number of processing threads
		 */
		LoadThread(ArrayList<Contig> contigs_, Streamer cris_, 
				int minlen_, int tid_, int threads_) {
			contigs=contigs_;
			cris=cris_;
			minlen=minlen_;
			tid=tid_;
			threads=threads_;
		}
		
		@Override
		public void run() {
			et=new EntropyTracker(entropyK, entropyWindow, false);
			if(call16S || call18S) {
				caller=GeneTools.makeGeneCaller();
			}
			synchronized(this) {
				runInner();
			}
		}
		
		/** Core processing logic that routes to appropriate processing method.
		 * Processes either pre-loaded contigs or streams new contigs from input. */
		private void runInner() {
			if(cris==null) {//Calculate data on existing contigs
				runOnContigs();
			}else {//Load contigs concurrently
				runOnCris();
			}
			success=true;
		}
		
		/** Processes pre-loaded contigs using thread-based work distribution.
		 * Each thread processes every nth contig based on thread ID and total thread count. */
		void runOnContigs() {
			for(int i=tid; i<contigs.size(); i+=threads) {
				Contig c=contigs.get(i);
				processContig(c);
			}
		}
		
		/**
		 * Processes contigs from concurrent read input stream.
		 * Continuously fetches read lists, converts reads to contigs, processes them,
		 * and adds processed contigs to the shared contig list.
		 */
		void runOnCris() {
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			
			//As long as there is a nonempty read list...
			while(ln!=null){
				ArrayList<Contig> localContigs=new ArrayList<Contig>(ln.size());
				for(Read r : ln) {
					Contig c=loadContig(r);
					if(c!=null) {
						processContig(c);
						localContigs.add(c);
					}
				}
				synchronized(contigs) {contigs.addAll(localContigs);}
				
				//Fetch a new list
				ln=cris.nextList();
			}
		}
		
		/**
		 * Converts a read to a contig and processes it completely.
		 * @param r The read to convert and process
		 * @return The processed contig, or null if the read was too short
		 */
		Contig processRead(Read r) {
			Contig c=loadContig(r);
			if(c==null) {return null;}
			processContig(c);
			return c;
		}
		
		/**
		 * Loads a read into a Contig object with metadata extraction.
		 * Parses taxonomy ID if enabled, filters by minimum length, detects ribosomal genes,
		 * and creates a fully initialized Contig with associated metadata.
		 *
		 * @param r The read to convert to a contig
		 * @return A new Contig object with metadata, or null if read is too short
		 */
		Contig loadContig(Read r) {
			contigsLoadedT++;
			basesLoadedT+=r.length();
			int tid=-1;
			if(parseTID) {
				tid=parseTaxID(r.name());
				if(tid>0) {
					synchronized(sizeMap) {
						sizeMap.increment(tid, r.length());
					}
				}
			}
			if(r.length()<minlen) {return null;}
			contigsRetainedT++;
			basesRetainedT+=r.length();
			Contig c=new Contig(r.name(), r.bases, (int)r.numericID);
			byte[][] ssu=callSSU(r);
			synchronized(c) {
				c.labelTaxid=tid;
				if(ssu!=null) {
					c.r16S=ssu[0];
					c.r18S=ssu[1];
				}
			}
			return c;
		}
		
		/**
		 * Detects and extracts 16S or 18S ribosomal RNA genes from a read.
		 * Uses gene calling to identify ribosomal sequences and returns the first
		 * 16S or 18S sequence found.
		 *
		 * @param r The read to analyze for ribosomal genes
		 * @return Array containing [16S_sequence, 18S_sequence] or null if none found
		 */
		byte[][] callSSU(Read r) {
			if(caller==null || r.length()<900) {return null;}
			assert(call16S || call18S);
			ArrayList<Orf> genes=caller.callGenes(r);
			if(genes==null || genes.isEmpty()) {return null;}
			byte[] r16s, r18s;
			for(Orf orf : genes) {
				if(orf.is16S()) {
					r16s=CallGenes.fetch(orf, r).bases;
					return new byte[][] {r16s, null};
				}else if(orf.is18S()) {
					r18s=CallGenes.fetch(orf, r).bases;
					return new byte[][] {null, r18s};
				}
			}
			return null;
		}
		
		/**
		 * Processes a contig by calculating various genomic metrics.
		 * Computes k-mer counts, normalized depth, entropy (if enabled), strandedness
		 * (if enabled), and parses depth information from headers. All processing
		 * is synchronized on the contig object.
		 *
		 * @param c The contig to process and analyze
		 */
		void processContig(Contig c) {
			synchronized(c) {
//				System.err.println("Thread "+tid+" got lock on "+c.name+", "+c.id()+", "+c.size());
				contigsProcessedT++;
				basesProcessedT+=c.size();
				c.loadCounts();
				if(c.numDepths()>1) {c.fillNormDepth();}
				if(calcEntropy) {
					if(!calcEntropyFast) {
						c.entropy=et.averageEntropy(c.bases, false);
					}else {
						//This would need regeneration of nns but uses much less CPU time loading.
						c.entropy=EntropyTracker.calcEntropyFromCounts(c.trimers);
					}
					c.entropy=AdjustEntropy.compensate(c.gc(), c.entropy);
				}
				if(calcStrandedness) {
					c.dimers=new int[16];
					c.strandedness=EntropyTracker.strandedness(c.bases, c.dimers, 2);
					c.hh=KmerTracker.HH(c.dimers);
					c.caga=KmerTracker.CAGA(c.dimers);
				}
				if(parseDepth) {
					boolean b=DataLoader.parseAndSetDepth(c, lps, lpt);
					if(!b) {
						KillSwitch.kill("Could not parse depth from header "+c.name+
								"\nThis program needs a sam file, a cov file, or labeled contigs.");
					}
					assert(b) : "Could not parse depth from "+c.name;
				}
				
				assert(c.tetramers!=null && c.numTetramers>0);
			}
		}
		
		/** Thread ID for work distribution among parallel threads */
		final int tid;
		/** Total number of processing threads for workload distribution */
		final int threads;
		/** Minimum contig length threshold for retention */
		final int minlen;
		/** Shared list of contigs to process or populate */
		final ArrayList<Contig> contigs;
		/** Streamer for loading contigs on demand */
		final Streamer cris;
		/** Entropy tracker for calculating sequence entropy metrics */
		private EntropyTracker et;
		/** Gene caller for detecting ribosomal RNA genes */
		private GeneCaller caller;
//		final int[] counts=(calcEntropy ? new int[1<<(entropyK*2)] : null);
		/** Thread processing success status */
		boolean success=false;
		/** Number of contigs processed by this thread */
		int contigsProcessedT=0;
		/** Total bases processed by this thread */
		long basesProcessedT=0;
		/** Line parser for single-character delimited strings */
		LineParserS1 lps=new LineParserS1('_');
		/** Line parser for four-character delimited strings */
		LineParserS4 lpt=new LineParserS4(",,=,");
		
		/** Number of contigs loaded by this thread */
		int contigsLoadedT=0;
		/** Total bases loaded by this thread */
		long basesLoadedT=0;
		/** Number of contigs retained after length filtering by this thread */
		int contigsRetainedT=0;
		/** Total bases retained after length filtering by this thread */
		long basesRetainedT=0;
	}
	
	/** Output stream for progress messages and error reporting */
	public PrintStream outstream=System.err;
	
	/** Whether to parse depth information from contig headers */
	public final boolean parseDepth;
	/** Whether to parse taxonomy IDs from contig headers */
	public final boolean parseTID;
	/** Map tracking contig sizes by taxonomy ID */
	public final IntLongHashMap sizeMap;

	/** Total number of contigs loaded across all threads */
	public int contigsLoaded=0;
	/** Total bases loaded across all threads */
	public long basesLoaded=0;
	/** Total number of contigs retained after filtering across all threads */
	public int contigsRetained=0;
	/** Total bases retained after filtering across all threads */
	public long basesRetained=0;
	
	/** Processing error status across all threads */
	public boolean errorState=false;
	/** Global flag to enable entropy calculation for processed contigs */
	public static boolean calcEntropy=true;
	/** Global flag to use fast entropy calculation method from k-mer counts */
	public static boolean calcEntropyFast=false;
	/** Global flag to enable strandedness calculation for processed contigs */
	public static boolean calcStrandedness=true;
	/** Global flag to enable 16S ribosomal RNA gene detection */
	public static boolean call16S=false;
	/** Global flag to enable 18S ribosomal RNA gene detection */
	public static boolean call18S=false;
	/**
	 * Override for the number of loading threads, -1 uses automatic determination
	 */
	public static int loadThreadsOverride=-1;
	
}
