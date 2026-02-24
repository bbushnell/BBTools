package clade;

import java.util.ArrayList;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import bin.AdjustEntropy;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.EntropyTracker;

/**
 * Loads a single file using multiple threads sharing one ConcurrentReadInputStream.
 * Designed for large single-file inputs (e.g., large genome assemblies or bins)
 * where multithreading within a single file provides speedup.
 * In non-perContig mode, each thread accumulates k-mer counts into a private
 * partial Clade; results are merged after all threads complete.
 * In perContig mode, each thread produces individual Clades per sequence;
 * results are collected and sorted by numericID to preserve input order.
 *
 * @author Chloe
 * @date February 23, 2026
 */
public class CladeLoaderSF extends CladeObject implements Accumulator<CladeLoaderSF.ProcessThread> {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Default constructor. Loads entropy model if needed.
	 */
	public CladeLoaderSF() {
		if(AdjustEntropy.kLoaded!=4 || AdjustEntropy.wLoaded!=150) {
			AdjustEntropy.load(4, 150);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Load clades from a single file using multiple threads.
	 * For clade-format files, delegates to CladeLoaderMF.loadCladesFromClade().
	 * For non-perContig mode, all reads are merged into one Clade.
	 * For perContig mode, returns one Clade per sequence, in input order.
	 * @param fname File to load
	 * @param perContig Whether to treat each contig as a separate clade
	 * @param minContig Minimum contig length to process
	 * @param maxReads Maximum reads to process (-1 means no limit)
	 * @param finish Whether to call finish() on completed Clades
	 * @return List of Clades loaded from the file
	 */
	public ArrayList<Clade> loadFile(String fname, boolean perContig,
			int minContig, long maxReads, boolean finish) {
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, false);
		if(ff.clade()) {
			return CladeLoaderMF.loadCladesFromClade(ff, maxReads);
		}
		return spawnThreads(ff, perContig, minContig, maxReads, finish);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Opens a ConcurrentReadInputStream for the file, spawns threads to process it,
	 * then merges results.
	 */
	private ArrayList<Clade> spawnThreads(FileFormat ff, boolean perContig,
			int minContig, long maxReads, boolean finish) {

		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff, null);
		cris.start();

		final int threads=Shared.threads();
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++) {
			alpt.add(new ProcessThread(cris, i, perContig, minContig, finish));
		}

		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		ReadWrite.closeStream(cris);

		ArrayList<Clade> result=new ArrayList<Clade>();

		if(!perContig) {
			//Merge all partial Clades into one.
			//Manual merge bypasses the taxID>0 assertion in Clade.add(Clade),
			//which is designed for merging known-organism Clades.
			Clade merged=new Clade(-1, -1, ff.simpleName());
			for(ProcessThread pt : alpt) {
				Clade partial=pt.partialClade;
				if(partial.bases==0) {continue;}
				Tools.add(merged.counts, partial.counts);
				if(merged.bases+partial.bases>0) {
					merged.entropy=(merged.entropy*merged.bases+partial.entropy*partial.bases)
							/(float)(merged.bases+partial.bases);
				}
				merged.bases+=partial.bases;
				merged.contigs+=partial.contigs;
			}
			if(finish && merged.bases>0) {merged.finish();}
			if(merged.bases>0) {result.add(merged);}
		} else {
			//Collect per-contig Clades from all threads and sort by numericID
			//to restore input order.
			ArrayList<long[]> entries=new ArrayList<long[]>(); // {numericID, threadIdx, cladeIdx}
			for(int t=0; t<alpt.size(); t++) {
				ProcessThread pt=alpt.get(t);
				for(int j=0; j<pt.perContigClades.size(); j++) {
					entries.add(new long[]{pt.numericIDs.get(j), t, j});
				}
			}
			entries.sort((a, b) -> Long.compare(a[0], b[0]));
			for(long[] entry : entries) {
				result.add(alpt.get((int)entry[1]).perContigClades.get((int)entry[2]));
			}
		}

		return result;
	}

	/**
	 * Accumulate statistics from a completed ProcessThread.
	 */
	@Override
	public final void accumulate(ProcessThread pt) {
		synchronized(pt) {
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			errorState|=(!pt.success);
		}
	}

	/**
	 * Check if processing was successful.
	 */
	@Override
	public final boolean success() {return !errorState;}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Thread that reads batches from a shared ConcurrentReadInputStream.
	 * In non-perContig mode, accumulates reads into a private partial Clade.
	 * In perContig mode, creates one Clade per sequence and tracks numericIDs.
	 */
	static class ProcessThread extends Thread {

		ProcessThread(ConcurrentReadInputStream cris_, int tid_,
				boolean perContig_, int minContig_, boolean finish_) {
			cris=cris_;
			tid=tid_;
			perContig=perContig_;
			minContig=minContig_;
			finish=finish_;
			partialClade=new Clade(-1, -1, "partial_"+tid_);
		}

		@Override
		public void run() {
			ListNum<Read> ln=cris.nextList();
			while(ln!=null && ln.size()>0) {
				for(Read r : ln) {
					if(r.bases==null || r.bases.length<minContig) {continue;}
					if(!perContig) {
						partialClade.add(r.bases, et);
					} else {
						int taxID=CladeObject.resolveTaxID(r.id);
						Clade c=new Clade(taxID<1 ? -1 : taxID, -1, r.id);
						c.add(r.bases, et);
						if(finish) {c.finish();}
						perContigClades.add(c);
						numericIDs.add(r.numericID);
					}
					readsProcessedT++;
					basesProcessedT+=r.bases.length;
				}
				cris.returnList(ln);
				ln=cris.nextList();
			}
			if(ln!=null) {cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());}
			success=true;
		}

		/** Number of reads processed by this thread */
		long readsProcessedT=0;
		/** Number of bases processed by this thread */
		long basesProcessedT=0;
		/** True only if this thread completed successfully */
		boolean success=false;

		/** Partial Clade accumulating k-mer counts in non-perContig mode */
		final Clade partialClade;
		/** Per-contig Clades produced in perContig mode */
		final ArrayList<Clade> perContigClades=new ArrayList<Clade>();
		/** numericIDs parallel to perContigClades, for sorting into input order */
		final ArrayList<Long> numericIDs=new ArrayList<Long>();

		/** Shared input stream - ConcurrentReadInputStream is thread-safe */
		private final ConcurrentReadInputStream cris;
		/** Thread ID */
		final int tid;
		private final boolean perContig;
		private final int minContig;
		private final boolean finish;
		/** Per-thread entropy tracker - not thread-safe, must not be shared */
		private final EntropyTracker et=new EntropyTracker(entropyK, entropyWindow, false);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	public long readsProcessed=0;
	/** Number of bases processed */
	public long basesProcessed=0;
	/** True if an error was encountered */
	public boolean errorState=false;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();

}
