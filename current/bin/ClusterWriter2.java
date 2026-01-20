package bin;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Multithreaded version of ClusterWriter.
 * Writes clusters to files using a fixed number of threads to improve I/O throughput.
 * Uses asynchronous job submission for the shared 'chaff' file to avoid lock contention.
 */
public class ClusterWriter2 {
	
	public ClusterWriter2(PrintStream outstream_, boolean writeChaff_, boolean loud_, 
			boolean overwrite_, boolean append_, int threads_, boolean sorted_){
		outstream=outstream_;
		overwrite=overwrite_;
		append=append_;
		writeChaff=writeChaff_;
		loud=loud_;
		threads=Math.max(1, threads_);
		sorted=sorted_;
	}

	/**
	 * Writes clustered contigs to output files using specified pattern.
	 * Supports individual bin files with % substitution or single combined output.
	 * Optionally writes small clusters to a separate chaff file.
	 * Uses multiple threads to write files in parallel.
	 *
	 * @param pattern Output file pattern with % for bin numbering
	 * @param clusters List of bins/clusters to output
	 * @param minBases Minimum base count threshold for main output
	 * @param minContigs Minimum contig count threshold for main output
	 */
	void outputClusters(String pattern, String sizeHist, ArrayList<? extends Bin> clusters,
			long minBases, int minContigs, Binner binner, DataLoader loader) {
		Timer t=new Timer();

		if(pattern!=null) {
			if(pattern.contains(".") && !pattern.contains("%") && pattern.contains("#")) {
				pattern=Tools.replaceLastInstanceOf(pattern, '#', '%');//A common mistake I make
				//Maybe prudent to make sure it is in a filename rather than dirname
			}
			if(!pattern.contains(".") && !pattern.contains("%")) {
				if(!pattern.endsWith("/")) {pattern=pattern+"/";}
				pattern=pattern+"bin_%.fa";
			}
			
			outstream.println("Writing clusters to "+pattern);
		}
		
		// Reset stats
		clustersWritten.set(0);
		contigsWritten.set(0);
		basesWritten.set(0);
		
		if(pattern!=null && pattern.indexOf('%')>=0) {
			// Multithreaded mode for multiple files
			
			ByteStreamWriter chaff=null;
			if(writeChaff) {
				chaff=ByteStreamWriter.makeBSW(pattern.replaceFirst("%", "chaff"), overwrite, append, true);
			}
			
			final AtomicInteger nextIndex=new AtomicInteger(0);
			ArrayList<WriteThread> threadList=new ArrayList<WriteThread>(threads);
			
			for(int i=0; i<threads; i++){
				WriteThread wt=new WriteThread(nextIndex, clusters, pattern, minBases, minContigs, chaff, i);
				threadList.add(wt);
				wt.start();
			}
			
			// Wait for completion
			for(WriteThread wt : threadList){
				try {
					wt.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			
			if(chaff!=null) {chaff.poisonAndWait();}
			
		}else {
			// Single file output - stick to single thread to maintain order and simplicity
			final ByteBuilder bb=new ByteBuilder(65536);
			final ByteStreamWriter bsw=ByteStreamWriter.makeBSW(pattern, overwrite, append, true);
			for(int i=0; i<clusters.size(); i++) {
				Bin a=clusters.get(i);
				printBin(a, bsw, bb, i+1);
				clustersWritten.incrementAndGet();
				contigsWritten.addAndGet(a.numContigs());
				basesWritten.addAndGet(a.size());
			}
			if(bsw!=null) {bsw.poisonAndWait();}
		}
		
		float cpct=contigsWritten.get()*100f/loader.contigsLoaded;
		float bpct=basesWritten.get()*100f/loader.basesLoaded;
		
		if(loud) {
			outstream.println("\nMetric   \t        In\t       Out\tPercent");
			outstream.println("Clusters \t"+Tools.padLeft(0, 10)+"\t"+Tools.padLeft(clustersWritten.get(), 10));
			outstream.println("Contigs  \t"+Tools.padLeft(loader.contigsLoaded, 10)+"\t"+
					Tools.padLeft(contigsWritten.get(), 10)+"\t"+Tools.format("%.2f%%", cpct));
			outstream.println("Bases    \t"+Tools.padLeft(loader.basesLoaded, 10)+"\t"+
					Tools.padLeft(basesWritten.get(), 10)+"\t"+Tools.format("%.2f%%", bpct));
		}
		
		if(sizeHist!=null) {
			final ByteStreamWriter bsw=ByteStreamWriter.makeBSW(sizeHist, overwrite, append, true);
			bsw.print("#Size\tGood\tBad\n");
			for(int i=0; i<binner.goodMergeSize.length; i++) {
				long size=1L<<i;
				bsw.print(size).tab().print(binner.goodMergeSize[i]).tab().print(binner.badMergeSize[i]).nl();
			}
			bsw.poisonAndWait();
		}
		t.stop("Time:");
	}
	
	/** Worker thread for processing clusters */
	private class WriteThread extends Thread {
		
		WriteThread(AtomicInteger nextIndex_, ArrayList<? extends Bin> clusters_, String pattern_, 
				long minBases_, int minContigs_, ByteStreamWriter chaff_, int tid_){
			nextIndex=nextIndex_;
			clusters=clusters_;
			pattern=pattern_;
			minBases=minBases_;
			minContigs=minContigs_;
			chaff=chaff_;
			tid=tid_;
		}
		
		@Override
		public void run(){
			// Buffer for individual bin files (reused and cleared)
			final ByteBuilder binBB=new ByteBuilder(65536);
			
			// Buffer for chaff (accumulated and swapped)
			final ByteBuilder chaffBB=(chaff==null ? null : new ByteBuilder(65536));
			
			int i=nextIndex.getAndIncrement();
			
			while(i<clusters.size()){
				Bin a=clusters.get(i);
				if(a.size()>=minBases && a.numContigs()>=minContigs) {
					// Case 1: Write to individual file
					String fname=pattern.replaceFirst("%", Integer.toString(i));
					// Handle variable substitutions
					if(fname.contains("%contam")) {fname=fname.replaceFirst("%contam", String.format("con%.4f", a.contam));}
					else if(fname.contains("contam%")) {fname=fname.replaceFirst("contam%", String.format("con%.4f", a.contam));}
					if(fname.contains("%comp")) {fname=fname.replaceFirst("%comp", String.format("com%.4f", a.completeness));}
					else if(fname.contains("comp%")) {fname=fname.replaceFirst("comp%", String.format("com%.4f", a.completeness));}
					
					// Synchronous write for exclusive file
					final ByteStreamWriter bsw=ByteStreamWriter.makeBSW(fname, overwrite, append, a.size()>=MIN_SIZE_FOR_SUBPROCESS);
					printBin(a, bsw, binBB, -1); // Uses standard print/clear logic
					bsw.poisonAndWait();
					
					clustersWritten.incrementAndGet();
					contigsWritten.addAndGet(a.numContigs());
					basesWritten.addAndGet(a.size());
				}else {
					// Case 2: Write to shared Chaff file
					if(chaff!=null){
						// Use async logic; chaffBB is updated (swapped) if flushed
						printBinChaff(a, chaff, chaffBB, i+1);
					}
					if(sorted && tid>0) {break;}
				}
				
				i=nextIndex.getAndIncrement();
			}
			
			// Final flush for any remaining chaff data
			if(chaff!=null && chaffBB!=null && chaffBB.length()>0) {
				assert(false); //Should not happen
				chaff.addJob(chaffBB);
			}
		}
		
		final AtomicInteger nextIndex;
		final ArrayList<? extends Bin> clusters;
		final String pattern;
		final long minBases;
		final int minContigs;
		final ByteStreamWriter chaff;
		final int tid;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Standard Methods      ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Standard synchronous printing for exclusive files.
	 * Clears the ByteBuilder after writing.
	 */
	private void printBin(Bin a, ByteStreamWriter bsw, ByteBuilder bb, int id) {
		if(a.isCluster()) {printCluster((Cluster)a, bsw, bb, id);}
		else {printContig((Contig)a, bsw, bb, id);}
	}
	
	private void printCluster(Cluster a, ByteStreamWriter bsw, ByteBuilder bb, int id) {
		ArrayList<Contig> contigs=a.contigs;
		Collections.sort(contigs);
		for(Contig c : contigs) {
			c.appendTo(bb, id);
			if(bb.length>=65536) {
				if(bsw!=null) {bsw.print(bb);}
				bb.clear();
			}
		}
		if(bsw!=null && !bb.isEmpty()) {bsw.print(bb);}
		bb.clear();
	}
	
	private void printContig(Contig c, ByteStreamWriter bsw, ByteBuilder bb, int id) {
		c.appendTo(bb, id);
		if(bsw!=null && !bb.isEmpty()) {bsw.print(bb);}
		bb.clear();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Chaff Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Asynchronous printing for shared Chaff file.
	 * Returns the active ByteBuilder (which may be a new instance if the old one was flushed).
	 * Does NOT clear the buffer implicitly unless flushed to BSW.
	 */
	private void printBinChaff(Bin a, ByteStreamWriter bsw, ByteBuilder bb, int id) {
		if(a.isCluster()) {printClusterChaff((Cluster)a, bsw, bb, id);}
		else {printContigChaff((Contig)a, bsw, bb, id);}
	}
	
	private void printClusterChaff(Cluster a, ByteStreamWriter bsw, ByteBuilder bb, int id) {
		ArrayList<Contig> contigs=a.contigs;
		Collections.sort(contigs);
		bb.clear();
		for(Contig c : contigs) {
			c.appendTo(bb, id);
		}
		ByteBuilder bb2=new ByteBuilder(bb.toBytes());
		bsw.addJob(bb2);
		bb.clear();
	}
	
	private void printContigChaff(Contig c, ByteStreamWriter bsw, ByteBuilder bb, int id) {
		bb.clear();
		c.appendTo(bb, id);
		ByteBuilder bb2=new ByteBuilder(bb.toBytes());
		bsw.addJob(bb2);
		bb.clear();
	}
	
	private PrintStream outstream;
	
	final AtomicLong clustersWritten=new AtomicLong(0);
	final AtomicLong contigsWritten=new AtomicLong(0);
	final AtomicLong basesWritten=new AtomicLong(0);
	
	final boolean overwrite;
	final boolean append;
	final boolean writeChaff;
	final boolean loud;
	final boolean sorted;
	final int threads;
	
	static int MIN_SIZE_FOR_SUBPROCESS=1000000;
}