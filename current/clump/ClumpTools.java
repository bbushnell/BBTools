package clump;

import java.util.ArrayList;

import bloom.KCountArray;
import fileIO.ReadWrite;
import shared.Shared;
import stream.ConcurrentCollectionReadInputStream;
import stream.Read;

/**
 * @author Brian Bushnell
 * @date Nov 12, 2015
 *
 */
public class ClumpTools {
	
	/** Returns the current cached k-mer count array table.
	 * @return The cached KCountArray table, or null if no table is currently loaded */
	public static KCountArray table(){
		return table;
	}
	
	/**
	 * Creates or retrieves a k-mer count array from a collection of reads.
	 * Processes the read collection through a concurrent input stream and generates
	 * a k-mer count array using PivotSet with the specified parameters.
	 * This method is thread-safe and clears any existing cached table.
	 *
	 * @param reads ArrayList of Read objects to process for k-mer counting
	 * @param k K-mer size for counting operations
	 * @param minCount Minimum occurrence threshold for k-mer inclusion
	 * @return A new KCountArray containing k-mer frequency data from the input reads
	 */
	public static synchronized KCountArray getTable(ArrayList<Read> reads, int k, int minCount){
		fname1=fname2=null;
		table=null;
		ConcurrentCollectionReadInputStream cris=new ConcurrentCollectionReadInputStream(reads, null, -1);
		cris.start();
		table=PivotSet.makeKcaStatic(cris, k, minCount, Shared.AMINO_IN);
		ReadWrite.closeStream(cris);
		return table;
	}
	
	/**
	 * Creates or retrieves a k-mer count array from input file names with caching.
	 * If the same file names have been processed before and a table exists, returns
	 * the cached table. Otherwise, creates a new table using PivotSet with the
	 * specified input files and parameters. This method is thread-safe.
	 *
	 * @param fname1_ Primary input file name (e.g., forward reads or single-end file)
	 * @param fname2_ Secondary input file name (e.g., reverse reads, may be null)
	 * @param k_ K-mer size for counting operations
	 * @param minCount_ Minimum occurrence threshold for k-mer inclusion
	 * @return KCountArray containing k-mer frequency data, either cached or newly created
	 */
	public static synchronized KCountArray getTable(String fname1_, String fname2_, int k_, int minCount_){
		if(fname1==null || !fname1.equals(fname1_) || table==null){
			fname1=fname1_;
			fname2=fname2_;
			String[] args=new String[] {"in1="+fname1, "in2="+fname2, "k="+k_, "minCount="+minCount_};
			table=PivotSet.makeSet(args);
		}
		return table;
	}
	
	/**
	 * Clears the cached k-mer count table and associated file names.
	 * Resets both the cached table reference and the file name cache,
	 * forcing the next getTable call to create a new table. This method is thread-safe.
	 */
	public static synchronized void clearTable() {
		fname1=fname2=null;
		table=null;
	}
	
	private static String fname1=null, fname2=null;
	/** Cached k-mer count array table for reuse across method calls */
	private static KCountArray table=null;
	
}
