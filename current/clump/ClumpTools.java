package clump;

import java.util.ArrayList;

import bloom.KCountArray;
import fileIO.ReadWrite;
import shared.Shared;
import stream.ConcurrentCollectionReadInputStream;
import stream.Read;

/**
 * Utility class for managing k-mer count arrays and retrieving statistical tables
 * for genomic read data.
 * Provides synchronized methods for creating, caching, and managing global k-mer count tables
 * from different input sources including read collections and input files.
 *
 * @author Brian Bushnell
 * @date Nov 12, 2015
 */
public class ClumpTools {
	
	/** Returns the cached table. NOTE: unsynchronized read of the static (build/clear ARE synchronized);
	 * fine for clump's build-then-use flow where the table is set before workers read it. */
	public static KCountArray table(){
		return table;
	}
	
	//Reads-overload: ALWAYS rebuilds (nulls fname1 first), so it never reuses the file-overload's cache and
	//leaves fname1==null -> a later file-overload call rebuilds too. No cross-overload stale-reuse.
	public static synchronized KCountArray getTable(ArrayList<Read> reads, int k, int minCount){
		fname1=fname2=null;
		table=null;
		ConcurrentCollectionReadInputStream cris=new ConcurrentCollectionReadInputStream(reads, null, -1);
		cris.start();
		table=PivotSet.makeKcaStatic(cris, k, minCount, Shared.AMINO_IN);
		ReadWrite.closeStream(cris);
		return table;
	}
	
	//[clump/ClumpTools#001] FIXED (Brian-approved 2026-06-24): the cache key now compares ALL params
	//(fname1_, fname2_, k_, minCount_), not just fname1_. The old key ignored k/minCount/fname2, so reusing
	//the same fname1 with different k/minCount/fname2 (no clearTable between - the table is static and persists
	//across Clumpify runs in one JVM, e.g. RQCFilter) returned the STALE table -> wrong dedup/error-correction.
	//Strictly safety-positive: a fuller key only ever forces a CORRECT rebuild, never a wrong reuse; it is a
	//no-op for normal single-run usage (where fname1=in1 already changes per pass). VERIFIED on mruber.fq.gz
	//(694MB/10.7M reads, groups=11 dedupe): identical read set (sorted-content md5 27c818b0 == baseline) and
	//same speed (24.9-25.7s, within noise). NB raw file md5 varies run-to-run from clumpify's parallel group
	//ordering even with unchanged code, so the order-independent (sorted) md5 is the correctness fingerprint.
	public static synchronized KCountArray getTable(String fname1_, String fname2_, int k_, int minCount_){
		final boolean fname2Same=(fname2==null ? fname2_==null : fname2.equals(fname2_));
		if(fname1==null || !fname1.equals(fname1_) || k!=k_ || minCount!=minCount_ || !fname2Same || table==null){
			fname1=fname1_;
			fname2=fname2_;
			k=k_;
			minCount=minCount_;
			String[] args=new String[] {"in1="+fname1, "in2="+fname2, "k="+k_, "minCount="+minCount_};
			table=PivotSet.makeSet(args);
		}
		return table;
	}

	public static synchronized void clearTable() {
		fname1=fname2=null;
		k=minCount=-1;
		table=null;
	}

	/**
	 * Cached secondary input file name (paired mate), if provided for the table.
	 */
	private static String fname1=null, fname2=null;
	/** Cached k and minCount the current table was built with (part of the cache key, see #001) */
	private static int k=-1, minCount=-1;
	private static KCountArray table=null;
	
}
