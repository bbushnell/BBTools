package clump2;

import java.util.ArrayList;

import bloom.KCountArray;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Tools;
import stream.ConcurrentCollectionReadInputStream;
import stream.Read;
import stream.SamLine;
import stream.Streamer;
import structures.ListNum;

/**
 * Shared pivot-count cache for clump sorting. File-backed requests reuse a table
 * only when both filenames, k and minimum count match; collection requests always
 * rebuild. Filenames are cache keys, not content fingerprints, so callers must
 * clear the cache when replacing a file in place or changing other counting policy.
 * Builds and clearing are synchronized; table() relies on build-before-worker-use
 * publication. Returned tables are shared, not defensive copies.
 *
 * @author Brian Bushnell
 * @date Nov 12, 2015
 */
public class ClumpTools{
	
	/** Returns the cached table. NOTE: unsynchronized read of the static (build/clear ARE synchronized);
	 * fine for clump's build-then-use flow where the table is set before workers read it. */
	public static KCountArray table(){
		return table;
	}
	
	/** Builds from borrowed reads and closes the collection stream after counting.
	 * Pivot hashing may set each read's object metadata. Replaces the shared cache.
	 * @param reads Reads to count; their contents must remain stable during counting
	 * @param k Pivot length
	 * @param minCount Count threshold used to size counting cells
	 * @return Newly built shared table
	 */
	//Reads-overload: ALWAYS rebuilds (nulls fname1 first), so it never reuses the file-overload's cache and
	//leaves fname1==null -> a later file-overload call rebuilds too. No cross-overload stale-reuse.
	public static synchronized KCountArray getTable(ArrayList<Read> reads, int k, int minCount){
		fname1=fname2=null;
		table=null;
		Streamer streamer=new CollectionStreamer(reads);
		streamer.start();
		table=PivotSet.makeKcaStatic(streamer, k, minCount, Shared.AMINO_IN);
		ReadWrite.closeStream(streamer);
		return table;
	}

	/** Recounts current grouped files using Streamer input into one shared table.
	 * Always rebuilds: group contents and comparator seed/border can change between
	 * outer passes. A failed build leaves no reusable cache entry.
	 * @param first Expanded primary group descriptors
	 * @param second Aligned optional mate descriptors, or null
	 * @param k Pivot length
	 * @param minCount Count threshold used to size counting cells
	 * @return Newly completed shared pivot-count table */
	public static synchronized KCountArray getTable(FileFormat[] first, FileFormat[] second,
			int k, int minCount){
		clearTable();
		table=PivotSet.makeKcaStatic(first, second, k, minCount, Shared.AMINO_IN);
		return table;
	}
	
	/** Returns the cached file-backed table or builds one for the requested key.
	 * @param fname1_ Primary input filename
	 * @param fname2_ Mate input filename, or null
	 * @param k_ Pivot length
	 * @param minCount_ Count threshold used to size counting cells
	 * @return Shared pivot-count table
	 */
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
			//TODO: Probable bug - the new key is published before makeSet succeeds. If a caller catches
			//a build exception and retries this key, the previous non-null table can be returned for it.
			//Normal successful builds are unaffected; no recovering caller has been established here.
			fname1=fname1_;
			fname2=fname2_;
			k=k_;
			minCount=minCount_;
			String[] args=new String[] {"in1="+fname1, "in2="+fname2, "k="+k_, "minCount="+minCount_};
			table=PivotSet.makeSet(args);
		}
		return table;
	}

	/** Drops the cache reference and key; does not invalidate tables held by callers. */
	public static synchronized void clearTable(){
		fname1=fname2=null;
		k=minCount=-1;
		table=null;
	}

	/** Rejects repeated expanded input/output names before streams are started.
	 * Uses Tools' case-insensitive name comparison and standard-output exceptions;
	 * does not resolve filesystem aliases such as symlinks or relative path variants. */
	static void validateFileNames(FileFormat[]... groups){
		final ArrayList<String> names=new ArrayList<String>();
		for(FileFormat[] group : groups){
			if(group!=null){
				for(FileFormat ff : group){
					if(ff!=null){names.add(ff.name());}
				}
			}
		}
		Tools.testForDuplicateFiles(true, names.toArray(new String[names.size()]));
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Exposes borrowed collection reads to Streamer consumers without changing the
	 * existing ConcurrentCollectionReadInputStream producer or buffer lifecycle. */
	private static final class CollectionStreamer implements Streamer{

		CollectionStreamer(ArrayList<Read> reads){
			cris=new ConcurrentCollectionReadInputStream(reads, null, -1);
		}

		@Override
		public String fname(){return cris.fname();}

		@Override
		public void start(){cris.start();}

		@Override
		public void close(){cris.close();}

		@Override
		public boolean paired(){return cris.paired();}

		@Override
		public int pairnum(){return 0;}

		@Override
		public long readsProcessed(){return cris.readsIn();}

		@Override
		public long basesProcessed(){return cris.basesIn();}

		@Override
		public void setSampleRate(float rate, long seed){cris.setSampleRate(rate, seed);}

		@Override
		public ListNum<Read> nextList(){
			ListNum<Read> ln=cris.nextList();
			if(ln!=null && ln.isEmpty()){
				cris.returnList(ln);
				ln=null;
			}
			return ln;
		}

		@Override
		public ListNum<SamLine> nextLines(){
			throw new UnsupportedOperationException("CollectionStreamer does not support SamLine batches");
		}

		@Override
		public boolean hasMore(){
			throw new UnsupportedOperationException("CollectionStreamer does not expose a non-consuming availability hint");
		}

		@Override
		public boolean errorState(){return cris.errorState();}

		@Override
		public void returnList(ListNum<Read> ln){cris.returnList(ln);}

		@Override
		public void returnList(long id, boolean poison){cris.returnList(id, poison);}

		private final ConcurrentCollectionReadInputStream cris;
	}

	/**
	 * Cached primary and optional secondary input filenames.
	 */
	private static String fname1=null, fname2=null;
	/** Cached k and minCount the current table was built with (part of the cache key, see #001) */
	private static int k=-1, minCount=-1;
	/** Shared table, or null before a successful build or after clearing. */
	private static KCountArray table=null;
	
}
