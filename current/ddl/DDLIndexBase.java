package ddl;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicLong;

import cardinality.DynamicDemiLog;

/**
 * Abstract base for the DDL sketch inverted index (maps a (bucket,value) pair to the
 * clade IDs holding that value in that bucket). Holds everything that is independent of
 * the storage layout: the query-ranking (topHits), the stats accounting, ANI, and the
 * bulk-build thread fan-out. Subclasses supply the storage plus its two hot primitives:
 *   accumulateCounts() -- the query inner loop, and
 *   buildBucketRange() -- the build inner loop.
 *
 * Three implementations, selectable via {@link #create(int)} / {@link #USE_CSR2} / {@link #USE_CSR}:
 *   CSRIndex2   -- IDs and offsets packed as 21-bit fields (3 per long); ~1/3 smaller than CSR,
 *                  bit-identical results, ceiling 2^21-1 sketches. THE DEFAULT.
 *   DDLIndexCSR -- packed int[] + int[] offsets per bucket; drops per-cell object overhead (csr2=f).
 *   DDLIndex    -- int[bucket][value][] matrix; simple, exact; the REFERENCE/fallback (csr=f).
 *
 * @author Brian Bushnell, Noire
 * @date August 3, 2026
 */
public abstract class DDLIndexBase {

	protected DDLIndexBase(int buckets, int values){
		this.buckets=buckets;
		this.values=values;
	}

	/** Selects storage: matrix reference if csr=f, else CSRIndex2 (default) or 32-bit CSR if csr2=f. */
	public static DDLIndexBase create(int buckets){return create(buckets, VALUES);}
	public static DDLIndexBase create(int buckets, int values){
		if(!USE_CSR){return new DDLIndex(buckets, values);}//csr=f -> matrix reference
		if(USE_CSR2){return new CSRIndex2(buckets, values);}//default
		return new DDLIndexCSR(buckets, values);//csr2=f -> 32-bit CSR
	}

	/*--------------------------------------------------------------*/
	/*----------------          Bulk Build          ----------------*/
	/*--------------------------------------------------------------*/

	/** Builds the whole index from records, fanning buckets across threads.
	 *  filledBuckets/numClades setup is shared; the per-bucket storage fill is delegated. */
	public void addAll(ArrayList<DDLRecord> records, int threads){
		if(records==null || records.isEmpty()){return;}
		final int numRecords=records.size();

		numClades=numRecords;
		if(numClades>filledBuckets.length){filledBuckets=new int[numClades];}
		for(int i=0; i<numRecords; i++){
			filledBuckets[i]=records.get(i).ddl.filledBuckets();
		}

		final int actualThreads=Math.min(threads, buckets);
		if(actualThreads<2){
			buildBucketRange(records, 0, buckets);
			return;
		}

		ExecutorService executor=Executors.newFixedThreadPool(actualThreads);
		ArrayList<Future<?>> futures=new ArrayList<>(actualThreads);
		final int bucketsPerThread=(buckets+actualThreads-1)/actualThreads;
		for(int t=0; t<actualThreads; t++){
			final int bStart=t*bucketsPerThread;
			final int bEnd=Math.min(bStart+bucketsPerThread, buckets);
			futures.add(executor.submit(()->buildBucketRange(records, bStart, bEnd)));
		}
		for(Future<?> f : futures){
			try{f.get();}catch(Exception e){throw new RuntimeException(e);}
		}
		executor.shutdown();
	}

	/** Populates this index's storage for buckets [bStart, bEnd) from all records. */
	protected abstract void buildBucketRange(ArrayList<DDLRecord> records, int bStart, int bEnd);

	/*--------------------------------------------------------------*/
	/*----------------         Query + Rank         ----------------*/
	/*--------------------------------------------------------------*/

	public int[] query(DynamicDemiLog query){return query(query.toAbsoluteArray());}

	/** Match count per cladeID for the query. The stats accounting is shared here; the
	 *  storage-specific walk is delegated to accumulateCounts(). */
	public final int[] query(char[] maxArray){
		final int[] counts=new int[numClades];
		final long work=accumulateCounts(maxArray, counts);
		totalIndexWork.addAndGet(work);
		long candidates=0;
		for(int i=0; i<numClades; i++){
			if(counts[i]>0){candidates++;}
		}
		totalCandidates.addAndGet(candidates);
		return counts;
	}

	/** For each nonzero query bucket value, increment counts[cladeID] for every clade holding
	 *  that value in that bucket. Returns the total number of clade-IDs visited (index work). */
	protected abstract long accumulateCounts(char[] maxArray, int[] counts);

	/** Finds the top N clades by match count.
	 * @return int[][2] of {cladeID, matchCount}, sorted by matchCount descending. */
	public int[][] topHits(DynamicDemiLog query, int maxHits){
		final int[] counts=query(query);
		final int[][] hits=new int[maxHits][2];
		int minVal=0, minIdx=0;
		int filled=0;
		for(int id=0; id<counts.length; id++){
			final int c=counts[id];
			if(c<=minVal && filled>=maxHits){continue;}
			if(filled<maxHits){
				hits[filled][0]=id;
				hits[filled][1]=c;
				filled++;
				if(filled==maxHits){
					minVal=Integer.MAX_VALUE;
					for(int i=0; i<filled; i++){
						if(hits[i][1]<minVal){minVal=hits[i][1]; minIdx=i;}
					}
				}
			}else{
				hits[minIdx][0]=id;
				hits[minIdx][1]=c;
				minVal=Integer.MAX_VALUE;
				for(int i=0; i<filled; i++){
					if(hits[i][1]<minVal){minVal=hits[i][1]; minIdx=i;}
				}
			}
		}
		Arrays.sort(hits, 0, filled, (a, b)->b[1]-a[1]);
		if(filled<maxHits){return Arrays.copyOf(hits, filled);}
		return hits;
	}

	/** ANI between the query and one indexed clade, from its match count. */
	public float ani(int queryMatches, int queryFilled, int cladeID, int k){
		if(queryMatches<=0){return 0;}
		final int refFilled=filledBuckets[cladeID];
		final int minDiv=Math.min(queryFilled, refFilled);
		if(minDiv<=0){return 0;}
		final double c=Math.min(1.0, (double)queryMatches/minDiv);
		return (float)Math.exp(Math.log(c)/k);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Meta             ----------------*/
	/*--------------------------------------------------------------*/

	public int numClades(){return numClades;}
	/** Number of non-empty (bucket,value) cells, for memory diagnostics. Storage-specific. */
	public abstract long populatedCells();
	public long totalIndexWork(){return totalIndexWork.get();}
	public long totalCandidates(){return totalCandidates.get();}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	protected final int buckets, values;
	protected int[] filledBuckets=new int[1024];
	protected int numClades=0;
	protected final AtomicLong totalIndexWork=new AtomicLong(0);
	protected final AtomicLong totalCandidates=new AtomicLong(0);

	/** When false (csr=f) the matrix reference {@link DDLIndex} is used, overriding csr2. True by
	 *  default; leave true for either packed index. */
	public static boolean USE_CSR=true;

	/** THE DEFAULT (true): the 21-bit bit-packed {@link CSRIndex2} -- validated bit-identical to the
	 *  matrix, ~33% smaller than CSR at 32k (24.6->16.3 GB) for ~1% slower compare; ceiling 2^21-1
	 *  sketches. Set csr2=f for the 32-bit {@link DDLIndexCSR}. */
	public static boolean USE_CSR2=true;

	public static final int BUCKETS=4096;
	public static final int VALUES=65536;
}
