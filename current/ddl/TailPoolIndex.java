package ddl;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.locks.ReadWriteLock;

import structures.IntList;
import template.Accumulator;
import template.ThreadWaiter;

/** Immutable reverse index of a caller-defined, frozen pool of protein-tail sketches.
 * An external family sketch retrieves every tail meeting an exact shared-register
 * floor. This class imposes no top-K, self exclusion, alignment or family policy.
 * Sparse per-bin CSR avoids a fixed width-by-65536 table and a single global
 * posting array. Input sketches are borrowed during construction, not retained.
 * @author Brian Bushnell, Fischl
 */
public final class TailPoolIndex{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Copies the ordinal mapping; callers must keep sketch rows unchanged until this constructor returns.
	 * Pool slots are input-array positions. Original ordinals must be unique and
	 * nonnegative, but need not be sorted or contiguous. Equal sequences remain
	 * separate tails when they have different original ordinals.
	 */
	public TailPoolIndex(final int width_, final int[] originalOrdinals,
			final char[][] tailSketches, final int buildThreads){
		require(width_>0 && width_<=DDLFamilyProfiles.MAX_BUCKETS && Integer.bitCount(width_)==1,
			"Tail pool sketches need a supported power-of-two width");
		require(originalOrdinals!=null && tailSketches!=null && originalOrdinals.length==tailSketches.length,
			"Every tail sketch needs one original ordinal");
		require(buildThreads>0 && tailSketches.length<=Integer.MAX_VALUE-8,
			"Tail pool construction needs positive threads and bounded per-bin arrays");
		width=width_; ordinals=originalOrdinals.clone();
		validateOrdinals(ordinals);
		final char[][] source=tailSketches.clone();
		for(char[] row:source){
			require(row!=null && row.length==width, "Every tail sketch must have the exact index width");
		}
		keys=new char[width][]; offsets=new int[width][]; postings=new int[width][];
		if(source.length==0){
			for(int bin=0; bin<width; bin++){
				keys[bin]=EMPTY_KEYS;
				offsets[bin]=EMPTY_OFFSETS;
				postings[bin]=EMPTY_IDS;
			}
		}else if(buildThreads==1){buildBins(source, 0, width);}
		else{new Construction(source, Math.min(width, buildThreads)).build();}
		long total=0;
		for(int[] bin:postings){total=Math.addExact(total, bin.length);}
		postingCount=total;
	}

	/** The normal ordered pool needs no temporary sort; arbitrary mappings are checked without reordering slots. */
	private static void validateOrdinals(final int[] mapping){
		boolean ordered=true;
		for(int i=0; i<mapping.length; i++){
			require(mapping[i]>=0, "Original tail ordinals cannot be negative");
			if(i>0 && mapping[i]<=mapping[i-1]){ordered=false;}
		}
		if(!ordered){
			final int[] sorted=mapping.clone(); Arrays.sort(sorted);
			for(int i=1; i<sorted.length; i++){require(sorted[i]!=sorted[i-1], "Original tail ordinals must be unique");}
		}
	}

	/** Counts occupied keys, prefixes their slices, then stores each pool slot once per nonzero bin. */
	private void buildBins(final char[][] source, final int first, final int last){
		assert(first>=0 && first<last && last<=width) : "Construction workers must own disjoint valid bin ranges";
		final int[] counts=new int[65536], touched=new int[Math.min(source.length, 65535)];
		for(int bin=first; bin<last; bin++){
			int distinct=0, total=0;
			for(char[] row:source){
				final int value=row[bin];
				if(value==0){continue;}
				if(counts[value]++==0){touched[distinct++]=value;}
				total++;
			}
			Arrays.sort(touched, 0, distinct);
			final char[] binKeys=new char[distinct];
			final int[] binOffsets=new int[distinct+1], binPostings=new int[total];
			for(int i=0; i<distinct; i++){
				final int value=touched[i]; binKeys[i]=(char)value;
				binOffsets[i+1]=binOffsets[i]+counts[value]; counts[value]=binOffsets[i];
			}
			for(int slot=0; slot<source.length; slot++){
				final int value=source[slot][bin];
				if(value!=0){binPostings[counts[value]++]=slot;}
			}
			for(int i=0; i<distinct; i++){
				assert(counts[touched[i]]==binOffsets[i+1]) : "Every counted tail key must fill exactly its CSR slice";
				counts[touched[i]]=0;
			}
			keys[bin]=binKeys; offsets[bin]=binOffsets; postings[bin]=binPostings;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Query            ----------------*/
	/*--------------------------------------------------------------*/

	/** Clears and fills both lists with all qualifying pool slots and their exact hit counts.
	 * Results are ordered by ascending pool slot, not score. minimumHits=0 includes
	 * every tail, including zero-hit tails; a floor above width returns no matches.
	 * No query identity is inferred, so a sketch equal to an indexed row is not
	 * excluded. Output suffixes beyond list.size are unspecified. Callers apply
	 * the family-specific floor, alignment direction, exclusions and deduplication.
	 * Each simultaneous query needs its own workspace and output lists.
	 */
	public int collectHits(final char[] query, final int minimumHits,
			final IntList outPoolSlots, final IntList outHits, final Workspace work){
		require(query!=null && query.length==width && minimumHits>=0,
			"Tail pool query requires the exact width and a nonnegative hit floor");
		require(outPoolSlots!=null && outHits!=null && outPoolSlots!=outHits && outPoolSlots.array!=outHits.array,
			"Tail pool IDs and counts need distinct reusable output storage");
		require(work!=null && !work.inUse, "Tail pool queries need exclusive, non-reentrant workspace");
		work.inUse=true;
		try{
			work.begin(ordinals.length); outPoolSlots.clear(); outHits.clear();
			if(minimumHits>width || ordinals.length==0){return 0;}
			for(int bin=0; bin<width; bin++){
				final char key=query[bin]; if(key==0){continue;}
				final int cell=Arrays.binarySearch(keys[bin], key); if(cell<0){continue;}
				final int end=offsets[bin][cell+1];
				work.postingsVisited+=end-offsets[bin][cell];
				for(int p=offsets[bin][cell]; p<end; p++){
					final int slot=postings[bin][p];
					if(work.counts[slot]++==0){work.touched[work.touchedSize++]=slot;}
				}
			}
			work.touchedTails=work.touchedSize;
			if(minimumHits==0){
				for(int slot=0; slot<ordinals.length; slot++){outPoolSlots.add(slot);}
			}else{
				for(int i=0; i<work.touchedSize; i++){
					final int slot=work.touched[i];
					if(work.counts[slot]>=minimumHits){outPoolSlots.add(slot);}
				}
				outPoolSlots.sort();
			}
			for(int i=0; i<outPoolSlots.size; i++){outHits.add(work.counts[outPoolSlots.array[i]]);}
			assert(outPoolSlots.size==outHits.size) : "Every returned pool slot must retain its exact shared-hit count";
			return outPoolSlots.size;
		}finally{work.end();}
	}

	/** Original identity is explicit and independent of family IDs and sketch equality. */
	public int originalOrdinal(final int poolSlot){
		if(poolSlot<0 || poolSlot>=ordinals.length){throw new IllegalArgumentException("Invalid tail pool slot: "+poolSlot);}
		return ordinals[poolSlot];
	}
	public int size(){return ordinals.length;}
	public int width(){return width;}
	public long postingCount(){return postingCount;}

	/** Retained primitive payload only; excludes object headers, caller rows and temporary/query scratch. */
	public long primitiveBytes(){
		if(ordinals.length==0){return 4;}// The single zero offset array is shared by all empty bins.
		long bytes=4L*ordinals.length;
		for(int bin=0; bin<width; bin++){
			bytes+=2L*keys[bin].length+4L*offsets[bin].length+4L*postings[bin].length;
		}
		return bytes;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Query state         ----------------*/
	/*--------------------------------------------------------------*/

	/** Scratch belongs to one query worker; snapshots are immutable and shareable.
	 * The inUse flag catches reentry, not concurrent misuse. A workspace is never
	 * safe to share between active workers. Counters describe its last valid query.
	 */
	public static final class Workspace{
		/** Grows both scratch arrays together and resets metrics before traversal. */
		private void begin(final int size){
			assert(inUse && touchedSize==0) : "Previous query must clear touched scores before workspace reuse";
			if(counts.length<size){
				final int[] nextCounts=new int[size], nextTouched=new int[size];
				counts=nextCounts; touched=nextTouched;
			}
			postingsVisited=0; touchedTails=0;
		}
		/** Clears only touched scores, including after an unsuccessful query. */
		private void end(){
			for(int i=0; i<touchedSize; i++){counts[touched[i]]=0;}
			touchedSize=0; inUse=false;
		}
		public long postingsVisited(){return postingsVisited;}
		public int touchedTails(){return touchedTails;}
		private int[] counts=EMPTY_IDS, touched=EMPTY_IDS;
		private int touchedSize, touchedTails;
		private long postingsVisited;
		private boolean inUse;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Construction         ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed rows and worker objects become unreachable after all construction workers join. */
	private final class Construction implements Accumulator<BuildWorker>{
		/** Assigns disjoint contiguous bin ranges; each worker has independent counters. */
		Construction(final char[][] source_, final int threads){
			source=source_;
			workers=new ArrayList<BuildWorker>(threads);
			for(int i=0; i<threads; i++){
				workers.add(new BuildWorker(this, (int)((long)width*i/threads), (int)((long)width*(i+1)/threads)));
			}
		}
		/** Joins every started worker before returning or propagating a construction failure. */
		void build(){
			try{
				if(!ThreadWaiter.startAndWait(workers, this)){fail(new IllegalStateException("Tail pool construction did not finish"));}
			}catch(Throwable cause){
				fail(cause);
				final ArrayList<BuildWorker> started=new ArrayList<BuildWorker>();
				for(BuildWorker worker:workers){if(worker.getState()!=Thread.State.NEW){started.add(worker);}}
				ThreadWaiter.waitForThreadsToFinish(started);
			}
			if(failure!=null){throw new IllegalStateException("Tail pool index failed before publication", failure);}
			if(done!=workers.size()){throw new IllegalStateException("Tail pool index omitted a construction range");}
		}
		synchronized void fail(final Throwable cause){if(failure==null){failure=cause;}}
		@Override public void accumulate(final BuildWorker worker){
			if(!worker.success){fail(new IllegalStateException("Tail pool bin worker failed"));}
			done++;
		}
		@Override public boolean success(){return failure==null && done==workers.size();}
		@Override public ReadWriteLock rwlock(){return null;}// Aggregation follows worker joins.
		final char[][] source;
		final ArrayList<BuildWorker> workers;
		volatile Throwable failure;
		int done;
	}

	/** Each worker fills a disjoint bin range and owns its counting scratch. */
	private final class BuildWorker extends Thread{
		/** Holds the source only through its short-lived construction coordinator. */
		BuildWorker(final Construction owner_, final int first_, final int last_){
			super("ddl-tail-pool-"+first_);
			owner=owner_; first=first_; last=last_;
		}
		@Override public void run(){
			try{buildBins(owner.source, first, last); success=true;}
			catch(Throwable cause){owner.fail(cause);}
		}
		final Construction owner;
		final int first, last;
		boolean success;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	private static void require(final boolean valid, final String message){if(!valid){throw new IllegalArgumentException(message);}}
	private final int width;
	private final int[] ordinals;
	private final char[][] keys;
	private final int[][] offsets, postings;
	private final long postingCount;
	private static final char[] EMPTY_KEYS=new char[0];
	private static final int[] EMPTY_IDS=new int[0], EMPTY_OFFSETS={0};
}
