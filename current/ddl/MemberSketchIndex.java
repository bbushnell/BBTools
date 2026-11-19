package ddl;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.locks.ReadWriteLock;

import template.Accumulator;
import template.ThreadWaiter;

/** Immutable sparse DDL postings for neighbors among one family's original members.
 * Only populated nonzero keys have offsets; small families do not allocate a full
 * width-by-65536 offset matrix. Returned IDs are local input positions, with ties
 * broken by original member ordinal. Queries share this index but never scratch.
 * The neighbor limit bounds output, not posting traversal in a repetitive family.
 * @author Brian, Fischl
 */
public final class MemberSketchIndex{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Copies ordinals and sketches before construction, so callers may later reuse inputs.
	 * Ordinals must be distinct nonnegative original-member IDs; input order is unrestricted.
	 * Construction workers own disjoint bins. The caller budgets them within its allocation.
	 * Empty families are valid; failed construction throws before an index is returned.
	 */
	public MemberSketchIndex(final int width_, final int[] memberOrdinals,
			final char[][] sketches, final int buildThreads){
		require(width_>0 && width_<=DDLFamilyProfiles.MAX_BUCKETS && Integer.bitCount(width_)==1,
			"Member sketches require a supported power-of-two width");
		require(memberOrdinals!=null && sketches!=null && memberOrdinals.length==sketches.length,
			"Every member sketch needs one original ordinal");
		require(buildThreads>0, "Member index construction needs positive worker count");
		width=width_; ordinals=memberOrdinals.clone(); rows=new char[sketches.length][];
		final int[] sorted=ordinals.clone();
		Arrays.sort(sorted);
		for(int i=0; i<sorted.length; i++){
			require(sorted[i]>=0 && (i==0 || sorted[i]!=sorted[i-1]),
				"Original member ordinals must be nonnegative and unique");
			require(sketches[i]!=null && sketches[i].length==width,
				"Every member must have an exact-width sketch");
			rows[i]=sketches[i].clone();
		}
		keys=new char[width][]; offsets=new int[width][]; postings=new int[width][];
		if(rows.length==0){
			for(int bin=0; bin<width; bin++){
				keys[bin]=EMPTY_KEYS; offsets[bin]=EMPTY_OFFSETS; postings[bin]=EMPTY_IDS;
			}
		}else if(buildThreads==1){
			buildBins(0,width);
		}else{
			new Construction(Math.min(width,buildThreads)).build();
		}
		long count=0;
		for(int bin=0; bin<width; bin++){count+=postings[bin].length;}
		postingCount=count;
	}

	/** Counts only touched keys, prefixes them in unsigned order, then fills local member IDs.
	 * Scratch is reused per worker and cleared only at touched values, not across all 65536.
	 */
	private void buildBins(final int first, final int last){
		assert(first>=0 && first<last && last<=width) : "Bin owners must cover disjoint valid index ranges";
		final int[] counts=new int[65536],touched=new int[Math.min(rows.length,65535)];
		for(int bin=first; bin<last; bin++){
			int distinct=0,total=0;
			for(char[] row:rows){
				final int value=row[bin];
				if(value==0){continue;}
				if(counts[value]++==0){touched[distinct++]=value;}
				total++;
			}
			Arrays.sort(touched,0,distinct);
			final char[] binKeys=new char[distinct];
			final int[] binOffsets=new int[distinct+1],binPostings=new int[total];
			for(int i=0; i<distinct; i++){
				final int value=touched[i];
				binKeys[i]=(char)value;
				binOffsets[i+1]=binOffsets[i]+counts[value];
				counts[value]=binOffsets[i];
			}
			for(int id=0; id<rows.length; id++){
				final int value=rows[id][bin];
				if(value!=0){binPostings[counts[value]++]=id;}
			}
			for(int i=0; i<distinct; i++){
				assert(counts[touched[i]]==binOffsets[i+1]) : "Every counted key must fill exactly its own CSR slice";
				counts[touched[i]]=0;
			}
			keys[bin]=binKeys; offsets[bin]=binOffsets; postings[bin]=binPostings;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Queries            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns up to limit nonself neighbors with at least minHits matching nonzero registers.
	 * Order is descending hits, then ascending original ordinal. Output slots beyond
	 * the returned count are untouched. No zero-hit completion or alignment is performed.
	 * Member IDs in both the request and output are local input positions, not ordinals.
	 * The caller exclusively owns the workspace and both output arrays during this call;
	 * the in-use flag detects misuse but does not synchronize competing callers.
	 * Scratch scores are cleared in finally after a query begins, including early returns.
	 */
	public int topNeighbors(final int member, final int limit, final int minHits,
			final int[] outMembers, final int[] outHits, final Workspace work){
		checkMember(member);
		require(limit>=0 && minHits>0, "Neighbor limit must be nonnegative and minimum shared hits positive");
		require(outMembers!=null && outHits!=null && outMembers.length>=limit && outHits.length>=limit,
			"Neighbor output arrays must hold the requested limit");
		require(outMembers!=outHits && work!=null && !work.inUse,
			"Neighbor queries require distinct output arrays and exclusive scratch");
		work.begin(rows.length,Math.min(limit,Math.max(0,rows.length-1)));
		try{
			if(limit==0 || rows.length<2 || minHits>width){return 0;}
			final char[] query=rows[member];
			for(int bin=0; bin<width; bin++){
				if(query[bin]==0){continue;}
				final int cell=Arrays.binarySearch(keys[bin],query[bin]);
				if(cell<0){continue;}
				final int end=offsets[bin][cell+1];
				work.postingsVisited+=end-offsets[bin][cell];
				for(int p=offsets[bin][cell]; p<end; p++){
					final int id=postings[bin][p];
					if(id==member){continue;}
					if(work.hits[id]++==0){work.touched[work.touchedSize++]=id;}
				}
			}
			final int cap=Math.min(limit,rows.length-1);
			for(int i=0; i<work.touchedSize; i++){
				final int id=work.touched[i];
				if(work.hits[id]>=minHits){offer(id,cap,work);}
			}
			final int found=work.heapSize;
			for(int i=found-1; i>=0; i--){
				final int id=work.heap[0];
				outMembers[i]=id; outHits[i]=work.hits[id];
				work.heap[0]=work.heap[--work.heapSize];
				if(work.heapSize>0){siftDown(work);}
			}
			return found;
		}finally{work.end();}
	}

	/** Keeps the worst retained candidate at the root of the bounded heap. */
	private void offer(final int id, final int cap, final Workspace work){
		assert(cap>0 && work.heapSize<=cap) : "Only a positive bounded neighbor request enters heap selection";
		if(work.heapSize<cap){
			int child=work.heapSize++;
			while(child>0){
				final int parent=(child-1)>>>1;
				if(!worse(id,work.heap[parent],work.hits)){break;}
				work.heap[child]=work.heap[parent]; child=parent;
			}
			work.heap[child]=id;
		}else if(worse(work.heap[0],id,work.hits)){
			work.heap[0]=id; siftDown(work);
		}
	}

	/** Restores the worst-first heap after replacing or removing its root. */
	private void siftDown(final Workspace work){
		assert(work.heapSize>0) : "Root repair requires a nonempty retained-neighbor heap";
		final int id=work.heap[0];
		int parent=0;
		for(int child=1; child<work.heapSize; child=(parent<<1)+1){
			if(child+1<work.heapSize && worse(work.heap[child+1],work.heap[child],work.hits)){child++;}
			if(!worse(work.heap[child],id,work.hits)){break;}
			work.heap[parent]=work.heap[child]; parent=child;
		}
		work.heap[parent]=id;
	}

	/** A larger original ordinal loses a tie independently of local member order. */
	private boolean worse(final int a, final int b, final int[] hits){
		return hits[a]<hits[b] || (hits[a]==hits[b] && ordinals[a]>ordinals[b]);
	}

	/** Direct score for a selected pair; zero values never count as shared evidence. */
	public int sharedHits(final int a, final int b){
		checkMember(a); checkMember(b);
		int hits=0;
		for(int bin=0; bin<width; bin++){if(rows[a][bin]!=0 && rows[a][bin]==rows[b][bin]){hits++;}}
		return hits;
	}
	/** Number of original members, including members whose sketches contain only zeros. */
	public int size(){return rows.length;}
	/** Number of registers in every copied member sketch. */
	public int width(){return width;}
	/** Maps a checked local input position to its copied original-member ordinal. */
	public int memberOrdinal(final int localId){checkMember(localId); return ordinals[localId];}
	/** Total nonzero registers across all members; each contributes one posting entry. */
	public long postingCount(){return postingCount;}

	/** Primitive payload only, excluding object headers, build scratch and query workspaces. */
	public long primitiveBytes(){
		long bytes=4L*ordinals.length+2L*rows.length*width;
		for(int bin=0; bin<width; bin++){
			bytes+=2L*keys[bin].length+4L*offsets[bin].length+4L*postings[bin].length;
		}
		return bytes;
	}
	/** Builds a diagnostic only on failure; normal neighbor queries allocate no message strings. */
	private void checkMember(final int id){
		if(id<0 || id>=rows.length){throw new IllegalArgumentException("Local member index outside this family: "+id);}
	}
	/** Rejects invalid public inputs even when JVM assertions are disabled. */
	private static void require(final boolean valid, final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Worker State           ----------------*/
	/*--------------------------------------------------------------*/

	/** Reusable exclusive query scratch. One instance per simultaneously active query worker.
	 * Buffers grow on demand and may be reused with another index after a query finishes.
	 * Sharing a workspace concurrently is unsupported; no internal lock grants ownership.
	 */
	public static final class Workspace{
		/** Reserves buffers, resets the visit counter and marks an already-cleared workspace active. */
		private void begin(final int members, final int limit){
			assert(!inUse && touchedSize==0 && heapSize==0) : "Previous neighbor query must release and clear its scratch";
			if(hits.length<members){hits=new int[members]; touched=new int[members];}
			if(heap.length<limit){heap=new int[limit];}
			postingsVisited=0; inUse=true;
		}
		/** Clears touched scores and heap state, preserving buffers and the last visit count. */
		private void end(){
			for(int i=0; i<touchedSize; i++){hits[touched[i]]=0;}
			touchedSize=0; heapSize=0; inUse=false;
		}
		/** Posting entries visited by the last query, including its own stored entries. */
		public long postingsVisited(){return postingsVisited;}
		/** Scores by local member ID, IDs needing score cleanup, and the worst-first ID heap.
		 * Arrays start empty and keep their grown capacities between exclusive queries.
		 */
		private int[] hits=EMPTY_IDS,
			touched=EMPTY_IDS,
			heap=EMPTY_IDS;
		/** Number of populated touched and heap entries, respectively. */
		private int touchedSize, heapSize;
		/** Misuse guard, not an atomic ownership claim or synchronization primitive. */
		private boolean inUse;
		/** Visit count retained after end(), reset at the next successful begin(). */
		private long postingsVisited;
	}

	/** Joins every started bin worker and rejects a partially built index on failure. */
	private final class Construction implements Accumulator<BuildWorker>{
		/** Divides bins into disjoint half-open ranges without starting workers. */
		Construction(final int count){
			workers=new ArrayList<BuildWorker>(count);
			for(int i=0; i<count; i++){
				workers.add(new BuildWorker(this,(int)((long)width*i/count),(int)((long)width*(i+1)/count)));
			}
		}
		/** Starts and joins workers; a startup failure still joins every worker already started. */
		void build(){
			try{
				if(!ThreadWaiter.startAndWait(workers,this)){fail(new IllegalStateException("Incomplete member-index build"));}
			}catch(Throwable cause){
				fail(cause);
				final ArrayList<BuildWorker> started=new ArrayList<BuildWorker>();
				for(BuildWorker worker:workers){if(worker.getState()!=Thread.State.NEW){started.add(worker);}}
				ThreadWaiter.waitForThreadsToFinish(started);
			}
			if(failure!=null){throw new IllegalStateException("Member index failed before publication",failure);}
			if(done!=workers.size()){throw new IllegalStateException("Member index omitted a bin range");}
		}
		/** Retains the first failure reported by a worker or by the construction caller. */
		synchronized void fail(final Throwable cause){if(failure==null){failure=cause;}}
		/** Called by ThreadWaiter on the construction caller after all workers are joined. */
		@Override
		public void accumulate(final BuildWorker worker){
			if(!worker.success){fail(new IllegalStateException("Member-index worker did not finish"));}
			done++;
		}
		/** Requires successful completion of every assigned bin range. */
		@Override
		public boolean success(){return failure==null && done==workers.size();}
		/** No accumulator lock is needed because accumulation follows worker joins. */
		@Override
		public ReadWriteLock rwlock(){return null;}// Accumulation follows worker joins.
		/** Fixed construction-worker roster, retained until all started workers terminate. */
		final ArrayList<BuildWorker> workers;
		/** First failure; synchronized writers preserve it and volatile readers observe it. */
		volatile Throwable failure;
		/** Joined workers accumulated by the construction caller, including failed workers. */
		int done;
	}

	/** A worker owns all arrays for its bin range and no shared per-posting counters. */
	private final class BuildWorker extends Thread{
		/** Binds the failure recipient and exclusive half-open bin range. */
		BuildWorker(final Construction owner_, final int first_, final int last_){
			super("ddl-member-index-"+first_);
			owner=owner_;
			first=first_;
			last=last_;
		}
		/** Publishes success only after all assigned bins are built; reports any thrown failure. */
		@Override
		public void run(){
			try{buildBins(first,last); success=true;}
			catch(Throwable cause){owner.fail(cause);}
		}
		/** Coordinator that retains worker failures for the constructing caller. */
		final Construction owner;
		/** Exclusive half-open bin range [first,last). */
		final int first, last;
		/** Read by the construction caller only after this worker is joined. */
		boolean success;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Supported register count shared by every member row and every per-bin array. */
	private final int width;
	/** Copied original ordinals indexed by local input position. */
	private final int[] ordinals;
	/** Copied member rows and sorted distinct nonzero keys for each bin, respectively. */
	private final char[][] rows, keys;
	/** Per-bin key slice boundaries and the local member IDs in those slices. */
	private final int[][] offsets, postings;
	/** Number of retained nonzero member-register postings. */
	private final long postingCount;
	/** Shared immutable empty key array for an empty family. */
	private static final char[] EMPTY_KEYS=new char[0];
	/** Shared empty postings and the single terminal zero required by empty key offsets. */
	private static final int[] EMPTY_IDS=new int[0], EMPTY_OFFSETS={0};
}
