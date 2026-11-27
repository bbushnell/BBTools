package ddl;

import java.util.Arrays;

import cardinality.DynamicDemiLog;

/** Immutable inverted lookup for a small collection of native DDL objects.
 * Legacy key-array inputs remain supported. Native queries borrow maxArray
 * without the copy made by DDLIndexBase.query(DynamicDemiLog).
 * Each nonempty reference bin contributes one packed key16/reference16 posting.
 * Per-bin postings are sorted, so lookup finds an equal-key run by binary search.
 * This avoids the value-wide offset tables used by DDLIndexBase implementations;
 * those indexes remain available for very large reference collections. The actual
 * crossover depends on the data and must be measured, not inferred from the name.
 * This class supports at most65536references and65536bins because of its packing.
 *
 * Constructor inputs are copied into postings, not retained. Queries borrow their
 * inputs read-only. Returned or supplied counters belong to the caller; separate
 * workers may query one index concurrently with separate scratch. Query methods
 * replace counts; accumulateCounts methods add and return postings visited, as in
 * DDLIndexBase. There are no record-loading, ANI, hashing or statistics side effects.
 * Empty/empty bins never match. All references must share one key/bin geometry.
 * @author Brian Bushnell, Raiden
 */
public final class DDLIndexSmall{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Copies one key array per reference, preserving its original ordinal. */
	public DDLIndexSmall(final char[][] references){
		this(references, null);
	}

	/** Copies native reference winners into the existing postings. Object order
	 * defines reference ordinals; later object reuse cannot alter this index. */
	public DDLIndexSmall(final DynamicDemiLog[] references){
		this(null, references);
	}

	private DDLIndexSmall(final char[][] arrays, final DynamicDemiLog[] objects){
		require((arrays==null)!=(objects==null), "Supply one complete reference array, either native DDL objects or legacy keys");
		numClades=objects==null ? arrays.length : objects.length;
		require(numClades>0 && numClades<=65536,
			"Small DDL postings require one to65536references whose ordinals fit sixteen bits");
		final char[] first=reference(arrays, objects, 0);
		require(first!=null, "Every indexed reference must supply its complete nonnull key array");
		buckets=first.length;
		require(buckets>0 && buckets<=65536, "Packed DDL query bins require a width of one to65536");
		final int[] sizes=new int[buckets];
		for(int id=0; id<numClades; id++){
			final char[] reference=reference(arrays, objects, id);
			require(reference!=null && reference.length==buckets, "All references must share the same complete DDL bin geometry");
			for(int bin=0; bin<buckets; bin++){
				if(reference[bin]!=0){sizes[bin]++;}
			}
		}
		postings=new int[buckets][];
		int populated=0;
		long entries=0;
		for(int bin=0; bin<buckets; bin++){
			postings[bin]=sizes[bin]==0 ? EMPTY : new int[sizes[bin]];
			if(sizes[bin]>0){populated++;}
			entries+=sizes[bin];
		}
		populatedBuckets=populated;
		postingEntries=entries;
		Arrays.fill(sizes, 0);
		for(int id=0; id<numClades; id++){
			final char[] reference=reference(arrays, objects, id);
			for(int bin=0; bin<buckets; bin++){
				final int key=reference[bin];
				if(key!=0){postings[bin][sizes[bin]++]=(key<<16)|id;}
			}
		}
		long cells=0;
		for(final int[] bin:postings){
			Arrays.sort(bin);
			int previous=0;
			for(final int posting:bin){
				final int key=posting>>>16;
				if(key!=previous){cells++;}
				previous=key;
			}
		}
		populatedCells=cells;
		assert(cells<=entries && cells>=populatedBuckets)
			: "Every populated bin/key cell has at least one posting, and every nonempty bin has a cell";
	}

	/** Borrows the original representation; no intermediate reference table is built. */
	private static char[] reference(final char[][] arrays, final DynamicDemiLog[] objects, final int id){
		if(objects==null){return arrays[id];}
		require(objects[id]!=null, "Every indexed native DDL object must be nonnull: ordinal="+id);
		return objects[id].maxArray();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Queries            ----------------*/
	/*--------------------------------------------------------------*/

	/** Whether an existing sorted posting run contains this nonempty key.
	 * This is the non-null posting-list test for this index representation;
	 * no result list, extra index or reference-key bitset is allocated.
	 */
	public boolean hasPosting(final int bin, final int key){
		assert(bin>=0 && bin<buckets && key>=0 && key<=65535)
			: shared.KillSwitch.assertDie("Posting membership requires a bound bin and unsigned16-bit key: bin="+bin+", key="+key);
		if(key==0){return false;}
		final int[] entries=postings[bin];
		final int target=key<<16;
		int low=0;
		int high=entries.length;
		while(low<high){
			final int mid=(low+high)>>>1;
			if(entries[mid]<target){low=mid+1;}
			else{high=mid;}
		}
		return low<entries.length && (entries[low]>>>16)==key;
	}

	/** Same native-object query and fresh result semantics as DDLIndexBase. */
	public int[] query(final DynamicDemiLog query){return query(objectKeys(query));}

	/** Native-object query with caller-owned replacement counters. */
	public int[] query(final DynamicDemiLog query, final int[] counts){return query(objectKeys(query), counts);}

	/** Replaces raw hits and sums of L/(4+L), where L is the exact bin/key
	 * posting run, never the entire bin. Optional lengths retain that run size
	 * for every query bin, including explicit zero for absent/empty keys. */
	public void querySpecificity(final DynamicDemiLog query, final int[] counts,
			final double[] sums, final int[] lengths){
		final char[] keys=objectKeys(query);
		check(keys, counts);
		require(sums!=null && sums.length==numClades && (lengths==null || lengths.length==buckets),
			"Specificity scratch must bind every reference and optional query bin");
		Arrays.fill(counts, 0);
		Arrays.fill(sums, 0);
		if(lengths!=null){Arrays.fill(lengths, 0);}
		for(int bin=0; bin<buckets; bin++){
			final int key=keys[bin];
			if(key==0){continue;}
			final int[] entries=postings[bin];
			final int target=key<<16;
			int low=0, high=entries.length;
			while(low<high){
				final int mid=(low+high)>>>1;
				if(entries[mid]<target){low=mid+1;}else{high=mid;}
			}
			final int start=low;
			while(low<entries.length && (entries[low]>>>16)==key){low++;}
			final int length=low-start;
			if(lengths!=null){lengths[bin]=length;}
			final double contribution=length/(4.0+length);
			for(int i=start; i<low; i++){
				final int id=entries[i]&65535;
				counts[id]++;
				sums[id]+=contribution;
			}
		}
	}

	/** Queries a frozen index, then REPLACES scores for explicitly changed rows.
	 * Current rows and IDs are caller-owned; the index remains immutable. A
	 * serial owner can publish exact greedy updates between occasional rebuilds.
	 * Duplicate changed IDs repeat work but never count a row twice. */
	public int[] queryReplacing(final DynamicDemiLog query, final int[] counts,
			final DynamicDemiLog[] current, final int[] changed, final int changedCount){
		require(current!=null && current.length==numClades && changed!=null && changed!=counts && changedCount>=0 && changedCount<=changed.length,
			"Replacement scoring needs the complete current roster and a bounded changed-ID prefix");
		final char[] keys=objectKeys(query);
		query(keys, counts);
		for(int i=0; i<changedCount; i++){
			final int id=changed[i];
			require(id>=0 && id<numClades && current[id]!=null && current[id].buckets==buckets,
				"A changed DDL row must retain its stable reference ID and matching geometry");
			final char[] reference=current[id].maxArray();
			int equal=0;
			for(int bin=0; bin<buckets; bin++){if(keys[bin]!=0 && keys[bin]==reference[bin]){equal++;}}
			counts[id]=equal;
			assert(counts[id]>=0 && counts[id]<=buckets)
				: "Replacing a frozen score must retain at most one matching observation per nonempty DDL bin";
		}
		return counts;
	}

	/** Adds native-object matches, following the existing query-first/scratch-second order. */
	public long accumulateCounts(final DynamicDemiLog query, final int[] counts){return accumulateCounts(objectKeys(query), counts);}

	private static char[] objectKeys(final DynamicDemiLog query){
		require(query!=null, "A native DDL query object must be nonnull");
		return query.maxArray();
	}

	/** Same input/output order and fresh-count semantics as DDLIndexBase.query. */
	public int[] query(final char[] keys){
		final int[] counts=new int[numClades];
		accumulateCounts(keys, counts);
		return counts;
	}

	/** Replaces all supplied counters without allocating a result array. */
	public int[] query(final char[] keys, final int[] counts){
		check(keys, counts);
		Arrays.fill(counts, 0);
		accumulateCounts(keys, counts);
		return counts;
	}

	/** Adds one match per equal nonempty bin; returns postings visited. Unlike the
	 * protected DDLIndexBase primitive, this scratch interface is public. */
	public long accumulateCounts(final char[] keys, final int[] counts){
		check(keys, counts);
		long work=0;
		for(int bin=0; bin<buckets; bin++){
			if(keys[bin]!=0){work+=addMatches(bin, keys[bin], counts);}
		}
		return work;
	}

	/** Fresh counts for a unique packed bin16/key16 prefix, rather than a DDL.
	 * Multiple distinct keys in one bin are allowed. The caller deduplicates exact
	 * bin/key pairs; a duplicate would count twice. Only the prefix is inspected. */
	public int[] query(final int[] packedKeys, final int count){
		final int[] counts=new int[numClades];
		accumulateCounts(packedKeys, count, counts);
		return counts;
	}

	/** Replaces supplied counts from the same deduplicated packed-prefix contract. */
	public int[] query(final int[] packedKeys, final int count, final int[] counts){
		check(packedKeys, count, counts);
		Arrays.fill(counts, 0);
		accumulateCounts(packedKeys, count, counts);
		return counts;
	}

	/** Adds matches from unique live bin/key pairs; unused array capacity is ignored. */
	public long accumulateCounts(final int[] packedKeys, final int count, final int[] counts){
		check(packedKeys, count, counts);
		long work=0;
		for(int i=0; i<count; i++){
			final int packed=packedKeys[i];
			checkKey(packed);
			work+=addMatches(packed>>>16, packed&65535, counts);
		}
		return work;
	}

	/** Replaces separate counts for contiguous bin partitions in one lookup pass.
	 * ends contains strictly increasing exclusive bin ends, finishing at buckets.
	 * Each counts row has one counter per reference; rows and input arrays must
	 * not alias. These are query scratch, not additional per-reference metadata. */
	public int[][] query(final int[] packedKeys, final int count, final int[] ends, final int[][] counts){
		checkPartitions(packedKeys, count, ends, counts);
		for(final int[] row:counts){Arrays.fill(row, 0);}
		accumulatePartitioned(packedKeys, count, ends, counts);
		return counts;
	}

	/** Additive partitioned variant; follows the same exclusive ends and prefix. */
	public long accumulateCounts(final int[] packedKeys, final int count, final int[] ends, final int[][] counts){
		checkPartitions(packedKeys, count, ends, counts);
		return accumulatePartitioned(packedKeys, count, ends, counts);
	}

	private long accumulatePartitioned(final int[] keys, final int count, final int[] ends, final int[][] counts){
		assert(ends.length==counts.length) : "Each contiguous bin partition owns one complete counter row";
		long work=0;
		for(int i=0; i<count; i++){
			final int packed=keys[i], bin=packed>>>16;
			checkKey(packed);
			int part=0;
			while(bin>=ends[part]){part++;}
			work+=addMatches(bin, packed&65535, counts[part]);
		}
		return work;
	}

	/** Signed sorting and lower bounds agree even for unsigned keys32768..65535. */
	private long addMatches(final int bin, final int key, final int[] counts){
		assert(bin>=0 && bin<buckets && key>0 && key<=65535 && counts.length==numClades)
			: "Posting lookup requires a valid bin/nonempty key and counters in original reference order";
		final int[] entries=postings[bin];
		final int target=key<<16;
		int low=0, high=entries.length;
		while(low<high){
			final int mid=(low+high)>>>1;
			if(entries[mid]<target){low=mid+1;}
			else{high=mid;}
		}
		final int start=low;
		for(; low<entries.length && (entries[low]>>>16)==key; low++){
			final int id=entries[low]&65535;
			assert(counts[id]>=0 && counts[id]<Integer.MAX_VALUE)
				: "Additive nonnegative match counters must not overflow their int storage: reference="+id;
			counts[id]++;
		}
		return low-start;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Validation          ----------------*/
	/*--------------------------------------------------------------*/

	private void check(final char[] keys, final int[] counts){
		require(keys!=null && keys.length==buckets && counts!=null && counts.length==numClades,
			"A DDL query and its counters must retain the bound bin/reference dimensions");
	}

	private void check(final int[] keys, final int count, final int[] counts){
		require(keys!=null && count>=0 && count<=keys.length && counts!=null && counts.length==numClades && keys!=counts,
			"Packed queries require a valid live prefix and separate complete reference counters");
	}

	private void checkKey(final int packed){
		require((packed>>>16)<buckets && (packed&65535)!=0,
			"Every live packed key must identify an existing bin and a nonzero DDL key");
	}

	private void checkPartitions(final int[] keys, final int count, final int[] ends, final int[][] counts){
		require(keys!=null && count>=0 && count<=keys.length && ends!=null && ends.length>0
			&& counts!=null && counts.length==ends.length, "Partitioned queries need an exact live prefix and one output per bin interval");
		int previous=0;
		for(int part=0; part<ends.length; part++){
			require(ends[part]>previous && ends[part]<=buckets, "Partition ends must increase strictly within the bound DDL width");
			require(counts[part]!=null && counts[part].length==numClades && counts[part]!=keys && counts[part]!=ends,
				"Partition counters must be complete, separate from query keys and bin ends");
			for(int other=0; other<part; other++){
				require(counts[other]!=counts[part], "Distinct bin partitions must not share mutable counters");
			}
			previous=ends[part];
		}
		require(previous==buckets, "Partition ends must cover every DDL bin exactly once");
	}

	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Meta             ----------------*/
	/*--------------------------------------------------------------*/

	public int numClades(){return numClades;}
	public int buckets(){return buckets;}
	/** Distinct nonempty bin/key cells, as in DDLIndexBase.populatedCells. */
	public long populatedCells(){return populatedCells;}
	/** Bins with at least one posting; distinct keys within a bin count only once. */
	public int populatedBuckets(){return populatedBuckets;}
	public int emptyBuckets(){return buckets-populatedBuckets;}
	/** One outer reference per bin, including shared empty-array references. */
	public int pointerSlots(){return buckets;}
	/** Packed integers only; excludes array headers/padding and outer references. */
	public long postingPayloadBytes(){return 4L*postingEntries;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int numClades, buckets, populatedBuckets;
	private final long postingEntries, populatedCells;
	private final int[][] postings;
	private static final int[] EMPTY=new int[0];
}
