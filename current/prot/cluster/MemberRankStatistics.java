package prot.cluster;

import java.util.Arrays;

/** Incremental retained-member statistics for the length/KD ranking pivot.
 * One publisher owns this mutable accumulator; it is not thread-safe. Readers
 * use phase snapshots or an externally coordinated published prefix. Adding a
 * member costs amortized O(log N); both statistics are O(1) reads, without member rescans.
 * Lengths retain duplicate observations in primitive heaps; each member KD gets
 * equal weight regardless of its peptide length. Existing model scalars are separate.
 * @author Brian Bushnell, Keqing
 */
public final class MemberRankStatistics{

	/** Appends one measured member; undefined KD must be handled explicitly by the caller. */
	public void add(final int length, final double kd){
		require(length>0, "Retained-member ranking requires a positive peptide length");
		validateKd(kd);
		require(count<MAX_MEMBERS, "Member statistics must fit FamilyMembers' supported population bound");
		if(count==0){singleLength=length;}
		else if(count==1){
			lower=capacity(lower,1); upper=capacity(upper,1);
			lower[0]=Math.min(singleLength,length); upper[0]=Math.max(singleLength,length);
			lowerSize=upperSize=1;
		}else{
			if(length<=lower[0]){lower=offer(lower,lowerSize,length,true); lowerSize++;}
			else{upper=offer(upper,upperSize,length,false); upperSize++;}
			if(lowerSize>upperSize+1){
				final int moved=lower[0]; removeRoot(lower,--lowerSize,true);
				upper=offer(upper,upperSize,moved,false); upperSize++;
			}else if(upperSize>lowerSize){
				final int moved=upper[0]; removeRoot(upper,--upperSize,false);
				lower=offer(lower,lowerSize,moved,true); lowerSize++;
			}
		}
		// Compensated summation keeps the mean of member means, not a residue-weighted pool.
		final double adjusted=kd-compensation, next=sumKd+adjusted;
		compensation=(next-sumKd)-adjusted; sumKd=next;
		count++;
		assert(count==1 || lowerSize+upperSize==count && lowerSize>=upperSize &&
			lowerSize<=upperSize+1 && lower[0]<=upper[0]) :
			"Median heaps must partition the complete member prefix around its central order statistics";
	}

	public int count(){return count;}

	/** Even populations average the two central lengths using overflow-safe arithmetic. */
	public double medianLength(){
		require(count>0, "An empty population has no measured median; use a retained nonempty snapshot explicitly");
		return count==1 ? singleLength : lowerSize>upperSize ? lower[0] : ((long)lower[0]+upper[0])/2.0;
	}

	public double meanKd(){
		require(count>0, "An empty population has no measured KD mean; use a retained nonempty snapshot explicitly");
		final double mean=sumKd/count;
		validateKd(mean);
		return mean;
	}

	/** Freezes this population, or preserves the exact last-nonempty provenance for an empty template.
	 * The caller supplies the measurement layer; clear() never invents a replacement descriptor.
	 */
	public Snapshot freeze(final int layer, final Snapshot lastNonempty){
		require(layer>=0, "Member-statistic provenance requires a nonnegative measurement layer");
		if(count==0){
			require(lastNonempty!=null, "Empty template has no last-nonempty member descriptors; consensus fallback is forbidden");
			return lastNonempty;
		}
		return new Snapshot(count,medianLength(),meanKd(),layer);
	}

	/** The owning publisher clears only after readers have captured their frozen descriptors. */
	public void clear(){
		count=lowerSize=upperSize=singleLength=0; sumKd=compensation=0;
	}

	/** Validated immutable descriptors, also suitable for explicit seed/codec hydration. */
	public static final class Snapshot{
		public Snapshot(final int members_, final double median_, final double meanKd_, final int layer_){
			require(members_>0 && members_<=MAX_MEMBERS && layer_>=0,
				"A ranking snapshot must name a real nonempty bounded population and measurement layer");
			require(Double.isFinite(median_) && median_>=1 && median_<=Integer.MAX_VALUE && median_*2==Math.rint(median_*2) &&
				((members_&1)==0 || median_==Math.rint(median_)),
				"Member median must be a finite positive integer or half-integer within peptide length bounds");
			validateKd(meanKd_);
			members=members_; medianLength=median_; meanKd=meanKd_; layer=layer_;
		}
		public final int members,layer;
		public final double medianLength,meanKd;
	}

	/** Reused by one query worker; publication fills all fields from one coherent prefix.
	 * Do not share this scratch object between threads or retain it across another read.
	 */
	public static final class View{
		public int members,layer;
		public double medianLength,meanKd;
	}

	/** Duplicate-safe heap insertion; generic LongHeap deliberately requires distinct payloads. */
	private static int[] offer(int[] heap, final int size, final int value, final boolean maximum){
		assert(size>=0 && size<=heap.length && value>0) : "Heap insertion uses an existing member prefix and a positive length";
		heap=capacity(heap,size+1);
		int at=size;
		while(at>0){
			final int parent=(at-1)>>>1;
			if(!before(value,heap[parent],maximum)){break;}
			heap[at]=heap[parent]; at=parent;
		}
		heap[at]=value;
		return heap;
	}

	/** Removes the old root after the caller decreases its heap size by one. */
	private static void removeRoot(final int[] heap, final int size, final boolean maximum){
		assert(size>=0 && size<heap.length) : "A median rebalance removes exactly one existing root";
		if(size==0){return;}
		final int value=heap[size]; int at=0;
		for(int child=1; child<size; child=at*2+1){
			if(child+1<size && before(heap[child+1],heap[child],maximum)){child++;}
			if(!before(heap[child],value,maximum)){break;}
			heap[at]=heap[child]; at=child;
		}
		heap[at]=value;
	}
	private static boolean before(final int a, final int b, final boolean maximum){return maximum ? a>b : a<b;}
	private static int[] capacity(final int[] heap, final int required){
		assert(required>0 && required<=MAX_MEMBERS) : "Median storage is bounded by the accepted FamilyMembers population";
		return required<=heap.length ? heap : Arrays.copyOf(heap,(int)Math.min(MAX_MEMBERS,Math.max(required,Math.max(2L,heap.length*2L))));
	}
	private static void validateKd(final double value){
		require(Double.isFinite(value) && value>=-4.5 && value<=4.5,
			"Ranking KD requires a measured mean within native Kyte-Doolittle bounds [-4.5,4.5]");
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	private static final int MAX_MEMBERS=1<<30;
	private static final int[] EMPTY=new int[0];
	private int[] lower=EMPTY,upper=EMPTY;
	private int lowerSize,upperSize,count,singleLength;
	private double sumKd,compensation;
}
