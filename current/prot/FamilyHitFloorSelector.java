package prot;

import java.util.Arrays;

/** Proposes one short-family admission vector from terminal-winner counts.
 * Counts describe the previous exact classifier replay, not independent family truth.
 * Removed winners only price a proposal: competition can expose other winners, so
 * every returned vector requires a new complete replay before selection.
 * @author Brian Bushnell, Keqing
 */
public final class FamilyHitFloorSelector{

	/** Borrows read-only histograms for this call. Rows have seven bins per family:
	 * zero through five hits, then six or more. Each accepted query contributes once
	 * to its terminal winner; rejected queries contribute no winner. Support is the
	 * frozen eligibility mask from the initial fitting-positive pass, not this round.
	 * Physical ranks must be unique; ties prefer the smaller original rank.
	 */
	public static Proposal propose(final byte[] current, final int[] physicalRanks,
			final boolean[] supported, final long positiveQueries, final long[] positiveWinners,
			final long[] negativeQueries, final long[][] negativeWinners){
		return new Selection(current, physicalRanks, supported, positiveQueries,
			positiveWinners, negativeQueries, negativeWinners).run();
	}

	/** Owns its output arrays. Predicted counts are old winners surviving the proposed
	 * floors, never certified new classifier outcomes or proof of an error ceiling. */
	public static final class Proposal{
		private Proposal(final byte[] floors_, final long lost_, final long[] remaining_, final int changed_){
			floors=floors_; attributedPositiveLoss=lost_;
			remainingNegativeWinners=remaining_; changedFamilies=changed_;
		}
		public final byte[] floors;
		public final long attributedPositiveLoss;
		public final long[] remainingNegativeWinners;
		public final int changedFamilies;
	}

	private static final class Selection{
		Selection(final byte[] current, final int[] physicalRanks, final boolean[] supported,
				final long positiveQueries, final long[] positives,
				final long[] negativeQueries, final long[][] negatives){
			require(current!=null && current.length>0 && current.length<=Integer.MAX_VALUE/BINS,
				"Floor selection needs a nonempty roster whose seven-bin offsets fit int");
			final int n=current.length;
			require(physicalRanks!=null && physicalRanks.length==n && supported!=null && supported.length==n,
				"Every floor needs its physical family rank and initial-positive support state");
			require(negativeQueries!=null && negativeQueries.length>0 && negativeQueries.length<=4 &&
				negatives!=null && negatives.length==negativeQueries.length,
				"Floor selection needs one denominator and histogram for each of at most four negative classes");
			original=current; floors=current.clone(); ranks=physicalRanks; eligible=supported;
			positive=positives; negative=negatives;
			final int[] ordered=ranks.clone(); Arrays.sort(ordered);
			for(int f=0; f<n; f++){
				require(floors[f]>=1 && floors[f]<=MAX_FLOOR,
					"Released floor range is1..6; ordinal="+f+", floor="+floors[f]);
				require(ordered[f]>=0 && (f==0 || ordered[f]>ordered[f-1]),
					"Physical ranks must be nonnegative and unique for deterministic proposal ties");
			}
			validateHistogram(positive, positiveQueries, "positive");
			remaining=new long[negative.length]; allowance=new long[negative.length];
			active=new boolean[negative.length];
			for(int p=0; p<negative.length; p++){
				remaining[p]=validateHistogram(negative[p], negativeQueries[p], "negative class "+p);
				// Integer equivalent of accepts*1,000,000<=queries*100, without overflow.
				allowance[p]=negativeQueries[p]/10000;
			}
			heap=new int[n]; boundary=new int[n]; loss=new long[n]; removal=new long[n];
		}

		Proposal run(){
			assert(size==0) : "A selector is single-use; its heap must start empty";
			for(int f=0; f<floors.length; f++){
				if(!eligible[f] && floors[f]<MAX_FLOOR){apply(f, MAX_FLOOR);}
			}
			refreshClasses(); rebuild();
			while(size>0){
				final int family=pop();
				apply(family, boundary[family]);
				// A completed class no longer prices proposals. Rebuild at most once per class.
				if(refreshClasses()){rebuild();}
				else if(price(family)){push(family);}
			}
			int changed=0;
			for(int f=0; f<floors.length; f++){
				require(floors[f]>=original[f] && floors[f]<=MAX_FLOOR && (eligible[f] || floors[f]==MAX_FLOOR),
					"A proposal may only raise floors and must retain the unsupported-family baseline");
				if(floors[f]!=original[f]){changed++;}
			}
			return new Proposal(floors, positiveLoss, remaining, changed);
		}

		/** Returns the accepted count and rejects censored or double-counted winner histograms. */
		private long validateHistogram(final long[] counts, final long queries, final String label){
			require(counts!=null && counts.length==floors.length*BINS && queries>0,
				label+" requires seven bins per family and a positive query denominator");
			long total=0;
			for(int f=0; f<floors.length; f++){
				for(int hit=0; hit<BINS; hit++){
					final long count=counts[f*BINS+hit];
					require(count>=0 && (hit>=floors[f] || count==0),
						label+" has a negative count or a winner below its measured floor; family="+f+", hit="+hit);
					require(count<=queries-total,
						label+" winner counts exceed query denominator; each accepted query must appear once");
					total+=count;
				}
			}
			return total;
		}

		/** Price the next observed negative-winning boundary in currently over-limit classes. */
		private boolean price(final int family){
			if(!eligible[family] || floors[family]==MAX_FLOOR){return false;}
			final int offset=family*BINS, floor=floors[family];
			int next=0;
			for(int hit=floor; hit<MAX_FLOOR && next==0; hit++){
				for(int p=0; p<negative.length; p++){
					if(active[p] && negative[p][offset+hit]>0){next=hit+1; break;}
				}
			}
			if(next==0){return false;}
			long pos=0, neg=0;
			for(int hit=floor; hit<next; hit++){
				pos=Math.addExact(pos, positive[offset+hit]);
				for(int p=0; p<negative.length; p++){
					if(active[p]){neg=Math.addExact(neg, negative[p][offset+hit]);}
				}
			}
			assert(neg>0) : "A next negative-winning boundary must remove at least one priced winner";
			boundary[family]=next; loss[family]=pos; removal[family]=neg;
			return true;
		}

		private void apply(final int family, final int next){
			final int floor=floors[family], offset=family*BINS;
			require(next>floor && next<=MAX_FLOOR,
				"Each proposal step must raise a floor within1..6; family="+family+", old="+floor+", next="+next);
			for(int hit=floor; hit<next; hit++){
				positiveLoss=Math.addExact(positiveLoss, positive[offset+hit]);
				for(int p=0; p<negative.length; p++){
					remaining[p]-=negative[p][offset+hit];
					assert(remaining[p]>=0) : "A terminal winner cannot be removed twice across disjoint floor intervals";
				}
			}
			floors[family]=(byte)next;
		}

		private boolean refreshClasses(){
			boolean changed=false;
			for(int p=0; p<active.length; p++){
				final boolean over=remaining[p]>allowance[p];
				changed|=active[p]!=over; active[p]=over;
			}
			return changed;
		}

		/** Primitive family-ordinal heap; no candidate object per family or per raise. */
		private void rebuild(){
			size=0;
			for(int f=0; f<floors.length; f++){if(price(f)){heap[size++]=f;}}
			for(int i=size/2-1; i>=0; i--){down(i);}
		}
		private boolean before(final int a, final int b){
			// Fail loudly on unsupported population sizes rather than round a price tie.
			final long left=Math.multiplyExact(loss[a], removal[b]);
			final long right=Math.multiplyExact(loss[b], removal[a]);
			return left<right || left==right && ranks[a]<ranks[b];
		}
		private int pop(){
			assert(size>0) : "Only a nonempty proposal heap has a next family";
			final int result=heap[0]; heap[0]=heap[--size];
			if(size>0){down(0);}
			return result;
		}
		private void push(final int family){
			assert(size<heap.length) : "A family is reinserted only after removing its preceding proposal";
			int i=size++;
			while(i>0){
				final int parent=(i-1)>>>1;
				if(!before(family, heap[parent])){break;}
				heap[i]=heap[parent]; i=parent;
			}
			heap[i]=family;
		}
		private void down(int i){
			assert(i>=0 && i<size) : "Heap repair needs an occupied family-ordinal slot";
			final int family=heap[i];
			while(i<size/2){
				int child=i*2+1;
				if(child+1<size && before(heap[child+1], heap[child])){child++;}
				if(!before(heap[child], family)){break;}
				heap[i]=heap[child]; i=child;
			}
			heap[i]=family;
		}

		final byte[] original, floors;
		final int[] ranks, heap, boundary;
		final boolean[] eligible, active;
		final long[] positive, remaining, allowance, loss, removal;
		final long[][] negative;
		int size;
		long positiveLoss;
	}

	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}
	public static final int BINS=7, MAX_FLOOR=6;
}
