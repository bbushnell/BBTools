package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

import map.LongHashSet;
import map.LongObjectMap;

/** Optional post-pruning preference against incompatible unsupported opposite-strand CDSs.
 * Native scores, support bonuses and candidate lists are never changed. Each graph starts
 * with cleared penalties, including disabled graphs. Annotation is not an input.
 * @author Sayu
 */
public final class TailCascadeOppositePenalty{
	/** Disabled production path: clear reused node state without allocating an index or report. */
	public static void clear(final ArrayList<Orf>[] frames){
		for(ArrayList<Orf> list:frames){
			for(Orf orf:list){orf.tailCascadeOppositePenalty=0;}
		}
	}

	/*--------------------------------------------------------------*/
	/* Retained-candidate prepass */
	/*--------------------------------------------------------------*/

	/** Reset the graph, select supported-stop representatives, then query interval indexes.
	 * Each unsupported candidate receives the configured penalty once, regardless of the
	 * number of conflicting representatives. All qualifying pairs remain observable.
	 */
	public static Result apply(final ArrayList<Orf>[] frames, final float penalty, final int maxOverlap){
		return apply(frames, penalty, maxOverlap, true);
	}

	/** Ordinary production callers do not retain conflict pairs when observation is disabled. */
	public static Result apply(final ArrayList<Orf>[] frames, final float penalty, final int maxOverlap,
			final boolean observe){
		require(frames!=null && frames.length==6, "Require six retained frame lists");
		require(Float.isFinite(penalty) && penalty>=0 && maxOverlap>=0, "Invalid opposite-strand penalty or overlap limit");
		final long start=System.nanoTime();
		final Result result=new Result(observe);
		String contig=null;
		final LongObjectMap<Orf> representatives=new LongObjectMap<Orf>(Orf.class);
		for(ArrayList<Orf> list:frames){
			require(list!=null, "Missing retained frame list");
			for(Orf orf:list){
				require(orf!=null, "Null retained candidate");
				orf.tailCascadeOppositePenalty=0;
				result.retainedCandidates++;
				if(penalty==0){continue;}
				if(contig==null){contig=orf.scafName;}
				require(contig.equals(orf.scafName), "A penalty graph must contain only one contig");
				if(!completeCds(orf)){continue;}
				require(Float.isFinite(orf.proteinBonus) && orf.proteinBonus>=0, "Support bonus must be finite and nonnegative");
				if(orf.proteinBonus==0){result.eligibleCandidates++; continue;}
				require(Float.isFinite(orf.tailCascadeNativeScore),
					"Supported representative needs a recorded finite native score");
				final long key=stopKey(orf);
				final Orf old=representatives.get(key);
				if(old==null || prefer(orf, old)){representatives.put(key, orf);}
			}
		}
		if(penalty>0){
			result.supportedStops=representatives.size();
			final ArrayList<Orf> plus=new ArrayList<Orf>(), minus=new ArrayList<Orf>();
			for(Orf representative:representatives.values()){
				if(representative!=null){(representative.strand==0 ? plus : minus).add(representative);}
			}
			final Index[] indexes={new Index(plus), new Index(minus)};
			for(ArrayList<Orf> list:frames){
				for(Orf candidate:list){
					if(completeCds(candidate) && candidate.proteinBonus==0){
						indexes[1-candidate.strand].visit(candidate, penalty, maxOverlap, result);
					}
				}
			}
		}
		result.prepassNanos=System.nanoTime()-start;
		return result;
	}

	/** Existing score comparator uses ascending positional order to break equal scores. */
	private static boolean prefer(final Orf candidate, final Orf old){
		return candidate.tailCascadeNativeScore>old.tailCascadeNativeScore ||
			(candidate.tailCascadeNativeScore==old.tailCascadeNativeScore && candidate.compareTo(old)<0);
	}

	/** Both biological boundaries must be complete; RNA and synthetic edges are excluded. */
	public static boolean completeCds(final Orf orf){
		return orf.type==ProkObject.CDS && !orf.startTruncated && !orf.stopTruncated;
	}

	/** The real-stop identity is strand plus the genomic leftmost base of its codon. */
	public static long stopKey(final Orf orf){
		final int left=orf.strand==0 ? orf.stop-2 : orf.start;
		require(left>=0 && (orf.strand==0 || orf.strand==1), "Malformed complete CDS stop");
		return ((long)left<<1)|orf.strand;
	}

	private static long coordinateKey(final Orf orf){
		return ((long)orf.start<<32)|((long)orf.stop<<1)|orf.strand;
	}

	/*--------------------------------------------------------------*/
	/* Balanced interval index */
	/*--------------------------------------------------------------*/

	/** A static median tree over sorted starts, with subtree maximum endpoints.
	 * Subtrees outside the query interval are pruned; only geometrically overlapping
	 * representatives reach the compatibility test. Dense overlap still requires reporting its pairs.
	 */
	private static final class Index{
		Index(final ArrayList<Orf> input){
			orfs=input.toArray(new Orf[input.size()]);
			Arrays.sort(orfs, POSITION);
			maximumEnd=new int[orfs.length];
			build(0, orfs.length);
		}

		private int build(final int from, final int to){
			if(from>=to){return -1;}
			final int mid=(from+to)>>>1;
			maximumEnd[mid]=Math.max(orfs[mid].stop, Math.max(build(from, mid), build(mid+1, to)));
			return maximumEnd[mid];
		}

		void visit(final Orf candidate, final float penalty, final int maxOverlap, final Result result){
			visit(0, orfs.length, candidate, penalty, maxOverlap, result);
		}

		private void visit(final int from, final int to, final Orf candidate, final float penalty,
				final int maxOverlap, final Result result){
			if(from>=to){return;}
			final int mid=(from+to)>>>1;
			if(maximumEnd[mid]<candidate.start || orfs[from].start>candidate.stop){return;}
			visit(from, mid, candidate, penalty, maxOverlap, result);
			final Orf supported=orfs[mid];
			if(supported.start<=candidate.stop && supported.stop>=candidate.start){
				result.overlapComparisons++;
				if(!candidate.isValidPrev(supported, maxOverlap) && !supported.isValidPrev(candidate, maxOverlap)){
					if(candidate.tailCascadeOppositePenalty==0){
						candidate.tailCascadeOppositePenalty=penalty;
						result.penalizedCount++;
						if(result.observe){result.penalized.add(candidate);}
					}
					result.pairCount++;
					if(result.observe){
						result.conflictCandidates.add(candidate);
						result.conflictSupported.add(supported);
					}
					result.affectedStops.add(stopKey(supported));
				}
			}
			visit(mid+1, to, candidate, penalty, maxOverlap, result);
		}

		final Orf[] orfs;
		final int[] maximumEnd;
	}

	private static final Comparator<Orf> POSITION=(a,b)->{
		int comparison=Integer.compare(a.start, b.start);
		if(comparison==0){comparison=Integer.compare(a.stop, b.stop);}
		if(comparison==0){comparison=Integer.compare(a.frame, b.frame);}
		return comparison;
	};

	/*--------------------------------------------------------------*/
	/* Per-caller observation */
	/*--------------------------------------------------------------*/

	/** Primitive observations are emitted after selection, with no mutable candidate exposure. */
	public interface Observer{
		void conflict(int start, int end, int strand, float nativeScore, float penalty, boolean selected,
			int supportedStart, int supportedEnd, int supportedStrand, float supportedNativeScore, boolean supportedSelected);
		void complete(int retainedCandidates, int supportedStops, int eligibleCandidates,
			int penalizedCandidates, int selectedPenalizedCandidates, int affectedSupportedStops,
			long overlapComparisons, long prepassNanos);
	}

	/** One graph's observations; counters distinguish candidate counts from conflict pairs. */
	public static final class Result{
		private Result(final boolean observe_){
			observe=observe_;
			penalized=observe ? new ArrayList<Orf>() : null;
			conflictCandidates=observe ? new ArrayList<Orf>() : null;
			conflictSupported=observe ? new ArrayList<Orf>() : null;
		}
		/** Emit only immutable primitive values after the caller has finalized its chosen path. */
		public void report(final ArrayList<Orf> selected, final Observer observer){
			if(observer==null){return;}
			require(observe, "Conflict observation was disabled for this graph");
			final LongHashSet selectedCoordinates=new LongHashSet();
			if(selected!=null){
				for(Orf orf:selected){if(orf.type==ProkObject.CDS){selectedCoordinates.add(coordinateKey(orf));}}
			}
			int selectedPenalized=0;
			for(Orf candidate:penalized){if(selectedCoordinates.contains(coordinateKey(candidate))){selectedPenalized++;}}
			for(int i=0; i<conflictCandidates.size(); i++){
				final Orf candidate=conflictCandidates.get(i), supported=conflictSupported.get(i);
				observer.conflict(candidate.start, candidate.stop, candidate.strand, candidate.tailCascadeNativeScore,
					candidate.tailCascadeOppositePenalty, selectedCoordinates.contains(coordinateKey(candidate)),
					supported.start, supported.stop, supported.strand, supported.tailCascadeNativeScore,
					selectedCoordinates.contains(coordinateKey(supported)));
			}
			observer.complete(retainedCandidates, supportedStops, eligibleCandidates, penalized.size(),
				selectedPenalized, affectedStops.size(), overlapComparisons, prepassNanos);
		}

		public int penalizedCandidates(){return penalizedCount;}
		public int affectedSupportedStops(){return affectedStops.size();}
		public int conflictPairs(){return pairCount;}
		public int supportedStops(){return supportedStops;}
		public long overlapComparisons(){return overlapComparisons;}
		public long prepassNanos(){return prepassNanos;}

		private int retainedCandidates, supportedStops, eligibleCandidates;
		private int penalizedCount, pairCount;
		private long overlapComparisons, prepassNanos;
		private final boolean observe;
		private final ArrayList<Orf> penalized, conflictCandidates, conflictSupported;
		private final LongHashSet affectedStops=new LongHashSet();
	}

	private static void require(final boolean valid, final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}
}
