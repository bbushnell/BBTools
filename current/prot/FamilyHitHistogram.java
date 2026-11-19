package prot;

import java.util.Arrays;

/** Worker-local terminal-winner and optional uncensored raw-hit counts.
 * The caller supplies one complete, unpruned, unweighted workspace per query.
 * This collector checks its shape, not the classifier/package provenance.
 * Bin6 contains six or more hits; bin0 is derived only for raw observations.
 * @author Sayu
 */
public final class FamilyHitHistogram{

	public FamilyHitHistogram(final int families_, final boolean captureRaw){
		require(families_>0 && families_<=Integer.MAX_VALUE/BINS,
			"Histogram needs a positive live family count with representable seven-bin offsets");
		families=families_;
		winners=new long[families*BINS]; raw=captureRaw ? new long[winners.length] : null;
		queryScores=new int[families]; generations=new int[families];
	}

	/** Borrows arrays only during this call. Untouched scores are never inspected.
	 * Validate the entire touched roster before changing any histogram counter.
	 * One collector belongs to one worker; concurrent calls are unsupported. */
	public void beginQuery(final int[] touched, final int[] scores, final int touchedCount){
		require(!pending,"A query is unfinished; finish it before beginning another");
		require(touched!=null && scores!=null && scores.length==families &&
			touchedCount>=0 && touchedCount<=families && touchedCount<=touched.length,
			"Completed workspace must match the live roster and touched count");
		// Every bin is bounded by the number of queries, so this also protects raw increments.
		Math.addExact(completed,1L);
		if(generation==Integer.MAX_VALUE){Arrays.fill(generations,0);generation=0;}
		generation++;
		for(int i=0;i<touchedCount;i++){
			final int family=touched[i];
			require(family>=0 && family<families,"Touched family is outside the live roster");
			require(generations[family]!=generation,"Touched families must be unique within a query");
			require(scores[family]>0,"Touched family must have a positive raw-hit count");
			generations[family]=generation;queryScores[family]=scores[family];
		}
		if(raw!=null){
			for(int i=0;i<touchedCount;i++){
				final int family=touched[i];
				raw[family*BINS+Math.min(queryScores[family],BINS-1)]++;
			}
		}
		pending=true;
	}

	/** Close the current query with one matching winner or the native -1/-1 rejection.
	 * Invalid terminal observations leave the query open and unexportable. */
	public void finishQuery(final int winnerOrdinal, final int rawHits){
		require(pending,"A terminal observation requires a begun unfinished query");
		require(winnerOrdinal>=-1 && winnerOrdinal<families,"Winner is outside the live roster");
		if(winnerOrdinal<0){require(rawHits==-1,"Explicit rejection requires winner=-1 and rawHits=-1");}
		else{
			require(rawHits>0 && generations[winnerOrdinal]==generation && queryScores[winnerOrdinal]==rawHits,
				"Winner and raw hits must match this query's completed workspace");
		}
		final long nextCompleted=Math.addExact(completed,1L);
		final long nextAccepted=winnerOrdinal<0 ? accepted : Math.addExact(accepted,1L);
		if(winnerOrdinal>=0){
			final int bin=winnerOrdinal*BINS+Math.min(rawHits,BINS-1);
			winners[bin]=Math.addExact(winners[bin],1L);
		}
		completed=nextCompleted;accepted=nextAccepted;pending=false;
	}

	/** Merge completed worker counts with the same ordinal roster and capture mode.
	 * Roster identity beyond its size remains the owning classifier's responsibility.
	 * No source arrays are retained or modified; self-merging is rejected. */
	public void mergeFrom(final FamilyHitHistogram other){
		require(other!=null && other!=this && other.families==families && (other.raw==null)==(raw==null),
			"Merge requires distinct collectors with matching roster size and raw capture mode");
		require(!pending && !other.pending,"Cannot merge an unfinished query on either collector");
		// Check denominators first: each bin and accepted total is bounded by its completed count.
		final long nextCompleted=Math.addExact(completed,other.completed);
		final long nextAccepted=Math.addExact(accepted,other.accepted);
		for(int i=0;i<winners.length;i++){winners[i]=Math.addExact(winners[i],other.winners[i]);}
		if(raw!=null){for(int i=0;i<raw.length;i++){raw[i]=Math.addExact(raw[i],other.raw[i]);}}
		completed=nextCompleted;accepted=nextAccepted;
	}

	/** Detached seven-bin rows in live-family ordinal order; rejected queries have no winner. */
	public long[] copyWinnerBins(){
		require(!pending,"Cannot export an unfinished query");
		return winners.clone();
	}

	/** Each raw row sums to completedQueries, including all-zero and rejected queries.
	 * Zero bins are filled once at export, never by a family-by-query zero loop. */
	public long[] copyRawBins(){
		require(!pending,"Cannot export an unfinished query");
		require(raw!=null,"Raw capture is disabled");
		final long[] copy=raw.clone();
		for(int family=0;family<families;family++){
			final int offset=family*BINS;long positive=0;
			for(int hit=1;hit<BINS;hit++){positive=Math.addExact(positive,copy[offset+hit]);}
			if(positive>completed){throw new IllegalStateException("Raw observations exceed completed query denominator");}
			copy[offset]=completed-positive;
		}
		return copy;
	}

	public long completedQueries(){return completed;}
	public long acceptedQueries(){return accepted;}

	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}
	public static final int BINS=7;
	private final int families;
	private final long[] winners,raw;
	private final int[] queryScores,generations;
	private long completed,accepted;
	private int generation;
	private boolean pending;
}
