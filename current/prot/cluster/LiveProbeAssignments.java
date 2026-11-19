package prot.cluster;

import java.util.Arrays;

/** Accepted assignments from one independent local probe. One worker owns each
 * collector; it is not thread-safe. Rejected rows consume no storage. Snapshot
 * creation seals the collector, sorts by original input ordinal, and rejects
 * duplicate accepted ordinals within or across the two origin groups.
 * @author Brian Bushnell, Sayu
 */
public final class LiveProbeAssignments implements LocalProbeEvaluator.AssignmentSink{

	/** The second family ID is -1 for a single-model probe. Family IDs need not
	 * be dense or ordered; original ordinals are bounded by sequenceCount.
	 */
	public LiveProbeAssignments(final int sequenceCount,final int firstFamilyId,final int secondFamilyId){
		if(sequenceCount<0){throw new IllegalArgumentException("Negative sequence count");}
		if(firstFamilyId<0 || secondFamilyId< -1 || firstFamilyId==secondFamilyId){
			throw new IllegalArgumentException("Expected one or two distinct nonnegative family IDs");
		}
		this.sequenceCount=sequenceCount;
		this.firstFamilyId=firstFamilyId;
		this.secondFamilyId=secondFamilyId;
	}

	/** Winner -1 is an omitted rejection; smaller values are invalid.
	 * Ordinal validity and one-shot state are checked even for rejected rows.
	 * Duplicate accepted ordinals are detected by snapshot(), not on insertion.
	 */
	@Override
	public void accept(final int originalOrdinal,final int winningFamilyId,final boolean fromPool){
		requireOpen();
		if(originalOrdinal<0 || originalOrdinal>=sequenceCount){
			throw new IllegalArgumentException("Original tail ordinal is outside the frozen input: "+originalOrdinal);
		}
		if(winningFamilyId< -1){throw new IllegalArgumentException("Invalid winning family ID: "+winningFamilyId);}
		if(winningFamilyId==-1){return;}
		if(winningFamilyId!=firstFamilyId && winningFamilyId!=secondFamilyId){
			throw new IllegalArgumentException("Winner is not a replacement family: "+winningFamilyId);
		}
		final long pair=((long)originalOrdinal<<32)|(winningFamilyId&0xffffffffL);
		if(fromPool){
			if(rescuedSize==rescued.length){rescued=grow(rescued);}
			rescued[rescuedSize++]=pair;
		}else{
			if(retainedSize==retained.length){retained=grow(retained);}
			retained[retainedSize++]=pair;
		}
	}

	/** Seals this collector even if duplicate validation fails. Each returned
	 * group is sorted independently by original ordinal; winner IDs stay paired.
	 * No backing array is exposed or retained by a mutable collector.
	 */
	public Snapshot snapshot(){
		requireOpen();
		sealed=true;
		final long[] pool=Arrays.copyOf(rescued,rescuedSize);
		final long[] affected=Arrays.copyOf(retained,retainedSize);
		rescued=EMPTY; retained=EMPTY;
		Arrays.sort(pool); Arrays.sort(affected);
		checkUnique(pool); checkUnique(affected);
		for(int p=0,a=0;p<pool.length && a<affected.length;){
			final int po=ordinal(pool[p]),ao=ordinal(affected[a]);
			if(po==ao){throw duplicate(po);}
			if(po<ao){p++;}else{a++;}
		}
		return new Snapshot(pool,affected);
	}

	private void requireOpen(){
		if(sealed){throw new IllegalStateException("Assignment collector is already sealed");}
	}

	private static long[] grow(final long[] values){
		if(values.length>=MAX_ARRAY_SIZE){throw new IllegalStateException("Too many accepted assignments");}
		final int capacity=(int)Math.min(MAX_ARRAY_SIZE,Math.max(8L,(long)values.length+(values.length>>1)+1L));
		return Arrays.copyOf(values,capacity);
	}

	private static void checkUnique(final long[] values){
		for(int i=1;i<values.length;i++){
			if(ordinal(values[i-1])==ordinal(values[i])){throw duplicate(ordinal(values[i]));}
		}
	}

	private static IllegalArgumentException duplicate(final int ordinal){
		return new IllegalArgumentException("Duplicate accepted tail ordinal: "+ordinal);
	}

	private static int ordinal(final long pair){return (int)(pair>>>32);}
	private static int winner(final long pair){return (int)pair;}

	/** Immutable accepted-only evidence; absent rows must not be interpreted as
	 * automatic retention. No publication, ownership change or gain calculation
	 * occurs here. Snapshot arrays are private and exclusively owned.
	 */
	public static final class Snapshot{
		private Snapshot(final long[] rescued,final long[] retained){
			this.rescued=rescued; this.retained=retained;
		}
		public int rescuedCount(){return rescued.length;}
		public int retainedCount(){return retained.length;}
		public int rescuedOrdinal(final int index){return ordinal(rescued[index]);}
		public int rescuedWinner(final int index){return winner(rescued[index]);}
		public int retainedOrdinal(final int index){return ordinal(retained[index]);}
		public int retainedWinner(final int index){return winner(retained[index]);}

		/** Defensive sorted original-ordinal copy for proposal overlap selection. */
		public int[] rescuedOrdinals(){
			final int[] ordinals=new int[rescued.length];
			for(int i=0;i<ordinals.length;i++){ordinals[i]=ordinal(rescued[i]);}
			return ordinals;
		}

		private final long[] rescued,retained;
	}

	private static final long[] EMPTY=new long[0];
	private static final int MAX_ARRAY_SIZE=Integer.MAX_VALUE-8;
	private final int sequenceCount,firstFamilyId,secondFamilyId;
	private long[] rescued=EMPTY,retained=EMPTY;
	private int rescuedSize,retainedSize;
	private boolean sealed;
}
