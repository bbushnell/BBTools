package prot.cluster;

import java.util.Arrays;
import java.util.BitSet;

import ddl.FamilySketchBuilder;
import ddl.TailPoolIndex;
import prot.ProteinSequence;
import structures.IntList;

/** Frozen residual-pool evaluation against only the replacements of one local move.
 * The caller binds the index and sketch builder to the same input/geometry recipe.
 * Encoded sequences and replacement models are borrowed read-only. No baseline
 * membership is cleared or rewritten: temporary unassignment is represented by
 * the explicit affected ordinals, which must be unique and outside the pool.
 * Independent workers share the baseline and index, but own all query scratch.
 * @author Brian Bushnell, Keqing
 */
public final class LocalProbeEvaluator{

	public LocalProbeEvaluator(final ProteinSequence[] sequences_,final TailPoolIndex pool_,
			final FamilySketchBuilder sketches_,final TailAlignmentGate gate_,
			final CandidateFilter filter_,final int minimumLength_){
		require(sequences_!=null && pool_!=null && sketches_!=null && gate_!=null && filter_!=null,
			"Local evaluation requires frozen sequences, pool, sketch recipe and native gates");
		require(pool_.width()==sketches_.width() && minimumLength_>=gate_.minSpan,
			"Pool and query geometry must agree, with input length at least the alignment span floor");
		filter_.validateWidth(pool_.width());
		sequences=sequences_; pool=pool_; sketches=sketches_; gate=gate_; filter=filter_;
		minimumLength=minimumLength_; poolOrdinals=new BitSet(sequences.length);
		for(int slot=0; slot<pool.size(); slot++){
			final int ordinal=pool.originalOrdinal(slot);
			checkOrdinal(ordinal);
			require(!poolOrdinals.get(ordinal),"Residual pool cannot repeat an original tail ordinal");
			poolOrdinals.set(ordinal);
		}
	}

	/** Each active thread needs its own worker; no pool-sized scratch is shared. */
	public Worker newWorker(){return new Worker();}
	public int poolSize(){return pool.size();}

	private void checkOrdinal(final int ordinal){
		requireOrdinal(ordinal>=0 && ordinal<sequences.length,"Original tail ordinal is outside the frozen input",ordinal);
		requireOrdinal(sequences[ordinal]!=null && sequences[ordinal].enc!=null &&
			sequences[ordinal].enc.length>=minimumLength,
			"Residual and affected tails must pass the baseline length filter",ordinal);
	}

	/** Optional diagnostic stream for visited tails. Pool tails absent from both hit
	 * lists are known rejects and omitted. Affected tails are always visited once.
	 * A negative winner means no replacement passed; sequence identity is ordinal,
	 * so equal FASTA names or equal residues never collapse different input rows.
	 */
	public interface AssignmentSink{
		void accept(int originalOrdinal,int winningFamilyId,boolean fromPool);
	}

	public final class Worker{
		private Worker(){
			querySketches=sketches.fork(); querySketch=new char[pool.width()];
			scorer=new LocalProbeScorer(gate,filter,pool.width(),minimumLength);
		}

		/** Evaluates a fusion (second=null) or fission (two distinct replacement IDs).
		 * All affected members are rechecked; none receive automatic credit. Results
		 * count winners, not passing pairs, and net gain subtracts affected.length.
		 * Discarding this result restores the baseline because evaluation never
		 * changes it. A failed diagnostic callback aborts the probe, not the baseline.
		 */
		public Result evaluate(final Family first,final Family second,final int[] affected,
				final AssignmentSink sink){
			require(!inUse,"Local probe workers must not be shared concurrently or reentered");
			inUse=true;
			try{
				final int firstFloor=scorer.requiredHits(first);
				final int secondFloor=second==null ? 0 : scorer.requiredHits(second);
				require(second==null || first.familyId!=second.familyId,
					"Two replacement models need distinct stable IDs for native tie breaking");
				scorer.beginProbe(first,second);
				require(affected!=null,"Local probes need the actual affected baseline membership");
				final int[] members=affected.clone(); Arrays.sort(members);
				for(int i=0; i<members.length; i++){
					checkOrdinal(members[i]);
					requireOrdinal(!poolOrdinals.get(members[i]),"An affected member cannot already belong to the residual pool",members[i]);
					requireOrdinal(i==0 || members[i]!=members[i-1],"Affected family memberships overlap at original ordinal",members[i]);
				}
				final Counts counts=new Counts();
				pool.collectHits(first.sketchUnsafe(),firstFloor,firstSlots,firstHits,indexWork);
				counts.postings=indexWork.postingsVisited();
				secondSlots.clear(); secondHits.clear();
				if(second!=null){
					pool.collectHits(second.sketchUnsafe(),secondFloor,secondSlots,secondHits,indexWork);
					counts.postings=Math.addExact(counts.postings,indexWork.postingsVisited());
				}
				// Sorted slot lists admit each original pool row once, even if both
				// children retrieve it. Absence means below the model-specific floor.
				int a=0,b=0;
				while(a<firstSlots.size || b<secondSlots.size){
					final int left=a<firstSlots.size ? firstSlots.array[a] : Integer.MAX_VALUE;
					final int right=b<secondSlots.size ? secondSlots.array[b] : Integer.MAX_VALUE;
					final int slot=Math.min(left,right),ordinal=pool.originalOrdinal(slot);
					prepare(ordinal);
					if(left==slot){scorer.consider(first,firstHits.array[a++]);}
					if(right==slot){scorer.consider(second,secondHits.array[b++]);}
					counts.poolVisited++;
					if(record(ordinal,true,first,sink,counts)){counts.poolAccepted++;}
				}
				for(int ordinal:members){
					prepare(ordinal);
					querySketches.sketchEncoded(sequences[ordinal].enc,querySketch);
					scorer.consider(first,shared(querySketch,first.sketchUnsafe()));
					if(second!=null){scorer.consider(second,shared(querySketch,second.sketchUnsafe()));}
					if(record(ordinal,false,first,sink,counts)){counts.affectedAccepted++;}
				}
				return new Result(pool.size(),members.length,first.familyId,
					second==null ? -1 : second.familyId,counts);
			}finally{scorer.endProbe(); inUse=false;}
		}

		private void prepare(final int ordinal){
			if(!scorer.prepare(sequences[ordinal].enc)){
				throw new IllegalStateException("Frozen eligible tail changed length during evaluation: "+ordinal);
			}
		}
		private boolean record(final int ordinal,final boolean fromPool,final Family first,
				final AssignmentSink sink,final Counts counts){
			counts.candidates+=scorer.candidates(); counts.alignments+=scorer.alignments();
			counts.lowHits+=scorer.lowHitRejections(); counts.hydro+=scorer.hydroRejections();
			final int winner=scorer.winner();
			if(winner>=0){if(winner==first.familyId){counts.first++;}else{counts.second++;}}
			if(sink!=null){sink.accept(ordinal,winner,fromPool);}
			return winner>=0;
		}
		private final FamilySketchBuilder querySketches;
		private final char[] querySketch;
		private final LocalProbeScorer scorer;
		private final TailPoolIndex.Workspace indexWork=new TailPoolIndex.Workspace();
		private final IntList firstSlots=new IntList(),firstHits=new IntList();
		private final IntList secondSlots=new IntList(),secondHits=new IntList();
		private boolean inUse;
	}

	private static int shared(final char[] a,final char[] b){
		assert(a.length==b.length) : "Exact hit counts require the previously validated frozen width";
		int hits=0;
		for(int i=0; i<a.length; i++){if(a[i]!=0 && a[i]==b[i]){hits++;}}
		return hits;
	}

	/** Counts are per probe. Skipped pool pairs are not included in lowHitRejections;
	 * poolVisited states how many distinct residual tails reached query processing.
	 */
	public static final class Result{
		private Result(final int poolSize_,final int affected_,final int firstId_,final int secondId_,final Counts c){
			poolSize=poolSize_; affected=affected_; firstId=firstId_; secondId=secondId_;
			poolVisited=c.poolVisited; poolAccepted=c.poolAccepted; affectedAccepted=c.affectedAccepted;
			accepted=poolAccepted+affectedAccepted; netGain=accepted-affected;
			firstAccepted=c.first; secondAccepted=c.second;
			candidates=c.candidates; alignments=c.alignments; lowHitRejections=c.lowHits;
			hydroRejections=c.hydro; postingsVisited=c.postings;
			if(poolAccepted>poolVisited || poolVisited>poolSize || affectedAccepted>affected ||
				accepted!=firstAccepted+secondAccepted || candidates!=alignments+lowHitRejections+hydroRejections){
				throw new IllegalStateException("Local probe counts violate unique-tail or disjoint-screen accounting");
			}
		}
		public final int poolSize,affected,firstId,secondId;
		public final long poolVisited,poolAccepted,affectedAccepted,accepted,netGain,firstAccepted,secondAccepted;
		public final long candidates,alignments,lowHitRejections,hydroRejections,postingsVisited;
	}
	private static final class Counts{
		long poolVisited,poolAccepted,affectedAccepted,first,second,candidates,alignments,lowHits,hydro,postings;
	}
	private static void require(final boolean valid,final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}
	private static void requireOrdinal(final boolean valid,final String message,final int ordinal){
		if(!valid){throw new IllegalArgumentException(message+": "+ordinal);}
	}
	private final ProteinSequence[] sequences;
	private final TailPoolIndex pool;
	private final FamilySketchBuilder sketches;
	private final TailAlignmentGate gate;
	private final CandidateFilter filter;
	private final int minimumLength;
	private final BitSet poolOrdinals;
}
