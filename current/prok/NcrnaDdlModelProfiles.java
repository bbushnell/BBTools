package prok;

import cardinality.CardinalityTracker;
import cardinality.DynamicDemiLog;
import kmer.AbstractKmerTable;

/** Constructs one geometry's member-voted profiles in current library order.
 * The member loader has already bound exact construction IDs and dispositions.
 * Sequence arrays are borrowed read-only; only the finished profiles survive.
 * A native sketch contributes its one winning nonzero key per member/bucket.
 * This is not a union of member sketches or a vote weighted by sequence depth.
 *
 * One owner constructs profiles before ranking workers start. Process-global
 * native encoding settings must remain fixed while these profiles are used.
 * One native sketch is cleared between members; Votes consumes its borrowed
 * array before the next clear. No member sketch survives profile construction.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlModelProfiles{

	/*--------------------------------------------------------------*/
	/*----------------          Construction        ----------------*/
	/*--------------------------------------------------------------*/

	/** Hashes each retained member with the experiment's fixed native DNA codec.
	 * Each width is constructed independently, before any plurality is chosen.
	 * @return Independently owned profiles with the member loader's model order */
	public static NcrnaDdlProfile[] build(final NcrnaDdlMembers members, final int k, final int buckets){
		return build(members, k, buckets, true, false);
	}

	/** The new grid uses forward words, optionally with BBDuk's single-middle-base mask. */
	public static NcrnaDdlProfile[] build(final NcrnaDdlMembers members, final int k, final int buckets,
			final boolean canonical, final boolean middle){
		if(members==null || k<1 || k>31 || buckets<32 || buckets>2048 || Integer.bitCount(buckets)!=1){
			throw new IllegalArgumentException("Member profiles require bound members, k1..31 and power-of-two width32..2048");
		}
		if(DynamicDemiLog.exponentBits()!=4 || DynamicDemiLog.blacklistExists()){
			throw new IllegalStateException("The ranking experiment requires exponent4 and no process-global DDL blacklist");
		}
		if(middle && (canonical || k<3)){
			throw new IllegalArgumentException("Middle-masked member profiles require forward words with k at least3");
		}
		final NcrnaDdlProfile[] result=new NcrnaDdlProfile[members.models()];
		final DynamicDemiLog sketch=DynamicDemiLog.create(buckets, k, 12345, 0, false, false);
		final long wordMask=AbstractKmerTable.makeMiddleMask(k, false, middle ? 1 : 0);
		final CardinalityTracker.KmerConsumer forward=canonical ? null : (word, length)->sketch.add(word&wordMask);
		for(int model=0; model<result.length; model++){
			final NcrnaDdlProfile.Votes votes=new NcrnaDdlProfile.Votes(buckets);
			final byte[][] sequences=members.members(model);
			for(final byte[] sequence:sequences){
				sketch.clear();
				if(canonical){sketch.hash(sequence, null);}
				else{CardinalityTracker.forEachKmer(sequence, k, k, false, false, forward);}
				votes.add(sketch.maxArray());
			}
			result[model]=votes.finish();
			assert(result[model].members==sequences.length)
				: "Every kept construction member must contribute exactly one unweighted native sketch: model="+model;
		}
		return result;
	}

	/** Both lengths compete inside each member before its one vote per bucket.
	 * A partitioned geometry instead reserves one half per length; plurality is
	 * still applied to actual per-member winners, never pooled raw word counts.
	 */
	public static NcrnaDdlProfile[] buildDual(final NcrnaDdlMembers members, final int shortK, final int longK,
			final int buckets, final boolean partitioned){
		if(members==null || members.models()<1){
			throw new IllegalArgumentException("Dual profiles require bound construction members in actual model order");
		}
		final NcrnaDdlDualKeys query=new NcrnaDdlDualKeys(shortK, longK, buckets, partitioned);
		final NcrnaDdlProfile[] result=new NcrnaDdlProfile[members.models()];
		for(int model=0; model<result.length; model++){
			final NcrnaDdlProfile.Votes votes=new NcrnaDdlProfile.Votes(buckets);
			final byte[][] sequences=members.members(model);
			for(final byte[] sequence:sequences){
				query.build(sequence);
				votes.add(query.memberKeys());
			}
			result[model]=votes.finish();
			assert(result[model].members==sequences.length)
				: "Every retained member supplies one combined sketch vote, irrespective of word count";
		}
		return result;
	}

	/** Unpartitioned forward masked10 profiles for the compact consensus study.
	 * Hash members once per width, then select every requested consensus weight
	 * from the same unchanged vote counts. Result order is weight, then model.
	 * All raw consensus keys are eligible; member contributions remain one native
	 * maximum per bucket. Zero weight is the highest-key member-only control.
	 */
	public static NcrnaDdlProfile[][] buildConsensus(final NcrnaDdlMembers members, final byte[][] consensuses,
			final int buckets, final int[] weights){
		if(members==null || consensuses==null || consensuses.length!=members.models()
				|| buckets<32 || buckets>2048 || Integer.bitCount(buckets)!=1 || weights==null || weights.length<1){
			throw new IllegalArgumentException("Compact consensus profiles need model-ordered sequences, weights and power-of-two width32..2048");
		}
		if(DynamicDemiLog.exponentBits()!=4 || DynamicDemiLog.blacklistExists()){
			throw new IllegalStateException("Compact profiles require fixed exponent4 and no process-global DDL blacklist");
		}
		for(final int weight:weights){
			if(weight<NcrnaDdlProfile.CONSENSUS_WINS){throw new IllegalArgumentException("Consensus weight must be nonnegative or CONSENSUS_WINS");}
		}
		for(int model=0; model<consensuses.length; model++){
			final byte[] sequence=consensuses[model];
			if(sequence==null){throw new IllegalArgumentException("Every model needs its bound consensus sequence");}
			for(int offset=0; offset<sequence.length; offset++){
				final int base=sequence[offset]&~32;
				if(base!='A' && base!='C' && base!='G' && base!='T' && base!='N'){
					throw new IllegalArgumentException("The compact native DNA codec requires ACGTN consensus bases: model="+model+", offset="+offset);
				}
			}
		}
		final NcrnaDdlProfile[][] result=new NcrnaDdlProfile[weights.length][members.models()];
		final DynamicDemiLog sketch=DynamicDemiLog.create(buckets, 10, 12345, 0, false, false);
		final long wordMask=AbstractKmerTable.makeMiddleMask(10, false, 1);
		final CardinalityTracker.KmerConsumer forward=(word, length)->sketch.add(word&wordMask);
		final NcrnaDdlQueryKeys consensus=new NcrnaDdlQueryKeys(buckets, 10);
		for(int model=0; model<members.models(); model++){
			final NcrnaDdlProfile.Votes votes=new NcrnaDdlProfile.Votes(buckets);
			final byte[][] sequences=members.members(model);
			for(final byte[] sequence:sequences){
				sketch.clear();
				CardinalityTracker.forEachKmer(sequence, 10, 10, false, false, forward);
				votes.add(sketch.maxArray());
			}
			consensus.build(consensuses[model], null, false, true);
			for(int variant=0; variant<weights.length; variant++){
				result[variant][model]=votes.finish(consensus, weights[variant]);
				assert(result[variant][model].members==sequences.length)
					: "Consensus pseudo-votes must not change the kept construction-member count";
			}
		}
		return result;
	}

	/** Member sketches vote only in their own partitions; the model consensus
	 * supplies an independent native maximum sketch in the final quarter. */
	public static NcrnaDdlProfile[] buildPartition(final NcrnaDdlMembers members, final byte[][] consensuses,
			final int shortK, final int buckets){
		return buildPartition(members, consensuses, shortK, buckets, shortK!=0);
	}

	/** Selects consensus content while retaining only its plain DDL keys. */
	public static NcrnaDdlProfile[] buildPartition(final NcrnaDdlMembers members, final byte[][] consensuses,
			final int shortK, final int buckets, final boolean consensusShort){
		if(members==null || members.models()<1 || consensuses==null || consensuses.length!=members.models()){
			throw new IllegalArgumentException("Partition profiles require one bound consensus per construction model");
		}
		final NcrnaDdlPartitionKeys codec=new NcrnaDdlPartitionKeys(shortK, buckets, consensusShort);
		final NcrnaDdlProfile[] result=new NcrnaDdlProfile[members.models()];
		for(int model=0; model<result.length; model++){
			if(consensuses[model]==null){throw new IllegalArgumentException("Partition model lacks its actual consensus sequence: "+model);}
			final NcrnaDdlProfile.Votes votes=new NcrnaDdlProfile.Votes(buckets);
			final byte[][] sequences=members.members(model);
			for(final byte[] sequence:sequences){votes.add(codec.memberKeys(sequence));}
			result[model]=votes.finishPartition(codec.consensusKeys(consensuses[model]));
			assert(result[model].members==sequences.length)
				: "A separate consensus partition cannot change the retained construction-member count";
		}
		return result;
	}

	private NcrnaDdlModelProfiles(){}
}
