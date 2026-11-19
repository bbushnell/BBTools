package prot.cluster;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import ddl.FamilySketchBuilder;
import dna.AminoAcid;
import prot.Blosum62;
import prot.ConsensusScalars;
import prot.ProteinSequence;

/** Family model and shared membership. Concurrent callers finish model/scalar
 * mutation before publication; append-only membership and its optional coherent
 * member-ranking descriptors may change afterward. */
public final class Family {

	public Family(int familyId_,String id_,int representativeIndex_,byte[] representative_){
		this(familyId_,id_,representativeIndex_,representative_,false);
	}

	/** A consensus input represents this many original tails during size-restricted merging. */
	public static Family mergeSeed(int familyId,String id,int input,byte[] sequence,int originalMembers){
		if(originalMembers<1){throw new IllegalArgumentException("Merge seeds need a positive original-tail count");}
		final Family family=new Family(familyId,id,input,sequence);
		family.mergeMembers=new AtomicLong(originalMembers);
		return family;
	}

	/** Seed-only models own no current input ordinals; an exported count is metadata. */
	private Family(int familyId_,String id_,int representativeIndex_,byte[] representative_,boolean seedOnly){
		if(familyId_<0){throw new IllegalArgumentException("Negative family ID");}
		familyId=familyId_;
		id=id_;sketchIdentity=id_;
		validateId(id);
		if(seedOnly ? representativeIndex_!=-1 : representativeIndex_<0){
			throw new IllegalArgumentException("An assigned family needs a representative input; a seed-only model uses -1");
		}
		representativeIndex=representativeIndex_;
		representative=copySequence(representative_, "representative");
		consensus=representative.clone();
		members=new FamilyMembers();
		if(!seedOnly){members.add(representativeIndex);}
		nextSketchSize=1;
		nextConsensusSize=5;
		scalars=new double[ConsensusScalars.SCALAR_COUNT];
		monomers=new float[MONOMER_WIDTH];
		monomerCounts=new long[MONOMER_WIDTH];
	}

	/** Loads a model for fixed-seed reassignment without inventing original member IDs.
	 * The first frozen pass may report the exported count; after membership clearing,
	 * every count comes solely from actual assignments in the new input population.
	 */
	public static Family fromSeed(ClusterSketchCodec.Entry entry){
		return fromSeed(entry, false);
	}

	/** Loads a model for fixed-seed reassignment, optionally accepting an empty template. */
	public static Family fromSeed(ClusterSketchCodec.Entry entry,boolean allowEmpty){
		final int min=allowEmpty ? 0 : 1;
		if(entry==null || entry.members<min){
			throw new IllegalArgumentException("A fixed seed needs a "+(allowEmpty ? "nonnegative" : "nonempty")+" source-family count");
		}
		final Family seed=new Family(entry.familyId,entry.name,-1,entry.consensus,true);
		seed.setSketch(entry.sketch); seed.sketchIdentity=entry.sketchIdentity; seed.setMonomers(entry.monomers);
		seed.setScalars(entry.scalars,entry.scalarMask);
		seed.seedReferenceMembers=entry.members;
		assert(seed.size()==0) : "Loading a seed count must never fabricate input-ordinal membership";
		return seed;
	}

	/** Regenerated versions share append-only membership but own model/statistic arrays. */
	public Family(Family source){
		familyId=source.familyId; id=source.id; sketchIdentity=source.sketchIdentity; representativeIndex=source.representativeIndex;
		representative=source.representative.clone(); consensus=source.consensus.clone();
		rawConsensus=null; members=source.members;
		nextSketchSize=source.nextSketchSize; nextConsensusSize=source.nextConsensusSize;
		scalarMask=source.scalarMask; occupiedBins=source.occupiedBins;
		ddl=source.ddl==null ? null : source.ddl.clone();
		monomers=source.monomers.clone();
		scalars=source.scalars.clone();
		monomerCounts=new long[MONOMER_WIDTH];
		seedReferenceMembers=source.seedReferenceMembers;
		mergeMembers=source.mergeMembers;
	}

	public int size(){return members.size();}
	/** Source count is usable only before the first assignment pass consumes the seed. */
	public int referenceSize(){return size()>0 ? size() : seedReferenceMembers;}
	public boolean isEmpty(){return size()==0;}
	public int memberIndex(int i){
		return members.get(i);
	}
	public int[] memberIndexes(){return members.toArray();}
	public MembershipView membershipView(){return new MembershipView(members);}
	boolean sharesMembers(Family other){return other!=null && members==other.members;}

	/** Attach descriptors before publication; model copies share this live member-prefix state. */
	public void enableMemberRanking(final MemberRankInputs inputs, final int layer, final MemberRankStatistics.Snapshot seed){
		if(mergeMembers!=null){
			throw new IllegalStateException("Consensus-input merge weights require original-tail descriptor binding, not consensus-as-member statistics");
		}
		members.enableRanking(inputs,layer,seed);
	}
	public boolean hasMemberRanking(){return members.hasRanking();}
	/** Publication may borrow descriptors only from the same original ordinal cache. */
	boolean memberRankingUses(final MemberRankInputs inputs){return members.rankingUses(inputs);}
	/** Capture once before a frozen assignment pass, never once per candidate. */
	public MemberRankStatistics.Snapshot memberRankingSnapshot(){return members.rankingSnapshot();}
	/** Fill a query worker's private view without allocating or rescanning members. */
	public void readMemberRanking(final MemberRankStatistics.View out){members.readRanking(out);}
	public void memberRankingLayer(final int layer){members.rankingLayer(layer);}

	public void addMember(int index){
		if(mergeMembers!=null){throw new IllegalStateException("Weighted merge families require atomic original-tail admission");}
		members.add(index);
	}

	/** Current original-tail weight, shared across every immutable version of a merge family. */
	public long mergeMemberCount(){
		if(mergeMembers==null){throw new IllegalStateException("Family has no original-tail merge weights: "+id);}
		return mergeMembers.get();
	}

	/** Reserves size eligibility atomically before publishing the consensus-input ordinal.
	 * A failed reserve changes neither count nor membership. An append failure aborts the run.
	 */
	public boolean tryAddMergeMember(int input,int originalMembers,int maximumSize){
		if(input<0 || input==Integer.MAX_VALUE || originalMembers<1 || maximumSize<1 || mergeMembers==null){
			throw new IllegalArgumentException("Weighted merge admission needs a valid ordinal, positive weight and size limit");
		}
		for(;;){
			final long before=mergeMembers.get();
			assert(before>0) : "Weighted merge groups begin with a positive seed count and only grow";
			if(originalMembers>=maximumSize && before>=maximumSize){return false;}
			final long after=Math.addExact(before,originalMembers);
			if(mergeMembers.compareAndSet(before,after)){
				members.add(input);
				return true;
			}
		}
	}

	/** Clears assignment for a frozen competitive pass; sequence payloads remain intact. */
	public void clearMembers(){
		if(mergeMembers!=null){throw new IllegalStateException("Weighted merge groups cannot enter reassignment");}
		members.clear(); seedReferenceMembers=0;
	}

	public boolean shouldRebuildSketch(){return size()>=nextSketchSize;}

	public void markSketchRebuilt(){
		markSketchRebuilt(size());
	}
	public void markSketchRebuilt(MembershipView view){
		validate(view);
		markSketchRebuilt(view.size());
	}
	public void markSketchRebuilt(int rebuiltSize){
		if(rebuiltSize<0){throw new IllegalArgumentException("Negative rebuilt size");}
		while(nextSketchSize<=rebuiltSize){
			nextSketchSize=Math.multiplyExact(nextSketchSize, 2);
		}
	}

	public boolean shouldRebuildConsensus(){return size()>=nextConsensusSize;}

	public void markConsensusRebuilt(){
		markConsensusRebuilt(size());
	}
	public void markConsensusRebuilt(MembershipView view){
		validate(view);
		markConsensusRebuilt(view.size());
	}
	public void markConsensusRebuilt(int rebuiltSize){
		if(rebuiltSize<0){throw new IllegalArgumentException("Negative rebuilt size");}
		while(nextConsensusSize<=rebuiltSize){
			nextConsensusSize=Math.addExact(Math.multiplyExact(nextConsensusSize-1, 2), 1);
		}
	}

	public void setRepresentative(byte[] sequence){representative=copySequence(sequence, "representative");}
	public void setConsensus(byte[] sequence){consensus=copySequence(sequence, "consensus");}
	public void setSketch(char[] sketch){
		ddl=sketch==null ? null : sketch.clone();
		occupiedBins=countOccupied(ddl);
	}
	public void setMonomers(float[] values){
		if(values==null || values.length!=monomers.length){throw new IllegalArgumentException("Need 20 monomer frequencies");}
		for(int i=0; i<values.length; i++){
			if(!Float.isFinite(values[i]) || values[i]<0){throw new IllegalArgumentException("Invalid monomer frequency");}
			monomers[i]=values[i];
		}
	}
	public void setScalars(double[] values,int mask){
		if(values==null || values.length!=scalars.length || mask<0 || mask>ConsensusScalars.ALL_DEFINED){
			throw new IllegalArgumentException("Need five scalars and a five-bit validity mask");
		}
		for(int i=0; i<values.length; i++){
			if(!Double.isFinite(values[i])){throw new IllegalArgumentException("Invalid scalar");}
			scalars[i]=values[i];
		}
		scalarMask=mask;
	}

	public void refreshStats(ProteinSequence[] sequences,FamilySketchBuilder sketchBuilder){
		refreshStats(sequences, sketchBuilder, membershipView());
	}

	public void refreshStats(ProteinSequence[] sequences,FamilySketchBuilder sketchBuilder,MembershipView view){
		if(sequences==null || sketchBuilder==null){throw new IllegalArgumentException("Missing sequence store or sketch builder");}
		validate(view);
		final int size=view.size();
		if(size<1){throw new IllegalStateException("Cannot refresh an empty family");}
		if(ddl==null || ddl.length!=sketchBuilder.width()){ddl=new char[sketchBuilder.width()];}
		sketchBuilder.reset();
		Arrays.fill(monomerCounts,0);
		long residues=0;
		for(int i=0; i<size; i++){
			final int member=view.memberIndex(i);
			if(member<0 || member>=sequences.length || sequences[member]==null){
				throw new IllegalStateException("Missing sequence for family "+id+" member index "+member);
			}
			final byte[] enc=sequences[member].enc;
			sketchBuilder.addEncoded(enc);
			for(byte residue : enc){
				if(residue>=0 && residue<MONOMER_WIDTH){
					monomerCounts[residue]++;
					residues++;
				}else if(residue!=Blosum62.X_CODE){
					throw new IllegalStateException("Invalid encoded residue in "+id+": "+residue);
				}
			}
		}
		sketchBuilder.writeProfile(ddl,id,consensus);
		sketchIdentity=id;
		occupiedBins=countOccupied(ddl);
		if(residues==0){
			Arrays.fill(monomers,0);
		}else{
			for(int i=0; i<monomers.length; i++){monomers[i]=(float)(monomerCounts[i]/(double)residues);}
		}
		scalarMask=ConsensusScalars.computeScalarMask(rawConsensus(), consensus.length, scalars, 0);
	}

	public String sketchIdentity(){return sketchIdentity;}

	public byte[] representative(){return representative.clone();}
	public byte[] consensus(){return consensus.clone();}
	byte[] consensusUnsafe(){return consensus;}
	public char[] sketch(){return ddl==null ? null : ddl.clone();}
	char[] sketchUnsafe(){return ddl;}
	/** Cached nonzero registers for this exact immutable model version. */
	public int occupiedBins(){return occupiedBins;}
	public float[] monomers(){return monomers.clone();}
	float[] monomersUnsafe(){return monomers;}
	public double[] scalars(){return scalars.clone();}
	double[] scalarsUnsafe(){return scalars;}
	public int scalarMask(){return scalarMask;}

	private static int countOccupied(char[] sketch){
		int count=0;
		if(sketch!=null){for(char value : sketch){if(value!=0){count++;}}}
		assert(sketch==null || count<=sketch.length) : "Flex hit thresholds count at most one occupied value per register";
		return count;
	}

	private static byte[] copySequence(byte[] sequence,String name){
		if(sequence==null || sequence.length<1){throw new IllegalArgumentException("Empty "+name);}
		final byte[] copy=sequence.clone();
		for(byte residue : copy){
			if(residue<0 || residue>21 || residue==20){
				throw new IllegalArgumentException("Invalid "+name+" residue: "+residue);
			}
		}
		return copy;
	}

	private static void validateId(String id){
		if(id==null || id.isEmpty()){throw new IllegalArgumentException("Missing family ID");}
		for(int i=0; i<id.length(); i++){
			final char c=id.charAt(i);
			if(c<=32 || c>=127){throw new IllegalArgumentException("Family ID must be a printable ASCII token: "+id);}
		}
	}

	private byte[] rawConsensus(){
		if(rawConsensus==null || rawConsensus.length<consensus.length){
			rawConsensus=new byte[consensus.length];
		}
		for(int i=0; i<consensus.length; i++){
			final byte residue=consensus[i];
			rawConsensus[i]=residue==Blosum62.X_CODE ? (byte)'X' : AminoAcid.numberToAcid[residue];
		}
		return rawConsensus;
	}

	private void validate(MembershipView view){
		if(view==null || view.members!=members){throw new IllegalArgumentException("Membership view belongs to a different family");}
	}

	/** Fixed published-prefix view for a single model/statistic rebuild. */
	public static final class MembershipView {
		private MembershipView(FamilyMembers members_){
			members=members_;
			size=members.size();
		}
		public int size(){return size;}
		public int memberIndex(int index){
			if(index<0 || index>=size){throw new ArrayIndexOutOfBoundsException(index);}
			return members.get(index);
		}
		public int[] memberIndexes(){
			final int[] copy=new int[size];
			for(int i=0; i<size; i++){copy[i]=memberIndex(i);}
			return copy;
		}
		boolean viewOf(Family family){return family!=null && members==family.members;}
		private final FamilyMembers members;
		private final int size;
	}

	public final int familyId;
	public final String id;
	public final int representativeIndex;
	private int nextSketchSize, nextConsensusSize, scalarMask,occupiedBins,seedReferenceMembers;
	private final FamilyMembers members;
	private AtomicLong mergeMembers;
	private byte[] representative, consensus, rawConsensus;
	private char[] ddl;
	private String sketchIdentity;
	private final float[] monomers;
	private final double[] scalars;
	private final long[] monomerCounts;
	private static final int MONOMER_WIDTH=20;
}
