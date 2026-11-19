package prot.cluster;

import java.util.Arrays;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Per-family identity and BLOSUM cohesion for the final published consensus. */
public final class ConsensusQualityReport{

	public ConsensusQualityReport(){}

	public void publish(Snapshot snapshot){
		if(snapshot==null){return;}
		ensure(snapshot.familyId+1);
		if(snapshots[snapshot.familyId]==null){published++;}
		snapshots[snapshot.familyId]=snapshot;
	}

	public Snapshot snapshot(int familyId){
		if(familyId<0){throw new IllegalArgumentException("Negative family ID: "+familyId);}
		final Snapshot snapshot=familyId<snapshots.length ? snapshots[familyId] : null;
		if(snapshot==null){throw new IllegalStateException("Missing consensus quality for family "+familyId);}
		return snapshot;
	}

	public int published(){return published;}

	public void write(String path,boolean overwrite,Family[] retained){
		if(path==null || path.isEmpty()){throw new IllegalArgumentException("Quality report output path is required");}
		if(retained==null){throw new IllegalArgumentException("Retained family list is required");}
		final ByteBuilder bb=new ByteBuilder();
		bb.append("family_id\tmembers\tvalue_available\tidentity_mean\tscore_ratio_mean\tunavailable_reason\n");
		for(Family family : retained){
			if(family==null){throw new IllegalArgumentException("Retained family list contains null");}
			final Snapshot snapshot=family.familyId<snapshots.length ? snapshots[family.familyId] : null;
			if(snapshot==null){throw new IllegalStateException("Missing consensus quality for family "+family.familyId);}
			if(snapshot.members!=family.size()){
				throw new IllegalStateException("Consensus quality member count differs for family "+family.familyId+
					": snapshot="+snapshot.members+" retained="+family.size());
			}
			bb.append(family.familyId).tab().append(snapshot.members).tab().append(snapshot.available ? 1 : 0).tab();
			if(snapshot.available){
				bb.appendSlow(snapshot.identityMean).tab().appendSlow(snapshot.scoreRatioMean).tab();
			}else{bb.append("NA\tNA\t").append(snapshot.reason);}
			bb.nl();
		}
		final ByteStreamWriter writer=new ByteStreamWriter(path,overwrite,false,false);
		writer.start();
		writer.print(bb);
		if(writer.poisonAndWait()){throw new IllegalStateException("Consensus quality report output failed: "+path);}
	}

	public static Snapshot available(int familyId,int members,double identitySum,double scoreRatioSum){
		validate(familyId,members);
		if(!Double.isFinite(identitySum) || identitySum<0 || identitySum>members){
			throw new IllegalArgumentException("Identity sum must be finite in [0,members]: "+identitySum);
		}
		if(!Double.isFinite(scoreRatioSum)){
			throw new IllegalArgumentException("Score-ratio sum must be finite: "+scoreRatioSum);
		}
		final double identityMean=identitySum/members, scoreRatioMean=scoreRatioSum/members;
		if(!Double.isFinite(identityMean) || !Double.isFinite(scoreRatioMean)){
			throw new IllegalArgumentException("Consensus quality means must be finite");
		}
		return new Snapshot(familyId,members,true,identityMean,scoreRatioMean,"");
	}

	public static Snapshot unavailable(int familyId,int members,String reason){
		validate(familyId,members);
		if(reason==null || reason.isEmpty()){throw new IllegalArgumentException("Unavailable reason is required");}
		for(int i=0; i<reason.length(); i++){
			final char c=reason.charAt(i);
			if(c<=32 || c>=127){throw new IllegalArgumentException("Unavailable reason must be a printable ASCII token: "+reason);}
		}
		return new Snapshot(familyId,members,false,-1,-1,reason);
	}

	private static void validate(int familyId,int members){
		if(familyId<0){throw new IllegalArgumentException("Negative family ID: "+familyId);}
		if(members<1){throw new IllegalArgumentException("Consensus quality needs at least one member: "+members);}
	}

	private void ensure(int size){
		if(size>snapshots.length){
			int next=snapshots.length;
			while(next<size){next=Math.multiplyExact(next,2);}
			snapshots=Arrays.copyOf(snapshots,next);
		}
	}

	public static final class Snapshot{
		private Snapshot(int familyId_,int members_,boolean available_,double identityMean_,
				double scoreRatioMean_,String reason_){
			familyId=familyId_; members=members_; available=available_;
			identityMean=identityMean_; scoreRatioMean=scoreRatioMean_; reason=reason_;
		}
		public int familyId(){return familyId;}
		public int members(){return members;}
		public boolean available(){return available;}
		public double identityMean(){return identityMean;}
		public double scoreRatioMean(){return scoreRatioMean;}
		public String reason(){return reason;}
		private final int familyId,members;
		private final boolean available;
		private final double identityMean,scoreRatioMean;
		private final String reason;
	}

	private Snapshot[] snapshots=new Snapshot[16];
	private int published;
}
