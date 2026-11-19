package prot;

import java.util.HashSet;

import structures.IntHashMap;

/** Immutable, live-roster-bound tail-GC ranking penalty and optional sample-SD interval.
 * The caller supplies rows joined to the accepted family report and collection map;
 * this class validates their identities/counts but does not establish file provenance.
 * Query availability is an argument to every score, never retained worker state.
 * @author Keqing, Qiqi
 */
public final class FamilyGcPenalty{

	/** Arrays follow live-family ordinal order; physical ranks and stable IDs may have gaps.
	 * Missing means have an explicit false flag and finite zero storage, never NaN.
	 * The coefficient is a finite, nonnegative multiplier of GC fractions in [0,1].
	 */
	public FamilyGcPenalty(int[] stableIds_,int[] ranks_,String[] ids_,long[] members_,
			boolean[] available_,double[] means_,double coefficient_){
		this(stableIds_,ranks_,ids_,members_,available_,means_,new double[ids_==null ? 0 : ids_.length],coefficient_,false,0,30);
	}

	/** Sample SD describes the same member nucleotide-tail interval as the mean.
	 * Filtering uses a literal closed interval; missing data and small families are neutral. */
	public FamilyGcPenalty(int[] stableIds_,int[] ranks_,String[] ids_,long[] members_,
			boolean[] available_,double[] means_,double[] sampleSd_,double coefficient_,
			boolean filter_,double z_,int minimumMembers_){
		require(stableIds_!=null && ranks_!=null && ids_!=null && members_!=null && available_!=null && means_!=null,
			"GC binding requires complete identity, count, availability and mean arrays");
		final int n=ids_.length;
		require(n>0 && stableIds_.length==n && ranks_.length==n && members_.length==n && available_.length==n && means_.length==n,
			"Every GC array must describe the same nonempty live roster");
		require(Double.isFinite(coefficient_) && coefficient_>=0,"GC coefficient must be finite and nonnegative");
		require(sampleSd_!=null && sampleSd_.length==n,"GC sample SD must cover the same roster");
		require(Double.isFinite(z_) && z_>=0 && minimumMembers_>=2,"GC filter requires finite nonnegative z and at least two members");
		stableIds=stableIds_.clone();ranks=ranks_.clone();ids=ids_.clone();members=members_.clone();
		available=available_.clone();means=means_.clone();coefficient=coefficient_;
		sampleSd=sampleSd_.clone();filter=filter_;z=z_;minimumMembers=minimumMembers_;
		final IntHashMap seenStable=new IntHashMap(),seenRanks=new IntHashMap();
		final HashSet<String> seenIds=new HashSet<String>();
		for(int i=0;i<n;i++){
			require(stableIds[i]>=0 && seenStable.get(stableIds[i])<0,"Duplicate or negative stable GC family ID at ordinal "+i);
			require(ranks[i]>=0 && seenRanks.get(ranks[i])<0,"Duplicate or negative physical GC rank at ordinal "+i);
			require(ids[i]!=null && !ids[i].isEmpty() && ids[i].indexOf('\t')<0 && ids[i].indexOf('\n')<0 && ids[i].indexOf('\r')<0 && seenIds.add(ids[i]),
				"GC family labels must be unique nonempty table tokens at ordinal "+i);
			require(members[i]>0,"GC binding contains an empty family in the live roster at ordinal "+i);
			require(Double.isFinite(means[i]) && means[i]>=0 && means[i]<=1,"GC mean must be finite in [0,1] at ordinal "+i);
			require(available[i] || means[i]==0,"Unavailable GC mean must use explicit zero storage at ordinal "+i);
			require(Double.isFinite(sampleSd[i]) && sampleSd[i]>=0 && sampleSd[i]<=1,"GC sample SD must be finite in [0,1]");
			require(available[i] || sampleSd[i]==0,"Unavailable GC must use zero SD storage");
			require(members[i]>1 || sampleSd[i]==0,"A singleton has no sample SD");
			seenStable.put(stableIds[i],i);seenRanks.put(ranks[i],i);
		}
	}

	/** Verify the actual classifier roster before attaching this immutable binding.
	 * Stable numeric IDs are validated by the report/map join; the classifier exposes
	 * display IDs, physical ranks and counts, which must also match exactly here.
	 */
	public void requireRoster(FamilyTailCascade.LoadedClassifier classifier){
		require(classifier!=null && classifier.familyCount()==ids.length,"GC binding and classifier must have the same live family count");
		for(int i=0;i<ids.length;i++){
			requireRow(i,classifier.familyRankAt(i),classifier.familyIdAt(i),classifier.familyMembersAt(i));
		}
	}

	/** Derived search views may replace gates, but not the identity/count order of GC means. */
	void requireRoster(FamilyReference[] families){
		require(families!=null && families.length==ids.length,"GC binding must cover the complete derived live roster");
		for(int i=0;i<families.length;i++){
			final FamilyReference family=families[i];
			require(family!=null,"GC binding cannot retain an absent live family");
			requireRow(i,family.rank,family.id,family.members);
		}
	}

	/** Array form supports a validated resident roster before classifier construction. */
	public void requireRoster(int[] actualRanks,String[] actualIds,long[] actualMembers){
		require(actualRanks!=null && actualIds!=null && actualMembers!=null && actualRanks.length==ids.length &&
			actualIds.length==ids.length && actualMembers.length==ids.length,"GC binding requires one matching row per live family");
		for(int i=0;i<ids.length;i++){requireRow(i,actualRanks[i],actualIds[i],actualMembers[i]);}
	}

	/** Preserve baseline bits exactly for zero weight or unavailable query/family GC.
	 * Finite negative adjusted scores are valid ranking values, never clipped or gated.
	 * Unavailable query GC uses finite zero storage; a measured zero remains distinct.
	 */
	public double score(double baseline,int ordinal,boolean queryAvailable,double queryGc){
		requireOrdinal(ordinal);
		require(Double.isFinite(baseline),"GC ranking requires an already finite native baseline score");
		require(Double.isFinite(queryGc) && queryGc>=0 && queryGc<=1,"Query GC must be a finite fraction in [0,1]");
		require(queryAvailable || queryGc==0,"Unavailable query GC must not carry a previous query's value");
		if(coefficient==0 || !queryAvailable || !available[ordinal]){return baseline;}
		final double score=baseline-coefficient*Math.abs(queryGc-means[ordinal]);
		require(Double.isFinite(score),"GC penalty overflowed the finite native ranking domain");
		return score;
	}

	public int families(){return ids.length;}
	/** No epsilon or special widening for zero SD. */
	public boolean passes(int ordinal,boolean queryAvailable,double queryGc){
		requireOrdinal(ordinal);
		require(Double.isFinite(queryGc) && queryGc>=0 && queryGc<=1 && (queryAvailable || queryGc==0),
			"GC filter requires explicit availability and a finite query fraction");
		if(!filter || !queryAvailable || !available[ordinal] || members[ordinal]<minimumMembers){return true;}
		final double radius=z*sampleSd[ordinal];
		return queryGc>=means[ordinal]-radius && queryGc<=means[ordinal]+radius;
	}
	public boolean filterEnabled(){return filter;}
	public double sampleSdAt(int ordinal){requireOrdinal(ordinal);return sampleSd[ordinal];}
	public double coefficient(){return coefficient;}
	public boolean meanAvailableAt(int ordinal){requireOrdinal(ordinal);return available[ordinal];}
	public double meanAt(int ordinal){
		requireOrdinal(ordinal);require(available[ordinal],"An unavailable GC mean is not a measured zero");return means[ordinal];
	}
	public int stableIdAt(int ordinal){requireOrdinal(ordinal);return stableIds[ordinal];}

	private void requireRow(int ordinal,int rank,String id,long count){
		require(rank==ranks[ordinal] && ids[ordinal].equals(id) && count==members[ordinal],
			"GC identity/rank/member count differs from the live roster at ordinal "+ordinal);
	}
	private void requireOrdinal(int ordinal){
		if(ordinal<0 || ordinal>=ids.length){throw new IllegalArgumentException("GC scoring requires a live family ordinal: "+ordinal);}
	}
	private static void require(boolean valid,String message){if(!valid){throw new IllegalArgumentException(message);}}

	private final int[] stableIds,ranks;
	private final String[] ids;
	private final long[] members;
	private final boolean[] available;
	private final double[] means;
	private final double coefficient;
	private final double[] sampleSd;
	private final boolean filter;
	private final double z;
	private final int minimumMembers;
}
