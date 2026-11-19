package prot;

import java.nio.file.Path;

/** Immutable caller-side fourth-root family rank score bound to one live roster.
 * This deliberately reuses the diagnostic FamilyHeuristicRanking arithmetic;
 * it adds only production binding, pinning and one selected arm.
 * @author Yelan
 */
public final class FamilyHeuristicScore{

	public FamilyHeuristicScore(FamilyTailCascade.LoadedClassifier base,Path medians_,String medianPin_,
			Path familyHydro_,String hydroPin_,int minAa,int maxAa,String armName_){
		this(base, medians_, medianPin_, familyHydro_, hydroPin_, minAa, maxAa,
			armName_, FamilyCalibration.VerificationMode.SHA80);
	}

	/** Compact family resources fold the selected live-order median and KD arrays into the family file. */
	public FamilyHeuristicScore(FamilyTailCascade.LoadedClassifier base,double[] medians_,
			double[] hydrophobicities_,int minAa,int maxAa,String armName_){
		require(base!=null && medians_!=null && armName_!=null,
			"Heuristic scoring needs a classifier, medians and one arm");
		require(minAa>=1 && maxAa>=minAa && maxAa<=200,
			"Heuristic query bounds must be an inclusive interval within1..200; member bounds are separate");
		config=base.config();minQueryAa=minAa;maxQueryAa=maxAa;
		require(boundTo(config),
			"Heuristic scoring requires combined25/ranking, all eligible families and no learned rank model");
		arm=arm(armName_);armName=FamilyHeuristicRanking.armName(arm);
		require(arm==6 || arm>=8 && arm<FamilyHeuristicRanking.ALL_ARMS,
			"Caller heuristic arm must be hitfourth or one of its released H/cosine variants: "+armName_);
		require(usesHydrophobicity()==(hydrophobicities_!=null),
			"Hydrophobicity arms require compact KD means; H0 must not receive them");
		metadata=new FamilyHeuristicMetadata(config, medians_, hydrophobicities_);
	}

	/** Explicit inherited-identity mode affects input verification, never ranking arithmetic. */
	public FamilyHeuristicScore(FamilyTailCascade.LoadedClassifier base,Path medians_,String medianPin_,
			Path familyHydro_,String hydroPin_,int minAa,int maxAa,String armName_,
			FamilyCalibration.VerificationMode verificationMode){
		require(base!=null && medians_!=null && armName_!=null,"Heuristic scoring needs a classifier, medians and one arm");
		require(minAa>=1 && maxAa>=minAa && maxAa<=200,
			"Heuristic query bounds must be an inclusive interval within1..200; member bounds are separate");
		config=base.config();minQueryAa=minAa;maxQueryAa=maxAa;
		require(boundTo(config),
			"Heuristic scoring requires combined25/ranking, all eligible families and no learned rank model");
		arm=arm(armName_);armName=FamilyHeuristicRanking.armName(arm);
		require(arm==6 || arm>=8 && arm<FamilyHeuristicRanking.ALL_ARMS,
			"Caller heuristic arm must be hitfourth or one of its released H/cosine variants: "+armName_);
		require(usesHydrophobicity() ? familyHydro_!=null && hydroPin_!=null : familyHydro_==null && hydroPin_==null,
			"Hydrophobicity arms require their pinned KD table; H0 must not receive one");
		metadata=new FamilyHeuristicMetadata(config, medians_, medianPin_, familyHydro_, hydroPin_, verificationMode);
	}

	/** Bind the existing immutable metadata to a validated profile-only index view. */
	FamilyHeuristicScore rebindIndex(FamilyTailCascade.SearchConfig target){
		return new FamilyHeuristicScore(this,target);
	}

	/** Profile-only views share already validated immutable member metadata. They
	 * cannot change its reference roster or ranking coordinates, and read no file. */
	private FamilyHeuristicScore(FamilyHeuristicScore source,FamilyTailCascade.SearchConfig target){
		require(source!=null && target!=null && target.families==source.config.families &&
			target.geometry==source.config.geometry && target.width==source.config.width && target.seed==source.config.seed &&
			target.scalarBounds==source.config.scalarBounds && target.memberScalarMeans==source.config.memberScalarMeans &&
			target.features==source.config.features && target.cosineMode==source.config.cosineMode,
			"Profile-only heuristic rebinding must preserve the exact roster, geometry and ranking coordinates");
		config=target;metadata=source.metadata;minQueryAa=source.minQueryAa;maxQueryAa=source.maxQueryAa;
		arm=source.arm;armName=source.armName;
		require(boundTo(target),"Replacement index must retain complete combined25 ranking and no learned model");
	}

	/** Score one already prepared candidate. H0 deliberately ignores queryMean. */
	public double score(int queryLength,int ordinal,int hits,double cosine,double queryMean){
		validateQueryLength(queryLength);
		require(ordinal>=0 && ordinal<config.families.length,"Heuristic score needs a live family ordinal");
		if(!usesHydrophobicity()){return FamilyHeuristicRanking.score(arm,hits,cosine,queryLength,metadata.median(ordinal));}
		final double similarity=FamilyHeuristicRanking.hydrophobicitySimilarity(queryMean,metadata.hydrophobicity(ordinal));
		return FamilyHeuristicRanking.score(arm,hits,cosine,queryLength,metadata.median(ordinal),similarity);
	}

	/** Apply the released query-length guard before any DDL early exit. */
	public void validateQueryLength(int length){
		if(length<minQueryAa || length>maxQueryAa){
			throw new IllegalArgumentException("Heuristic query length is outside its released interval: "+
				length+" not in "+minQueryAa+".."+maxQueryAa);
		}
	}

	public boolean usesHydrophobicity(){return arm>=8;}
	public String armName(){return armName;}
	/** Read-only diagnostic access to the same validated member median used by score(). */
	public double familyMedianAt(int ordinal){
		return metadata.median(ordinal);
	}
	/** Read-only diagnostic access; absent KD binding is not a measured zero. */
	public double familyHydrophobicityAt(int ordinal){
		return metadata.hydrophobicity(ordinal);
	}

	/** Reject changed contents before and after a reproducible caller run. */
	public void verify(){
		metadata.verify();
	}

	/** The classifier may change gates and immutable member means; the search index and live roster must not. */
	boolean boundTo(FamilyTailCascade.SearchConfig candidate){
		if(candidate==null || candidate.index!=config.index || candidate.families.length!=config.families.length ||
				candidate.index.usesPruning() || candidate.index.usesDirectCompletion() ||
				candidate.scalarBounds==null || candidate.cosineMode!=FamilyTailCascade.CosineMode.RANKING ||
				candidate.features!=FamilyTailCascade.RankingFeatures.COMBINED25 ||
				candidate.lengthRanking!=null){
			return false;
		}
		for(int i=0;i<config.families.length;i++){
			final FamilyReference a=config.families[i],b=candidate.families[i];
			if(a.rank!=b.rank || a.members!=b.members || !a.id.equals(b.id)){return false;}
		}
		return true;
	}

	private static int arm(String name){
		for(int i=0;i<FamilyHeuristicRanking.ALL_ARMS;i++){
			if(FamilyHeuristicRanking.armName(i).equals(name)){return i;}
		}
		throw new IllegalArgumentException("Unknown heuristic arm: "+name);
	}

	private static void require(boolean valid,String why){if(!valid){throw new IllegalArgumentException(why);}}

	private final FamilyTailCascade.SearchConfig config;
	private final int minQueryAa,maxQueryAa;
	private final String armName;
	private final FamilyHeuristicMetadata metadata;
	private final int arm;
}
