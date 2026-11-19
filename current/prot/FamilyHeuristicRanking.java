package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;

import fileIO.ByteFile;
import parse.LineParser1;
import prot.FamilyTailCascade.Comparison;
import prot.FamilyTailCascade.SearchConfig;
import structures.IntHashMap;

/** Rank-only length comparisons over one immutable classifier roster.
 * Provided targets are memberships or explicitly frozen earlier oracle labels;
 * the alternative is the best calibrated-metric alignment in the active strict
 * top50 lists' union. The caller retains the provided label's provenance.
 * No classifier threshold, family, index or production default is modified.
 * @author Brian, Keqing
 */
public final class FamilyHeuristicRanking{

	public FamilyHeuristicRanking(FamilyTailCascade.LoadedClassifier base,Path medians_,
			String pin_,int minQueryAa,int maxQueryAa)throws Exception{
		this(base,medians_,pin_,minQueryAa,maxQueryAa,"six");
	}
	/** Explicit arm selection preserves the original six-arm API and output by default. */
	public FamilyHeuristicRanking(FamilyTailCascade.LoadedClassifier base,Path medians_,
			String pin_,int minQueryAa,int maxQueryAa,String mode)throws Exception{
		this(base,medians_,pin_,minQueryAa,maxQueryAa,mode,null,null);
	}
	/** Hydrophobicity arms require an explicitly pinned retained-member mean table; never use consensus as fallback. */
	public FamilyHeuristicRanking(FamilyTailCascade.LoadedClassifier base,Path medians_,
			String pin_,int minQueryAa,int maxQueryAa,String mode,Path familyHydro,
			String familyHydroPin)throws Exception{
		require(base!=null && medians_!=null,"Ranking requires a loaded package and pinned member medians");
		activeArms=selectArms(mode);
		cosineSweep=isCosineMode(mode);rawAdmission=isPivotMode(mode);
		range=new MemberLengthRange(minQueryAa,maxQueryAa);
		require(minQueryAa==40 && maxQueryAa==99 || minQueryAa==100 && maxQueryAa==149,
			"Released ranking query intervals are40..99 or100..149, without clipping");
		original=base.config();
		require(original.top==50 && original.width==1024 && original.geometry.isExact(5) &&
			original.features==FamilyTailCascade.RankingFeatures.COMBINED25 &&
			original.cosineMode==FamilyTailCascade.CosineMode.RANKING && original.scalarBounds!=null &&
			!original.index.usesPruning() && original.lengthRanking==null &&
			original.heuristicRanking==null && original.ddlAccept>1024,
			"Rank comparisons require the unchanged combined25 top50 package without "+
			"learned ranking, pruning or sketch acceptance");
		families=original.families;require(families.length>0,"Rank comparisons need a nonempty live roster");
		require(isHydrophobicityMode(mode) ? familyHydro!=null && familyHydroPin!=null :
			familyHydro==null && familyHydroPin==null,
			"Hydro mode requires its pinned retained-member means; other modes must not silently ignore a mean table");
		hydroMeans=familyHydro==null ? null : familyHydro.toAbsolutePath().normalize();hydroPin=familyHydroPin;
		medians=medians_.toAbsolutePath().normalize();pin=pin_;verify();lengths=readMedians();
		familyHydrophobicity=hydroMeans==null ? null : readHydrophobicity();verify();
		physicalToLive=new IntHashMap(families.length*2+1);
		for(int i=0;i<families.length;i++){
			require(physicalToLive.get(families[i].rank)<0,"Physical ranks must uniquely identify live families");
			physicalToLive.put(families[i].rank,i);
		}
		allHits=new SearchConfig(original.index,families,families.length,
			rawAdmission ? original.ddlReject : 1,1025,original.geometry,original.seed,1024,
			original.cosineMode,original.features,null,original.scalarBounds,false,false,original.identityFloor,
			original.blosumFloor,false,false,original.positionMode,original.partitions,original.overlap,original.blacklist,
			families[0].spectrum.length==FamilyReferenceCodec.MONOMER_WIDTH,original.memberScalarMeans);
	}

	public int familyCount(){return families.length;}
	public int familyRankAt(int ordinal){return families[ordinal].rank;}
	public String familyIdAt(int ordinal){return families[ordinal].id;}
	public int familyOrdinalAtRank(int rank){return physicalToLive.get(rank);}
	public Worker newWorker(){return new Worker();}
	public int armCount(){return activeArms.length;}
	public String armNameAt(int arm){
		require(arm>=0 && arm<activeArms.length,"Unknown active heuristic arm");
		return NAMES[activeArms[arm]];
	}
	public String[] armNames(){
		final String[] names=new String[activeArms.length];
		for(int i=0;i<names.length;i++){names[i]=armNameAt(i);}
		return names;
	}
	/** Every H mode requires member means and provided targets; higher weights never enable a new oracle. */
	public static boolean isHydrophobicityMode(String mode){
		return "hydro".equals(mode) || "hydrohigh".equals(mode) ||
			"hydroextra".equals(mode) || "hydro020".equals(mode) || isCosineMode(mode);
	}
	public static boolean isCosineMode(String mode){
		return "cosine".equals(mode) || "hydrocosine".equals(mode) || isPivotMode(mode);
	}
	/** The new coefficient comparison shares production raw-hit admission; old diagnostics stay reproducible. */
	public static boolean isPivotMode(String mode){return "pivotcompare".equals(mode);}
	public boolean usesRawAdmission(){return rawAdmission;}
	public int minimumHits(){return rawAdmission ? original.ddlReject : 1;}
	/** Exact active weight order for run provenance, including the unchanged zero-weight control. */
	public String hydrophobicityWeights(){
		require(familyHydrophobicity!=null,"Weight provenance requires a bound member-H comparison");
		final StringBuilder weights=new StringBuilder();
		for(int form:activeArms){
			if(weights.length()>0){weights.append(',');}
			if(form==6){
				weights.append('0');
			}else{
				weights.append(form<15 ? HYDRO_WEIGHTS[form-8] : COSINE_HYDRO_WEIGHTS[form-15]);
			}
		}
		return weights.toString();
	}
	/** One explicit cosine coefficient per active arm; the original contribution remains coefficient1. */
	public String cosineWeights(){
		require(cosineSweep,"Cosine coefficient provenance requires an explicit cosine-weight comparison");
		final StringBuilder weights=new StringBuilder();
		for(int form:activeArms){
			if(weights.length()>0){weights.append(',');}
			weights.append(form<15 ? 1.0 : COSINE_WEIGHTS[form-15]);
		}
		return weights.toString();
	}
	public static String armName(int arm){require(arm>=0 && arm<ALL_ARMS,"Unknown heuristic arm");return NAMES[arm];}
	private static int[] selectArms(String mode){
		if("six".equals(mode)){return new int[]{0,1,2,3,4,5};}
		if("eight".equals(mode)){return new int[]{0,1,2,3,4,5,6,7};}
		if("quarter".equals(mode)){return new int[]{6,7};}
		if("hydro".equals(mode)){return new int[]{6,8,9,10};}
		if("hydrohigh".equals(mode)){return new int[]{6,11,12,13};}
		if("hydroextra".equals(mode)){return new int[]{6,14};}
		if("hydro020".equals(mode)){return new int[]{6,10};}
		if("cosine".equals(mode)){return new int[]{10,15,16,17};}
		if("hydrocosine".equals(mode)){return new int[]{6,10,15,16,17};}
		if(isPivotMode(mode)){return new int[]{16,18,19};}
		throw new IllegalArgumentException("Rank arms must be six, eight, quarter, hydro, "+
			"hydrohigh, hydroextra, hydro020, cosine, hydrocosine or pivotcompare: "+mode);
	}
	public void verify()throws Exception{
		if(original.memberScalarMeans!=null){original.memberScalarMeans.verify();}
		require(pin!=null && pin.matches("[0-9a-f]{20}") && Files.isRegularFile(medians) && !Files.isSymbolicLink(medians) &&
			pin.equals(DigestSuffix.file(medians.toString())),"Member median binding changed: "+medians);
		if(hydroMeans!=null){
			require(hydroPin!=null && hydroPin.matches("[0-9a-f]{20}") &&
				Files.isRegularFile(hydroMeans) && !Files.isSymbolicLink(hydroMeans) &&
				hydroPin.equals(DigestSuffix.file(hydroMeans.toString())),"Member hydrophobicity binding changed: "+hydroMeans);
		}
	}

	/** The existing combined25 contribution is additive in raw-hit units. */
	public static double score(int arm,int hits,double cosine,int queryLength,double familyMedian){
		require(arm>=0 && arm<8 && hits>0 && hits<=1024 && Double.isFinite(cosine) && cosine>=0 && cosine<=1 &&
			queryLength>0 && Double.isFinite(familyMedian) && familyMedian>0,
			"Rank scores need finite hits, cosine and positive lengths");
		final double ratio=Math.min(queryLength,familyMedian)/Math.max(queryLength,familyMedian);
		switch(arm){
			case 0:return (hits+cosine)*ratio;
			case 1:return hits*ratio+cosine;
			case 2:return (hits+cosine)*Math.sqrt(ratio);
			case 3:return hits*Math.sqrt(ratio)+cosine;
			case 4:return hits;
			case 5:return FamilyTailCascade.combinedRankingScore(hits,cosine);
			case 6:return hits*Math.sqrt(Math.sqrt(ratio))+cosine;
			case 7:return (hits+cosine)*Math.sqrt(Math.sqrt(ratio));
			default:throw new IllegalArgumentException("Unknown heuristic score: "+arm);
		}
	}

	/** New arms add a bounded similarity to the selected fourth-root score; old forms stay exact. */
	public static double score(int arm,int hits,double cosine,int queryLength,double familyMedian,double hydroSimilarity){
		if(arm<8){return score(arm,hits,cosine,queryLength,familyMedian);}
		require(arm<ALL_ARMS && Double.isFinite(hydroSimilarity) && hydroSimilarity>=0 && hydroSimilarity<=1,
			"Hydrophobicity arms require a measured similarity in0..1, never an unavailable-value substitute");
		if(arm>=15){
			require(Double.isFinite(cosine) && cosine>=0 && cosine<=1,
				"Cosine multipliers require the measured combined25 cosine in0..1");
			return cosineScore(arm,score(6,hits,0,queryLength,familyMedian),cosine,hydroSimilarity);
		}
		return score(6,hits,cosine,queryLength,familyMedian)+HYDRO_WEIGHTS[arm-8]*hydroSimilarity;
	}
	/** Shared literal coefficient arithmetic for the public oracle and cached worker path. */
	private static double cosineScore(int arm,double lengthScore,double cosine,double hydroSimilarity){
		assert(arm>=15 && arm<ALL_ARMS && Double.isFinite(lengthScore) &&
			Double.isFinite(cosine) && Double.isFinite(hydroSimilarity)) :
			"Cosine comparisons use appended arms and already measured finite query/candidate coordinates";
		return lengthScore+COSINE_WEIGHTS[arm-15]*cosine+COSINE_HYDRO_WEIGHTS[arm-15]*hydroSimilarity;
	}
	/** ConsensusScalars.KD spans -4.5..4.5; dividing the mean difference by9 fixes the scale without fitting. */
	public static double hydrophobicitySimilarity(double queryMean,double familyMean){
		require(Double.isFinite(queryMean) && Double.isFinite(familyMean) &&
			queryMean>=-4.5 && queryMean<=4.5 && familyMean>=-4.5 && familyMean<=4.5,
			"Hydrophobicity similarity needs defined Kyte-Doolittle means within the native residue scale[-4.5,4.5]");
		return 1-Math.abs(queryMean-familyMean)/9;
	}
	/** Binds every live family's mean by physical rank, identity and retained count.
	 * Empty roster slots remain explicit. */
	private double[] readHydrophobicity(){
		final HashMap<String,Integer> ordinals=new HashMap<String,Integer>();
		for(int i=0;i<families.length;i++){
			require(ordinals.put(families[i].id,i)==null,
				"Live family identities must be unique for member-mean binding");
		}
		final double[] means=new double[families.length];
		final boolean[] seen=new boolean[families.length];
		final java.util.HashSet<String> ids=new java.util.HashSet<String>();
		final IntHashMap rows=new IntHashMap(families.length*2+1);
		final ByteFile input=ByteFile.makeByteFile1(hydroMeans.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		int count=0;
		try{
			require(Arrays.equals(input.nextLine(),HYDRO_HEADER),"Unexpected retained-member hydrophobicity table schema");
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line);
				require(lp.terms()==4,"Member hydrophobicity rows require rank, identity, count and mean");
				final int rank=lp.parseInt(0);
				final String id=lp.parseString(1);
				final long n=lp.parseLong(2);
				require(rank>=0 && n>=0 && !id.isEmpty() && ids.add(id) && rows.get(rank)<0,
					"Member mean rows need unique nonnegative ranks/identities and valid counts");
				rows.put(rank,0);
				final Integer at=ordinals.get(id);
				if(n==0){
					require(at==null && lp.termEquals("NA",3),"Only a non-live zero-member roster slot may have an undefined mean");
					continue;
				}
				require(at!=null && !seen[at] && families[at].rank==rank && families[at].members==n,
					"Member hydrophobicity rank, identity and count must match the exact live reference roster");
				means[at]=lp.parseDouble(3);
				hydrophobicitySimilarity(means[at],means[at]);
				seen[at]=true;
				count++;
			}
		}finally{require(!input.close(),"Member hydrophobicity table read failed");}
		require(count==families.length,"Member hydrophobicity table is incomplete; no consensus or zero fallback");
		return means;
	}

	/** Arrays are borrowed synchronously. Rank0 means unavailable, never first place.
	 * Empty oracle targets have ordinal-1/median0; a missing known target retains its real median.
	 */
	@FunctionalInterface
	public interface Observer{
		void accept(long queryOrdinal,String queryId,int queryLength,int targetOrdinal,boolean oracle,int targetHits,
			double targetMedian,boolean targetPassesDdlReject,boolean targetInRawWindow,int allHitCandidates,
			int ordinaryEligibleCandidates,int unionSize,long alignments,int[] ranks,int[] listSizes);
	}

	/** One instance belongs to one consuming thread; callbacks may not reenter it. */
	public final class Worker{
		private Worker(){
			search=new FamilyTailCascade.Searcher(allHits);
			seen=new long[families.length];
			hydroSimilarities=familyHydrophobicity==null ? null : new double[families.length];
			hydroBaseScores=familyHydrophobicity==null ? null : new double[families.length];
			cosineLengthScores=cosineSweep ? new double[families.length] : null;
		}
		public void rankObserve(long ordinal,String id,byte[] peptide,int truth,
				Observer observer){
			require(!busy && ordinal>=0 && id!=null && peptide!=null && observer!=null && truth>=-1 && truth<families.length,
				"Rank observation needs a stable ordinal, valid target and non-reentrant worker");
			require(range.contains(peptide.length),
				"Query length outside the declared rank cohort: "+peptide.length+" expected="+range);
			require(familyHydrophobicity==null || truth>=0,
				"Hydrophobicity follow-ons require a supplied frozen target; new oracle alignments are not released");
			busy=true;
			try{rank(ordinal,id,peptide,truth,observer);}finally{busy=false;}
		}
		private void rank(long ordinal,String id,byte[] raw,int truth,Observer observer){
			require(generation<Long.MAX_VALUE,"Rank union generation exhausted; stale dedup stamps must never be reused");
			generation++;
			final ProteinSequence query=new ProteinSequence(id,raw);
			require(query.length()==raw.length,"Rank queries must not be silently clipped or trimmed");
			search.prepare(raw,query,result,counters,null,false);
			if(hydroSimilarities!=null && !search.shortlist.isEmpty()){
				require(search.queryScalarsReady && (search.queryScalarMask&1)!=0,
					"Combined25 preparation must supply this query's measured hydrophobicity "+
					"before ranking its positive-hit families");
				for(Comparison candidate:search.shortlist){
					hydroSimilarities[candidate.ordinal]=
						hydrophobicitySimilarity(search.queryScalars[0],familyHydrophobicity[candidate.ordinal]);
					if(cosineLengthScores==null){
						hydroBaseScores[candidate.ordinal]=
							score(6,candidate.ddlHits,candidate.combinedCosine,raw.length,lengths[candidate.ordinal]);
					}else{
						require(candidate.combinedMeasured && Double.isFinite(candidate.combinedCosine) &&
							candidate.combinedCosine>=0 && candidate.combinedCosine<=1,
							"Cosine-weight arms must use measured combined25 scores from this query's all-hit preparation");
						final double lengthScore=score(6,candidate.ddlHits,0,raw.length,lengths[candidate.ordinal]);
						cosineLengthScores[candidate.ordinal]=lengthScore;
						hydroBaseScores[candidate.ordinal]=lengthScore+candidate.combinedCosine;
					}
				}
			}
			Arrays.fill(sizes,0);Arrays.fill(ranks,0);
			int eligible=0,cutoff=original.ddlReject;
			Comparison target=null;
			for(Comparison candidate:search.shortlist){
				require(candidate.combinedMeasured && candidate.ddlHits>0,
					"All-hit preparation must provide combined25 for each positive-hit candidate");
				if(candidate.ddlHits>=original.ddlReject){eligible++;}
				if(candidate.ddlPosition==50){cutoff=Math.max(cutoff,candidate.ddlHits);}
				if(candidate.ordinal==truth){target=candidate;}
			}
			final boolean oracle=truth<0;
			for(Comparison candidate:search.shortlist){
				for(int arm=0;arm<activeArms.length;arm++){
					if(activeArms[arm]==CURRENT && candidate.ddlHits<cutoff){continue;}
					if(oracle){offer(arm,candidate,raw.length);}
					else if(sizes[arm]<TOP){sizes[arm]++;}
				}
			}
			int unionSize=0,targetOrdinal=truth;
			if(oracle){
				for(int arm=0;arm<activeArms.length;arm++){
					for(int i=0;i<sizes[arm];i++){
						final int family=tops[arm][i];
						if(seen[family]!=generation){
							require(unionSize<union.length,"The active strict top50 lists cannot exceed their total candidate capacity");
							seen[family]=generation;union[unionSize++]=family;
						}
					}
				}
				final int self=ProteinSearcher.selfScore(query);double best=0;
				for(int i=0;i<unionSize;i++){
					final int family=union[i];
					FamilyCalibration.align(query,families[family].consensus,self,metrics);
					require(Double.isFinite(metrics.blosum),"Oracle alignment must return a finite full-query score");
					if(targetOrdinal<0 || metrics.blosum>best ||
							metrics.blosum==best && families[family].rank<families[targetOrdinal].rank){
						targetOrdinal=family;best=metrics.blosum;
					}
				}
				if(targetOrdinal>=0){target=search.comparisons[targetOrdinal];}
			}
			final int hits=target==null ? 0 : target.ddlHits;
			final boolean window=hits>0 && hits>=cutoff;
			if(target!=null){
				for(int arm=0;arm<activeArms.length;arm++){
					if(activeArms[arm]==CURRENT && !window){continue;}
					ranks[arm]=1;targetScores[arm]=value(arm,target,raw.length);
				}
				for(Comparison candidate:search.shortlist){
					if(candidate.ordinal==targetOrdinal){continue;}
					for(int arm=0;arm<activeArms.length;arm++){
						if(ranks[arm]==0 || activeArms[arm]==CURRENT && candidate.ddlHits<cutoff){continue;}
						if(better(value(arm,candidate,raw.length),candidate.ordinal,
								targetScores[arm],targetOrdinal)){ranks[arm]++;}
					}
				}
			}
			observer.accept(ordinal,id,raw.length,targetOrdinal,oracle,hits,
				targetOrdinal<0 ? 0 : lengths[targetOrdinal],hits>0 && hits>=original.ddlReject,window,
				search.shortlist.size(),eligible,unionSize,oracle ? unionSize : 0,ranks,sizes);
		}
		private double value(int arm,Comparison candidate,int length){
			final int form=activeArms[arm];
			if(hydroBaseScores!=null){
				final double base=hydroBaseScores[candidate.ordinal];
				if(form>=15){
					return cosineScore(form,cosineLengthScores[candidate.ordinal],
						candidate.combinedCosine,hydroSimilarities[candidate.ordinal]);
				}
				return form==6 ? base : base+HYDRO_WEIGHTS[form-8]*hydroSimilarities[candidate.ordinal];
			}
			return score(form,candidate.ddlHits,candidate.combinedCosine,length,lengths[candidate.ordinal]);
		}
		/** Worst-first fixed heap: O(log50) insertion and no per-candidate allocation. */
		private void offer(int arm,Comparison candidate,int length){
			final int family=candidate.ordinal;final double value=value(arm,candidate,length);
			final int[] ids=tops[arm];final double[] values=topScores[arm];
			int size=sizes[arm];
			if(size<TOP){
				int child=size++;sizes[arm]=size;
				while(child>0){
					final int parent=(child-1)>>>1;
					if(!better(values[parent],ids[parent],value,family)){break;}
					ids[child]=ids[parent];values[child]=values[parent];child=parent;
				}
				ids[child]=family;values[child]=value;return;
			}
			if(!better(value,family,values[0],ids[0])){return;}
			int parent=0;
			while(parent*2+1<size){
				int child=parent*2+1;
				if(child+1<size && better(values[child],ids[child],values[child+1],ids[child+1])){child++;}
				if(!better(value,family,values[child],ids[child])){break;}
				ids[parent]=ids[child];values[parent]=values[child];parent=child;
			}
			ids[parent]=family;values[parent]=value;
		}
		private boolean better(double a,int ai,double b,int bi){return a>b || a==b && families[ai].rank<families[bi].rank;}
		private final FamilyTailCascade.Searcher search;
		private final FamilyTailCascade.Result result=new FamilyTailCascade.Result();
		private final FamilyTailCascade.Counters counters=new FamilyTailCascade.Counters();
		private final FamilyCalibration.Metrics metrics=new FamilyCalibration.Metrics();
		private final int[][] tops=new int[activeArms.length][TOP];
		private final double[][] topScores=new double[activeArms.length][TOP];
		private final int[] sizes=new int[activeArms.length],ranks=new int[activeArms.length],
			union=new int[activeArms.length*TOP];
		private final double[] targetScores=new double[activeArms.length];
		private final long[] seen;
		/** Caches are populated only for this query's positive-hit candidates before any arm reads them. */
		private final double[] hydroSimilarities,hydroBaseScores,cosineLengthScores;
		private long generation;private boolean busy;
	}

	private double[] readMedians(){
		final HashMap<String,Integer> ordinals=new HashMap<String,Integer>();
		for(int i=0;i<families.length;i++){
			require(ordinals.put(families[i].id,i)==null,"Live family IDs must be unique for median binding");
		}
		final double[] values=new double[families.length];
		final boolean[] seen=new boolean[values.length];
		int count=0;
		final ByteFile input=ByteFile.makeByteFile1(medians.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			require(Arrays.equals(input.nextLine(),MEDIAN_HEADER),"Unexpected member median table schema");
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line);
				require(lp.terms()==6,"Member median rows need six fields");
				final Integer ordinal=ordinals.get(lp.parseString(0));
				require(ordinal!=null && !seen[ordinal] && lp.parseLong(1)==families[ordinal].members,
					"Median IDs/counts must match the exact live family roster once");
				final double value=lp.parseDouble(3);
				require(Double.isFinite(value) && value>=40 && value<=200 && value*2==Math.rint(value*2),
					"Member median must be a finite integer or half integer in40..200; no imputation");
				values[ordinal]=value;seen[ordinal]=true;count++;
			}
		}finally{require(!input.close(),"Member median table read failed");}
		require(count==families.length,"Member median table is incomplete");
		return values;
	}
	private static void require(boolean valid,String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}
	private final SearchConfig original,allHits;
	private final FamilyReference[] families;
	private final Path medians,hydroMeans;
	private final String pin,hydroPin;
	private final double[] lengths,familyHydrophobicity;
	private final MemberLengthRange range;private final IntHashMap physicalToLive;
	private final int[] activeArms;
	private final boolean cosineSweep,rawAdmission;
	/** Historical default count; explicit eight/quarter instances use armCount(). */
	public static final int ARMS=6,ALL_ARMS=20;
	private static final int TOP=50,CURRENT=5;
	private static final double[] HYDRO_WEIGHTS={.05,.1,.2,.4,.8,1.6,3.2};
	private static final double[] COSINE_WEIGHTS={.5,1.5,2,1.25,1.5};
	private static final double[] COSINE_HYDRO_WEIGHTS={.2,.2,.2,.1,.1};
	private static final String[] NAMES={"sumratio","hitratio","sumsqrt","hitsqrt","raw","current","hitfourth","sumfourth",
		"hitfourth_h005","hitfourth_h010","hitfourth_h020",
		"hitfourth_h040","hitfourth_h080","hitfourth_h160","hitfourth_h320",
		"hitfourth_h020_c050","hitfourth_h020_c150","hitfourth_h020_c200",
		"hitfourth_h010_c125","hitfourth_h010_c150"};
	private static final byte[] MEDIAN_HEADER=
		"family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] HYDRO_HEADER=
		"family_rank\tfamily_id\tN\thydrophobicity_mean".getBytes(StandardCharsets.US_ASCII);
}
