package prot;

import java.util.Arrays;

/** Experimental stop-tail search over existing compatible family references.
 * EXHAUSTIVE preserves the original full50 search; FAST preselects before top50
 * and bounds rank-order search after strict supported score improvements.
 * FAST competitors describe visited families only, not global confidence.
 * @author Keqing, Sayu
 */
public final class ProteinTailAssigner{
	public ProteinTailAssigner(final HbmTailGate.Policy policy_, final ProteinFamilyAssigner... packs_){
		this(policy_, new SearchPolicy(Mode.EXHAUSTIVE, 0, 0, 0), packs_);
	}
	public ProteinTailAssigner(final HbmTailGate.Policy policy_, final SearchPolicy search_, final ProteinFamilyAssigner... packs_){
		if(policy_==null || search_==null || packs_==null){throw new IllegalArgumentException("Tail search needs explicit acceptance/search policies and reference packs");}
		packs=packs_.clone(); policy=policy_; search=search_;
		// Reuse existing cross-pack identity/background validation before sharing scores.
		new ProteinFamilyCompetition(packs);
	}
	public Worker newWorker(){return new Worker();}
	public String rule(){return policy.profileEvidence ? "hbm_stop_tail_profile_joint_v1" : RULE;}
	public String competitorScope(){return search.mode==Mode.FAST ? "visited" : "full_top50";}
	public enum Mode{EXHAUSTIVE, FAST}
	/** Seed thresholds are experimental global filters, not per-family calibration. */
	public static final class SearchPolicy{
		public SearchPolicy(Mode mode_, int minSeeds_, double minSeedDensity_, int lookahead_){
			if(mode_==null || minSeeds_<0 || !Double.isFinite(minSeedDensity_) || minSeedDensity_<0 || minSeedDensity_>1 || lookahead_<0){
				throw new IllegalArgumentException("Invalid tail search policy");
			}
			if(mode_==Mode.EXHAUSTIVE && (minSeeds_!=0 || minSeedDensity_!=0 || lookahead_!=0)){
				throw new IllegalArgumentException("EXHAUSTIVE requires original seed0/density0/lookahead0");
			}
			mode=mode_; minSeeds=minSeeds_; minSeedDensity=minSeedDensity_; lookahead=lookahead_;
		}
		public final Mode mode;
		public final int minSeeds, lookahead;
		public final double minSeedDensity;
	}

	/** One worker owns all query scratch; no synchronization or cross-query evidence sharing. */
	public final class Worker{
		private Worker(){
			final int n=packs.length;
			shortlists=new ProteinSearcher.ShortlistScratch[n]; hits=new int[n][];
			composition=new float[n][]; scores=new float[n][];
			medians=new double[n][]; spreads=new double[n][];
			final HbmCompetitiveGate[] gates=new HbmCompetitiveGate[n];
			final String[][] ids=new String[n][]; final int[][] familyIds=new int[n][];
			qSizes=new int[n];
			for(int p=0; p<n; p++){
				final FamilyShortlistSidecar sidecar=packs[p].library().sidecar;
				final ProteinSearcher.ShortlistScratch s=new ProteinSearcher.ShortlistScratch(sidecar.dims, sidecar.nFamilies, TOP);
				shortlists[p]=s; hits[p]=s.sk; composition[p]=s.compScores; scores[p]=s.f4;
				medians[p]=sidecar.medLen; spreads[p]=sidecar.sdLnLen;
				gates[p]=packs[p].library().gate; ids[p]=packs[p].library().ids; familyIds[p]=packs[p].library().familyIds;
			}
			work=new SearchScratch(policy, search, gates, ids, familyIds);
		}

		/** Input is bounded translated ASCII amino acids; no initiator-M rewrite. */
		public Hit assign(final String id, final byte[] aminoAcids){
			final ProteinSequence query=new ProteinSequence(id, aminoAcids);
			if(query.length()>policy.maxLength){throw new IllegalArgumentException("Stop query exceeds tail cap: "+query.length());}
			if(search.mode==Mode.FAST){
				//All-short queries fail length first, so no seed feature is read.
				if(query.length()>=policy.minPaired){
					for(int p=0; p<packs.length; p++){
						ProteinSearcher.scoreFamiliesKmers(query, packs[p].library().sidecar, shortlists[p]);
						qSizes[p]=shortlists[p].qSize;
					}
				}
				if(!work.hasCandidate(query.enc, hits, qSizes)){return work.emptyHit();}
			}
			for(int p=0; p<packs.length; p++){
				if(search.mode==Mode.FAST){ProteinSearcher.scoreFamiliesF4FromSeeds(query, packs[p].library().sidecar, shortlists[p]);}
				else{ProteinSearcher.scoreFamiliesF4(query, packs[p].library().sidecar, shortlists[p]);}
				qSizes[p]=shortlists[p].qSize;
			}
			FamilyShortlistScores.scoreTail(hits, composition, medians, spreads, query.length(), scores);
			return work.assign(query.enc, scores, hits, qSizes);
		}
		private final ProteinSearcher.ShortlistScratch[] shortlists;
		private final int[][] hits;
		private final int[] qSizes;
		private final float[][] composition, scores;
		private final double[][] medians, spreads;
		private final SearchScratch work;
	}

	/** Production selection scratch, also used by small bound-profile fixtures. */
	public static final class SearchScratch{
		public SearchScratch(HbmTailGate.Policy policy_, SearchPolicy search_, HbmCompetitiveGate[] gates_, String[][] ids_, int[][] families_){
			assert(policy_!=null && search_!=null && gates_.length==ids_.length && ids_.length==families_.length) : "Tail scratch must bind the same packs as retrieval";
			policy=policy_; search=search_; gates=gates_; ids=ids_; families=families_;
		}
		public Hit assign(byte[] query, float[][] scores, int[][] seeds, int[] qSizes){
			assert(query!=null && query.length>0 && query.length<=policy.maxLength) : "Tail scratch needs an encoded bounded query";
			Arrays.fill(rejected, 0); winner.clear();
			return search.mode==Mode.EXHAUSTIVE ? exhaustive(query, scores) : fast(query, scores, seeds, qSizes);
		}
		/** No composition/ranking reads. Only a complete zero-survivor scan has usable rejection totals. */
		public boolean hasCandidate(byte[] query, int[][] seeds, int[] qSizes){
			assert(search.mode==Mode.FAST && query.length>0 && query.length<=policy.maxLength) : "Only FAST uses the bounded pre-composition eligibility scan";
			emptyLengthRejected=emptySeedRejected=0;
			for(int p=0; p<gates.length; p++){
				for(int r=0; r<ids[p].length; r++){
					if(!lengthCanPass(query.length, p, r)){emptyLengthRejected++; continue;}
					if(!seedsCanPass(seeds[p][r], qSizes[p])){emptySeedRejected++; continue;}
					return true;
				}
			}
			return false;
		}
		/** Called only by Worker after a complete empty scan and no composition work. */
		private Hit emptyHit(){
			Arrays.fill(rejected, 0); winner.clear();
			return new Hit(-1, -1, -1, null, 0, -1, -1, null, 0, rejected,
				0, emptyLengthRejected, emptySeedRejected, 0, 0, 0, 0, true);
		}
		private boolean lengthCanPass(int length, int p, int r){
			return length>=policy.minPaired && gates[p].tailLengthCanPass(length, r, policy);
		}
		private boolean seedsCanPass(int seeds, int qSize){
			assert(seeds>=0 && qSize>=0) : "Seed eligibility needs current nonnegative covering counts and distinct-query denominator";
			return seeds>=search.minSeeds && (search.minSeedDensity<=0 || (qSize>0 && seeds/(double)qSize>=search.minSeedDensity));
		}
		private Hit exhaustive(byte[] query, float[][] scores){
			FamilyShortlistScores.top(scores, topPacks, topRanks);
			int count=0;
			for(int k=0; k<TOP; k++){
				final int p=topPacks[k], rank=topRanks[k];
				final HbmPositionModel.TailResult a=gates[p].alignTail(query, rank, policy.maxLength, false);
				if(a==null){rejected[HbmTailGate.NO_ALIGNMENT]++; continue;}
				int at=count;
				while(at>0 && HbmCompetitiveGate.better(a.score, ids[p][rank],
					orderedScores[at-1], ids[orderedPacks[at-1]][orderedRanks[at-1]])){
					orderedScores[at]=orderedScores[at-1]; orderedPacks[at]=orderedPacks[at-1]; orderedRanks[at]=orderedRanks[at-1]; at--;
				}
				orderedScores[at]=a.score; orderedPacks[at]=p; orderedRanks[at]=rank; count++;
			}
			int best=-1, runner=-1, recorded=0;
			for(int k=0; k<count; k++){
				final int p=orderedPacks[k], rank=orderedRanks[k];
				final HbmCompetitiveGate gate=gates[p];
				final HbmPositionModel.TailResult a=gate.alignTail(query, rank, policy.maxLength, true); recorded++;
				assert(a!=null && a.score==orderedScores[k]) : "Recording a tail path must preserve its optimum score";
				final int status=gate.evaluateTail(query, rank, a, policy, scratch);
				if(status!=HbmTailGate.SUPPORTED){rejected[status]++; continue;}
				if(best<0){best=k; winner.copyFrom(scratch);}
				else{runner=k; break;}
			}
			final int p=best<0 ? -1 : orderedPacks[best], rank=best<0 ? -1 : orderedRanks[best];
			final int other=best<0 || count<2 ? -1 : orderedScores[best==0 ? 1 : 0];
			return new Hit(p, rank, best<0 ? -1 : families[p][rank], best<0 ? null : ids[p][rank],
				best<0 ? 0 : orderedScores[best], runner<0 ? -1 : orderedScores[runner], other, best<0 ? null : winner, recorded, rejected,
				TOP, 0, 0, TOP, TOP+recorded, best<0 ? 0 : runner<0 ? 1 : 2, 0);
		}
		private Hit fast(byte[] query, float[][] scores, int[][] seeds, int[] qSizes){
			int count=0, lengthRejected=0, seedRejected=0;
			for(int p=0; p<gates.length; p++){
				for(int r=0; r<ids[p].length; r++){
					if(!lengthCanPass(query.length, p, r)){lengthRejected++; continue;}
					assert(seeds[p][r]>=0 && qSizes[p]>=0 && Float.isFinite(scores[p][r])) : "F4 and query seed features must be valid before selection";
					if(!seedsCanPass(seeds[p][r], qSizes[p])){seedRejected++; continue;}
					final float score=scores[p][r];
					if(count==TOP && score<=scores[topPacks[count-1]][topRanks[count-1]]){continue;}
					int at=Math.min(count, TOP-1);
					while(at>0 && score>scores[topPacks[at-1]][topRanks[at-1]]){at--;}
					if(count<TOP){count++;}
					for(int k=count-1; k>at; k--){topPacks[k]=topPacks[k-1]; topRanks[k]=topRanks[k-1];}
					topPacks[at]=p; topRanks[at]=r;
				}
			}
			int best=-1, bestScore=0, runner=-1, positive=-1, positiveScore=-1, secondPositive=-1;
			int visited=0, recorded=0, accepted=0, stop=count;
			for(int k=0; k<stop; k++){
				final int p=topPacks[k], r=topRanks[k]; visited++;
				final HbmPositionModel.TailResult a=gates[p].alignTail(query, r, policy.maxLength, true);
				if(a==null){rejected[HbmTailGate.NO_ALIGNMENT]++; continue;}
				recorded++;
				if(a.score>positiveScore){secondPositive=positiveScore; positiveScore=a.score; positive=k;}
				else if(a.score>secondPositive){secondPositive=a.score;}
				final int status=gates[p].evaluateTail(query, r, a, policy, scratch);
				if(status!=HbmTailGate.SUPPORTED){rejected[status]++; continue;}
				accepted++;
				if(best<0 || HbmCompetitiveGate.better(a.score, ids[p][r], bestScore, ids[topPacks[best]][topRanks[best]])){
					if(a.score>bestScore){stop=k+1+Math.min(search.lookahead, count-k-1);}
					runner=best<0 ? -1 : bestScore; best=k; bestScore=a.score; winner.copyFrom(scratch);
				}else if(a.score>runner){runner=a.score;}
			}
			int rejects=0; for(int n : rejected){rejects+=n;}
			assert(visited==accepted+rejects && recorded<=visited) : "Every visited tail must be supported or have exactly one rejection";
			final int p=best<0 ? -1 : topPacks[best], r=best<0 ? -1 : topRanks[best];
			final int other=best<0 ? -1 : best==positive ? secondPositive : positiveScore;
			return new Hit(p, r, best<0 ? -1 : families[p][r], best<0 ? null : ids[p][r], bestScore, runner, other,
				best<0 ? null : winner, recorded, rejected, count, lengthRejected, seedRejected, visited, visited, accepted, count-visited);
		}
		private final HbmTailGate.Policy policy;
		private final SearchPolicy search;
		private final HbmCompetitiveGate[] gates;
		private final String[][] ids;
		private final int[][] families;
		private int emptyLengthRejected, emptySeedRejected;
		private final int[] topPacks=new int[TOP], topRanks=new int[TOP];
		private final int[] orderedPacks=new int[TOP], orderedRanks=new int[TOP], orderedScores=new int[TOP];
		private final int[] rejected=new int[HbmTailGate.STATUS_COUNT];
		private final HbmTailGate.Metrics scratch=new HbmTailGate.Metrics(), winner=new HbmTailGate.Metrics();
	}

	/** Detached evidence. Missing competitors are -1, not infinite confidence margins. */
	public static final class Hit{
		public Hit(int p, int r, int f, String id, int score, int runner, int other, HbmTailGate.Metrics m, int recorded_, int[] rejected_){
			this(p, r, f, id, score, runner, other, m, recorded_, rejected_, TOP, 0, 0, TOP, TOP+recorded_, r<0 ? 0 : runner<0 ? 1 : 2, 0);
		}
		Hit(int p, int r, int f, String id, int score, int runner, int other, HbmTailGate.Metrics m, int recorded_, int[] rejected_,
				int shortlisted_, int lengthRejected_, int seedRejected_, int visited_, int scored_, int accepted_, int unvisited_){
			this(p, r, f, id, score, runner, other, m, recorded_, rejected_, shortlisted_, lengthRejected_, seedRejected_,
				visited_, scored_, accepted_, unvisited_, false);
		}
		private Hit(int p, int r, int f, String id, int score, int runner, int other, HbmTailGate.Metrics m, int recorded_, int[] rejected_,
				int shortlisted_, int lengthRejected_, int seedRejected_, int visited_, int scored_, int accepted_, int unvisited_, boolean compositionSkipped_){
			compositionSkipped=compositionSkipped_;
			assert(!compositionSkipped || (r<0 && shortlisted_==0 && scored_==0)) : "Composition bypass is exact only for zero eligible families and zero DP work";
			shortlisted=shortlisted_; lengthRejected=lengthRejected_; seedRejected=seedRejected_; visited=visited_;
			scored=scored_; accepted=accepted_; unvisited=unvisited_; earlyStopped=unvisited_>0;
			packIndex=p; rank=r; familyId=f; representative=id; score64=score; runnerUpScore64=runner; bestOtherScore64=other; recorded=recorded_;
			paired=m==null ? 0 : m.paired; pairedCore=m==null ? 0 : m.pairedCore; terminalPaired=m==null ? 0 : m.terminalPaired;
			terminalScore64=m==null ? Integer.MIN_VALUE : m.terminalScore64;
			queryStart=m==null ? -1 : m.queryStart; referenceStart=m==null ? -1 : m.referenceStart;
			querySuffix=m==null ? -1 : m.querySuffix; referenceSuffix=m==null ? -1 : m.referenceSuffix;
			identity=m==null ? Float.NaN : m.identity; pairedFraction=m==null ? Float.NaN : m.pairedFraction;
			coreCoverage=m==null ? Float.NaN : m.coreCoverage; rejected=rejected_.clone();
		}
		public boolean assigned(){return rank>=0;}
		public final int packIndex, rank, familyId, score64, recorded;
		/** scored counts every DP fill, including EXHAUSTIVE recorded reruns; recorded counts nonnull recorded paths. */
		public final int shortlisted, lengthRejected, seedRejected, visited, scored, accepted, unvisited;
		public final boolean earlyStopped;
		/** True only when Worker returned before composition because no family could pass the existing filters. */
		public final boolean compositionSkipped;
		/** Second supported hit versus strongest other positive score, even when it fails the gate. */
		public final int runnerUpScore64, bestOtherScore64;
		public final int paired, pairedCore, terminalPaired, queryStart, referenceStart, querySuffix, referenceSuffix;
		public final int terminalScore64;
		public final String representative;
		public final float identity, pairedFraction, coreCoverage;
		/** Rejections among scored-null or traced candidates; untraced lower scores are not classified. */
		public final int[] rejected;
	}
	public static final String RULE="hbm_stop_tail_experimental_v1";
	private static final int TOP=50;
	private final ProteinFamilyAssigner[] packs;
	private final HbmTailGate.Policy policy;
	private final SearchPolicy search;
}
