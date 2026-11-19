package prot;

import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import map.LongHashSet;
import map.LongLongListHashMap;
import structures.IntList;
import structures.LongList;

/** Shared stop-tail detector used by the benchmark and opt-in callers.
 * Ties keep deterministic shortlist order. Returned Result/Hit instances copy
 * every scalar they expose, so later Scratch reuse cannot mutate decisions.
 * Stop-tail representatives are evidence for terminal stops, not functional labels.
 * @author Yelan
 */
public final class StopTailDetector {

	/** Loads a pinned legacy k-mer shortlist or exhaustive detector reference. */
	public static Reference load(String ref,String refPin,int shortlist,int k,HbmInput hbm){
		verify(ref, refPin, hbm);
		final Reference out=new Reference(ref, shortlist, k, null, hbm);
		verify(ref, refPin, hbm);
		return out;
	}

	/** Loads a pinned paired DDL/covering detector reference. */
	public static Reference load(String ref,String refPin,PairedSettings paired,HbmInput hbm) throws Exception{
		verify(ref, refPin, hbm);
		paired.verify();
		final Reference out=new Reference(ref, 0, 0, paired.load(ref), hbm);
		verify(ref, refPin, hbm);
		paired.verify();
		return out;
	}

	/** Loads a detector from an already-instantiated paired ranker in the prot package. */
	static Reference load(String ref,int shortlist,int k,HbmTailPairedShortlist paired,HbmInput hbm){
		return new Reference(ref, shortlist, k, paired, hbm);
	}

	public enum Method {
		AAI("aai"), BLOSUM("blosum"), HBM("hbm");
		Method(String name_){name=name_;}
		public final String name;
	}

	public static final class Policy {
		public Policy(double coverage_,int minPaired_,int terminalPaired_,boolean endPaired_){
			require(Double.isFinite(coverage_) && coverage_>0 && coverage_<=1 && minPaired_>0
				&& terminalPaired_>0, "Invalid geometry policy");
			coverage=coverage_; minPaired=minPaired_; terminalPaired=terminalPaired_; endPaired=endPaired_;
		}
		boolean reject(int qLen,int rLen,HbmPositionModel.TailResult tr,HbmTailEvidence.Metrics m){
			return (endPaired && tr.path[tr.path.length-1]!='m') ||
				m.paired<Math.max(minPaired, (int)Math.ceil(coverage*Math.max(qLen, rLen))) ||
				m.terminalPaired<terminalPaired;
		}
		final double coverage;
		final int minPaired, terminalPaired;
		final boolean endPaired;
	}

	/** Path/pin bundle for the optional stop-tail HBM scorer. */
	public static final class HbmInput {
		static HbmInput parse(HashMap<String,String> opt){
			final String[] keys={"hbm", "hbmsha80", "provenance", "provenancesha80", "background", "backgroundsha80", "hbmcutoffs"};
			int present=0;
			for(String key : keys){if(opt.containsKey(key)){present++;}}
			require(present==0 || present==keys.length,
				"HBM scoring requires hbm/hbmsha80/provenance/provenancesha80/background/backgroundsha80/hbmcutoffs");
			return present==0 ? null : new HbmInput(opt.get("hbm"), opt.get("hbmsha80"), opt.get("provenance"),
				opt.get("provenancesha80"), opt.get("background"), opt.get("backgroundsha80"));
		}
		public HbmInput(String hbm_,String hbmPin_,String provenance_,String provenancePin_,String background_,String backgroundPin_){
			hbm=hbm_; hbmPin=hbmPin_; provenance=provenance_; provenancePin=provenancePin_;
			background=background_; backgroundPin=backgroundPin_;
		}
		void verify(){
			HbmProfilePilot.requireHash(hbm, hbmPin);
			HbmProfilePilot.requireHash(provenance, provenancePin);
			HbmProfilePilot.requireHash(background, backgroundPin);
		}
		final String hbm, hbmPin, provenance, provenancePin, background, backgroundPin;
	}

	/** Public pinned file-based settings for paired DDL plus covering shortlists. */
	public static final class PairedSettings {
		public PairedSettings(String profilesList_,String profilesListPin_,String sidecar_,String sidecarPin_,
				String kmers_,String kmersPin_,int buckets_,String mode_,int k_,long seed_,boolean indexed_,
				String namespace_,int ddlLimit_,int coverLimit_){
			profilesList=profilesList_; profilesListPin=profilesListPin_; sidecar=sidecar_;
			sidecarPin=sidecarPin_; kmers=kmers_; kmersPin=kmersPin_;
			buckets=buckets_; mode=mode_; k=k_; seed=seed_; indexed=indexed_; namespace=namespace_;
			ddlLimit=ddlLimit_; coverLimit=coverLimit_;
		}
		HbmTailPairedShortlist load(String ref) throws Exception{
			return HbmTailPairedShortlist.load(profilesList, sidecar, ref, kmers, buckets, mode, k, seed,
				indexed, namespace, ddlLimit, coverLimit);
		}
		void verify(){
			HbmProfilePilot.requireHash(profilesList, profilesListPin);
			HbmProfilePilot.requireHash(sidecar, sidecarPin);
			HbmProfilePilot.requireHash(kmers, kmersPin);
		}
		private final String profilesList, profilesListPin, sidecar, sidecarPin, kmers, kmersPin, mode, namespace;
		private final int buckets, k, ddlLimit, coverLimit;
		private final long seed;
		private final boolean indexed;
	}

	/** Raw ASCII plus encoded amino acids for a single candidate stop-tail. */
	public static final class Query {
		public Query(String id_,byte[] raw_){
			require(raw_!=null, "Stop-tail query bases may not be null: "+id_);
			raw=Arrays.copyOf(raw_, raw_.length);
			protein=new ProteinSequence(id_, raw);
			id=protein.id;
			require(protein.length()==raw.length && protein.length()>=100 && protein.length()<=200,
				"Expected unmodified 100..200-AA tail: "+id);
			for(byte b : protein.enc){require(b>=0 && b<20, "Stop-tail queries reject ambiguity and stop codons: "+id);}
		}
		public int length(){return raw.length;}
		public final String id;
		final byte[] raw;
		final ProteinSequence protein;
	}

	/** Immutable detector reference. Scratch and Result objects hold all mutable state. */
	public static final class Reference {
		private Reference(String ref,int shortlist,int k_,HbmTailPairedShortlist paired_,HbmInput hbm){
			require(paired_!=null || shortlist>=0 && k_>=0 && (shortlist==0 || k_>0), "Invalid stop-tail shortlist");
			paired=paired_;
			final java.util.List<ProteinSequence> seqs=HbmTailLibrary.readTails(ref);
			require(!seqs.isEmpty(), "Empty stop-tail reference");
			ids=new String[seqs.size()]; refs=new byte[seqs.size()][]; profiles=new HbmPositionModel[seqs.size()];
			index=paired==null && shortlist>0 ? new LongLongListHashMap(seqs.size()*16+1) : null;
			exhaustive=paired==null && shortlist==0;
			alphabet=paired==null && shortlist>0 ? ReducedAlphabetSeedAssay.Alphabet.named("aa20") : null;
			k=k_;
			final HashSet<String> seen=new HashSet<String>();
			for(int i=0; i<seqs.size(); i++){
				final ProteinSequence s=seqs.get(i);
				require(s.length()>=100 && s.length()<=200 && safe(s.id) && seen.add(s.id), "Invalid consensus tail: "+s.id);
				for(byte b : s.enc){require(b>=0 && b<20, "Consensus tails reject ambiguity and stops: "+s.id);}
				ids[i]=s.id; refs[i]=s.enc; profiles[i]=HbmPositionModel.fromConsensus(s.enc);
				if(index!=null){
					final LongHashSet words=ReducedAlphabetSeedAssay.kmerSet(s.enc, alphabet, k, ReducedAlphabetSeedAssay.Boundary.NONE);
					for(long word : words.toArray()){index.put(word, i);}
				}
			}
			if(paired!=null){paired.validateRoster(ids);}
			hbmProfiles=hbm==null ? null : loadHbm(hbm);
			rankLimit=paired==null ? (shortlist==0 ? ids.length : Math.min(shortlist, ids.length)) : paired.limit();
		}

		public Scratch makeScratch(){return new Scratch(this);}
		public int rankLimit(){return rankLimit;}

		public Result assign(Query query,Policy policy,int minHits,Scratch scratch){
			require(query!=null && policy!=null && scratch!=null, "Stop-tail assignment requires query, policy and scratch");
			require(scratch.owner==this, "Stop-tail scratch belongs to a different reference");
			require(minHits>=0, "minkmerhits must be nonnegative");
			require(!exhaustive || minHits==0, "Exhaustive stop-tail assignment requires minkmerhits=0");
			final int candidates=shortlist(query, scratch.ranks, minHits, scratch);
			scratch.clear();
			final byte[] queryEnc=query.protein.enc;
			int blosumScored=0, blosumTraces=0, blosumGeometryFail=0, blosumAccepted=0;
			int hbmScored=0, hbmTraces=0, hbmGeometryFail=0, hbmAccepted=0;
			for(int i=0; i<candidates; i++){
				final int rank=scratch.ranks[i];
				blosumScored++;
				final HbmPositionModel.TailResult tr=profiles[rank].alignTail(queryEnc, HbmTailProfileRefiner.CAP, true);
				if(tr!=null){
					blosumTraces++;
					final HbmTailEvidence.Metrics m=HbmTailEvidence.measure(queryEnc, refs[rank], profiles[rank], tr);
					if(policy.reject(query.length(), refs[rank].length, tr, m)){blosumGeometryFail++;}
					else{
						blosumAccepted++;
						scratch.aai.update(Method.AAI, m.matches/(double)tr.path.length, tr.score, tr.path.length, rank, ids[rank], tr, m);
						scratch.blosum.update(Method.BLOSUM, tr.score/(64.0*tr.path.length), tr.score, tr.path.length, rank, ids[rank], tr, m);
					}
				}
				if(hbmProfiles!=null){
					hbmScored++;
					final HbmPositionModel.TailResult htr=hbmProfiles[rank].alignTail(queryEnc, HbmTailProfileRefiner.CAP, true);
					if(htr!=null){
						hbmTraces++;
						final HbmTailEvidence.Metrics hm=HbmTailEvidence.measure(queryEnc, refs[rank], hbmProfiles[rank], htr);
						if(policy.reject(query.length(), refs[rank].length, htr, hm)){hbmGeometryFail++;}
						else{
							hbmAccepted++;
							scratch.hbm.update(Method.HBM, htr.score/(64.0*htr.path.length), htr.score,
								htr.path.length, rank, ids[rank], htr, hm);
						}
					}
				}
			}
			return new Result(candidates, noKmers(scratch), scratch.aai, scratch.blosum,
				hbmProfiles==null ? null : scratch.hbm,
				new Work(blosumScored, blosumTraces, blosumGeometryFail, blosumAccepted),
				new Work(hbmScored, hbmTraces, hbmGeometryFail, hbmAccepted));
		}

		HbmPositionModel[] loadHbm(HbmInput hbm){
			try{
				final HashMap<String,byte[]> byId=new HashMap<String,byte[]>();
				for(int i=0; i<ids.length; i++){byId.put(ids[i], refs[i]);}
				final byte[][] provenance=HbmBundleLoader.loadSemanticProvenance(hbm.provenance);
				final double[] background=NewFamilyHbmBuilder.readBackground(hbm.background);
				final HbmBundleLoader.Loaded loaded=HbmBundleLoader.load(Paths.get(hbm.hbm), Arrays.asList(ids), byId::get, provenance);
				require(loaded.familyCount()==ids.length, "HBM family count differs from consensus roster");
				for(int i=0; i<ids.length; i++){require(ids[i].equals(loaded.repId(i)), "HBM roster differs at rank "+i);}
				final HbmPositionModel[] out=loaded.positionModels("logodds", .01, true, -4, background);
				for(int i=0; i<out.length; i++){require(out[i].length()==refs[i].length, "HBM length differs at rank "+i+": "+ids[i]);}
				return out;
			}catch(Exception e){throw new RuntimeException(e);}
		}

		int shortlist(Query query,int[] ranks,int minHits,Scratch scratch){
			if(paired!=null){
				final int count=paired.rank(query.raw, query.protein, scratch.paired);
				System.arraycopy(scratch.paired.union, 0, ranks, 0, count);
				return count;
			}
			if(exhaustive){
				for(int i=0; i<ranks.length; i++){ranks[i]=i;}
				return ranks.length;
			}
			final LongHashSet set=scratch.set;
			final LongList words=scratch.words;
			final IntList touched=scratch.touched;
			final int[] hits=scratch.hits;
			ReducedAlphabetSeedAssay.fillKmerSet(query.protein.enc, alphabet, k, ReducedAlphabetSeedAssay.Boundary.NONE, set, words);
			for(int i=0; i<touched.size; i++){hits[touched.get(i)]=0;}
			touched.clear();
			for(int i=0; i<words.size; i++){
				final LongList list=index.get(words.get(i));
				if(list==null){continue;}
				for(int j=0; j<list.size; j++){
					final int rank=(int)list.get(j);
					if(hits[rank]++==0){touched.add(rank);}
				}
			}
			int count=0;
			for(int i=0; i<touched.size; i++){
				final int rank=touched.get(i);
				if(hits[rank]<minHits || (count>=ranks.length && !better(rank, ranks[count-1], hits))){continue;}
				int j=(count<ranks.length ? count++ : count-1);
				while(j>0 && better(rank, ranks[j-1], hits)){
					ranks[j]=ranks[j-1]; j--;
				}
				ranks[j]=rank;
			}
			return count;
		}
		boolean noKmers(Scratch scratch){
			return paired==null ? !exhaustive && scratch.words.size==0 :
				scratch.paired.queryRegisters==0 && scratch.paired.coverQueryRegisters==0;
		}
		boolean better(int a,int b,int[] hits){
			return hits[a]>hits[b] || (hits[a]==hits[b] && ids[a].compareTo(ids[b])<0);
		}
		final String[] ids;
		final byte[][] refs;
		final HbmPositionModel[] profiles;
		final HbmPositionModel[] hbmProfiles;
		final LongLongListHashMap index;
		final ReducedAlphabetSeedAssay.Alphabet alphabet;
		final HbmTailPairedShortlist paired;
		final boolean exhaustive;
		final int k, rankLimit;
	}

	public static final class Scratch {
		private Scratch(Reference ref){
			owner=ref;
			ranks=new int[ref.rankLimit];
			if(ref.paired==null){
				set=new LongHashSet(); words=new LongList(); touched=new IntList();
				hits=new int[ref.ids.length]; paired=null;
			}else{
				set=null; words=null; touched=null; hits=null; paired=ref.paired.makeScratch();
			}
		}
		void clear(){aai.clear(); blosum.clear(); hbm.clear();}
		final Reference owner;
		final int[] ranks;
		final LongHashSet set;
		final LongList words;
		final IntList touched;
		final int[] hits;
		final HbmTailPairedShortlist.Scratch paired;
		final MutableHit aai=new MutableHit(Method.AAI);
		final MutableHit blosum=new MutableHit(Method.BLOSUM);
		final MutableHit hbm=new MutableHit(Method.HBM);
	}

	public static final class Result {
		private Result(int candidates_,boolean noKmers_,MutableHit aai_,MutableHit blosum_,MutableHit hbm_,
				Work blosumWork_,Work hbmWork_){
			candidates=candidates_; noKmers=noKmers_;
			hits[METHOD_AAI]=new Hit(aai_);
			hits[METHOD_BLOSUM]=new Hit(blosum_);
			if(hbm_!=null){hits[METHOD_HBM]=new Hit(hbm_);}
			blosumWork=blosumWork_; hbmWork=hbmWork_;
		}
		public Hit hit(Method method){
			final Hit hit=hits[method.ordinal()];
			if(hit==null){throw new IllegalStateException("Stop-tail method was not loaded: "+method.name);}
			return hit;
		}
		public double score(Method method){return hit(method).score;}
		public final int candidates;
		public final boolean noKmers;
		public final Work blosumWork, hbmWork;
		private final Hit[] hits=new Hit[3];
	}

	public static final class Hit {
		private Hit(MutableHit source){
			method=source.method; geometryPassed=source.hasAlignment();
			score=source.score; score64=source.score64; columns=source.columns;
			rank=source.rank; representative=source.representative;
			paired=source.paired; terminalPaired=source.terminalPaired;
			queryStart=source.queryStart; queryEnd=source.queryEnd;
			referenceStart=source.referenceStart; referenceEnd=source.referenceEnd;
		}
		public boolean passes(double cutoff){
			require(Double.isFinite(cutoff), "Stop-tail score cutoff must be finite: "+cutoff);
			return geometryPassed && score>=cutoff;
		}
		public final Method method;
		/** True only when this method found a finite alignment passing geometry gates. */
		public final boolean geometryPassed;
		public final double score;
		/** Raw score64 for the selected alignment; AAI normalizes matches from the same trace. */
		public final int score64, columns, rank;
		public final String representative;
		public final int paired, terminalPaired, queryStart, queryEnd, referenceStart, referenceEnd;
	}

	public static final class Work {
		private Work(int scored_,int traces_,int geometryFail_,int accepted_){
			scored=scored_; traces=traces_; geometryFail=geometryFail_; accepted=accepted_;
		}
		public final int scored, traces, geometryFail, accepted;
	}

	private static final class MutableHit {
		MutableHit(Method method_){method=method_; clear();}
		void clear(){
			score=Double.NEGATIVE_INFINITY; score64=Integer.MIN_VALUE; columns=rank=-1; representative=null;
			paired=terminalPaired=0; queryStart=queryEnd=referenceStart=referenceEnd=-1;
		}
		boolean hasAlignment(){return score!=Double.NEGATIVE_INFINITY;}
		void update(Method method_,double score_,int score64_,int columns_,int rank_,String representative_,
				HbmPositionModel.TailResult tr,HbmTailEvidence.Metrics m){
			assert(method_==method && Double.isFinite(score_) && columns_>0 && tr!=null && m!=null)
				: "Geometry-passed detector hit must keep finite score and measured geometry";
			if(score_<=score){return;}
			score=score_; score64=score64_; columns=columns_; rank=rank_; representative=representative_;
			paired=m.paired; terminalPaired=m.terminalPaired;
			queryStart=tr.queryStart; queryEnd=tr.queryEnd;
			referenceStart=tr.start; referenceEnd=tr.end;
		}
		final Method method;
		double score;
		int score64, columns, rank, paired, terminalPaired, queryStart, queryEnd, referenceStart, referenceEnd;
		String representative;
	}

	private static boolean safe(String s){
		if(s==null || s.isEmpty()){return false;}
		for(int i=0; i<s.length(); i++){final char c=s.charAt(i); if(c<=32 || c>=127 || c=='/' || c=='\\'){return false;}}
		return true;
	}

	private static void verify(String ref,String refPin,HbmInput hbm){
		HbmProfilePilot.requireHash(ref, refPin);
		if(hbm!=null){hbm.verify();}
	}

	private static void require(boolean ok,String message){if(!ok){throw new IllegalArgumentException(message);}}

	public static final int METHOD_AAI=Method.AAI.ordinal(), METHOD_BLOSUM=Method.BLOSUM.ordinal(), METHOD_HBM=Method.HBM.ordinal();
	private StopTailDetector(){}
}
