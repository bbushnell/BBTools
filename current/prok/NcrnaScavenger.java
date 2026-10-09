package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.IdentityHashMap;

import consensus.BaseGraph;
import dna.AminoAcid;
import idaligner.AlignmentStats;
import idaligner.EndClippedAligner;
import idaligner.QuantumAligner;
import idaligner.ScrabbleAligner;
import map.LongHashSet;
import ml.CellNet;
import shared.Tools;
import structures.IntList;

/**
 * Generic conserved-ncRNA scavenger: the family-agnostic core of TrnaCaller's
 * scavenger pass (conserved-kmer seed -> candidate window -> kmer-index
 * shortlist -> ScrabbleAligner/QuantumAligner verification -> HBM fallback),
 * factored out per Noire's read of TrnaCaller (2026-08-23) so any conserved
 * ncRNA family with a consensus library, HBM model set, and conserved
 * long-kmer set can be called the same way -- without tRNA's PGM-based
 * candidate generation (callTrnas/scanInner/findRegions/extractTrnas),
 * anticodon logic, acceptor-stem trimming, or intron handling, all of which
 * stay tRNA-only in TrnaCaller.
 *
 * <p>Boundary trimming here is deliberately simplified from TrnaCaller's:
 * an extended window is aligned to the winning model with traceback and the
 * ORF is snapped to the alignment's own extent (rStart/rStop) -- no
 * acceptor-stem search, no anticodon extraction. That structural refinement
 * is tRNA-specific and stays in TrnaCaller.
 *
 * <p>Per-family tunables (window pad, min length, kmer set/k) are
 * constructor parameters since they genuinely differ per family (tRNA
 * windowPad=83; RNase P/SRP need their own). The remaining knobs (index k,
 * identity thresholds, etc.) are mutable instance fields defaulted to
 * TrnaCaller's measured-best tRNA values, per Noire's "start with tRNA's
 * defaults, tune per family later" -- not yet wired to per-family values.
 * @author Neptune, Noire, Brian Bushnell
 */
public class NcrnaScavenger {

	public NcrnaScavenger(byte[][] library_, BaseGraph[] models_, String[] modelNames_,
			LongHashSet kmerSet_, int kLong_, int minLen_, int windowPad_){
		this(library_, models_, modelNames_, kmerSet_, kLong_, minLen_, windowPad_,
				7, 60, true, 11f, 0.48f, 0.072f, 12);
	}

	public NcrnaScavenger(byte[][] library_, BaseGraph[] models_, String[] modelNames_,
			LongHashSet kmerSet_, int kLong_, int minLen_, int windowPad_,
			int indexK_, int indexTopN_, boolean adaptive_,
			float adaptFloor_, float adaptTopFrac_, float adaptQFrac_, int fixedMinHits_){
		this(library_, models_, modelNames_, kmerSet_, kLong_, minLen_, windowPad_,
				indexK_, indexTopN_, adaptive_, adaptFloor_, adaptTopFrac_, adaptQFrac_, fixedMinHits_,
				0f, 20f, 0.75f, 0.65f);
	}

	/** Forward-ported from Noire's ncRNA-family-loading tree (2026-08-28, C3 merge) -- see
	 * NcrnaFamily's matching constructor javadoc for the scoreA/scoreB/idPass/idBorderline
	 * rationale (Noire's #1 recall lever, +11.5pp rnasep). idPass/idBorderline are now
	 * FINAL, set here -- previously mutable instance fields defaulted to 0.75f/0.65f
	 * (unchanged as the fallback defaults for this overload's own delegation below). */
	public NcrnaScavenger(byte[][] library_, BaseGraph[] models_, String[] modelNames_,
			LongHashSet kmerSet_, int kLong_, int minLen_, int windowPad_,
			int indexK_, int indexTopN_, boolean adaptive_,
			float adaptFloor_, float adaptTopFrac_, float adaptQFrac_, int fixedMinHits_,
			float scoreA_, float scoreB_, float idPass_, float idBorderline_){
		this(library_, models_, modelNames_, kmerSet_, kLong_, minLen_, windowPad_,
				indexK_, indexTopN_, adaptive_, adaptFloor_, adaptTopFrac_, adaptQFrac_, fixedMinHits_,
				scoreA_, scoreB_, idPass_, idBorderline_,
				null, null, null, null, -1, -1, -1, -1, 0f);
	}

	/** Full constructor, adding C3's boundary-precision-NN resources (Noire's spec,
	 * plans/c3_ncrnaboundaryscorer_spec.md; G11, 2026-08-28). boundary5NetTemplate==null
	 * means OFF -- structural, matches NcrnaFamily's own default. Per-instance CLONES its
	 * own net copies from the shared read-only templates (mirrors TrnaCaller's
	 * boundary5Net/boundary3Net constructor-clone pattern exactly, same reason: CellNet.
	 * feedForward mutates per-Cell state, so concurrent use of ONE net object across
	 * per-thread NcrnaScavenger instances would corrupt it). meanLen is asserted >0
	 * whenever a template is given -- there is no safe shared default across families of
	 * very different length (rnasep ~380bp vs srp_small ~95bp), so a caller passing a real
	 * net without a real meanLen is a construction-time bug, not a runtime one -- fail
	 * loud immediately rather than silently score every candidate with lengthRatio=0. */
	public NcrnaScavenger(byte[][] library_, BaseGraph[] models_, String[] modelNames_,
			LongHashSet kmerSet_, int kLong_, int minLen_, int windowPad_,
			int indexK_, int indexTopN_, boolean adaptive_,
			float adaptFloor_, float adaptTopFrac_, float adaptQFrac_, int fixedMinHits_,
			float scoreA_, float scoreB_, float idPass_, float idBorderline_,
			CellNet boundary5NetTemplate, CellNet boundary3NetTemplate,
			TrnaBoundaryFeatures.NinemerTable boundaryStartTable_, TrnaBoundaryFeatures.NinemerTable boundaryStopTable_,
			int boundaryStartInside_, int boundaryStartOutside_, int boundaryStopInside_, int boundaryStopOutside_,
			float boundaryMeanLen_){
		this(library_, models_, modelNames_, kmerSet_, kLong_, minLen_, windowPad_,
				indexK_, indexTopN_, adaptive_, adaptFloor_, adaptTopFrac_, adaptQFrac_, fixedMinHits_,
				scoreA_, scoreB_, idPass_, idBorderline_, boundary5NetTemplate, boundary3NetTemplate,
				boundaryStartTable_, boundaryStopTable_, boundaryStartInside_, boundaryStartOutside_,
				boundaryStopInside_, boundaryStopOutside_, boundaryMeanLen_,
				NcrnaFamily.LEGACY_START_OFFSETS, NcrnaFamily.LEGACY_STOP_OFFSETS);
	}

	public NcrnaScavenger(byte[][] library_, BaseGraph[] models_, String[] modelNames_,
			LongHashSet kmerSet_, int kLong_, int minLen_, int windowPad_,
			int indexK_, int indexTopN_, boolean adaptive_,
			float adaptFloor_, float adaptTopFrac_, float adaptQFrac_, int fixedMinHits_,
			float scoreA_, float scoreB_, float idPass_, float idBorderline_,
			CellNet boundary5NetTemplate, CellNet boundary3NetTemplate,
			TrnaBoundaryFeatures.NinemerTable boundaryStartTable_, TrnaBoundaryFeatures.NinemerTable boundaryStopTable_,
			int boundaryStartInside_, int boundaryStartOutside_, int boundaryStopInside_, int boundaryStopOutside_,
			float boundaryMeanLen_, int[] boundaryStartOffsets_, int[] boundaryStopOffsets_){
		this(library_, models_, modelNames_, kmerSet_, kLong_, minLen_, windowPad_,
				indexK_, indexTopN_, adaptive_, adaptFloor_, adaptTopFrac_, adaptQFrac_, fixedMinHits_,
				scoreA_, scoreB_, idPass_, idBorderline_, boundary5NetTemplate, boundary3NetTemplate,
				boundaryStartTable_, boundaryStopTable_, boundaryStartInside_, boundaryStartOutside_,
				boundaryStopInside_, boundaryStopOutside_, boundaryMeanLen_,
				boundaryStartOffsets_, boundaryStopOffsets_, 0f, 0f);
	}

	public NcrnaScavenger(byte[][] library_, BaseGraph[] models_, String[] modelNames_,
			LongHashSet kmerSet_, int kLong_, int minLen_, int windowPad_,
			int indexK_, int indexTopN_, boolean adaptive_,
			float adaptFloor_, float adaptTopFrac_, float adaptQFrac_, int fixedMinHits_,
			float scoreA_, float scoreB_, float idPass_, float idBorderline_,
			CellNet boundary5NetTemplate, CellNet boundary3NetTemplate,
			TrnaBoundaryFeatures.NinemerTable boundaryStartTable_, TrnaBoundaryFeatures.NinemerTable boundaryStopTable_,
			int boundaryStartInside_, int boundaryStartOutside_, int boundaryStopInside_, int boundaryStopOutside_,
			float boundaryMeanLen_, int[] boundaryStartOffsets_, int[] boundaryStopOffsets_,
			float boundaryMarginStart_, float boundaryMarginStop_){
		TrnaKmerIndex.validateConfiguration(indexK_, adaptFloor_, adaptTopFrac_, adaptQFrac_, fixedMinHits_);
		TrnaKmerIndex.validateTopN(indexTopN_);
		if(!Float.isFinite(boundaryMarginStart_) || boundaryMarginStart_<0f ||
				!Float.isFinite(boundaryMarginStop_) || boundaryMarginStop_<0f){
			throw new IllegalArgumentException("ncRNA boundary margins must be finite and >=0: "
				+boundaryMarginStart_+", "+boundaryMarginStop_);
		}
		library=library_;
		models=models_;
		modelNames=modelNames_;
		annotate=(modelNames!=null);
		kmerSet=kmerSet_;
		kLong=kLong_;
		minLen=minLen_;
		windowPad=windowPad_;
		indexK=indexK_;
		indexTopN=indexTopN_;
		adaptiveMinHits=adaptive_;
		adaptFloor=adaptFloor_;
		adaptTopFrac=adaptTopFrac_;
		adaptQFrac=adaptQFrac_;
		indexMinHitsDefault=fixedMinHits_;
		scoreA=scoreA_;
		scoreB=scoreB_;
		idPass=idPass_;
		idBorderline=idBorderline_;
		kmerIndex=(library!=null ? new TrnaKmerIndex(library, indexK, adaptiveMinHits,
			adaptFloor, adaptTopFrac, adaptQFrac, indexMinHitsDefault) : null);
		boundary5Net=(boundary5NetTemplate!=null ? boundary5NetTemplate.copy(false) : null);
		boundary3Net=(boundary3NetTemplate!=null ? boundary3NetTemplate.copy(false) : null);
		boundaryStartTable=boundaryStartTable_;
		boundaryStopTable=boundaryStopTable_;
		boundaryStartInside=boundaryStartInside_;
		boundaryStartOutside=boundaryStartOutside_;
		boundaryStopInside=boundaryStopInside_;
		boundaryStopOutside=boundaryStopOutside_;
		boundaryMeanLen=boundaryMeanLen_;
		boundaryStartOffsets=boundaryStartOffsets_.clone();
		boundaryStopOffsets=boundaryStopOffsets_.clone();
		boundaryMarginStart=boundaryMarginStart_;
		boundaryMarginStop=boundaryMarginStop_;
		if(boundary5Net!=null && !(boundaryMeanLen>0)){
			throw new IllegalArgumentException("NcrnaScavenger built with a boundary-precision net but boundaryMeanLen="
				+boundaryMeanLen_+" (must be >0).");
		}
	}

	public long alignmentCount(){return alignmentCount;}
	/** Sum of candidate input spans at actual alignment calls, including endpoint trimming. */
	public long alignedBases(){return alignedBases;}
	/** HBM evaluations count separately; reusing a traceback is not another alignment. */
	public long hbmScoreCalls(){return hbmScoreCalls;}
	public long hbmBasesScored(){return hbmBasesScored;}
	/** Seed occurrences admitted to this family's scavenger after its input/minimum-length guards. */
	public long kmerHitCount(){return kmerHitCount;}
	/** Windows scheduled across both passes, before the seed-count and model filters. */
	public long windowCount(){return windowCount;}
	public long boundaryNanos(){return boundaryNanos;}
	boolean hasPerModelBoundaryResources(){return boundaryNetsByModel!=null;}
	boolean endpointUsesPerSiteScrabble(){return NcrnaBoundaryScorer.usesPerSiteScrabble(boundaryFeatureVersion);}

	/** Development-harness support for model-dispatched V3 endpoint resources.
	 * The ordinary family-wide resources remain the fallback and all production
	 * construction paths remain unchanged. Each net is cloned because feedForward
	 * mutates Cell state. Consensus tables are read-only after loading; the existing
	 * family-wide tables supply the second half of the V3 profile blend. */
	void setPerModelBoundaryResources(CellNet[] nets,
			TrnaBoundaryFeatures.NinemerTable[] startTables,
			TrnaBoundaryFeatures.NinemerTable[] stopTables){
		if(nets==null || startTables==null || stopTables==null || nets.length!=library.length
				|| startTables.length!=library.length || stopTables.length!=library.length){
			throw new IllegalArgumentException("Per-model boundary resources must exactly match consensus count "+library.length);
		}
		boundaryNetsByModel=new CellNet[library.length];
		boundaryStartTablesByModel=startTables.clone();
		boundaryStopTablesByModel=stopTables.clone();
		for(int i=0; i<library.length; i++){
			if(nets[i]==null || startTables[i]==null || stopTables[i]==null){
				throw new IllegalArgumentException("Missing per-model boundary resources for consensus index "+i);
			}
			boundaryNetsByModel[i]=nets[i].copy(false);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Scavenger Pass       ----------------*/
	/*--------------------------------------------------------------*/

	public ArrayList<Orf> scavenge(String name, byte[] bases, int strand, ArrayList<int[]> called){
		if(bases==null || bases.length<minLen || library==null || kmerSet==null){return new ArrayList<Orf>();}
		final int[] hits=findKmerHitPositions(bases);
		return scavenge(name, bases, strand, called, hits, votingEnabled() || joinedPositions!=null ? findKmerHitKeys(bases, hits) : null);
	}

	/** Runs the unchanged family-local scavenger using a hit stream produced by
	 * the shared 17-mer front end. */
	ArrayList<Orf> scavenge(String name, byte[] bases, int strand, ArrayList<int[]> called,
			int[] hitPositions){
		return scavenge(name, bases, strand, called, hitPositions, null);
	}

	ArrayList<Orf> scavenge(String name, byte[] bases, int strand, ArrayList<int[]> called,
			int[] hitPositions, long[] hitKeys){
		ArrayList<Orf> results=new ArrayList<>();
		if(bases==null || bases.length<minLen || library==null || kmerSet==null){return results;}
		if(hitPositions==null){throw new IllegalArgumentException("Null ncRNA seed-hit stream for "+family);}
		if((votingEnabled() || joinedPositions!=null) && (hitKeys==null || hitKeys.length!=hitPositions.length)){
			throw new IllegalArgumentException("Seed voting requires one key per hit for "+family);
		}
		kmerHitCount+=hitPositions.length;
		if(workloadSink!=null){workloadSink.seedHits(name, strand, Arrays.copyOf(hitPositions, hitPositions.length));}
		if(diagSink!=null){diagSink.seed(diagFamily,name,strand,bases.length,Arrays.copyOf(hitPositions,hitPositions.length));}
		if(DEBUG){System.err.println("DEBUG scavenge name="+name+" strand="+strand+" bases.length="+bases.length
			+" hitPositions="+Arrays.toString(hitPositions));}
		if(hitPositions.length==0){return results;}
		final ArrayList<int[]> joinedPrior=hasJoinedRoute()?new ArrayList<int[]>(called):null;
		ArrayList<int[]> windows=prepareCandidateWindows(name, bases, strand, 1, hitPositions, hitKeys);
		if(DEBUG){System.err.println("DEBUG windows after collapse: "+dumpWindows(windows));}
		windows=subtractClaimed(windows, called);
		if(DEBUG){System.err.println("DEBUG windows after subtractClaimed: "+dumpWindows(windows));}
		if(refreshClaimedWindows){
			alignPassAgainstSnapshot(name,bases,strand,1,windows,new ArrayList<int[]>(called),hitPositions,hitKeys,called,results);
		}else{
			for(int[] original : windows){
				windowCount++;
				if(workloadSink!=null){workloadSink.scheduledWindow(name, strand, 1, original[0], original[1]);}
				Orf orf=alignWindow(name, bases, strand, 1, original[0], original[1], windowSource(original), hitPositions, hitKeys);
				if(DEBUG){System.err.println("DEBUG alignWindow("+original[0]+","+original[1]+") -> "+(orf==null ? "null" : (orf.start+"-"+orf.stop+" score="+orf.orfScore)));}
				if(orf!=null){called.add(new int[]{orf.start, orf.stop});results.add(orf);emitInstrumentation(orf);}
			}
			if(hasJoinedRoute()){
				final ArrayList<Orf> joined=new ArrayList<Orf>();appendJoinedCandidates(name,bases,strand,joinedPrior,hitPositions,hitKeys,joined);
				for(Orf orf:joined){called.add(new int[]{orf.start,orf.stop});results.add(orf);emitInstrumentation(orf);}
			}
		}
		if(scavengePass2){
			int[] nearHits=findNearbyUnclaimed(hitPositions, called, bases.length);
			long[] nearKeys=(votingEnabled() ? subsetKeys(hitPositions, hitKeys, nearHits) : null);
			if(nearHits.length>0){
				ArrayList<int[]> pass2Windows=prepareCandidateWindows(name, bases, strand, 2, nearHits, nearKeys);
				pass2Windows=subtractClaimed(pass2Windows, called);
				if(refreshClaimedWindows){
					alignPassAgainstSnapshot(name,bases,strand,2,pass2Windows,new ArrayList<int[]>(called),nearHits,nearKeys,called,results);
				}else{
					for(int[] original : pass2Windows){
						windowCount++;
						if(workloadSink!=null){workloadSink.scheduledWindow(name, strand, 2, original[0], original[1]);}
						Orf orf=alignWindow(name, bases, strand, 2, original[0], original[1], windowSource(original), nearHits, nearKeys);
						if(orf!=null){called.add(new int[]{orf.start, orf.stop});results.add(orf);emitInstrumentation(orf);}
					}
				}
			}
		}
		return results;
	}

	/** Evaluates one pass against a stable snapshot of claims that existed before the
	 * pass.  New candidates cannot starve later, higher-scoring candidates merely by
	 * appearing first; mutually overlapping same-family candidates are resolved after
	 * all windows have been evaluated, then the winners become claims for the next pass. */
	private void alignPassAgainstSnapshot(String name,byte[] bases,int strand,int pass,ArrayList<int[]> windows,
			ArrayList<int[]> claimSnapshot,int[] hitPositions,long[] hitKeys,ArrayList<int[]> called,ArrayList<Orf> results){
		final ArrayList<Orf> candidates=new ArrayList<Orf>();
		for(int[] original:windows){
			final ArrayList<int[]> pending=new ArrayList<int[]>(1);pending.add(original);
			for(int[] w:subtractClaimed(pending,claimSnapshot)){
				windowCount++;
				if(workloadSink!=null){workloadSink.scheduledWindow(name,strand,pass,w[0],w[1]);}
				final Orf orf=alignWindow(name,bases,strand,pass,w[0],w[1],windowSource(w),hitPositions,hitKeys);
				if(DEBUG){System.err.println("DEBUG alignWindow("+w[0]+","+w[1]+") -> "+(orf==null ? "null" : (orf.start+"-"+orf.stop+" score="+orf.orfScore)));}
				if(orf!=null){candidates.add(orf);}
			}
		}
		if(pass==1 && hasJoinedRoute()){appendJoinedCandidates(name,bases,strand,claimSnapshot,hitPositions,hitKeys,candidates);}
		commitSnapshotCandidates(candidates,called,results);
	}
	private void appendJoinedCandidates(String name,byte[] bases,int strand,ArrayList<int[]> claimSnapshot,int[] hitPositions,long[] hitKeys,ArrayList<Orf> candidates){
		assert(hasJoinedRoute() && claimSnapshot!=null):"Joined windows must see claims from before ordinary candidates are committed";
			if(cmVerifier!=null || !pacBioConsensusAlignment || !pacBioRolling){throw new IllegalStateException("Joined fast replay must use rolling MSA and no CM");}
			final Euk18sJoinedProposals.Group group=joinedPositions==null?joinedProposals.get(name,strand):liveJoinedGroup(bases.length,strand,hitPositions,hitKeys);
			if(group!=null){for(Euk18sJoinedProposals.Proposal proposal:group.rows){
				if(proposal.stop>bases.length){throw new IllegalArgumentException("Joined proposal exceeds its original shred: "+name);}
				final int lo=strand==0?proposal.start-1:bases.length-proposal.stop,hi=strand==0?proposal.stop-1:bases.length-proposal.start;
				final ArrayList<int[]> pending=new ArrayList<int[]>(1);pending.add(new int[]{lo,hi,WINDOW_JOINED});
				for(int[] w:subtractClaimed(pending,claimSnapshot)){
					windowCount++;if(workloadSink!=null){workloadSink.scheduledWindow(name,strand,1,w[0],w[1]);}
					final Orf orf=alignJoinedWindow(name,bases,strand,w[0],w[1],proposal.model,hitPositions,hitKeys);
					if(orf!=null){candidates.add(orf);}
				}
			}}
	}
	private boolean hasJoinedRoute(){return joinedProposals!=null || joinedPositions!=null;}
	private Euk18sJoinedProposals.Group liveJoinedGroup(int length,int strand,int[] centers,long[] keys){
		assert(joinedPositions!=null && centers!=null && keys!=null):"Live proposals consume the existing family hit stream and its bound positions";
		if(joinedBuilder==null){joinedBuilder=new SeedInsertionWindowBuilder();}
		final SeedInsertionWindowBuilder.Result built=joinedBuilder.build(centers,keys,length,joinedPositions,joinedSideSupport,60,60,256,4096);
		final Euk18sJoinedProposals.Group group=new Euk18sJoinedProposals.Group();
		final map.LongHashSet[] seen=new map.LongHashSet[joinedPositions.modelCount()];
		for(SeedInsertionWindowBuilder.Proposal p:built.proposals){
			if(!p.eligible()){continue;}if(seen[p.model]==null){seen[p.model]=new map.LongHashSet();}
			if(!seen[p.model].add(((long)p.start<<32)|(p.stop&0xffffffffL))){continue;}
			final int start=strand==0?p.start+1:length-p.stop,stop=strand==0?p.stop+1:length-p.start;
			group.rows.add(new Euk18sJoinedProposals.Proposal(start,stop,p.model));
		}
		return group;
	}
	/** Same acceptance and score as the fast caller, with tagged span/core handling. */
	Orf alignJoinedWindow(String name,byte[] bases,int strand,int start,int stop,int model,int[] hitPositions,long[] hitKeys){
		if(cmVerifier!=null || !pacBioConsensusAlignment || !pacBioRolling || start<0 || stop>=bases.length || start>stop){throw new IllegalArgumentException("Joined fast window must be physical, rolling and CM-free");}
		final byte[] seq=copyRegionUpper(bases,start,stop+1);
		if(modelAttemptSink!=null){kmerIndex.shortlist(seq,indexTopN,indexScoreMargin);}// Diagnostic counts must not reuse the preceding ordinary window's counts.
		return alignWindowOnce(name,bases,strand,1,start,stop,WINDOW_JOINED,hitPositions,hitKeys,seq,new int[]{model},countHits(hitPositions,start+8,stop-8),false,true);
	}

	/** Commits only snapshot-resolver winners and emits capture for those returned calls. */
	void commitSnapshotCandidates(ArrayList<Orf> candidates,ArrayList<int[]> called,ArrayList<Orf> results){
		final ArrayList<Orf> selected=selectBestNonOverlapping(candidates);
		for(Orf orf:selected){
			called.add(new int[]{orf.start,orf.stop});results.add(orf);emitInstrumentation(orf);
		}
		emitResolvedInstrumentation(candidates,selected);
	}

	/** Keeps the maximum-score call in each reciprocal-half-overlap conflict component.
	 * Low-overlap adjacent same-family genes remain distinct; this is required for the
	 * Pneumocystis LSU repeat pair that shares only 277 bp of roughly 3.4 kb. */
	static ArrayList<Orf> selectBestNonOverlapping(ArrayList<Orf> candidates){
		final ArrayList<Orf> sorted=new ArrayList<Orf>(candidates);
		sorted.sort((a,b)->compareSnapshotCandidates(a,b));
		final int n=sorted.size();if(n<2){return sorted;}
		final int[] parent=new int[n],winner=new int[n];for(int i=0;i<n;i++){parent[i]=i;winner[i]=-1;}
		for(int i=0;i<n;i++){for(int j=i+1;j<n;j++){if(snapshotConflict(sorted.get(i),sorted.get(j))){unionSnapshot(parent,i,j);}}}
		for(int i=0;i<n;i++){final int root=findSnapshot(parent,i),old=winner[root];if(old<0||sorted.get(i).orfScore>sorted.get(old).orfScore||(sorted.get(i).orfScore==sorted.get(old).orfScore&&compareSnapshotCandidates(sorted.get(i),sorted.get(old))<0)){winner[root]=i;}}
		final ArrayList<Orf> selected=new ArrayList<Orf>();
		for(int i=0;i<n;i++){if(winner[findSnapshot(parent,i)]==i){selected.add(sorted.get(i));}}
		Collections.sort(selected);
		for(int i=0;i<selected.size();i++){for(int j=i+1;j<selected.size();j++){assert(!snapshotConflict(selected.get(i),selected.get(j))) : "Snapshot resolver retained a reciprocal-half-overlap same-family conflict";}}
		return selected;
	}

	static boolean snapshotConflict(Orf a,Orf b){
		if(a.strand!=b.strand||!a.scafName.equals(b.scafName)||compareNullable(a.ncrnaFamily,b.ncrnaFamily)!=0){return false;}
		final int overlap=Tools.min(a.stop,b.stop)-Tools.max(a.start,b.start)+1;
		return overlap>0&&2L*overlap>=a.stop-a.start+1L&&2L*overlap>=b.stop-b.start+1L;
	}

	private static int findSnapshot(int[] parent,int x){int root=x;while(parent[root]!=root){root=parent[root];}while(parent[x]!=x){final int next=parent[x];parent[x]=root;x=next;}return root;}
	private static void unionSnapshot(int[] parent,int a,int b){final int x=findSnapshot(parent,a),y=findSnapshot(parent,b);if(x!=y){parent[Tools.max(x,y)]=Tools.min(x,y);}}

	/** Total biological/provenance order used before DP.  Exact weight ties prefer
	 * the earlier item in this order because the recurrence keeps the existing path. */
	static int compareSnapshotCandidates(Orf a,Orf b){
		int x=Integer.compare(a.stop,b.stop);if(x!=0){return x;}
		x=Integer.compare(a.start,b.start);if(x!=0){return x;}
		x=-Float.compare(a.orfScore,b.orfScore);if(x!=0){return x;}
		x=compareNullable(a.trnaModel,b.trnaModel);if(x!=0){return x;}
		x=compareNullable(a.ncrnaFamily,b.ncrnaFamily);if(x!=0){return x;}
		x=Integer.compare(a.type,b.type);if(x!=0){return x;}
		x=Integer.compare(a.strand,b.strand);if(x!=0){return x;}
		return Integer.compare(a.frame,b.frame);
	}

	private static int compareNullable(String a,String b){
		if(a==b){return 0;}if(a==null){return -1;}if(b==null){return 1;}return a.compareTo(b);
	}

	/** Temporary diagnostic flag (Neptune, 2026-08-23): traces kmer hits, candidate window
	 * construction, and per-model alignment scores for FN root-cause diagnosis. Default false,
	 * zero cost when off. Not wired to a command-line flag -- toggle+recompile for one-off use. */
	static boolean DEBUG=Boolean.getBoolean("ncrna.debug");

	private static String dumpWindows(ArrayList<int[]> windows){
		StringBuilder sb=new StringBuilder("[");
		for(int[] w : windows){sb.append("(").append(w[0]).append(",").append(w[1]).append(") ");}
		sb.append("]");
		return sb.toString();
	}

	int[] findKmerHitPositions(byte[] bases){
		if(kmerSet==null || kLong<=0 || kLong>31 || bases.length<kLong){return EMPTY;}
		final long kmask=~((-1L)<<(2*kLong));
		final byte[] bton=AminoAcid.baseToNumber;
		IntList hits=new IntList();
		long kmer=0; int len=0;
		for(int i=0; i<bases.length; i++){
			final int x=bton[bases[i]];
			if(x>=0){
				kmer=((kmer<<2)|x)&kmask; len++;
				if(len>=kLong && kmerSet.contains(kmer)){hits.add(i-kLong/2);}
			}else{len=0; kmer=0;}
		}
		return hits.toArray();
	}

	private long[] findKmerHitKeys(byte[] bases, int[] hits){
		final long[] keys=new long[hits.length];
		final byte[] bton=AminoAcid.baseToNumber;
		for(int i=0; i<hits.length; i++){
			final int start=hits[i]-kLong/2; long key=0;
			for(int j=0; j<kLong; j++){key=(key<<2)|bton[bases[start+j]];}
			keys[i]=key;
		}
		return keys;
	}

	/** Restricts vote observation to generated proposals, before claims or endpoint voting.
	 * Diagnostic-only keys never replace the science stream used by computeVotes. */
	private ArrayList<int[]> prepareCandidateWindows(String name, byte[] bases, int strand, int pass, int[] positions, long[] keys){
		if(voteDiagSink==null || !voteWindows){return prepareCandidateWindows(positions, keys, bases.length);}
		assert(voteTrace==null) : shared.KillSwitch.assertDie("Window trace context must not leak between scavenger passes or endpoint voting");
		try{
			final long[] observedKeys=keys!=null ? keys : diagnosticKeys(bases, positions);
			voteTrace=new NcrnaVoteDiagnostics(voteDiagSink, voteDiagFamily, name, strand, bases.length, pass,
				voteSlack, windowPad, positions, observedKeys, voteTable);
			return prepareCandidateWindows(positions, keys, bases.length);
		}catch(RuntimeException | AssertionError e){
			shared.KillSwitch.assertDie("Generated-vote diagnostics failed on caller worker: "+e);throw e;
		}finally{voteTrace=null;}
	}

	/** Unknown or invalid diagnostic-only reconstruction remains unavailable. */
	private long[] diagnosticKeys(byte[] bases, int[] positions){
		final long[] keys=new long[positions.length];Arrays.fill(keys, NcrnaVoteDiagnostics.UNKNOWN);
		if(kLong<1 || kLong>31){return keys;}
		for(int i=0; i<positions.length; i++){
			final long start=(long)positions[i]-kLong/2;
			if(start<0 || start+kLong>bases.length){continue;}
			long key=0;boolean valid=true;
			for(int j=0; j<kLong; j++){
				final int base=bases[(int)start+j];
				final int x=base<0 || base>=AminoAcid.baseToNumber.length ? -1 : AminoAcid.baseToNumber[base];
				if(x<0){valid=false;break;}key=(key<<2)|x;
			}
			if(valid){keys[i]=key;}
		}
		return keys;
	}

	/** One window preparation seam for both passes; opt-in changes no legacy family. */
	ArrayList<int[]> prepareCandidateWindows(int[] positions, long[] keys, int seqLen){
		if(voteBeforePadding){
			if(!voteWindows || voteTable==null){throw new IllegalStateException("Vote-before-padding requires enabled window votes and a bound table: "+family);}
			if(voteWindowBuilder==null){
				voteWindowBuilder=new SeedVoteWindowBuilder();final int[] lengths=new int[library.length];
				for(int i=0; i<lengths.length; i++){lengths[i]=library[i].length;}Arrays.sort(lengths);fallbackCoreLength=lengths[lengths.length/2];
			}
			final ArrayList<int[]> windows;
			try{windows=voteWindowBuilder.build(positions, keys, seqLen, voteTable, kLong, minLen, windowPad, voteSlack, collapseFrac, fallbackCoreLength, voteTrace);}
			catch(RuntimeException | AssertionError e){shared.KillSwitch.assertDie("Seed-voted window construction failed on caller worker: "+e);throw e;}
			votedWindows+=voteWindowBuilder.voted;voteWindowFallbacks+=voteWindowBuilder.fallbackWindows;return windows;
		}
		final ArrayList<int[]> windows=collapseByIntersection(buildCandidateWindows(positions, seqLen));
		if(voteWindows){applyVoteWindows(windows, positions, keys, seqLen);}return windows;
	}

	private ArrayList<int[]> buildCandidateWindows(int[] hitPositions, int seqLen){
		ArrayList<int[]> windows=new ArrayList<>();
		for(int center : hitPositions){
			int start=(int)Math.max(0L, (long)center-windowPad);
			int stop=(int)Math.min(seqLen-1L, (long)center+windowPad+kLong);
			windows.add(new int[]{start, stop, WINDOW_COLLAPSE});
		}
		return windows;
	}

	/** No explicit upper bound on the merged window size -- matches TrnaCaller's own
	 * collapseByIntersection, which has never capped this either (SCAV_QUANTUM_THRESH
	 * only switches the aligner for large windows, it doesn't reject them). */
	private ArrayList<int[]> collapseByIntersection(ArrayList<int[]> windows){
		if(windows.size()<2){return windows;}
		windows.sort((a, b)->a[0]-b[0]);
		ArrayList<int[]> result=new ArrayList<>();
		int[] current=windows.get(0);
		for(int i=1; i<windows.size(); i++){
			int[] next=windows.get(i);
			int overlap=Tools.min(current[1], next[1])-Tools.max(current[0], next[0]);
			int shorter=Tools.min(current[1]-current[0], next[1]-next[0]);
			if(overlap>0 && overlap>=shorter*collapseFrac){
				current=new int[]{Tools.max(current[0], next[0]), Tools.min(current[1], next[1]), WINDOW_COLLAPSE};
			}else{
				if(current[1]-current[0]>=minLen){result.add(current);}
				else{diagnoseCollapsedDrop(current);}
				current=next;
			}
		}
		if(current[1]-current[0]>=minLen){result.add(current);}
		else{diagnoseCollapsedDrop(current);}
		return result;
	}

	private void diagnoseCollapsedDrop(int[] window){
		if(voteTrace==null){return;}
		voteTrace.resetVote();
		voteTrace.legacy(window[0], window[1], window[0], window[1], "dropped", "SHORT_COLLAPSED_WINDOW");
	}

	/** Subtracts every claimed interval from every window, keeping BOTH surviving remainders when a
	 * claim lands strictly inside a window (2026-08-27 fix, mirrors the identical fix in
	 * TrnaCaller.subtractClaimed -- see that javadoc for the root-cause trace; found via the tRNA path
	 * first, applied here because this class shares the same buggy pattern). Carries a list of
	 * surviving segments through every claim so a claim strictly inside a window correctly produces
	 * TWO output windows (left + right), not zero or one. */
	private ArrayList<int[]> subtractClaimed(ArrayList<int[]> windows, ArrayList<int[]> claimed){
		ArrayList<int[]> result=new ArrayList<>();
		for(int[] w : windows){
			ArrayList<int[]> segments=new ArrayList<>();
			segments.add(new int[]{w[0], w[1], windowSource(w)});
			for(int[] c : claimed){
				ArrayList<int[]> next=new ArrayList<>();
				for(int[] seg : segments){
					final int lo=seg[0], hi=seg[1];
					if(hi<c[0] || lo>c[1]){
						next.add(seg);//no overlap with this claim -- unchanged
						continue;
					}
					if(lo<c[0]){next.add(new int[]{lo, c[0]-1, windowSource(seg)});}//left remainder survives
					if(hi>c[1]){next.add(new int[]{c[1]+1, hi, windowSource(seg)});}//right remainder survives
					//neither branch fires -> claim fully covers this segment, it's consumed
				}
				segments=next;
			}
			//TODO: Probable bug (pre-existing, not fixed here -- see the matching TODO in
			//TrnaCaller.subtractClaimed) -- hi-lo>=minLen treats an inclusive [lo,hi] range as if it
			//had hi-lo bases, when it actually has hi-lo+1. Not changed here; flagging only.
			for(int[] seg : segments){
				if(seg[1]-seg[0]>=minLen){result.add(seg);}
			}
		}
		return result;
	}

	private int[] findNearbyUnclaimed(int[] allHits, ArrayList<int[]> claimed, int seqLen){
		IntList result=new IntList();
		for(int pos : allHits){
			boolean inside=false, nearby=false;
			for(int[] c : claimed){
				if(pos>=c[0] && pos<=c[1]){inside=true; break;}
				if(pos>=c[0]-nearbyPad && pos<=c[1]+nearbyPad){nearby=true;}
			}
			if(!inside && nearby){result.add(pos);}
		}
		return result.toArray();
	}

	private Orf alignWindow(String name, byte[] bases, int strand, int pass, int wStart, int wStop, int windowSource,
			int[] hitPositions, long[] hitKeys){
		if(!Float.isNaN(modelClipRescueId) && (!modelClipRescue || modelEndClipping
				|| !Float.isFinite(modelClipRescueId) || modelClipRescueId<0 || modelClipRescueId>1)){
			throw new IllegalStateException("Local-rescue identity override requires rescue mode and a finite cutoff in [0,1]: "+family);
		}
		// Length caps currently belong to the one-alignment rRNA path, where every
		// model has its own measured span before ranking or early acceptance.
		if(maxLen<Integer.MAX_VALUE && !reuseConsensusAlignment){
			throw new IllegalStateException("A finite ncRNA maxLen requires reuseConsensusAlignment: "+family);
		}
		if(pacBioConsensusAlignment && (!Euk18sRuntimeConfig.FAMILY.equals(family) || !reuseConsensusAlignment
				|| modelEndClipping || models!=null || voteEnds || rrnaEndpointFeatures!=null
				|| boundaryNetsByModel!=null || boundary5Net!=null || boundary3Net!=null || cmVerifier!=null)){
			throw new IllegalStateException("Experimental PacBio primary requires euk18S reused alignment, no HBM/NN/CM, and alignment ends: "+family);
		}
		if((modelEndClipping || modelClipRescue) && (!reuseConsensusAlignment || models!=null || voteEnds
				|| rrnaEndpointFeatures!=null || boundaryNetsByModel!=null || boundary5Net!=null || boundary3Net!=null)){
			throw new IllegalStateException("Experimental model clipping requires reused alignment, no HBM/endpoint nets, and alignment ends: "+family);
		}
		final int wLen=wStop-wStart+1;
		if(wLen<minLen){
			diagnoseWindow(name,bases,strand,pass,wStart,wStop,"REJECT_WLEN",-1,-1,-1,0f,0f,0f,-1,-1,false);
			return null;
		}
		byte[] seq=copyRegionUpper(bases, wStart, wStop+1);
		final int gateHits=kmerHits(seq), khits=lastSeedOccurrences;
		// Diagnostic khits/postSeedWindow retain occurrence units. In distinct mode a
		// REJECT_KHITS row may therefore show raw khits>=minKmerHits; controls report the mode.
		if(DEBUG){System.err.println("DEBUG alignWindow wLen="+wLen+" khits="+khits+" minKmerHits="+minKmerHits
			+" quantumThresh="+quantumThresh+" usingQuantum="+(wLen>quantumThresh));}
		if(gateHits<minKmerHits){
			diagnoseWindow(name,bases,strand,pass,wStart,wStop,"REJECT_KHITS",khits,-1,-1,0f,0f,0f,-1,-1,false);
			return null;
		}
		if(workloadSink!=null){workloadSink.postSeedWindow(name, strand, pass, wStart, wStop, khits);}
		int[] shortlist=shortlistByKmer(seq, indexTopN);
		if(DEBUG){System.err.println("DEBUG shortlist.length="+shortlist.length+" shortlist="+Arrays.toString(shortlist));}
		if(!mappingOracleExhaustive && strictIndexCutoff && kmerIndex!=null && kmerIndex.lastMaxShared()<indexMinHitsDefault){
			diagnoseWindow(name,bases,strand,pass,wStart,wStop,"REJECT_INDEX",khits,shortlist.length,-1,0f,0f,0f,-1,-1,false);
			return null;
		}
		if(reuseConsensusAlignment){
			if(modelClipRescue){
				assert(!modelEndClipping && !suppressRejectedWindow) : "Rescue must begin with unchanged Quantum and no leaked diagnostic suppression";
				final Orf accepted;
				// A scheduled window has one terminal diagnostic outcome, although
				// MODEL rows count every real alignment in both verifier phases.
				suppressRejectedWindow=true;
				try{accepted=alignWindowOnce(name,bases,strand,pass,wStart,wStop,windowSource,hitPositions,hitKeys,seq,shortlist,khits,false);}
				finally{suppressRejectedWindow=false;}
				if(accepted!=null){return accepted;}
				return alignWindowOnce(name,bases,strand,pass,wStart,wStop,windowSource,hitPositions,hitKeys,seq,shortlist,khits,true);
			}
			return alignWindowOnce(name,bases,strand,pass,wStart,wStop,windowSource,hitPositions,hitKeys,seq,shortlist,khits,modelEndClipping);
		}
		float bestId=0; int bestModel=-1;
		int bestStart=0, bestStop=wLen-1;
		final boolean quantumUsed=wLen>quantumThresh;
		final boolean quantumFeatures=(boundaryFeatureVersion==NcrnaBoundaryScorer.FEATURES_V2
			|| boundaryFeatureVersion==NcrnaBoundaryScorer.FEATURES_V3);
		final float[] shortlistIdentity=(quantumFeatures ? new float[shortlist.length] : null);
		if(quantumUsed){
			int[] pos=new int[4];
			for(int j=0; j<shortlist.length; j++){
				int m=shortlist[j];
				float id=QuantumAligner.alignStatic(library[m], seq, pos);
				diagnoseModel(name,bases.length,strand,pass,wStart,wStop,m,id);
				alignedBases+=seq.length;
				if(shortlistIdentity!=null){shortlistIdentity[j]=id;}
				alignmentCount++;
				if(modelAttemptSink!=null){
					final int alignedLength=quantumTraceAlignedLength(modelAttemptSink,pos,wLen);
					modelAttemptSink.modelAttempt(name, strand, pass, wStart, wStop,
						j+1, m, kmerIndex.lastSharedCount(m), id, alignedLength, id>=idPass);
				}
				if(DEBUG){System.err.println("DEBUG   quantum model="+m+" modelLen="+library[m].length+" id="+id+" pos="+Arrays.toString(pos));}
				if(id>bestId){bestId=id; bestModel=m; bestStart=pos[0]; bestStop=pos[1];}
				if(!mappingOracleExhaustive && rankedModelFallback && id>=idPass){break;}
			}
		}else{
			for(int j=0; j<shortlist.length; j++){
				int m=shortlist[j];
				float id=ScrabbleAligner.alignStatic(seq, library[m], null);
				diagnoseModel(name,bases.length,strand,pass,wStart,wStop,m,id);
				alignedBases+=seq.length;
				if(shortlistIdentity!=null){shortlistIdentity[j]=id;}
				alignmentCount++;
				if(modelAttemptSink!=null){modelAttemptSink.modelAttempt(name, strand, pass, wStart, wStop,
					j+1, m, kmerIndex.lastSharedCount(m), id, -1, id>=idPass);}
				if(DEBUG){System.err.println("DEBUG   scrabble model="+m+" modelLen="+library[m].length+" id="+id);}
				if(id>bestId){bestId=id; bestModel=m;}
				if(!mappingOracleExhaustive && rankedModelFallback && id>=idPass){break;}
			}
		}
		if(DEBUG){System.err.println("DEBUG best: bestId="+bestId+" bestModel="+bestModel+" idBorderline="+idBorderline+" idPass="+idPass);}
		if(bestId<idBorderline || bestModel<0){
			diagnoseWindow(name,bases,strand,pass,wStart,wStop,"REJECT_ID",khits,shortlist.length,bestModel,bestId,0f,0f,-1,-1,false);
			return null;
		}
		final int orfStart=wStart+bestStart;
		final int orfStop=wStart+bestStop;
		//Forward-ported from Noire's tree (2026-08-28, C3 merge): inclusive-length check (+1)
		//alongside the scoreA/scoreB formula below -- the OLD `orfStop-orfStart<minLen` undercounted
		//an inclusive [orfStart,orfStop] span by one base, the same off-by-one class flagged (but
		//deliberately NOT fixed) in subtractClaimed's TODO comment; this one Noire did fix, as part
		//of the same scoreA/scoreB commit that introduced orfLen.
		if(orfStop-orfStart+1<minLen){
			diagnoseWindow(name,bases,strand,pass,wStart,wStop,"REJECT_LEN",khits,shortlist.length,bestModel,bestId,0f,0f,orfStart,orfStop,false);
			return null;
		}
		Orf orf=new Orf(name, orfStart, orfStop, strand, 0, bases, false, outputType);
		orf.ncrnaFamily=family;//set once, covers the IDPASS/RESCUE/HBM accept branches below -- all reuse this same Orf instance
		final int orfLen=orfStop-orfStart+1;
		//Forward-ported from Noire's tree (2026-08-28, C3 merge): scoreA/scoreB per-family score
		//formula replaces the flat bestId*100 -- Noire's #1 recall lever (+11.5pp rnasep).
		//NcrnaScavenger-only; TrnaCaller's tRNA orfScore (bestId*100) is untouched.
		orf.orfScore=scoreA+scoreB*orfLen*bestId*bestId;

		if(cmVerifier!=null){
			boolean trimmed=false;
			if(annotate && modelNames!=null && bestModel<modelNames.length){
				orf.trnaModel=modelNames[bestModel];
				trimmed=finishAcceptedBoundary(orf, bases, bestModel, wStart, wStop, windowSource, bestId, orfLen,
					bestId, quantumUsed, hitPositions, hitKeys);
			}else{applyVotedEnds(orf, hitPositions, hitKeys, wStart, wStop, bases.length);}
			final boolean accepted=cmVerifier.verify(orf, bases, wStart, wStop, "ncrna_seed");
			if(!accepted && pendingInstrumentation!=null){pendingInstrumentation.remove(orf);}
			diagnoseWindow(name,bases,strand,pass,wStart,wStop,accepted ? "ACCEPT_CM" : "REJECT_CM",
				khits,shortlist.length,bestModel,bestId,0f,0f,orf.start,orf.stop,trimmed);
			return accepted ? orf : null;
		}
		if(bestId>=idPass){
			boolean trimmed=false;
			if(annotate && modelNames!=null && bestModel<modelNames.length){
				orf.trnaModel=modelNames[bestModel];
				trimmed=finishAcceptedBoundary(orf, bases, bestModel, wStart, wStop, windowSource, bestId, orfLen,
					bestId, quantumUsed, hitPositions, hitKeys);
			}else{applyVotedEnds(orf, hitPositions, hitKeys, wStart, wStop, bases.length);}
			diagnoseWindow(name,bases,strand,pass,wStart,wStop,"ACCEPT_IDPASS",khits,shortlist.length,bestModel,bestId,0f,0f,orf.start,orf.stop,trimmed);
			return orf;
		}else{
			byte[] orfSeq=copyRegionUpper(bases, orfStart, orfStop+1);
			if(DEBUG){System.err.println("DEBUG rescue orfStart="+orfStart+" orfStop="+orfStop
				+" orfSeq.length="+orfSeq.length+" (from Quantum pos bestStart="+bestStart+" bestStop="+bestStop+")");}
			float bestReId=0; int bestReModel=-1, bestReRank=-1;
			float bestHbm=-999, bestHbmIdentity=Float.NaN; int bestHbmModel=-1, bestHbmRank=-1;
			for(int j=0; j<shortlist.length; j++){
				final int m=shortlist[j];
				final float id=ScrabbleAligner.alignStatic(orfSeq, library[m], null);
				alignedBases+=orfSeq.length;
				alignmentCount++;
				if(DEBUG){System.err.println("DEBUG   rescue-scrabble model="+m+" modelLen="+library[m].length+" id="+id);}
				if(id>bestReId){bestReId=id; bestReModel=m; bestReRank=j;}
				if(!mappingOracleExhaustive && id>=idPass){break;}
				if(models!=null && m<models.length && id>=idBorderline){
					final float hbm=TrnaConsensusBuilder.scoreAgainstModel(orfSeq, models[m]);
					hbmScoreCalls++; hbmBasesScored+=orfSeq.length;
					if(DEBUG){System.err.println("DEBUG   rescue-hbm model="+m+" hbm="+hbm);}
					if(hbm>bestHbm){bestHbm=hbm; bestHbmModel=m; bestHbmIdentity=id; bestHbmRank=j;}
				}
			}
			if(DEBUG){System.err.println("DEBUG rescue result: bestReId="+bestReId+" bestReModel="+bestReModel
				+" bestHbm="+bestHbm+" bestHbmModel="+bestHbmModel+" idPass="+idPass+" hbmPass="+hbmPass);}
			if(bestReId>=idPass && bestReModel>=0){
				boolean trimmed=false;
				if(annotate && modelNames!=null && bestReModel<modelNames.length){
					orf.trnaModel=modelNames[bestReModel];
					trimmed=finishAcceptedBoundary(orf, bases, bestReModel, wStart, wStop, windowSource, bestReId, orfLen,
						(shortlistIdentity==null ? bestReId : shortlistIdentity[bestReRank]), quantumUsed, hitPositions, hitKeys);
				}else{applyVotedEnds(orf, hitPositions, hitKeys, wStart, wStop, bases.length);}
				diagnoseWindow(name,bases,strand,pass,wStart,wStop,"ACCEPT_RESCUE",khits,shortlist.length,bestReModel,bestId,bestReId,bestHbm,orf.start,orf.stop,trimmed);
				return orf;
			}else if(bestHbm>=hbmPass && bestHbmModel>=0){
				boolean trimmed=false;
				if(annotate && modelNames!=null && bestHbmModel<modelNames.length){
					orf.trnaModel=modelNames[bestHbmModel];
					trimmed=finishAcceptedBoundary(orf, bases, bestHbmModel, wStart, wStop, windowSource, bestHbmIdentity, orfLen,
						(shortlistIdentity==null ? bestHbmIdentity : shortlistIdentity[bestHbmRank]), quantumUsed, hitPositions, hitKeys);
				}else{applyVotedEnds(orf, hitPositions, hitKeys, wStart, wStop, bases.length);}
				diagnoseWindow(name,bases,strand,pass,wStart,wStop,"ACCEPT_HBM",khits,shortlist.length,bestHbmModel,bestId,bestReId,bestHbm,orf.start,orf.stop,trimmed);
				return orf;
			}
			diagnoseWindow(name,bases,strand,pass,wStart,wStop,"REJECT_RESCUE",khits,shortlist.length,bestReModel,bestId,bestReId,bestHbm,orfStart,orfStop,false);
		}
		return null;
	}

	/** Production rRNA path: align each shortlisted consensus exactly once and
	 * retain its Quantum traceback for detection, coordinates, endpoint features,
	 * and optional HBM rescue. */
	private Orf alignWindowOnce(String name, byte[] bases, int strand, int pass, int wStart, int wStop,
			int windowSource, int[] hitPositions, long[] hitKeys, byte[] seq, int[] shortlist, int khits, boolean clipModel){
		return alignWindowOnce(name,bases,strand,pass,wStart,wStop,windowSource,hitPositions,hitKeys,seq,shortlist,khits,clipModel,false);
	}
	private Orf alignWindowOnce(String name, byte[] bases, int strand, int pass, int wStart, int wStop,
			int windowSource, int[] hitPositions, long[] hitKeys, byte[] seq, int[] shortlist, int khits, boolean clipModel,boolean joined){
		if(trimAlignmentExtent){throw new IllegalStateException("Alignment-reuse path cannot perform a second endpoint alignment: "+family);}
		final boolean[] aligned=new boolean[library.length];
		float bestId=0, bestHbm=-999, bestHbmIdentity=Float.NaN;
		float bestPassId=0;
		int bestPassModel=-1;
		AlignmentStats bestPassStats=null;
		int bestModel=-1, bestHbmModel=-1;
		AlignmentStats bestHbmStats=null;
		for(int j=0; j<shortlist.length; j++){
			final int m=shortlist[j];
			if(m<0 || m>=library.length){throw new IllegalStateException("Invalid shortlisted model "+m+" for "+family);}
			final float modelPass=clipModel && modelClipRescue && !Float.isNaN(modelClipRescueId)
				? modelClipRescueId : modelThresholds==null ? idPass : modelThresholds.pass(m);
			final float modelBorderline=modelThresholds==null?idBorderline:modelThresholds.borderline(m);
			if(aligned[m]){throw new IllegalStateException("Consensus aligned more than once at one locus: family="+family+", model="+m);}
			aligned[m]=true;
			final AlignmentStats stats;
			//TODO: Probable endpoint bug - Quantum can lose a better 3-prime path even when the
			//window contains the full locus: euk18S NW_026953564.1 and NC_009907.1 were 155/199bp
			//short; same-model full GlocalAligner on the same containing windows restores both
			//CM ends. Reproducer: dev/euk18s_saved_window_exact_20261001.sh, source-window trace
			//seal e324e60215c630244495. No fallback here until its accuracy/cost is measured;
			//window truncation is a separate loss and must not be attributed to the aligner.
			final float id;
			if(clipModel){
				if(endClippedAligner==null){endClippedAligner=new EndClippedAligner();}
				final EndClippedAligner.Result clipped=endClippedAligner.align(library[m], seq);
				stats=clipped;
				if(diagSink!=null){diagSink.modelClip(family,name,strand,pass,wStart,wStop,m,clipped,minLen);}
				id=stats==null ? 0f : stats.identity;
			}else if(pacBioConsensusAlignment){
				if(pacBioAligner==null){pacBioAligner=new Euk18sPacBioAligner(pacBioCosts, pacBioRolling);}
				if(!pacBioAligner.costs.sameValues(pacBioCosts) || pacBioAligner.rolling!=pacBioRolling){
					throw new IllegalStateException("PacBio costs or engine changed after worker initialization: "+family);
				}
				final Euk18sPacBioAligner.Result r=pacBioAligner.align(library[m], seq, diagSink!=null);
				id=r.identity;stats=r.validBounds ? r : null;
				if(r.comparedQuantum){alignmentCount++;alignedBases+=seq.length;}
				if(diagSink!=null){diagSink.modelPacBio(family,name,strand,pass,wStart,wStop,m,r,modelPass,joined?seq.length:maxLen);}
			}else{
				stats=new AlignmentStats(true);
				id=QuantumAligner.alignAndTraceStatic(library[m], seq, stats);
			}
			diagnoseModel(name,bases.length,strand,pass,wStart,wStop,m,id);
			alignedBases+=seq.length;
			alignmentCount++;
			final int alignedLength=stats==null ? 0 : quantumAlignedLength(stats,seq.length);
			final boolean lengthAllowed=joined ? joinedLengthAllowed(stats,seq.length) : stats!=null && alignedLength<=maxLen && (!clipModel || alignedLength>=minLen);
			if(modelAttemptSink!=null){
				modelAttemptSink.modelAttempt(name,strand,pass,wStart,wStop,j+1,m,
					kmerIndex.lastSharedCount(m),id,alignedLength,id>=modelPass && lengthAllowed);
			}
			// An overlong high-identity model must not suppress a later valid model.
			if(!lengthAllowed){continue;}
			if(DEBUG){System.err.println("DEBUG   single-quantum model="+m+" modelLen="+library[m].length
				+" id="+id+" rStart="+stats.rStart+" rStop="+stats.rStop);}
			if(id>bestId){bestId=id;bestModel=m;}
			final float candidateFloor=cmVerifier==null ? modelPass : modelBorderline;
			if(id>=candidateFloor && id>bestPassId){bestPassId=id;bestPassModel=m;bestPassStats=stats;}
			if(cmVerifier==null && models!=null && m<models.length && id>=modelBorderline && hbmPass<=1f){
				final float hbm=TrnaConsensusBuilder.scoreAlignedAgainstModel(seq,stats,models[m]);
				hbmScoreCalls++; hbmBasesScored+=alignedLength;
				if(DEBUG){System.err.println("DEBUG   reused-quantum-hbm model="+m+" hbm="+hbm);}
				if(hbm>bestHbm){bestHbm=hbm;bestHbmModel=m;bestHbmIdentity=id;bestHbmStats=stats;}
			}
			if(!mappingOracleExhaustive && rankedModelFallback && id>=modelPass){break;}
		}
		final int model;
		final float identity;
		final AlignmentStats stats;
		if(bestPassModel>=0){model=bestPassModel;identity=bestPassId;stats=bestPassStats;}
		else if(bestHbm>=hbmPass && bestHbmModel>=0){model=bestHbmModel;identity=bestHbmIdentity;stats=bestHbmStats;}
		else{
			diagnoseWindow(name,bases,strand,pass,wStart,wStop,"REJECT_VERIFY",khits,shortlist.length,bestModel,bestId,0f,bestHbm,-1,-1,false);
			return null;
		}
		final int orfStart=wStart+stats.rStart, orfStop=wStart+stats.rStop;
		if(orfStop-orfStart+1<minLen){
			diagnoseWindow(name,bases,strand,pass,wStart,wStop,"REJECT_LEN",khits,shortlist.length,model,bestId,0f,bestHbm,orfStart,orfStop,false);
			return null;
		}
		final Orf orf=new Orf(name,orfStart,orfStop,strand,0,bases,false,outputType);
		orf.ncrnaFamily=family;
		final int orfLen=orfStop-orfStart+1;
		orf.orfScore=scoreA+scoreB*orfLen*identity*identity;
		boolean trimmed=false;
		if(annotate && modelNames!=null && model<modelNames.length){
			orf.trnaModel=modelNames[model];
			trimmed=finishAcceptedBoundary(orf,bases,model,wStart,wStop,windowSource,identity,orfLen,
				identity,true,hitPositions,hitKeys);
		}else{applyVotedEnds(orf,hitPositions,hitKeys,wStart,wStop,bases.length);}
		if(!joined && orf.stop-orf.start+1>maxLen){
			if(pendingInstrumentation!=null){pendingInstrumentation.remove(orf);}
			diagnoseWindow(name,bases,strand,pass,wStart,wStop,"REJECT_MAXLEN",khits,shortlist.length,
				model,bestId,0f,bestHbm,orf.start,orf.stop,trimmed);
			return null;
		}
		if(cmVerifier!=null){
			final boolean accepted=cmVerifier.verify(orf, bases, wStart, wStop, "ncrna_reused",
				modelNames==null ? null : modelNames[model], stats);
			if(!accepted && pendingInstrumentation!=null){pendingInstrumentation.remove(orf);}
			diagnoseWindow(name,bases,strand,pass,wStart,wStop,accepted ? "ACCEPT_CM" : "REJECT_CM",
				khits,shortlist.length,model,bestId,0f,bestHbm,orf.start,orf.stop,trimmed);
			return accepted ? orf : null;
		}
		diagnoseWindow(name,bases,strand,pass,wStart,wStop,bestPassModel>=0 ? "ACCEPT_IDPASS" : "ACCEPT_HBM",
			khits,shortlist.length,model,bestId,0f,bestHbm,orf.start,orf.stop,trimmed);
		return orf;
	}
	static boolean joinedLengthAllowed(AlignmentStats stats,int windowLength){
		if(!(stats instanceof Euk18sPacBioAligner.Result)){return false;}
		return JoinedRnaTrace.valid((Euk18sPacBioAligner.Result)stats,windowLength);
	}

	/** Initial model observations do not retain sequences or request extra alignments. */
	private void diagnoseModel(String name,int length,int strand,int pass,int start,int stop,int model,float identity){
		if(diagSink==null){return;}
		final String modelName=modelNames==null ? null : modelNames[model];
		diagSink.model(diagFamily,name,strand,length,pass,start,stop,model,modelName,identity);
	}

	/** Read-only reporting helper: no model-name lookup or output work when the sink is off. */
	private void diagnoseWindow(String name,byte[] bases,int strand,int pass,int start,int stop,
			String outcome,int hits,int shortlistSize,int model,float id,float reId,float hbm,
			int orfStart,int orfStop,boolean trimmed){
		if(diagSink==null || (suppressRejectedWindow && outcome.startsWith("REJECT"))){return;}
		final String modelName=(modelNames!=null && model>=0 && model<modelNames.length ? modelNames[model] : null);
		diagSink.window(diagFamily,name,strand,bases.length,pass,start,stop,outcome,hits,shortlistSize,
			model,modelName,id,reId,hbm,orfStart,orfStop,trimmed);
	}

	/** Sink gate kept explicit so trace-only coordinate validation cannot alter the
	 * default caller path.  The null-sink return is a fixture-visible proof that
	 * malformed diagnostic coordinates remain irrelevant while tracing is off. */
	static int quantumTraceAlignedLength(NcrnaModelAttemptInstrumentSink sink,int[] pos,int windowLength){
		return sink==null ? -1 : quantumAlignedLength(pos,windowLength);
	}

	/** Converts QuantumAligner's inclusive reference coordinates into the diagnostic
	 * model-attempt span. Invalid coordinates would corrupt the trace used to compare
	 * accepted and rejected attempts, so fail loudly rather than emit a plausible value. */
	static int quantumAlignedLength(int[] pos,int windowLength){
		if(pos==null || pos.length<2 || pos[0]<0 || pos[1]<pos[0] || pos[1]>=windowLength){
			throw new IllegalStateException("Invalid Quantum alignment coordinates for model-attempt trace: "
				+Arrays.toString(pos)+", windowLength="+windowLength);
		}
		return pos[1]-pos[0]+1;
	}

	static int quantumAlignedLength(AlignmentStats stats,int windowLength){
		if(stats==null || stats.rStart<0 || stats.rStop<stats.rStart || stats.rStop>=windowLength){
			throw new IllegalStateException("Invalid Quantum traceback coordinates for model attempt: "
				+(stats==null ? "null" : stats.rStart+"-"+stats.rStop)+", windowLength="+windowLength);
		}
		return stats.rStop-stats.rStart+1;
	}

	/** Preserves the selected anchor order explicitly: alignment trim, optional
	 * hybrid5 voting, then NN refinement around the resulting per-endpoint anchor. */
	private boolean finishAcceptedBoundary(Orf orf, byte[] bases, int model, int wStart, int wStop,
			int windowSource, float acceptedIdentity, int alignedLength, float locusAni, boolean locusAniFromQuantum,
			int[] hitPositions, long[] hitKeys){
		if(rrnaEndpointFeatures!=null){
			// Resource-sized endpoint nets observe the selected raw anchor. Voting
			// must happen before features/NN, never overwrite an NN-refined result.
			applyVotedEnds(orf,hitPositions,hitKeys,wStart,wStop,bases.length);
			captureRrnaEndpointFeatures(orf,bases,model);
			if(instrumentSink!=null){captureInstrumentation(orf,bases,model,wStart,wStop,false,false,windowSource,acceptedIdentity,alignedLength);}
			return false;
		}
		final boolean trimSucceeded;
		if(trimAlignmentExtent){trimSucceeded=trimToAlignmentExtent(orf, bases, model, wStart, wStop);}
		else{trimSucceeded=false;}
		final boolean hasBoundaryNets=(boundaryNetsByModel!=null || (boundary5Net!=null && boundary3Net!=null));
		final boolean nnInvoked=(trimSucceeded || boundaryOnRawEndpoints) && hasBoundaryNets;
		if(instrumentSink!=null){captureInstrumentation(orf,bases,model,wStart,wStop,trimSucceeded,nnInvoked,windowSource,acceptedIdentity,alignedLength);}
		final int alignedStart=orf.start, alignedStop=orf.stop;
		if(nnInvoked && !nnCentre5Vote && !nnCentre3Vote){refineBoundaryNN(orf, bases, model, locusAni, locusAniFromQuantum, 0);}
		final int votedMask=applyVotedEnds(orf, hitPositions, hitKeys, wStart, wStop, bases.length);
		if(nnInvoked && nnCentre5Vote && nnCentre3Vote){refineBoundaryNN(orf, bases, model, locusAni, locusAniFromQuantum, votedMask);}
		else if(nnInvoked && nnCentre5Vote!=nnCentre3Vote){
			final int centreMask=selectMixedNnAnchors(orf,alignedStart,alignedStop,nnCentre5Vote,nnCentre3Vote,votedMask);
			refineBoundaryNN(orf, bases, model, locusAni, locusAniFromQuantum, centreMask);
		}
		return trimSucceeded;
	}

	/** Restores the alignment-centred endpoint of a mixed-centre experiment after
	 * voting has supplied the other endpoint. Returns the per-end vote provenance
	 * mask consumed by boundary tracing. */
	static int selectMixedNnAnchors(Orf orf,int alignedStart,int alignedStop,
			boolean centre5Vote,boolean centre3Vote,int votedMask){
		assert(centre5Vote!=centre3Vote) : "Mixed-anchor helper requires exactly one vote-centred endpoint";
		if(!centre5Vote){orf.start=alignedStart;}if(!centre3Vote){orf.stop=alignedStop;}
		return (centre5Vote ? votedMask&1 : 0)|(centre3Vote ? votedMask&2 : 0);
	}

	/**
	 * Generic boundary trim: aligns an extended window around the candidate to
	 * the winning model with traceback, and snaps the ORF boundaries to the
	 * alignment's own extent (rStart/rStop) directly -- no acceptor-stem
	 * search, no anticodon extraction (both tRNA-specific; stay in
	 * TrnaCaller). Mirrors TrnaCaller.trimOrf's "window within consensus
	 * span" guard.
	 */
	/** Package-visible (was private) so NcrnaBoundaryInstrumentSinkTest can directly construct
	 * the trim-no-op case (extended window shorter than the model) without needing to coax it
	 * out of the full scavenge() pipeline's alignment dynamics -- deterministic and fast versus
	 * fragile fixture engineering for the same real code path. */
	boolean trimToAlignmentExtent(Orf orf, byte[] bases, int model, int wStart, int wStop){
		final int xFrom=Tools.max(0, orf.start-trimExt);
		final int xTo=Tools.min(bases.length-1, orf.stop+trimExt);
		byte[] seqX=copyRegionUpper(bases, xFrom, xTo+1);
		final byte[] cons=library[model];
		AlignmentStats stats=new AlignmentStats(true);
		stats.doTrace=true;
		if(seqX.length>=cons.length){ScrabbleAligner.alignAndTraceStatic(cons, seqX, stats);}
		else{ScrabbleAligner.alignAndTraceStatic(seqX, cons, stats);}
		alignedBases+=seqX.length;
		alignmentCount++;
		//Citan, 2026-08-28: the original code `return`ed here on either guard, which ALSO
		//skipped instrumentation capture for these loci -- silently undercounting the accepted
		//denominator for exactly the accepted-but-unrefined cases. Restructured to a boolean so
		//capture always fires once per accepted locus (see below). trimSucceeded is the
		//STRUCTURAL eligibility signal for downstream bootstrap work: true means orf.start/
		//orf.stop below are a real post-trim position (usable, pending the driver's own
		//sweep-reachability check against truth) -- false means they are the untouched raw
		//alignWindow span and must not be used as if they were a trim result.
		boolean trimSucceeded=false;
		if(stats.matchString!=null && seqX.length>=cons.length){
			final int rStart=stats.rStart, rStop=stats.rStop;
			if(rStart>=0 && rStop>rStart && rStop<seqX.length){
				orf.start=xFrom+rStart;
				orf.stop=xFrom+rStop;
				trimSucceeded=true;
			}
		}
		return trimSucceeded;
	}

	/** Builds the private window copy (never a live bases[] reference -- see
	 * NcrnaBoundaryInstrumentSink's thread-safety javadoc) and invokes the sink. Split out of
	 * trimToAlignmentExtent so the hot (instrumentation-off) path is a single null-check with
	 * no other cost. */
	void captureInstrumentation(Orf orf, byte[] bases, int model, int wStart, int wStop,
			boolean trimSucceeded, boolean nnInvoked, int windowSource, float acceptedIdentity, int alignedLength){
		final int copyFrom=Tools.max(0, orf.start-INSTRUMENT_CAPTURE_PAD);
		final int copyTo=Tools.min(bases.length-1, orf.stop+INSTRUMENT_CAPTURE_PAD);
		final byte[] windowCopy=Arrays.copyOfRange(bases, copyFrom, copyTo+1);
		final InstrumentationCapture capture=new InstrumentationCapture(orf.scafName,orf.strand,model,wStart,wStop,
			orf.start,orf.stop,windowCopy,copyFrom,trimSucceeded,nnInvoked,windowSource,acceptedIdentity,alignedLength);
		final InstrumentationCapture old=pendingInstrumentation.put(orf,capture);
		assert(old==null) : "Each accepted candidate may queue exactly one boundary capture; guards duplicate provenance for "+orf.scafName;
	}

	/** Emits capture only after the candidate becomes a returned call. Snapshot refresh may
	 * align several mutually overlapping accepted candidates, so firing inside alignWindow
	 * overcounts the final call denominator and breaks trace/call conservation. */
	private void emitInstrumentation(Orf orf){
		if(instrumentSink==null){return;}
		final InstrumentationCapture c=pendingInstrumentation.remove(orf);
		assert(c!=null) : "Every returned accepted call must retain its queued boundary capture; guards trace/call conservation for "+orf.scafName;
		instrumentSink.capture(c.contig,c.strand,c.model,c.wStart,c.wStop,c.postStart,c.postStop,
			c.window,c.windowOffset,c.trimSucceeded,c.nnInvoked,c.windowSource,c.identity,c.alignedLength);
	}

	/** Emits evidence for accepted candidates rejected by the snapshot resolver, then drops
	 * their queued state. The committed overlapping winner is chosen by maximum overlap,
	 * then candidate score/order, solely to make the resolution link deterministic. */
	private void emitResolvedInstrumentation(ArrayList<Orf> candidates,ArrayList<Orf> selected){
		if(pendingInstrumentation==null){return;}
		final IdentityHashMap<Orf,Boolean> kept=new IdentityHashMap<Orf,Boolean>();for(Orf orf:selected){kept.put(orf,Boolean.TRUE);}
		for(Orf loser:candidates){if(kept.containsKey(loser)){continue;}final InstrumentationCapture c=pendingInstrumentation.remove(loser);assert(c!=null) : "Every resolved-out candidate must retain queued capture evidence for "+loser.scafName;
			final Orf winner=resolvedComponentWinner(loser,candidates,kept);
			if(winner==null){throw new IllegalStateException("Resolved-out candidate lacks a committed winner in its reciprocal-overlap component: "+loser.scafName+":"+loser.start+"-"+loser.stop);}
			instrumentSink.resolvedOut(loser,winner,c.model,c.wStart,c.wStop,c.postStart,c.postStop,c.window,c.windowOffset,c.trimSucceeded,c.nnInvoked,c.windowSource,c.identity,c.alignedLength);
		}
	}

	/** Finds the one committed winner in loser's reciprocal-overlap component.
	 * Component traversal preserves evidence linkage even for an A-B-C conflict
	 * chain whose maximum-score winner does not directly overlap every loser. */
	private static Orf resolvedComponentWinner(Orf loser,ArrayList<Orf> candidates,IdentityHashMap<Orf,Boolean> kept){
		final IdentityHashMap<Orf,Boolean> seen=new IdentityHashMap<Orf,Boolean>();final ArrayList<Orf> queue=new ArrayList<Orf>();queue.add(loser);seen.put(loser,Boolean.TRUE);
		for(int q=0;q<queue.size();q++){final Orf current=queue.get(q);if(kept.containsKey(current)){return current;}for(Orf candidate:candidates){if(!seen.containsKey(candidate)&&snapshotConflict(current,candidate)){seen.put(candidate,Boolean.TRUE);queue.add(candidate);}}}
		return null;
	}

	/** Applies the boundary-precision NN's refinement to an already-trimmed, already-verified
	 * Orf. Builds a window large enough for the configured sweep and tip-feature profile,
	 * then defers to NcrnaBoundaryScorer.refineBoundariesWithCounts.
	 * No-op (leaves orf untouched) if the padded window can't hold a valid base candidate. */
	private void refineBoundaryNN(Orf orf, byte[] bases, int model, float locusAni, boolean locusAniFromQuantum, int votedMask){
		final boolean dispatched=(boundaryNetsByModel!=null);
		final CellNet startNet=(dispatched ? boundaryNetsByModel[model] : boundary5Net);
		final CellNet stopNet=(dispatched ? boundaryNetsByModel[model] : boundary3Net);
		final TrnaBoundaryFeatures.NinemerTable startTable=(dispatched ? boundaryStartTablesByModel[model] : boundaryStartTable);
		final TrnaBoundaryFeatures.NinemerTable stopTable=(dispatched ? boundaryStopTablesByModel[model] : boundaryStopTable);
		final int PAD=boundaryRefinementPad(boundaryStartOffsets, boundaryStopOffsets,
			boundaryStartInside, boundaryStartOutside, boundaryStopInside, boundaryStopOutside);
		final int winStart=Tools.max(0, orf.start-PAD);
		final int winStop=Tools.min(bases.length-1, orf.stop+PAD);
		final byte[] window=copyRegionUpper(bases, winStart, winStop+1);
		final int s=orf.start-winStart, e=orf.stop-winStart;
		if(s<0 || e>=window.length || e-s<15){return;}
		final float contigGC=contigGC(bases);
		final BaseGraph modelGraph=(models!=null && model<models.length ? models[model] : null);
		final float[] startFuzz,stopFuzz;
		if(boundaryFeatureVersion==NcrnaBoundaryScorer.FEATURES_V2 || boundaryFeatureVersion==NcrnaBoundaryScorer.FEATURES_V3){
			if(!locusAniFromQuantum){throw new IllegalStateException("boundaryfeatures=v"+boundaryFeatureVersion
				+" requires a QuantumAligner locus identity; family="+family);}
		}
		if(boundaryFeatureVersion==NcrnaBoundaryScorer.FEATURES_V2){
			if(boundaryFuzzConstants==null){throw new IllegalStateException("boundaryfeatures=v2 requires loaded model/end fuzz constants; family="+family);}
			startFuzz=boundaryFuzzConstants.values(model,true);
			stopFuzz=boundaryFuzzConstants.values(model,false);
		}else{startFuzz=null;stopFuzz=null;}
		final NcrnaBoundaryScorer.RefinementResult result;
		final long boundaryStartNanos=System.nanoTime();
		try{
			result=NcrnaBoundaryScorer.refineBoundariesWithCounts(
				startNet, stopNet, window, s, e,
				library[model], modelGraph, startTable, stopTable,
				boundaryStartInside, boundaryStartOutside, boundaryStopInside, boundaryStopOutside,
				contigGC, boundaryMeanLen, boundaryStartOffsets, boundaryStopOffsets,
				boundaryMarginStart, boundaryMarginStop, winStart==0, winStop==bases.length-1,
				boundaryFeatureVersion,locusAni,startFuzz,stopFuzz,
				(dispatched ? boundaryStartTable : null), (dispatched ? boundaryStopTable : null));
		}finally{boundaryNanos+=System.nanoTime()-boundaryStartNanos;}
		final int anchorStart=orf.start, anchorStop=orf.stop;
		final int startOffset=applyBoundaryScoreCutoff(result.startOffset,result.startChosenScore,nnCutoff);
		final int stopOffset=applyBoundaryScoreCutoff(result.stopOffset,result.stopChosenScore,nnCutoff);
		orf.start+=startOffset;
		orf.stop+=stopOffset;
		if(boundaryScoreSink!=null){boundaryScoreSink.scored(orf.scafName, orf.strand, model,
			anchorStart, anchorStop, startOffset, stopOffset,
			result.startSitesScored, result.stopSitesScored,
			result.startChosenScore, result.stopChosenScore,
			NcrnaBoundaryScorer.centeredRadius(boundaryStartOffsets),
			NcrnaBoundaryScorer.centeredRadius(boundaryStopOffsets),
			(votedMask&1)!=0, (votedMask&2)!=0);}
	}

	static int applyBoundaryScoreCutoff(int offset,float score,float cutoff){return Float.isNaN(cutoff)||score>=cutoff ? offset : 0;}

	/** Enough genomic context for every configured candidate and every five-position
	 * enrichment profile.  The prior fixed PAD=10 clipped the measured G16 radii
	 * 12--25 and could silently replace real tip features with out-of-bounds zeros. */
	static int boundaryRefinementPad(int[] startOffsets, int[] stopOffsets,
			int startInside, int startOutside, int stopInside, int stopOutside){
		if(startOffsets==null || startOffsets.length<1 || stopOffsets==null || stopOffsets.length<1){
			throw new IllegalArgumentException("Boundary refinement requires nonempty start and stop offset arrays");
		}
		if(startInside<0 || startOutside<0 || stopInside<0 || stopOutside<0){
			throw new IllegalArgumentException("Boundary refinement requires nonnegative table geometry");
		}
		int maxOffset=0;
		for(int x : startOffsets){maxOffset=Tools.max(maxOffset, Math.abs(x));}
		for(int x : stopOffsets){maxOffset=Tools.max(maxOffset, Math.abs(x));}
		final int maxK=Tools.max(startInside+startOutside, stopInside+stopOutside);
		return Tools.max(10, maxOffset+2+maxK);
	}

	/** Returns an uppercase private copy for ncRNA matching and refinement without
	 * modifying the shared, potentially soft-masked genome sequence. */
	private static byte[] copyRegionUpper(byte[] bases, int from, int to){
		final byte[] copy=Arrays.copyOfRange(bases, from, to);
		Tools.toUpperCase(copy);
		return copy;
	}

	/** Per-contig GC cache (identity-keyed on the bases[] reference), mirrors TrnaCaller's own
	 * contigGC cache exactly -- a NcrnaScavenger instance processes one contig/strand's bases[]
	 * across many calls within one scavenge() invocation, so recomputing GC from scratch per
	 * call would rescan the whole contig once per locus. */
	private float contigGC(byte[] bases){
		if(bases!=gcCacheBases){gcCacheValue=shared.Tools.calcGC(bases); gcCacheBases=bases;}
		return gcCacheValue;
	}
	private byte[] gcCacheBases=null;
	private float gcCacheValue=0;

	/** One scan records raw occurrences for diagnostics and optionally distinct matching keys for the gate. */
	private int kmerHits(byte[] seq){
		lastSeedOccurrences=0;
		if(seedDistinct && distinctSeedKeys!=null){distinctSeedKeys.clear();}
		if(kmerSet==null){lastSeedOccurrences=Integer.MAX_VALUE;return Integer.MAX_VALUE;}
		if(kLong<=0 || kLong>31 || seq.length<kLong){return 0;}
		if(seedDistinct && distinctSeedKeys==null){distinctSeedKeys=new LongHashSet(16);}
		final long kmask=~((-1L)<<(2*kLong));
		final byte[] bton=AminoAcid.baseToNumber;
		long kmer=0; int len=0, hits=0, distinct=0;
		for(int i=0; i<seq.length; i++){
			final int x=bton[seq[i]];
			if(x>=0){
				kmer=((kmer<<2)|x)&kmask; len++;
				if(len>=kLong && kmerSet.contains(kmer)){
					hits++;
					if(seedDistinct && distinctSeedKeys.add(kmer)){distinct++;}
				}
			}else{len=0; kmer=0;}
		}
		lastSeedOccurrences=hits;
		assert(distinct<=hits) : shared.KillSwitch.assertDie("Distinct matching seed keys cannot outnumber raw occurrences from the same window scan");
		return seedDistinct ? distinct : hits;
	}

	private int[] shortlistByKmer(byte[] seq, int topN){
		if(mappingOracleExhaustive || kmerIndex==null || library==null){
			//Refresh diagnostic shared-count scratch for this query, but discard its filtering.
			if(mappingOracleExhaustive && kmerIndex!=null){kmerIndex.shortlist(seq,library.length);}
			int[] all=new int[library.length];
			for(int i=0; i<all.length; i++){all[i]=i;}
			return all;
		}
		return kmerIndex.shortlist(seq, topN, indexScoreMargin);
	}

	private boolean votingEnabled(){return voteTable!=null && (voteWindows || voteEnds);}

	private void applyVoteWindows(ArrayList<int[]> windows, int[] positions, long[] keys, int seqLen){
		for(int[] window : windows){
			final int from=window[0], to=window[1];
			if(voteTrace!=null){voteTrace.resetVote();}
			final boolean tooFew=countHits(positions, from, to)<minKmerHits;
			if(tooFew || !computeVotes(positions, keys, from, to)){
				if(voteTrace!=null){voteTrace.legacy(from, to, from, to, "retained_fallback",
					tooFew ? "TOO_FEW_SEEDS" : voteTable==null ? "NO_VOTE_TABLE"
					: voteTrace.hasWeightedVote() ? "INVALID_VOTE_SPAN" : "NO_TRAINED_SEEDS");}
				window[2]=WINDOW_FALLBACK; voteWindowFallbacks++; continue;
			}
			final int start=Tools.max(0, (int)Math.round(voteScratch[0])-voteSlack);
			final int stop=Tools.min(seqLen-1, (int)Math.round(voteScratch[1])+voteSlack);
			if(stop-start+1<minLen){
				if(voteTrace!=null){voteTrace.legacy(from, to, from, to, "retained_fallback", "SHORT_VOTE");}
				window[2]=WINDOW_FALLBACK; voteWindowFallbacks++; continue;
			}
			if(voteTrace!=null){voteTrace.legacy(from, to, start, stop, "retained_vote", "VOTE_RETAINED");}
			window[0]=start; window[1]=stop; window[2]=WINDOW_VOTED; votedWindows++;
		}
	}

	static String windowSourceName(int source){return source==WINDOW_VOTED ? "voted" : source==WINDOW_FALLBACK ? "fallback" : source==WINDOW_JOINED ? "joined" : "collapse";}
	private static int windowSource(int[] window){return window.length>2 ? window[2] : WINDOW_COLLAPSE;}

	private int applyVotedEnds(Orf orf, int[] positions, long[] keys, int from, int to, int seqLen){
		if(!voteEnds){return 0;}
		if(!computeVotes(positions, keys, from, to)){
			voteEndFallbackStart++; voteEndFallbackStop++; return 0;
		}
		int start=orf.start, stop=orf.stop;
		final boolean useStart=voteScratch[2]<=voteEndsMaxSd, useStop=voteScratch[3]<=voteEndsMaxSd;
		if(useStart){start=(int)Math.round(voteScratch[0]);}else{voteEndFallbackStart++;}
		if(useStop){stop=(int)Math.round(voteScratch[1]);}else{voteEndFallbackStop++;}
		if(start>=0 && stop>=start && stop<seqLen){
			orf.start=start; orf.stop=stop; if(useStart){votedStarts++;} if(useStop){votedStops++;}
			return (useStart ? 1 : 0)|(useStop ? 2 : 0);
		}else{
			voteEndInvalidSpans++;
			if(useStart){voteEndFallbackStart++;}
			if(useStop){voteEndFallbackStop++;}
			return 0;
		}
	}

	/** Task 16-selected W2-D25 policy, isolated so the combiner remains swappable. */
	static double voteWeight(double sd, double distance){
		return 1.0/(Math.max(1.0, sd*sd)*(1.0+distance/25.0));
	}

	private boolean computeVotes(int[] positions, long[] keys, int from, int to){
		if(voteTable==null || positions==null || keys==null || positions.length!=keys.length){return false;}
		double leftSum=0, leftWeight=0, rightSum=0, rightWeight=0; int topLeft=0, topRight=0;
		Arrays.fill(topLeftWeights, -1); Arrays.fill(topRightWeights, -1);
		for(int i=0; i<positions.length; i++){
			final int center=positions[i]; if(center<from || center>to){continue;}
			final SeedOffsetTable.KmerInfo info=voteTable.get(keys[i]); if(info==null || !info.trained()){continue;}
			final double predLeft=center-kLong/2-info.leftOffset;
			final double predRight=center+kLong/2+info.rightOffset;
			final double wl=voteWeight(info.leftSD, info.leftOffset), wr=voteWeight(info.rightSD, info.rightOffset);
			leftSum+=predLeft*wl; leftWeight+=wl; rightSum+=predRight*wr; rightWeight+=wr;
			if(voteTrace!=null){voteTrace.contribution(predLeft, predRight);}
			topLeft=insertTop(predLeft, wl, topLeftValues, topLeftWeights, topLeft);
			topRight=insertTop(predRight, wr, topRightValues, topRightWeights, topRight);
		}
		if(!(leftWeight>0) || !(rightWeight>0)){return false;}
		voteScratch[0]=leftSum/leftWeight; voteScratch[1]=rightSum/rightWeight;
		voteScratch[2]=populationSd(topLeftValues, topLeft); voteScratch[3]=populationSd(topRightValues, topRight);
		if(voteTrace!=null){voteTrace.weighted(voteScratch[0], voteScratch[1], leftWeight, rightWeight, voteScratch[2], voteScratch[3]);}
		return voteScratch[1]>=voteScratch[0];
	}

	private static int insertTop(double value, double weight, double[] values, double[] weights, int size){
		int pos=Math.min(size, 10); while(pos>0 && weight>weights[pos-1]){pos--;}
		if(pos<10){for(int j=Math.min(size,9); j>pos; j--){weights[j]=weights[j-1]; values[j]=values[j-1];}
			weights[pos]=weight; values[pos]=value; if(size<10){size++;}}
		return size;
	}

	private static double populationSd(double[] values, int size){
		if(size<1){return Double.POSITIVE_INFINITY;} double mean=0,m2=0;
		for(int i=0;i<size;i++){double d=values[i]-mean;mean+=d/(i+1);m2+=d*(values[i]-mean);}
		return Math.sqrt(Math.max(0,m2/size));
	}

	private static int countHits(int[] positions, int from, int to){int n=0;for(int p:positions){if(p>=from&&p<=to){n++;}}return n;}

	static long[] subsetKeys(int[] allPositions, long[] allKeys, int[] subset){
		final long[] out=new long[subset.length]; int j=0;
		for(int i=0;i<allPositions.length&&j<subset.length;i++){if(allPositions[i]==subset[j]){out[j++]=allKeys[i];}}
		if(j!=subset.length){throw new IllegalStateException("Could not match pass-2 seed keys");}
		return out;
	}

	/*--------------------------------------------------------------*/

	private final byte[][] library;
	private final BaseGraph[] models;
	private final String[] modelNames;
	private final boolean annotate;
	private final LongHashSet kmerSet;
	private final int kLong;
	private final int minLen;
	private final int windowPad;
	private final TrnaKmerIndex kmerIndex;

	private static final int[] EMPTY=new int[0];
	static final int WINDOW_COLLAPSE=0, WINDOW_VOTED=1, WINDOW_FALLBACK=2, WINDOW_JOINED=3;
	Euk18sJoinedProposals joinedProposals=null;
	SeedModelPositionTable joinedPositions=null;
	int joinedSideSupport=1;
	private SeedInsertionWindowBuilder joinedBuilder;

	//Per-family tunables not yet wired to per-family values (Noire, 2026-08-23: "start
	//with tRNA's defaults, tune per family later") -- default to TrnaCaller's measured-best
	//tRNA config. Mutable instance fields (not static): unlike windowPad/minLen/kLong, these
	//don't yet have measured per-family values, but must stay per-instance since coexisting
	//scavengers for different families will eventually need different ones.
	int indexK=7;
	int indexTopN=60;
	int indexMinHitsDefault=12;
	boolean adaptiveMinHits=true;
	float adaptFloor=11;
	float adaptTopFrac=0.48f;
	float adaptQFrac=0.072f;
	int minKmerHits=1;
	/** Opt-in gate policy; raw occurrence diagnostics and window construction are unchanged. */
	boolean seedDistinct=false;
	/** Private scratch reused across candidate windows; absent while distinct mode is unused. */
	private LongHashSet distinctSeedKeys=null;
	private int lastSeedOccurrences=0;
	int maxLen=Integer.MAX_VALUE;
	/** Production controls; claimed-window refresh prevents duplicate calls from stale window geometry. */
	boolean refreshClaimedWindows=true;
	boolean rankedModelFallback=false;
	boolean strictIndexCutoff=false;
	/** -1 preserves topN; otherwise retain every positive score within this gap. */
	int indexScoreMargin=-1;
	boolean trimAlignmentExtent=true;
	int outputType=ProkObject.RNA;
	/** Explicit experiment-only permission to refine QuantumAligner endpoints when
	 * traceback trimming is disabled. Default false preserves every existing caller. */
	boolean boundaryOnRawEndpoints=false;
	boolean nnCentre5Vote=false, nnCentre3Vote=false;
	float nnCutoff=Float.NaN;
	SeedOffsetTable voteTable=null;
	boolean voteWindows=false, voteEnds=false;
	boolean voteBeforePadding=false;
	private SeedVoteWindowBuilder voteWindowBuilder=null;
	private int fallbackCoreLength;
	int voteSlack=60;
	float voteEndsMaxSd=5f;
	long votedWindows=0, voteWindowFallbacks=0, votedStarts=0, votedStops=0;
	long voteEndFallbackStart=0, voteEndFallbackStop=0, voteEndInvalidSpans=0;
	private final double[] voteScratch=new double[4], topLeftValues=new double[10], topLeftWeights=new double[10];
	private final double[] topRightValues=new double[10], topRightWeights=new double[10];
	//Forward-ported from Noire's tree (2026-08-28, C3 merge): idPass/idBorderline are now
	//FINAL, set unconditionally by every constructor (0.75f/0.65f remain the fallback values
	//for the 8-arg and 15-arg delegating overloads) -- previously mutable with the same
	//defaults but never actually reassigned outside a constructor, so this is a safety
	//tightening, not a behavior change for any existing caller.
	final float idPass;
	final float idBorderline;
	private NcrnaModelThresholds modelThresholds=null;
	void setModelThresholds(NcrnaModelThresholds thresholds){
		if(thresholds==null || !reuseConsensusAlignment){throw new IllegalArgumentException("Model cutoffs require the single-alignment path");}
		thresholds.validate(modelNames);modelThresholds=thresholds;
	}
	//Forward-ported from Noire's tree (2026-08-28, C3 merge): per-family orfScore formula
	//constants -- see the 15-arg constructor's javadoc for the full rationale.
	final float scoreA;
	final float scoreB;
	float hbmPass=0.75f;
	int quantumThresh=120;
	int nearbyPad=200;
	float collapseFrac=0.9f;
	int trimExt=10;
	boolean scavengePass2=true;
	boolean reuseConsensusAlignment=false;
	boolean pacBioConsensusAlignment=false;
	boolean pacBioRolling=false;
	align2.PacBioScoreParameters pacBioCosts=align2.PacBioScoreParameters.DEFAULT;
	private Euk18sPacBioAligner pacBioAligner;
	boolean modelEndClipping=false;
	boolean modelClipRescue=false;
	/** NaN preserves the per-model pass cutoff; Quantum never consumes this override. */
	float modelClipRescueId=Float.NaN;
	private boolean suppressRejectedWindow=false;
	private EndClippedAligner endClippedAligner;
	/** Explicit family identity (e.g. "s18", "r58", "lsu"), set post-construction by
	 * GeneCaller.makeRnas() from NcrnaFamily.name -- mirrors the existing hbmPass/
	 * collapseFrac post-construction assignment pattern above. Propagated onto every
	 * generic ncRNA Orf this scavenger creates (Ganyu/Qiqi design, 2026-09-17). */
	String family;
	CmRnaVerifier cmVerifier;

	private long alignmentCount=0;
	private long alignedBases=0, hbmScoreCalls=0, hbmBasesScored=0;
	private long kmerHitCount=0, windowCount=0;
	private long boundaryNanos=0;

	//C3 boundary-precision-NN resources (G11, 2026-08-28) -- boundary5Net==null (this
	//instance's default unless the full constructor is used with real templates) means OFF,
	//structurally: refineBoundaryNN is never called (see trimToAlignmentExtent's pre-call
	//guard) and these fields are never read. Per-instance CLONES of the family's shared
	//read-only templates (thread safety -- see the full constructor's javadoc).
	private final CellNet boundary5Net, boundary3Net;
	private final TrnaBoundaryFeatures.NinemerTable boundaryStartTable, boundaryStopTable;
	private CellNet[] boundaryNetsByModel=null;
	private TrnaBoundaryFeatures.NinemerTable[] boundaryStartTablesByModel=null, boundaryStopTablesByModel=null;
	private final int boundaryStartInside, boundaryStartOutside, boundaryStopInside, boundaryStopOutside;
	private final float boundaryMeanLen;
	private final int[] boundaryStartOffsets, boundaryStopOffsets;
	private final float boundaryMarginStart, boundaryMarginStop;
	int boundaryFeatureVersion=NcrnaBoundaryScorer.FEATURES_V1;
	NcrnaBoundaryFuzzConstants boundaryFuzzConstants=null;

	//Boundary-NN instrumentation (Citan/Brian, 2026-08-28) -- opt-in, off by default. A setter
	//rather than another constructor param: NcrnaScavenger already has 3 overloaded
	//constructors with 8-16 params each; threading one more opt-in field through all of them
	//would touch every call site (CallGenes/NcrnaFamily/GeneCaller) for a feature that's off in
	//every production run. Null-checked once per ACCEPTED locus only (see
	//trimToAlignmentExtent) -- not the per-candidate-window hot path, so this costs nothing on
	//the scan loop either way.
	private NcrnaBoundaryInstrumentSink instrumentSink=null;
	private IdentityHashMap<Orf,InstrumentationCapture> pendingInstrumentation=null;
	private BoundaryScoreSink boundaryScoreSink=null;
	private RrnaEndpointCallerFeatures rrnaEndpointFeatures=null;
	void setRrnaEndpointFeatures(RrnaEndpointCallerFeatures.Resources resources,RrnaEndpointCallerFeatures.Sink sink){
		if(resources==null || !Arrays.equals(modelNames,resources.names) || library.length!=resources.refs.length){throw new IllegalArgumentException("Endpoint feature library must match the caller model order");}
		for(int i=0;i<library.length;i++){if(!Arrays.equals(library[i],resources.refs[i])){throw new IllegalArgumentException("Endpoint consensus bases differ from caller model "+i);}}
		if(!family.equals(resources.family) || !reuseConsensusAlignment || trimAlignmentExtent || ((voteEnds || voteWindows) && voteTable==null)
				|| boundaryNetsByModel!=null || boundary5Net!=null || boundary3Net!=null){
			throw new IllegalArgumentException("Family-matched endpoint processing requires reused Quantum alignment, optional bound votes, and no legacy endpoint net");
		}
		rrnaEndpointFeatures=new RrnaEndpointCallerFeatures(resources,sink);
	}
	/** Test-visible caller seam; no work or validation occurs while observation is off. */
	void captureRrnaEndpointFeatures(Orf orf,byte[] bases,int model){
		if(rrnaEndpointFeatures==null){return;}
		final long start=System.nanoTime(),rawLength=orf.stop-orf.start+1L;
		try{rrnaEndpointFeatures.capture(orf,bases,model,contigGC(bases),minLen,maxLen);alignmentCount++;alignedBases+=rawLength;}
		catch(RuntimeException | AssertionError e){shared.KillSwitch.assertDie("Endpoint inference/observation failed on caller worker: "+e);throw e;}
		finally{boundaryNanos+=System.nanoTime()-start;}
	}

	//B4 seed-trigger/workload instrumentation (Citan/G11, 2026-08-29) -- opt-in, off by default,
	//same cost model as instrumentSink above (single null-check per call site, zero cost when
	//off). Separate field from instrumentSink: different question (workload/gate-metric capture
	//on EVERY scanned strand and EVERY scheduled window, vs accepted-locus-only boundary state),
	//different consumer (a B4/B5 family evaluator driver, not boundary-NN training).
	private NcrnaWorkloadInstrumentSink workloadSink=null;
	private NcrnaModelAttemptInstrumentSink modelAttemptSink=null;
	private NcrnaStageDiagSink diagSink=null;
	private String diagFamily=null;
	private boolean mappingOracleExhaustive=false;

	void setStageDiagSink(NcrnaStageDiagSink sink,String familyName){diagSink=sink;diagFamily=familyName;}
	/** Independent generated-window observer; null leaves the scientific path unchanged. */
	void setVoteDiagSink(NcrnaVoteDiagSink sink, String familyName){voteDiagSink=sink;voteDiagFamily=familyName;}
	private NcrnaVoteDiagSink voteDiagSink;
	private String voteDiagFamily;
	private NcrnaVoteDiagnostics voteTrace;
	/** Calibration-only exhaustive verifier arm: bypass mapping filters and early exits,
	 * retaining seed/window/identity/HBM gates. Model traversal is library order, not rank. */
	void setMappingOracleExhaustive(boolean value){mappingOracleExhaustive=value;}

	/** Arms (or disarms, via null) B4 workload instrumentation capture. Package-visible: only an
	 * evaluation driver in this package should call this, mirroring setInstrumentSink's scoping
	 * rationale (never wired to a production CallGenes flag without an explicit opt-in gate). No
	 * arm-time validation needed here (unlike setInstrumentSink's modelNames-length assert) --
	 * seedHits/scheduledWindow never index into modelNames, so there is no equivalent silent-skip
	 * hazard to guard against at arm time. */
	void setWorkloadSink(NcrnaWorkloadInstrumentSink sink){workloadSink=sink;}
	void setModelAttemptSink(NcrnaModelAttemptInstrumentSink sink){modelAttemptSink=sink;}
	void setBoundaryScoreSink(BoundaryScoreSink sink){boundaryScoreSink=sink;}

	/** Development-only endpoint-sweep evidence, two rows per refined locus. */
	interface BoundaryScoreSink{
		void scored(String contig, int strand, int model, int anchorStart, int anchorStop,
			int startOffset, int stopOffset, int startSitesScored, int stopSitesScored,
			float startChosenScore, float stopChosenScore,
			int startRadius, int stopRadius, boolean startVotedAnchor, boolean stopVotedAnchor);
	}

	/** Arms (or disarms, via null) boundary-NN instrumentation capture. Package-visible: only
	 * an instrumentation driver in this package should call this, never production CallGenes
	 * flag plumbing without an explicit opt-in flag gating it. */
	/** Fail-loud arming, per Citan (2026-08-28): all 3 accepted branches in alignWindow call
	 * trimToAlignmentExtent only inside `if(annotate && modelNames!=null && bestModel&lt;
	 * modelNames.length)` -- if modelNames were null or shorter than library, an accepted locus
	 * would silently skip trim AND capture entirely, undercounting the accepted denominator
	 * again, exactly the class of bug just fixed for the trim-guard case. Validating at ARM time
	 * (not silently at run time) guarantees that once instrumentation is successfully armed,
	 * every accepted bestModel is GUARANTEED to be a valid modelNames index, so this specific
	 * gate can never again be the reason a capture is skipped. A plain assert (not
	 * KillSwitch.assertDie): this runs on the single calling thread that builds and arms an
	 * instrumentation driver, never a producer/consumer worker thread, so an AssertionError here
	 * cannot leave anything else silently hung -- the textbook case for a plain assert per the
	 * assertions skill. Sink-off (sink==null) never runs this check -- zero change to production
	 * arming behavior (there is none) or cost. */
	void setInstrumentSink(NcrnaBoundaryInstrumentSink sink){
		if(sink!=null){
			assert(modelNames!=null && modelNames.length==library.length) : "Cannot arm boundary-NN "
				+"instrumentation: modelNames is "+(modelNames==null ? "null" : "length "+modelNames.length)
				+" but library has "+(library==null ? "null" : ""+library.length)+" models -- they must be "
				+"non-null and equal length, or alignWindow's annotate/modelNames guard would silently skip "
				+"trim+capture for some accepted loci (NcrnaScavenger.java, the 3 trimToAlignmentExtent call "
				+"sites), undercounting the accepted denominator just like the trim-guard bug this replaces.";
		}
		instrumentSink=sink;
		pendingInstrumentation=(sink==null ? null : new IdentityHashMap<Orf,InstrumentationCapture>());
	}

	private static final class InstrumentationCapture{
		InstrumentationCapture(String c,int s,int m,int ws,int we,int ps,int pe,byte[] w,int wo,boolean t,boolean n,int source,float id,int al){contig=c;strand=s;model=m;wStart=ws;wStop=we;postStart=ps;postStop=pe;window=w;windowOffset=wo;trimSucceeded=t;nnInvoked=n;windowSource=source;identity=id;alignedLength=al;}
		final String contig;final int strand,model,wStart,wStop,postStart,postStop,windowOffset,windowSource,alignedLength;final byte[] window;final boolean trimSucceeded,nnInvoked;final float identity;
	}

	/** Padding beyond the post-trim [start,stop] captured into the instrumentation window copy
	 * -- must cover both the family-configured boundary candidate arrays and the enrichment profile's own local
	 * +-2 radius plus the widest currently-staged k-mer window (k=11, srp_small) -- 4(sweep)+
	 * 2(local radius)+11(k)=17 is the true minimum reach past the boundary; 30 leaves real
	 * margin without meaningfully growing the copy. */
	static final int INSTRUMENT_CAPTURE_PAD=30;
}
