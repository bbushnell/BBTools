package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import aligner.SingleStateAlignerFlat2;
import aligner.SingleStateAlignerFlat3;
import aligner.SingleStateAlignerFlatFloat;
import dna.AminoAcid;
import dna.GeneticCode;
import idaligner.IDAligner;
import map.LongHashSet;
import ml.CellNet;
import ml.CellNetParser;
import shared.KillSwitch;
import shared.Tools;
import simd.Vector;
import stream.Read;
import structures.FloatList;
import structures.IntList;


/**
 * This class calls genes within a single thread.
 * @author Brian
 * @date Sep 24, 2018
 *
 */
public class GeneCaller extends ProkObject {
	
	/*--------------------------------------------------------------*/
	/*----------------             Init             ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructs a GeneCaller with specified thresholds and gene model.
	 *
	 * @param minLen_ Minimum gene length in bases
	 * @param maxOverlapSameStrand_ Maximum overlap for same-strand genes
	 * @param maxOverlapOppositeStrand_ Maximum overlap for opposite-strand genes
	 * @param minStartScore_ Minimum start codon score threshold
	 * @param minStopScore_ Minimum stop codon score threshold
	 * @param minInnerScore_ Minimum average inner k-mer score
	 * @param minOrfScore_ Minimum overall ORF score
	 * @param minAvgScore_ Minimum score per base
	 * @param pgm_ Gene model containing probability matrices
	 */
	public GeneCaller(int minLen_, int maxOverlapSameStrand_, int maxOverlapOppositeStrand_,
			float minStartScore_, float minStopScore_, float minInnerScore_,
			float minOrfScore_, float minAvgScore_, GeneModel pgm_){
		this(minLen_, maxOverlapSameStrand_, maxOverlapOppositeStrand_, minStartScore_, minStopScore_,
			minInnerScore_, minOrfScore_, minAvgScore_, pgm_, null);
	}

	/** The run's code is independent of replaceable PGM statistics; null preserves legacy starts/stops. */
	public GeneCaller(int minLen_, int maxOverlapSameStrand_, int maxOverlapOppositeStrand_,
			float minStartScore_, float minStopScore_, float minInnerScore_,
			float minOrfScore_, float minAvgScore_, GeneModel pgm_, GeneticCode geneticCode_){
		geneticCode=geneticCode_;
		minLen=minLen_;
		maxOverlapSameStrand=maxOverlapSameStrand_;
		maxOverlapOppositeStrand=maxOverlapOppositeStrand_;
		pgm=pgm_;

		minStartScore=minStartScore_;
		minStopScore=minStopScore_;
		minInnerScore=minInnerScore_;
		minOrfScore=minOrfScore_;
		minAvgScore=minAvgScore_;
		assert(pgm!=null);
		orfNet=(orfNetTemplate!=null ? orfNetTemplate.copy(false) : null);
		orfNetLow=(orfNetTemplateLow!=null ? orfNetTemplateLow.copy(false) : null);
		orfNetMid=(orfNetTemplateMid!=null ? orfNetTemplateMid.copy(false) : null);
		orfNetHigh=(orfNetTemplateHigh!=null ? orfNetTemplateHigh.copy(false) : null);
		if(orfNetTemplateGC!=null){
			orfNetsGC=new CellNet[orfNetTemplateGC.length];
			for(int i=0; i<orfNetTemplateGC.length; i++){
				orfNetsGC[i]=(orfNetTemplateGC[i]!=null ? orfNetTemplateGC[i].copy(false) : null);
			}
		}else{orfNetsGC=null;}
		if(stopNetTemplateGC!=null){
			stopNetsGC=new CellNet[stopNetTemplateGC.length];
			for(int i=0; i<stopNetTemplateGC.length; i++){
				stopNetsGC[i]=(stopNetTemplateGC[i]!=null ? stopNetTemplateGC[i].copy(false) : null);
			}
		}else{stopNetsGC=null;}
		if(metaNetTemplateGC!=null){
			metaNetsGC=new CellNet[metaNetTemplateGC.length];
			for(int i=0; i<metaNetTemplateGC.length; i++){
				metaNetsGC[i]=(metaNetTemplateGC[i]!=null ? metaNetTemplateGC[i].copy(false) : null);
			}
		}else{metaNetsGC=null;}
		metaNetInstance=(metaNetSingle!=null ? metaNetSingle.copy(false) : null);
		gc=pgm.gc();
	}

	public void setPGM(GeneModel pgm_){
		pgm=pgm_;
		gc=pgm.gc();
		trnaCaller=null;
	}

	/** Explicit full-family opt-in; each caller owns independent alignment scratch. */
	public void setProteinSupport(final ProteinSupport reference){
		rejectProteinSupportConflict(reference!=null, stopTailSupport!=null || tailCascade!=null);
		proteinSupport=reference==null ? null : reference.newWorker();
	}
	/** Explicit stop-tail opt-in; mutually exclusive with full-family evidence. */
	public void setStopTailSupport(final StopTailSupport reference){
		rejectProteinSupportConflict(proteinSupport!=null || tailCascade!=null, reference!=null);
		stopTailSupport=reference==null ? null : reference.newWorker();
	}
	private static void rejectProteinSupportConflict(final boolean oldSupport, final boolean tailSupport){
		if(oldSupport && tailSupport){
			throw new IllegalStateException("Use either full-family or stop-tail support, not both");
		}
	}
	ProteinSupport.Worker proteinSupport;
	/** Run-owned optional experiment; never stored on the shared family roster. */
	public void setCobalaminDdl(final NcrnaDdlProfileResource profiles, final int cap){
		if(ncrnaScavengers!=null || cobalaminDdlProfiles!=null){
			throw new IllegalArgumentException("Install one RF00174 profile on a fresh caller before lazy worker construction");
		}
		NcrnaDdlShortlist.validateCap(profiles, cap);
		cobalaminDdlProfiles=profiles;
		cobalaminDdlCap=cap;
	}
	private NcrnaDdlProfileResource cobalaminDdlProfiles;
	/** Published with the complete family roster; all mutable search state stays on workers. */
	public static CobalaminCallerResources cobalaminCallerResources;
	private CobalaminCallerResources cobalaminServing;
	private int cobalaminDdlCap;
	/** Explicit three-member experiment. Profiles are either absent (native), one
	 * joint resource, or three independent member resources in separate mode. */
	public void setCobalaminClan(final boolean joint, final NcrnaDdlProfileResource[] profiles){
		setCobalaminClan(joint, profiles, true);
	}
	/** Native ranking can be omitted only for the explicit joint DDL experiment. */
	public void setCobalaminClan(final boolean joint, final NcrnaDdlProfileResource[] profiles, final boolean nativeRanking){
		if(ncrnaScavengers!=null || cobalaminDdlProfiles!=null || cobalaminClanConfigured){
			throw new IllegalStateException("Configure clan execution once on a fresh caller before any scavenger exists");
		}
		if(profiles!=null && profiles.length!=(joint ? 1 : 3)){
			throw new IllegalArgumentException("Joint execution needs one combined profile; separate execution needs three member profiles");
		}
		if(!nativeRanking && (!joint || profiles==null)){
			throw new IllegalArgumentException("An index-free clan requires joint execution and its explicit DDL profiles");
		}
		cobalaminClanConfigured=true;
		cobalaminJoint=joint;
		cobalaminNativeRanking=nativeRanking;
		cobalaminMemberProfiles=profiles==null ? null : profiles.clone();
	}
	private boolean cobalaminClanConfigured, cobalaminJoint;
	private boolean cobalaminNativeRanking=true;
	private NcrnaDdlProfileResource[] cobalaminMemberProfiles;
	private int[] cobalaminMemberSlots;
	private boolean inactiveClanMember(final int index){
		return cobalaminJoint && cobalaminMemberSlots!=null
			&& (index==cobalaminMemberSlots[1] || index==cobalaminMemberSlots[2]);
	}
	StopTailSupport.Worker stopTailSupport;
	/** Shared bounded search pool; this caller owns only its contig evidence and counters. */
	public void setTailCascadeSupport(final StopTailCascadeSupport reference, final float bonus, final boolean protect){
		setTailCascadeSupport(reference, bonus, protect, null);
	}
	/** Optional manifest identity accompanies emitted evidence; direct library fixtures may omit it. */
	public void setTailCascadeSupport(final StopTailCascadeSupport reference, final float bonus,
			final boolean protect, final String pin){
		setTailCascadeSupport(reference, bonus, protect, pin, false);
	}
	/** Opt-in graded evidence reuses cached winning-family margins; binary behavior remains the default. */
	public void setTailCascadeSupport(final StopTailCascadeSupport reference, final float bonus,
			final boolean protect, final String pin, final boolean marginMode){
		configureTailCascadeSupport(reference, bonus, protect, pin, marginMode, false);
	}
	/** Bind explicitly structural resources, preserving an inherited identity or literal unmeasured
	 * in each emitted ORF. This entry point does not establish payload byte identity.
	 */
	public void setTailCascadeSupportStructural(final StopTailCascadeSupport reference, final float bonus,
			final boolean protect, final String pin, final boolean marginMode){
		configureTailCascadeSupport(reference, bonus, protect, pin, marginMode, true);
	}
	/** Explicit information mode shares the same ordinary/replay score insertion. */
	public void setTailCascadeInformationSupport(final StopTailCascadeSupport reference,final float bonus,
			final boolean protect,final String pin,final boolean structural){
		if(reference!=null && bonus!=0){reference.requireMatchInformation();}
		configureTailCascadeSupport(reference,bonus,protect,pin,false,structural);
		tailCascadeInformationMode=true;
	}
	/** Formula constants belong to the immutable formula, not an additional tailbonus multiplier. */
	public void setTailCascadeFormulaSupport(final StopTailCascadeSupport reference,final prot.FamilyBonusFormula formula,
			final boolean protect,final String pin,final boolean structural){
		setTailCascadeFormulaSupport(reference,formula,protect,pin,structural,false);
	}
	/** Preserve the explicit raw-divisor override through caller-worker configuration. */
	public void setTailCascadeFormulaSupport(final StopTailCascadeSupport reference,final prot.FamilyBonusFormula formula,
			final boolean protect,final String pin,final boolean structural,final boolean rawScoreI0Override){
		if(reference==null || formula==null || tailCascadeOppositePenalty!=0){
			throw new IllegalArgumentException("Formula scoring requires an explicit classifier/formula "+
				"and no opposite-strand penalty");
		}
		final prot.FamilyBonusFormula bound=reference.bindFormulaEvidence(formula,rawScoreI0Override);
		configureTailCascadeSupport(reference,0,protect,pin,false,structural);
		tailCascadeFormula=bound;
	}
	/** Share scoring/protection setup while keeping the ordinary identity contract strict. */
	private void configureTailCascadeSupport(final StopTailCascadeSupport reference, final float bonus,
			final boolean protect, final String pin, final boolean marginMode, final boolean structural){
		rejectProteinSupportConflict(proteinSupport!=null || stopTailSupport!=null, reference!=null);
		if(!Float.isFinite(bonus) || bonus<0){
			throw new IllegalArgumentException("Tail bonus must be finite and nonnegative");
		}
		if(structural){
			if(pin==null || (!pin.matches("[0-9a-f]{20}") && !pin.equals("unmeasured"))){
				throw new IllegalArgumentException("Structural tail resource identity must be inherited sha80 or unmeasured");
			}
		}else if(pin!=null && !pin.matches("[0-9a-f]{20}")){
			throw new IllegalArgumentException("Tail resource identity must be sha80");
		}
		if(reference!=null){cascadeCode(geneticCode); if(marginMode){reference.requireAlignmentMargins();}}
		tailCascade=reference; tailCascadeBonus=bonus; tailCascadeProtect=protect; tailCascadePin=pin;
		tailCascadeMarginMode=marginMode;
		tailCascadeInformationMode=false;
		tailCascadeFormula=null;
	}
	/** Observe legacy stops without replacing the caller's native start mask or output code. */
	public static GeneticCode cascadeCode(final GeneticCode explicit){
		if(explicit!=null){return explicit;}
		final GeneticCode code=GeneticCode.forTable(11);
		for(int i=0; i<64; i++){
			if(GeneModel.isStopCodon(i)!=code.isStop(i)){
				throw new IllegalArgumentException("Tail classification requires an explicit genetic code "+
					"for a nonstandard legacy stop mask");
			}
		}
		return code;
	}
	private StopTailCascadeSupport tailCascade;
	private StopTailCascadeSupport.ContigEvidence tailCascadeEvidence;
	private StopTailCascadeSupport.Replay preparedTailReplay;
	/** Active only inside one targeted diagnostic replay; never part of ordinary caller policy. */
	private long targetedTailStop=-1;
	private float targetedTailBonus;
	private boolean[] tailCascadeParents;
	private float tailCascadeBonus;
	private boolean tailCascadeProtect;
	private boolean tailCascadeMarginMode;
	private boolean tailCascadeInformationMode;
	private prot.FamilyBonusFormula tailCascadeFormula;
	private String tailCascadePin;
	public final TailCascadeStats tailCascadeStats=new TailCascadeStats();
	/** Null by default; receives primitive snapshots without access to mutable candidates. */
	private TailCascadeObserver tailCascadeObserver;
	private GeneCallerPathTrace.Observer pathObserver;
	private NcrnaPathScoreBasis pathScoreBasis;
	private GeneCallerPathTrace pathTrace;
	private boolean pathActive;
	/** Private experiment only; absent in ordinary callers. */
	private NcrnaCdsStartSelection rnaCdsStartSelection;
	private float tailCascadeOppositePenalty;
	private TailCascadeOppositePenalty.Observer oppositePenaltyObserver;
	/** Optional post-pruning preference. Native and zero-bonus callers remain disabled. */
	public void setTailCascadeOppositePenalty(final float penalty){
		if(pathActive || tailCascadeEvidence!=null || preparedTailReplay!=null){
			throw new IllegalStateException("Change opposite-strand policy only between completed contigs");
		}
		if(!Float.isFinite(penalty) || penalty<0){throw new IllegalArgumentException("Invalid opposite-strand penalty");}
		if(tailCascadeFormula!=null && penalty!=0){
			throw new IllegalArgumentException("Formula comparison cannot combine with opposite-strand penalties");
		}
		tailCascadeOppositePenalty=penalty;
	}
	/** Observation is local to this caller and receives only primitive values. */
	public void setOppositePenaltyObserver(final TailCascadeOppositePenalty.Observer observer){
		if(pathActive || tailCascadeEvidence!=null || preparedTailReplay!=null){
			throw new IllegalStateException("Change opposite-strand observer only between completed contigs");
		}
		oppositePenaltyObserver=observer;
	}
	/** Install between contigs; callbacks receive snapshots and may not reenter this caller. */
	public void setPathObserver(final GeneCallerPathTrace.Observer observer){
		if(pathActive || tailCascadeEvidence!=null || preparedTailReplay!=null){
			throw new IllegalStateException("Change the path observer only between completed contigs");
		}
		pathObserver=observer;
	}
	/** Install between contigs; the table is private to this native input worker. */
	public void setPathScoreBasis(final NcrnaPathScoreBasis basis){
		if(pathActive){throw new IllegalStateException("Change score provenance only between completed contigs");}
		pathScoreBasis=basis;
		installPathScoreBasis();
	}
	private void installPathScoreBasis(){
		if(ncrnaScavengers==null){return;}
		for(int i=0; i<ncrnaScavengers.size(); i++){
			final NcrnaScavenger scorer=ncrnaScavengers.get(i);
			scorer.setPathScoreObserver(pathScoreBasis==null || inactiveClanMember(i) ? null : pathScoreBasis.bind(scorer));
		}
	}
	/** Read-only policy capture is allowed at the path observation seam. */
	public GeneCallerPathPolicy pathPolicy(){
		verifyPathReplayPolicy();
		return new GeneCallerPathPolicy(this, geneticCode);
	}
	void verifyPathReplayPolicy(){
		if(proteinSupport!=null || stopTailSupport!=null || tailCascade!=null || ssuCompetition!=null
			|| rnaCdsStartSelection!=null || iterationCandidates!=null || cmVerifier!=null
			|| CallGenes.NCRNA_REFRESH_CLAIMED_WINDOWS){
			throw new IllegalStateException("Path replay excludes protein/SSU/start-reselection/candidate-reuse/CM/refresh policies");
		}
		if(ncrnaScavengers!=null){
			for(final NcrnaScavenger scorer:ncrnaScavengers){
				if(scorer.refreshClaimedWindows){throw new IllegalStateException("Path replay requires effective frozen scavenger windows");}
			}
		}
	}
	static float pathMinorStartMultiplier(){return MINOR_START_MULTIPLIER;}
	GeneticCode pathGeneticCode(){return geneticCode;}
	/** Option-off, clan-only experiment. The selector itself has no family-specific policy. */
	public void setRnaAwareCdsStarts(final boolean enabled){
		if(pathActive){throw new IllegalStateException("Change RNA-aware start selection only between contigs");}
		rnaCdsStartSelection=enabled ? new NcrnaCdsStartSelection() : null;
	}
	/** The ordinary pre-merge seam; exposed for finite option-off/family-scope fixtures. */
	public int selectRnaAwareCdsStarts(final ArrayList<Orf>[] frames, final ArrayList<Orf>[] rnas){
		if(rnaCdsStartSelection==null || rnas==null){return 0;}
		int removed=0;
		for(final ArrayList<Orf> list:rnas){
			for(final Orf rna:list){
				if(NcrnaFamilyRegistry.isCobalamin(rna.ncrnaFamily)){
					removed+=rnaCdsStartSelection.apply(frames, rna, maxOverlapSameStrand, maxOverlapOppositeStrand);
				}
			}
		}
		return removed;
	}
	public void setTailCascadeObserver(final TailCascadeObserver observer){
		if(tailCascadeEvidence!=null || preparedTailReplay!=null){
			throw new IllegalStateException("Change the tail observer only between completed contigs");
		}
		tailCascadeObserver=observer;
	}
	/** Observation only, on this caller's thread. Candidate coordinates are zero-based
	 * inclusive genomic coordinates. ordinal=-1 means no classified real stop; realStop
	 * distinguishes an unscanned real stop from a truncated boundary. Each admitted
	 * candidate is reported before pruning, then again if retained. Neither callback
	 * may change scoring, admission, pruning or the selected path. */
	public interface TailCascadeObserver{
		void begin(StopTailCascadeSupport.ContigEvidence evidence);
		void parent(int ordinal);
		void candidate(int start, int stop, int strand, boolean realStop, int ordinal,
			int candidateAaLength, float nativeScore, float bonus, float totalScore, boolean retained);
		void end();
	}
	/** Optional observer; never changes candidate admission, scores or path selection. */
	ProteinSupportFate proteinFateSink;
	private ArrayList<ProteinSupportFate.Entry> proteinFates;

	/** Apply the separate post-pruning preference exactly once to a node contribution. */
	public static float pathRootScore(final Orf orf){
		return subtractOppositePenalty(orf, orf.orfScore);
	}
	private static float subtractOppositePenalty(final Orf orf, final float score){
		return orf.tailCascadeOppositePenalty==0 ? score : score-orf.tailCascadeOppositePenalty;
	}

	/** Add conserved-family evidence once after overlap/minor-start terms, then the separate penalty. */
	public static float proteinAwareScore(final Orf orf, final int overlap){
		assert(overlap>=0) : "DP overlap is a nonnegative number of bases";
		if(overlap==0 || orf.type==ProkObject.RNA){return pathRootScore(orf);}
		float score=orf.calcOrfScore(overlap);
		if(MINOR_START_MULTIPLIER!=1f && orf.type==CDS){
			score=minorStartScore(score, orf.startCodon, orf.startTruncated, MINOR_START_MULTIPLIER);
		}
		// Formula evidence may be a signed contribution; overlap DP must preserve it.
		return subtractOppositePenalty(orf, orf.proteinBonus!=0 ? score+orf.proteinBonus : score);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Calls genes on a collection of reads.
	 * Processes each read individually and combines results.
	 * @param reads Collection of reads to process
	 * @return List of all identified ORFs across all reads
	 */
	public ArrayList<Orf> callGenes(ArrayList<Read> reads){
		ArrayList<Orf> orfs=new ArrayList<Orf>();
		for(Read r : reads) {
			ArrayList<Orf> list=callGenes(r);
			if(list!=null) {orfs.addAll(list);}
		}
		return orfs;
	}
	
	/**
	 * Calls genes on a single read using default settings.
	 * @param r Read to process
	 * @return List of identified ORFs
	 */
	public ArrayList<Orf> callGenes(Read r){
		return callGenes(r, pgm, true);
	}
	
	/**
	 * Calls genes on a single read with configurable ORF breaking.
	 * @param r Read to process
	 * @param breakOrfs Whether to break long ORFs into multiple starts
	 * @return List of identified ORFs
	 */
	public ArrayList<Orf> callGenes(Read r, boolean breakOrfs){
		return callGenes(r, pgm, breakOrfs);
	}

	/** Returns ALL scored ORF candidates without DP path selection.
	 * Used for generating training vectors for the NN scorer. */
	public ArrayList<Orf> getAllScoredOrfs(Read r, GeneModel pgm_){
		pgm=pgm_;
		final byte[] bases=r.bases;
		final String name=r.id;
		ArrayList<Orf>[] frameLists=makeOrfs(name, bases, minLen, geneticCode);
		ArrayList<Orf>[] brokenLists=breakOrfs(frameLists, bases);
		ArrayList<Orf> all=new ArrayList<Orf>();
		for(ArrayList<Orf> list : brokenLists){
			if(list!=null){all.addAll(list);}
		}
		Collections.sort(all);
		return all;
	}
	
	/**
	 * Main gene calling method with full parameter control.
	 * Identifies ORFs on both strands, calls RNA features, and finds optimal path.
	 *
	 * @param r Read to process
	 * @param pgm_ Gene model to use for scoring
	 * @param breakOrfs Whether to generate multiple starts per stop
	 * @return List of non-overlapping ORFs in optimal arrangement
	 */
	public ArrayList<Orf> callGenes(Read r, GeneModel pgm_, boolean breakOrfs){
		if(pathActive){throw new IllegalStateException("A path observer must not reenter its caller");}
		if(proteinFateSink==null && tailCascade==null && pathObserver==null && boundaryTrace==null){return callGenesInner(r, pgm_, breakOrfs);}
		if(boundaryTrace!=null){
			if(iterationCandidates!=null || iterationSeeds!=null){
				throw new IllegalArgumentException("Boundary tracing does not support iteration-candidate or seed reuse");
			}
		}
		try{
			if(boundaryTrace!=null){boundaryTrace.begin(r.id, r.length(), r.numericID);}
			return callGenesInner(r, pgm_, breakOrfs);
		}
		catch(RuntimeException | Error failure){
			if(boundaryTrace!=null){boundaryTrace.clear();}
			KillSwitch.kill((boundaryTrace==null ? "Protein fate diagnostics" : "Boundary trace diagnostics")+" failed on "+r.id+": "+failure);
			throw failure;
		}finally{if(boundaryTrace!=null){boundaryTrace.clear();}}
	}
	/** Only Replay owns the private input copy; cached evidence cannot be supplied
	 * separately from its contig. Native controls share the identical input bytes. */
	ArrayList<Orf> callGenesWithPreparedTail(final Read read, final GeneModel model,
			final StopTailCascadeSupport.Replay replay){
		if(preparedTailReplay!=null || proteinSupport!=null || stopTailSupport!=null){
			throw new IllegalStateException("Tail replay requires a non-reentrant native or cascade caller");
		}
		if(cascadeCode(geneticCode)!=replay.evidence.code || tailCascade!=null && tailCascade!=replay.owner){
			throw new IllegalArgumentException("Replay and caller must use the same genetic code and classifier");
		}
		if(tailCascade==null){return callGenes(read, model, true);}
		preparedTailReplay=replay;
		try{return callGenes(read, model, true);}
		finally{preparedTailReplay=null;}
	}
	/** Raise one accepted stop's binary bonus while retaining normal admission, protection,
	 * pruning and DP. Operon/lookback changes can follow a different path; they are outcomes,
	 * not independently changed parameters. Replay owns the private contig and evidence. */
	ArrayList<Orf> callGenesWithTargetedTailBonus(final Read read, final GeneModel model,
			final StopTailCascadeSupport.Replay replay, final int ordinal, final float bonus){
		if(preparedTailReplay!=null || pathActive || targetedTailStop!=-1){
			throw new IllegalStateException("Targeted bonus replay requires an idle caller");
		}
		if(tailCascade!=replay.owner || tailCascadeMarginMode || tailCascadeInformationMode ||
			tailCascadeFormula!=null || ordinal<0 || ordinal>=replay.evidence.size() ||
			!replay.evidence.assignment(ordinal).assigned()){
			throw new IllegalArgumentException("Targeted replay requires the same binary classifier and an accepted stop");
		}
		if(!Float.isFinite(bonus) || bonus<tailCascadeBonus){
			throw new IllegalArgumentException("Targeted stop bonus must be finite and at least the ordinary binary bonus");
		}
		targetedTailStop=replay.evidence.key(ordinal);targetedTailBonus=bonus;
		try{return callGenesWithPreparedTail(read,model,replay);}
		finally{targetedTailStop=-1;targetedTailBonus=0;}
	}
	/** Same caller path in both modes; diagnostic failures must not strand output consumers. */
	private ArrayList<Orf> callGenesInner(Read r, GeneModel pgm_, boolean breakOrfs){
		if((proteinSupport!=null || stopTailSupport!=null || tailCascade!=null) && !breakOrfs){
			throw new IllegalArgumentException("Protein support requires candidate start enumeration");
		}
		if(proteinFateSink!=null && proteinSupport==null){
			throw new IllegalStateException("Protein fate diagnostics require the actual protein-support search");
		}
		proteinFates=proteinFateSink==null ? null : new ArrayList<ProteinSupportFate.Entry>();
		if(ncrnaCpuProfile==null && NcrnaCpuProfile.enabled()){ncrnaCpuProfile=NcrnaCpuProfile.create(ncrnaFamilies.size());}
		pgm=pgm_;
		assert(pgm!=null);

		final String name=r.id;
		final byte[] bases=r.bases;
		assert(iterationCandidates==null || (iterationBases==bases && breakOrfs))
			: KillSwitch.assertDie("Candidate reuse is confined to one unchanged contig with ordinary CDS start enumeration");
		contigGC=calcContigGC(bases);
		activeNet=selectNet(contigGC);
		activeStopNet=selectStopNet(contigGC);
		activeMetaNet=selectMetaNet(contigGC);
		if(tailCascade!=null){
			if(preparedTailReplay==null){
				final StopTailCandidateScan.Scan scan=tailCascade.scan(name, bases, cascadeCode(geneticCode));
				tailCascadeEvidence=tailCascade.classify(scan);
				tailCascadeStats.under100+=scan.under100;
				tailCascadeStats.from100To199+=scan.from100To199;
			}else{
				tailCascadeEvidence=preparedTailReplay.evidence;
				tailCascadeStats.under100+=preparedTailReplay.under100;
				tailCascadeStats.from100To199+=preparedTailReplay.from100To199;
				tailCascadeStats.reused+=tailCascadeEvidence.size();
			}
			tailCascadeParents=new boolean[tailCascadeEvidence.size()];
			if(tailCascadeObserver!=null){tailCascadeObserver.begin(tailCascadeEvidence);}
			tailCascadeStats.searched+=tailCascadeEvidence.size();
			tailCascadeStats.accepted+=tailCascadeEvidence.accepted;
			for(int i=0; i<tailCascadeEvidence.size(); i++){
				final prot.FamilyTailCascade.Assignment assignment=tailCascadeEvidence.assignment(i);
				if(assignment.assigned()){tailCascadeStats.acceptedStages[TailCascadeStats.stage(assignment)]++;}
			}
		}

		ArrayList<Orf>[] brokenLists=iterationCandidates==null ? null : iterationCandidates.restoreCds(pgm_);
		if(brokenLists==null){
			// All scored pre-path starts, including alternatives sharing a stop.
			final ArrayList<Orf>[] frameLists=makeOrfs(name, bases, minLen, geneticCode);
			brokenLists=(breakOrfs ? breakOrfs(frameLists, bases) : frameLists);
			if(iterationCandidates!=null){
				if(!breakOrfs){throw new IllegalArgumentException("Candidate reuse requires ordinary CDS start enumeration");}
				iterationCandidates.captureCds(brokenLists, pgm_);
			}
		}
		if(!breakOrfs) {
			for(int i=0; i<brokenLists.length; i++) {
				if(brokenLists[i]==null) {brokenLists[i]=new ArrayList<Orf>(1);}
			}
		}
		
		ArrayList<Orf>[] rnaLists=null;
		final int rlen=r.length();
		if(calltRNA || (call16S && rlen>800) || (call23S && rlen>1500) || call5S ||
			(call18S && rlen>1000) || !ncrnaFamilies.isEmpty()){
			rnaLists=makeRnas(name, bases, brokenLists);
			if(ncrnaCpuProfile!=null){ncrnaCpuProfile.begin(ncrnaFamilies.size(), NcrnaCpuProfile.SSU_APPLY);}
			try{applySsuCompetition(rnaLists, bases);}
			finally{if(ncrnaCpuProfile!=null){ncrnaCpuProfile.end();}}
			if(NcrnaScavenger.DEBUG){
				for(ArrayList<Orf> list : rnaLists){for(Orf o : list){
				System.err.println("DEBUG rnaList candidate: "+o.start+"-"+o.stop+" strand="+o.strand+
					" type="+o.type+" score="+o.orfScore);
				}}
				for(Orf o : brokenLists[0]){
					if(o.stop>=100 && o.start<=400){
						System.err.println("DEBUG brokenLists[0] near locus: "+o.start+"-"+o.stop+
							" score="+o.orfScore+" type="+o.type);
					}
				}
			}
			// Private experiment: reuse retained, scored starts before shared path search.
			// No score, codon, stop or RNA report is changed. The accepted RNA may still
			// lose the path, so this option stays off pending whole-genome CDS evaluation.
			selectRnaAwareCdsStarts(brokenLists, rnaLists);
			if(rnaLists[0]!=null && !rnaLists[0].isEmpty()) {
				brokenLists[0].addAll(rnaLists[0]);
				Collections.sort(brokenLists[0]);
			}
			if(rnaLists[1]!=null && !rnaLists[1].isEmpty()) {
				brokenLists[3].addAll(rnaLists[1]);
				Collections.sort(brokenLists[3]);
			}
			if(NcrnaScavenger.DEBUG){
				System.err.println("DEBUG brokenLists[0] AFTER merge, near locus:");
				for(Orf o : brokenLists[0]){
					if(o.stop>=100 && o.start<=400){
						System.err.println("DEBUG   "+o.start+"-"+o.stop+" score="+o.orfScore+
							" type="+o.type);
					}
				}
			}
		}
		
		boolean printAllOrfs=false;
		boolean printRnas=false;
		if(printAllOrfs){
			ArrayList<Orf> temp=new ArrayList<Orf>();
			for(ArrayList<Orf> broken : brokenLists){
				temp.addAll(broken);
			}
			Collections.sort(temp);
//			System.err.print(temp.size()+"; ");
			return temp;
		}
		
		if(printRnas && rnaLists!=null){
			ArrayList<Orf> temp=new ArrayList<Orf>();
			for(ArrayList<Orf> list : rnaLists){
				temp.addAll(list);
			}
			Collections.sort(temp);
			return temp;
		}
		
		stCds2.add(brokenLists);
		
		//Find the optimal path through Orfs
		if(ncrnaDiagSink!=null){
			for(ArrayList<Orf> list : brokenLists){recordNcrnaPath(list, name, rlen, false);}
		}
		final ArrayList<Orf> path;
		if(ncrnaCpuProfile!=null){ncrnaCpuProfile.begin(ncrnaFamilies.size(), NcrnaCpuProfile.PATH);}
		try{path=findPath(brokenLists, bases);}
		finally{if(ncrnaCpuProfile!=null){ncrnaCpuProfile.end();}}
		if(tailCascade!=null){finishTailCascade(path);}
		if(proteinFateSink!=null){
			proteinFateSink.write(proteinFates, path, maxOverlapSameStrand, maxOverlapOppositeStrand);
			proteinFates=null;
		}
		if(ssuCompetition!=null){
			if(ncrnaCpuProfile!=null){ncrnaCpuProfile.begin(ncrnaFamilies.size(), NcrnaCpuProfile.SSU_FINISH);}
			try{ssuCompetition.finish(path);}
			catch(RuntimeException | AssertionError e){KillSwitch.assertDie("SSU final-path verification failed: "+e);throw e;}
			finally{if(ncrnaCpuProfile!=null){ncrnaCpuProfile.end();}}
		}
		if(ncrnaDiagSink!=null){recordNcrnaPath(path, name, rlen, true);}
//		geneStartsOut+=path.size();

		if(callCDS){stCdsPass.add(path);}
		if(calltRNA){sttRNA.add(path);}
		if(call16S){st16s.add(path);}
		if(call23S){st23s.add(path);}
		if(call5S){st5s.add(path);}
		if(call18S){st18s.add(path);}
		
		// Selection and its diagnostics use raw discovery bounds. Downstream GFF,
		// sequence extraction and length reports consume the refined output bounds.
		if(boundaryTrace!=null){boundaryTrace.beforeReport(path, maxOverlapSameStrand, maxOverlapOppositeStrand);}
		NcrnaReportedOrf.applyTo(path);
		if(boundaryTrace!=null){boundaryTrace.finish(path);}
		return path;
	}

	/** Compare only admitted RNA candidates, in restored genomic coordinates, before shared DP. */
	private void applySsuCompetition(ArrayList<Orf>[] rnas, byte[] bases){
		if(!call16S || !SsuOverlapResolver.hasClippedMito(rnas)){return;}
		try{
			if(ssuCompetition==null){ssuCompetition=new SsuOverlapResolver(
				Boolean.getBoolean("bbtools.diagnostic.ssuVeto"), Boolean.getBoolean("bbtools.verify.ssuVeto"));}
			ssuCompetition.apply(rnas, bases);
		}catch(RuntimeException | AssertionError e){KillSwitch.assertDie("SSU full-model competition failed: "+e);throw e;}
	}

	/** Observe the exact shared-DP boundary, distinct from verifier accepts that
	 * scavenger snapshot resolution may discard. Never exposes mutable Orfs. */
	private void recordNcrnaPath(ArrayList<Orf> list,String name,int length,boolean retained){
		assert(ncrnaDiagSink!=null && list!=null) : "Path observation is reached only when explicitly enabled";
		for(Orf orf : list){
			if(orf.ncrnaFamily==null){continue;}
			//makeRnas flips minus-strand calls to genomic coordinates before the DP.
			final int start=orf.strand==0 ? orf.start : length-1-orf.stop;
			final int stop=orf.strand==0 ? orf.stop : length-1-orf.start;
			if(orf.acceptedRrnaIdentity>=0){
				ncrnaDiagSink.path(orf.ncrnaFamily, name, orf.strand, length,
					orf.trnaModel, start, stop, retained, orf.acceptedRrnaIdentity);
			}else{
				ncrnaDiagSink.path(orf.ncrnaFamily, name, orf.strand, length,
					orf.trnaModel, start, stop, retained);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Generates lists of all max-length non-overlapping Orfs per frame.
	 * There IS overlap between frames.
	 * All Orfs come out flipped to + orientation. 
	 * */
	public ArrayList<Orf>[] makeRnas(String name, byte[] bases){
		return makeRnas(name, bases, null);
	}

	/** Supplies only pre-path CDS candidates; isolated RNA calls have no CDS context. */
	private ArrayList<Orf>[] makeRnas(String name, byte[] bases, ArrayList<Orf>[] cdsCandidates){
		if(ncrnaCpuProfile==null && NcrnaCpuProfile.enabled()){ncrnaCpuProfile=NcrnaCpuProfile.create(ncrnaFamilies.size());}
		@SuppressWarnings("unchecked")
		ArrayList<Orf>[] array=new ArrayList[2];
		array[0]=new ArrayList<Orf>();
		array[1]=new ArrayList<Orf>();
		//Only the legacy PGM branch below consumes these whole-contig arrays.
		//Allocate lazily there so a conserved-seed-only call needs neither one.
		float[] scores=null;
		int[] kmersSeen=null;
		initializeNcrnaScavengers();
		if(ncrnaScavengers!=null){
			for(NcrnaScavenger scavenger : ncrnaScavengers){
				scavenger.prepareCdsStopConstraint(cdsCandidates, bases.length);
			}
		}
		//The existing in-place flip at the end of each iteration makes bases
		//strand-oriented here and restores its original orientation before return.
		//Generic scavengers therefore share that buffer with the legacy callers;
		//a second reverse-complement copy would retain another whole contig per worker.
		for(int strand=0; strand<2; strand++){
			final byte[] strandBases=bases;
			final ConservedRnaSeedIndex.ScanResult sharedSeedHits;
			if(conservedRnaSeedIndex==null){sharedSeedHits=null;}
			else if(iterationSeeds!=null && iterationSeeds[strand]!=null){
				assert(iterationBases==bases) : "Seed reuse is confined to the same borrowed contig and strand";
				sharedSeedHits=iterationSeeds[strand];
			}else{
				if(ncrnaCpuProfile!=null){ncrnaCpuProfile.begin(ncrnaFamilies.size(), NcrnaCpuProfile.SHARED_SCAN);}
				try{sharedSeedHits=conservedRnaSeedIndex.scan(strandBases);}
				finally{if(ncrnaCpuProfile!=null){ncrnaCpuProfile.end();}}
				if(iterationSeeds!=null){iterationSeeds[strand]=sharedSeedHits;}
				if(strandBases.length>=ConservedRnaSeedIndex.K){sharedSweepPasses++;}
			}
			//Generic ncRNA scavenger pass (Noire's TrnaCaller factoring, 2026-08-23): one
			//independent claimed-region list per family per strand -- these families don't
			//overlap each other or tRNA biologically, so there's no reason to cross-claim.
			if(ncrnaScavengers!=null){
				for(int fi=0; fi<ncrnaScavengers.size(); fi++){
					if(inactiveClanMember(fi)){continue;}
					final int reuseSlot=2*fi+strand;
					if(iterationCandidates!=null && !iterationCandidates.isTarget(fi)){
						final ArrayList<Orf> retained=iterationCandidates.restoreRna(reuseSlot);
						if(retained!=null){
							array[strand].addAll(retained);
							continue;
						}
					}
					ArrayList<int[]> calledPos=new ArrayList<>();
					final NcrnaScavenger scavenger=ncrnaScavengers.get(fi);
					final int seedSlot=ncrnaSeedSlots[fi];
					ArrayList<Orf> found;
					if(ncrnaCpuProfile!=null){ncrnaCpuProfile.begin(fi, NcrnaCpuProfile.SCAVENGE);}
					try{
						if(seedSlot>=0){
							if(sharedSeedHits==null){
								throw new IllegalStateException("Missing shared 17-mer scan for "+
									ncrnaFamilies.get(fi).name);
							}
							final int[] hits=sharedSeedHits.hits(seedSlot);
							if(VERIFY_UNIFIED_SEEDS){
								verifySeedHits(name, ncrnaFamilies.get(fi).name, strand, hits,
									scavenger.findKmerHitPositions(strandBases));
							}
							found=scavenger.scavenge(name, strandBases, strand, calledPos, hits, sharedSeedHits.keys(seedSlot));
						}else{
							found=scavenger.scavenge(name, strandBases, strand, calledPos);
						}
					}finally{if(ncrnaCpuProfile!=null){ncrnaCpuProfile.end();}}
					if(strand==1 && found!=null){
						for(Orf orf : found){orf.flip();}
					}
					if(iterationCandidates!=null && !iterationCandidates.isTarget(fi)){
						iterationCandidates.captureRna(reuseSlot, found);
					}
					if(found!=null){array[strand].addAll(found);}
				}
			}
			for(StatsContainer sc : pgm.rnaContainers){
				if(ProkObject.callType(sc.type)){
					if(sc.type==tRNA && trnaLibrary!=null && trnaLibrary.length>0){
						final int reuseSlot=2*ncrnaFamilies.size()+strand;
						if(iterationCandidates!=null){
							final ArrayList<Orf> retained=iterationCandidates.restoreRna(reuseSlot);
							if(retained!=null){
								array[strand].addAll(retained);
								continue;
							}
						}
						if(ncrnaCpuProfile!=null){ncrnaCpuProfile.begin(ncrnaFamilies.size(), NcrnaCpuProfile.TRNA);}
						try{
							if(trnaCaller==null){
								trnaCaller=new TrnaCaller(pgm, trnaLibrary, trnaModels, trnaModelNames);
								trnaCaller.cmVerifier=cmVerifier;
							}
							ArrayList<Orf> list;
							if(TrnaCaller.SCAVENGE_ONLY){
								list=new ArrayList<>();
							}else{
								list=trnaCaller.callTrnas(name, bases, strand);
							}
							//Scavenger pass: find tRNAs at kmer-hit positions missed by PGM-based candidate generation
							if((TrnaCaller.SCAVENGE || TrnaCaller.SCAVENGE_ONLY) && list!=null){
								ArrayList<int[]> calledPos=new ArrayList<>();
								for(Orf orf : list){calledPos.add(new int[]{orf.start, orf.stop});}
								ArrayList<Orf> scavenged;
								if(trnaSeedSlot>=0){
									if(sharedSeedHits==null){
										throw new IllegalStateException("Missing shared 17-mer scan for tRNA");
									}
									final int[] hits=sharedSeedHits.hits(trnaSeedSlot);
									if(VERIFY_UNIFIED_SEEDS){
										verifySeedHits(name, "tRNA", strand, hits,
											trnaCaller.findKmerHitPositions(bases));
									}
									scavenged=trnaCaller.scavengeTrnas(name, bases, strand, calledPos, hits);
								}else{
									scavenged=trnaCaller.scavengeTrnas(name, bases, strand, calledPos);
								}
								if(scavenged!=null){list.addAll(scavenged);}
							}
							if(strand==1 && list!=null){
								for(Orf orf : list){orf.flip();}
							}
							if(iterationCandidates!=null){iterationCandidates.captureRna(reuseSlot, list);}
							if(list!=null){array[strand].addAll(list);}
						}finally{if(ncrnaCpuProfile!=null){ncrnaCpuProfile.end();}}
					}else{
						if(scores==null){
							scores=new float[bases.length];
							kmersSeen=(lsuKmers==null && ssuKmers==null && trnaKmers==null &&
								r5SKmers==null) ? null : new int[bases.length];
						}
						ArrayList<Orf> list=makeRnasForStrand(name, bases, strand, sc, scores,
							(sc.kmerSet()==null ? null : kmersSeen), false, -1);//TODO: Make this loop through all RNA types
						if(strand==1 && list!=null){
							for(Orf orf : list){
								assert(orf.strand==strand);
								orf.flip();
							}
						}
						if(list!=null){array[strand].addAll(list);}
					}
				}
			}
			Collections.sort(array[strand]);
			Vector.reverseComplementInPlaceFast(bases);
		}
		return array;
	}

	/** Returns this worker's live family caller after native resources have loaded.
	 * Policy iteration changes this instance, never the shared family descriptor.
	 */
	public NcrnaScavenger ncrnaScavenger(final int family){
		if(family<0 || family>=ncrnaFamilies.size()){
			throw new IllegalArgumentException("Unknown loaded ncRNA family index: "+family);
		}
		initializeNcrnaScavengers();
		return ncrnaScavengers.get(family);
	}

	/** Optional target counters for the owning iteration worker between calls. */
	public NcrnaAlignmentWork ncrnaAlignmentWork(final int family){
		return ncrnaCpuProfile==null ? null : ncrnaCpuProfile.alignmentWork(family);
	}

	/** Shares only immutable seed hits while a driver repeats policies on one
	 * unchanged contig. Every policy still generates its own candidates and path.
	 * Pair with endNcrnaIteration in finally; no sequence survives that scope.
	 */
	public void beginNcrnaIteration(final Read read){
		if(read==null || read.bases==null || iterationSeeds!=null){
			throw new IllegalStateException("A seed-reuse scope needs one intact contig and cannot nest");
		}
		iterationBases=read.bases;
		iterationSeeds=new ConservedRnaSeedIndex.ScanResult[2];
	}

	/** Optional reuse of candidates independent of the target shortlist. Limited
	 * to the measured prokaryotic callers; ordinary callers never enter this scope.
	 * Family claims remain independent in makeRnas. TPP still consumes the exact
	 * original pre-path CDS starts, and every policy reruns shared path selection.
	 */
	public void beginNcrnaIteration(final Read read, final int target){
		beginNcrnaIteration(read, new int[]{target});
	}

	/** Recompute every selected member while reusing unchanged callers' candidates.
	 * Copies target membership and preserves the single-target dependency guards.
	 */
	public void beginNcrnaIteration(final Read read, final int[] targets){
		if(proteinSupport!=null || proteinFateSink!=null || stopTailSupport!=null || tailCascade!=null
				|| cmVerifier!=null || ssuCompetition!=null
			|| call16S || call23S || call5S || call18S
			|| (calltRNA && (trnaLibrary==null || trnaLibrary.length==0))){
			throw new IllegalArgumentException("Candidate reuse requires seeded prok RNA callers without "+
				"protein/CM/SSU-competition state");
		}
		for(final NcrnaFamily family:ncrnaFamilies){
			if(family.outputType==r18S || family.name.equals("r58") || family.name.equals("lsu") ||
					family.name.equalsIgnoreCase("euk5S")){
				throw new IllegalArgumentException("Eukaryotic candidate reuse needs separate dependency "+
					"verification: "+family.name);
			}
		}
		final NcrnaCandidateReuse candidates=new NcrnaCandidateReuse(targets, ncrnaFamilies.size());
		beginNcrnaIteration(read);
		iterationCandidates=candidates;
	}

	public void endNcrnaIteration(){
		iterationBases=null;
		iterationSeeds=null;
		iterationCandidates=null;
	}

	/** Lazy construction stays on the owning worker, including CPU-clock spans. */
	private void initializeNcrnaScavengers(){
		if(ncrnaScavengers==null && cobalaminCallerResources!=null && !cobalaminCallerResources.matches(ncrnaFamilies)){
			cobalaminCallerResources=null;
		}
		if(ncrnaScavengers==null && !cobalaminClanConfigured && cobalaminCallerResources!=null && CallGenes.COBALAMIN_ENABLED){
			cobalaminServing=cobalaminCallerResources;
			setCobalaminClan(true, new NcrnaDdlProfileResource[]{cobalaminServing.memberProfiles()}, false);
		}
		if(ncrnaCpuProfile==null && NcrnaCpuProfile.enabled()){
			ncrnaCpuProfile=NcrnaCpuProfile.create(ncrnaFamilies.size());
		}
		if(ncrnaScavengers==null && !ncrnaFamilies.isEmpty()){
			if(cobalaminClanConfigured){
				cobalaminMemberSlots=new int[]{-1, -1, -1};
				final String[] names={"cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689"};
				for(int i=0; i<ncrnaFamilies.size(); i++){
					for(int member=0; member<3; member++){
						if(names[member].equals(ncrnaFamilies.get(i).name)){
							if(cobalaminMemberSlots[member]>=0){throw new IllegalArgumentException("Repeated clan member registration");}
							cobalaminMemberSlots[member]=i;
						}
					}
				}
				for(final int slot:cobalaminMemberSlots){
					if(slot<0){throw new IllegalArgumentException("The clan experiment requires all three loaded cobalamin members");}
				}
				if(cobalaminJoint && (iterationCandidates!=null || cmVerifier!=null
					|| boundaryTrace!=null && !boundaryTrace.invocationsOnly()
					|| ncrnaDiagSink!=null || ncrnaVoteDiagSink!=null || CallGenes.NCRNA_REFRESH_CLAIMED_WINDOWS)){
					throw new IllegalArgumentException("Joint cobalamin requires frozen passes without candidate reuse or CM verification");
				}
			}
			ncrnaScavengers=new ArrayList<>(ncrnaFamilies.size());
			for(NcrnaFamily fam : ncrnaFamilies){
				if(ncrnaCpuProfile!=null){ncrnaCpuProfile.begin(ncrnaScavengers.size(), NcrnaCpuProfile.INIT);}
				try{
				final boolean sketchOnly=cobalaminDdlProfiles!=null && cobalaminDdlProfiles.dual() && !cobalaminDdlProfiles.orderOnly()
					&& "cobalamin_rf00174".equals(fam.name);
				NcrnaScavenger scavenger=new NcrnaScavenger(fam.library, fam.models, fam.modelNames,
					fam.kmerSet, fam.kLong, fam.minLen, fam.windowPad,
					fam.indexK, fam.indexTopN, fam.adaptive,
					fam.adaptFloor, fam.adaptTopFrac, fam.adaptQFrac, fam.fixedMinHits,
					fam.scoreA, fam.scoreB, fam.idPass, fam.idBorderline,
					fam.boundary5NetTemplate, fam.boundary3NetTemplate,
					fam.boundaryStartTable, fam.boundaryStopTable,
					fam.boundaryStartInside, fam.boundaryStartOutside, fam.boundaryStopInside, fam.boundaryStopOutside,
					fam.boundaryMeanLen, fam.boundaryStartOffsets, fam.boundaryStopOffsets,
					fam.boundaryMarginStart, fam.boundaryMarginStop, !sketchOnly
						&& !(cobalaminJoint && isCobalaminMember(fam.name)));
				applyNcrnaFamilyControls(scavenger, fam);
				if(boundaryTrace!=null && boundaryTrace.accepts(fam)){
					scavenger.setBoundaryTrace(boundaryTrace.family(fam));
				}
				// D50 keeps18S on its fast path even when small-family CM is enabled.
				scavenger.cmVerifier=(fam.outputType!=r18S && cmVerifier!=null && cmVerifier.supports(fam.outputType) ? cmVerifier : null);
				if(cobalaminDdlProfiles!=null && "cobalamin_rf00174".equals(fam.name)){
					scavenger.useDdlProfiles(cobalaminDdlProfiles, cobalaminDdlCap);
				}
				scavenger.setStageDiagSink(ncrnaDiagSink,fam.name);
				scavenger.setVoteDiagSink(ncrnaVoteDiagSink, fam.name);
				scavenger.setAttemptLog(Euk18sAttemptLog.forFamily(fam.name));
				ncrnaScavengers.add(scavenger);
				}finally{if(ncrnaCpuProfile!=null){ncrnaCpuProfile.end();}}
			}
			if(cobalaminClanConfigured){initializeCobalaminClan();}
			if(pathScoreBasis!=null){installPathScoreBasis();}
		}
	}

	private static boolean isCobalaminMember(final String name){
		return "cobalamin_rf00174".equals(name) || "cobalamin_rf01482".equals(name) || "cobalamin_rf01689".equals(name);
	}

	/** Three registry identities remain intact; only the main execution slot owns
	 * the shared walk. Resource-only minor contexts never enter the caller loop. */
	private void initializeCobalaminClan(){
		assert(cobalaminMemberSlots!=null && cobalaminMemberSlots.length==3)
			: "Clan initialization follows complete loaded-member binding";
		if(cobalaminJoint){
			final NcrnaFamily[] families=new NcrnaFamily[3];
			final NcrnaScavenger[] refiners=new NcrnaScavenger[3];
			for(int member=0; member<3; member++){
				families[member]=ncrnaFamilies.get(cobalaminMemberSlots[member]);
				refiners[member]=ncrnaScavengers.get(cobalaminMemberSlots[member]);
			}
			final int slot=cobalaminMemberSlots[0];
			if(ncrnaCpuProfile!=null){ncrnaCpuProfile.begin(slot, NcrnaCpuProfile.INIT);}
			try{
				final NcrnaClanScavenger clan=new NcrnaClanScavenger(families, refiners, cobalaminNativeRanking);
				if(cobalaminMemberProfiles!=null){
					if(cobalaminNativeRanking){clan.useMemberDdlProfiles(cobalaminMemberProfiles[0]);}
					else{clan.useDdlOnlyProfiles(cobalaminMemberProfiles[0]);}
				}
				if(cobalaminServing!=null){cobalaminServing.configure(clan);}
				ncrnaScavengers.set(slot, clan);
			}finally{if(ncrnaCpuProfile!=null){ncrnaCpuProfile.end();}}
			if(ncrnaCpuProfile!=null){ncrnaCpuProfile.executionName(slot, "cobalamin_clan");}
		}else if(cobalaminMemberProfiles!=null){
			for(int member=0; member<3; member++){
				final NcrnaScavenger worker=ncrnaScavengers.get(cobalaminMemberSlots[member]);
				if(member==0){worker.useDdlProfiles(cobalaminMemberProfiles[member], worker.indexTopN);}
				else{worker.useMemberDdlProfiles(cobalaminMemberProfiles[member]);}
			}
		}
	}

	public static void applyNcrnaProductionControls(NcrnaScavenger scavenger){
		scavenger.refreshClaimedWindows=CallGenes.NCRNA_REFRESH_CLAIMED_WINDOWS;
	}

	/** Attach only to a fresh worker; source ordinals belong to the actual input scheduler. */
	public void setBoundaryTrace(final NcrnaBoundaryTrace.Worker trace){
		if(ncrnaScavengers!=null || boundaryTrace!=null){throw new IllegalStateException("Configure boundary tracing before lazy caller initialization");}
		boundaryTrace=trace;
	}

	public static void applyNcrnaFamilyControls(NcrnaScavenger scavenger, NcrnaFamily family){
		if(family.modelPatience!=-1){
			if(scavenger.hasNativeRankingIndex()){
				scavenger.setShortlistPolicy(family.indexTopN, family.adaptive, family.adaptFloor,
					family.adaptTopFrac, family.adaptQFrac, family.fixedMinHits, family.modelPatience);
			}else{scavenger.setModelPatience(family.modelPatience);}
		}
		scavenger.hbmPass=family.hbmPass;
		scavenger.collapseFrac=family.collapseFrac;
		scavenger.family=family.name;
		if("tmrna".equals(family.name)){
			scavenger.setTmrnaRescueMargins(CallGenes.TMRNA_RESCUE_ID_MARGIN,CallGenes.TMRNA_RESCUE_HBM_MARGIN);
		}
		scavenger.voteTable=family.voteTable; scavenger.voteWindows=family.voteWindows; scavenger.voteEnds=family.voteEnds;
		scavenger.voteBeforePadding=family.voteBeforePadding;
		scavenger.preserveContainingWindows=family.preserveContainingWindows;
		scavenger.collapsedWindowPad=family.collapsedWindowPad;
		scavenger.voteSlack=family.voteSlack; scavenger.voteEndsMaxSd=family.voteEndsMaxSd;
		scavenger.outputType=family.outputType;
		scavenger.minKmerHits=family.seedMinHits;
		scavenger.seedDistinct=family.seedDistinct;
		scavenger.maxLen=family.maxLen;
		scavenger.quantumThresh=family.quantumThresh;
		scavenger.scavengePass2=family.scavengePass2;
		scavenger.rankedModelFallback=family.rankedModelFallback;
		scavenger.strictIndexCutoff=family.strictIndexCutoff;
		scavenger.indexScoreMargin=family.indexScoreMargin;
		scavenger.trimAlignmentExtent=family.trimAlignmentExtent;
		scavenger.reuseConsensusAlignment=family.reuseConsensusAlignment;
		if(family.quantumCapMultiplier!=1 && (!family.reuseConsensusAlignment
				|| !"23S".equals(family.name) || family.quantumCapMultiplier!=2)){
			throw new IllegalArgumentException("Only reused23S verification may select Quantum cap2: "+family.name);
		}
		scavenger.quantumCapMultiplier=family.quantumCapMultiplier;
		scavenger.pacBioConsensusAlignment=family.pacBioConsensusAlignment;
		scavenger.pacBioCosts=family.pacBioCosts;
		scavenger.pacBioRolling=family.pacBioRolling;
		scavenger.joinedProposals=family.joinedProposals;
		scavenger.joinedPositions=family.joinedPositions;
		scavenger.joinedSideSupport=family.joinedSideSupport;
		scavenger.modelEndClipping=family.modelEndClipping;
		scavenger.modelClipRescue=family.modelClipRescue;
		scavenger.prok16sPartialMinLen=family.prok16sPartialMinLen;
		scavenger.physicalEdgeMinLen=family.physicalEdgeMinLen;
		scavenger.modelClipRescueId=family.modelClipRescueId;
		if(family.modelThresholds!=null){scavenger.setModelThresholds(family.modelThresholds);}
		scavenger.boundaryFeatureVersion=family.boundaryFeatureVersion;
		scavenger.boundarySharedConfidence=family.boundarySharedConfidence;
		scavenger.boundaryOnRawEndpoints=family.boundaryOnRawEndpoints;
		scavenger.deferBoundaryReport=family.deferBoundaryReport;
		scavenger.nnCentre5Vote=family.nnCentre5Vote;
		scavenger.nnCentre3Vote=family.nnCentre3Vote;
		scavenger.nnCutoff=family.nnCutoff;
		scavenger.nnCutoffEnabled=family.nnCutoffEnabled;
		if(family.boundaryNetsByModel!=null){
			scavenger.setPerModelBoundaryResources(family.boundaryNetsByModel,
				family.boundaryStartTablesByModel, family.boundaryStopTablesByModel);
		}
		if(family.rrnaEndpointFeatures!=null){
			scavenger.setRrnaEndpointFeatures(family.rrnaEndpointFeatures,
				family.rrnaEndpointFeatureSink);
		}
		if(family.boundaryRefinerFactory!=null){
			scavenger.setBoundaryRefiner(family.boundaryRefinerFactory.create(family));
		}
		if(family.cdsStopConstraintFactory!=null){
			if(family.boundaryRefinerFactory!=null || family.rrnaEndpointFeatures!=null){
				throw new IllegalArgumentException("CDS stop constraint requires the ordinary ncRNA boundary search: "+family.name);
			}
			scavenger.setCdsStopConstraint(family.cdsStopConstraintFactory.create());
		}
		applyNcrnaProductionControls(scavenger);
	}

	/** Builds the single shared 17-mer index after all requested resources have loaded and
	 * before worker threads start. Non-17-mer experimental/legacy paths retain their existing
	 * scanners. */
	public static synchronized void initializeConservedRnaSeedIndex(){
		final ArrayList<LongHashSet> seedSets=new ArrayList<>();
		final ArrayList<Boolean> retainKeys=new ArrayList<>();
		ncrnaSeedSlots=new int[ncrnaFamilies.size()];
		Arrays.fill(ncrnaSeedSlots, -1);
		for(int i=0; i<ncrnaFamilies.size(); i++){
			final NcrnaFamily family=ncrnaFamilies.get(i);
			if(family.kLong==ConservedRnaSeedIndex.K && family.kmerSet!=null){
				ncrnaSeedSlots[i]=seedSets.size();
				seedSets.add(family.kmerSet);
				retainKeys.add(family.joinedPositions!=null || family.voteTable!=null && (family.voteWindows || family.voteEnds));
			}
		}
		trnaSeedSlot=-1;
		if(calltRNA && (TrnaCaller.SCAVENGE || TrnaCaller.SCAVENGE_ONLY)
				&& kLongTRna==ConservedRnaSeedIndex.K && trnaKmers!=null){
			trnaSeedSlot=seedSets.size();
			seedSets.add(trnaKmers);
			retainKeys.add(false);
		}
		if(seedSets.isEmpty()){conservedRnaSeedIndex=null;}
		else{
			final boolean[] retain=new boolean[retainKeys.size()]; for(int i=0;i<retain.length;i++){retain[i]=retainKeys.get(i);}
			conservedRnaSeedIndex=new ConservedRnaSeedIndex(seedSets.toArray(new LongHashSet[seedSets.size()]), retain);
		}
	}

	public static int conservedRnaSeedSlotCount(){
		return conservedRnaSeedIndex==null ? 0 : conservedRnaSeedIndex.slotCount();
	}

	private static void verifySeedHits(String contig, String family, int strand, int[] unified, int[] legacy){
		if(!Arrays.equals(unified, legacy)){
			throw new IllegalStateException("Unified 17-mer hit-stream mismatch: contig="+contig
				+", family="+family+", strand="+strand+", unified="+Arrays.toString(unified)
				+", legacy="+Arrays.toString(legacy));
		}
	}
	
	/** Designed for quickly calling a single SSU */
	public Orf makeRna(String name, byte[] bases, int type){
		final float[] scores=new float[bases.length];//TODO: Big and slow; make a FloatList?
		StatsContainer sc=pgm.allContainers[type];
		final int[] kmersSeen=(sc.kmerSet()==null ? null : new int[bases.length]);//TODO: IntList?
		
		int strand=0;
		ArrayList<Orf> list=makeRnasForStrand(name, bases, strand, sc, scores, kmersSeen, true, -1);
		final Orf best1=pickBest(list);
		assert(best1==null || best1.start>=0 && best1.stop<bases.length) : bases.length+"\n"+best1;
		if(best1!=null && best1.orfScore>-999){return best1;}
		
		strand++;
		Vector.reverseComplementInPlaceFast(bases);
		list=makeRnasForStrand(name, bases, strand, sc, scores, kmersSeen, true, -1);
		Vector.reverseComplementInPlaceFast(bases);
		if(strand==1 && list!=null){
			for(Orf orf : list){
				assert(orf.strand==strand);
				orf.flip();
			}
		}
		final Orf best2=pickBest(list);
		assert(best2==null || best2.start>=0 && best2.stop<bases.length) : bases.length+"\n"+best2;
		if(best2!=null && best2.orfScore>-999){return best2;}
		return best1!=null ? best1 : best2;
	}
	
	/**
	 * Selects the highest-scoring ORF from a list.
	 * @param list List of ORFs to evaluate
	 * @return ORF with highest score, or null if list is empty
	 */
	final Orf pickBest(ArrayList<Orf> list){
		if(list==null){return null;}
		Orf best=null;
		for(Orf orf : list){
			if(best==null || orf.orfScore>best.orfScore){
				best=orf;
			}
		}
		return best;
	}
	
	/** 
	 * Generates lists of all max-length non-overlapping Orfs per frame.
	 * There IS overlap between frames.
	 * All Orfs come out flipped to + orientation. 
	 * */
	static ArrayList<Orf>[] makeOrfs(String name, byte[] bases, int minlen){
		return makeOrfs(name, bases, minlen, null);
	}

	/** Uses one explicit code on both strands while preserving the caller's input bases. */
	static ArrayList<Orf>[] makeOrfs(String name, byte[] bases, int minlen, GeneticCode geneticCode){
		@SuppressWarnings("unchecked")
		int orfs=0;
		ArrayList<Orf>[] array=new ArrayList[6];
		for(int strand=0; strand<2; strand++){
			for(int frame=0; frame<3; frame++){
				ArrayList<Orf> list=makeOrfsForFrame(name, bases, frame, strand, minlen, geneticCode);
				orfs+=(list==null ? 0 : list.size());
				array[frame+3*strand]=list;
				if(strand==1 && list!=null){
					for(Orf orf : list){
						assert(orf.frame==frame);
						assert(orf.strand==strand);
						orf.flip();
					}
				}
			}
			Vector.reverseComplementInPlaceFast(bases);
		}
//		System.err.print("orfs="+orfs+", ");
		return array;
	}
	
	/** Run the actual bounded path search on six caller-ordered, fresh retained lists.
	 * Mutates candidate path scores/links and output counters. Exposed for diagnostic
	 * fixtures; it does not generate, score or prune candidates. With an observer,
	 * all candidates must belong to this contig and have no previous path links. */
	public ArrayList<Orf> findPath(ArrayList<Orf>[] frameLists, byte[] bases){
		// Sequence is used only for the observation length. Preserve the legacy
		// unobserved null-bases case while requiring real sequence when observed.
		return findPath(frameLists, pathObserver!=null ? bases.length : 0);
	}

	/** Searches saved candidate metadata with the native DP and no dummy sequence.
	 * The caller and all effective scoring settings must match the captured graph.
	 * Candidates are mutable, caller-owned and fresh for this trial, as above.
	 */
	public ArrayList<Orf> findPath(ArrayList<Orf>[] frameLists, int contigLength){
		if(contigLength<0){throw new IllegalArgumentException("Path observations require a nonnegative contig length: "+contigLength);}
		if(pathActive){throw new IllegalStateException("The path search is not reentrant");}
		pathActive=true;
		try{
			final float penalty=tailCascade!=null && tailCascadeBonus>0 ? tailCascadeOppositePenalty : 0;
			final TailCascadeOppositePenalty.Result opposite;
			if(penalty>0){
				opposite=TailCascadeOppositePenalty.apply(frameLists, penalty, maxOverlapOppositeStrand,
					oppositePenaltyObserver!=null);
			}
			else{TailCascadeOppositePenalty.clear(frameLists); opposite=null;}
			pathTrace=pathObserver==null ? null :
				new GeneCallerPathTrace(pathObserver, frameLists, contigLength, tailCascade!=null, pathScoreBasis);
			final ArrayList<Orf> selected=findPathInner(frameLists);
			if(opposite!=null){opposite.report(selected, oppositePenaltyObserver);}
			return selected;
		}finally{pathTrace=null; pathActive=false;}
	}

	private ArrayList<Orf> findPathInner(ArrayList<Orf>[] frameLists){
		ArrayList<Orf> all=new ArrayList<Orf>();
		for(ArrayList<Orf> list : frameLists){all.addAll(list);}
		if(all.isEmpty()){
			if(pathTrace!=null){pathTrace.finish(null, 0);}
			return all;
		}
		Collections.sort(all);
		
		for(Orf orf : all){
			orf.pathScorePlus=-999999;
			orf.pathScoreMinus=-999999;
		}
		
		int[] lastPositionScored=KillSwitch.allocInt1D(6);
		Arrays.fill(lastPositionScored, -1);

		//Index of highest-scoring ORF in this frame, with prev on the plus strand
		int[] bestIndexPlus=KillSwitch.allocInt1D(6);
		//Index of highest-scoring ORF in this frame, with prev on the minus strand
		int[] bestIndexMinus=KillSwitch.allocInt1D(6);
		//Highest-scoring ORF in this frame, with prev on the plus strand
		Orf[] bestOrfPlus=new Orf[6];
		//Highest-scoring ORF in this frame, with prev on the minus strand
		Orf[] bestOrfMinus=new Orf[6];

		int[][] bestIndex=new int[][] {bestIndexPlus, bestIndexMinus};
		Orf[][] bestOrf=new Orf[][] {bestOrfPlus, bestOrfMinus};
		
		for(Orf orf : all){
			final int myListNum=3*orf.strand+orf.frame;
			if(pathTrace!=null){pathTrace.start(orf, lastPositionScored[myListNum]+1);}
			calcPathScore(orf, frameLists, lastPositionScored, bestIndex);
			if(pathTrace!=null){pathTrace.scored(orf);}
			if(bestOrfPlus[myListNum]==null || orf.pathScorePlus>=bestOrfPlus[myListNum].pathScorePlus){
				bestOrfPlus[myListNum]=orf;
				bestIndexPlus[myListNum]=lastPositionScored[myListNum];
				assert(frameLists[myListNum].get(lastPositionScored[myListNum])==orf);
			}
			if(bestOrfMinus[myListNum]==null || orf.pathScoreMinus>=bestOrfMinus[myListNum].pathScoreMinus){
				bestOrfMinus[myListNum]=orf;
				bestIndexMinus[myListNum]=lastPositionScored[myListNum];
				assert(frameLists[myListNum].get(lastPositionScored[myListNum])==orf);
			}
		}
		
		Orf best=bestOrf[0][0];
		for(Orf[] array : bestOrf){
			for(Orf orf : array){
				if(best==null || (orf!=null && orf.pathScore()>best.pathScore())){
					best=orf;
				}
			}
		}
		//best = global max of pathScore() over both prev-strand tables
		//(bestOrfPlus/bestOrfMinus) across all 6 frames; the chain is reconstructed
		//by walking orf.prev() backward. pathScorePlus/pathScoreMinus are tracked
		//per-Orf separately precisely so the DP can pick a predecessor on either
		//strand without clobbering the other strand's running chain.
		ArrayList<Orf> bestPath=new ArrayList<Orf>();
		for(Orf orf=best; orf!=null; orf=orf.prev()){
			bestPath.add(orf);
			if(orf.type==CDS){geneStartsOut++;}
			else if(orf.type==tRNA){tRNAOut++;}
			else if(orf.type==r16S){r16SOut++;}
			else if(orf.type==r23S){r23SOut++;}
			else if(orf.type==r5S){r5SOut++;}
			else if(orf.type==r18S){r18SOut++;}
			else if(orf.type==RNA){rnaOut++;}
			if(orf.ncrnaFamily!=null){
				if(ncrnaOutByFamily==null){ncrnaOutByFamily=new long[ncrnaFamilies.size()];}
				int familyIndex=-1;
				for(int i=0; i<ncrnaFamilies.size(); i++){
					if(orf.ncrnaFamily.equals(ncrnaFamilies.get(i).name)){familyIndex=i; break;}
				}
				if(familyIndex<0){throw new IllegalStateException("Unknown committed ncRNA family: "+orf.ncrnaFamily);}
				ncrnaOutByFamily[familyIndex]++;
			}
		}
		Collections.sort(bestPath);
		if(pathTrace!=null){pathTrace.finish(best, bestPath.size());}
		return bestPath;
	}
	
	/**
	 * Calculate the best path to this ORF.
	 * @param orf
	 * @param frameLists
	 * @param lastPositionScored
	 * @param bestIndex
	 */
	private void calcPathScore(Orf orf, ArrayList<Orf>[] frameLists, int[] lastPositionScored, int[][] bestIndex){
		final int myListNum=3*orf.strand+orf.frame;

//		System.err.println("* "+orf);
//		System.err.println("* "+Arrays.toString(lastPositionScored));
//		System.err.println();
		
		for(int listStrand=0; listStrand<2; listStrand++){
			for(int listFrame=0; listFrame<3; listFrame++){
				int listNum=listFrame+3*listStrand;
				ArrayList<Orf> list=frameLists[listNum];
				int lastPos=lastPositionScored[listNum];
				int bestPos=bestIndex[listStrand][listNum];
				if(listStrand==0){
					calcPathScorePlus(orf, list, listStrand, lastPos, bestPos, listNum);
				}else{
					calcPathScoreMinus(orf, list, listStrand, lastPos, bestPos, listNum);
				}
			}
		}
		
//		System.err.println(myListNum+", "+Arrays.toString(lastPositionScored)+", "+frameLists[myListNum].size());
		
		lastPositionScored[myListNum]++;
		assert(frameLists[myListNum].get(lastPositionScored[myListNum])==orf) :
			myListNum+"\n"+orf+"\n"+frameLists[myListNum].get(lastPositionScored[myListNum])+"\n"
			+Arrays.toString(lastPositionScored)+"\n"+frameLists[myListNum].get(lastPositionScored[myListNum]+1);
		
		//These are sanity checks to make sure that the path did not break in the middle.
		//Safe to disable.
//		assert(orf.prevPlus!=null || orf.stop<100000);
//		assert(orf.prevMinus!=null || orf.stop<100000);
//		assert(orf.pathScore>-10) : orf.pathScore+"\n"+orf+"\n"+orf.prev+"\n";
	}
	
	/**
	 * Calculate the best path to this ORF from a plus-strand previous ORF.
	 * @param orf
	 * @param list
	 * @param listStrand
	 * @param lastPos
	 * @param bestPos
	 */
	private void calcPathScorePlus(final Orf orf, final ArrayList<Orf> list, final int listStrand,
			final int lastPos, final int bestPos, final int listNum){
		assert(listStrand==0);
		if(lastPos<0){
			if(pathTrace!=null){pathTrace.window(listNum, lastPos, bestPos, 0, -1);}
			if(orf.prevPlus==null){
				orf.pathScorePlus=pathRootScore(orf);
				orf.pathLengthPlus=1;
			}
			return;
		}
		if(list.isEmpty()){return;}
		
//		System.err.println("\nExamining   \t"+orf+"\nlastPos="+lastPos+", bestPos="+bestPos+", sameFrame="+sameFrame);
		boolean found=false;
		final boolean sameStrand=(orf.strand==listStrand);
		final int maxOverlap=(sameStrand ? maxOverlapSameStrand : maxOverlapOppositeStrand);
		final int min=Tools.max(0, bestPos-lookbackPlus);
		int lastVisited=-1;
		for(int i=lastPos; i>=min || (i>0 && !found); i--){
			lastVisited=i;
			Orf prev=list.get(i);
			assert(prev!=orf) : prev;
//			System.err.println("Comparing to \t"+prev);
			if(orf.isValidPrev(prev, maxOverlap)){
				int overlap=Tools.max(0, prev.stop-orf.start+1);
				//TODO: Probable bug - seeded rRNA uses r5S/r16S/r23S output types,
				//not generic RNA, so positive overlap discards NcrnaScavenger's
				//A+B*length*identity^2 score and applies the CDS heuristic to its
				//unpopulated start/stop/kmer fields. Native5S trace (2026-10-13,
				//Thermococcus celer) measured69.70996->4.5397615 at12bp overlap.
				//Only annotated here: the current assignment is observation-only.
				float orfScore=proteinAwareScore(orf, overlap);
				
				final float prevScore=prev.pathScore();
				final int prevLength=prev.pathLength();
				
				float pathScore;
				final int pathLength;
				if(sameStrand){
					pathLength=prevLength+1;
					pathScore=prevScore+orfScore;
					pathScore+=(FLAT_RNA_TRANSITIONS && orf.type!=CDS ? p0 :
						p0+p1*(Tools.mid(p5*(p2+pathLength), p6*(p3-pathLength), p4)));
					if(EUK5S_ARRAY_BONUS && euk5sArrayEdge(prev,orf,EUK5S_ARRAY_MAX_GAP)){pathScore+=30f;}
				}else{
					pathLength=1;
					pathScore=prev.pathScore()+orfScore;
					pathScore+=(FLAT_RNA_TRANSITIONS && orf.type!=CDS ? p0 :
						q1+Tools.mid(q2*prevLength, q3+q4*prevLength, q5));
				}
				
				if(overlap<1 && prevScore>0){found=true;}
				if(pathTrace!=null){pathTrace.edge(listNum, i, overlap, prevScore, prevLength,
					orfScore, pathScore, pathLength, pathScore>=orf.pathScorePlus);}
				if(pathScore>=orf.pathScorePlus){
					orf.pathScorePlus=pathScore;
					orf.prevPlus=prev;
					orf.pathLengthPlus=pathLength;
//					System.err.println("Set as best");
				}
			}
//			if(found && prev.stop<maxOverlap-2000 && orf.prevPlus!=null){
//				System.err.println("Breaking");
//				break;
//			}
		}
		if(pathTrace!=null){pathTrace.window(listNum, lastPos, bestPos, min, lastVisited);}
	}
	
	/**
	 * Calculate the best path to this ORF from a minus-strand previous ORF.
	 * @param orf
	 * @param list
	 * @param listStrand
	 * @param lastPos
	 * @param bestPos
	 */
	private void calcPathScoreMinus(final Orf orf, final ArrayList<Orf> list, final int listStrand,
			final int lastPos, final int bestPos, final int listNum){
		assert(listStrand==1);
		if(lastPos<0){
			if(pathTrace!=null){pathTrace.window(listNum, lastPos, bestPos, 0, -1);}
			if(orf.prevMinus==null){
				orf.pathScoreMinus=pathRootScore(orf);
				orf.pathLengthMinus=1;
			}
			return;
		}
		if(list.isEmpty()){return;}
		
//		System.err.println("\nExamining   \t"+orf+"\nlastPos="+lastPos+", bestPos="+bestPos+", sameFrame="+sameFrame);
		boolean found=false;
		final boolean sameStrand=(orf.strand==listStrand);
		final int maxOverlap=(sameStrand ? maxOverlapSameStrand : maxOverlapOppositeStrand);
		final int min=Tools.max(0, bestPos-lookbackMinus);
		int lastVisited=-1;
		for(int i=lastPos; i>=min || (i>0 && !found); i--){
			lastVisited=i;
			Orf prev=list.get(i);
			assert(prev!=orf) : prev;
//			System.err.println("Comparing to \t"+prev);
			if(orf.isValidPrev(prev, maxOverlap)){
				int overlap=Tools.max(0, prev.stop-orf.start+1);
				//TODO: Same seeded-rRNA/CDS overlap-rescoring issue documented in
				//calcPathScorePlus; the native12bp5S example reaches this branch.
				float orfScore=proteinAwareScore(orf, overlap);
				
				final float prevScore=prev.pathScore();
				final int prevLength=prev.pathLength();
				
				float pathScore;
				final int pathLength;
				if(sameStrand){
					pathLength=prevLength+1;
					pathScore=prevScore+orfScore;
					pathScore+=(FLAT_RNA_TRANSITIONS && orf.type!=CDS ? p0 :
						p0+p1*(Tools.mid(p5*(p2+pathLength), p6*(p3-pathLength), p4)));
					if(EUK5S_ARRAY_BONUS && euk5sArrayEdge(prev,orf,EUK5S_ARRAY_MAX_GAP)){pathScore+=30f;}
				}else{
					pathLength=1;
					pathScore=prev.pathScore()+orfScore;
					pathScore+=(FLAT_RNA_TRANSITIONS && orf.type!=CDS ? p0 :
						q1+Tools.mid(q2*prevLength, q3+q4*prevLength, q5));
				}
				if(overlap<1 && prevScore>0){found=true;}
				if(pathTrace!=null){pathTrace.edge(listNum, i, overlap, prevScore, prevLength,
					orfScore, pathScore, pathLength, pathScore>=orf.pathScoreMinus);}
				if(pathScore>=orf.pathScoreMinus){
					orf.pathScoreMinus=pathScore;
					orf.prevMinus=prev;
					orf.pathLengthMinus=pathLength;
//					System.err.println("Set as best");
				}
			}
//			if(found && prev.stop<maxOverlap-2000 && orf.prevMinus!=null){
//				System.err.println("Breaking");
//				break;
//			}
		}
		if(pathTrace!=null){pathTrace.window(listNum, lastPos, bestPos, min, lastVisited);}
	}
	
	/** Eligibility of this proposed immediate DP predecessor edge, not a search
	 * for nearby RNA across other retained genes. Genomic inclusive coordinates
	 * apply on both strands; overlapping calls never receive the diagnostic bonus. */
	public static boolean euk5sArrayEdge(Orf prev,Orf next,int maxGap){
		assert(prev!=null && next!=null && maxGap>=0) :
			"A proposed DP edge needs two real candidates and a nonnegative genomic gap limit";
		// addEuk5sFamily configures outputType=r5S, shared with legacy5S;
		// the explicit family tag, not the numeric output type, defines eligibility.
		if(prev.strand!=next.strand
			|| !"euk5S".equals(prev.ncrnaFamily) || !"euk5S".equals(next.ncrnaFamily)){return false;}
		assert(prev.scafName.equals(next.scafName)) : "GeneCaller DP predecessors must belong to the same contig";
		final long gap=(long)next.start-prev.stop-1;
		return gap>=0 && gap<=maxGap;
	}

	/** 
	 * Generates a list of maximal-length Orfs only (non-overlapping).
	 * All Orfs come out in native orientation (unflipped). 
	 * */
	static ArrayList<Orf> makeOrfsForFrame(String name, byte[] bases, int startFrame, int strand, int minlen){
		return makeOrfsForFrame(name, bases, startFrame, strand, minlen, null);
	}

	/** Enumerates selected-code stops and starts; edge/gap provenance survives later start selection. */
	public static ArrayList<Orf> makeOrfsForFrame(String name, byte[] bases, int startFrame,
			int strand, int minlen, GeneticCode geneticCode){
//		assert(false) : "TODO";
		// TODO: Probable CLI-validation bug - CallGenes accepts minlen=1 or 2,
		// then worker threads fail here instead of rejecting the parameter up front.
		assert(minlen>=3);
		if(bases==null || bases.length<minlen){return null;}
		ArrayList<Orf> orfs=new ArrayList<Orf>();
		if(!ProkObject.callCDS){return orfs;}
//		int mask=63;
		int code=0;
		int start=-2;
		boolean startTruncated=true;
		int frame=0;
		int pos=startFrame;
		
		
		for(; pos<bases.length; pos++){
			byte b=bases[pos];
			int x=AminoAcid.baseToNumber[b];
//			code=((code<<2)|x)&mask;
			code=((code<<2)|x);
			frame++;
			if(frame==3){
				frame=0;
				if(start>=0){
					if(isStopCodon(code, geneticCode) || code<0){//NOTE: This adds a stop codon wherever there are Ns.
						int len=pos-start+1;
						if(len>=minlen){
							Orf f=new Orf(name, start, pos, strand, startFrame, bases, true, CDS, startTruncated, code<0);
							orfs.add(f);
						}
						start=-1;
					}
				}else{
					if(start==-2){
						//At a contig edge, retain a complete first sense codon as a
						//synthetic truncated start.  A stop or ambiguous codon cannot
						//start a CDS; including it produces a leading '*' or X.
						start=(validSyntheticEdgeStart(code, geneticCode) ? pos-2 : -1);
						startTruncated=true;
					}else if(start<0 && isStartCodon(code, geneticCode)){
						start=pos-2;
						startTruncated=false;
					}
				}
				code=0;
			}
		}

		//Close a truncated ORF at its last complete in-frame codon, in either code path.
		if(start>=0){
			//pos is the sequence length; frame counts its unconsumed trailing bases.
			//start is a codon boundary in this frame, so the retained span is whole codons.
			//Removing one more base leaves a two-base codon and shifts stop scoring.
			pos-=frame+1;
			int len=pos-start+1;
			//A bare worker AssertionError can strand ordered output or input queues.
			assert(len>=3 && len%3==0) : KillSwitch.assertDie("Terminal ORF must retain complete codons "+
				"for Orf.getCodon and breakOrf scoring: "
					+start+", "+pos+", frame="+startFrame+", remainder="+frame);
			if(len>=minlen){
				assert(pos<bases.length) : start+", "+pos+", "+bases.length;
				Orf f=new Orf(name, start, pos, strand, startFrame, bases, true, CDS, startTruncated, true);
				orfs.add(f);
			}
		}

//		System.err.print(" *"+orfs.size()+"* ");
		return orfs;
	}
	
	/** 
	 * Generates a list of maximal-length RNAs (non-overlapping).
	 * All RNAs come out in native orientation (unflipped). 
	 * */
	ArrayList<Orf> makeRnasForStrand(String name, byte[] bases, int strand, StatsContainer sc,
			float[] scores, int[] kmersSeen, boolean quitEarly, float bias){
		final int window=sc.lengthAvg;
		if(bases==null || bases.length*2<window){return null;}
		ArrayList<Orf> orfs=new ArrayList<Orf>(sc.type==tRNA ? 32 : 8);
		
		final FrameStats inner=sc.inner;
//		final FrameStats start=sc.start;
//		final FrameStats stop=sc.stop;
		
		final int k=inner.k;
		final int mask=inner.mask;
//		final float invLen=sc.invLengthAvg;
		final int halfWindow=window/2;
		final int maxWindow=(int)(window*1.5f);
		final int maxWindow2=(int)(window*2.5f);
//		final int slop=Tools.max(50, window/8);
		int len=0;
		int kmer=0;
		float currentScore=0;
		float currentScoreAbs=0;
		bias=(bias>-1 ? bias : biases[sc.type]);
		final float maxBias=biases[sc.type]*1.45f;
		
		float thresh=cutoff1[sc.type];
		float prevScore=0;
		int lastStart=0;
		
		float max=0;
		int maxPos=0;
		
		for(int pos=0; pos<bases.length; pos++){
			final byte b=bases[pos];
			assert(b>=0 && b<128) : "Invalid base b="+((int)b)+"; pos="+pos+"\n"+new String(bases)+"\n";
			final int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			
			if(x>=0){
				len++;
				if(len>=k){
					//RNA scoring is frameless: inner.probs[0] only (a single frame, unlike CDS's 3).
					//dif=prob-bias is the per-kmer above/below-average increment; currentScore is a
					//Smith-Waterman-style local accumulator clamped at 0 -- a run of above-bias kmers
					//grows it, a below-bias stretch resets it to 0. That reset is exactly where a
					//candidate region [lastStart,maxPos] is closed in the else-branch below.
					float prob=inner.probs[0][kmer];
					float dif=prob-bias;//Prob above 1 is more likely than average
					currentScoreAbs+=prob;
					currentScore=Tools.max(0, currentScore+dif);
				}
				
				if(currentScore>0){
					if(currentScore>max){
						max=currentScore;
						maxPos=pos;
					}
					if(prevScore<=0){
						lastStart=pos;
					}
				}else{
					int rnaLen=maxPos-lastStart;
					if(max>thresh && rnaLen>=halfWindow){
						if(rnaLen>maxWindow){
							if(bias<=maxBias){
								orfs=null;
								float biasMult=(rnaLen>8*window ? 1.2f : rnaLen>4*window ? 1.1f : 1.05f);
								return makeRnasForStrand(name, bases, strand, sc, scores, kmersSeen, quitEarly, bias*biasMult);
							}
						}
						if(rnaLen<=maxWindow2){
							Orf orf=new Orf(name, lastStart, maxPos, strand, 0, bases, false, sc.type);
							orfs.add(orf);
							orf.orfScore=max;
							if(verbose){
								final int start2=(strand==0 ? lastStart : bases.length-maxPos-1);
								final int stop2=(strand==0 ? maxPos : bases.length-lastStart-1);
								System.err.println("Made Orf "+start2+"\t"+stop2+"\t"+max);
							}
						}
					}
					max=0;
					lastStart=pos;
				}
//				System.err.println("i="+i+"\tscore="+score+"\tmax="+max+"\tmaxPos="+maxPos+"\tprevScore="+prevScore+"\tlastStart="+lastStart);
				prevScore=currentScore;
				
//				if(pos>=223000 && pos<232000){
////					System.err.println("i="+i+"\tscore="+score+"\tmax="+max+"\tmaxPos="+maxPos+"\tprevScore="+prevScore+"\tlastStart="+lastStart);
//					System.out.println(pos+"\t"+currentScore);
//				}
				
			}else{
				len=0;
				kmer=0;
			}
			scores[pos]=currentScoreAbs;
		}
		
//		System.err.println("size="+orfs.size()+", type="+Orf.typeStrings[sc.type]);
		
		
		{
			int rnaLen=maxPos-lastStart;
			if(max>thresh && rnaLen>=halfWindow){
				if(rnaLen>maxWindow){
					if(bias<=maxBias){
						orfs=null;
						float biasMult=(rnaLen>8*window ? 1.2f : rnaLen>4*window ? 1.1f : 1.05f);
						return makeRnasForStrand(name, bases, strand, sc, scores, kmersSeen, quitEarly, bias*biasMult);
					}
				}
				if(rnaLen<=maxWindow2){
					Orf orf=new Orf(name, lastStart, maxPos, strand, 0, bases, false, sc.type);
					orfs.add(orf);
					orf.orfScore=max;
					if(verbose){
						final int start2=(strand==0 ? lastStart : bases.length-maxPos-1);
						final int stop2=(strand==0 ? maxPos : bases.length-lastStart-1);
						System.err.println(start2+"\t"+stop2+"\t"+max);
					}
				}
			}
		}
		
		if(kmersSeen!=null && orfs.size()>0 && sc.kmerSet()!=null){
			fillKmersSeen(bases, kmersSeen, sc.kmerSet(), sc.kLongLen());
		}
		
		float cutoff=cutoff2[sc.type];
		
		for(int i=0; i<orfs.size(); i++){
			Orf orf=orfs.get(i);
//			System.err.println(orf.orfScore);
			boolean good=refineRna(orf, bases, strand, sc, scores, kmersSeen);
			final boolean cm=(cmVerifier!=null && cmVerifier.supports(sc.type));
			final boolean admitted=cm ? good && cmVerifier.verify(orf, bases, cmWindowStart, cmWindowStop, "legacy_rna")
				: admitRnaCandidate(orf.orfScore, good, cutoff, sc.type);
			if(!admitted){
				if(verbose){System.err.println("REJECT: "+orf.toStringFlipped());}
				orfs.set(i, null);
			}else{
				if(verbose){System.err.println("ACCEPT: "+orf.toStringFlipped());}
				if(quitEarly){
					orfs.clear();
					orfs.add(orf);
					return orfs;
				}
			}
		}
		Tools.condenseStrict(orfs);
		
//		assert(false);
		
//		for(int pos=0; pos<bases.length; pos++){
//			final byte b=bases[pos];
//			final int x=AminoAcid.baseToNumber[b];
//			kmer=((kmer<<2)|x)&mask;
//			
//			if(x>=0){
//				len++;
//				if(len>=k){
//					float prob=inner.probs[0][kmer];
//					float dif=prob-1.2f;//Prob above 1 is more likely than average
//					currentScore=Tools.max(0, currentScore+dif);
//					if(currentScore>0){
//						currentStreak++;
//					}else{
//						currentStreak=0;
//					}
//					if(currentScore>200 && currentStreak>1500 && currentStreak<1700){
//						Orf orf=new Orf(name, pos-currentStreak-1, pos, strand, 0, bases, false);
//						orfs.add(orf);
//						orf.orfScore=currentScore;
//						orf.startScore=start.scorePoint(orf.start, bases);
//						orf.stopScore=stop.scorePoint(orf.stop, bases);
//						currentStreak=0;
//						currentScore=0;
//					}
//				}
//			}else{
//				len=0;
//				kmer=0;
//			}
//		}
		
		return orfs;
	}

	boolean admitRnaCandidate(float orfScore, boolean good, float cutoff, int type){
		final boolean scoreFail=orfScore<cutoff, goodFail=!good;
		if(scoreFail){rrnaGateScoreFail[type]++;}
		if(goodFail){rrnaGateGoodFail[type]++;}
		if(scoreFail && goodFail){rrnaGateBothFail[type]++;}
		if(!scoreFail && !goodFail){rrnaGateAdmitted[type]++; return true;}
		return false;
	}
	
	/**
	 * Fills array with cumulative count of observed k-mers from reference set.
	 * Used to validate RNA predictions against known k-mer signatures.
	 *
	 * @param bases Sequence bases to analyze
	 * @param kmersSeen Output array of cumulative k-mer counts
	 * @param set Reference k-mer set for this feature type
	 * @param k K-mer length for rolling hash
	 */
	void fillKmersSeen(byte[] bases, int[] kmersSeen, LongHashSet set, int k){
		final long mask=~((-1L)<<(2*k));
		long kmer=0;
		int len=0;
		int seen=0;
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			final int num=AminoAcid.baseToNumber[b];
			if(num>=0){
				len++;
				kmer=((kmer<<2)|num)&mask;
				if(len>=k && set.contains(kmer)){seen++;}
			}else{
				len=0;
			}
			kmersSeen[i]=seen;
		}
	}
	
	/**
	 * Refines RNA boundaries using start/stop models and alignment validation.
	 * Tests all combinations of potential start and stop positions within slop range.
	 * Validates refined ORF against consensus sequences if available.
	 *
	 * @param orf RNA candidate to refine (modified in place)
	 * @param bases Sequence bases
	 * @param strand Current strand
	 * @param sc Statistics container for scoring
	 * @param scores Cumulative k-mer scores array
	 * @param kmersSeen K-mer observation counts (may be null)
	 * @return true if refinement successful, false if ORF should be rejected
	 */
	boolean refineRna(Orf orf, byte[] bases, int strand, StatsContainer sc, float[] scores, int[] kmersSeen){
		rrnaCandidatesGenerated[sc.type]++;
		if(orf==null){return false;}
		if(verbose){System.err.println("REFINE: "+orf.toStringFlipped());}
		final int window=sc.lengthAvg;
//		final int halfWindow=window/2;
//		final int maxWindow=(int)(window*1.5f);
		final int slop=Tools.max(60, 10+window/8);
		final float invWindow=sc.invLengthAvg;
		final float oldScore=orf.orfScore;
		IntList starts=new IntList();
		IntList stops=new IntList();
		FloatList startScores=new FloatList();
		FloatList stopScores=new FloatList();

		final int leftmost=Tools.max(0, orf.start-slop);
		final int rightmost=Tools.min(bases.length-1, orf.stop+slop);
		if(cmVerifier!=null){cmWindowStart=leftmost;cmWindowStop=rightmost;}
		if(kmersSeen!=null){
			if(kmersSeen[leftmost]>=kmersSeen[rightmost]){
//				System.err.println("Bad: "+oldScore);
				rrnaG5Reject[sc.type]++;
				orf.orfScore=-999;
				return false;
			}else{
//				System.err.println("Good: "+oldScore);
			}
		}
		if(verbose){System.err.println("REFINE2");}
		
		{//starts
			final int left=leftmost;
			final int right=Tools.min(bases.length-1, orf.stop+slop-window);
			final float thresh=cutoff3[sc.type];
			fillPoints(left, right, bases, sc.start, thresh, starts, startScores);
		}
		if(verbose){System.err.println("starts: "+starts.size);}
//		if((orf.start+"").startsWith("146") || true){System.err.println(starts);}
		if(verbose){System.err.println(startScores);}

		{//stops
			final int left=Tools.max(0, orf.start-slop+window);
			final int right=rightmost;
			final float thresh=cutoff4[sc.type];
			fillPoints(left, right, bases, sc.stop, thresh, stops, stopScores);
		}
		if(verbose){System.err.println("stops: "+stops.size);}
//		if((orf.start+"").startsWith("146") || true){System.err.println(stops);}
		if(verbose){System.err.println(stopScores);}

		final float innerThresh=cutoff5[sc.type];
		
		final int minlen=Tools.max(window/2, window-slop);
		final int maxlen=window+slop;
		
		orf.orfScore=0;
		int lastStop=0;
		for(int i=0; i<starts.size; i++){
			final int start=starts.get(i);
			final int startSeen=kmersSeen==null ? 0 : kmersSeen[start];
			final float startScore=startScores.get(i);
//			System.err.println("start="+start);
			for(int j=lastStop; j<stops.size; j++){
				final int stop=stops.get(j);
				final int rnalen=stop-start+1;
//				System.err.println("stop="+stop);
				if(rnalen<minlen){
					lastStop=j+1;
				}else if(rnalen>maxlen){
//					System.err.println("broke");
					break;
				}else if(kmersSeen!=null && kmersSeen[stop]<=startSeen){//TODO: Test this
//					//skip
				}else{
					final int len=stop-start+1;
					final float stopScore=stopScores.get(j);
					final float innerScore=(scores[stop]-scores[start])/len;
//					System.err.println("innerScore="+innerScore);
					assert(rnalen<=maxlen) : "start="+start+", stop="+stop+", rnalen="+rnalen+", minlen="+minlen+", maxlen="+maxlen;
					if(innerScore>=innerThresh){
					//refined score = inner^2 (c) * length-fit (d) * geomean of start/stop
					//point scores (sqrt(a*b)). d is a triangular fit peaking at len==window (==1.0)
					//and falling 2.4*invWindow per base off-length; a,b are floored at 0.25 so
					//a weak/negative start or stop degrades but never zeroes or sign-flips the
					//product. The best (start,stop) pair within [minlen,maxlen] wins.
						final float a=Tools.max(startScore+0.25f, 0.25f);
						final float b=Tools.max(stopScore+0.25f, 0.25f);
						final float c=innerScore*innerScore;
						final float d=(window-(2.4f*Tools.absdif(len, window)))*invWindow;
						final float score=c*d*(float)Math.sqrt(a*b);
						if(verbose && score>=0.2f*orf.orfScore){
							final int start2=(strand==0 ? start : bases.length-stop-1);
							final int stop2=(strand==0 ? stop : bases.length-start-1);
						System.err.println(start2+"-"+stop2+", "+startScore+", "+stopScore+", "+
							innerScore+", "+(score*scoreMult[sc.type])+", "+oldScore);
						}
						if(score>orf.orfScore){
							orf.start=start;
							orf.stop=stop;
							orf.kmerScore=innerScore*len;
//							if(verbose){
//								final int start2=(strand==0 ? start : bases.length-stop-1);
//								final int stop2=(strand==0 ? stop : bases.length-start-1);
//								System.err.println(start2+"-"+stop2+", "+startScore+", "+stopScore+", "+innerScore+", "+(score*scoreMult[sc.type])+", "+oldScore);
//							}
							orf.startScore=startScore;
							orf.stopScore=stopScore;
							orf.orfScore=score;
						}
					}else{
						assert(true);
					}
				}
			}
		}
		orf.orfScore*=scoreMult[sc.type];
		
		if(starts.isEmpty() || stops.isEmpty()){
			if(verbose){System.err.println("No starts or stops.");}
			rrnaG6Reject[sc.type]++;
			orf.orfScore=Tools.min(-999, orf.orfScore-9999);
			return false;
		}
		return refineByAlignment(orf, bases, strand, sc);//Returns true if no consensus is present
	}
	
	/**
	 * Refines ORF boundaries by alignment to consensus sequences.
	 * Attempts alignment with all available consensus sequences for the feature type.
	 *
	 * @param orf ORF to refine
	 * @param bases Sequence bases
	 * @param strand Current strand
	 * @param sc Statistics container with consensus data
	 * @return true if successful alignment found
	 */
	public boolean refineByAlignment(Orf orf, byte[] bases, int strand, StatsContainer sc){
		if(sc.type==r16S){orf.ssuFullModel=null;orf.ssuFullModelName=null;}
		rrnaRefineByAlignmentCalls[sc.type]++;
		if(verbose){System.err.println("ALIGN");}
		Read[] consensus=ProkObject.consensusReads(sc.type);
		if(consensus==null || consensus.length==0){return true;}
		//Counting seam starts here, deliberately -- an empty/no-consensus candidate is
		//explicitly NOT an alignment candidate (Citan's decision, 2026-09-02) and must
		//return before every counter below, including rejected.
		assert(sc.type>=0 && sc.type<rrnaCandidatesEnteringRefinement.length) : sc.type;
		rrnaCandidatesEnteringRefinement[sc.type]++;
		boolean refined=false;
		int acceptedIndex=-1;
		int attemptsMade=0;
//		System.err.println("Initial: "+orf.start+", "+orf.stop);
		for(int i=0; i<consensus.length; i++){
			Read r=consensus[i];
			final int start0=orf.start, stop0=orf.stop;
			final float score0=orf.orfScore;
//			refined=refineByAlignment(orf, bases, strand, sc, r.bases, 15, 15, 2);
			refined=refineByAlignment(orf, bases, strand, sc, r.bases, sc.startSlop(), sc.stopSlop(), 2);
			if(!refined){
				orf.start=start0;
				orf.stop=stop0;
				orf.orfScore=score0;
			}
			attemptsMade++;
			if(i==0){rrnaUniversalAttempts[sc.type]++;}else{rrnaFallbackAttempts[sc.type]++;}
			if(refined){
				acceptedIndex=i;
				if(sc.type==r16S){
					assert(orf.length()>0 && Float.isFinite(lastAttemptIdentity) &&
						lastAttemptIdentity>=0f && lastAttemptIdentity<=1f)
						: KillSwitch.assertDie("Successful16S refinement must retain its accepted identity "+
							"for length-times-identity-squared scoring; identity="+lastAttemptIdentity+
							" length="+orf.length());
					// Match NcrnaScavenger's SSU score scale after actual endpoint refinement.
					// Candidate generation and the existing identity acceptance remain unchanged.
					orf.orfScore=orf.length()*lastAttemptIdentity*lastAttemptIdentity;
					orf.ssuFullModel=r.bases;
					orf.ssuFullModelName=r.id;
				}
			}
			if(attemptSink!=null){
				assert(!"UNKNOWN".equals(lastAttemptReason));
				// The strand-1 pass operates on an in-place reverse complement.  Emit
				// forward-reference coordinates so diagnostic consumers can join them to
				// the original input sequence without duplicating caller internals.
				final int preStart=(strand==0 ? start0 : bases.length-stop0-1);
				final int preStop=(strand==0 ? stop0 : bases.length-start0-1);
				final int postStart=(strand==0 ? orf.start : bases.length-orf.stop-1);
				final int postStop=(strand==0 ? orf.stop : bases.length-orf.start-1);
				final String candidateId=System.identityHashCode(bases)+"_"+System.identityHashCode(this)+
					"_"+strand+"_"+preStart+"_"+preStop;
				attemptSink.onAttempt(i, consensus.length, preStart, preStop, score0, refined, postStart, postStop, strand,
						orf.scafName, r.id, r.bases.length, lastAttemptIdentity, lastAttemptReason, candidateId);
			}
			if(refined || !rrnaFallback){break;}
		}
		rrnaAttemptCountHist[sc.type]=growHistRow(rrnaAttemptCountHist[sc.type], attemptsMade);
		rrnaAttemptCountHist[sc.type][attemptsMade]++;
		if(refined){
			if(acceptedIndex>0){rrnaFallbackOnlyRescues[sc.type]++;}
			if(verbose){System.err.println("Aligned to: "+orf.start+", "+orf.stop);}
		}else{
			rrnaRejectedCandidates[sc.type]++;
			if(verbose){System.err.println("Alignment failed.");}
			orf.orfScore=Tools.min(-999, orf.orfScore-9999);
		}
		return refined;
	}
	
	/**
	 * Refines ORF by alignment to specific consensus sequence.
	 * Dispatches to either IDAligner or SingleStateAligner based on configuration.
	 *
	 * @param orf ORF to refine
	 * @param bases Sequence bases
	 * @param strand Current strand
	 * @param sc Statistics container
	 * @param consensus Reference sequence for alignment
	 * @param startSlop Maximum start position adjustment allowed
	 * @param stopSlop Maximum stop position adjustment allowed
	 * @param recurLimit Recursion limit for boundary expansion
	 * @return true if successful alignment with sufficient identity
	 */
	boolean refineByAlignment(Orf orf, byte[] bases, int strand, StatsContainer sc, byte[] consensus, 
			final int startSlop, final int stopSlop, int recurLimit){
		lastAttemptIdentity=-1f;
		if(attemptSink!=null){lastAttemptReason="UNKNOWN";}
		if(useIDAligner) {
			return refineByAlignment_IDA(orf, bases, strand, 
					sc, consensus, startSlop, stopSlop, recurLimit);
		}else {
			return refineByAlignment_SSA(orf, bases, strand, 
					sc, consensus, startSlop, stopSlop, recurLimit);
		}
	}
	
	/**
	 * Refines ORF boundaries using IDAligner for fast identity-based alignment.
	 * Pads ORF region and aligns consensus to find best boundaries within identity threshold.
	 *
	 * @param orf ORF to refine (modified in place)
	 * @param bases Sequence bases
	 * @param strand Current strand
	 * @param sc Statistics container with identity threshold
	 * @param consensus Reference sequence
	 * @param startSlop Maximum start adjustment
	 * @param stopSlop Maximum stop adjustment
	 * @param recurLimit Unused for IDAligner
	 * @return true if alignment meets identity threshold
	 */
	boolean refineByAlignment_IDA(Orf orf, byte[] bases, int strand, StatsContainer sc, byte[] consensus, 
			final int startSlop, final int stopSlop, int recurLimit){
		final int start0=orf.start;
		final int stop0=orf.stop;
		
		assert(start0>=0 && start0<bases.length) : start0+", "+stop0;
		assert(stop0>=0 && stop0<bases.length) : start0+", "+stop0;
		
		final float minID=sc.minIdentity();
		final int padding=Tools.min(alignmentPaddingIDA, 30+sc.lengthAvg/4);
		final int a=Tools.max(0, orf.start-padding);
		final int b=Tools.min(bases.length-1, orf.stop+padding);
		final int reflen=b-a+1;
		assert(reflen>0) : reflen;
		if(reflen>10*sc.lengthAvg && reflen>20000){
			System.err.println("Skipped reflen "+reflen+"/"+sc.lengthAvg+" for "
					+ "seqlen="+bases.length+", orf="+orf.toString());
			assert(false);
			//TODO: Possibly change return to -1, 0, 1 ("can't align")
			//Should be a limit on window size...
			//Also consider shrinking matrix after jumbo alignments
			if(attemptSink!=null){lastAttemptIdentity=-1f; lastAttemptReason="REFLEN_TOO_LARGE";}
			return false;
		}
		assert(a>=0 && b<bases.length) : a+", "+b;
		IDAligner ida=idaligner.Factory.makeIDAligner();
		int[] pos=new int[2];
		float id=ida.align(consensus, bases, pos, a, b);
		if(id<minID){
			if(attemptSink!=null){lastAttemptIdentity=id; lastAttemptReason="IDENTITY_BELOW_THRESHOLD";}
			return false;
		}
		if(cmVerifier!=null){cmWindowStart=a;cmWindowStop=b;}
		
		final int rstart=Tools.max(pos[0], 0);
		final int rstop=Tools.min(pos[1], bases.length-1);
		
		if(Tools.absdif(rstart, start0)>startSlop){orf.start=rstart;}
		if(Tools.absdif(rstop, stop0)>stopSlop){orf.stop=rstop;}
		assert(orf.length()>0) : "\nreflen="+reflen+", orflen="+orf.length()+
			", id"+id+", minid="+minID+", strand="+orf.strand+"\n"+
			"qlen="+consensus.length+", rlen="+bases.length+", a="+a+", b="+b+"\n"+
			"pos="+Arrays.toString(pos)+", start0="+start0+", stop0="+stop0+"\n"+orf;
			//+new String(consensus)+"\n"+new String(bases, a, b-a+1)+"\n";
		lastAttemptIdentity=id;
		if(attemptSink!=null){lastAttemptReason=(orf.length()>0 ? "ACCEPTED" : "EMPTY_RESULT");}
		return orf.length()>0;
	}
	
	/**
	 * Refines ORF boundaries using SingleStateAligner for detailed alignment.
	 * Performs dynamic programming alignment with optional boundary expansion.
	 * Can recursively expand search region if initial alignment extends beyond padding.
	 *
	 * @param orf ORF to refine (modified in place)
	 * @param bases Sequence bases
	 * @param strand Current strand
	 * @param sc Statistics container with identity threshold
	 * @param consensus Reference sequence
	 * @param startSlop Maximum start adjustment
	 * @param stopSlop Maximum stop adjustment
	 * @param recurLimit Recursion limit for boundary expansion
	 * @return true if alignment meets identity threshold
	 */
	boolean refineByAlignment_SSA(Orf orf, byte[] bases, int strand, StatsContainer sc, byte[] consensus, 
			final int startSlop, final int stopSlop, int recurLimit){
		final int start0=orf.start;
		final int stop0=orf.stop;
		
		assert(start0>=0 && start0<bases.length) : start0+", "+stop0;
		assert(stop0>=0 && stop0<bases.length) : start0+", "+stop0;
		
		final float minID=sc.minIdentity();
		final int padding=Tools.min(alignmentPadding, 30+sc.lengthAvg/4);
		final int a=Tools.max(0, orf.start-padding);
		final int b=Tools.min(bases.length-1, orf.stop+padding);
		final int reflen=b-a+1;
		if(reflen>10*sc.lengthAvg && reflen>20000){
			System.err.println("Skipped reflen "+reflen+"/"+sc.lengthAvg+" for "
					+ "seqlen="+bases.length+", orf="+orf.toString());
			assert(false);
			//TODO: Possibly change return to -1, 0, 1 ("can't align")
			//Should be a limit on window size...
			//Also consider shrinking matrix after jumbo alignments
			if(attemptSink!=null){lastAttemptIdentity=-1f; lastAttemptReason="REFLEN_TOO_LARGE";}
			return false;
		}
		assert(a>=0 && b<bases.length) : a+", "+b;
		SingleStateAlignerFlat2 ssa=getSSA();
		final int minScore=ssa.minScoreByIdentity(consensus.length, minID);
		int[] max=ssa.fillUnlimited(consensus, bases, a, b, minScore);
		if(max==null){
			if(attemptSink!=null){lastAttemptIdentity=-1f; lastAttemptReason="NO_ALIGNMENT_FOUND";}
			return false;
		}
		
		final int rows=max[0];
		final int maxCol=max[1];
		final int maxState=max[2];
//		final int maxScore=max[3];
//		final int maxStart=max[4];
		
		//returns {score, bestRefStart, bestRefStop} 
		//padded: {score, bestRefStart, bestRefStop, padLeft, padRight};
		final int[] score=ssa.score(consensus, bases, a, b, rows, maxCol, maxState);
		final int rstart=Tools.max(score[1], 0);
		final int rstop=Tools.min(score[2], bases.length-1);
		final float id=ssa.tracebackIdentity(consensus, bases, a, b, rows, maxCol, maxState, null);
		
		assert(rstart>=0 && rstart<bases.length) : "bases="+bases.length+
			", rstart="+rstart+", rstop="+rstop+", a="+a+", b="+b+"\n"+"score="+Arrays.toString(score)+", id="+id;
		assert(rstop>=0 && rstop<bases.length) : rstart+", "+rstop;
		
		if(score.length>3 && recurLimit>0){
			final int padLeft=score[3];
			final int padRight=score[4];
			//TODO:  This takes extra memory.  It may be better to cap the width or ignore start0/stop0 here.
			int rstart2=Tools.max(0, Tools.min(start0, rstart)-10);
			int rstop2=Tools.min(bases.length-1, Tools.max(stop0, rstop)+10);
			assert(rstart2>=0 && rstart2<bases.length) : rstart2+", "+rstop2;
			assert(rstop2>=0 && rstop2<bases.length) : rstart2+", "+rstop2;
			orf.start=rstart2;
			orf.stop=rstop2;
			if(orf.start<a || orf.stop>b){
				boolean ret=refineByAlignment(orf, bases, strand, sc, consensus, 0, 0, recurLimit-1);
				if(ret){return true;}
			}
			orf.start=start0;
			orf.stop=stop0;
//			assert(false) : "Don't traceback after recur.";
		}
		if(verbose){
			System.err.println("Identity: "+id);
		}
//		assert(score.length==3) : "TODO: Handle padding requests.";
		
//		System.err.println("Identity: "+String.format("%.2f", 100*id)+"; location: "+rstart+"-"+rstop);
		if(id<minID){
			if(attemptSink!=null){lastAttemptIdentity=id; lastAttemptReason="IDENTITY_BELOW_THRESHOLD";}
			return false;
		}
		
		
		if(cmVerifier!=null){cmWindowStart=a;cmWindowStop=b;}
		if(Tools.absdif(rstart, start0)>startSlop){orf.start=rstart;}
		if(Tools.absdif(rstop, stop0)>stopSlop){orf.stop=rstop;}
		lastAttemptIdentity=id;
		if(attemptSink!=null){lastAttemptReason="ACCEPTED";}
		return true;
	}
	
	/**
	 * Finds high-scoring positions within a sequence region using iterative threshold reduction.
	 * Continues lowering threshold until at least 8 points found or minimum threshold reached.
	 *
	 * @param left Start position of search region
	 * @param right End position of search region
	 * @param bases Sequence bases
	 * @param fs Frame statistics for position scoring
	 * @param thresh Initial score threshold
	 * @param points Output list of qualifying positions
	 * @param scores Output list of corresponding scores
	 */
	void fillPoints(final int left, final int right, final byte[] bases, final FrameStats fs, float thresh,
			final IntList points, final FloatList scores){
		//TODO: Confirmed bug CG01 - nonpositive thresholds can loop forever with fewer
		//than 8 qualifying positions; trnastop=0 reproduces through the legacy tRNA CLI.
		//Keep its bounded-loop repair separate from genetic-code integration.
		points.clear();
		scores.clear();
		final float minThresh=thresh;//thresh*0.05f;
//		System.err.println(left+", "+right+", "+thresh+", "+minThresh);
		while(points.size()<8 && thresh>=minThresh){
//			System.err.println("Running with thresh="+thresh);
			points.clear();
			scores.clear();
			for(int i=left; i<right; i++){
				float score=fs.scorePoint(i, bases);
//				System.err.println(i+", "+score);
				if(score>=thresh){
					points.add(i);
					scores.add(score);
				}
			}
			thresh=thresh*0.75f;
		}
	}
	
	/**
	 * Generate all possible genes from each Orf, and return them in a new set of lists.
	 * @param frameLists
	 * @param bases
	 * @return Lists of orfs.
	 */
	private ArrayList<Orf>[] breakOrfs(ArrayList<Orf>[] frameLists, byte[] bases){

		@SuppressWarnings("unchecked")
		ArrayList<Orf>[] brokenLists=new ArrayList[6];
		for(int strand=0; strand<2; strand++){
			for(int frame=0; frame<3; frame++){
				int fnum=frame+3*strand;
				ArrayList<Orf> longest=frameLists[fnum]; //Longest Orf per stop
				// All high scoring Orfs, including multiple per stop, with different starts.
				ArrayList<Orf> broken=new ArrayList<Orf>();
				if(longest!=null) {
					final Object[] evidence=prepareProteinStops(longest, bases);
					for(int i=0; i<longest.size(); i++){
						final Orf orf=longest.get(i);
						assert(orf.frame==frame);
						assert(orf.strand==strand);
						ArrayList<Orf> temp=breakOrf(orf, bases, evidence==null ? null : evidence[i]);
						if(temp!=null){
							broken.addAll(temp);
						}
					}
				}
				Collections.sort(broken);
				brokenLists[fnum]=broken;
			}
			//Reverse-complement bases after processing each strand
			Vector.reverseComplementInPlaceFast(bases);
		}
		return brokenLists;
	}
	
	/**
	 * Generate an Orf for each possible start codon.
	 * Retain only the high-scoring ones.
	 * @param longest Longest open reading frame for a given stop.
	 * @param bases Bases, oriented for this Orf.
	 * @return List of Orfs.
	 */
	private ArrayList<Orf> breakOrf(Orf longest, byte[] bases, Object stopEvidence){
		assert(longest.start<longest.stop);
		final int flipped=longest.flipped();
		if(flipped==1){longest.flip();}//Now the orf is aligned to its native strand
		
		geneStopsMade++;
		
		final FrameStats innerStats=pgm.statsCDS.inner;
		final FrameStats startStats=pgm.statsCDS.start;
		final FrameStats stopStats=pgm.statsCDS.stop;
		//TODO: Confirmed compatibility bug CG02 - an expanded start mask can admit
		//codons whose negatives were only noise-sampled by GeneModel.findStartCodons
		//when this PGM was trained. In the shipped model, central CTG enrichment is
		//27.93 versus ATG 2.51 (FrameStats.calculate); this preferentially selects
		//wrong starts. Training and inference must use compatible candidate sets.
		
		final String name=longest.scafName;
		final int start=longest.start;
		final int stop=longest.stop;
		final int strand=longest.strand;
		final int max=Tools.min(longest.stop-2, longest.stop-minLen+4);

		assert(pgm!=null) : pgm;
		assert(pgm.statsCDS!=null) : pgm;
		assert(pgm.statsCDS.inner!=null) : pgm.statsCDS;
		assert(pgm.statsCDS.inner.k>0) : pgm.statsCDS.inner;
		
		final int k=innerStats.k;
		final int mask=~((-1)<<(2*k));

		final float stopScore=stopStats.scorePoint(longest.stop, bases);
		stCds.geneStopScoreSum+=stopScore;
		stCds.geneStopScoreCount++;
		final ProteinSupportFate.Entry fate=proteinFates!=null &&
				((ProteinSupport.StopEvidence)stopEvidence).supported()
			? new ProteinSupportFate.Entry(longest, (ProteinSupport.StopEvidence)stopEvidence,
				stopScore, minLen, minStartScore) : null;
		if(fate!=null){proteinFates.add(fate);}
		
		ArrayList<Orf> broken=new ArrayList<Orf>();
		int created=0;
		
		int codon=0;
		int kmer=0;
		int len=0;
		int numKmers=0;
		float currentScore=0;
		for(int pos=start, currentFrame=0; pos<=stop; pos++){
			final byte b=bases[pos];
			final int x=AminoAcid.baseToNumber[b];
			codon=((codon<<2)|x);
			kmer=((kmer<<2)|x)&mask;
			
			if(x>=0){
				len++;
				if(len>=k){
					float prob=innerStats.probs[currentFrame][kmer];
					float dif=prob-0.99f;//Prob above 1 is more likely than average
					currentScore+=dif;
//					outstream.println("pos="+pos+"\tdif="+Tools.format("%.4f", dif)+",\tscore="+Tools.format("%.4f", currentScore)+
//							"\tasStart="+Tools.format("%.4f", pgm.calcStartScore(pos-2, bases))+"\tasStop="+Tools.format("%.4f", stopStats.scorePoint(pos, bases))+
//							"\tcodon="+AminoAcid.kmerToString(kmer, 3)+" frame="+(currentFrame));
				}else{
//					outstream.println("pos="+pos+"\tdif="+Tools.format("%.4f", 0.0)+",\tscore="+Tools.format("%.4f", 0.0)+
//							"\tasStart="+Tools.format("%.4f", pgm.calcStartScore(pos-2, bases))+"\tasStop="+Tools.format("%.4f", stopStats.scorePoint(pos, bases))+
//							"\tcodon="+AminoAcid.kmerToString(kmer, 3)+" frame="+(currentFrame));
				}
			}else{
				len=0;
				kmer=0;
			}
			
			currentFrame++;
//			outstream.println("pos="+pos+", codon="+AminoAcid.kmerToString(kmer, 3)+", frame="+currentFrame+", start="+start+", isStartCodon="+pgm.isStartCodon(codon));
			if(currentFrame>2){
				currentFrame=0;
				final boolean syntheticEdgeStart=(pos==start+2 && validSyntheticEdgeStart(codon, geneticCode));
				if(fate!=null && (syntheticEdgeStart || isStartCodon(codon, geneticCode))){
					if(pos>=max){fate.rejectedLength++;}
					else if(created>=breakLimit){fate.rejectedLimit++;}
				}
				if(pos<max && created<breakLimit && (syntheticEdgeStart || isStartCodon(codon, geneticCode))){
//					outstream.println(x);
					int glen=stop-pos+3;
					assert(glen>=minLen) : "glen="+glen+", minLen="+minLen+", pos="+pos+", max="+max+", start="+start;
					
					int oStart=pos-2;
					float startScore=startStats.scorePoint(oStart, bases);

					stCds.geneStartScoreSum+=startScore;
					stCds.geneStartScoreCount++;
					
					stCds.lengthSum+=(stop-(pos-2)+1);
					stCds.lengthCount++;
					
					final boolean nativeEligible=created<breakLimit && (startScore>=minStartScore || pos<6);
					if(fate!=null && !nativeEligible){fate.rejectedStartScore++;}
					if(nativeEligible){
						Orf orf=new Orf(name, pos-2, stop, strand, longest.frame, bases, false, longest.type,
							pos-2==start && longest.startTruncated, longest.stopTruncated);
						
						geneStartsMade++;
						orf.kmerScore=currentScore;
						orf.startScore=startScore;
						orf.stopScore=stopScore;

						assert(orf.frame==longest.frame);
						assert(orf.strand==strand);

						if(strand==1){orf.flip();}
						broken.add(orf);
						if(nativeEligible){created++;}
					}
				}
				codon=0;
			}
		}
		
		final int size=broken.size();
		if(size<1){
			if(fate!=null){fate.retained(broken);}
			if(supported(stopEvidence)){withoutStarts(stopEvidence);}
			return broken;
		}
		Orf best=null;
		Orf bestStart=null;
		for(int i=0; i<size; i++){
			Orf orf=broken.get(i);
			
			//This fixes scores because they were generated together, from start to stop, to make this O(N) instead of O(N^2).
			orf.kmerScore=currentScore-orf.kmerScore;
			if(dumpOrfScores){
				outstream.println("ORFDUMP\t"+orf.scafName+"\t"+(orf.start+1)+"\t"+(orf.stop+1)+
					"\t"+orf.strand+"\t"+orf.startScore+"\t"+orf.stopScore+"\t"+
					orf.averageKmerScore()+"\t"+orf.length());
			}
			final boolean nearEdge=(orf.start<NN_EDGE_MARGIN || orf.stop+NN_EDGE_MARGIN>=bases.length);
			if(!nearEdge && activeMetaNet!=null && activeNet!=null && activeStopNet!=null){
				orf.orfScore=orf.calcOrfScoreHybridMeta(activeNet, activeStopNet, activeMetaNet, contigGC, bases, longest.length());
			}else if(!nearEdge && activeNet!=null && activeStopNet!=null){
				orf.orfScore=orf.calcOrfScoreHybridBoth(activeNet, activeStopNet, contigGC, bases);
			}else if(!nearEdge && activeNet!=null){
				orf.orfScore=orf.calcOrfScoreHybridStart(activeNet, contigGC, bases);
			}else if(!nearEdge && activeStopNet!=null){
				orf.orfScore=orf.calcOrfScoreHybridStop(activeStopNet, contigGC, bases);
			}else{
				orf.orfScore=orf.calcOrfScore();
			}
			if(MINOR_START_MULTIPLIER!=1f){
				orf.orfScore=minorStartScore(orf.orfScore, orf.startCodon, orf.startTruncated, MINOR_START_MULTIPLIER);
			}
			annotateProtein(orf, stopEvidence);
			if(fate!=null){fate.admitted(orf);}
			if(best==null || orf.orfScore>=best.orfScore){best=orf;}
			if(bestStart==null || orf.startScore>=bestStart.startScore){bestStart=orf;}
			
			stCds.geneInnerScoreSum+=orf.averageKmerScore();
			stCds.geneInnerScoreCount++;
		}
		if(broken.isEmpty()){return broken;}
		final int retained=broken.size(), sizeCutoff=Tools.max(5, retained/2);
		if(tailCascade!=null){
			// Snapshot/prune native scores before adding any evidence. The first
			// sorted candidate wins score ties by the existing positional comparator.
			Collections.sort(broken, PFeature.featureComparatorScore);
			final Orf nativeWinner=broken.get(0);
			final prot.FamilyTailCascade.Assignment assignment=(prot.FamilyTailCascade.Assignment)stopEvidence;
			final boolean targeted=targetedTailStop>=0 && assignment!=null && assignment.assigned() &&
				!nativeWinner.stopTruncated &&
				StopTailCandidateScan.stopKey(tailStopLeft(nativeWinner),nativeWinner.strand)==targetedTailStop;
			prot.FamilyMatchInformation.View formulaView=null;
			if(tailCascadeFormula!=null && tailCascadeFormula.requiresFamilySetView() && assignment!=null &&
					assignment.assigned()){
				final int ordinal=tailCascadeEvidence.ordinalAtStop(tailStopLeft(nativeWinner),nativeWinner.strand);
				if(ordinal<0 || tailCascadeEvidence.assignment(ordinal)!=assignment){
					throw new IllegalStateException("Formula view must come from this candidate's exact accepted stop evidence");
				}
				formulaView=tailCascadeEvidence.primaryView(ordinal);
			}
			final float appliedBonus=assignment==null || tailCascadeFormula!=null ? 0 :
				tailCascadeInformationMode ? assignment.informationBonus(tailCascadeBonus) :
				assignment.stopBonus(targeted ? targetedTailBonus : tailCascadeBonus,tailCascadeMarginMode);
			for(int i=0; i<retained; i++){
				final Orf orf=broken.get(i);
				if(!Float.isFinite(orf.orfScore)){
					throw new IllegalStateException("Nonfinite native ORF score at "+
						orf.scafName+":"+orf.start);
				}
				orf.tailCascadeNativeScore=orf.orfScore;
				if(assignment!=null && assignment.assigned() && !(keepAtLeastOneOrf && i==0) &&
					softPruned(orf, i, sizeCutoff, best, bestStart)){
					tailCascadeStats.nativePruned++;
				}
			}
			for(Orf orf:broken){
				orf.tailCascadeAssignment=assignment;
				orf.tailCascadePin=tailCascadePin;
				if(assignment!=null && assignment.assigned()){
					orf.proteinBonus=tailCascadeFormula==null ? appliedBonus :
						tailCascadeFormula.score(assignment,formulaView,candidateAminoAcidLength(orf));
					if(appliedBonus!=0){orf.orfScore+=appliedBonus;}
					else if(orf.proteinBonus!=0){orf.orfScore+=orf.proteinBonus;}
					if(!Float.isFinite(orf.orfScore)){throw new IllegalArgumentException("Tail bonus overflows the native ORF score");}
					orf.tailCascadeProtected=tailCascadeProtect && orf==nativeWinner;
				}
				if(tailCascadeObserver!=null){observeTailCandidate(orf, false);}
			}
		}
		
		//Sort to by score descending to eliminate low-scoring copies. 
		if(broken.size()>1){Collections.sort(broken, PFeature.featureComparatorScore);}
		if(tailCascade!=null && broken.get(0).orfScore>best.orfScore){best=broken.get(0);}
		
		int removed=0;
		for(int i=0; i<retained; i++){
			Orf orf=broken.get(i);
			if(tailCascade!=null && orf.tailCascadeAssignment!=null && orf.tailCascadeAssignment.assigned()
					&& !(keepAtLeastOneOrf && i==0) && softPruned(orf, i, sizeCutoff, best, bestStart)){
				tailCascadeStats.afterBonusPruned++;
				if(orf.tailCascadeProtected){tailCascadeStats.guaranteeRetained++;}
			}
			if(protectsFromPruning(orf) || (keepAtLeastOneOrf && i==0)) {
				//skip
			}else if(orf.averageKmerScore()<minInnerScore || orf.orfScore<minOrfScore || 
					orf.orfScore/orf.length()<minAvgScore || 
					orf.orfScore<0.5f*best.orfScore-10 ||
					(orf.startScore<bestStart.startScore-0.55f && orf.kmerScore<best.kmerScore*1.1f &&
					orf!=best)){
				if(fate!=null){
					final String reason=orf.averageKmerScore()<minInnerScore ? "middle_score"
						: orf.orfScore<minOrfScore ? "orf_score"
						: orf.orfScore/orf.length()<minAvgScore ? "average_score"
						: orf.orfScore<0.5f*best.orfScore-10 ? "relative_orf_score" : "relative_start_score";
					fate.rejected(orf, reason);
				}
				broken.set(i, null);
				removed++;
			}else if(i>sizeCutoff){
				if(fate!=null){fate.rejected(orf, "rank_limit");}
				broken.set(i, null);
				removed++;
			}
		}
		if(removed>0){
			Tools.condenseStrict(broken);
		}
		if(tailCascade!=null && supported(stopEvidence) && !broken.isEmpty()){
			tailCascadeStats.retainedStops++;
		}
		if(tailCascadeObserver!=null && tailCascade!=null){
			for(Orf orf:broken){observeTailCandidate(orf, true);}
		}
		if(fate!=null){fate.retained(broken);}
		
		geneStartsRetained+=broken.size();
		geneStopsRetained+=(broken.size()>0 ? 1 : 0);
		
		if(flipped==1){longest.flip();}
		return broken;
	}

	private Object[] prepareProteinStops(final ArrayList<Orf> longest, final byte[] bases){
		if(tailCascade!=null){
			final Object[] evidence=new Object[longest.size()];
			for(int i=0; i<longest.size(); i++){
				final Orf orf=longest.get(i);
				if(orf.stopTruncated){continue;}
				final int ordinal=tailCascadeEvidence.ordinalAtStop(tailStopLeft(orf), orf.strand);
				if(ordinal<0){continue;}
				if(tailCascadeParents[ordinal]){
					throw new IllegalStateException("Duplicate native parent for one biological stop");
				}
				tailCascadeParents[ordinal]=true;
				if(tailCascadeObserver!=null){tailCascadeObserver.parent(ordinal);}
				evidence[i]=tailCascadeEvidence.assignment(ordinal);
			}
			return evidence;
		}
		return proteinSupport!=null ? proteinSupport.prepareStops(longest, bases, geneticCode)
			: stopTailSupport!=null ? stopTailSupport.prepareStops(longest, bases, geneticCode) : null;
	}
	private void annotateProtein(final Orf orf, final Object evidence){
		if(proteinSupport!=null){proteinSupport.annotate(orf, (ProteinSupport.StopEvidence)evidence);}
		else if(stopTailSupport!=null){stopTailSupport.annotate(orf, (StopTailSupport.StopEvidence)evidence);}
	}
	private static boolean supported(final Object evidence){
		if(evidence instanceof prot.FamilyTailCascade.Assignment){
			return ((prot.FamilyTailCascade.Assignment)evidence).assigned();
		}
		return evidence instanceof ProteinSupport.StopEvidence ? ((ProteinSupport.StopEvidence)evidence).supported()
			: evidence instanceof StopTailSupport.StopEvidence && ((StopTailSupport.StopEvidence)evidence).supported();
	}
	private void withoutStarts(final Object evidence){
		if(evidence instanceof prot.FamilyTailCascade.Assignment){tailCascadeStats.noStarts++;}
		if(evidence instanceof ProteinSupport.StopEvidence){proteinSupport.withoutStarts++;}
		else if(evidence instanceof StopTailSupport.StopEvidence){stopTailSupport.withoutStarts++;}
	}
	public static boolean protectsFromPruning(final Orf orf){
		if(orf.tailCascadeProtected){
			assert(orf.tailCascadeAssignment!=null && orf.tailCascadeAssignment.assigned()) :
				"Only an accepted biological stop may protect its native winner";
			return true;
		}
		// Cascade protection is explicit; the legacy no-metadata rule treats any
		// positive proteinBonus as protection and would bypass tailprotect=false.
		if(orf.tailCascadeAssignment!=null){return false;}
		return orf.stopTailEvidence!=null ? StopTailSupport.protectsFromPruning(orf) :
			ProteinSupport.protectsFromPruning(orf);
	}

	/** The existing soft filters, evaluated separately from homology protection. */
	private boolean softPruned(final Orf orf, final int rank, final int sizeCutoff, final Orf best,
			final Orf bestStart){
		assert(best!=null && bestStart!=null && rank>=0) :
			"Pruning comparisons require a nonempty native candidate group";
		return orf.averageKmerScore()<minInnerScore || orf.orfScore<minOrfScore ||
			orf.orfScore/orf.length()<minAvgScore || orf.orfScore<0.5f*best.orfScore-10 ||
			(orf.startScore<bestStart.startScore-0.55f && orf.kmerScore<best.kmerScore*1.1f && orf!=best)
			|| rank>sizeCutoff;
	}

	/** Original genomic left coordinate of the biological stop codon in either stored orientation.
	 * makeOrfs flips minus parents again after their constructor flip; broken candidates
	 * instead retain their constructor-relative flip. Strand alone cannot identify the frame. */
	public static int tailStopLeft(final Orf orf){
		assert(orf!=null && orf.type==CDS && !orf.stopTruncated) :
			"Only a real CDS terminal codon can identify classifier evidence";
		return orf.strand==0 ? (orf.flipped()==0 ? orf.stop-2 : orf.scaflen-orf.start-3)
			: (orf.flipped()==0 ? orf.scaflen-orf.stop-1 : orf.start);
	}

	/** Translated CDS length for the formula endpoint, excluding the terminal stop
	 * codon when the candidate has one. Legacy null-code terminal clipping can
	 * leave one or two extra edge bases; only complete truncated codons translate. */
	static int candidateAminoAcidLength(final Orf orf){
		final int codingBases=orf.stopTruncated ? orf.length()/3*3 : orf.length()-3;
		if(codingBases<=0 || codingBases%3!=0){
			throw new IllegalArgumentException("Candidate ORF length must encode whole translated codons: length="+
				orf.length()+" stopTruncated="+orf.stopTruncated);
		}
		return codingBases/3;
	}

	/** Snapshot only; broken CDS candidates already have genomic coordinates on both strands. */
	private void observeTailCandidate(final Orf orf, final boolean retained){
		assert(orf.type==CDS && orf.flipped()==orf.strand && tailCascadeEvidence!=null)
			: "breakOrf constructs genomic CDS candidates before the tail observation hook";
		final boolean realStop=!orf.stopTruncated;
		final int ordinal=realStop ? tailCascadeEvidence.ordinalAtStop(tailStopLeft(orf), orf.strand) : -1;
		tailCascadeObserver.candidate(orf.start, orf.stop, orf.strand, realStop, ordinal,
			candidateAminoAcidLength(orf), orf.tailCascadeNativeScore, orf.proteinBonus, orf.orfScore, retained);
	}

	/** Count unique selected biological stops, then release this contig's immutable results. */
	private void finishTailCascade(final ArrayList<Orf> path){
		assert(tailCascadeEvidence!=null && tailCascadeParents!=null) :
			"Every classified contig must finish its own evidence map";
		for(int i=0; i<tailCascadeParents.length; i++){
			if(!tailCascadeParents[i] && tailCascadeEvidence.assignment(i).assigned()){tailCascadeStats.noParent++;}
		}
		final map.LongHashSet selected=new map.LongHashSet(16);
		for(Orf orf:path){
			if(orf.tailCascadeAssignment!=null && orf.tailCascadeAssignment.assigned()){
				if(selected.add(StopTailCandidateScan.stopKey(tailStopLeft(orf), orf.strand))){
					tailCascadeStats.selected++;
					tailCascadeStats.selectedStages[TailCascadeStats.stage(orf.tailCascadeAssignment)]++;
				}
			}
		}
		if(tailCascadeObserver!=null){tailCascadeObserver.end();}
		tailCascadeEvidence=null; tailCascadeParents=null;
	}

	/** Per-caller totals; merge only after its processing thread has joined. */
	public static final class TailCascadeStats{
		public void add(final TailCascadeStats other){
			assert(other!=null && other!=this) : "Only completed independent caller totals may be merged";
			searched+=other.searched; reused+=other.reused; accepted+=other.accepted;
			under100+=other.under100; from100To199+=other.from100To199;
			noParent+=other.noParent; noStarts+=other.noStarts; nativePruned+=other.nativePruned;
			afterBonusPruned+=other.afterBonusPruned; guaranteeRetained+=other.guaranteeRetained; selected+=other.selected;
			retainedStops+=other.retainedStops;
			for(int i=0; i<STAGES.length; i++){
				acceptedStages[i]+=other.acceptedStages[i];
				selectedStages[i]+=other.selectedStages[i];
			}
		}
		private static int stage(final prot.FamilyTailCascade.Assignment assignment){
			for(int i=0; i<STAGES.length; i++){if(STAGES[i].equals(assignment.terminal)){return i;}}
			throw new IllegalStateException("Unknown accepted tail-classifier stage: "+assignment.terminal);
		}
		public static int stageCount(){return STAGES.length;}
		public static String stageName(final int index){return STAGES[index];}
		/** Logical evidence coverage; reused is the subset not classified again in this arm. */
		public long searched, reused, accepted, under100, from100To199, noParent, noStarts,
			nativePruned, afterBonusPruned, guaranteeRetained, retainedStops, selected;
		public final long[] acceptedStages=new long[STAGES.length], selectedStages=new long[STAGES.length];
		private static final String[] STAGES={"ddl_accept", "cosine_accept", "identity_accept", "blosum_accept"};
	}

	/** Measurement-only whole-score penalty before best-candidate selection, pruning and DP.
	 * Codons are captured in biological orientation by Orf; unknown5prime edges are exempt.
	 * A nonpositive score is preserved because multiplication would otherwise improve it. */
	static float minorStartScore(float score, int codon, boolean startTruncated, float multiplier){
		assert(multiplier>0f && multiplier<=1f && Float.isFinite(multiplier))
			: "The diagnostic multiplier must attenuate positive scores without changing their sign";
		if(multiplier==1f || !(score>0f) || startTruncated){return score;}
		// Packed A/C/G/T codons: CTG=30, ATT=15, ATC=13, ATA=12 (GeneticCode.codon).
		return codon==30 || codon==15 || codon==13 || codon==12 ? score*multiplier : score;
	}

	/** Invalid JVM diagnostic policy fails before any caller worker or output is started. */
	static float parseMinorStartMultiplier(String value){
		final float factor;
		try{factor=Float.parseFloat(value);}
		catch(NumberFormatException e){
			throw new IllegalArgumentException("Invalid bbtools.diagnostic.minorStartMultiplier: "+value, e);
		}
		if(!Float.isFinite(factor) || factor<=0f || factor>1f){
			throw new IllegalArgumentException("bbtools.diagnostic.minorStartMultiplier must be finite in (0,1]: "+value);
		}
		return factor;
	}

	/** True when a complete first edge codon may seed a truncated CDS. */
	static boolean validSyntheticEdgeStart(final int code){
		return validSyntheticEdgeStart(code, null);
	}

	static boolean validSyntheticEdgeStart(final int code, final GeneticCode geneticCode){
		return code>=0 && !isStopCodon(code, geneticCode);
	}

	/** Negative rolling codons may be below -1; never pass them into GeneticCode's strict index API. */
	private static boolean isStartCodon(final int code, final GeneticCode geneticCode){
		return code>=0 && (geneticCode==null ? GeneModel.isStartCodon(code) : geneticCode.isStart(code));
	}

	private static boolean isStopCodon(final int code, final GeneticCode geneticCode){
		return code>=0 && (geneticCode==null ? GeneModel.isStopCodon(code) : geneticCode.isStop(code));
	}
	
	/**
	 * Return a score of how well this sequence matches a feature based on inner kmer frequencies.
	 * Intended to see whether a short read matches a feature, or which side matches better.
	 * Generally this would be called twice, with bases being reverse-complemented in between.
	 * @param bases Bases, oriented on the intended strand.
	 * @param type Type of feature, e.g. ProkObject.CDS.
	 * @return Average score of best-scoring frame for this strand.
	 */
	public float scoreFeature(byte[] bases, int type){
		assert(pgm!=null) : pgm;
		StatsContainer sc=pgm.allContainers[type];
		assert(sc!=null) : type+", "+pgm;
		final FrameStats innerStats=pgm.allContainers[type].inner;
		assert(innerStats!=null) : sc;
		assert(innerStats.k>0) : innerStats;
		final int frames=(type==CDS ? 3 : 1);
		
		final int k=innerStats.k;
		final int mask=~((-1)<<(2*k));
		
		int codon=0;
		int kmer=0;
		int len=0;
		int numKmers=0;
		float[] frameScores=new float[frames];
		for(int pos=0, currentFrame=0; pos<bases.length; pos++){
			final byte b=bases[pos];
			final int x=AminoAcid.baseToNumber[b];
			codon=((codon<<2)|x);
			kmer=((kmer<<2)|x)&mask;
			
			if(x>=0){
				len++;
				if(len>=k){
					for(int frame2=0; frame2<frames; frame2++) {
						int frame=(currentFrame+frame2)%frames;
						float prob=innerStats.probs[frame][kmer];
//						float dif=prob-0.99f;//Prob above 1 is more likely than average
						frameScores[frame]+=prob;
						numKmers++;
					}
				}
			}else{
				len=0;
				kmer=0;
			}
			
			currentFrame++;
			if(currentFrame>2){currentFrame=0;}
		}
//		System.err.println(Arrays.toString(frameScores));
		return (Tools.max(frameScores)*frames)/Tools.max(1, numKmers);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** 
	 * Current gene model.
	 * TODO: Dynamically swap this as needed for contigs with varying GC.
	 */
	GeneModel pgm;
	/** Immutable run policy; changing PGM statistics never changes the selected codon assignments. */
	private final GeneticCode geneticCode;
	private TrnaCaller trnaCaller;
	/** Run-local, per-worker verifier; null is the unchanged fast caller. */
	CmRnaVerifier cmVerifier;
	/** Actual padded window of the successful legacy refinement, per worker. */
	private int cmWindowStart=-1, cmWindowStop=-1;
	public static byte[][] trnaLibrary;
	public static consensus.BaseGraph[] trnaModels;
	public static String[] trnaModelNames;

	/** Shared (read-only after CallGenes.loadNcrnaResources), one entry per generic
	 * ncRNA family with staged resources. Empty if none are staged -- makeRnas skips
	 * the ncRNA loop entirely in that case. */
	public static final ArrayList<NcrnaFamily> ncrnaFamilies=new ArrayList<>();
	/** One shared, read-only-after-build primitive map for every active 17-mer scavenger. */
	private static ConservedRnaSeedIndex conservedRnaSeedIndex;
	/** Generic-family index -> shared 17-mer slot; -1 retains the legacy non-17-mer scan. */
	private static int[] ncrnaSeedSlots=new int[0];
	/** Shared 17-mer slot for the tRNA scavenger, or -1 when inactive/non-17-mer. */
	private static int trnaSeedSlot=-1;
	/** Opt-in production equivalence assertion against the pre-unification scanners. */
	private static final boolean VERIFY_UNIFIED_SEEDS=Boolean.getBoolean("ncrna.unifiedSeedVerify");
	/** Per-instance scavengers, one per ncrnaFamilies entry, lazily built (mirrors
	 * trnaCaller's own per-instance-per-thread construction -- NcrnaScavenger's
	 * TrnaKmerIndex is per-thread-mutable, so each GeneCaller instance needs its own,
	 * built from the SHARED read-only family bundle). Null until first use. */
	private ArrayList<NcrnaScavenger> ncrnaScavengers;
	private byte[] iterationBases;
	private ConservedRnaSeedIndex.ScanResult[] iterationSeeds;
	private NcrnaCandidateReuse iterationCandidates;
	private NcrnaCpuProfile ncrnaCpuProfile;
	/** Read only after this worker joins; return scalar counters, never retain its models. */
	NcrnaCpuProfile completedCpuProfile(){
		if(ncrnaCpuProfile!=null && ncrnaScavengers!=null){
			for(int i=0; i<ncrnaScavengers.size(); i++){
				if(inactiveClanMember(i)){continue;}
				final NcrnaScavenger s=ncrnaScavengers.get(i);
				if(s instanceof NcrnaClanScavenger){
					final NcrnaClanScavenger clan=(NcrnaClanScavenger)s;
					for(int member=0; member<3; member++){
						ncrnaCpuProfile.cmEndpoints(cobalaminMemberSlots[member], clan.cmEndpointCounters(member));
					}
				}else{ncrnaCpuProfile.cmEndpoints(i, s.cmEndpointCounters());}
				ncrnaCpuProfile.activity(i, s.alignmentCount(), s.alignedBases(), s.hbmScoreCalls(),
					s.hbmBasesScored(), s.kmerHitCount(), s.windowCount());
				ncrnaCpuProfile.verification(i, s.alignmentCpuNanos(), s.shortlistCpuNanos(),
					s.hbmCpuNanos(), s.modelLengthChecks(), s.modelLengthSkips(), s.modelSpanRejects());
				if(s instanceof NcrnaClanScavenger){ncrnaCpuProfile.clanComponents((NcrnaClanScavenger)s);}
			}
		}
		return ncrnaCpuProfile;
	}
	private SsuOverlapResolver ssuCompetition;
	private long[] ncrnaOutByFamily;
	/** Per-worker shared 17-mer scans, including the optional tRNA slot; never per family. */
	private long sharedSweepPasses=0;
	public long sharedSweepPasses(){return sharedSweepPasses;}

	/** Returns completed alignment counts in the same order as ncrnaFamilies.
	 * Called only after this GeneCaller's worker thread has terminated. */
	public long[] ncrnaAlignmentCounts(){
		long[] counts=new long[ncrnaFamilies.size()];
		if(ncrnaScavengers!=null){
			assert(ncrnaScavengers.size()==counts.length);
			for(int i=0; i<ncrnaScavengers.size(); i++){
				if(inactiveClanMember(i)){continue;}
				counts[i]=ncrnaScavengers.get(i).alignmentCount();
				if(ssuCompetition!=null && Euk18sRuntimeConfig.FAMILY.equals(ncrnaFamilies.get(i).name)){
					counts[i]+=ssuCompetition.alignmentCount;
				}
			}
		}
		return counts;
	}

	long[] ncrnaOutputCounts(){return ncrnaOutByFamily==null ? new long[ncrnaFamilies.size()] : ncrnaOutByFamily;}

	/** Merge only after worker termination. Slots are original16S alignments,
	 * mitochondrial alignments and prealignment mitochondrial length skips.
	 * These are additional counters; ordinary family alignment totals are untouched. */
	public void accumulateMitoTierCounts(long[] totals){
		assert(totals!=null && totals.length==3) : "Mitochondrial work totals require original/mito/length-skip slots";
		if(ncrnaScavengers==null){return;}
		assert(ncrnaScavengers.size()==ncrnaFamilies.size()) :
			"Family order must remain stable until joined-worker work aggregation";
		for(int i=0; i<ncrnaScavengers.size(); i++){
			final NcrnaFamily family=ncrnaFamilies.get(i);
			if(!MitoSsuModelPolicy.coverage(family.name,family.modelNames)){continue;}
			final NcrnaScavenger s=ncrnaScavengers.get(i);
			totals[0]=Math.addExact(totals[0],s.original16sAlignmentCount());
			totals[1]=Math.addExact(totals[1],s.mitoAlignmentCount());
			totals[2]=Math.addExact(totals[2],s.mitoLengthSkips());
		}
	}

	public long[] ncrnaAlignedBases(){return ncrnaCostCounts(0);}
	public long[] ncrnaHbmScoreCalls(){return ncrnaCostCounts(1);}
	public long[] ncrnaHbmBasesScored(){return ncrnaCostCounts(2);}
	/** Read only after this worker terminates, in the same family order as alignment counts. */
	private long[] ncrnaCostCounts(int kind){
		assert(kind>=0 && kind<=2) : "Cost selector must name aligned bases, HBM calls, or HBM bases";
		final long[] counts=new long[ncrnaFamilies.size()];
		if(ncrnaScavengers!=null){
			assert(ncrnaScavengers.size()==counts.length) : "Family registry changed before cost aggregation";
			for(int i=0;i<counts.length;i++){
				if(inactiveClanMember(i)){continue;}
				final NcrnaScavenger s=ncrnaScavengers.get(i);
				counts[i]=(kind==0 ? s.alignedBases() : kind==1 ? s.hbmScoreCalls() : s.hbmBasesScored());
				if(kind==0 && ssuCompetition!=null && Euk18sRuntimeConfig.FAMILY.equals(ncrnaFamilies.get(i).name)){
					counts[i]+=ssuCompetition.alignedBases;
				}
			}
		}
		return counts;
	}

	/** Family-local seed occurrences after input/minLen guards; read after worker termination. */
	public long[] ncrnaKmerHitCounts(){
		final long[] counts=new long[ncrnaFamilies.size()];
		if(ncrnaScavengers!=null){
			assert(ncrnaScavengers.size()==counts.length) :
				"Family registry changed before seed-count aggregation; per-family statistics would be misattributed";
			for(int i=0; i<counts.length; i++){counts[i]=ncrnaScavengers.get(i).kmerHitCount();}
		}
		return counts;
	}

	/** Scheduled windows across both passes, including windows rejected before alignment. */
	public long[] ncrnaWindowCounts(){
		final long[] counts=new long[ncrnaFamilies.size()];
		if(ncrnaScavengers!=null){
			assert(ncrnaScavengers.size()==counts.length) :
				"Family registry changed before window-count aggregation; per-family statistics would be misattributed";
			for(int i=0; i<counts.length; i++){counts[i]=ncrnaScavengers.get(i).windowCount();}
		}
		return counts;
	}

	/** Generated voted windows across both passes, before claimed-window subtraction; read after worker termination. */
	long[] ncrnaGeneratedVotedWindowCounts(){return ncrnaGeneratedVoteCounts(false);}
	/** Generated fallback windows across both passes, before claimed-window subtraction; read after worker termination. */
	long[] ncrnaGeneratedFallbackWindowCounts(){return ncrnaGeneratedVoteCounts(true);}
	private long[] ncrnaGeneratedVoteCounts(boolean fallback){
		final long[] counts=new long[ncrnaFamilies.size()];
		if(ncrnaScavengers!=null){
			assert(ncrnaScavengers.size()==counts.length) :
				"Family registry changed before generated-window aggregation; per-family counts would be misattributed";
			for(int i=0; i<counts.length; i++){
				final NcrnaScavenger scavenger=ncrnaScavengers.get(i);
				counts[i]=fallback ? scavenger.voteWindowFallbacks : scavenger.votedWindows;
			}
		}
		return counts;
	}

	/** Per-instance NNs for ORF scoring (thread-safe copies). Null = use heuristic. */
	final CellNet orfNet;
	final CellNet orfNetLow;
	final CellNet orfNetMid;
	final CellNet orfNetHigh;
	final CellNet[] orfNetsGC;
	final CellNet[] stopNetsGC;
	final CellNet[] metaNetsGC;
	final CellNet metaNetInstance;
	float gc;

	/** Shared templates loaded once; each GeneCaller copies for thread safety. */
	public static CellNet orfNetTemplate=null;
	public static CellNet orfNetTemplateLow=null;
	public static CellNet orfNetTemplateMid=null;
	public static CellNet orfNetTemplateHigh=null;
	public static CellNet[] orfNetTemplateGC=null;
	public static float[] gcMeans=null;

	static final float GC_LOW_THRESHOLD=0.42f;
	static final float GC_HIGH_THRESHOLD=0.58f;

	public static void loadOrfNet(String path){
		if(path!=null){
			orfNetTemplate=CellNetParser.load(path);
		}
	}

	public static void loadOrfNetsGC(String lowPath, String midPath, String highPath){
		if(lowPath!=null){orfNetTemplateLow=CellNetParser.load(lowPath);}
		if(midPath!=null){orfNetTemplateMid=CellNetParser.load(midPath);}
		if(highPath!=null){orfNetTemplateHigh=CellNetParser.load(highPath);}
	}

	public static void loadOrfNetsGCN(String[] paths, float[] means){
		gcMeans=means;
		orfNetTemplateGC=new CellNet[paths.length];
		for(int i=0; i<paths.length; i++){
			if(paths[i]!=null){orfNetTemplateGC[i]=CellNetParser.load(paths[i]);}
		}
	}

	public static CellNet[] stopNetTemplateGC=null;
	public static CellNet[] metaNetTemplateGC=null;
	public static CellNet metaNetSingle=null;

	/** Skip NN scoring near contig ends; edge ORFs lack upstream/downstream NN context. */
	static final int NN_EDGE_MARGIN=21;

	public static void loadStopNetsGCN(String[] paths){
		stopNetTemplateGC=new CellNet[paths.length];
		for(int i=0; i<paths.length; i++){
			if(paths[i]!=null){stopNetTemplateGC[i]=CellNetParser.load(paths[i]);}
		}
	}

	public static void loadMetaNetsGCN(String[] paths){
		metaNetTemplateGC=new CellNet[paths.length];
		for(int i=0; i<paths.length; i++){
			if(paths[i]!=null){metaNetTemplateGC[i]=CellNetParser.load(paths[i]);}
		}
	}

	public static void loadMetaNetSingle(String path){
		metaNetSingle=CellNetParser.load(path);
	}

	/** Select the appropriate net based on contig GC content — picks the net whose training mean GC is closest. */
	CellNet selectNet(float contigGC){
		if(orfNetsGC!=null && gcMeans!=null){
			int best=0;
			float bestDist=Math.abs(contigGC-gcMeans[0]);
			for(int i=1; i<gcMeans.length; i++){
				float dist=Math.abs(contigGC-gcMeans[i]);
				if(dist<bestDist && orfNetsGC[i]!=null){
					bestDist=dist;
					best=i;
				}
			}
			if(orfNetsGC[best]!=null){return orfNetsGC[best];}
		}
		if(orfNetLow!=null && contigGC<GC_LOW_THRESHOLD){return orfNetLow;}
		if(orfNetHigh!=null && contigGC>GC_HIGH_THRESHOLD){return orfNetHigh;}
		if(orfNetMid!=null){return orfNetMid;}
		return orfNet;
	}

	/** Shared with annotated-corpus extraction so boundary training uses inference's
	 * complete-contig GC fraction, ambiguous-base convention and float rounding.
	 * Borrows a nonnull sequence without changing it; zero defined bases returns zero.
	 */
	public static float calcContigGC(byte[] bases){
		long gc=0, at=0;
		for(byte b : bases){
			int x=AminoAcid.baseToNumber[b];
			if(x==1 || x==2){gc++;}
			else if(x==0 || x==3){at++;}
		}
		return (float)(gc/(double)Tools.max(1, gc+at));
	}

	float contigGC;
	CellNet activeNet;
	CellNet activeStopNet;
	CellNet activeMetaNet;

	CellNet selectMetaNet(float contigGC){
		if(metaNetsGC!=null && gcMeans!=null){
			int best=0;
			float bestDist=Math.abs(contigGC-gcMeans[0]);
			for(int i=1; i<gcMeans.length; i++){
				float dist=Math.abs(contigGC-gcMeans[i]);
				if(dist<bestDist && metaNetsGC[i]!=null){bestDist=dist; best=i;}
			}
			if(metaNetsGC[best]!=null){return metaNetsGC[best];}
		}
		if(metaNetInstance!=null){return metaNetInstance;}
		return null;
	}

	CellNet selectStopNet(float contigGC){
		if(stopNetsGC!=null && gcMeans!=null){
			int best=0;
			float bestDist=Math.abs(contigGC-gcMeans[0]);
			for(int i=1; i<gcMeans.length; i++){
				float dist=Math.abs(contigGC-gcMeans[i]);
				if(dist<bestDist && stopNetsGC[i]!=null){bestDist=dist; best=i;}
			}
			if(stopNetsGC[best]!=null){return stopNetsGC[best];}
		}
		return null;
	}

	//Gene-calling cutoffs
	/** Minimum gene length in bases to consider */
	final int minLen;
	/** Maximum overlap allowed between genes on the same strand */
	final int maxOverlapSameStrand;
	/** Maximum overlap allowed between genes on opposite strands */
	final int maxOverlapOppositeStrand;
	/** Minimum start codon score threshold */
	final float minStartScore;
	/** Minimum stop codon score threshold */
	final float minStopScore;
	/** Minimum average inner k-mer score threshold */
	final float minInnerScore;
	/** Minimum overall ORF score threshold */
	final float minOrfScore;
	/** Minimum average score per base threshold */
	final float minAvgScore;

//	static float[] cutoff1=new float[] {0, 40f, 200f, 150f, 45f};
//	static float[] cutoff2=new float[] {0, 44f, 300f, 150f, 60f};
//	static float[] cutoff3=new float[] {0, 1.7f, 1.5f, 1.4f, 1.8f};
//	static float[] cutoff4=new float[] {0, 1.6f, 1.5f, 0.6f, 1.1f};
//	static float[] cutoff5=new float[] {0, 1.0f, 1.0f, 1.0f, 1.55f};//inner score
//	static float[] scoreMult=new float[] {1f, 8f, 200f, 175f, 12f};
//	static float[] biases=new float[] {1f, 1.17f, 1.2f, 1.2f, 1.15f};

//	static float[] cutoff1=new float[] {0, 40f, 140f, 150f, 45f};
//	static float[] cutoff2=new float[] {0, 44f, 300f, 150f, 45f};
//	static float[] cutoff3=new float[] {0, 1.7f, 1.1f, 1.3f, 1.8f};
//	static float[] cutoff4=new float[] {0, 1.6f, 1.3f, 0.5f, 1.1f};
//	static float[] cutoff5=new float[] {0, 1.0f, 1.0f, 1.0f, 1.3f};//inner score
//	static float[] scoreMult=new float[] {1f, 8f, 200f, 175f, 12f};
//	static float[] biases=new float[] {1f, 1.17f, 1.2f, 1.2f, 1.15f};

	//for k=4,2,2 
//	static float[] cutoff1=new float[] {0, 40f, 140f, 150f, 40f};
//	static float[] cutoff2=new float[] {0, 44f, 300f, 150f, 45f};
//	static float[] cutoff3=new float[] {0, 1.7f, 1.1f, 1.1f, 1.8f};
//	static float[] cutoff4=new float[] {0, 1.6f, 1.3f, 0.45f, 1.1f};
//	static float[] cutoff5=new float[] {0, 1.0f, 1.0f, 1.0f, 1.3f};//inner score
//	static float[] scoreMult=new float[] {1f, 8f, 200f, 175f, 12f};
//	static float[] biases=new float[] {1f, 1.17f, 1.2f, 1.2f, 1.22f};

//	//for k=5,3,3
//	static float[] cutoff1=new float[] {0, 40f, 300f, 300f, 135f};//sum score
//	static float[] cutoff2=new float[] {0, 44f, 300f, 400f, 100f};//orf score
//	static float[] cutoff3=new float[] {0, 1.7f, 1.5f, 1.5f, 3.5f};//start score
//	static float[] cutoff4=new float[] {0, 1.6f, 1.5f, 1.4f, 1.8f};//stop score
//	static float[] cutoff5=new float[] {0, 1.0f, 1.1f, 1.15f, 1.4f};//inner score
//	static float[] scoreMult=new float[] {1f, 8f, 80f, 120f, 5f};//score mult
//	static float[] biases=new float[] {1f, 1.17f, 1.2f, 1.2f, 1.22f};//bias for sum score

	//for k=6,3,3 - this has low scores for correct 23s stops; may try k=2 or k=4 for that end 
	//Also, 5S seems to score low very for archaea
//	public static int CDS=0, tRNA=1, /*r16S=2,*/ r23S=3, r5S=4, r18S=5, r28S=6, RNA=7;
//	static float[] cutoff1=new float[] {0, 100f, 600f, 400f, 300f};//sum score
//	static float[] cutoff2=new float[] {0, 48f, 300f, 300f, 32f};//orf score
//	static float[] cutoff3=new float[] {0, 3.0f, 1.8f, 1.6f, 4.0f};//start score
//	static float[] cutoff4=new float[] {0, 2.8f, 2.0f, 0.90f, 1.9f};//stop score
//	static float[] cutoff5=new float[] {0, 2.8f, 1.1f, 1.55f, 1.4f};//inner score
//	static float[] scoreMult=new float[] {1f, 1f, 35f, 80f, 1.0f};//score mult
//	static float[] biases=new float[] {1f, 1.25f, 1.30f, 1.30f, 1.22f};//16S bias for sum score: 1.25 best for archs, 1.4 best for bacts

////	{"CDS", "tRNA", "16S", "23S", "5S", "18S", "28S", "RNA"};
//	static float[] cutoff1=new float[] {0, 90f, 300f, 400f, 270f, 300f};//sum score
//	static float[] cutoff2=new float[] {0, 48f, 300f, 300f, 32f, 300f};//orf score
//	static float[] cutoff3=new float[] {0, 2.8f, 1.65f, 1.6f, 3.8f, 1.65f};//start score
//	static float[] cutoff4=new float[] {0, 2.6f, 1.70f, 0.90f, 1.8f, 1.70f};//stop score
//	static float[] cutoff5=new float[] {0, 2.8f, 1.70f, 1.55f, 1.4f, 1.70f};//inner score
//	static float[] scoreMult=new float[] {1f, 1f, 35f, 80f, 1.0f, 35f};//score mult
//	static float[] biases=new float[] {1f, 1.50f, 1.30f, 1.30f, 1.30f, 1.30f};

////	{"CDS", "tRNA", "16S", "23S", "5S", "18S", "28S", "RNA"};
//	static float[] cutoff1=new float[] {0, 90f, 300f, 400f, 90f, 300f};//sum score
//	static float[] cutoff2=new float[] {0, 48f, 300f, 300f, 24f, 300f};//orf score
//	static float[] cutoff3=new float[] {0, 2.8f, 1.65f, 1.6f, 2.0f, 1.65f};//start score
//	static float[] cutoff4=new float[] {0, 2.6f, 1.70f, 0.90f, 1.0f, 1.70f};//stop score
//	static float[] cutoff5=new float[] {0, 2.8f, 1.70f, 1.55f, 2.6f, 1.70f};//inner score
//	static float[] scoreMult=new float[] {1f, 1f, 35f, 80f, 1.25f, 35f};//score mult
//	static float[] biases=new float[] {1f, 1.50f, 1.30f, 1.30f, 1.50f, 1.30f};

////	{"CDS", "tRNA", "16S", "23S", "5S", "18S", "28S", "RNA"};
//	static float[] cutoff1=new float[] {0, 90f, 300f, 400f, 100f, 300f};//sum score
//	static float[] cutoff2=new float[] {0, 48f, 300f, 300f, 32f, 300f};//orf score
//	static float[] cutoff3=new float[] {0, 2.8f, 1.65f, 1.6f, 1.8f, 1.65f};//start score
//	static float[] cutoff4=new float[] {0, 2.6f, 1.70f, 0.90f, 1.0f, 1.70f};//stop score
//	static float[] cutoff5=new float[] {0, 2.8f, 1.70f, 1.55f, 2.6f, 1.70f};//inner score
//	static float[] scoreMult=new float[] {1f, 1f, 35f, 80f, 1.25f, 35f};//score mult
//	static float[] biases=new float[] {1f, 1.50f, 1.30f, 1.30f, 1.55f, 1.30f};

//	{"CDS", "tRNA", "16S", "23S", "5S", "18S", "28S", "RNA"};
	//NOTE (feature forward-lead, Eru 2026-06-19): cutoff1-5/scoreMult/biases
	//are type-indexed and sized 6 = {CDS0,tRNA1,16S2,23S3,5S4,18S5}. Adding
	//5.8S/28S/25S/ITS (the commented header above shows 28S=6,RNA=7 as the
	//intended next slots) requires extending ALL six of these arrays in lockstep,
	//new ProkObject type constants, and trained rnaContainers. See
	//plans/prok_ribo_trna_its_calling.plan.
	/** Sum score cutoffs for initial RNA detection by feature type */
	static float[] cutoff1=new float[] {0, 20f, 300f, 400f, 100f, 400f};//sum score
	/** ORF score cutoffs for refined RNA evaluation by feature type */
	static float[] cutoff2=new float[] {0, 36f, 300f, 300f, 32f, 300f};//orf score //tRNA is very sensitive here
	/** Start score thresholds for RNA boundary refinement by feature type */
	static float[] cutoff3=new float[] {0, 2.4f, 1.65f, 1.6f, 1.8f, 1.5f};//start score
	/** Stop score thresholds for RNA boundary refinement by feature type */
	static float[] cutoff4=new float[] {0, 1.5f, 1.70f, 0.90f, 1.0f, 1.1f};//stop score
	/** Inner k-mer score thresholds for RNA validation by feature type */
	static float[] cutoff5=new float[] {0, 2.2f, 1.70f, 1.55f, 2.6f, 1.5f};//inner score
	/** Score multipliers applied to final scores by feature type */
	static float[] scoreMult=new float[] {1f, 1.0f, 35f, 80f, 1.25f, 35f};//score mult
	/** Bias values for k-mer probability calculations by feature type */
	static float[] biases=new float[] {1f, 1.45f, 1.30f, 1.30f, 1.55f, 1.50f};
	
	/** Counter for total gene stops created during processing */
	long geneStopsMade=0;
	/** Counter for total gene starts created during processing */
	long geneStartsMade=0;
	/** Counter for gene starts retained after filtering */
	long geneStartsRetained=0;
	/** Counter for gene stops retained after filtering */
	long geneStopsRetained=0;
	/** Counter for final gene starts in output path */
	long geneStartsOut=0;

	/** Counter for tRNA features in final output */
	long tRNAOut=0;
	/** Counter for 16S rRNA features in final output */
	long r16SOut=0;
	/** Counter for 23S rRNA features in final output */
	long r23SOut=0;
	/** Counter for 5S rRNA features in final output */
	long r5SOut=0;
	/** Counter for 18S rRNA features in final output */
	long r18SOut=0;
	/** Counter for generic ncRNA (RNase P, SRP, ...) features in final output */
	long rnaOut=0;

	/** rrnafallback attempt counters (Citan-authorized, 2026-09-02), indexed by
	 * ProkObject type constant (CDS=0..RNA=7). Unconditional, always collected (same cost
	 * class as r16SOut etc. above, none of which are opt-in gated) -- per the 203-genome
	 * A/B benchmark contract's stratify-by-5S/16S/23S requirement. Incremented ONLY inside
	 * the OUTER refineByAlignment(Orf,byte[],int,StatsContainer) loop -- the seam is
	 * deliberately outer-loop-only, never inside refineByAlignment_IDA/_SSA, so a single
	 * counting site covers both aligner modes identically.
	 * <p>A candidate counts as "entering refinement" only once a NONEMPTY consensus array
	 * is confirmed present (Citan's decision, 2026-09-02) -- the empty/no-consensus early
	 * return is explicitly not an alignment candidate. Consequently
	 * rrnaCandidatesEnteringRefinement and rrnaUniversalAttempts are always numerically
	 * IDENTICAL under current control flow (every nonempty-consensus candidate runs
	 * exactly one i==0 iteration, unconditionally) -- both fields are retained anyway as
	 * explicitly equal-under-current-control-flow contract fields (Citan's decision,
	 * 2026-09-02), matching the batch-run plan's own two separately-named counter columns,
	 * not collapsed into one. */
	/** Per-type candidate Orfs reaching refineRna before its internal gates. Diagnostic only. */
	long[] rrnaCandidatesGenerated=new long[8];
	long[] rrnaRefineByAlignmentCalls=new long[8];
	long[] rrnaGateScoreFail=new long[8], rrnaGateGoodFail=new long[8],
		rrnaGateBothFail=new long[8], rrnaGateAdmitted=new long[8];
	long[] rrnaG5Reject=new long[8], rrnaG6Reject=new long[8];
	public long[] rrnaCandidatesEnteringRefinement=new long[8];
	/** The universal (first, i==0) consensus entry is attempted exactly once per candidate
	 * that reaches the loop, in BOTH rrnafallback=f and =t modes. See field-group javadoc
	 * above for why this is currently always equal to rrnaCandidatesEnteringRefinement. */
	public long[] rrnaUniversalAttempts=new long[8];
	/** A fallback (i>0) attempt is only possible when rrnafallback=t AND the universal
	 * attempt failed. */
	public long[] rrnaFallbackAttempts=new long[8];
	/** Incremented once per candidate whose FINAL accepted attempt had index>0 -- i.e. the
	 * universal entry failed but a later fallback entry rescued it. */
	public long[] rrnaFallbackOnlyRescues=new long[8];
	/** Incremented once per candidate where every attempted consensus entry failed (the
	 * loop exited with refined==false). */
	public long[] rrnaRejectedCandidates=new long[8];
	/** Per-type histogram of attempts-made-per-candidate (index = attempt count) -- lets a
	 * downstream wrapper compute EXACT mean/median/p95 attempts-per-candidate per type
	 * from the merged histogram, not an estimate. Each row is GROWABLE, not fixed-size
	 * (Citan's correction, 2026-09-02): a row starts null and grows via
	 * growHistRow/Arrays.copyOf only as far as an actually-observed attempt count
	 * requires, so a future resource with a longer consensus shortlist than any type
	 * ships today still gets exact (not clamped/lossy) p95 data. */
	public long[][] rrnaAttemptCountHist=new long[8][];

	/** Grows (never shrinks) a histogram row so index `idx` is valid, preserving existing
	 * counts and zero-filling new slots. `row` may be null (treated as length 0). */
	static long[] growHistRow(long[] row, int idx){
		if(row==null){return new long[idx+1];}
		if(row.length<=idx){return Arrays.copyOf(row, idx+1);}
		return row;
	}

	/** Score tracker for CDS features during processing */
	ScoreTracker stCds=new ScoreTracker(CDS);
	/** Secondary score tracker for CDS features */
	ScoreTracker stCds2=new ScoreTracker(CDS);
	/** Score tracker for CDS features passing final filters */
	ScoreTracker stCdsPass=new ScoreTracker(CDS);
	/** Score tracker for tRNA features */
	ScoreTracker sttRNA=new ScoreTracker(tRNA);
	/** Score tracker for 16S rRNA features */
	ScoreTracker st16s=new ScoreTracker(r16S);
	/** Score tracker for 23S rRNA features */
	ScoreTracker st23s=new ScoreTracker(r23S);
	/** Score tracker for 5S rRNA features */
	ScoreTracker st5s=new ScoreTracker(r5S);
	/** Score tracker for 18S rRNA features */
	ScoreTracker st18s=new ScoreTracker(r18S);
	
	/** Array of all score trackers for different feature types */
	ScoreTracker[] trackers=new ScoreTracker[] {stCdsPass, sttRNA, st16s, st23s, st5s, st18s};
	
	/**
	 * Gets thread-local SingleStateAlignerFlat2 instance.
	 * Creates new instance if none exists for current thread.
	 * @return SingleStateAlignerFlat2 for current thread
	 */
	public static SingleStateAlignerFlat2 getSSA(){
		SingleStateAlignerFlat2 ssa=localSSA.get();
		if(ssa==null){
			synchronized(localSSA){
				ssa=localSSA.get();
				if(ssa==null){
					ssa=new SingleStateAlignerFlat2();
					localSSA.set(ssa);
				}
			}
		}
		return ssa;
	}
	
	/**
	 * Gets thread-local SingleStateAlignerFlat3 instance.
	 * Creates new instance if none exists for current thread.
	 * @return SingleStateAlignerFlat3 for current thread
	 */
	public static SingleStateAlignerFlat3 getSSA3(){
		SingleStateAlignerFlat3 ssa=localSSA3.get();
		if(ssa==null){
			synchronized(localSSA3){
				ssa=localSSA3.get();
				if(ssa==null){
					ssa=new SingleStateAlignerFlat3();
					localSSA3.set(ssa);
				}
			}
		}
		return ssa;
	}
	
	/**
	 * Gets thread-local SingleStateAlignerFlatFloat instance.
	 * Creates new instance if none exists for current thread.
	 * @return SingleStateAlignerFlatFloat for current thread
	 */
	public static SingleStateAlignerFlatFloat getSSAF(){
		SingleStateAlignerFlatFloat ssa=localSSAF.get();
		if(ssa==null){
			synchronized(localSSAF){
				ssa=localSSAF.get();
				if(ssa==null){
					ssa=new SingleStateAlignerFlatFloat();
					localSSAF.set(ssa);
				}
			}
		}
		return ssa;
	}
	
	/** Whether to retain at least one ORF even if it fails quality filters */
	public boolean keepAtLeastOneOrf=false;

	/** Test-only, null-by-default observer for ordered rRNA consensus attempts. */
	private NcrnaBoundaryTrace.Worker boundaryTrace=null;
	private RefinementAttemptSink attemptSink=null;
	public void setAttemptSink(RefinementAttemptSink sink){attemptSink=sink;}
	private NcrnaStageDiagSink ncrnaDiagSink=null;
	void setNcrnaDiagSink(NcrnaStageDiagSink sink){
		ncrnaDiagSink=sink;
		if(ncrnaScavengers!=null){
			for(int i=0;i<ncrnaScavengers.size();i++){
				ncrnaScavengers.get(i).setStageDiagSink(sink,ncrnaFamilies.get(i).name);
			}
		}
	}
	/** Accepted refinement identity is needed for16S scoring even without a diagnostic sink. */
	private float lastAttemptIdentity=-1f;
	/** Separate observer for generated vote proposals; existing stage rows stay unchanged. */
	private NcrnaVoteDiagSink ncrnaVoteDiagSink;
	void setNcrnaVoteDiagSink(NcrnaVoteDiagSink sink){
		ncrnaVoteDiagSink=sink;
		if(ncrnaScavengers!=null){
			assert(ncrnaScavengers.size()==ncrnaFamilies.size()) :
				"Stable family ordering is required to label existing scavenger vote observers";
			for(int i=0; i<ncrnaScavengers.size(); i++){
				ncrnaScavengers.get(i).setVoteDiagSink(sink, ncrnaFamilies.get(i).name);
			}
		}
	}
	private String lastAttemptReason="UNKNOWN";

	/** Thread-local storage for SingleStateAlignerFlat2 instances */
	private static ThreadLocal<SingleStateAlignerFlat2> localSSA=new ThreadLocal<SingleStateAlignerFlat2>();
	/** Thread-local storage for SingleStateAlignerFlat3 instances */
	private static ThreadLocal<SingleStateAlignerFlat3> localSSA3=new ThreadLocal<SingleStateAlignerFlat3>();
	/** Thread-local storage for SingleStateAlignerFlatFloat instances */
	private static ThreadLocal<SingleStateAlignerFlatFloat> localSSAF=new ThreadLocal<SingleStateAlignerFlatFloat>();
//	public static int maxAlignmentEndpointDifference=15;
	/** Padding bases added around ORF boundaries for alignment */
	public static int alignmentPadding=300;
	/** Padding bases for IDAligner-specific alignments */
	public static int alignmentPaddingIDA=80;//TODO: Recur for edge alignments, and add start hint
	/** Whether to use IDAligner instead of SingleStateAligner for refinement */
	public static boolean useIDAligner=false;//Quantum/Drifting are fastest
	
	/** Maximum number of alternative starts to generate per stop codon */
	public static int breakLimit=12;
	/** Immutable experiment policy shared by all refinement/final callers in this JVM.
	 * Default1.0 bypasses the hook; it changes no genetic-code permissions or model counts. */
	private static final float MINOR_START_MULTIPLIER=parseMinorStartMultiplier(
		System.getProperty("bbtools.diagnostic.minorStartMultiplier", "1.0"));
	/** Lookback distance for plus-strand path scoring */
	public static int lookbackPlus=70;
	/** Lookback distance for minus-strand path scoring */
	public static int lookbackMinus=25;

	/** Measurement-only JVM switch. Default false preserves production path scores.
	 * Applies to incoming RNA edges only; incoming CDS rules, predecessor selection,
	 * run lengths, overlaps and candidate generation are unchanged. */
	public static final boolean FLAT_RNA_TRANSITIONS=Boolean.getBoolean("bbtools.diagnostic.flatRnaTransitions");
	/** G11-approved measurement only: +30 on adjacent same-strand euk5S edges.
	 * Static false removes the predicate from ordinary production scoring. */
	public static final boolean EUK5S_ARRAY_BONUS=Boolean.getBoolean("bbtools.diagnostic.euk5sArrayBonus");
	public static final int EUK5S_ARRAY_MAX_GAP=Integer.parseInt(
		System.getProperty("bbtools.diagnostic.euk5sArrayMaxGap","1000"));
	static{
		if(EUK5S_ARRAY_MAX_GAP<0 || (EUK5S_ARRAY_BONUS && FLAT_RNA_TRANSITIONS)){
			throw new IllegalArgumentException("Array-gap limit must be nonnegative; array and global "+
				"flat-RNA experiments cannot be combined");
		}
	}
	
//	pathScore+=p0+p1*(Tools.mid(p5*(p2+pathLength), p6*(p3-pathLength), p4));
	
	//Operon length modifiers for same strand
	/** Base penalty constant for same-strand operon scoring */
	public static float p0=-30f;
	/** Length-dependent coefficient for same-strand operon scoring */
	public static float p1=-0.35f; //Sensitive - higher decreases FP and TP
	/** Positive length bias parameter for same-strand scoring */
	public static float p2=4.0f;//Insensitive - higher positive decreases FP and TP
	/** Negative length bias parameter for same-strand scoring */
	public static float p3=12f; //Higher decreases FP (substantially) and TP (slightly)
	/** Minimum penalty floor for same-strand scoring */
	public static float p4=-10f; //Lower decreases FP (weakly) and TP (greatly)
	/** Positive length scaling factor for same-strand scoring */
	public static float p5=2.0f; //Insensitive - lower increases FP and TP
	/** Negative length scaling factor for same-strand scoring */
	public static float p6=2f; //Greater increases FP and TP
	
//	pathScore+=q1+Tools.mid(q2*prevLength, q3+q4*prevLength, q5);
	
	//Operon length modifiers for opposite strand
	/** Base penalty constant for opposite-strand transitions */
	public static float q1=-36f;
	/** Previous length coefficient for opposite-strand scoring */
	public static float q2=-1.6f; //q2 and q4 must have opposite signs
	/** Base value for opposite-strand length adjustment */
	public static float q3=-12f;
	/** Previous length multiplier for opposite-strand scoring */
	public static float q4=3.0f; //(Lower [even negative] decreases FP and TP)
	/** Maximum penalty ceiling for opposite-strand scoring */
	public static float q5=-40f;
	
	/** Output stream for debugging messages */
	private static PrintStream outstream=System.err;
	/** Whether to print verbose debugging output */
	static boolean verbose;
	static boolean dumpOrfScores;
	
}
