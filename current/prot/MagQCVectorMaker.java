package prot;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.TreeSet;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.IntHashSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.IntList;
import structures.IntLongHashMap;
import tax.TaxNode;
import tax.TaxTree;
import tracker.KmerTracker;
import stream.bam.BgzfSettings;

/**
 * Generates MAG-QC training vectors by synthesizing bins from the per-contig
 * precompute cache (Stage 1) and emitting each bin's feature vector with its
 * ACHIEVED completeness and contamination as continuous regression targets.
 *
 * <p>A synthetic bin = a random subset of one TARGET organism's contigs (to hit a
 * sampled completeness against the target's true genome size, plasmids included)
 * plus, optionally, contigs from one or more CONTAMINANT organisms (to hit a sampled
 * contamination). Labels are computed from the explicit selected target, not a
 * re-inferred winner: completeness = target_bp / genomeSize[target];
 * contamination = foreign_bp / (target_bp + foreign_bp). Foreign target bases are
 * c/(1-c)*cleanBases. The generator emits the ACHIEVED labels after whole-contig
 * selection.
 *
 * <p>Organisms are split globally BEFORE sampling; a held-out organism appears as
 * neither target nor contaminant in the training file (Barbara design review #3).
 * The feature vector uses the same cache the deployment path will produce, so the
 * aggregation here is the train==serve contract.
 *
 * <p>Vector layout: [family features] + [13 global stats] + [phylum one-hot], then
 * the two outputs completeness, contamination. The family-count encoding is selected
 * by enc= : ratio (default) N/(1+N); raw min(N,32)/32; log log2(1+N)/log2(65); two =
 * two columns per family (presence 0/1, then excess-copies min(N-1,16)/16); norm = 0 if
 * absent else N/avgCopyWhenPresent[family] (per-family baseline so present-at-typical is
 * ~1.0 and duplication is &gt;1.0, calibrated to each family's natural copy number). The
 * two/norm schemes keep the duplication signal that flags contamination uncompressed
 * (CheckM2 uses raw counts for exactly this reason).
 *
 * <p>Idiom rewrite (2026-08): I/O runs on {@link ByteFile}/{@link ByteStreamWriter}
 * with a reused {@link LineParser1} (byte-range parsing, zero per-field String
 * allocation except where a value must persist as a String or is on a load-once,
 * bounded-frequency path). Every tid-keyed lookup used inside the per-bin hot path
 * ({@code makeBin}) is a primitive {@link structures.IntHashMap}/{@link IntLongHashMap}
 * or a dense array (byTid's contig lists, since IntHashMap can't hold list values) -
 * no boxed {@code HashMap<Integer,...>} on the hot path. {@code writeSet}'s per-attempt
 * buffers ({@code fam}/{@code glob}/{@code labels}) and the per-bin {@link Agg}
 * accumulator are now hoisted out of the attempt loop and reused (cleared, not
 * reallocated). {@code fmt()}'s exact string representation (whole numbers print
 * without a decimal point, everything else at 6 fixed decimals) is replicated by
 * {@link #appendFmt(ByteBuilder, double)} for the hot output paths; the String-
 * returning {@link #fmt(double)} is KEPT for the aggregator's roundtrip-through-string
 * rounding (low frequency, and an actual String is what {@code Float.parseFloat} needs).
 * NOT threaded this pass (Brian's explicit constraint): {@code makeBin} draws from
 * one sequential {@link Random} stream per output set, and parallelizing bin synthesis
 * would reorder those draws and change the output. Algorithm, RNG call sequence, and
 * every output byte are UNCHANGED - see {@code MagQCVectorMakerTest}/
 * {@code MagQCAggVectorTest} for the pinned behavioral contract, verified against the
 * original by UMP45's differential (byte-identical across all fixture modes).
 */
public class MagQCVectorMaker implements Cloneable {

	/** Resource-only construction is restricted to the validated prepared-input factory. */
	private MagQCVectorMaker(){}

	/**
	 * Makes a replay-only worker view. Loaded caches, vocabularies, indexes and bundle
	 * masters are immutable after process() staging and remain shared; every scratch object
	 * touched by replayBin/format* is private to the worker. This is deliberately not a
	 * general-purpose copy operation.
	 */
	private MagQCVectorMaker copyForReplayThread(){return copyForReplayThread(false);}

	/** Prepared locked workers retain subnet bindings, but never share formatter scratch. */
	private MagQCVectorMaker copyForReplayThread(boolean lockedSubnets){
		try{
			final MagQCVectorMaker x=(MagQCVectorMaker)super.clone();
			x.replayDiag=null;
			x.sharedPreparedSubnets=lockedSubnets;
			x.preparedWorker=preparedDenseInference;
			x.lastContributionCounts=null;
			x.anticodonRows=0; x.anticodonMissingRows=0;
			x.refCdsNativeIdx=new IntList(); x.refCdsForeignIdx=new IntList();
			x.benchNative=new ArrayList<Contig>(); x.benchForeign=new ArrayList<Contig>();
			x.lastNcServe=lastNcServe.clone(); x.lastAntiObs=lastAntiObs.clone(); x.lastAntiServe=lastAntiServe.clone();
			x.lastNcObs=lastNcObs.clone();
			x.lastFamObs=lastFamObs==null ? null : lastFamObs.clone();
			x.cleanFamBuf=cleanFamBuf==null ? null : cleanFamBuf.clone();
			x.rowCtx=rowCtx==null ? null : rowCtx.clone();
			x.quantizedRowCtx=quantizedRowCtx==null ? null : quantizedRowCtx.clone();
			x.subnetSharedInputs=null; x.subnetObservedInputs=null;
			x.preparedGlob=preparedGlob==null ? null : preparedGlob.clone();
			x.observedByItem=observedByItem==null ? null : observedByItem.clone();
			if(aggSubnets!=null){
				x.aggSubnets=new ArrayList<AggSubnet>(aggSubnets.size());
				for(AggSubnet s : aggSubnets){
					final AggSubnet c=new AggSubnet(); c.name=s.name; c.type=s.type; c.numObs=s.numObs;
					c.ranks=s.ranks;
					if(lockedSubnets){
						assert(preparedDenseInference && s.bundleSubnet!=null) : "Locked subnet sharing requires the prepared resource master";
						c.bundleSubnet=s.bundleSubnet; c.preparedOutput=new float[s.bundleSubnet.expectedOutputs];
					}else if(s.bundleSubnet!=null){
						// Replay workers own explicit InferenceNet instances.  Detaching the
						// worker view from bundleSubnet is intentional: it prevents the replay
						// hot path from consulting the bundle's ThreadLocal on every subnet call.
						c.net=s.bundleSubnet.newReplayWorker(); c.bundleSubnet=null;
					}else{
						c.net=s.net==null ? null : s.net.copy(false); c.bundleSubnet=null;
					}
					c.buf=s.buf==null ? null : s.buf.clone(); x.aggSubnets.add(c);
				}
			}
			return x;
		}catch(CloneNotSupportedException e){throw new AssertionError(e);}
	}

	public static void main(String[] args){
		MagQCVectorMaker x=new MagQCVectorMaker(args);
		Throwable failure=null;
		try{x.process();}
		catch(RuntimeException e){failure=e; throw e;}
		catch(Error e){failure=e; throw e;}
		finally{x.closeContributionIndex(failure);}
	}

	public MagQCVectorMaker(String[] args){
		// A null display class suppresses command-line echo: legacy exact-source
		// bindings may contain full digests, while operator-facing records use sha80.
		args=new PreParser(args,null,false).args;
		for(String arg : args){
			int eq=arg.indexOf('=');
			if(eq<0){continue;}
			String a=arg.substring(0, eq).toLowerCase(), b=arg.substring(eq+1);
			if(a.equals("cache")){cacheFile=b;}
			else if(a.equals("sizemap")){sizemapFile=b;}
			else if(a.equals("familylist")){familyFile=b;}
			else if(a.equals("features")){featuresFile=b;}
			else if(a.equals("taxpgm")){taxpgmFile=b;}
			else if(a.equals("labels")){labelsFile=b;}
			else if(a.equals("tree")){treeFile=b;}
			else if(a.equals("out")){out=b;}
			else if(a.equals("outval")){outval=b;}
			else if(a.equals("emitdense")){emitDense=parseBool(b);}
			else if(a.equals("n")){n=Long.parseLong(b);}
			else if(a.equals("valn")){valn=Long.parseLong(b);}
			else if(a.equals("valfrac")){valfrac=Double.parseDouble(b);}
			else if(a.equals("seed")){seed=Long.parseLong(b);}
			else if(a.equals("splitseed")){splitSeed=Long.parseLong(b); splitSeedSet=true;}
			else if(a.equals("minlen")){minlen=Integer.parseInt(b);}
			else if(a.equals("mixcomp")){mixComp=Double.parseDouble(b);}
			else if(a.equals("mixcont")){mixCont=Double.parseDouble(b);}
			else if(a.equals("cleanspike")){cleanSpike=Double.parseDouble(b);}
			else if(a.equals("multicontamprob")){multiContamProb=Double.parseDouble(b);}
			else if(a.equals("perfectfrac")){perfectFrac=Double.parseDouble(b);}
			else if(a.equals("nearperfectfrac")){nearPerfectFrac=Double.parseDouble(b);}
			else if(a.equals("extraperfectfrac")){extraPerfectFrac=Double.parseDouble(b);}
			else if(a.equals("extranearperfectfrac")){extraNearPerfectFrac=Double.parseDouble(b);}
			else if(a.equals("extraspikeset")){extraSpikeSet=b.toLowerCase();}
			else if(a.equals("unshreddedcache")){unshreddedCacheFile=b;}
			else if(a.equals("extraunshreddedperfectfrac")){extraUnshreddedPerfectFrac=Double.parseDouble(b);}
			else if(a.equals("benchtruth")){benchTruthFile=b; benchMode=true;}
			else if(a.equals("benchmanifest")){benchManifestFile=b; benchMode=true;}
			else if(a.equals("benchvec")){benchVecFile=b; benchMode=true;}
			else if(a.equals("benchvecmono")){benchVecMonoFile=b; benchMode=true;}
			else if(a.equals("benchbins")){benchBins=Long.parseLong(b);}
			else if(a.equals("samefamprob")){sameFamProb=Double.parseDouble(b);}
			else if(a.equals("enc")){enc=parseEnc(b);}
			else if(a.equals("subnet")){subnetName=b.toLowerCase();}
			else if(a.equals("subnetout")){subnetOut=b;}
			else if(a.equals("subnetvalout")){subnetValOut=b;}
			else if(a.equals("subnetlabels")){subnetLabels=parseSubnetLabels(b);}
			else if(a.equals("subnetfeatures")){subnetFeatures=parseSubnetFeatures(b);}
			else if(a.equals("expectedcopytable")){expectedCopyTableFile=b;}
			else if(a.equals("expectedcopytablesha256") || a.equals("expectedcopytablesha80")){expectedCopyTableSha256=b;}
			else if(a.equals("subnetpopulations")){subnetPopulationsFile=b;}
			else if(a.equals("subnetpopulationssha256") || a.equals("subnetpopulationssha80")){subnetPopulationsSha256=b;}
			else if(a.equals("subsetfile")){subsetFile=b;}
			else if(a.equals("aggmanifest")){aggManifestFile=b;}
			else if(a.equals("bundle")){aggBundleFile=b;}
			else if(a.equals("subnetmanifest")){releaseManifestFile=b;}
			else if(a.equals("subnetmanifestsha80")){releaseManifestSha80=b;}
			else if(a.equals("aggout")){aggOut=b;}
			else if(a.equals("aggvalout")){aggValOut=b;}
			else if(a.equals("aggids")){
				if(!("t".equalsIgnoreCase(b) || "true".equalsIgnoreCase(b) || "1".equals(b) ||
					"f".equalsIgnoreCase(b) || "false".equalsIgnoreCase(b) || "0".equals(b))){
					throw new IllegalArgumentException("aggids requires t or f");
				}
				aggIds=parseBool(b);
			}
			else if(a.equals("densehead")){denseHead=Integer.parseInt(b);}
			else if(a.equals("aggobs")){aggObsServe=parseAggObs(b);}
			else if(a.equals("poolmode")){poolMode=parsePoolMode(b);}
			else if(a.equals("paneltids")){panelTidsFile=b;}
			else if(a.equals("panelout")){panelOutFile=b;}
			else if(a.equals("binmanifestout")){binManifestOutFile=b;}
			else if(a.equals("binmanifestin")){binManifestInFile=b;}
			else if(a.equals("binmodel")){binModel=b;}
			else if(a.equals("replaytmp")){replayTmpDir=b;}
			else if(a.equals("replayseed")){replaySeed=Long.parseLong(b); replaySeedSet=true;}
			else if(a.equals("replayshuffle")){replayShuffle=parseBool(b);}
			else if(a.equals("replaymemory")){replayMemory=parseBool(b);}
			else if(a.equals("replaydiag")){replayDiagnostics=parseBool(b);}
			else if(a.equals("replaythreads")){replayThreads=Integer.parseInt(b);}
			else if(a.equals("replayblockrows")){replayBlockRows=Integer.parseInt(b);}
			else if(a.equals("replayparts")){replayParts=Integer.parseInt(b);}
			else if(a.equals("replaytrainstart")){replayTrainStart=Parse.parseKMG(b); replayTrainStartSet=true;}
			else if(a.equals("replaytrainend")){replayTrainEnd=Parse.parseKMG(b); replayTrainEndSet=true;}
			else if(a.equals("replayvalstart")){replayValStart=Parse.parseKMG(b); replayValStartSet=true;}
			else if(a.equals("replayvalend")){replayValEnd=Parse.parseKMG(b); replayValEndSet=true;}
			else if(a.equals("d39")){d39Sampling=parseBool(b);}
			else if(a.equals("samplemodel")){sampleModel=b;}
			else if(a.equals("splitmod")){splitModulus=Integer.parseInt(b); splitModulusSet=true;}
			else if(a.equals("samplingtaxonomy")){samplingTaxonomyFile=b;}
			else if(a.equals("samplingunknown")){
				if(!"error".equals(b) && !"ordinary".equals(b)){throw new IllegalArgumentException("samplingunknown must be error or ordinary");}
				samplingUnknownOrdinary="ordinary".equals(b);
			}
			else if(a.equals("samplephylum")){samplePhylum=b;}
			else if(a.equals("excludetids")){excludedTidsFile=b;}
			else if(a.equals("referencecdsindex")){referenceCdsIndexFile=b;}
			else if(a.equals("referencecdscontributionindex")){referenceCdsContributionIndexFile=b;}
			else if(a.equals("referencecdscontributionindexsha256")){referenceCdsContributionIndexSha256=b;}
			else if(a.equals("referencecdscontributionmanifestsha256")){referenceCdsContributionManifestSha256=b;}
			else if(a.equals("referencecdscontributioncachesha256")){referenceCdsContributionCacheSha256=b;}
			else if(a.equals("refcdstable")){refCdsTableFile=b;}
			else if(a.equals("refcdstablesha256")){refCdsTableSha256=b;}
			else if(a.equals("noorganismholdout")){noOrganismHoldout=parseBool(b);}
			else if(a.equals("deterministic")){if(Parse.parseBoolean(b)){shared.Shared.SIMD_FMA=false; shared.Shared.SIMD_FEED_FORWARD=false;}}
			else{System.err.println("Warning: unknown arg "+arg);}
		}
		if(cacheFile==null || sizemapFile==null || (emitDense && out==null)){
			throw new RuntimeException("Required: cache= sizemap= (out=|emitdense=f) (taxpgm=|labels=) [outval= familylist= tree= n= valn= ...]");
		}
		if(!emitDense){
			if(!d39Sampling && binManifestInFile==null){throw new RuntimeException("emitdense=f requires d39=t generation or manifest replay");}
			if(!((subnetName!=null && subnetOut!=null) || (aggOut!=null && (aggManifestFile!=null || aggBundleFile!=null)))){
				throw new RuntimeException("emitdense=f requires a subnet or aggregator training output");
			}
		}
		// Step 4 of the rebuild path (UMP45, 2026-09-03): labels=<C1 labels.tsv> REPLACES taxpgm=
		// as the domain/phylum source (deployment-faithful QuickClade predictions, not TID/GTDB
		// ground truth) -- exactly one of the two is required, never both (ambiguous which wins)
		// and never neither (no domain/phylum source at all).
		if((taxpgmFile==null)==(labelsFile==null)){
			throw new RuntimeException("Exactly one of taxpgm= or labels= is required (not both, not neither).");
		}
		if(aggManifestFile!=null && aggBundleFile!=null){
			throw new RuntimeException("aggmanifest= and bundle= are mutually exclusive");
		}
		if((releaseManifestFile==null)!=(releaseManifestSha80==null) ||
				(releaseManifestFile!=null && aggBundleFile==null)){
			throw new IllegalArgumentException("subnetmanifest and subnetmanifestsha80 require each other and bundle=");
		}
		// subnetfeatures=six (opt-in, UMP45 design 2026-09-09, revision2): four learned outputs of a
		// four-output subnet plus two table-derived expected-copy ratios, replacing the legacy
		// [ratio,log_obs,log_pred,zero_flag]+pooled-ratio block. Every check below fires here so a
		// misconfigured six run fails before any writer opens (Part A of the design's rejection matrix).
		if(subnetFeatures==SUBNETFEATURES_SIX && aggBundleFile==null){
			throw new RuntimeException("subnetfeatures=six requires bundle=");
		}
		if(subnetFeatures==SUBNETFEATURES_SIX){
			if(expectedCopyTableFile==null){throw new RuntimeException("subnetfeatures=six requires expectedcopytable=");}
			if(!MagQCExpectedCopyTableTyped.validDigest(expectedCopyTableSha256)){
				throw new RuntimeException("subnetfeatures=six requires expectedcopytablesha80= (20 hex), or legacy sha256");}
			if(subnetPopulationsFile==null){throw new RuntimeException("subnetfeatures=six requires subnetpopulations=");}
			if(!MagQCExpectedCopyTableTyped.validDigest(subnetPopulationsSha256)){
				throw new RuntimeException("subnetfeatures=six requires subnetpopulationssha80= (20 hex), or legacy sha256");}
		}else{
			if(expectedCopyTableFile!=null || expectedCopyTableSha256!=null
					|| subnetPopulationsFile!=null || subnetPopulationsSha256!=null){
				throw new RuntimeException("subnetfeatures=legacy does not consume expected-copy inputs; "
					+"remove them or set subnetfeatures=six");
			}
		}
		if(subnetFeatures==SUBNETFEATURES_SIX && !aggObsServe){
			throw new RuntimeException("subnetfeatures=six requires aggobs=serve (whole-bin); the four-output "
				+"subnets are trained on whole-bin observations, so aggobs=clean would be a train/serve mismatch");
		}
		taxonomyFromLabels=(labelsFile!=null);
		if(!extraSpikeSet.equals("train") && !extraSpikeSet.equals("val") && !extraSpikeSet.equals("both")){
			throw new RuntimeException("Unknown extraspikeset="+extraSpikeSet+" (train|val|both)");
		}
		if(extraUnshreddedPerfectFrac>0 && unshreddedCacheFile==null){
			throw new RuntimeException("extraunshreddedperfectfrac>0 requires unshreddedcache=");
		}
		if(binManifestOutFile!=null && binManifestOutFile.equals(cacheFile)){
			throw new RuntimeException("binmanifestout must differ from cache=");
		}
		if(binManifestOutFile!=null && binManifestInFile!=null){throw new RuntimeException("binmanifestout= and binmanifestin= are mutually exclusive");}
		if(aggIds && (binManifestInFile==null || (aggOut==null && aggValOut==null))){
			throw new IllegalArgumentException("aggids=t requires binmanifestin= and an aggregator output");
		}
		if(binModel!=null){
			MagQCBinManifest.validateModelName(binModel);
			if(binManifestInFile==null){throw new RuntimeException("binmodel= requires binmanifestin=");}
			if(benchMode || panelTidsFile!=null){throw new RuntimeException("binmodel= is supported only for manifest replay");}
		}
		if(replayShuffle && !replaySeedSet){throw new IllegalArgumentException("replayshuffle=t requires replayseed=");}
		if(replaySeedSet && !replayShuffle){throw new IllegalArgumentException("replayseed= requires replayshuffle=t");}
		if((replayShuffle || replayMemory || replayParts>1) && binModel==null){
			throw new IllegalArgumentException("Replay shuffle, memory and part options require binmanifestin= with binmodel=");
		}
		if(replayTmpDir!=null && binModel==null){throw new RuntimeException("replaytmp= requires binmodel= replay");}
		if(replayThreads<1 || replayThreads>64){throw new IllegalArgumentException("replaythreads must be 1..64: "+replayThreads);}
		if(replayBlockRows<0 || replayBlockRows>64){throw new IllegalArgumentException("replayblockrows must be 0..64: "+replayBlockRows);}
		if(replayBlockRows>0 && binManifestInFile==null){throw new RuntimeException("replayblockrows= requires binmanifestin= replay");}
		if(replayParts<1 || replayParts>64){throw new IllegalArgumentException("replayparts must be 1..64: "+replayParts);}
		if(hasReplayRange() && (binManifestInFile==null || binModel==null)){
			throw new RuntimeException("replaytrain*/replayval* bounds require binmanifestin= with binmodel=");
		}
		if(splitModulusSet && !d39Sampling){throw new RuntimeException("splitmod= applies only to d39=t generation; replay preserves stored splits");}
		if(samplingTaxonomyFile!=null || samplePhylum!=null){
			if(!d39Sampling || samplingTaxonomyFile==null || samplePhylum==null || samplePhylum.isEmpty() || !taxonomyFromLabels){
				throw new RuntimeException("Phylum sampling requires d39=t, samplingtaxonomy=, samplephylum=, and C1 labels= input features");
			}
		}
		if(d39Sampling){
			MagQCBinSplit.validate(splitModulus);
			if(binManifestInFile!=null || binManifestOutFile==null || benchMode || panelTidsFile!=null){
				throw new RuntimeException("d39=t requires binmanifestout= generation, not replay/bench/panel");
			}
			MagQCBinManifest.validateModelName(sampleModel);
			if(!noOrganismHoldout || subnetLabels!=SUBNETLABELS_GENE ||
					refCdsTableFile==null || referenceCdsContributionIndexFile==null){
				throw new RuntimeException("d39=t requires noorganismholdout=t, subnetlabels=gene, and both exact gene sources");
			}
			if(n<1 || n>Integer.MAX_VALUE || valn<0 || valn>Integer.MAX_VALUE){
				throw new RuntimeException("D39 row counts must fit positive int train / nonnegative int validation sizes");
			}
			if(subnetOut==null || (valn>0 && ((emitDense && outval==null) || subnetValOut==null))){
				throw new RuntimeException("D39 requires subnetout= and, when valn>0, both validation outputs (outval= when emitdense=t, subnetvalout= always)");
			}
			if(aggOut!=null && valn>0 && aggValOut==null){
				throw new RuntimeException("D39 aggregator generation with valn>0 also requires aggvalout=");
			}
			if(perfectFrac!=0 || nearPerfectFrac!=0 || extraPerfectFrac!=0 || extraNearPerfectFrac!=0 || extraUnshreddedPerfectFrac!=0){
				throw new RuntimeException("d39=t replaces legacy spike flags; do not combine them");
			}
		}
		if(referenceCdsIndexFile!=null && binManifestInFile==null){throw new RuntimeException("referencecdsindex= requires binmanifestin=");}
		if(referenceCdsIndexFile!=null && (benchMode || panelTidsFile!=null)){
			throw new RuntimeException("referencecdsindex= is supported only for manifest replay; bench/panel labels are not overridden");
		}
		// This rejection is LEGACY-ONLY (subnetlabels=gene, opt-in, is the point of the combination it
		// forbids here): gene mode's own validation block below requires exactly the flags this would
		// reject, and computes real all-gene subnet targets instead of falling back to legacy semantics.
		if((referenceCdsIndexFile!=null || referenceCdsContributionIndexFile!=null) &&
				(subnetOut!=null || subnetValOut!=null) && subnetLabels==SUBNETLABELS_LEGACY){
			throw new RuntimeException("referencecdsindex=/referencecdscontributionindex= do not support "
				+"subnetout=/subnetvalout= under subnetlabels=legacy"
				+": subnet labels would retain legacy semantics (use subnetlabels=gene)");
		}
		// refcdstable= (opt-in, UMP45's reusable producer, brief 2026-09-09 §4): the SAME all-gene
		// reference-CDS labels as referencecdsindex=, but from a precomputed per-assembly/per-shred
		// table instead of a per-bin sidecar. The legacy referencecdsindex= remains exclusive;
		// refcdstable= and referencecdscontributionindex= may be paired in replay so one shared
		// production manifest can label ordinary/HQ shred rows and whole-genome breakpoint rows,
		// unlike referencecdsindex= this one is not replay-only (see the makeBin/writeReplaySet call
		// sites). binmanifestout= is mandatory whenever it drives GENERATION (no binmanifestin=), per
		// Yoimiya 2026-09-09 (v): an audited run must be able to replay/verify what it labeled.
		final int geneLabelSources=(referenceCdsIndexFile==null ? 0 : 1)+
			(refCdsTableFile==null ? 0 : 1)+(referenceCdsContributionIndexFile==null ? 0 : 1);
		final boolean dualExactSources=(refCdsTableFile!=null && referenceCdsContributionIndexFile!=null);
		if((referenceCdsIndexFile!=null && geneLabelSources>1) || geneLabelSources>2 ||
				(geneLabelSources==2 && !dualExactSources)){
			throw new RuntimeException("referencecdsindex=, refcdstable=, and referencecdscontributionindex= "
				+"may not be combined except for the exact refcdstable= + referencecdscontributionindex= replay pair");
		}
		if(referenceCdsContributionIndexFile!=null){
			if(binManifestInFile==null && !d39Sampling){throw new RuntimeException("referencecdscontributionindex= requires binmanifestin= or d39=t");}
			if(unshreddedCacheFile==null){throw new RuntimeException("referencecdscontributionindex= requires unshreddedcache=");}
			if(minlen!=0){throw new RuntimeException("referencecdscontributionindex= requires minlen=0 so the whole assembly is retained");}
			requireLowerHex64(referenceCdsContributionIndexSha256,
				"referencecdscontributionindexsha256=");
			requireLowerHex64(referenceCdsContributionManifestSha256,
				"referencecdscontributionmanifestsha256=");
			requireLowerHex64(referenceCdsContributionCacheSha256,
				"referencecdscontributioncachesha256=");
			if(benchMode || panelTidsFile!=null){
				throw new RuntimeException("referencecdscontributionindex= is supported only for manifest replay");
			}
		}
		if(refCdsTableFile!=null){
			if(subnetLabels!=SUBNETLABELS_GENE){throw new RuntimeException("refcdstable= requires subnetlabels=gene");}
			// Caller-pinned manifest hash (Yoimiya's 2026-09-09 provenance-gap finding + UMP45's
			// three-argument loader): pinning only the cache/shred-list identity still lets a
			// changed table+manifest silently relabel unchanged observations, so this is required,
			// not optional -- same convention as expectedcopytablesha256= above.
			if(refCdsTableSha256==null || !refCdsTableSha256.matches("[0-9a-f]{64}")){
				throw new RuntimeException("refcdstable= requires refcdstablesha256= (64 hex)");
			}
			if(binManifestInFile==null && binManifestOutFile==null){
				throw new RuntimeException("refcdstable= in generation mode (binmanifestin= unset) requires "
					+"binmanifestout= for an audited run (Yoimiya 2026-09-09)");
			}
			if(benchMode || panelTidsFile!=null){
				throw new RuntimeException("refcdstable= is not supported for bench/panel modes");
			}
			if(extraUnshreddedPerfectFrac>0){
				throw new RuntimeException("refcdstable= does not support extraunshreddedperfectfrac>0 "
					+"(unshredded whole-genome contigs are not shreds in the survival table)");
			}
		}
		// subnetlabels=gene (opt-in, UMP45 2026-09-08 design + 2026-09-09 refcdstable= extension):
		// all-gene reference-CDS labels on the subnet row instead of the legacy native-total regression
		// target, from either label source above. Every required flag is named explicitly in its own
		// message -- see subnet_gene_replay_design_20260908.md "the change" §1.
		if(subnetLabels==SUBNETLABELS_GENE){
			if(geneLabelSources!=1 && !dualExactSources){
				throw new RuntimeException("subnetlabels=gene requires one exact reference-CDS source, or the "
					+"refcdstable= + referencecdscontributionindex= replay pair");
			}
			if(subnetName==null){throw new RuntimeException("subnetlabels=gene requires subnet=ncrna|rrna|trna_anticodon|famset");}
			if(subnetOut==null && subnetValOut==null){throw new RuntimeException("subnetlabels=gene requires subnetout= or subnetvalout=");}
		}
		// A4 explicit-panel emission (PATH_TO_PRODUCTION_v1, 2026-09-02): a deterministic,
		// non-random alternative to the writeSet/writeBench sampling paths, added for the A4
		// end-to-end differential (Java-hits cache vs mmseqs-hits cache, same panel tids, same
		// composite net) -- see writePanel()'s own javadoc for what it does and why the labels
		// are hardcoded. Off by default (panelTidsFile==null): zero effect on any existing path.
		if(panelTidsFile!=null){
			if(panelOutFile==null){throw new RuntimeException("paneltids= requires panelout=");}
			if(unshreddedCacheFile==null){throw new RuntimeException("paneltids= requires unshreddedcache= (the panel's real whole-genome contigs).");}
			if(aggManifestFile==null && aggBundleFile==null){throw new RuntimeException("paneltids= requires aggmanifest= or bundle= (formatAggRow needs the loaded subnets).");}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Per-contig          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * One frozen subnet participating in aggregator-vector emission: its explicit
	 * observation type and subset definition (family ranks for famset, or the
	 * special observation accessor for ncrna/rrna/trna_anticodon),
	 * its loaded CellNet, and an exact-width input buffer (CellNet.applyInput asserts
	 * the width, which doubles as a wiring guard).
	 */
	static final class AggSubnet {
		String name;
		String type;
		int numObs;
		int[] ranks;//null for the special observation types
		ml.CellNet net;
		MagQCNetBundle.Subnet bundleSubnet;
		float[] buf;
		/** Prepared-input serving cache; this AggSubnet is confined to its owning prepared maker thread. */
		private ml.CellNet preparedDenseNet;
		/** Caller-owned outputs from optional shared locked inference. */
		private float[] preparedOutput;
		ml.CellNet netForCurrentThread(){return bundleSubnet==null ? net : bundleSubnet.netForCurrentThread();}
		ml.CellNet netForPreparedThread(){
			if(bundleSubnet==null){return net;}
			if(preparedDenseNet==null){preparedDenseNet=bundleSubnet.newReplayWorker();}
			return preparedDenseNet;
		}
	}

	/** One cached contig's sufficient statistics; family counts stored sparsely. */
	static final class Contig {
		int tid, length, gc, acgt, cds, mapped, coding, r16, r23, r5, rother, trna;
		long glenSum, glenSq;
		int[] famRank, famCount;
		int[] dimer;//16 dinucleotide counts (KmerTracker k=2 native index), cache field 18; null on pre-rebuild caches
		int[] antiRank, antiCount;//sparse tRNA anticodon code(0-63)->count, cache field 17; EMPTY if absent
		String name;//contig_id (cache field 0); populated only in benchmark mode, for the FASTA manifest
		int tableIdx=-1;//refcdstable= only: this shred's row index in ReferenceCdsShredSurvivalTableReader
	}

	/*--------------------------------------------------------------*/
	/*----------------            Loaders           ----------------*/
	/*--------------------------------------------------------------*/

	private void loadAux(){
		if(aggBundleFile!=null){loadConfiguredBundle();}
		if(frozenInputs!=null && !taxonomyFromLabels){
			throw new IllegalArgumentException("Frozen gene subnets require QuickClade labels= input");
		}
		// family list -> count of family columns
		if(familyFile!=null){loadFamilyCount();}
		if(taxonomyFromLabels){
			loadLabelsTaxonomy();
		}else{
			// taxpgm: tid <tab> phylum <tab> pgm -> tid->phylum (staged), phylum vocabulary.
			// phylumIndex needs every distinct phylum name collected first (TreeSet, sorted
			// deterministic order), so tid->phylum is staged in parallel primitive/String
			// lists here and resolved to tid2phylumIdx (int->int) in a second pass below.
			final IntList stagedTid=new IntList();
			final ArrayList<String> stagedPhylum=new ArrayList<String>();
			{
				final ByteFile bf=ByteFile.makeByteFile(taxpgmFile, true);
				final LineParser1 lp=new LineParser1((byte)'\t');
				final TreeSet<String> phyla=new TreeSet<String>();
				for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
					if(line.length==0){continue;}
					lp.set(line);
					if(lp.terms()>=2){
						final int tid=lp.parseInt(0);
						final String phy=lp.parseString(1);
						stagedTid.add(tid); stagedPhylum.add(phy);
						phyla.add(phy);
					}
				}
				bf.close();
				phylumList=new ArrayList<String>(phyla);
				phylumList.add("other");
				for(int i=0; i<phylumList.size(); i++){phylumIndex.put(phylumList.get(i), i);}
				numPhyla=phylumList.size();
			}
			for(int i=0; i<stagedTid.size(); i++){
				final Integer pi=phylumIndex.get(stagedPhylum.get(i));
				tid2phylumIdx.put(stagedTid.get(i), pi==null ? phylumIndex.get("other") : pi);
			}
		}
		loadSizemap();
	}

	/**
	 * Step 4 of the rebuild path (UMP45, 2026-09-03): builds the domain/phylum one-hot sources
	 * from a REAL C1 QuickClade labels.tsv instead of taxpgm's TID/GTDB ground truth --
	 * deployment-faithful, since QuickClade is exactly what production runs. Keyed by the
	 * FILENAME tid parsed from each row's source_rel (labels.tsv column 3), per UMP45's spec --
	 * the training cache's usable tid set does not include the header-tid exceptions this would
	 * otherwise need to reconcile (those tids are excluded from the training cache upstream).
	 *
	 * <p>Three real status values, three different one-hot outcomes -- NOT a two-way
	 * classified/unknown split:
	 * <ul>
	 * <li><b>classified</b>: domain AND phylum both set to their real one-hot columns.</li>
	 * <li><b>partial</b> (a real, expected deployment case -- QuickClade confident on domain,
	 *     genuinely unable to resolve phylum for a novel organism): domain one-hot SET, phylum
	 *     one-hot ALL ZERO. Deliberately NOT the "other" bucket -- "other" means "a real phylum
	 *     value that just isn't in the vocabulary", a different thing from "the server could not
	 *     call the phylum at all" (UMP45's explicit design note, 2026-09-03).</li>
	 * <li><b>unknown</b> (no hit at all): BOTH one-hots all zero.</li>
	 * </ul>
	 *
	 * <p>"No column" is represented as index numPhyla (phylum) / DOMAINS (domain) -- an
	 * out-of-range but NON-NEGATIVE sentinel, stored directly, never routed through the "other"
	 * fallback. Every one-hot write loop in this class already does {@code i==idx ? '1' : '0'}
	 * for i in the real [0,numPhyla)/[0,DOMAINS) column range, so this sentinel naturally
	 * produces an all-zero one-hot with no change needed to any writer. MUST NOT be -1: the
	 * phylum index returned here flows straight out of {@code makeBin}, whose return value does
	 * double duty as the sampling-success/failure signal ({@code targetPhylumIdx<0} means
	 * "reject and retry") -- an earlier version of this method used -1 and silently discarded
	 * every partial-status bin as a failed sampling attempt, producing zero output rows for any
	 * such tid (found and fixed the same session, 2026-09-03).
	 *
	 * <p>A cache tid with NO matching labels.tsv row is a real data problem (a training tid
	 * this pipeline forgot to classify or exclude) and crashes loud -- see {@link #tid2phylumIdx0}
	 * and {@link #domainIdxOf}, which throw rather than silently defaulting.
	 */
	private void loadLabelsTaxonomy(){
		final ByteFile bf=ByteFile.makeByteFile(labelsFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		final TreeSet<String> phyla=new TreeSet<String>();
		final IntList stagedTid=new IntList();
		final ArrayList<String> stagedPhylum=new ArrayList<String>();//null = explicit no-column
		final ArrayList<String> stagedDomain=new ArrayList<String>();//null = explicit no-column
		long rows=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			lp.set(line);
			if(lp.terms()<14){continue;}
			final String sourceRel=lp.parseString(2);
			final String status=lp.parseString(8);
			final String domain=lp.parseString(10);
			final String phylum=lp.parseString(11);
			final int tid=parseTidFromBasename(basenameOf(sourceRel));
			if(tid<0){throw new RuntimeException("labels row source_rel has no parseable filename tid: "+sourceRel);}
			stagedTid.add(tid);
			if(status.equals("classified")){
				stagedPhylum.add(phylum); phyla.add(phylum);
			}else if(status.equals("partial") || status.equals("unknown")){
				stagedPhylum.add(null);
			}else{
				throw new RuntimeException("labels row has an unrecognized status '"+status+"' for "+sourceRel);
			}
			stagedDomain.add(status.equals("unknown") ? null : domain);
			rows++;
		}
		bf.close();
		if(rows==0){throw new RuntimeException("labels file has no data rows: "+labelsFile);}
		phylumList=frozenInputs==null ? new ArrayList<String>(phyla) : new ArrayList<String>(frozenInputs.phyla());
		if(frozenInputs==null){phylumList.add("other");}
		phylumIndex.clear();
		for(int i=0; i<phylumList.size(); i++){phylumIndex.put(phylumList.get(i), i);}
		numPhyla=phylumList.size();
		for(int i=0; i<stagedTid.size(); i++){
			final int tid=stagedTid.get(i);
			final String phy=stagedPhylum.get(i);
			// "No column" MUST be a non-negative sentinel outside the real [0,numPhyla) range,
			// never -1: makeBin's return value (ultimately this same phylum index) does DOUBLE
			// DUTY as the sampling-success/failure signal ("targetPhylumIdx<0" -> reject and
			// retry) -- a real bug found here (Eru, 2026-09-03): using -1 for "valid bin, no
			// phylum column" made every partial-status bin look like a FAILED sampling attempt,
			// silently discarded and retried until the cap, producing zero output rows for any
			// tid with a partial label. numPhyla is always non-negative and never a real column
			// index, so every one-hot write loop (i<numPhyla) naturally still emits all-zero.
			tid2phylumIdxLabels.put(tid, phylumColumn(phy));
			final String dom=stagedDomain.get(i);
			// Same reasoning applies to domain, for consistency and future-proofing, even though
			// no current caller treats domainIdxOf's return as a rejection signal today.
			tid2domainIdxLabels.put(tid, domainColumn(dom));
		}
		System.err.println("labels taxonomy: "+rows+" rows, "+numPhyla+" phyla (incl. other), from "+labelsFile);
	}

	private static String basenameOf(String path){
		final int i=path.lastIndexOf('/');
		return i<0 ? path : path.substring(i+1);
	}

	/** Shared family cardinality; the configured bundle independently verifies the exact family file. */
	private void loadFamilyCount(){
		final ByteFile input=ByteFile.makeByteFile(familyFile, true);
		int count=0;
		try{
			input.nextLine();
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				if(line.length>0){count++;}
			}
		}finally{if(input.close()){throw new RuntimeException("I/O error reading family list "+familyFile);}}
		numFam=count;
	}

	/** Null means QuickClade could not call a phylum; a real unseen name means the other column. */
	private int phylumColumn(String phylum){
		if(phylum==null){return numPhyla;}
		final Integer column=phylumIndex.get(phylum);
		return column==null ? phylumIndex.get("other") : column;
	}

	/** Null is the unknown-status all-zero domain block, distinct from a named other domain. */
	private static int domainColumn(String domain){return domain==null ? DOMAINS : domainIndex(domain);}

	/** Parses a leading tid_NNN_ or tid|NNN| prefix, matching every bash-side script's own
	 *  extract_tid_from_name() convention exactly (quickclade_label_prep_v1.sh and friends). */
	private static int parseTidFromBasename(String base){
		if(base.startsWith("tid_")){
			final int end=base.indexOf('_', 4);
			if(end>4){try{return Integer.parseInt(base.substring(4, end));}catch(NumberFormatException e){return -1;}}
		}else if(base.startsWith("tid|")){
			final int end=base.indexOf('|', 4);
			if(end>4){try{return Integer.parseInt(base.substring(4, end));}catch(NumberFormatException e){return -1;}}
		}
		return -1;
	}

	/** sizemap: tid <tab> bp. Last row for a tid wins (matches HashMap.put's overwrite
	 *  semantics) - IntLongHashMap.put() does NOT overwrite, so remove-then-put. */
	private void loadSizemap(){
		{
			final ByteFile bf=ByteFile.makeByteFile(sizemapFile, true);
			final LineParser1 lp=new LineParser1((byte)'\t');
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0){continue;}
				lp.set(line);
				if(lp.terms()>=2){
					final int tid=lp.parseInt(0);
					final long bp=lp.parseLong(1);
					genomeSize.remove(tid);
					genomeSize.put(tid, bp);
				}
			}
			bf.close();
		}
	}

	/**
	 * Loads the per-contig cache into per-tid contig lists (dense tid->index +
	 * parallel ArrayList<Contig>[]).  Package visibility is intentional: the
	 * FASTA-in differential harness uses this exact training-path parser rather
	 * than maintaining a second cache reader.  Normal production execution
	 * still reaches it only through {@link #process()}.
	 */
	void loadCache(){
		final ByteFile bf=ByteFile.makeByteFile(cacheFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		final IntList rankBuf=new IntList(), countBuf=new IntList();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			lp.set(line);
			final int terms=lp.terms();
			assert(terms>=17) : "Malformed cache row: "+terms+" fields (need >=17): "+new String(line);
			final Contig c=new Contig();
			c.tid=lp.parseInt(1);
			c.length=lp.parseInt(3);
			if(c.length<minlen){continue;}
			if(benchMode || binManifestOutFile!=null || binManifestInFile!=null){c.name=lp.parseString(0);}//retain IDs for benchmark/shared selection manifests
			if(refCdsTable!=null){
				assert(c.name!=null) : "refcdstable= requires shred names, guaranteed by constructor validation "
					+"(binmanifestin=/binmanifestout= required whenever refcdstable= is set)";
				c.tableIdx=refCdsTable.requireShred(c.name);
			}
			c.gc=lp.parseInt(4);
			c.acgt=lp.parseInt(5);
			c.cds=lp.parseInt(6);
			c.mapped=lp.parseInt(7);
			c.glenSum=lp.parseLong(8);
			c.glenSq=lp.parseLong(9);
			c.coding=lp.parseInt(10);
			c.r16=lp.parseInt(11);
			c.r23=lp.parseInt(12);
			c.r5=lp.parseInt(13);
			c.rother=lp.parseInt(14);
			c.trna=lp.parseInt(15);
			final int flen=lp.length(16);
			if(flen==0){c.famRank=EMPTY; c.famCount=EMPTY;}
			else{
				rankBuf.clear(); countBuf.clear();
				parseFamCounts(lp.line(), lp.a(), lp.b(), rankBuf, countBuf);
				c.famRank=rankBuf.toArray();
				c.famCount=countBuf.toArray();
			}
			// Rebuild fields (backward-compatible): field 17 = tRNA anticodon sparse code:count; field 18 = 16
			// dinucleotide counts (dense CSV). Absent on pre-rebuild 17-field caches -> anti EMPTY, dimer null.
			if(terms>=18){
				final int alen=lp.length(17);
				if(alen==0){c.antiRank=EMPTY; c.antiCount=EMPTY;}
				else{
					rankBuf.clear(); countBuf.clear();
					parseFamCounts(lp.line(), lp.a(), lp.b(), rankBuf, countBuf);
					c.antiRank=rankBuf.toArray();
					c.antiCount=countBuf.toArray();
				}
			}else{c.antiRank=EMPTY; c.antiCount=EMPTY;}
			if(terms>=19){
				lp.length(18);
				c.dimer=parseDimers(lp.line(), lp.a(), lp.b());
			}
			int idx=tidToIdx.get(c.tid);
			if(idx<0){
				idx=contigLists.size();
				tidToIdx.put(c.tid, idx);
				contigLists.add(new ArrayList<Contig>());
				if(idx>=domainIdxArr.length){domainIdxArr=Arrays.copyOf(domainIdxArr, Math.max(idx+1, domainIdxArr.length*2));}
				domainIdxArr[idx]=domainIndex(lp.parseString(2));
			}
			contigLists.get(idx).add(c);
		}
		bf.close();
	}

	/**
	 * Reconstructs the family-count component for one loaded organism using the
	 * same {@link Agg#add(Contig)} accumulation used by training-bin creation.
	 * This is a package-visible differential seam, not a second vectorization
	 * implementation; the caller must have invoked {@link #loadCache()} first.
	 * @param tid Organism id whose loaded cache rows should be summed.
	 * @param familyCount Number of family-rank columns to return.
	 * @return Dense rank-indexed family copy counts.
	 */
	int[] familyCountsForTid(final int tid, final int familyCount){
		if(familyCount<0){throw new RuntimeException("Negative family count: "+familyCount);}
		final ArrayList<Contig> contigs=getContigs(tid);
		if(contigs==null){throw new RuntimeException("No loaded cache rows for tid "+tid);}
		final int[] fam=new int[familyCount];
		final Agg agg=new Agg();
		agg.setFam(fam);
		for(final Contig c : contigs){
			for(final int rank : c.famRank){
				if(rank<0 || rank>=familyCount){
					throw new RuntimeException("Cache family rank "+rank+
						" outside requested vector width "+familyCount+" for tid "+tid);
				}
			}
			agg.add(c);
		}
		return fam;
	}

	/** Parses "rank:count;rank:count;..." within [a,b) of line into the (cleared) output
	 *  lists, in order. Zero allocation beyond the two IntList's own backing-array growth. */
	private static void parseFamCounts(byte[] line, int a, int b, IntList outRank, IntList outCount){
		int start=a;
		for(int i=a; i<=b; i++){
			if(i==b || line[i]==';'){
				if(i>start){
					int colon=-1;
					for(int j=start; j<i; j++){if(line[j]==':'){colon=j; break;}}
					assert(colon>start) : "Malformed famcounts field (missing ':' in a rank:count pair).";
					outRank.add(Parse.parseInt(line, start, colon));
					outCount.add(Parse.parseInt(line, colon+1, i));
				}
				start=i+1;
			}
		}
	}

	/** Parses exactly 16 comma-separated ints within [a,b) into a new int[16] (dense dinucleotide
	 *  counts, KmerTracker native index order). Zero-alloc beyond the returned array. */
	private static int[] parseDimers(byte[] line, int a, int b){
		final int[] out=new int[16];
		int idx=0, start=a;
		for(int i=a; i<=b; i++){
			if(i==b || line[i]==','){
				if(i>start && idx<16){out[idx]=Parse.parseInt(line, start, i);}
				idx++;
				start=i+1;
			}
		}
		assert(idx==16) : "dimer field expected 16 comma-separated counts, got "+idx+": "+new String(line, a, b-a);
		return out;
	}

	/** Returns tid's contig list, or null if tid was never seen in the cache (mirrors
	 *  the original {@code byTid.get(tid)}'s null-on-absent contract exactly). */
	private ArrayList<Contig> getContigs(int tid){
		final int idx=tidToIdx.get(tid);
		return idx<0 ? null : contigLists.get(idx);
	}

	/** Loads the SECOND (unshredded-genome) per-contig cache into a parallel tid->contig-list
	 *  map, for the unshredded-isolate-spike append path (Brian 2026-08-26: a "perfect" bin
	 *  built from real, whole-genome CallGenes-called contigs instead of 20kb shred fragments --
	 *  zero shredding artifacts). Same 17-19 field format as the primary cache (produced by the
	 *  SAME unmodified CacheBuilder, just fed whole-genome GFFs/proteins/sequences instead of
	 *  shred ones) -- this is a deliberate parse duplication of loadCache(), not a shared helper,
	 *  because it targets different fields (tidToIdx2/contigLists2) and does not need domain
	 *  resolution (the organism's domain is already known via the primary cache's tidToIdx/
	 *  domainIdxArr for the same tid; nothing in makeBin's unshredded path reads domain). */
	private void loadCache2(){
		final ByteFile bf=ByteFile.makeByteFile(unshreddedCacheFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		final IntList rankBuf=new IntList(), countBuf=new IntList();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			lp.set(line);
			final int terms=lp.terms();
			assert(terms>=17) : "Malformed unshredded cache row: "+terms+" fields (need >=17): "+new String(line);
			final Contig c=new Contig();
			//Same benchMode-gated name retention as loadCache() (Eru, 2026-08-29): loadCache2's
			//own doc comment says it targets different fields and skips domain resolution, but it
			//missed this one -- every isolate bin's manifest rows came out with contig name
			//literally "null" (3171 rows across 50 bins, caught inspecting bench_v4b_manifest.tsv),
			//since Contig.name defaults null and nothing here ever set it. BenchmarkBinWriter would
			//crash-loud on the first such lookup (no shred is named "null"), so this was a real,
			//not-yet-hit blocker for scoring the isolate regime.
			if(benchMode || binManifestOutFile!=null || binManifestInFile!=null){c.name=lp.parseString(0);}
			c.tid=lp.parseInt(1);
			c.length=lp.parseInt(3);
			if(c.length<minlen){continue;}
			c.gc=lp.parseInt(4);
			c.acgt=lp.parseInt(5);
			c.cds=lp.parseInt(6);
			c.mapped=lp.parseInt(7);
			c.glenSum=lp.parseLong(8);
			c.glenSq=lp.parseLong(9);
			c.coding=lp.parseInt(10);
			c.r16=lp.parseInt(11);
			c.r23=lp.parseInt(12);
			c.r5=lp.parseInt(13);
			c.rother=lp.parseInt(14);
			c.trna=lp.parseInt(15);
			final int flen=lp.length(16);
			if(flen==0){c.famRank=EMPTY; c.famCount=EMPTY;}
			else{
				rankBuf.clear(); countBuf.clear();
				parseFamCounts(lp.line(), lp.a(), lp.b(), rankBuf, countBuf);
				c.famRank=rankBuf.toArray();
				c.famCount=countBuf.toArray();
			}
			if(terms>=18){
				final int alen=lp.length(17);
				if(alen==0){c.antiRank=EMPTY; c.antiCount=EMPTY;}
				else{
					rankBuf.clear(); countBuf.clear();
					parseFamCounts(lp.line(), lp.a(), lp.b(), rankBuf, countBuf);
					c.antiRank=rankBuf.toArray();
					c.antiCount=countBuf.toArray();
				}
			}else{c.antiRank=EMPTY; c.antiCount=EMPTY;}
			if(terms>=19){
				lp.length(18);
				c.dimer=parseDimers(lp.line(), lp.a(), lp.b());
			}
			int idx=tidToIdx2.get(c.tid);
			if(idx<0){
				idx=contigLists2.size();
				tidToIdx2.put(c.tid, idx);
				contigLists2.add(new ArrayList<Contig>());
			}
			contigLists2.get(idx).add(c);
		}
		bf.close();
	}

	/** Returns tid's UNSHREDDED (real whole-genome) contig list, or null if tid was never seen
	 *  in unshreddedCacheFile. Mirrors getContigs exactly, against the second cache. */
	private ArrayList<Contig> getContigsUnshredded(int tid){
		final int idx=tidToIdx2.get(tid);
		return idx<0 ? null : contigLists2.get(idx);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Process            ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads a family-feature subset (one rank per line) for reduced-width vectors. */
	private int[] loadRanks(String file){
		final IntList l=new IntList();
		final ByteFile bf=ByteFile.makeByteFile(file, true);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			int a=0, b=line.length;
			while(a<b && line[a]<=' '){a++;}
			while(b>a && line[b-1]<=' '){b--;}
			if(b>a){l.add(Parse.parseInt(line, a, b));}
		}
		bf.close();
		final int[] a=l.toArray();
		System.err.println("feature subset: "+a.length+" family ranks kept");
		return a;
	}

	private static boolean parseAggObs(String s){
		s=s.toLowerCase();
		if(s.equals("serve") || s.equals("whole") || s.equals("wholebin")){return true;}
		if(s.equals("clean") || s.equals("target") || s.equals("targetonly")){return false;}
		throw new RuntimeException("Unknown aggobs="+s+" (serve|clean)");
	}

	private static int parsePoolMode(String s){
		s=s.toLowerCase();
		if(s.equals("trainval")){return POOL_TRAINVAL;}
		if(s.equals("valsplit")){return POOL_VALSPLIT;}
		if(s.equals("allbutc")){return POOL_ALLBUTC;}
		throw new RuntimeException("Unknown poolmode="+s+" (trainval|valsplit|allbutc)");
	}

	/** subnetlabels=legacy (default): native-only obs, single native-total target (unchanged).
	 *  subnetlabels=gene (opt-in): whole-bin obs, two all-gene reference-CDS targets -- see
	 *  {@link #formatFamsetGeneRow}/{@link #formatNcrnaGeneRow} and subnet_gene_replay_design_20260908.md. */
	private static int parseSubnetLabels(String s){
		s=s.toLowerCase();
		if(s.equals("legacy")){return SUBNETLABELS_LEGACY;}
		if(s.equals("gene")){return SUBNETLABELS_GENE;}
		throw new RuntimeException("Unknown subnetlabels="+s+" (legacy|gene)");
	}

	/** subnetfeatures=legacy (default): byte-identical to today (4 learned-transform cols + 1 pooled
	 *  ratio per subnet). subnetfeatures=six (opt-in, UMP45 design 2026-09-09 revision2): 4 learned
	 *  outputs of a four-output subnet + 2 table-derived expected-copy ratios, no pooled column. */
	private static int parseSubnetFeatures(String s){
		s=s.toLowerCase();
		if(s.equals("legacy")){return SUBNETFEATURES_LEGACY;}
		if(s.equals("six")){return SUBNETFEATURES_SIX;}
		throw new RuntimeException("Unknown subnetfeatures="+s+" (legacy|six)");
	}

	/**
	 * Loads the aggregator manifest: tab-separated rows
	 * base, numobs, numIn, subset_path, val_path, net_path (val_path unused here),
	 * optionally type. The original six-column form remains valid; the seventh type
	 * column is required for new rrna and trna_anticodon rows. Special rows use
	 * subset_path "-" and are observed from Agg's frozen special counters. Every net's input width is asserted against
	 * obs+phylum+SHARED_CONTEXT_WIDTH - the FROZEN 2026-08-24 vector-layout rebuild
	 * retired the old per-net representation flags (sndomain/snhhcaga/sngenelen/
	 * snbinscaled/sncodingaffine): every subnet now gets the identical standard context,
	 * so there is nothing left to override per row. #-lines are comments.
	 */
	private void loadAggManifest(String file){
		aggSubnets=new ArrayList<AggSubnet>();
		final ByteFile bf=ByteFile.makeByteFile(file, true);
		for(byte[] lineB=bf.nextLine(); lineB!=null; lineB=bf.nextLine()){
			String line=new String(lineB).trim();
			if(line.length()==0 || line.charAt(0)=='#'){continue;}
			String[] p=line.split("\t",-1);
			if(p.length!=6 && p.length!=7){throw new RuntimeException("Manifest row needs 6 or 7 columns "
				+"(base numobs numIn subset_path val_path net_path [type]): "+line);}
			AggSubnet s=new AggSubnet();
			s.name=p[0];
			if(p.length==6 && ("rrna".equals(s.name) || "trna_anticodon".equals(s.name))){
				throw new RuntimeException("Special manifest row "+s.name+" requires explicit seventh type column");
			}
			s.type=p.length==7 ? p[6] : ("ncrna".equals(s.name) ? "ncrna" : "famset");
			noteAggObservationType(s.type);
			s.numObs=Integer.parseInt(p[1]);
			final int specialCount=MagQCObservationLayout.specialCount(s.type);
			if(specialCount<0){
				s.ranks=loadRanks(p[3]);
				if(s.ranks.length!=s.numObs){throw new RuntimeException("Subset "+s.name
					+": "+s.ranks.length+" ranks but numobs="+s.numObs);}
			}else{
				if(!"-".equals(p[3])){throw new RuntimeException(s.type+" row requires subset_path=-: "+s.name);}
				if(s.numObs!=specialCount){throw new RuntimeException(s.type+" row numobs must be "+specialCount
					+": "+s.name);}
			}
			s.net=ml.CellNetParser.load(p[5]);
			if(s.net==null){throw new RuntimeException("Failed to load net "+p[5]);}
			final int expected=s.numObs+numPhyla+SHARED_CONTEXT_WIDTH;
			final int manifestIn=Integer.parseInt(p[2]);
			if(manifestIn>=0 && manifestIn!=expected){throw new RuntimeException("Subset "+s.name
				+": manifest numIn="+manifestIn+" but expected "+expected);}
			if(s.net.numInputs()!=expected){throw new RuntimeException("Subset "+s.name
				+": net takes "+s.net.numInputs()+" inputs but expected "+expected
				+" (this net predates the FROZEN vector-layout rebuild and must be retrained)");}
			s.buf=new float[expected];
			aggSubnets.add(s);
		}
		bf.close();
		if(aggSubnets.isEmpty()){throw new RuntimeException("Empty aggregator manifest: "+file);}
		System.err.println("aggregator manifest: "+aggSubnets.size()+" subnets loaded");
	}

	/** Loads the self-contained release bundle into the same AggSubnet shape used by loose nets.
	 *  Dispatches on subnetFeatures: legacy (default, byte-identical to before this flag existed)
	 *  or six (opt-in, requires every subnet to declare four outputs -- see loadAggBundleSix). */
	private void loadAggBundle(String file){
		try{
			if(familyFile==null){throw new RuntimeException("bundle= requires familylist= for the release hash gate");}
			if(subnetFeatures==SUBNETFEATURES_SIX){loadAggBundleSix(file);}
			else{loadAggBundleLegacy(file);}
		}catch(Exception e){if(e instanceof RuntimeException){throw (RuntimeException)e;} throw new RuntimeException("Failed to load bundle "+file,e);}
	}

	/** Loads once, before taxonomy columns are constructed, and validates independent release inputs. */
	private void loadConfiguredBundle(){
		if(loadedBundle!=null){return;}
		assert(aggBundleFile!=null) : "Bundle initialization requires the configured bundle path";
		if(familyFile==null){throw new IllegalArgumentException("bundle= requires familylist=");}
		try{
			final java.nio.file.Path path=java.nio.file.Paths.get(aggBundleFile);
			final MagQCNetBundle candidate=subnetFeatures==SUBNETFEATURES_SIX ?
					MagQCNetBundle.loadMultiOutput(path) : MagQCNetBundle.load(path);
			candidate.requireConsumerRelease(releaseManifestFile==null ? null : java.nio.file.Paths.get(releaseManifestFile),
					releaseManifestSha80);
			candidate.requireFamilyList(java.nio.file.Paths.get(familyFile));
			final MagQCNetBundle.FrozenInputs inputs=candidate.metadata("input_definition")==null ? null : candidate.requireFrozenInputs();
			if(inputs!=null && subnetFeatures!=SUBNETFEATURES_SIX){
				throw new IllegalArgumentException("Frozen gene subnets require subnetfeatures=six");
			}
			loadedBundle=candidate;
			frozenInputs=inputs;
			if(verbose){System.err.println("replay SIMD: Shared.SIMD="+shared.Shared.SIMD
				+" SIMD_FMA="+shared.Shared.SIMD_FMA
				+" SIMD_FEED_FORWARD="+shared.Shared.SIMD_FEED_FORWARD
				+" Vector.SIMD_FMA_SPARSE="+simd.Vector.SIMD_FMA_SPARSE);}
		}catch(Exception e){throw new IllegalArgumentException("Invalid configured subnet bundle",e);}
	}

	/** Records whether an aggregator source needs structural anticodon snapshots. */
	private void noteAggObservationType(String type){
		if("trna_anticodon".equals(type)){aggNeedsAnticodons=true;}
	}

	private void loadAggBundleLegacy(String file) throws Exception {
		loadConfiguredBundle();
		final MagQCNetBundle bundle=loadedBundle;
		bundle.requireFamilyList(java.nio.file.Paths.get(familyFile));
		aggSubnets=new ArrayList<AggSubnet>(bundle.size());
		for(int i=0; i<bundle.size(); i++){
			final MagQCNetBundle.Subnet bs=bundle.subnet(i);
			final int expected=bs.numObs+numPhyla+SHARED_CONTEXT_WIDTH;
			if(bs.order!=i){throw new RuntimeException("bundle ordering mismatch at "+i+" for "+bs.id);}
			if(bs.expectedInputs!=expected){throw new RuntimeException("Subset "+bs.id+": bundle expectedInputs="
				+bs.expectedInputs+" but vectorizer computes "+expected);}
			if(bs.expectedOutputs!=1){throw new RuntimeException("Subset "+bs.id+": expected one output");}
			noteAggObservationType(bs.type);
			final int specialCount=MagQCObservationLayout.specialCount(bs.type);
			if(specialCount<0){
				if(bs.familyRanks==null || bs.familyRanks.length!=bs.numObs){throw new RuntimeException("Subset "+bs.id
					+": bundle rank count does not match numObs="+bs.numObs);}
				for(int r:bs.familyRanks){if(r<0 || r>=numFam){throw new RuntimeException("Subset "+bs.id
						+": rank "+r+" outside vectorizer family list [0,"+numFam+")");}}
			}else{
				if(bs.familyRanks!=null){throw new RuntimeException("Special subset "+bs.id+" must not carry family ranks");}
				if(bs.numObs!=specialCount){throw new RuntimeException("Subset "+bs.id+" type "+bs.type
					+" requires numObs="+specialCount+", got "+bs.numObs);}
			}
			final AggSubnet s=new AggSubnet(); s.name=bs.id; s.type=bs.type; s.numObs=bs.numObs;
			s.ranks=bs.familyRanks==null ? null : bs.familyRanks.clone(); s.bundleSubnet=bs; s.buf=new float[expected];
			aggSubnets.add(s);
		}
		if(aggSubnets.size()!=bundle.size()){throw new RuntimeException("bundle subnet count changed during load");}
		subnetBlockWidth=4; pooledCols=1;
		System.err.println("aggregator bundle (legacy): "+aggSubnets.size()+" subnets loaded");
	}

	/** subnetfeatures=six loader (design revision2 §4.1): loadMultiOutput (accepts schema 1 or 2),
	 *  every subnet must declare four outputs (bounded-fixture constraint, not a topology selection),
	 *  names/units validated via the shared {@link #requireFourOutputContract} helper (also used by
	 *  MagQCTool's phase-0 check), and the expected-copy binding compiled once here. */
	private void loadAggBundleSix(String file) throws Exception {
		loadConfiguredBundle();
		final MagQCNetBundle bundle=loadedBundle;
		bundle.requireFamilyList(java.nio.file.Paths.get(familyFile));
		if(frozenInputs!=null){frozenInputs.requireVocabulary(phylumList);}
		aggSubnets=new ArrayList<AggSubnet>(bundle.size());
		for(int i=0; i<bundle.size(); i++){
			final MagQCNetBundle.Subnet bs=bundle.subnet(i);
			final int expected=bs.numObs+numPhyla+SHARED_CONTEXT_WIDTH;
			if(bs.order!=i){throw new RuntimeException("bundle ordering mismatch at "+i+" for "+bs.id);}
			if(bs.expectedInputs!=expected){throw new RuntimeException("Subset "+bs.id+": bundle expectedInputs="
				+bs.expectedInputs+" but vectorizer computes "+expected);}
			if(bs.expectedOutputs!=4){throw new RuntimeException("subnetfeatures=six requires four-output subnets; "
				+bs.id+" declares "+bs.expectedOutputs);}
			noteAggObservationType(bs.type);
			final int specialCount=MagQCObservationLayout.specialCount(bs.type);
			if(specialCount<0){
				if(bs.familyRanks==null || bs.familyRanks.length!=bs.numObs){throw new RuntimeException("Subset "+bs.id
					+": bundle rank count does not match numObs="+bs.numObs);}
				for(int r:bs.familyRanks){if(r<0 || r>=numFam){throw new RuntimeException("Subset "+bs.id
						+": rank "+r+" outside vectorizer family list [0,"+numFam+")");}}
			}else{
				if(bs.familyRanks!=null){throw new RuntimeException("Special subset "+bs.id+" must not carry family ranks");}
				if(bs.numObs!=specialCount){throw new RuntimeException("Subset "+bs.id+" type "+bs.type
					+" requires numObs="+specialCount+", got "+bs.numObs);}
			}
			final AggSubnet s=new AggSubnet(); s.name=bs.id; s.type=bs.type; s.numObs=bs.numObs;
			s.ranks=bs.familyRanks==null ? null : bs.familyRanks.clone(); s.bundleSubnet=bs; s.buf=new float[expected];
			aggSubnets.add(s);
		}
		if(aggSubnets.size()!=bundle.size()){throw new RuntimeException("bundle subnet count changed during load");}
		final String[][] contract=requireFourOutputContract(bundle);
		sixNames=contract[0]; sixUnits=contract[1];
		expectedCopy=MagQCExpectedCopyBinding.bind(expectedCopyTableFile,expectedCopyTableSha256,bundle,
			subnetPopulationsFile,subnetPopulationsSha256,numFam);
		observedByItem=new int[expectedCopy.itemWidth()];
		if(observedByItem.length!=numFam+NCRNA_OBS && observedByItem.length!=numFam+NCRNA_OBS+TRNA_ANTICODON_OBS){
			throw new RuntimeException("Unsupported expected-copy item width "+observedByItem.length
				+"; expected "+(numFam+NCRNA_OBS)+" or "+(numFam+NCRNA_OBS+TRNA_ANTICODON_OBS));
		}
		if(observedByItem.length==numFam+NCRNA_OBS){
			for(AggSubnet s:aggSubnets){if("trna_anticodon".equals(s.type)){
				throw new RuntimeException("trna_anticodon requires extended expected-copy items (64 structural + unknown)");
			}}
		}
		subnetBlockWidth=6; pooledCols=0;
		if(verbose){System.err.println("aggregator bundle (six): "+aggSubnets.size()+" subnets loaded, "
			+"expected-copy table="+expectedCopyTableFile+" declaration="+subnetPopulationsFile);}
	}

	/** Validates the four-output contract shared by MVM's six-feature loader and MagQCTool's phase-0
	 *  check (design revision2 §4.1 step 3/§4.3 step 2): every subnet in the bundle must declare
	 *  exactly 4 outputs, and every subnet's validated outputNames()/outputUnits() (already checked
	 *  against the bundle's own frozen private constants at load time, :526-527) must equal subnet
	 *  0's exactly -- MVM never retypes the contract, it only requires internal agreement.
	 *  @return {sixNames, sixUnits}: the four learned names/units followed by the two
	 *          {@link MagQCExpectedCopyFeatures#featureNames()} derived-feature names/units. */
	static String[][] requireFourOutputContract(MagQCNetBundle bundle){
		if(bundle.size()==0){throw new RuntimeException("subnetfeatures=six: empty bundle");}
		final MagQCNetBundle.Subnet s0=bundle.subnet(0);
		if(s0.expectedOutputs!=4){throw new RuntimeException("subnetfeatures=six requires four-output subnets; "
			+s0.id+" declares "+s0.expectedOutputs);}
		final String[] names0=s0.outputNames(), units0=s0.outputUnits();
		if(names0==null || names0.length!=4 || units0==null || units0.length!=4){
			throw new RuntimeException("subnetfeatures=six: subnet 0 output contract malformed");
		}
		for(int i=0; i<bundle.size(); i++){
			final MagQCNetBundle.Subnet bs=bundle.subnet(i);
			if(bs.expectedOutputs!=4){throw new RuntimeException("subnetfeatures=six requires four-output subnets; "
				+bs.id+" declares "+bs.expectedOutputs);}
			if(!Arrays.equals(bs.outputNames(),names0) || !Arrays.equals(bs.outputUnits(),units0)){
				throw new RuntimeException("subnetfeatures=six: subnet "+bs.id+" output names/units differ from subnet 0");
			}
		}
		final String[] fn=MagQCExpectedCopyFeatures.featureNames();
		final String[] sixNames=new String[names0.length+fn.length];
		System.arraycopy(names0, 0, sixNames, 0, names0.length);
		System.arraycopy(fn, 0, sixNames, names0.length, fn.length);
		final String[] sixUnits=new String[units0.length+fn.length];
		System.arraycopy(units0, 0, sixUnits, 0, units0.length);
		for(int i=0; i<fn.length; i++){sixUnits[units0.length+i]=MagQCExpectedCopyFeatures.FEATURE_UNIT;}
		return new String[][]{sixNames,sixUnits};
	}

	/**
	 * Selects the raw dense head: the K most-prevalent families (organism presence
	 * count over all usable orgs, ties broken by rank), a reference-DB constant like
	 * avgCopyWhenPresent. Their whole-bin counts feed the aggregator raw (enc=two
	 * presence+excess), giving it a direct low-missingness signal alongside the
	 * subnet summaries.
	 */
	private int[] computeDenseHead(IntList usable, int k){
		if(frozenInputs!=null){
			if(k<numFam){throw new IllegalArgumentException("Frozen composite inputs require densehead>=family count "+numFam);}
			final int[] ranks=new int[numFam];
			for(int i=0; i<numFam; i++){ranks[i]=i;}
			return ranks;
		}
		final int[] present=new int[numFam];
		final int[] orgCount=new int[numFam];
		for(int i=0; i<usable.size(); i++){
			Arrays.fill(orgCount, 0);
			for(Contig c : getContigs(usable.get(i))){
				for(int j=0; j<c.famRank.length; j++){orgCount[c.famRank[j]]+=c.famCount[j];}
			}
			for(int f=0; f<numFam; f++){if(orgCount[f]>0){present[f]++;}}
		}
		Integer[] order=new Integer[numFam];
		for(int i=0; i<numFam; i++){order[i]=i;}
		Arrays.sort(order, (a, b) -> (present[a]!=present[b] ? present[b]-present[a] : a-b));
		final int[] head=new int[Math.min(k, numFam)];
		for(int i=0; i<head.length; i++){head[i]=order[i];}
		return head;
	}

	void process(){
		loadAux();
		loadExcludedTids();
		if(referenceCdsContributionIndexFile!=null){
			final String observedCacheSha256=MagQCBinManifest.sha256File(unshreddedCacheFile);
			if(!referenceCdsContributionCacheSha256.equals(observedCacheSha256)){
				throw new RuntimeException("reference-CDS contribution cache SHA-256 mismatch: expected "+
					referenceCdsContributionCacheSha256+" observed "+observedCacheSha256+
					" for "+unshreddedCacheFile);
			}
			try{
				referenceCdsContributionIndex=ReferenceCdsGeneContributionIndex.Reader.load(
					referenceCdsContributionIndexFile,referenceCdsContributionIndexSha256,
					referenceCdsContributionManifestSha256,referenceCdsContributionCacheSha256,numFam);
			}catch(java.io.IOException e){
				throw new RuntimeException("Could not load reference-CDS contribution index "+
					referenceCdsContributionIndexFile,e);
			}
		}
		// Loaded BEFORE loadCache() so loadCache/loadCache2 can resolve every shred's table row
		// (crash-loud on any cache shred absent from the table) while parsing, not in a second pass.
		// Three-argument loader (UMP45 + Yoimiya's provenance-gap fix, 2026-09-09): the manifest
		// itself is hash-pinned by refcdstablesha256=, and expectedShredsSha256 binds every table's
		// declared shred list to THIS run's actual cache -- MagQCBinManifest.sha256File is the same
		// hashing convention the bin-manifest replay path already trusts.
		if(refCdsTableFile!=null){
			refCdsTable=ReferenceCdsShredSurvivalTableReader.load(refCdsTableFile, refCdsTableSha256,
				MagQCBinManifest.sha256File(cacheFile));
		}
		loadCache();

		// usable tids: present in cache + sizemap, with contigs
		final IntList usable=new IntList();
		final int[] allTids=tidToIdx.toArray();
		for(int tid : allTids){
			if(genomeSize.contains(tid) && !excludedTids.contains(tid) && !getContigs(tid).isEmpty()){usable.add(tid);}
		}
		usable.sort();
		if(d39Sampling && usable.size()<2){throw new RuntimeException("D39 needs at least two eligible organisms");}
		if(samplePhylum!=null){loadSamplingPhylum(usable);}
		System.err.println("usable orgs="+usable.size()+", numFam="+numFam+", numPhyla="+numPhyla);

		// avg~0.5 normalization scales (Brian 2026-08-24): corpus reference for the size/#genes/gene-length
		// context features. Computed now (layout-independent); WIRED INTO the emit methods in the one-pass
		// vector-layout rewrite. Currently compute-and-log only, so all existing output is byte-identical.
		computeNormScales(usable);

		// optional taxonomy for same-family contaminant bias
		if(treeFile!=null){
			tree=TaxTree.loadTaxTree(treeFile, System.err, false, false);
			for(int i=0; i<usable.size(); i++){
				final int tid=usable.get(i);
				TaxNode fn=(tree==null ? null : tree.getNodeAtLevel(tid, TaxTree.FAMILY));
				tid2family.put(tid, fn==null ? -1 : fn.id);
			}
		}

		// global organism split BEFORE sampling. shuffleInPlace replicates
		// java.util.Collections.shuffle(List,Random)'s EXACT algorithm/call sequence -
		// IntList's own shuffle() draws from Shared.threadLocalRandom(), NOT this seeded
		// stream, and would silently desync the whole train/val split from the original.
		// splitseed= (Eru 2026-09-04, plan v2 E1a): default unset -> uses seed, byte-identical
		// to pre-splitseed behavior. When set, decouples the holdout draw from the sample draw
		// (line ~924-926, seed*2+1/seed*2+2) so shard k can run seed=k splitseed=1 for an
		// IDENTICAL holdout C across all shards while sampling different rows per shard.
		IntList trainTids, valTids;
		if(noOrganismHoldout){
			// Opt-in redesign path: validation is a random ROW sample from the same organism pool,
			// never an organism holdout.  Keep the legacy split untouched when this flag is absent.
			trainTids=subList(usable, 0, usable.size());
			valTids=subList(usable, 0, usable.size());
			System.err.println("noorganismholdout=true: train and validation draw from all "+usable.size()+" organisms");
		}else{
			final Random split=new Random(splitSeedSet ? splitSeed : seed);
			shuffleInPlace(usable, split);
			final int nVal=(int)Math.round(usable.size()*valfrac);
			valTids=subList(usable, 0, nVal);
			trainTids=subList(usable, nVal, usable.size());
		}
		System.err.println("train orgs="+trainTids.size()+", val orgs="+valTids.size());
		{//val tid list (sorted, cheap at this scale): lets a test/operator verify holdout identity
			//across seed/splitseed combinations without any behavior change to file output.
			IntList sortedVal=new IntList(valTids.size());
			for(int i=0; i<valTids.size(); i++){sortedVal.add(valTids.get(i));}
			sortedVal.sort();
			System.err.println("val tids="+sortedVal);
		}

		// precompute recoverable bp per tid
		for(int i=0; i<usable.size(); i++){
			final int tid=usable.get(i);
			long sum=0; for(Contig c : getContigs(tid)){sum+=c.length;}
			recoverable.put(tid, sum);
		}

		// unshredded-genome cache (isolate spike, Brian 2026-08-26): optional second cache of
		// real whole-genome contigs, loaded only when unshreddedcache= is set (extraunshreddedperfectfrac
		// then draws forced-perfect bins from THIS instead of the shredded cache -- zero shredding
		// artifacts). Off by default: unshreddedCacheFile==null skips this entirely, no behavior change.
		if(unshreddedCacheFile!=null){
			loadCache2();
			int withRecov=0;
			for(int i=0; i<usable.size(); i++){
				final int tid=usable.get(i);
				final ArrayList<Contig> cl=getContigsUnshredded(tid);
				if(cl==null || cl.isEmpty()){continue;}
				long sum=0; for(Contig c : cl){sum+=c.length;}
				recoverable2.put(tid, sum);
				withRecov++;
			}
			System.err.println("unshredded cache: "+tidToIdx2.size()+" orgs loaded, "
				+withRecov+"/"+usable.size()+" usable orgs have unshredded contigs");
		}

		// precompute each organism's NATIVE ncRNA complement (the subnet denominator):
		// {r16,r23,r5,rother,trna} summed over ALL of the tid's contigs.
		for(int i=0; i<usable.size(); i++){
			final int tid=usable.get(i);
			int nR16=0, nR23=0, nR5=0, nRother=0, nTrna=0;
			for(Contig c : getContigs(tid)){
				nR16+=c.r16; nR23+=c.r23; nR5+=c.r5; nRother+=c.rother; nTrna+=c.trna;
			}
			nativeNcR16.put(tid, nR16); nativeNcR23.put(tid, nR23); nativeNcR5.put(tid, nR5);
			nativeNcRother.put(tid, nRother); nativeNcTrna.put(tid, nTrna);
		}
		// SAME, but from the unshredded cache -- a FORCE_PERFECT_UNSHREDDED bin's observed ncRNA
		// must be checked against the UNSHREDDED native complement, not the shredded one (the two
		// differ: the unshredded/whole-genome GFF is the more complete, boundary-artifact-free
		// ground truth). Keyed by the same tid; formatNcrnaRow picks the right map via lastUsedUnshredded.
		if(unshreddedCacheFile!=null){
			for(int i=0; i<usable.size(); i++){
				final int tid=usable.get(i);
				final ArrayList<Contig> cl=getContigsUnshredded(tid);
				if(cl==null){continue;}
				int nR16=0, nR23=0, nR5=0, nRother=0, nTrna=0;
				for(Contig c : cl){
					nR16+=c.r16; nR23+=c.r23; nR5+=c.r5; nRother+=c.rother; nTrna+=c.trna;
				}
				nativeNcR16_2.put(tid, nR16); nativeNcR23_2.put(tid, nR23); nativeNcR5_2.put(tid, nR5);
				nativeNcRother_2.put(tid, nRother); nativeNcTrna_2.put(tid, nTrna);
			}
		}
		// Structural anticodon expected-copy target: sum the 64 structural codes over the
		// target's whole genome.  The residual unknown count is emitted as the explicit 65th
		// model input and is checked against the all-genome tRNA total below.
		if("trna_anticodon".equals(subnetName)){
			for(int i=0; i<usable.size(); i++){
				final int tid=usable.get(i);
				nativeAntiTotal.put(tid, structuralAnticodonTotal(getContigs(tid), tid));
			}
			if(unshreddedCacheFile!=null){
				for(int i=0; i<usable.size(); i++){
					final int tid=usable.get(i);
					final ArrayList<Contig> cl=getContigsUnshredded(tid);
					if(cl==null){continue;}
					nativeAntiTotal_2.put(tid, structuralAnticodonTotal(cl, tid));
				}
			}
		}

		if(enc==ENC_NORM){
			avgCopy=computeAvgCopy(usable);
			System.err.println("enc=norm: computed avgCopyWhenPresent over "+usable.size()+" orgs");
		}
		if(featuresFile!=null){keptRanks=loadRanks(featuresFile);}
		precomputeNStrings();

		final int baseFam=(keptRanks!=null ? keptRanks.length : numFam);
		final int famCols=baseFam*(enc==ENC_TWO ? 2 : 1);
		numInputs=famCols+SHARED_CONTEXT_WIDTH+numPhyla;//formatRow's new context block, not the old raw NUM_GLOBALS dump
		final boolean subnetNcrna=("ncrna".equals(subnetName));
		subnetAnticodon=("trna_anticodon".equals(subnetName));
		subnetRrna=("rrna".equals(subnetName));
		subnetFamset=("famset".equals(subnetName));
		if(subnetName!=null && !subnetNcrna && !subnetFamset && !subnetAnticodon && !subnetRrna){
			throw new RuntimeException("Unknown subnet="+subnetName+" (ncrna|rrna|trna_anticodon|famset)");
		}
		final boolean subnet=subnetNcrna || subnetFamset || subnetAnticodon || subnetRrna;
		if(subnetFamset){
			if(subsetFile==null){throw new RuntimeException("subnet=famset requires subsetfile=<one rank per line>");}
			subsetRanks=loadRanks(subsetFile);
			subsetMask=new boolean[numFam];
			for(int r : subsetRanks){subsetMask[r]=true;}
			lastFamObs=new int[subsetRanks.length];
			// Native subset complement per organism: the subset families' counts summed over
			// ALL the tid's contigs (the famset subnet's denominator target).
			for(int i=0; i<usable.size(); i++){
				final int tid=usable.get(i);
				int sum=0;
				for(Contig c : getContigs(tid)){
					for(int j=0; j<c.famRank.length; j++){if(subsetMask[c.famRank[j]]){sum+=c.famCount[j];}}
				}
				nativeFamTotal.put(tid, sum);
			}
			// SAME, from the unshredded cache (see the ncRNA_2 block above for why a separate map
			// is required rather than reusing nativeFamTotal).
			if(unshreddedCacheFile!=null){
				for(int i=0; i<usable.size(); i++){
					final int tid=usable.get(i);
					final ArrayList<Contig> cl=getContigsUnshredded(tid);
					if(cl==null){continue;}
					int sum=0;
					for(Contig c : cl){
						for(int j=0; j<c.famRank.length; j++){if(subsetMask[c.famRank[j]]){sum+=c.famCount[j];}}
					}
					nativeFamTotal_2.put(tid, sum);
				}
			}
		}
		// Subnet input width: obs block + phylum one-hot + shared context (FROZEN, standard on
		// every net - the old optional domain/hhcaga/genelen blocks are retired, see CTX_N above).
		final int obsCols=(subnetFamset ? subsetRanks.length : subnetAnticodon ? TRNA_ANTICODON_OBS :
			subnetRrna ? RRNA_OBS : NCRNA_OBS);
		subnetInputs=obsCols+numPhyla+SHARED_CONTEXT_WIDTH;
		// Aggregator mode: load manifest, then dense head + buffers.
		if(aggManifestFile!=null || aggBundleFile!=null){
			if(aggManifestFile!=null){loadAggManifest(aggManifestFile);}else{loadAggBundle(aggBundleFile);}
			if(aggOut==null){throw new RuntimeException("aggmanifest= or bundle= requires aggout=");}
			denseRanks=computeDenseHead(usable, denseHead);
			cleanFamBuf=new int[numFam];
			rowCtx=new double[CTX_N];
			numAggInputs=aggInputWidth();
			System.err.println("agg: "+aggSubnets.size()+" subnets, dense head "+denseRanks.length
				+", numAggInputs="+numAggInputs+", obs="+(aggObsServe ? "serve" : "clean"));
		}
		if(rowCtx==null){rowCtx=new double[CTX_N];}//non-aggregator runs still need it for formatRow/subnet rows
		// A4 explicit-panel emission: bypasses the poolMode/writeSet/writeBench machinery entirely
		// (constructor already asserted the aggregator input and unshredded cache are set when this is
		// non-null, so denseRanks/cleanFamBuf/rowCtx/numAggInputs/aggSubnets above are all ready).
		if(panelTidsFile!=null){writePanel(); return;}
		// poolmode=valsplit: both output sets come from the ORIGINAL val orgs (never seen
		// by any subnet trained on the seed-matched train side): first half = aggregator-train
		// (B), second half = final-test (C). Stacking discipline for free (Barbara).
		// poolmode=allbutc: train pool = EVERY usable org except C (Brian 2026-08-11: "hold out
		// vectors, never organisms" - 49 orgs starved the aggregator; the same C stays out so the
		// novel-org readout remains comparable across modes). C is IDENTICAL to valsplit's C.
		if(!noOrganismHoldout && (poolMode==POOL_VALSPLIT || poolMode==POOL_ALLBUTC)){
			final int half=valTids.size()/2;
			if(half<1){throw new RuntimeException("poolmode needs >=2 val orgs, have "+valTids.size());}
			final IntList b=subList(valTids, 0, half);
			final IntList c=subList(valTids, half, valTids.size());
			if(poolMode==POOL_VALSPLIT){
				trainTids=b;
			}else{
				trainTids.addAll(b);//all usable orgs except C
			}
			valTids=c;
			System.err.println("poolmode="+(poolMode==POOL_VALSPLIT ? "valsplit" : "allbutc")
				+": aggregator-train orgs="+trainTids.size()+", final-test orgs="+c.size());
		}
		if(benchMode){
			if(noOrganismHoldout){throw new RuntimeException("benchmark mode requires an explicit held-out pool; disable noorganismholdout");}
			if(benchTruthFile==null || benchManifestFile==null){
				throw new RuntimeException("benchmark mode needs benchtruth= and benchmanifest=");
			}
			if(poolMode!=POOL_ALLBUTC && poolMode!=POOL_VALSPLIT){
				throw new RuntimeException("benchmark mode requires poolmode=allbutc (a defined held-out pool C)");
			}
			writeBench(valTids, benchBins, new Random(seed*3+7));
			System.err.println("done (benchmark).");
			return;
		}
		final MagQCBinManifest.Manifest replay=(binManifestInFile==null ? null : loadReplayManifest(binManifestInFile));
		try(MagQCBinReplayStore store=replayStore){
			if(store!=null){validateReplayRanges(store);}
			if(referenceCdsIndexFile!=null){referenceCdsIndex=MagQCReferenceCdsIndex.load(referenceCdsIndexFile, binManifestInFile);}
			openBinManifest();
			try{
				if(replay==null){
					if(d39Sampling){writeD39Set(out,subnetOut,aggOut,trainTids,(int)n,new Random(seed*2+1),"train");}
					else{writeSet(out, (subnet ? subnetOut : null), aggOut, trainTids, n, new Random(seed*2+1), "train");}
					if(validationOutputRequested() && valn>0 && !valTids.isEmpty()){
						if(d39Sampling){writeD39Set(outval,subnetValOut,aggValOut,valTids,(int)valn,new Random(seed*2+2),"val");}
						else{writeSet(outval, (subnet ? subnetValOut : null), aggValOut, valTids, valn, new Random(seed*2+2), "val");}
					}
				}else{
					writeReplaySet(out, (subnet ? subnetOut : null), aggOut, replay, "train",store,
						replayTrainStartResolved,replayTrainEndResolved);
					if(validationOutputRequested()){
						writeReplaySet(outval, (subnet ? subnetValOut : null), aggValOut, replay, "val",store,
							replayValStartResolved,replayValEndResolved);
					}
				}
			}finally{
				closeBinManifest();
			}
		}
		if(subnetAnticodon){
			System.err.println("structural anticodon diagnostic: rows="+anticodonRows
				+" rows_with_unknown="+anticodonMissingRows+" (unknown is an explicit 65th input; vector width="+TRNA_ANTICODON_OBS+")");
		}
		System.err.println("done.");
	}

	/** Dense suppression must not accidentally suppress requested subnet validation vectors. */
	private boolean validationOutputRequested(){
		return emitDense ? outval!=null : subnetValOut!=null || aggValOut!=null;
	}

	/** Reads the explicit approved exclusion ledger, never inferring exclusions from missing labels. */
	private void loadExcludedTids(){
		if(excludedTidsFile==null){return;}
		final ByteFile input=ByteFile.makeByteFile(excludedTidsFile,false);
		final LineParser1 parser=new LineParser1((byte)'\t');
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				parser.set(line);
				final int tid=parser.parseInt(0);
				if(tid<1 || !excludedTids.add(tid)){
					throw new RuntimeException("Invalid or duplicate excluded tid: "+tid);
				}
			}
		}finally{if(input.close()){throw new RuntimeException("I/O error reading "+excludedTidsFile);}}
	}

	/**
	 * Loads D141 reference groups solely for source sampling. C1 feature maps are untouched.
	 * Every eligible organism must have exactly one reference row; excluded organisms cannot
	 * enter either sampling pool. At least two tracked organisms permit same-phylum mixtures.
	 */
	private void loadSamplingPhylum(IntList eligible){
		assert(samplePhylum!=null && samplingTaxonomyFile!=null) : "Constructor requires paired reference sampling arguments";
		final IntList[] pools=samplingPhylumPools(samplingTaxonomyFile,samplePhylum,eligible,samplingUnknownOrdinary);
		samplingInside=pools[0]; samplingOutside=pools[1];
		System.err.println("D39 reference sampling phylum="+samplePhylum+" inside="+samplingInside.size()+
			" outside="+samplingOutside.size()+" ordinary_only="+(eligible.size()-samplingInside.size()-samplingOutside.size()));
	}

	/** D145: unknown reference rows stay in the caller's ordinary pool, never either forced pool. */
	static IntList[] samplingPhylumPools(String taxonomy,String phylum,IntList eligible,boolean allowUnknown){
		assert(taxonomy!=null && phylum!=null && eligible!=null) : "Reference sampling needs explicit taxonomy, phylum and eligibility";
		final IntHashSet seen=new IntHashSet(eligible.size()*2+1),inside=new IntHashSet(eligible.size()+1);
		final IntHashSet unknown=new IntHashSet(32);
		final ByteFile input=ByteFile.makeByteFile(taxonomy,false);
		final LineParser1 parser=new LineParser1((byte)'\t');
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				parser.set(line);
				if(parser.terms()<2){throw new RuntimeException("Sampling taxonomy requires tid and phylum columns");}
				final int tid=parser.parseInt(0);
				if(tid<1 || !seen.add(tid)){throw new RuntimeException("Invalid or duplicate sampling-taxonomy tid: "+tid);}
				if(parser.length(1)==0){throw new RuntimeException("Empty reference sampling phylum for tid: "+tid);}
				if(parser.termEquals("-",1)){
					if(!allowUnknown){throw new RuntimeException("Unresolved reference phylum requires samplingunknown=ordinary: "+tid);}
					unknown.add(tid);
				}else if(parser.termEquals(phylum,1)){inside.add(tid);}
			}
		}finally{if(input.close()){throw new RuntimeException("I/O error reading "+taxonomy);}}
		final IntList samplingInside=new IntList(),samplingOutside=new IntList();
		for(int i=0; i<eligible.size(); i++){
			final int tid=eligible.get(i);
			if(!seen.contains(tid)){throw new RuntimeException("Eligible tid lacks reference sampling taxonomy: "+tid);}
			if(unknown.contains(tid)){continue;}
			(inside.contains(tid) ? samplingInside : samplingOutside).add(tid);
		}
		if(samplingInside.size()<2 || samplingOutside.size()<1){
			throw new RuntimeException("Phylum sampling needs at least two eligible tracked organisms and one outside: "+phylum);
		}
		return new IntList[]{samplingInside,samplingOutside};
	}

	/** Closes the optional memory-mapped contribution index after every normal or failed run. */
	private void closeContributionIndex(Throwable priorFailure){
		if(referenceCdsContributionIndex==null){return;}
		try{referenceCdsContributionIndex.close();}
		catch(java.io.IOException e){
			final RuntimeException closeFailure=new RuntimeException(
				"Could not close reference-CDS contribution index",e);
			if(priorFailure!=null){priorFailure.addSuppressed(closeFailure);}
			else{throw closeFailure;}
		}
		finally{referenceCdsContributionIndex=null;}
	}

	/** Initializes only the deterministic cache/vectorizer state needed by the deployment seam.
	 * No training or sampled-bin output is performed. */
	void initializeSingleBin(){
		if(aggBundleFile==null){throw new RuntimeException("single-bin initialization requires bundle=");}
		assert(enc!=ENC_NORM) : "single-bin aggregator deployment does not support enc=norm";
		loadAux(); loadCache();
		final IntList usable=new IntList();
		for(int tid : tidToIdx.toArray()){if(genomeSize.contains(tid) && !getContigs(tid).isEmpty()){usable.add(tid);}}
		usable.sort(); computeNormScales(usable); precomputeNStrings();
		final int famCols=(enc==ENC_TWO ? numFam*2 : numFam);
		numInputs=famCols+SHARED_CONTEXT_WIDTH+numPhyla;
		loadAggBundle(aggBundleFile); denseRanks=computeDenseHead(usable,denseHead);
		cleanFamBuf=new int[numFam]; rowCtx=new double[CTX_N];
		numAggInputs=aggInputWidth();
		if(numAggInputs<=0){throw new RuntimeException("Invalid single-bin aggregator width "+numAggInputs);}
	}

	/**
	 * Initializes inference from frozen resources alone, without a reference cache, TIDs or labels.
	 * Prepared bins supply whole-bin observations; scales and vocabulary always come from the bundle.
	 */
	static MagQCVectorMaker initializePrepared(String bundle, String family, String release, String releasePin,
			String table, String tablePin, String populations, String populationsPin){
		return initializePrepared(bundle, family, release, releasePin, table, tablePin, populations, populationsPin, true);
	}

	/** The public client can silence startup diagnostics without changing replay-tool defaults. */
	static MagQCVectorMaker initializePrepared(String bundle, String family, String release, String releasePin,
			String table, String tablePin, String populations, String populationsPin, boolean verbose){
		final MagQCVectorMaker vm=new MagQCVectorMaker();
		vm.verbose=verbose;
		vm.aggBundleFile=bundle; vm.familyFile=family;
		vm.releaseManifestFile=release; vm.releaseManifestSha80=releasePin;
		vm.expectedCopyTableFile=table; vm.expectedCopyTableSha256=tablePin;
		vm.subnetPopulationsFile=populations; vm.subnetPopulationsSha256=populationsPin;
		// Production composite recipe: run_composite_131_canary_20260928.sh uses enc=raw.
		vm.subnetFeatures=SUBNETFEATURES_SIX; vm.aggObsServe=true; vm.enc=ENC_RAW;
		vm.loadConfiguredBundle();
		if(vm.frozenInputs==null){throw new IllegalArgumentException("Prepared inference requires a frozen subnet input definition");}
		vm.loadFamilyCount();
		vm.phylumList=new ArrayList<String>(vm.frozenInputs.phyla());
		vm.numPhyla=vm.phylumList.size();
		for(int i=0; i<vm.numPhyla; i++){vm.phylumIndex.put(vm.phylumList.get(i), i);}
		vm.computeNormScales(new IntList());
		vm.precomputeNStrings();
		vm.loadAggBundle(bundle);
		// D196 contract: prepared-input BBNet evaluation uses the same private dense inference
		// representation as replay workers, avoiding sparse/dense accumulation drift while leaving
		// legacy bundle-serving APIs and generation paths unchanged.
		vm.preparedDenseInference=true;
		vm.denseHead=vm.numFam; vm.denseRanks=vm.computeDenseHead(new IntList(), vm.numFam);
		vm.rowCtx=new double[CTX_N]; vm.numAggInputs=vm.aggInputWidth();
		vm.preparedGlob=new double[NUM_GLOBALS];
		return vm;
	}

	/** Number of family ranks accepted by the prepared-bin parser. */
	int preparedFamilyCount(){return numFam;}
	/**
	 * Shares validated resources while giving one batch worker its own formatter and subnet state.
	 * The replay copy already detaches every scratch buffer used by formatAggRowCore and expands
	 * each subnet through the same dense inference factory used by prepared serving.
	 */
	MagQCVectorMaker newPreparedWorker(){return newPreparedWorker(false);}

	/** Selects existing private copies or one locked dense instance per subnet. */
	MagQCVectorMaker newPreparedWorker(boolean lockedSubnets){
		if(preparedWorker){throw new IllegalStateException("Create prepared workers from the initialized resource master");}
		if(!preparedDenseInference || preparedGlob==null || frozenInputs==null){
			throw new IllegalStateException("Prepared workers require initializePrepared resources");
		}
		for(AggSubnet subnet:aggSubnets){
			if(subnet.bundleSubnet==null){
				throw new IllegalStateException("Create prepared workers from the initialized resource master, not another worker");
			}
		}
		final MagQCVectorMaker worker=copyForReplayThread(lockedSubnets);
		assert(worker.preparedGlob!=preparedGlob && worker.rowCtx!=rowCtx && worker.aggSubnets!=aggSubnets) :
			"Concurrent prepared bins must not share formatter scratch or mutable subnet workers";
		return worker;
	}
	/** Input-only output width, derived from the bound release rather than a caller-supplied number. */
	int preparedInputWidth(){return numAggInputs;}

	/** Appends exactly the replay formatter's input fields for one validated whole-bin observation. */
	void formatPreparedBin(ByteBuilder out, MagQCPreparedBin bin){
		if(preparedGlob==null || bin.families.length!=numFam){throw new IllegalStateException("Prepared inference state/family width mismatch");}
		bin.validate();
		fillFeatureGlobals(bin.length, bin.agg, preparedGlob);
		snapshotServeSpecial(bin.agg);
		final boolean unknown=bin.status.equals("unknown"), classified=bin.status.equals("classified");
		final int phylum=phylumColumn(classified ? bin.phylum : null);
		final int domain=domainColumn(unknown ? null : bin.domain);
		formatAggRowCore(out, bin.families, preparedGlob, phylum, null, bin.agg, false, domain, bin.id);
	}

	/** Aggregator row width: subnetBlockWidth*subnets + pooledCols + dense head(2x) + phylum +
	 *  shared context + raw ncRNA two-channel(2x5) + mapped-fraction(1). subnetBlockWidth/pooledCols
	 *  are 4/1 under subnetfeatures=legacy (byte-identical to the pre-six formula) and 6/0 under six
	 *  (design revision2 §4.1, "root: one shared calculation") -- set by loadAggBundleLegacy/Six. */
	private int aggInputWidth(){
		if(subnetFeatures==SUBNETFEATURES_SIX){
			return sixCompositeInputWidth(aggSubnets.size(),denseRanks.length,numPhyla);
		}
		return aggSubnets.size()*subnetBlockWidth+pooledCols+2*denseRanks.length+numPhyla
			+SHARED_CONTEXT_WIDTH+2*NCRNA_OBS+1;
	}

	/** Shared width calculation for six-feature generation, packaging and deployment preflight. */
	static int sixCompositeInputWidth(int subnetCount,int familyCount,int phylumCount){
		if(subnetCount<1 || familyCount<1 || phylumCount<1){
			throw new IllegalArgumentException("Composite inputs require subnets, families and phylum columns");
		}
		final long width=6L*subnetCount+2L*familyCount+phylumCount+SHARED_CONTEXT_WIDTH+2*NCRNA_OBS+1;
		if(width>Integer.MAX_VALUE){throw new IllegalArgumentException("Composite input width exceeds integer range");}
		return (int)width;
	}

	/**
	 * Complete immutable serving contract: ordered six-value subnet blocks (four learned,
	 * O/E, X/E), full rank-ordered family presence/excess-cap16 pairs, frozen phylum/context,
	 * five ncRNA presence and N/(1+N) pairs, then mapped/CDS. No legacy pooled ratio or
	 * output clamp. Bump this ID when any column order, transform or numerical convention changes.
	 */
	static final String JOINT_COMPOSITE_LAYOUT="magqc_six_full_families_ratio_v1";

	/** Package-visible for MagQCTool's phase-0 configuration validation (design revision2 §4.3 step
	 *  2): derives numFam/numPhyla through the SAME read-only vocabulary loading MVM itself uses,
	 *  rather than a second, divergence-prone reimplementation in MagQCTool. Touches no cache row,
	 *  writer, or output path. */
	void loadAuxForPhase0(){loadAux();}

	/** Uses the already-validated release vocabulary when inspecting deployment labels. */
	void loadAuxForPhase0(MagQCNetBundle.FrozenInputs inputs){
		assert(inputs!=null) : "Frozen phase-0 loading requires the validated bundle input definition";
		frozenInputs=inputs;
		loadAux();
	}
	int numFam(){return numFam;}
	int numPhyla(){return numPhyla;}
	/** Package-visible, read-only provenance snapshot (Phase1ContextScaleRecovery): the exact ordered
	 *  phylum vocabulary loadAux() built from labels= or taxpgm=, copied so later loader mutation
	 *  cannot change a recorded result. Requires loadAux() to have completed. */
	java.util.List<String> phylumVocabulary(){
		if(phylumList==null || phylumList.size()!=numPhyla || numPhyla<1){
			throw new IllegalStateException("phylumVocabulary() requires a completed loadAux(); numPhyla="+numPhyla
				+" list="+(phylumList==null ? "null" : Integer.toString(phylumList.size())));
		}
		return java.util.Collections.unmodifiableList(new ArrayList<String>(phylumList));
	}

	/** Replicates java.util.Collections.shuffle(List,Random)'s exact algorithm and RNG
	 *  call sequence (Fisher-Yates, i from size down to 2, swap(i-1, rnd.nextInt(i))) -
	 *  determinism-critical: the global train/val organism split depends on drawing the
	 *  SAME sequence of rnd.nextInt() calls the original produced. */
	private static void shuffleInPlace(IntList list, Random rnd){
		for(int i=list.size(); i>1; i--){
			final int j=rnd.nextInt(i);
			final int tmp=list.get(i-1);
			list.set(i-1, list.get(j));
			list.set(j, tmp);
		}
	}

	private static IntList subList(IntList src, int from, int to){
		final IntList out=new IntList(Math.max(1, to-from));
		for(int i=from; i<to; i++){out.add(src.get(i));}
		return out;
	}

	/** Builds family->tids index within a pool (for same-family contaminant selection). */
	private HashMap<Integer,IntList> familyIndex(IntList pool){
		final HashMap<Integer,IntList> m=new HashMap<Integer,IntList>();
		for(int i=0; i<pool.size(); i++){
			final int tid=pool.get(i);
			final int fam=tid2family.get(tid);//IntHashMap.get returns -1 when absent, matching getOrDefault(tid,-1)
			IntList l=m.get(fam);
			if(l==null){m.put(fam, l=new IntList());}
			l.add(tid);
		}
		return m;
	}

	private void writeSet(String file, String subnetFile, String aggFile, IntList pool, long count, Random rnd, String tag){
		final HashMap<Integer,IntList> famIdx=familyIndex(pool);
		final ByteStreamWriter bsw=new ByteStreamWriter(file, true, false, true);
		bsw.start();
		final ByteBuilder bb=new ByteBuilder(numInputs*4+64);
		bb.append("#dims\t").append(numInputs).append("\t2\t0").nl();
		bsw.print(bb); bb.clear();

		final ByteStreamWriter sbsw=(subnetFile==null ? null : new ByteStreamWriter(subnetFile, true, false, true));
		final ByteBuilder bb2=(sbsw==null ? null : new ByteBuilder(subnetInputs*4+64));
		if(sbsw!=null){
			sbsw.start();
			bb2.append("#dims\t").append(subnetInputs)
				.append(subnetLabels==SUBNETLABELS_GENE ? "\t2\t0" : "\t1\t0").nl();
			sbsw.print(bb2); bb2.clear();
		}

		final ByteStreamWriter absw=(aggFile==null ? null : new ByteStreamWriter(aggFile, true, false, true));
		final ByteBuilder bb3=(absw==null ? null : new ByteBuilder(numAggInputs*4+64));
		if(absw!=null){
			absw.start();
			bb3.append("#dims\t").append(numAggInputs).append("\t2\t0").nl();
			absw.print(bb3); bb3.clear();
		}

		// Per-attempt buffers hoisted OUT of the loop and reused (S3): fam[]/glob[]/labels[]
		// are cleared, never reallocated. fam[] must be explicitly zeroed every attempt -
		// Agg.add() only INCREMENTS specific indices, it never zeroes the whole array, so a
		// stale value from a REJECTED prior attempt would otherwise bleed into the next one.
		// Agg itself (S4) is built ONCE and reset() between attempts, so its scratchArr/lens
		// buffers (already designed to reuse-if-large-enough) actually get to amortize.
		final double[] labels=new double[2];
		final int[] fam=new int[numFam];
		final double[] glob=new double[NUM_GLOBALS];
		final Agg agg=new Agg();
		agg.setFam(fam);

		long made=0, tries=0;
		// try/finally (2026-09-03): makeBin can now throw a real RuntimeException mid-loop (labels=
		// mode, a cache tid absent from labels.tsv) -- bsw/sbsw/absw are non-daemon ByteStreamWriter
		// threads parked on their own job queue, and without this the poisonAndWait() calls below
		// never ran on that path, leaking a live thread forever and hanging the whole JVM at exit
		// (found via a real hang in MagQCVectorMakerLabelsTest's missing-tid crash-loud case).
		try{
			while(made<count && tries<count*20+1000){
				tries++;
				Arrays.fill(fam, 0, numFam, 0);
				agg.reset();
				if(binManifestOutFile!=null){benchNative.clear(); benchForeign.clear();}
				final int targetPhylumIdx=makeBin(pool, famIdx, rnd, fam, glob, labels, agg, FORCE_NONE);
				if(targetPhylumIdx<0){continue;}
				writeBinManifestRow(tag, labels);
				emitBin(bb, bb2, bb3, bsw, sbsw, absw, fam, glob, targetPhylumIdx, labels, agg);
				made++;
				if((made%50000)==0){System.err.println(tag+": "+made+"/"+count);}
			}
			// Additional isolate spike (Brian 2026-08-26): APPENDED on top of the base `count` bins above,
			// never replacing them -- distinct from perfectFrac/nearPerfectFrac, which shrink the normal
			// population instead. Off (0.0) by default, so this loop runs zero iterations and the base
			// draw's output is completely unaffected. Scope controlled by extraSpikeSet (train/val/both);
			// `tag` here is literally "train" or "val" (writeSet's own caller-supplied label), so the check
			// is a direct string match, no separate train/val flag needed.
			if(extraSpikeSet.equals("both") || extraSpikeSet.equals(tag)){
				final long extraPerfect=Math.round(extraPerfectFrac*count);
				final long extraNear=Math.round(extraNearPerfectFrac*count);
				long madeExtra=0, extraTries=0;
				while(madeExtra<extraPerfect && extraTries<extraPerfect*20+1000){
					extraTries++;
					Arrays.fill(fam, 0, numFam, 0);
					agg.reset();
					if(binManifestOutFile!=null){benchNative.clear(); benchForeign.clear();}
					final int targetPhylumIdx=makeBin(pool, famIdx, rnd, fam, glob, labels, agg, FORCE_PERFECT);
					if(targetPhylumIdx<0){continue;}
					writeBinManifestRow(tag, labels);
					emitBin(bb, bb2, bb3, bsw, sbsw, absw, fam, glob, targetPhylumIdx, labels, agg);
					madeExtra++;
				}
				long madeExtraNear=0; long extraNearTries=0;
				while(madeExtraNear<extraNear && extraNearTries<extraNear*20+1000){
					extraNearTries++;
					Arrays.fill(fam, 0, numFam, 0);
					agg.reset();
					if(binManifestOutFile!=null){benchNative.clear(); benchForeign.clear();}
					final int targetPhylumIdx=makeBin(pool, famIdx, rnd, fam, glob, labels, agg, FORCE_NEARPERFECT);
					if(targetPhylumIdx<0){continue;}
					writeBinManifestRow(tag, labels);
					emitBin(bb, bb2, bb3, bsw, sbsw, absw, fam, glob, targetPhylumIdx, labels, agg);
					madeExtraNear++;
				}
				System.err.println(tag+": appended "+madeExtra+" extra-perfect + "+madeExtraNear+" extra-near-perfect bins"
					+" (extraperfectfrac="+extraPerfectFrac+", extranearperfectfrac="+extraNearPerfectFrac+")");
				made+=madeExtra+madeExtraNear;

				// Unshredded-perfect append (Brian 2026-08-26): ADDITIONAL to the shredded-perfect append
				// above, not a replacement -- both spikes are additive on top of the base draw. Skipped
				// entirely (zero iterations, zero RNG draws) when extraUnshreddedPerfectFrac==0, so a run
				// without unshreddedcache= is completely unaffected.
				final long extraUnshreddedPerfect=Math.round(extraUnshreddedPerfectFrac*count);
				long madeExtraUnshredded=0, extraUnshreddedTries=0;
				while(madeExtraUnshredded<extraUnshreddedPerfect && extraUnshreddedTries<extraUnshreddedPerfect*20+1000){
					extraUnshreddedTries++;
					Arrays.fill(fam, 0, numFam, 0);
					agg.reset();
					if(binManifestOutFile!=null){benchNative.clear(); benchForeign.clear();}
					final int targetPhylumIdx=makeBin(pool, famIdx, rnd, fam, glob, labels, agg, FORCE_PERFECT_UNSHREDDED);
					if(targetPhylumIdx<0){continue;}
					writeBinManifestRow(tag, labels);
					emitBin(bb, bb2, bb3, bsw, sbsw, absw, fam, glob, targetPhylumIdx, labels, agg);
					madeExtraUnshredded++;
				}
				System.err.println(tag+": appended "+madeExtraUnshredded+" extra-UNSHREDDED-perfect bins"
					+" (extraunshreddedperfectfrac="+extraUnshreddedPerfectFrac+")");
				made+=madeExtraUnshredded;
			}
		}finally{
			bsw.poisonAndWait();
			if(sbsw!=null){sbsw.poisonAndWait(); System.err.println(tag+": wrote "+made+" ncRNA-subnet rows to "+subnetFile);}
			if(absw!=null){absw.poisonAndWait(); System.err.println(tag+": wrote "+made+" aggregator rows to "+aggFile);}
			System.err.println(tag+": wrote "+made+" rows to "+file+" (tries="+tries+")");
		}
	}

	/**
	 * Generates the corpus-wide D39 mixture with exact integer category counts.
	 * Every accepted draw is replayed through the same exact-label path used by
	 * binmanifestin before emission. Realization fingerprints span both output
	 * splits, preventing the same selected sources/cuts from entering both sets.
	 * When samplephylum is set, the same path also enforces the per-origin phylum strata.
	 */
	private void writeD39Set(String file,String subnetFile,String aggFile,IntList pool,
			int count,Random rnd,String tag){
		if(pool.size()<2){throw new RuntimeException("D39 contaminated categories require at least two source organisms");}
		for(int i=0; i<pool.size(); i++){
			final int tid=pool.get(i);
			if(taxonomyFromLabels && !tid2phylumIdxLabels.containsKey(tid)){
				throw new RuntimeException("Eligible D39 tid has no C1 label: "+tid+"; check the explicit exclusion ledger");
			}
			if(getContigsUnshredded(tid)==null || referenceCdsContributionIndex.nativeTotal(tid)<1){
				throw new RuntimeException("D39 needs intact reference genes and contigs for every eligible tid: "+tid);
			}
		}
		final int[] remaining=MagQCD39Quotas.allocate(count,samplePhylum!=null);
		final int[] requested=remaining.clone();
		final int[] forced={FORCE_PERFECT_UNSHREDDED,FORCE_D39_HQ,FORCE_D39_ISOLATE,FORCE_D39_ORDINARY};
		final HashMap<Integer,IntList> famIdx=familyIndex(pool);
		final IntList[] targetPools=samplePhylum==null ? new IntList[]{pool} :
			new IntList[]{samplingInside,samplingOutside,samplingInside,pool};
		final IntList[] foreignPools=samplePhylum==null ? new IntList[]{pool} :
			new IntList[]{samplingOutside,samplingInside,samplingInside,pool};
		final ArrayList<HashMap<Integer,IntList>> foreignIndexes=new ArrayList<HashMap<Integer,IntList>>();
		for(IntList foreignPool : foreignPools){foreignIndexes.add(foreignPool==pool ? famIdx : familyIndex(foreignPool));}
		final double[] labels=new double[2],glob=new double[NUM_GLOBALS];
		final int[] fam=new int[numFam];
		final Agg agg=new Agg(); agg.setFam(fam);
		final ByteBuilder dense=(emitDense ? new ByteBuilder(numInputs*4+64) : null),subnet=new ByteBuilder(subnetInputs*4+64);
		final ByteBuilder composite=new ByteBuilder(numAggInputs*4+64);
		final IntList broken=new IntList();
		final java.security.MessageDigest digest;
		try{digest=java.security.MessageDigest.getInstance("SHA-256");}
		catch(java.security.NoSuchAlgorithmException e){throw new IllegalStateException(e);}
		ByteStreamWriter dw=null,sw=null,aw=null;
		Throwable failure=null;
		try{
			if(emitDense){
				dw=new ByteStreamWriter(file,true,false,true); dw.start();
				dense.append("#dims\t").append(numInputs).append("\t2\t0\n"); dw.print(dense); dense.clear();
			}
			if(subnetFile!=null){
				sw=new ByteStreamWriter(subnetFile,true,false,true); sw.start();
				subnet.append("#dims\t").append(subnetInputs).append("\t2\t0\n"); sw.print(subnet); subnet.clear();
			}
			if(aggFile!=null){
				aw=new ByteStreamWriter(aggFile,true,false,true); aw.start();
				composite.append("#dims\t").append(numAggInputs).append("\t2\t0\n"); aw.print(composite); composite.clear();
			}
			int made=0;
			long attempts=0;
			while(made<count){
				if(++attempts>count*200L*splitModulus+1000){
					throw new RuntimeException("Cannot fill unique D39 quotas for "+sampleModel+" "+tag+
						": remaining="+Arrays.toString(remaining)+" attempts="+attempts);
				}
				int roll=rnd.nextInt(count-made),slot=0;
				while(roll>=remaining[slot]){roll-=remaining[slot++];}
				final int category=slot%4,role=slot/4;
				Arrays.fill(fam,0); agg.reset(); benchNative.clear(); benchForeign.clear(); broken.clear();
				final boolean requireForeign=samplePhylum!=null && (role==1 || role==2);
				if(makeBin(targetPools[role],foreignPools[role],foreignIndexes.get(role),rnd,fam,glob,labels,agg,
						forced[category],requireForeign)<0){continue;}
				java.util.Collections.sort(benchNative,D39_CONTIG_ORDER);
				java.util.Collections.sort(benchForeign,D39_CONTIG_ORDER);
				String cuts="-";
				if(category==0){
					final int genes=referenceCdsContributionIndex.nativeTotal(lastTarget);
					// Preserve the specified zero-based min(rand(500),rand(400),rand(400)) draw.
					// Zero cuts are intact near-perfect examples; fingerprints deduplicate them.
					final int remove=Math.min(genes,Math.min(rnd.nextInt(500),Math.min(rnd.nextInt(400),rnd.nextInt(400))));
					final IntHashSet seen=new IntHashSet(16);
					while(broken.size()<remove){final int ordinal=rnd.nextInt(genes); if(seen.add(ordinal)){broken.add(ordinal);}}
					broken.sort();
					final ByteBuilder cutText=new ByteBuilder();
					for(int j=0; j<broken.size(); j++){if(j>0){cutText.append(';');} cutText.append('g').append(broken.get(j));}
					if(broken.size()>0){cuts=cutText.toString();}
				}
				final boolean whole=(category==0 || category==2);
				final D39Fingerprint fingerprint=d39Fingerprint(digest,whole,broken);
				// Complete clean versions of one genome stay together even across fragmentation.
				// Primary-source labels already came from the exact survival table in makeBin;
				// an intact-source bin retains all reference genes exactly when it has no cuts.
				final boolean completeClean=benchForeign.isEmpty() && (whole ? broken.size()==0 : labels[0]==1.0);
				final byte[] splitKey;
				if(completeClean){
					digest.reset(); digest.update((byte)2); updateFingerprintInt(digest,lastTarget);
					splitKey=digest.digest();
				}else{splitKey=fingerprint.bytes;}
				if(MagQCBinSplit.validation(splitKey,splitModulus)!=tag.equals("val")){continue;}
				if(!d39Realizations.add(fingerprint)){continue;}
				final java.util.TreeSet<Integer> foreignTids=new java.util.TreeSet<Integer>();
				for(Contig c : benchForeign){foreignTids.add(c.tid);}
				final StringBuilder tids=new StringBuilder();
				for(int tid : foreignTids){if(tids.length()>0){tids.append(';');} tids.append(tid);}
				final String roleTag=samplePhylum==null ? "" : "_"+MagQCD39Quotas.ROLES[role];
				final MagQCBinManifest.Bin bin=new MagQCBinManifest.Bin(sampleModel+"_"+seed+"_"+tag+roleTag+"_"+made,
					tag,lastTarget,tids.length()==0 ? "-" : tids.toString(),lastMixtureComp,lastMixtureCont,
					lastSpikeClass,cuts,contigListText(benchNative,'n'),contigListText(benchForeign,'f'),
					sampleModel+":"+tag+":"+made,MagQCBinManifest.SCHEMA_VERSION_ENCODED);
				// Sampling merely chooses sources. Canonical replay owns features and exact gene labels.
				Arrays.fill(fam,0); agg.reset();
				final int phylum=replayBin(bin,fam,glob,labels,agg);
				if(phylum<0){throw new RuntimeException("D39 accepted draw failed canonical replay: "+bin.id);}
				final ReferenceCdsSurvivalLabelReader.Counts reference=exactReplayCounts(bin);
				if(reference==null){throw new AssertionError("D39 constructor requires both exact reference-gene sources");}
				labels[0]=reference.completeness; labels[1]=reference.contamination;
				appendManifestBin(bin);
				manifestOrdinal++;
				emitBin(dense,subnet,composite,dw,sw,aw,fam,glob,phylum,labels,agg);
				remaining[slot]--; made++;
				if(made%50000==0){System.err.println(tag+": D39 "+made+"/"+count);}
			}
			System.err.println(tag+": D39 exact categories breakpoint,HQ,isolate,ordinary"+
				(samplePhylum==null ? "" : " grouped by main,foreign,both,remaining")+"="+
				Arrays.toString(requested)+" rows="+made+" attempts="+attempts);
		}catch(RuntimeException e){failure=e; throw e;}
		catch(Error e){failure=e; throw e;}
		finally{
			boolean ioError=false;
			if(dw!=null){ioError|=dw.poisonAndWait();}
			if(sw!=null){ioError|=sw.poisonAndWait();}
			if(aw!=null){ioError|=aw.poisonAndWait();}
			if(ioError){
				final RuntimeException e=new RuntimeException("I/O error writing D39 "+tag+" vectors");
				if(failure!=null){failure.addSuppressed(e);}else{throw e;}
			}
		}
	}

	/** Canonical source-index fingerprint, shared across train and validation in this run. */
	private D39Fingerprint d39Fingerprint(java.security.MessageDigest digest,boolean whole,IntList broken){
		digest.reset(); digest.update((byte)(whole ? 1 : 0));
		updateFingerprintInt(digest,lastTarget);
		if(!whole){
			updateFingerprintInt(digest,benchNative.size());
			for(Contig c : benchNative){
				assert(c.tableIdx>=0) : "Primary cache contigs must be bound to the survival table before D39 sampling";
				updateFingerprintInt(digest,c.tableIdx);
			}
		}
		updateFingerprintInt(digest,benchForeign.size());
		for(Contig c : benchForeign){updateFingerprintInt(digest,c.tableIdx);}
		updateFingerprintInt(digest,broken.size());
		for(int i=0; i<broken.size(); i++){updateFingerprintInt(digest,broken.get(i));}
		return new D39Fingerprint(digest.digest());
	}

	/** Fixed-width integer encoding avoids per-contig Strings or byte buffers while hashing. */
	private static void updateFingerprintInt(java.security.MessageDigest digest,int value){
		assert(value>=0) : "D39 fingerprints encode nonnegative tids, source indices, counts and ordinals";
		digest.update((byte)(value>>>24)); digest.update((byte)(value>>>16));
		digest.update((byte)(value>>>8)); digest.update((byte)value);
	}

	/** A retained digest key bounds deduplication memory independently of contig-name lengths. */
	private static final class D39Fingerprint{
		D39Fingerprint(byte[] bytes_){bytes=bytes_; hash=Arrays.hashCode(bytes);}
		@Override public int hashCode(){return hash;}
		@Override public boolean equals(Object other){return other instanceof D39Fingerprint && Arrays.equals(bytes,((D39Fingerprint)other).bytes);}
		final byte[] bytes;
		final int hash;
	}

	private static final java.util.Comparator<Contig> D39_CONTIG_ORDER=new java.util.Comparator<Contig>(){
		@Override public int compare(Contig a,Contig b){
			int c=Integer.compare(a.tableIdx,b.tableIdx);
			if(c==0){c=Integer.compare(a.tid,b.tid);}
			return c==0 ? a.name.compareTo(b.name) : c;
		}
	};

	/** Loads and verifies an explicit selection manifest before replay. */
	private MagQCBinManifest.Manifest loadReplayManifest(String file){
		if(binModel!=null){
			replayStore=MagQCBinReplayStore.build(file,binModel,
				replayTmpDir==null ? null : java.nio.file.Paths.get(replayTmpDir),new MagQCBinReplayStore.Validator(){
					@Override public void header(MagQCBinManifest.Manifest m){validateReplayHeader(m);}
					@Override public void bin(MagQCBinManifest.Bin b,MagQCBinManifest.Manifest m){validateReplayBin(b,m);}
				},replayShuffle ? Long.valueOf(replaySeed) : null,replayMemory);
			System.err.println("Replay indexed "+replayStore.rowsRead()+" source rows; selected train="+
				replayStore.size("train")+" val="+replayStore.size("val")+"; selected payload bytes on disk="+replayStore.spoolBytes());
			return replayStore.header();
		}
		final MagQCBinManifest.Manifest m=MagQCBinManifest.load(file);
		validateReplayHeader(m);
		for(MagQCBinManifest.Bin b : m.bins){validateReplayBin(b,m);}
		return m;
	}

	/** Returns whether any split-local named-replay range bound was supplied. */
	private boolean hasReplayRange(){return replayTrainStartSet || replayTrainEndSet || replayValStartSet || replayValEndSet;}

	/** Validates explicit split-local half-open replay ranges before any replay writer opens. */
	private void validateReplayRanges(MagQCBinReplayStore store){
		final int trainSize=store.size("train"),valSize=store.size("val");
		replayTrainStartResolved=validateReplayRange("train",replayTrainStartSet ? replayTrainStart : 0,
			replayTrainEndSet ? replayTrainEnd : trainSize,trainSize);
		replayTrainEndResolved=(int)(replayTrainEndSet ? replayTrainEnd : trainSize);
		replayValStartResolved=validateReplayRange("val",replayValStartSet ? replayValStart : 0,
			replayValEndSet ? replayValEnd : valSize,valSize);
		replayValEndResolved=(int)(replayValEndSet ? replayValEnd : valSize);
	}

	/** Validates one split-local half-open range without clamping explicit bounds. */
	private static int validateReplayRange(String split,long start,long end,int size){
		if(start<0 || end<0){throw new IllegalArgumentException("Replay "+split+" range must be nonnegative: ["+start+","+end+")");}
		if(start>end){throw new IllegalArgumentException("Replay "+split+" range is reversed: ["+start+","+end+")");}
		if(start>size || end>size){throw new IllegalArgumentException("Replay "+split+" range ["+start+","+end+") exceeds selected model size "+size);}
		return (int)start;
	}

	/** Verifies source bindings before any selected payload is accepted for replay. */
	private void validateReplayHeader(MagQCBinManifest.Manifest m){
		if(!cacheFile.equals(m.sourcePath) || "-".equals(m.sourceSha256)
			|| !m.sourceSha256.equals(MagQCBinManifest.sha256File(cacheFile))){
			throw new RuntimeException("Manifest primary source hash does not match cache="+cacheFile);
		}
		if(!"-".equals(m.secondarySourcePath)){
			if(unshreddedCacheFile==null || !m.secondarySourcePath.equals(unshreddedCacheFile)){
				throw new RuntimeException("Manifest secondary source does not match unshreddedcache=");
			}
			if(!m.secondarySourceSha256.equals(MagQCBinManifest.sha256File(unshreddedCacheFile))){
				throw new RuntimeException("Manifest secondary source hash does not match unshreddedcache=");
			}
		}
	}

	/** Applies exclusion and exact-label capability checks even to nonselected model rows. */
	private void validateReplayBin(MagQCBinManifest.Bin b,MagQCBinManifest.Manifest m){
		if(excludedTids.contains(b.targetTid)){
			throw new RuntimeException("Replay selects an explicitly excluded target tid: "+b.targetTid);
		}
		for(MagQCBinManifest.Selection s : MagQCBinManifest.decodeSelections(b,true)){
			if(excludedTids.contains(s.tid)){
				throw new RuntimeException("Replay selects an explicitly excluded foreign tid: "+s.tid);
			}
		}
		// The legacy per-bin index and the shred table alone cannot apply individual gene cuts.
		// A paired contribution index is the exact breakpoint source; its dispatch is row-local so
		// ordinary/HQ rows in the SAME manifest continue through refcdstable=.
		if(referenceCdsIndexFile!=null && !"-".equals(b.breakpoints)){
			throw new RuntimeException("referencecdsindex= does not support breakpoint bins: "+b.id+
				" breakpoints="+b.breakpoints);
		}
		if(referenceCdsContributionIndexFile==null && !"-".equals(b.breakpoints)){
			throw new RuntimeException("Breakpoint bins require referencecdscontributionindex=: "+b.id+
				" breakpoints="+b.breakpoints);
		}
		if(!"-".equals(b.breakpoints) && !"perfect_unshredded".equals(b.spikeClass)){
			throw new RuntimeException("Breakpoint bins require spike_class=perfect_unshredded: "+b.id+
				" spike_class="+b.spikeClass);
		}
		if(!"-".equals(b.breakpoints) &&
				(!"-".equals(b.foreignContigs) || !"-".equals(b.contaminantTids))){
			throw new RuntimeException("Breakpoint clone may not contain foreign contigs: "+b.id);
		}
		if("perfect_unshredded".equals(b.spikeClass) && "-".equals(m.secondarySourcePath)){
			throw new RuntimeException("Manifest has perfect_unshredded row but no secondary source provenance: "+b.id);
		}
		if("perfect_unshredded".equals(b.spikeClass) && referenceCdsContributionIndexFile==null &&
				subnetLabels==SUBNETLABELS_GENE){
			throw new RuntimeException("perfect_unshredded gene labels require "+
				"referencecdscontributionindex=: "+b.id);
		}
		if(referenceCdsContributionIndexFile!=null && refCdsTableFile==null &&
				!"perfect_unshredded".equals(b.spikeClass) && subnetLabels==SUBNETLABELS_GENE){
			throw new RuntimeException("Non-unshredded gene-label replay requires refcdstable= in addition "+
				"to referencecdscontributionindex=: "+b.id);
		}
	}

	/** One fully formatted replay row, produced by exactly one worker. */
	private static final class ReplayOutput {
		final ByteBuilder dense,subnet,agg,ids,sidecar;
		ReplayOutput(ByteBuilder dense_,ByteBuilder subnet_,ByteBuilder agg_,ByteBuilder ids_,ByteBuilder sidecar_){
			dense=dense_; subnet=subnet_; agg=agg_; ids=ids_; sidecar=sidecar_;
		}
	}

	/** Optional per-worker replay measurements.  All scopes are nanoseconds and are
	 * accumulated only when replaydiag=t; the normal replay path has no timer calls. */
	private static final class ReplayDiagnostics {
		long sourceRowGetNanos;
		long replayReconstructionNanos;
		long contextNanos;
		long expectedCopyNanos;
		long subnetInputNanos,applyInputNanos;
		long subnetInferenceNanos;
		long formatAggRowNanos;
		long formatAggRowExclusiveNanos;
		long feedForwardCalls;
	}

	/** Per-worker row scratch; buffers are reused only after the ordered publisher copies them. */
	private static final class ReplayScratch {
		ReplayScratch(int numFam,int denseWidth,int subnetWidth,int aggWidth,boolean dense,boolean subnet,boolean agg,
				boolean ids,boolean sidecar){
			labels=new double[2]; fam=new int[numFam]; glob=new double[NUM_GLOBALS]; accumulator=new Agg(); accumulator.setFam(fam);
			denseOut=dense ? new ByteBuilder(denseWidth*4+64) : null;
			subnetOut=subnet ? new ByteBuilder(subnetWidth*4+64) : null;
			aggOut=agg ? new ByteBuilder(aggWidth*4+64) : null;
			idsOut=ids ? new ByteBuilder(64) : null;
			sidecarOut=sidecar ? new ByteBuilder(192) : null;
		}
		void reset(){
			Arrays.fill(labels,0); Arrays.fill(fam,0); Arrays.fill(glob,0); accumulator.reset();
			if(denseOut!=null){denseOut.clear();} if(subnetOut!=null){subnetOut.clear();}
			if(aggOut!=null){aggOut.clear();} if(idsOut!=null){idsOut.clear();} if(sidecarOut!=null){sidecarOut.clear();}
		}
		final double[] labels;
		final int[] fam;
		final double[] glob;
		final Agg accumulator;
		final ByteBuilder denseOut,subnetOut,aggOut,idsOut,sidecarOut;
	}

	/**
	 * Candidate-only state for one row in a subnet-major replay block.  The row's
	 * source reconstruction and all output buffers remain private until the
	 * ordered coordinator publishes them, so block evaluation cannot alter
	 * output order or alias another row's scratch.
	 */
	private static final class ReplayBlockRow {
		final ReplayScratch output;
		final int[] cleanFam;
		final int[] ncServe=new int[NCRNA_OBS], ncObs=new int[NCRNA_OBS];
		final int[] antiServe=new int[TRNA_ANTICODON_OBS], antiObs=new int[TRNA_ANTICODON_OBS];
		final double[] rowCtx=new double[CTX_N];
		final float[] quantizedRowCtx=new float[CTX_N];
		final int[] observedByItem;
		final float[] subnetObservedInputs, subnetSharedInputs;
		final double[] subnetOutputs, subnetDerived;
		final long[] subnetObs;
		final double[] subnetPred;
		String binId;
		long contextNanos, expectedCopyNanos, subnetInputNanos, subnetInferenceNanos;
		int phylumIdx, domainIdx;
		long sumObs;
		double sumPred;

		ReplayBlockRow(MagQCVectorMaker owner,int subnetCount,boolean hasSubnet,boolean hasAgg){
			output=new ReplayScratch(owner.numFam,owner.numInputs,owner.subnetInputs,owner.numAggInputs,
				owner.emitDense,hasSubnet,hasAgg,owner.aggIds,
				owner.subnetLabels==SUBNETLABELS_GENE && hasSubnet);
			cleanFam=new int[owner.numFam];
			final int expectedWidth=owner.expectedCopy==null ? 0 : owner.expectedCopy.itemWidth();
			observedByItem=expectedWidth>0 ? new int[expectedWidth] : null;
			subnetObservedInputs=new float[owner.numFam+NCRNA_OBS+TRNA_ANTICODON_OBS];
			subnetSharedInputs=new float[owner.numPhyla+CTX_N+DOMAINS];
			subnetOutputs=new double[subnetCount*4]; subnetDerived=new double[subnetCount*2];
			subnetObs=new long[subnetCount]; subnetPred=new double[subnetCount];
		}

		void reset(){
			output.reset(); sumObs=0; sumPred=0;
			Arrays.fill(subnetOutputs,0); Arrays.fill(subnetDerived,0);
			Arrays.fill(subnetObs,0); Arrays.fill(subnetPred,0);
			if(observedByItem!=null){Arrays.fill(observedByItem,0);}
			binId=null; contextNanos=0; expectedCopyNanos=0;
			subnetInputNanos=0; subnetInferenceNanos=0;
		}
	}

	/** Reusable bounded candidate block; rows are emitted in their original ordinal order. */
	private static final class ReplayBlockScratch {
		final ReplayBlockRow[] rows;
		ReplayBlockScratch(MagQCVectorMaker owner,int capacity,boolean hasSubnet,boolean hasAgg){
			rows=new ReplayBlockRow[capacity];
			final int subnetCount=owner.aggSubnets==null ? 0 : owner.aggSubnets.size();
			for(int i=0; i<capacity; i++){rows[i]=new ReplayBlockRow(owner,subnetCount,hasSubnet,hasAgg);}
		}
	}

	/** Ordered publication gate: workers compute concurrently, but output remains manifest order. */
	private static final class ReplayCoordinator {
		private final ByteStreamWriter[] dense,subnet,agg,ids,sidecar;
		private final int partCount;
		private final String tag;
		private final boolean diagnostics;
		private final long startedNano=System.nanoTime();
		private final Object lock=new Object();
		private int next;
		private Throwable failure;
		private long anticodonRows,anticodonMissingRows;
		private long sourceRowGetNanos,replayReconstructionNanos,contextNanos,expectedCopyNanos,subnetInferenceNanos;
		private long subnetInputNanos,applyInputNanos;
		private long formatAggRowNanos,formatAggRowExclusiveNanos,writerEnqueueNanos;
		private long feedForwardCalls;
		ReplayCoordinator(ByteStreamWriter[] dense_,ByteStreamWriter[] subnet_,ByteStreamWriter[] agg_,
				ByteStreamWriter[] ids_,ByteStreamWriter[] sidecar_,int partCount_,String tag_,boolean diagnostics_){
			dense=dense_; subnet=subnet_; agg=agg_; ids=ids_; sidecar=sidecar_; partCount=partCount_; tag=tag_; diagnostics=diagnostics_;
		}
		/** Publishes after the ordered wait; writer_enqueue excludes lock wait and final
		 * compressor/drain time, which occurs after this method returns. */
		void publish(int ordinal,ReplayOutput out){
			synchronized(lock){
				while(failure==null && ordinal!=next){
					try{lock.wait();}
					catch(InterruptedException e){Thread.currentThread().interrupt(); throw new RuntimeException("Replay worker interrupted",e);}
				}
				if(failure!=null){throw new RuntimeException("Replay publication aborted",failure);}
				final int part=ordinal%partCount;
				final long writerStart=(diagnosticsEnabled() ? System.nanoTime() : 0);
				if(dense!=null){dense[part].print(out.dense);}
				if(subnet!=null){subnet[part].print(out.subnet);}
				if(agg!=null){agg[part].print(out.agg);}
				if(ids!=null){ids[part].print(out.ids);}
				if(sidecar!=null){sidecar[part].print(out.sidecar);}
				if(diagnosticsEnabled()){writerEnqueueNanos+=System.nanoTime()-writerStart;}
				next++;
				if(next%10000==0){
					final double seconds=(System.nanoTime()-startedNano)/1e9;
					System.err.println(tag+": replay progress rows="+next+" elapsed_seconds="+seconds+" rows_per_second="+(next/Math.max(seconds,1e-9)));
				}
				lock.notifyAll();
			}
		}
		void fail(Throwable t){
			synchronized(lock){if(failure==null){failure=t;} lock.notifyAll();}
		}
		Throwable failure(){synchronized(lock){return failure;}}
		int published(){synchronized(lock){return next;}}
		void addDiagnostics(long rows,long missing){synchronized(lock){anticodonRows+=rows; anticodonMissingRows+=missing;}}
		void addDiagnostics(ReplayDiagnostics d){
			if(d==null){return;}
			synchronized(lock){
				sourceRowGetNanos+=d.sourceRowGetNanos; replayReconstructionNanos+=d.replayReconstructionNanos;
				contextNanos+=d.contextNanos;
				expectedCopyNanos+=d.expectedCopyNanos; subnetInferenceNanos+=d.subnetInferenceNanos;
				subnetInputNanos+=d.subnetInputNanos; applyInputNanos+=d.applyInputNanos;
				formatAggRowNanos+=d.formatAggRowNanos; formatAggRowExclusiveNanos+=d.formatAggRowExclusiveNanos;
				feedForwardCalls+=d.feedForwardCalls;
			}
		}
		boolean diagnosticsEnabled(){return diagnostics;}
		long sourceRowGetNanos(){synchronized(lock){return sourceRowGetNanos;}}
		long replayReconstructionNanos(){synchronized(lock){return replayReconstructionNanos;}}
		long contextNanos(){synchronized(lock){return contextNanos;}}
		long expectedCopyNanos(){synchronized(lock){return expectedCopyNanos;}}
		long subnetInputNanos(){synchronized(lock){return subnetInputNanos;}}
		long applyInputNanos(){synchronized(lock){return applyInputNanos;}}
		long subnetInferenceNanos(){synchronized(lock){return subnetInferenceNanos;}}
		long formatAggRowNanos(){synchronized(lock){return formatAggRowNanos;}}
		long formatAggRowExclusiveNanos(){synchronized(lock){return formatAggRowExclusiveNanos;}}
		long writerEnqueueNanos(){synchronized(lock){return writerEnqueueNanos;}}
		long feedForwardCalls(){synchronized(lock){return feedForwardCalls;}}
		long anticodonRows(){synchronized(lock){return anticodonRows;}}
		long anticodonMissingRows(){synchronized(lock){return anticodonMissingRows;}}
		boolean hasSubnet(){return subnet!=null;}
		boolean hasAgg(){return agg!=null;}
	}

	/** Native ProcessThread-shaped replay worker with private maker scratch and file cursor. */
	private final class ReplayProcessThread extends Thread {
		private MagQCVectorMaker worker;
		private ReplayScratch scratch;
		private ReplayDiagnostics diagnostics;
		private final MagQCBinReplayStore store;
		private MagQCBinReplayStore.Reader reader;
		private final ArrayList<MagQCBinManifest.Bin> rows;
		private final String tag;
		private final int first,last,ordinal,threadCount;
		private final ReplayCoordinator coordinator;
		ReplayProcessThread(MagQCBinReplayStore store_,ArrayList<MagQCBinManifest.Bin> rows_,String tag_,
				int first_,int last_,int ordinal_,int threadCount_,ReplayCoordinator coordinator_){
			store=store_; rows=rows_; tag=tag_; first=first_; last=last_;
			ordinal=ordinal_; threadCount=threadCount_; coordinator=coordinator_;
			setName("MagQC-replay-"+tag+"-"+ordinal);
		}
		@Override public void run(){
			try{
				// Keep ownership and allocation inside the worker thread.  This includes
				// every bundle InferenceNet clone, per-worker arrays/IntLists, and all
				// reusable output scratch; the constructor only records immutable inputs.
				worker=copyForReplayThread();
				diagnostics=(coordinator.diagnosticsEnabled() ? new ReplayDiagnostics() : null);
				worker.replayDiag=diagnostics;
				final ReplayBlockScratch block=(replayBlockRows>0 ?
					new ReplayBlockScratch(worker,replayBlockRows,coordinator.hasSubnet(),coordinator.hasAgg()) : null);
				if(block==null){
					scratch=new ReplayScratch(worker.numFam,worker.numInputs,worker.subnetInputs,worker.numAggInputs,
						worker.emitDense,coordinator.hasSubnet(),coordinator.hasAgg(),worker.aggIds,
						worker.subnetLabels==SUBNETLABELS_GENE && coordinator.hasSubnet());
				}
				reader=(store==null ? null : store.openReader());
				if(block==null){
					for(int row=first+ordinal; row<last; row+=threadCount){
						if(coordinator.failure()!=null){return;}
						final long sourceRowStart=(diagnostics==null ? 0 : System.nanoTime());
						final MagQCBinManifest.Bin b=reader==null ? rows.get(row) : reader.get(tag,row);
						if(diagnostics!=null){diagnostics.sourceRowGetNanos+=System.nanoTime()-sourceRowStart;}
						coordinator.publish(row-first,worker.formatReplayOutput(b,tag,row,first,scratch,coordinator.hasSubnet(),coordinator.hasAgg()));
					}
				}else{
					final int ownedFirst=first+ordinal;
					for(int blockFirst=ownedFirst; blockFirst<last; blockFirst+=threadCount*replayBlockRows){
						int count=0;
						while(count<replayBlockRows && blockFirst+count*threadCount<last){count++;}
						for(int i=0; i<count; i++){
							if(coordinator.failure()!=null){return;}
							final int rowIndex=blockFirst+i*threadCount;
							final long sourceRowStart=(diagnostics==null ? 0 : System.nanoTime());
							final MagQCBinManifest.Bin b=reader==null ? rows.get(rowIndex) : reader.get(tag,rowIndex);
							if(diagnostics!=null){diagnostics.sourceRowGetNanos+=System.nanoTime()-sourceRowStart;}
							worker.prepareReplayBlockRow(block.rows[i],b,tag,rowIndex,first);
						}
						if(coordinator.hasAgg()){worker.evaluateReplayBlock(block,count);}
						for(int i=0; i<count; i++){
							final int rowIndex=blockFirst+i*threadCount;
							coordinator.publish(rowIndex-first,worker.finishReplayBlockRow(block.rows[i]));
						}
					}
				}
			}catch(Throwable t){coordinator.fail(t);}
			finally{
				if(diagnostics!=null){coordinator.addDiagnostics(diagnostics);}
				if(worker!=null){coordinator.addDiagnostics(worker.anticodonRows,worker.anticodonMissingRows);}
				if(reader!=null){try{reader.close();}catch(Throwable t){coordinator.fail(t);}}
			}
		}
	}

	private ReplayOutput formatReplayOutput(MagQCBinManifest.Bin b,String tag,int rowIndex,int firstRow,ReplayScratch scratch,
			boolean hasSubnet,boolean hasAgg){
		if(binModel==null && !b.modelRows.equals(tag+":"+rowIndex)){
			throw new RuntimeException("Manifest model_rows/order mismatch for "+b.id+": expected "+tag+":"+rowIndex+" got "+b.modelRows);
		}
		scratch.reset();
		final double[] labels=scratch.labels;
		final int[] fam=scratch.fam;
		final double[] glob=scratch.glob;
		final Agg agg=scratch.accumulator;
		final ReplayDiagnostics diag=replayDiag;
		final long reconstructionStart=(diag==null ? 0 : System.nanoTime());
		final int targetPhylumIdx=replayBin(b,fam,glob,labels,agg);
		if(targetPhylumIdx<0){throw new RuntimeException("Manifest row rejected during replay: "+b.id);}
		final ReferenceCdsSurvivalLabelReader.Counts reference=exactReplayCounts(b);
		if(reference!=null){labels[0]=reference.completeness; labels[1]=reference.contamination;}
		if(diag!=null){diag.replayReconstructionNanos+=System.nanoTime()-reconstructionStart;}
		final ByteBuilder dense=scratch.denseOut;
		final ByteBuilder subnet=scratch.subnetOut;
		final ByteBuilder aggOut=scratch.aggOut;
		if(dense!=null){formatRow(dense,fam,glob,targetPhylumIdx,labels,agg);}
		if(subnet!=null){
			if(subnetLabels==SUBNETLABELS_GENE){
				if(subnetFamset){formatFamsetGeneRow(subnet,fam,glob,targetPhylumIdx,labels,agg);}
				else if(subnetAnticodon || subnetRrna){formatSpecialGeneRow(subnet,glob,targetPhylumIdx,labels,agg);}
				else{formatNcrnaGeneRow(subnet,glob,targetPhylumIdx,labels,agg);}
			}else{
				if(subnetFamset){formatFamsetRow(subnet,glob,targetPhylumIdx,agg);}
				else if(subnetAnticodon || subnetRrna){formatSpecialRow(subnet,glob,targetPhylumIdx,agg);}
				else{formatNcrnaRow(subnet,glob,targetPhylumIdx,agg);}
			}
		}
		if(aggOut!=null){
			final long formatStart=(diag==null ? 0 : System.nanoTime());
			final long nestedStart=(diag==null ? 0 : diag.contextNanos+diag.expectedCopyNanos+diag.subnetInputNanos+diag.subnetInferenceNanos);
			formatAggRow(aggOut,fam,glob,targetPhylumIdx,labels,agg);
			if(diag!=null){
				final long formatElapsed=System.nanoTime()-formatStart;
				final long nestedElapsed=(diag.contextNanos+diag.expectedCopyNanos+diag.subnetInputNanos+diag.subnetInferenceNanos)-nestedStart;
				diag.formatAggRowNanos+=formatElapsed;
				diag.formatAggRowExclusiveNanos+=Math.max(0,formatElapsed-nestedElapsed);
			}
		}
		final ByteBuilder ids=scratch.idsOut;
		if(ids!=null){ids.append((rowIndex-firstRow)/replayParts).tab().append(b.id).nl();}
		final ByteBuilder sidecar=scratch.sidecarOut;
		if(sidecar!=null){
			if(reference==null){throw new IllegalStateException("Gene replay sidecar requires exact reference labels for "+b.id);}
			appendSubnetGeneReplaySidecarRow(sidecar,rowIndex,b,reference);
		}
		return new ReplayOutput(dense,subnet,aggOut,ids,sidecar);
	}

	/** Prepares one row for the unreleased bounded subnet-major candidate. */
	private void prepareReplayBlockRow(ReplayBlockRow row, MagQCBinManifest.Bin b, String tag, int rowIndex,
			int firstRow){
		row.reset();
		row.binId=b.id;
		if(binModel==null && !b.modelRows.equals(tag+":"+rowIndex)){
			throw new RuntimeException("Manifest model_rows/order mismatch for "+b.id+": expected "+tag+":"+rowIndex+" got "+b.modelRows);
		}
		final double[] labels=row.output.labels;
		final int[] fam=row.output.fam;
		final double[] glob=row.output.glob;
		final Agg agg=row.output.accumulator;
		final ReplayDiagnostics diag=replayDiag;
		final long reconstructionStart=(diag==null ? 0 : System.nanoTime());
		final int targetPhylumIdx=replayBin(b,fam,glob,labels,agg);
		if(targetPhylumIdx<0){throw new RuntimeException("Manifest row rejected during replay: "+b.id);}
		final ReferenceCdsSurvivalLabelReader.Counts reference=exactReplayCounts(b);
		if(reference!=null){labels[0]=reference.completeness; labels[1]=reference.contamination;}
		if(diag!=null){diag.replayReconstructionNanos+=System.nanoTime()-reconstructionStart;}
		if(row.output.denseOut!=null){formatRow(row.output.denseOut,fam,glob,targetPhylumIdx,labels,agg);}
		if(row.output.subnetOut!=null){
			if(subnetLabels==SUBNETLABELS_GENE){
				if(subnetFamset){formatFamsetGeneRow(row.output.subnetOut,fam,glob,targetPhylumIdx,labels,agg);}
				else if(subnetAnticodon || subnetRrna){formatSpecialGeneRow(row.output.subnetOut,glob,targetPhylumIdx,labels,agg);}
				else{formatNcrnaGeneRow(row.output.subnetOut,glob,targetPhylumIdx,labels,agg);}
			}else{
				if(subnetFamset){formatFamsetRow(row.output.subnetOut,glob,targetPhylumIdx,agg);}
				else if(subnetAnticodon || subnetRrna){formatSpecialRow(row.output.subnetOut,glob,targetPhylumIdx,agg);}
				else{formatNcrnaRow(row.output.subnetOut,glob,targetPhylumIdx,agg);}
			}
		}
		row.phylumIdx=targetPhylumIdx;
		row.domainIdx=domainIdxOf(lastTarget);
		System.arraycopy(lastNcServe,0,row.ncServe,0,NCRNA_OBS);
		System.arraycopy(lastNcObs,0,row.ncObs,0,NCRNA_OBS);
		System.arraycopy(lastAntiServe,0,row.antiServe,0,TRNA_ANTICODON_OBS);
		System.arraycopy(lastAntiObs,0,row.antiObs,0,TRNA_ANTICODON_OBS);
		if(cleanFamBuf!=null){System.arraycopy(cleanFamBuf,0,row.cleanFam,0,numFam);}
		else{Arrays.fill(row.cleanFam,0);}
		final long contextStart=(diag==null ? 0 : System.nanoTime());
		computeContext(row.rowCtx,glob,agg);
		for(int i=0; i<CTX_N; i++){row.quantizedRowCtx[i]=Float.parseFloat(fmt(row.rowCtx[i]));}
		final int[] inputFam=(aggObsServe ? fam : row.cleanFam);
		final int[] inputNc=(aggObsServe ? row.ncServe : row.ncObs);
		final int[] inputAnti=(aggObsServe ? row.antiServe : row.antiObs);
		for(int i=0; i<row.subnetObservedInputs.length; i++){row.subnetObservedInputs[i]=0;}
		for(int r=0; r<numFam; r++){row.subnetObservedInputs[r]=inputFam[r];}
		for(int i=0; i<NCRNA_OBS; i++){row.subnetObservedInputs[numFam+i]=inputNc[i];}
		for(int i=0; i<TRNA_ANTICODON_OBS; i++){row.subnetObservedInputs[numFam+NCRNA_OBS+i]=inputAnti[i];}
		Arrays.fill(row.subnetSharedInputs,0);
		if(row.phylumIdx>=0 && row.phylumIdx<numPhyla){row.subnetSharedInputs[row.phylumIdx]=1;}
		System.arraycopy(row.quantizedRowCtx,0,row.subnetSharedInputs,numPhyla,CTX_N);
		if(row.domainIdx>=0 && row.domainIdx<DOMAINS){row.subnetSharedInputs[numPhyla+CTX_N+row.domainIdx]=1;}
		if(row.observedByItem!=null){
			for(int r=0; r<numFam; r++){row.observedByItem[r]=fam[r];}
			for(int i=0; i<NCRNA_OBS; i++){row.observedByItem[numFam+i]=row.ncServe[i];}
			if(row.observedByItem.length==numFam+NCRNA_OBS+TRNA_ANTICODON_OBS){
				for(int i=0; i<TRNA_ANTICODON_OBS; i++){
					row.observedByItem[numFam+NCRNA_OBS+i]=row.antiServe[i];
				}
			}
		}
		if(diag!=null){
			row.contextNanos=System.nanoTime()-contextStart;
			diag.contextNanos+=row.contextNanos;
		}
		if(row.output.idsOut!=null){row.output.idsOut.append((rowIndex-firstRow)/replayParts).tab().append(b.id).nl();}
		if(row.output.sidecarOut!=null){
			if(reference==null){throw new IllegalStateException("Gene replay sidecar requires exact reference labels for "+b.id);}
			appendSubnetGeneReplaySidecarRow(row.output.sidecarOut,rowIndex,b,reference);
		}
	}

	/** Evaluates each frozen subnet across the bounded row block, without changing row order. */
	private void evaluateReplayBlock(ReplayBlockScratch block, int count){
		if(aggSubnets==null || count<1){return;}
		final boolean six=(subnetFeatures==SUBNETFEATURES_SIX);
		for(int si=0; si<aggSubnets.size(); si++){
			final AggSubnet s=aggSubnets.get(si);
			for(int ri=0; ri<count; ri++){
				final ReplayBlockRow row=block.rows[ri];
				final ReplayDiagnostics diag=replayDiag;
				final long inputStart=(diag==null ? 0 : System.nanoTime());
				final long obsTotal=prepareReplayBlockSubnetInput(row,s);
				final ml.CellNet net=s.netForCurrentThread();
				final long applyStart=(diag==null ? 0 : System.nanoTime());
				net.applyInput(s.buf);
				if(diag!=null){
					final long applied=System.nanoTime();
					diag.applyInputNanos+=applied-applyStart;
					row.subnetInputNanos+=applied-inputStart;
					diag.subnetInputNanos+=applied-inputStart;
				}
				final long inferenceStart=(diag==null ? 0 : System.nanoTime());
				net.feedForward();
				if(diag!=null){
					final long inferenceElapsed=System.nanoTime()-inferenceStart;
					row.subnetInferenceNanos+=inferenceElapsed;
					diag.feedForwardCalls++; diag.subnetInferenceNanos+=inferenceElapsed;
				}
				if(six){
					for(int k=0; k<4; k++){
						final float o=net.getOutput(k);
						if(!Float.isFinite(o)){throw new RuntimeException("subnetfeatures=six: non-finite learned output "+sixNames[k]
							+" from subnet "+s.name+" (order "+si+") for bin "+row.binId);}
						row.subnetOutputs[si*4+k]=o;
					}
					final long expectedCopyStart=(diag==null ? 0 : System.nanoTime());
					final MagQCExpectedCopyFeatures.Result derived=expectedCopy.compute(si,row.observedByItem);
					if(diag!=null){
						final long expectedCopyElapsed=System.nanoTime()-expectedCopyStart;
						row.expectedCopyNanos+=expectedCopyElapsed;
						diag.expectedCopyNanos+=expectedCopyElapsed;
					}
					row.subnetDerived[si*2]=derived.observedExpected;
					row.subnetDerived[si*2+1]=derived.excessExpected;
				}else{
					final double pred=net.getOutput(0);
					row.subnetObs[si]=obsTotal;
					row.subnetPred[si]=pred;
					row.sumObs+=obsTotal;
					row.sumPred+=Math.max(0,pred);
				}
			}
		}
		for(int ri=0; ri<count; ri++){formatReplayBlockAggRow(block.rows[ri],six);}
	}

	/** Builds exactly one subnet input from a prepared row, reusing the worker's bounded subnet buffer. */
	private long prepareReplayBlockSubnetInput(ReplayBlockRow row, AggSubnet s){
		final int[] famArr=(aggObsServe ? row.output.fam : row.cleanFam);
		final int[] ncArr=(aggObsServe ? row.ncServe : row.ncObs);
		final int[] antiArr=(aggObsServe ? row.antiServe : row.antiObs);
		final float[] observed=row.subnetObservedInputs;
		int p=0; long obsTotal=0;
		if("ncrna".equals(s.type)){
			System.arraycopy(observed,numFam,s.buf,0,NCRNA_OBS); p=NCRNA_OBS;
			if(subnetFeatures!=SUBNETFEATURES_SIX){for(int i=0; i<NCRNA_OBS; i++){obsTotal+=ncArr[i];}}
		}else if("rrna".equals(s.type)){
			System.arraycopy(observed,numFam,s.buf,0,RRNA_OBS); p=RRNA_OBS;
			if(subnetFeatures!=SUBNETFEATURES_SIX){for(int i=0; i<RRNA_OBS; i++){obsTotal+=ncArr[i];}}
		}else if("trna_anticodon".equals(s.type)){
			System.arraycopy(observed,numFam+NCRNA_OBS,s.buf,0,TRNA_ANTICODON_OBS); p=TRNA_ANTICODON_OBS;
			if(subnetFeatures!=SUBNETFEATURES_SIX){for(int i=0; i<TRNA_ANTICODON_OBS; i++){obsTotal+=antiArr[i];}}
		}else if("famset".equals(s.type)){
			for(int r : s.ranks){s.buf[p++]=observed[r]; if(subnetFeatures!=SUBNETFEATURES_SIX){obsTotal+=famArr[r];}}
		}else{throw new RuntimeException("Unknown aggregator subnet type "+s.type+" for "+s.name);}
		System.arraycopy(row.subnetSharedInputs,0,s.buf,p,row.subnetSharedInputs.length); p+=row.subnetSharedInputs.length;
		assert(p==s.buf.length) : s.name+": filled "+p+" of "+s.buf.length;
		return obsTotal;
	}

	/** Emits one prepared row after all subnet evaluations for its block are complete. */
	private void formatReplayBlockAggRow(ReplayBlockRow row, boolean six){
		final ByteBuilder bb=row.output.aggOut;
		if(bb==null){return;}
		final ReplayDiagnostics diag=replayDiag;
		final long serializationStart=(diag==null ? 0 : System.nanoTime());
		if(six){
			for(int si=0; si<aggSubnets.size(); si++){
				for(int k=0; k<4; k++){appendFmt(bb,row.subnetOutputs[si*4+k]); bb.tab();}
				appendFmt(bb,row.subnetDerived[si*2]); bb.tab();
				appendFmt(bb,row.subnetDerived[si*2+1]); bb.tab();
			}
		}else{
			for(int si=0; si<aggSubnets.size(); si++){
				final long obsTotal=row.subnetObs[si]; final double pred=row.subnetPred[si];
				appendFmt(bb,Math.min(RATIO_CAP,obsTotal/Math.max(0.5,pred))); bb.tab();
				appendFmt(bb,log2(1+obsTotal)); bb.tab();
				appendFmt(bb,log2(1+Math.max(0,pred))); bb.tab();
				bb.append(obsTotal==0 ? '1' : '0'); bb.tab();
			}
			appendFmt(bb,Math.min(RATIO_CAP,row.sumObs/Math.max(1,row.sumPred))); bb.tab();
		}
		for(int r : denseRanks){appendTwo(bb,row.output.fam[r]);}
		for(int i=0; i<numPhyla; i++){bb.append(i==row.phylumIdx ? '1' : '0'); bb.tab();}
		appendContextDomain(bb,row.rowCtx,row.domainIdx);
		final int[] ncArr=(aggObsServe ? row.ncServe : row.ncObs);
		for(int i=0; i<NCRNA_OBS; i++){appendCountTwo(bb,ncArr[i]);}
		appendFmt(bb,row.output.glob[7]);
		bb.tab(); appendFmt(bb,row.output.labels[0]); bb.tab(); appendFmt(bb,row.output.labels[1]); bb.nl();
		if(diag!=null){
			// Match row-major formatAggRowNanos with disjoint measured scopes.  Preparation and
			// inference happen subnet-major before this method, so inclusive time is synthesized
			// from this row's context/input/inference/expected-copy scopes plus local serialization;
			// never use a wall timer spanning preparation of other rows in the block.
			final long serialization=System.nanoTime()-serializationStart;
			final long nested=row.contextNanos+row.expectedCopyNanos+row.subnetInputNanos+row.subnetInferenceNanos;
			final long inclusive=nested+serialization;
			diag.formatAggRowNanos+=inclusive;
			diag.formatAggRowExclusiveNanos+=serialization;
		}
	}

	private ReplayOutput finishReplayBlockRow(ReplayBlockRow row){
		return new ReplayOutput(row.output.denseOut,row.output.subnetOut,row.output.aggOut,row.output.idsOut,row.output.sidecarOut);
	}

	/** Names one native output part without creating an uncompressed aggregate intermediate. */
	private static String replayPartPath(String file,int part,int parts){
		if(parts==1){return file;}
		final String marker=String.format(".part%02dof%02d",part,parts);
		if(file.endsWith(".vec.gz")){return file.substring(0,file.length()-7)+marker+".vec.gz";}
		if(file.endsWith(".ids.tsv.gz")){return file.substring(0,file.length()-11)+marker+".ids.tsv.gz";}
		if(file.endsWith(".gz")){return file.substring(0,file.length()-3)+marker+".gz";}
		return file+marker;
	}

	/** Identity companions use the stable D195 .ids.tsv.gz suffix for vector .vec.gz bases. */
	private static String replayIdentityPath(String aggFile){
		if(aggFile.endsWith(".vec.gz")){return aggFile.substring(0,aggFile.length()-7)+".ids.tsv.gz";}
		return aggFile+".ids.tsv";
	}

	private String replayManifestSha256(){
		if(binManifestInFile==null){throw new IllegalStateException("Replay manifest digest requested without binmanifestin=");}
		if(replayManifestSha256==null){replayManifestSha256=MagQCBinManifest.sha256File(binManifestInFile);}
		return replayManifestSha256;
	}

	private static ByteStreamWriter[] openReplayWriters(String file,int parts){
		final ByteStreamWriter[] writers=new ByteStreamWriter[parts];
		try{
			for(int i=0; i<parts; i++){
				writers[i]=new ByteStreamWriter(replayPartPath(file,i,parts),true,false,true);
				writers[i].start();
			}
			return writers;
		}catch(RuntimeException|Error failure){
			if(writers!=null){for(ByteStreamWriter writer:writers){if(writer!=null){writer.poisonAndWait();}}}
			throw failure;
		}
	}

	private static boolean poisonReplayWriters(ByteStreamWriter[] writers){
		boolean failed=false;
		if(writers!=null){for(ByteStreamWriter writer:writers){if(writer!=null && writer.poisonAndWait()){failed=true;}}}
		return failed;
	}

	private static void appendReplayPartHeader(ByteBuilder header,String tag,int part,int parts,String manifestSha80,
			boolean includeSplit,boolean includeManifestHash,boolean shuffled,long seed){
		if(parts>1 || shuffled){
			if(includeManifestHash){header.append("#bin_manifest_sha80\t").append(manifestSha80).nl();}
			if(shuffled){header.append("#replay_shuffle\tseeded_fisher_yates\n").append("#replay_seed\t").append(seed).nl();}
			if(includeSplit){header.append("#split\t").append(tag).nl();}
		}
		if(parts>1){
			header.append("#part_index\t").append(part).nl();
			header.append("#part_count\t").append(parts).nl();
			header.append("#row_partition\temitted_row_index_modulo\t").append(parts).nl();
		}
	}

	/** Replays exact manifest selections; no sampler/RNG calls occur. */
	private void writeReplaySet(String file, String subnetFile, String aggFile,
			MagQCBinManifest.Manifest manifest, String tag,MagQCBinReplayStore store,int start,int end){
		final ArrayList<MagQCBinManifest.Bin> rows=new ArrayList<MagQCBinManifest.Bin>();
		if(store==null){for(MagQCBinManifest.Bin b : manifest.bins){if(tag.equals(b.split)){rows.add(b);}}}
		final int rowCount=(store==null ? rows.size() : store.size(tag));
		final int first=(store==null ? 0 : start),last=(store==null ? rowCount : end);
		final String replayManifestSha80=(binManifestInFile==null ? "-" : replayManifestSha256().substring(44));
		if(replayParts>1){
			// D195 parts are required to be BGZF, with bounded native compression rather than one
			// full thread pool per sink.  The mixer uses the same explicit policy for its outputs.
			ReadWrite.USE_UNBGZIP=false; ReadWrite.ALLOW_NATIVE_BGZF=true; ReadWrite.USE_BGZF=true;
			ReadWrite.FORCE_BGZIP=true; ReadWrite.PREFER_NATIVE_BGZF_OUT=true; ReadWrite.ZIPLEVEL=4;
			ReadWrite.setZipThreads(1); BgzfSettings.USE_MULTITHREADED_BGZF=false;
		}
		final ByteBuilder binding=aggIds && aggFile!=null ?
			compositeReplayBinding(aggFile, file, subnetFile, tag, first, last) : null;
		// Setup phase (Yoimiya 2026-09-09 review): a LATER writer's own construction/header can fail
		// after EARLIER writers already started their non-daemon I/O threads (same class of hazard as
		// the makeBin try/finally fix, 2026-09-03) -- openSubnetGeneReplaySidecar in particular reads
		// external files to hash into its header and can throw. Every writer that DID open must be
		// closed before the failure propagates, or it leaks a live thread and hangs the JVM at exit.
		ByteStreamWriter[] bsw=null, sbsw=null, absw=null, gsw=null, isw=null;
		final ByteBuilder bb, bb2, bb3;
		try{
			bsw=(emitDense ? openReplayWriters(file,replayParts) : null);
			bb=(bsw==null ? null : new ByteBuilder(numInputs*4+128));
			if(bsw!=null){
				for(int i=0; i<replayParts; i++){
					bb.append("#dims\t").append(numInputs).append("\t2\t0").nl();
					appendReplayPartHeader(bb,tag,i,replayParts,replayManifestSha80,true,true,replayShuffle,replaySeed); bsw[i].print(bb); bb.clear();
				}
			}
			sbsw=(subnetFile==null ? null : openReplayWriters(subnetFile,replayParts));
			bb2=(sbsw==null ? null : new ByteBuilder(subnetInputs*4+128));
			if(sbsw!=null){
				for(int i=0; i<replayParts; i++){
					bb2.append("#dims\t").append(subnetInputs).append(subnetLabels==SUBNETLABELS_GENE ? "\t2\t0" : "\t1\t0").nl();
					appendReplayPartHeader(bb2,tag,i,replayParts,replayManifestSha80,true,true,replayShuffle,replaySeed); sbsw[i].print(bb2); bb2.clear();
				}
			}
			absw=(aggFile==null ? null : openReplayWriters(aggFile,replayParts));
			bb3=(absw==null ? null : new ByteBuilder(numAggInputs*4+256));
			if(absw!=null){
				for(int i=0; i<replayParts; i++){
					bb3.append("#dims\t").append(numAggInputs).append("\t2\t0").nl();
					if(binding!=null){bb3.append(binding);}
					appendReplayPartHeader(bb3,tag,i,replayParts,replayManifestSha80,false,false,replayShuffle,replaySeed); absw[i].print(bb3); bb3.clear();
				}
			}
			// Gene-label replay sidecar (opt-in, subnetlabels=gene only): manifest-identity/order
			// companion to THIS subnet file, one row per emitted subnet row -- see "the change" §4 in
			// subnet_gene_replay_design_20260908.md. Null (zero I/O, zero behavior change) under legacy
			// mode or when no subnet file is being written at all.
			if(sbsw==null || subnetLabels!=SUBNETLABELS_GENE){gsw=null;}
			else{
				gsw=new ByteStreamWriter[replayParts];
				for(int i=0; i<replayParts; i++){
					gsw[i]=openSubnetGeneReplaySidecar(replayPartPath(subnetFile+".bins.tsv",i,replayParts),tag,i,replayParts);
				}
			}
			if(binding!=null){
				isw=openReplayWriters(replayIdentityPath(aggFile),replayParts);
				for(int i=0; i<replayParts; i++){
					final ByteBuilder idsHeader=new ByteBuilder("#schema_version\tmagqc_composite_rows_v1\n")
						.append("#input_width\t").append(numAggInputs).nl().append("#target_width\t2\n")
						.append(binding);
					appendReplayPartHeader(idsHeader,tag,i,replayParts,replayManifestSha80,false,false,replayShuffle,replaySeed);
					idsHeader.append("#columns\trow_index\tbin_id\n"); isw[i].print(idsHeader);
				}
			}
		}catch(RuntimeException|Error setupFailure){
			poisonReplayWriters(bsw); poisonReplayWriters(sbsw); poisonReplayWriters(absw);
			poisonReplayWriters(gsw); poisonReplayWriters(isw);
			throw setupFailure;
		}
		final StringBuilder ioErrors=new StringBuilder();//names of every writer whose poisonAndWait reported an error
		long made=0, identityRows=0;
		try{
			final int count=last-first;
			final int threadCount=Math.min(replayThreads,count);
			final ReplayCoordinator coordinator=new ReplayCoordinator(bsw,sbsw,absw,isw,gsw,replayParts,tag,replayDiagnostics);
			final ReplayProcessThread[] workers=new ReplayProcessThread[threadCount];
			final boolean[] started=new boolean[threadCount];
			Throwable lifecycleFailure=null;
			for(int i=0; i<threadCount && lifecycleFailure==null; i++){
				try{
					workers[i]=new ReplayProcessThread(store,rows,tag,first,last,i,threadCount,coordinator);
				}
				catch(Throwable failure){lifecycleFailure=failure;}
			}
			if(lifecycleFailure==null){
				for(int i=0; i<threadCount; i++){
					try{workers[i].start(); started[i]=true;}
					catch(Throwable failure){lifecycleFailure=failure; coordinator.fail(failure); break;}
				}
			}
			boolean interrupted=false;
			for(int i=0; i<threadCount; i++){
				if(started[i]){
					boolean joined=false;
					while(!joined){
						try{workers[i].join(); joined=true;}
						catch(InterruptedException e){interrupted=true;}
					}
				}
			}
			if(interrupted){Thread.currentThread().interrupt();}
			if(lifecycleFailure!=null){coordinator.fail(lifecycleFailure);}
			final Throwable workerFailure=coordinator.failure();
			if(workerFailure instanceof RuntimeException){throw (RuntimeException)workerFailure;}
			if(workerFailure instanceof Error){throw (Error)workerFailure;}
			if(workerFailure!=null){throw new RuntimeException("Parallel replay failed ("+tag+")",workerFailure);}
			made=coordinator.published();
			if(replayDiagnostics){
				final long denom=Math.max(1,made);
				final long formatExclusive=coordinator.formatAggRowExclusiveNanos();
				final long formatPlusWriter=coordinator.formatAggRowNanos()+coordinator.writerEnqueueNanos();
				final long formatExclusivePlusWriter=formatExclusive+coordinator.writerEnqueueNanos();
				System.err.println(tag+": replay diagnostics rows="+made
					+" subnet_feedForward_calls="+coordinator.feedForwardCalls()
					+" subnet_feedForward_calls_per_emitted_row="+(coordinator.feedForwardCalls()/(double)denom));
				System.err.println(tag+": replay diagnostics nanos_total"
					+" source_row_get="+coordinator.sourceRowGetNanos()
					+" replay_reconstruction="+coordinator.replayReconstructionNanos()
					+" context="+coordinator.contextNanos()
					+" expected_copy="+coordinator.expectedCopyNanos()
					+" subnet_input_assembly="+coordinator.subnetInputNanos()
					+" apply_input_subset="+coordinator.applyInputNanos()
					+" subnet_feedforward="+coordinator.subnetInferenceNanos()
					+" formatAggRow_inclusive="+coordinator.formatAggRowNanos()
					+" formatAggRow_exclusive="+formatExclusive
					+" writer_enqueue="+coordinator.writerEnqueueNanos()
					+" formatAggRow_plus_writer_inclusive="+formatPlusWriter
					+" formatAggRow_plus_writer_exclusive="+formatExclusivePlusWriter);
				System.err.println(tag+": replay diagnostics nanos_per_emitted_row"
					+" source_row_get="+(coordinator.sourceRowGetNanos()/(double)denom)
					+" replay_reconstruction="+(coordinator.replayReconstructionNanos()/(double)denom)
					+" context="+(coordinator.contextNanos()/(double)denom)
					+" expected_copy="+(coordinator.expectedCopyNanos()/(double)denom)
					+" subnet_input_assembly="+(coordinator.subnetInputNanos()/(double)denom)
					+" apply_input_subset="+(coordinator.applyInputNanos()/(double)denom)
					+" subnet_feedforward="+(coordinator.subnetInferenceNanos()/(double)denom)
					+" formatAggRow_inclusive="+(coordinator.formatAggRowNanos()/(double)denom)
					+" formatAggRow_exclusive="+(formatExclusive/(double)denom)
					+" writer_enqueue="+(coordinator.writerEnqueueNanos()/(double)denom)
					+" formatAggRow_plus_writer_inclusive="+(formatPlusWriter/(double)denom)
					+" formatAggRow_plus_writer_exclusive="+(formatExclusivePlusWriter/(double)denom));
			}
			anticodonRows+=coordinator.anticodonRows(); anticodonMissingRows+=coordinator.anticodonMissingRows();
			if(made!=count){throw new IllegalStateException("Replay row conservation failed ("+tag+"): published="+made+" expected="+count);}
			identityRows=(isw==null ? 0 : made);
		}finally{
			// Close EVERY writer and accumulate their error flags (Yoimiya review 2026-09-08: the dense/subnet/
			// aggregator returns were ignored, so gene replay could report success with incomplete vectors).
			// Never throw inside finally: it would mask an in-flight exception; reject after normal completion.
			// ByteStreamWriter KillSwitch-terminates the JVM on a write IOException, so these flags carry the
			// close/flush-time failures -- which the pinned BBTools DROPS (ReadWrite.finishWriting discards
			// close()'s return; records/BBTOOLS_BUGS_FOUND_v1.md row 13). The fixture wrapper proves this path
			// works against an isolated one-line-patched ReadWrite; production correctness needs Brian's upstream fix.
			if(poisonReplayWriters(bsw)){ioErrors.append(" out=").append(file);}
			if(poisonReplayWriters(sbsw)){ioErrors.append(" subnetout=").append(subnetFile);}
			if(poisonReplayWriters(absw)){ioErrors.append(" aggout=").append(aggFile);}
			if(poisonReplayWriters(gsw)){ioErrors.append(" sidecar=").append(subnetFile).append(".bins.tsv");}
			if(poisonReplayWriters(isw)){ioErrors.append(" identities=").append(replayIdentityPath(aggFile));}
		}
		if(ioErrors.length()>0){throw new RuntimeException("I/O error writing replay vectors ("+tag+"):"+ioErrors);}
		if(isw!=null && identityRows!=made){throw new IllegalStateException("Aggregator replay/identity row conservation failed");}
		System.err.println(tag+": replayed "+made+" rows from "+binManifestInFile);
	}

	/**
	 * Builds shared vector/identity provenance before any replay writer opens.
	 * The manifest is hashed once; both opt-in headers receive identical binding bytes.
	 * Row indices in its body are output-local even when named source ranges start later.
	 */
	private ByteBuilder compositeReplayBinding(String aggFile, String denseFile, String subnetFile,
			String split, int first, int last){
		if(FileFormat.isStdout(aggFile)){throw new IllegalArgumentException("aggids=t requires a file aggregator output");}
		final String identityFile=replayIdentityPath(aggFile);
		final HashSet<java.nio.file.Path> outputs=new HashSet<java.nio.file.Path>();
		for(String name:new String[]{aggFile, identityFile, emitDense ? denseFile : null, subnetFile,
			subnetFile!=null && subnetLabels==SUBNETLABELS_GENE ? subnetFile+".bins.tsv" : null}){
			if(name!=null){
				for(int part=0; part<replayParts; part++){
					final java.nio.file.Path path=java.nio.file.Paths.get(replayPartPath(name,part,replayParts)).toAbsolutePath().normalize();
					if(!outputs.add(path)){throw new IllegalArgumentException("Replay vector/identity output paths overlap");}
					if(java.nio.file.Files.exists(path,java.nio.file.LinkOption.NOFOLLOW_LINKS)){
						throw new IllegalArgumentException("Replay vector/identity output must be fresh: "+path);
					}
				}
			}
		}
		assert(first>=0 && last>=first) : "Replay range validation must precede identity ledger construction";
		final String fullHash=replayManifestSha256();
		return new ByteBuilder("#split\t").append(split).nl()
			.append("#bin_model\t").append(binModel==null ? "-" : binModel).nl()
			.append("#bin_manifest_sha80\t").append(fullHash.substring(44)).nl()
			.append("#source_row_start\t").append(first).nl().append("#source_row_end\t").append(last).nl();
	}

	/** Counts the manifest's selected primary-cache shreds. Whole-genome native selections are
	 * deliberately excluded when {@code includeNative=false}; that hybrid mode takes its exact
	 * native count from the contribution index and uses this table only for foreign shreds. */
	private ReferenceCdsSurvivalLabelReader.Counts shredTableCounts(
			MagQCBinManifest.Bin b, boolean includeNative){
		if(refCdsTable==null){throw new RuntimeException("Reference-CDS shred table is not loaded for "+b.id);}
		refCdsNativeIdx.clear();
		if(includeNative){
			for(MagQCBinManifest.Selection s : MagQCBinManifest.decodeSelections(b,false)){
				refCdsNativeIdx.add(refCdsTable.requireShred(s.contigId));
			}
		}
		refCdsForeignIdx.clear();
		for(MagQCBinManifest.Selection s : MagQCBinManifest.decodeSelections(b,true)){
			refCdsForeignIdx.add(refCdsTable.requireShred(s.contigId));
		}
		return refCdsTable.counts(b.targetTid,refCdsNativeIdx,refCdsForeignIdx);
	}

	/** One exact gene-label dispatch, shared by saved-manifest replay and D39 generation. */
	private ReferenceCdsSurvivalLabelReader.Counts exactReplayCounts(MagQCBinManifest.Bin bin){
		if(referenceCdsIndex!=null){return referenceCdsIndex.labelsFor(bin);}
		if(lastContributionCounts!=null){return lastContributionCounts;}
		if(refCdsTable!=null){return shredTableCounts(bin,true);}
		if(referenceCdsContributionIndex!=null){throw new RuntimeException("Contribution-index replay produced no counts for "+bin.id);}
		return null;
	}

	/** Opens the gene-label replay sidecar for one subnet output file (subnetlabels=gene only):
	 *  header identifies the exact bin manifest + reference-CDS index by path and sha256, so a
	 *  reader can re-verify provenance without trusting the file's own claims. See "the change"
	 *  §4 in subnet_gene_replay_design_20260908.md for the schema. */
	private ByteStreamWriter openSubnetGeneReplaySidecar(String subnetFile){
		return openSubnetGeneReplaySidecar(subnetFile+".bins.tsv",null,-1,1);
	}

	private ByteStreamWriter openSubnetGeneReplaySidecar(String sidecarFile,String split,int part,int parts){
		// Header built BEFORE the writer opens (Yoimiya 2026-09-09 review): a hash-source failure
		// here must never leave an opened, unpoisoned writer thread behind. Branches on whichever
		// label source is actually active -- the constructor guarantees exactly one of
		// referenceCdsIndexFile/refCdsTableFile is set whenever subnetlabels=gene reaches this path.
		final ByteBuilder header=new ByteBuilder(320);
		header.append("#schema_version\tsubnet_gene_replay_bins_v1\n");
		if(binModel!=null){header.append("#bin_model\t").append(binModel).nl();}
		header.append("#bin_manifest\t").append(binManifestInFile).append("\t#\tsha256\t")
			.append(replayManifestSha256()).append('\n');
		if(referenceCdsIndexFile!=null){
			header.append("#reference_cds_index\t").append(referenceCdsIndexFile).append("\t#\tsha256\t")
				.append(MagQCBinManifest.sha256File(referenceCdsIndexFile)).append('\n');
		}else{
			if(referenceCdsContributionIndexFile!=null){
				header.append("#reference_cds_contribution_index\t").append(referenceCdsContributionIndexFile)
					.append("\t#\tsha256\t").append(referenceCdsContributionIndexSha256).append('\n');
			}
			if(refCdsTableFile!=null){
				header.append("#reference_cds_table\t").append(refCdsTableFile).append("\t#\tsha256\t")
					.append(ReferenceCdsSurvivalLabelReader.sha256File(refCdsTableFile)).append('\n');
			}
			assert(referenceCdsContributionIndexFile!=null || refCdsTableFile!=null) :
				"subnetlabels=gene requires an exact reference-CDS source (constructor validation)";
		}
		header.append("#subnet\t").append(subnetName).append("\tsubnet_inputs\t").append(subnetInputs).append("\ttargets\t2\n");
		if(parts>1){
			final String manifestSha80=replayManifestSha256().substring(44);
			appendReplayPartHeader(header,split,part,parts,manifestSha80,true,true,replayShuffle,replaySeed);
		}
		header.append("#columns\trow_index\tbin_id\tsplit\tmodel_rows\ttarget_tid\tnative_total\tnative_retained\t")
			.append("foreign_retained\tbin_total_genes\tcompleteness\tcontamination\n");
		final FileFormat ff=FileFormat.testOutput(sidecarFile, FileFormat.TXT, null, false, true, false, false);
		final ByteStreamWriter w=new ByteStreamWriter(ff); w.start();
		w.print(header);
		return w;
	}

	/** One row of the gene-label replay sidecar: every field copied verbatim from the manifest
	 *  Bin the subnet row was just emitted from, or from the SAME reference-CDS Counts that filled
	 *  labels[] for that row -- no invented identities (design §4). */
	private void writeSubnetGeneReplaySidecarRow(ByteStreamWriter w, ByteBuilder row, int rowIndex, MagQCBinManifest.Bin b,
			ReferenceCdsSurvivalLabelReader.Counts reference){
		row.clear();
		appendSubnetGeneReplaySidecarRow(row,rowIndex,b,reference);
		w.print(row);
	}

	/** Formats one sidecar row without binding it to a writer thread. */
	private void appendSubnetGeneReplaySidecarRow(ByteBuilder row, int rowIndex, MagQCBinManifest.Bin b,
			ReferenceCdsSurvivalLabelReader.Counts reference){
		row.clear();
		row.append(rowIndex).tab().append(b.id).tab().append(b.split).tab().append(b.modelRows).tab()
			.append(b.targetTid).tab().append(reference.nativeTotal).tab().append(reference.nativeRetained).tab()
			.append(reference.foreignRetained).tab().append(reference.binTotalGenes).tab();
		appendFmt(row, reference.completeness); row.tab(); appendFmt(row, reference.contamination); row.nl();
	}

	private int replayBin(MagQCBinManifest.Bin b,int[] fam,double[] glob,double[] labels,Agg agg){
		lastContributionCounts=null;
		final boolean useUnshredded="perfect_unshredded".equals(b.spikeClass);
		final long gsize=genomeSize.get(b.targetTid);
		if(gsize<=0){throw new RuntimeException("Manifest target tid has no positive genome size: "+b.targetTid);}
		final ArrayList<Contig> nativeList=new ArrayList<Contig>();
		long cleanBp=0;
		// Bin.schemaVersion (not a separately threaded version) gates decoding -- correct
		// regardless of which Manifest, if any, is in scope at this call site.
		for(MagQCBinManifest.Selection s : MagQCBinManifest.decodeSelections(b,false)){
			if(s.tid!=b.targetTid){throw new RuntimeException("Native selection tid mismatch in "+b.id);}
			final Contig c=findReplayContig(s,useUnshredded); if(c==null){throw new RuntimeException("Native contig not found in "+b.id+": "+s.contigId);}
			nativeList.add(c); agg.add(c); cleanBp+=c.length;
		}
		if(cleanBp<=0){return -1;}
		final boolean contributionRow=(referenceCdsContributionIndex!=null && useUnshredded);
		if(contributionRow){
			lastContributionCounts=applyBreakpointMutation(b,nativeList,fam,agg,useUnshredded);
		}
		lastTarget=b.targetTid; lastUsedUnshredded=useUnshredded; lastMixtureComp=b.compFraction; lastMixtureCont=b.contFraction;
		snapshotTargetSpecial(agg);
		if(subsetRanks!=null){for(int i=0;i<subsetRanks.length;i++){lastFamObs[i]=fam[subsetRanks[i]];}}
		if(cleanFamBuf!=null){System.arraycopy(fam,0,cleanFamBuf,0,numFam);}
		long foreignBp=0;
		for(MagQCBinManifest.Selection s : MagQCBinManifest.decodeSelections(b,true)){
			final Contig c=findReplayContig(s,false); if(c==null){throw new RuntimeException("Foreign contig not found in "+b.id+": "+s.contigId);}
			agg.add(c); foreignBp+=c.length;
		}
		if(contributionRow && foreignBp>0){
			if(!"-".equals(b.breakpoints)){
				throw new RuntimeException("Breakpoint clone may not contain foreign contigs: "+b.id);
			}
			if(refCdsTable==null){
				throw new RuntimeException("Unshredded isolate contamination requires refcdstable= for foreign genes: "+b.id);
			}
			final ReferenceCdsSurvivalLabelReader.Counts foreign=shredTableCounts(b,false);
			if(foreign.nativeTotal!=lastContributionCounts.nativeTotal){
				throw new RuntimeException("Reference-CDS table/index native-total mismatch for tid "+b.targetTid+
					": table="+foreign.nativeTotal+" index="+lastContributionCounts.nativeTotal);
			}
			lastContributionCounts=new ReferenceCdsSurvivalLabelReader.Counts(
				lastContributionCounts.nativeTotal,lastContributionCounts.nativeRetained,foreign.foreignRetained);
		}
		benchNative.clear(); benchForeign.clear(); benchNative.addAll(nativeList);
		benchForeign.addAll(resolveReplayList(b));
		snapshotServeSpecial(agg);
		final long totalBp=cleanBp+foreignBp; if(totalBp<=0){return -1;} lastTotalBp=totalBp;
		final double sampledFraction=Math.min(1.0,cleanBp/(double)gsize);
		final double meanContigLen=agg.contigsGE500>0 ? agg.bpGE500/(double)agg.contigsGE500 : 0;
		labels[0]=sampledFraction*fragmentationFactor(meanContigLen); labels[1]=foreignBp/(double)totalBp;
		// Every exact reference-CDS source replaces these provisional base labels in
		// writeReplaySet. A long foreign intergenic span must not reject a valid
		// gene-labeled bin before that replacement. Keep the cap for legacy replay.
		if(referenceCdsIndex==null && refCdsTable==null && referenceCdsContributionIndex==null &&
				labels[1]>CONT_MAX){return -1;}
		fillFeatureGlobals(totalBp, agg, glob);
		glob[1]=log2(agg.contigs); glob[2]=log2(agg.l50());
		glob[9]=agg.r16; glob[10]=agg.r23; glob[11]=agg.r5; glob[12]=agg.trna;
		final Integer pi=tid2phylumIdx0(b.targetTid); return pi==null ? phylumIndex.get("other") : pi;
	}

	/**
	 * Applies D39's exact whole-genome breakpoint mutation. A breakpoint token is
	 * the canonical {@code g<zero-based-gene-ordinal>} form; duplicate or missing
	 * genes fail rather than being silently counted twice. The first breakpoint
	 * regime is deliberately narrow: an intact unshredded target assembly and no
	 * foreign contigs. DNA, GC, ACGT, dimer, contig, and RNA fields are untouched.
	 */
	private ReferenceCdsSurvivalLabelReader.Counts applyBreakpointMutation(
			MagQCBinManifest.Bin b, ArrayList<Contig> nativeList, int[] fam,
			Agg agg, boolean useUnshredded){
		if(!useUnshredded){
			throw new RuntimeException("Contribution-index breakpoint replay requires spike_class=perfect_unshredded: "+b.id);
		}
		final boolean hasForeign=!"-".equals(b.foreignContigs) || !"-".equals(b.contaminantTids);
		if(hasForeign && !"-".equals(b.breakpoints)){
			throw new RuntimeException("Breakpoint clone may not contain foreign contigs: "+b.id);
		}
		if(hasForeign && refCdsTable==null){
			throw new RuntimeException("Unshredded isolate contamination requires refcdstable= for foreign genes: "+b.id);
		}
		final ArrayList<Contig> complete=getContigsUnshredded(b.targetTid);
		if(complete==null || complete.isEmpty() || nativeList.size()!=complete.size()){
			throw new RuntimeException("Breakpoint replay requires the complete unshredded target assembly for tid "+
				b.targetTid+": selected="+nativeList.size()+" available="+(complete==null ? 0 : complete.size()));
		}
		final HashSet<Contig> selected=new HashSet<Contig>();
		for(Contig c : nativeList){
			if(!selected.add(c)){throw new RuntimeException("Breakpoint replay selected one native contig twice in "+b.id);}
		}
		for(Contig c : complete){
			if(!selected.contains(c)){throw new RuntimeException("Breakpoint replay omitted native contig "+c.name+" in "+b.id);}
		}

		final int nativeTotal=referenceCdsContributionIndex.nativeTotal(b.targetTid);
		if(agg.cds!=nativeTotal){
			throw new RuntimeException("Whole-genome cache/index CDS mismatch for tid "+b.targetTid+
				": cache="+agg.cds+" index="+nativeTotal);
		}
		int familyTotal=0;
		for(int count : fam){familyTotal+=count;}
		if(familyTotal!=agg.mapped){
			throw new RuntimeException("Whole-genome cache mapped/family mismatch for tid "+b.targetTid+
				": mapped="+agg.mapped+" family_copies="+familyTotal);
		}

		final IntHashSet seen=new IntHashSet(16);
		int broken=0;
		if(!"-".equals(b.breakpoints)){
			for(String token : b.breakpoints.split(";",-1)){
				final int ordinal=parseBreakpointOrdinal(token,b.id);
				if(!seen.add(ordinal)){
					throw new RuntimeException("Duplicate breakpoint token "+token+" in "+b.id);
				}
				final int length=referenceCdsContributionIndex.length(b.targetTid,ordinal);
				final int family=referenceCdsContributionIndex.familyRank(b.targetTid,ordinal);
				final long square=(long)length*length;
				if(agg.cds<1 || agg.glenSum<length || agg.glenSq<square || agg.coding<length){
					throw new RuntimeException("Breakpoint mutation would underflow gene statistics for "+token+" in "+b.id);
				}
				agg.cds--; agg.glenSum-=length; agg.glenSq-=square; agg.coding-=length;
				if(family>=0){
					if(family>=fam.length || agg.mapped<1 || fam[family]<1){
						throw new RuntimeException("Breakpoint mutation would underflow assigned family "+family+
							" for "+token+" in "+b.id);
					}
					agg.mapped--; fam[family]--;
				}
				broken++;
			}
		}
		if(agg.cds!=nativeTotal-broken){
			throw new AssertionError("Breakpoint CDS conservation failed for "+b.id);
		}
		return new ReferenceCdsSurvivalLabelReader.Counts(nativeTotal,nativeTotal-broken,0);
	}

	/** Parses one canonical zero-based breakpoint token such as {@code g7}. */
	private static int parseBreakpointOrdinal(String token, String binId){
		if(token==null || token.length()<2 || token.charAt(0)!='g' ||
				(token.length()>2 && token.charAt(1)=='0')){
			throw new RuntimeException("Malformed breakpoint token '"+token+"' in "+binId+
				"; expected canonical g<ordinal>");
		}
		long value=0;
		for(int i=1; i<token.length(); i++){
			final char c=token.charAt(i);
			if(c<'0' || c>'9'){
				throw new RuntimeException("Malformed breakpoint token '"+token+"' in "+binId+
					"; expected canonical g<ordinal>");
			}
			value=value*10+(c-'0');
			if(value>Integer.MAX_VALUE){
				throw new RuntimeException("Breakpoint ordinal overflow in token '"+token+"' of "+binId);
			}
		}
		return (int)value;
	}

	private Contig findReplayContig(MagQCBinManifest.Selection s, boolean unshredded){
		final ArrayList<Contig> list=unshredded ? getContigsUnshredded(s.tid) : getContigs(s.tid); if(list==null){return null;}
		Contig found=null;
		for(Contig c : list){
			if(s.contigId.equals(c.name)){
				if(found!=null){throw new RuntimeException("Ambiguous replay contig tid="+s.tid+" name="+s.contigId);}
				found=c;
			}
		}
		return found;
	}
	private ArrayList<Contig> resolveReplayList(MagQCBinManifest.Bin bin){
		final ArrayList<Contig> out=new ArrayList<Contig>(); for(MagQCBinManifest.Selection s : MagQCBinManifest.decodeSelections(bin,true)){final Contig c=findReplayContig(s,false); if(c==null){throw new RuntimeException("Replay contig not found: "+s.contigId);} out.add(c);} return out;
	}

	/** Opens the shared manifest only for the opt-in implementation path. */
	private void openBinManifest(){
		if(binManifestOutFile==null){return;}
		final String cacheHash=MagQCBinManifest.sha256File(cacheFile);
		final FileFormat ff=FileFormat.testOutput(binManifestOutFile, FileFormat.TXT, null, false, true, false, false);
		binManifestWriter=new ByteStreamWriter(ff); binManifestWriter.start();
		// v2 (Ady/UMP45, 2026-09-09 format repair): every identifier this writer emits is
		// IdentifierCodec-encoded (see contigListText) -- real contig headers can contain commas/
		// semicolons/pipes that v1 could not represent at all (UMP45's real-corpus gate: 547,044
		// comma-bearing ids). The schema-version header is the only encoding gate.
		binManifestWriter.print(new ByteBuilder().append("#schema_version\t").append(MagQCBinManifest.SCHEMA_VERSION_ENCODED).append('\n')
			.append("#source\t").append(cacheFile).append("\t#\tsha256\t").append(cacheHash).append('\n'));
		if(unshreddedCacheFile!=null){
			binManifestWriter.print(new ByteBuilder().append("#source_secondary\t").append(unshreddedCacheFile)
				.append("\t#\tsha256\t").append(MagQCBinManifest.sha256File(unshreddedCacheFile)).append('\n'));
		}
		if(d39Sampling){
			binManifestWriter.print(new ByteBuilder().append("#realization_split\t").append(MagQCBinSplit.VERSION)
				.append("\tmodulus=").append(splitModulus).append("\tvalidation_bucket=0\n"));
		}
		binManifestWriter.print(new ByteBuilder().append("#bin_id\tsplit\ttarget_tid\tcontaminant_tids\tcomp_requested\tcont_requested\tspike_class\tbreakpoints\tnative_contigs\tforeign_contigs\tmodel_rows\n"));
	}

	private void closeBinManifest(){
		if(binManifestWriter!=null){
			if(binManifestWriter.poisonAndWait()){throw new RuntimeException("I/O error writing "+binManifestOutFile);}
			final String h=MagQCBinManifest.sha256File(binManifestOutFile);
			final FileFormat ff=FileFormat.testOutput(binManifestOutFile+".sha256", FileFormat.TXT, null, false, true, false, false);
			final ByteStreamWriter bsw=new ByteStreamWriter(ff); bsw.start(); bsw.print(new ByteBuilder().append(h).append("  ").append(binManifestOutFile).nl());
			if(bsw.poisonAndWait()){throw new RuntimeException("I/O error writing "+binManifestOutFile+".sha256");}
			System.err.println("bin manifest: wrote "+manifestOrdinal+" rows to "+binManifestOutFile+" (sha256="+h+")");
			binManifestWriter=null;
		}
	}

	/** Emits one explicit selection row after makeBin has populated benchNative/benchForeign. */
	private void writeBinManifestRow(String split, double[] labels){
		if(binManifestWriter==null){return;}
		final String id=String.format("bin%012d", manifestOrdinal++);
		final java.util.TreeSet<Integer> foreignTids=new java.util.TreeSet<Integer>();
		for(Contig c : benchForeign){foreignTids.add(c.tid);}
		final long row=("train".equals(split) ? trainOutputOrdinal++ : valOutputOrdinal++);
		final String nativeText=contigListText(benchNative, 'n'), foreignText=contigListText(benchForeign, 'f');
		final StringBuilder tids=new StringBuilder();
		for(int tid : foreignTids){if(tids.length()>0){tids.append(';');} tids.append(tid);}
		// Constructing the immutable row here makes the MVM emission path itself subject to
		// exactly the same role/instance validation as a loaded or written manifest.
		final MagQCBinManifest.Bin validated=new MagQCBinManifest.Bin(id, split, lastTarget,
			tids.length()==0 ? "-" : tids.toString(), lastMixtureComp, lastMixtureCont,
			lastSpikeClass, "-", nativeText, foreignText, split+":"+row, MagQCBinManifest.SCHEMA_VERSION_ENCODED);
		appendManifestBin(validated);
	}

	/** Writes an already validated immutable bin without changing its IDs or memberships. */
	private void appendManifestBin(MagQCBinManifest.Bin validated){
		assert(binManifestWriter!=null) : "Explicit sampled bins require an open provenance manifest";
		final ByteBuilder b=new ByteBuilder(512);
		b.append(validated.id).tab().append(validated.split).tab().append(validated.targetTid).tab()
			.append(validated.contaminantTids).tab().append(Double.toString(validated.compFraction)).tab()
			.append(Double.toString(validated.contFraction)).tab().append(validated.spikeClass).tab()
			.append(validated.breakpoints).tab().append(validated.nativeContigs).tab()
			.append(validated.foreignContigs).tab().append(validated.modelRows).nl();
		binManifestWriter.print(b);
	}

	/** Emits IdentifierCodec-encoded contig_id components (v2 -- see openBinManifest): a real
	 *  contig header can contain commas/semicolons/pipes that would otherwise collide with this
	 *  list's own structural delimiters. The instance-id suffix (role char + zero-padded ordinal)
	 *  is always plain ASCII and never needs escaping, but is encoded too for uniformity -- cheap,
	 *  since IdentifierCodec.encode's fast path returns the same String instance unchanged when
	 *  nothing needs escaping. */
	private static String contigListText(ArrayList<Contig> list, char role){
		if(list==null || list.isEmpty()){return "-";}
		final StringBuilder out=new StringBuilder();
		for(int i=0; i<list.size(); i++){
			if(i>0){out.append(';');}
			final Contig c=list.get(i);
			if(c.name==null){throw new RuntimeException("Manifest selection requires contig names; enable cache name retention");}
			final String ordinal=role+String.format(java.util.Locale.ROOT, "%06d", i);
			out.append(c.tid).append('|').append(IdentifierCodec.encode(c.name)).append('|').append(IdentifierCodec.encode(ordinal));
		}
		return out.toString();
	}

	/** Formats and writes one bin's row to the main output, and (if enabled) the subnet and
	 *  aggregator side outputs. Factored out of writeSet's main loop so the additional-spike
	 *  append loop (same method, same file handles) doesn't duplicate the emission logic. */
	private void emitBin(ByteBuilder bb, ByteBuilder bb2, ByteBuilder bb3, ByteStreamWriter bsw,
			ByteStreamWriter sbsw, ByteStreamWriter absw, int[] fam, double[] glob, int targetPhylumIdx,
			double[] labels, Agg agg){
		if(bsw!=null){
			formatRow(bb, fam, glob, targetPhylumIdx, labels, agg);
			bsw.print(bb); bb.clear();
		}
		if(sbsw!=null){
			if(subnetLabels==SUBNETLABELS_GENE){
				if(subnetFamset){formatFamsetGeneRow(bb2, fam, glob, targetPhylumIdx, labels, agg);}
				else if(subnetAnticodon || subnetRrna){formatSpecialGeneRow(bb2, glob, targetPhylumIdx, labels, agg);}
				else{formatNcrnaGeneRow(bb2, glob, targetPhylumIdx, labels, agg);}
			}else{
				if(subnetFamset){formatFamsetRow(bb2, glob, targetPhylumIdx, agg);}
				else if(subnetAnticodon || subnetRrna){formatSpecialRow(bb2, glob, targetPhylumIdx, agg);}
				else{formatNcrnaRow(bb2, glob, targetPhylumIdx, agg);}
			}
			sbsw.print(bb2); bb2.clear();
		}
		if(absw!=null){
			formatAggRow(bb3, fam, glob, targetPhylumIdx, labels, agg);
			absw.print(bb3); bb3.clear();
		}
	}

	/** Benchmark generation: draws {@code count} synthetic bins from the held-out pool (C in allbutc
	 *  mode) using the SAME makeBin sampler as training, and emits (1) a truth table
	 *  (binID, tid, completeness, contamination, totalBp, nContigs, nForeign), (2) a contig-name
	 *  manifest (binID, contig, native|foreign) for splitting the shred FASTA into per-bin FASTAs,
	 *  and (3) optional aggregator vectors so our net scores the IDENTICAL bins CheckM does. The
	 *  truth labels are the makeBin achieved labels (completeness=cleanBp/gsize, contamination=
	 *  foreignBp/totalBp) - exact ground truth by construction. When extraunshreddedperfectfrac>0,
	 *  ADDITIONALLY appends round(extraUnshreddedPerfectFrac*count) pristine, non-shredded,
	 *  real-genome bins (FORCE_PERFECT_UNSHREDDED) on top of the base draw -- same sizing as
	 *  writeSet's equivalent train/val append, so the benchmark actually contains the
	 *  CheckM2-competitive regime the spike was built to test (UMP45, 2026-08-29). When
	 *  benchvecmono= is also set, additionally emits a monolith-format (numInputs-in, formatRow)
	 *  vector per bin alongside the aggregator vector, so the monolith control (trained on raw
	 *  fam/glob features, no subnet predictions) can be scored on the SAME bins as the composite
	 *  without a dimension mismatch (UMP45, 2026-08-29). */
	private void writeBench(IntList pool, long count, Random rnd){
		assert(benchMode) : "writeBench called outside benchmark mode";
		assert(pool!=null && !pool.isEmpty()) : "benchmark pool empty (need held-out orgs; use poolmode=allbutc)";
		final HashMap<Integer,IntList> famIdx=familyIndex(pool);
		final ByteStreamWriter tsw=new ByteStreamWriter(benchTruthFile, true, false, true); tsw.start();
		final ByteStreamWriter msw=new ByteStreamWriter(benchManifestFile, true, false, true); msw.start();
		final ByteStreamWriter vsw=(benchVecFile==null ? null : new ByteStreamWriter(benchVecFile, true, false, true));
		//Monolith-format vec, alongside the aggregator vec (UMP45 2026-08-29): the composite's
		//benchvec is 9659-in agg vectors, but the monolith control is 8943-in/2-out (raw
		//fam/glob features, no subnet predictions) -- scoring it against the 9659-wide vec would
		//be a dimension mismatch. formatRow is the SAME writer writeSet's primary (non-agg)
		//output already uses, so this is byte-for-byte the monolith's native training format,
		//just for the benchmark's bins instead of a fresh random draw.
		final ByteStreamWriter vmsw=(benchVecMonoFile==null ? null : new ByteStreamWriter(benchVecMonoFile, true, false, true));
		final ByteBuilder tb=new ByteBuilder(256), mb=new ByteBuilder(256);
		final ByteBuilder vb=(vsw==null ? null : new ByteBuilder(numAggInputs*4+64));
		final ByteBuilder vmb=(vmsw==null ? null : new ByteBuilder(numInputs*4+64));
		tb.append("#binID\ttid\tcompleteness\tcontamination\ttotalBp\tnContigs\tnForeign").nl(); tsw.print(tb); tb.clear();
		mb.append("#binID\tcontig\trole").nl(); msw.print(mb); mb.clear();
		if(vsw!=null){vsw.start(); vb.append("#dims\t").append(numAggInputs).append("\t2\t0").nl(); vsw.print(vb); vb.clear();}
		if(vmsw!=null){vmsw.start(); vmb.append("#dims\t").append(numInputs).append("\t2\t0").nl(); vmsw.print(vmb); vmb.clear();}

		final double[] labels=new double[2];
		final int[] fam=new int[numFam];
		final double[] glob=new double[NUM_GLOBALS];
		final Agg agg=new Agg(); agg.setFam(fam);
		long made=0, tries=0;
		while(made<count && tries<count*20+1000){
			tries++;
			Arrays.fill(fam, 0, numFam, 0); agg.reset();
			benchNative.clear(); benchForeign.clear();
			final int targetPhylumIdx=makeBin(pool, famIdx, rnd, fam, glob, labels, agg, FORCE_NONE);
			if(targetPhylumIdx<0){continue;}
			emitBenchBin(tb, mb, vb, vmb, tsw, msw, vsw, vmsw, fam, glob, targetPhylumIdx, labels, agg, made);
			made++;
		}

		// Unshredded-isolate append (UMP45 2026-08-29): writeSet's equivalent block (its
		// unshredded-perfect append) already lets training/validation see pristine,
		// non-shredded, real-genome bins alongside the shredded majority. writeBench had no
		// counterpart, so the benchmark never contained the CheckM2-competitive regime the
		// spike was built to test in the first place -- makeBin was only ever called here
		// with FORCE_NONE. ADDITIONAL to the base `count` draw above, never a replacement;
		// gated the same way as writeSet's version, so extraunshreddedperfectfrac=0 (or no
		// unshreddedcache=) leaves this at zero iterations and bench output unchanged.
		if(extraUnshreddedPerfectFrac>0){
			final long extraUnshreddedPerfect=Math.round(extraUnshreddedPerfectFrac*count);
			long madeExtraUnshredded=0, extraUnshreddedTries=0;
			while(madeExtraUnshredded<extraUnshreddedPerfect
					&& extraUnshreddedTries<extraUnshreddedPerfect*20+1000){
				extraUnshreddedTries++;
				Arrays.fill(fam, 0, numFam, 0); agg.reset();
				benchNative.clear(); benchForeign.clear();
				final int targetPhylumIdx=makeBin(pool, famIdx, rnd, fam, glob, labels, agg,
						FORCE_PERFECT_UNSHREDDED);
				if(targetPhylumIdx<0){continue;}
				emitBenchBin(tb, mb, vb, vmb, tsw, msw, vsw, vmsw, fam, glob, targetPhylumIdx, labels, agg, made);
				made++;
				madeExtraUnshredded++;
			}
			System.err.println("bench: appended "+madeExtraUnshredded+" extra-UNSHREDDED-perfect bins"
				+" (extraunshreddedperfectfrac="+extraUnshreddedPerfectFrac+")");
		}

		tsw.poisonAndWait(); msw.poisonAndWait(); if(vsw!=null){vsw.poisonAndWait();} if(vmsw!=null){vmsw.poisonAndWait();}
		System.err.println("bench: wrote "+made+" bins (tries="+tries+") to "+benchTruthFile+" + "+benchManifestFile
			+(benchVecFile==null ? "" : " + "+benchVecFile)+(benchVecMonoFile==null ? "" : " + "+benchVecMonoFile));
	}

	/**
	 * Deterministic explicit-TID emission for a fixed, frozen panel (PATH_TO_PRODUCTION_v1 A4,
	 * 2026-09-02): for EACH tid listed in {@code panelTidsFile} (one tid per line, manifest
	 * order), emits exactly ONE aggregator row built from ALL of that tid's real contigs in
	 * {@code unshreddedCacheFile} -- no random sampling (makeBin/selectContigs are not called
	 * by this method at all), no injected contamination. labels are hardcoded 1.0 completeness
	 * / 0.0 contamination: this is a definitional fact about how the row was constructed (100%
	 * of the target's own data, 0% foreign), not a claim about the underlying assembly's
	 * biological completeness. Crashes loud if a listed tid has no contigs in
	 * unshreddedCacheFile -- a panel tid absent from its own claimed input is a broken
	 * manifest/cache pairing, not a case to skip silently.
	 */
	private void writePanel(){
		assert(unshreddedCacheFile!=null) : "writePanel requires unshreddedcache=";
		assert(aggManifestFile!=null || aggBundleFile!=null) : "writePanel requires aggmanifest= or bundle= (formatAggRow needs the loaded subnets)";
		final IntList tids=new IntList();
		{
			final ByteFile bf=ByteFile.makeByteFile(panelTidsFile, true);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				tids.add(Parse.parseInt(line, 0, line.length));
			}
			bf.close();
		}
		if(tids.isEmpty()){throw new RuntimeException("paneltids="+panelTidsFile+" contained no tids.");}
		System.err.println("panel: "+tids.size()+" tids from "+panelTidsFile);

		// Validate EVERYTHING before opening the writer (Elly's review, 2026-09-02): a
		// mid-loop RuntimeException after bsw.start() would skip poisonAndWait() and leave the
		// writer thread hanging forever waiting for more input -- exactly the producer-crashes-
		// consumer-hangs class this project has hit before. Fail fast on the whole manifest at
		// once (duplicates + missing contigs), before any writer exists to hang.
		{
			final java.util.HashSet<Integer> seen=new java.util.HashSet<Integer>();
			final ArrayList<Integer> dups=new ArrayList<Integer>();
			final ArrayList<Integer> missing=new ArrayList<Integer>();
			for(int i=0; i<tids.size(); i++){
				final int tid=tids.get(i);
				if(!seen.add(tid)){dups.add(tid);}
				final ArrayList<Contig> contigs=getContigsUnshredded(tid);
				if(contigs==null || contigs.isEmpty()){missing.add(tid);}
			}
			if(!dups.isEmpty() || !missing.isEmpty()){
				throw new RuntimeException("paneltids="+panelTidsFile+" is invalid -- A4 requires exactly "+
					"one row per selected genome. Duplicate tids: "+dups+". Tids with no contigs in "+
					unshreddedCacheFile+": "+missing+".");
			}
		}

		final ByteStreamWriter bsw=new ByteStreamWriter(panelOutFile, true, false, true);
		bsw.start();
		final ByteBuilder bb=new ByteBuilder(numAggInputs*4+64);
		bb.append("#dims\t").append(numAggInputs).append("\t2\t0").nl();
		bsw.print(bb); bb.clear();

		final int[] fam=new int[numFam];
		final double[] glob=new double[NUM_GLOBALS];
		final double[] labels=new double[2];
		final Agg agg=new Agg(); agg.setFam(fam);

		// try/finally (Elly's review, 2026-09-02): the pre-loop validation above rules out bad
		// tids, but it cannot rule out a failure INSIDE formatAggRow/net.feedForward (a malformed
		// net, an assertion, an OOM) -- poisonAndWait() must run on every exit path once the
		// writer is open, or the writer thread hangs forever waiting for more input.
		try{
			for(int i=0; i<tids.size(); i++){
				final int tid=tids.get(i);
				final int phylumIdx=prepareCacheBin(tid, fam, glob, agg, true);

				labels[0]=1.0;//completeness: this row uses 100% of tid's own data, by construction
				labels[1]=0.0;//contamination: zero foreign data injected, by construction

				formatAggRow(bb, fam, glob, phylumIdx, labels, agg);
				bsw.print(bb); bb.clear();
			}
		}finally{
			bsw.poisonAndWait();
		}
		System.err.println("panel: wrote "+tids.size()+" rows to "+panelOutFile);
	}

	/** Formats and writes one bin's truth-table row, contig-manifest rows, and (if enabled)
	 *  aggregator vector row and/or monolith vector row -- factored out of writeBench's main
	 *  loop so the isolate-append loop (same method, same file handles) doesn't duplicate the
	 *  emission logic, matching writeSet/emitBin's existing pattern. Reads benchNative/
	 *  benchForeign/lastTarget/lastTotalBp, which makeBin just populated for the bin at index
	 *  `made`. The monolith vector (vmsw/vmb) is the SAME fam/glob features as the aggregator
	 *  vector, just formatted with formatRow instead of formatAggRow -- one bin, two vector
	 *  representations, so the composite and monolith control can each be scored on the format
	 *  they were actually trained on, against the SAME truth table (UMP45 2026-08-29). */
	private void emitBenchBin(ByteBuilder tb, ByteBuilder mb, ByteBuilder vb, ByteBuilder vmb,
			ByteStreamWriter tsw, ByteStreamWriter msw, ByteStreamWriter vsw, ByteStreamWriter vmsw,
			int[] fam, double[] glob, int targetPhylumIdx, double[] labels, Agg agg, long made){
		final int nContigs=benchNative.size()+benchForeign.size();
		assert(nContigs>0) : "bin passed makeBin with 0 collected contigs (bin"+made+", tid "+lastTarget+")";
		assert(labels[0]>=0 && labels[0]<=1.0001) : "completeness out of range: "+labels[0];
		assert(labels[1]>=0 && labels[1]<=CONT_MAX+1e-9) : "contamination "+labels[1]+" > CONT_MAX "+CONT_MAX;
		final String binID="bin"+made;
		tb.append(binID).tab().append(lastTarget).tab();
		appendFmt(tb, labels[0]); tb.tab(); appendFmt(tb, labels[1]); tb.tab();
		tb.append(lastTotalBp).tab().append(nContigs).tab().append(benchForeign.size()).nl();
		tsw.print(tb); tb.clear();
		for(final Contig c : benchNative){mb.append(binID).tab().append(c.name).append("\tnative").nl();}
		for(final Contig c : benchForeign){mb.append(binID).tab().append(c.name).append("\tforeign").nl();}
		msw.print(mb); mb.clear();
		if(vsw!=null){formatAggRow(vb, fam, glob, targetPhylumIdx, labels, agg); vsw.print(vb); vb.clear();}
		if(vmsw!=null){formatRow(vmb, fam, glob, targetPhylumIdx, labels, agg); vmsw.print(vmb); vmb.clear();}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Bin sampler         ----------------*/
	/*--------------------------------------------------------------*/

	/** Sentinel forcedType values for {@link #makeBin}: roll normally (existing behavior, every
	 *  pre-existing call site), or force the bin-type decision instead of drawing it. */
	private static final int FORCE_NONE=-1, FORCE_PERFECT=0, FORCE_NEARPERFECT=1, FORCE_PERFECT_UNSHREDDED=2;

	/**
	 * Synthesizes one bin into the provided fam[]/glob[] accumulators and labels[].
	 * @param forcedType FORCE_NONE (roll perfect/near-perfect normally, per perfectFrac/nearPerfectFrac --
	 *        the ONLY behavior of every pre-existing call site, RNG-sequence-identical to before this
	 *        parameter existed) or FORCE_PERFECT/FORCE_NEARPERFECT (skip the roll, used only by the NEW
	 *        additional-spike append loop in {@link #writeSet} so the base draw's RNG sequence is
	 *        completely untouched).
	 * @return the target's phylum index, or -1 if the attempt failed (retry).
	 */
	private int makeBin(IntList pool, HashMap<Integer,IntList> famIdx,
			Random rnd, int[] fam, double[] glob, double[] labels, Agg agg, int forcedType){
		return makeBin(pool,pool,famIdx,rnd,fam,glob,labels,agg,forcedType,false);
	}

	/**
	 * Draws native and foreign sources from independently constrained pools. The foreign
	 * family index must derive from foreignPool, so same-family bias cannot escape a quota.
	 * Foreign-required strata suppress clean draws but retain the existing positive-mixture distribution.
	 */
	private int makeBin(IntList pool,IntList foreignPool,HashMap<Integer,IntList> famIdx,
			Random rnd,int[] fam,double[] glob,double[] labels,Agg agg,int forcedType,boolean requireForeign){
		assert(pool.size()>0 && foreignPool.size()>0) : "Preflight must provide nonempty native and foreign source pools";
		final int target=pool.get(rnd.nextInt(pool.size()));
		final long gsize=genomeSize.get(target);
		// unshredded-spike (Brian 2026-08-26): FORCE_PERFECT_UNSHREDDED draws its clean contigs from
		// the REAL whole-genome cache (getContigsUnshredded/recoverable2) instead of the shredded one
		// -- a genuinely unfragmented "perfect" bin, zero shred-boundary artifacts. Every other
		// forcedType (including plain FORCE_PERFECT, which stays shredded per the additive-not-
		// replacing design) is unaffected; selectContigs/Agg below are identical either way, only
		// the SOURCE list differs.
		final boolean d39Ordinary=(forcedType==FORCE_D39_ORDINARY),d39Hq=(forcedType==FORCE_D39_HQ);
		final boolean d39Isolate=(forcedType==FORCE_D39_ISOLATE);
		final boolean useUnshredded=(forcedType==FORCE_PERFECT_UNSHREDDED || d39Isolate);
		final ArrayList<Contig> targetContigs=useUnshredded ? getContigsUnshredded(target) : getContigs(target);
		final long recov=useUnshredded ? recoverable2.get(target) : recoverable.get(target);
		if(gsize<=0 || recov<=0 || targetContigs==null || targetContigs.isEmpty()){return -1;}

		// Bin-type mixture (isolate/high-quality spike, Brian 2026-08-24): a fraction of bins are drawn
		// PERFECT (comp=1, contam=0) or NEAR-PERFECT (comp 0.90-1.0, contam 0-0.05) so the net gets sharp at
		// the low-contamination/high-completeness end (where CheckM2 beats us; the isolate-QC regime). When
		// both fracs are 0 the roll is skipped, so the RNG sequence - and all existing output - is unchanged.
		// forcedType!=FORCE_NONE (the additional-spike append loop, Brian 2026-08-26) skips the roll
		// entirely and decides the type directly -- draws zero extra Random values either way.
		final boolean perfectBin, nearPerfectBin;
		if(forcedType==FORCE_PERFECT || useUnshredded){perfectBin=true; nearPerfectBin=false;}
		else if(d39Ordinary || d39Hq){perfectBin=false; nearPerfectBin=false;}
		else if(forcedType==FORCE_NEARPERFECT){perfectBin=false; nearPerfectBin=true;}
		else{
			final boolean spike=(perfectFrac>0 || nearPerfectFrac>0);
			final double binTypeRoll=(spike ? rnd.nextDouble() : 1.0);
			perfectBin=(binTypeRoll<perfectFrac);
			nearPerfectBin=(!perfectBin && binTypeRoll<perfectFrac+nearPerfectFrac);
		}
		lastSpikeClass=d39Hq ? "hq" : perfectBin ? (useUnshredded ? "perfect_unshredded" : "perfect")
			: (nearPerfectBin ? "nearperfect" : "ordinary");
		// sampled completeness (flat + sqrt-high mixture), or the spiked high-completeness box
		final double comp=perfectBin ? 1.0 : d39Hq ? 0.90+0.10*Math.sqrt(rnd.nextDouble()) :
			(nearPerfectBin ? 0.90+0.10*rnd.nextDouble() : sampleComp(rnd));
		lastMixtureComp=comp;
		// UMP45's catch (2026-08-26): for the unshredded path, comp is always 1.0 (perfectBin
		// forced), but gsize is the SHREDDED sizemap value -- if shredding dropped any remainder/
		// partial fragment, gsize can be LESS than the real whole-genome content (recov), which
		// would cap targetBp below the true total and let selectContigs randomly omit a whole
		// contig from a "perfect" bin. Use recov (the exact sum of cache2's own contigs) instead,
		// so selectContigs's loop (added<targetBp) is guaranteed to add every contig. The
		// completeness LABEL below still divides by the true gsize, unaffected by this.
		// Extended for explicit parity (2026-08-27, UMP45: no numeric effect -- recov<=gsize always
		// for the shredded corpus, verified on real tids 9/24/25, so this already drew every
		// available shred either way): ANY perfectBin (shredded FORCE_PERFECT included, not just
		// FORCE_PERFECT_UNSHREDDED) targets its own recoverable content, not a possibly-larger
		// gsize -- "draw everything available for a perfect bin" now reads the same for both flavors.
		final long targetBp=perfectBin ? recov : (long)(comp*gsize);
		final long cleanBp=selectContigs(targetContigs, targetBp, rnd, agg, collectSelections()?benchNative:null);
		if(cleanBp<=0){return -1;}
		// Snapshot the TARGET's observed ncRNA (before contaminants are added) and its id, so a
		// subnet emitter can pair this bin's observed ncRNA with the target's native complement.
		// Read-only w.r.t. the global vector (agg is only inspected).
		lastTarget=target;
		lastUsedUnshredded=useUnshredded;
		snapshotTargetSpecial(agg);
		// Same snapshot for a famset subnet: the subset families' observed counts, target-only
		// (fam[] holds only target contributions here; contaminants are added below).
		if(subsetRanks!=null){
			for(int i=0; i<subsetRanks.length; i++){lastFamObs[i]=fam[subsetRanks[i]];}
		}
		// Full target-only fam snapshot for aggregator aggobs=clean (fam[] gains
		// contaminant counts below; this preserves the pre-contaminant state).
		if(cleanFamBuf!=null){System.arraycopy(fam, 0, cleanFamBuf, 0, numFam);}

		// sampled contamination: 0 for perfect bins, U[0,0.05] for near-perfect, else clean-spike/square-low
		final double cont=(d39Isolate || requireForeign) ? sampleCont(rnd) : (perfectBin || d39Hq) ? 0.0 : (nearPerfectBin ? 0.05*rnd.nextDouble()
			: (rnd.nextDouble()<cleanSpike ? 0.0 : sampleCont(rnd)));
		lastMixtureCont=cont;
		long foreignBp=0;
		if(cont>0){
			final long foreignTarget=(long)((cont/(1.0-cont))*cleanBp);
			if(foreignTarget>0){
				int nContam=1;
				if(rnd.nextDouble()<multiContamProb){nContam=2+rnd.nextInt(2);}//2 or 3
				final long per=Math.max(1, foreignTarget/nContam);
				for(int k=0; k<nContam; k++){
					final int ctid=pickContaminant(foreignPool, famIdx, target, rnd);
					if(ctid<0){break;}
					foreignBp+=selectContigs(getContigs(ctid), per, rnd, agg, collectSelections()?benchForeign:null);
				}
			}
		}

		final long totalBp=cleanBp+foreignBp;
		if((d39Isolate || requireForeign) && foreignBp==0){return -1;}
		if(totalBp<=0 || agg.contigs<=0){return -1;}
		lastTotalBp=totalBp;
		// Whole-bin (serve-faithful) ncRNA observed, contaminants included - what a
		// deployed bin actually shows. glob[] lacks rother, so capture all 5 here.
		snapshotServeSpecial(agg);

		// achieved labels from the explicit target
		// Completeness relabel (Brian's DECIDED direction via UMP45, 2026-08-27, calibrated in
		// records/COMPLETENESS_LABEL_FORMULA.md): a fragmentation PENALTY on top of the existing
		// basepair-coverage fraction, not a replacement for it. sampledFraction (the ORIGINAL
		// formula) already carries any real basepair shortfall (e.g. multishred's ~1.2%
		// remainder loss, verified this session on real tids); f() additionally penalizes
		// gene-BREAK loss from fragmentation, using ONLY contigs >=500bp for the mean (junk
		// sub-500bp contigs can't hold a full gene and would drag the mean down without
		// representing real fragmentation loss -- see Agg.bpGE500/contigsGE500).
		final double sampledFraction=Math.min(1.0, cleanBp/(double)gsize);
		final double meanContigLen=agg.contigsGE500>0 ? agg.bpGE500/(double)agg.contigsGE500 : 0;
		labels[0]=sampledFraction*fragmentationFactor(meanContigLen);   // completeness
		labels[1]=foreignBp/(double)totalBp;                 // contamination
		// Whole-contig overshoot can push foreign past clean, flipping the dominant organism
		// (Barbara #4). Such a bin's contamination relative to the chosen target is out of the
		// 0-50% spec and ambiguous vs a bin-grader, so reject and retry.
		if(labels[1]>CONT_MAX){return -1;}
		// refcdstable= (opt-in): overwrite the basepair-derived labels above with all-gene
		// reference-CDS ones, exactly as writeReplaySet does via referenceCdsIndex.labelsFor() for
		// the replay path -- see the constructor validation (subnetlabels=gene is required, and
		// FORCE_PERFECT_UNSHREDDED is rejected there too). benchNative holds only target-tid
		// Contigs and benchForeign only non-target-tid Contigs (selectContigs' own call sites,
		// above); collectSelections()==true is guaranteed here because refcdstable= in generation
		// mode requires binmanifestout=.
		if(refCdsTable!=null && !(d39Sampling && useUnshredded)){
			assert(collectSelections()) : "refcdstable= generation requires binmanifestout= (constructor validation)";
			refCdsNativeIdx.clear();
			for(int i=0; i<benchNative.size(); i++){refCdsNativeIdx.add(benchNative.get(i).tableIdx);}
			refCdsForeignIdx.clear();
			for(int i=0; i<benchForeign.size(); i++){refCdsForeignIdx.add(benchForeign.get(i).tableIdx);}
			final ReferenceCdsSurvivalLabelReader.Counts counts=refCdsTable.counts(target, refCdsNativeIdx, refCdsForeignIdx);
			labels[0]=counts.completeness; labels[1]=counts.contamination;
		}

		// global stats
		glob[0]=log2(totalBp);
		glob[1]=log2(agg.contigs);
		glob[2]=log2(agg.l50());
		glob[3]=agg.acgt>0 ? agg.gc/(double)agg.acgt : 0;
		glob[4]=totalBp>0 ? agg.coding/(double)totalBp : 0;
		glob[5]=agg.cds>0 ? agg.glenSum/(double)agg.cds : 0;
		glob[6]=agg.geneStd();
		glob[7]=agg.cds>0 ? agg.mapped/(double)agg.cds : 0;
		glob[8]=log2(agg.richness());
		glob[9]=agg.r16; glob[10]=agg.r23; glob[11]=agg.r5; glob[12]=agg.trna;

		final Integer pi=tid2phylumIdx0(target);
		return pi==null ? phylumIndex.get("other") : pi;
	}

	/** tid2phylumIdx.get() returns int (-1 sentinel); wraps as Integer only at this one
	 *  low-frequency (once-per-bin) call site to keep makeBin's return contract unchanged.
	 *  In labels-mode, the explicit "no column" case (partial/unknown) is returned as the
	 *  numPhyla sentinel (see loadLabelsTaxonomy -- NEVER -1, which would collide with
	 *  makeBin's own "reject and retry" convention for its return value), never as null --
	 *  null here means "route to the other bucket", which only the legacy taxpgm path (an
	 *  unmapped-but-real phylum value) should ever trigger. A tid genuinely absent from
	 *  labels.tsv is a real data problem and crashes loud rather than silently defaulting,
	 *  matching this whole pipeline's crash-loud discipline. */
	private Integer tid2phylumIdx0(int target){
		if(taxonomyFromLabels){
			final Integer v=tid2phylumIdxLabels.get(target);
			if(v==null){throw new RuntimeException("No QuickClade label for tid "+target+" in "+labelsFile);}
			return v;
		}
		final int v=tid2phylumIdx.get(target);
		return v<0 ? null : Integer.valueOf(v);
	}

	/** Picks a contaminant tid != target, favoring the target's family. */
	private int pickContaminant(IntList pool, HashMap<Integer,IntList> famIdx,
			int target, Random rnd){
		if(!tid2family.isEmpty() && rnd.nextDouble()<sameFamProb){
			final int fam=tid2family.get(target);
			final IntList l=famIdx.get(fam);
			if(l!=null && l.size()>1){
				for(int t=0; t<8; t++){final int c=l.get(rnd.nextInt(l.size())); if(c!=target){return c;}}
			}
		}
		for(int t=0; t<8; t++){final int c=pool.get(rnd.nextInt(pool.size())); if(c!=target){return c;}}
		return -1;
	}

	/** Randomly adds contigs (by shuffled draw) until reaching targetBp; returns bp added.
	 *  When {@code collect!=null} (benchmark mode), the selected Contigs are appended to it so a
	 *  FASTA manifest can be emitted - the sampling itself is unchanged (collect is inspect-only). */
	private long selectContigs(ArrayList<Contig> contigs, long targetBp, Random rnd, Agg agg, ArrayList<Contig> collect){
		final int nc=contigs.size();
		final int[] order=agg.scratch(nc);
		for(int i=0; i<nc; i++){order[i]=i;}
		// partial Fisher-Yates: shuffle enough to draw without replacement
		long added=0;
		for(int i=0; i<nc && added<targetBp; i++){
			final int j=i+rnd.nextInt(nc-i);
			final int tmp=order[i]; order[i]=order[j]; order[j]=tmp;
			final Contig c=contigs.get(order[i]);
			agg.add(c);
			if(collect!=null){collect.add(c);}
			added+=c.length;
		}
		return added;
	}

	private boolean collectSelections(){return benchMode || binManifestOutFile!=null;}

	/*--------------------------------------------------------------*/
	/*----------------          Aggregator          ----------------*/
	/*--------------------------------------------------------------*/

	/** Accumulates per-bin sufficient statistics over selected contigs. Reused (S4) across
	 *  attempts via reset() instead of being constructed fresh per makeBin() call - this is
	 *  what lets scratch()'s reuse-if-large-enough buffer actually amortize as designed. */
	static final class Agg {
		int[] fam;
		final int[] anti=new int[STRUCTURAL_ANTICODON_OBS];
		long gc, acgt, coding, glenSum, glenSq;
		int contigs, cds, mapped, r16, r23, r5, rother, trna;
		// Completeness relabel (Brian's direction via UMP45, 2026-08-27,
		// records/COMPLETENESS_LABEL_FORMULA.md): >=500bp-filtered accumulators feeding the
		// fragmentation term's meanContigLen. Junk contigs under 500bp can't hold a full
		// bacterial gene (~1kb avg) and would drag the mean down without representing real
		// gene-fragmentation loss -- kept SEPARATE from the unfiltered contigs/cleanBp above,
		// which still drive the coverage term (min(1.0,cleanBp/gsize)) unchanged.
		long bpGE500; int contigsGE500;
		int[] lens=new int[64]; int nlens=0;
		final int[] dimer=new int[16];//summed dinucleotide counts -> bin-faithful HH/CAGA (null-dimer contigs skipped)
		int[] scratchArr;
		void setFam(int[] fam_){fam=fam_;}
		void reset(){
			gc=acgt=coding=glenSum=glenSq=0;
			contigs=cds=mapped=r16=r23=r5=rother=trna=0;
			Arrays.fill(anti, 0);
			bpGE500=0; contigsGE500=0;
			nlens=0;
			Arrays.fill(dimer, 0);
		}
		int[] scratch(int need){
			if(scratchArr==null || scratchArr.length<need){scratchArr=new int[need];}
			return scratchArr;
		}
		void add(Contig c){
			contigs++; gc+=c.gc; acgt+=c.acgt; coding+=c.coding; cds+=c.cds; mapped+=c.mapped;
			glenSum+=c.glenSum; glenSq+=c.glenSq; r16+=c.r16; r23+=c.r23; r5+=c.r5; rother+=c.rother; trna+=c.trna;
			for(int i=0; i<c.antiRank.length; i++){
				final int rank=c.antiRank[i], count=c.antiCount[i];
				if(rank<0 || rank>=STRUCTURAL_ANTICODON_OBS || count<0){
					throw new RuntimeException("Invalid structural anticodon pair rank="+rank+" count="+count+" tid="+c.tid);
				}
				anti[rank]+=count;
			}
			if(c.length>=500){bpGE500+=c.length; contigsGE500++;}
			for(int i=0; i<c.famRank.length; i++){fam[c.famRank[i]]+=c.famCount[i];}
			if(c.dimer!=null){for(int i=0; i<16; i++){dimer[i]+=c.dimer[i];}}
			if(nlens>=lens.length){lens=Arrays.copyOf(lens, lens.length*2);}
			lens[nlens++]=c.length;
		}
		double geneStd(){
			if(cds<=0){return 0;}
			double mean=glenSum/(double)cds;
			double var=glenSq/(double)cds-mean*mean;
			return var>0 ? Math.sqrt(var) : 0;
		}
		/** Bin-faithful homopolymer-dimer fraction from the SUMMED dinucleotide counts (additive; never
		 *  average per-contig ratios). Zero when no dimer data is present (KmerTracker guards the denominator). */
		float hh(){return KmerTracker.HH(dimer);}
		/** Bin-faithful CAGA strand-skew metric from the summed dinucleotide counts. */
		float caga(){return KmerTracker.CAGA(dimer);}
		int richness(){int r=0; for(int v : fam){if(v>0){r++;}} return r;}
		long l50(){
			if(nlens<=0){return 1;}
			int[] a=Arrays.copyOf(lens, nlens);
			Arrays.sort(a);
			long total=0; for(int v : a){total+=v;}
			long half=total/2, cum=0;
			for(int i=a.length-1; i>=0; i--){cum+=a[i]; if(cum>=half){return a[i];}}
			return a[0];
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Distributions        ----------------*/
	/*--------------------------------------------------------------*/

	// completeness in [COMP_MIN,1]: flat + sqrt(flat) (extra density high)
	private double sampleComp(Random rnd){
		double u=rnd.nextDouble();
		double base=(rnd.nextDouble()<mixComp) ? Math.sqrt(u) : u;
		return COMP_MIN+(1.0-COMP_MIN)*base;
	}
	// contamination in [0,CONT_MAX]: flat + square(flat) (extra density low)
	private double sampleCont(Random rnd){
		double u=rnd.nextDouble();
		double base=(rnd.nextDouble()<mixCont) ? u*u : u;
		return CONT_MAX*base;
	}

	// Fragmentation-BREAK penalty multiplying the existing basepair-coverage fraction
	// (labels[0]=sampledFraction*fragmentationFactor(meanContigLen)). ★ DECIDED (Brian,
	// 2026-08-27), calibrated by UMP45 on primary bytes -- full derivation in
	// records/COMPLETENESS_LABEL_FORMULA.md, mag-qc repo. Curve: f(L)=max(0,1-K/L), K=345 (bp),
	// meanContigLen computed from >=500bp contigs only (Agg.bpGE500/contigsGE500).
	// Calibration: gene-recovery of 20kb-shredding = 0.9711 aggregate (Sigma-fam ratio, shredded
	// vs unshredded cache, 36,760 orgs); all-shreds basepair coverage C=Sigma(recov)/Sigma(gsize)
	// =0.988 (multishred's ~1.2% remainder loss, already carried by sampledFraction). Pin:
	// f(20kb)=0.9711/0.988=0.9828 -> K=20000*(1-0.9828)=345. Check: sampledFraction(0.988) *
	// f(20kb)(0.9828) = 0.971 exactly, Brian's anchor. A closed (~1Mb+) genome -> f~1.0.
	// Best-effort heuristic, not rigor (Brian): estimating an assembly's gene-completeness from
	// that same assembly is inherently circular; some corpus references are JGI-trimmed
	// (500bp/end, <2kb dropped) so the absolute anchor is noisier than the curve SHAPE.
	private static final double FRAGMENTATION_K=345.0;
	private static double fragmentationFactor(double meanContigLen){
		if(meanContigLen<=0){return 1.0;}//degenerate bin (contigsGE500==0): no penalty, not a divide-by-zero
		return Math.max(0.0, 1.0-FRAGMENTATION_K/meanContigLen);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Formatting          ----------------*/
	/*--------------------------------------------------------------*/

	private void precomputeNStrings(){
		nOverN1=new String[N_STR_MAX+1];
		//The agg dense head uses the two-channel encoding regardless of the global enc mode.
		if(enc==ENC_TWO || aggManifestFile!=null || aggBundleFile!=null){excessArr=new String[N_STR_MAX+1];}
		double logCap=Math.log(1+LOG_CAP)/LOG2;
		for(int i=0; i<=N_STR_MAX; i++){
			nOverN1[i]=fmt(encodeCount(i, logCap));
			if(excessArr!=null){int ex=Math.min(Math.max(i-1, 0), EXC_CAP); excessArr[i]=fmt(ex/(double)EXC_CAP);}
		}
	}
	/** Single-column encoding of a family's summed count under the active enc mode. */
	private double encodeCount(int count, double logCap){
		if(enc==ENC_LOG){double v=(Math.log(1+count)/LOG2)/logCap; return v>1 ? 1 : v;}
		if(enc==ENC_RAW){return Math.min(count, RAW_CAP)/(double)RAW_CAP;}
		//TODO: Probable bug - 1+count overflows at Integer.MAX_VALUE; retained for numerical compatibility in the D188 extraction.
		return count/(double)(1+count);//ENC_RATIO (ENC_TWO handles presence separately)
	}
	private void appendFamStr(ByteBuilder bb, int count){
		if(count<=N_STR_MAX){bb.append(nOverN1[count]);}
		else{appendFmt(bb, encodeCount(count, Math.log(1+LOG_CAP)/LOG2));}
	}
	/** Two-channel family: presence (0/1) then excess-copies min(count-1,cap)/cap. */
	private void appendTwo(ByteBuilder bb, int count){
		bb.append(count>0 ? '1' : '0'); bb.tab();
		if(count<=N_STR_MAX){bb.append(excessArr[count]);}
		else{final int ex=Math.min(count-1, EXC_CAP); appendFmt(bb, ex/(double)EXC_CAP);}
		bb.tab();
	}
	/** Fixed-notation float, no exponent (RegressionTrainer's fast parser rejects 'e').
	 *  String-returning form, KEPT for the aggregator's per-subnet context roundtrip (Float.parseFloat(fmt(v)))
	 *  which genuinely needs a String - low frequency (once per bin, not per family). */
	private static String fmt(double v){
		if(v==(long)v){return Long.toString((long)v);}
		return String.format("%.6f", v);
	}
	/** Zero-allocation equivalent of fmt(), appended directly - MUST replicate fmt()'s
	 *  whole-number shortcut exactly: ByteBuilder's own fast append(double,decimals) always
	 *  prints the requested decimals (String.format-style), but a naive append(v,6) would
	 *  differ from fmt() on whole numbers, which appendSlow's precise-but-slow path does NOT
	 *  special-case either - so the whole-number branch is replicated here explicitly. */
	private static void appendFmt(ByteBuilder bb, double v){
		if(v==(long)v){bb.append((long)v);}
		else{bb.appendSlow(v, 6);}
	}
	private static int parseEnc(String s){
		s=s.toLowerCase();
		if(s.equals("ratio")){return ENC_RATIO;}
		if(s.equals("log") || s.equals("log1p")){return ENC_LOG;}
		if(s.equals("raw") || s.equals("linear")){return ENC_RAW;}
		if(s.equals("two") || s.equals("twochannel")){return ENC_TWO;}
		if(s.equals("norm") || s.equals("avgnorm")){return ENC_NORM;}
		throw new RuntimeException("Unknown enc="+s+" (ratio|raw|log|two|norm)");
	}

	private static boolean parseBool(String s){
		return s==null || s.equals("t") || s.equals("true") || s.equals("1") || s.equals("yes");
	}

	/** Requires one caller-pinned, canonical lowercase SHA-256 value. */
	private static void requireLowerHex64(String value, String flag){
		if(value==null || !value.matches("[0-9a-f]{64}")){
			throw new RuntimeException(flag+" requires 64 lowercase hexadecimal characters");
		}
	}

	/** Maps a domain string to the 8-way one-hot index [bact,arch,fungi,plant,animal,protist,virus,other]. */
	private static int domainIndex(String d){
		if(d==null){return DOMAIN_OTHER;}
		final String s=d.toLowerCase();
		if(s.startsWith("bacteri")){return 0;}
		if(s.startsWith("archae")){return 1;}
		if(s.startsWith("fung")){return 2;}
		if(s.contains("viridiplant") || s.startsWith("plant")){return 3;}
		if(s.startsWith("metazoa") || s.startsWith("animal")){return 4;}
		if(s.startsWith("protist")){return 5;}
		if(s.startsWith("vir")){return 6;}
		return DOMAIN_OTHER;//incl. bare "eukaryota" (subkingdom needs tax lineage)
	}


	/**
	 * Per-family expected copy number WHEN PRESENT, over the reference organisms:
	 * for each family, the mean of the organism-level count across the organisms where
	 * that family appears (count &gt; 0). A per-family baseline so enc=norm can encode a
	 * present-at-typical-copies family as ~1.0 and a duplicated family as &gt;1.0,
	 * calibrated to each family's own natural copy number (Brian 2026-08-06). This is a
	 * reference-DB constant (not a label), computed over ALL usable orgs for stability.
	 */
	private double[] computeAvgCopy(IntList usable){
		final long[] sum=new long[numFam];
		final int[] present=new int[numFam];
		final int[] orgCount=new int[numFam];
		for(int i=0; i<usable.size(); i++){
			Arrays.fill(orgCount, 0);
			for(Contig c : getContigs(usable.get(i))){
				for(int j=0; j<c.famRank.length; j++){orgCount[c.famRank[j]]+=c.famCount[j];}
			}
			for(int f=0; f<numFam; f++){if(orgCount[f]>0){sum[f]+=orgCount[f]; present[f]++;}}
		}
		final double[] avg=new double[numFam];
		for(int f=0; f<numFam; f++){avg[f]=present[f]>0 ? sum[f]/(double)present[f] : 1.0;}
		return avg;
	}
	/** enc=norm: 0 if absent, else count / avgCopyWhenPresent[rank]. */
	private void appendNormStr(ByteBuilder bb, int count, int rank){
		if(count==0){bb.append('0'); return;}
		final double a=avgCopy[rank];
		appendFmt(bb, a>0 ? count/a : count);
	}

	/**
	 * Corpus-scale reference for the avg~0.5 normalized context features (Brian 2026-08-24). Over the
	 * usable orgs, computes the mean of log2(genome bp), log2(1+total CDS), mean gene length, and
	 * gene-length std; each such feature is later emitted as {@code 0.5*value/mean} so the population
	 * centers near 0.5 (bounded, well-scaled net inputs instead of huge raw magnitudes). This is an
	 * ORG-level reference (stable and deterministic); partial bins center a little below 0.5, which is
	 * fine - the goal is bounded scaling, not exact 0.5. The four scales are printed to stderr using
	 * round-trip double strings. With a frozen-input bundle, uses its training scales directly,
	 * independently of the current corpus. Size/#genes are scaled on log2
	 * (wide dynamic range); gene-length mean/std on the raw value.
	 */
	private void computeNormScales(IntList usable){
		if(frozenInputs!=null){
			scaleLogBp=frozenInputs.meanLog2Bp;
			scaleLogCds=frozenInputs.meanLog2Cds;
			scaleGlen=frozenInputs.meanGeneLength;
			scaleGlenStd=frozenInputs.meanGeneLengthStd;
			return;
		}
		double sumLogBp=0, sumLogCds=0, sumGlen=0, sumGlenStd=0;
		int nBp=0, nCds=0, nGlen=0;
		for(int i=0; i<usable.size(); i++){
			final int tid=usable.get(i);
			final long bp=genomeSize.get(tid);
			if(bp>0){sumLogBp+=log2(bp); nBp++;}
			long cds=0, glenSum=0, glenSq=0;
			for(Contig c : getContigs(tid)){cds+=c.cds; glenSum+=c.glenSum; glenSq+=c.glenSq;}
			sumLogCds+=log2(1+cds); nCds++;
			if(cds>0){
				final double mean=glenSum/(double)cds;
				final double var=glenSq/(double)cds-mean*mean;
				sumGlen+=mean; sumGlenStd+=(var>0 ? Math.sqrt(var) : 0); nGlen++;
			}
		}
		scaleLogBp=(nBp>0 && sumLogBp>0 ? sumLogBp/nBp : 1);
		scaleLogCds=(nCds>0 && sumLogCds>0 ? sumLogCds/nCds : 1);
		scaleGlen=(nGlen>0 && sumGlen>0 ? sumGlen/nGlen : 1);
		scaleGlenStd=(nGlen>0 && sumGlenStd>0 ? sumGlenStd/nGlen : 1);
		System.err.println("normScales (avg~0.5 reference): mean_log2bp="+scaleLogBp+" mean_log2cds="+scaleLogCds
			+" mean_glen="+scaleGlen+" mean_glenStd="+scaleGlenStd);
	}

	private void formatRow(ByteBuilder bb, int[] fam, double[] glob, int phylumIdx, double[] labels, Agg agg){
		if(enc==ENC_TWO){
			if(keptRanks!=null){
				for(int i=0; i<keptRanks.length; i++){appendTwo(bb, fam[keptRanks[i]]);}
			}else{
				for(int i=0; i<numFam; i++){appendTwo(bb, fam[i]);}
			}
		}else if(enc==ENC_NORM){
			if(keptRanks!=null){
				for(int i=0; i<keptRanks.length; i++){final int r=keptRanks[i]; appendNormStr(bb, fam[r], r); bb.tab();}
			}else{
				for(int i=0; i<numFam; i++){appendNormStr(bb, fam[i], i); bb.tab();}
			}
		}else{
			if(keptRanks!=null){
				for(int i=0; i<keptRanks.length; i++){appendFamStr(bb, fam[keptRanks[i]]); bb.tab();}
			}else{
				for(int i=0; i<numFam; i++){appendFamStr(bb, fam[i]); bb.tab();}
			}
		}
		computeContext(rowCtx, glob, agg);
		appendContext(bb, rowCtx, lastTarget);
		for(int i=0; i<numPhyla; i++){bb.append(i==phylumIdx ? '1' : '0'); bb.tab();}
		appendFmt(bb, labels[0]); bb.tab();
		appendFmt(bb, labels[1]); bb.nl();
	}

	/**
	 * Fills ctx[0..CTX_N) with the FROZEN shared-context scalar block for the current bin
	 * (magqc_rebuild_20260824.plan "FROZEN VECTOR LAYOUT", Brian-confirmed 2026-08-24):
	 * size/#genes/gene-length mean+std at the persisted avg~0.5 corpus scale (computeNormScales),
	 * GC/coding-fraction/log2(richness) unchanged, HH/CAGA BIN-FAITHFUL from agg's SUMMED
	 * dinucleotide counts (additive discipline - never averaged per-contig; KmerTracker.HH/CAGA
	 * guard the zero-denominator case). Shared by every emit method so every net - subnets,
	 * aggregator, and the Tier-A monolith - sees byte-identical context.
	 */
	private void computeContext(double[] ctx, double[] glob, Agg agg){
		ctx[CTX_SIZE]=0.5*glob[0]/scaleLogBp;
		ctx[CTX_GC]=glob[3];
		ctx[CTX_CODING]=glob[4];
		ctx[CTX_GENES]=0.5*log2(1+agg.cds)/scaleLogCds;
		ctx[CTX_L2RICH]=glob[8];
		ctx[CTX_GLEN]=0.5*glob[5]/scaleGlen;
		ctx[CTX_GLENSTD]=0.5*glob[6]/scaleGlenStd;
		ctx[CTX_HH]=agg.hh();
		ctx[CTX_CAGA]=agg.caga();
	}

	/** Appends the FROZEN shared-context block for tid: the CTX_N scalars (already computed
	 *  into ctx by {@link #computeContext}) followed by the domain one-hot(8) - shared-context
	 *  item 10, categorical so it isn't part of the scalar array. */
	private void appendContext(ByteBuilder bb, double[] ctx, int tid){
		appendContextDomain(bb, ctx, domainIdxOf(tid));
	}

	/** Same context encoding when taxonomy is already resolved from a prepared QuickClade row. */
	private void appendContextDomain(ByteBuilder bb, double[] ctx, int domain){
		for(int i=0; i<CTX_N; i++){appendFmt(bb, ctx[i]); bb.tab();}
		for(int i=0; i<DOMAINS; i++){bb.append(i==domain ? '1' : '0'); bb.tab();}
	}

	/** Exact shared arithmetic for the globals consumed by computeContext and the composite tail. */
	private static void fillFeatureGlobals(long totalBp, Agg agg, double[] glob){
		assert(totalBp>0 && agg.fam!=null) : "Replay and prepared bins require positive length and whole-bin family counts";
		glob[0]=log2(totalBp);
		glob[3]=agg.acgt>0 ? agg.gc/(double)agg.acgt : 0;
		glob[4]=totalBp>0 ? agg.coding/(double)totalBp : 0;
		glob[5]=agg.cds>0 ? agg.glenSum/(double)agg.cds : 0;
		glob[6]=agg.geneStd(); glob[7]=agg.cds>0 ? agg.mapped/(double)agg.cds : 0;
		glob[8]=log2(agg.richness());
	}

	/** Two-channel encoding for an aggregator-direct raw count (ncRNA types): presence(0/1)
	 *  then encodeCount(count) under the active enc mode - distinct from the family two-channel
	 *  {@link #appendTwo} (presence+excess-copies), since these aren't duplication-flag counts. */
	private void appendCountTwo(ByteBuilder bb, int count){
		bb.append(count>0 ? '1' : '0'); bb.tab();
		appendFmt(bb, encodeCount(count, Math.log(1+LOG_CAP)/LOG2)); bb.tab();
	}

	/** Snapshots the target-only special observations before contaminants are added. */
	private void snapshotTargetSpecial(Agg agg){
		lastNcObs[0]=agg.r16; lastNcObs[1]=agg.r23; lastNcObs[2]=agg.r5;
		lastNcObs[3]=agg.rother; lastNcObs[4]=agg.trna;
		copyAnticodonSnapshot(lastAntiObs, agg);
	}

	/** Snapshots whole-bin special observations after contaminants are added. */
	private void snapshotServeSpecial(Agg agg){
		lastNcServe[0]=agg.r16; lastNcServe[1]=agg.r23; lastNcServe[2]=agg.r5;
		lastNcServe[3]=agg.rother; lastNcServe[4]=agg.trna;
		copyAnticodonSnapshot(lastAntiServe, agg);
	}

	/**
	 * Returns the whole-genome sum of structurally classified tRNA anticodons.
	 * The returned value is deliberately the classified portion only; callers pair
	 * it with the native tRNA total to form the explicit unknown residual.  Cache
	 * ranks are structural codes in [0,64), and malformed ranks/counts are fatal
	 * because silently folding them into the unknown column would change the model
	 * semantics.
	 */
	private static int structuralAnticodonTotal(ArrayList<Contig> contigs, int tid){
		if(contigs==null){throw new RuntimeException("Missing contigs for structural anticodon total tid="+tid);}
		int total=0;
		for(Contig c : contigs){
			if(c.antiRank==null || c.antiCount==null || c.antiRank.length!=c.antiCount.length){
				throw new RuntimeException("Malformed structural anticodon arrays tid="+tid);
			}
			for(int i=0; i<c.antiRank.length; i++){
				final int rank=c.antiRank[i], count=c.antiCount[i];
				if(rank<0 || rank>=STRUCTURAL_ANTICODON_OBS || count<0){
					throw new RuntimeException("Invalid structural anticodon pair rank="+rank+" count="+count+" tid="+tid);
				}
				total+=count;
				assert(total>=0) : "Structural anticodon total overflow tid="+tid;
			}
		}
		return total;
	}

	/** Copies classified anticodon counts and appends the nonnegative tRNA residual.
	 * This is used for both target-only and serve-faithful snapshots, so the
	 * destination always has exactly 65 cells and sums to the aggregate tRNA count.
	 */
	private void copyAnticodonSnapshot(int[] destination, Agg agg){
		if(!subnetAnticodon && !aggNeedsAnticodons){return;}
		if(destination==null || destination.length<TRNA_ANTICODON_OBS){
			throw new RuntimeException("Structural anticodon snapshot destination width="
				+(destination==null ? -1 : destination.length)+" expected="+TRNA_ANTICODON_OBS);
		}
		if(agg==null || agg.anti==null || agg.anti.length!=STRUCTURAL_ANTICODON_OBS){
			throw new RuntimeException("Malformed structural anticodon aggregate");
		}
		System.arraycopy(agg.anti, 0, destination, 0, STRUCTURAL_ANTICODON_OBS);
		int classified=0;
		for(int i=0; i<STRUCTURAL_ANTICODON_OBS; i++){
			assert(agg.anti[i]>=0) : "Negative structural anticodon aggregate rank="+i;
			classified+=agg.anti[i];
		}
		final int unknown=agg.trna-classified;
		if(unknown<0){throw new RuntimeException("Structural anticodon counts exceed total tRNA: classified="
			+classified+" trna="+agg.trna+" tid="+lastTarget);}
		destination[STRUCTURAL_ANTICODON_OBS]=unknown;
		assert(classified+unknown==agg.trna) : "Structural anticodon snapshot does not conserve tRNA tid="+lastTarget;
	}

	/** Records how often a generated row has an explicit unknown residual. */
	private void noteAnticodonDiagnostic(int nativeUnknown){
		if(subnetAnticodon){
			anticodonRows++;
			if(nativeUnknown<0){throw new RuntimeException("Native structural anticodon counts exceed total tRNA: unknown="+nativeUnknown
				+" tid="+lastTarget);}
			if(nativeUnknown>0){anticodonMissingRows++;}
			assert(anticodonMissingRows<=anticodonRows) : "anticodon diagnostic counters crossed";
		}
	}

	/** Emits the D175 four-rRNA or 64-structural-plus-unknown-anticodon subnet row. */
	private void formatSpecialRow(ByteBuilder bb, double[] glob, int phylumIdx, Agg agg){
		long observed=0;
		if(subnetAnticodon){
			for(int i=0; i<TRNA_ANTICODON_OBS; i++){bb.append(lastAntiObs[i]); bb.tab(); observed+=lastAntiObs[i];}
		}else if(subnetRrna){
			for(int i=0; i<RRNA_OBS; i++){bb.append(lastNcObs[i]); bb.tab(); observed+=lastNcObs[i];}
		}else{throw new RuntimeException("formatSpecialRow called for subnet="+subnetName);}
		for(int i=0; i<numPhyla; i++){bb.append(i==phylumIdx ? '1' : '0'); bb.tab();}
		computeContext(rowCtx, glob, agg);
		appendContext(bb, rowCtx, lastTarget);
		final int nativeTotal;
		if(subnetAnticodon){
			final int structural=lastUsedUnshredded ? nativeAntiTotal_2.get(lastTarget) : nativeAntiTotal.get(lastTarget);
			nativeTotal=lastUsedUnshredded ? nativeNcTrna_2.get(lastTarget) : nativeNcTrna.get(lastTarget);
			noteAnticodonDiagnostic(nativeTotal-structural);
		}else{
				nativeTotal=(lastUsedUnshredded ? nativeNcR16_2.get(lastTarget)+nativeNcR23_2.get(lastTarget)
				+nativeNcR5_2.get(lastTarget)+nativeNcRother_2.get(lastTarget)
				: nativeNcR16.get(lastTarget)+nativeNcR23.get(lastTarget)
				+nativeNcR5.get(lastTarget)+nativeNcRother.get(lastTarget));
		}
		assert(observed<=nativeTotal) : "special subnet observed "+observed+" > native "+nativeTotal
			+" (subnet "+subnetName+", tid "+lastTarget+")";
		bb.append(nativeTotal).nl();
	}

	/**
	 * Emits one ncRNA-subnet training row for the current bin: the target organism's
	 * observed ncRNA counts (r16,r23,r5,rother,trna) plus the FROZEN shared context (phylum
	 * one-hot + the standard context block) as inputs, and the target's NATIVE ncRNA
	 * complement (summed over its whole genome) as the single regression target. The
	 * subnet thus learns the EXPECTED denominator; completeness = observed/expected is
	 * derived downstream (Barbara's refinement of the subset-relative-label design).
	 */
	private void formatNcrnaRow(ByteBuilder bb, double[] glob, int phylumIdx, Agg agg){
		int obsTotal=0;
		for(int i=0; i<NCRNA_OBS; i++){bb.append(lastNcObs[i]); bb.tab(); obsTotal+=lastNcObs[i];}
		for(int i=0; i<numPhyla; i++){bb.append(i==phylumIdx ? '1' : '0'); bb.tab();}
		computeContext(rowCtx, glob, agg);
		appendContext(bb, rowCtx, lastTarget);
		// A FORCE_PERFECT_UNSHREDDED bin's clean contigs came from cache2, whose native ncRNA
		// complement can legitimately exceed the shredded cache's (shredding loses/truncates
		// features at contig boundaries) -- must compare against the SAME source's native total.
		final int nativeTotal=lastUsedUnshredded
			? nativeNcR16_2.get(lastTarget)+nativeNcR23_2.get(lastTarget)+nativeNcR5_2.get(lastTarget)
				+nativeNcRother_2.get(lastTarget)+nativeNcTrna_2.get(lastTarget)
			: nativeNcR16.get(lastTarget)+nativeNcR23.get(lastTarget)+nativeNcR5.get(lastTarget)
				+nativeNcRother.get(lastTarget)+nativeNcTrna.get(lastTarget);
		// The bin's target contigs are a subset of the genome, so observed<=native always.
		assert(obsTotal<=nativeTotal) : "ncRNA observed "+obsTotal+" > native "+nativeTotal+" (tid "+lastTarget+")";
		bb.append(nativeTotal).nl();
	}

	/**
	 * Emits one famset-subnet training row for the current bin: the subset families'
	 * observed counts (target organism only, snapshotted before contaminants) plus the
	 * SAME FROZEN shared context as the ncRNA subnet, and the target's NATIVE subset
	 * total (the subset families' counts over its whole genome) as the single
	 * regression target. The subset is defined by subsetfile= (one family rank per
	 * line), so per-phylum marker sets and co-occurrence modules train with identical
	 * machinery; evaluate with SubnetRatioScore numobs=(subset size).
	 */
	private void formatFamsetRow(ByteBuilder bb, double[] glob, int phylumIdx, Agg agg){
		int obsTotal=0;
		for(int i=0; i<lastFamObs.length; i++){bb.append(lastFamObs[i]); bb.tab(); obsTotal+=lastFamObs[i];}
		for(int i=0; i<numPhyla; i++){bb.append(i==phylumIdx ? '1' : '0'); bb.tab();}
		computeContext(rowCtx, glob, agg);
		appendContext(bb, rowCtx, lastTarget);
		// see formatNcrnaRow's comment: an unshredded-sourced bin must be checked against the
		// unshredded native subset total, not the shredded one.
		final int nativeTotal=lastUsedUnshredded ? nativeFamTotal_2.get(lastTarget) : nativeFamTotal.get(lastTarget);
		assert(nativeTotal>=0) : "nativeFamTotal missing for tid "+lastTarget+" (unshredded="+lastUsedUnshredded+")";
		// The bin's target contigs are a subset of the genome, so observed<=native always.
		assert(obsTotal<=nativeTotal) : "famset observed "+obsTotal+" > native "+nativeTotal+" (tid "+lastTarget+")";
		bb.append(nativeTotal).nl();
	}

	/**
	 * Gene-label counterpart to {@link #formatFamsetRow} (subnetlabels=gene, opt-in, replay-only).
	 * Two differences from the legacy formatter, both from subnet_gene_replay_design_20260908.md
	 * "the change" §2:
	 * <ul>
	 * <li>Observation block is WHOLE-BIN, serve-faithful, contaminants included -- {@code fam[]} at
	 *     emit time already includes both native and foreign contigs (see {@link Agg#add}), unlike
	 *     the native-only {@code lastFamObs} snapshot the legacy formatter reads.</li>
	 * <li>The row carries the TWO all-gene reference-CDS labels ({@code labels[0]}=completeness,
	 *     {@code labels[1]}=contamination, from {@link ReferenceCdsSurvivalLabelReader.Counts}) instead
	 *     of the single native-total regression target. {@code labels[]} is the SAME array
	 *     {@code writeReplaySet} filled for the dense {@code out=} row, so the two agree byte-for-byte
	 *     by construction -- never re-derived from {@code nativeFamTotal}/{@code lastFamObs} sums.</li>
	 * </ul>
	 * Phylum one-hot + shared context are otherwise IDENTICAL to {@code formatFamsetRow} (same column
	 * order), so {@code subnetInputs} (the obs-block width) is unchanged; only the observation source
	 * and target count/values differ.
	 */
	private void formatFamsetGeneRow(ByteBuilder bb, int[] fam, double[] glob, int phylumIdx, double[] labels, Agg agg){
		for(int i=0; i<subsetRanks.length; i++){bb.append(fam[subsetRanks[i]]); bb.tab();}
		for(int i=0; i<numPhyla; i++){bb.append(i==phylumIdx ? '1' : '0'); bb.tab();}
		computeContext(rowCtx, glob, agg);
		appendContext(bb, rowCtx, lastTarget);
		// Reference-CDS Counts always yields ratios in [0,1] (ReferenceCdsSurvivalLabelReader.Counts:
		// completeness=nativeRetained/nativeTotal, contamination=foreignRetained/binTotalGenes, both
		// non-negative numerator<=denominator by construction) -- a value outside [0,1] means labels[]
		// was populated from something other than that Counts object (design §2's explicit ban on
		// substituting tracked-family counts).
		assert(labels[0]>=0 && labels[0]<=1 && labels[1]>=0 && labels[1]<=1) :
			"subnetlabels=gene requires all-gene reference-CDS labels in [0,1]: completeness="
			+labels[0]+" contamination="+labels[1]+" tid="+lastTarget;
		appendFmt(bb, labels[0]); bb.tab(); appendFmt(bb, labels[1]); bb.nl();
	}

	/** Gene-label counterpart to {@link #formatNcrnaRow}: see {@link #formatFamsetGeneRow}'s javadoc
	 *  for the two differences from the legacy formatter. Observation block is {@code lastNcServe}
	 *  (whole-bin ncRNA, serve-faithful) rather than the native-only {@code lastNcObs}. */
	private void formatNcrnaGeneRow(ByteBuilder bb, double[] glob, int phylumIdx, double[] labels, Agg agg){
		for(int i=0; i<NCRNA_OBS; i++){bb.append(lastNcServe[i]); bb.tab();}
		for(int i=0; i<numPhyla; i++){bb.append(i==phylumIdx ? '1' : '0'); bb.tab();}
		computeContext(rowCtx, glob, agg);
		appendContext(bb, rowCtx, lastTarget);
		assert(labels[0]>=0 && labels[0]<=1 && labels[1]>=0 && labels[1]<=1) :
			"subnetlabels=gene requires all-gene reference-CDS labels in [0,1]: completeness="
			+labels[0]+" contamination="+labels[1]+" tid="+lastTarget;
		appendFmt(bb, labels[0]); bb.tab(); appendFmt(bb, labels[1]); bb.nl();
	}

	/** Gene-label counterpart for the D175 special subnet modes.  The observation block is
	 * whole-bin/serve-faithful, while the two all-gene labels remain the canonical targets. */
	private void formatSpecialGeneRow(ByteBuilder bb, double[] glob, int phylumIdx, double[] labels, Agg agg){
		if(subnetAnticodon){
			final int structural=lastUsedUnshredded ? nativeAntiTotal_2.get(lastTarget) : nativeAntiTotal.get(lastTarget);
			final int nativeTrna=lastUsedUnshredded ? nativeNcTrna_2.get(lastTarget) : nativeNcTrna.get(lastTarget);
			noteAnticodonDiagnostic(nativeTrna-structural);
		}
		if(subnetAnticodon){
			for(int i=0; i<TRNA_ANTICODON_OBS; i++){bb.append(lastAntiServe[i]); bb.tab();}
		}else if(subnetRrna){
			for(int i=0; i<RRNA_OBS; i++){bb.append(lastNcServe[i]); bb.tab();}
		}else{throw new RuntimeException("formatSpecialGeneRow called for subnet="+subnetName);}
		for(int i=0; i<numPhyla; i++){bb.append(i==phylumIdx ? '1' : '0'); bb.tab();}
		computeContext(rowCtx, glob, agg);
		appendContext(bb, rowCtx, lastTarget);
		assert(labels[0]>=0 && labels[0]<=1 && labels[1]>=0 && labels[1]<=1) :
			"subnetlabels=gene requires all-gene reference-CDS labels in [0,1]: completeness="
			+labels[0]+" contamination="+labels[1]+" tid="+lastTarget;
		appendFmt(bb, labels[0]); bb.tab(); appendFmt(bb, labels[1]); bb.nl();
	}

	/**
	 * Emits one aggregator training row for the current bin. For every manifest
	 * subnet: gathers its observed counts (whole-bin by default - serve-faithful,
	 * contaminants included; target-only under aggobs=clean), builds the subnet's
	 * input exactly as its training rows were built (same column order, same fmt()
	 * rounding so in-process values match what a file round-trip would deliver),
	 * feed-forwards the frozen net, and emits [ratio, log-obs, log-pred, zero-flag].
	 * Then the pooled ratio baseline, the raw dense head (enc=two presence+excess of
	 * the top-K prevalent families), phylum one-hot, the FROZEN shared context, the
	 * raw ncRNA two-channel(5x2) direct aggregator input, mapped-fraction, and the
	 * GLOBAL comp/contam targets.
	 */
	private void formatAggRow(ByteBuilder bb, int[] fam, double[] glob, int phylumIdx, double[] labels, Agg agg){
		formatAggRowCore(bb, fam, glob, phylumIdx, labels, agg, true);
	}

	/** Pure extraction of the established aggregator formatter.  The false-label mode emits
	 * only the input feature row for deployment; the training formatter above remains the
	 * sole target-bearing path. */
	private void formatAggRowCore(ByteBuilder bb, int[] fam, double[] glob, int phylumIdx, double[] labels, Agg agg, boolean includeTargets){
		formatAggRowCore(bb, fam, glob, phylumIdx, labels, agg, includeTargets, domainIdxOf(lastTarget), null);
	}

	/** Shared numerical fan-out; explicit taxonomy removes the serving path's need for a reference TID. */
	private void formatAggRowCore(ByteBuilder bb, int[] fam, double[] glob, int phylumIdx, double[] labels,
			Agg agg, boolean includeTargets, int domainIdx, String preparedId){
		final int[] famArr=(aggObsServe ? fam : cleanFamBuf);
		final int[] ncArr=(aggObsServe ? lastNcServe : lastNcObs);
		final ReplayDiagnostics diag=replayDiag;
		// Per-bin context, computed ONCE (shared by every subnet AND the aggregator tail - every
		// net now sees the identical FROZEN context, the old per-net rep flags are retired).
		final long contextStart=(diag==null ? 0 : System.nanoTime());
		computeContext(rowCtx, glob, agg);
		if(quantizedRowCtx==null || quantizedRowCtx.length!=CTX_N){quantizedRowCtx=new float[CTX_N];}
		for(int i=0; i<CTX_N; i++){quantizedRowCtx[i]=Float.parseFloat(fmt(rowCtx[i]));}
		final int[] antiArr=aggObsServe ? lastAntiServe : lastAntiObs;
		prepareSubnetInputs(famArr,ncArr,antiArr,phylumIdx,domainIdx);
		if(diag!=null){diag.contextNanos+=System.nanoTime()-contextStart;}

		final boolean six=(subnetFeatures==SUBNETFEATURES_SIX);
		if(six){
			// Whole-bin observations for the table-derived features, stated once per row (design
			// revision2 §4.1 step "Formatter"): NOT famArr/ncArr -- those alias fam/lastNcServe
			// under six anyway (aggobs=clean is rejected at construction), but reading fam/
			// lastNcServe directly matches the spec's explicit "NOT famArr, NOT ncArr" wording.
			Arrays.fill(observedByItem, 0);
			for(int r=0; r<numFam; r++){observedByItem[r]=fam[r];}
			for(int i=0; i<NCRNA_OBS; i++){observedByItem[numFam+i]=lastNcServe[i];}
			if(observedByItem.length==numFam+NCRNA_OBS+TRNA_ANTICODON_OBS){
				for(int i=0; i<TRNA_ANTICODON_OBS; i++){
					observedByItem[numFam+NCRNA_OBS+i]=lastAntiServe[i];
				}
			}
		}
		long sumObs=0;
		double sumPred=0;
		for(int si=0; si<aggSubnets.size(); si++){
			final long inputStart=(diag==null ? 0 : System.nanoTime());
			final AggSubnet s=aggSubnets.get(si);
			// Input layout mirrors formatFamsetRow/formatNcrnaRow exactly: obs, phylum one-hot,
			// shared context (CTX_N scalars + domain one-hot). Identical under legacy and six --
			// only what happens AFTER feedForward differs.
			int p=0;
			long obsTotal=0;
			if("ncrna".equals(s.type)){
				System.arraycopy(subnetObservedInputs,numFam,s.buf,0,NCRNA_OBS); p=NCRNA_OBS;
				if(!six){for(int i=0; i<NCRNA_OBS; i++){obsTotal+=ncArr[i];}}
			}else if("rrna".equals(s.type)){
				System.arraycopy(subnetObservedInputs,numFam,s.buf,0,RRNA_OBS); p=RRNA_OBS;
				if(!six){for(int i=0; i<RRNA_OBS; i++){obsTotal+=ncArr[i];}}
			}else if("trna_anticodon".equals(s.type)){
				System.arraycopy(subnetObservedInputs,numFam+NCRNA_OBS,s.buf,0,TRNA_ANTICODON_OBS); p=TRNA_ANTICODON_OBS;
				if(!six){for(int i=0; i<TRNA_ANTICODON_OBS; i++){obsTotal+=antiArr[i];}}
			}else if("famset".equals(s.type)){
				for(int r : s.ranks){s.buf[p++]=subnetObservedInputs[r]; if(!six){obsTotal+=famArr[r];}}
			}else{
				throw new RuntimeException("Unknown aggregator subnet type "+s.type+" for "+s.name);
			}
			// The complete shared suffix already has training-file rounding and taxonomy.
			System.arraycopy(subnetSharedInputs,0,s.buf,p,subnetSharedInputs.length);
			p+=subnetSharedInputs.length;
			assert(p==s.buf.length) : s.name+": filled "+p+" of "+s.buf.length;
			final ml.CellNet net;
			if(sharedPreparedSubnets){
				net=null;
				final long inferenceStart=(diag==null ? 0 : System.nanoTime());
				s.bundleSubnet.scorePreparedSynced(s.buf, s.preparedOutput);
				if(diag!=null){diag.feedForwardCalls++; diag.subnetInferenceNanos+=System.nanoTime()-inferenceStart;}
			}else{
				net=(preparedDenseInference ? s.netForPreparedThread() : s.netForCurrentThread());
				final long applyStart=(diag==null ? 0 : System.nanoTime());
				net.applyInput(s.buf);
				if(diag!=null){
					final long applied=System.nanoTime();
					diag.applyInputNanos+=applied-applyStart;
					diag.subnetInputNanos+=applied-inputStart;
				}
				final long inferenceStart=(diag==null ? 0 : System.nanoTime());
				net.feedForward();
				if(diag!=null){diag.feedForwardCalls++;}
				if(diag!=null){diag.subnetInferenceNanos+=System.nanoTime()-inferenceStart;}
			}
			if(six){
				for(int k=0; k<4; k++){
					final float o=sharedPreparedSubnets ? s.preparedOutput[k] : net.getOutput(k);
					if(!Float.isFinite(o)){throw new RuntimeException("subnetfeatures=six: non-finite learned output "
						+sixNames[k]+" from subnet "+s.name+" (order "+si+") for bin "+
						(preparedId==null ? "tid="+lastTarget : preparedId));}
					appendFmt(bb, o); bb.tab();
				}
				final long expectedCopyStart=(diag==null ? 0 : System.nanoTime());
				final MagQCExpectedCopyFeatures.Result r=expectedCopy.compute(si, observedByItem);
				if(diag!=null){diag.expectedCopyNanos+=System.nanoTime()-expectedCopyStart;}
				appendFmt(bb, r.observedExpected); bb.tab();
				appendFmt(bb, r.excessExpected); bb.tab();
			}else{
				final double pred=sharedPreparedSubnets ? s.preparedOutput[0] : net.getOutput(0);
				final double ratio=Math.min(RATIO_CAP, obsTotal/Math.max(0.5, pred));
				appendFmt(bb, ratio); bb.tab();
				appendFmt(bb, log2(1+obsTotal)); bb.tab();
				appendFmt(bb, log2(1+Math.max(0, pred))); bb.tab();
				bb.append(obsTotal==0 ? '1' : '0'); bb.tab();
				sumObs+=obsTotal;
				sumPred+=Math.max(0, pred);
			}
		}
		if(!six){
			// Pooled ratio -- DROPPED under six (design revision2 §4.2: a gene-fraction output has
			// no defined denominator to pool; no replacement pooled statistic is proposed).
			appendFmt(bb, Math.min(RATIO_CAP, sumObs/Math.max(1, sumPred))); bb.tab();
		}
		// Dense head is ALWAYS whole-bin (the aggregator's direct deployment signal, not a
		// subnet input; aggobs= only A/Bs the subnet-obs question - settled with Eru 2026-08-11).
		for(int r : denseRanks){appendTwo(bb, fam[r]);}
		for(int i=0; i<numPhyla; i++){bb.append(i==phylumIdx ? '1' : '0'); bb.tab();}
		appendContextDomain(bb, rowCtx, domainIdx);
		// Aggregator-only direct inputs (FROZEN VECTOR LAYOUT): raw ncRNA two-channel (5 types x
		// [presence, encodeCount]) then mapped-fraction. Whole-bin/target-only follows the same
		// aggobs= switch as the subnet inputs above (ncArr), for serve-faithfulness consistency.
		for(int i=0; i<NCRNA_OBS; i++){appendCountTwo(bb, ncArr[i]);}
		appendFmt(bb, glob[7]);
		if(includeTargets){
			bb.tab(); appendFmt(bb, labels[0]); bb.tab(); appendFmt(bb, labels[1]);
		}
		bb.nl();
	}

	/**
	 * Prepares the observed counts and shared taxonomy/context suffix once per row.
	 * Subnets gather their observation ranks and copy the suffix without formatting
	 * or recomputing features. The arrays belong to this maker/worker and are reused.
	 * Context values were quantized through the unchanged training-file formatter.
	 */
	private void prepareSubnetInputs(final int[] fam,final int[] nc,final int[] anti,
			final int phylumIdx,final int domainIdx){
		assert(fam.length>=numFam && nc.length>=NCRNA_OBS && anti.length>=TRNA_ANTICODON_OBS)
			: "Subnet input layout requires all family, ncRNA and anticodon observations";
		final int observedWidth=numFam+NCRNA_OBS+TRNA_ANTICODON_OBS;
		final int sharedWidth=numPhyla+CTX_N+DOMAINS;
		if(subnetObservedInputs==null || subnetObservedInputs.length!=observedWidth){subnetObservedInputs=new float[observedWidth];}
		if(subnetSharedInputs==null || subnetSharedInputs.length!=sharedWidth){subnetSharedInputs=new float[sharedWidth];}
		for(int i=0; i<numFam; i++){subnetObservedInputs[i]=fam[i];}
		for(int i=0; i<NCRNA_OBS; i++){subnetObservedInputs[numFam+i]=nc[i];}
		for(int i=0; i<TRNA_ANTICODON_OBS; i++){subnetObservedInputs[numFam+NCRNA_OBS+i]=anti[i];}
		Arrays.fill(subnetSharedInputs,0);
		if(phylumIdx>=0 && phylumIdx<numPhyla){subnetSharedInputs[phylumIdx]=1;}
		System.arraycopy(quantizedRowCtx,0,subnetSharedInputs,numPhyla,CTX_N);
		if(domainIdx>=0 && domainIdx<DOMAINS){subnetSharedInputs[numPhyla+CTX_N+domainIdx]=1;}
	}

	/** Builds one bin's Agg and global fields from the already-loaded 19-field cache rows.
	 * This is the same accumulation/state setup used by paneltids=; it exists as a shared seam
	 * so deployment cannot grow a second cache-row interpretation. */
	private int prepareCacheBin(final int tid, final int[] fam, final double[] glob, final Agg agg, final boolean unshredded){
		final ArrayList<Contig> contigs=unshredded ? getContigsUnshredded(tid) : getContigs(tid);
		if(contigs==null || contigs.isEmpty()){throw new RuntimeException("No loaded cache rows for tid "+tid);}
		Arrays.fill(fam, 0, numFam, 0); agg.reset();
		long totalBp=0;
		for(final Contig c : contigs){agg.add(c); totalBp+=c.length;}
		lastTarget=tid; lastUsedUnshredded=unshredded;
		snapshotTargetSpecial(agg);
		snapshotServeSpecial(agg);
		if(cleanFamBuf!=null){System.arraycopy(fam, 0, cleanFamBuf, 0, numFam);}
		fillFeatureGlobals(totalBp, agg, glob);
		glob[1]=log2(agg.contigs); glob[2]=log2(agg.l50());
		glob[9]=agg.r16; glob[10]=agg.r23;
		glob[11]=agg.r5; glob[12]=agg.trna;
		final Integer pi=tid2phylumIdx0(tid);
		return pi==null ? phylumIndex.get("other") : pi;
	}

	/** Returns one input-only aggregator row for a single tid represented in the loaded cache.
	 * Callers must initialize this instance through process() so bundle, vocabulary, scales, and
	 * dense-head state are the exact production state. Targets are deliberately never emitted. */
	float[] formatSingleBinAggRow(final int tid){
		final float[] row=parseSingleBinAggRow(formatSingleBinAggRowBytes(tid));
		if(row.length!=numAggInputs){throw new RuntimeException("single-bin aggregator row width="+row.length+" expected "+numAggInputs);}
		return row;
	}

	/** Returns the input-only row in the exact ByteBuilder representation used by the formatter. */
	ByteBuilder formatSingleBinAggRowBytes(final int tid){
		if(aggSubnets==null || numAggInputs<=0){throw new IllegalStateException("aggregator state is not initialized");}
		final int[] fam=new int[numFam]; final double[] glob=new double[NUM_GLOBALS]; final Agg agg=new Agg(); agg.setFam(fam);
		final int phylumIdx=prepareCacheBin(tid, fam, glob, agg, false);
		final ByteBuilder bb=new ByteBuilder(numAggInputs*8+64);
		formatAggRowCore(bb, fam, glob, phylumIdx, null, agg, false);
		return bb;
	}

	static float[] parseSingleBinAggRow(final ByteBuilder bb){
		final String[] fields=bb.toString().trim().split("\\t",-1);
		final float[] row=new float[fields.length];
		for(int i=0; i<row.length; i++){row[i]=Float.parseFloat(fields[i]);}
		return row;
	}

	/** In labels-mode, returns the labels-derived domain index (real column, or the DOMAINS
	 *  sentinel for the unknown-status explicit no-column case -- see
	 *  {@link #loadLabelsTaxonomy}), throwing on a tid genuinely absent from labels.tsv.
	 *  Otherwise, unchanged: the cache-field-2-derived domain (filename/membership ground
	 *  truth), defaulting an unrecognized tid to DOMAIN_OTHER. */
	private int domainIdxOf(int tid){
		if(taxonomyFromLabels){
			final Integer di=tid2domainIdxLabels.get(tid);
			if(di==null){throw new RuntimeException("No QuickClade label for tid "+tid+" in "+labelsFile);}
			return di;
		}
		final int v=tidToIdx.get(tid);
		return v<0 || v>=domainIdxArr.length ? DOMAIN_OTHER : domainIdxArr[v];
	}

	private static double log2(double v){return v<=1 ? 0 : Math.log(v)/LOG2;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String cacheFile, sizemapFile, familyFile, taxpgmFile, treeFile, out, outval, featuresFile;
	/** Disable unused dense row construction for subnet-only D39 generation and replay. */
	private boolean emitDense=true;
	private String binManifestOutFile;
	private String binManifestInFile;
	private String replayManifestSha256;
	/** Optional named membership list in the shared bin manifest's model_rows column. */
	private String binModel;
	/** Optional scratch directory for named-model replay; default is the JVM temporary directory. */
	private String replayTmpDir;
	private MagQCBinReplayStore replayStore;
	private long replaySeed;
	private boolean replaySeedSet;
	private boolean replayShuffle;
	private boolean replayMemory;
	/** Opt-in replay stage counters/timers; deliberately false for the production path. */
	private boolean replayDiagnostics=false;
	/** Set only on a replay worker when replayDiagnostics is enabled. */
	private ReplayDiagnostics replayDiag;
	private int replayThreads=8;
	/** Unreleased row-blocking candidate; zero retains the frozen row-major path. */
	private int replayBlockRows=0;
	private int replayParts=1;
	private long replayTrainStart=-1,replayTrainEnd=-1,replayValStart=-1,replayValEnd=-1;
	private boolean replayTrainStartSet,replayTrainEndSet,replayValStartSet,replayValEndSet;
	private int replayTrainStartResolved,replayTrainEndResolved,replayValStartResolved,replayValEndResolved;
	/** Opt-in D39 sampler with optional reference-phylum quotas independent of input taxonomy. */
	private boolean d39Sampling=false;
	private String sampleModel;
	private int splitModulus=MagQCBinSplit.DEFAULT_MODULUS;
	private boolean splitModulusSet=false;
	private String samplingTaxonomyFile,samplePhylum;
	/** Explicit D145 opt-in; unknown references remain eligible for unstratified ordinary draws. */
	private boolean samplingUnknownOrdinary=false;
	private IntList samplingInside,samplingOutside;
	private String excludedTidsFile;
	private final IntHashSet excludedTids=new IntHashSet(32);
	private final HashSet<D39Fingerprint> d39Realizations=new HashSet<D39Fingerprint>();
	private static final int FORCE_D39_ORDINARY=4,FORCE_D39_HQ=5,FORCE_D39_ISOLATE=6;
	private String referenceCdsIndexFile;
	private MagQCReferenceCdsIndex referenceCdsIndex;
	private String referenceCdsContributionIndexFile,referenceCdsContributionIndexSha256;
	private String referenceCdsContributionManifestSha256,referenceCdsContributionCacheSha256;
	private ReferenceCdsGeneContributionIndex.Reader referenceCdsContributionIndex;
	private ReferenceCdsSurvivalLabelReader.Counts lastContributionCounts;
	private String refCdsTableFile, refCdsTableSha256;
	private ReferenceCdsShredSurvivalTableReader refCdsTable;
	// Reusable per-bin scratch (refcdstable= only; hoisted so makeBin/writeReplaySet never allocate
	// per-attempt, matching the fam[]/glob[]/labels[] hoisting elsewhere in this class).
	private IntList refCdsNativeIdx=new IntList();
	private IntList refCdsForeignIdx=new IntList();
	private ByteStreamWriter binManifestWriter;
	private long manifestOrdinal=0;
	private long trainOutputOrdinal=0, valOutputOrdinal=0;
	// Step 4 of the rebuild path (UMP45, 2026-09-03): labels=<C1 labels.tsv> replaces taxpgm= as
	// the domain/phylum one-hot source. Exactly one of taxpgmFile/labelsFile is non-null after the
	// constructor's check; taxonomyFromLabels caches that decision for loadAux/tid2phylumIdx0/
	// domainIdxOf. tid2phylumIdxLabels/tid2domainIdxLabels are keyed by FILENAME tid (from
	// labels.tsv's source_rel), storing a real one-hot column index, or -1 for an explicit
	// no-column case (partial's phylum, unknown's domain+phylum) -- see loadLabelsTaxonomy().
	private String labelsFile=null;
	private boolean taxonomyFromLabels=false;
	private final HashMap<Integer,Integer> tid2phylumIdxLabels=new HashMap<Integer,Integer>();
	private final HashMap<Integer,Integer> tid2domainIdxLabels=new HashMap<Integer,Integer>();
	private String subnetName, subnetOut, subnetValOut, subsetFile;
	// subnetlabels= (opt-in, default legacy): see parseSubnetLabels's javadoc.
	private int subnetLabels=SUBNETLABELS_LEGACY;
	// subnetfeatures= (opt-in, default legacy): see parseSubnetFeatures's javadoc. Six-mode-only
	// fields: expectedCopy/observedByItem/sixNames/sixUnits are null/unset under legacy.
	private int subnetFeatures=SUBNETFEATURES_LEGACY;
	private static final int SUBNETFEATURES_LEGACY=0, SUBNETFEATURES_SIX=1;
	private String expectedCopyTableFile, expectedCopyTableSha256, subnetPopulationsFile, subnetPopulationsSha256;
	private MagQCExpectedCopyBinding expectedCopy;
	private int[] observedByItem;
	private String[] sixNames, sixUnits;
	// Aggregator row width knobs (design revision2 §4.1, "root: one shared calculation"):
	// 4/1 under legacy (byte-identical to the pre-six literal formula), 6/0 under six.
	private int subnetBlockWidth=4, pooledCols=1;
	private String aggManifestFile, aggBundleFile, aggOut, aggValOut;
	private boolean aggIds=false;
	private String releaseManifestFile,releaseManifestSha80;
	private MagQCNetBundle loadedBundle;
	private MagQCNetBundle.FrozenInputs frozenInputs;
	// A4 explicit-panel emission (2026-09-02): see writePanel(). Both null unless paneltids= is given.
	private String panelTidsFile=null, panelOutFile=null;
	private int denseHead=100;
	private boolean aggObsServe=true;
	private int poolMode=POOL_TRAINVAL;
	private boolean noOrganismHoldout=false;
	private static final int POOL_TRAINVAL=0, POOL_VALSPLIT=1, POOL_ALLBUTC=2;
	private static final int SUBNETLABELS_LEGACY=0, SUBNETLABELS_GENE=1;
	private ArrayList<AggSubnet> aggSubnets;
	private int[] denseRanks;
	private int[] cleanFamBuf;
	// FROZEN shared-context scratch buffer (CTX_N scalars), reused across every format*Row call
	// for the current bin - computeContext() overwrites it fully each time, no partial-update risk.
	private double[] rowCtx;
	/** Per-row exact formatter quantization, shared by all aggregator subnets in that row. */
	private float[] quantizedRowCtx;
	/** Per-worker observation cache and common suffix, overwritten once per bin. */
	private float[] subnetObservedInputs,subnetSharedInputs;
	private double[] preparedGlob;
	/** Optional prepared-worker mode and clone guard; never enabled for replay generation. */
	private boolean sharedPreparedSubnets, preparedWorker;
	/** True only for initializePrepared(); keeps D196 dense serving scoped to prepared input. */
	private boolean preparedDenseInference=false;
	/** Startup diagnostics remain enabled for developer replay tools, optional for public serving. */
	private boolean verbose=true;
	private int numAggInputs;
	private int[] lastNcServe=new int[5];
	private int[] lastAntiObs=new int[TRNA_ANTICODON_OBS];
	private int[] lastAntiServe=new int[TRNA_ANTICODON_OBS];
	private long anticodonRows, anticodonMissingRows;
	private boolean subnetFamset=false;
	private boolean subnetAnticodon=false, subnetRrna=false;
	/** True when an aggregator manifest/bundle contains a tRNA anticodon subnet. */
	private boolean aggNeedsAnticodons=false;
	private int[] subsetRanks;
	private boolean[] subsetMask;
	private int[] lastFamObs;
	private final IntHashMap nativeFamTotal=new IntHashMap();
	private final IntHashMap nativeFamTotal_2=new IntHashMap();//unshredded-cache counterpart
	private final IntHashMap nativeAntiTotal=new IntHashMap();
	private final IntHashMap nativeAntiTotal_2=new IntHashMap();//unshredded-cache counterpart
	private int subnetInputs;
	private int lastTarget;
	private boolean lastUsedUnshredded;//true when lastTarget's bin drew clean contigs from cache2
	private long lastTotalBp;
	private double lastMixtureComp, lastMixtureCont;
	private String lastSpikeClass="ordinary";
	private int[] lastNcObs=new int[5];
	private int[] keptRanks;
	private long n=400000, valn=40000, seed=1, splitSeed=1;
	private boolean splitSeedSet=false;//when false, the organism split uses seed (byte-identical to pre-splitseed behavior)
	private double valfrac=0.10, mixComp=0.5, mixCont=0.5, cleanSpike=0.15, multiContamProb=0.15, sameFamProb=0.70;
	private int minlen=0;
	private int enc=ENC_RATIO;

	private int numFam=8000, numPhyla, numInputs;
	private TaxTree tree;

	// byTid replacement: dense tid->index (IntHashMap can't hold list values) + a
	// parallel growable list-of-lists indexed by that dense int.
	private final IntHashMap tidToIdx=new IntHashMap();
	private final ArrayList<ArrayList<Contig>> contigLists=new ArrayList<ArrayList<Contig>>();
	private int[] domainIdxArr=new int[64];//parallel to tidToIdx's dense index (was tid2domainIdx)

	private final IntLongHashMap genomeSize=new IntLongHashMap();
	// Benchmark mode (off unless benchtruth=/benchmanifest= given): emit held-out synthetic bins as a
	// truth table + contig-name manifest (+ optional aggregator vectors) so CheckM1/CheckM2 and our net
	// score the IDENTICAL bins. Additive - existing output paths and the differential gate are untouched.
	private boolean benchMode=false;
	private String benchManifestFile=null, benchTruthFile=null, benchVecFile=null, benchVecMonoFile=null;
	private long benchBins=500;
	private ArrayList<Contig> benchNative=new ArrayList<Contig>();
	private ArrayList<Contig> benchForeign=new ArrayList<Contig>();
	// Isolate/high-quality training spike (Brian 2026-08-24): fraction of bins drawn perfect / near-perfect.
	private double perfectFrac=0.0, nearPerfectFrac=0.0;
	// ADDITIONAL isolate spike (Brian 2026-08-26, MAG-QC v4b): fraction of the BASE count N, appended
	// on top after the normal draw completes -- NOT a replacement fraction like perfectFrac/nearPerfectFrac
	// above (which shrink the normal-distribution population). Default 0 (off): the append loop in
	// writeSet then runs zero iterations and draws zero extra Random values, so existing output is
	// byte-identical when these are unset -- this is purely additive.
	private double extraPerfectFrac=0.0, extraNearPerfectFrac=0.0;
	// Which output set(s) the additional spike applies to: train, val, or both. UMP45's call
	// (2026-08-26): both -- val draws its perfect/near-perfect bins from the VAL org pool (writeSet
	// runs per-set, so this is leakage-safe by construction), mirroring the ~9% perfect fraction in
	// the 550-bin benchmark so epoch/model selection tracks the real deployment target. Brian can
	// override to train-only via extraspikeset=train, a one-flag change.
	private String extraSpikeSet="both";
	private final IntLongHashMap recoverable=new IntLongHashMap();
	// UNSHREDDED-genome isolate spike (Brian 2026-08-26): a second, optional per-contig cache
	// (unshreddedcache=) of REAL whole-genome CallGenes output, loaded only when set. Additive
	// exactly like extraPerfectFrac above (default 0 -> zero iterations, byte-identical output
	// when unset), but draws FORCE_PERFECT_UNSHREDDED bins from THIS cache's contigs instead of
	// the shredded one -- a "perfect" bin with zero shred-boundary artifacts.
	private String unshreddedCacheFile=null;
	private double extraUnshreddedPerfectFrac=0.0;
	private final IntHashMap tidToIdx2=new IntHashMap();
	private final ArrayList<ArrayList<Contig>> contigLists2=new ArrayList<ArrayList<Contig>>();
	private final IntLongHashMap recoverable2=new IntLongHashMap();
	private final IntHashMap nativeNcR16=new IntHashMap(), nativeNcR23=new IntHashMap(), nativeNcR5=new IntHashMap(),
		nativeNcRother=new IntHashMap(), nativeNcTrna=new IntHashMap();
	private final IntHashMap nativeNcR16_2=new IntHashMap(), nativeNcR23_2=new IntHashMap(),
		nativeNcR5_2=new IntHashMap(), nativeNcRother_2=new IntHashMap(), nativeNcTrna_2=new IntHashMap();
	private final IntHashMap tid2phylumIdx=new IntHashMap();
	private final IntHashMap tid2family=new IntHashMap();
	private final HashMap<String,Integer> phylumIndex=new HashMap<String,Integer>();
	private ArrayList<String> phylumList;
	private String[] nOverN1;
	private String[] excessArr;
	private double[] avgCopy;
	// avg~0.5 normalization scales (corpus reference over usable orgs; persisted with the nets). Wired into
	// the emit methods in the one-pass vector-layout rewrite; computed by computeNormScales().
	private double scaleLogBp=1, scaleLogCds=1, scaleGlen=1, scaleGlenStd=1;

	private static final int NUM_GLOBALS=13;
	private static final int NCRNA_OBS=5;
	private static final int RRNA_OBS=4;
	private static final int STRUCTURAL_ANTICODON_OBS=64;
	private static final int TRNA_ANTICODON_OBS=STRUCTURAL_ANTICODON_OBS+1;
	private static final int DOMAINS=8, DOMAIN_OTHER=7;
	private static final int N_STR_MAX=4096;
	private static final int ENC_RATIO=0, ENC_LOG=1, ENC_RAW=2, ENC_TWO=3, ENC_NORM=4;
	private static final int RAW_CAP=32, LOG_CAP=64, EXC_CAP=16;
	private static final double COMP_MIN=0.10, CONT_MAX=0.50, LOG2=Math.log(2);
	/** Cap for the per-subnet and pooled obs/pred ratios (duplication saturates at 2x). */
	private static final double RATIO_CAP=2.0;
	/** Indices into the FROZEN shared-context scalar block (Brian 2026-08-24 vector-layout
	 *  rebuild, magqc_rebuild_20260824.plan "FROZEN VECTOR LAYOUT"): 9 scalars, standard on
	 *  EVERY net (subnets + aggregator + monolith) - the old per-net optional sndomain/
	 *  snhhcaga/sngenelen/snbinscaled/sncodingaffine representation flags are RETIRED, there
	 *  is now exactly one encoding per feature. Domain one-hot(8) is the 10th shared-context
	 *  item but is categorical, not scalar, so it is appended separately right after this
	 *  block (see appendContext) rather than living in this array. */
	private static final int CTX_SIZE=0, CTX_GC=1, CTX_CODING=2, CTX_GENES=3, CTX_L2RICH=4,
		CTX_GLEN=5, CTX_GLENSTD=6, CTX_HH=7, CTX_CAGA=8, CTX_N=9;
	/** Shared-context wire width: the 9 scalars above plus the domain one-hot(8). Package-visible:
	 *  MagQCTool's phase-0 validation needs it to check a six-mode bundle's expectedInputs without
	 *  reimplementing the formula (design revision2 §4.3 step 2). */
	static final int SHARED_CONTEXT_WIDTH=CTX_N+DOMAINS;
	private static final int[] EMPTY=new int[0];
}
