package prot.cluster;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.concurrent.locks.ReadWriteLock;

import ddl.DDLFamilySearch.StructuralInputs;
import ddl.DDLKeyTieBreak;
import ddl.DDLPentamerBlacklist;
import ddl.FamilySketchBuilder;
import ddl.TailPoolIndex;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prot.FamilyKeyPolicyInputs;
import prot.FamilyCalibration.VerificationMode;
import prot.ProteinSubstitutionMatrix;
import shared.Resources;
import shared.Shared;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.IntList;
import template.Accumulator;
import template.ThreadWaiter;

/** Native Stage-B collection and original-member union rebuilding. Mutual mode
 * replaces consensus clustering with repeated frozen reciprocal pairs, retaining
 * every original member. The legacy union mode remains available for old recipes.
 * A merged seed does not claim that every original member passed it: reports
 * distinguish accepted consensus proposals from retained provisional seeds.
 */
public final class DdlClustStageB{

	/** Reports expected resource absence after normal cleanup; malformed inputs still propagate. */
	public static void main(String[] args) throws Exception{
		try{
			final DdlClustStageB tool=new DdlClustStageB(args);
			try{tool.process();}
			finally{Shared.closeStream(tool.outstream);}
		}catch(Resources.MissingResourceException e){
			System.err.println(e.getMessage());
			System.exit(1);
		}
	}

	static LayerResult runResidentLayer(final String[] args,final LiveProbeState state,final ClusteringInput input,
			final FamilySketchBuilder builder,final int[] fissionParentIds)throws Exception{
		return runResidentLayer(args,state,input,builder,fissionParentIds,null);
	}

	static LayerResult runResidentLayer(final String[] args,final LiveProbeState state,final ClusteringInput input,
			final FamilySketchBuilder builder,final int[] fissionParentIds,final ProteinSubstitutionMatrix matrix)throws Exception{
		final DdlClustStageB tool=new DdlClustStageB(args,true,matrix);
		return tool.runResidentLayer(state,input,builder,fissionParentIds);
	}

	public DdlClustStageB(final String[] args){
		this(args,false,null);
	}

	private DdlClustStageB(String[] args,final boolean residentLayer_){
		this(args,residentLayer_,null);
	}

	private DdlClustStageB(String[] args,final boolean residentLayer_,final ProteinSubstitutionMatrix residentMatrix){
		residentLayer=residentLayer_;
		final PreParser pre=new PreParser(args, getClass(), false);
		args=pre.args; outstream=pre.outstream;
		final Parser parser=new Parser(); parser.overwrite=true;
		AlignmentAnchor parsedAlignmentAnchor=null;
		Boolean parsedPinRight=null, parsedPinLeft=null;
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String key=pair[0].toLowerCase(Locale.ROOT);
			final String value=pair.length>1 && !pair[1].equalsIgnoreCase("null") ? pair[1] : null;
			if(key.equals("mode")){mode=value==null ? null : value.toLowerCase(Locale.ROOT);}
			else if(key.equals("phase")){livePhase=value==null ? null : value.toLowerCase(Locale.ROOT);livePhaseSet=true;}
			else if(key.equals("selection")){
				require("predicted".equalsIgnoreCase(value) || "measured".equalsIgnoreCase(value),
					"selection must be predicted or measured");
				rankByPrediction="predicted".equalsIgnoreCase(value);liveSelectionSet=true;
			}
			else if(key.equals("codecs")){codecs=value;}
			else if(key.equals("initialfraction")){initialFraction=Double.parseDouble(value);initialFractionSet=true;}
			else if(key.equals("pricedfraction")){pricedFraction=Double.parseDouble(value);pricedFractionSet=true;}
			else if(key.equals("pricemode")){
				require("floor".equalsIgnoreCase(value) || "fraction".equalsIgnoreCase(value),"pricemode must be floor or fraction");
				priceMode=LiveProbeLayer.PriceMode.valueOf(value.toUpperCase(Locale.ROOT));priceModeSet=true;
			}
			else if(key.equals("setupslotcost")){setupSlotCost=Parse.parseIntKMG(value);setupSlotCostSet=true;}
			else if(key.equals("belowmode")){
				require("constant".equalsIgnoreCase(value) || "linear".equalsIgnoreCase(value) || "sqrt".equalsIgnoreCase(value) ||
					"prior".equalsIgnoreCase(value) || "priorscaled".equalsIgnoreCase(value),
					"belowmode must be constant, linear, sqrt, prior or priorscaled");
				belowMode=LiveProbeLayer.BelowMode.valueOf(value.toUpperCase(Locale.ROOT));belowModeSet=true;
			}
			else if(key.equals("edgefraction")){edgeFraction=Double.parseDouble(value);edgeFractionSet=true;}
			else if(key.equals("sizelambda")){sizeLambda=Double.parseDouble(value);sizeLambdaSet=true;}
			else if(key.equals("fissiontopfraction")){fissionTopFraction=Double.parseDouble(value);eligibilityFractionSet=true;}
			else if(key.equals("fusionbottomfraction")){fusionBottomFraction=Double.parseDouble(value);eligibilityFractionSet=true;}
			else if(key.equals("stageacodecs")){stageACodecs=value;}
			else if(key.equals("stageamembers")){stageAMembers=value;}
			else if(key.equals("stagebcodec")){stageBCodec=value;}
			else if(key.equals("stagebmembers")){stageBMembers=value;}
			else if(key.equals("partitionmap")){partitionMap=value;}
			else if(key.equals("roster")){roster=value;}
			else if(key.equals("report")){report=value;}
			else if(key.equals("passesreport")){passesReport=value;}
			else if(key.equals("fusionreport")){fusionReport=value;}
			else if(key.equals("fissionparents")){fissionParents=value;}
			else if(key.equals("fusionproductsout")){fusionProductsOut=value;}
			else if(key.equals("metrics")){liveMetrics=value;}
			else if(key.equals("sizehist")){liveSizeHist=value;}
			else if(key.equals("stateout")){stateOut=value;}
			else if(key.equals("currentstate")){currentState=Parse.parseBoolean(value);currentStateSet=true;}
			else if(key.equals("nextfamilyid")){nextFamilyId=Parse.parseIntKMG(value);nextFamilyIdSet=true;}
			else if(key.equals("slotfraction")){slotFraction=Double.parseDouble(value);slotFractionSet=true;}
			else if(key.equals("members")){membersOut=value;}
			else if(key.equals("probes")){probePlan=value;probePlanSet=true;}
			else if(key.equals("rebuildreport")){rebuildReport=value;}
			else if(key.equals("auditreport")){auditReport=value;}
			else if(key.equals("featuresreport")){featuresReport=value;}
			else if(key.equals("widefeatures")){wideFeatures=value;}
			else if(key.equals("featurecontext")){featureContext=value;}
			else if(key.equals("featurelayer")){featureLayer=Parse.parseIntKMG(value);featureLayerSet=true;}
			else if(key.equals("fitreport")){fitReport=value;}
			else if(key.equals("refinementreport")){refinementReport=value;}
			else if(key.equals("localrefinerounds")){localRefineRounds=Parse.parseIntKMG(value);localRefineRoundsSet=true;}
			else if(key.equals("selectedseeds")){selectedSeeds=value;localConfigured=true;}
			else if(key.equals("baselinemembers")){baselineMembers=value;localConfigured=true;}
			else if(key.equals("memberranking")){memberRanking=Parse.parseBoolean(value);memberRankingSet=true;}
			else if(key.equals("memberhydroweight")){memberHydroWeight=MemberRankInputs.checkedHydroWeight(Double.parseDouble(value));}
			else if(key.equals("selectedcount")){selectedCount=Parse.parseIntKMG(value);localConfigured=true;}
			else if(key.equals("fissionprobes")){fissionProbes=Parse.parseIntKMG(value);localPlanConfigured=true;}
			else if(key.equals("fusionprobes")){fusionProbes=Parse.parseIntKMG(value);localPlanConfigured=true;}
			else if(key.equals("probeseed")){probeSeed=Long.parseLong(value);probeSeedSet=true;localPlanConfigured=true;}
			else if(key.equals("probestratified")){probeStratified=Parse.parseBoolean(value);localPlanConfigured=true;}
			else if(key.equals("minlen")){localMinLength=Parse.parseIntKMG(value);localConfigured=true;}
			else if(key.equals("minhits")){localMinimumHits=Parse.parseIntKMG(value);localFlex=false;localConfigured=true;}
			else if(key.equals("passhits")){
				require(value!=null,"passhits requires flex or a nonnegative integer");
				localFlex="flex".equalsIgnoreCase(value);
				if(!localFlex){localMinimumHits=Parse.parseIntKMG(value);}
				localConfigured=true;
			}
			else if(key.equals("maxhydrodelta")){localHydroDelta=Double.parseDouble(value);localConfigured=true;}
			else if(key.equals("mergemaxsize")){mergeMaxSize=Parse.parseIntKMG(value); mergeMaxSizeSet=true;}
			else if(parseStructural(key,value)){}
			else if(key.equals("shortlist")){mergeShortlist=Parse.parseIntKMG(value);}
			else if(key.equals("frozenindexthreads")){frozenIndexThreads=Parse.parseIntKMG(value);}
			else if(key.equals("mergeminid")){mergeIdentity=Double.parseDouble(value);}
			else if(key.equals("mergeminratio")){mergeRatio=Double.parseDouble(value);}
			else if(key.equals("blacklist")){blacklist=value;}
			else if(key.equals("blacklistsha80")){blacklistSha80=value;}
			else if(key.equals("blacklistk")){blacklistK=Parse.parseIntKMG(value);}
			else if(key.equals("blacklistmode")){blacklistMode=value;}
			else if(key.equals("blacklistprovenance")){blacklistProvenance=value;}
			else if(key.equals("blacklistwords") || key.equals("blacklist_words")){blacklistWords=Parse.parseIntKMG(value);}
			else if(key.equals("verification")){verification=parseVerification(value);}
			else if(key.equals("keytiebreak")){keyPolicy=DDLKeyTieBreak.Policy.parse(value); keyPolicySpecified=true;}
			else if(key.equals("consensusbonus")){consensusMemberBonus=Parse.parseKMG(value);consensusBonusSpecified=true;}
			else if(key.equals("rawmanifest")){rawManifest=value;}
			else if(key.equals("expectedrawtails")){expectedRawTails=Parse.parseKMG(value);}
			else if(key.equals("rawminlen")){rawMinLength=Parse.parseIntKMG(value);}
			else if(key.equals("rawmaxlen")){rawMaxLength=Parse.parseIntKMG(value);}
			else if(key.equals("rawtailcrop")){rawTailCrop=Parse.parseIntKMG(value);}
			else if(key.equals("minspan")){minSpan=Parse.parseIntKMG(value);}
			else if(key.equals("minpaired")){minPaired=Parse.parseIntKMG(value);}
			else if(key.equals("maxalignlen")){maxAlignLength=Parse.parseIntKMG(value);}
			else if(key.equals("minid")){minIdentity=Double.parseDouble(value);}
			else if(key.equals("minblosumratio")){minRatio=Double.parseDouble(value);}
			else if(key.equals("matrix")){matrixName=value;}
			else if(key.equals("matrixfile")){matrixFile=value;}
			else if(key.equals("pinright")){parsedPinRight=Parse.parseBoolean(value);}
			else if(key.equals("pinleft")){parsedPinLeft=Parse.parseBoolean(value);}
			else if(key.equals("anchorterminus")){
				if(value==null){throw new IllegalArgumentException("anchorterminus requires left, right or none");}
				parsedAlignmentAnchor=parseAlignmentAnchor(value);
			}
			else if(key.equals("trimleft")){trimLeft=Parse.parseBoolean(value);}
			else if(key.equals("alignmentmode")){explicitAlignmentMode=TailAlignmentGate.parseRightMode(value);}
			else if(key.equals("consensustrimleft")){consensusTrimLeft=Parse.parseBoolean(value);}
			else if(key.equals("trimdepth")){trimDepth=Float.parseFloat(value);}
			else if(key.equals("denominator")){
				if("aligned".equalsIgnoreCase(value)){denominator=TailAlignmentGate.Denominator.ALIGNED_QUERY;}
				else if("full".equalsIgnoreCase(value)){denominator=TailAlignmentGate.Denominator.FULL_QUERY;}
				else{throw new IllegalArgumentException("denominator must be aligned or full");}
			}else if(parser.parse(arg, key, value)){}
			else if(Tools.looksLikeInputSequenceStream(arg)){parser.in1=arg;}
			else{throw new IllegalArgumentException("Unknown parameter "+arg);}
		}
		in=parser.in1; out=parser.out1; overwrite=parser.overwrite;
		alignmentAnchor=parsedAlignmentAnchor==null ?
			resolveAlignmentAnchor(parsedPinRight,parsedPinLeft) :
			resolveAlignmentAnchor(parsedAlignmentAnchor,parsedPinRight,parsedPinLeft);
		// These entrypoints reconstruct actual original memberships before ranking.
		if(!memberRankingSet){memberRanking=supportsMemberRanking(mode);}
		require(memberRanking || memberHydroWeight==MemberRankInputs.DEFAULT_HYDRO_WEIGHT,
			"Nondefault memberhydroweight requires memberranking=t");
		workerThreads=Shared.threads();
		require(!initialFractionSet || "select".equals(mode),"initialfraction requires mode=select");
		require(frozenIndexThreads>=1 && frozenIndexThreads<=workerThreads &&
			(frozenIndexThreads==1 || "mutual".equals(mode) || "layercake".equals(mode) || "probe".equals(mode) || "localprobe".equals(mode) || "localplan".equals(mode) || "liveprobe".equals(mode)),
			"frozenindexthreads must lie in [1,t]; parallel index construction requires mode=mutual, layercake, probe, localprobe, localplan or liveprobe");
		require(parser.in2==null && parser.out2==null && parser.outsingle==null,
			"Stage-B accepts one complete input stream and no paired/singleton outputs");
		require(parser.qfin1==null && parser.qfin2==null && parser.qfout1==null && parser.qfout2==null,
			"Stage-B does not support quality streams");
		require(parser.extin==null && parser.extout==null, "Stage-B does not support renamed paired extensions");
		require(!parser.append, "Stage-B outputs must be freshly written, not appended");
		require(parser.maxReads<0, "Stage-B must see every parent ordinal; reads= is unsupported");
		if(!Float.isFinite(trimDepth) || trimDepth<0 || trimDepth>1){throw new IllegalArgumentException("trimdepth must be in [0,1]");}
		if(maxAlignLength<1 || minSpan<1 || minSpan>maxAlignLength){throw new IllegalArgumentException("Require 1 <= minspan <= maxalignlen");}
		if(!Double.isFinite(minIdentity) || minIdentity<0 || minIdentity>1){throw new IllegalArgumentException("minid must be in [0,1]");}
		if(!Double.isFinite(minRatio)){throw new IllegalArgumentException("minblosumratio must be finite");}
		if(alignmentAnchor==AlignmentAnchor.NONE){throw new IllegalArgumentException("Stage-B union rebuilds require pinned alignment");}
		if(mergeMaxSize<1 || mergeShortlist<1 || mergeShortlist==Integer.MAX_VALUE){
			throw new IllegalArgumentException("Mutual merging needs positive mergemaxsize and a shortlist with room to exclude self");
		}
		if(!Double.isFinite(mergeIdentity) || mergeIdentity<0 || mergeIdentity>1 || !Double.isFinite(mergeRatio)){
			throw new IllegalArgumentException("Mutual merge identity must be in [0,1] and ratio must be finite");
		}
		if("select".equals(mode)){
			require(!nonBudgetStructuralConfigured && !probePlanSet,"Global selection accepts familybudget, not structural operations");
			ClusterSeedSelector.requestedCount(familyBudget,initialFraction);
		}else if("layercake".equals(mode)){
			require(scheduleMode!=null,"mode=layercake requires explicit schedule=none, fission, fusion or layercake");
			require(!probePlanSet,"probes requires mode=probe");
			require(!mergeMaxSizeSet,"mode=layercake uses fusionsizeratio and the shared objective, not mergemaxsize");
			structuralSettings(0); structuralOptions(); structuralObjective();
		}else if("probe".equals(mode) || "localprobe".equals(mode)){
			require(probePlanSet && probePlan!=null,"Probe modes require probes= with a persisted candidate table");
			require(!scheduleControlsSet && !objectiveWeightsSet && !mergeMaxSizeSet,
				"Raw probes do not accept schedule controls, fitted objective weights/margin or mergemaxsize");
			require(Double.isFinite(temperature) && temperature>=0,"Probe temperature must be finite and nonnegative");
			structuralOptions();structuralObjective();
		}else if("liveprobe".equals(mode)){
			final boolean combined="combined".equals(livePhase),layer="layer".equals(livePhase);
			if((combined || layer) && !fissionMinSizeSet){fissionMinSize=1;}
			if(!localRefineRoundsSet){localRefineRounds=2;}
			require(localRefineRounds>=1 && localRefineRounds<=2,"Live probes require localrefinerounds=1 or 2");
			require("fission".equals(livePhase) || "fusion".equals(livePhase) || "cycle".equals(livePhase) || combined || layer,
				"mode=liveprobe requires explicit phase=fission, fusion, cycle, combined or layer");
			require(!gainMarginSet || !liveMarginSlot && Double.isFinite(gainMargin) && gainMargin==0,
				"Live probes require symmetric slotfraction thresholds; omit gainmargin or use gainmargin=0 (not slot/nonzero offsets)");
			require(!combined || fissionMinSize==1,"Combined probes require fissionminsize=1 so every active family is requested");
			require(!layer || fissionMinSize==1 && !mergeMaxSizeSet,
				"Live layers require fissionminsize=1 and no mergemaxsize cap; feasibility is decided by the native proposals");
			require(!probePlanSet && !scheduleControlsSet && !fittedWeightsSet && !temperatureSet,
				"Live probes do not accept persisted requests, schedule controls, fitted coefficients or temperature");
			structuralOptions();structuralObjective();
		}else if("localplan".equals(mode)){
			require(!probePlanSet && !scheduleControlsSet && !objectiveWeightsSet && !mergeMaxSizeSet,
				"Local request generation uses explicit random counts, not persisted requests, schedule controls or fitted weights");
			require(probeSeedSet && fissionProbes>=0 && fusionProbes>=0 && (long)fissionProbes+fusionProbes>0 &&
				(long)fissionProbes+fusionProbes<=Integer.MAX_VALUE-8,"Local request generation requires probeseed and nonnegative fissionprobes/fusionprobes with positive bounded total");
			structuralOptions();
		}else{require(!structuralConfigured && !probePlanSet,"Structural schedule controls require mode=layercake; probe requests require mode=probe");}
		if("localprobe".equals(mode) || "localplan".equals(mode) || "liveprobe".equals(mode)){
			require(selectedCount>0 && localMinLength>=minSpan,"Local probes require positive selectedcount and minlen at least minspan");
			new CandidateFilter(localMinimumHits,0,localHydroDelta,localHydroDelta,localFlex);
		}else{require(!localConfigured,"Selected-baseline and local screen flags require mode=localprobe, localplan or liveprobe");}
		require(rebuildReport==null || "localprobe".equals(mode),"rebuildreport requires mode=localprobe");
		require(auditReport==null || "localprobe".equals(mode),"auditreport requires mode=localprobe");
		final boolean combined="liveprobe".equals(mode) && "combined".equals(livePhase);
		final boolean layer="liveprobe".equals(mode) && "layer".equals(livePhase);
		require(!pricedFractionSet || layer,"pricedfraction requires liveprobe phase=layer");
		require(!currentStateSet && !nextFamilyIdSet && stateOut==null || layer,
			"currentstate, nextfamilyid and stateout require liveprobe phase=layer");
		require(fissionParents==null && fusionProductsOut==null || layer,
			"fissionparents and fusionproductsout require liveprobe phase=layer");
		require(!layer || stateOut!=null,"Live layers require stateout to preserve the fresh-family ID cursor");
		require(!layer || residentLayer || (currentState ? selectedSeeds==null && nextFamilyIdSet && nextFamilyId>=0 :
			selectedSeeds!=null && !nextFamilyIdSet),
			"A current layer requires currentstate=t, nextfamilyid and no selectedseeds; an initial layer requires selectedseeds and no nextfamilyid");
		require(!slotFractionSet || "liveprobe".equals(mode),"slotfraction requires mode=liveprobe");
		require(!memberRanking || supportsMemberRanking(mode),
			"memberranking requires mode=liveprobe, localplan, mutual, probe, localprobe or layercake");
		require(Double.isFinite(slotFraction) && slotFraction>=0,"slotfraction must be finite and nonnegative");
		require(wideFeatures==null || combined || layer,"widefeatures requires liveprobe phase=combined or layer");
		require(wideFeatures!=null || featureContext==null && !featureLayerSet,"featurecontext/featurelayer require widefeatures");
		if(wideFeatures!=null){wideOutput=new WideProbeFeatures.Output(wideFeatures,featureContext,featureLayer);}
		require(!liveSelectionSet || combined,"selection requires liveprobe phase=combined");
		require(featuresReport==null || "localprobe".equals(mode) && auditReport!=null || combined,
			"featuresreport requires localprobe with auditreport, or liveprobe phase=combined");
		require(fitReport==null || combined,"fitreport requires liveprobe phase=combined");
		require(!liveMarginSlot || combined,"gainmargin=slot requires liveprobe phase=combined");
		require(localRefineRounds>=0 && localRefineRounds<=2,"localrefinerounds must be0,1or2");
		require(localRefineRounds==0 || "liveprobe".equals(mode) || "localprobe".equals(mode) && auditReport!=null && refinementReport!=null,
			"Local refinement requires mode=liveprobe or mode=localprobe with auditreport and refinementreport");
		require(refinementReport==null || "localprobe".equals(mode),"refinementreport requires mode=localprobe");
		require(refinementReport==null || localRefineRounds>0,"refinementreport requires positive localrefinerounds");
		require(!livePhaseSet || "liveprobe".equals(mode),"phase requires mode=liveprobe");
		require(fusionReport==null || "liveprobe".equals(mode) && ("cycle".equals(livePhase) || layer),
			"fusionreport requires mode=liveprobe phase=cycle or layer");
		require(liveMetrics==null && liveSizeHist==null || "liveprobe".equals(mode) && ("cycle".equals(livePhase) || combined || layer),
			"fusionreport, metrics and sizehist require a supported live phase; metrics and sizehist support cycle, combined or layer");
		require(!localPlanConfigured || "localplan".equals(mode),"Random request counts/seed require mode=localplan");
		require(!priceModeSet && !setupSlotCostSet && !edgeFractionSet && !belowModeSet && !sizeLambdaSet || layer,
			"pricemode, setupslotcost, edgefraction, belowmode and sizelambda require liveprobe phase=layer");
		require(!eligibilityFractionSet || layer,"fissiontopfraction and fusionbottomfraction require liveprobe phase=layer");
		require(Double.isFinite(sizeLambda) && sizeLambda>=0,"sizelambda must be finite and nonnegative");
		require(Double.isFinite(fissionTopFraction) && fissionTopFraction>=0 && fissionTopFraction<=1 &&
			Double.isFinite(fusionBottomFraction) && fusionBottomFraction>=0 && fusionBottomFraction<=1,
			"eligibility fractions must be finite and in [0,1]");
		require(fissionParents==null || fissionTopFraction==1,
			"fissionparents requires fissiontopfraction=1; tagged fusion products are already the eligible fission frontier");
		require(Double.isFinite(edgeFraction) && edgeFraction>0,"edgefraction must be finite and positive");
		require(!setupSlotCostSet || setupSlotCost>=0,"setupslotcost must be nonnegative");
		require(!layer || priceMode!=LiveProbeLayer.PriceMode.FLOOR || setupSlotCostSet,
			"Floor pricing requires explicit setupslotcost from the preserved step3 continuation");
		final boolean defaultMatrixFile=matrixFile==null && ProteinSubstitutionMatrix.DEFAULT_NAME.equals(matrixName);
		if(defaultMatrixFile){matrixFile=ProteinSubstitutionMatrix.defaultResourcePath();}
		if(matrixFile!=null && !Files.isRegularFile(Paths.get(matrixFile))){
			throw new IllegalArgumentException("Matrix resource must be a regular disk file: "+matrixFile);
		}
		matrix=residentMatrix==null ? ProteinSubstitutionMatrix.load(matrixName,defaultMatrixFile ? null : matrixFile) : residentMatrix;
		if(!matrix.name.equals(matrixName)){
			throw new IllegalArgumentException("Resident Stage-B matrix differs from the configured selector: "+matrixName);
		}
	}

	/** Only modes that choose candidates use the member-derived ranking pivot. */
	private static boolean supportsMemberRanking(final String mode){
		return "liveprobe".equals(mode) || "localplan".equals(mode) || "mutual".equals(mode) ||
			"probe".equals(mode) || "localprobe".equals(mode) || "layercake".equals(mode);
	}

	/** Structural controls cannot silently alter an existing collection or merge recipe. */
	private boolean parseStructural(final String key,final String value){
		if(key.equals("schedule")){scheduleMode=LayerCakeSchedule.Mode.valueOf(value.toUpperCase(Locale.ROOT));}
		else if(key.equals("temperature")){temperature=Double.parseDouble(value);temperatureSet=true;}
		else if(key.equals("mintemperature")){minimumTemperature=Double.parseDouble(value);}
		else if(key.equals("cooling")){cooling=Double.parseDouble(value);}
		else if(key.equals("familybudget")){familyBudget=Parse.parseIntKMG(value);}
		else if(key.equals("structuraltarget")){structuralTarget=Double.parseDouble(value);}
		else if(key.equals("maxcycles")){maxCycles=Parse.parseIntKMG(value);}
		else if(key.equals("largesize")){largeSize=Parse.parseIntKMG(value);}
		else if(key.equals("gainmargin")){
			liveMarginSlot="slot".equalsIgnoreCase(value);
			if(!liveMarginSlot){gainMargin=Double.parseDouble(value);}
			gainMarginSet=true;
		}
		else if(key.equals("cohesionscale")){cohesionScale=Double.parseDouble(value);}
		else if(key.equals("separationscale")){separationScale=Double.parseDouble(value);}
		else if(key.equals("sizescale")){sizeScale=Double.parseDouble(value);}
		else if(key.equals("sizeform")){sizeForm=LocalFamilyObjective.SizeForm.valueOf(value.toUpperCase(Locale.ROOT));}
		else if(key.equals("fissionminsize")){fissionMinSize=Parse.parseIntKMG(value);fissionMinSizeSet=true;}
		else if(key.equals("minchildsize")){minChildSize=Parse.parseIntKMG(value);}
		else if(key.equals("memberneighbors")){memberNeighbors=Parse.parseIntKMG(value);}
		else if(key.equals("memberminhits")){memberMinHits=Parse.parseIntKMG(value);}
		else if(key.equals("evidenceneighbors")){evidenceNeighbors=Parse.parseIntKMG(value);}
		else if(key.equals("fusionsizeratio")){fusionSizeRatio=Double.parseDouble(value);}
		else{return false;}
		scheduleControlsSet|=key.equals("schedule") || key.equals("mintemperature") || key.equals("cooling") || key.equals("structuraltarget") || key.equals("maxcycles");
		objectiveWeightsSet|=key.equals("gainmargin") || key.equals("cohesionscale") || key.equals("separationscale") || key.equals("sizescale");
		fittedWeightsSet|=key.equals("cohesionscale") || key.equals("separationscale") || key.equals("sizescale");
		structuralConfigured=true;
		nonBudgetStructuralConfigured|=!key.equals("familybudget");
		return true;
	}

	private LayerCakeSchedule.Settings structuralSettings(final long inputs){
		return new LayerCakeSchedule.Settings(scheduleMode,temperature,minimumTemperature,cooling,
			familyBudget,structuralTarget,inputs,maxCycles);
	}
	private LayerCakeClustering.Options structuralOptions(){
		return new LayerCakeClustering.Options(workerThreads,frozenIndexThreads,fissionMinSize,minChildSize,
			memberNeighbors,memberMinHits,evidenceNeighbors,mergeShortlist,fusionSizeRatio,consensusTrimLeft,trimDepth,familyBudget,localRefineRounds,pricedFraction);
	}
	private LocalFamilyObjective structuralObjective(){
		if("liveprobe".equals(mode)){
			// Construction requires a legacy objective object, whose margin must be positive.
			// It is NOT the live decision function: run() receives gainMargin separately
			// and compares measured coverage gains against the phase-signed slot cost.
			return new LocalFamilyObjective(new LocalFamilyObjective.Parameters(largeSize,1,0,0,0,sizeForm));
		}
		return new LocalFamilyObjective(new LocalFamilyObjective.Parameters(largeSize,gainMargin,
			cohesionScale,separationScale,sizeScale,sizeForm));
	}

	public void process() throws Exception{
		configureProteinInput();
		if("collect".equals(mode)){collect();}
		else if("select".equals(mode)){ClusterSeedSelector.write(split(codecs,"codecs"),familyBudget,initialFraction,out,roster,report,overwrite);}
		else if("union".equals(mode)){union();}
		else if("localprobe".equals(mode)){localProbe();}
		else if("localplan".equals(mode)){localPlan();}
		else if("liveprobe".equals(mode)){liveProbe();}
		else if("mutual".equals(mode) || "layercake".equals(mode) || "probe".equals(mode)){mutual();}
		else{throw new IllegalArgumentException("mode must be collect, select, union, mutual, layercake, probe, localprobe, localplan or liveprobe");}
	}

	/** One phase or one ordered fission/fusion cycle publishes current models and exact winners. */
	private void liveProbe() throws Exception{
		require(in!=null && out!=null && codecs!=null && (currentState || selectedSeeds!=null) && baselineMembers!=null &&
			membersOut!=null && report!=null,
			"mode=liveprobe requires in, out, codecs, selectedseeds (unless currentstate=t), baselinemembers, members and report");
		require(stageACodecs==null && stageAMembers==null && stageBCodec==null && stageBMembers==null &&
			partitionMap==null && roster==null,
			"Live probes use one selected baseline, not legacy Stage-A/union inputs");
		final boolean layer="layer".equals(livePhase),cycle="cycle".equals(livePhase),combined="combined".equals(livePhase);
		require(combined ? featuresReport!=null && fitReport!=null && passesReport!=null && liveMetrics!=null && liveSizeHist!=null :
			cycle || layer ? fusionReport!=null && passesReport!=null && liveMetrics!=null && liveSizeHist!=null : passesReport==null,
			combined ? "Combined probes require featuresreport, fitreport, passesreport, metrics and sizehist" :
			layer ? "Live layer requires fusionreport, passesreport, metrics and sizehist" :
			"Live cycle requires fusionreport, passesreport, metrics and sizehist; single phases do not accept passesreport");
		final String[] paths=split(codecs,"codecs");
		require(paths.length==1,"Live probes require one model codec matching the selected baseline or complete current state");
		final String[] ordinaryOutputs=combined ? new String[]{out,membersOut,report,featuresReport,fitReport,passesReport,liveMetrics,liveSizeHist} :
			layer ? new String[]{out,membersOut,report,fusionReport,passesReport,liveMetrics,liveSizeHist,stateOut} :
			cycle ? new String[]{out,membersOut,report,fusionReport,passesReport,liveMetrics,liveSizeHist} :
			new String[]{out,membersOut,report};
		final int extra=(wideFeatures==null ? 0 : 1)+(fusionProductsOut==null ? 0 : 1);
		final String[] outputs=extra==0 ? ordinaryOutputs : Arrays.copyOf(ordinaryOutputs,ordinaryOutputs.length+extra);
		int outAt=ordinaryOutputs.length;
		if(wideFeatures!=null){outputs[outAt++]=wideFeatures;}
		if(fusionProductsOut!=null){outputs[outAt++]=fusionProductsOut;}
		PathRegistry.prepare().inputs(in,paths[0],selectedSeeds,baselineMembers,fissionParents,blacklist,matrixFile).outputs(outputs);
		require(Tools.testOutputFiles(overwrite,false,false,outputs),"Live probe outputs are unavailable");
		final ClusterSketchCodec.Header h;
		final FamilySketchBuilder builder;
		final CandidateFilter filter;
		final ClusteringInput input;
		LiveProbeState state;
		final TailAlignmentGate memberGate,mergeGate;
		try(PhaseCostLog.Span phase=PhaseCostLog.begin("input_loading_and_state")){
			final ClusterSketchCodec.Loaded source=ClusterSketchCodec.read(paths[0]);
			h=source.header;
			require(!currentState || selectedCount==source.families.length,
				"Current-state selectedcount must equal the entire serialized roster; continuation does not select or prune models");
			builder=familySketchBuilder(h,loadBlacklist(h,outputs),outputs);
			filter=new CandidateFilter(localMinimumHits,0,localHydroDelta,localHydroDelta,localFlex);
			filter.validateWidth(h.width);
			input=ClusteringInput.load(in,-1,false);
			final MemberRankInputs memberInputs=memberRanking ? new MemberRankInputs(input.sequences,memberHydroWeight) : null;
			final LocalProbeBaseline baseline=currentState ?
				memberRanking ? LocalProbeBaseline.loadCurrent(source,baselineMembers,input.sequences,localMinLength,nextFamilyId,memberInputs) :
				LocalProbeBaseline.loadCurrent(source,baselineMembers,input.sequences,localMinLength,nextFamilyId) :
				memberRanking ? LocalProbeBaseline.load(source,selectedSeeds,baselineMembers,input.sequences,localMinLength,selectedCount,memberInputs) :
				LocalProbeBaseline.load(source,selectedSeeds,baselineMembers,input.sequences,localMinLength,selectedCount);
			state=LiveProbeState.fromBaseline(baseline,input.sequences.length);
			if(memberRanking){state=state.withMemberRanking(input.sequences,memberInputs,baseline.memberRankingLayer);}
			if(fissionParents!=null){state=state.withFusionProductIds(readFamilyIds(fissionParents));}
			memberGate=gate();
			mergeGate=new TailAlignmentGate(memberGate.mode,denominator,maxAlignLength,minSpan,mergeIdentity,mergeRatio,null,memberGate.matrix,minPaired);
		}
		if(layer){liveLayer(state,input,builder,memberGate,mergeGate,filter);return;}
		if(combined){liveCombined(state,input,builder,memberGate,mergeGate,filter);return;}
		if(cycle){liveCycle(state,input,builder,memberGate,mergeGate,filter);return;}
		final LayerCakeSchedule.Phase phase=LayerCakeSchedule.Phase.valueOf(livePhase.toUpperCase(Locale.ROOT));
		final LiveProbeClustering.Result result=LiveProbeClustering.runWithFraction(state,input.sequences,builder,memberGate,mergeGate,
			filter,localMinLength,structuralObjective(),structuralOptions(),phase,slotFraction,mergeMaxSizeSet ? mergeMaxSize : Integer.MAX_VALUE);
		final LiveProbeState next=result.state();
		keyInputs.verify();
		writeLiveMembers(next,input);
		ClusterSketchCodec.write(out,overwrite,next.models(),h);
		result.writeReport(report,overwrite);
		outstream.println("Live "+livePhase+" pass: candidates="+result.candidates()+", available="+result.available()+
			", accepted="+result.accepted()+", slot_cost="+result.slotCost()+", threshold="+result.threshold()+
			", covered="+next.assigned()+", residual="+next.residualOrdinals().length+", filtered="+next.filtered()+
			", next_family_id="+next.firstNewId());
		outstream.println("Current-state outputs retain all models and original ordinals; they are not selected-source loader inputs.");
	}

	/** One layer preserves every current template and uses one layer-start B for all waves. */
	private void liveLayer(final LiveProbeState initial,final ClusteringInput input,final FamilySketchBuilder builder,
			final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,final CandidateFilter filter)throws Exception{
		final LayerResult layer=runResidentLayer(initial,input,builder,memberGate,mergeGate,filter,
			fissionParents==null ? null : initial.fusionProductIds());
		final LiveProbeState next=layer.result.state();
		try(PhaseCostLog.Span phase=PhaseCostLog.begin("layer_output")){
			keyInputs.verify();
			writeResidentLayerReports(layer);
			writeLiveMembers(next,input);
			final ClusterSketchCodec.Header h=next.sourceHeader();
			ClusterSketchCodec.write(out,overwrite,next.models(),h);
			if(fusionProductsOut!=null){writeFamilyIds(fusionProductsOut,layer.fusionProductIds);}
			writeLiveState(layer);
		}
		printLiveLayer(layer,initial,next);
	}

	private LayerResult runResidentLayer(final LiveProbeState initial,final ClusteringInput input,
			final FamilySketchBuilder builder,final int[] fissionParentIds)throws Exception{
		require("liveprobe".equals(mode) && "layer".equals(livePhase),
			"Resident Stage-B layers require mode=liveprobe phase=layer");
		require(initial!=null && input!=null && memberRanking==(initial.memberRankingInputs()!=null),
			"Resident Stage-B layer must preserve the effective ranking mode and original cache binding");
		require(!memberRanking || initial.memberRankingInputs().belongsTo(input.sequences),
			"Resident Stage-B ranking requires the identical original population/cache");
		require(!memberRanking || initial.memberRankingInputs().hydroWeight==memberHydroWeight,
			"Resident Stage-B cannot change the population-bound memberhydroweight");
		for(Family family:initial.models()){
			require(memberRanking ? family.memberRankingUses(initial.memberRankingInputs()) : !family.hasMemberRanking(),
				"Resident Stage-B models must agree with the effective ranking mode and original cache binding");
		}
		final ClusterSketchCodec.Header h=initial.sourceHeader();
		final DDLPentamerBlacklist blocked=loadBlacklist(h,passesReport,liveMetrics,liveSizeHist,stateOut);
		validateResidentKeyPolicy(h,builder,blocked);
		require(h.sameKeyPolicy(h.withKeyPolicy(builder.policyContext())),
			"Resident Stage-B builder must preserve the source key policy and raw count domain");
		require(passesReport!=null && liveMetrics!=null && liveSizeHist!=null && stateOut!=null,
			"Resident Stage-B layers require passesreport, metrics, sizehist and stateout");
		require(fissionParentIds==null || fissionTopFraction==1,
			"Explicit resident fission parents require fissiontopfraction=1");
		final TailAlignmentGate memberGate=gate();
		final TailAlignmentGate mergeGate=new TailAlignmentGate(memberGate.mode,denominator,maxAlignLength,minSpan,mergeIdentity,mergeRatio,null,memberGate.matrix,minPaired);
		final CandidateFilter filter=new CandidateFilter(localMinimumHits,0,localHydroDelta,localHydroDelta,localFlex);
		filter.validateWidth(builder.width());
		final LayerResult layer=runResidentLayer(initial,input,builder,memberGate,mergeGate,filter,fissionParentIds);
		try(PhaseCostLog.Span phase=PhaseCostLog.begin("layer_output")){
			keyInputs.verify();
			writeResidentLayerReports(layer);
			writeLiveState(layer);
		}
		printLiveLayer(layer,initial,layer.result.state());
		return layer;
	}

	private LayerResult runResidentLayer(final LiveProbeState initial,final ClusteringInput input,
			final FamilySketchBuilder builder,final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,
			final CandidateFilter filter,final int[] fissionParentIds){
		final long started=System.nanoTime();
		final LiveProbeMetrics metrics=new LiveProbeMetrics(input.sequences,memberGate,workerThreads,familyBudget);
		final LiveProbeMetrics.Snapshot before;
		try(PhaseCostLog.Span phase=PhaseCostLog.begin("metrics_before")){
			before=measureLiveMetrics(metrics,initial.models(),"before");
		}
		final LiveProbeLayer.Outputs reports=new LiveProbeLayer.Outputs(report,fusionReport,passesReport,wideOutput,overwrite);
		final LiveProbeLayer.Result result=LiveProbeLayer.run(initial,input.sequences,builder,memberGate,mergeGate,
			filter,localMinLength,structuralObjective(),structuralOptions(),slotFraction,priceMode,setupSlotCost,edgeFraction,
			belowMode,sizeLambda,fissionTopFraction,fusionBottomFraction,reports,
			fissionParentIds);
		final LiveProbeState next=result.state();
		final LiveProbeMetrics.Snapshot afterFission,afterFusion;
		try(PhaseCostLog.Span phase=PhaseCostLog.begin("metrics_after")){
			afterFission=measureLiveMetrics(metrics,result.afterFissionState().models(),"after_fission");
			afterFusion=measureLiveMetrics(metrics,next.models(),"after_fusion");
		}
		require(result.fissionNetGain()==afterFission.covered()-before.covered() &&
			result.fusionNetGain()==afterFusion.covered()-afterFission.covered(),
			"Measured layer gains must equal the coverage changes in the two published states");
		require(next.sequenceCount()==initial.sequenceCount() && next.filtered()==initial.filtered(),
			"A layer preserves the original input population and length filter across every wave");
		return new LayerResult(result,before,afterFission,afterFusion,next.fusionProductIds(),started);
	}

	private void writeResidentLayerReports(final LayerResult layer){
		writeLiveMetrics(layer.before,layer.afterFission,layer.afterFusion);
	}

	private void writeLiveState(final LayerResult layer){
		final LiveProbeLayer.Result result=layer.result;
		final LiveProbeState next=result.state();
		// Written last: downstream continuation carries this cursor unchanged through
		// seed-only refinement, even if fusion retired the largest published family ID.
		final ByteStreamWriter writer=open(stateOut);
		try{
			writer.println("key\tvalue");
			writer.println("schema\tlive_layer_state_v1");
			writer.println("input_records\t"+next.sequenceCount());
			writer.println("models\t"+next.models().length);
			writer.println("active_models\t"+next.active().length);
			writer.println("assigned\t"+next.assigned());
			writer.println("residual\t"+next.residualOrdinals().length);
			writer.println("filtered\t"+next.filtered());
			writer.println("next_family_id\t"+next.firstNewId());
			writer.println("family_budget\t"+familyBudget);
			writer.println("price_mode\t"+result.priceMode().name().toLowerCase(Locale.ROOT));
			writer.println("below_mode\t"+result.belowMode().name().toLowerCase(Locale.ROOT));
			writer.println("size_lambda\t"+result.sizeLambda());
			writer.println("setup_slot_cost\t"+(result.setupSlotCost()<0 ? "NA" : Integer.toString(result.setupSlotCost())));
			writer.println("budget_slot_cost\t"+result.initialBudgetSlotCost());
			writer.println("edge_fraction\t"+result.edgeFraction());
			writer.println("priced_fraction\t"+result.initialPricedFraction());
			writer.println("fission_top_fraction\t"+fissionTopFraction);
			writer.println("fusion_bottom_fraction\t"+fusionBottomFraction);
			writer.println("priced_rank\t"+result.initialPricedRank());
			writer.println("initial_positive_families\t"+result.initialPositiveCount());
			writer.println(LiveProbeLayer.appendSlotCost(new ByteBuilder().append("initial_slot_cost\t"),result.initialSlotCost()));
			writer.println("slot_fraction\t"+result.slotFraction());
			writer.println("fission_waves\t"+result.fissionWaves());
			writer.println("accepted_fissions\t"+result.acceptedFissions());
			writer.println("accepted_fusions\t"+result.acceptedFusions());
			writer.println("fusion_product_families\t"+layer.fusionProductIds.length);
			writer.println("engine_elapsed_ns\t"+result.elapsedNanos());
			layer.elapsedNanos=System.nanoTime()-layer.started;
			writer.println("total_elapsed_ns\t"+layer.elapsedNanos);
		}finally{close(writer,stateOut);}
	}

	private void printLiveLayer(final LayerResult layer,final LiveProbeState initial,final LiveProbeState next){
		final LiveProbeLayer.Result result=layer.result;
		outstream.println("Live layer: waves="+result.fissionWaves()+", fissions="+result.acceptedFissions()+
			", fusions="+result.acceptedFusions()+", initial slot cost="+result.initialSlotCost()+
			", covered="+initial.assigned()+" -> "+next.assigned()+", models="+next.models().length+
			", next_family_id="+next.firstNewId());
	}

	/** Measures one frozen mixed batch; selection can use a fitted prediction or exact measured gains. */
	private void liveCombined(final LiveProbeState initial,final ClusteringInput input,final FamilySketchBuilder builder,
			final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,final CandidateFilter filter)throws Exception{
		final LiveProbeMetrics metrics=new LiveProbeMetrics(input.sequences,memberGate,workerThreads,familyBudget);
		final LiveProbeMetrics.Snapshot before=measureLiveMetrics(metrics,initial.models(),"before");
		final double splitThreshold=Math.max(0,slotFraction*before.slotCost());
		final double fusionThreshold=splitThreshold==0 ? 0 : -splitThreshold;
		require(Double.isFinite(splitThreshold),"slotfraction times current rank-budget size must remain finite");
		outstream.println("Combined frozen probes starting; actual slot cost="+before.slotCost()+
			", selection="+(rankByPrediction ? "predicted" : "measured")+
			", fission threshold="+splitThreshold+", fusion threshold="+fusionThreshold);
		final LiveProbeClustering.CombinedResult result=LiveProbeClustering.runCombined(initial,input.sequences,builder,memberGate,mergeGate,filter,localMinLength,
				structuralObjective(),structuralOptions(),mergeMaxSizeSet ? mergeMaxSize : Integer.MAX_VALUE,featuresReport,overwrite,
				rankByPrediction,splitThreshold,fusionThreshold,wideOutput);
		final LiveProbeState next=result.state();
		final LiveProbeMetrics.Snapshot after=measureLiveMetrics(metrics,next.models(),"after_combined");
		require(result.slotCost()==before.slotCost() && result.topMembersBefore()==before.topMembers() &&
			result.topMembersAfter()==after.topMembers(),
			"Combined selection and independent published-state census must agree on diagnostic top-budget coverage");
		require(result.netGain()==after.covered()-before.covered() && next.filtered()==initial.filtered(),
			"Combined measured gains must equal the published coverage change without changing the length filter");
		keyInputs.verify();
		result.writeReport(report,overwrite);result.writeFit(fitReport,overwrite);
		writeLiveMetrics(new LiveProbeMetrics.Snapshot[]{before,after},new String[]{"before","after_combined"});
		final ByteStreamWriter writer=open(passesReport);
		try{
			writer.println("phase\tcandidates\tfission_candidates\tfusion_candidates\tmeasured\tfitted\taccepted\tslot_cost"+
				"\tfission_threshold\tfusion_threshold\tnet_gain\telapsed_ns\tcovered_before\tcovered_after\ttop_members_before\ttop_members_after"+
				"\tfission_measured\tfusion_measured\tfission_fitted\tfusion_fitted\tfission_accepted\tfusion_accepted\ttop_budget_decreased");
			final ByteBuilder row=new ByteBuilder().append("combined").tab().append(result.candidates()).tab()
				.append(result.fissionCandidates()).tab().append(result.fusionCandidates()).tab().append(result.measured()).tab()
				.append(result.fitted()).tab().append(result.accepted()).tab().append(result.slotCost()).tab();
			row.appendSlow(result.fissionThreshold()).tab().appendSlow(result.fusionThreshold());
			row.tab().append(result.netGain()).tab().append(result.elapsedNanos()).tab()
				.append(before.covered()).tab().append(after.covered()).tab().append(before.topMembers()).tab().append(after.topMembers()).tab()
				.append(result.fissionMeasured()).tab().append(result.fusionMeasured()).tab().append(result.fissionFitted()).tab()
				.append(result.fusionFitted()).tab().append(result.fissionAccepted()).tab().append(result.fusionAccepted()).tab()
				.append(result.topBudgetDecreased());
			writer.println(row);
		}finally{close(writer,passesReport);}
		writeLiveMembers(next,input);
		final ClusterSketchCodec.Header h=next.sourceHeader();
		ClusterSketchCodec.write(out,overwrite,next.models(),h);
		outstream.println("Combined probes: requested="+result.candidates()+", measured="+result.measured()+", fitted="+result.fitted()+
			", accepted="+result.accepted()+", covered="+before.covered()+" -> "+after.covered()+
			", top-"+familyBudget+" members="+before.topMembers()+" -> "+after.topMembers());
	}

	/** Reuses loaded inputs and publishes the actual fission state before proposing any fusion. */
	private void liveCycle(final LiveProbeState initial,final ClusteringInput input,final FamilySketchBuilder builder,
			final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,final CandidateFilter filter)throws Exception{
		final LiveProbeMetrics metrics=new LiveProbeMetrics(input.sequences,memberGate,workerThreads,familyBudget);
		final LiveProbeMetrics.Snapshot before=measureLiveMetrics(metrics,initial.models(),"before");
		final LayerCakeClustering.Options options=structuralOptions();
		final LocalFamilyObjective objective=structuralObjective();
		final int maximum=mergeMaxSizeSet ? mergeMaxSize : Integer.MAX_VALUE;
		outstream.println("Live cycle fission starting");
		final LiveProbeClustering.Result fission=LiveProbeClustering.runWithFraction(initial,input.sequences,builder,memberGate,mergeGate,
			filter,localMinLength,objective,options,LayerCakeSchedule.Phase.FISSION,slotFraction,maximum);
		livePhaseCompleted("fission",fission);
		final LiveProbeMetrics.Snapshot afterFission=measureLiveMetrics(metrics,fission.state().models(),"after_fission");
		outstream.println("Live cycle fusion starting");
		final LiveProbeClustering.Result fusion=LiveProbeClustering.runWithFraction(fission.state(),input.sequences,builder,memberGate,mergeGate,
			filter,localMinLength,objective,options,LayerCakeSchedule.Phase.FUSION,slotFraction,maximum);
		livePhaseCompleted("fusion",fusion);
		final LiveProbeMetrics.Snapshot afterFusion=measureLiveMetrics(metrics,fusion.state().models(),"after_fusion");
		final LiveProbeState next=fusion.state();
		assert(next.sequenceCount()==initial.sequenceCount() && next.filtered()==initial.filtered()) :
			"A cycle reuses the same original input population and length filter for both phases";
		keyInputs.verify();
		writeLiveMembers(next,input);
		final ClusterSketchCodec.Header h=next.sourceHeader();
		ClusterSketchCodec.write(out,overwrite,next.models(),h);
		fission.writeReport(report,overwrite);fusion.writeReport(fusionReport,overwrite);
		writeLiveMetrics(before,afterFission,afterFusion);
		writeLivePasses(initial,fission,fusion);
		outstream.println("Live cycle: fission accepted="+fission.accepted()+", fusion accepted="+fusion.accepted()+
			", covered="+initial.assigned()+" -> "+next.assigned()+", residual="+next.residualOrdinals().length+
			", next_family_id="+next.firstNewId());
	}

	private LiveProbeMetrics.Snapshot measureLiveMetrics(final LiveProbeMetrics metrics,final Family[] models,final String stage){
		outstream.println("Live cycle metrics "+stage+" starting; models="+models.length);
		final LiveProbeMetrics.Snapshot measured=metrics.measure(models);
		outstream.println("Live cycle metrics "+stage+" completed; seconds="+(measured.elapsedNanos()*1e-9)+
			", covered="+measured.covered()+", rescored_families="+measured.measuredFamilies()+", reused_families="+measured.reusedFamilies());
		return measured;
	}

	private void livePhaseCompleted(final String phase,final LiveProbeClustering.Result result){
		outstream.println("Live cycle "+phase+" completed; seconds="+(result.elapsedNanos()*1e-9)+
			", candidates="+result.candidates()+", probed="+result.available()+", accepted="+result.accepted()+", net_gain="+result.netGain());
	}

	private void writeLivePasses(final LiveProbeState initial,final LiveProbeClustering.Result fission,
			final LiveProbeClustering.Result fusion){
		final ByteStreamWriter writer=open(passesReport);
		try{
			writer.println("phase\tcandidates\tprobed\tunavailable\taccepted\tslot_cost\tthreshold\tnet_gain\telapsed_ns\tcovered_before\tcovered_after\tresidual_before\tresidual_after");
			livePassRow(writer,"fission",initial,fission);
			livePassRow(writer,"fusion",fission.state(),fusion);
		}finally{close(writer,passesReport);}
	}

	/** Phase time excludes the separately measured before/after consensus scoring. */
	private static void livePassRow(final ByteStreamWriter writer,final String phase,final LiveProbeState before,
			final LiveProbeClustering.Result result){
		final LiveProbeState after=result.state();
		assert(result.netGain()==(long)after.assigned()-before.assigned()) : "Phase counters must match the published coverage delta";
		final ByteBuilder row=new ByteBuilder();
		row.append(phase).tab().append(result.candidates()).tab().append(result.available()).tab()
			.append(result.candidates()-result.available()).tab().append(result.accepted()).tab();
		LiveProbeLayer.appendSlotCost(row,result.slotCost()).tab()
			.appendSlow(result.threshold()).tab().append(result.netGain()).tab().append(result.elapsedNanos()).tab()
			.append(before.assigned()).tab().append(after.assigned()).tab().append(before.residualOrdinals().length).tab()
			.append(after.residualOrdinals().length);
		writer.println(row);
	}

	private void writeLiveMetrics(final LiveProbeMetrics.Snapshot before,final LiveProbeMetrics.Snapshot afterFission,
			final LiveProbeMetrics.Snapshot afterFusion){
		writeLiveMetrics(new LiveProbeMetrics.Snapshot[]{before,afterFission,afterFusion},new String[]{"before","after_fission","after_fusion"});
	}

	/** Labels actual committed states; a mixed pass has no intermediate single-operation state. */
	private void writeLiveMetrics(final LiveProbeMetrics.Snapshot[] snapshots,final String[] stages){
		assert(snapshots.length==stages.length) : "Each census belongs to exactly one published state";
		ByteStreamWriter summary=null,histogram=null;
		try{
			summary=open(liveMetrics);histogram=open(liveSizeHist);
			summary.println(LiveProbeMetrics.SUMMARY_HEADER);histogram.println(LiveProbeMetrics.HISTOGRAM_HEADER);
			for(int i=0;i<snapshots.length;i++){
				snapshots[i].writeSummary(summary,stages[i]);snapshots[i].writeHistogram(histogram,stages[i]);
			}
		}finally{close(summary,liveMetrics,histogram,liveSizeHist);}
	}

	/** Validate the complete published partition before opening a membership writer. No payload is recomputed. */
	private void writeLiveMembers(final LiveProbeState state,final ClusteringInput input){
		require(state.sequenceCount()==input.sequences.length,"Live state and membership output must share original ordinals");
		final int[] owners=new int[input.sequences.length]; Arrays.fill(owners,-2);
		for(int ordinal:state.residualOrdinals()){
			require(ordinal>=0 && ordinal<owners.length && owners[ordinal]==-2,
				"Live residual output requires unique valid original ordinals");
			owners[ordinal]=-1;
		}
		int assigned=0,filtered=0;
		for(Family family:state.models()){
			for(int j=0;j<family.size();j++){
				final int ordinal=family.memberIndex(j);
				require(ordinal>=0 && ordinal<owners.length && owners[ordinal]==-2,
					"Live members and residuals must partition original ordinals without overlap");
				owners[ordinal]=family.familyId;assigned++;
			}
		}
		for(int i=0;i<owners.length;i++){
			final boolean shortInput=input.sequences[i].length()<localMinLength;
			require((owners[i]==-2)==shortInput,"Only length-filtered original inputs may be absent from live coverage and residuals");
			if(shortInput){filtered++;}
		}
		require(assigned==state.assigned() && filtered==state.filtered(),"Live output census must match published state");
		final ByteStreamWriter writer=open(membersOut);
		try{
			writer.println("input_index\tinput_id\tfamily_id\tstatus");
			final ByteBuilder row=new ByteBuilder();
			for(int i=0;i<owners.length;i++){
				final int owner=owners[i];
				row.clear().append(i).tab().append(input.sequences[i].id).tab().append(owner>=0 ? owner : -1).tab()
					.append(owner>=0 ? "assigned" : owner==-1 ? "unassigned" : "filtered_length");
				writer.println(row);
			}
		}finally{close(writer,membersOut);}
	}

	/** Generates a bounded uniform pilot from the same unchanged selected baseline.
	 * The full size-filtered neighbor table is built once; no residual index or
	 * response scoring is needed to choose requests. Scientific run controllers
	 * retain the input/configuration pins alongside the emitted sampling metadata.
	 */
	private void localPlan() throws Exception{
		require(in!=null && out!=null && codecs!=null && selectedSeeds!=null && baselineMembers!=null && report!=null,
			"mode=localplan requires in, out, codecs, selectedseeds, baselinemembers and report");
		require(stageACodecs==null && stageAMembers==null && stageBCodec==null && stageBMembers==null &&
			partitionMap==null && roster==null && passesReport==null && membersOut==null,
			"Local request generation does not consume legacy Stage-A/union inputs or member outputs");
		final String[] paths=split(codecs,"codecs");
		require(paths.length==1,"Local request generation requires the single unchanged original codec");
		PathRegistry.prepare().inputs(in,paths[0],selectedSeeds,baselineMembers,blacklist,matrixFile).outputs(out,report);
		require(Tools.testOutputFiles(overwrite,false,false,out,report),"Local request/metadata outputs are unavailable");
		final ClusterSketchCodec.Loaded source=ClusterSketchCodec.read(paths[0]);
		final ClusterSketchCodec.Header h=source.header;
		final FamilySketchBuilder builder=familySketchBuilder(h,loadBlacklist(h,out,report),out,report);
		final ClusteringInput input=ClusteringInput.load(in,-1,false);
		final MemberRankInputs memberInputs=memberRanking ? new MemberRankInputs(input.sequences,memberHydroWeight) : null;
		final LocalProbeBaseline baseline=memberRanking ?
			LocalProbeBaseline.load(source,selectedSeeds,baselineMembers,input.sequences,localMinLength,selectedCount,memberInputs) :
			LocalProbeBaseline.load(source,selectedSeeds,baselineMembers,input.sequences,localMinLength,selectedCount);
		final TailAlignmentGate geometry=gate();
		final TailAlignmentGate mergeGate=new TailAlignmentGate(geometry.mode,denominator,maxAlignLength,minSpan,mergeIdentity,mergeRatio,null,geometry.matrix,minPaired);
		final LayerCakeClustering.Options options=structuralOptions();
		final int[] counts=new int[baseline.selected.length];
		for(int i=0;i<counts.length;i++){counts[i]=baseline.selected[i].size();}
		final FamilyNeighborTable neighbors=memberRanking ?
			FamilyNeighborTable.buildMemberRanked(baseline.selected,counts,builder,mergeGate,
				mergeShortlist,fusionSizeRatio,frozenIndexThreads,workerThreads,memberInputs) :
			FamilyNeighborTable.build(baseline.selected,counts,builder,mergeGate,
				mergeShortlist,fusionSizeRatio,frozenIndexThreads,workerThreads);
		final LocalProbeRequests.Selection selection=LocalProbeRequests.sample(baseline,neighbors,options,mergeGate,
			fissionProbes,fusionProbes,probeSeed,probeStratified);
		keyInputs.verify();
		selection.writeRequests(out,overwrite);
		final ByteStreamWriter writer=open(report);
		try{
			writer.println("key\tvalue");
			if(memberRanking){
				writer.println("member_ranking\ttrue");
				writer.println("memberhydroweight\t"+memberHydroWeight);
				writer.println("member_ranking_layer\t"+baseline.memberRankingLayer);
			}
			writer.println("sampling_method\t"+selection.samplingMethod());
			writer.println("scope\t"+(probeStratified ? "five_size_rank_bands_within_operation" : "uniform_within_operation_pilot"));
			writer.println("master_seed\t"+selection.seed);
			writer.println("fission_seed\t"+selection.fissionSeed());
			writer.println("fusion_seed\t"+selection.fusionSeed());
			writer.println("selected_models\t"+baseline.selected.length);
			writer.println("positive_models\t"+baseline.active.length);
			writer.println("residual_tails\t"+baseline.residualOrdinals.length);
			writer.println("fission_requests\t"+selection.fissionCount);
			writer.println("fusion_requests\t"+selection.fusionCount);
			writer.println("eligible_fissions\t"+selection.fissionPoolSize);
			writer.println("reciprocal_fusions\t"+selection.reciprocalFusionPoolSize);
			writer.println("eligible_fusions\t"+selection.fusionPoolSize);
			writer.println("fission_floor\t"+selection.fissionFloor);
			writer.println("fusion_size_ratio\t"+selection.fusionSizeRatio);
			writer.println("merge_shortlist\t"+selection.mergeShortlist);
			writer.println("merge_min_identity\t"+selection.mergeMinIdentity);
			writer.println("merge_min_ratio\t"+selection.mergeMinScoreRatio);
			writer.println("neighbor_alignments\t"+neighbors.attemptedAlignments);
			writer.println("neighbor_primitive_bytes\t"+neighbors.primitiveBytes());
			if(selection.stratified()){
				for(int op=0;op<2;op++){
					final boolean split=op==0;final String prefix=split ? "fission_band_" : "fusion_band_";
					for(int band=0;band<selection.strata();band++){
						final String key=prefix+band;
						writer.println(key+"_population\t"+selection.stratumPopulation(split,band));
						writer.println(key+"_quota\t"+selection.stratumQuota(split,band));
						writer.println(key+"_min_size\t"+selection.stratumMinSize(split,band));
						writer.println(key+"_max_size\t"+selection.stratumMaxSize(split,band));
						writer.println(key+"_seed\t"+selection.stratumSeed(split,band));
					}
				}
			}
		}finally{close(writer,report);}
		outstream.println("Wrote "+selection.size()+(probeStratified ? " size-stratified" : " uniform")+
			" pilot requests without replacement; later infeasible probes retain their rows.");
	}

	/** One immutable setup serves every persisted local request. The original
	 * codec supplies unchanged payloads; the selected roster and current assignment
	 * table supply identities/memberships. No Stage-A partition map is inferred.
	 */
	private void localProbe() throws Exception{
		final long startup=System.nanoTime();
		require(in!=null && out!=null && codecs!=null && selectedSeeds!=null && baselineMembers!=null,
			"mode=localprobe requires in, out, codecs, selectedseeds, baselinemembers and probes");
		require(stageACodecs==null && stageAMembers==null && stageBCodec==null && stageBMembers==null &&
			partitionMap==null && roster==null && report==null && passesReport==null && membersOut==null,
			"Local responses use one selected baseline and one out table; legacy Stage-A/union/output flags are unsupported");
		final String[] paths=split(codecs,"codecs");
		require(paths.length==1,"Local probes require the single original model codec used by selected-baseline setup");
		final PathRegistry registry=PathRegistry.prepare().inputs(in,paths[0],selectedSeeds,baselineMembers,probePlan,blacklist,matrixFile);
		registry.outputs(out);
		if(rebuildReport!=null){registry.outputs(rebuildReport);}
		if(auditReport!=null){registry.outputs(auditReport);}
		if(featuresReport!=null){registry.outputs(featuresReport);}
		if(refinementReport!=null){registry.outputs(refinementReport);}
		require(Tools.testOutputFiles(overwrite,false,false,out,rebuildReport,auditReport,featuresReport,refinementReport),"Local response output is unavailable");
		final StructuralProbeBatch.Plan requests=StructuralProbeBatch.loadRequests(probePlan);
		final ClusterSketchCodec.Loaded source=ClusterSketchCodec.read(paths[0]);
		final ClusterSketchCodec.Header h=source.header;
		final FamilySketchBuilder builder=familySketchBuilder(h,loadBlacklist(h,out,rebuildReport,auditReport,featuresReport,refinementReport),
			out,rebuildReport,auditReport,featuresReport,refinementReport);
		final CandidateFilter filter=new CandidateFilter(localMinimumHits,0,localHydroDelta,localHydroDelta,localFlex);
		filter.validateWidth(h.width);
		final ClusteringInput input=ClusteringInput.load(in,-1,false);
		final MemberRankInputs memberInputs=memberRanking ? new MemberRankInputs(input.sequences,memberHydroWeight) : null;
		final LocalProbeBaseline baseline=memberRanking ?
			LocalProbeBaseline.load(source,selectedSeeds,baselineMembers,input.sequences,localMinLength,selectedCount,memberInputs) :
			LocalProbeBaseline.load(source,selectedSeeds,baselineMembers,input.sequences,localMinLength,selectedCount);
		LocalProbeBatch.validate(requests,baseline);
		final TailAlignmentGate memberGate=gate();
		final TailAlignmentGate mergeGate=new TailAlignmentGate(memberGate.mode,denominator,maxAlignLength,minSpan,mergeIdentity,mergeRatio,null,memberGate.matrix,minPaired);
		final TailPoolIndex pool=LocalProbePool.build(baseline,input.sequences,builder,workerThreads);
		final LocalProbeEvaluator evaluator=new LocalProbeEvaluator(input.sequences,pool,builder,memberGate,filter,localMinLength);
		final LayerCakeClustering engine=LayerCakeClustering.localProbes(baseline.active,baseline.selected,input.sequences,
			builder,memberGate,mergeGate,structuralObjective(),structuralOptions(),baseline.firstNewId);
		if(memberInputs!=null){engine.withMemberRanking(memberInputs,baseline.memberRankingLayer);}
		outstream.println("Local probe input loading and shared pool/evidence indexes seconds: "+((System.nanoTime()-startup)*1e-9));
		LocalProbeBatch.run(requests,baseline,engine,evaluator,workerThreads,temperature,out,rebuildReport,auditReport,featuresReport,refinementReport,overwrite);
		keyInputs.verify();
		outstream.println("Completed "+requests.size()+" independent local responses against "+baseline.selected.length+
			" selected models and "+pool.size()+" residual tails; affected members were rechecked and the baseline was unchanged.");
	}

	/** Repeated disjoint reciprocal pairs, using actual original members throughout.
	 * Unpaired models are unchanged. A rejected consensus proposal keeps the smaller-ID
	 * provisional model and the complete union, matching the existing union fallback.
	 */
	private void mutual() throws Exception{
		require(in!=null && out!=null && roster!=null && partitionMap!=null && stageACodecs!=null &&
			stageAMembers!=null && report!=null && passesReport!=null && membersOut!=null,
			"mode="+mode+" requires in, out, roster, partitionmap, stageacodecs, stageamembers, report, passesreport and members");
		require(stageBCodec==null && stageBMembers==null,
			"Mutual merging never consumes transitive Stage-B consensus clusters");
		final String[] paths=split(stageACodecs, "stageacodecs"), memberPaths=split(stageAMembers, "stageamembers");
		require(paths.length==memberPaths.length, "Stage-A codecs and membership lists must have equal lengths");
		final PathRegistry registry=PathRegistry.prepare().inputs(paths).inputs(memberPaths).inputs(in,partitionMap,roster,blacklist,probePlan,matrixFile);
		final StructuralProbeBatch.Plan probes="probe".equals(mode) ? StructuralProbeBatch.load(probePlan,out) : null;
		final String[] modelOutputs=probes==null ? new String[]{out} : probes.outputs();
		final String[] outputs=Arrays.copyOf(modelOutputs,modelOutputs.length+3);
		outputs[modelOutputs.length]=report;outputs[modelOutputs.length+1]=passesReport;outputs[modelOutputs.length+2]=membersOut;
		require(Tools.testOutputFiles(overwrite,false,false,outputs),"Structural outputs are unavailable");
		registry.outputs(outputs);
		final CodecSet loaded=loadCodecs(paths);
		final Roster source=Roster.load(roster, loaded);
		final ClusteringInput input=ClusteringInput.load(in, -1, false);
		final int[][] parentByLocal=loadParentMap(partitionMap, paths.length, input);
		final IntList[] originals=loadStageAMembers(memberPaths, parentByLocal, input, source);
		final MemberRankInputs memberInputs=memberRanking ? new MemberRankInputs(input.sequences,memberHydroWeight) : null;
		int measurementLayer=memberRanking ? initialMutualRankingLayer(loaded) : -1;
		final ClusterSketchCodec.Header header=loaded.header;
		final DDLPentamerBlacklist blocked=loadBlacklist(header,outputs);
		final FamilySketchBuilder builder=familySketchBuilder(header,blocked,outputs);
		Family[] families=new Family[source.rows.length];
		int ordinal=0;
		long originalCount=0;
		for(ClusterSketchCodec.Loaded file : loaded.files){
			for(ClusterSketchCodec.Entry entry : file.families){
				final IntList members=originals[ordinal]; members.sort();
				final Family family=new Family(ordinal, source.rows[ordinal].id, members.array[0], entry.consensus);
				for(int j=1; j<members.size; j++){family.addMember(members.array[j]);}
				family.setSketch(entry.sketch); family.setMonomers(entry.monomers);
				family.setScalars(entry.scalars, entry.scalarMask);
				require(family.size()==entry.members, "Original memberships must agree before reciprocal selection");
				if(memberInputs!=null){family.enableMemberRanking(memberInputs,measurementLayer,null);}
				families[ordinal++]=family; originalCount+=family.size();
			}
		}
		final TailAlignmentGate memberGate=gate();
		final TailAlignmentGate mergeGate=new TailAlignmentGate(memberGate.mode, denominator,
			maxAlignLength, minSpan, mergeIdentity, mergeRatio,null,memberGate.matrix,minPaired);
		if(probes!=null){
			probes.validate(families.length);
			final LayerCakeClustering engine=new LayerCakeClustering(families,input.sequences,builder,memberGate,mergeGate,
				structuralObjective(),structuralOptions(),null);
			if(memberInputs!=null){engine.withMemberRanking(memberInputs,measurementLayer);}
			StructuralProbeBatch.run(probes,engine,families,input.sequences,header,temperature,report,passesReport,membersOut,overwrite);
			keyInputs.verify();
			outstream.println("Produced "+probes.size()+" independent probe rows; rejected rows have empty codecs and explicit status. No reassignment was run.");
			return;
		}
		if("layercake".equals(mode)){
			structural(families,input,builder,memberGate,mergeGate,header,originalCount,memberInputs,measurementLayer);
			return;
		}
		ByteStreamWriter pairs=null, passes=null;
		try{
			pairs=open(report); passes=open(passesReport);
			pairs.println("pass\tleft_family_id\tright_family_id\tleft_members\tright_members\tmerged_family_id\tmerged_members\tproposals\taccepted\trejected\tunchanged\tincompatible_snapshots\tconstruction_alignments\tleft_identity\tleft_ratio\tright_identity\tright_ratio");
			passes.println("pass\tfamilies_before\teligible_queries\tattempted_alignments\tspan_valid_alignments\treciprocal_pairs\tgate_rejected_pairs\tmerged_pairs\tfamilies_after\toriginal_members\telapsed_ns");
			final ByteBuilder row=new ByteBuilder();
			for(int pass=1; ; pass++){
				final long started=System.nanoTime();
				final int[] weights=new int[families.length];
				for(int i=0; i<weights.length; i++){weights[i]=families[i].size();}
				final MutualNearestMerge.Result selected=memberRanking ?
					MutualNearestMerge.selectMemberRanked(families,weights,builder,mergeGate,
						mergeMaxSize,mergeShortlist,workerThreads,frozenIndexThreads,memberInputs) :
					MutualNearestMerge.select(families,weights,builder,mergeGate,
						mergeMaxSize,mergeShortlist,workerThreads,frozenIndexThreads);
				if(memberRanking && selected.pairs.length>0){
					require(measurementLayer<Integer.MAX_VALUE,"Mutual union cannot advance beyond the maximum member-ranking layer");
					measurementLayer=Math.addExact(measurementLayer,1);
				}
				final MutualRebuild rebuild=new MutualRebuild(families,selected.pairs,input,builder,memberGate,
					memberInputs,measurementLayer);
				rebuild.run();
				final Family[] next=new Family[families.length-selected.pairs.length];
				final Family[] replacements=new Family[families.length];
				final boolean[] removed=new boolean[families.length];
				for(int i=0; i<selected.pairs.length; i++){
					final MutualNearestMerge.Pair pair=selected.pairs[i];
					final int left=pair.leftIndex, right=pair.rightIndex;
					require(left!=right && !removed[left] && !removed[right] && replacements[left]==null && replacements[right]==null,
						"Frozen reciprocal pairs must be disjoint before membership publication");
					final Family merged=rebuild.merged[i];
					replacements[left]=merged; removed[right]=true;
					final ConsensusRefiner.Update update=rebuild.updates[i];
					row.clear().append(pass).tab().append(families[left].familyId).tab().append(families[right].familyId)
						.tab().append(weights[left]).tab().append(weights[right]).tab().append(merged.familyId)
						.tab().append(merged.size()).tab().append(update.proposals()).tab().append(update.accepted())
						.tab().append(update.rejected()).tab().append(update.unchanged()).tab()
						.append(update.incompatibleSnapshots()).tab().append(update.alignments())
						.tab().appendSlow(selected.nearestIdentityFor(left)).tab().appendSlow(selected.nearestRatioFor(left))
						.tab().appendSlow(selected.nearestIdentityFor(right)).tab().appendSlow(selected.nearestRatioFor(right));
					pairs.println(row);
				}
				int kept=0; long members=0;
				for(int i=0; i<families.length; i++){
					if(removed[i]){continue;}
					final Family family=replacements[i]==null ? families[i] : replacements[i];
					next[kept++]=family; members+=family.size();
				}
				require(kept==next.length && members==originalCount,
					"Reciprocal unions must conserve every original member and remove exactly one family per pair");
				row.clear().append(pass).tab().append(families.length).tab().append(selected.eligibleQueries)
					.tab().append(selected.attemptedAlignments).tab().append(selected.spanValidAlignments)
					.tab().append(selected.reciprocalPairs).tab().append(selected.gateRejectedPairs)
					.tab().append(selected.pairs.length).tab().append(next.length).tab().append(members)
					.tab().append(System.nanoTime()-started);
				passes.println(row);
				outstream.println("Mutual merge pass "+pass+": "+families.length+" -> "+next.length+
					" families; "+selected.pairs.length+" pairs; "+members+" original members");
				families=next;
				if(selected.pairs.length==0){break;}
				assert(pass<source.rows.length) : "Each nonterminal reciprocal pass removes at least one source family";
			}
		}finally{close(pairs, report, passes, passesReport);}
		keyInputs.verify();
		writeMutualMembers(families, input, originalCount);
		ClusterSketchCodec.write(out, overwrite, families, header);
	}

	/** Fresh measurements must follow imported provenance; legacy sources have no prior layer. */
	private static int initialMutualRankingLayer(final CodecSet source){
		assert(source!=null && source.files!=null) : "Mutual descriptor hydration follows validated Stage-A codec loading";
		int layer=0;
		for(ClusterSketchCodec.Loaded file:source.files){
			for(ClusterSketchCodec.Entry entry:file.families){
				if(entry.memberRanking!=null){
					require(entry.memberRanking.layer<Integer.MAX_VALUE,
						"Mutual hydration cannot advance beyond the maximum member-ranking layer");
					layer=Math.max(layer,Math.addExact(entry.memberRanking.layer,1));
				}
			}
		}
		return layer;
	}

	/** Reuses the validated original-member loader; NONE exports unchanged globally numbered seeds. */
	private void structural(final Family[] initial,final ClusteringInput input,final FamilySketchBuilder builder,
			final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,final ClusterSketchCodec.Header header,
			final long originalCount,final MemberRankInputs memberInputs,final int measurementLayer)throws Exception{
		ByteStreamWriter decisions=null,passes=null;
		final LayerCakeSchedule.Result result;
		try{
			decisions=open(report); passes=open(passesReport);
			final ByteStreamWriter phaseWriter=passes;
			passes.println("cycle\tphase\ttemperature\taccepted_operations\tfamilies_before\tfamilies_after\toriginal_members\ttop_budget_members\tslot_cost\telapsed_ns\tstop");
			final LayerCakeClustering engine=new LayerCakeClustering(initial,input.sequences,builder,memberGate,mergeGate,
				structuralObjective(),structuralOptions(),decisions);
			if(memberInputs!=null){engine.withMemberRanking(memberInputs,measurementLayer);}
			final ByteBuilder row=new ByteBuilder();
			final LayerCakeSchedule.Summary start=LayerCakeSchedule.summarize(initial,familyBudget,input.sequences.length);
			structuralRow(row,0,"initial",temperature,0,start,start,0,"CONTINUE"); passes.println(row);
			result=new LayerCakeSchedule(structuralSettings(input.sequences.length)).run(initial,engine,
				new LayerCakeSchedule.Listener(){
					@Override public void completed(final LayerCakeSchedule.PhaseSummary phase){
						structuralRow(row,phase.cycle,phase.phase.name(),phase.temperature,phase.acceptedOperations,
							phase.before,phase.after,phase.elapsedNanos,"CONTINUE"); phaseWriter.println(row);
						outstream.println("Structural "+phase.phase+" cycle "+phase.cycle+": "+phase.before.families+
							" -> "+phase.after.families+" families; "+phase.acceptedOperations+" accepted operations");
					}
				});
			final double lastTemperature=result.history.length==0 ? temperature : result.history[result.history.length-1].temperature;
			structuralRow(row,result.cycles,"final",lastTemperature,0,result.summary,result.summary,0,result.stop.name()); passes.println(row);
			outstream.println("Structural stop: "+result.stop+"; raw top-"+familyBudget+" membership "+result.summary.topMembers+
				" / "+input.sequences.length+" inputs. This is a structural census, not calibrated retention.");
			outstream.println("Structural work: fission proposals="+engine.fissionProposals+", fusion proposals="+engine.fusionProposals+
				", quality alignments="+engine.qualityAlignments+", neighbor alignments="+engine.neighborAlignments+
				", rebuild alignments="+engine.rebuildAlignments+", selection alignments="+engine.selectionAlignments+
				", member postings="+engine.memberPostings+", neighbor postings="+engine.neighborPostings+
				", retained edges="+engine.fissionRetainedEdges+", union edges="+engine.fissionUnionEdges);
		}finally{close(decisions,report,passes,passesReport);}
		keyInputs.verify();
		writeMutualMembers(result.families,input,originalCount);
		ClusterSketchCodec.write(out,overwrite,result.families,header);
	}

	private static void structuralRow(final ByteBuilder row,final int cycle,final String phase,final double temperature,
			final long changes,final LayerCakeSchedule.Summary before,final LayerCakeSchedule.Summary after,
			final long nanos,final String stop){
		row.clear().append(cycle).tab().append(phase).tab().appendSlow(temperature).tab().append(changes)
			.tab().append(before.families).tab().append(after.families).tab().append(after.assignedTails)
			.tab().append(after.topMembers).tab().append(after.slotCost).tab().append(nanos).tab().append(stop);
	}

	private void writeMutualMembers(Family[] families, ClusteringInput input, long expected){
		final int[] owners=new int[input.sequences.length]; Arrays.fill(owners, -1);
		long assigned=0;
		for(Family family : families){
			for(int j=0; j<family.size(); j++){
				final int member=family.memberIndex(j);
				require(member>=0 && member<owners.length && owners[member]<0,
					"Final reciprocal unions must have unique valid original ordinals");
				owners[member]=family.familyId; assigned++;
			}
		}
		require(assigned==expected, "Final reciprocal membership output must conserve original assigned tails");
		final ByteStreamWriter writer=open(membersOut);
		try{
			writer.println("parent_input_index\tinput_id\tfamily_id");
			final ByteBuilder row=new ByteBuilder();
			for(int i=0; i<owners.length; i++){
				if(owners[i]<0){continue;}
				row.clear().append(i).tab().append(input.sequences[i].id).tab().append(owners[i]); writer.println(row);
			}
		}finally{close(writer, membersOut);}
	}

	private final class MutualRebuild implements Accumulator<MutualWorker>{
		MutualRebuild(Family[] families_, MutualNearestMerge.Pair[] pairs_, ClusteringInput input_,
				FamilySketchBuilder builder_, TailAlignmentGate gate_,MemberRankInputs memberInputs_,int measurementLayer_){
			families=families_; pairs=pairs_; input=input_; builder=builder_; gate=gate_;
			memberInputs=memberInputs_;measurementLayer=measurementLayer_;
			require(memberInputs==null ? measurementLayer==-1 : memberInputs.belongsTo(input.sequences) && measurementLayer>=0,
				"Mutual union descriptors require the same original population/cache and an explicit measurement layer");
			merged=new Family[pairs.length]; updates=new ConsensusRefiner.Update[pairs.length];
		}
		void run(){
			if(pairs.length==0){return;}
			final MutualWorker[] workers=new MutualWorker[Math.min(workerThreads, pairs.length)];
			for(int i=0; i<workers.length; i++){workers[i]=new MutualWorker(this);}
			try{
				if(!ThreadWaiter.startAndWait(Arrays.asList(workers), this)){
					if(failure==null){failure=new IllegalStateException("Reciprocal union worker completion failed");}
				}
			}catch(Throwable cause){
				if(failure==null){failure=cause;}
				final ArrayList<MutualWorker> started=new ArrayList<MutualWorker>();
				for(MutualWorker worker : workers){if(worker.getState()!=Thread.State.NEW){started.add(worker);}}
				ThreadWaiter.waitForThreadsToFinish(started);
			}
			if(failure!=null || done!=pairs.length){
				throw new IllegalStateException("Reciprocal union workers did not complete every frozen pair", failure);
			}
		}
		@Override public void accumulate(MutualWorker worker){
			done+=worker.done; if(failure==null){failure=worker.failure;}
		}
		@Override public ReadWriteLock rwlock(){return null;}
		@Override public boolean success(){return failure==null && done==pairs.length;}
		final Family[] families, merged;
		final MutualNearestMerge.Pair[] pairs;
		final ClusteringInput input;
		final FamilySketchBuilder builder;
		final TailAlignmentGate gate;
		final MemberRankInputs memberInputs;
		final int measurementLayer;
		final ConsensusRefiner.Update[] updates;
		final AtomicInteger next=new AtomicInteger();
		int done;
		Throwable failure;
	}

	private final class MutualWorker extends Thread{
		MutualWorker(MutualRebuild parent_){parent=parent_;}
		@Override public void run(){
			try(ConsensusRefiner driver=new ConsensusRefiner(1, consensusTrimLeft, trimDepth, null, true);
					ConsensusRefiner refiner=driver.fork();
					ClusteringMetrics metrics=new ClusteringMetrics(null, null, true, 1, Long.MAX_VALUE)){
				final FamilySketchBuilder builder=parent.builder.fork();
				for(int i=parent.next.getAndIncrement(); i<parent.pairs.length; i=parent.next.getAndIncrement()){
					final MutualNearestMerge.Pair pair=parent.pairs[i];
					final Family a=parent.families[pair.leftIndex], b=parent.families[pair.rightIndex];
					require(a.size()<mergeMaxSize && b.size()<mergeMaxSize,
						"Both frozen pair families must be below the original-member threshold");
					final Family seed=a.familyId<b.familyId ? a : b;
					final IntList members=new IntList(Math.addExact(a.size(), b.size()));
					for(int j=0; j<a.size(); j++){members.add(a.memberIndex(j));}
					for(int j=0; j<b.size(); j++){members.add(b.memberIndex(j));}
					members.sort();
					final Family merged=new Family(seed.familyId, seed.id, members.array[0], seed.consensusUnsafe());
					for(int j=1; j<members.size; j++){
						require(members.array[j]!=members.array[j-1], "Reciprocal pair originals must be disjoint");
						merged.addMember(members.array[j]);
					}
					metrics.begin("consensus", 0, merged.size());
					parent.updates[i]=refiner.updateAfterConcurrentGreedy(merged, parent.input.sequences, parent.gate, metrics);
					metrics.end(merged.size(), merged.size(), 0, 1, 1);
					merged.refreshStats(parent.input.sequences, builder);
					// Mutual unions retain every original, including members of a rejected consensus proposal.
					if(parent.memberInputs!=null){merged.enableMemberRanking(parent.memberInputs,parent.measurementLayer,null);}
					parent.merged[i]=merged; done++;
				}
			}catch(Throwable t){failure=t;}
		}
		final MutualRebuild parent;
		int done;
		Throwable failure;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Collect           ----------------*/
	/*--------------------------------------------------------------*/

	private void collect(){
		require(out!=null && roster!=null, "mode=collect requires out= and roster=");
		require(Tools.testOutputFiles(overwrite, false, false, new String[]{out, roster}), "Stage-B collect outputs are unavailable");
		final String[] paths=split(codecs, "codecs");
		PathRegistry.prepare().inputs(paths).inputs(matrixFile).outputs(out, roster);
		final CodecSet loaded=loadCodecs(paths);
		ByteStreamWriter fasta=null, table=null;
		try{
			fasta=open(out); table=open(roster);
			final ByteBuilder row=new ByteBuilder(), seq=new ByteBuilder();
			table.println(ROSTER_HEADER);
			int ordinal=0;
			for(int shard=0; shard<loaded.files.length; shard++){
				for(ClusterSketchCodec.Entry entry : loaded.files[shard].families){
					require(entry.members>0, "Stage-A codec contains an empty family: "+entry.name);
					final String id=stageAId(shard, entry.familyId);
					row.clear().append(ordinal).tab().append(id).tab().append(shard).tab()
						.append(entry.familyId).tab().append(entry.members);
					table.println(row);
					seq.clear().append('>').append(id).append(" shard=").append(shard)
						.append(" family=").append(entry.familyId).append(" members=").append(entry.members).nl();
					appendConsensus(entry, seq);
					fasta.print(seq);
					ordinal++;
				}
			}
		}finally{
			close(table, roster, fasta, out);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Union            ----------------*/
	/*--------------------------------------------------------------*/

	private void union() throws Exception{
		require(in!=null && out!=null && roster!=null && partitionMap!=null && stageBCodec!=null &&
			stageBMembers!=null && stageACodecs!=null && stageAMembers!=null && report!=null,
			"mode=union requires in=, out=, roster=, partitionmap=, stageacodecs=, stageamembers=, stagebcodec=, stagebmembers= and report=");
		require(Tools.testOutputFiles(overwrite, false, false, new String[]{out, report}), "Stage-B union outputs are unavailable");
		final String[] stageAPaths=split(stageACodecs, "stageacodecs");
		final String[] memberPaths=split(stageAMembers, "stageamembers");
		require(stageAPaths.length==memberPaths.length, "Stage-A codec and member file lists must have the same length");
		PathRegistry.prepare().inputs(stageAPaths).inputs(memberPaths).inputs(in, partitionMap, roster,
			stageBCodec, stageBMembers, blacklist, matrixFile).outputs(out, report);
		final CodecSet stageA=loadCodecs(stageAPaths);
		final ClusterSketchCodec.Loaded stageB=ClusterSketchCodec.read(stageBCodec);
		validateGeometry(stageA.header, stageB.header, stageBCodec);
		final Roster stageARoster=Roster.load(roster, stageA);
		final StageB stageBClusters=StageB.load(stageB);
		final ClusteringInput input=ClusteringInput.load(in, -1, false);
		final int[][] parentByLocal=loadParentMap(partitionMap, stageAPaths.length, input);
		final IntList[] originalsByStageA=loadStageAMembers(memberPaths, parentByLocal, input, stageARoster);
		final int[] ownerByStageA=loadStageBMembers(stageBMembers, stageBClusters, stageARoster);
		final IntList[] stageAByStageB=groupStageAByStageB(ownerByStageA, stageBClusters);
		final DDLPentamerBlacklist blacklist=loadBlacklist(stageA.header,out,report);
		final FamilySketchBuilder builder=familySketchBuilder(stageA.header,blacklist,out,report);
		final UnionResults results=buildUnions(stageB, stageAByStageB,
			originalsByStageA, stageARoster, input, builder);
		final ByteStreamWriter writer=open(report);
		try{
			writer.println(REPORT_HEADER);
			for(int i=0; i<results.unions.length; i++){
				writeReport(writer, stageB.families[i], results.unions[i], results.updates[i], stageA.header);
			}
		}finally{close(writer, report);}
		keyInputs.verify();
		ClusterSketchCodec.write(out, overwrite, results.unions, stageA.header.withFamilies(results.unions.length));
	}

	/** Rebuilds exactly one Stage-B group: stage-A consensuses choose members,
	 * but final union membership is original sample ordinals sorted for determinism.
	 */
	private Family buildUnion(ClusterSketchCodec.Entry seed, IntList stageAOrdinals, IntList[] originalsByStageA,
			Roster stageARoster, ClusteringInput input){
		final Family family=Family.fromSeed(seed);
		final IntList members=new IntList(seed.members*4);
		int sourceMembers=0, consensuses=0;
		for(int i=0; i<stageAOrdinals.size; i++){
			final int ordinal=stageAOrdinals.array[i];
			final IntList originals=originalsByStageA[ordinal];
			members.addAll(originals);
			sourceMembers=Math.addExact(sourceMembers, stageARoster.rows[ordinal].sourceMembers);
			consensuses++;
		}
		require(consensuses==seed.members, "Stage-B membership count disagrees with codec for "+seed.name);
		require(sourceMembers==members.size, "Stage-A source-member counts disagree with assigned rows for "+seed.name);
		members.sort();
		int previous=-1;
		for(int i=0; i<members.size; i++){
			final int member=members.array[i];
			require(member>=0 && member<input.sequences.length && member!=previous,
				"Original-member union contains an invalid or duplicate parent ordinal for "+seed.name);
			family.addMember(member); previous=member;
		}
		require(family.size()==sourceMembers && family.size()>0, "Stage-B union did not retain exact source members for "+seed.name);
		return family;
	}

	private void writeReport(ByteStreamWriter writer, ClusterSketchCodec.Entry seed, Family family,
			ConsensusRefiner.Update update, ClusterSketchCodec.Header header){
		final String status=update.accepted()>0 ? "union_accepted" :
			update.unchanged()>0 ? "union_unchanged" : "union_rejected";
		final String reason=update.accepted()>0 ? "all_members_pass" :
			update.unchanged()>0 ? "same_sequence" :
			update.incompatibleSnapshots()>0 ? "snapshot_no_positive_placement" :
			update.rejected()>0 ? "strict_new_validation_failed" : "not_run";
		final ByteBuilder row=new ByteBuilder();
		row.append(seed.familyId).tab().append(seed.name).tab().append(seed.members).tab()
			.append(family.size()).tab().append(family.size()).tab().append(status).tab().append(reason).tab()
			.append(update.proposals()).tab().append(update.accepted()).tab().append(update.rejected()).tab()
			.append(update.unchanged()).tab().append(update.incompatibleSnapshots()).tab()
			.append(update.alignments()).tab().append(header.k).tab().append(header.width).tab()
			.append(header.sketchSeed).tab().append(header.positional ? 1 : 0).tab()
			.append(header.positionMode.flag).tab().append(header.partitions).tab().append(header.overlap).tab()
			.append(header.blacklistSha80==null ? "." : header.blacklistSha80).tab()
			.append(minSpan).tab().append(maxAlignLength).tab().appendSlow(minIdentity).tab()
			.appendSlow(minRatio).tab().append(denominator.name()).tab().append(pinRight() ? 1 : 0).tab()
			.append(trimLeft ? 1 : 0).tab().append(consensusTrimLeft ? 1 : 0).tab().appendSlow(trimDepth);
		writer.println(row);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Stage-A TSVs        ----------------*/
	/*--------------------------------------------------------------*/

	private IntList[] loadStageAMembers(String[] paths, int[][] parentByLocal, ClusteringInput input,
			Roster roster){
		final IntList[] originals=new IntList[roster.rows.length];
		for(int i=0; i<originals.length; i++){originals[i]=new IntList(roster.rows[i].sourceMembers);}
		for(int shard=0; shard<paths.length; shard++){
			final BitSet seen=new BitSet(parentByLocal[shard].length);
			final ByteFile bf=ByteFile.makeByteFile(paths[shard], false);
			try{
				final LineParser1 parser=new LineParser1('\t');
				requireHeader(bf.nextLine(), "input_index\tinput_id\tfamily_id\tstatus", paths[shard]);
				for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
					parser.set(line);
					require(parser.terms()==4, "Stage-A membership rows need four columns");
					final int local=parser.parseInt(0);
					require(local>=0 && local<parentByLocal[shard].length && !seen.get(local),
						"Duplicate or out-of-range shard-local membership row");
					seen.set(local);
					final int parent=parentByLocal[shard][local];
					require(parent>=0, "Stage-A membership row lacks a partition-map parent");
					require(parser.termEquals(input.sequences[parent].id, 1), "Stage-A input_id disagrees with parent input");
					final int familyId=parser.parseInt(2);
					if(parser.termEquals("assigned", 3)){
						final int ordinal=roster.ordinalFor(shard, familyId);
						require(familyId>=0 && ordinal>=0, "Stage-A assigned row targets an unknown family");
						originals[ordinal].add(parent);
					}else{
						require(familyId==-1 && (parser.termEquals("unassigned", 3) || parser.termEquals("filtered_length", 3)),
							"Stage-A nonassigned rows need family -1 and a known status");
					}
				}
				require(seen.cardinality()==parentByLocal[shard].length, "Stage-A membership file omitted shard-local rows");
			}finally{close(bf, paths[shard]);}
		}
		for(int i=0; i<originals.length; i++){
			require(originals[i].size==roster.rows[i].sourceMembers,
				"Stage-A codec member count disagrees with assigned rows for "+roster.rows[i].id);
		}
		return originals;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Stage-B TSVs        ----------------*/
	/*--------------------------------------------------------------*/

	private int[] loadStageBMembers(String path, StageB stageB, Roster roster){
		final int[] owner=new int[roster.rows.length];
		Arrays.fill(owner, -1);
		final IntHashMap counts=new IntHashMap(stageB.entries.length*2+1);
		final ByteFile bf=ByteFile.makeByteFile(path, false);
		try{
			final LineParser1 parser=new LineParser1('\t');
			requireHeader(bf.nextLine(), "input_index\tinput_id\tfamily_id\tstatus", path);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				parser.set(line);
				require(parser.terms()==4, "Stage-B membership rows need four columns");
				final int ordinal=parser.parseInt(0), familyId=parser.parseInt(2);
				require(ordinal>=0 && ordinal<owner.length && owner[ordinal]<0, "Duplicate or out-of-range Stage-A ordinal in Stage-B members");
				require(parser.termEquals(roster.rows[ordinal].id, 1), "Stage-B input_id disagrees with roster");
				require(parser.termEquals("assigned", 3) && familyId>=0 && stageB.indexFor(familyId)>=0,
					"Every Stage-A consensus must be assigned to a known Stage-B family");
				owner[ordinal]=familyId;
				final int old=counts.get(familyId);
				counts.put(familyId, old<0 ? 1 : old+1);
			}
		}finally{close(bf, path);}
		for(int ordinal=0; ordinal<owner.length; ordinal++){
			require(owner[ordinal]>=0, "Stage-B omitted Stage-A ordinal "+ordinal);
		}
		for(ClusterSketchCodec.Entry entry : stageB.entries){
			require(counts.get(entry.familyId)==entry.members, "Stage-B codec member count disagrees with assigned rows for "+entry.name);
		}
		return owner;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Parent Map          ----------------*/
	/*--------------------------------------------------------------*/

	private int[][] loadParentMap(String path, int shards, ClusteringInput input){
		final int[] maxLocal=new int[shards];
		Arrays.fill(maxLocal, -1);
		final BitSet parents=new BitSet(input.sequences.length);
		int rows=0;
		ByteFile bf=ByteFile.makeByteFile(path, false);
		try{
			final LineParser1 parser=new LineParser1('\t');
			requireHeader(bf.nextLine(), ClusteringShardWriter.MAPPING_HEADER, path);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				parser.set(line);
				require(parser.terms()==4, "Partition mapping rows need four columns");
				final int parent=parser.parseInt(0), shard=parser.parseInt(1), local=parser.parseInt(2);
				require(parent>=0 && parent<input.sequences.length && !parents.get(parent), "Duplicate or out-of-range parent ordinal");
				require(shard>=0 && shard<shards && local>=0 && local<input.sequences.length,
					"Out-of-range shard or local ordinal");
				require(parser.termEquals(input.sequences[parent].id, 3), "Partition input_id disagrees with parent input");
				maxLocal[shard]=Math.max(maxLocal[shard], local);
				parents.set(parent); rows++;
			}
		}finally{close(bf, path);}
		require(rows==input.sequences.length, "Partition mapping must cover every input record exactly once");
		final int[][] parentByLocal=new int[shards][];
		for(int shard=0; shard<shards; shard++){
			require(maxLocal[shard]>=0, "Partition mapping has no rows for shard "+shard);
			parentByLocal[shard]=new int[maxLocal[shard]+1];
			Arrays.fill(parentByLocal[shard], -1);
		}
		bf=ByteFile.makeByteFile(path, false);
		try{
			final LineParser1 parser=new LineParser1('\t');
			requireHeader(bf.nextLine(), ClusteringShardWriter.MAPPING_HEADER, path);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				parser.set(line);
				final int parent=parser.parseInt(0), shard=parser.parseInt(1), local=parser.parseInt(2);
				require(parentByLocal[shard][local]<0, "Duplicate shard-local mapping row");
				parentByLocal[shard][local]=parent;
			}
		}finally{close(bf, path);}
		for(int shard=0; shard<parentByLocal.length; shard++){
			for(int local=0; local<parentByLocal[shard].length; local++){
				require(parentByLocal[shard][local]>=0, "Partition mapping omits a shard-local ordinal");
			}
		}
		return parentByLocal;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Codecs/Roster       ----------------*/
	/*--------------------------------------------------------------*/

	private CodecSet loadCodecs(String[] paths){
		final ClusterSketchCodec.Loaded[] loaded=new ClusterSketchCodec.Loaded[paths.length];
		for(int i=0; i<loaded.length; i++){
			loaded[i]=ClusterSketchCodec.read(paths[i]);
			if(i>0){validateGeometry(loaded[0].header, loaded[i].header, paths[i]);}
		}
		return new CodecSet(loaded, loaded[0].header);
	}

	private DDLPentamerBlacklist loadBlacklist(ClusterSketchCodec.Header header,String...outputs) throws Exception{
		if(blacklistSha80!=null && !blacklistSha80.equals(header.blacklistSha80)){
			throw new IllegalArgumentException("blacklistsha80 does not match the codec header");
		}
		if(blacklistK>=0 || blacklistMode!=null || blacklistProvenance!=null){
			return ddl.ProteinBlacklistPolicy.loadExplicit(header.geometry,blacklist,header.blacklistSha80,blacklistWords,
				blacklistK,blacklistMode,blacklistProvenance,verification==VerificationMode.STRUCTURAL);
		}
		final ddl.ProteinBlacklistPolicy.Descriptor policy=ddl.ProteinBlacklistPolicy.policy(header.geometry);
		if(policy!=null){
			ddl.ProteinBlacklistPolicy.requireStored(header.geometry,header.blacklistSha80);
			if(blacklistWords>=0 && blacklistWords!=policy.words){
				throw new IllegalArgumentException("Continuation requires "+policy.words+" words for "+policy.label);
			}
			if(blacklist==null){
				blacklist=policy.resource().toString();
				// The implicit resource was not available to the initial input/output alias registry.
				final PathRegistry registry=PathRegistry.prepare().inputs(blacklist);
				for(String output:outputs){if(output!=null){registry.outputs(output);}}
			}
			final DDLPentamerBlacklist selected=policy.load(Paths.get(blacklist), verification==VerificationMode.STRUCTURAL);
			return ddl.ProteinBlacklistPolicy.resolve(header.geometry,selected);
		}
		if(header.blacklistSha80==null){
			if(blacklist!=null){throw new IllegalArgumentException("Plain codec headers do not declare a blacklist");}
			if(blacklistWords>=0){throw new IllegalArgumentException("blacklistwords requires a codec with blacklist_sha80");}
			return null;
		}
		if(blacklist==null){throw new IllegalArgumentException("A codec with blacklist_sha80 needs blacklist=");}
		return verification==VerificationMode.STRUCTURAL ?
			DDLPentamerBlacklist.loadStructural(Paths.get(blacklist), header.blacklistSha80, blacklistWords) :
			DDLPentamerBlacklist.load(Paths.get(blacklist), header.blacklistSha80);
	}

	/** Builds a fresh sketcher whose key policy exactly matches a serialized codec header. */
	private FamilySketchBuilder familySketchBuilder(ClusterSketchCodec.Header header,DDLPentamerBlacklist blocked,
			String...outputs)throws Exception{
		final DDLKeyTieBreak.PolicyContext context=keyPolicyContext(header,blocked,outputs);
		return new FamilySketchBuilder(header.width,header.geometry,header.sketchSeed,header.positionMode,
			header.partitions,header.overlap,blocked,context,keyPolicySpecified || context.policy!=DDLKeyTieBreak.Policy.SMALLEST);
	}

	/** Loads standalone raw-tail counts once and binds any serialized nondefault key policy. */
	private DDLKeyTieBreak.PolicyContext keyPolicyContext(ClusterSketchCodec.Header header,
			DDLPentamerBlacklist blocked,String...outputs)throws Exception{
		final DDLKeyTieBreak.Policy policy=effectiveKeyPolicy(header);
		final Path manifest=rawManifest==null ? null : Paths.get(rawManifest);
		protectRawPolicyOutputs(policy,manifest,outputs);
		final DDLKeyTieBreak.PolicyContext context=FamilyKeyPolicyInputs.load(policy,manifest,expectedRawTails,
			rawMinLength,rawMaxLength,rawTailCrop,header.width,header.geometry,header.sketchSeed,header.positionMode,
			header.partitions,header.overlap,blocked,keyInputs,outstream,null,effectiveConsensusBonus(header));
		require(header.sameKeyPolicy(header.withKeyPolicy(context)),
			"Stage-B builder key policy differs from the source codec");
		return context;
	}

	/** Checks Cake-forwarded resident metadata against borrowed counts without recounting. */
	private void validateResidentKeyPolicy(ClusterSketchCodec.Header header,FamilySketchBuilder builder,
			DDLPentamerBlacklist blocked)throws Exception{
		final DDLKeyTieBreak.Policy policy=effectiveKeyPolicy(header);
		final Path manifest=rawManifest==null ? null : Paths.get(rawManifest);
		protectRawPolicyOutputs(policy,manifest,report,fusionReport,passesReport,liveMetrics,liveSizeHist,stateOut,wideFeatures);
		FamilyKeyPolicyInputs.validateResident(builder.policyContext(),policy,manifest,expectedRawTails,
			rawMinLength,rawMaxLength,rawTailCrop,header.width,header.geometry,header.sketchSeed,
			header.positionMode,header.partitions,header.overlap,blocked,keyInputs,effectiveConsensusBonus(header));
		require(header.sameKeyPolicy(header.withKeyPolicy(builder.policyContext())),
			"Resident Stage-B builder key policy differs from the source codec");
	}

	private long effectiveConsensusBonus(ClusterSketchCodec.Header header){
		return consensusBonusSpecified ? consensusMemberBonus : header.consensusMemberBonus;
	}

	private DDLKeyTieBreak.Policy effectiveKeyPolicy(ClusterSketchCodec.Header header){
		return keyPolicySpecified ? keyPolicy : header.keyPolicy;
	}

	/** Raw full-view inputs may equal ordinary inputs, but must never alias outputs. */
	private void protectRawPolicyOutputs(DDLKeyTieBreak.Policy policy,Path manifest,String...outputs)throws Exception{
		if(policy!=DDLKeyTieBreak.Policy.RAW_TAIL_RARITY){return;}
		final PathRegistry registry=PathRegistry.prepare();
		for(Path input:FamilyKeyPolicyInputs.declaredInputs(manifest,keyInputs)){registry.inputs(input.toString());}
		for(String output:outputs){if(output!=null){registry.outputs(output);}}
	}

	private static VerificationMode parseVerification(final String value){
		if("sha80".equalsIgnoreCase(value)){return VerificationMode.SHA80;}
		if("structural".equalsIgnoreCase(value)){return VerificationMode.STRUCTURAL;}
		throw new IllegalArgumentException("verification must be sha80 or structural");
	}

	private static void validateGeometry(ClusterSketchCodec.Header a, ClusterSketchCodec.Header b, String path){
		require(a.sameModelBinding(b),
			"Cluster model binding disagrees in "+path);
	}

	private UnionResults buildUnions(ClusterSketchCodec.Loaded stageB, IntList[] stageAByStageB,
			IntList[] originalsByStageA, Roster stageARoster, ClusteringInput input,
			FamilySketchBuilder prototype){
		final Family[] unions=new Family[stageB.families.length];
		final ConsensusRefiner.Update[] updates=new ConsensusRefiner.Update[stageB.families.length];
		final int threads=Tools.mid(1, workerThreads, Math.max(1, unions.length));
		final AtomicInteger next=new AtomicInteger(0);
		final AtomicReference<Throwable> error=new AtomicReference<Throwable>();
		final TailAlignmentGate gate=gate();
		try(ConsensusRefiner driver=new ConsensusRefiner(1, consensusTrimLeft, trimDepth, null, true)){
			final StageBWorker[] workers=new StageBWorker[threads];
			try{
				for(int i=0; i<workers.length; i++){
					workers[i]=new StageBWorker(next, error, driver.fork(), stageB, stageAByStageB,
						originalsByStageA, stageARoster, input, prototype, gate, unions, updates);
				}
				for(StageBWorker worker : workers){worker.start();}
			}catch(Throwable t){error.compareAndSet(null, t);}
			joinStarted(workers);
			final Throwable thrown=error.get();
			if(thrown!=null){
				if(thrown instanceof RuntimeException){throw (RuntimeException)thrown;}
				if(thrown instanceof Error){throw (Error)thrown;}
				throw new RuntimeException(thrown);
			}
		}
		for(int i=0; i<unions.length; i++){
			require(unions[i]!=null && updates[i]!=null,
				"Stage-B worker omitted family "+stageB.families[i].name);
		}
		return new UnionResults(unions, updates);
	}

	private static void joinStarted(StageBWorker[] workers){
		boolean interrupted=false;
		for(StageBWorker worker : workers){
			if(worker==null || worker.getState()==Thread.State.NEW){continue;}
			while(worker.isAlive()){
				try{worker.join();}
				catch(InterruptedException e){
					interrupted=true; worker.error.compareAndSet(null, e);
				}
			}
		}
		for(StageBWorker worker : workers){
			if(worker!=null && worker.getState()==Thread.State.NEW){worker.closeUnstarted();}
		}
		if(interrupted){Thread.currentThread().interrupt();}
	}

	private IntList[] groupStageAByStageB(int[] ownerByStageA, StageB stageB){
		final IntList[] lists=new IntList[stageB.entries.length];
		for(int i=0; i<lists.length; i++){lists[i]=new IntList();}
		for(int ordinal=0; ordinal<ownerByStageA.length; ordinal++){
			final int dense=stageB.indexFor(ownerByStageA[ordinal]);
			require(dense>=0 && dense<lists.length, "Stage-A row assigned to an unknown Stage-B family");
			lists[dense].add(ordinal);
		}
		for(int dense=0; dense<lists.length; dense++){
			require(lists[dense].size==stageB.entries[dense].members,
				"Stage-B dense membership count disagrees with codec for "+stageB.entries[dense].name);
		}
		return lists;
	}

	private final class StageBWorker extends Thread{
		StageBWorker(AtomicInteger next_, AtomicReference<Throwable> error_,
				ConsensusRefiner refiner_, ClusterSketchCodec.Loaded stageB_,
				IntList[] stageAByStageB_, IntList[] originalsByStageA_, Roster stageARoster_,
				ClusteringInput input_, FamilySketchBuilder prototype_,
				TailAlignmentGate gate_, Family[] unions_,
				ConsensusRefiner.Update[] updates_){
			next=next_; error=error_; refiner=refiner_; stageB=stageB_;
			stageAByStageB=stageAByStageB_; originalsByStageA=originalsByStageA_;
			stageARoster=stageARoster_; input=input_; prototype=prototype_;
			gate=gate_; unions=unions_; updates=updates_;
		}

		@Override
		public void run(){
			try(ConsensusRefiner localRefiner=refiner;
					ClusteringMetrics metrics=new ClusteringMetrics(null, null, true, 1, Long.MAX_VALUE)){
				final FamilySketchBuilder builder=prototype.fork();
				for(int index=next.getAndIncrement(); error.get()==null && index<unions.length;
						index=next.getAndIncrement()){
					final Family family=buildUnion(stageB.families[index], stageAByStageB[index],
						originalsByStageA, stageARoster, input);
					metrics.begin("consensus", 0, family.size());
					final ConsensusRefiner.Update update=localRefiner.updateAfterConcurrentGreedy(family,
						input.sequences, gate, metrics);
					metrics.end(family.size(), family.size(), 0, 1, 1);
					family.refreshStats(input.sequences, builder);
					unions[index]=family; updates[index]=update;
				}
			}catch(Throwable t){error.compareAndSet(null, t);}
		}

		private final AtomicInteger next;
		private final AtomicReference<Throwable> error;
		private final ConsensusRefiner refiner;
		private final ClusterSketchCodec.Loaded stageB;
		private final IntList[] stageAByStageB, originalsByStageA;
		private final Roster stageARoster;
		private final ClusteringInput input;
		private final FamilySketchBuilder prototype;
		private final TailAlignmentGate gate;
		private final Family[] unions;
		private final ConsensusRefiner.Update[] updates;

		void closeUnstarted(){refiner.close();}
	}

	private static final class UnionResults{
		UnionResults(Family[] unions_, ConsensusRefiner.Update[] updates_){
			unions=unions_; updates=updates_;
		}
		final Family[] unions;
		final ConsensusRefiner.Update[] updates;
	}

	private static final class CodecSet{
		CodecSet(ClusterSketchCodec.Loaded[] files_, ClusterSketchCodec.Header header_){
			files=files_; header=header_;
		}
		final ClusterSketchCodec.Loaded[] files;
		final ClusterSketchCodec.Header header;
	}

	private static final class Roster{
		static Roster load(String path, CodecSet stageA){
			final int count=countFamilies(stageA.files), shards=stageA.files.length;
			final Row[] rows=new Row[count];
			final int[] expectedShard=new int[count];
			final ClusterSketchCodec.Entry[] expectedEntry=new ClusterSketchCodec.Entry[count];
			for(int shard=0, ordinal=0; shard<stageA.files.length; shard++){
				for(ClusterSketchCodec.Entry entry : stageA.files[shard].families){
					expectedShard[ordinal]=shard;
					expectedEntry[ordinal]=entry;
					ordinal++;
				}
			}
			final IntHashMap[] ordinals=new IntHashMap[shards];
			for(int i=0; i<ordinals.length; i++){ordinals[i]=new IntHashMap(stageA.files[i].families.length*2+1);}
			final ByteFile bf=ByteFile.makeByteFile(path, false);
			try{
				final LineParser1 parser=new LineParser1('\t');
				requireHeader(bf.nextLine(), ROSTER_HEADER, path);
				int seen=0;
				for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine(), seen++){
					require(seen<rows.length, "Roster has more rows than the Stage-A codec set");
					parser.set(line);
					require(parser.terms()==5, "Stage-A roster rows need five columns");
					final int ordinal=parser.parseInt(0), shard=parser.parseInt(2), familyId=parser.parseInt(3);
					require(ordinal==seen && shard>=0 && shard<shards, "Stage-A roster is out of order");
					final ClusterSketchCodec.Entry expected=expectedEntry[seen];
					require(shard==expectedShard[seen], "Roster no longer matches the ordered Stage-A codec shard");
					require(familyId==expected.familyId && parser.parseInt(4)==expected.members,
						"Roster no longer matches the ordered Stage-A codec list");
					require(expected.members>0, "Stage-A roster contains an empty family");
					final String id=parser.parseString(1);
					require(id.equals(stageAId(shard, familyId)), "Roster Stage-A ID is not canonical");
					rows[seen]=new Row(id, shard, familyId, expected.members);
					require(ordinals[shard].get(familyId)<0, "Roster repeats Stage-A family "+id);
					ordinals[shard].put(familyId, ordinal);
				}
				require(seen==rows.length, "Roster has fewer rows than the Stage-A codec set");
			}finally{close(bf, path);}
			return new Roster(rows, ordinals);
		}

		private Roster(Row[] rows_, IntHashMap[] ordinals_){
			rows=rows_; ordinals=ordinals_;
		}

		int ordinalFor(int shard,int familyId){
			return shard>=0 && shard<ordinals.length ? ordinals[shard].get(familyId) : -1;
		}

		final Row[] rows;
		private final IntHashMap[] ordinals;
	}

	private static int countFamilies(ClusterSketchCodec.Loaded[] loaded){
		int count=0;
		for(ClusterSketchCodec.Loaded file : loaded){count=Math.addExact(count, file.families.length);}
		return count;
	}

	private static final class Row{
		Row(String id_, int shard_, int familyId_, int sourceMembers_){
			id=id_; shard=shard_; familyId=familyId_; sourceMembers=sourceMembers_;
		}
		final String id;
		final int shard, familyId, sourceMembers;
	}

	private static final class StageB{
		static StageB load(ClusterSketchCodec.Loaded loaded){
			final IntHashMap index=new IntHashMap(loaded.families.length*2+1);
			for(int i=0; i<loaded.families.length; i++){
				require(index.get(loaded.families[i].familyId)<0,
					"Stage-B codec repeats family "+loaded.families[i].familyId);
				index.put(loaded.families[i].familyId, i);
			}
			return new StageB(loaded.families, index);
		}
		private StageB(ClusterSketchCodec.Entry[] entries_, IntHashMap index_){
			entries=entries_; index=index_;
		}
		int indexFor(int familyId){return index.get(familyId);}
		final ClusterSketchCodec.Entry[] entries;
		private final IntHashMap index;
	}

	/** Shared package guard for Stage-B and global selection input/output aliases. */
	static final class PathRegistry{
		static PathRegistry prepare(){return new PathRegistry();}

		PathRegistry inputs(String...paths){
			if(paths!=null){for(String path : paths){input(path);}}
			return this;
		}

		PathRegistry outputs(String...paths){
			if(paths!=null){for(String path : paths){output(path);}}
			return this;
		}

		private void input(String name){
			if(name==null){return;}
			requireDisk(name);
			final Path path=Paths.get(name).toAbsolutePath().normalize();
			final Path real=real(path, name);
			require(Files.isRegularFile(real), "Stage-B input is not a regular file: "+name);
			add(real, name);
		}

		private void output(String name){
			require(name!=null, "Stage-B output path is missing");
			requireDisk(name);
			final Path path=Paths.get(name).toAbsolutePath().normalize();
			final Path parent=path.getParent();
			require(parent!=null, "Stage-B output lacks a parent directory: "+name);
			final Path realParent=real(parent, name);
			require(Files.isDirectory(realParent), "Stage-B output parent is not a directory: "+name);
			final Path resolved=resolveOutput(path, realParent, name);
			add(resolved, name);
		}

		private static Path resolveOutput(Path path, Path realParent, String name){
			try{
				if(Files.exists(path, LinkOption.NOFOLLOW_LINKS)){
					final Path real=path.toRealPath();
					require(Files.isRegularFile(real), "Stage-B output target is not a regular file: "+name);
					return real;
				}
			}catch(NoSuchFileException e){throw new IllegalArgumentException("Dangling Stage-B output symlink: "+name, e);}
			catch(IOException e){throw new IllegalArgumentException("Could not resolve Stage-B output: "+name, e);}
			return realParent.resolve(path.getFileName().toString()).normalize();
		}

		private void add(Path path,String name){
			for(Path old : paths){
				require(!sameFile(old, path), "Stage-B input/output paths alias: "+name);
			}
			paths.add(path);
		}

		private static Path real(Path path,String name){
			try{return path.toRealPath();}
			catch(IOException e){throw new IllegalArgumentException("Could not resolve Stage-B path: "+name, e);}
		}

		private static boolean sameFile(Path a,Path b){
			if(a.equals(b)){return true;}
			try{return Files.isSameFile(a, b);}
			catch(IOException e){return false;}
		}

		private static void requireDisk(String name){
			final Path fileName=Paths.get(name).getFileName();
			final String base=fileName==null ? name : fileName.toString();
			require(!FileFormat.isStdio(base) && !"/dev/null".equals(name) &&
				!"devnull".equalsIgnoreCase(base) && !"null".equalsIgnoreCase(base),
				"Stage-B requires on-disk files: "+name);
		}

		private final ArrayList<Path> paths=new ArrayList<Path>();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Helpers           ----------------*/
	/*--------------------------------------------------------------*/

	private TailAlignmentGate gate(){
		return new TailAlignmentGate(alignmentMode(),denominator,maxAlignLength,minSpan,minIdentity,minRatio,null,matrix,minPaired);
	}
	private TailAlignmentGate.Mode alignmentMode(){
		if(explicitAlignmentMode!=null){
			require(alignmentAnchor==AlignmentAnchor.RIGHT && trimLeft && denominator==TailAlignmentGate.Denominator.ALIGNED_QUERY,
				"Explicit right-pinned alignmentmode requires right anchor, trimleft=t and denominator=aligned");
			return explicitAlignmentMode;
		}
		if(alignmentAnchor==AlignmentAnchor.LEFT){return TailAlignmentGate.Mode.PIN_LEFT_TRIM_RIGHT;}
		if(alignmentAnchor==AlignmentAnchor.RIGHT){
			return trimLeft ? TailAlignmentGate.Mode.PIN_RIGHT_TRIM_LEFT :
				TailAlignmentGate.Mode.PIN_RIGHT_FULL_QUERY;
		}
		assert(alignmentAnchor==AlignmentAnchor.NONE) : "Unhandled anchor terminus: "+alignmentAnchor;
		return TailAlignmentGate.Mode.FULL_QUERY_UNPINNED;
	}
	private boolean pinRight(){return alignmentAnchor==AlignmentAnchor.RIGHT;}
	private static AlignmentAnchor parseAlignmentAnchor(final String value){
		if(value.equalsIgnoreCase("left") || value.equalsIgnoreCase("start") ||
				value.equalsIgnoreCase("head")){return AlignmentAnchor.LEFT;}
		if(value.equalsIgnoreCase("right") || value.equalsIgnoreCase("end") ||
				value.equalsIgnoreCase("tail")){return AlignmentAnchor.RIGHT;}
		if(value.equalsIgnoreCase("none") || value.equalsIgnoreCase("unpinned")){
			return AlignmentAnchor.NONE;
		}
		throw new IllegalArgumentException("anchorterminus requires left/start/head, right/end/tail, or none/unpinned");
	}
	private static AlignmentAnchor resolveAlignmentAnchor(final Boolean parsedPinRight,final Boolean parsedPinLeft){
		if(Boolean.TRUE.equals(parsedPinLeft) && Boolean.TRUE.equals(parsedPinRight)){
			throw new IllegalArgumentException("pinleft=t and pinright=t are contradictory");
		}
		return Boolean.TRUE.equals(parsedPinLeft) ? AlignmentAnchor.LEFT :
			Boolean.FALSE.equals(parsedPinRight) ? AlignmentAnchor.NONE : AlignmentAnchor.RIGHT;
	}
	private static AlignmentAnchor resolveAlignmentAnchor(final AlignmentAnchor anchor,
			final Boolean parsedPinRight,final Boolean parsedPinLeft){
		if(Boolean.TRUE.equals(parsedPinLeft) && anchor!=AlignmentAnchor.LEFT){
			throw new IllegalArgumentException("pinleft=t contradicts anchorterminus="+anchor);
		}
		if(Boolean.TRUE.equals(parsedPinRight) && anchor!=AlignmentAnchor.RIGHT){
			throw new IllegalArgumentException("pinright=t contradicts anchorterminus="+anchor);
		}
		if(Boolean.FALSE.equals(parsedPinLeft) && anchor==AlignmentAnchor.LEFT){
			throw new IllegalArgumentException("pinleft=f contradicts anchorterminus=LEFT");
		}
		if(Boolean.FALSE.equals(parsedPinRight) && anchor==AlignmentAnchor.RIGHT){
			throw new IllegalArgumentException("pinright=f contradicts anchorterminus=RIGHT");
		}
		return anchor;
	}

	private static void appendConsensus(ClusterSketchCodec.Entry entry, ByteBuilder out){
		out.append(ClusteringInput.decode(entry.consensus)).nl();
	}

	private static String[] split(String csv, String name){
		require(csv!=null && !csv.isEmpty(), "Missing "+name);
		final String[] split=csv.split(",");
		for(String s : split){require(s!=null && !s.isEmpty(), "Empty path in "+name);}
		return split;
	}

	private static void configureProteinInput(){
		Shared.AMINO_IN=true;
		Shared.TRIM_READ_DESCRIPTION=false;
		Read.TO_UPPER_CASE=false;
		Read.LOWER_CASE_TO_N=false;
		Read.DOT_DASH_X_TO_N=false;
		Read.FIX_HEADER=false;
		// Compression settings are caller policy parsed before process() starts.
		Read.JUNK_MODE=Read.IGNORE_JUNK;
	}

	private static String stageAId(int shard,int familyId){
		return "s"+shard+"_f"+familyId;
	}

	private ByteStreamWriter open(String path){
		final ByteStreamWriter writer=new ByteStreamWriter(path, overwrite, false, true);
		writer.start();
		return writer;
	}

	private static void close(ByteStreamWriter writer, String path){
		if(writer!=null && writer.poisonAndWait()){throw new IllegalStateException("Failed Stage-B output: "+path);}
	}

	private static void close(ByteStreamWriter a, String aPath, ByteStreamWriter b, String bPath){
		RuntimeException error=null;
		try{close(a, aPath);}
		catch(RuntimeException e){error=e;}
		try{close(b, bPath);}
		catch(RuntimeException e){
			if(error==null){error=e;}
			else{error.addSuppressed(e);}
		}
		if(error!=null){throw error;}
	}

	private static void close(ByteFile file, String path){
		if(file.close()){throw new IllegalStateException("Failed Stage-B input: "+path);}
	}

	private static void requireHeader(byte[] line, String expected, String path){
		require(line!=null && expected.equals(new String(line)), "Unexpected header in "+path);
	}

	private static int[] readFamilyIds(final String path){
		final ByteFile file=ByteFile.makeByteFile(path,true);
		final LineParser1 parser=new LineParser1('\t');
		final IntList ids=new IntList();
		boolean header=false;
		long lineNumber=0;
		try{
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				lineNumber++;
				if(line.length<1 || line[0]=='#'){continue;}
				parser.set(line);
				require(parser.terms()==1, path+":"+lineNumber+": expected exactly one family_id column");
				if(!header){
					require(parser.termEquals("family_id",0), path+":"+lineNumber+": expected family_id header");
					header=true; continue;
				}
				final int id=parseNonnegativeInt(line, path, lineNumber);
				ids.add(id);
			}
		}finally{close(file,path);}
		require(header, path+": expected family_id header");
		return ids.toArray();
	}

	private static int parseNonnegativeInt(final byte[] line, final String path, final long lineNumber){
		long value=0;
		for(byte b : line){
			require(Tools.isDigit(b), path+":"+lineNumber+": expected a nonnegative integer family_id");
			value=value*10+(b-'0');
			require(value<=Integer.MAX_VALUE, path+":"+lineNumber+": family_id exceeds Integer.MAX_VALUE");
		}
		return (int)value;
	}

	private void writeFamilyIds(final String path, final int[] ids){
		final ByteStreamWriter writer=open(path);
		try{
			writer.println("family_id");
			final ByteBuilder row=new ByteBuilder();
			for(int id:ids){row.clear().append(id);writer.println(row);}
		}finally{close(writer,path);}
	}

	private static void require(boolean ok,String message){
		if(!ok){throw new IllegalArgumentException(message);}
	}

	private String mode, in, out, livePhase;
	private String codecs, stageACodecs, stageAMembers, stageBCodec, stageBMembers;
	private double initialFraction=1.2,pricedFraction=LayerCakeClustering.Options.DEFAULT_PRICED_FRACTION;
	private LiveProbeLayer.PriceMode priceMode=LiveProbeLayer.PriceMode.FLOOR;
	private LiveProbeLayer.BelowMode belowMode=LiveProbeLayer.BelowMode.PRIORSCALED;
	private int setupSlotCost=-1;
	private double edgeFraction=LiveProbeLayer.DEFAULT_EDGE_FRACTION;
	private double sizeLambda=LiveProbeLayer.DEFAULT_SIZE_LAMBDA;
	private boolean priceModeSet,setupSlotCostSet,edgeFractionSet,belowModeSet,sizeLambdaSet;
	private boolean initialFractionSet,pricedFractionSet,nonBudgetStructuralConfigured;
	private String partitionMap, roster, report, blacklist, blacklistSha80;
	private String matrixName=ProteinSubstitutionMatrix.DEFAULT_NAME, matrixFile;
	private DDLKeyTieBreak.Policy keyPolicy=DDLKeyTieBreak.Policy.SMALLEST;
	private boolean keyPolicySpecified;
	private long consensusMemberBonus;
	private boolean consensusBonusSpecified;
	private String rawManifest;
	private long expectedRawTails=-1;
	private int rawMinLength=1, rawMaxLength=Integer.MAX_VALUE, rawTailCrop=0;
	private final StructuralInputs keyInputs=new StructuralInputs();
	private int blacklistWords=-1,blacklistK=-1;
	private String blacklistMode,blacklistProvenance;
	private VerificationMode verification=VerificationMode.SHA80;
	private String passesReport, membersOut, rebuildReport, auditReport, featuresReport, refinementReport;
	private String wideFeatures,featureContext;
	private int featureLayer;
	private boolean featureLayerSet;
	private WideProbeFeatures.Output wideOutput;
	private String fusionReport,liveMetrics,liveSizeHist,fitReport,stateOut,fissionParents,fusionProductsOut;
	private boolean currentState,currentStateSet,nextFamilyIdSet;
	private final boolean residentLayer;
	private int nextFamilyId=-1;
	private int localRefineRounds;
	private boolean probeStratified;
	private String probePlan;
	private String selectedSeeds,baselineMembers;
	private boolean localConfigured,localFlex=true,localPlanConfigured,probeSeedSet,memberRanking,memberRankingSet;
	private double memberHydroWeight=MemberRankInputs.DEFAULT_HYDRO_WEIGHT;
	private int selectedCount=80000,localMinLength=100,localMinimumHits=4;
	private int fissionProbes=-1,fusionProbes=-1;
	private long probeSeed;
	private double localHydroDelta=.55;
	private int mergeMaxSize=20, mergeShortlist=50;
	private double mergeIdentity=0.50, mergeRatio=0.14;
	private AlignmentAnchor alignmentAnchor=AlignmentAnchor.RIGHT;
	private boolean overwrite;
	private boolean trimLeft=true, consensusTrimLeft=true;
	private int minSpan=70, minPaired=0, maxAlignLength=Integer.MAX_VALUE;
	private double minIdentity=0.35, minRatio=0.14;
	private float trimDepth=0;
	private TailAlignmentGate.Denominator denominator=TailAlignmentGate.Denominator.ALIGNED_QUERY;
	private TailAlignmentGate.Mode explicitAlignmentMode;
	private enum AlignmentAnchor{RIGHT, LEFT, NONE}
	private int workerThreads=1;
	private int frozenIndexThreads=1;
	private boolean mergeMaxSizeSet,structuralConfigured;
	private boolean probePlanSet,scheduleControlsSet,objectiveWeightsSet;
	private boolean livePhaseSet,gainMarginSet,fittedWeightsSet,temperatureSet,localRefineRoundsSet;
	private boolean liveMarginSlot,fissionMinSizeSet,liveSelectionSet;
	private boolean eligibilityFractionSet;
	private boolean rankByPrediction;
	private boolean slotFractionSet;
	private double slotFraction=1;
	private double fissionTopFraction=1,fusionBottomFraction=1;
	private LayerCakeSchedule.Mode scheduleMode;
	private LocalFamilyObjective.SizeForm sizeForm=LocalFamilyObjective.SizeForm.LOG2;
	private double temperature=1,minimumTemperature=0,cooling=.5,structuralTarget=.8;
	private double gainMargin=1,cohesionScale=1,separationScale=1,sizeScale=1,fusionSizeRatio=8;
	private int familyBudget=80000,maxCycles=0,largeSize=100,fissionMinSize=100,minChildSize=5;
	private int memberNeighbors=8,memberMinHits=1,evidenceNeighbors=10;
	private PrintStream outstream;
	private final ProteinSubstitutionMatrix matrix;

	private static final String ROSTER_HEADER="stagea_ordinal\tstagea_id\tshard\tstagea_family_id\tsource_members";
	private static final String REPORT_HEADER="stageb_family_id\tstageb_id\tstagea_consensuses\t"+
		"source_union_members\tunique_source_members\tstatus\treason\tproposals\taccepted\trejected\t"+
		"unchanged\tincompatible_snapshots\tconstruction_alignments\tk\twidth\tsketchseed\t"+
		"positional\tpositionmode\tpartitions\toverlap\tblacklist_sha80\tminspan\tmaxalignlen\tminid\t"+
		"minblosumratio\tdenominator\tpinright\ttrimleft\tconsensustrimleft\ttrimdepth";

	static final class LayerResult{
		LiveProbeLayer.Result result(){return result;}
		long elapsedNanos(){return elapsedNanos;}

		private LayerResult(final LiveProbeLayer.Result result_,final LiveProbeMetrics.Snapshot before_,
				final LiveProbeMetrics.Snapshot afterFission_,final LiveProbeMetrics.Snapshot afterFusion_,
				final int[] fusionProductIds_,final long started_){
			result=result_;before=before_;afterFission=afterFission_;afterFusion=afterFusion_;
			fusionProductIds=fusionProductIds_;started=started_;
			require(result!=null && before!=null && afterFission!=null && afterFusion!=null &&
				fusionProductIds!=null,"Resident layer result needs complete states and fusion-product IDs");
		}

		private final LiveProbeLayer.Result result;
		private final LiveProbeMetrics.Snapshot before,afterFission,afterFusion;
		private final int[] fusionProductIds;
		private final long started;
		private long elapsedNanos=-1;
	}
}
