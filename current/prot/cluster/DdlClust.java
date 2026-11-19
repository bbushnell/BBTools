package prot.cluster;

import java.io.PrintStream;
import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.NoSuchFileException;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;

import ddl.FamilySketchBuilder;
import ddl.DDLKeyTieBreak;
import ddl.DDLFamilySearch.StructuralInputs;
import ddl.DDLPentamerBlacklist;
import ddl.ProteinKmerGeometry;
import ddl.ProteinBlacklistPolicy;
import ddl.ProteinPositionMode;
import fileIO.ReadWrite;
import fileIO.ByteStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prot.DigestSuffix;
import prot.FamilyCalibration.VerificationMode;
import prot.ProteinSubstitutionMatrix;
import prot.FamilyKeyPolicyInputs;
import shared.Resources;
import shared.Shared;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import structures.IntHashMap;

/** Parameterized greedy protein clustering with parallel frozen passes.
 * Native argument parsing and FASTA streaming follow template.A_SampleStreamerMT;
 * greedy admission is serial by default, with explicit asynchronous workers available.
 * Frozen reassignment and family rebuilds also use private workers.
 * Ordered in=file1,file2 lists share one ordinal space; tailcrop=N optionally
 * retains the last N biological residues during loading, before all clustering.
 * @author Keqing
 */
public final class DdlClust {

	/** Reports expected resource absence after normal cleanup; malformed inputs still propagate. */
	public static void main(String[] args){
		try{
			final PreParser pre=new PreParser(args, DdlClust.class, false);
			try{
				if(DdlClustCake.requested(pre.args)){new DdlClustCake(pre.args, pre.outstream).process();}
				else{new DdlClust(pre).process();}
			}finally{Shared.closeStream(pre.outstream);}
		}catch(Resources.MissingResourceException e){
			System.err.println(e.getMessage());
			System.exit(1);
		}
	}

	public DdlClust(String[] args){
		this(new PreParser(args, DdlClust.class, false));
	}

	private DdlClust(final PreParser pre){
		this(pre,null,null);
	}

	/** Runs the ordinary assignment engine on borrowed seeds/input; no input is reread.
	 * Null input is allowed only for the first assignment and loads the population once.
	 */
	static Assignment assignResident(final String[] args,final ClusterSketchCodec.Loaded seeds,final ClusteringInput input){
		return assignResident(args,seeds,input,null,0);
	}

	/** Ranked continuations borrow the identical cache and population; they never hydrate again. */
	static Assignment assignResident(final String[] args,final ClusterSketchCodec.Loaded seeds,
			final ClusteringInput input,final MemberRankInputs inputs,final int layer){
		return assignResident(args,seeds,input,inputs,layer,null);
	}

	/** A Cake continuation borrows its once-loaded immutable alignment matrix. */
	static Assignment assignResident(final String[] args,final ClusterSketchCodec.Loaded seeds,
			final ClusteringInput input,final MemberRankInputs inputs,final int layer,final ProteinSubstitutionMatrix matrix){
		return assignResident(args,seeds,input,inputs,layer,matrix,null,null);
	}

	/** Reuses the original immutable key counts and their complete input-stability tracker. */
	static Assignment assignResident(final String[] args,final ClusterSketchCodec.Loaded seeds,
			final ClusteringInput input,final MemberRankInputs inputs,final int layer,final ProteinSubstitutionMatrix matrix,
			final DDLKeyTieBreak.PolicyContext keyContext,final StructuralInputs keyInputs){
		if(seeds==null){throw new IllegalArgumentException("Resident assignment needs validated fixed seeds");}
		final DdlClust tool=new DdlClust(new PreParser(args,DdlClust.class,false),seeds,null,matrix,keyContext,keyInputs);
		try{return tool.process(input,inputs,layer);}finally{Shared.closeStream(tool.outstream);}
	}

	/** Delay original-population loading and ranked selection until constructor normalization is complete. */
	static Assignment assignInitialRanked(final String[] args,final String codec,final int budget,final double fraction){
		return assignInitialRanked(args,new String[]{codec},budget,fraction);
	}

	/** Delay original-population loading and ranked selection until constructor normalization is complete. */
	static Assignment assignInitialRanked(final String[] args,final String[] codecs,final int budget,final double fraction){
		return assignInitialRanked(args,codecs,budget,fraction,null);
	}

	/** The initial resident pass uses the same score/gap binding as all later layers. */
	static Assignment assignInitialRanked(final String[] args,final String[] codecs,final int budget,
			final double fraction,final ProteinSubstitutionMatrix matrix){
		final SeedSelection selection=new SeedSelection(codecs,budget,fraction);
		final DdlClust tool=new DdlClust(new PreParser(args,DdlClust.class,false),null,selection,matrix);
		try{return tool.process(null,null,0);}finally{Shared.closeStream(tool.outstream);}
	}

	private DdlClust(final PreParser pre,final ClusterSketchCodec.Loaded residentSeeds_,final SeedSelection selection_){
		this(pre,residentSeeds_,selection_,null);
	}

	private DdlClust(final PreParser pre,final ClusterSketchCodec.Loaded residentSeeds_,
			final SeedSelection selection_,final ProteinSubstitutionMatrix residentMatrix){
		this(pre,residentSeeds_,selection_,residentMatrix,null,null);
	}

	private DdlClust(final PreParser pre,final ClusterSketchCodec.Loaded residentSeeds_,
			final SeedSelection selection_,final ProteinSubstitutionMatrix residentMatrix,
			final DDLKeyTieBreak.PolicyContext residentKeyContext,final StructuralInputs residentKeyInputs){
		if((residentKeyContext==null)!=(residentKeyInputs==null)){
			throw new IllegalArgumentException("Resident key context and its original input tracker must travel together");
		}
		keyInputs=residentKeyInputs==null ? new StructuralInputs() : residentKeyInputs;
		residentSeeds=residentSeeds_;seedSelection=selection_;
		final String[] args=pre.args; outstream=pre.outstream;
		final Parser parser=new Parser(); parser.overwrite=true;
		ProteinPositionMode explicitPositionMode=null;
		AlignmentAnchor parsedAlignmentAnchor=null;
		Boolean parsedPinRight=null, parsedPinLeft=null;
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String key=pair[0].toLowerCase(Locale.ROOT);
			final String value=pair.length>1 && !pair[1].equalsIgnoreCase("null") ? pair[1] : null;
			if(key.equals("memberranking")){memberRanking=Parse.parseBoolean(value);}
			else if(key.equals("memberhydroweight")){memberHydroWeight=MemberRankInputs.checkedHydroWeight(Double.parseDouble(value));}
			else if(key.equals("seedmembers")){seedMembers=value;}
			else if(key.equals("seedmemberssha80")){seedMembersPin=value;}
			else if(key.equals("keytiebreak")){keyPolicy=DDLKeyTieBreak.Policy.parse(value); keyPolicySpecified=true;}
			else if(key.equals("consensusbonus")){consensusMemberBonus=Parse.parseKMG(value);consensusBonusSpecified=true;}
			else if(key.equals("rawmanifest")){rawManifest=value;}
			else if(key.equals("expectedrawtails")){expectedRawTails=Parse.parseKMG(value);}
			else if(key.equals("rawminlen")){rawMinLength=Parse.parseIntKMG(value);}
			else if(key.equals("rawmaxlen")){rawMaxLength=Parse.parseIntKMG(value);}
			else if(key.equals("rawtailcrop")){rawTailCrop=Parse.parseIntKMG(value);}
			else if(key.equals("minlen")){minLength=Parse.parseIntKMG(value);}
			else if(key.equals("minspan")){minSpan=Parse.parseIntKMG(value);}
			else if(key.equals("minpaired")){minPaired=Parse.parseIntKMG(value);}
			else if(key.equals("maxalignlen")){maxAlignLength=Parse.parseIntKMG(value);}
			else if(key.equals("tailcrop")){tailCrop=Parse.parseIntKMG(value);}
			else if(key.equals("mincluster")){minCluster=Parse.parseIntKMG(value);}
			else if(key.equals("prunefinal")){pruneFinal=Parse.parseBoolean(value);}
			else if(key.equals("preserveallseeds")){preserveAllSeeds=Parse.parseBoolean(value);}
			else if(key.equals("greedythreads")){greedyThreads=Parse.parseIntKMG(value); greedyThreadsSet=true;}
			else if(key.equals("frozenindexthreads")){frozenIndexThreads=Parse.parseIntKMG(value);}
			else if(key.equals("shortlist")){shortlist=Parse.parseIntKMG(value);}
			else if(key.equals("passes")){passes=Parse.parseIntKMG(value);}
			else if(key.equals("consensuspasses")){consensusPasses=Parse.parseIntKMG(value);}
			else if(key.equals("topn")){topN=Parse.parseIntKMG(value);}
			else if(key.equals("seedlimit")){seedLimit=Parse.parseIntKMG(value);}
			else if(key.equals("passchanges")){passChanges=value;}
			else if(key.equals("stopchangefraction")){stopChangeFraction=Double.parseDouble(value);}
			else if(key.equals("affectedfamilies")){affectedFamilies=value; affectedFamiliesSet=true;}
			else if(key.equals("affectedcounts")){affectedCounts=value; affectedCountsSet=true;}
			else if(key.equals("skipfinalrebuild")){skipFinalRebuild=Parse.parseBoolean(value);}
			else if(key.equals("minid")){minIdentity=Double.parseDouble(value);}
			else if(key.equals("minblosumratio")){minRatio=Double.parseDouble(value);}
			else if(key.equals("matrix")){matrixName=value;}
			else if(key.equals("matrixfile")){matrixFile=value;}
			else if(key.equals("mergeroster")){mergeRoster=value;}
			else if(key.equals("mergemaxsize")){mergeMaximumSize=Parse.parseIntKMG(value);}
			else if(key.equals("mergeminid")){
				mergeMinIdentity=Double.parseDouble(value); mergeMinIdentitySet=true;
				if(!Double.isFinite(mergeMinIdentity)){
					throw new IllegalArgumentException("mergeminid must be finite: "+mergeMinIdentity);
				}
			}
			else if(key.equals("mergeminratio")){
				mergeMinRatio=Double.parseDouble(value); mergeMinRatioSet=true;
				if(!Double.isFinite(mergeMinRatio)){
					throw new IllegalArgumentException("mergeminratio must be finite: "+mergeMinRatio);
				}
			}
			else if(key.equals("minhits")){minimumHits=Parse.parseIntKMG(value); passHitsFlex=false;}
			else if(key.equals("passhits")){
				if(value==null){throw new IllegalArgumentException("passhits requires flex or a nonnegative integer");}
				passHitsFlex="flex".equalsIgnoreCase(value);
				if(!passHitsFlex){minimumHits=Parse.parseIntKMG(value);}
			}
			else if(key.equals("greedyhits")){
				if(value==null){throw new IllegalArgumentException("greedyhits requires flex, ramp or a nonnegative integer");}
				greedyRamp="ramp".equalsIgnoreCase(value);
				greedyHits=greedyRamp || "flex".equalsIgnoreCase(value) ? -1 : Parse.parseIntKMG(value);
			}
			else if(key.equals("rampfamilies")){rampFamilies=Parse.parseIntKMG(value);}
			else if(key.equals("rampminmembers")){rampMinMembers=Parse.parseIntKMG(value);}
			else if(key.equals("maxhydrodelta")){maximumHydroDelta=Double.parseDouble(value);}
			else if(key.equals("greedymaxhydrodelta")){greedyHydroDelta=Double.parseDouble(value);}
			else if(key.equals("pinright")){parsedPinRight=Parse.parseBoolean(value);}
			else if(key.equals("alignmentmode")){explicitAlignmentMode=TailAlignmentGate.parseRightMode(value);}
			else if(key.equals("pinleft")){parsedPinLeft=Parse.parseBoolean(value);}
			else if(key.equals("anchorterminus")){
				if(value==null){throw new IllegalArgumentException("anchorterminus requires left, right or none");}
				parsedAlignmentAnchor=parseAlignmentAnchor(value);
			}
			else if(key.equals("indexed")){indexed=Parse.parseBoolean(value);}
			else if(key.equals("concurrentgreedy")){concurrentGreedy=Parse.parseBoolean(value);}
			else if(key.equals("merge")){merge=Parse.parseBoolean(value);}
			else if(key.equals("scavenge")){
				if(value!=null && !value.equalsIgnoreCase("t") && !value.equalsIgnoreCase("true") && !value.equals("1") &&
					!value.equalsIgnoreCase("f") && !value.equalsIgnoreCase("false") && !value.equals("0")){
					throw new IllegalArgumentException("scavenge supports only boolean singleton replay, not size-two replay");
				}
				scavenge=Parse.parseBoolean(value);
			}
			else if(key.equals("snapshotfraction")){snapshotFraction=Double.parseDouble(value);}
			else if(key.equals("progress")){progress=Parse.parseBoolean(value);}
			else if(key.equals("trimleft")){leftClip=Parse.parseBoolean(value);}
			else if(key.equals("consensustrimleft")){consensusTrimLeft=Parse.parseBoolean(value);}
			else if(key.equals("trimdepth")){trimDepth=Float.parseFloat(value);}
			else if(key.equals("denominator")){
				if("aligned".equalsIgnoreCase(value)){denominator=TailAlignmentGate.Denominator.ALIGNED_QUERY;}
				else if("full".equalsIgnoreCase(value)){denominator=TailAlignmentGate.Denominator.FULL_QUERY;}
				else{throw new IllegalArgumentException("denominator must be aligned or full");}
			}else if(key.equals("seed")){seed=Long.parseLong(value);}
			else if(key.equals("order")){
				if(value==null){throw new IllegalArgumentException("order requires shuffle, hydro or hydrobin");}
				inputOrder=ClusteringOrder.Mode.valueOf(value.toUpperCase(Locale.ROOT));
			}
			else if(key.equals("sketchseed")){sketchSeed=Long.parseLong(value);}
			else if(key.equals("k")){k=Parse.parseIntKMG(value);}
			else if(key.equals("aabits") || key.equals("aapartition")){
				if(aaBits!=null){throw new IllegalArgumentException("Use either aabits or aapartition, not both");}
				aaBits=value;
			}
			else if(key.equals("width")){width=Parse.parseIntKMG(value);}
			else if(key.equals("positional")){positional=Parse.parseBoolean(value);}
			else if(key.equals("positionmode") || key.equals("positionalmode")){explicitPositionMode=ProteinPositionMode.parse(value);}
			else if(key.equals("partitions")){partitions=Parse.parseIntKMG(value);}
			else if(key.equals("overlap")){overlap=Parse.parseIntKMG(value);}
			else if(key.equals("blacklist")){blacklistPath=value;}
			else if(key.equals("blacklistsha80")){blacklistPin=value;}
			else if(key.equals("blacklistk")){blacklistK=Parse.parseIntKMG(value);}
			else if(key.equals("blacklistmode")){blacklistMode=value;}
			else if(key.equals("blacklistprovenance")){blacklistProvenance=value;}
			else if(key.equals("blacklistwords") || key.equals("blacklist_words")){blacklistWords=Parse.parseIntKMG(value);}
			else if(key.equals("verification")){verification=parseVerification(value);}
			else if(key.equals("members")){members=value;}
			else if(key.equals("ddlout")){ddlout=value;}
			else if(key.equals("consensus")){consensus=value;}
			else if(key.equals("stats")){stats=value;}
			else if(key.equals("refinement")){refinement=value;}
			else if(key.equals("checkpoints")){checkpoints=value;}
			else if(key.equals("hist")){hist=value;}
			else if(key.equals("pairhist")){pairhist=value;}
			else if(key.equals("assignmetrics")){assignmetrics=value;}
			else if(key.equals("quality")){quality=value;}
			else if(key.equals("seeds")){seedPath=value;}
			else if(key.equals("selectedseeds")){selectedSeeds=value;}
			else if(key.equals("seedoverlay")){seedOverlayPath=value; seedOverlaySet=true;}
			else if(key.equals("seedremove")){seedRemove=value; seedRemoveSet=true;}
			else if(key.equals("seedexclude")){seedExclude=value; seedExcludeSet=true;}
			else if(key.equals("decisionaudit")){decisionaudit=value;}
			else if(key.equals("greedyaudit")){greedyaudit=value;}
			else if(key.equals("auditoldmembers")){auditoldmembers=value;}
			else if(key.equals("auditqueries")){auditqueries=value;}
			else if(key.equals("checkpointevery")){checkpointEvery=Parse.parseKMG(value);}
			else if(parser.parse(arg, key, value)){}
			else if(i==0 && Tools.looksLikeInputSequenceStream(arg)){parser.in1=arg;}
			else{throw new IllegalArgumentException("Unknown clustering parameter: "+arg);}
		}
		input=parser.in1; output=parser.out1; overwrite=parser.overwrite; maxReads=parser.maxReads;
		workerThreads=Shared.threads();
		alignmentAnchor=parsedAlignmentAnchor==null ?
			resolveAlignmentAnchor(parsedPinRight,parsedPinLeft) :
			resolveAlignmentAnchor(parsedAlignmentAnchor,parsedPinRight,parsedPinLeft);
		if(!memberRanking && memberHydroWeight!=MemberRankInputs.DEFAULT_HYDRO_WEIGHT){
			throw new IllegalArgumentException("Nondefault memberhydroweight requires memberranking=t");
		}
		if(!memberRanking && (seedMembers!=null || seedMembersPin!=null || seedSelection!=null)){
			throw new IllegalArgumentException("Seed hydration requires explicit memberranking=t");
		}
		if(memberRanking){
			if(maxReads>=0 || merge || mergeRoster!=null || seedOverlaySet || seedRemoveSet || seedExcludeSet ||
				seedSelection!=null && (seedLimit!=0 || affectedFamiliesSet || seedPath!=null) ||
				residentSeeds==null && seedSelection==null && seedLimit!=0){
				throw new IllegalArgumentException("Member ranking requires complete original inputs; overlays, exclusions, merging and file-seed limits are not bound");
			}
			if(!hasSeeds()){
				if(seedMembers!=null || seedMembersPin!=null){
					throw new IllegalArgumentException("Raw member ranking has no seeds to hydrate; seedmembers and seedmemberssha80 require seeds");
				}
			}else if(residentSeeds==null){
				final boolean hydrating=seedMembers!=null || seedMembersPin!=null;
				if(hydrating){
					if(seedSelection!=null && seedSelection.codecs.length!=1 ||
						seedMembers==null || !Files.isRegularFile(Paths.get(seedMembers)) ||
						seedMembersPin==null || !seedMembersPin.matches("[0-9a-f]{20}")){
						throw new IllegalArgumentException("Initial member ranking needs seedmembers and seedmemberssha80 for the original complete population");
					}
				}
			}else if(seedMembers!=null || seedMembersPin!=null){
				throw new IllegalArgumentException("Resident continuation must reuse its original cache, not rehydrate seedmembers");
			}
		}
		if(workerThreads<1){throw new IllegalArgumentException("Clustering needs at least one worker");}
		if(frozenIndexThreads<1 || frozenIndexThreads>workerThreads ||
			frozenIndexThreads>1 && (!hasSeeds() || !indexed || merge)){
			throw new IllegalArgumentException("frozenindexthreads must lie in [1,t]; parallel construction requires seeds, indexed=t and merge=f");
		}
		if(greedyThreadsSet && (!concurrentGreedy || greedyThreads<1 || greedyThreads>workerThreads)){
			throw new IllegalArgumentException("greedythreads requires concurrentgreedy=t and a count in [1,t]");
		}
		if(seedOverlaySet || seedRemoveSet){
			if(!seedOverlaySet || !seedRemoveSet || seedOverlayPath==null || seedOverlayPath.isEmpty() ||
				seedRemove==null || seedRemove.isEmpty() || seedPath==null){
				throw new IllegalArgumentException("seedoverlay and seedremove must be supplied together with seeds");
			}
			removedSeedIds=parseRemovedSeeds(seedRemove);
		}
		if(seedLimit<0 || seedLimit>0 && !hasSeeds()){
			throw new IllegalArgumentException("seedlimit must be nonnegative; a positive limit requires seeds");
		}
		if(seedExcludeSet){
			if(!hasSeeds() || seedLimit<=0 || seedExclude==null || seedExclude.isEmpty()){
				throw new IllegalArgumentException("seedexclude requires seeds, positive seedlimit and a nonempty ID list");
			}
			excludedSeedIds=parseFamilyIds(seedExclude,"seedexclude");
		}
		if(selectedSeeds!=null && !hasSeeds()){throw new IllegalArgumentException("selectedseeds requires seeds");}
		if(affectedFamiliesSet || affectedCountsSet){
			if(!affectedFamiliesSet || !affectedCountsSet || affectedFamilies==null || affectedFamilies.isEmpty() ||
				affectedCounts==null || affectedCounts.isEmpty()){
				throw new IllegalArgumentException("affectedfamilies and affectedcounts must be supplied together");
			}
			affectedIds=parseFamilyIds(affectedFamilies,"affectedfamilies");
			java.util.Arrays.sort(affectedIds);
			affectedOffsets=new IntHashMap(affectedIds.length);
			for(int i=0;i<affectedIds.length;i++){affectedOffsets.put(affectedIds[i],i);}
		}
		if((affectedCountsSet || skipFinalRebuild) &&
			(!hasSeeds() || merge || minCluster!=1 || topN!=0)){
			throw new IllegalArgumentException("affectedcounts/skipfinalrebuild require seeds, merge=f, mincluster=1 and topn=0");
		}
		if(affectedCountsSet && seedLimit!=0){
			throw new IllegalArgumentException("affectedcounts requires seedlimit=0");
		}
		if(skipFinalRebuild && (passes!=1 || quality!=null || consensus!=null || ddlout!=null && !preserveAllSeeds)){
			throw new IllegalArgumentException("skipfinalrebuild requires passes=1 and quality/consensus=null; ddlout requires preserveallseeds");
		}
		if(!Double.isFinite(stopChangeFraction) || stopChangeFraction<0 || stopChangeFraction>=1){
			throw new IllegalArgumentException("stopchangefraction must be finite and in [0,1): "+stopChangeFraction);
		}
		if(stopChangeFraction>0 && (!hasSeeds() || merge || topN!=0 || pruneFinal || skipFinalRebuild)){
			throw new IllegalArgumentException("stopchangefraction requires fixed seeds, ordinary rebuilds, merge=false, topn=0 and prunefinal=false");
		}
		if(passChanges!=null && (!hasSeeds() || merge || topN!=0)){
			throw new IllegalArgumentException("passchanges requires seeds, merge=f and topn=0");
		}
		if(pruneFinal && (!hasSeeds() || passes<1 || merge)){
			throw new IllegalArgumentException("prunefinal requires seeds, passes>=1 and merge=false");
		}
		if(preserveAllSeeds && (!hasSeeds() || passes!=1 || consensusPasses!=2 ||
			minCluster!=1 || topN!=0 || seedLimit!=0 || seedExcludeSet || seedOverlaySet ||
			merge || pruneFinal || quality!=null)){
			throw new IllegalArgumentException("preserveallseeds requires all fixed seeds, passes=1, consensuspasses=2, mincluster=1, topn=0, and no exclusion, overlay, merge, pruning or quality report");
		}
		if(input==null || parser.in2!=null || parser.out2!=null || parser.append){
			throw new IllegalArgumentException("Supply protein FASTA input or an ordered comma-separated list; paired files and append mode are unsupported");
		}
		inputFiles=ClusteringInput.inputFiles(input);
		if(tailCrop<0){throw new IllegalArgumentException("tailcrop must be nonnegative; zero preserves complete input proteins");}
		if(minLength<1 || minSpan<1 || minSpan>minLength || minCluster<1 || shortlist<1 ||
			passes<0 || passes==Integer.MAX_VALUE || consensusPasses<1 || topN<0 || checkpointEvery<1){
			throw new IllegalArgumentException("Invalid clustering lengths, cluster count, pass count or checkpoint interval");
		}
		if(output!=null && (output.indexOf('#')<0 || output.indexOf('#')!=output.lastIndexOf('#'))){
			throw new IllegalArgumentException("out= requires exactly one # placeholder for stable family ID");
		}
		if(mergeMaximumSize<0 || (mergeRoster!=null)!=(mergeMaximumSize>0)){
			throw new IllegalArgumentException("mergeroster and positive mergemaxsize require each other");
		}
		if(mergeRoster!=null && (passes!=0 || minCluster!=1 || topN!=0 || maxReads>=0 || seedPath!=null ||
			merge || scavenge || greedyRamp || pairhist!=null || decisionaudit!=null || greedyaudit!=null || assignmetrics!=null)){
			throw new IllegalArgumentException("mergeroster requires all consensus inputs, passes=0, mincluster=1, topn=0 and no seeds, merge, scavenge, ramp or decision observers");
		}
		if(assignmetrics!=null && passes<1){throw new IllegalArgumentException("assignmetrics requires passes>=1 to observe a frozen assignment pass");}
		if(merge && (passes<1 || !pinned() || pairhist!=null || decisionaudit!=null ||
			greedyaudit!=null || assignmetrics!=null || greedyRamp)){
			throw new IllegalArgumentException("merge requires passes>=1, a pinned alignment and no pairhist, audits, assignmetrics or ramp");
		}
		if(hasSeeds() && (passes<1 || concurrentGreedy || scavenge || greedyRamp ||
			pairhist!=null || decisionaudit!=null || greedyaudit!=null)){
			throw new IllegalArgumentException("seeds requires passes>=1 and no greedy mode, scavenging, ramp, pairhist or decision audit");
		}
		if(scavenge && (pairhist!=null || greedyaudit!=null)){
			throw new IllegalArgumentException("scavenge does not support pairhist or greedyaudit; use ordinary stage counters");
		}
		if((decisionaudit!=null || greedyaudit!=null)!=(auditqueries!=null) ||
			(decisionaudit!=null && (workerThreads!=1 || passes<1))){
			throw new IllegalArgumentException("An audit prefix and auditqueries require each other; decisionaudit also requires t=1 and passes>=1");
		}
		if((greedyaudit==null)!=(auditoldmembers==null) ||
			(greedyaudit!=null && (concurrentGreedy || width!=1024 || greedyHits>4 || shortlist>500))){
			throw new IllegalArgumentException("greedyaudit and auditoldmembers require each other, serial greedy, width=1024, shortlist<=500 and a greedy hit floor in0-4 or flex/ramp");
		}
		if(decisionaudit!=null && passHitsFlex){
			throw new IllegalArgumentException("The existing decision audit records one fixed hit floor; use numeric passhits with decisionaudit");
		}
		if(rampFamilies<0 || rampMinMembers<1 || (greedyRamp && rampFamilies==0)){
			throw new IllegalArgumentException("greedyhits=ramp requires explicit rampfamilies>0 and rampminmembers>=1");
		}
		if(!Double.isFinite(snapshotFraction) || snapshotFraction<=0 ||
			concurrentGreedy && (!indexed || !pinned() || decisionaudit!=null)){
			throw new IllegalArgumentException("snapshotfraction must be positive and finite; concurrentgreedy requires indexed=t, a pinned alignment and no serial decision audit");
		}
		// Gate construction validates finite numeric thresholds before any output is opened.
		final boolean defaultMatrixFile=matrixFile==null && ProteinSubstitutionMatrix.DEFAULT_NAME.equals(matrixName);
		if(defaultMatrixFile){matrixFile=ProteinSubstitutionMatrix.defaultResourcePath();}
		if(matrixFile!=null && !Files.isRegularFile(Paths.get(matrixFile))){
			throw new IllegalArgumentException("Matrix resource must be a regular disk file: "+matrixFile);
		}
		matrix=residentMatrix==null ? ProteinSubstitutionMatrix.load(matrixName,defaultMatrixFile ? null : matrixFile) : residentMatrix;
		if(!matrix.name.equals(matrixName)){
			throw new IllegalArgumentException("Resident clustering matrix differs from the configured selector: "+matrixName);
		}
		gate=makeGate();
		mergeGate=makeGate(mergeIdentity(), mergeRatio());
		if(explicitPositionMode!=null && positional!=null && positional.booleanValue()!=explicitPositionMode.positional()){
			throw new IllegalArgumentException("positional="+positional+" contradicts positionmode="+explicitPositionMode.flag);
		}
		positionMode=explicitPositionMode!=null ? explicitPositionMode :
			(positional==null ? (pinned() ? ProteinPositionMode.FIXED : ProteinPositionMode.PLAIN) :
				(positional.booleanValue() ? ProteinPositionMode.FIXED : ProteinPositionMode.PLAIN));
		Shared.AMINO_IN=true;
		geometry=ProteinKmerGeometry.make(k,aaBits);
		try{
			// Legacy policy constants describe12K6 only; the descriptor also routes16K5 raw filtering.
			final ProteinBlacklistPolicy.Descriptor policy=ProteinBlacklistPolicy.policy(geometry);
			if(blacklistK>=0 || blacklistMode!=null || blacklistProvenance!=null){
				blacklist=ProteinBlacklistPolicy.loadExplicit(geometry,blacklistPath,blacklistPin,blacklistWords,
					blacklistK,blacklistMode,blacklistProvenance,verification==VerificationMode.STRUCTURAL);
			}else if(policy!=null){
				if(blacklistPin!=null && !blacklistPin.equals(policy.sha80) ||
					blacklistWords>=0 && blacklistWords!=policy.words){
					throw new IllegalArgumentException("Clustering requires the "+policy.label+" blacklist");
				}
				final DDLPentamerBlacklist supplied;
				if(blacklistPath==null){blacklistPath=policy.resource().toString();supplied=null;}
				else{supplied=policy.load(Paths.get(blacklistPath), verification==VerificationMode.STRUCTURAL);}
				blacklist=ProteinBlacklistPolicy.resolve(geometry,supplied);
				blacklistPin=blacklist.sha80();blacklistWords=blacklist.size();
			}else if(blacklistPath==null){
				if(blacklistPin!=null || blacklistWords>=0){throw new IllegalArgumentException("blacklistsha80/blacklistwords require blacklist=");}
				blacklist=null;
			}else if(verification==VerificationMode.STRUCTURAL){
				blacklist=DDLPentamerBlacklist.loadStructural(Paths.get(blacklistPath), blacklistPin, blacklistWords);
			}else{
				blacklist=DDLPentamerBlacklist.load(Paths.get(blacklistPath),
					blacklistPin==null ? DigestSuffix.file(blacklistPath) : blacklistPin);
			}
		}catch(Exception failure){throw new IllegalArgumentException("Could not load the protein blacklist", failure);}
		if((residentSeeds!=null || seedSelection!=null) && (seedPath!=null || seedOverlaySet || maxReads>=0 || output!=null || mergeRoster!=null)){
			throw new IllegalArgumentException("Resident assignment requires one complete shared input, resident seeds and no file-seed overlay/cluster FASTA export");
		}
		final ClusterSketchCodec.Loaded originalSeeds=residentSeeds!=null ? residentSeeds : seedPath==null ? null : ClusterSketchCodec.read(seedPath);
		loadedSeedFamilies=originalSeeds==null ? 0 : originalSeeds.families.length;
		seedModels=seedOverlaySet ? ClusterSketchCodec.applyOverlay(originalSeeds,
			ClusterSketchCodec.read(seedOverlayPath),removedSeedIds) : originalSeeds;
		if(memberRanking && seedModels!=null && seedMembers==null){
			for(ClusterSketchCodec.Entry entry:seedModels.families){
				if(entry.memberRanking==null){throw new IllegalArgumentException("Ranked assignment requires seed member descriptors for family "+entry.familyId);}
			}
		}
		if(affectedIds!=null){
			final int[] found=new int[affectedIds.length];
			for(ClusterSketchCodec.Entry entry:seedModels.families){
				final int offset=affectedOffsets.get(entry.familyId);
				if(offset>=0){found[offset]++;}
			}
			for(int count:found){if(count!=1){throw new IllegalArgumentException("Affected family ID is absent or duplicated in the effective seed roster");}}
		}
		filter=new CandidateFilter(minimumHits,greedyHits,maximumHydroDelta,greedyHydroDelta,passHitsFlex);
		filter.validateWidth(width);
		if(!Float.isFinite(trimDepth) || trimDepth<0 || trimDepth>1){throw new IllegalArgumentException("trimdepth must be in [0,1]");}
		final ArrayList<String> fixed=new ArrayList<String>();
		fixed.add(members); fixed.add(ddlout); fixed.add(consensus); fixed.add(stats);
		fixed.add(refinement); fixed.add(checkpoints); fixed.add(hist); fixed.add(assignmetrics); fixed.add(passChanges); fixed.add(affectedCounts); fixed.add(selectedSeeds);
		if(quality!=null){
			fixed.add(quality+".cohesion.tsv");
			for(String path:ConsensusNeighborQualityReport.outputPaths(quality)){fixed.add(path);}
		}
		if(pairhist!=null){for(String path : AssignmentHistograms.outputPaths(pairhist)){fixed.add(path);}}
		if(decisionaudit!=null){for(String path : AssignmentDecisionAudit.outputPaths(decisionaudit)){fixed.add(path);}}
		if(greedyaudit!=null){for(String path : GreedyDecisionAudit.outputPaths(greedyaudit)){fixed.add(path);}}
		final String[] outputs=fixed.toArray(new String[fixed.size()]);
		if(!Tools.testOutputFiles(overwrite, false, false, outputs)){throw new IllegalArgumentException("Output files are unavailable");}
		final ArrayList<String> paths=new ArrayList<String>();
		for(String file:inputFiles){paths.add(file);}
		if(seedPath!=null){paths.add(seedPath);}
		if(seedMembers!=null){paths.add(seedMembers);}
		if(seedSelection!=null){for(String codec:seedSelection.codecs){paths.add(codec);}}
		if(seedOverlayPath!=null){paths.add(seedOverlayPath);}
		if(mergeRoster!=null){paths.add(mergeRoster);}
		if(blacklistPath!=null){paths.add(blacklistPath);}
		if(matrixFile!=null){paths.add(matrixFile);}
		if(auditqueries!=null){paths.add(auditqueries);}
		if(auditoldmembers!=null){paths.add(auditoldmembers);}
		for(String path : outputs){if(path!=null){paths.add(path);}}
		checkDistinct(paths);
		checkInputAliases(paths);
		checkQualityAliases(paths);
		checkPassChangesAliases(paths);
		checkAffectedCountAliases(paths);
		checkSelectedSeedAliases(paths);
		checkSeedAliases(paths);
		checkMatrixAliases(paths);
		checkMergeRosterAliases(paths);
		try{
			rawPolicyInputs=keyPolicy==DDLKeyTieBreak.Policy.RAW_TAIL_RARITY ?
				FamilyKeyPolicyInputs.declaredInputs(rawManifest==null ? null : Paths.get(rawManifest),keyInputs) : new Path[0];
		}catch(Exception failure){throw new IllegalArgumentException("Could not validate full-view key-policy inputs",failure);}
		checkRawPolicyOutputAliases(fixed);
		// Preserve original FASTA bytes; validation/encoding is performed by ProteinSequence.
		Shared.setThreads(1); Shared.AMINO_IN=true; Shared.TRIM_READ_DESCRIPTION=false;
		Read.TO_UPPER_CASE=Read.LOWER_CASE_TO_N=Read.DOT_DASH_X_TO_N=Read.FIX_HEADER=false;
		// Native Read validation must not fix/toss residues before ProteinSequence sees them.
		// ProteinSequence supplies strict alphabet/stop validation on the untouched bytes.
		Read.JUNK_MODE=Read.IGNORE_JUNK;
		ReadWrite.setZipThreads(1);
		if(!consensusBonusSpecified && residentKeyContext!=null){consensusMemberBonus=residentKeyContext.consensusMemberBonus;}
		try{
			final Path manifest=rawManifest==null ? null : Paths.get(rawManifest);
			if(residentKeyContext==null){
				keyContext=FamilyKeyPolicyInputs.load(keyPolicy,manifest,expectedRawTails,rawMinLength,rawMaxLength,rawTailCrop,
					width,geometry,sketchSeed,positionMode,partitions,overlap,blacklist,keyInputs,outstream,null,consensusMemberBonus);
			}else{
				FamilyKeyPolicyInputs.validateResident(residentKeyContext,keyPolicy,manifest,expectedRawTails,
					rawMinLength,rawMaxLength,rawTailCrop,width,geometry,sketchSeed,positionMode,partitions,overlap,blacklist,keyInputs,consensusMemberBonus);
				keyContext=residentKeyContext;
			}
		}catch(Exception failure){throw new IllegalArgumentException("Could not bind the full-view family key policy",failure);}
		// The study explicitly requests packed selection even for its smallest arm.
		// An ordinary implicit default retains the existing direct modal finalizer.
		sketches=new FamilySketchBuilder(width,geometry,sketchSeed,positionMode,partitions,overlap,blacklist,keyContext,
			keyPolicySpecified || keyPolicy!=DDLKeyTieBreak.Policy.SMALLEST || consensusMemberBonus>0);
		checkSeedGeometry(seedModels);
	}

	private boolean hasSeeds(){return seedPath!=null || residentSeeds!=null || seedSelection!=null;}

	private void checkSeedGeometry(final ClusterSketchCodec.Loaded models){
		if(models==null){return;}
		final ClusterSketchCodec.Header h=models.header;
		if(!sketches.matchesGeometry(h.width,h.geometry,h.sketchSeed,h.positionMode,h.partitions,h.overlap,h.blacklistSha80)){
			throw new IllegalArgumentException("Seed model geometry or blacklist differs from the configured query sketches");
		}
		if(!h.sameKeyPolicy(h.withKeyPolicy(keyContext))){
			throw new IllegalArgumentException("Seed model key policy, encoding or full-view count binding differs from the assignment");
		}
	}

	private TailAlignmentGate makeGate(){
		return makeGate(minIdentity, minRatio);
	}

	private TailAlignmentGate makeGate(double identity,double ratio){
		return new TailAlignmentGate(alignmentMode(), denominator, maxAlignLength, minSpan, identity, ratio,null,matrix,minPaired);
	}

	private static VerificationMode parseVerification(final String value){
		if("sha80".equalsIgnoreCase(value)){return VerificationMode.SHA80;}
		if("structural".equalsIgnoreCase(value)){return VerificationMode.STRUCTURAL;}
		throw new IllegalArgumentException("verification must be sha80 or structural");
	}

	private TailAlignmentGate.Mode alignmentMode(){
		if(explicitAlignmentMode!=null){
			if(alignmentAnchor!=AlignmentAnchor.RIGHT || !leftClip || denominator!=TailAlignmentGate.Denominator.ALIGNED_QUERY){
				throw new IllegalArgumentException("Explicit right-pinned alignmentmode requires pinright=t, trimleft=t and denominator=aligned");
			}
			return explicitAlignmentMode;
		}
		if(alignmentAnchor==AlignmentAnchor.LEFT){return TailAlignmentGate.Mode.PIN_LEFT_TRIM_RIGHT;}
		if(alignmentAnchor==AlignmentAnchor.RIGHT){
			return leftClip ? TailAlignmentGate.Mode.PIN_RIGHT_TRIM_LEFT :
				TailAlignmentGate.Mode.PIN_RIGHT_FULL_QUERY;
		}
		assert(alignmentAnchor==AlignmentAnchor.NONE) : "Unhandled anchor terminus: "+alignmentAnchor;
		return TailAlignmentGate.Mode.FULL_QUERY_UNPINNED;
	}
	private boolean pinRight(){return alignmentAnchor==AlignmentAnchor.RIGHT;}
	private boolean pinLeft(){return alignmentAnchor==AlignmentAnchor.LEFT;}
	private boolean pinned(){return alignmentAnchor!=AlignmentAnchor.NONE;}
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
	private static AlignmentAnchor resolveAlignmentAnchor(final AlignmentAnchor parsedAlignmentAnchor,
			final Boolean parsedPinRight,final Boolean parsedPinLeft){
		if(Boolean.TRUE.equals(parsedPinLeft) && parsedAlignmentAnchor!=AlignmentAnchor.LEFT){
			throw new IllegalArgumentException("pinleft=t contradicts anchorterminus="+parsedAlignmentAnchor);
		}
		if(Boolean.TRUE.equals(parsedPinRight) && parsedAlignmentAnchor!=AlignmentAnchor.RIGHT){
			throw new IllegalArgumentException("pinright=t contradicts anchorterminus="+parsedAlignmentAnchor);
		}
		if(Boolean.FALSE.equals(parsedPinLeft) && parsedAlignmentAnchor==AlignmentAnchor.LEFT){
			throw new IllegalArgumentException("pinleft=f contradicts anchorterminus=LEFT");
		}
		if(Boolean.FALSE.equals(parsedPinRight) && parsedAlignmentAnchor==AlignmentAnchor.RIGHT){
			throw new IllegalArgumentException("pinright=f contradicts anchorterminus=RIGHT");
		}
		return parsedAlignmentAnchor;
	}

	private double mergeIdentity(){return mergeMinIdentitySet ? mergeMinIdentity : minIdentity;}
	private double mergeRatio(){return mergeMinRatioSet ? mergeMinRatio : minRatio;}

	private void process(){process(null,null,1);}

	private Assignment process(final ClusteringInput residentInput,final MemberRankInputs borrowedInputs,final int layer){
		final long start=System.nanoTime();
		verifyKeyInputs();
		if(layer<0 || !memberRanking && borrowedInputs!=null ||
			memberRanking && residentSeeds!=null && (residentInput==null || borrowedInputs==null) ||
			borrowedInputs!=null && (residentInput==null || !borrowedInputs.belongsTo(residentInput.sequences))){
			throw new IllegalArgumentException("Resident ranking requires the identical original population/cache and an explicit nonnegative layer");
		}
		if(borrowedInputs!=null && borrowedInputs.hydroWeight!=memberHydroWeight){
			throw new IllegalArgumentException("Resident assignment cannot change memberhydroweight from "+borrowedInputs.hydroWeight+" to "+memberHydroWeight);
		}
		outstream.println(hasSeeds() ? "Fixed seeds skip greedy; all-input reassignment uses "+workerThreads+" worker(s)." : concurrentGreedy ? "Greedy admission uses "+(greedyThreadsSet ? greedyThreads : Math.min(workerThreads,16))+
			" workers and one background index builder; family results may vary with scheduling." :
			"Greedy clustering is serial; frozen reassignment and family rebuilds use "+workerThreads+" worker(s).");
		final ClusteringInput loaded=residentInput==null ? ClusteringInput.load(inputFiles,maxReads,output!=null,tailCrop) : residentInput;
		if(residentInput!=null && residentSeeds==null){throw new IllegalArgumentException("Borrowed input requires resident fixed-seed assignment");}
		ClusterSketchCodec.Loaded seedModels=this.seedModels;
		MemberRankInputs rankInputs=borrowedInputs;
		int memberLayer=layer;
		int nextFamilyId=-1;
		if(seedSelection!=null){
			final ClusterSeedSelector.RankedSelection selected;
			try(PhaseCostLog.Span phase=PhaseCostLog.begin("seed_selection")){
				selected=ClusterSeedSelector.selectRanked(seedSelection.codecs,seedSelection.budget,seedSelection.fraction,
					seedMembers==null ? null : Paths.get(seedMembers),seedMembersPin,loaded.sequences,0);
			}
			seedModels=selected.models;rankInputs=selected.inputs;nextFamilyId=selected.nextFamilyId;
			memberLayer=nextMemberLayer(seedModels.families,memberLayer);
			checkSeedGeometry(seedModels);
		}else if(memberRanking && !hasSeeds()){
			// Bind original ordinals once; greedy families publish only their admitted members.
			rankInputs=new MemberRankInputs(loaded.sequences,memberHydroWeight);
		}else if(memberRanking && residentSeeds==null && seedMembers==null){
			// Imported file seeds keep their original sparse IDs, order and immutable model payloads.
			memberLayer=nextMemberLayer(seedModels.families,memberLayer);
			rankInputs=new MemberRankInputs(loaded.sequences,memberHydroWeight);
		}else if(memberRanking && residentSeeds==null){
			final int n=seedModels.families.length;
			final int[] ids=new int[n],counts=new int[n];
			for(int i=0;i<n;i++){ids[i]=seedModels.families[i].familyId;counts[i]=seedModels.families[i].members;}
			final SeedMemberRankings hydration=SeedMemberRankings.read(Paths.get(seedMembers),seedMembersPin,loaded.sequences,ids,counts,0);
			rankInputs=hydration.inputs;
			final ClusterSketchCodec.Entry[] entries=new ClusterSketchCodec.Entry[n];
			for(int i=0;i<n;i++){entries[i]=ClusterSketchCodec.withMemberRanking(seedModels.families[i],hydration.family(ids[i]));}
			seedModels=ClusterSketchCodec.resident(seedModels.header,entries);
		}
		// Selection/hydration returns snapshots, not bound Family objects; share its
		// observations under the chosen coefficient before the engine attaches models.
		if(rankInputs!=null && borrowedInputs==null){rankInputs=rankInputs.withHydroWeight(memberHydroWeight);}
		final int loadedSeedFamilies=seedSelection==null ? this.loadedSeedFamilies : seedModels.families.length;
		final int[] mergeWeights=mergeRoster==null ? null : MergeRoster.load(mergeRoster,loaded.sequences);
		final AssignmentHistograms pairs=pairhist==null ? null : new AssignmentHistograms(passes,width,shortlist);
		final ConsensusQualityReport cohesion=quality==null ? null : new ConsensusQualityReport();
		final GreedyClusterer.Result result;
		try(ClusteringMetrics metrics=new ClusteringMetrics(checkpoints, hist, overwrite, shortlist, checkpointEvery, progress ? outstream : null);
			ConsensusRefiner refiner=new ConsensusRefiner(consensusPasses, consensusTrimLeft, trimDepth, refinement, overwrite);
			AssignmentDecisionAudit audit=decisionaudit==null ? null : new AssignmentDecisionAudit(
				decisionaudit,auditqueries,loaded.sequences,passes,shortlist,overwrite,gate,filter);
			GreedyDecisionAudit greedyAudit=greedyaudit==null ? null : new GreedyDecisionAudit(
				greedyaudit,auditoldmembers,auditqueries,loaded.sequences,shortlist,width,overwrite)){
			final GreedyClusterer engine=new GreedyClusterer(loaded.sequences, minLength, shortlist, minCluster,
				passes, topN, seed, sketches, gate, refiner, metrics, indexed,workerThreads,pairs,filter,audit,frozenIndexThreads);
			if(concurrentGreedy){engine.concurrentGreedy(snapshotFraction);}
			if(greedyThreadsSet){engine.greedyThreads(greedyThreads);}
			if(pruneFinal){engine.pruneFinal(true);}
			if(scavenge){engine.scavenge(true);}
			if(greedyRamp){engine.greedyRamp(rampFamilies,rampMinMembers);}
			if(inputOrder!=ClusteringOrder.Mode.SHUFFLE){engine.inputOrder(inputOrder);}
			if(greedyAudit!=null){engine.greedyAudit(greedyAudit);}
			if(assignmetrics!=null){engine.assignmentMetrics(new AssignmentMetricSummary(loaded.sequences.length));}
			if(cohesion!=null){engine.consensusQuality(cohesion);}
			if(memberRanking){engine.memberRanking(rankInputs,memberLayer);}
			if(seedModels!=null){engine.fixedSeeds(seedModels,seedLimit,excludedSeedIds,preserveAllSeeds);}
			if(skipFinalRebuild){engine.skipFinalRebuild(true);}
			if(passChanges!=null){engine.observePassChanges();}
			if(stopChangeFraction>0){engine.stopChangeFraction(stopChangeFraction);}
			if(merge){engine.merge(true,mergeGate);}
			if(mergeWeights!=null){engine.mergeSizes(mergeWeights,mergeMaximumSize);}
			if(selectedSeeds!=null){engine.writeSelectedSeeds(selectedSeeds,overwrite);}
			result=engine.run();
			verifyKeyInputs();
			checkClusterPaths(result);
			if(passChanges!=null){engine.writePassChanges(passChanges,overwrite);}
			if(affectedCounts!=null){writeAffectedCounts(result);}
			if(assignmetrics!=null){engine.writeAssignmentMetrics(assignmetrics,overwrite);}
			if(cohesion!=null){
				cohesion.write(quality+".cohesion.tsv",overwrite,result.families);
				ConsensusNeighborQualityReport.write(result.families,cohesion,sketches,gate,workerThreads,quality,overwrite);
			}
			ClusteringOutput.members(members, overwrite, loaded, result, minLength);
			ClusteringOutput.consensus(consensus, overwrite, result);
			ClusteringOutput.clusters(output, overwrite, loaded, result);
			if(ddlout!=null){ClusterSketchCodec.write(ddlout,overwrite,result.families,modelHeader(result.families.length));}
			String[] settings=configuration();
			if(minPaired>0){
				settings=java.util.Arrays.copyOf(settings,settings.length+1);
				settings[settings.length-1]="minpaired\t"+minPaired;
			}
			if(explicitAlignmentMode!=null){
				final int offset=settings.length;
				settings=java.util.Arrays.copyOf(settings,offset+1);
				settings[offset]="alignmentmode\t"+(explicitAlignmentMode==TailAlignmentGate.Mode.PIN_RIGHT_OVERHANG_ONLY ?
					"pin_right_overhang_only_aligned_query" : "pin_right_trim_left_aligned_query");
			}
			if(keyPolicySpecified || keyPolicy!=DDLKeyTieBreak.Policy.SMALLEST || consensusMemberBonus>0){
				final int offset=settings.length;
				settings=java.util.Arrays.copyOf(settings,offset+(keyPolicy==DDLKeyTieBreak.Policy.RAW_TAIL_RARITY ? 7 : 2));
				settings[offset]="keytiebreak\t"+keyPolicy.label;
				settings[offset+1]="keyencoding\t"+keyContext.encodingVersion;
				if(keyPolicy==DDLKeyTieBreak.Policy.RAW_TAIL_RARITY){
					settings[offset+2]="rawmanifest\t"+rawManifest;
					settings[offset+3]="expectedrawtails\t"+expectedRawTails;
					settings[offset+4]="rawminlen\t"+rawMinLength;
					settings[offset+5]="rawmaxlen\t"+rawMaxLength;
					settings[offset+6]="rawtailcrop\t"+rawTailCrop;
				}
			}
			if(consensusMemberBonus>0){
				final int offset=settings.length;settings=java.util.Arrays.copyOf(settings,offset+2);
				settings[offset]="consensusbonus\t"+consensusMemberBonus;
				settings[offset+1]="consensusbonusrule\tall-distinct-consensus-bin-key-v1";
			}
			if(!matrix.same(ProteinSubstitutionMatrix.BLOSUM62)){
				final int offset=settings.length;
				settings=java.util.Arrays.copyOf(settings,offset+4);
				settings[offset]="alignment_matrix\t"+matrix.name;
				settings[offset+1]="alignment_matrix_sha80\t"+matrix.sha80;
				settings[offset+2]="alignment_matrix_scale\t"+matrix.scale;
				settings[offset+3]="alignment_matrix_gap\t"+matrix.linearGap;
			}
			final int rankingOffset=settings.length;
			settings=java.util.Arrays.copyOf(settings,rankingOffset+(memberRanking ? 5 : 1));
			settings[rankingOffset]="memberranking\t"+memberRanking;
			if(memberRanking){
				settings[rankingOffset+4]="memberhydroweight\t"+memberHydroWeight;
				settings[rankingOffset+1]="member_ranking_layer\t"+memberLayer;
				settings[rankingOffset+2]="seedmembers_sha80\t"+(!hasSeeds() ? "none" : seedMembersPin==null ?
					residentSeeds!=null ? "resident" : "serialized" : seedMembersPin);
				settings[rankingOffset+3]="member_ranking_source\t"+(!hasSeeds() ? "raw" : residentSeeds!=null ? "resident" :
					seedMembersPin==null ? "serialized" : "file");
			}
			if(skipFinalRebuild){
				settings=java.util.Arrays.copyOf(settings,settings.length+1);
				settings[settings.length-1]="skipfinalrebuild\ttrue";
			}
			if(frozenIndexThreads>1){
				settings=java.util.Arrays.copyOf(settings,settings.length+1);
				settings[settings.length-1]="frozenindexthreads\t"+frozenIndexThreads;
			}
			if(mergeWeights!=null){
				final String[] extra=engine.mergeSizeStatistics();
				final int offset=settings.length;
				settings=java.util.Arrays.copyOf(settings,offset+2+extra.length);
				settings[offset]="mergemaxsize\t"+mergeMaximumSize;
				settings[offset+1]="merge_size_units\toriginal_tail_members";
				System.arraycopy(extra,0,settings,offset+2,extra.length);
			}
			if(pruneFinal){
				settings=java.util.Arrays.copyOf(settings,settings.length+1);
				settings[settings.length-1]="prunefinal\ttrue";
			}
			if(preserveAllSeeds){
				settings=java.util.Arrays.copyOf(settings,settings.length+1);
				settings[settings.length-1]="preserveallseeds\ttrue";
			}
			if(merge){
				final String[] extra=engine.mergeStatistics();
				final int offset=settings.length;
				settings=java.util.Arrays.copyOf(settings,offset+extra.length+2);
				settings[offset]="mergeminid\t"+mergeIdentity();
				settings[offset+1]="mergeminratio\t"+mergeRatio();
				System.arraycopy(extra,0,settings,offset+2,extra.length);
			}
			if(seedModels!=null){
				final int offset=settings.length;
				settings=java.util.Arrays.copyOf(settings,offset+2);
				settings[offset]="fixed_seeds\ttrue";
				settings[offset+1]="loaded_seed_families\t"+loadedSeedFamilies;
			}
			if(seedOverlaySet){
				final int offset=settings.length;
				settings=java.util.Arrays.copyOf(settings,offset+4);
				settings[offset]="seed_overlay_operation\t"+(removedSeedIds.length==1 ? "fission" : "fusion");
				settings[offset+1]="removed_seed_families\t"+removedSeedIds.length;
				settings[offset+2]="replacement_seed_families\t"+(removedSeedIds.length==1 ? 2 : 1);
				settings[offset+3]="effective_seed_families\t"+seedModels.families.length;
			}
			if(seedLimit>0){
				final int offset=settings.length;
				settings=java.util.Arrays.copyOf(settings,offset+2);
				settings[offset]="seedlimit\t"+seedLimit;
				settings[offset+1]="selected_seed_families\t"+Math.min(seedLimit,seedModels.families.length);
			}
			if(scavenge){
				final String[] extra=engine.scavengeStatistics();
				final int offset=settings.length;
				settings=java.util.Arrays.copyOf(settings,offset+extra.length);
				System.arraycopy(extra,0,settings,offset,extra.length);
			}
			if(inputOrder!=ClusteringOrder.Mode.SHUFFLE){
				final int offset=settings.length;
				settings=java.util.Arrays.copyOf(settings,offset+2);
				settings[offset]="order\t"+inputOrder.toString().toLowerCase(Locale.ROOT);
				settings[offset+1]="hydrobinwidth\t"+ClusteringOrder.BIN_WIDTH;
			}
			if(tailCrop>0){
				settings=java.util.Arrays.copyOf(settings,settings.length+1);
				settings[settings.length-1]="tailcrop\t"+tailCrop;
			}
			if(greedyRamp){
				final int offset=settings.length;
				settings=java.util.Arrays.copyOf(settings,offset+2);
				settings[offset]="rampfamilies\t"+rampFamilies;
				settings[offset+1]="rampminmembers\t"+rampMinMembers;
			}
			if(concurrentGreedy){
				final String[] extra=engine.concurrentStatistics();
				final int offset=settings.length;
				settings=java.util.Arrays.copyOf(settings,offset+extra.length);
				System.arraycopy(extra,0,settings,offset,extra.length);
			}
			ClusteringOutput.summary(stats, overwrite, result, seed, passes, (System.nanoTime()-start)*1e-9, metrics, refiner, settings);
			if(pairs!=null){pairs.write(pairhist,overwrite);}
			outstream.println("Assigned "+result.assigned+" of "+result.eligible+" eligible proteins to "+result.families.length+" families.");
			outstream.println("Alignments: "+metrics.assignmentAlignments()+" assignment, "+metrics.constructionAlignments()+" construction/validation.");
			outstream.println("Candidates screened before alignment: "+metrics.hitFiltered()+" by hits, "+metrics.hydroFiltered()+" by hydrophobicity.");
			outstream.println("Consensus proposals: "+refiner.accepted()+" accepted, "+refiner.rejected()+" rejected, "+refiner.unchanged()+" unchanged.");
		}
		verifyKeyInputs();
		return new Assignment(loaded,result,sketches,modelHeader(result.families.length),rankInputs,nextFamilyId,keyInputs);
	}

	/** Every exported or resident model carries the same policy as its actual builder. */
	private ClusterSketchCodec.Header modelHeader(final int families){
		return new ClusterSketchCodec.Header(geometry,width,sketchSeed,positionMode,partitions,overlap,families,
			blacklist==null ? null : blacklist.sha80()).withKeyPolicy(keyContext);
	}

	/** Detects changes to the original counting inputs before accepting dependent output. */
	private void verifyKeyInputs(){
		try{keyInputs.verify();}
		catch(IOException failure){throw new IllegalStateException("Full-view key-policy inputs became unavailable",failure);}
	}

	static final class Assignment{
		Assignment(final ClusteringInput input_,final GreedyClusterer.Result result_,final FamilySketchBuilder builder_,
				final ClusterSketchCodec.Header header_,final MemberRankInputs inputs_,final int nextFamilyId_,
				final StructuralInputs keyInputs_){
			assert(input_!=null && result_!=null && input_.sequences.length==result_.inputs()) :
				"Resident assignment returns the same original ordinal population and completed result";
			assert(inputs_==null || inputs_.belongsTo(input_.sequences)) : "Assignment must return its original cache by identity";
			input=input_;result=result_;builder=builder_;header=header_;memberInputs=inputs_;nextFamilyId=nextFamilyId_;keyInputs=keyInputs_;
		}
		final ClusteringInput input;
		final GreedyClusterer.Result result;
		final FamilySketchBuilder builder;
		final ClusterSketchCodec.Header header;
		final MemberRankInputs memberInputs;
		final int nextFamilyId;
		final StructuralInputs keyInputs;
	}

	/** Continue after imported or newly hydrated seed descriptors without wrapping. */
	static int nextMemberLayer(final ClusterSketchCodec.Entry[] families,final int floor){
		if(families==null){throw new IllegalArgumentException("Member-ranking layer needs resident entries");}
		if(floor<0){throw new IllegalArgumentException("Member-ranking layer floor must be nonnegative");}
		int layer=floor;
		for(ClusterSketchCodec.Entry family:families){
			if(family.memberRanking==null){throw new IllegalArgumentException("Ranked publication needs seed member descriptors");}
			if(family.memberRanking.layer==Integer.MAX_VALUE){
				throw new IllegalArgumentException("Member-ranking publication layer space is exhausted");
			}
			layer=Math.max(layer,family.memberRanking.layer+1);
		}
		return layer;
	}

	/** Ranked Cake either hydrates one pinned source or imports serialized descriptors from every selected source. */
	private static final class SeedSelection{
		SeedSelection(final String[] codecs_,final int budget_,final double fraction_){
			if(codecs_==null || codecs_.length<1){
				throw new IllegalArgumentException("Ranked initial selection requires at least one original codec");
			}
			codecs=codecs_.clone();
			for(String codec:codecs){
				if(codec==null || codec.isEmpty() || codec.indexOf(',')>=0 || !Files.isRegularFile(Paths.get(codec))){
					throw new IllegalArgumentException("Ranked initial selection requires existing explicit codec paths");
				}
			}
			ClusterSeedSelector.requestedCount(budget_,fraction_);
			budget=budget_;fraction=fraction_;
		}
		final String[] codecs;
		final int budget;
		final double fraction;
	}

	/** Explicit known seed IDs, including pruned empty references, get exact compact counts. */
	private void writeAffectedCounts(final GreedyClusterer.Result result){
		final long[] counts=new long[affectedIds.length];
		for(Family family:result.families){
			final int offset=affectedOffsets.get(family.familyId);
			if(offset>=0){counts[offset]=family.size();}
		}
		final ByteStreamWriter writer=new ByteStreamWriter(affectedCounts,overwrite,false,true);
		writer.start();
		try{
			writer.println("family_id\tassigned_members\tinput_records\teligible_records\ttotal_assigned\ttotal_unassigned\tfiltered_records");
			final ByteBuilder row=new ByteBuilder();
			for(int i=0;i<affectedIds.length;i++){
				row.clear().append(affectedIds[i]).tab().append(counts[i]).tab().append(result.inputs())
					.tab().append(result.eligible).tab().append(result.assigned).tab().append(result.inputs()-result.assigned)
					.tab().append(result.inputs()-result.eligible); writer.println(row);
			}
		}finally{if(writer.poisonAndWait()){throw new IllegalStateException("Affected-count output failed: "+affectedCounts);}}
	}

	/** Effective algorithm settings accompany results, independently of shell-history provenance. */
	private String[] configuration(){
		return new String[]{"minlen\t"+minLength, "minspan\t"+minSpan, "mincluster\t"+minCluster,
			"minhits\t"+(passHitsFlex ? "flex" : Integer.toString(minimumHits)), "greedyhits\t"+(greedyRamp ? "ramp" : greedyHits<0 ? "flex" : Integer.toString(greedyHits)),
			"maxhydrodelta\t"+maximumHydroDelta, "greedymaxhydrodelta\t"+greedyHydroDelta,
			"topn\t"+topN, "shortlist\t"+shortlist, "minid\t"+minIdentity, "minblosumratio\t"+minRatio,
			"stopchangefraction\t"+stopChangeFraction,
			"pinright\t"+pinRight(), "pinleft\t"+pinLeft(), "trimleft\t"+leftClip, "consensustrimleft\t"+consensusTrimLeft,
			"trimdepth\t"+trimDepth, "denominator\t"+denominator.name(), "maxalignlen\t"+maxAlignLength,
			"consensuspasses\t"+consensusPasses, "k\t"+k, "aabits\t"+(aaBits==null ? "none" : aaBits),
			"width\t"+width, "sketchseed\t"+sketchSeed,
			"positional\t"+positionMode.positional(), "positionmode\t"+positionMode.flag,
			"partitions\t"+partitions, "overlap\t"+overlap,
			"ranking\t"+(memberRanking ? "member_length_kd_hits_plus_monomer_cosine" : concurrentGreedy ? "snapshot_prefix_tail_plus_monomer_cosine" :
				indexed ? "exact_indexed_hits_plus_monomer_cosine" : "exhaustive_hits_plus_monomer_cosine"),
			"verification\t"+verification.name().toLowerCase(Locale.ROOT),
			"blacklist_sha80\t"+(blacklist==null ? "none" : blacklist.sha80()),
			"blacklist_words\t"+(blacklist==null ? 0 : blacklist.size()), "clustering_threads\t"+workerThreads};
	}

	/** Validates concrete family filenames against every fixed output and the input before writing them. */
	private void checkClusterPaths(GreedyClusterer.Result result){
		if(output==null){return;}
		final ArrayList<String> all=new ArrayList<String>();
		for(String file:inputFiles){all.add(file);}
		all.add(blacklistPath); all.add(members); all.add(ddlout); all.add(consensus); all.add(stats);
		all.add(matrixFile);
		all.add(seedPath); all.add(seedOverlayPath); all.add(seedMembers);
		if(seedSelection!=null){for(String codec:seedSelection.codecs){all.add(codec);}}
		all.add(mergeRoster);
		all.add(refinement); all.add(checkpoints); all.add(hist); all.add(assignmetrics); all.add(passChanges); all.add(affectedCounts); all.add(selectedSeeds);
		if(quality!=null){
			all.add(quality+".cohesion.tsv");
			for(String path:ConsensusNeighborQualityReport.outputPaths(quality)){all.add(path);}
		}
		all.add(auditqueries);
		all.add(auditoldmembers);
		if(pairhist!=null){for(String path : AssignmentHistograms.outputPaths(pairhist)){all.add(path);}}
		if(decisionaudit!=null){for(String path : AssignmentDecisionAudit.outputPaths(decisionaudit)){all.add(path);}}
		if(greedyaudit!=null){for(String path : GreedyDecisionAudit.outputPaths(greedyaudit)){all.add(path);}}
		final String[] clusters=new String[result.families.length];
		for(int i=0; i<clusters.length; i++){
			clusters[i]=output.replace("#", Integer.toString(result.families[i].familyId)); all.add(clusters[i]);
		}
		checkDistinct(all);
		checkInputAliases(all);
		checkQualityAliases(all);
		checkPassChangesAliases(all);
		checkAffectedCountAliases(all);
		checkSelectedSeedAliases(all);
		checkSeedAliases(all);
		checkMatrixAliases(all);
		checkMergeRosterAliases(all);
		checkRawPolicyOutputAliases(java.util.Arrays.asList(clusters));
		if(!Tools.testOutputFiles(overwrite, false, false, clusters)){
			throw new IllegalArgumentException("Cluster FASTA paths conflict with another input/output or existing file");
		}
	}

	/** Protects every ordered source, including hard links to a later file in the list. */
	private void checkInputAliases(ArrayList<String> paths){
		assert(inputFiles!=null && inputFiles.length>0 && paths!=null) :
			"Input parsing must establish the complete source list before output alias preflight";
		final Path[] sources=new Path[inputFiles.length];
		for(int i=0; i<sources.length; i++){sources[i]=Paths.get(inputFiles[i]);}
		final ArrayList<Path> targets=new ArrayList<Path>(paths.size());
		for(String name:paths){if(name!=null){targets.add(Paths.get(name));}}
		checkInputAliases(sources,targets);
	}

	/** Checks aliases without an all-inputs squared filesystem walk. This is a
	 * preflight over stable declared paths, not protection against concurrent replacement.
	 * BasicFileAttributes.fileKey guarantees equal nonnull keys for the same file;
	 * keys select candidates only, and isSameFile still confirms every rejection.
	 * Providers without keys retain the portable pairwise check. Missing outputs
	 * and native stdin names are skipped exactly as in the original preflight.
	 */
	public static void checkInputAliases(Path[] sources,Iterable<Path> targets){
		checkInputAliases(sources,targets,true);
	}

	/** Disabling keys exercises the conservative fallback without weakening alias checks. */
	public static void checkInputAliases(Path[] sources,Iterable<Path> targets,boolean useFileKeys){
		if(sources==null || sources.length==0 || targets==null){
			throw new IllegalArgumentException("Alias preflight requires sources and declared target paths");
		}
		try{
			final HashMap<Path,InputPathIdentity> cache=new HashMap<Path,InputPathIdentity>();
			final HashMap<Object,ArrayList<Path>> keyed=new HashMap<Object,ArrayList<Path>>();
			final ArrayList<Path> unkeyed=new ArrayList<Path>(),all=new ArrayList<Path>();
			for(Path target:targets){
				if(target==null){continue;}
				final InputPathIdentity identity=inputIdentity(target,cache);
				if(identity==null){continue;}
				all.add(target);
				if(!useFileKeys || identity.key==null){unkeyed.add(target);}
				else{
					ArrayList<Path> bucket=keyed.get(identity.key);
					if(bucket==null){bucket=new ArrayList<Path>();keyed.put(identity.key,bucket);}
					bucket.add(target);
				}
			}
			for(Path source:sources){
				if(source==null){throw new IllegalArgumentException("Null protein input in alias preflight");}
				final InputPathIdentity identity=inputIdentity(source,cache);
				if(identity==null){continue;}
				if(!useFileKeys || identity.key==null){rejectInputAliases(source,all);}
				else{
					final ArrayList<Path> candidates=keyed.get(identity.key);
					if(candidates!=null){rejectInputAliases(source,candidates);}
					rejectInputAliases(source,unkeyed);
				}
			}
		}catch(IOException failure){throw new IllegalArgumentException("Could not check protein input file identities",failure);}
	}

	/** Reads each declared path once; absence differs from an existing file with no key. */
	private static InputPathIdentity inputIdentity(Path path,HashMap<Path,InputPathIdentity> cache)throws IOException{
		assert(path!=null && cache!=null) : "Alias preflight must name an actual path before reading its identity";
		if(cache.containsKey(path)){return cache.get(path);}
		InputPathIdentity result=null;
		try{result=new InputPathIdentity(Files.readAttributes(path,BasicFileAttributes.class).fileKey());}
		catch(NoSuchFileException missing){/* New outputs and native stdin have no filesystem identity. */}
		cache.put(path,result);
		return result;
	}

	/** A candidate key is never sufficient by itself to reject a path. */
	private static void rejectInputAliases(Path source,ArrayList<Path> candidates)throws IOException{
		assert(source!=null && candidates!=null) : "Only existing input identities may reach alias confirmation";
		for(Path target:candidates){
			if(!source.equals(target) && Files.isSameFile(source,target)){
				throw new IllegalArgumentException("Protein input aliases another input or output: "+source+" and "+target);
			}
		}
	}

	private static final class InputPathIdentity{
		InputPathIdentity(Object key_){key=key_;}
		final Object key;
	}

	/** The full counting view may equal an ordinary input, but must never alias an output. */
	private void checkRawPolicyOutputAliases(final Iterable<String> outputs){
		try{
			for(Path source:rawPolicyInputs){
				for(String name:outputs){
					if(name==null){continue;}
					final Path target=Paths.get(name).toAbsolutePath().normalize();
					if(source.toAbsolutePath().normalize().equals(target) ||
						Files.exists(target) && Files.isSameFile(source,target)){
						throw new IllegalArgumentException("Full-view key-policy input aliases an output: "+source+" and "+name);
					}
				}
			}
		}catch(IOException failure){throw new IllegalArgumentException("Could not check full-view key-policy output aliases",failure);}
	}

	/** Matrix resources remain immutable even when an output uses a different hard-link path. */
	private void checkMatrixAliases(ArrayList<String> paths){
		if(matrixFile==null){return;}
		try{
			final Path source=Paths.get(matrixFile);
			if(!Files.isRegularFile(source)){throw new IllegalArgumentException("Matrix resource must be a regular disk file: "+matrixFile);}
			for(String other:paths){
				if(other==null || other.equals(matrixFile)){continue;}
				final Path candidate=Paths.get(other);
				if(Files.exists(candidate) && Files.isSameFile(source,candidate)){
					throw new IllegalArgumentException("Matrix resource aliases another input or output: "+other);
				}
			}
		}catch(IOException failure){throw new IllegalArgumentException("Could not check matrix resource file identity",failure);}
	}

	/** Original-tail counts are inputs and must survive even output hard-link aliases. */
	private void checkMergeRosterAliases(ArrayList<String> paths){
		if(mergeRoster==null){return;}
		try{
			final Path roster=Paths.get(mergeRoster);
			if(!Files.isRegularFile(roster)){throw new IllegalArgumentException("Merge roster must be an existing disk file: "+mergeRoster);}
			for(String other:paths){
				if(other==null || other.equals(mergeRoster)){continue;}
				final Path candidate=Paths.get(other);
				if(Files.exists(candidate) && Files.isSameFile(roster,candidate)){
					throw new IllegalArgumentException("Merge roster aliases another input or output: "+other);
				}
			}
		}catch(IOException failure){throw new IllegalArgumentException("Could not check merge roster file identity",failure);}
	}

	/** Structural edits still remove only one or two original models. */
	private static int[] parseRemovedSeeds(final String text){
		final int comma=text.indexOf(',');
		if(comma>=0 && text.indexOf(',',comma+1)>=0){throw new IllegalArgumentException("seedremove accepts one or two IDs");}
		return parseFamilyIds(text,"seedremove");
	}

	/** Parses a nonempty distinct-ID CSV, never KMG counts or empty list terms. */
	private static int[] parseFamilyIds(final String text,final String option){
		int count=1;
		for(int i=0;i<text.length();i++){if(text.charAt(i)==','){count++;}}
		final int[] ids=new int[count];
		final IntHashMap seen=new IntHashMap(count);
		for(int i=0,start=0;i<ids.length;i++){
			final int comma=text.indexOf(',',start);
			final int end=comma<0 ? text.length() : comma;
			if(start==end){throw new IllegalArgumentException("Empty "+option+" ID");}
			long value=0;
			for(int j=start;j<end;j++){
				final char c=text.charAt(j);
				if(c<'0' || c>'9'){throw new IllegalArgumentException(option+" IDs must be decimal integers");}
				value=value*10+c-'0';
				if(value>=Integer.MAX_VALUE){throw new IllegalArgumentException(option+" ID exceeds the supported range");}
			}
			final int id=(int)value;
			if(seen.get(id)>=0){throw new IllegalArgumentException("Duplicate "+option+" ID");}
			seen.put(id,i);ids[i]=id;start=end+1;
		}
		return ids;
	}

	/** Serialized seed inputs and original hydration assignments must survive output hard-link aliases. */
	private void checkSeedAliases(ArrayList<String> paths){
		final ArrayList<String> inputSeeds=new ArrayList<String>();
		inputSeeds.add(seedPath);inputSeeds.add(seedOverlayPath);inputSeeds.add(seedMembers);
		if(seedSelection!=null){for(String codec:seedSelection.codecs){inputSeeds.add(codec);}}
		for(String inputSeed:inputSeeds){
			if(inputSeed==null){continue;}
			try{
				final Path seed=Paths.get(inputSeed);
				for(String other:paths){
					if(other==null || other.equals(inputSeed)){continue;}
					final Path candidate=Paths.get(other);
					if(Files.exists(candidate) && Files.isSameFile(seed,candidate)){
						throw new IllegalArgumentException("Seed model paths alias an input or output: "+other);
					}
				}
			}catch(IOException failure){throw new IllegalArgumentException("Could not check seed model file identity",failure);}
		}
	}

	/** Selected-model metadata must not overwrite any input or output through a hard link. */
	private void checkSelectedSeedAliases(ArrayList<String> paths){
		if(selectedSeeds==null){return;}
		try{
			final Path target=Paths.get(selectedSeeds);
			if(!Files.exists(target)){return;}
			for(String other:paths){
				if(other==null || other.equals(selectedSeeds)){continue;}
				final Path candidate=Paths.get(other);
				if(Files.exists(candidate) && Files.isSameFile(target,candidate)){
					throw new IllegalArgumentException("Selected-seed output aliases an input or another output: "+other);
				}
			}
		}catch(IOException failure){throw new IllegalArgumentException("Could not check selected-seed file identity",failure);}
	}

	/** Compact count output has the same hard-link protection as the optional pass report. */
	private void checkAffectedCountAliases(ArrayList<String> paths){
		if(affectedCounts==null){return;}
		try{
			final Path target=Paths.get(affectedCounts);
			if(!Files.exists(target)){return;}
			for(String other:paths){
				if(other==null || other.equals(affectedCounts)){continue;}
				final Path candidate=Paths.get(other);
				if(Files.exists(candidate) && Files.isSameFile(target,candidate)){
					throw new IllegalArgumentException("Affected-count output aliases an input or another output: "+other);
				}
			}
		}catch(IOException failure){throw new IllegalArgumentException("Could not check affected-count file identity",failure);}
	}

	/** The optional observer must not overwrite an input or another output via a hard link. */
	private void checkPassChangesAliases(ArrayList<String> paths){
		if(passChanges==null){return;}
		try{
			final Path target=Paths.get(passChanges);
			if(!Files.exists(target)){return;}
			for(String other:paths){
				if(other==null || other.equals(passChanges)){continue;}
				final Path candidate=Paths.get(other);
				if(Files.exists(candidate) && Files.isSameFile(target,candidate)){
					throw new IllegalArgumentException("Pass-change output aliases an input or another output: "+other);
				}
			}
		}catch(IOException failure){throw new IllegalArgumentException("Could not check pass-change file identity",failure);}
	}

	/** Only new quality outputs need these bounded inode checks against existing inputs/outputs. */
	private void checkQualityAliases(ArrayList<String> paths){
		if(quality==null){return;}
		final ArrayList<String> qualityPaths=new ArrayList<String>();
		qualityPaths.add(quality+".cohesion.tsv");
		for(String path:ConsensusNeighborQualityReport.outputPaths(quality)){qualityPaths.add(path);}
		try{
			for(String name:qualityPaths){
				final Path target=Paths.get(name);
				if(!Files.exists(target)){continue;}
				for(String other:paths){
					if(other==null || name.equals(other)){continue;}
					final Path candidate=Paths.get(other);
					if(Files.exists(candidate) && Files.isSameFile(target,candidate)){
						throw new IllegalArgumentException("Quality output paths alias an input or another output: "+name+" and "+other);
					}
				}
			}
		}catch(IOException failure){throw new IllegalArgumentException("Could not resolve quality output file identities",failure);}
	}

	/** Detects relative/symlink aliases before any output can overwrite an input or another output. */
	private static void checkDistinct(ArrayList<String> paths){
		final HashSet<String> seen=new HashSet<String>();
		try{
			for(String path : paths){
				if(path==null){continue;}
				final String lower=path.toLowerCase(Locale.ROOT);
				final String key=lower.equals("stdin") || lower.startsWith("stdin.") ? "stream:stdin" :
					lower.equals("stdout") || lower.startsWith("stdout.") ? "stream:stdout" :
					lower.equals("stderr") || lower.startsWith("stderr.") ? "stream:stderr" : new File(path).getCanonicalPath();
				if(!seen.add(key)){throw new IllegalArgumentException("Input/output paths alias the same file or stream: "+path);}
			}
		}catch(IOException failure){throw new IllegalArgumentException("Could not resolve clustering input/output paths", failure);}
	}

	private final PrintStream outstream;
	private final String input, output;
	private final String[] inputFiles;
	private String members, ddlout, consensus, checkpoints, hist, pairhist, decisionaudit, auditqueries;
	private String greedyaudit,auditoldmembers,assignmetrics,quality;
	private String seedPath,seedOverlayPath,seedRemove,seedExclude,passChanges,selectedSeeds;
	private boolean seedOverlaySet,seedRemoveSet,seedExcludeSet;
	private String affectedFamilies,affectedCounts;
	private boolean affectedFamiliesSet,affectedCountsSet,skipFinalRebuild;
	private int[] affectedIds;
	private IntHashMap affectedOffsets;
	private int[] removedSeedIds,excludedSeedIds;
	private final int loadedSeedFamilies;
	private int seedLimit;
	private String mergeRoster;
	private int mergeMaximumSize;
	private final ClusterSketchCodec.Loaded seedModels;
	private final ClusterSketchCodec.Loaded residentSeeds;
	private final SeedSelection seedSelection;
	private boolean memberRanking=true;
	private double memberHydroWeight=MemberRankInputs.DEFAULT_HYDRO_WEIGHT;
	private String seedMembers,seedMembersPin;
	private String stats="ddlclust.stats.tsv", refinement="ddlclust.refinement.tsv";
	private final boolean overwrite;
	private final long maxReads;
	private final int workerThreads;
	private int greedyThreads;
	private int frozenIndexThreads=1;
	private boolean greedyThreadsSet,pruneFinal,preserveAllSeeds;
	private final TailAlignmentGate gate,mergeGate;
	private String matrixName=ProteinSubstitutionMatrix.DEFAULT_NAME,matrixFile;
	private final ProteinSubstitutionMatrix matrix;
	private DDLKeyTieBreak.Policy keyPolicy=DDLKeyTieBreak.Policy.SMALLEST;
	private boolean keyPolicySpecified;
	private long consensusMemberBonus;
	private boolean consensusBonusSpecified;
	private String rawManifest;
	private long expectedRawTails=-1;
	private int rawMinLength=1, rawMaxLength=Integer.MAX_VALUE, rawTailCrop=0;
	private final DDLKeyTieBreak.PolicyContext keyContext;
	private final StructuralInputs keyInputs;
	private final Path[] rawPolicyInputs;
	private final FamilySketchBuilder sketches;
	private final ProteinKmerGeometry geometry;
	private final CandidateFilter filter;
	private final DDLPentamerBlacklist blacklist;
	private final ProteinPositionMode positionMode;
	private String blacklistPath, blacklistPin, aaBits, blacklistMode, blacklistProvenance;
	private int blacklistK=-1;
	private int blacklistWords=-1;
	private VerificationMode verification=VerificationMode.SHA80;
	private int minLength=100, minSpan=70, minPaired=0, minCluster=5, shortlist=50, passes=3, consensusPasses=2, topN=80000;
	private int maxAlignLength=Integer.MAX_VALUE, k=5, width=1024, partitions=4, overlap=2;
	private int tailCrop=0;
	private long seed=1, sketchSeed=0, checkpointEvery=10000;
	private double minIdentity=0.35, minRatio=0.14;
	private double mergeMinIdentity,mergeMinRatio;
	private int minimumHits=4,greedyHits=-1;
	private boolean passHitsFlex;
	private boolean greedyRamp;
	private int rampFamilies,rampMinMembers=1;
	private ClusteringOrder.Mode inputOrder=ClusteringOrder.Mode.SHUFFLE;
	private double maximumHydroDelta=0.55,greedyHydroDelta=0.53;
	private double stopChangeFraction;
	private AlignmentAnchor alignmentAnchor=AlignmentAnchor.RIGHT;
	private TailAlignmentGate.Mode explicitAlignmentMode;
	private boolean leftClip=true, consensusTrimLeft=true;
	private boolean indexed=true, progress=true;
	private boolean concurrentGreedy=false;
	private boolean scavenge=false;
	private boolean merge=false;
	private boolean mergeMinIdentitySet,mergeMinRatioSet;
	private double snapshotFraction=0.02;
	private Boolean positional;
	private float trimDepth=0;
	private TailAlignmentGate.Denominator denominator=TailAlignmentGate.Denominator.ALIGNED_QUERY;

	private enum AlignmentAnchor{RIGHT, LEFT, NONE}
}
