package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicReference;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import ddl.DDLFamilyProfiles;
import ddl.DDLFamilySearch;
import ddl.DDLPentamerBlacklist;
import ddl.ProteinBlacklistPolicy;
import ddl.ProteinKmerGeometry;
import ddl.ProteinPositionMode;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;

/** Standalone calibrated DDL/cosine/identity/BLOSUM family cascade.
 * @author Yelan, Sayu
 */
public final class FamilyTailCascade{

	private FamilyTailCascade(){}

	/** Runtime floors only raise frozen gates. Directed rounding preserves strict float acceptance. */
	private static float acceptFloor(HashMap<String,String> options,String name,double fallback,double min,double max){
		final double value=options.containsKey(name) ? Double.parseDouble(required(options,name)) : fallback;
		if(!Double.isFinite(value) || value<min || value>max){
			throw new IllegalArgumentException(name+" must be finite in ["+min+","+max+"]");
		}
		return FamilyGateCodec.roundUpper(value);
	}

	private static FamilyAlignmentMetric.Mode parseAlignmentMode(String label){
		for(FamilyAlignmentMetric.Mode mode:FamilyAlignmentMetric.Mode.values()){
			if(mode.label.equals(label)){return mode;}
		}
		throw new IllegalArgumentException("Unknown family alignment mode: "+label);
	}

	public static void main(String[] args){
		try{
			final PreParser pp=new PreParser(args,FamilyTailCascade.class,false);
			try{run(pp.args);}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	public static void run(String[] args) throws Exception{
		final HashMap<String,String> options=new HashMap<String,String>();
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg,'=');
			if(pair.length!=2){throw new IllegalArgumentException("Require flag=value arguments");}
			final String key=pair[0].toLowerCase(Locale.ROOT);
			if(!OPTIONS.contains(key) || options.put(key,pair[1])!=null){
				throw new IllegalArgumentException("Unknown or duplicate option: "+key);
			}
		}
		MemberLengthRange.from(options);// Validate member provenance; it is not a query-length gate.
		FamilyBlacklistBinding.applyDefaults(options);
		final int top=integer(options,"top"), k=integer(options,"k"), width=integer(options,"width");
		final ProteinKmerGeometry geometry=geometry(options,k);
		final int ddlReject=options.containsKey("ddlreject") ? integer(options,"ddlreject") : 8;
		final int ddlAccept=options.containsKey("ddlaccept") ? integer(options,"ddlaccept") : 26;
		final int workers=threads(options);
		final long seed=number(options,"seed"), maxReads=number(options,"maxreads");
		final CosineMode cosineMode=cosineMode(options);
		final RankingFeatures features=rankingFeatures(options,cosineMode);
		final ScalarBounds scalarBounds=scalarBounds(options,features);
		final String mode=options.containsKey("mode") ? options.get("mode") : "search";
		if(!mode.equals("search") && !mode.equals("union") && !mode.equals("grid")){
			throw new IllegalArgumentException("mode must be search, union or grid");
		}
		final boolean union=mode.equals("union");
		final boolean grid=mode.equals("grid");
		final FamilyThresholdGrid.GridDefinition gridDefinition=grid ?
			FamilyThresholdGrid.gridDefinition(options.get("gridaxes")) : null;
		final boolean stageCpu=optionalBoolean(options,"stagecpu",false);
		if(stageCpu && (!mode.equals("search") || options.containsKey("rankout"))){
			throw new IllegalArgumentException("stagecpu requires ordinary search without rank diagnostics");
		}
		final boolean candidateOnly=optionalBoolean(options,"candidateonly",false);
		if(candidateOnly && (!union || options.containsKey("supplementalout") || options.containsKey("auditout"))){
			throw new IllegalArgumentException("candidateonly requires union mode and no supplementalout/auditout");
		}
		final float identityAcceptFloor=acceptFloor(options,"identityacceptfloor",30,30,100);
		final float blosumAcceptFloor=acceptFloor(options,"blosumacceptfloor",0,0,Float.MAX_VALUE);
		final boolean floorAudit=optionalBoolean(options,"flooraudit",false);
		if(union && (identityAcceptFloor!=30 || blosumAcceptFloor!=0 || floorAudit) ||
				floorAudit && !options.containsKey("auditout")){
			throw new IllegalArgumentException("Accept floors require search mode; flooraudit requires auditout");
		}
		final boolean scalarGates=optionalBoolean(options,"scalargates",!union && options.containsKey("familyreferences"));
		final double scalarZ=options.containsKey("scalarz") ?
			Double.parseDouble(required(options,"scalarz")) : scalarGates ? 3.8906 : 0;
		if(!Double.isFinite(scalarZ) || (scalarZ!=0 && scalarZ!=3.2905 && scalarZ!=3.8906) ||
				(scalarGates && (scalarZ==0 || !options.containsKey("familyreferences") || union))){
			throw new IllegalArgumentException("Scalar autofail requires familyreferences, search mode, "+
				"and scalarz=3.2905 or3.8906");
		}
		if(union && (cosineMode!=CosineMode.RANKING || features!=RankingFeatures.COMBINED25 ||
				!options.containsKey("familyreferences") || options.containsKey("rankout") ||
				maxReads>(candidateOnly ? 1000000 : 100000))){
			throw new IllegalArgumentException("Union requires familyreferences, combined25 ranking, no rankout; "+
				"maxreads<=100000 (candidateonly<=1000000)");
		}
		if(!union && options.containsKey("supplementalout")){
			throw new IllegalArgumentException("supplementalout requires mode=union");
		}
		if(top<1 || workers<1 || geometry.k()<4 || geometry.k()>6 || (width!=512 && width!=1024) ||
				seed!=0 || maxReads<1 || ddlReject<0 || ddlAccept<=ddlReject){
			throw new IllegalArgumentException("Require top>=1, k=4..6, width=512/1024, seed=0, maxreads>=1, "+
				"and 0<=ddlreject<ddlaccept");
		}
		boolean collection=false;
		for(String key : COLLECTION_OPTIONS){collection|=options.containsKey(key);}
		final boolean familyReferences=options.containsKey("familyreferences");
		final ProteinPositionMode positionMode=positionMode(options,familyReferences && !union);
		final boolean positional=positionMode.positional();
		final int partitions=options.containsKey("partitions") ? integer(options,"partitions") : 4;
		final int overlap=options.containsKey("overlap") ? integer(options,"overlap") : 8;
		final boolean replacement=positional || options.containsKey("ddlprofiles");
		if(!geometry.isExact() && !replacement){
			throw new IllegalArgumentException("Irregular amino-acid geometry requires replacement profiles");
		}
		final boolean indexEarlyReject=indexEarlyReject(options);
		final boolean indexDirect=optionalBoolean(options,"indexdirect",false);
		final boolean indexPruned=indexPruned(indexEarlyReject,indexDirect,options);
		if(grid && indexDirect){throw new IllegalArgumentException("Grid does not support direct-completion indexing");}
		if(indexPruned && (!replacement || options.containsKey("rankout") ||
				!(mode.equals("search") || grid && indexEarlyReject))){
			throw new IllegalArgumentException("Experimental index pruning requires ordinary replacement-profile search "+
				"or explicit early-rejection grid mode without rank diagnostics");
		}
		if(replacement && (!familyReferences || union) ||
				!positional && containsAny(options,new String[]{"partitions","overlap"}) ||
				options.containsKey("ddlprofiles")!=options.containsKey("ddlprofilessha80")){
			throw new IllegalArgumentException("DDL replacement requires packaged search and path+sha80; "+
				"partitions/overlap require positional=t");
		}
		if(positional){DDLFamilyProfiles.validatePositionalLayout(partitions,overlap,geometry);}
		if(options.containsKey("blacklist")!=options.containsKey("blacklistsha80") ||
				(options.containsKey("blacklist") && !replacement)){
			throw new IllegalArgumentException("Experimental pentamer blacklist requires path+sha80 "+
				"and matching profile replacement");
		}
		final Path ddlProfiles=replacement ? input(options,"ddlprofiles") : null;
		if(grid && (!familyReferences || union || top!=50 || cosineMode!=CosineMode.RANKING ||
				features!=RankingFeatures.COMBINED25 || !scalarGates || ddlReject!=gridDefinition.R[0] || ddlAccept!=1025 ||
				containsAny(options,new String[]{"rankout","auditout","supplementalout","flooraudit","candidateonly",
					"identityacceptfloor","blosumacceptfloor"}))){
			throw new IllegalArgumentException("Grid requires packaged combined25/scalar search, top50, ddlreject"+
				gridDefinition.R[0]+"/ddlaccept1025 and no per-query ordinary outputs");
		}
		if(!grid && containsAny(options,new String[]{"gridaxes","gridcohort","gridwitnessout",
				"gridoffset","gridwitnessselected"})){
			throw new IllegalArgumentException("grid options require mode=grid");
		}
		if(collection && familyReferences){
			throw new IllegalArgumentException("Family reference files and legacy collections are mutually exclusive");
		}
		final long loadStart=System.nanoTime();
		final Path models=familyReferences ? null : input(options,"models"), in=input(options,"in");
		final Path out=Paths.get(required(options,"out")).toAbsolutePath().normalize();
		if(Files.exists(out)){throw new IllegalArgumentException("Output exists: "+out);}
		final Path rankOut=options.containsKey("rankout") ? path(required(options,"rankout")) : null;
		final Path auditOut=options.containsKey("auditout") ? path(required(options,"auditout")) : null;
		final Path supplementalOut=union && !candidateOnly ? path(required(options,"supplementalout")) : null;
		if(supplementalOut!=null && (Files.exists(supplementalOut) || Files.isSymbolicLink(supplementalOut) ||
				supplementalOut.toFile().getCanonicalFile().equals(out.toFile().getCanonicalFile()) ||
				(auditOut!=null && supplementalOut.toFile().getCanonicalFile().equals(auditOut.toFile().getCanonicalFile())))){
			throw new IllegalArgumentException("Supplemental output exists or aliases another output: "+supplementalOut);
		}
		if(auditOut!=null && (Files.exists(auditOut) || Files.isSymbolicLink(auditOut) ||
				auditOut.toFile().getCanonicalFile().equals(out.toFile().getCanonicalFile()) ||
				(rankOut!=null && auditOut.toFile().getCanonicalFile().equals(rankOut.toFile().getCanonicalFile())))){
			throw new IllegalArgumentException("Audit output exists or aliases another output: "+auditOut);
		}
		if(rankOut!=null && (rankOut.toFile().getCanonicalFile().equals(out.toFile().getCanonicalFile()) ||
				Files.exists(rankOut) || Files.isSymbolicLink(rankOut))){
			throw new IllegalArgumentException("Rank output exists or aliases primary output: "+rankOut);
		}
		if((options.containsKey("truth") || options.containsKey("truthsha80")) && rankOut==null && !grid){
			throw new IllegalArgumentException("Truth input requires rankout");
		}
		if(options.containsKey("truthsha80") && !options.containsKey("truth")){
			throw new IllegalArgumentException("truthsha80 requires truth");
		}
		final Path truthPath=optionalInput(options,"truth");
		Shared.AMINO_IN=true; Shared.setThreads(workers); ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		final DDLPentamerBlacklist blacklist=FamilyBlacklistBinding.officialGeometry(options) ?
			FamilyBlacklistBinding.load(options,FamilyCalibration.VerificationMode.SHA80) : options.containsKey("blacklist") ?
			DDLPentamerBlacklist.load(input(options,"blacklist"),required(options,"blacklistsha80")) : null;
		final DDLFamilySearch.Index index;
		final FamilyReference[] families;
		final boolean compactFamilyReferences;
		long familyReferenceLoadNanos=0, familyReferenceIndexNanos=0;
		FamilyReferenceLoader.Timing referenceTiming=null;
		FamilyAlignmentMetric alignmentMetric=FamilyAlignmentMetric.LEGACY;
		if(familyReferences){
			if(options.containsKey("alignmentmode") || options.containsKey("alignmentminspan") ||
					options.containsKey("alignmentminpaired")){
				throw new IllegalArgumentException("alignmentmode/alignmentminspan are supported only by callable "+
					"classifier packages");
			}
			if(containsAny(options,SINGLE_SHARD_OPTIONS) || containsAny(options,COLLECTION_OPTIONS) ||
					containsAny(options,LEGACY_MODEL_OPTIONS)){
				throw new IllegalArgumentException("Legacy reference flags cannot be combined with familyreferences");
			}
			final Integer maxFamilies=optionalInteger(options,"maxfamilies");
			final Path familyReferencePath=input(options,"familyreferences");
			final long referenceLoadStart=System.nanoTime();
			final FamilyReferenceLoader.Loaded loaded=FamilyReferenceLoader.loadMetricBound(
				familyReferencePath.toString(),workers,maxFamilies==null ? 1000000 : maxFamilies.intValue());
			final long referenceLoadStop=System.nanoTime();
			referenceTiming=loaded.timing;
			if(loaded.header.k!=k || loaded.header.width!=width || loaded.header.seed!=seed){
				throw new IllegalArgumentException("Family reference header does not match k/width/seed");
			}
			if(scalarGates && Double.compare(loaded.header.scalarZ,scalarZ)!=0){
				throw new IllegalArgumentException("Scalar gate z differs from the pinned reference package");
			}
			alignmentMetric=loaded.header.alignmentMetric;
			if(!alignmentMetric.same(FamilyAlignmentMetric.LEGACY)){
				throw new IllegalArgumentException("Direct familyreferences CLI scoring is released only for legacy "+
					"alignment packages");
			}
			compactFamilyReferences=loaded.header.spectrumWidth==FamilyReferenceCodec.MONOMER_WIDTH;
			if(compactFamilyReferences){validateCompactReferenceUse(cosineMode,features,rankOut);}
			families=loaded.live;
			final long indexStart=System.nanoTime();
			index=replacement ?
				loadReplacementIndex(families,ddlProfiles,geometry,workers,positionMode,
					partitions,overlap,blacklist,options) :
				loadIndex(families,workers);
			final long indexStop=System.nanoTime();
			familyReferenceLoadNanos=referenceLoadStop-referenceLoadStart;
			familyReferenceIndexNanos=indexStop-indexStart;
		}else if(collection){
			final CollectionFiles files=loadCollection(options);
			index=loadIndex(files.references,k,seed,width,out.getParent());
			validateProfileCensus(files.references,files.census,width);
			if(index.families()!=files.census.nonempty){
				throw new IllegalArgumentException("Global DDL roster does not match live census families");
			}
			final Models fixed=loadModels(models,input(options,"consensus"),index,files.census,
				integer(options,"maxreferenceaa"));
			final Gates gates=files.a48 ? loadGatesA48(files.lowers,files.uppers,index.names,files.census) :
				loadGates(files.lowers,files.uppers,index.names,files.census);
			final LinkedSpectra spectra=loadSpectra(files.references,index.names);
			families=legacyFamilies(fixed.values,gates.values,spectra.values,
				metadata(files.census,index.names),width);
			compactFamilyReferences=false;
		}else{
			final Path references=input(options,"references"), lower=input(options,"lower"), upper=input(options,"upper");
			final ReferenceFiles referenceFiles=loadReferenceManifest(references);
			index=loadIndex(new ReferenceFiles[] {referenceFiles},k,seed,width,out.getParent());
			final Models fixed=loadModels(models,index,integer(options,"maxreferenceaa"));
			final Gates gates=loadGates(lower,upper,index.names);
			final LinkedSpectra spectra=loadSpectra(new ReferenceFiles[] {referenceFiles},index.names);
			families=legacyFamilies(fixed.values,gates.values,spectra.values,
				metadata(referenceFiles.members,index.names),width);
			compactFamilyReferences=false;
		}
		final HashMap<String,Truth> truth=rankOut==null ? null : loadTruth(truthPath,families);
		if(union && families.length==0){throw new IllegalArgumentException("Union requires at least one live family");}
		final SearchConfig config=new SearchConfig(index,families,top,ddlReject,ddlAccept,geometry,seed,
			width,cosineMode,features,truth,scalarBounds,union,scalarGates,identityAcceptFloor,
			blosumAcceptFloor,floorAudit,candidateOnly,positionMode,partitions,overlap,blacklist,
			compactFamilyReferences,null,alignmentMetric);
		if(grid){
			final FamilyThresholdGrid.Summary summary=FamilyThresholdGrid.run(config,in,out,truthPath,options,workers,maxReads);
			pinConfiguredInputs(options,collection,familyReferences,models,in);
			if(ddlProfiles!=null){pin(ddlProfiles,required(options,"ddlprofilessha80"));}
			if(blacklist!=null){pin(path(required(options,"blacklist")),blacklist.sha80());}
			if(truthPath!=null){pin(truthPath,required(options,"truthsha80"));}
			pin(path(required(options,"scalarbounds")),required(options,"scalarboundssha80"));
			System.out.println("FAMILY_THRESHOLD_GRID_PASS reads="+summary.reads+" tuples="+summary.tuples+
				" actual_alignments="+summary.alignments+" first_ordinal="+summary.offset+
				" end_ordinal="+(summary.offset+summary.reads));
			return;
		}
		final ByteStreamWriter writer=writer(out);
		final ByteStreamWriter rankWriter=rankOut==null ? null : writer(rankOut);
		final ByteStreamWriter auditWriter=auditOut==null ? null : writer(auditOut);
		final FileFormat format=FileFormat.testInput(in.toString(),FileFormat.FASTA,null,false,true);
		final Streamer input=StreamerFactory.makeStreamer(format,null,true,maxReads,false,true,1);
		final long loadStop=System.nanoTime(), queryStart=loadStop;
		final CascadeAccumulator accumulator=new CascadeAccumulator(union ? families.length : 0);
		boolean writerStarted=false, rankWriterStarted=false, auditWriterStarted=false, inputStarted=false;
		try{
			writer.forcePrint(candidateOnly ? CANDIDATE_HEADER : union ? UNION_PAIR_HEADER : OUTPUT_HEADER);
			writer.start();
			writerStarted=true;
			if(rankWriter!=null){
				rankWriter.forcePrint(scalarBounds==null ? RANK_HEADER : COMBINED_RANK_HEADER);
				rankWriter.start(); rankWriterStarted=true;
			}
			if(auditWriter!=null){
				auditWriter.forcePrint(floorAudit ? FLOOR_AUDIT_HEADER : AUDIT_HEADER);
				auditWriter.start(); auditWriterStarted=true;
			}
			final ArrayList<ProcessThread> processThreads=new ArrayList<ProcessThread>(workers);
			for(int i=0; i<workers; i++){
				processThreads.add(new ProcessThread(i,input,writer,rankWriter,auditWriter,config,
					accumulator.failure,stageCpu));
			}
			try{
				final boolean success=startInputAndWait(input,processThreads,accumulator);
				inputStarted=true;
				if(!success){
					throw new IllegalStateException("FamilyTailCascade worker failed", accumulator.failure.get());
				}
			}finally{
				if(inputStarted && ReadWrite.closeStream(input)){
					throw new IllegalStateException("Query FASTA read failed");
				}
			}
		}finally{
			boolean failed=false;
			try{if(writerStarted){failed=writer.poisonAndWait();}}
			finally{
				try{if(rankWriterStarted){failed=rankWriter.poisonAndWait() || failed;}}
				finally{if(auditWriterStarted){failed=auditWriter.poisonAndWait() || failed;}}
			}
			if(failed){throw new IllegalStateException("Cascade output writer failed");}
		}
		final long queryStop=System.nanoTime();
		final long repinStart=System.nanoTime();
		pinConfiguredInputs(options,collection,familyReferences,models,in);
		if(ddlProfiles!=null){pin(ddlProfiles,required(options,"ddlprofilessha80"));}
		if(blacklist!=null){pin(path(required(options,"blacklist")),blacklist.sha80());}
		if(truthPath!=null){pin(truthPath,required(options,"truthsha80"));}
		if(scalarBounds!=null){pin(path(required(options,"scalarbounds")),required(options,"scalarboundssha80"));}
		final long repinStop=System.nanoTime();
		accumulator.counters.validate();
		if(candidateOnly){
			final Counters c=accumulator.counters;
			if(c.candidateDp!=0 || c.candidateUnion!=0 || c.candidateRows>c.rankFamilies ||
				c.candidateRows<c.queryCandidates || c.candidateRows>(long)top*c.queryCandidates){
				throw new IllegalStateException("Candidate-only selection performed DP or violated top cap");
			}
			System.out.println("FAMILY_TAIL_CANDIDATES_PASS reads="+c.reads+" candidates="+c.candidateRows+
				" candidate_queries="+c.queryCandidates+" ddl_reject="+c.queryDdlReject+" ddl_accept="+c.queryDdlAccept+
				" ranked_window_pairs="+c.rankFamilies+" dp=0");
			System.out.println("FAMILY_TAIL_CANDIDATES_RECIPE scope=deployed_top_after_full_ddl_top_plus_ties_ranking"+
				" ddlreject="+ddlReject+" ddlaccept="+ddlAccept+" top="+top+" no_alignment=true no_scalar_autofail=true");
			return;
		}
		if(union){
			if(accumulator.counters.candidateUnion!=accumulator.counters.rankFamilies){
				throw new IllegalStateException("Union did not cover every ranked window pair");
			}
			final UnionStats stats=accumulator.counters.unionStats;
			stats.write(supplementalOut,families,accumulator.counters.reads,accumulator.counters.candidateUnion);
			System.out.println("FAMILY_TAIL_UNION_PASS reads="+accumulator.counters.reads+
				" pairs="+accumulator.counters.candidateUnion+" ddl_reject="+accumulator.counters.queryDdlReject+
				" ddl_accept="+accumulator.counters.queryDdlAccept+" reaching_dp="+accumulator.counters.queryReachingDp+
				" ranked_window_pairs="+accumulator.counters.rankFamilies);
			System.out.println("FAMILY_TAIL_UNION_RECIPE scope=full_ddl_top_plus_ties_window "+
				"raw_metrics=double baseline_merge=exporter"+
				" ddlreject="+ddlReject+" ddlaccept="+ddlAccept+" deployment_top="+top);
			System.out.println("FAMILY_TAIL_UNION_TIMING load_seconds="+seconds(loadStop-loadStart)+
				" query_writer_seconds="+seconds(queryStop-queryStart)+" repin_seconds="+seconds(repinStop-repinStart));
			return;
		}
		System.out.println("FAMILY_TAIL_CASCADE_PASS reads="+accumulator.counters.reads+
			" residues="+accumulator.counters.residues+
			" assigned="+accumulator.counters.assigned+" rejected="+accumulator.counters.rejected);
		if(replacement){
			System.out.println("FAMILY_TAIL_CASCADE_SKETCH_RECIPE sketch="+positionMode.sketch+
				" profiles_sha80="+required(options,"ddlprofilessha80")+" family_gates_and_spectra=unchanged");
		}
		if(positional){
			System.out.println("FAMILY_TAIL_CASCADE_POSITIONAL_LAYOUT partitions="+partitions+" overlap="+overlap+
				" tag_bits="+DDLFamilyProfiles.positionalTagBits(partitions)+(positionMode==ProteinPositionMode.RELATIVE ?
				" position_basis=residue-length" : " segment_words="+DDLFamilyProfiles.positionalSegmentWords(partitions)));
		}
		System.out.println("FAMILY_TAIL_CASCADE_QUERY_COUNTS ddl_reject="+accumulator.counters.queryDdlReject+
			" ddl_accept="+accumulator.counters.queryDdlAccept+" cosine_accept="+accumulator.counters.queryCosineAccept+
			" identity_accept="+accumulator.counters.queryIdentityAccept+
			" blosum_upper_accept="+accumulator.counters.queryBlosumUpperAccept+
			" blosum_lower_accept="+accumulator.counters.queryBlosumLowerAccept+
			" exhausted="+accumulator.counters.queryExhausted);
		System.out.println("FAMILY_TAIL_CASCADE_STAGE_COUNTS cosine_all_pruned="+
			accumulator.counters.queryCosineAllPruned+" reaching_dp="+accumulator.counters.queryReachingDp);
		System.out.println("FAMILY_TAIL_CASCADE_CANDIDATE_COUNTS visits="+accumulator.counters.candidateVisits()+
			" ddl_low="+accumulator.counters.candidateDdlLow+" ddl_high="+accumulator.counters.candidateDdlHigh+
			" cosine_upper="+accumulator.counters.candidateCosineUpper+
			" cosine_low="+accumulator.counters.candidateCosineLow+
			" cosine_survive="+accumulator.counters.candidateCosineSurvive+
			" dp="+accumulator.counters.candidateDp+
			" identity_upper="+accumulator.counters.candidateIdentityUpper+
			" identity_low="+accumulator.counters.candidateIdentityLow+
			" blosum_upper="+accumulator.counters.candidateBlosumUpper+
			" blosum_lower="+accumulator.counters.candidateBlosumLower+
			" blosum_reject="+accumulator.counters.candidateBlosumReject);
		if(accumulator.counters.candidateAlignmentUnavailable>0 || accumulator.counters.candidateAlignmentFailed>0){
			System.out.println("FAMILY_TAIL_CASCADE_ALIGNMENT_FAILURES"+
				" alignment_unavailable="+accumulator.counters.candidateAlignmentUnavailable+
				" alignment_failed="+accumulator.counters.candidateAlignmentFailed);
		}
		System.out.println("FAMILY_TAIL_CASCADE_WORK_COUNTS batches="+accumulator.counters.batches+
			" rank_families="+accumulator.counters.rankFamilies+
			" rank_cosines="+accumulator.counters.rankCosines+
			" rank_zero_norm="+accumulator.counters.rankZeroNorm+
			" rank_bound_skips="+accumulator.counters.rankBoundSkips+
			" gate_cosines="+accumulator.counters.gateCosines);
		System.out.println("FAMILY_TAIL_CASCADE_TIMING load_validation_seconds="+seconds(loadStop-loadStart)+
			" family_reference_load_seconds="+seconds(familyReferenceLoadNanos)+
			" family_reference_index_seconds="+seconds(familyReferenceIndexNanos)+
			" query_writer_seconds="+seconds(queryStop-queryStart)+
			" final_repin_seconds="+seconds(repinStop-repinStart)+
			" worker_process_wall_seconds="+seconds(accumulator.counters.processNanos));
		accumulator.counters.printAccounting(ddlReject);
		accumulator.counters.printScalarAccounting(scalarGates,scalarZ);
		if(rankWriter!=null){accumulator.counters.printTruth(scalarBounds!=null);}
		if(referenceTiming!=null){
			final FamilyReferenceLoader.Timing lt=referenceTiming;
			System.out.println("FAMILY_TAIL_CASCADE_LOADER_TIMING total_wall_seconds="+seconds(lt.totalWallNanos)+
				" open_wall_seconds="+seconds(lt.openWallNanos)+
				" worker_phase_wall_seconds="+seconds(lt.workerPhaseWallNanos)+
				" merge_wall_seconds="+seconds(lt.mergeWallNanos)+
				" census_wall_seconds="+seconds(lt.censusWallNanos)+
				" close_wall_seconds="+seconds(lt.closeWallNanos)+
				" read_batch_worker_seconds="+seconds(lt.readBatchWorkerNanos)+
				" parse_worker_seconds="+seconds(lt.parseWorkerNanos)+
				" parsed_records="+lt.parsedRecords+
				" metadata_gates_worker_seconds="+seconds(lt.metadataGatesWorkerNanos)+
				" consensus_ddl_worker_seconds="+seconds(lt.consensusDdlWorkerNanos)+
				" spectrum_worker_seconds="+seconds(lt.spectrumWorkerNanos)+
				" dto_norm_worker_seconds="+seconds(lt.dtoNormWorkerNanos));
		}
		if(stageCpu){accumulator.counters.printStageCpu();}
		final DDLFamilySearch.Index.BuildTiming bt=index.buildTiming;
		if(bt!=null){
			System.out.println("FAMILY_TAIL_CASCADE_INDEX_TIMING threads="+bt.threads+
				" setup_wall_seconds="+seconds(bt.setupWallNanos)+
				" count_prefix_wall_seconds="+seconds(bt.countPrefixWallNanos)+
				" allocation_wall_seconds="+seconds(bt.allocationWallNanos)+
				" fill_wall_seconds="+seconds(bt.fillWallNanos)+
				" count_worker_seconds="+seconds(bt.countWorkerNanos)+
				" prefix_worker_seconds="+seconds(bt.prefixWorkerNanos)+
				" cursor_worker_seconds="+seconds(bt.cursorWorkerNanos)+
				" fill_worker_seconds="+seconds(bt.fillWorkerNanos));
		}
	}

	/** Load the packaged cascade for in-process callers that already own protein extraction. */
	public static LoadedClassifier loadCallable(HashMap<String,String> options,int threads) throws Exception{
		return loadCallable(options, threads, null, null);
	}

	/** Package wrapper owns this explicit file-stability context through post-evaluation verify(). */
	public static LoadedClassifier loadCallableStructural(HashMap<String,String> options,int threads,
			DDLFamilySearch.StructuralInputs inputs)throws Exception{
		require(inputs!=null,"Structural file loading needs a retained input-stability context");
		return loadCallable(options,threads,null,null,FamilyCalibration.VerificationMode.STRUCTURAL,null,false,inputs);
	}

	/** Separate fragment package: reuses only validated reference payloads and the index.
	 * Historical reference gates/metric stay intact and are never fragment acceptance gates. */
	static FragmentFamilySearch loadFragments(final HashMap<String, String> options, final int threads)throws Exception{
		FragmentFamilySearch.validatePolicy(options);
		require(threads>0, "Fragment reference loading requires a positive worker count");
		final Path references=input(options, "familyreferences"), profiles=input(options, "ddlprofiles");
		pin(references, required(options, "familyreferencessha80")); pin(profiles, required(options, "ddlprofilessha80"));
		final FamilyReferenceLoader.Loaded loaded=FamilyReferenceLoader.loadMetricBound(references.toString(),
			threads,1000000);
		require(loaded.header.k==6 && loaded.header.width==512 && loaded.header.seed==0 &&
			loaded.header.spectrumWidth==FamilyReferenceCodec.MONOMER_WIDTH,
			"Fragment policy binds K6/512/seed0 references with stored monomer spectra");
		for(FamilyReference f:loaded.live){
			require(f.consensus.length()>0 && f.consensus.length()<=100,
				"Fragment reference coverage is limited to the declared retained100-aa cap: "+f.id);
		}
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(6, required(options, "aabits"));
		final DDLFamilySearch.Index index=DDLFamilySearch.Index.parallel(replacementNames(loaded.live),
			loadReplacementRows(loaded.live, profiles, geometry, ProteinPositionMode.FIXED, 4, 2, null, 512), null, threads);
		final SearchConfig config=new SearchConfig(index, loaded.live, 50, 5, 513, geometry, 0, 512,
			CosineMode.RANKING, RankingFeatures.MONOMER, null, null, false, false, 0, 0, false, false,
			ProteinPositionMode.FIXED, 4, 2, null, true, null, FamilyAlignmentMetric.LEGACY);
		pin(references, required(options, "familyreferencessha80")); pin(profiles, required(options, "ddlprofilessha80"));
		return new FragmentFamilySearch(config, options);
	}

	/** Preparation transfers immutable references and modal rows directly to calibration.
	 * File-backed and resident packages share policy, census and index validation below.
	 */
	static LoadedClassifier loadPrepared(HashMap<String,String> options, int threads,
			FamilyReferenceCodec.Header header, FamilyReference[] references,
			FamilyMemberPreparation.Prepared[] prepared) throws Exception{
		return loadPrepared(options,threads,header,references,prepared,FamilyCalibration.VerificationMode.SHA80);
	}

	/** Structural preparation keeps the validated in-memory roster and never enters a file-backed loader. */
	static LoadedClassifier loadPrepared(HashMap<String,String> options, int threads,
			FamilyReferenceCodec.Header header, FamilyReference[] references, FamilyMemberPreparation.Prepared[] prepared,
			FamilyCalibration.VerificationMode verificationMode) throws Exception{
		require(prepared!=null && references!=null && prepared.length==references.length,
			"Resident preparation requires one prepared payload per physical reference");
		return loadCallable(options, threads, FamilyReferenceLoader.fromMemoryMetricBound(header, references),
			prepared, verificationMode);
	}

	/** Transfer the exact immutable filter used to prepare members; never reload an unmeasured list. */
	static LoadedClassifier loadPrepared(HashMap<String,String> options, int threads,
			FamilyReferenceCodec.Header header, FamilyReference[] references,
			FamilyMemberPreparation.Prepared[] prepared,
			FamilyCalibration.VerificationMode verificationMode,DDLPentamerBlacklist blacklist) throws Exception{
		require(prepared!=null && references!=null && prepared.length==references.length,
			"Resident preparation requires one prepared payload per physical reference");
		return loadCallable(options,threads,FamilyReferenceLoader.fromMemoryMetricBound(header,references),
			prepared,verificationMode,blacklist,true);
	}

	/** A resident roster replaces only file loading; all deployed policy checks remain common. */
	private static LoadedClassifier loadCallable(HashMap<String,String> options, int threads,
			FamilyReferenceLoader.Loaded resident, FamilyMemberPreparation.Prepared[] prepared) throws Exception{
		return loadCallable(options,threads,resident,prepared,FamilyCalibration.VerificationMode.SHA80);
	}

	private static LoadedClassifier loadCallable(HashMap<String,String> options, int threads,
			FamilyReferenceLoader.Loaded resident, FamilyMemberPreparation.Prepared[] prepared,
			FamilyCalibration.VerificationMode verificationMode) throws Exception{
		return loadCallable(options,threads,resident,prepared,verificationMode,null,false);
	}

	private static LoadedClassifier loadCallable(HashMap<String,String> options, int threads,
			FamilyReferenceLoader.Loaded resident, FamilyMemberPreparation.Prepared[] prepared,
			FamilyCalibration.VerificationMode verificationMode,
			DDLPentamerBlacklist residentBlacklist,boolean directBlacklist) throws Exception{
		return loadCallable(options,threads,resident,prepared,verificationMode,residentBlacklist,directBlacklist,null);
	}

	private static LoadedClassifier loadCallable(HashMap<String,String> options, int threads,
			FamilyReferenceLoader.Loaded resident, FamilyMemberPreparation.Prepared[] prepared,
			FamilyCalibration.VerificationMode verificationMode,
			DDLPentamerBlacklist residentBlacklist,boolean directBlacklist,
			DDLFamilySearch.StructuralInputs structuralInputs) throws Exception{
		require(structuralInputs==null || verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL &&
			resident==null && prepared==null &&
			(!options.containsKey("blacklist") || FamilyBlacklistBinding.structuralFileGeometry(options)),
			"Structural file loading requires an explicit file-backed input context");
		require(verificationMode!=null && (verificationMode==FamilyCalibration.VerificationMode.SHA80 ||
			structuralInputs!=null || resident!=null && prepared!=null &&
			(!options.containsKey("blacklist") || directBlacklist)),
			"Structural loading needs a file-stability context or resident references/modal rows and a transferred filter");
		require(!directBlacklist || resident!=null && prepared!=null,
			"A directly transferred filter is valid only for resident preparation");
		requireOnly(options,CALLABLE_OPTIONS);
		if(containsAny(options,COMPACT_EMBEDDED_OPTIONS)){
			throw new IllegalArgumentException("Embedded ranking sidecars require compact family references");
		}
		FamilyBlacklistBinding.requireOfficialOptions(options);
		FamilyBlacklistBinding.validateOptions(options,verificationMode);
		final MemberLengthRange memberRange=MemberLengthRange.from(options);
		final int top=integer(options,"top"), k=integer(options,"k"), width=integer(options,"width");
		final ProteinKmerGeometry geometry=geometry(options,k);
		final int ddlReject=options.containsKey("ddlreject") ? integer(options,"ddlreject") : 8;
		final int ddlAccept=options.containsKey("ddlaccept") ? integer(options,"ddlaccept") : 26;
		final long seed=number(options,"seed");
		if(top<1 || threads<1 || geometry.k()<4 || geometry.k()>6 || (width!=512 && width!=1024) ||
				seed!=0 || ddlReject<0 || ddlAccept<=ddlReject){
			throw new IllegalArgumentException("Callable cascade requires top>=1, threads>=1, k=4..6, "+
				"width=512/1024, seed=0, and 0<=ddlreject<ddlaccept");
		}
		final CosineMode cosineMode=cosineMode(options);
		final RankingFeatures features=rankingFeatures(options,cosineMode);
		final ScalarBounds scalarBounds=scalarBounds(options,features,verificationMode);
		final boolean scalarGates=optionalBoolean(options,"scalargates",true);
		final double scalarZ=options.containsKey("scalarz") ?
			Double.parseDouble(required(options,"scalarz")) : scalarGates ? 3.8906 : 0;
		final boolean replacement=options.containsKey("ddlprofiles");
		if(!geometry.isExact() && !replacement){
			throw new IllegalArgumentException("Irregular amino-acid geometry requires replacement profiles");
		}
		if(indexPruned(options) && !replacement){
			throw new IllegalArgumentException("Experimental index pruning requires replacement profiles");
		}
		if(options.containsKey("ddlprofiles")!=options.containsKey("ddlprofilessha80")){
			throw new IllegalArgumentException("Replacement profiles require path and sha80");
		}
		final ProteinPositionMode positionMode=positionMode(options,replacement);
		final boolean positional=positionMode.positional();
		final int partitions=options.containsKey("partitions") ? integer(options,"partitions") : 4;
		final int overlap=options.containsKey("overlap") ? integer(options,"overlap") : 2;
		if(!replacement && positional || !positional && containsAny(options,new String[]{"partitions","overlap"})){
			throw new IllegalArgumentException("DDL replacement requires path+sha80; partitions/overlap require positional=t");
		}
		if(positional){DDLFamilyProfiles.validatePositionalLayout(partitions,overlap,geometry);}
		final DDLPentamerBlacklist blacklist=directBlacklist ?
			residentBlacklist : FamilyBlacklistBinding.load(options,verificationMode);
		if(structuralInputs!=null && blacklist!=null){input(options,"blacklist",structuralInputs);}
		FamilyBlacklistBinding.verifyLoaded(options,blacklist,verificationMode);
		FamilyBlacklistBinding.validateGeometry(blacklist,geometry);
		if(options.containsKey("alignmentminpaired") && !options.containsKey("alignmentmode")){
			throw new IllegalArgumentException("Paired minimum requires explicit alignmentmode and alignmentminspan");
		}
		if(options.containsKey("alignmentmode")!=options.containsKey("alignmentminspan")){
			throw new IllegalArgumentException("Alignment metric overrides require both alignmentmode and alignmentminspan");
		}
		if(blacklist!=null && !replacement){
			throw new IllegalArgumentException("Pentamer blacklist requires replacement profiles");
		}
		if(!Double.isFinite(scalarZ) || (scalarZ!=0 && scalarZ!=3.2905 && scalarZ!=3.8906) ||
				scalarGates && scalarZ==0){
			throw new IllegalArgumentException("Scalar autofail requires scalarz=3.2905 or3.8906");
		}
		final FamilyReferenceLoader.Loaded loaded=resident==null ? FamilyReferenceLoader.loadMetricBound(
			input(options,"familyreferences",structuralInputs).toString(),threads,1000000) : resident;
		if(options.containsKey(MemberLengthRange.MIN_KEY)){
			for(FamilyReference family:loaded.all){
				if(!memberRange.containsTotal(family.members,family.residues)){
					throw new IllegalArgumentException("Reference census does not fit declared member range "+
						memberRange+": "+family.id);
				}
			}
		}
		if(loaded.header.k!=k || loaded.header.width!=width || loaded.header.seed!=seed){
			throw new IllegalArgumentException("Family reference header does not match k/width/seed");
		}
		if(scalarGates && Double.compare(loaded.header.scalarZ,scalarZ)!=0){
			throw new IllegalArgumentException("Scalar gate z differs from the pinned reference package");
		}
		if(options.containsKey("matrixfile")){input(options,"matrixfile",structuralInputs);}
		final ProteinSubstitutionMatrix selectedMatrix=FamilyScoringPackage.resolveMatrix(options);
		if(options.containsKey("alignmentmode")){
			final int alignmentMinSpan=integer(options,"alignmentminspan");
			if(alignmentMinSpan>1000){throw new IllegalArgumentException("alignmentminspan must be at most1000");}
			final FamilyAlignmentMetric optionMetric=new FamilyAlignmentMetric(
				parseAlignmentMode(required(options,"alignmentmode")),alignmentMinSpan,selectedMatrix,
				options.containsKey("alignmentminpaired") ? integer(options,"alignmentminpaired") : 0);
			if(!optionMetric.same(loaded.header.alignmentMetric)){
				throw new IllegalArgumentException("Reference alignment metric differs from the requested callable package metric");
			}
		}else if(resident==null && !loaded.header.alignmentMetric.same(FamilyAlignmentMetric.LEGACY)){
			throw new IllegalArgumentException("Explicit new-metric references require alignmentmode and alignmentminspan");
		}
		if(!selectedMatrix.same(loaded.header.alignmentMetric.matrix)){
			throw new IllegalArgumentException("Reference matrix differs from the requested callable package scores");
		}
		final boolean compact=loaded.header.spectrumWidth==FamilyReferenceCodec.MONOMER_WIDTH;
		if(compact){validateCompactReferenceUse(cosineMode,features,null);}
		final Path ddlProfiles=replacement && resident==null ? input(options,"ddlprofiles",structuralInputs) : null;
		final DDLFamilySearch.Index index;
		if(resident!=null){
			require(replacement, "Resident preparation must retain its final replacement-profile binding");
			final char[][] rows=preparedRows(loaded, prepared, geometry, positionMode, partitions, overlap,
				blacklist,directBlacklist);
			index=replacementIndex(loaded.live, rows, threads, positionMode, partitions, options);
		}else{
			if(replacement && structuralInputs!=null){
				final char[][] rows=loadReplacementRows(loaded.live,ddlProfiles,geometry,positionMode,
					partitions,overlap,blacklist,width,structuralInputs);
				index=replacementIndex(loaded.live,rows,threads,positionMode,partitions,options);
			}else{
				index=replacement ? loadReplacementIndex(loaded.live,ddlProfiles,geometry,threads,
					positionMode,partitions,overlap,blacklist,options) :
					loadIndex(loaded.live,threads);
			}
		}
		final SearchConfig config=new SearchConfig(index,loaded.live,top,ddlReject,ddlAccept,geometry,seed,
			width,cosineMode,features,null,scalarBounds,false,scalarGates,
			acceptFloor(options,"identityacceptfloor",30,30,100),
			acceptFloor(options,"blosumacceptfloor",0,0,Float.MAX_VALUE),
			false,false,positionMode,partitions,overlap,blacklist,compact,null,
			loaded.header.alignmentMetric);
		return new LoadedClassifier(config,loaded.header,loaded.all);
	}

	/** Compact v1 family resources fold final positional rows and all live side tables together. */
	public static LoadedClassifier loadCompactCallable(HashMap<String,String> options,int threads,
			FamilyCompactLoader.Loaded loaded,DDLPentamerBlacklist blacklist) throws Exception{
		require(loaded!=null && loaded.header!=null && loaded.all!=null && loaded.live!=null,
			"Compact callable loading needs a fully parsed family resource");
		require(threads>0,"Compact callable loading needs positive worker count");
		requireOnly(options,CALLABLE_OPTIONS);
		final MemberLengthRange memberRange=MemberLengthRange.from(options);
		final int top=integer(options,"top"),k=integer(options,"k"),width=integer(options,"width");
		final ProteinKmerGeometry geometry=geometry(options,k);
		validateCompactBlacklist(options,blacklist,geometry);
		final int ddlReject=integer(options,"ddlreject"),ddlAccept=integer(options,"ddlaccept");
		final long seed=number(options,"seed");
		if(top<1 || k!=6 || width!=FamilyCompactCodec.DDL_WIDTH || seed!=0 ||
				ddlReject<0 || ddlAccept<=ddlReject || ddlAccept!=1025){
			throw new IllegalArgumentException(
				"Compact callable resources require top>=1, k=6, width=512, seed=0 and ddlaccept=1025");
		}
		final CosineMode cosineMode=cosineMode(options);
		final RankingFeatures features=rankingFeatures(options,cosineMode);
		validateCompactReferenceUse(cosineMode,features,null);
		final ScalarBounds scalarBounds=compactScalarBounds(loaded.header,options,features);
		final boolean scalarGates=optionalBoolean(options,"scalargates",true);
		final double scalarZ=options.containsKey("scalarz") ?
			Double.parseDouble(required(options,"scalarz")) : scalarGates ? 3.8906 : 0;
		final ProteinPositionMode positionMode=positionMode(options,true);
		final int partitions=options.containsKey("partitions") ? integer(options,"partitions") : 4;
		final int overlap=options.containsKey("overlap") ? integer(options,"overlap") : 2;
		if(!Double.isFinite(scalarZ) || (scalarZ!=3.2905 && scalarZ!=3.8906) ||
				positionMode!=ProteinPositionMode.FIXED || partitions!=4 || overlap!=2 ||
				options.containsKey("familyreferences") || options.containsKey("ddlprofiles")){
			throw new IllegalArgumentException(
				"Compact callable resources bind scalar gates, fixed P4/O2 profiles and embedded DDL rows");
		}
		require(loaded.all.length==compactInt(loaded.header,"families") &&
			loaded.live.length==compactInt(loaded.header,"live") &&
			k==compactInt(loaded.header,"k") && width==compactInt(loaded.header,"width") &&
			seed==compactLong(loaded.header,"seed") && ddlReject==compactInt(loaded.header,"ddlreject") &&
			ddlAccept==compactInt(loaded.header,"ddlaccept") &&
			scalarGates==compactBoolean(loaded.header,"scalargates") &&
			!compactBoolean(loaded.header,"indexearlyreject") &&
			"tail".equals(loaded.header.get("namespace")) &&
			"nonempty".equals(loaded.header.get("mode")) &&
			ProteinPositionMode.FIXED.sketch.equals(loaded.header.get("sketch")) &&
			geometry.name().equals(loaded.header.get("geometry")) &&
			"4".equals(loaded.header.get("partitions")) &&
			"2".equals(loaded.header.get("overlap")) &&
			required(options,"blacklistk").equals(loaded.header.get("blacklistk")) &&
			required(options,"blacklistmode").equals(loaded.header.get("blacklistmode")) &&
			required(options,"blacklistwords").equals(loaded.header.get("blacklistwords")) &&
			required(options,"blacklistprovenance").equals(loaded.header.get("blacklistprovenance")) &&
			"2".equals(loaded.header.get("tag_bits")) &&
			"49".equals(loaded.header.get("segment_words")) &&
			FamilyReferenceCodec.FORMULA_VERSION.equals(loaded.header.get("formulaVersion")) &&
			"combined25".equals(loaded.header.get("rankingfeatures")) &&
			"ranking".equals(loaded.header.get("cosine")),
			"Compact header does not match the callable resource policy");
		if(options.containsKey(MemberLengthRange.MIN_KEY)){
			for(FamilyReference family:loaded.all){
				if(!memberRange.containsTotal(family.members,family.residues)){
					throw new IllegalArgumentException("Compact census does not fit declared member range "+
						memberRange+": "+family.id);
				}
			}
		}
		final ProteinSubstitutionMatrix selectedMatrix=FamilyScoringPackage.resolveMatrix(options);
		require(selectedMatrix.name.equals(loaded.header.get("matrix")) &&
			selectedMatrix.binding().equals(loaded.header.get("matrixbinding")),
			"Compact header matrix binding differs from the callable package");
		final int alignmentMinSpan=integer(options,"alignmentminspan");
		if(alignmentMinSpan>1000){throw new IllegalArgumentException("alignmentminspan must be at most1000");}
		final FamilyAlignmentMetric metric=new FamilyAlignmentMetric(
			parseAlignmentMode(required(options,"alignmentmode")),alignmentMinSpan,selectedMatrix,
			options.containsKey("alignmentminpaired") ? integer(options,"alignmentminpaired") : 0);
		final DDLFamilySearch.Index index=replacementIndex(loaded.live,embeddedRows(loaded.live,width),
			threads,positionMode,partitions,options);
		final SearchConfig config=new SearchConfig(index,loaded.live,top,ddlReject,ddlAccept,geometry,seed,
			width,cosineMode,features,null,scalarBounds,false,scalarGates,
			acceptFloor(options,"identityacceptfloor",30,30,100),
			acceptFloor(options,"blosumacceptfloor",0,0,Float.MAX_VALUE),
			false,false,positionMode,partitions,overlap,blacklist,true,null,metric);
		return bindCompactSidecars(new LoadedClassifier(config,compactHeader(loaded,k,width,seed,scalarZ,metric),
			loaded.all),options,loaded);
	}

	/** Compact family resources fold the released live-order ranking side tables into each record. */
	private static LoadedClassifier bindCompactSidecars(LoadedClassifier base,HashMap<String,String> options,
			FamilyCompactLoader.Loaded loaded)throws Exception{
		final boolean heuristic=options.containsKey("heuristic") || options.containsKey("minqueryaa") ||
			options.containsKey("maxqueryaa");
		final boolean gc=options.containsKey("familygclambda") || options.containsKey("gcfilter") ||
			options.containsKey("gc_z");
		require((options.containsKey("heuristic") && options.containsKey("minqueryaa") &&
			options.containsKey("maxqueryaa"))==heuristic,
			"Compact heuristic ranking requires heuristic, minqueryaa and maxqueryaa together");
		require(!gc || heuristic,"Compact GC ranking requires the released heuristic pivot");
		if(!heuristic){return base;}
		require("hitfourth_h020_c150".equals(required(options,"heuristic")),
			"Compact callable resources bind the released common heuristic pivot");
		final LoadedClassifier ranked=base.withHeuristicRanking(new FamilyHeuristicScore(base,
			loaded.lengthMedians,loaded.kdMeans,integer(options,"minqueryaa"),
			integer(options,"maxqueryaa"),required(options,"heuristic")));
		if(!gc){return ranked;}
		require(options.containsKey("familygclambda") && options.containsKey("gcfilter") &&
			options.containsKey("gc_z"),
			"Compact GC ranking requires familygclambda, gcfilter and gc_z together");
		final FamilyGcPenalty penalty=FamilyGcInputs.fromCompact(loaded.live,loaded.liveStableIds,
			loaded.gcMeasured,loaded.gcMeans,loaded.gcSampleSds,Double.parseDouble(required(options,"familygclambda")),
			optionalBoolean(options,"gcfilter",true),Double.parseDouble(required(options,"gc_z")));
		return ranked.withPackagedFamilyGcPenalty(penalty);
	}

	/** Compact callable packages retain blacklist policy metadata but do not carry content digests. */
	private static void validateCompactBlacklist(HashMap<String,String> options,DDLPentamerBlacklist filter,
			ProteinKmerGeometry geometry){
		require(options.containsKey("blacklist") && filter!=null,"Compact callable resources require a loaded blacklist");
		for(String key:new String[]{"blacklistk","blacklistmode","blacklistwords","blacklistprovenance"}){
			require(options.containsKey(key),"Compact callable blacklist is missing "+key);
		}
		require(!options.containsKey("blacklistsha80"),
			"Compact callable packages bind blacklist provenance, not a stored blacklist sha80");
		try{filter.validate(geometry);}
		catch(IllegalArgumentException e){
			throw new IllegalArgumentException("Compact callable blacklist metadata differs from the loaded filter",e);
		}
		require(compactBlacklistGeometry(geometry) && integer(options,"blacklistk")==6 &&
			DDLPentamerBlacklist.MatchMode.parse(required(options,"blacklistmode"))==
				DDLPentamerBlacklist.MatchMode.REDUCED_LITERAL &&
			integer(options,"blacklistwords")==2048 &&
			filter.size()==2048 && filter.wordLength()==6 &&
			filter.matchMode()==DDLPentamerBlacklist.MatchMode.REDUCED_LITERAL &&
			required(options,"blacklistprovenance").equals(filter.provenance()),
			"Compact callable blacklist metadata differs from the loaded filter");
	}

	/** Compact packages currently accept only nonoverlapping K6 dimers with 120 ordered-pair classes. */
	private static boolean compactBlacklistGeometry(ProteinKmerGeometry geometry){
		final String prefix="dimer6+1,3,5+";
		final String name=geometry==null ? "" : geometry.name();
		return name.startsWith(prefix) && name.substring(prefix.length()).split("/",-1).length==120;
	}

	/** Compact packages pin the shared combined25 scalar normalization in the ten-line header. */
	private static ScalarBounds compactScalarBounds(FamilyCompactCodec.Header header,HashMap<String,String> options,
			RankingFeatures features){
		if(features!=RankingFeatures.COMBINED25){
			if(containsAny(options,new String[]{"scalarbounds","scalarboundssha80","scalarscale"})){
				throw new IllegalArgumentException("Scalar bounds/scale require rankingfeatures=combined25");
			}
			return null;
		}
		require(!options.containsKey("scalarbounds") && !options.containsKey("scalarboundssha80"),
			"Compact callable resources embed scalar bounds in the family header");
		final double[] min=new double[ConsensusScalars.SCALAR_COUNT],
			max=new double[ConsensusScalars.SCALAR_COUNT];
		for(int i=0;i<SCALAR_NAMES.length;i++){
			min[i]=Double.parseDouble(compactHeader(header,SCALAR_NAMES[i]+"_min"));
			max[i]=Double.parseDouble(compactHeader(header,SCALAR_NAMES[i]+"_max"));
		}
		return new ScalarBounds(min,max,Double.parseDouble(required(options,"scalarscale")));
	}

	private static FamilyReferenceCodec.Header compactHeader(FamilyCompactLoader.Loaded loaded,
			int k,int width,long seed,double scalarZ,FamilyAlignmentMetric metric){
		final String schema=metric.minimumPaired>0 ? FamilyReferenceCodec.PAIRED_MONOMER_SCHEMA :
			metric.matrix.same(ProteinSubstitutionMatrix.BLOSUM62) ?
			FamilyReferenceCodec.ALIGNMENT_MONOMER_SCHEMA : FamilyReferenceCodec.MATRIX_MONOMER_SCHEMA;
		return new FamilyReferenceCodec.Header(schema,k,width,seed,loaded.all.length,loaded.live.length,
			loaded.all.length-loaded.live.length,FamilyReferenceCodec.FORMULA_VERSION,
			FamilyReferenceCodec.SPECTRUM_ENCODING,scalarZ,metric);
	}

	private static char[][] embeddedRows(FamilyReference[] families,int width){
		require(families!=null,"Compact embedded rows require a live family roster");
		final char[][] rows=new char[families.length][];
		for(int i=0; i<families.length; i++){
			final FamilyReference family=families[i];
			require(family!=null && family.members>0 && family.ddl!=null && family.ddl.length==width,
				"Every compact live family must carry one embedded DDL row");
			rows[i]=family.ddl;
		}
		return rows;
	}

	private static int compactInt(FamilyCompactCodec.Header header,String key){
		return Math.toIntExact(compactLong(header,key));
	}

	private static long compactLong(FamilyCompactCodec.Header header,String key){
		final String value=compactHeader(header,key);
		long x=0;
		for(int i=0; i<value.length(); i++){
			final int digit=value.charAt(i)-'0';
			if(digit<0 || digit>9 || x>(Long.MAX_VALUE-digit)/10){
				throw new IllegalArgumentException("Compact header "+key+" must be a nonnegative integer");
			}
			x=x*10+digit;
		}
		return x;
	}

	private static boolean compactBoolean(FamilyCompactCodec.Header header,String key){
		final String value=compactHeader(header,key);
		if("true".equals(value)){return true;}
		if("false".equals(value)){return false;}
		throw new IllegalArgumentException("Compact header "+key+" must be true or false");
	}

	private static String compactHeader(FamilyCompactCodec.Header header,String key){
		final String value=header.get(key);
		if(value==null || value.length()==0){throw new IllegalArgumentException("Compact header missing "+key);}
		return value;
	}

	/** Match modal rows to the physical reference roster before making the live index view. */
	private static char[][] preparedRows(FamilyReferenceLoader.Loaded loaded, FamilyMemberPreparation.Prepared[] prepared,
			ProteinKmerGeometry geometry, ProteinPositionMode mode, int partitions, int overlap,
			DDLPentamerBlacklist blacklist,boolean directBlacklist){
		require(prepared!=null && prepared.length==loaded.all.length,
			"Prepared profiles must cover the complete physical roster");
		final char[][] rows=new char[loaded.live.length][];
		final String blacklistPin=blacklist==null ? "none" : blacklist.sha80();
		int live=0;
		for(int i=0; i<prepared.length; i++){
			final FamilyMemberPreparation.Prepared p=prepared[i];
			final FamilyReference f=loaded.all[i];
			require(p!=null && p.rank==f.rank && p.id.equals(f.id) && p.members==f.members &&
				p.residues==f.residues && p.validWords==f.validWords,
				"Prepared modal profile identity/census differs from its reference at ordinal "+i);
			p.verifyBlacklist(blacklist,directBlacklist);
			require(p.width==loaded.header.width && p.geometry.equals(geometry.name()) &&
				p.positionMode==mode && p.blacklistSha80.equals(blacklistPin) &&
				(!mode.positional() || p.partitions==partitions && p.overlap==overlap),
				"Prepared modal profile width or geometry differs from the deployed query policy: "+p.id);
			if(p.members>0){
				require(p.plain!=null && p.plain.length==p.width && p.positional!=null &&
					p.positional.length==p.width && live<rows.length,
					"Every live prepared family needs both complete modal rows at the declared width: "+p.id);
				rows[live++]=p.positional;
			}
		}
		require(live==rows.length, "Resident modal rows must match the validated live reference census");
		return rows;
	}

	private static ProteinKmerGeometry geometry(HashMap<String,String> options,int k){
		if(options.containsKey("aabits") && options.containsKey("aapartition")){
			throw new IllegalArgumentException("Use either aabits or aapartition, not both");
		}
		return ProteinKmerGeometry.make(k, options.containsKey("aabits") ?
			options.get("aabits") : options.get("aapartition"));
	}

	/** Preserve legacy defaults while rejecting an explicit boolean that contradicts the stored mode. */
	private static ProteinPositionMode positionMode(HashMap<String,String> options,boolean positionalDefault){
		return ProteinPositionMode.resolve(options.get("positionmode"),null,
			options.containsKey("positional") ? Boolean.valueOf(optionalBoolean(options,"positional",false)) : null,
			positionalDefault ? ProteinPositionMode.FIXED : ProteinPositionMode.PLAIN);
	}

	/** Evaluation-only interpolation/extrapolation of the original acceptance margin, in double precision. */
	public static double evaluationGate(float memberLower,float upper,float floor,double scale){
		require(Float.isFinite(memberLower) && Float.isFinite(upper) && Float.isFinite(floor) && upper>=memberLower,
			"Evaluation margins need finite gates with upper at least member lower");
		require(Double.isFinite(scale) && scale>=0,"Evaluation margin scale must be finite and nonnegative");
		final float current=Math.max(upper,floor);
		if(scale==1){return current;}
		final double base=Math.max(memberLower,floor),gate=base+scale*((double)current-base);
		require(Double.isFinite(gate),"Scaled evaluation threshold overflowed; use a smaller margin scale");
		return gate;
	}

	public static double memberMultiplierGate(float memberLower,float floor,double multiplier){
		require(Float.isFinite(memberLower) && Float.isFinite(floor),
			"Member-multiplier gates need finite stored lowers and runtime floors");
		require(Double.isFinite(multiplier) && multiplier>=0,"Member cutoff multiplier must be finite and nonnegative");
		final double member=memberLower*multiplier,gate=Math.max(floor,member);
		require(Double.isFinite(gate),"Member cutoff multiplier overflowed; use a smaller multiplier");
		return gate;
	}

	/** Lower-bound override is a search view only; the frozen family gates remain immutable. */
	private static float effectiveIdentityLower(FamilyReference family,float floor,float reduction){
		final float lower=Math.min(family.identityLower,Math.max(family.identityLower-reduction,floor));
		assert(Float.isFinite(lower) && lower>=0 && lower<=100) : "Identity lower reduction must clamp to a valid percent";
		return lower;
	}

	/** Immutable alternate thresholds; reference records and all ranking data remain unchanged. */
	private static final class EvaluationMargins{
		EvaluationMargins(SearchConfig config,double scale_){
			this(config,scale_,false);
		}
		EvaluationMargins(SearchConfig config,double value,boolean memberMultiplier_){
			scale=memberMultiplier_ ? 1 : value; memberMultiplier=memberMultiplier_ ? value : Double.NaN;
			memberMultiplierMode=memberMultiplier_;
			identity=new double[config.families.length];blosum=new double[identity.length];
			for(int i=0;i<identity.length;i++){
				final FamilyReference f=config.families[i];
				if(f.members==0){zeroMemberFamilies++;}
				if(f.alignmentCalibration!=FamilyReference.AlignmentCalibration.ENABLED){alignmentUnavailableFamilies++;}
				identity[i]=memberMultiplier_ ? memberMultiplierGate(f.identityLower,config.identityFloor,value) :
					evaluationGate(f.identityLower,f.identityUpper,config.identityFloor,scale);
				blosum[i]=memberMultiplier_ ? memberMultiplierGate(f.blosumLower,config.blosumFloor,value) :
					evaluationGate(f.blosumLower,f.blosumFinal,config.blosumFloor,scale);
				if(memberMultiplier_){
					identityFloorBinding+=binding(f.identityLower,config.identityFloor,value);
					blosumFloorBinding+=binding(f.blosumLower,config.blosumFloor,value);
				}
			}
		}
		private static int binding(float lower,float floor,double multiplier){
			final double member=lower*multiplier;
			return Double.compare(floor,member)>=0 ? 1 : 0;
		}
		final double scale,memberMultiplier;
		final double[] identity,blosum;
		final boolean memberMultiplierMode;
		int identityFloorBinding,blosumFloorBinding,zeroMemberFamilies,alignmentUnavailableFamilies;
	}

	public static final class LoadedClassifier{
		LoadedClassifier(SearchConfig config_,FamilyReferenceCodec.Header header_,FamilyReference[] all_){
			require(config_!=null && header_!=null && all_!=null,"Callable classifier requires a loaded reference package");
			require(config_.width==header_.width && config_.k==header_.k && config_.seed==header_.seed,
				"Callable query geometry must agree with the reference header before workers allocate their sketches");
			require(all_.length==header_.totalFamilies && config_.families.length==header_.liveFamilies,
				"Callable live and physical family counts must match the reference header");
			config=config_; header=header_; all=all_;
			liveToPhysical=liveToPhysical(config.families,all);
		}
		/** Immutable loaded package; callers may create one worker per consuming thread. */
		public Worker newWorker(){return new Worker(config);}
		/** Caller-only evidence on the unchanged pinned alignment path. */
		public Worker newWorker(FamilyMatchInformation.Binding information){
			if(information!=null){requireMatchInformation();}
			return new Worker(config,false,information);
		}
		public void requireMatchInformation(){
			require(hasAlignmentOnlyAcceptance() && !config.alignmentMetric.same(FamilyAlignmentMetric.LEGACY) &&
				!config.union && !config.candidateOnly,
				"Information bonuses require pinned alignment-only classification; shortcuts and legacy/provider "+
					"metrics cannot supply exact residue evidence");
		}
		/** Opt-in timed worker for benchmark diagnostics; scientific scoring stays identical. */
		public Worker newWorker(boolean stageCpu){return new Worker(config,stageCpu);}
		public int familyCount(){return config.families.length;}
		/** Read-only statistics from this classifier's actual live-family index. */
		public DDLFamilySearch.Index.PostingStatistics postingStatistics(){return config.index.postingStatistics();}
		/** Roster metadata uses live-family ordinals; original ranks may contain empty-model gaps. */
		public int familyRankAt(int ordinal){return config.families[ordinal].rank;}
		public String familyIdAt(int ordinal){return config.families[ordinal].id;}
		public long familyMembersAt(int ordinal){return config.families[ordinal].members;}
		/** Immutable consensus sequence for diagnostics that remeasure the deployed ranked window. */
		public ProteinSequence familyConsensusAt(int ordinal){return config.families[ordinal].consensus;}
		/** True only when identity/BLOSUM gates are calibrated for this live family. */
		public boolean familyAlignmentEnabledAt(int ordinal){
			return config.families[ordinal].alignmentCalibration==FamilyReference.AlignmentCalibration.ENABLED;
		}
		/** Effective strict identity-acceptance gate after global floors and margin scaling. */
		public double familyIdentityGateAt(int ordinal){
			return config.margins==null ? Math.max(config.families[ordinal].identityUpper,config.identityFloor) :
				config.margins.identity[ordinal];
		}
		/** Effective strict BLOSUM-acceptance gate after global floors and margin scaling. */
		public double familyBlosumGateAt(int ordinal){
			return config.margins==null ? Math.max(config.families[ordinal].blosumFinal,config.blosumFloor) :
				config.margins.blosum[ordinal];
		}
		/** Effective inclusive lower gate used before the BLOSUM route. */
		public float familyEffectiveIdentityLowerAt(int ordinal){
			return effectiveIdentityLower(config.families[ordinal],config.identityFloor,config.identityLowerReduction);
		}
		/** Pinned alignment mode and span of the loaded package. */
		public String alignmentModeName(){return config.alignmentMetric.name();}
		public int alignmentMinimumSpan(){return config.alignmentMetric.minimumSpan;}
		public int alignmentMinimumPaired(){return config.alignmentMetric.minimumPaired;}
		/** Effective score identity used by every worker; diagnostic swaps retain separate calibrated gates. */
		public ProteinSubstitutionMatrix alignmentMatrix(){return config.alignmentMetric.matrix;}
		/** Original score identity retained by the immutable reference/calibration header. */
		public ProteinSubstitutionMatrix calibrationMatrix(){return header.alignmentMetric.matrix;}
		/** Explicit diagnostic swap after ordinary package validation; never rewrites empirical gates.
		 * Fresh workers rebuild their scorer and selected-matrix normalization scratch.
		 * Only explicitly released band/blend names are admitted; the original
		 * matrix is also allowed as an exact-parity control. No fitted gate is changed.
		 */
		public LoadedClassifier withUncalibratedMatrix(final ProteinSubstitutionMatrix matrix){
			require(matrix!=null && "pairwise40".equals(header.alignmentMetric.matrix.name),
				"Uncalibrated matrix comparison requires the original pairwise40 reference binding");
			require(config.alignmentMetric.same(header.alignmentMetric),
				"Apply a diagnostic matrix once to the original validated alignment metric");
			require(matrix.same(header.alignmentMetric.matrix) || "pairwise30to50".equals(matrix.name) ||
				"blend25".equals(matrix.name) || "blend50".equals(matrix.name) || "blend75".equals(matrix.name) ||
				"pairwise25to40".equals(matrix.name) || "pairwise30to40".equals(matrix.name) ||
				"pairwise35".equals(matrix.name) || "pairwise25to45".equals(matrix.name),
				"Only released uncalibrated band/blend matrices or exact original-matrix parity are supported");
			require(matrix.scale==header.alignmentMetric.matrix.scale &&
				matrix.linearGap==header.alignmentMetric.matrix.linearGap,
				"Uncalibrated swap must preserve the frozen caller formula's score and gap units");
			require(hasAlignmentOnlyAcceptance() &&
				config.alignmentMetric.mode==FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY &&
				config.lengthPolicy==null && config.hitFloors==null && config.familyAdmission==null &&
				!config.lengthScaledShortAdmission &&
				!config.union && !config.candidateOnly,
				"Matrix diagnostic preserves fixed admission and pinned alignment; fitted/union/length-policy "+
					"variants are outside this comparison");
			if(matrix.same(config.alignmentMetric.matrix)){return this;}
			final FamilyAlignmentMetric metric=new FamilyAlignmentMetric(config.alignmentMetric.mode,
				config.alignmentMetric.minimumSpan,matrix,config.alignmentMetric.minimumPaired);
			return new LoadedClassifier(new SearchConfig(config,metric),header,all);
		}
		/** Exact amino-acid geometry and effective accepted margin loaded for this view. */
		public String geometryName(){return config.geometry.name();}
		public double evaluationMarginScale(){return config.margins==null ? 1 : config.margins.scale;}
		public boolean evaluationMemberMultiplierMode(){return config.margins!=null && config.margins.memberMultiplierMode;}
		/** Returns the active multiplier; callers must check the mode before requesting this quantity. */
		public double evaluationMemberMultiplier(){
			assert(evaluationMemberMultiplierMode()) :
				"Member multiplier is meaningless when its mode is off; check evaluationMemberMultiplierMode() before access";
			return evaluationMemberMultiplierMode() ? config.margins.memberMultiplier : Double.NaN;
		}
		public int memberMultiplierIdentityFloorBindings(){
			return evaluationMemberMultiplierMode() ? config.margins.identityFloorBinding : 0;
		}
		public int memberMultiplierBlosumFloorBindings(){
			return evaluationMemberMultiplierMode() ? config.margins.blosumFloorBinding : 0;
		}
		public int memberMultiplierZeroMemberFamilies(){
			return evaluationMemberMultiplierMode() ? config.margins.zeroMemberFamilies : 0;
		}
		public int memberMultiplierAlignmentUnavailableFamilies(){
			return evaluationMemberMultiplierMode() ? config.margins.alignmentUnavailableFamilies : 0;
		}
		/** Exact immutable raw-word filter used by query sketching; null means unfiltered. */
		public DDLPentamerBlacklist blacklist(){return config.blacklist;}
		/** Compare actual indexed registers without exposing mutable index storage. */
		public boolean matchesProfiles(char[][] rows){return config.index.matchesProfiles(rows);}
		/** Experimental profile/filter view with frozen references, alignment gates and
		 * ranking metadata. Rows must have been generated using this exact immutable
		 * raw-word filter and the parent's geometry; register bytes alone cannot prove
		 * that provenance. IDs/ranks bind every row to the live roster. Input rows are
		 * copied before parallel index construction and are not retained afterward. */
		public LoadedClassifier withReplacementProfiles(String[] ids,int[] ranks,char[][] rows,
				DDLPentamerBlacklist filter,int threads){
			require(ids!=null && ranks!=null && rows!=null && threads>0 &&
				ids.length==config.families.length && ranks.length==ids.length && rows.length==ids.length,
				"Replacement profiles require one named/ranked row per live family and positive workers");
			require(!config.index.usesPruning() && !config.index.usesDirectCompletion() &&
				config.lengthRanking==null && config.hitFloors==null && config.familyAdmission==null &&
				config.lengthPolicy==null && !config.lengthScaledShortAdmission && config.memberScalarMeans==null &&
				config.inverseRarityRanking==null && !config.union && config.truth==null && !config.candidateOnly,
				"Profile-only comparison preserves ordinary complete-index policy; fitted admission and diagnostic "+
					"policies need separate bindings");
			FamilyBlacklistBinding.validateGeometry(filter,config.geometry);
			final char[][] copy=new char[rows.length][];
			for(int i=0; i<rows.length; i++){
				require(config.families[i].id.equals(ids[i]) && config.families[i].rank==ranks[i] &&
					rows[i]!=null && rows[i].length==config.width,
					"Replacement row identity/rank/width differs at live ordinal "+i);
				copy[i]=rows[i].clone();
			}
			final DDLFamilySearch.Index replacement=DDLFamilySearch.Index.parallel(
				replacementNames(config.families),copy,null,threads);
			final SearchConfig derived=new SearchConfig(config,config.families,config.margins,null,config.lengthMode,
				null,config.heuristicRanking,false,null,null,replacement,filter);
			return new LoadedClassifier(derived,header,all);
		}
		/** Expose the resolved production policy so callers can record the configuration they actually loaded. */
		public boolean usesEarlyRejection(){return config.index.usesPruning() && !config.index.usesDirectCompletion();}
		/** Margin-based caller evidence requires every accepted family to have alignment measurements. */
		public boolean hasAlignmentOnlyAcceptance(){
			return config.ddlAccept>config.width && config.cosineMode!=CosineMode.FULL;
		}
		/** Absolute scale relative to the frozen parent gates, never compounded with a previous evaluation view. */
		public LoadedClassifier withEvaluationMarginScale(double scale){
			require(Double.isFinite(scale) && scale>=0,"Evaluation margin scale must be finite and nonnegative");
			if(scale==1 && config.margins==null){return this;}
			require(config.lengthPolicy==null,
				"Length policy is a terminal deployment view and cannot change its frozen effective gates");
			final EvaluationMargins margins=scale==1 ? null : new EvaluationMargins(config,scale);
			return new LoadedClassifier(new SearchConfig(config,config.families,margins),header,all);
		}
		/** Optional replay view using stored member-derived lowers scaled by a finite nonnegative multiplier. */
		public LoadedClassifier withMemberCutoffMultiplier(double multiplier){
			require(Double.isFinite(multiplier) && multiplier>=0,"Member cutoff multiplier must be finite and nonnegative");
			require(config.lengthPolicy==null,
				"Length policy is a terminal deployment view and cannot change its frozen effective gates");
			final EvaluationMargins margins=new EvaluationMargins(config,multiplier,true);
			return new LoadedClassifier(new SearchConfig(config,config.families,margins),header,all);
		}
		/** Absolute member-lower identity relaxation in percentage points; direct upper gates stay unchanged. */
		public LoadedClassifier withIdentityLowerReduction(float percentagePoints){
			require(Float.isFinite(percentagePoints) && percentagePoints>=0,
				"Identity lower reduction must be finite and nonnegative");
			if(percentagePoints==config.identityLowerReduction){return this;}
			require(config.lengthPolicy==null,
				"Length policy is a terminal deployment view and cannot change its frozen effective gates");
			return new LoadedClassifier(new SearchConfig(config,config.families,config.margins,config.lengthRanking,
				config.lengthMode,config.memberScalarMeans,config.heuristicRanking,config.lengthScaledShortAdmission,
				config.familyAdmission,config.hitFloors,config.index,config.blacklist,percentagePoints),header,all);
		}
		/** Explicit short-query experiment; the parent index, arrays and reference gates remain immutable.
		 * Query lengths are checked after terminal-stop normalization, before either admission decision. */
		public LoadedClassifier withLengthScaledShortAdmission(){
			if(config.lengthScaledShortAdmission){return this;}
			return new LoadedClassifier(new SearchConfig(config,config.families,config.margins,
				config.lengthRanking,config.lengthMode,config.memberScalarMeans,config.heuristicRanking,true),header,all);
		}
		public boolean usesLengthScaledShortAdmission(){return config.lengthScaledShortAdmission;}
		/** Uses retained family means before ranking; the parent and calibrated gates remain unchanged. */
		public LoadedClassifier withFamilyMeanAdmission(){
			if(config.familyAdmission!=null){return this;}
			return new LoadedClassifier(new SearchConfig(config,config.families,config.margins,
				config.lengthRanking,config.lengthMode,config.memberScalarMeans,config.heuristicRanking,
				config.lengthScaledShortAdmission,new FamilyLengthAdmission(config.families)),header,all);
		}
		public boolean usesFamilyMeanAdmission(){return config.familyAdmission!=null;}
		/** Initial fitting probe, bound to the exact live references and their frozen gates.
		 * This is not a selected policy; measure and freeze support before selecting a vector. */
		public FamilyHitFloorVector newFamilyHitFloorProbe(){
			return FamilyHitFloorVector.forClassifier(config);
		}
		/** Explicit family-specific raw admission. Attach ranking and margin views FIRST;
		 * later copies must preserve those choices and the exact reference/index binding. */
		public LoadedClassifier withFamilyHitFloors(final FamilyHitFloorVector floors){
			require(floors!=null && floors.boundToPolicy(config),
				"Family-specific floors must bind the exact loaded references, index, ranking, margins and gates");
			if(config.hitFloors==floors){return this;}
			return new LoadedClassifier(new SearchConfig(config,config.families,config.margins,
				config.lengthRanking,config.lengthMode,config.memberScalarMeans,config.heuristicRanking,
				config.lengthScaledShortAdmission,config.familyAdmission,floors),header,all);
		}
		public boolean usesFamilyHitFloors(){return config.hitFloors!=null;}
		/** Binds one accepted shared span-penalty coefficient without modifying serialized package bytes. */
		public LoadedClassifier withLengthPolicy(final FamilyLengthPolicy policy){
			require(policy!=null,"Length policy must be loaded and reference-validated by the caller resource");
			require(config.lengthPolicy==null,"Length policy is a terminal deployment view");
			require(config.alignmentMetric.mode==FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY,
				"Length policy requires the pinned right/trim-left aligned-query metric");
			return new LoadedClassifier(new SearchConfig(config,policy),header,all);
		}
		/** Uses an accepted policy certificate but applies a caller-supplied amplitude for tradeoff measurement. */
		public LoadedClassifier withExperimentalLengthAmplitude(final FamilyLengthPolicy policy,final float amplitude){
			require(policy!=null,"Length policy must be loaded and reference-validated by the caller resource");
			FamilyLengthAcceptance.validateAmplitude(amplitude);
			require(config.lengthPolicy==null,"Length policy is a terminal deployment view");
			require(config.alignmentMetric.mode==FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY,
				"Length policy requires the pinned right/trim-left aligned-query metric");
			return new LoadedClassifier(new SearchConfig(config,policy,amplitude),header,all);
		}
		public boolean usesLengthPolicy(){return config.lengthPolicy!=null;}
		/** Explicit opt-in; the index, references, calibrated metrics and gates stay shared. */
		public LoadedClassifier withLengthRanking(FamilyLengthScore binding,FamilyLengthScore.Mode mode){
			require(binding!=null && mode!=null && binding.boundTo(config),
				"Length model and medians must bind this exact loaded family/index roster");
			require(config.memberScalarMeans==null,
				"A learned ranking model fitted on consensus scalars cannot silently receive member-mean features");
			require(config.heuristicRanking==null && config.inverseRarityRanking==null,
				"Learned and explicit heuristic ranking policies cannot be combined silently");
			return new LoadedClassifier(new SearchConfig(config,binding,mode),header,all);
		}
		/** Explicit pre-top ranking over the existing DDL-eligible pool, without changing admission or calibrated gates. */
		public LoadedClassifier withHeuristicRanking(FamilyHeuristicScore binding)throws Exception{
			require(config.lengthPolicy==null,
				"Length policy is a terminal deployment view and cannot change its frozen ranking");
			require(binding!=null && binding.boundTo(config),"Heuristic inputs must bind this exact loaded family/index roster");
			require(config.lengthRanking==null && config.features==RankingFeatures.COMBINED25 &&
				config.cosineMode==CosineMode.RANKING && config.scalarBounds!=null && config.ddlAccept>config.width &&
				config.inverseRarityRanking==null && !config.union && config.truth==null,
				"Caller heuristics require combined25 ranking, alignment-only acceptance and no learned/union/truth mode");
			binding.verify();
			return new LoadedClassifier(new SearchConfig(config,config.families,config.margins,null,
				config.lengthMode,config.memberScalarMeans,binding),header,all);
		}
		public boolean usesHeuristicRanking(){return config.heuristicRanking!=null;}
		/** Attach immutable family means to the released pivot. Every worker call must
		 * then supply query GC availability explicitly, including neutral missing values.
		 * This transient API does not authorize exporting or recalibrating its binding. */
		public LoadedClassifier withFamilyGcPenalty(FamilyGcPenalty binding){
			return withFamilyGcPenalty(binding,false);
		}
		/** Only the package loader supplies the retained resource recipe needed for export. */
		LoadedClassifier withPackagedFamilyGcPenalty(FamilyGcPenalty binding){
			return withFamilyGcPenalty(binding,true);
		}
		private LoadedClassifier withFamilyGcPenalty(FamilyGcPenalty binding,boolean recipeBound){
			require(binding!=null && config.familyGcPenalty==null,"GC ranking needs one explicit immutable family binding");
			// Caller attachment preserves frozen evaluation margins; it does not refit or export gates.
			requireCalibrationUpdatesAllowed(true);
			require(config.heuristicRanking!=null && !config.alignmentMetric.same(FamilyAlignmentMetric.LEGACY),
				"GC ranking requires the released metric-bound common-pivot classifier");
			binding.requireRoster(this);
			return new LoadedClassifier(new SearchConfig(config,binding,recipeBound),header,all);
		}
		public boolean usesFamilyGcPenalty(){return config.familyGcPenalty!=null;}
		public void verifyHeuristicRanking()throws Exception{
			if(config.heuristicRanking!=null){config.heuristicRanking.verify();}
		}
		/** Explicit inverse-posting-rarity replacement for the raw-hit term in the bound heuristic pivot. */
		public LoadedClassifier withInverseRarityRanking(FamilyInverseRarityScore binding)throws Exception{
			require(config.lengthPolicy==null,
				"Length policy is a terminal deployment view and cannot change its frozen ranking");
			require(binding!=null && binding.boundTo(config),
				"Inverse-rarity inputs must bind this exact loaded family/index roster");
			require(config.lengthRanking==null && config.inverseRarityRanking==null &&
				config.features==RankingFeatures.COMBINED25 && config.cosineMode==CosineMode.RANKING &&
				config.scalarBounds!=null && config.ddlAccept>config.width && !config.union && config.truth==null,
				"Inverse rarity requires the released explicit heuristic pivot without learned/union/truth mode");
			binding.verify();
			return new LoadedClassifier(new SearchConfig(config,config.families,config.margins,null,
				config.lengthMode,config.memberScalarMeans,config.heuristicRanking,binding),header,all);
		}
		public boolean usesInverseRarityRanking(){return config.inverseRarityRanking!=null;}
		public void verifyInverseRarityRanking(){
			if(config.inverseRarityRanking!=null){config.inverseRarityRanking.verify();}
		}
		/** Explicit alternative ranking array; the original classifier and calibrated gates remain unchanged. */
		public LoadedClassifier withMemberScalarMeans(FamilyMemberScalarMeans binding){
			require(binding!=null && binding.boundTo(config),"Member means must bind this exact loaded family/index roster");
			require(config.scalarBounds!=null && config.features==RankingFeatures.COMBINED25 &&
				config.cosineMode==CosineMode.RANKING,
				"Member-mean arrays replace combined25 ranking coordinates only, never dimer-calibrated gates");
			require(config.lengthRanking==null && config.inverseRarityRanking==null,
				"Member-mean features require their own comparison before reuse by a fitted ranking model");
			binding.verify();
			return new LoadedClassifier(new SearchConfig(config,config.families,config.margins,null,
				config.lengthMode,binding),header,all);
		}
		/** Identify and recheck the optional external means before and after a reproducible evaluation. */
		public boolean usesMemberScalarMeans(){return config.memberScalarMeans!=null;}
		public void verifyMemberScalarMeans(){if(config.memberScalarMeans!=null){config.memberScalarMeans.verify();}}
		/** Creates a shared three-arm configuration; scales are absolute, never compounded. */
		public LengthPairedClassifier lengthPaired(FamilyLengthScore binding,double[] scales){
			require(config.alignmentMetric.same(FamilyAlignmentMetric.LEGACY),
				"Length-paired evidence currently reuses legacy glocal alignments only");
			return new LengthPairedClassifier(this,binding,scales);
		}
		FamilyCalibrationUpdates newCalibrationUpdates(){
			requireCalibrationUpdatesAllowed();
			return new FamilyCalibrationUpdates(config.families,config.ddlReject,config.ddlAccept,
				config.identityFloor,config.blosumFloor,config.alignmentMetric);
		}
		LoadedClassifier withCalibrationUpdates(FamilyCalibrationUpdates updates){
			requireCalibrationUpdatesAllowed();
			require(config.margins==null && updates!=null && updates.boundTo(config),
				"Calibration maxima must be bound to this classifier's live roster and deployed policy");
			final FamilyReference[] live=updates.updated();
			final FamilyReference[] physical=all.clone();
			require(live.length==config.families.length,"Calibration update changed the live roster width");
			for(int i=0;i<live.length;i++){
				final FamilyReference family=live[i], parent=config.families[i];
				require(family!=null && family.rank==parent.rank && family.id.equals(parent.id),
					"Calibration update changed the live roster order or family identity");
				physical[liveToPhysical[i]]=family;
			}
			return new LoadedClassifier(config.withFamilies(live),header,physical);
		}
		void writeReferences(Path out){
			writeReferences(out, Long.MAX_VALUE);
		}
		/** Preparation's decoded-output bound also covers the final, possibly updated references. */
		long writeReferences(Path out, long maxBytes){
			requireCalibrationUpdatesAllowed();
			require(config.familyAdmission==null,
				"Reference serialization alone does not retain the family-mean admission policy");
			require(!config.lengthScaledShortAdmission,
				"Reference serialization does not retain the length-scaled query admission policy");
			require(config.margins==null,"Evaluation-only thresholds are not a recalibrated reference package");
			require(config.memberScalarMeans==null,
				"Reference serialization does not retain external member means; do not silently export the old array");
			require(config.lengthPolicy==null,
				"Reference serialization does not retain a transient length-penalty deployment policy");
			require(config.inverseRarityRanking==null,
				"Reference serialization does not retain a transient inverse-rarity ranking binding");
			require(!config.alignmentMetric.same(FamilyAlignmentMetric.LEGACY) || config.heuristicRanking==null,
				"Reference serialization does not retain an external caller ranking binding");
			require(maxBytes>0, "Decoded preparation output exceeds maxoutputbytes before final references");
			require(out!=null && !Files.exists(out) && !Files.isSymbolicLink(out),
				"Calibration reference output exists or is a symbolic link: "+out);
			final ByteStreamWriter writer=new ByteStreamWriter(out.toString(),false,false,false);
			final ByteBuilder buffer=new ByteBuilder();
			long bytes=0;
			writer.start();
			try{
				FamilyReferenceCodec.appendHeader(header,buffer);
				bytes=buffer.length();
				require(bytes<=maxBytes, "Decoded preparation output exceeds maxoutputbytes in final reference header");
				writer.print(buffer);
				for(FamilyReference family:all){
					FamilyReferenceCodec.appendFiveLines(family,buffer.clear(),header);
					bytes=Math.addExact(bytes, buffer.length());
					require(bytes<=maxBytes, "Decoded preparation output exceeds maxoutputbytes in final references");
					writer.print(buffer);
				}
			}finally{require(!writer.poisonAndWait(),"Calibration reference writer failed; no complete references");}
			return bytes;
		}
		private void requireCalibrationUpdatesAllowed(){
			requireCalibrationUpdatesAllowed(false);
		}
		/** Only GC attachment may retain an existing evaluation margin. All calibration
		 * and export callers use the no-argument guard and still reject that view. */
		private void requireCalibrationUpdatesAllowed(boolean allowEvaluationMargins){
			require(config.alignmentMetric.matrix.same(header.alignmentMetric.matrix),
				"An uncalibrated evaluation-matrix swap cannot be refitted or exported as calibrated references");
			require(config.familyGcPenalty==null || config.familyGcRecipeBound,
				"Transient GC ranking cannot be exported or recalibrated until its family resources "+
					"and coefficient have a reproducible package binding");
			require(config.hitFloors==null,
				"The family-floor experiment freezes alignment gates; reference export or recalibration must "+
					"not discard its vector");
			require(config.lengthPolicy==null,
				"The length-penalty deployment policy is transient and must not be serialized or recalibrated");
			if(config.alignmentMetric.same(FamilyAlignmentMetric.LEGACY)){
				require(config.familyAdmission==null,
					"Family-mean admission must not reuse calibration measured with another candidate policy");
				require(!config.lengthScaledShortAdmission,
					"The short-query admission comparison freezes family gates; recalibration is not part of this policy");
				require(config.margins==null,"Evaluation-only margins cannot be used for calibration updates");
				require(config.memberScalarMeans==null,"Experimental member-mean arrays have no released recalibration recipe");
				require(config.heuristicRanking==null,"An explicit caller ranking policy has no released recalibration recipe");
				require(config.inverseRarityRanking==null,
					"An inverse-rarity caller ranking policy has no released recalibration recipe");
				return;
			}
			require(config.heuristicRanking!=null && config.heuristicRanking.armName().equals("hitfourth_h020_c150") &&
				(config.width==512 || config.width==1024) && config.seed==0 &&
				config.ddlReject>=0 && config.ddlReject<config.ddlAccept && config.ddlAccept==1025 &&
				config.releasedMetricProfileGeometry() &&
				!config.index.usesPruning() && !config.index.usesDirectCompletion() &&
				config.features==RankingFeatures.COMBINED25 && config.cosineMode==CosineMode.RANKING &&
				config.scalarBounds!=null && config.familyAdmission==null && !config.lengthScaledShortAdmission &&
				(allowEvaluationMargins || config.margins==null) && config.memberScalarMeans==null && config.lengthRanking==null &&
				config.lengthMode==FamilyLengthScore.Mode.OFF && config.inverseRarityRanking==null &&
				!config.union && config.truth==null && !config.candidateOnly,
				"Metric-bound calibration updates are released only for the fixed P4/O2 or plain common-pivot classifier");
		}
		SearchConfig config(){return config;}
		private static int[] liveToPhysical(FamilyReference[] live,FamilyReference[] physical){
			final int[] map=new int[live.length];
			int next=0;
			for(int slot=0;slot<physical.length;slot++){
				final FamilyReference family=physical[slot];
				require(family!=null,"Physical reference roster must be complete");
				if(family.members>0){
					require(next<live.length && live[next]==family,"Live reference roster must follow physical rank order");
					map[next++]=slot;
				}
			}
			require(next==live.length,"Live reference roster has entries missing from the physical roster");
			return map;
		}
		private final SearchConfig config;
		private final FamilyReferenceCodec.Header header;
		private final FamilyReference[] all;
		private final int[] liveToPhysical;
	}

	/** One index/model binding shared by all panel workers and all three ranking modes. */
	public static final class LengthPairedClassifier{
		private LengthPairedClassifier(LoadedClassifier base,FamilyLengthScore binding,double[] scales){
			require(base!=null && binding!=null && scales!=null && scales.length==3,
				"Paired classification requires OFF/multiply/add absolute margin scales");
			for(int arm=0;arm<3;arm++){
				configs[arm]=base.withLengthRanking(binding,LENGTH_MODES[arm]).withEvaluationMarginScale(scales[arm]).config;
			}
		}
		public LengthPairedWorker newWorker(){return new LengthPairedWorker(configs);}
		public int familyCount(){return configs[0].families.length;}
		public int familyRankAt(int ordinal){return configs[0].families[ordinal].rank;}
		public String familyIdAt(int ordinal){return configs[0].families[ordinal].id;}
		private final SearchConfig[] configs=new SearchConfig[3];
	}

	/** Receives borrowed arrays in OFF/multiply/add1 order, synchronously; never retain them or re-enter the worker.
	 * Oracle scores are the existing glocal raw-score/full-query-self metric. A mask bit
	 * denotes strict pre-scalar top50 membership, not acceptance. Empty union is ordinal-1/mask0.
	 */
	@FunctionalInterface
	public interface LengthPairedObserver{
		void accept(int queryLength,String[] terminals,int[] familyOrdinals,int[] ddlMatches,
			float[] identities,float[] blosumScores,
			int oracleOrdinal,int oracleMask,int unionSize);
	}

	/** Worker-local union alignment cache; no duplicate DP across ranking arms. */
	public static final class LengthPairedWorker{
		private LengthPairedWorker(SearchConfig[] configs_){
			configs=configs_;searcher=new Searcher(configs[0]);
			final int size=configs[0].families.length;
			seen=new long[size];identityCache=new double[size];blosumCache=new double[size];
			provider=(query,candidate,self,out)->{
				require(seen[candidate.ordinal]==generation,
					"Native classifier may only reuse metrics from this query's strict-list union");
				out.identityPct=identityCache[candidate.ordinal];out.blosum=blosumCache[candidate.ordinal];
			};
		}
		/** Reuses index scores and NN predictions; each union member is DP-aligned exactly once. */
		public void assignObserve(String id,byte[] raw,LengthPairedObserver observer){
			require(observer!=null && generation<Long.MAX_VALUE,
				"Paired observer and a fresh finite query generation are required");
			generation++;
			final ProteinSequence query=new ProteinSequence(id,raw);
			final byte[] clipped=query.length()==raw.length ? raw : Arrays.copyOf(raw,query.length());
			int unionSize=0;
			for(int arm=0;arm<3;arm++){
				searcher.prepare(clipped,query,result,counters,null,arm>0,LENGTH_MODES[arm]);
				topSizes[arm]=searcher.rankedTopSize;
				System.arraycopy(searcher.rankedTop,0,tops[arm],0,topSizes[arm]);
				for(int i=0;i<topSizes[arm];i++){
					final int ordinal=tops[arm][i];
					if(seen[ordinal]!=generation){
						require(unionSize<union.length,"Three strict top50 lists cannot exceed150 distinct families");
						seen[ordinal]=generation;union[unionSize++]=ordinal;
					}
				}
			}
			final int self=ProteinSearcher.selfScore(query);int oracle=-1;double best=0;
			for(int i=0;i<unionSize;i++){
				final int ordinal=union[i];final FamilyReference family=configs[0].families[ordinal];
				FamilyCalibration.align(query,family.consensus,self,metrics);
				identityCache[ordinal]=metrics.identityPct;blosumCache[ordinal]=metrics.blosum;
				if(oracle<0 || metrics.blosum>best ||
						metrics.blosum==best && family.rank<configs[0].families[oracle].rank){
					oracle=ordinal;best=metrics.blosum;
				}
			}
			int mask=0;
			for(int arm=0;arm<3;arm++){
				for(int i=0;i<topSizes[arm];i++){if(tops[arm][i]==oracle){mask|=1<<arm;}}
				if(searcher.prepare(clipped,query,result,counters,null,true,LENGTH_MODES[arm])){
					searcher.alignPrepared(query,result,counters,null,provider,configs[arm].margins);
				}
				terminals[arm]=result.terminal;winners[arm]=result.family;ddl[arm]=result.ddl;
				identities[arm]=result.identity;blosum[arm]=result.blosum;
				require(terminals[arm]!=null && winners[arm]>=-1 && winners[arm]<configs[0].families.length &&
					Float.isFinite(identities[arm]) && Float.isFinite(blosum[arm]),
					"Paired native outcomes must contain finite metrics and a valid terminal");
			}
			require((unionSize==0)==(oracle<0) && (oracle<0 ? mask==0 : mask>0),
				"The shared oracle must belong to at least one nonempty strict list");
			observer.accept(query.length(),terminals,winners,ddl,identities,blosum,oracle,mask,unionSize);
		}
		private final SearchConfig[] configs;
		private final Searcher searcher;
		private final Result result=new Result();private final Counters counters=new Counters();
		private final FamilyCalibration.Metrics metrics=new FamilyCalibration.Metrics();
		private final AlignmentProvider provider;
		private final int[][] tops=new int[3][50];private final int[] topSizes=new int[3],union=new int[150];
		private final long[] seen;private long generation;
		private final double[] identityCache,blosumCache;
		private final String[] terminals=new String[3];private final int[] winners=new int[3],ddl=new int[3];
		private final float[] identities=new float[3],blosum=new float[3];
	}
	private static final FamilyLengthScore.Mode[] LENGTH_MODES={FamilyLengthScore.Mode.OFF,
		FamilyLengthScore.Mode.MULTIPLY,FamilyLengthScore.Mode.ADD};
	@FunctionalInterface
	private interface AlignmentProvider{
		void align(ProteinSequence query,Comparison candidate,int self,FamilyCalibration.Metrics out);
	}

	public static final class Worker{
		Worker(SearchConfig config_){this(config_,false);}
		Worker(SearchConfig config_,boolean stageCpu_){this(config_,stageCpu_,null);}
		Worker(SearchConfig config_,boolean stageCpu_,FamilyMatchInformation.Binding information){
			require(!stageCpu_ || information==null,
				"Information timing needs an explicitly bound worker; do not silently discard its binding");
			config=config_;
			timing=stageCpu_ ? new StageCpu() : null;
			searcher=timing==null ? new Searcher(config_,false,information) : new TimedSearcher(config_,timing);
		}
		/** Single-threaded scratch owner; do not share one worker across caller threads. */
		public Assignment assign(String id,byte[] raw){
			score(id,raw);
			return new Assignment(result,config);
		}
		/** Explicit measured/neutral GC for this query only; unavailable GC uses zero storage. */
		public Assignment assign(String id,byte[] raw,boolean gcAvailable,double queryGc){
			score(id,raw,gcAvailable,queryGc);
			return new Assignment(result,config);
		}
		/** Development-only primary alignment counts without another alignment or path clone.
		 * The callback runs once per query, with measured counts only for an accepted winner.
		 * Counts are borrowed scratch: read during this call, never mutate or retain them.
		 * Callback exceptions propagate; scratch detaches before the next query.
		 */
		public void assignResidues(String id,byte[] raw,ResidueObserver observer){
			requireNoQueryGc();assignResidues(id,raw,false,0,observer);
		}
		/** Explicit query GC with the same borrowed residue-observer lifetime. */
		public void assignResidues(String id,byte[] raw,boolean gcAvailable,double queryGc,ResidueObserver observer){
			require(observer!=null && !active && searcher.residueCounts==null,
				"Residue observation cannot reenter or replace an active classifier query");
			require(config.ddlAccept>config.width && config.cosineMode!=CosineMode.FULL &&
				!config.alignmentMetric.same(FamilyAlignmentMetric.LEGACY) && !config.union &&
				!config.candidateOnly && config.truth==null,
				"Development counts require ordinary pinned alignment-only classification");
			if(residueScratch==null){residueScratch=new FamilyAlignmentCounts();}
			residueScratch.clear();searcher.residueCounts=residueScratch;
			try{
				score(id,raw,gcAvailable,queryGc);
				final Assignment assignment=new Assignment(result,config);
				if(assignment.assigned()){
					require(residueScratch.measured(),"Accepted winner lacks its exact reference/match traceback counts");
				}else{residueScratch.clear();}
				active=true;
				observer.accept(assignment,residueScratch);
			}finally{active=false;searcher.residueCounts=null;residueScratch.clear();}
		}
		/** Observe the existing raw, ranking, scalar and alignment decisions for one query.
		 * No candidate is added and no extra alignment is performed. Complete raw counts
		 * require exhaustive indexing; cosine-gated and reused diagnostic paths are separate.
		 * Callbacks run synchronously while this worker rejects reentry. A thrown callback
		 * propagates, detaches the observer and leaves the next query free to reset scratch.
		 * Returns the same immutable terminal assignment as assign().
		 */
		public Assignment assignTrace(String id,byte[] raw,CandidateObserver observer){
			requireNoQueryGc();return assignTrace(id,raw,false,0,observer);
		}
		/** Explicit query GC with ordinary candidate/terminal tracing and reentry guards. */
		public Assignment assignTrace(String id,byte[] raw,boolean gcAvailable,double queryGc,CandidateObserver observer){
			require(observer!=null && !active && searcher.candidateObserver==null,
				"Candidate tracing needs an observer and cannot reenter or replace an active query");
			require(!config.index.usesPruning() && !config.index.usesDirectCompletion() &&
				!config.union && config.truth==null && !config.candidateOnly &&
				config.cosineMode!=CosineMode.FULL && config.ddlAccept>config.width,
				"Candidate tracing requires exhaustive ordinary alignment-only search without cosine gates or diagnostic reuse");
			searcher.candidateObserver=observer;
			try{assignObserve(id,raw,gcAvailable,queryGc,observer,null);return new Assignment(result,config);}
			finally{searcher.candidateObserver=null;}
		}
		/** Observe the deployed shortlist and scalar gates without performing any alignments.
		 * Uses the same exhaustive preparation and strict top cap as assignTrace(), including
		 * terminal-stop clipping and no backfill after scalar rejection. All ranking/scalar
		 * callbacks complete before the terminal callback. No alignment callback is emitted.
		 * The terminal is candidates_selected, ddl_reject, scalar_reject or no_hit. Its family
		 * and DDL match fields are -1, scores are unmeasured zero placeholders, and state is
		 * NO_ALIGNMENT. Selection is not an assignment or a biological acceptance decision.
		 * Callback exceptions propagate; reentry is forbidden and scratch detaches in finally.
		 */
		public Selection selectTrace(String id,byte[] raw,CandidateObserver observer){
			requireNoQueryGc();return selectTrace(id,raw,false,0,observer);
		}
		/** Explicit query GC for selection-only tracing; clears state even when a callback fails. */
		public Selection selectTrace(String id,byte[] raw,boolean gcAvailable,double queryGc,CandidateObserver observer){
			require(observer!=null && !active && searcher.candidateObserver==null,
				"Candidate selection needs an observer and cannot reenter or replace an active query");
			require(!config.index.usesPruning() && !config.index.usesDirectCompletion() &&
				!config.union && config.truth==null && !config.candidateOnly &&
				config.cosineMode!=CosineMode.FULL && config.ddlAccept>config.width,
				"Candidate selection requires exhaustive ordinary alignment-only search without cosine gates or diagnostic reuse");
			final long start=timing==null ? 0 : timing.now();
			active=true;searcher.candidateObserver=observer;
			try{
				searcher.setQueryGc(gcAvailable,queryGc);
				final ProteinSequence query=new ProteinSequence(id,raw);
				final byte[] clipped=query.length()==raw.length ? raw : Arrays.copyOf(raw,query.length());
				final boolean selected=searcher.prepare(clipped,query,result,counters,null,false);
				require(result.family==-1 && result.dpAttempts==0 &&
					result.alignmentState==FamilyAlignmentMetric.State.NO_ALIGNMENT,
					"Selection-only preparation cannot assign a family or perform an alignment");
				if(selected){result.terminal="candidates_selected";}
				final Selection selection=new Selection(query.length(),result);
				observer.admission(query.length(),result.admissionFloor,result.admissionFloorMax,result.ddlActive,0);
				observer.accept(query.length(),selection.terminal,-1,-1,0,0,FamilyAlignmentMetric.State.NO_ALIGNMENT);
				return selection;
			}finally{
				searcher.clearQueryGc();searcher.candidateObserver=null;active=false;
				if(timing!=null){counters.stageWorkerNanos+=timing.now()-start;}
			}
		}
		/** Scores once and calls the observer synchronously without allocating an Assignment.
		 * The observer receives values, never mutable model or worker state. It must not
		 * re-enter this worker. A family ordinal of -1 means no accepted family.
		 */
		public void assignObserve(String id,byte[] raw,AssignmentObserver observer){
			assignObserve(id,raw,observer,null);
		}
		/** Explicit measured/missing query GC and one synchronous terminal observation. */
		public void assignObserve(String id,byte[] raw,boolean gcAvailable,double queryGc,AssignmentObserver observer){
			assignObserve(id,raw,gcAvailable,queryGc,observer,null);
		}
		/** Opt-in raw observation runs once after exhaustive index scoring, including
		 * queries rejected by every floor. Borrowed workspace arrays expire on return;
		 * the callback must neither mutate/retain them nor reenter this worker. */
		public void assignObserve(String id,byte[] raw,AssignmentObserver observer,RawHitObserver rawObserver){
			requireNoQueryGc();assignObserve(id,raw,false,0,observer,rawObserver);
		}
		/** Explicit query GC with optional uncensored raw-hit observation from the same index pass. */
		public void assignObserve(String id,byte[] raw,boolean gcAvailable,double queryGc,
				AssignmentObserver observer,RawHitObserver rawObserver){
			require(observer!=null && !active && searcher.rawHitObserver==null,
				"An assignment observer is required and callbacks cannot reenter their worker");
			if(rawObserver!=null){
				require(!config.index.usesPruning() && !config.index.usesDirectCompletion() &&
					!config.union && config.truth==null && !config.candidateOnly,
					"Uncensored raw-hit observation requires one ordinary exhaustive-index query, "+
						"not a diagnostic or reused-score grid");
			}
			final int length;
			searcher.rawHitObserver=rawObserver;
			try{length=score(id,raw,gcAvailable,queryGc);}finally{searcher.rawHitObserver=null;}
			active=true;
			try{
				if(!Float.isFinite(result.identity) || result.identity<0 || result.identity>100 ||
						!Float.isFinite(result.blosum) || result.ddl< -1 || result.ddl>config.width ||
						result.family< -1 || result.family>=config.families.length){
					throw new IllegalStateException("Observed assignment contains invalid family or numeric state");
				}
				observer.admission(length,result.admissionFloor,result.admissionFloorMax,result.ddlActive,result.dpAttempts);
				if(result.assigned()){
					require(result.rankedPosition>=1 && result.rankedPosition<=config.top ||
						result.rankedPosition==-1 && (result.terminal.equals("ddl_accept") || result.terminal.equals("cosine_accept")),
						"An accepting rank must come from the final capped shortlist, or be unavailable for a pre-ranking shortcut");
					observer.acceptingRank(result.family,result.rankedPosition);
				}
				final boolean measuredWinner=result.assigned() && result.alignmentState==FamilyAlignmentMetric.State.MEASURED;
				observer.accept(length,result.terminal,result.family,result.ddl,result.identity,result.blosum,result.alignmentState,
					measuredWinner ? result.querySpan : 0,measuredWinner ? result.referenceSpan : 0);
			}finally{active=false;}
		}
		private int score(String id,byte[] raw){
			requireNoQueryGc();return score(id,raw,false,0);
		}
		private void requireNoQueryGc(){
			require(config.familyGcPenalty==null,"A GC-bound worker requires explicit measured/missing query GC on every call");
		}
		private int score(String id,byte[] raw,boolean gcAvailable,double queryGc){
			require(!active,"Classifier callbacks cannot reenter their worker or replace its current query workspace");
			final long start=timing==null ? 0 : timing.now();
			active=true;
			try{
				searcher.setQueryGc(gcAvailable,queryGc);
				final ProteinSequence query=new ProteinSequence(id,raw);
				final byte[] clipped=query.length()==raw.length ? raw : Arrays.copyOf(raw,query.length());
				searcher.search(clipped,query,result,counters,null);
				if(result.assigned()){counters.assigned++;}else{counters.rejected++;}
				return query.length();
			}finally{
				searcher.clearQueryGc();active=false;
				if(timing!=null){counters.stageWorkerNanos+=timing.now()-start;}
			}
		}
		/** Cumulative ordinary-search candidate visits to the alignment stage. */
		public long alignmentAttempts(){return counters.candidateDp;}
		/** Visits skipped because the family lacks enabled alignment calibration. */
		public long alignmentUnavailable(){return counters.candidateAlignmentUnavailable;}
		/** Alignment calls which returned no usable measurement, including short spans. */
		public long alignmentFailures(){return counters.candidateAlignmentFailed;}
		/** Opt-in frozen-panel diagnostics add one baseline DP per actual overhang DP.
		 * Enable before scoring; never enable during fitting or for cached providers. */
		public void enableOverhangDiagnostics(){
			require(!active && counters.candidateDp==0 && searcher.overhangOnly && !config.union,
				"Overhang diagnostics require a fresh non-union overhang worker");
			searcher.overhangDiagnostics=true;
		}
		public long overhangCompared(){return searcher.overhangCompared;}
		public long overhangTrimPositive(){return searcher.overhangTrimPositive;}
		public long overhangTrimPositiveNonpositive(){return searcher.overhangTrimPositiveNonpositive;}
		public long overhangMeasured(){return searcher.overhangMeasured;}
		public long overhangSpanFailed(){return searcher.overhangSpanFailed;}
		public long overhangIdentityEligible(){return searcher.overhangIdentityEligible;}
		public long overhangIdentityRemoved(){return searcher.overhangIdentityRemoved;}
		public StageCpuCounts stageCpuCounts(){
			if(timing==null){return null;}
			counters.stageIndexNanos=timing.index;
			counters.stageShortlistNanos=timing.shortlist;
			counters.stageAlignmentNanos=timing.alignment;
			return new StageCpuCounts(counters);
		}
		private final SearchConfig config;
		private final Searcher searcher;
		private final StageCpu timing;
		private boolean active;
		private FamilyAlignmentCounts residueScratch;
		private final Result result=new Result();
		private final Counters counters=new Counters();
	}

	/** Immutable selection outcome. Selected families passed the strict cap and scalar gates;
	 * their ordinals arrive through CandidateObserver. No family has been aligned or assigned. */
	public static final class Selection{
		private Selection(int queryLength_,Result result){
			queryLength=queryLength_;terminal=result.terminal;
			rawHitFloor=result.admissionFloor;rawHitFloorMax=result.admissionFloorMax;
			admittedFamilies=result.ddlActive;selectedFamilies=result.dpShortlist;
		}
		public final int queryLength,rawHitFloor,rawHitFloorMax,admittedFamilies,selectedFamilies;
		public final String terminal;
	}

	/** Reduced timing counters for callers that use the packaged Worker API. */
	public static final class StageCpuCounts{
		public StageCpuCounts(){}
		private StageCpuCounts(Counters counters){
			workerNanos=counters.stageWorkerNanos;
			indexNanos=counters.stageIndexNanos;
			shortlistNanos=counters.stageShortlistNanos;
			alignmentNanos=counters.stageAlignmentNanos;
			validate();
		}
		public void add(StageCpuCounts other){
			if(other==null){return;}
			workerNanos+=other.workerNanos;
			indexNanos+=other.indexNanos;
			shortlistNanos+=other.shortlistNanos;
			alignmentNanos+=other.alignmentNanos;
			validate();
		}
		private void validate(){
			if(workerNanos<0 || indexNanos<0 || shortlistNanos<0 || alignmentNanos<0 ||
					workerNanos<indexNanos+shortlistNanos+alignmentNanos){
				throw new IllegalStateException("Stage CPU counters do not partition worker time");
			}
		}
		public long workerNanos,indexNanos,shortlistNanos,alignmentNanos;
	}

	/** One synchronous callback per completed full index score. Only touched[0..size)
	 * is populated; scores is indexed by live family ordinal. Untouched scores are zero.
	 * Arrays are borrowed and must not be retained, mutated, or used to reenter the worker. */
	@FunctionalInterface
	public interface RawHitObserver{void accept(int[] touched,int[] scores,int size);}
	@FunctionalInterface
	public interface ResidueObserver{void accept(Assignment assignment,FamilyAlignmentCounts counts);}

	/** Optional observations of actual candidate decisions, followed by the ordinary terminal.
	 * Ordinals identify the loaded live roster; clients can obtain IDs/ranks from LoadedClassifier.
	 * Only primitive values and immutable labels escape. Callbacks must not reenter the worker.
	 * Missing events mean the stage was not visited, never a measured zero or a passing gate.
	 */
	public interface CandidateObserver extends AssignmentObserver{
		/** One touched family after complete index scoring; all untouched families have zero hits. */
		default void rawCandidate(int familyOrdinal,int rawHits,int admissionFloor){}
		/** Complete raw-admitted inverse-rarity feature pool before strict top truncation.
		 * Arrays are borrowed worker scratch; read during this callback only. */
		default void fullAdmittedFeatures(String queryId,int queryLength,int size,int[] ordinals,int[] physicalRanks,
				int[] hits,double[] rarity,double[] lengthFactor,double[] combinedCosine,double[] hydrophobicitySimilarity,
				int[] scalarFailureMask){}
		/** Positions are one-based within the admitted raw order and the ranked window.
		 * The latter precedes scalar filtering; selected marks the strict top cap without backfill.
		 * Below-floor families have only raw events, not fabricated ranking positions. */
		default void rankedCandidate(int familyOrdinal,int rawHits,int rawPosition,int rankedPosition,boolean selected){}
		/** Actual scalar check for a selected candidate. Zero means no enabled/defined gate failed.
		 * No event means unchecked; bits use the same five-scalar ordering as scalarFailureMask. */
		default void scalarCandidate(int familyOrdinal,int rankedPosition,int failureMask){}
		/** Rounded production scores and the exact effective gates before acceptance/rejection.
		 * Identity acceptance precedes identity-lower rejection and strict BLOSUM acceptance.
		 * A score minus its upper gate may be negative; neither margin is clamped here. */
		default void alignedCandidate(int familyOrdinal,int attempt,float identity,float blosum,
				double identityUpper,double blosumUpper,float identityLower){}
		/** Measured attempt with the span operands used by optional length-penalty policy. */
		default void alignedCandidate(int familyOrdinal,int attempt,float identity,float blosum,
				double identityUpper,double blosumUpper,float identityLower,int querySpan,int referenceSpan){
			alignedCandidate(familyOrdinal,attempt,identity,blosum,identityUpper,blosumUpper,identityLower);
		}
		/** Attempt without measured scores: calibration unavailable or an actual failed alignment.
		 * The reason distinguishes unavailable calibration from a measured attempt's failure state. */
		default void unalignedCandidate(int familyOrdinal,int attempt,String reason,FamilyAlignmentMetric.State state){}
		/** Failed alignment with traceback coordinates. Meaningful only after pairwise DP ran. */
		default void unalignedCandidate(int familyOrdinal,int attempt,String reason,FamilyAlignmentMetric.State state,
				int queryStart,int queryEnd,int referenceStart,int referenceEnd){
			unalignedCandidate(familyOrdinal,attempt,reason,state);
		}
	}

	/** Receives one terminal result. Alignment scores are meaningful only on aligned paths. */
	@FunctionalInterface
	public interface AssignmentObserver{
		void accept(int queryLength,String terminal,int liveFamilyOrdinal,int ddlMatches,float identity,float blosum);
		/** Once per accepted query, before its terminal callback. Rank is one-based in
		 * the final sorted shortlist BEFORE scalar/GC filtering; -1 explicitly means
		 * a shortcut accepted before ranking. Rejections have no event. This is neither
		 * raw DDL position nor compressed alignment-attempt number. */
		default void acceptingRank(int liveFamilyOrdinal,int rankedPosition){}
		default void accept(int queryLength,String terminal,int liveFamilyOrdinal,int ddlMatches,float identity,float blosum,
				FamilyAlignmentMetric.State alignmentState){
			accept(queryLength,terminal,liveFamilyOrdinal,ddlMatches,identity,blosum);
		}
		/** Actual winning traceback spans, not query length. Both are zero for an
		 * unassigned query or a winner without a measured alignment. The default
		 * preserves the existing state-aware and six-argument terminal callbacks.
		 */
		default void accept(int queryLength,String terminal,int liveFamilyOrdinal,int ddlMatches,float identity,float blosum,
				FamilyAlignmentMetric.State alignmentState,int querySpan,int referenceSpan){
			accept(queryLength,terminal,liveFamilyOrdinal,ddlMatches,identity,blosum,alignmentState);
		}
		/** Untrimmed raw-hit-admitted pool and actual DP attempts, including rejected queries.
		 * Called once before accept; the default preserves existing observers and lambdas. */
		default void admission(int queryLength,int rawHitFloor,int admittedFamilies,int dpAttempts){}
		/** Configured floor range; both endpoints are equal for fixed/query policies.
		 * Legacy observers receive the minimum. Range-aware reporters override this overload. */
		default void admission(int queryLength,int minimumFloor,int maximumFloor,int admittedFamilies,int dpAttempts){
			admission(queryLength,minimumFloor,admittedFamilies,dpAttempts);
		}
	}

	public static final class Assignment{
		Assignment(Result result,SearchConfig config){
			rawHitFloor=result.admissionFloor; rawHitFloorMax=result.admissionFloorMax;
			familyMeanAdmission=config.familyAdmission!=null; admittedFamilies=result.ddlActive;
			familyHitFloorAdmission=config.hitFloors!=null;
			terminal=result.terminal; ddlMatches=result.ddl; ddlPosition=result.ddlPosition;
			dpAttempt=result.dpAttempt; dpAttempts=result.dpAttempts;
			shortlistCandidates=result.shortlistCandidates;alignmentCalls=result.alignmentCalls;
			alignmentAccepts=result.alignmentAccepts;alignmentInputSkips=result.alignmentInputSkips;
			alignmentUnavailable=result.alignmentUnavailable;alignmentProviderCalls=result.alignmentProviderCalls;
			identity=result.identity; blosum=result.blosum; alignmentState=result.alignmentState;
			scoreMeasured=result.combinedMeasured || result.cosineMeasured;
			identityMeasured=blosumMeasured=alignmentState==FamilyAlignmentMetric.State.MEASURED;
			score=result.combinedMeasured ? result.combinedCosine : result.cosineMeasured ? result.cosine : 0;
			if(result.assigned()){
				final FamilyReference family=config.families[result.family];
				familyId=family.id; rosterRank=family.rank;
			}else{familyId=null; rosterRank=-1;}
			alignmentMarginMeasured=result.assigned() && identityMeasured && blosumMeasured;
			informationMeasured=result.assigned() && result.informationMeasured;
			informationRatio=informationMeasured ? result.informationRatio : 0;
			rawScore=alignmentMarginMeasured ? result.rawScore : 0;
			denominatorSelf=alignmentMarginMeasured ? result.denominatorSelf : 0;
			querySpan=alignmentMarginMeasured ? result.querySpan : 0;
			referenceSpan=alignmentMarginMeasured ? result.referenceSpan : 0;
			if(alignmentMarginMeasured){
				final FamilyReference family=config.families[result.family];
				final double identityGate=config.margins==null ?
					Math.max(family.identityUpper,config.identityFloor) : config.margins.identity[result.family];
				final double blosumGate=config.margins==null ?
					Math.max(family.blosumFinal,config.blosumFloor) : config.margins.blosum[result.family];
				final double penalty=config.lengthPolicy==null ? 0 :
					FamilyLengthAcceptance.penalty(querySpan,referenceSpan,config.lengthAmplitude);
				alignmentMargin=normalizedAlignmentMargin(identity,blosum,identityGate+penalty,
					blosumGate+FamilyLengthAcceptance.BLOSUM_UNIT*penalty);
				require(alignmentMargin>0,"An alignment-accepted winner must strictly exceed an applied acceptance gate");
			}else{alignmentMargin=0;}
			validate();
		}
		/** Dimensionless winning-gate excess; uses the same rounded metrics and double gates as acceptance. */
		public static double normalizedAlignmentMargin(double identity,double blosum,double identityGate,double blosumGate){
			require(Double.isFinite(identity) && identity>=0 && identity<=100 && Double.isFinite(blosum) &&
				Double.isFinite(identityGate) && Double.isFinite(blosumGate),
				"Alignment margins need finite measured scores and applied gates");
			final double margin=Math.max(0,Math.max((identity-identityGate)/100,blosum-blosumGate));
			require(Double.isFinite(margin),"Alignment gate subtraction overflowed the finite margin");
			return margin;
		}
		public boolean assigned(){return familyId!=null;}
		/** Shared production/replay score contribution; a rejected query contributes zero. */
		public float stopBonus(float scale,boolean marginMode){
			require(Float.isFinite(scale) && scale>=0,"Stop bonus scale must be finite and nonnegative");
			if(!assigned()){return 0;}
			if(!marginMode){return scale;}
			require(alignmentMarginMeasured,
				"Accepted stop lacks measured alignment margin; shortcuts cannot supply graded evidence");
			final double value=scale*alignmentMargin;
			require(Double.isFinite(value) && value>=0 && value<=Float.MAX_VALUE,
				"Scaled alignment margin exceeds the finite caller score range");
			return (float)value;
		}
		/** Independent information rule shared by ordinary calling and replay. */
		public float informationBonus(float scale){
			require(Float.isFinite(scale) && scale>=0,"Information bonus scale must be finite and nonnegative");
			if(!assigned() || scale==0){return 0;}
			require(informationMeasured,"Accepted stop lacks measured information; no flat-bonus or zero fallback is allowed");
			final double value=scale*Math.min(1,informationRatio);
			require(Double.isFinite(value) && value>=0 && value<=Float.MAX_VALUE,
				"Information bonus exceeded the finite caller score range");
			return (float)value;
		}
		private void validate(){
			if(rawHitFloor<0 || rawHitFloorMax<rawHitFloor || admittedFamilies<0 || terminal==null ||
					ddlMatches<-1 || ddlMatches>1024 || ddlPosition<0 || dpAttempt<0 || dpAttempts<0 ||
					!Double.isFinite(score) || score<0 || score>1025 || !Float.isFinite(identity) ||
					identity<0 || identity>100 || !Float.isFinite(blosum) || (familyId==null)!=(rosterRank<0) ||
					denominatorSelf<0 || (!alignmentMarginMeasured && (rawScore!=0 || denominatorSelf!=0)) ||
					(alignmentMarginMeasured && denominatorSelf<=0) ||
					querySpan<0 || referenceSpan<0 || (querySpan==0)!=(referenceSpan==0) ||
					!Double.isFinite(alignmentMargin) || alignmentMargin<0 || (!alignmentMarginMeasured && alignmentMargin!=0) ||
					!Double.isFinite(informationRatio) || informationRatio<0 || (!informationMeasured && informationRatio!=0)){
				throw new IllegalStateException("Callable assignment exposed invalid numeric state");
			}
		}
		public final String terminal, familyId;
		public final int rosterRank, ddlMatches, ddlPosition, dpAttempt, dpAttempts;
		/** Ordinary-search work only; copied once before worker scratch is reused.
		 * shortlistCandidates is the ranked top cap before scalar filtering.
		 * Calls exclude unavailable calibration, short-input skips and cached providers.
		 * Search returns at the first accepted alignment, hence accepts is zero or one.
		 */
		public final int shortlistCandidates,alignmentCalls,alignmentAccepts,alignmentInputSkips,
			alignmentUnavailable,alignmentProviderCalls;
		/** Configured minimum/maximum floors and exact admitted pool before ranking or scalar gates.
		 * A family-dependent policy has no single query floor; each candidate uses its own bound. */
		public final int rawHitFloor,rawHitFloorMax,admittedFamilies;
		public final boolean familyMeanAdmission;
		public final boolean familyHitFloorAdmission;
		public final double score;
		public final float identity, blosum;
		public final FamilyAlignmentMetric.State alignmentState;
		public final boolean scoreMeasured, identityMeasured, blosumMeasured;
		/** Zero is inert when unmeasured; shortcuts and rejected queries never fabricate alignment evidence. */
		public final double alignmentMargin;
		public final boolean alignmentMarginMeasured;
		/** Uncapped I/I0 of the primary winner; measured zero is usable evidence. */
		public final double informationRatio;
		public final boolean informationMeasured;
		public final int rawScore, denominatorSelf;
		/** Actual measured traceback span; zero only when no real DP alignment was assigned. */
		public final int querySpan, referenceSpan;
	}

	private static final class CascadeAccumulator implements Accumulator<ProcessThread>{
		CascadeAccumulator(int unionFamilies){if(unionFamilies>0){counters.unionStats=new UnionStats(unionFamilies);}}
		@Override
		public void accumulate(ProcessThread pt){
			counters.add(pt.counters);
			success&=pt.success;
		}
		@Override
		public boolean success(){return success && failure.get()==null;}
		@Override
		public ReadWriteLock rwlock(){return rwlock;}
		final Counters counters=new Counters();
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final ReadWriteLock rwlock=new ReentrantReadWriteLock();
		boolean success=true;
	}

	private static final class ProcessThread extends Thread{
		ProcessThread(final int tid_,final Streamer input_,final ByteStreamWriter writer_,final ByteStreamWriter rankWriter_,
				final ByteStreamWriter auditWriter_,final SearchConfig config_,
				final AtomicReference<Throwable> sharedFailure_,final boolean stageCpu_){
			super("FamilyTailCascade-"+tid_);
			input=input_; writer=writer_; rankWriter=rankWriter_; auditWriter=auditWriter_; config=config_;
			sharedFailure=sharedFailure_; stageCpu=stageCpu_;
			if(config.union){counters.unionStats=new UnionStats(config.families.length);}
		}
		@Override
		public void run(){
			try{
				final StageCpu timing=stageCpu ? new StageCpu() : null;
				final Searcher searcher=timing==null ? new Searcher(config) : new TimedSearcher(config,timing);
				final Result result=new Result();
				final long workerStart=timing==null ? 0 : timing.now();
				for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
					final ByteBuilder out=new ByteBuilder();
					final ByteBuilder ranks=rankWriter==null ? null : new ByteBuilder();
					final ByteBuilder audit=auditWriter==null ? null : new ByteBuilder();
					try{
						if(sharedFailure.get()==null){
							final long start=System.nanoTime();
							try{process(batch,out,ranks,audit,searcher,result);}
							finally{counters.processNanos+=System.nanoTime()-start;}
						}
					}catch(Throwable t){
						out.clear();
						if(ranks!=null){ranks.clear();}
						if(audit!=null){audit.clear();}
						sharedFailure.compareAndSet(null,t);
					}finally{
						writer.add(out,batch.id);
						if(rankWriter!=null){rankWriter.add(ranks,batch.id);}
						if(auditWriter!=null){auditWriter.add(audit,batch.id);}
					}
				}
				if(timing!=null){
					counters.stageWorkerNanos=timing.now()-workerStart;
					counters.stageIndexNanos=timing.index;
					counters.stageShortlistNanos=timing.shortlist;
					counters.stageAlignmentNanos=timing.alignment;
				}
				success=(sharedFailure.get()==null);
			}catch(Throwable t){
				t.printStackTrace();
				KillSwitch.kill("FamilyTailCascade unrecoverable worker dispatch failure");
			}
		}
		private void process(final ListNum<Read> batch,final ByteBuilder out,final ByteBuilder ranks,final ByteBuilder audit,
				final Searcher searcher,final Result result){
			final ByteBuilder line=new ByteBuilder();
			counters.batches++;
			for(Read read : batch.list){
				if(read==null || read.bases==null || read.mate!=null){
					throw new IllegalArgumentException("FamilyTailCascade requires unpaired protein FASTA input");
				}
				final ProteinSequence query=new ProteinSequence(firstToken(read.id),read.bases);
				final byte[] raw=Arrays.copyOf(read.bases, query.length());
				counters.reads++; counters.residues+=query.length();
				if(config.union && config.unionQueries.putIfAbsent(query.id,Boolean.TRUE)!=null){
					throw new IllegalArgumentException("Duplicate union query ID: "+query.id);
				}
				searcher.search(raw,query,result,counters,config.union ? out : null);
				if(result.assigned()){counters.assigned++;}else{counters.rejected++;}
				if(!config.union){writeResult(out,line,result,query,config.families);}
				if(ranks!=null){writeRanks(ranks,result,query.id,config.top,config.families,config.scalarBounds!=null);}
				if(audit!=null){
					audit.append(query.id).tab().append(result.terminal).tab().append(result.ddlActive).tab()
						.append(result.cosineWindow).tab().append(result.dpShortlist).tab().append(result.dpAttempts).tab()
						.append(result.identityRejects).tab().append(result.blosumRejects).tab()
						.append(result.cosineWindow>0 && result.dpAttempts>0 ? "true" : "false");
					if(config.floorAudit){
						if(result.deferredFamily<0){audit.append("\tNA\tNA\tNA\tNA");}
						else{audit.tab().append(config.families[result.deferredFamily].id).tab()
							.append(Float.toString(result.deferredIdentity)).tab().append(Float.toString(result.deferredBlosum)).tab()
							.append(result.deferredBlosumPass ? "true" : "false");}
					}
					audit.nl();
				}
			}
		}
		final Streamer input;
		final ByteStreamWriter writer,rankWriter,auditWriter;
		final SearchConfig config;
		final boolean stageCpu;
		final AtomicReference<Throwable> sharedFailure;
		final Counters counters=new Counters();
		boolean success=false;
	}

	private static void writeResult(final ByteBuilder out,final ByteBuilder line,final Result result,
			final ProteinSequence query,final FamilyReference[] families){
		line.clear().append(query.id).tab().append(result.assigned() ? "assigned" : "no_hit").tab();
		if(result.assigned()){
			final FamilyReference family=families[result.family];
			line.append(family.id).tab().append(family.rank).tab().append(result.ddl).tab();
			if(result.combinedMeasured){line.append(Double.toString(result.combinedCosine));}
			else if(result.cosineMeasured){line.append(Float.toString(result.cosine));}
			else{line.append("NA");}
			line.tab();
			if(result.dpAttempt>0){line.append(Float.toString(result.identity)).tab().append(Float.toString(result.blosum));}
			else{line.append("NA\tNA");}
		}else{line.append("NA\tNA\tNA\tNA\tNA\tNA");}
		out.append(line).nl();
	}

	private static void writeRanks(ByteBuilder out,Result result,String query,int top,FamilyReference[] families,
			boolean combined){
		out.append(query).tab().append(result.truth==null ? "NA" : result.truth.family).tab()
			.append(TRUTH_STATUS[result.boundary]).tab().append(result.windowSize).tab().append(top);
		// Keep the original15 fields in place; the fourth mode is appended after DP attribution.
		for(int mode=0; mode<3; mode++){
			out.tab();
			if(result.truthRanks[mode]==0){out.append("NA");}else{out.append(result.truthRanks[mode]);}
		}
		for(int mode=0; mode<3; mode++){out.tab().append(TRUTH_STATUS[truthStatus(result,mode,top)]);}
		out.tab().append(result.assigned() ? families[result.family].id : "NA").tab();
		if(result.assigned()){
			out.append(result.ddlPosition).tab();
			if(result.dpAttempt==0){out.append("NA");}else{out.append(result.dpAttempt);}
		}else{out.append("NA\tNA");}
		out.tab().append(result.dpAttempts).tab();
		if(result.truthRanks[3]==0){out.append("NA");}else{out.append(result.truthRanks[3]);}
		out.tab().append(TRUTH_STATUS[truthStatus(result,3,top)]);
		if(combined){
			out.tab(); if(result.truthRanks[4]==0){out.append("NA");}else{out.append(result.truthRanks[4]);}
			out.tab().append(TRUTH_STATUS[truthStatus(result,4,top)]);
		}
		out.nl();
	}

	private static int truthStatus(Result result,int mode,int top){
		return result.boundary<4 ? result.boundary : result.truthRanks[mode]<=top ? 5 : 4;
	}

	private static boolean precedes(float score,int rank,float targetScore,int targetRank){
		final int order=Float.compare(score,targetScore);
		return order>0 || (order==0 && rank<targetRank);
	}

	private static boolean precedes(double score,int rank,double targetScore,int targetRank){
		final int order=Double.compare(score,targetScore);
		return order>0 || (order==0 && rank<targetRank);
	}

	static final class Truth{
		Truth(String family_,int ordinal_){family=family_; ordinal=ordinal_;}
		final String family;
		final int ordinal;
	}

	/** Retain exact tokens once; the map is immutable while worker threads read it. */
	static HashMap<String,Truth> loadTruth(Path file,FamilyReference[] families){
		final HashMap<String,Truth> result=new HashMap<String,Truth>();
		if(file==null){return result;}
		final HashMap<String,Integer> ordinals=new HashMap<String,Integer>();
		for(int i=0; i<families.length; i++){ordinals.put(families[i].id,i);}
		final ByteFile reader=ByteFile.makeByteFile(file.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			header(reader,"query_id\texpected_family_id");
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				lp.set(line);
				if(lp.terms()!=2){throw new IllegalArgumentException("Malformed truth row: require two tokens");}
				final String query=lp.parseString(0),family=lp.parseString(1);
				truthToken(query); truthToken(family);
				final Integer ordinal=ordinals.get(family);
				if(result.put(query,new Truth(family,ordinal==null ? -1 : ordinal.intValue()))!=null){
					throw new IllegalArgumentException("Duplicate truth query ID: "+query);
				}
			}
		}finally{if(reader.close()){throw new IllegalStateException("Truth map read failed");}}
		return result;
	}

	private static void truthToken(String token){
		if(token==null || token.isEmpty()){throw new IllegalArgumentException("Empty truth ID");}
		for(int i=0; i<token.length(); i++){
			if(token.charAt(i)<=32 || token.charAt(i)>=127){throw new IllegalArgumentException("Malformed truth ID: "+token);}
		}
	}

	static final class SearchConfig{
		SearchConfig(DDLFamilySearch.Index index_,FamilyReference[] families_,int top_,
				int ddlReject_,int ddlAccept_,ProteinKmerGeometry geometry_,long seed_,int width_,
				CosineMode cosineMode_,RankingFeatures features_,HashMap<String,Truth> truth_,
				ScalarBounds scalarBounds_,boolean union_,boolean scalarGates_,float identityFloor_,
				float blosumFloor_,boolean floorAudit_,boolean candidateOnly_,boolean positional_,
				int partitions_,int overlap_,DDLPentamerBlacklist blacklist_,boolean compactFamilyReferences_){
			this(index_,families_,top_,ddlReject_,ddlAccept_,geometry_,seed_,width_,cosineMode_,features_,truth_,scalarBounds_,
				union_,scalarGates_,identityFloor_,blosumFloor_,floorAudit_,candidateOnly_,
				positional_ ? ProteinPositionMode.FIXED : ProteinPositionMode.PLAIN,partitions_,overlap_,
				blacklist_,compactFamilyReferences_);
		}
		SearchConfig(DDLFamilySearch.Index index_,FamilyReference[] families_,int top_,
				int ddlReject_,int ddlAccept_,ProteinKmerGeometry geometry_,long seed_,int width_,
				CosineMode cosineMode_,RankingFeatures features_,HashMap<String,Truth> truth_,
				ScalarBounds scalarBounds_,boolean union_,boolean scalarGates_,float identityFloor_,
				float blosumFloor_,boolean floorAudit_,boolean candidateOnly_,
				ProteinPositionMode positionMode_,int partitions_,int overlap_,DDLPentamerBlacklist blacklist_,
				boolean compactFamilyReferences_){
			this(index_,families_,top_,ddlReject_,ddlAccept_,geometry_,seed_,width_,cosineMode_,features_,truth_,scalarBounds_,
				union_,scalarGates_,identityFloor_,blosumFloor_,floorAudit_,candidateOnly_,positionMode_,partitions_,overlap_,
				blacklist_,compactFamilyReferences_,null);
		}
		/** All-hit diagnostic views retain the same optional scalar binding as their source classifier. */
		SearchConfig(DDLFamilySearch.Index index_,FamilyReference[] families_,int top_,
				int ddlReject_,int ddlAccept_,ProteinKmerGeometry geometry_,long seed_,int width_,
				CosineMode cosineMode_,RankingFeatures features_,HashMap<String,Truth> truth_,
				ScalarBounds scalarBounds_,boolean union_,boolean scalarGates_,float identityFloor_,
				float blosumFloor_,boolean floorAudit_,boolean candidateOnly_,
				ProteinPositionMode positionMode_,int partitions_,int overlap_,DDLPentamerBlacklist blacklist_,
				boolean compactFamilyReferences_,FamilyMemberScalarMeans memberScalarMeans_){
			this(index_,families_,top_,ddlReject_,ddlAccept_,geometry_,seed_,width_,cosineMode_,features_,truth_,scalarBounds_,
				union_,scalarGates_,identityFloor_,blosumFloor_,floorAudit_,candidateOnly_,positionMode_,partitions_,overlap_,
				blacklist_,compactFamilyReferences_,memberScalarMeans_,FamilyAlignmentMetric.LEGACY);
		}
		SearchConfig(DDLFamilySearch.Index index_,FamilyReference[] families_,int top_,
				int ddlReject_,int ddlAccept_,ProteinKmerGeometry geometry_,long seed_,int width_,
				CosineMode cosineMode_,RankingFeatures features_,HashMap<String,Truth> truth_,
				ScalarBounds scalarBounds_,boolean union_,boolean scalarGates_,float identityFloor_,
				float blosumFloor_,boolean floorAudit_,boolean candidateOnly_,
				ProteinPositionMode positionMode_,int partitions_,int overlap_,DDLPentamerBlacklist blacklist_,
				boolean compactFamilyReferences_,FamilyMemberScalarMeans memberScalarMeans_,
				FamilyAlignmentMetric alignmentMetric_){
			require(Float.isFinite(identityFloor_) && Float.isFinite(blosumFloor_),
				"Accept floors must be finite");
			require(positionMode_!=null,"Classifier mode must match the validated replacement profiles and query sketches");
			require(alignmentMetric_!=null,"Classifier must retain its calibrated alignment metric");
			index=index_; families=families_; top=top_; ddlReject=ddlReject_;
			ddlAccept=ddlAccept_; geometry=geometry_; k=geometry.k(); seed=seed_; width=width_; cosineMode=cosineMode_;
			features=features_; truth=truth_;
			scalarBounds=scalarBounds_;
			memberScalarMeans=memberScalarMeans_;
			union=union_; unionQueries=union ? new java.util.concurrent.ConcurrentHashMap<String,Boolean>() : null;
			scalarGates=scalarGates_;
			candidateOnly=candidateOnly_;
			positionMode=positionMode_; positional=positionMode.positional();
			partitions=partitions_; overlap=overlap_;
			ddl.ProteinBlacklistPolicy.requireStored(geometry,blacklist_==null ? null : blacklist_.sha80());
			blacklist=ddl.ProteinBlacklistPolicy.resolve(geometry,blacklist_);
			alignmentMetric=alignmentMetric_;
			identityFloor=identityFloor_; blosumFloor=blosumFloor_; floorAudit=floorAudit_;
			identityLowerReduction=0;
			margins=null;
			lengthRanking=null;lengthMode=FamilyLengthScore.Mode.OFF;heuristicRanking=null;inverseRarityRanking=null;
			familyGcPenalty=null;
			familyGcRecipeBound=false;
			lengthScaledShortAdmission=false;
			familyAdmission=null;
			hitFloors=null;
			lengthPolicy=null;lengthAmplitude=0;
			monomerSpectra=features!=RankingFeatures.DIMER || truth!=null ?
				(compactFamilyReferences_ ? new float[families.length][] : new float[families.length][MONOMER_WIDTH]) : null;
			monomerNorms=monomerSpectra==null ? null : new float[families.length];
			if(monomerSpectra!=null){
				// Shared immutable cache: do not repeat family row sums in each worker or per query.
				for(int i=0; i<families.length; i++){
					if(compactFamilyReferences_){
						if(families[i].spectrum.length!=MONOMER_WIDTH){
							throw new IllegalArgumentException("Compact family reference has non-projected spectrum: "+families[i].id);
						}
						monomerSpectra[i]=families[i].spectrum;
						monomerNorms[i]=families[i].spectrumNorm;
					}else{
						dimerRowSums(families[i].spectrum,monomerSpectra[i]);
						monomerNorms[i]=norm(monomerSpectra[i]);
					}
				}
			}
			centeredSpectra=features==RankingFeatures.CENTERED_MONOMER || truth!=null ?
				new float[families.length][MONOMER_WIDTH] : null;
			background=centeredSpectra==null ? null : monomerBackground(families,monomerSpectra);
			centeredNorms=centeredSpectra==null ? null : new float[families.length];
			if(centeredSpectra!=null){
				for(int i=0; i<families.length; i++){
					centerMonomers(monomerSpectra[i],background,centeredSpectra[i]);
					centeredNorms[i]=norm(centeredSpectra[i]);
				}
			}
			combinedSpectra=scalarBounds==null ? null : new double[families.length][COMBINED_WIDTH];
			combinedNorms=scalarBounds==null ? null : new double[families.length];
			fillCombined();
		}
		private SearchConfig(SearchConfig source,FamilyReference[] families_,EvaluationMargins margins_){
			this(source,families_,margins_,source.lengthRanking,source.lengthMode);
		}
		private SearchConfig(SearchConfig source,FamilyGcPenalty penalty,boolean recipeBound){
			this(source,source.families,source.margins,source.lengthRanking,source.lengthMode,source.memberScalarMeans,
				source.heuristicRanking,source.inverseRarityRanking,source.lengthScaledShortAdmission,source.familyAdmission,
				source.hitFloors,source.index,source.blacklist,source.identityLowerReduction,source.lengthPolicy,
				source.alignmentMetric,source.lengthAmplitude,penalty,recipeBound);
		}
		private SearchConfig(SearchConfig source,FamilyLengthScore binding,FamilyLengthScore.Mode mode){
			this(source,source.families,source.margins,binding,mode);
		}
		private SearchConfig(SearchConfig source,FamilyReference[] families_,EvaluationMargins margins_,
				FamilyLengthScore binding,FamilyLengthScore.Mode mode){
			this(source,families_,margins_,binding,mode,source.memberScalarMeans);
		}
		private SearchConfig(SearchConfig source,FamilyReference[] families_,EvaluationMargins margins_,
				FamilyLengthScore binding,FamilyLengthScore.Mode mode,FamilyMemberScalarMeans memberScalarMeans_){
			this(source,families_,margins_,binding,mode,memberScalarMeans_,source.heuristicRanking);
		}
		private SearchConfig(SearchConfig source,FamilyReference[] families_,EvaluationMargins margins_,
				FamilyLengthScore binding,FamilyLengthScore.Mode mode,FamilyMemberScalarMeans memberScalarMeans_,
				FamilyHeuristicScore heuristicRanking_){
			this(source,families_,margins_,binding,mode,memberScalarMeans_,heuristicRanking_,source.inverseRarityRanking);
		}
		private SearchConfig(SearchConfig source,FamilyReference[] families_,EvaluationMargins margins_,
				FamilyLengthScore binding,FamilyLengthScore.Mode mode,FamilyMemberScalarMeans memberScalarMeans_,
				FamilyHeuristicScore heuristicRanking_,
				FamilyInverseRarityScore inverseRarityRanking_){
			this(source,families_,margins_,binding,mode,memberScalarMeans_,heuristicRanking_,inverseRarityRanking_,
				source.lengthScaledShortAdmission,source.familyAdmission,source.hitFloors,source.index,source.blacklist);
		}
		private SearchConfig(SearchConfig source,FamilyReference[] families_,EvaluationMargins margins_,
				FamilyLengthScore binding,FamilyLengthScore.Mode mode,FamilyMemberScalarMeans memberScalarMeans_,
				FamilyHeuristicScore heuristicRanking_,
				boolean lengthScaledShortAdmission_){
			this(source,families_,margins_,binding,mode,memberScalarMeans_,heuristicRanking_,
				lengthScaledShortAdmission_,source.familyAdmission);
		}
		private SearchConfig(SearchConfig source,FamilyReference[] families_,EvaluationMargins margins_,
				FamilyLengthScore binding,FamilyLengthScore.Mode mode,FamilyMemberScalarMeans memberScalarMeans_,
				FamilyHeuristicScore heuristicRanking_,
				boolean lengthScaledShortAdmission_,FamilyLengthAdmission familyAdmission_){
			this(source,families_,margins_,binding,mode,memberScalarMeans_,heuristicRanking_,
				lengthScaledShortAdmission_,familyAdmission_,source.hitFloors);
		}
		private SearchConfig(SearchConfig source,FamilyReference[] families_,EvaluationMargins margins_,
				FamilyLengthScore binding,FamilyLengthScore.Mode mode,FamilyMemberScalarMeans memberScalarMeans_,
				FamilyHeuristicScore heuristicRanking_,
				boolean lengthScaledShortAdmission_,FamilyLengthAdmission familyAdmission_,FamilyHitFloorVector hitFloors_){
			this(source,families_,margins_,binding,mode,memberScalarMeans_,heuristicRanking_,
				lengthScaledShortAdmission_,familyAdmission_,hitFloors_,source.index,source.blacklist);
		}
		private SearchConfig(SearchConfig source,FamilyReference[] families_,EvaluationMargins margins_,
				FamilyLengthScore binding,FamilyLengthScore.Mode mode,FamilyMemberScalarMeans memberScalarMeans_,
				FamilyHeuristicScore heuristicRanking_,
				boolean lengthScaledShortAdmission_,FamilyLengthAdmission familyAdmission_,FamilyHitFloorVector hitFloors_,
				DDLFamilySearch.Index index_,DDLPentamerBlacklist blacklist_){
			this(source,families_,margins_,binding,mode,memberScalarMeans_,heuristicRanking_,source.inverseRarityRanking,
				lengthScaledShortAdmission_,familyAdmission_,hitFloors_,index_,blacklist_,source.identityLowerReduction);
		}
		private SearchConfig(SearchConfig source,FamilyReference[] families_,EvaluationMargins margins_,
				FamilyLengthScore binding,FamilyLengthScore.Mode mode,FamilyMemberScalarMeans memberScalarMeans_,
				FamilyHeuristicScore heuristicRanking_,
				boolean lengthScaledShortAdmission_,FamilyLengthAdmission familyAdmission_,FamilyHitFloorVector hitFloors_,
				DDLFamilySearch.Index index_,DDLPentamerBlacklist blacklist_,float identityLowerReduction_){
			this(source,families_,margins_,binding,mode,memberScalarMeans_,heuristicRanking_,source.inverseRarityRanking,
				lengthScaledShortAdmission_,familyAdmission_,hitFloors_,index_,blacklist_,identityLowerReduction_,
				source.lengthPolicy,source.alignmentMetric,source.lengthAmplitude);
		}
		/** Reuse matrix-independent index/ranking arrays, retaining every original gate and policy. */
		private SearchConfig(SearchConfig source,FamilyAlignmentMetric metric){
			this(source,source.families,source.margins,source.lengthRanking,source.lengthMode,source.memberScalarMeans,
				source.heuristicRanking,source.inverseRarityRanking,source.lengthScaledShortAdmission,
				source.familyAdmission,source.hitFloors,
				source.index,source.blacklist,source.identityLowerReduction,source.lengthPolicy,metric,source.lengthAmplitude);
		}
		private SearchConfig(SearchConfig source,FamilyLengthPolicy lengthPolicy_){
			this(source,lengthPolicy_,lengthPolicy_.amplitude());
		}
		private SearchConfig(SearchConfig source,FamilyLengthPolicy lengthPolicy_,float lengthAmplitude_){
			this(source,source.families,source.margins,source.lengthRanking,source.lengthMode,source.memberScalarMeans,
				source.heuristicRanking,source.inverseRarityRanking,source.lengthScaledShortAdmission,
				source.familyAdmission,source.hitFloors,source.index,
				source.blacklist,source.identityLowerReduction,lengthPolicy_,
				new FamilyAlignmentMetric(source.alignmentMetric.mode,FamilyLengthAcceptance.MIN_OVERLAP,
					source.alignmentMetric.matrix,source.alignmentMetric.minimumPaired),lengthAmplitude_);
		}
		private SearchConfig(SearchConfig source,FamilyReference[] families_,EvaluationMargins margins_,
				FamilyLengthScore binding,FamilyLengthScore.Mode mode,FamilyMemberScalarMeans memberScalarMeans_,
				FamilyHeuristicScore heuristicRanking_,
				FamilyInverseRarityScore inverseRarityRanking_,
				boolean lengthScaledShortAdmission_,FamilyLengthAdmission familyAdmission_,FamilyHitFloorVector hitFloors_,
				DDLFamilySearch.Index index_,DDLPentamerBlacklist blacklist_){
			this(source,families_,margins_,binding,mode,memberScalarMeans_,heuristicRanking_,inverseRarityRanking_,
				lengthScaledShortAdmission_,familyAdmission_,hitFloors_,index_,blacklist_,source.identityLowerReduction);
		}
		private SearchConfig(SearchConfig source,FamilyReference[] families_,EvaluationMargins margins_,
				FamilyLengthScore binding,FamilyLengthScore.Mode mode,FamilyMemberScalarMeans memberScalarMeans_,
				FamilyHeuristicScore heuristicRanking_,
				FamilyInverseRarityScore inverseRarityRanking_,
				boolean lengthScaledShortAdmission_,FamilyLengthAdmission familyAdmission_,FamilyHitFloorVector hitFloors_,
				DDLFamilySearch.Index index_,DDLPentamerBlacklist blacklist_,float identityLowerReduction_){
			this(source,families_,margins_,binding,mode,memberScalarMeans_,heuristicRanking_,inverseRarityRanking_,
				lengthScaledShortAdmission_,familyAdmission_,hitFloors_,index_,blacklist_,identityLowerReduction_,
				source.lengthPolicy,source.alignmentMetric,source.lengthAmplitude);
		}
		private SearchConfig(SearchConfig source,FamilyReference[] families_,EvaluationMargins margins_,
				FamilyLengthScore binding,FamilyLengthScore.Mode mode,FamilyMemberScalarMeans memberScalarMeans_,
				FamilyHeuristicScore heuristicRanking_,
				FamilyInverseRarityScore inverseRarityRanking_,
				boolean lengthScaledShortAdmission_,FamilyLengthAdmission familyAdmission_,FamilyHitFloorVector hitFloors_,
				DDLFamilySearch.Index index_,DDLPentamerBlacklist blacklist_,float identityLowerReduction_,
				FamilyLengthPolicy lengthPolicy_,FamilyAlignmentMetric alignmentMetric_,float lengthAmplitude_){
			this(source,families_,margins_,binding,mode,memberScalarMeans_,heuristicRanking_,inverseRarityRanking_,
				lengthScaledShortAdmission_,familyAdmission_,hitFloors_,index_,blacklist_,identityLowerReduction_,
				lengthPolicy_,alignmentMetric_,lengthAmplitude_,source.familyGcPenalty,source.familyGcRecipeBound);
		}
		private SearchConfig(SearchConfig source,FamilyReference[] families_,EvaluationMargins margins_,
				FamilyLengthScore binding,FamilyLengthScore.Mode mode,FamilyMemberScalarMeans memberScalarMeans_,
				FamilyHeuristicScore heuristicRanking_,FamilyInverseRarityScore inverseRarityRanking_,
				boolean lengthScaledShortAdmission_,FamilyLengthAdmission familyAdmission_,
				FamilyHitFloorVector hitFloors_,DDLFamilySearch.Index index_,DDLPentamerBlacklist blacklist_,
				float identityLowerReduction_,FamilyLengthPolicy lengthPolicy_,FamilyAlignmentMetric alignmentMetric_,
				float lengthAmplitude_,FamilyGcPenalty familyGcPenalty_,
				boolean familyGcRecipeBound_){
			require(Float.isFinite(identityLowerReduction_) && identityLowerReduction_>=0,
				"Identity lower reduction must be finite and nonnegative");
			index=index_; families=families_; top=source.top; ddlReject=source.ddlReject;
			ddlAccept=source.ddlAccept; geometry=source.geometry; k=source.k; seed=source.seed; width=source.width;
			cosineMode=source.cosineMode; features=source.features; truth=source.truth;
			scalarBounds=source.scalarBounds;
			memberScalarMeans=memberScalarMeans_;
			union=source.union; unionQueries=union ? new java.util.concurrent.ConcurrentHashMap<String,Boolean>() : null;
			scalarGates=source.scalarGates;
			candidateOnly=source.candidateOnly;
			positionMode=source.positionMode; positional=source.positional;
			partitions=source.partitions; overlap=source.overlap;
			ddl.ProteinBlacklistPolicy.requireStored(geometry,blacklist_==null ? null : blacklist_.sha80());
			blacklist=ddl.ProteinBlacklistPolicy.resolve(geometry,blacklist_);
			alignmentMetric=alignmentMetric_;
			identityFloor=source.identityFloor; blosumFloor=source.blosumFloor; floorAudit=source.floorAudit;
			identityLowerReduction=identityLowerReduction_;
			margins=margins_;
			lengthRanking=binding;lengthMode=mode;
			heuristicRanking=heuristicRanking_==null || index==source.index ?
				heuristicRanking_ : heuristicRanking_.rebindIndex(this);
			inverseRarityRanking=inverseRarityRanking_;
			familyGcPenalty=familyGcPenalty_;
			familyGcRecipeBound=familyGcRecipeBound_;
			require(!familyGcRecipeBound || familyGcPenalty!=null,
				"A GC recipe certificate requires the exact loaded penalty binding");
			if(familyGcPenalty!=null){
				familyGcPenalty.requireRoster(families);
				require(heuristicRanking!=null && heuristicRanking.armName().equals("hitfourth_h020_c150") &&
					binding==null && mode==FamilyLengthScore.Mode.OFF && inverseRarityRanking==null &&
					memberScalarMeans==null && lengthPolicy_==null && hitFloors_==null && familyAdmission_==null &&
					!lengthScaledShortAdmission_ && !union && truth==null && !candidateOnly,
					"Derived GC views preserve the released heuristic and cannot silently combine experimental "+
						"ranking/admission policies");
			}
			lengthScaledShortAdmission=lengthScaledShortAdmission_;
			familyAdmission=familyAdmission_;
			hitFloors=hitFloors_;
			lengthPolicy=lengthPolicy_;lengthAmplitude=lengthPolicy==null ? 0 : lengthAmplitude_;
			if(hitFloors!=null){
				require(hitFloors.boundToPolicy(source) && hitFloors.boundTo(families) && families==source.families &&
					margins==source.margins && heuristicRanking==source.heuristicRanking &&
					inverseRarityRanking==source.inverseRarityRanking,
					"Family-floor views preserve the exact roster, ranking and margins; set those policies before attaching floors");
				require(familyAdmission==null && !lengthScaledShortAdmission &&
					geometry.isExact(5) && width==1024 && seed==0 && ddlReject==6 && ddlAccept==1025 && top==50 &&
					positionMode==ProteinPositionMode.FIXED && partitions==4 && overlap==2 &&
					!index.usesPruning() && !index.usesDirectCompletion(),
					"Family floors require the fixed-six K5/P4/O2/1024 parent with exhaustive indexing, top50 and disabled shortcuts");
				require(binding==null && mode==FamilyLengthScore.Mode.OFF && memberScalarMeans==null &&
					heuristicRanking!=null && heuristicRanking.armName().equals("hitfourth_h020_c150") &&
					features==RankingFeatures.COMBINED25 && cosineMode==CosineMode.RANKING && scalarBounds!=null && scalarGates &&
					alignmentMetric.mode==FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY &&
					alignmentMetric.minimumSpan==28 &&
					inverseRarityRanking==null && !union && truth==null && !candidateOnly && !floorAudit,
					"Family floors preserve the released short-family metric, scalars and common pivot; "+
						"diagnostic and fitted-ranking modes are separate");
				for(FamilyReference family:families){
					require(family.members>0 && family.members<=Long.MAX_VALUE/99 &&
						family.residues>=40*family.members && family.residues<=99*family.members,
						"Family floors are released only for retained-member totals within forty through "+
							"ninety-nine residues per member");
				}
			}
			require(familyAdmission==null || !lengthScaledShortAdmission && familyAdmission.boundTo(families),
				"Family-mean and query-length admission cannot be combined; family counts must retain their bound roster");
			if(lengthScaledShortAdmission || familyAdmission!=null){
				require(geometry.isExact(5) && width==1024 && seed==0 && ddlReject==6 && ddlAccept==1025 &&
					top==50 && positionMode==ProteinPositionMode.FIXED && partitions==4 && overlap==2 &&
					!index.usesPruning() && !index.usesDirectCompletion(),
					"Length-based admission requires the released exact-K5/P4/O2/1024 policy, fixed6 anchor, "+
						"top50, disabled shortcut and exhaustive index");
				require(binding==null && mode==FamilyLengthScore.Mode.OFF &&
					(heuristicRanking==null || heuristicRanking.armName().equals("hitfourth_h020_c150")) && memberScalarMeans==null &&
					features==RankingFeatures.COMBINED25 && cosineMode==CosineMode.RANKING && scalarBounds!=null &&
					inverseRarityRanking==null &&
					!union && truth==null && !candidateOnly,
					"Length-based admission preserves old combined25 arrays and permits only the released pivot; "+
						"no calibration, diagnostic or learned policies");
			}
			if(lengthPolicy!=null){
				FamilyLengthAcceptance.validateAmplitude(lengthAmplitude);
				require(alignmentMetric.mode==FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY &&
					alignmentMetric.minimumSpan==FamilyLengthAcceptance.MIN_OVERLAP,
					"Length policy requires explicit transient span28 aligned-query measurement");
				require(binding==null && mode==FamilyLengthScore.Mode.OFF && memberScalarMeans==null &&
					familyAdmission==null && hitFloors==null && !lengthScaledShortAdmission &&
					(margins==null || Double.compare(margins.scale,.985)==0) &&
					(heuristicRanking==null || heuristicRanking.armName().equals("hitfourth_h020_c150")) &&
					identityLowerReduction==0 && inverseRarityRanking==null &&
					!index.usesPruning() && !index.usesDirectCompletion() &&
					ddlAccept>width && cosineMode!=CosineMode.FULL && !floorAudit && !union && truth==null && !candidateOnly,
					"Length policy preserves frozen gates and permits only the released common-pivot ranking");
			}
			require(heuristicRanking==null || binding==null && heuristicRanking.boundTo(this),
				"Derived caller views must preserve the exact heuristic roster and cannot add a learned ranking model");
			monomerSpectra=source.monomerSpectra; monomerNorms=source.monomerNorms;
			centeredSpectra=source.centeredSpectra; centeredNorms=source.centeredNorms; background=source.background;
			if(memberScalarMeans==source.memberScalarMeans){
				combinedSpectra=source.combinedSpectra; combinedNorms=source.combinedNorms;
				require(memberScalarMeans==null || memberScalarMeans.boundTo(this),
					"Derived member-mean views must retain the bound rank/identity/count roster");
			}else{
				combinedSpectra=new double[families.length][COMBINED_WIDTH]; combinedNorms=new double[families.length];
				fillCombined();
			}
			require(inverseRarityRanking==null || binding==null && inverseRarityRanking.boundTo(this),
				"Derived caller views must preserve the exact inverse-rarity roster and cannot add a learned ranking model");
		}
		/** The fitted metric contract permits fixed-coordinate P4/O2 profiles and their plain sibling. */
		boolean releasedMetricProfileGeometry(){
			return positionMode==ProteinPositionMode.PLAIN ||
				positionMode.fixedCoordinates() && partitions==4 && overlap==2;
		}
		/** Build immutable normalized family vectors once; query scalar arithmetic is unchanged. */
		private void fillCombined(){
			require(memberScalarMeans==null || scalarBounds!=null && features==RankingFeatures.COMBINED25 &&
				cosineMode==CosineMode.RANKING && memberScalarMeans.boundTo(this),
				"Member means require combined25 ranking and the exact bound family roster; "+
					"dimer gates must not consume these vectors");
			if(combinedSpectra==null){return;}
			assert(features==RankingFeatures.COMBINED25 && cosineMode==CosineMode.RANKING) :
				"Combined vectors must never reach dimer-calibrated cosine gates";
			final double[] rawMeans=memberScalarMeans==null ? null : new double[ConsensusScalars.SCALAR_COUNT];
			for(int i=0;i<families.length;i++){
				if(memberScalarMeans==null){
					// Historical comparison/default path; explicit member-mean views repair the scalar coordinates.
					scalarBounds.combine(monomerSpectra[i],families[i].consensusRaw,false,combinedSpectra[i]);
				}else{
					memberScalarMeans.copyMeans(i,rawMeans);
					scalarBounds.combineComputed(monomerSpectra[i],false,rawMeans,ConsensusScalars.ALL_DEFINED,combinedSpectra[i]);
				}
				combinedNorms[i]=combinedNorm(combinedSpectra[i]);
			}
		}
		SearchConfig withFamilies(FamilyReference[] families_){
			require(hitFloors==null,"Family-floor views cannot replace their frozen reference gates or membership");
			require(lengthPolicy==null,"A length-penalty deployment view cannot replace its frozen reference gates");
			require(familyAdmission==null,
				"Family-mean admission must not silently retain floors while replacing reference membership");
			require(!lengthScaledShortAdmission,"The length-scaled admission study does not authorize family recalibration");
			require(families_!=null && families_.length==families.length,"Updated classifier roster width changed");
			require(margins==null,"Calibration cannot silently discard an evaluation-only margin policy");
			require(lengthRanking==null,"Frozen learned ranking does not authorize per-family calibration updates");
			require(inverseRarityRanking==null,"Inverse-rarity ranking does not authorize per-family calibration updates");
			require(heuristicRanking==null || !alignmentMetric.same(FamilyAlignmentMetric.LEGACY) &&
				heuristicRanking.armName().equals("hitfourth_h020_c150"),
				"Only metric-bound common-pivot calibration may retain a frozen caller heuristic across gate updates");
			return new SearchConfig(this,families_,null);
		}
		final DDLFamilySearch.Index index;
		final FamilyReference[] families;
		final ProteinKmerGeometry geometry;
		final int top, ddlReject, ddlAccept, k, width, partitions, overlap;
		final DDLPentamerBlacklist blacklist;
		final FamilyAlignmentMetric alignmentMetric;
		final long seed;
		final CosineMode cosineMode;
		final float[][] monomerSpectra,centeredSpectra;
		final float[] monomerNorms,centeredNorms,background;
		final ScalarBounds scalarBounds;
		final FamilyMemberScalarMeans memberScalarMeans;
		final double[][] combinedSpectra;
		final double[] combinedNorms;
		final RankingFeatures features;
		final HashMap<String,Truth> truth;
		final boolean union,scalarGates,floorAudit,candidateOnly,positional;
		final ProteinPositionMode positionMode;
		final float identityFloor,blosumFloor,identityLowerReduction;
		final EvaluationMargins margins;
		final FamilyLengthScore lengthRanking;
		final FamilyLengthScore.Mode lengthMode;
		final FamilyHeuristicScore heuristicRanking;
		final FamilyGcPenalty familyGcPenalty;
		final boolean familyGcRecipeBound;
		final FamilyInverseRarityScore inverseRarityRanking;
		final boolean lengthScaledShortAdmission;
		final FamilyLengthAdmission familyAdmission;
		final FamilyHitFloorVector hitFloors;
		final FamilyLengthPolicy lengthPolicy;
		final float lengthAmplitude;
		final java.util.concurrent.ConcurrentHashMap<String,Boolean> unionQueries;
	}

	/** Worker-local raw maxima; merge is order-independent. No baseline counts or margins here. */
	private static final class UnionStats{
		UnionStats(int size){
			counts=new long[size]; identity=new double[size]; blosum=new double[size];
		}
		void observe(int ordinal,FamilyCalibration.Metrics metrics){
			if(!Double.isFinite(metrics.identityPct) || metrics.identityPct<0 || metrics.identityPct>100 ||
					!Double.isFinite(metrics.blosum)){
				throw new IllegalArgumentException("Nonfinite or invalid raw union alignment metric");
			}
			if(counts[ordinal]==0){identity[ordinal]=metrics.identityPct; blosum[ordinal]=metrics.blosum;}
			else{
				identity[ordinal]=Math.max(identity[ordinal],metrics.identityPct);
				blosum[ordinal]=Math.max(blosum[ordinal],metrics.blosum);
			}
			counts[ordinal]++;
		}
		void add(UnionStats other){
			for(int i=0; i<counts.length; i++){
				if(other.counts[i]==0){continue;}
				if(counts[i]==0){identity[i]=other.identity[i]; blosum[i]=other.blosum[i];}
				else{identity[i]=Math.max(identity[i],other.identity[i]); blosum[i]=Math.max(blosum[i],other.blosum[i]);}
				counts[i]+=other.counts[i];
			}
		}
		void write(Path out,FamilyReference[] families,long reads,long pairs){
			long sum=0;
			for(long count : counts){
				if(count<0 || count>reads || count>100000){throw new IllegalStateException("Invalid per-family union count");}
				sum+=count;
			}
			if(sum!=pairs){throw new IllegalStateException("Union pair count does not reconcile");}
			// This final table is serial; the ordered batch writer requires add(buffer,batchId), not print.
			final ByteStreamWriter output=new ByteStreamWriter(out.toString(),false,false,false); output.start();
			try{
				output.print("roster_rank\tfamily_id\taligned_count\tidentity_pct_upper\tblosum_a_sq_upper\n");
				for(int i=0; i<families.length; i++){
					if(counts[i]==0){continue;}
					output.print(new ByteBuilder().append(families[i].rank).tab().append(families[i].id).tab().append(counts[i]).tab()
						.append(Double.toString(identity[i])).tab().append(Double.toString(blosum[i])).nl());
				}
			}finally{if(output.poisonAndWait()){throw new IllegalStateException("Union supplemental writer failed");}}
		}
		final long[] counts;
		final double[] identity,blosum;
	}

	static final class Counters{
		void add(final Counters c){
			reads+=c.reads; residues+=c.residues; assigned+=c.assigned; rejected+=c.rejected;
			if(c.stageWorkerNanos!=0){
				stageWorkerNanos+=c.stageWorkerNanos; stageIndexNanos+=c.stageIndexNanos;
				stageShortlistNanos+=c.stageShortlistNanos; stageAlignmentNanos+=c.stageAlignmentNanos;
			}
			queryDdlReject+=c.queryDdlReject; queryDdlAccept+=c.queryDdlAccept;
			queryCosineAccept+=c.queryCosineAccept; queryIdentityAccept+=c.queryIdentityAccept;
			queryBlosumUpperAccept+=c.queryBlosumUpperAccept; queryBlosumLowerAccept+=c.queryBlosumLowerAccept;
			queryExhausted+=c.queryExhausted;
			queryExhaustedAfterDp+=c.queryExhaustedAfterDp;
			queryDdlEmpty+=c.queryDdlEmpty;
			ddlActiveAllSum+=c.ddlActiveAllSum; ddlActivePostDdlSum+=c.ddlActivePostDdlSum;
			postDdlQueries+=c.postDdlQueries;
			for(int i=0; i<winnerCount.length; i++){
				winnerCount[i]+=c.winnerCount[i]; winnerDdlPositionSum[i]+=c.winnerDdlPositionSum[i];
				winnerDpPositionSum[i]+=c.winnerDpPositionSum[i];
			}
			queryCosineAllPruned+=c.queryCosineAllPruned; queryReachingDp+=c.queryReachingDp;
			queryScalarAllPruned+=c.queryScalarAllPruned; scalarComputations+=c.scalarComputations;
			queryGcAllPruned+=c.queryGcAllPruned;gcPairsChecked+=c.gcPairsChecked;gcPairsRemoved+=c.gcPairsRemoved;
			scalarPairsChecked+=c.scalarPairsChecked;
			scalarPairsRemoved+=c.scalarPairsRemoved; scalarPairsKept+=c.scalarPairsKept;
			for(int i=0;i<5;i++){
				scalarFailed[i]+=c.scalarFailed[i]; scalarDisabled[i]+=c.scalarDisabled[i];
				scalarQueryUndefined[i]+=c.scalarQueryUndefined[i]; scalarTested[i]+=c.scalarTested[i];
			}
			candidateDdlLow+=c.candidateDdlLow; candidateDdlHigh+=c.candidateDdlHigh;
			candidateCosineUpper+=c.candidateCosineUpper; candidateCosineLow+=c.candidateCosineLow;
			candidateCosineSurvive+=c.candidateCosineSurvive; candidateDp+=c.candidateDp;
			candidateUnion+=c.candidateUnion;
			candidateRows+=c.candidateRows; queryCandidates+=c.queryCandidates;
			if(c.unionStats!=null){unionStats.add(c.unionStats);}
			candidateIdentityUpper+=c.candidateIdentityUpper;
			candidateIdentityLow+=c.candidateIdentityLow; candidateBlosumUpper+=c.candidateBlosumUpper;
			candidateBlosumLower+=c.candidateBlosumLower; candidateBlosumReject+=c.candidateBlosumReject;
			candidateAlignmentUnavailable+=c.candidateAlignmentUnavailable; candidateAlignmentFailed+=c.candidateAlignmentFailed;
			batches+=c.batches; processNanos+=c.processNanos;
			rankFamilies+=c.rankFamilies; rankCosines+=c.rankCosines;
			rankZeroNorm+=c.rankZeroNorm; rankBoundSkips+=c.rankBoundSkips;
			gateCosines+=c.gateCosines;
			diagnosticQueries+=c.diagnosticQueries; diagnosticFamilies+=c.diagnosticFamilies;
			for(int i=0; i<RANK_MODES.length; i++){
				for(int j=0; j<TRUTH_STATUS.length; j++){truthStatus[i][j]+=c.truthStatus[i][j];}
			}
		}
		long queryTerminals(){return queryDdlReject+queryDdlAccept+queryCosineAccept+queryIdentityAccept+
			queryBlosumUpperAccept+queryBlosumLowerAccept+queryExhausted+queryCandidates;}
		long candidateVisits(){
			return candidateDdlLow+candidateDdlHigh+candidateCosineUpper+candidateCosineLow+
				scalarPairsRemoved+candidateDp;
		}
		long dpPartition(){return candidateIdentityUpper+candidateIdentityLow+candidateBlosumUpper+
			candidateBlosumLower+candidateBlosumReject+candidateAlignmentUnavailable+candidateAlignmentFailed+candidateUnion;}
		/** Original DDL order is descending matches, then roster rank, before top/tie trimming.
		 * DP order is the 1-based actual attempt after cosine filtering/reranking; zero means no DP. */
		void winner(final int stage,final int ddlPosition,final int dpPosition){
			assert(stage>=0 && stage<winnerCount.length && ddlPosition>=1 &&
				(stage<2 ? dpPosition==0 : dpPosition>=1)) :
				"Winning positions must describe an actual stage visit, not the family roster rank: "+
				stage+"/"+ddlPosition+"/"+dpPosition;
			winnerCount[stage]++; winnerDdlPositionSum[stage]+=ddlPosition;
			winnerDpPositionSum[stage]+=dpPosition;
		}
		void printAccounting(int ddlReject){
			validate();
			if(gcPairsChecked>0){System.out.println("FAMILY_TAIL_CASCADE_GC_COUNTS pairs_checked="+gcPairsChecked+
				" pairs_removed="+gcPairsRemoved+" all_pruned_queries="+queryGcAllPruned);}
			System.out.println("FAMILY_TAIL_CASCADE_TERMINALS ddl_accept="+queryDdlAccept+
				" ddl_reject="+(queryDdlReject+queryDdlEmpty)+" cosine_accept="+queryCosineAccept+
				" identity_accept="+queryIdentityAccept+" blosum_accept="+(queryBlosumUpperAccept+queryBlosumLowerAccept)+
				" exhausted_after_dp="+queryExhaustedAfterDp+" cosine_exhausted="+queryCosineAllPruned+
				" scalar_exhausted="+queryScalarAllPruned+
				(gcPairsChecked>0 ? " gc_exhausted="+queryGcAllPruned : ""));
			System.out.println("FAMILY_TAIL_CASCADE_DDL_COUNTS threshold="+ddlReject+" untrimmed_all_sum="+ddlActiveAllSum+
				" all_queries="+reads+" untrimmed_post_ddl_sum="+ddlActivePostDdlSum+" post_ddl_queries="+postDdlQueries);
			for(int i=0; i<WINNER_STAGES.length; i++){
				System.out.println("FAMILY_TAIL_CASCADE_WINNER_"+WINNER_STAGES[i]+" count="+winnerCount[i]+
					" ddl_position_sum="+winnerDdlPositionSum[i]+" ddl_position_count="+winnerCount[i]+
					" dp_attempt_sum="+winnerDpPositionSum[i]+" dp_attempt_count="+(i<2 ? 0 : winnerCount[i]));
			}
		}
		void printTruth(boolean combined){
			for(int mode=0; mode<(combined ? 5 : 4); mode++){
				long total=0; for(long count : truthStatus[mode]){total+=count;}
				if(total!=reads || total!=diagnosticQueries){
					throw new IllegalStateException("Truth statuses must partition every query");
				}
				final StringBuilder line=new StringBuilder("FAMILY_TAIL_CASCADE_TRUTH_COUNTS mode=").append(RANK_MODES[mode]);
				line.append(" queries=").append(total).append(" mapped_roster_queries=")
					.append(total-truthStatus[mode][0]-truthStatus[mode][1]);
				for(int status=0; status<TRUTH_STATUS.length; status++){
					line.append(' ').append(TRUTH_STATUS[status]).append('=').append(truthStatus[mode][status]);
				}
				System.out.println(line);
			}
			System.out.println("FAMILY_TAIL_CASCADE_DIAGNOSTIC_WORK queries="+diagnosticQueries+
				" window_families="+diagnosticFamilies+" dimer_evaluations="+diagnosticFamilies+
				" monomer_evaluations="+diagnosticFamilies+
				" centered_monomer_evaluations="+diagnosticFamilies+
				(combined ? " combined25_evaluations="+diagnosticFamilies : ""));
		}
		void validate(){
			if(assigned+rejected!=reads || queryTerminals()!=reads || dpPartition()!=candidateDp ||
					queryExhausted!=queryExhaustedAfterDp+queryCosineAllPruned+queryScalarAllPruned+queryGcAllPruned+queryDdlEmpty ||
					rejected!=queryDdlReject+queryExhausted+queryCandidates ||
					postDdlQueries!=reads-queryDdlReject-queryDdlEmpty-queryDdlAccept ||
					queryReachingDp!=queryIdentityAccept+queryBlosumUpperAccept+queryBlosumLowerAccept+queryExhaustedAfterDp ||
					ddlActiveAllSum<ddlActivePostDdlSum || ddlActivePostDdlSum<0 ||
					winnerCount[0]!=queryDdlAccept || winnerCount[1]!=queryCosineAccept ||
					winnerCount[2]!=queryIdentityAccept || winnerCount[3]!=queryBlosumUpperAccept+queryBlosumLowerAccept ||
					queryCosineAllPruned>queryExhausted || queryReachingDp>candidateDp ||
					queryDdlReject+queryDdlAccept+queryCosineAccept+queryCosineAllPruned+
					queryScalarAllPruned+queryGcAllPruned+queryReachingDp>reads ||
					gcPairsRemoved>gcPairsChecked || queryGcAllPruned>gcPairsRemoved ||
					scalarPairsKept+scalarPairsRemoved!=scalarPairsChecked || scalarComputations>reads ||
					candidateCosineUpper+candidateCosineLow+candidateCosineSurvive!=gateCosines ||
					rankCosines+rankZeroNorm+rankBoundSkips!=rankFamilies){
				throw new IllegalStateException("FamilyTailCascade counters do not reconcile");
			}
			long failures=0;
			for(int i=0;i<5;i++){
				if(scalarDisabled[i]+scalarQueryUndefined[i]+scalarTested[i]!=scalarPairsChecked ||
						scalarFailed[i]>scalarTested[i]){
					throw new IllegalStateException("Scalar gate applicability counters do not reconcile");
				}
				failures+=scalarFailed[i];
			}
			if(failures<scalarPairsRemoved || failures>5*scalarPairsRemoved){
				throw new IllegalStateException("Scalar failure overlap does not reconcile");
			}
			for(int i=0; i<winnerCount.length; i++){
				if(winnerDdlPositionSum[i]<winnerCount[i] ||
						(i<2 ? winnerDpPositionSum[i]!=0 : winnerDpPositionSum[i]<winnerCount[i])){
					throw new IllegalStateException("Winning positions do not reconcile for "+WINNER_STAGES[i]);
				}
			}
		}
		void printScalarAccounting(boolean enabled,double z){
			System.out.println("FAMILY_TAIL_CASCADE_SCALAR_COUNTS enabled="+(enabled ? 1 : 0)+
				" pairs_checked="+scalarPairsChecked+
				" pairs_removed="+scalarPairsRemoved+" pairs_kept="+scalarPairsKept+" all_pruned_queries="+queryScalarAllPruned+
				" query_computations="+scalarComputations);
			System.out.println("FAMILY_TAIL_CASCADE_SCALAR_RECIPE z="+z+
				" placement=ranked_top_cap_before_dp backfill=false failure_counts=overlapping_per_scalar");
			for(int i=0;i<5;i++){
				System.out.println("FAMILY_TAIL_CASCADE_SCALAR_"+SCALAR_NAMES[i]+
					" failed="+scalarFailed[i]+" tested="+scalarTested[i]+
					" disabled="+scalarDisabled[i]+" query_undefined="+scalarQueryUndefined[i]);
			}
		}
		void printStageCpu(){
			final long search=stageIndexNanos+stageShortlistNanos+stageAlignmentNanos;
			if(stageIndexNanos<0 || stageShortlistNanos<0 || stageAlignmentNanos<0 || stageWorkerNanos<search){
				throw new IllegalStateException("Stage CPU intervals do not partition worker CPU");
			}
			final String[] names={"index","shortlist","alignment","other_worker"};
			final long[] values={stageIndexNanos,stageShortlistNanos,stageAlignmentNanos,stageWorkerNanos-search};
			for(int i=0;i<names.length;i++){
				System.out.println("FAMILY_TAIL_STAGE_CPU stage="+names[i]+" nanos="+values[i]+
					" seconds="+seconds(values[i])+" percent_worker_cpu="+(stageWorkerNanos==0 ? 0 : 100.0*values[i]/stageWorkerNanos)+
					" percent_search_cpu="+(i==3 || search==0 ? 0 : 100.0*values[i]/search));
			}
			System.out.println("FAMILY_TAIL_STAGE_CPU_TOTAL worker_nanos="+stageWorkerNanos+
				" search_nanos="+search+
				" clock=current_thread_cpu scope=search_workers_only separate_loader_writer_jit_gc_threads_excluded=true"+
				" worker_allocation_and_gc_assist_included=true"+
				" index=sketch_and_index_scan shortlist=prepare_except_index alignment=DP_and_acceptance_loop");
		}
		long stageWorkerNanos,stageIndexNanos,stageShortlistNanos,stageAlignmentNanos;
		long reads,residues,assigned,rejected;
		long queryDdlReject,queryDdlAccept,queryCosineAccept,queryIdentityAccept,queryCosineAllPruned,queryReachingDp;
		long queryBlosumUpperAccept,queryBlosumLowerAccept,queryExhausted;
		long queryDdlEmpty,queryExhaustedAfterDp,postDdlQueries;
		long candidateRows,queryCandidates;
		long queryScalarAllPruned,scalarComputations,scalarPairsChecked,scalarPairsRemoved,scalarPairsKept;
		long queryGcAllPruned,gcPairsChecked,gcPairsRemoved;
		final long[] scalarFailed=new long[5],scalarDisabled=new long[5],scalarQueryUndefined=new long[5],
			scalarTested=new long[5];
		long ddlActiveAllSum,ddlActivePostDdlSum;
		final long[] winnerCount=new long[4],winnerDdlPositionSum=new long[4],winnerDpPositionSum=new long[4];
		static final String[] WINNER_STAGES={"DDL","COSINE","IDENTITY","BLOSUM"};
		long candidateDdlLow,candidateDdlHigh,candidateCosineUpper,candidateCosineLow,candidateCosineSurvive,candidateDp;
		long candidateIdentityUpper,candidateIdentityLow,candidateBlosumUpper,candidateBlosumLower,candidateBlosumReject;
		long candidateAlignmentUnavailable,candidateAlignmentFailed;
		long candidateUnion;
		UnionStats unionStats;
		long batches,processNanos,rankFamilies,rankCosines,rankZeroNorm,rankBoundSkips,gateCosines;
		long diagnosticQueries,diagnosticFamilies;
		final long[][] truthStatus=new long[RANK_MODES.length][TRUTH_STATUS.length];
	}

	/** Opt-in thread CPU clock. Never instantiated or queried on the normal path. */
	private static final class StageCpu{
		StageCpu(){
			bean=java.lang.management.ManagementFactory.getThreadMXBean();
			if(!bean.isCurrentThreadCpuTimeSupported()){throw new IllegalStateException("Thread CPU timing unsupported");}
			if(!bean.isThreadCpuTimeEnabled()){bean.setThreadCpuTimeEnabled(true);}
			now();
		}
		long now(){
			final long value=bean.getCurrentThreadCpuTime();
			if(value<0){throw new IllegalStateException("Thread CPU timing unavailable");}
			return value;
		}
		final java.lang.management.ThreadMXBean bean;
		long index,shortlist,alignment;
	}
	/** Timing wrappers preserve every production branch, including early returns. */
	private static final class TimedSearcher extends Searcher{
		TimedSearcher(SearchConfig config,StageCpu timing_){super(config); timing=timing_;}
		@Override void scoreIndex(byte[] raw){
			final long start=timing.now();
			try{super.scoreIndex(raw);}finally{timing.index+=timing.now()-start;}
		}
		@Override boolean prepare(byte[] raw,ProteinSequence query,Result result,Counters counters,
				ByteBuilder pairs,boolean reuse){
			final long start=timing.now(),beforeIndex=timing.index;
			try{return super.prepare(raw,query,result,counters,pairs,reuse);}
			finally{timing.shortlist+=timing.now()-start-(timing.index-beforeIndex);}
		}
		@Override void alignPrepared(ProteinSequence query,Result result,Counters counters,ByteBuilder pairs){
			final long start=timing.now();
			try{super.alignPrepared(query,result,counters,pairs);}finally{timing.alignment+=timing.now()-start;}
		}
		final StageCpu timing;
	}

	static class Searcher{
		Searcher(final SearchConfig config){
			this(config,false);
		}
		/** Diagnostic workers may retain the strict ranked cap before scalar filtering.
		 * Ordinary workers keep their existing allocation and search behavior. */
		Searcher(final SearchConfig config,final boolean captureTop){
			this(config,captureTop,null);
		}
		Searcher(final SearchConfig config,final boolean captureTop,final FamilyMatchInformation.Binding information_){
			information=information_;
			index=config.index; families=config.families;
			indexPruned=index.usesPruning(); indexDirect=index.usesDirectCompletion(); indexOrder=index.scanOrder();
			keyOrder=indexPruned && indexOrder==0 ? new int[config.width] : null;
			occupied=indexPruned && indexOrder==0 ? new int[config.width] : null;
			keyOrigins=indexPruned && (indexOrder==1 || indexOrder==2) ? new ddl.DDLQueryOrigins(config.width) : null;
			truth=config.truth; monomerSpectra=config.monomerSpectra; monomerNorms=config.monomerNorms;
			centeredSpectra=config.centeredSpectra; centeredNorms=config.centeredNorms; background=config.background;
			centered=config.features==RankingFeatures.CENTERED_MONOMER;
			scalarBounds=config.scalarBounds; combinedSpectra=config.combinedSpectra; combinedNorms=config.combinedNorms;
			combined=scalarBounds!=null;
			union=config.union;
			candidateOnly=config.candidateOnly;
			positionMode=config.positionMode; positional=config.positional;
			partitions=config.partitions; overlap=config.overlap;
			blacklist=config.blacklist;
			legacyAlignment=config.alignmentMetric.same(FamilyAlignmentMetric.LEGACY);
			overhangOnly=config.alignmentMetric.mode==FamilyAlignmentMetric.Mode.PIN_RIGHT_OVERHANG_ONLY_ALIGNED_QUERY;
			scalarGates=config.scalarGates;
			identityFloor=config.identityFloor; blosumFloor=config.blosumFloor;
			identityLowerReduction=config.identityLowerReduction;
			margins=config.margins;
			lengthMode=config.lengthMode;lengthRanker=config.lengthRanking==null ? null : config.lengthRanking.newWorker();
			heuristicRanker=config.heuristicRanking;
			familyGcPenalty=config.familyGcPenalty;
			inverseRarityRanker=config.inverseRarityRanking;
			lengthScaledShortAdmission=config.lengthScaledShortAdmission;
			familyAdmission=config.familyAdmission;
			hitFloors=config.hitFloors;
			lengthPolicy=config.lengthPolicy;
			lengthAmplitude=config.lengthAmplitude;
			rankedTop=lengthRanker==null && !captureTop ? null : new int[config.top];
			geometry=config.geometry;
			final boolean monomer=config.features==RankingFeatures.MONOMER;
			useCosine=config.cosineMode!=CosineMode.OFF; gateCosine=config.cosineMode==CosineMode.FULL;
			top=config.top; ddlReject=config.ddlReject; ddlAccept=config.ddlAccept; k=config.k; seed=config.seed;
			sketch=new char[config.width]; comparisons=new Comparison[families.length];
			for(int i=0; i<comparisons.length; i++){
				comparisons[i]=new Comparison(i,families[i],
					centered ? centeredSpectra[i] : monomer ? monomerSpectra[i] : families[i].spectrum,
					centered ? centeredNorms[i] : monomer ? monomerNorms[i] : families[i].spectrumNorm);
			}
			work=new DDLFamilySearch.Workspace(index.families()); queryCounts=new ProteinSpectrum.Spectrum(2);
			inverseRarityScores=inverseRarityRanker==null ? null : new double[index.families()];
			inverseRarityScratch=inverseRarityRanker!=null && inverseRarityRanker.needsQueryMedian() ?
				new double[index.families()] : null;
			querySpectrum=new float[DIMER_WIDTH];
			queryMonomer=config.monomerSpectra==null ? null : new float[MONOMER_WIDTH];
			queryCentered=centeredSpectra==null ? null : new float[MONOMER_WIDTH];
			queryCombined=combined ? new double[COMBINED_WIDTH] : null;
			queryRankingSpectrum=centered ? queryCentered : monomer ? queryMonomer : querySpectrum;
			metrics=new FamilyCalibration.Metrics();
			alignment=config.alignmentMetric.newWorker(information==null ? null : information.background);
			observation=new FamilyAlignmentMetric.Observation();
		}
		/** Shared production candidate preparation. Grid reuse changes only the reject threshold;
		 * the scored DDL workspace is still for this exact query. No post-rerank pool is reused. */
		boolean prepare(final byte[] raw, final ProteinSequence query, final Result result, final Counters counters,
				final ByteBuilder pairs, final boolean reuseScores){
			return prepare(raw,query,result,counters,pairs,reuseScores,lengthMode);
		}
		boolean prepare(final byte[] raw,final ProteinSequence query,final Result result,final Counters counters,
				final ByteBuilder pairs,final boolean reuseScores,final FamilyLengthScore.Mode mode){
			require(familyGcPenalty==null || queryGcInitialized,
				"GC ranking requires query-specific availability before preparation");
			assert(!reuseScores || !union && truth==null) : "Grid score reuse must not invoke union or diagnostic side effects";
			require(mode!=null && (mode==FamilyLengthScore.Mode.OFF || lengthRanker!=null),
				"An enabled ranking mode needs its bound worker model");
			if(hitFloors!=null){
				require(raw.length>=40 && raw.length<=99,
					"Family-specific short floors require an actual normalized query length within forty through "+
						"ninety-nine residues");
			}
			if(heuristicRanker!=null){
				require(mode==FamilyLengthScore.Mode.OFF,"A caller heuristic cannot be mixed with a learned ranking arm");
				heuristicRanker.validateQueryLength(raw.length);
			}
			if(inverseRarityRanker!=null){
				require(mode==FamilyLengthScore.Mode.OFF,
					"An inverse-rarity caller heuristic cannot be mixed with a learned ranking arm");
				inverseRarityRanker.validateQueryLength(raw.length);
			}
			if(lengthRanker!=null){
				require(raw.length>=100 && raw.length<=199,
					"The frozen length model supports100..199-aa queries only, including rejected queries");
				if(!reuseScores){
					require(lengthGeneration<Long.MAX_VALUE,"Prediction cache query generation overflow");
					lengthGeneration++;
				}
			}
			rankedTopSize=0;
			result.clear();
			final int admissionFloor=hitFloors!=null ? hitFloors.minimum : familyAdmission!=null ? familyAdmission.minimum :
				lengthScaledShortAdmission ? ShortQueryAdmission.rejectFloor(query.length()) : ddlReject;
			result.admissionFloor=admissionFloor;
			result.admissionFloorMax=hitFloors!=null ? hitFloors.maximum :
				familyAdmission==null ? admissionFloor : familyAdmission.maximum;
			queryScalarsReady=false;
			for(Comparison c : dirty){c.reset();}
			dirty.clear();
			shortlist.clear();
			if(!reuseScores){
				scoreIndex(raw);
			}
			if(rawHitObserver!=null){
				require(!reuseScores,"A raw-hit observation counts each query once, never a reused grid score");
				rawHitObserver.accept(work.touched,work.scores,work.size);
			}
			int maxDdl=0, ddlActive=0;
			for(int i=0; i<work.size; i++){
				final int family=work.touched[i], score=work.scores[family];
				final int floor=rawHitFloor(family,admissionFloor);
				if(candidateObserver!=null){candidateObserver.rawCandidate(family,score,floor);}
				maxDdl=Math.max(maxDdl,score);
				if(score>=floor){ddlActive++;}
			}
			counters.ddlActiveAllSum+=ddlActive;
			result.ddlActive=ddlActive;
			final boolean belowAdmission=familyAdmission==null && hitFloors==null ? maxDdl<admissionFloor : ddlActive==0;
			if(belowAdmission && inverseRarityRanker!=null && candidateObserver!=null){fullAdmittedFeatures(query,scalarGates);}
			if(belowAdmission && truth==null){result.terminal="ddl_reject"; counters.queryDdlReject++; return false;}
			for(int i=0; i<work.size; i++){
				final int family=work.touched[i], ddl=work.scores[family];
				if(ddl>=rawHitFloor(family,admissionFloor)){
					final Comparison c=comparisons[family];
					c.ddlHits=ddl; shortlist.add(c); dirty.add(c);
				}
			}
			Collections.sort(shortlist,DDL_ORDER);
			if(truth!=null){diagnose(raw,query.id,result,counters);}
			if(belowAdmission){result.terminal="ddl_reject"; counters.queryDdlReject++; return false;}
			// With ddlreject=0, an empty touched roster still means no DDL candidate, not DP exhaustion.
			if(shortlist.isEmpty()){
				if(inverseRarityRanker!=null && candidateObserver!=null){fullAdmittedFeatures(query,scalarGates);}
				counters.queryDdlEmpty++; counters.queryExhausted++; return false;
			}
			Comparison best=shortlist.get(0);
			if(best.ddlHits>=ddlAccept){
				best.ddlPosition=1;
				result.set(best); counters.candidateDdlHigh++; counters.queryDdlAccept++;
				result.terminal="ddl_accept";
				counters.winner(0,1,0); return false;
			}
			counters.postDdlQueries++;
			counters.ddlActivePostDdlSum+=ddlActive;
			final boolean learned=mode!=FamilyLengthScore.Mode.OFF;
			final boolean heuristic=heuristicRanker!=null || inverseRarityRanker!=null;
			int limit=learned || heuristic ? shortlist.size() : Math.min(top,shortlist.size());
			if(useCosine){
				if(shortlist.size()>limit){
					final int cutoff=shortlist.get(limit-1).ddlHits;
					while(limit<shortlist.size() && shortlist.get(limit).ddlHits==cutoff){limit++;}
				}
				queryCounts.clear(); queryCounts.add(raw); l1(queryCounts,querySpectrum);
				result.cosineWindow=limit;
				if(combined || queryRankingSpectrum!=querySpectrum){dimerRowSums(querySpectrum,queryMonomer);}
				if(centered){centerMonomers(queryMonomer,background,queryCentered);}
				final float queryNorm=norm(queryRankingSpectrum);
				if(combined && truth==null){
					prepareQueryScalars(raw,counters);
					scalarBounds.combineComputed(queryMonomer,true,queryScalars,queryScalarMask,queryCombined);
				}
				final double queryCombinedNorm=combined ? combinedNorm(queryCombined) : 0;
				if((heuristicRanker!=null && heuristicRanker.usesHydrophobicity()) || inverseRarityRanker!=null){
					require((queryScalarMask&1)!=0,
						"An enabled caller H term requires a measured query KD mean; undefined coordinates cannot be imputed");
				}
				int survivors=0;
				for(int i=0; i<limit; i++){
					final Comparison c=shortlist.get(i); counters.rankFamilies++;
					c.ddlPosition=i+1;
					if(combined){
						c.combinedCosine=combinedCosine(queryCombined,combinedSpectra[c.ordinal],
							queryCombinedNorm,combinedNorms[c.ordinal]);
						c.combinedMeasured=true;
						if(learned){
							final float f=lengthRanker.predict(raw.length,c.ordinal,c.ddlHits,c.combinedCosine,lengthGeneration);
							c.lengthScore=mode==FamilyLengthScore.Mode.MULTIPLY ? c.ddlHits*(double)f : c.ddlHits+(double)f;
						}
						if(inverseRarityRanker!=null){
							c.inverseLengthFactor=inverseRarityRanker.lengthFactor(raw.length,c.ordinal);
							c.inverseHydrophobicity=inverseRarityRanker.hydrophobicitySimilarity(queryScalars[0],c.ordinal);
						}else if(heuristicRanker!=null){
							c.heuristicScore=heuristicRanker.score(raw.length,c.ordinal,c.ddlHits,c.combinedCosine,queryScalars[0]);
							if(familyGcPenalty!=null){
								c.heuristicScore=familyGcPenalty.score(c.heuristicScore,c.ordinal,queryGcAvailable,queryGc);
							}
							c.heuristicMeasured=true;
						}
					}else{
						c.cosine=cosine(queryRankingSpectrum,c.rankingSpectrum,queryNorm,c.rankingNorm);
						c.cosineMeasured=true;
					}
					if(combined ? queryCombinedNorm==0 || combinedNorms[c.ordinal]==0 :
							queryNorm==0 || c.rankingNorm==0){counters.rankZeroNorm++;}
					else{counters.rankCosines++;}
					// Ranking-only retains the whole DDL top+ties window; no gate was evaluated.
					if(!gateCosine){shortlist.set(survivors++,c); continue;}
					counters.gateCosines++;
					if(c.cosine>c.family.cosineAccept){
						result.set(c); counters.candidateCosineUpper++; counters.queryCosineAccept++;
						result.terminal="cosine_accept";
						counters.winner(1,c.ddlPosition,0); return false;
					}
					if(c.cosine<c.family.cosineLower){counters.candidateCosineLow++;}
					else{
						counters.candidateCosineSurvive++;
						shortlist.set(survivors++,c);
					}
				}
				while(shortlist.size()>survivors){shortlist.remove(shortlist.size()-1);}
				if(survivors==0){
					result.terminal="cosine_reject";
					counters.queryCosineAllPruned++;
					counters.queryExhausted++;
					return false;
				}
			}else{
				while(shortlist.size()>limit){shortlist.remove(shortlist.size()-1);}
				for(int i=0; i<shortlist.size(); i++){shortlist.get(i).ddlPosition=i+1;}
			}
			if(inverseRarityRanker!=null){scoreInverseRarityRanking(shortlist.size());}
			if(inverseRarityRanker!=null && candidateObserver!=null){
				fullAdmittedFeatures(query,scalarGates);
			}
			Collections.sort(shortlist,heuristic ? HEURISTIC_RANKING_ORDER :
				learned ? LENGTH_RANKING_ORDER : combined ? COMBINED_RANKING_ORDER : RANKING_ORDER);
			result.dpShortlist=Math.min(top,shortlist.size());
			result.shortlistCandidates=result.dpShortlist;
			// Stamp only the strict top cap; scalar/GC compaction preserves these objects.
			for(int i=0;i<result.dpShortlist;i++){shortlist.get(i).rankedPosition=i+1;}
			if(candidateObserver!=null){
				for(int i=0;i<shortlist.size();i++){
					final Comparison c=shortlist.get(i);
					candidateObserver.rankedCandidate(c.ordinal,c.ddlHits,c.ddlPosition,i+1,i<result.dpShortlist);
				}
			}
			if(rankedTop!=null){
				rankedTopSize=result.dpShortlist;
				for(int i=0;i<rankedTopSize;i++){rankedTop[i]=shortlist.get(i).ordinal;}
			}
			if(candidateOnly){
				for(int i=0;i<result.dpShortlist;i++){
					final Comparison c=shortlist.get(i);
					pairs.append(query.id).tab().append(c.family.rank).tab().append(c.family.id).tab()
						.append(c.ddlHits).tab().append(c.ddlPosition).tab().append(i+1).tab().append(result.cosineWindow)
						.tab().append("true").tab().append(Double.toString(c.combinedCosine)).nl();
				}
				counters.candidateRows+=result.dpShortlist; counters.queryCandidates++; result.terminal="candidates_selected";
				return false;
			}
			if(scalarGates){
				prepareQueryScalars(raw,counters);
				int retained=0;
				// Freeze the original top cap first. Filtering never backfills lower-ranked families.
				for(int i=0;i<result.dpShortlist;i++){
					final Comparison c=shortlist.get(i); final FamilyReference f=c.family;
					final int failures=scalarFailureMask(f,queryScalars,queryScalarMask);
					if(candidateObserver!=null){candidateObserver.scalarCandidate(c.ordinal,i+1,failures);}
					counters.scalarPairsChecked++;
					for(int s=0;s<5;s++){
						final int bit=1<<s;
						if((f.scalarMask&bit)==0){counters.scalarDisabled[s]++;}
						else if((queryScalarMask&bit)==0){counters.scalarQueryUndefined[s]++;}
						else{counters.scalarTested[s]++; if((failures&bit)!=0){counters.scalarFailed[s]++;}}
					}
					if(failures!=0){counters.scalarPairsRemoved++;}
					else{counters.scalarPairsKept++; shortlist.set(retained++,c);}
				}
				result.dpShortlist=retained;
				if(retained==0){
					result.clearAssignment(); result.terminal="scalar_reject";
					counters.queryScalarAllPruned++; counters.queryExhausted++; return false;
				}
			}
			if(familyGcPenalty!=null && familyGcPenalty.filterEnabled()){
				int retained=0;
				// Same frozen top cap; preserve the five scalar masks and do not backfill.
				for(int i=0;i<result.dpShortlist;i++){
					final Comparison c=shortlist.get(i);
					counters.gcPairsChecked++;
					if(familyGcPenalty.passes(c.ordinal,queryGcAvailable,queryGc)){shortlist.set(retained++,c);}
					else{counters.gcPairsRemoved++;}
				}
				result.dpShortlist=retained;
				if(retained==0){
					result.clearAssignment();result.terminal="gc_reject";
					counters.queryGcAllPruned++;counters.queryExhausted++;
					return false;
				}
			}
			return true;
		}
		private void scoreInverseRarityRanking(final int size){
			Comparison baselineBest=null;
			for(int i=0;i<size;i++){
				final Comparison c=shortlist.get(i);
				if(inverseRarityScratch!=null){
					inverseRarityScratch[i]=FamilyInverseRarityScore.rarityRatio(c.ddlHits,inverseRarityScores[c.ordinal]);
				}
				final double baseline=FamilyInverseRarityScore.scorePrepared(c.ddlHits, c.inverseLengthFactor,
					c.combinedCosine, c.inverseHydrophobicity);
				if(baselineBest==null || baseline>baselineBest.heuristicScore ||
						(baseline==baselineBest.heuristicScore && c.family.rank<baselineBest.family.rank)){
					baselineBest=c;
				}
				c.heuristicScore=baseline;
			}
			final double median=inverseRarityScratch==null ? 1 : FamilyInverseRarityScore.median(inverseRarityScratch, size);
			for(int i=0;i<size;i++){
				final Comparison c=shortlist.get(i);
				c.heuristicScore=inverseRarityRanker.scorePrepared(c.ddlHits, inverseRarityScores[c.ordinal],
					c.inverseLengthFactor, c.combinedCosine, c.inverseHydrophobicity, median);
				c.heuristicProtected=inverseRarityRanker.protectsBaselineRankOne() && c==baselineBest;
				c.heuristicMeasured=true;
			}
		}
		/** Both pool counting and shortlist insertion must consult the same immutable vector. */
		private int rawHitFloor(final int family,final int queryFloor){
			return hitFloors!=null ? hitFloors.floor(family) :
				familyAdmission!=null ? familyAdmission.floor(family) : queryFloor;
		}
		void scoreIndex(final byte[] raw){
			if(!indexPruned){
				sketchDdl(raw,geometry,seed,positionMode,partitions,overlap,blacklist,sketch);
				work.rarityScores=inverseRarityScores;
				index.score(sketch,work);
				return;
			}
			if(indexOrder==3){
				sketchDdl(raw,geometry,seed,positionMode,partitions,overlap,blacklist,sketch);
				index.scorePruned(sketch,work,null,sketch.length,sketch.length/2,indexDirect);
			}else if(indexOrder==0){
				sketchDdl(raw,geometry,seed,positionMode,partitions,overlap,blacklist,sketch);
				int keys=0;
				for(int bin=0; bin<sketch.length; bin++){if(sketch[bin]!=0){occupied[keys++]=bin;}}
				final int half=(keys+1)/2;
				for(int i=0; i<keys; i++){keyOrder[i]=occupied[i<half ? 2*i : 2*(i-half)+1];}
				index.scorePruned(sketch,work,keyOrder,keys,half,indexDirect);
			}else{
				DDLFamilyProfiles.sketchPositional(raw,geometry,seed,partitions,overlap,sketch,blacklist,keyOrigins,positionMode);
				final int[] order=keyOrigins.orderQuadrants(sketch,indexOrder==2);
				index.scorePruned(sketch,work,order,keyOrigins.size,keyOrigins.keysAt[1],indexDirect);
			}
		}
		private void fullAdmittedFeatures(ProteinSequence query,boolean scalarGates){
			final int size=shortlist.size();
			if(featureOrdinals==null || featureOrdinals.length<families.length){
				featureOrdinals=new int[families.length]; featureRanks=new int[families.length];
				featureHits=new int[families.length]; featureScalarMasks=new int[families.length];
				featureRarity=new double[families.length]; featureLength=new double[families.length];
				featureCosine=new double[families.length]; featureHydrophobicity=new double[families.length];
			}
			for(int i=0; i<size; i++){
				final Comparison c=shortlist.get(i);
				featureOrdinals[i]=c.ordinal; featureRanks[i]=c.family.rank; featureHits[i]=c.ddlHits;
				featureRarity[i]=inverseRarityScores[c.ordinal];
				featureLength[i]=c.inverseLengthFactor;
				featureCosine[i]=c.combinedCosine;
				featureHydrophobicity[i]=c.inverseHydrophobicity;
				featureScalarMasks[i]=scalarGates ? scalarFailureMask(c.family,queryScalars,queryScalarMask) : 0;
			}
			candidateObserver.fullAdmittedFeatures(query.id,query.length(),size,featureOrdinals,featureRanks,featureHits,
				featureRarity,featureLength,featureCosine,featureHydrophobicity,featureScalarMasks);
		}
		void search(final byte[] raw, final ProteinSequence query, final Result result,
				final Counters counters, final ByteBuilder pairs){
			if(!prepare(raw,query,result,counters,pairs,false)){return;}
			alignPrepared(query,result,counters,pairs);
		}
		void alignPrepared(final ProteinSequence query,final Result result,final Counters counters,final ByteBuilder pairs){
			alignPrepared(query,result,counters,pairs,null,margins);
		}
		void alignPrepared(final ProteinSequence query,final Result result,final Counters counters,final ByteBuilder pairs,
				final AlignmentProvider provider,final EvaluationMargins chosenMargins){
			require(provider==null || information==null,
				"An aggregate-only alignment provider cannot supply exact-match information");
			require(provider==null || residueCounts==null,"An aggregate-only alignment provider cannot supply residue counts");
			require(provider==null || lengthPolicy==null,
				"Length policy requires real traceback spans, not aggregate cached metrics");
			require(provider==null || !overhangDiagnostics,"Overhang diagnostics cannot compare aggregate cached providers");
			// Calibration visits the full SAME ranked DDL top+ties window. Deployment cap remains top.
			final int alignmentLimit=union ? shortlist.size() : result.dpShortlist;
			final int self=ProteinSearcher.selfScore(query);
			for(int i=0; i<alignmentLimit; i++){
				final Comparison c=shortlist.get(i);
				result.set(c);
				result.dpAttempt=i+1;
				result.dpAttempts++;
				if(i==0){counters.queryReachingDp++;}
				counters.candidateDp++;
				if(c.family.alignmentCalibration!=FamilyReference.AlignmentCalibration.ENABLED){
					result.terminal="alignment_unavailable";
					c.alignmentState=result.alignmentState=FamilyAlignmentMetric.State.NO_ALIGNMENT;
					counters.candidateAlignmentUnavailable++;result.alignmentUnavailable++;
					if(candidateObserver!=null){
						candidateObserver.unalignedCandidate(c.ordinal,i+1,result.terminal,result.alignmentState);
					}
					continue;
				}
				if(provider==null){
					alignment.measure(query,c.family.consensus,self,observation,residueCounts);
					// INPUT_TOO_SHORT returns before the aligner; failed/short-span DP still counts.
					if(observation.state==FamilyAlignmentMetric.State.INPUT_TOO_SHORT){result.alignmentInputSkips++;}
					else{result.alignmentCalls++;}
					if(overhangDiagnostics && observation.state!=FamilyAlignmentMetric.State.INPUT_TOO_SHORT){
						overhangCompared++;
						if(alignment.maximalTrimScore(query,c.family.consensus)>0){
							overhangTrimPositive++;
							if(observation.rawScore<=0){overhangTrimPositiveNonpositive++;}
						}
						if(observation.measured()){overhangMeasured++;}
						if(observation.state==FamilyAlignmentMetric.State.SPAN_TOO_SHORT){overhangSpanFailed++;}
					}
					c.alignmentState=result.alignmentState=observation.state;
					if(!observation.measured()){
						result.terminal=alignmentTerminal(observation.state);
						counters.candidateAlignmentFailed++;
						if(candidateObserver!=null){
							candidateObserver.unalignedCandidate(c.ordinal,i+1,result.terminal,result.alignmentState,
								observation.queryStart,observation.queryEnd,observation.referenceStart,observation.referenceEnd);
						}
						continue;
					}
					metrics.identityPct=observation.identityPct();
					metrics.blosum=observation.blosum();
					result.querySpan=observation.queryEnd-observation.queryStart+1;
					result.referenceSpan=observation.referenceEnd-observation.referenceStart+1;
					result.rawScore=observation.rawScore;
					result.denominatorSelf=observation.denominatorSelf;
					if(information!=null){
						require(observation.informationMeasured,"Pinned alignment did not measure requested exact-match information");
						result.informationRatio=information.ratio(observation.information);
						result.informationMeasured=true;
					}
				}else{
					provider.align(query,c,self,metrics);result.alignmentProviderCalls++;
					c.alignmentState=result.alignmentState=FamilyAlignmentMetric.State.MEASURED;
				}
				if(union){
					counters.unionStats.observe(c.ordinal,metrics); counters.candidateUnion++;
					pairs.append(query.id).tab().append(c.family.rank).tab().append(c.family.id).tab()
						.append(c.ddlHits).tab().append(c.ddlPosition).tab().append(i+1).tab().append(result.cosineWindow).tab()
						.append(i<top ? "true" : "false").tab().append(Double.toString(c.combinedCosine)).tab()
						.append(Double.toString(metrics.identityPct)).tab().append(Double.toString(metrics.blosum)).nl();
					result.terminal="union_measured"; continue;
				}
				result.identity=c.identity=(float)metrics.identityPct;
				result.blosum=c.blosum=(float)metrics.blosum;
				final double identityGate=chosenMargins==null ?
					Math.max(c.family.identityUpper,identityFloor) : chosenMargins.identity[c.ordinal];
				final double blosumGate=chosenMargins==null ?
					Math.max(c.family.blosumFinal,blosumFloor) : chosenMargins.blosum[c.ordinal];
				final float identityLower=effectiveIdentityLower(c.family,identityFloor,identityLowerReduction);
				if(candidateObserver!=null){
					candidateObserver.alignedCandidate(c.ordinal,i+1,result.identity,result.blosum,
						identityGate,blosumGate,identityLower,result.querySpan,result.referenceSpan);
				}
				final int lengthRoute=lengthPolicy==null ? FamilyLengthAcceptance.REJECT :
					FamilyLengthAcceptance.route(result.identity,result.blosum,identityLower,identityGate,blosumGate,
						result.querySpan,result.referenceSpan,lengthAmplitude);
				final boolean identityAccept=lengthPolicy==null ?
					result.identity>identityGate : lengthRoute==FamilyLengthAcceptance.IDENTITY;
				if(overhangDiagnostics && identityAccept){overhangIdentityEligible++;}
				// Signed calibration is retained above; biological evidence must be positive.
				// Cached providers retain the ratio with a strictly positive self denominator,
				// so its sign is exactly the raw-score sign even without a stored traceback.
				if(overhangOnly && (provider==null ? result.rawScore<=0 : metrics.blosum<=0)){
					if(overhangDiagnostics && identityAccept){overhangIdentityRemoved++;}
					result.terminal="blosum_reject";result.blosumRejects++;counters.candidateBlosumReject++;continue;
				}
				if(identityAccept){
					result.terminal="identity_accept";
					if(provider==null){result.alignmentAccepts++;}
					counters.candidateIdentityUpper++; counters.queryIdentityAccept++;
					counters.winner(2,c.ddlPosition,i+1); return;
				}
				if(result.deferredFamily<0 && result.identity>c.family.identityUpper){
					// First suppressed old identity winner; preserve it even if later candidates recover.
					assert(result.identity>=identityLower) :
						"Frozen identity upper includes member lower; a deferred accept must reach BLOSUM";
					result.deferredFamily=c.ordinal; result.deferredIdentity=result.identity; result.deferredBlosum=result.blosum;
					result.deferredBlosumPass=result.blosum>blosumGate;
				}
				if(lengthPolicy!=null){
					final double penalty=FamilyLengthAcceptance.penalty(result.querySpan,result.referenceSpan,lengthAmplitude);
					if(result.identity<identityLower+penalty){
						result.terminal="identity_reject"; result.identityRejects++; counters.candidateIdentityLow++; continue;
					}
				}else if(result.identity<identityLower){
					result.terminal="identity_reject"; result.identityRejects++; counters.candidateIdentityLow++; continue;
				}
				if(lengthPolicy==null ? result.blosum>blosumGate : lengthRoute==FamilyLengthAcceptance.BLOSUM){
					result.terminal="blosum_accept";
					if(provider==null){result.alignmentAccepts++;}
					counters.candidateBlosumUpper++; counters.queryBlosumUpperAccept++;
					counters.winner(3,c.ddlPosition,i+1); return;
				}
				result.terminal="blosum_reject"; result.blosumRejects++; counters.candidateBlosumReject++;
			}
			final FamilyAlignmentMetric.State finalState=result.alignmentState;
			result.clearAssignment();
			if(finalState!=FamilyAlignmentMetric.State.MEASURED){result.alignmentState=finalState;}
			counters.queryExhausted++;
			counters.queryExhaustedAfterDp++;
		}
		void measurePrepared(final ProteinSequence query,final Result result,final Counters counters,
				final int limit,final MetricSink sink){
			require(legacyAlignment,"NN feature export currently supports legacy full-query alignment metrics only");
			final int alignmentLimit=Math.min(Math.min(limit,result.dpShortlist),shortlist.size());
			final int self=ProteinSearcher.selfScore(query);
			for(int i=0; i<alignmentLimit; i++){
				final Comparison c=shortlist.get(i);
				result.set(c);
				result.dpAttempt=i+1;
				result.dpAttempts++;
				if(i==0){counters.queryReachingDp++;}
				counters.candidateDp++;
				FamilyCalibration.align(query,c.family.consensus,self,metrics);
				result.identity=c.identity=(float)metrics.identityPct;
				result.blosum=c.blosum=(float)metrics.blosum;
				sink.accept(i+1,c,metrics);
			}
		}
		interface MetricSink{
			void accept(int rank,Comparison candidate,FamilyCalibration.Metrics metrics);
		}
		/** Inspect the common DDL top+ties window without changing comparison scores or scoring counters. */
		void diagnose(byte[] raw,String id,Result result,Counters counters){
			counters.diagnosticQueries++;
			final Truth expected=truth.get(id); result.truth=expected;
			if(expected==null){result.boundary=0;}
			else if(expected.ordinal<0){result.boundary=1;}
			else{result.boundary=2;}
			int limit=Math.min(top,shortlist.size());
			if(limit>0){
				final int cutoff=shortlist.get(limit-1).ddlHits;
				while(limit<shortlist.size() && shortlist.get(limit).ddlHits==cutoff){limit++;}
			}
			result.windowSize=limit;
			int position=-1;
			if(expected!=null && expected.ordinal>=0){
				for(int i=0; i<shortlist.size(); i++){
					if(shortlist.get(i).ordinal==expected.ordinal){position=i; result.boundary=i<limit ? 4 : 3; break;}
				}
			}
			// Evaluate the entire common window, including DDL-autoaccepted queries and missing truth.
			queryCounts.clear(); queryCounts.add(raw); l1(queryCounts,querySpectrum);
			dimerRowSums(querySpectrum,queryMonomer);
			centerMonomers(queryMonomer,background,queryCentered);
			final float dimerNorm=norm(querySpectrum),monomerNorm=norm(queryMonomer),centeredNorm=norm(queryCentered);
			if(combined){
				prepareQueryScalars(raw,counters);
				scalarBounds.combineComputed(queryMonomer,true,queryScalars,queryScalarMask,queryCombined);
			}
			final double queryCombinedNorm=combined ? combinedNorm(queryCombined) : 0;
			for(int i=0; i<limit; i++){
				final Comparison c=shortlist.get(i);
				c.diagnosticDimer=cosine(querySpectrum,c.family.spectrum,dimerNorm,c.family.spectrumNorm);
				c.diagnosticMonomer=cosine(queryMonomer,monomerSpectra[c.ordinal],monomerNorm,monomerNorms[c.ordinal]);
				c.diagnosticCentered=cosine(queryCentered,centeredSpectra[c.ordinal],centeredNorm,centeredNorms[c.ordinal]);
				if(combined){
					c.diagnosticCombined=combinedCosine(queryCombined,combinedSpectra[c.ordinal],
						queryCombinedNorm,combinedNorms[c.ordinal]);
				}
				counters.diagnosticFamilies++;
			}
			if(position>=0 && position<limit){
				final Comparison target=shortlist.get(position);
				result.truthRanks[0]=position+1; result.truthRanks[1]=result.truthRanks[2]=result.truthRanks[3]=1;
				final float targetDimer=target.ddlHits+target.diagnosticDimer,targetMonomer=target.ddlHits+target.diagnosticMonomer;
				final float targetCentered=target.ddlHits+target.diagnosticCentered;
				if(combined){result.truthRanks[4]=1;}
				for(int i=0; i<limit; i++){
					final Comparison c=shortlist.get(i);
					if(precedes(c.ddlHits+c.diagnosticDimer,c.family.rank,targetDimer,target.family.rank)){result.truthRanks[1]++;}
					if(precedes(c.ddlHits+c.diagnosticMonomer,c.family.rank,targetMonomer,target.family.rank)){result.truthRanks[2]++;}
					if(precedes(c.ddlHits+c.diagnosticCentered,c.family.rank,targetCentered,target.family.rank)){
						result.truthRanks[3]++;
					}
					if(combined && precedes(c.ddlHits+c.diagnosticCombined,c.family.rank,
							target.ddlHits+target.diagnosticCombined,target.family.rank)){result.truthRanks[4]++;}
				}
			}
			for(int mode=0; mode<(combined ? 5 : 4); mode++){counters.truthStatus[mode][truthStatus(result,mode,top)]++;}
		}
		final DDLFamilySearch.Index index;
		final FamilyReference[] families;
		final ProteinKmerGeometry geometry;
		final int top, ddlAccept, k, partitions, overlap;
		int ddlReject;
		final DDLPentamerBlacklist blacklist;
		final float identityFloor,blosumFloor,identityLowerReduction;
		final EvaluationMargins margins;
		final FamilyLengthScore.Worker lengthRanker;
		final FamilyLengthScore.Mode lengthMode;
		final FamilyHeuristicScore heuristicRanker;
		final FamilyGcPenalty familyGcPenalty;
		private boolean queryGcInitialized,queryGcAvailable;
		private double queryGc;
		/** Reject stale/malformed values before any early-return scoring branch. */
		void setQueryGc(boolean available,double value){
			require(!queryGcInitialized,"Query GC must be cleared before worker reuse");
			require(Double.isFinite(value) && value>=0 && value<=1 && (available || value==0),
				"Query GC must be a finite measured fraction, or explicit missing with zero storage");
			queryGcAvailable=available;queryGc=value;queryGcInitialized=true;
		}
		void clearQueryGc(){queryGcInitialized=false;queryGcAvailable=false;queryGc=0;}
		final FamilyInverseRarityScore inverseRarityRanker;
		final boolean lengthScaledShortAdmission;
		final FamilyLengthAdmission familyAdmission;
		final FamilyHitFloorVector hitFloors;
		final FamilyLengthPolicy lengthPolicy;
		final float lengthAmplitude;
		final int[] rankedTop;int rankedTopSize;long lengthGeneration;
		final long seed;
		void prepareQueryScalars(byte[] raw,Counters counters){
			if(queryScalarsReady){return;}
			queryScalarMask=ConsensusScalars.computeScalarMask(raw,raw.length,queryScalars,0);
			queryScalarsReady=true; counters.scalarComputations++;
		}
		final boolean useCosine,gateCosine,centered,combined,union,scalarGates,candidateOnly,positional,legacyAlignment;
		final ProteinPositionMode positionMode;
		final double[] queryScalars=new double[5];
		int queryScalarMask;
		boolean queryScalarsReady;
		final ScalarBounds scalarBounds;
		final double[][] combinedSpectra;
		final double[] combinedNorms,queryCombined;
		final char[] sketch;
		final Comparison[] comparisons;
		final ArrayList<Comparison> shortlist=new ArrayList<Comparison>();
		final ArrayList<Comparison> dirty=new ArrayList<Comparison>();
		final DDLFamilySearch.Workspace work;
		final double[] inverseRarityScores, inverseRarityScratch;
		int[] featureOrdinals,featureRanks,featureHits,featureScalarMasks;
		double[] featureRarity,featureLength,featureCosine,featureHydrophobicity;
		RawHitObserver rawHitObserver;
		CandidateObserver candidateObserver;
		FamilyAlignmentCounts residueCounts;
		final boolean overhangOnly;
		boolean overhangDiagnostics;
		long overhangCompared,overhangTrimPositive,overhangTrimPositiveNonpositive;
		long overhangMeasured,overhangSpanFailed,overhangIdentityEligible,overhangIdentityRemoved;
		final boolean indexPruned,indexDirect;
		final int indexOrder;
		final int[] keyOrder,occupied;
		final ddl.DDLQueryOrigins keyOrigins;
		final ProteinSpectrum.Spectrum queryCounts;
		final float[] querySpectrum,queryRankingSpectrum,queryMonomer,queryCentered;
		final float[][] monomerSpectra,centeredSpectra;
		final float[] monomerNorms,centeredNorms,background;
		final HashMap<String,Truth> truth;
		final FamilyCalibration.Metrics metrics;
		final FamilyAlignmentMetric.Worker alignment;
		final FamilyAlignmentMetric.Observation observation;
		final FamilyMatchInformation.Binding information;
	}

	private static String alignmentTerminal(FamilyAlignmentMetric.State state){
		if(state==FamilyAlignmentMetric.State.INPUT_TOO_SHORT){return "input_too_short";}
		if(state==FamilyAlignmentMetric.State.NO_ALIGNMENT){return "no_alignment";}
		if(state==FamilyAlignmentMetric.State.SPAN_TOO_SHORT){return "span_too_short";}
		throw new IllegalArgumentException("Measured alignments have no failed terminal");
	}

	static final class Result{
		void clear(){clearAssignment(); dpAttempts=windowSize=0; boundary=0; truth=null; Arrays.fill(truthRanks,0);
			shortlistCandidates=alignmentCalls=alignmentAccepts=alignmentInputSkips=alignmentUnavailable=
				alignmentProviderCalls=0;
			deferredFamily=-1; deferredIdentity=deferredBlosum=0; deferredBlosumPass=false;
			terminal="no_hit";
			admissionFloor=admissionFloorMax=ddlActive=cosineWindow=dpShortlist=identityRejects=blosumRejects=0;
		}
		void clearAssignment(){
			informationRatio=0;informationMeasured=false;
			rawScore=denominatorSelf=0;
			querySpan=referenceSpan=0;
			family=-1; ddl=-1; cosine=identity=blosum=0; combinedCosine=0;
			cosineMeasured=combinedMeasured=false; ddlPosition=dpAttempt=0; rankedPosition=-1;
			alignmentState=FamilyAlignmentMetric.State.NO_ALIGNMENT;
		}
		void set(Comparison c){
			informationRatio=0;informationMeasured=false;
			rawScore=denominatorSelf=0;
			querySpan=referenceSpan=0;
			family=c.ordinal; ddl=c.ddlHits; cosine=c.cosine; identity=blosum=0;
			combinedCosine=c.combinedCosine;
			cosineMeasured=c.cosineMeasured; combinedMeasured=c.combinedMeasured;
			ddlPosition=c.ddlPosition; rankedPosition=c.rankedPosition; dpAttempt=0; alignmentState=c.alignmentState;
		}
		boolean assigned(){return family>=0;}
		double informationRatio;
		boolean informationMeasured;
		int rawScore,denominatorSelf;
		int querySpan,referenceSpan;
		int family=-1, ddl=-1;
		int deferredFamily=-1;
		float deferredIdentity,deferredBlosum;
		boolean deferredBlosumPass;
		float cosine,identity,blosum;
		double combinedCosine;
		boolean cosineMeasured,combinedMeasured;
		FamilyAlignmentMetric.State alignmentState=FamilyAlignmentMetric.State.NO_ALIGNMENT;
		int ddlPosition,dpAttempt,dpAttempts,windowSize,boundary;
		int rankedPosition=-1;
		int shortlistCandidates,alignmentCalls,alignmentAccepts,alignmentInputSkips,alignmentUnavailable,
			alignmentProviderCalls;
		String terminal="no_hit";
		int ddlActive,cosineWindow,dpShortlist,identityRejects,blosumRejects;
		int admissionFloor,admissionFloorMax;
		Truth truth;
		final int[] truthRanks=new int[RANK_MODES.length];
	}

	static final class Comparison{
		Comparison(int ordinal_,FamilyReference family_,float[] rankingSpectrum_,float rankingNorm_){
			ordinal=ordinal_; family=family_; rankingSpectrum=rankingSpectrum_; rankingNorm=rankingNorm_;
		}
		void reset(){
			ddlHits=-1; ddlPosition=0; rankedPosition=-1; cosine=identity=blosum=0; combinedCosine=heuristicScore=0;
			inverseLengthFactor=inverseHydrophobicity=0;
			cosineMeasured=combinedMeasured=heuristicMeasured=heuristicProtected=false;
			alignmentState=FamilyAlignmentMetric.State.NO_ALIGNMENT;
		}
		float rank(){return ddlHits+cosine;}
		final int ordinal;
		final FamilyReference family;
		final float[] rankingSpectrum;
		final float rankingNorm;
		int ddlHits=-1,ddlPosition=0;
		int rankedPosition=-1;
		float cosine,identity,blosum;
		float diagnosticDimer,diagnosticMonomer,diagnosticCentered;
		double combinedCosine,diagnosticCombined,lengthScore,heuristicScore,inverseLengthFactor,inverseHydrophobicity;
		boolean cosineMeasured,combinedMeasured,heuristicMeasured,heuristicProtected;
		FamilyAlignmentMetric.State alignmentState=FamilyAlignmentMetric.State.NO_ALIGNMENT;
	}

	private static FamilyReference[] legacyFamilies(Model[] models,Gate[] gates,
			ProteinSpectrum.Spectrum[] spectra,FamilyCensus.Entry[] metadata,int width){
		final FamilyReference[] families=new FamilyReference[models.length];
		for(int i=0; i<families.length; i++){
			final Model model=models[i];
			final Gate gate=gates[i];
			final ProteinSpectrum.Spectrum spectrum=spectra[i];
			final float[] frequencies=new float[DIMER_WIDTH];
			l1(spectrum,frequencies);
			final FamilyCensus.Entry entry=metadata[i];
			families[i]=new FamilyReference(model.id,model.rank,entry.members,entry.residues,
				entry.validWords,model.sequence,new char[width],frequencies,
				FamilyGateCodec.lower(gate.cosineLower),FamilyGateCodec.lower(gate.identityLower),
				FamilyGateCodec.lower(gate.blosumLower),
				(float)gate.cosineUpper,(float)gate.identityUpper,
				(float)gate.blosumFinal);
		}
		return families;
	}

	private static void l1(ProteinSpectrum.Spectrum counts,float[] frequencies){
		if(counts.total==0){
			Arrays.fill(frequencies,0);
		}else{
			for(int i=0; i<frequencies.length; i++){frequencies[i]=(float)(counts.counts[i]/(double)counts.total);}
		}
	}

	/** Raw row sums in ProteinSpectrum's base20 left/right order; not direct residue counts.
	 * Terminal residues and canonical runs of length1 contribute no outgoing dimer. */
	public static void dimerRowSums(float[] dimers,float[] monomers){
		assert(dimers!=null && dimers.length==DIMER_WIDTH && monomers!=null && monomers.length==MONOMER_WIDTH) :
			"ProteinSpectrum k2 uses400 row-major entries; ranking projection requires20 distinct output entries";
		for(int row=0; row<MONOMER_WIDTH; row++){
			float sum=0;
			for(int col=0; col<MONOMER_WIDTH; col++){sum+=dimers[row*MONOMER_WIDTH+col];}
			monomers[row]=sum;
		}
	}

	/** Frozen population bounds; not fitted to the loaded roster or query panel. */
	public static final class ScalarBounds{
		public ScalarBounds(double[] min_,double[] max_,double scale_){
			if(min_==null || max_==null || min_.length!=ConsensusScalars.SCALAR_COUNT ||
					max_.length!=ConsensusScalars.SCALAR_COUNT || (scale_!=1 && scale_!=0.1)){
				throw new IllegalArgumentException("Combined25 requires five bounds and scalarscale=1 or 0.1");
			}
			min=min_.clone(); max=max_.clone(); scale=scale_;
			for(int i=0; i<min.length; i++){
				if(!finite(min[i]) || !finite(max[i]) || !(min[i]<max[i]) || !finite(max[i]-min[i])){
					throw new IllegalArgumentException("Combined25 scalar bounds must have finite positive range: "+SCALAR_NAMES[i]);
				}
			}
		}
		/** Copy raw outgoing row sums unchanged; append full-double consensus/query scalars.
		 * Only queries have the literal-zero override. No per-sequence scratch allocation. */
		public void combine(float[] monomers,byte[] raw,boolean query,double[] out){
			assert(monomers!=null && monomers.length==MONOMER_WIDTH && raw!=null && out!=null && out.length==COMBINED_WIDTH) :
				"Combined25 needs20 raw rows plus5 scalars in a separate reusable output; no residue-frequency substitution";
			boolean nonzero=false;
			for(int i=0; i<MONOMER_WIDTH; i++){
				assert(Float.isFinite(monomers[i]) && monomers[i]>=0) :
					"Nonnegative finite dimer row sums are required for [0,1] cosine";
				out[i]=monomers[i]; nonzero|=monomers[i]!=0;
			}
			if(query && !nonzero){Arrays.fill(out,0); return;}
			final int mask=ConsensusScalars.computeScalarMask(raw,raw.length,out,MONOMER_WIDTH);
			if(mask==0){
				throw new IllegalArgumentException("Combined25 reference has no canonical residues for consensus scalars");
			}
			scaleScalars(out,mask);
		}
		/** Reuse the once-per-query raw scalars for ranking. Undefined coordinates rank as zero;
		 * the autofail gate separately skips them through the validity mask. */
		public void combineComputed(float[] monomers,boolean query,double[] scalars,int mask,double[] out){
			if(monomers==null || monomers.length!=MONOMER_WIDTH || scalars==null || scalars.length!=5 ||
					out==null || out.length!=COMBINED_WIDTH || (mask&~ConsensusScalars.ALL_DEFINED)!=0){
				throw new IllegalArgumentException(
					"Combined25 requires20 monomers,5 scalars,25 output coordinates and a five-bit mask");
			}
			boolean nonzero=false;
			for(int i=0;i<MONOMER_WIDTH;i++){
				assert(Float.isFinite(monomers[i]) && monomers[i]>=0) : "Nonnegative finite monomer rows required";
				out[i]=monomers[i]; nonzero|=monomers[i]!=0;
			}
			if(query && !nonzero){Arrays.fill(out,0); return;}
			if(mask==0){throw new IllegalArgumentException("Combined25 nonzero spectrum has no defined scalars");}
			for(int i=0;i<5;i++){out[MONOMER_WIDTH+i]=scalars[i];}
			scaleScalars(out,mask);
		}
		private void scaleScalars(double[] out,int mask){
			for(int i=0; i<ConsensusScalars.SCALAR_COUNT; i++){
				final double value=out[MONOMER_WIDTH+i];
				if(!finite(value)){throw new IllegalArgumentException("Nonfinite combined25 scalar: "+SCALAR_NAMES[i]);}
				out[MONOMER_WIDTH+i]=(mask&(1<<i))==0 ? 0 : scale*Math.max(0,Math.min(1,(value-min[i])/(max[i]-min[i])));
			}
		}
		private final double[] min,max;
		private final double scale;
	}

	/** All violated enabled/defined scalar intervals, not merely the first failure.
	 * Bounds are inclusive; a passing pair still needs the existing alignment gates. */
	public static int scalarFailureMask(FamilyReference family,double[] values,int queryMask){
		if(family==null || values==null || values.length!=5 || (queryMask&~ConsensusScalars.ALL_DEFINED)!=0){
			throw new IllegalArgumentException("Scalar gate needs a family, five finite query scalars and a five-bit mask");
		}
		int failures=0;
		final int enabled=family.scalarMask&queryMask;
		for(int i=0;i<5;i++){
			if(!Double.isFinite(values[i])){throw new IllegalArgumentException("Nonfinite query scalar: "+SCALAR_NAMES[i]);}
			if((enabled&(1<<i))!=0 && (values[i]<family.scalarBounds[2*i] ||
					values[i]>family.scalarBounds[2*i+1])){failures|=1<<i;}
		}
		return failures;
	}

	private static ScalarBounds scalarBounds(HashMap<String,String> options,RankingFeatures features) throws Exception{
		return scalarBounds(options,features,FamilyCalibration.VerificationMode.SHA80);
	}

	/** Explicit resident input policy; the scalar schema and arithmetic are unchanged. */
	private static ScalarBounds scalarBounds(HashMap<String,String> options,RankingFeatures features,
			FamilyCalibration.VerificationMode verificationMode) throws Exception{
		require(verificationMode!=null,"Scalar bounds need an explicit verification policy");
		if(features!=RankingFeatures.COMBINED25){
			if(containsAny(options,new String[]{"scalarbounds","scalarboundssha80","scalarscale"})){
				throw new IllegalArgumentException("Scalar bounds/scale require rankingfeatures=combined25");
			}
			return null;
		}
		final double scale=Double.parseDouble(required(options,"scalarscale"));
		final Path file=Paths.get(required(options,"scalarbounds")).toAbsolutePath().normalize();
		verifyScalarInput(file,required(options,"scalarboundssha80"),verificationMode);
		final double[] min=new double[ConsensusScalars.SCALAR_COUNT],max=new double[ConsensusScalars.SCALAR_COUNT];
		final ByteFile reader=ByteFile.makeByteFile(file.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		int rows=0;
		try{
			header(reader,"scalar\trows\tmin\tmax\tmin_family\tmax_family");
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				lp.set(line);
				if(rows>=min.length || lp.terms()!=6 || !lp.termEquals(SCALAR_NAMES[rows],0) || lp.parseLong(1)!=79565){
					throw new IllegalArgumentException("Combined25 bounds require five ordered scalar rows from79565 families");
				}
				truthToken(lp.parseString(4)); truthToken(lp.parseString(5));
				min[rows]=Double.parseDouble(lp.parseString(2)); max[rows]=Double.parseDouble(lp.parseString(3)); rows++;
			}
		}finally{if(reader.close()){throw new IllegalStateException("Scalar bounds read failed");}}
		if(rows!=min.length){throw new IllegalArgumentException("Combined25 bounds are incomplete");}
		verifyScalarInput(file,required(options,"scalarboundssha80"),verificationMode);
		return new ScalarBounds(min,max,scale);
	}
	/** Shared combined25 parser for family-side calibration; identical pinned bounds and scaling. */
	public static ScalarBounds combinedBounds(HashMap<String,String> options) throws Exception{
		return scalarBounds(options,RankingFeatures.COMBINED25);
	}
	/** Preparation supplies inherited identity explicitly while retaining the canonical scalar parser. */
	static ScalarBounds combinedBounds(HashMap<String,String> options,
			FamilyCalibration.VerificationMode verificationMode) throws Exception{
		return scalarBounds(options,RankingFeatures.COMBINED25,verificationMode);
	}
	private static void verifyScalarInput(Path file,String identity,
			FamilyCalibration.VerificationMode verificationMode) throws Exception{
		if(verificationMode==FamilyCalibration.VerificationMode.SHA80){pin(file,identity);return;}
		if(!"unmeasured".equals(identity)){DigestSuffix.requireSuffix(identity,"inherited scalar identity");}
		require(Files.isRegularFile(file) && !Files.isSymbolicLink(file),
			"Structural scalar input must be a regular file: "+file);
	}

	public static double combinedNorm(double[] values){
		assert(values!=null && values.length==COMBINED_WIDTH) :
			"Combined cosine requires25 coordinates, not legacy float spectra";
		double sum=0;
		for(double value : values){sum+=value*value;}
		return Math.sqrt(sum);
	}
	/** Shared deployed summed similarity; query/family ordinal breaks exact summed-score ties. */
	public static double combinedRankingScore(int ddl,double cosine){
		assert(ddl>=0 && ddl<=1024 && Double.isFinite(cosine) && cosine>=0 && cosine<=1) :
			"DDL registers and clamped combined cosine bound the ranking score";
		return ddl+cosine;
	}

	public static double combinedCosine(double[] query,double[] reference,double qNorm,double rNorm){
		assert(query!=null && reference!=null && query.length==COMBINED_WIDTH && reference.length==COMBINED_WIDTH &&
			finite(qNorm) && finite(rNorm) && qNorm>=0 && rNorm>=0) :
			"Combined cosine must use finite norms from the same25D vectors";
		if(qNorm==0 || rNorm==0){return 0;}
		double dot=0;
		for(int i=0; i<COMBINED_WIDTH; i++){dot+=query[i]*reference[i];}
		final double value=dot/(qNorm*rNorm);
		if(!Double.isFinite(value)){
			throw new IllegalArgumentException("Combined cosine calculation produced a nonfinite score");
		}
		return Math.max(0,Math.min(1,value));
	}

	/** Member-weighted mean of the existing float row sums, not renormalized residue frequencies.
	 * Visit the loaded roster in order; multiply/add each coordinate and the member total in double,
	 * divide once, then cast each background coordinate once to float. Empty/zero-dimer families
	 * contribute neither numerator nor denominator. validWords counts configured protein words, not dimers.
	 * No eligible members means no defined background; fail during SearchConfig construction, before writers open. */
	public static float[] monomerBackground(FamilyReference[] families,float[][] monomers){
		assert(families!=null && monomers!=null && families.length==monomers.length) :
			"Background rows must match the loaded family roster so member weights apply to the right spectra";
		final double[] sums=new double[MONOMER_WIDTH];
		double members=0;
		for(int i=0; i<families.length; i++){
			final FamilyReference family=families[i];
			if(family.members==0 || family.spectrumNorm==0){continue;}
			assert(monomers[i]!=null && monomers[i].length==MONOMER_WIDTH) :
				"Eligible background family needs20 outgoing-dimer row sums: "+family.id;
			final double weight=family.members;
			members+=weight;
			for(int a=0; a<MONOMER_WIDTH; a++){sums[a]+=weight*monomers[i][a];}
		}
		if(members==0){
			throw new IllegalArgumentException(
				"Centered monomer background requires at least one nonempty family with nonzero dimer spectrum");
		}
		final float[] background=new float[MONOMER_WIDTH];
		for(int a=0; a<MONOMER_WIDTH; a++){background[a]=(float)(sums[a]/members);}
		return background;
	}

	/** Shared reference/query transform. All-zero dimer rows stay literal zero, not -background.
	 * Input and background are immutable; the caller supplies a distinct, reusable output. */
	public static void centerMonomers(float[] monomers,float[] background,float[] centered){
		assert(monomers!=null && background!=null && centered!=null && monomers.length==MONOMER_WIDTH &&
			background.length==MONOMER_WIDTH && centered.length==MONOMER_WIDTH && centered!=monomers && centered!=background) :
			"Centering requires20-entry vectors and a distinct output to preserve raw ranks and the shared background";
		boolean nonzero=false;
		for(float value : monomers){nonzero|=value!=0;}
		if(!nonzero){Arrays.fill(centered,0); return;}
		for(int a=0; a<MONOMER_WIDTH; a++){centered[a]=monomers[a]-background[a];}
	}

	private static float norm(float[] values){
		float sum=0;
		for(float value : values){sum+=value*value;}
		return (float)Math.sqrt(sum);
	}

	private static float cosine(float[] query,float[] reference,float qNorm,float rNorm){
		if(qNorm==0 || rNorm==0){return 0;}
		final float value=simd.Vector.fma(query,reference)/(qNorm*rNorm);
		if(!Float.isFinite(value)){throw new IllegalArgumentException("Cosine calculation produced a nonfinite score");}
		return value>1 ? 1 : value< -1 ? -1 : value;
	}

	private static final Comparator<Comparison> DDL_ORDER=new Comparator<Comparison>(){
		@Override
		public int compare(Comparison a,Comparison b){
			final int ddl=Integer.compare(b.ddlHits,a.ddlHits);
			return ddl==0 ? Integer.compare(a.family.rank,b.family.rank) : ddl;
		}
	};

	private static final Comparator<Comparison> RANKING_ORDER=new Comparator<Comparison>(){
		@Override
		public int compare(Comparison a,Comparison b){
			final int score=Float.compare(b.rank(),a.rank());
			return score==0 ? Integer.compare(a.family.rank,b.family.rank) : score;
		}
	};

	private static final Comparator<Comparison> COMBINED_RANKING_ORDER=new Comparator<Comparison>(){
		@Override
		public int compare(Comparison a,Comparison b){
			assert(a.combinedMeasured && b.combinedMeasured) :
				"Combined25 ranking must score the whole top+ties window before ordering candidates";
			final int score=Double.compare(combinedRankingScore(b.ddlHits,b.combinedCosine),
				combinedRankingScore(a.ddlHits,a.combinedCosine));
			return score==0 ? Integer.compare(a.family.rank,b.family.rank) : score;
		}
	};
	private static final Comparator<Comparison> LENGTH_RANKING_ORDER=new Comparator<Comparison>(){
		@Override public int compare(Comparison a,Comparison b){
			assert(a.combinedMeasured && b.combinedMeasured &&
				Double.isFinite(a.lengthScore) && Double.isFinite(b.lengthScore)) :
				"Every eligible raw-hit candidate must have finite learned features before pre-truncation ordering";
			final int score=Double.compare(b.lengthScore,a.lengthScore);
			return score==0 ? Integer.compare(a.family.rank,b.family.rank) : score;
		}
	};
	private static final Comparator<Comparison> HEURISTIC_RANKING_ORDER=new Comparator<Comparison>(){
		@Override public int compare(Comparison a,Comparison b){
			assert(a.combinedMeasured && b.combinedMeasured && a.heuristicMeasured && b.heuristicMeasured &&
				Double.isFinite(a.heuristicScore) && Double.isFinite(b.heuristicScore)) :
				"Every DDL-eligible family must have a finite measured heuristic score before strict-top ordering";
			if(a.heuristicProtected!=b.heuristicProtected){return a.heuristicProtected ? -1 : 1;}
			final int score=Double.compare(b.heuristicScore,a.heuristicScore);
			return score==0 ? Integer.compare(a.family.rank,b.family.rank) : score;
		}
	};

	private static Gates loadGates(Path lower,Path upper,String[] names){
		final HashMap<String,Gate> gates=new HashMap<String,Gate>();
		final ByteFile low=ByteFile.makeByteFile(lower.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			header(low,FamilyCalibration.CUTOFF_HEADER);
			for(byte[] line=low.nextLine();line!=null;line=low.nextLine()){
				lp.set(line);
				if(lp.terms()!=10){throw new IllegalArgumentException("Wrong lower cutoff field count");}
				if(!lp.parseString(3).equals("CALIBRATED")){continue;}
				final String id=lp.parseString(1);
				if(gates.put(id,new Gate(gateDouble(lp,4),gateDouble(lp,6),gateDouble(lp,8)))!=null){
					throw new IllegalArgumentException("Duplicate lower gate for family: "+id);
				}
			}
		}finally{if(low.close()){throw new IllegalStateException("Lower cutoff read failed");}}
		final ByteFile up=ByteFile.makeByteFile(upper.toString(),false);
		try{
			header(up,"roster_rank\tfamily_id\tnegative_count\tcosine_upper\taligned_count"+
				"\tidentity_pct_upper\tblosum_a_sq_upper");
			for(byte[] line=up.nextLine();line!=null;line=up.nextLine()){
				lp.set(line);
				if(lp.terms()!=7){throw new IllegalArgumentException("Wrong upper cutoff field count");}
				final Gate gate=gates.get(lp.parseString(1));
				if(gate==null){continue;}
				if(gate.upperSet){throw new IllegalArgumentException("Duplicate upper gate for family: "+lp.parseString(1));}
				gate.setRawUpper(gateDouble(lp,3),gateDouble(lp,5),gateDouble(lp,6));
			}
		}finally{if(up.close()){throw new IllegalStateException("Upper cutoff read failed");}}
		final Gate[] linked=new Gate[names.length];
		for(int i=0; i<names.length; i++){
			final Gate gate=gates.remove(names[i].substring(5));
			if(gate==null || !gate.complete()){throw new IllegalArgumentException("Missing calibrated gates: "+names[i]);}
			linked[i]=gate;
		}
		if(!gates.isEmpty()){throw new IllegalArgumentException("Extra calibrated gates: "+gates.keySet());}
		return new Gates(linked);
	}

	private static Gates loadGates(Path[] lowers,Path[] uppers,String[] names,FamilyCensus census){
		final HashMap<String,Gate> gates=new HashMap<String,Gate>();
		final HashMap<String,Boolean> zeros=new HashMap<String,Boolean>();
		final LineParser1 lp=new LineParser1('\t');
		for(Path lower : lowers){
			final ByteFile low=ByteFile.makeByteFile(lower.toString(),false);
			try{
				header(low,FamilyCalibration.CUTOFF_HEADER);
				for(byte[] line=low.nextLine();line!=null;line=low.nextLine()){
					lp.set(line);
					if(lp.terms()!=10){throw new IllegalArgumentException("Wrong lower cutoff field count");}
					final int rank=lp.parseInt(0);
					final String id=lp.parseString(1);
					final int members=lp.parseInt(2);
					final FamilyCensus.Entry entry=census.get(rank,id);
					if(entry.members!=members){throw new IllegalArgumentException("Lower/census member mismatch: "+id);}
					final String status=lp.parseString(3);
					if(members<1){
						if(!status.equals("NO_MEMBERS") || !lp.parseString(4).equals("NA") || lp.parseInt(5)!=0 ||
								!lp.parseString(6).equals("NA") || lp.parseInt(7)!=0 ||
								!lp.parseString(8).equals("NA") || lp.parseInt(9)!=0){
							throw new IllegalArgumentException("Malformed zero-member lower row: "+id);
						}
						if(zeros.put(id,Boolean.TRUE)!=null){throw new IllegalArgumentException("Duplicate zero lower row: "+id);}
					}else{
						if(!status.equals("CALIBRATED")){throw new IllegalArgumentException("Uncalibrated live lower row: "+id);}
						final int cosineRetained=lp.parseInt(5), identityRetained=lp.parseInt(7), blosumRetained=lp.parseInt(9);
						final double cosine=gateDouble(lp,4), identity=gateDouble(lp,6), blosum=gateDouble(lp,8);
						if(!fraction(cosine) || !percent(identity) || !finite(blosum) ||
									!retained(members,cosineRetained) || !retained(members,identityRetained) ||
									!retained(members,blosumRetained)){
							throw new IllegalArgumentException("Malformed live lower row: "+id);
						}
						if(gates.put(id,new Gate(cosine,identity,blosum))!=null){
							throw new IllegalArgumentException("Duplicate lower gate for family: "+id);
						}
					}
				}
			}finally{if(low.close()){throw new IllegalStateException("Lower cutoff read failed");}}
		}
		final HashMap<String,Boolean> uppersSeen=new HashMap<String,Boolean>();
		for(Path upper : uppers){
			final ByteFile up=ByteFile.makeByteFile(upper.toString(),false);
			try{
				header(up,"roster_rank\tfamily_id\tnegative_count\tcosine_upper\taligned_count"+
				"\tidentity_pct_upper\tblosum_a_sq_upper");
				for(byte[] line=up.nextLine();line!=null;line=up.nextLine()){
					lp.set(line);
					if(lp.terms()!=7){throw new IllegalArgumentException("Wrong upper cutoff field count");}
					final int rank=lp.parseInt(0);
					final String id=lp.parseString(1);
					final FamilyCensus.Entry entry=census.get(rank,id);
					if(entry.members<1){throw new IllegalArgumentException("Upper row for zero-member family: "+id);}
					if(lp.parseInt(2)!=100000 || lp.parseInt(4)!=500){
						throw new IllegalArgumentException("Unexpected upper counts: "+id);
					}
					final double cosine=gateDouble(lp,3), identity=gateDouble(lp,5), blosum=gateDouble(lp,6);
					if(!fraction(cosine) || !percent(identity) || !finite(blosum)){
						throw new IllegalArgumentException("Malformed upper row: "+id);
					}
					final Gate gate=gates.get(id);
					if(gate==null){throw new IllegalArgumentException("Upper row lacks lower gate: "+id);}
					if(uppersSeen.put(id,Boolean.TRUE)!=null){
						throw new IllegalArgumentException("Duplicate upper gate for family: "+id);
					}
					gate.setRawUpper(cosine,identity,blosum);
				}
			}finally{if(up.close()){throw new IllegalStateException("Upper cutoff read failed");}}
		}
		final Gate[] linked=new Gate[names.length];
		for(int i=0; i<names.length; i++){
			final String id=names[i].substring(5);
			final Gate gate=gates.remove(id);
			if(gate==null || !gate.complete()){throw new IllegalArgumentException("Missing collection gates: "+id);}
			linked[i]=gate;
		}
		if(!gates.isEmpty() || zeros.size()!=census.zeros){
			throw new IllegalArgumentException("Collection gates do not match census");
		}
		return new Gates(linked);
	}

	private static Gates loadGatesA48(Path[] lowers,Path[] uppers,String[] names,FamilyCensus census){
		final HashMap<String,Gate> gates=new HashMap<String,Gate>();
		final HashMap<String,Boolean> zeros=new HashMap<String,Boolean>();
		final LineParser1 lp=new LineParser1('\t');
		for(Path lower : lowers){
			final ByteFile low=ByteFile.makeByteFile(lower.toString(),false);
			try{
				header(low,FamilyGateCodec.LOWER_HEADER);
				for(byte[] line=low.nextLine();line!=null;line=low.nextLine()){
					lp.set(line);
					if(lp.terms()!=10){throw new IllegalArgumentException("Wrong A48 lower cutoff field count");}
					final int rank=nonnegativeInt(lp,0);
					final String id=lp.parseString(1);
					final int members=nonnegativeInt(lp,2);
					final FamilyCensus.Entry entry=census.get(rank,id);
					if(entry.members!=members){throw new IllegalArgumentException("A48 lower/census member mismatch: "+id);}
					final String status=lp.parseString(3);
					if(members<1){
						if(!status.equals("NO_MEMBERS") || !lp.parseString(4).equals("NA") || nonnegativeInt(lp,5)!=0 ||
								!lp.parseString(6).equals("NA") || nonnegativeInt(lp,7)!=0 ||
								!lp.parseString(8).equals("NA") || nonnegativeInt(lp,9)!=0){
							throw new IllegalArgumentException("Malformed zero-member A48 lower row: "+id);
						}
						if(zeros.put(id,Boolean.TRUE)!=null){throw new IllegalArgumentException("Duplicate zero A48 lower row: "+id);}
					}else{
						if(!status.equals("CALIBRATED")){throw new IllegalArgumentException("Uncalibrated live A48 lower row: "+id);}
						final int cosineRetained=nonnegativeInt(lp,5);
						final int identityRetained=nonnegativeInt(lp,7);
						final int blosumRetained=nonnegativeInt(lp,9);
						final double cosine=FamilyGateCodec.readLower(lp,4);
						final double identity=FamilyGateCodec.readLower(lp,6);
						final double blosum=FamilyGateCodec.readLower(lp,8);
						if(!fraction(cosine) || !percent(identity) || !finite(blosum) ||
									!retained(members,cosineRetained) || !retained(members,identityRetained) ||
									!retained(members,blosumRetained)){
							throw new IllegalArgumentException("Malformed live A48 lower row: "+id);
						}
						if(gates.put(id,new Gate(cosine,identity,blosum))!=null){
							throw new IllegalArgumentException("Duplicate A48 lower gate for family: "+id);
						}
					}
				}
			}finally{if(low.close()){throw new IllegalStateException("A48 lower cutoff read failed");}}
		}
		final HashMap<String,Boolean> uppersSeen=new HashMap<String,Boolean>();
		for(Path upper : uppers){
			final ByteFile up=ByteFile.makeByteFile(upper.toString(),false);
			try{
				header(up,FamilyGateCodec.UPPER_HEADER);
				for(byte[] line=up.nextLine();line!=null;line=up.nextLine()){
					lp.set(line);
					if(lp.terms()!=7){throw new IllegalArgumentException("Wrong A48 upper cutoff field count");}
					final int rank=nonnegativeInt(lp,0);
					final String id=lp.parseString(1);
					final FamilyCensus.Entry entry=census.get(rank,id);
					if(entry.members<1){throw new IllegalArgumentException("A48 upper row for zero-member family: "+id);}
					if(nonnegativeInt(lp,2)!=100000 || nonnegativeInt(lp,4)!=500){
						throw new IllegalArgumentException("Unexpected A48 upper counts: "+id);
					}
					final double cosine=FamilyGateCodec.readCosineAccept(lp,3);
					final double identity=FamilyGateCodec.readUpper(lp,5);
					final double blosum=FamilyGateCodec.readFinalBlosum(lp,6);
					if(!finite(cosine) || cosine<0 ||
							!finite(identity) || identity<0 || identity>Math.nextUp(110f) || !finite(blosum)){
						throw new IllegalArgumentException("Malformed A48 upper row: "+id);
					}
					final Gate gate=gates.get(id);
					if(gate==null){throw new IllegalArgumentException("A48 upper row lacks lower gate: "+id);}
					if(uppersSeen.put(id,Boolean.TRUE)!=null){
						throw new IllegalArgumentException("Duplicate A48 upper gate for family: "+id);
					}
					gate.setFinalUpper(cosine,identity,blosum);
				}
			}finally{if(up.close()){throw new IllegalStateException("A48 upper cutoff read failed");}}
		}
		final Gate[] linked=new Gate[names.length];
		for(int i=0; i<names.length; i++){
			final String id=names[i].substring(5);
			final Gate gate=gates.remove(id);
			if(gate==null || !gate.complete()){throw new IllegalArgumentException("Missing A48 collection gates: "+id);}
			linked[i]=gate;
		}
		if(!gates.isEmpty() || zeros.size()!=census.zeros){
			throw new IllegalArgumentException("A48 collection gates do not match census");
		}
		return new Gates(linked);
	}

	private static double gateDouble(LineParser1 lp,int term){
		return Double.parseDouble(lp.parseString(term));
	}

	private static boolean fraction(final double value){return value>=0 && value<=1 && value<Double.POSITIVE_INFINITY;}
	private static boolean percent(final double value){return value>=0 && value<=100 && value<Double.POSITIVE_INFINITY;}
	private static boolean finite(final double value){
		return value<Double.POSITIVE_INFINITY && value>Double.NEGATIVE_INFINITY;
	}
	private static boolean retained(final int members, final int retained){
		return retained>=members-members/100 && retained<=members;
	}

	private static DDLFamilySearch.Index loadIndex(ReferenceFiles[] references,int k,long seed,int width,
			Path work) throws Exception{
		final Path list=Files.createTempFile(work,"family_tail_cascade_profiles",".tsv");
		try{
			final StringBuilder sb=new StringBuilder();
			for(ReferenceFiles r : references){
				sb.append(r.profiles).append('\t').append(r.profilesHash).append('\t')
					.append(r.summaryHash).append('\n');
			}
			write(list,sb.toString());
			return DDLFamilySearch.loadForSearch(list,width,"nonempty",k,seed,true,false,false);
		}finally{Files.deleteIfExists(list);}
	}

	private static DDLFamilySearch.Index loadIndex(FamilyReference[] families,int threads){
		final String[] names=new String[families.length];
		final char[][] rows=new char[families.length][];
		for(int i=0; i<families.length; i++){
			final FamilyReference family=families[i];
			names[i]="tail:"+family.id;
			rows[i]=family.ddl;
		}
		return DDLFamilySearch.Index.parallel(names,rows,null,threads);
	}

	/** Replace only DDL rows, preserving the pinned package's family order and immutable payloads. */
	static DDLFamilySearch.Index loadReplacementIndex(FamilyReference[] families,Path list,int threads,
			boolean positional,int partitions,int overlap,DDLPentamerBlacklist blacklist) throws Exception{
		return loadReplacementIndex(families,list,ProteinKmerGeometry.exact(5),threads,
			positional,partitions,overlap,blacklist);
	}

	/** Replace only DDL rows, preserving the pinned package's family order and immutable payloads. */
	static DDLFamilySearch.Index loadReplacementIndex(FamilyReference[] families,Path list,
			ProteinKmerGeometry geometry,int threads,boolean positional,int partitions,int overlap,
			DDLPentamerBlacklist blacklist) throws Exception{
		return loadReplacementIndex(families,list,geometry,threads,
			positional ? ProteinPositionMode.FIXED : ProteinPositionMode.PLAIN,partitions,overlap,blacklist);
	}

	static DDLFamilySearch.Index loadReplacementIndex(FamilyReference[] families,Path list,
			ProteinKmerGeometry geometry,int threads,ProteinPositionMode positionMode,int partitions,
			int overlap,DDLPentamerBlacklist blacklist) throws Exception{
		return DDLFamilySearch.Index.parallel(replacementNames(families),
			loadReplacementRows(families,list,geometry,positionMode,partitions,overlap,blacklist),null,threads);
	}

	/** Default to halfway rejection for ordinary four-partition positional search.
	 * Calibration and rank diagnostics retain exhaustive indexing. Explicit direct
	 * completion also suppresses the default, while an explicit early flag wins.
	 */
	private static boolean indexEarlyReject(HashMap<String,String> options){
		assert(options!=null) : "Index policy must resolve from the same validated options used by the loader";
		final boolean production=options.containsKey("familyreferences") && options.containsKey("ddlprofiles") &&
			(!options.containsKey("mode") || options.get("mode").equals("search")) && !options.containsKey("rankout") &&
			positionMode(options,true).positional() &&
			(!options.containsKey("partitions") || integer(options,"partitions")==4) &&
			!optionalBoolean(options,"indexdirect",false);
		return optionalBoolean(options,"indexearlyreject",production);
	}

	/** Validate resolved search policy before reference loading. */
	private static boolean indexPruned(HashMap<String,String> options){
		final boolean early=indexEarlyReject(options),direct=optionalBoolean(options,"indexdirect",false);
		return indexPruned(early,direct,options);
	}

	/** Validate an already resolved search policy before reference loading. */
	private static boolean indexPruned(boolean early,boolean direct,HashMap<String,String> options){
		if(early && direct){
			throw new IllegalArgumentException("indexearlyreject and indexdirect are mutually exclusive experimental modes");
		}
		DDLFamilySearch.Index.parseScanOrder(options.containsKey("indexorder") ? options.get("indexorder") : "plain");
		return early || direct;
	}

	/** The default early-rejection index retains no supplementary packed family rows. */
	static DDLFamilySearch.Index loadReplacementIndex(FamilyReference[] families,Path list,
			int threads,boolean positional,int partitions,int overlap,DDLPentamerBlacklist blacklist,
			HashMap<String,String> options) throws Exception{
		return loadReplacementIndex(families,list,ProteinKmerGeometry.exact(5),threads,
			positional,partitions,overlap,blacklist,options);
	}

	static DDLFamilySearch.Index loadReplacementIndex(FamilyReference[] families,Path list,
			ProteinKmerGeometry geometry,int threads,boolean positional,int partitions,int overlap,
			DDLPentamerBlacklist blacklist,HashMap<String,String> options) throws Exception{
		return loadReplacementIndex(families,list,geometry,threads,
			positional ? ProteinPositionMode.FIXED : ProteinPositionMode.PLAIN,partitions,overlap,blacklist,options);
	}

	static DDLFamilySearch.Index loadReplacementIndex(FamilyReference[] families,Path list,
			ProteinKmerGeometry geometry,int threads,ProteinPositionMode positionMode,
			int partitions,int overlap,DDLPentamerBlacklist blacklist,
			HashMap<String,String> options) throws Exception{
		if(indexPruned(options) && (!positionMode.positional() || partitions!=4)){
			throw new IllegalArgumentException("Experimental index policy requires four positional partitions");
		}
		final int width=options.containsKey("width") ? integer(options,"width") : 1024;
		return replacementIndex(families,
			loadReplacementRows(families,list,geometry,positionMode,partitions,overlap,blacklist,width),
			threads,positionMode,partitions,options);
	}

	/** File-backed and resident modal rows use the same exhaustive/pruned index policy. */
	private static DDLFamilySearch.Index replacementIndex(FamilyReference[] families, char[][] rows, int threads,
			ProteinPositionMode positionMode, int partitions, HashMap<String,String> options){
		assert(families!=null && rows!=null && rows.length==families.length) :
			"One modal row is required per validated live family";
		if(!indexPruned(options)){return DDLFamilySearch.Index.parallel(replacementNames(families), rows, null, threads);}
		if(!positionMode.positional() || partitions!=4){
			throw new IllegalArgumentException("Experimental index policy requires four positional partitions");
		}
		return DDLFamilySearch.Index.parallelPruned(replacementNames(families), rows, threads,
			indexEarlyReject(options),optionalBoolean(options,"indexdirect",false),
			options.containsKey("indexorder") ? options.get("indexorder") : "plain");
	}

	/** Validate the whole live package roster before callers shard or rebuild a local index. */
	static char[][] loadReplacementRows(FamilyReference[] families,Path list,boolean positional,
			int partitions,int overlap,DDLPentamerBlacklist blacklist) throws Exception{
		return loadReplacementRows(families,list,ProteinKmerGeometry.exact(5),
			positional,partitions,overlap,blacklist);
	}

	/** Validate the whole live package roster before callers shard or rebuild a local index. */
	static char[][] loadReplacementRows(FamilyReference[] families,Path list,ProteinKmerGeometry geometry,
			boolean positional,int partitions,int overlap,DDLPentamerBlacklist blacklist) throws Exception{
		return loadReplacementRows(families,list,geometry,
			positional ? ProteinPositionMode.FIXED : ProteinPositionMode.PLAIN,partitions,overlap,blacklist);
	}

	/** Require the exact stored coordinate system, even when both modes use positional tags. */
	static char[][] loadReplacementRows(FamilyReference[] families,Path list,ProteinKmerGeometry geometry,
			ProteinPositionMode positionMode,int partitions,int overlap,
			DDLPentamerBlacklist blacklist) throws Exception{
		return loadReplacementRows(families,list,geometry,positionMode,partitions,overlap,blacklist,1024);
	}

	/** Explicit width is supplied by the validated classifier; historical helper callers retain1024. */
	private static char[][] loadReplacementRows(FamilyReference[] families,Path list,ProteinKmerGeometry geometry,
			ProteinPositionMode positionMode,int partitions,int overlap,DDLPentamerBlacklist blacklist,
			int width) throws Exception{
		return loadReplacementRows(families,list,geometry,positionMode,partitions,overlap,blacklist,width,null);
	}

	static char[][] loadReplacementRows(FamilyReference[] families,Path list,ProteinKmerGeometry geometry,
			ProteinPositionMode positionMode,int partitions,int overlap,DDLPentamerBlacklist blacklist,
			DDLFamilySearch.StructuralInputs structuralInputs) throws Exception{
		return loadReplacementRows(families,list,geometry,positionMode,partitions,overlap,blacklist,1024,
			structuralInputs);
	}

	private static char[][] loadReplacementRows(FamilyReference[] families,Path list,ProteinKmerGeometry geometry,
			ProteinPositionMode positionMode,int partitions,int overlap,DDLPentamerBlacklist blacklist,int width,
			DDLFamilySearch.StructuralInputs structuralInputs) throws Exception{
		assert(families!=null && list!=null) : "Replacement needs the accepted live roster and a pinned profile list";
		require(width==512 || width==1024,"Replacement profiles require validated width512 or1024");
		require(positionMode!=null,"Replacement profile identity requires an explicit coordinate mode");
		for(FamilyReference family:families){
			require(family!=null && family.members>0 && family.ddl.length==width,
				"Every live reference row must match the configured replacement width before index construction");
		}
		final DDLFamilySearch.Index source=structuralInputs==null ?
			DDLFamilySearch.loadForAttribution(list,width,geometry,0,positionMode.sketch,
				positionMode.positional() ? partitions : 1,positionMode.positional() ? overlap : 0,blacklist) :
			DDLFamilySearch.loadForAttributionStructural(list,width,geometry,0,positionMode.sketch,
				positionMode.positional() ? partitions : 1,positionMode.positional() ? overlap : 0,blacklist,structuralInputs);
		if(source.families()!=families.length){
			throw new IllegalArgumentException("Positional and packaged live-family counts differ");
		}
		final HashMap<String,Integer> ordinals=new HashMap<String,Integer>(families.length*2);
		for(int i=0;i<source.names.length;i++){
			if(ordinals.put(source.names[i],i)!=null){
				throw new IllegalArgumentException("Duplicate positional family: "+source.names[i]);
			}
		}
		final char[][] rows=new char[families.length][];
		for(int i=0;i<families.length;i++){
			final String name="tail:"+families[i].id;
			final Integer ordinal=ordinals.remove(name);
			if(ordinal==null){throw new IllegalArgumentException("Missing positional profile for packaged family: "+name);}
			rows[i]=source.profileCopy(ordinal.intValue());
		}
		if(!ordinals.isEmpty()){throw new IllegalArgumentException("Extra positional families outside the accepted package");}
		return rows;
	}

	private static String[] replacementNames(FamilyReference[] families){
		final String[] names=new String[families.length];
		for(int i=0;i<families.length;i++){names[i]="tail:"+families[i].id;}
		return names;
	}

	static long sketchDdl(byte[] raw,int k,long seed,boolean positional,int partitions,int overlap,
			DDLPentamerBlacklist blacklist,char[] sketch){
		return sketchDdl(raw,ProteinKmerGeometry.exact(k),seed,positional,partitions,overlap,blacklist,sketch);
	}

	static long sketchDdl(byte[] raw,ProteinKmerGeometry geometry,long seed,boolean positional,int partitions,int overlap,
			DDLPentamerBlacklist blacklist,char[] sketch){
		return sketchDdl(raw,geometry,seed,positional ? ProteinPositionMode.FIXED : ProteinPositionMode.PLAIN,
			partitions,overlap,blacklist,sketch);
	}

	static long sketchDdl(byte[] raw,ProteinKmerGeometry geometry,long seed,
			ProteinPositionMode positionMode,int partitions,int overlap,DDLPentamerBlacklist blacklist,
			char[] sketch){
		require(positionMode!=null,"Query sketch mode must match the loaded family coordinate system");
		return positionMode.positional() ?
			DDLFamilyProfiles.sketchPositional(raw,geometry,seed,partitions,overlap,positionMode,sketch,blacklist) :
			DDLFamilyProfiles.sketch(raw,geometry,seed,sketch,blacklist);
	}

	private static LinkedSpectra loadSpectra(ReferenceFiles[] references,String[] names) throws Exception{
		final HashMap<String,ProteinSpectrum.Spectrum> loaded=new HashMap<String,ProteinSpectrum.Spectrum>();
		for(ReferenceFiles r : references){
			final HashMap<String,ProteinSpectrum.Spectrum> map=ProteinSpectrum.load(r.spectra.toString());
			for(String key : map.keySet()){
				if(loaded.put(key,map.get(key))!=null){throw new IllegalArgumentException("Duplicate dimer spectrum: "+key);}
			}
		}
		final ProteinSpectrum.Spectrum[] spectra=new ProteinSpectrum.Spectrum[names.length];
		for(int i=0; i<names.length; i++){
			spectra[i]=loaded.remove(names[i]);
			if(spectra[i]==null || spectra[i].k!=2){
				throw new IllegalArgumentException("Missing or invalid dimer spectrum: "+names[i]);
			}
		}
		if(!loaded.isEmpty()){
			throw new IllegalArgumentException("Unmatched spectra outside DDL profile roster: "+loaded.keySet());
		}
		return new LinkedSpectra(spectra);
	}

	private static FamilyCensus.Entry[] metadata(FamilyCensus census,String[] names){
		final FamilyCensus.Entry[] metadata=new FamilyCensus.Entry[names.length];
		for(int i=0; i<names.length; i++){metadata[i]=census.get(names[i].substring(5));}
		return metadata;
	}

	private static FamilyCensus.Entry[] metadata(Path members,String[] names){
		final HashMap<String,FamilyCensus.Entry> loaded=new HashMap<String,FamilyCensus.Entry>(names.length*2);
		final HashMap<String,Boolean> seen=new HashMap<String,Boolean>(names.length*2);
		final LineParser1 lp=new LineParser1('\t');
		final ByteFile input=ByteFile.makeByteFile(members.toString(),false);
		try{
			header(input,"roster_rank\tfamily_id\tmembers\tresidues\tvalid_words");
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line);
				if(lp.terms()!=5){throw new IllegalArgumentException("Wrong members field count");}
				final String id=lp.parseString(1);
				final int count=lp.parseInt(2);
				final long residues=Long.parseLong(lp.parseString(3));
				final long words=Long.parseLong(lp.parseString(4));
				if(count<0 || residues<0 || words<0 || (count==0)!=(residues==0 && words==0) ||
						seen.put(id,Boolean.TRUE)!=null){
					throw new IllegalArgumentException("Malformed members row: "+id);
				}
				if(count==0){continue;}
				if(loaded.put(id,new FamilyCensus.Entry(id,count,residues,words))!=null){
					throw new IllegalArgumentException("Duplicate live members row: "+id);
				}
			}
		}finally{if(input.close()){throw new IllegalStateException("Members read failed");}}
		final FamilyCensus.Entry[] metadata=new FamilyCensus.Entry[names.length];
		for(int i=0; i<names.length; i++){
			final String id=names[i].substring(5);
			final FamilyCensus.Entry entry=loaded.remove(id);
			if(entry==null){throw new IllegalArgumentException("Missing members row: "+id);}
			metadata[i]=entry;
		}
		if(!loaded.isEmpty()){
			throw new IllegalArgumentException("Unmatched members outside DDL profile roster: "+loaded.keySet());
		}
		return metadata;
	}

	private static void validateProfileCensus(ReferenceFiles[] references,FamilyCensus census,int width){
		final HashMap<String,Boolean> seen=new HashMap<String,Boolean>(census.nonempty*2);
		final LineParser1 lp=new LineParser1('\t');
		for(ReferenceFiles r : references){
			final ByteFile input=ByteFile.makeByteFile(r.profiles.toString(),false);
			try{
				header(input,DDLFamilyProfiles.HEADER16);
				for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
					lp.set(line);
					if(lp.terms()!=10){throw new IllegalArgumentException("Wrong profile census field count");}
					final int buckets=lp.parseInt(2), members=lp.parseInt(3);
					final long residues=Long.parseLong(lp.parseString(4)), words=Long.parseLong(lp.parseString(5));
					final String namespace=lp.parseString(0), id=lp.parseString(1), mode=lp.parseString(6);
					if(!namespace.equals("tail")){
						throw new IllegalArgumentException("Unexpected profile row for collection: "+namespace+":"+id);
					}
					if(buckets!=width || !mode.equals("nonempty")){continue;}
					final FamilyCensus.Entry entry=census.get(id);
					if(entry.members<1 || entry.members!=members || entry.residues!=residues || entry.validWords!=words){
						throw new IllegalArgumentException("Profile/census mismatch: "+id);
					}
					if(seen.put(id,Boolean.TRUE)!=null){throw new IllegalArgumentException("Duplicate profile census row: "+id);}
				}
			}finally{if(input.close()){throw new IllegalStateException("Profile census read failed");}}
		}
		if(seen.size()!=census.nonempty){
			throw new IllegalArgumentException("Profile census roster does not cover every live family");
		}
	}

	private static Models loadModels(Path file,DDLFamilySearch.Index index,int maxReference) throws Exception{
		final Model[] models=new Model[index.families()];
		final ArrayList<String> canonicalIds=new ArrayList<String>();
		final HashMap<String,Integer> storage=new HashMap<String,Integer>(index.families()*2);
		for(int i=0; i<index.names.length; i++){
			if(!index.names[i].startsWith("tail:") || storage.put(index.names[i].substring(5),i)!=null){
				throw new IllegalArgumentException("Duplicate or non-tail DDL family: "+index.names[i]);
			}
		}
		final ByteFile input=ByteFile.makeByteFile1(file.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		int rank=0,found=0;
		try{
			header(input,HbmTailLibrary.MODELS_HEADER);
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line);
				if(lp.terms()!=9 || lp.parseInt(0)!=rank){
					throw new IllegalArgumentException("Fixed model roster is not canonical");
				}
				final String id=lp.parseString(1);
				final Integer indexSlot=storage.remove(id);
				if(indexSlot!=null){
					final Path consensus=Paths.get(lp.parseString(2)).resolve("consensus.faa");
					final String hash=lp.parseString(6);
					pin(consensus,hash);
					models[indexSlot.intValue()]=new Model(rank,id,consensus(consensus,maxReference,id));
					pin(consensus,hash);
					canonicalIds.add("tail:"+id); found++;
				}
				rank++;
			}
		}finally{if(input.close()){throw new IllegalStateException("Fixed models read failed");}}
		if(!storage.isEmpty() || found!=models.length){
			throw new IllegalArgumentException("Fixed models missing profile families: "+storage.keySet());
		}
		return new Models(models,canonicalIds.toArray(new String[canonicalIds.size()]));
	}

	private static Models loadModels(Path file,Path consensus,DDLFamilySearch.Index index,
			FamilyCensus census,int maxReference) throws Exception{
		final Model[] models=new Model[index.families()];
		final ArrayList<String> canonicalIds=new ArrayList<String>();
		final HashMap<String,Integer> storage=new HashMap<String,Integer>(index.families()*2);
		for(int i=0; i<index.names.length; i++){
			final String name=index.names[i];
			if(!name.startsWith("tail:") || storage.put(name.substring(5),i)!=null){
				throw new IllegalArgumentException("Duplicate or non-tail DDL family: "+index.names[i]);
			}
			if(census.members(name.substring(5))<1){
				throw new IllegalArgumentException("DDL profile does not identify a live census family: "+name);
			}
		}
		final HashMap<String,ProteinSequence> sequences=loadCombinedConsensus(consensus,maxReference);
		final HashMap<String,Boolean> seen=new HashMap<String,Boolean>(census.size()*2);
		final ByteFile input=ByteFile.makeByteFile1(file.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		int rank=0,found=0;
		try{
			header(input,HbmTailLibrary.MODELS_HEADER);
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line);
				if(lp.terms()!=9 || lp.parseInt(0)!=rank){
					throw new IllegalArgumentException("Fixed model roster is not canonical");
				}
				final String id=lp.parseString(1);
				if(rank>=census.size() || !id.equals(census.id(rank)) || seen.put(id,Boolean.TRUE)!=null){
					throw new IllegalArgumentException("Fixed model roster does not match current census: "+id);
				}
				final ProteinSequence sequence=sequences.remove(id);
				if(sequence==null){throw new IllegalArgumentException("Missing combined consensus: "+id);}
				final Integer indexSlot=storage.remove(id);
				if(indexSlot!=null){
					models[indexSlot.intValue()]=new Model(rank,id,sequence);
					canonicalIds.add("tail:"+id); found++;
				}
				rank++;
			}
		}finally{if(input.close()){throw new IllegalStateException("Fixed models read failed");}}
		if(rank!=census.size() || !storage.isEmpty() || found!=models.length || !sequences.isEmpty()){
			throw new IllegalArgumentException("Fixed models/consensuses missing collection profile families");
		}
		return new Models(models,canonicalIds.toArray(new String[canonicalIds.size()]));
	}

	private static HashMap<String,ProteinSequence> loadCombinedConsensus(Path path,int cap){
		//TODO: Possible bug - large early-rejected FASTA inputs can strand FastaStreamerST's producer.
		//FastaStreamerST.close() does not cancel/join; its terminal outputQueue.put at FastaStreamerST.java:232
		//can block after validation exits this loop early; this was not identified as the leaked object.
		final HashMap<String,ProteinSequence> map=new HashMap<String,ProteinSequence>();
		final FileFormat format=FileFormat.testInput(path.toString(),FileFormat.FASTA,null,false,true);
		final Streamer input=StreamerFactory.makeStreamer(format,null,true,-1,false,true,1);
		input.start();
		try{
			for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
				for(Read read : batch.list){
					if(read==null || read.bases==null || read.mate!=null || read.bases.length<1 || read.bases.length>cap){
						throw new IllegalArgumentException("Invalid combined consensus record");
					}
					final String id=firstToken(read.id);
					if(map.put(id,new ProteinSequence(id,read.bases))!=null){
						throw new IllegalArgumentException("Duplicate combined consensus: "+id);
					}
				}
			}
		}finally{if(ReadWrite.closeStream(input)){throw new IllegalStateException("Combined consensus read failed");}}
		return map;
	}

	private static ProteinSequence consensus(Path path,int cap,String family){
		//TODO: Possible bug - large early-rejected FASTA inputs can strand FastaStreamerST's producer.
		//FastaStreamerST.close() does not cancel/join; its terminal outputQueue.put at FastaStreamerST.java:232
		//can block after validation exits this loop early; this was not identified as the leaked object.
		final FileFormat format=FileFormat.testInput(path.toString(),FileFormat.FASTA,null,false,true);
		final Streamer input=StreamerFactory.makeStreamer(format,null,true,-1,false,true,1);
		input.start();
		ProteinSequence found=null;
		try{
			for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
				for(Read read : batch.list){
					if(read==null || read.bases==null || read.mate!=null || read.bases.length<1 ||
							read.bases.length>cap || found!=null){
						throw new IllegalArgumentException("Fixed consensus must contain exactly one bounded protein");
					}
					final String id=firstToken(read.id);
					if(!id.equals(family)){
						throw new IllegalArgumentException("Fixed consensus ID "+id+" differs from models family "+family);
					}
					found=new ProteinSequence(id,read.bases);
				}
			}
		}finally{if(ReadWrite.closeStream(input)){throw new IllegalStateException("Fixed consensus read failed");}}
		if(found==null){throw new IllegalArgumentException("Missing fixed consensus: "+path);}
		return found;
	}

	private static ReferenceFiles loadReferenceManifest(Path manifest) throws Exception{
		final Path dir=manifest.getParent();
		final java.util.List<String> lines=Files.readAllLines(manifest,StandardCharsets.US_ASCII);
		if(lines.size()!=8 || !lines.get(0).equals("file\tsha80")){
			throw new IllegalArgumentException("Require exact full-reference result manifest");
		}
		final HashMap<String,Path> files=new HashMap<String,Path>();
		final HashMap<String,String> hashes=new HashMap<String,String>();
		for(int i=1;i<lines.size();i++){
			final String[] fields=lines.get(i).split("\t",-1);
			if(fields.length!=2){throw new IllegalArgumentException("Malformed reference manifest row");}
			final String logical=logicalReferenceName(fields[0]);
			if(files.put(logical,dir.resolve(fields[0]))!=null || hashes.put(logical,fields[1])!=null){
				throw new IllegalArgumentException("Duplicate reference manifest logical file: "+logical);
			}
			pin(dir.resolve(fields[0]),fields[1]);
		}
		for(String name : REFERENCE_FILES){
			if(!files.containsKey(name)){throw new IllegalArgumentException("Missing reference file: "+name);}
		}
		if(!new String(Files.readAllBytes(files.get("PASS")),StandardCharsets.US_ASCII)
				.equals("FAMILY_FULL_REFERENCES_PASS\n")){
			throw new IllegalArgumentException("Full-member references have not passed");
		}
		return new ReferenceFiles(files.get("profiles.tsv"),hashes.get("profiles.tsv"),
			files.get("spectra.tsv"),hashes.get("summary.tsv"),files.get("members.tsv"));
	}

	private static LowerFiles loadLowerManifest(Path manifest) throws Exception{
		final Path dir=manifest.getParent();
		final java.util.List<String> lines=Files.readAllLines(manifest,StandardCharsets.US_ASCII);
		if(lines.size()!=7 || !lines.get(0).equals("file\tsha80")){
			throw new IllegalArgumentException("Require exact lower result manifest");
		}
		final HashMap<String,Path> files=new HashMap<String,Path>();
		for(int i=1;i<lines.size();i++){
			final String[] fields=lines.get(i).split("\t",-1);
			if(fields.length!=2){throw new IllegalArgumentException("Malformed lower manifest row");}
			final String logical=logicalLowerName(fields[0]);
			if(files.put(logical,dir.resolve(fields[0]))!=null){
				throw new IllegalArgumentException("Duplicate lower manifest logical file: "+logical);
			}
			pin(dir.resolve(fields[0]),fields[1]);
		}
		for(String name : LOWER_FILES){
			if(!files.containsKey(name)){throw new IllegalArgumentException("Missing lower file: "+name);}
		}
		if(!new String(Files.readAllBytes(files.get("PASS")),StandardCharsets.US_ASCII)
				.equals("FAMILY_ALL_MEMBER_CALIBRATION_PASS\n")){
			throw new IllegalArgumentException("Full lower calibration has not passed");
		}
		return new LowerFiles(files.get("cutoffs.tsv"));
	}

	private static UpperFiles loadUpperManifest(Path manifest) throws Exception{
		final Path dir=manifest.getParent();
		final java.util.List<String> lines=Files.readAllLines(manifest,StandardCharsets.US_ASCII);
		if(lines.size()!=6 || !lines.get(0).equals("file\tsha80")){
			throw new IllegalArgumentException("Require exact upper result manifest");
		}
		final HashMap<String,Path> files=new HashMap<String,Path>();
		for(int i=1;i<lines.size();i++){
			final String[] fields=lines.get(i).split("\t",-1);
			if(fields.length!=2){throw new IllegalArgumentException("Malformed upper manifest row");}
			if(!contains(UPPER_FILES,fields[0]) || fields[0].endsWith(".gz") ||
					files.put(fields[0],dir.resolve(fields[0]))!=null){
				throw new IllegalArgumentException("Unexpected upper manifest file: "+fields[0]);
			}
			pin(dir.resolve(fields[0]),fields[1]);
		}
		for(String name : UPPER_FILES){
			if(!files.containsKey(name)){throw new IllegalArgumentException("Missing upper file: "+name);}
		}
		if(!new String(Files.readAllBytes(files.get("PASS")),StandardCharsets.US_ASCII)
				.equals("FAMILY_NEGATIVE_UPPER_CALIBRATION_PASS\n")){
			throw new IllegalArgumentException("Upper calibration has not passed");
		}
		return new UpperFiles(files.get("uppers.tsv"));
	}

	private static GateFiles loadGateManifest(Path manifest) throws Exception{
		final Path dir=manifest.getParent();
		final java.util.List<String> lines=Files.readAllLines(manifest,StandardCharsets.US_ASCII);
		final HashMap<String,String> values=new HashMap<String,String>();
		if(lines.size()!=GATE_KEYS.length+1 || !lines.get(0).equals("key\tvalue")){
			throw new IllegalArgumentException("Require exact A48 gate result manifest");
		}
		for(int i=1; i<lines.size(); i++){
			final String[] fields=lines.get(i).split("\t",-1);
			if(fields.length!=2 || !contains(GATE_KEYS,fields[0]) || values.put(fields[0],fields[1])!=null){
				throw new IllegalArgumentException("Malformed A48 gate manifest row");
			}
		}
		if(!FamilyGateCodec.SCHEMA.equals(values.get("schema"))){
			throw new IllegalArgumentException("A48 gate schema mismatch: "+values.get("schema"));
		}
		DigestSuffix.requireSuffix(values.get("raw_lower_sha80"),"raw lower calibration");
		DigestSuffix.requireSuffix(values.get("raw_upper_sha80"),"raw upper calibration");
		final Path cutoffs=dir.resolve("cutoffs.tsv"), uppers=dir.resolve("uppers.tsv");
		final Path recipe=dir.resolve("recipe.txt"), summary=dir.resolve("summary.tsv"), pass=dir.resolve("PASS");
		pin(cutoffs,values.get("cutoffs.tsv")); pin(uppers,values.get("uppers.tsv"));
		pin(recipe,values.get("recipe.txt")); pin(summary,values.get("summary.tsv"));
		pin(pass,values.get("PASS"));
		if(!new String(Files.readAllBytes(pass),StandardCharsets.US_ASCII).equals(FamilyGateExport.PASS+"\n")){
			throw new IllegalArgumentException("A48 gate export has not passed");
		}
		return new GateFiles(cutoffs,uppers);
	}

	private static CollectionFiles loadCollection(HashMap<String,String> options) throws Exception{
		for(String key : SINGLE_SHARD_OPTIONS){
			if(options.containsKey(key)){
				throw new IllegalArgumentException("Single-shard and collection flags are mutually exclusive: "+key);
			}
		}
		final Path referenceIndex=input(options,"referenceindex"), upperIndex=optionalInput(options,"upperindex");
		final Path gateIndex=optionalInput(options,"gateindex");
		if((upperIndex==null)==(gateIndex==null)){
			throw new IllegalArgumentException("Require exactly one of upperindex or gateindex");
		}
		final HashMap<Integer,IndexRow> references=loadReferenceIndex(referenceIndex,gateIndex==null);
		final HashMap<Integer,IndexRow> gates=gateIndex==null ? loadUpperIndex(upperIndex) : loadGateIndex(gateIndex);
		if(references.size()!=gates.size()){
			throw new IllegalArgumentException("Reference and gate indices cover different bucket counts");
		}
		final ArrayList<ReferenceFiles> refList=new ArrayList<ReferenceFiles>(references.size());
		final ArrayList<Path> lowerList=new ArrayList<Path>(references.size());
		final ArrayList<Path> upperList=new ArrayList<Path>(gates.size());
		final FamilyCensus census=new FamilyCensus();
		for(Integer bucket : references.keySet()){
			final IndexRow refRow=references.get(bucket), gateRow=gates.remove(bucket);
			if(gateRow==null || !refRow.sameCensus(gateRow)){
				throw new IllegalArgumentException("Reference/gate index mismatch for bucket: "+bucket);
			}
			final ReferenceFiles refs=loadReferenceManifest(refRow.referenceManifest);
			final LowerFiles lower;
			final UpperFiles upper;
			if(gateIndex==null){
				lower=loadLowerManifest(refRow.lowerManifest);
				upper=loadUpperManifest(gateRow.upperManifest);
			}else{
				final GateFiles gate=loadGateManifest(gateRow.gateManifest);
				lower=new LowerFiles(gate.cutoffs); upper=new UpperFiles(gate.uppers);
			}
			addCensus(refs.members,census,refRow);
			refList.add(refs); lowerList.add(lower.cutoffs); upperList.add(upper.uppers);
		}
		if(!gates.isEmpty()){
			throw new IllegalArgumentException("Gate index has buckets missing from reference index: "+gates.keySet());
		}
		census.finish();
		return new CollectionFiles(refList.toArray(new ReferenceFiles[refList.size()]),
			lowerList.toArray(new Path[lowerList.size()]),upperList.toArray(new Path[upperList.size()]),census,gateIndex!=null);
	}

	private static void addCensus(Path file,FamilyCensus census,IndexRow expected){
		final ByteFile input=ByteFile.makeByteFile(file.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		int families=0, nonempty=0, zeros=0, members=0;
		long residues=0;
		try{
			header(input,"roster_rank\tfamily_id\tmembers\tresidues\tvalid_words");
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line);
				if(lp.terms()!=5){throw new IllegalArgumentException("Wrong census field count");}
				final int rank=lp.parseInt(0), count=lp.parseInt(2);
				final long aa=Long.parseLong(lp.parseString(3));
				final long valid=Long.parseLong(lp.parseString(4));
				final String id=lp.parseString(1);
				if(count<0 || aa<0 || valid<0 || (count==0)!=(aa==0 && valid==0)){
					throw new IllegalArgumentException("Malformed census row: "+id);
				}
				census.add(rank,id,count,aa,valid);
				families++; members+=count; residues+=aa;
				if(count==0){zeros++;}else{nonempty++;}
			}
		}finally{if(input.close()){throw new IllegalStateException("Census read failed");}}
		if(families!=expected.families || nonempty!=expected.nonempty || zeros!=expected.zeros ||
				members!=expected.members || residues!=expected.residues){
			throw new IllegalArgumentException("Census totals do not match index");
		}
	}

	private static HashMap<Integer,IndexRow> loadReferenceIndex(Path index, boolean loadLower) throws Exception{
		final HashMap<Integer,IndexRow> map=new HashMap<Integer,IndexRow>();
		final ByteFile input=ByteFile.makeByteFile(index.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			header(input,"bucket\tfamilies\tnonempty\tzeros\tmembers\tresidues\treferences_manifest"+
				"\treferences_sha80\tlower_manifest\tlower_sha80\tlineage");
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line);
				if(lp.terms()!=11){throw new IllegalArgumentException("Wrong reference index field count");}
				final IndexRow row=new IndexRow(lp,10);
				row.referenceManifest=path(lp.parseString(6)); pin(row.referenceManifest,lp.parseString(7));
				if(loadLower){
					row.lowerManifest=path(lp.parseString(8)); pin(row.lowerManifest,lp.parseString(9));
				}
				if(map.put(Integer.valueOf(row.bucket),row)!=null){
					throw new IllegalArgumentException("Duplicate reference bucket: "+row.bucket);
				}
			}
		}finally{if(input.close()){throw new IllegalStateException("Reference index read failed");}}
		return map;
	}

	private static HashMap<Integer,IndexRow> loadUpperIndex(Path index) throws Exception{
		final HashMap<Integer,IndexRow> map=new HashMap<Integer,IndexRow>();
		final ByteFile input=ByteFile.makeByteFile(index.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			header(input,"bucket\tfamilies\tnonempty\tzeros\tmembers\tresidues\tupper_rows"+
				"\tupper_manifest\tupper_sha80\tlineage");
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line);
				if(lp.terms()!=10){throw new IllegalArgumentException("Wrong upper index field count");}
				final IndexRow row=new IndexRow(lp,9);
				if(lp.parseInt(6)!=row.nonempty){
					throw new IllegalArgumentException("Upper row count does not match nonempty count");
				}
				row.upperManifest=path(lp.parseString(7)); pin(row.upperManifest,lp.parseString(8));
				if(map.put(Integer.valueOf(row.bucket),row)!=null){
					throw new IllegalArgumentException("Duplicate upper bucket: "+row.bucket);
				}
			}
		}finally{if(input.close()){throw new IllegalStateException("Upper index read failed");}}
		return map;
	}

	private static HashMap<Integer,IndexRow> loadGateIndex(Path index) throws Exception{
		final HashMap<Integer,IndexRow> map=new HashMap<Integer,IndexRow>();
		final ByteFile input=ByteFile.makeByteFile(index.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			header(input,"bucket\tfamilies\tnonempty\tzeros\tmembers\tresidues\tgate_manifest\tgate_sha80\tlineage");
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line);
				if(lp.terms()!=9){throw new IllegalArgumentException("Wrong A48 gate index field count");}
				final IndexRow row=new IndexRow(lp,8,true);
				row.gateManifest=path(lp.parseString(6)); pin(row.gateManifest,lp.parseString(7));
				if(map.put(Integer.valueOf(row.bucket),row)!=null){
					throw new IllegalArgumentException("Duplicate A48 gate bucket: "+row.bucket);
				}
			}
		}finally{if(input.close()){throw new IllegalStateException("A48 gate index read failed");}}
		return map;
	}

	private static Path path(String value){return Paths.get(value).toAbsolutePath().normalize();}

	private static String logicalReferenceName(String name){
		if(name.equals("profiles.tsv.gz")){return "profiles.tsv";}
		if(name.equals("spectra.tsv.gz")){return "spectra.tsv";}
		for(String allowed : REFERENCE_FILES){if(name.equals(allowed)){return name;}}
		throw new IllegalArgumentException("Unexpected reference manifest file: "+name);
	}

	private static String logicalLowerName(String name){
		if(name.equals("cutoffs.tsv.gz")){return "cutoffs.tsv";}
		if(name.equals("spectra.tsv.gz")){return "spectra.tsv";}
		for(String allowed : LOWER_FILES){if(name.equals(allowed)){return name;}}
		throw new IllegalArgumentException("Unexpected lower manifest file: "+name);
	}

	private static boolean contains(String[] array,String value){
		for(String s : array){if(s.equals(value)){return true;}}
		return false;
	}

	private static void header(ByteFile input,String value){
		if(!Arrays.equals(input.nextLine(),value.getBytes(StandardCharsets.US_ASCII))){
			throw new IllegalArgumentException("Unexpected input header: "+value);
		}
	}

	private static Path input(HashMap<String,String> options,String key) throws Exception{
		final Path p=Paths.get(required(options,key)).toAbsolutePath().normalize();
		pin(p,required(options,key+"sha80"));
		return p;
	}

	/** Structural identities may be inherited pins or explicit unmeasured markers; neither is recomputed. */
	private static Path input(HashMap<String,String> options,String key,
			DDLFamilySearch.StructuralInputs structuralInputs)throws Exception{
		if(structuralInputs==null){return input(options,key);}
		final String identity=required(options,key+"sha80");
		if(!"unmeasured".equals(identity)){DigestSuffix.requireSuffix(identity,"Inherited input "+key);}
		final Path p=Paths.get(required(options,key)).toAbsolutePath().normalize();
		structuralInputs.track(p);return p;
	}

	private static Path optionalInput(HashMap<String,String> options,String key) throws Exception{
		final String value=options.get(key);
		if(value==null){return null;}
		final String hash=required(options,key+"sha80");
		final Path p=Paths.get(value).toAbsolutePath().normalize();
		pin(p,hash);
		return p;
	}

	private static void pin(Path path,String hash) throws Exception{
		if(!Files.isRegularFile(path) || Files.isSymbolicLink(path)){
			throw new IllegalArgumentException("Require a regular input: "+path);
		}
		if(!DigestSuffix.file(path.toString()).equals(DigestSuffix.requireSuffix(hash,"input sha80"))){
			throw new IllegalArgumentException("Input pin mismatch: "+path);
		}
	}

	private static void pinConfiguredInputs(HashMap<String,String> options,boolean collection,
			boolean familyReferences,Path models,Path in) throws Exception{
		pin(in,required(options,"insha80"));
		if(familyReferences){
			pin(Paths.get(required(options,"familyreferences")).toAbsolutePath().normalize(),
				required(options,"familyreferencessha80"));
		}else if(collection){
			pin(models,required(options,"modelssha80"));
			pin(Paths.get(required(options,"referenceindex")).toAbsolutePath().normalize(),
				required(options,"referenceindexsha80"));
			if(options.containsKey("gateindex")){
				pin(Paths.get(required(options,"gateindex")).toAbsolutePath().normalize(),required(options,"gateindexsha80"));
			}else{
				pin(Paths.get(required(options,"upperindex")).toAbsolutePath().normalize(),required(options,"upperindexsha80"));
			}
			pin(Paths.get(required(options,"consensus")).toAbsolutePath().normalize(),required(options,"consensussha80"));
		}else{
			pin(models,required(options,"modelssha80"));
			pin(Paths.get(required(options,"references")).toAbsolutePath().normalize(),required(options,"referencessha80"));
			pin(Paths.get(required(options,"lower")).toAbsolutePath().normalize(),required(options,"lowersha80"));
			pin(Paths.get(required(options,"upper")).toAbsolutePath().normalize(),required(options,"uppersha80"));
		}
	}

	private static long number(HashMap<String,String> o,String key){return Long.parseLong(required(o,key));}

	private static int integer(HashMap<String,String> o,String key){return Math.toIntExact(number(o,key));}

	private static int nonnegativeInt(LineParser1 lp,int term){
		final int length=lp.setBounds(term);
		if(length<1 || length>10){throw new IllegalArgumentException("A48 count/rank must fit a nonnegative int");}
		for(int i=lp.a(); i<lp.b(); i++){
			if(lp.line()[i]<'0' || lp.line()[i]>'9'){
				throw new IllegalArgumentException("A48 count/rank must be decimal");
			}
		}
		final long value=lp.parseLong(term);
		if(value>Integer.MAX_VALUE){throw new IllegalArgumentException("A48 count/rank exceeds signed32 range");}
		return (int)value;
	}

	private static boolean startInputAndWait(final Streamer input,final ArrayList<ProcessThread> workers,
			final CascadeAccumulator accumulator){
		try{
			input.start();
			final ArrayList<ProcessThread> started=new ArrayList<ProcessThread>(workers.size());
			if(SIMULATE_FIRST_START_FAILURE){throw new IllegalStateException("Synthetic FamilyTailCascade first-start failure");}
			for(ProcessThread pt : workers){
				pt.start();
				started.add(pt);
			}
			return ThreadWaiter.waitForThreadsToFinish(started,accumulator);
		}catch(Throwable t){
			t.printStackTrace();
			accumulator.failure.compareAndSet(null,t);
			KillSwitch.kill("FamilyTailCascade failed during worker startup");
		}
		return false;
	}

	private static int threads(HashMap<String,String> o){
		final Integer t=optionalInteger(o,"t"), threads=optionalInteger(o,"threads");
		if(t!=null && threads!=null && t.intValue()!=threads.intValue()){
			throw new IllegalArgumentException("Conflicting t and threads values");
		}
		final int value=(t!=null ? t.intValue() : threads!=null ? threads.intValue() : 1);
		if(value<1){throw new IllegalArgumentException("Threads must be positive: "+value);}
		return value;
	}

	private static Integer optionalInteger(HashMap<String,String> o,String key){
		final String value=o.get(key);
		return value==null ? null : Integer.valueOf(Math.toIntExact(Long.parseLong(value)));
	}

	private static boolean optionalBoolean(HashMap<String,String> o,String key,boolean defaultValue){
		final String value=o.get(key);
		return value==null ? defaultValue : Parse.parseBoolean(value);
	}

	private static void validateCompactReferenceUse(CosineMode cosineMode,RankingFeatures features,Path rankOut){
		if(cosineMode!=CosineMode.RANKING || features!=RankingFeatures.COMBINED25 || rankOut!=null){
			throw new IllegalArgumentException(
				"Compact 20-coordinate familyreferences support only combined25 ranking without rankout/truth diagnostics");
		}
	}

	/** Omitted cosine uses ranking only; explicit boolean aliases retain their historical gates. */
	private static CosineMode cosineMode(HashMap<String,String> options){
		final String value=options.get("cosine");
		if(value==null || value.equalsIgnoreCase("ranking")){return CosineMode.RANKING;}
		if(value.equalsIgnoreCase("full")){return CosineMode.FULL;}
		if(value.equalsIgnoreCase("off")){return CosineMode.OFF;}
		return Parse.parseBoolean(value) ? CosineMode.FULL : CosineMode.OFF;
	}

	enum CosineMode{FULL,RANKING,OFF}
	enum RankingFeatures{DIMER,MONOMER,CENTERED_MONOMER,COMBINED25}

	private static RankingFeatures rankingFeatures(HashMap<String,String> options,CosineMode mode){
		final String value=options.get("rankingfeatures");
		if(value==null || value.equalsIgnoreCase("dimer")){return RankingFeatures.DIMER;}
		final boolean centered=value.equalsIgnoreCase("centered_monomer");
		final boolean combined=value.equalsIgnoreCase("combined25");
		if(!centered && !combined && !value.equalsIgnoreCase("monomer")){
			throw new IllegalArgumentException("rankingfeatures requires dimer, monomer, centered_monomer or combined25");
		}
		if(mode!=CosineMode.RANKING){
			throw new IllegalArgumentException("rankingfeatures="+
				(combined ? "combined25" : centered ? "centered_monomer" : "monomer")+
				" requires cosine=ranking");
		}
		return combined ? RankingFeatures.COMBINED25 : centered ? RankingFeatures.CENTERED_MONOMER : RankingFeatures.MONOMER;
	}

	private static boolean containsAny(HashMap<String,String> options,String[] keys){
		for(String key : keys){if(options.containsKey(key)){return true;}}
		return false;
	}

	private static void require(boolean condition,String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	private static void requireOnly(HashMap<String,String> options,String[] keys){
		if(options==null){throw new IllegalArgumentException("Missing options");}
		for(String key : options.keySet()){
			boolean found=false;
			for(String allowed : keys){if(allowed.equals(key)){found=true; break;}}
			if(!found){throw new IllegalArgumentException("Unknown callable option: "+key);}
		}
	}

	private static String required(HashMap<String,String> o,String key){
		final String value=o.get(key);
		if(value==null || value.isEmpty()){throw new IllegalArgumentException("Missing option: "+key);}
		return value;
	}

	static ByteStreamWriter writer(Path path){
		final FileFormat format=FileFormat.testOutput(path.toString(),FileFormat.TEXT,null,false,false,false,true);
		return new ByteStreamWriter(format);
	}

	private static void write(Path path,String value) throws Exception{
		Files.write(path,value.getBytes(StandardCharsets.US_ASCII));
	}

	private static String seconds(long nanos){return String.format(Locale.ROOT,"%.6f",nanos*1e-9);}

	static String firstToken(String id){
		if(id==null || id.isEmpty()){throw new IllegalArgumentException("Missing FASTA ID");}
		int end=0; while(end<id.length() && !Character.isWhitespace(id.charAt(end))){end++;}
		if(end==0){throw new IllegalArgumentException("Missing FASTA ID");}
		return id.substring(0,end);
	}

	private static final class ReferenceFiles{
		ReferenceFiles(Path profiles_,String profilesHash_,Path spectra_,String summaryHash_,Path members_){
			profiles=profiles_; profilesHash=profilesHash_; spectra=spectra_; summaryHash=summaryHash_; members=members_;
		}
		final Path profiles,spectra,members;
		final String profilesHash,summaryHash;
	}

	private static final class LowerFiles{
		LowerFiles(Path cutoffs_){cutoffs=cutoffs_;}
		final Path cutoffs;
	}

	private static final class UpperFiles{
		UpperFiles(Path uppers_){uppers=uppers_;}
		final Path uppers;
	}

	private static final class GateFiles{
		GateFiles(Path cutoffs_,Path uppers_){cutoffs=cutoffs_; uppers=uppers_;}
		final Path cutoffs,uppers;
	}

	private static final class CollectionFiles{
		CollectionFiles(ReferenceFiles[] references_,Path[] lowers_,Path[] uppers_,FamilyCensus census_,boolean a48_){
			references=references_; lowers=lowers_; uppers=uppers_; census=census_; a48=a48_;
		}
		final ReferenceFiles[] references;
		final Path[] lowers,uppers;
		final FamilyCensus census;
		final boolean a48;
	}

	private static final class LinkedSpectra{
		LinkedSpectra(ProteinSpectrum.Spectrum[] values_){values=values_;}
		final ProteinSpectrum.Spectrum[] values;
	}

	private static final class Models{
		Models(Model[] values_,String[] canonicalIds_){values=values_; canonicalIds=canonicalIds_;}
		final Model[] values;
		final String[] canonicalIds;
	}

	private static final class Model{
		Model(int rank_,String id_,ProteinSequence sequence_){rank=rank_; id=id_; sequence=sequence_;}
		final int rank;
		final String id;
		final ProteinSequence sequence;
	}

	private static final class FamilyCensus{
		void add(int rank,String id,int members,long residues,long validWords){
			if(rank<0 || id==null || id.isEmpty() || members<0 || residues<0 || validWords<0){
				throw new IllegalArgumentException("Bad census row: "+id);
			}
			while(ids.size()<=rank){ids.add(null); counts.add(null);}
			if(ids.set(rank,id)!=null || byId.put(id,Integer.valueOf(rank))!=null){
				throw new IllegalArgumentException("Duplicate census family: "+id);
			}
			counts.set(rank,new Entry(id,members,residues,validWords));
			if(members==0){zeros++;}else{nonempty++;}
		}
		void finish(){
			for(int i=0; i<ids.size(); i++){
				if(ids.get(i)==null){throw new IllegalArgumentException("Missing census rank: "+i);}
			}
		}
		Entry get(int rank,String id){
			if(rank<0 || rank>=ids.size() || !ids.get(rank).equals(id)){
				throw new IllegalArgumentException("Family/rank not in current census: "+rank+"/"+id);
			}
			return counts.get(rank);
		}
		Entry get(String id){
			final Integer rank=byId.get(id);
			if(rank==null){throw new IllegalArgumentException("Family not in current census: "+id);}
			return counts.get(rank.intValue());
		}
		int members(String id){
			return get(id).members;
		}
		String id(int rank){return ids.get(rank);}
		int size(){return ids.size();}
		final ArrayList<String> ids=new ArrayList<String>();
		final ArrayList<Entry> counts=new ArrayList<Entry>();
		final HashMap<String,Integer> byId=new HashMap<String,Integer>();
		int nonempty,zeros;
		private static final class Entry{
			Entry(String id_,int members_,long residues_,long validWords_){
				id=id_; members=members_; residues=residues_; validWords=validWords_;
			}
			final String id;
			final int members;
			final long residues,validWords;
		}
	}

	private static final class IndexRow{
		IndexRow(LineParser1 lp,int lineageColumn){
			this(lp,lineageColumn,false);
		}
		IndexRow(LineParser1 lp,int lineageColumn,boolean checked){
			bucket=checked ? nonnegativeInt(lp,0) : lp.parseInt(0);
			families=checked ? nonnegativeInt(lp,1) : lp.parseInt(1);
			nonempty=checked ? nonnegativeInt(lp,2) : lp.parseInt(2);
			zeros=checked ? nonnegativeInt(lp,3) : lp.parseInt(3);
			members=checked ? nonnegativeInt(lp,4) : lp.parseInt(4);
			residues=Long.parseLong(lp.parseString(5));
			if(bucket<0 || families<0 || nonempty<0 || zeros<0 || members<0 || residues<0 ||
					families!=nonempty+zeros || (!lp.parseString(lineageColumn).equals("NEW") &&
					!lp.parseString(lineageColumn).equals("REUSED_ACCEPTED"))){
				throw new IllegalArgumentException("Malformed collection index row");
			}
		}
		boolean sameCensus(IndexRow b){
			return bucket==b.bucket && families==b.families && nonempty==b.nonempty &&
				zeros==b.zeros && members==b.members && residues==b.residues;
		}
		final int bucket,families,nonempty,zeros,members;
		final long residues;
		Path referenceManifest,lowerManifest,upperManifest,gateManifest;
	}

	private static final class Gates{
		Gates(Gate[] values_){values=values_;}
		final Gate[] values;
	}

	private static final class Gate{
		Gate(double cosineLower_,double identityLower_,double blosumLower_){
			if(!fraction(cosineLower_) || !percent(identityLower_) || !finite(blosumLower_)){
				throw new IllegalArgumentException("Lower gates must be finite and within their metric domains");
			}
			cosineLower=cosineLower_; identityLower=identityLower_; blosumLower=blosumLower_;
		}
		void setRawUpper(double cosineUpper_,double identityUpper_,double blosumUpper_){
			cosineUpperRaw=cosineUpper_; identityUpperRaw=identityUpper_;
			cosineUpper=FamilyGateCodec.cosineAccept(cosineLower,cosineUpper_);
			identityUpper=FamilyGateCodec.upper(identityLower,identityUpper_);
			blosumFinal=FamilyGateCodec.finalBlosum(blosumLower,blosumUpper_);
			upperSet=true;
		}
		void setFinalUpper(double cosineUpper_,double identityUpper_,double blosumFinal_){
			cosineUpperRaw=cosineUpper_; identityUpperRaw=identityUpper_;
			if(!finite(cosineUpper_) || cosineUpper_<0.72 || cosineUpper_<cosineLower){
				throw new IllegalArgumentException("Invalid final cosine gate");
			}
			cosineUpper=cosineUpper_;
			if(!finite(identityUpper_) || identityUpper_<Math.max(identityLower,30) ||
					identityUpper_>Math.nextUp(110f)){
				throw new IllegalArgumentException("Final identity gate must preserve30percent floor and member lower");
			}
			identityUpper=identityUpper_;
			if(!finite(blosumFinal_) || blosumFinal_<blosumLower){
				throw new IllegalArgumentException("Invalid final BLOSUM gate");
			}
			blosumFinal=blosumFinal_;
			upperSet=true;
		}
		boolean complete(){return upperSet;}
		final double cosineLower, identityLower, blosumLower;
		double cosineUpperRaw,identityUpperRaw;
		double cosineUpper,identityUpper,blosumFinal;
		boolean upperSet;
	}

	private static final String[] REFERENCE_FILES={
		"profiles.tsv","spectra.tsv","members.tsv","summary.tsv","recipe.txt","timing.tsv","PASS"
	};
	private static final String[] LOWER_FILES={"cutoffs.tsv","spectra.tsv","summary.tsv","recipe.txt","timing.tsv","PASS"};
	private static final String[] UPPER_FILES={"uppers.tsv","summary.tsv","timing.tsv","recipe.txt","PASS"};
	private static final String[] SINGLE_SHARD_OPTIONS={
		"references","referencessha80","lower","lowersha80","upper","uppersha80"
	};
	private static final String[] COLLECTION_OPTIONS={
		"referenceindex","referenceindexsha80","upperindex","upperindexsha80",
		"gateindex","gateindexsha80","consensus","consensussha80"};
	private static final String[] LEGACY_MODEL_OPTIONS={"models","modelssha80","maxreferenceaa"};
	private static final byte[] OUTPUT_HEADER=
		"query_id\tstatus\tfamily_id\troster_rank\tddl_matches\tcosine\tidentity_pct\tblosum_a_sq\n"
		.getBytes(StandardCharsets.US_ASCII);
	private static final byte[] UNION_PAIR_HEADER=("query_id\troster_rank\tfamily_id"+
		"\tddl_matches\tddl_position\twindow_position\twindow_size\t"+
		"in_deployed_top\tcombined_cosine\tidentity_pct\tblosum_a_sq\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] CANDIDATE_HEADER=("query_id\troster_rank\tfamily_id"+
		"\tddl_matches\tddl_position\twindow_position\twindow_size\t"+
		"in_deployed_top\tcombined_cosine\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] AUDIT_HEADER=("query_id\tterminal_stage\tddl_active_count"+
		"\tcosine_window_size\tdp_shortlist_size\tdp_attempts\t"+
		"identity_reject_count\tblosum_reject_count\tcosine_ranked_to_aligned\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FLOOR_AUDIT_HEADER=(new String(AUDIT_HEADER,StandardCharsets.US_ASCII).trim()+
		"\tdeferred_identity_family\tdeferred_identity_pct\tdeferred_blosum_a_sq"+
		"\tdeferred_blosum_pass\n").getBytes(StandardCharsets.US_ASCII);
	private static final String[] TRUTH_STATUS={
		"absent_mapping","absent_family","ddl_rejected","outside_window","in_window","in_top"
	};
	private static final String[] RANK_MODES={"ddl","dimer","monomer","centered_monomer","combined25"};
	private static final byte[] RANK_HEADER=("query_id\texpected_family_id\ttruth_boundary"+
		"\twindow_size\ttop\tddl_rank\tdimer_rank\tmonomer_rank\t"+
		"ddl_status\tdimer_status\tmonomer_status\tassigned_family_id\twinner_ddl_position\twinner_dp_attempt\tdp_attempts\t"+
		"centered_monomer_rank\tcentered_monomer_status\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] COMBINED_RANK_HEADER=(new String(RANK_HEADER,StandardCharsets.US_ASCII).trim()+
		"\tcombined25_rank\tcombined25_status\n").getBytes(StandardCharsets.US_ASCII);
	private static final java.util.List<String> OPTIONS=Arrays.asList("references","referencessha80","lower","lowersha80",
		"upper","uppersha80","models","modelssha80","in","insha80","out","top","k","seed","width","maxreads",
		"maxreferenceaa","ddlreject","ddlaccept","referenceindex","referenceindexsha80","upperindex","upperindexsha80",
		"gateindex","gateindexsha80","consensus","consensussha80","familyreferences","familyreferencessha80",
		"maxfamilies","cosine","rankingfeatures","rankout","auditout","truth","truthsha80","t","threads",
		"aabits","aapartition","scalarbounds","scalarboundssha80","scalarscale",
		"scalargates","scalarz","mode","supplementalout","stagecpu","identityacceptfloor",
		"blosumacceptfloor","flooraudit","candidateonly","positional","positionmode",
		"ddlprofiles","ddlprofilessha80","partitions","overlap","blacklist","blacklistsha80",
		"alignmentmode","alignmentminspan","alignmentminpaired","blacklistk","blacklistmode",
		"blacklistwords","blacklistprovenance","gridaxes","gridcohort","gridwitnessout",
		"gridoffset","gridwitnessselected","indexearlyreject","indexdirect","indexorder",
		"minmemberaa","maxmemberaa");
	private static final String[] CALLABLE_OPTIONS={"familyreferences","familyreferencessha80","top","k","seed","width",
		"ddlreject","ddlaccept","cosine","rankingfeatures","aabits","aapartition",
		"scalarbounds","scalarboundssha80","scalarscale","scalargates","scalarz",
		"identityacceptfloor","blosumacceptfloor","positional","positionmode","ddlprofiles",
		"ddlprofilessha80","partitions","overlap","blacklist","blacklistsha80","blacklistk",
		"blacklistmode","blacklistwords","blacklistprovenance","indexearlyreject",
		"indexdirect","indexorder","alignmentmode","alignmentminspan","alignmentminpaired",
		"minmemberaa","maxmemberaa","matrix","matrixfile","matrixfilesha80","matrixbinding",
		"heuristic","minqueryaa","maxqueryaa","familygclambda","gcfilter","gc_z"};
	private static final String[] COMPACT_EMBEDDED_OPTIONS={
		"heuristic","minqueryaa","maxqueryaa","familygclambda","gcfilter","gc_z"};
	private static final int DIMER_WIDTH=400,MONOMER_WIDTH=20,COMBINED_WIDTH=25;
	private static final String[] SCALAR_NAMES={"hydrophobicity","net_charge","aromaticity","fraction_small","pI_tail"};
	private static final String[] GATE_KEYS={"schema","raw_lower_sha80","raw_upper_sha80","cutoffs.tsv",
		"uppers.tsv","recipe.txt","summary.tsv","PASS"};
	private static final boolean SIMULATE_FIRST_START_FAILURE=Boolean.getBoolean("bbgap.familytailcascade.failfirststart");
}
