package prot;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import java.util.TreeSet;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicReference;

import ddl.DDLFamilyProfiles;
import ddl.DDLFamilySearch;
import ddl.DDLKeyRarity;
import ddl.DDLKeyTieBreak;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;
import ddl.DDLPentamerBlacklist;
import ddl.ProteinKmerGeometry;
import ddl.ProteinPositionMode;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.Parse;
import prot.FamilyCalibration.MemberInput;
import prot.FamilyCalibration.MemberInputs;
import prot.FamilyMemberPreparation.Prepared;
import prot.FamilyNegativeUpperCalibration.NegativePanel;
import prot.FamilyNegativeUpperCalibration.PackagedCalibration;
import prot.FamilyNegativeUpperCalibration.PackagedWorker;
import prot.FamilyNegativeUpperCalibration.MetricWorker;
import prot.FamilyNegativeUpperCalibration.MetricPairObserver;
import prot.FamilyNegativeUpperCalibration.MetricResult;
import prot.FamilyPanelRunner.Panel;
import structures.ByteBuilder;

/** Builds initial family gates from pinned members and explicitly named negative panels.
 * One coordinator owns input/output and barriers; a bounded pool owns independent family scratch.
 * Prepared payloads outlive member buckets, and each negative panel is loaded once for all families.
 * No positives, new calibration policy or provisional success markers enter this path.
 * @author Brian Bushnell, Keqing
 */
final class FamilyPreparation{

	FamilyPreparation(Path manifest){
		this(manifest, FamilyCalibration.VerificationMode.SHA80);
	}

	/** Policy belongs to this preparation; inherited pins are never new output digests. */
	FamilyPreparation(Path manifest,FamilyCalibration.VerificationMode verificationMode_){
		require(verificationMode_!=null, "Preparation requires an explicit verification policy");
		verificationMode=verificationMode_;
		path=manifest.toAbsolutePath().normalize();
		if(verificationMode==FamilyCalibration.VerificationMode.SHA80){pin(path, DigestSuffix.file(path.toString()));}
		else{
			require(Files.isRegularFile(path) && !Files.isSymbolicLink(path),
				"Structural preparation requires a regular recipe: "+path);
			pins.put(path, "unmeasured");
			rememberStructuralFile(path);
		}
		final ByteFile input=ByteFile.makeByteFile(path.toString(), false);
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				final String text=new String(line, StandardCharsets.UTF_8).trim();
				if(text.isEmpty() || text.startsWith("#")){continue;}
				final String[] pair=Parse.splitOnFirst(text, '=');
				require(pair.length==2, "Preparation requires key=value: "+text);
				final String key=pair[0].trim().toLowerCase(Locale.ROOT), value=pair[1].trim();
				require(!value.isEmpty() && ALLOWED.contains(key) && options.put(key, value)==null,
					"Unknown, empty or duplicate preparation option: "+key);
			}
		}finally{require(!input.close(), "Preparation manifest read failed: "+path);}
		for(String key:REQUIRED){required(key);}
		require(required("schema").equals("family_preparation_v1"), "Unsupported preparation schema");
		FamilyBlacklistBinding.applyDefaults(options);
		FamilyBlacklistBinding.validateOptions(options,verificationMode);
		for(String key:INPUTS){resolveInput(key, true);}
		resolveInput("blacklist", false);
		resolveInput("lengthmedians", false); resolveInput("familyhydro", false);
		resolveInput("matrixfile",false);
		if(!options.containsKey("matrixfile") && !options.containsKey("matrixbinding") &&
			(!options.containsKey("matrix") || ProteinSubstitutionMatrix.DEFAULT_NAME.equals(options.get("matrix")))){
			pin(java.nio.file.Paths.get(ProteinSubstitutionMatrix.defaultResourcePath()),
				ProteinSubstitutionMatrix.POSSUM55_FILE_SHA80);
		}
		final ProteinSubstitutionMatrix matrix=FamilyScoringPackage.resolveNewWorkMatrix(options);
		FamilyScoringPackage.validateRankingOptions(options,verificationMode);
		require(options.containsKey("alignmentmode")==options.containsKey("alignmentminspan"),
			"Alignment mode and minimum span must be explicitly paired");
		require(!options.containsKey("alignmentminpaired") || options.containsKey("alignmentmode"),
			"Paired minimum requires an explicit alignment mode and span");
		FamilyAlignmentMetric selectedMetric=FamilyAlignmentMetric.LEGACY;
		if(options.containsKey("alignmentmode")){
			FamilyAlignmentMetric.Mode selectedMode=null;
			for(FamilyAlignmentMetric.Mode mode:FamilyAlignmentMetric.Mode.values()){
				if(mode.label.equals(required("alignmentmode"))){selectedMode=mode;}
			}
			require(selectedMode!=null,"Unknown preparation alignment metric");
			selectedMetric=new FamilyAlignmentMetric(selectedMode,integer("alignmentminspan",0,1000),matrix,
				options.containsKey("alignmentminpaired") ? integer("alignmentminpaired",0,1000) : 0);
		}
		require(matrix.same(ProteinSubstitutionMatrix.BLOSUM62) || options.containsKey("alignmentmode"),
			"Custom matrix preparation requires a pinned alignment mode and span");
		alignmentMetric=selectedMetric; metricAware=!alignmentMetric.same(FamilyAlignmentMetric.LEGACY);
		require(metricAware==options.containsKey("heuristic"),
			"Explicit-metric preparation requires the common pivot; legacy preparation retains its original ordering");
		// Explicit fuller-coverage experiment: align the union of query-side raw
		// top100 and family-side raw top1000 negatives. The accepted November 2026
		// run selected 414,311,738 SHORT and 442,229,284 LONG pairs, then stayed
		// study-only after the fuller/+1 replay. See
		// results/fuller_pair_and_gate_diagnostics_20261116.
		rawPairUnion=options.containsKey("negativeselection");
		if(rawPairUnion){
			require(metricAware && required("negativeselection").equals("raw_pair_union"),
				"Raw-pair union selection is an explicit metric-aware preparation mode");
		}
		families=integer("expectedfamilies", 1, 1000000); members=number("expectedmembers", 0);
		buckets=integer("buckets", 1, Integer.MAX_VALUE); batches=integer("batches", 1, Integer.MAX_VALUE);
		require((long)buckets*batches<=Integer.MAX_VALUE, "Preparation chunk grid exceeds integer indexing");
		maxReference=integer("maxreferenceaa", 1, 1000);
		maxBucketMembers=number("maxbucketmembers", 1); maxBucketBytes=number("maxbucketbytes", 1);
		maxPreparedBytes=number("maxpreparedbytes", 1); maxOutputBytes=number("maxoutputbytes", 1);
		final String policy=required("retention");
		require(policy.equals("metric_specific") || policy.equals("uniform99"), "Unknown member retention policy: "+policy);
		retention=policy.equals("metric_specific") ? FamilyCalibration.MemberRetention.METRIC_SPECIFIC :
			FamilyCalibration.MemberRetention.UNIFORM_99;
		upperTop=integer("uppertop", 1, 1000000);
		memberRange=MemberLengthRange.from(options);
		minNegative=integer("minnegativeaa", 40, 200); maxNegative=integer("maxnegativeaa", minNegative, 200);
		require(!(options.containsKey("aabits") && options.containsKey("aapartition")),
			"Use either aabits or aapartition, not both");
		geometry=ProteinKmerGeometry.make(integer("k", 4, 6),
			options.containsKey("aabits") ? options.get("aabits") : options.get("aapartition"));
		width=integer("width", 512, 1024);
		require(number("seed", 0)==0 && (width==512 || width==1024),
			"Preparation uses the accepted K4..K6/seed0 kernels at width512 or1024");
		positionMode=ProteinPositionMode.resolve(options.get("positionmode"), null,
			Boolean.valueOf(bool("positional")), ProteinPositionMode.FIXED);
		require(bool("monomer"), "This preparation format requires monomer=t");
		partitions=integer("partitions", 1, Integer.MAX_VALUE); overlap=integer("overlap", 0, Integer.MAX_VALUE);
		if(positionMode.positional()){DDLFamilyProfiles.validatePositionalLayout(partitions, overlap, geometry);}
		else{require(partitions==1 && overlap==0, "Position-free preparation requires partitions=1 and overlap=0");}
		require(required("cosine").equals("ranking") && required("rankingfeatures").equals("combined25"),
			"Compact preparation requires cosine=ranking rankingfeatures=combined25");
		integer("top", 1, Integer.MAX_VALUE);
		final int reject=integer("ddlreject", 0, Integer.MAX_VALUE), accept=integer("ddlaccept", 1, Integer.MAX_VALUE);
		require(accept>reject, "Deployment requires ddlaccept greater than ddlreject");
		if(metricAware){
			require(reject>=1 && reject<=width,
				"Metric-aware preparation requires a positive deployed raw-hit floor within the sketch width");
		}
		decimal("identityacceptfloor", 30, 100); decimal("blosumacceptfloor", 0, Float.MAX_VALUE);
		blosumFloor=decimal("blosumfloor", 0, Float.MAX_VALUE); scalarZ=decimal("scalarz", 0, Double.MAX_VALUE);
		identityUpperMultiplier=identityPadding("identityuppermultiplier",
			FamilyGateCodec.DEFAULT_IDENTITY_UPPER_MULTIPLIER,1);
		identityUpperAdditive=identityPadding("identityupperadditive",FamilyGateCodec.DEFAULT_IDENTITY_UPPER_ADDITIVE,0);
		require(scalarZ==0 || scalarZ==3.2905 || scalarZ==3.8906, "Unsupported scalar z; use0,3.2905 or3.8906");
		require(!bool("scalargates") || scalarZ>0, "Enabled scalar gates require positive scalarz");
		final double scale=decimal("scalarscale", 0, Double.MAX_VALUE);
		require(scale==1 || scale==0.1, "Combined25 uses scalarscale=1 or0.1");
		final boolean early=bool("indexearlyreject"), direct=options.containsKey("indexdirect") && bool("indexdirect");
		require(!early || !direct, "Direct and early-reject index experiments are mutually exclusive");
		require(positionMode.positional() || !early && !direct, "Experimental index policies require positional profiles");
		if(metricAware){
			require(positionMode.fixedCoordinates() && partitions==4 && overlap==2 && accept==1025 && !early && !direct,
				"Released metric preparation requires P4/O2, exhaustive index and alignment-only acceptance");
			require(Integer.parseInt(required("minqueryaa"))<=minNegative &&
				Integer.parseInt(required("maxqueryaa"))>=maxNegative,
				"Every calibration negative must fit the packaged ranking query bounds before any scoring begins");
		}
		DDLFamilySearch.Index.parseScanOrder(options.containsKey("indexorder") ? options.get("indexorder") : "plain");
		final ArrayList<Panel> list=new ArrayList<Panel>(); final ArrayList<Integer> counts=new ArrayList<Integer>();
		readPanels(list, counts);
		panels=list.toArray(new Panel[list.size()]); queryCounts=new int[panels.length];
		int selected=-1;
		for(int i=0; i<panels.length; i++){
			queryCounts[i]=counts.get(i);
			if(panels[i].name.equals(required("cosinepanel"))){selected=i;}
		}
		require(selected>=0, "cosinepanel must name one declared negative panel"); cosinePanel=selected;
		verify();
	}

	/** Builds the initial classifier in memory; retained profile/audit outputs are final artifacts, not phase inputs. */
	FamilyTailCascade.LoadedClassifier build(Path output, int threads, PrintStream log) throws Exception{
		require(!built && !Files.exists(output) && threads>0,
			"Preparation requires a fresh output and positive worker count, once per instance");
		built=true; verify();
		final long started=System.nanoTime();
		final DDLPentamerBlacklist blacklist=FamilyBlacklistBinding.load(options,verificationMode);
		FamilyBlacklistBinding.validateGeometry(blacklist,geometry);
		final long consensusMemberBonus=options.containsKey("consensusbonus") ? number("consensusbonus",0) : 0;
		if(options.containsKey("keytiebreak") || consensusMemberBonus>0){
			final DDLKeyTieBreak.Policy policy=options.containsKey("keytiebreak") ?
				DDLKeyTieBreak.Policy.parse(required("keytiebreak")) : DDLKeyTieBreak.Policy.SMALLEST;
			final boolean raw=policy==DDLKeyTieBreak.Policy.RAW_TAIL_RARITY;
			final Path manifest=raw ? inputPath("rawmanifest") : null;
			if(raw){keyRawInputs=FamilyKeyPolicyInputs.declaredInputs(manifest,keyInputs);}
			keyPolicy=FamilyKeyPolicyInputs.load(policy,manifest,raw ? number("expectedrawtails",1) : 0,
				raw ? integer("rawminlen",1,Integer.MAX_VALUE) : 1,
				raw ? integer("rawmaxlen",1,Integer.MAX_VALUE) : 1,
				raw ? integer("rawtailcrop",0,Integer.MAX_VALUE) : 0,
				width,geometry,0,positionMode,partitions,overlap,blacklist,keyInputs,log,null,consensusMemberBonus);
		}
		final FamilyTailCascade.ScalarBounds bounds=FamilyTailCascade.combinedBounds(classifierOptions(), verificationMode);
		final MemberInputs source=new MemberInputs(inputPath("families"), required("familiessha80"),
			inputPath("chunks"), required("chunkssha80"), inputPath("consensus"), required("consensussha80"),
			families, members, buckets, batches, maxReference, memberRange, verificationMode);
		final Prepared[] prepared=new Prepared[families];
		final int workers=Math.min(threads, families);
		final FamilyMemberPreparation[] memberWorkers=new FamilyMemberPreparation[workers];
		for(int i=0; i<workers; i++){
			memberWorkers[i]=new FamilyMemberPreparation(geometry, partitions, overlap, blacklist, retention, true,
				memberRange, positionMode, alignmentMetric, width,keyPolicy);
		}
		// Reserve panel references and detached result/count storage before allocating it.
		// This is a conservative payload budget, not a claim about total JVM heap usage.
		final long evidenceBytes=metricAware ? Math.multiplyExact(Math.multiplyExact((long)panels.length,families),160L) : 0;
		require(evidenceBytes<=maxPreparedBytes,"Retained alignment evidence exceeds maxpreparedbytes");
		final AtomicLong retainedBytes=new AtomicLong(evidenceBytes);
		final ExecutorService pool=Executors.newFixedThreadPool(workers);
		final double[] cosine=new double[families], identity=new double[families], blosum=new double[families];
		final MetricResult[][] evidence=metricAware ? new MetricResult[panels.length][families] : null;
		final RawPairStats[][] rawEvidence=rawPairUnion ? new RawPairStats[panels.length][families] : null;
		final WitnessPair[] witnesses=rawPairUnion ? new WitnessPair[families] : null;
		final ByteBuilder timings=new ByteBuilder("phase\tpanel\tfamilies\tqueries\twall_seconds\n");
		RawPairSelection rawPairs=null;
		int liveCount=0; long memberTotal=0, residueTotal=0, wordTotal=0;
		try{
			final long memberStart=System.nanoTime();
			for(int bucket=0; bucket<buckets; bucket++){
				final MemberInput[] rows=source.loadBucket(bucket, maxBucketMembers, maxBucketBytes);
				parallel(pool, workers, rows.length, new Work(){
					@Override public void run(int worker, int first, int end){
						for(int f=first; f<end; f++){
							final MemberInput row=rows[f];
							// Reserve detached primitive arrays before allocating, not after an oversized bucket.
							final long bytes=row.members.length==0 ? 0 : 2L*row.consensus.length+12L*width+8L*400+8L*20;
							require(retainedBytes.addAndGet(bytes)<=maxPreparedBytes, "Prepared primitive arrays exceed maxpreparedbytes");
							prepared[row.rank]=memberWorkers[worker].prepare(row.rank, row.id, row.members, row.consensus);
							rows[f]=null;
						}
					}
				});
				log.println("Preparation read member bucket "+(bucket+1)+" of "+buckets+".");
			}
			source.verifyComplete();
			for(Prepared f:prepared){
				require(f!=null, "Every physical family rank must be prepared before negative calibration");
				if(f.members>0){liveCount++;}
				memberTotal=Math.addExact(memberTotal, f.members); residueTotal=Math.addExact(residueTotal, f.residues);
				wordTotal=Math.addExact(wordTotal, f.validWords);
			}
			require(liveCount>0 && memberTotal==members,
				"Preparation needs live families and the declared complete member population");
			final Prepared[] live=new Prepared[liveCount]; int slot=0;
			for(Prepared f:prepared){if(f.members>0){live[slot++]=f;}}
			final FamilyHeuristicMetadata rankMetadata=metricAware ?
				new FamilyHeuristicMetadata(live, inputPath("lengthmedians"), required("lengthmedianssha80"),
					inputPath("familyhydro"), required("familyhydrosha80"), verificationMode) : null;
			final int[] ordinals=metricAware ? new int[liveCount] : null, hitFloors=metricAware ? new int[liveCount] : null;
			if(metricAware){
				final int floor=integer("ddlreject",1,width);
				for(int f=0;f<liveCount;f++){ordinals[f]=f;hitFloors[f]=floor;}
			}
			timings.append("members\tNA\t").append(families).append("\t0\t")
				.appendSlow((System.nanoTime()-memberStart)*1e-9).nl();
			final NegativePanel[] metricPanels=rawPairUnion ? new NegativePanel[panels.length] : null;
			final PackagedCalibration[] metricCalibrations=rawPairUnion ? new PackagedCalibration[panels.length] : null;
			if(rawPairUnion){for(Prepared f:live){witnesses[f.rank]=new WitnessPair();}}
			for(int p=0; p<panels.length; p++){
				final Panel declared=panels[p]; declared.verify(); final long panelStart=System.nanoTime();
				final boolean raw=p==cosinePanel, firstPanel=p==0;
				final NegativePanel panel=raw ? FamilyNegativeUpperCalibration.loadPreparationPanel(declared.path,
					queryCounts[p], minNegative, maxNegative, bounds, geometry, positionMode, partitions, overlap, blacklist, width) :
					FamilyNegativeUpperCalibration.loadPackagedPanel(declared.path, queryCounts[p], minNegative, maxNegative,
						bounds, geometry, positionMode, partitions, overlap, blacklist, width);
				validateIds(panel, declared.name);
				final PackagedCalibration calibration=metricAware ?
					new PackagedCalibration(live,panel,alignmentMetric,rankMetadata,ordinals,hitFloors) :
					new PackagedCalibration(live, panel);
				if(rawPairUnion){metricPanels[p]=panel; metricCalibrations[p]=calibration;}
				if(rawPairUnion){
					if(raw){
						parallel(pool, workers, live.length, new Work(){
							@Override public void run(int worker, int first, int end){
								for(int f=first; f<end; f++){
									final Prepared family=live[f]; cosine[family.rank]=panel.rawCosineMaximum(family);
								}
							}
						});
					}
				}else{
					final PackagedWorker[] scorers=metricAware ? null : new PackagedWorker[workers];
					final MetricWorker[] metricScorers=metricAware ? new MetricWorker[workers] : null;
					for(int w=0; w<workers; w++){
						if(metricAware){metricScorers[w]=calibration.newMetricWorker(upperTop);}
						else{scorers[w]=calibration.newWorker(upperTop);}
					}
					final MetricResult[] panelEvidence=metricAware ? evidence[p] : null;
					parallel(pool, workers, live.length, new Work(){
						@Override public void run(int worker, int first, int end){
							for(int f=first; f<end; f++){
								final Prepared family=live[f]; final int rank=family.rank;
								if(metricAware){panelEvidence[rank]=metricScorers[worker].scoreFamilyForPreparation(f,null);}
								else{
									final FamilyNegativeUpperCalibration.FamilyResult result=scorers[worker].scoreFamily(f, null);
									identity[rank]=firstPanel ? result.identity : Math.max(identity[rank], result.identity);
									blosum[rank]=firstPanel ? result.blosum : Math.max(blosum[rank], result.blosum);
								}
								if(raw){cosine[rank]=panel.rawCosineMaximum(family);}
							}
						}
					});
				}
				declared.verify();
				timings.append("negative\t").append(declared.name).tab().append(live.length).tab().append(panel.size())
					.tab().appendSlow((System.nanoTime()-panelStart)*1e-9).nl();
				log.println("Preparation scored "+declared.name+": "+panel.size()+" negatives against "+live.length+" families.");
			}
			if(rawPairUnion){
				final long rawStart=System.nanoTime();
				rawPairs=scoreMetricPairs(pool, workers, live, metricPanels, metricCalibrations, evidence, rawEvidence, witnesses);
				timings.append("raw_pairs\tpooled\t").append(live.length).tab().append(queryTotal(metricPanels))
					.tab().appendSlow((System.nanoTime()-rawStart)*1e-9).nl();
				log.println("Preparation raw-pair union selected "+rawPairs.unionPairs+" family/query pairs from "
					+rawPairs.positivePairs+" positive hits; query frontier="+rawPairs.queryPairs
					+", family frontier="+rawPairs.familyPairs+".");
				log.println("Preparation aligned raw-pair union across "+metricPanels.length+" negative panels.");
			}
			if(rankMetadata!=null){rankMetadata.verify();}
		}finally{shutdown(pool);}
		verify(); Files.createDirectories(output.getParent()); Files.createDirectory(output);
		final FamilyTailCascade.LoadedClassifier classifier=writePayloads(output, prepared, liveCount, memberTotal,
			residueTotal, wordTotal, cosine, identity, blosum, evidence, rawEvidence, rawPairs, witnesses, blacklist, threads);
		timings.append("total\tNA\t").append(families).append("\t0\t").appendSlow((System.nanoTime()-started)*1e-9).nl();
		write(output.resolve("preparation.tsv"), timings);
		verify();
		return classifier;
	}

	/** Writes the last evaluated references once, then binds their final fingerprint in the package recipe. */
	void writeFinal(Path out, FamilyTailCascade.LoadedClassifier classifier){
		require(built && outputOptions!=null && !published && classifier!=null,
			"Final preparation output requires its resident classifier and may be published only once");
		final Path references=out.resolve("references.tsv.gz");
		require(references.toAbsolutePath().normalize().toString().equals(outputOptions.get("familyreferences")),
			"Final references must remain in the preparation output directory");
		verify();
		outputBytes=Math.addExact(outputBytes, classifier.writeReferences(references, maxOutputBytes-outputBytes));
		outputOptions.put("familyreferencessha80", recordOutput(references));
		final ByteBuilder recipe=new ByteBuilder();
		if(keyPolicy!=null){
			recipe.append("# keytiebreak=").append(keyPolicy.policy.label).nl()
				.append("# keyencoding=").append(keyPolicy.encodingVersion).nl();
			if(keyPolicy.consensusMemberBonus>0){
				recipe.append("# consensusbonus=").append(keyPolicy.consensusMemberBonus).nl()
					.append("# consensusbonusrule=all-distinct-consensus-bin-key-v1").nl();
			}
			if(keyPolicy.rawInputBinding!=null){
				recipe.append("# keyraw_base64url=").append(java.util.Base64.getUrlEncoder().withoutPadding().encodeToString(
					keyPolicy.rawInputBinding.getBytes(StandardCharsets.UTF_8))).nl();
			}
		}
		if(verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL){
			recipe.append("# Structural resident result; serialized reload has not been validated.\nverification=structural\n");
		}
		if(identityUpperMultiplier!=FamilyGateCodec.DEFAULT_IDENTITY_UPPER_MULTIPLIER ||
				identityUpperAdditive!=FamilyGateCodec.DEFAULT_IDENTITY_UPPER_ADDITIVE){
			recipe.append("# identityuppermultiplier=").appendSlow(identityUpperMultiplier).nl()
				.append("# identityupperadditive=").appendSlow(identityUpperAdditive).nl();
		}
		for(String key:new TreeSet<String>(outputOptions.keySet())){
			recipe.append(key).append('=').append(outputOptions.get(key)).nl();
		}
		final Path config=out.resolve("classifier.config");
		write(config, recipe); recordOutput(config);
		write(out.resolve("preparation_inputs.tsv"), inputTable());
		verify(); published=true;
	}

	/** Rebuilds only representative keys, preserving the supplied package's reference payload
	 * and numeric gates. Every default row must reproduce its baseline before publication.
	 * Structural verification is explicit: no input or output file digests are computed.
	 */
	void writeKeyRarityProfiles(FamilyScoringPackage baseline,Path output,DDLKeyRarity.Policy policy,
			int threads,PrintStream log)throws Exception{
		require(policy!=null,"An explicit rarity policy is required");
		writeProfileChanges(baseline,output,policy,null,0,threads,log);
	}

	/** Released fixed-family first look: families not rebuilt; lookup-only.
	 * results/family_sketch_phase_a_20261116 lost ALL/40-99 raw objective in
	 * every arm: W1 -415/-242, W10 -267/-196, W1000/W10000 -344/-218,
	 * fixed-quadrant -1271/-1890. None was adopted.
	 * Reconstruct every accepted C/W0 row before selecting the target strategy.
	 * Numeric gates remain unchanged here; each resulting arm MUST subsequently fit its own negatives.
	 */
	void writeSketchStrategyProfiles(FamilyScoringPackage baseline,Path output,ProteinPositionMode targetMode,
			long consensusBonus,int threads,PrintStream log)throws Exception{
		require(targetMode==ProteinPositionMode.FIXED || targetMode==ProteinPositionMode.FIXED_QUADRANT,
			"Sketch strategy requires fixed or fixed-quadrant geometry");
		require(consensusBonus==0 || consensusBonus==1 || consensusBonus==10 || consensusBonus==1000 || consensusBonus==10000,
			"Sketch strategy bonus must be one of the six released arms");
		require(targetMode==ProteinPositionMode.FIXED || consensusBonus==0,
			"Do not combine quadrant placement with a consensus bonus");
		require(positionMode==ProteinPositionMode.FIXED && partitions==4 && overlap==2 && width==512 && geometry.k()==6 &&
			geometry.name().startsWith("dimer6+1,3,5+") &&
			geometry.name().substring("dimer6+1,3,5+".length()).split("/",-1).length==120,
			"Released sketch strategies require unchanged d120-K6/512/P4/O2 accepted fixed baseline");
		writeProfileChanges(baseline,output,null,targetMode,consensusBonus,threads,log);
	}

	/** One member-bucket traversal/export path for legacy rarity and explicit sketch-strategy changes. */
	private void writeProfileChanges(FamilyScoringPackage baseline,Path output,DDLKeyRarity.Policy policy,
			ProteinPositionMode targetMode,long consensusBonus,int threads,PrintStream log)throws Exception{

		require(!built && baseline!=null && output!=null && (policy!=null || targetMode!=null) &&
			threads>0 && !Files.exists(output),
			"Profile rebuilding requires a fresh output, baseline, explicit policy and positive workers");
		require(verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL &&
			baseline.verificationMode()==verificationMode,"Profile-only comparison requires explicit structural verification");
		built=true; verify(); baseline.verify();
		final FamilyTailCascade.LoadedClassifier parent=baseline.load(threads);
		final FamilyTailCascade.SearchConfig config=parent.config();
		final boolean strategy=targetMode!=null;
		final ProteinPositionMode outputMode=strategy ? targetMode : positionMode;
		final String policyLabel=strategy ? DDLKeyTieBreak.Policy.CONSENSUS_DIVERSITY.label : policy.label;
		final HashMap<String,String> selected=baseline.resolvedOptions();
		require(selected.containsKey("ddlprofiles"),"Rarity rebuilding requires an explicit stored-profile baseline");
		require(parent.familyCount()==families && memberRange.same(MemberLengthRange.from(selected)) &&
			config.width==width && config.seed==0 && config.geometry.equals(geometry) && config.positionMode==positionMode &&
			config.partitions==partitions && config.overlap==overlap && config.alignmentMetric.same(alignmentMetric),
			"Baseline and member recipe must name the same complete live roster, geometry, member range and alignment contract");
		for(String key:new String[]{"blacklist","blacklistsha80","blacklistk","blacklistmode",
				"blacklistwords","blacklistprovenance"}){
			require(java.util.Objects.equals(options.get(key),selected.get(key)),
				"Baseline/member blacklist binding differs: "+key);
		}
		final int[] ordinalByRank=new int[families]; Arrays.fill(ordinalByRank,-1);
		final DDLFamilySearch.StructuralInputs rowInputs=new DDLFamilySearch.StructuralInputs();
		final char[][] baselineRows=loadKeyRarityBaseline(parent,
			java.nio.file.Paths.get(selected.get("ddlprofiles")),rowInputs);
		final char[][] resultRows=new char[families][];
		long baselineMembers=0;
		for(int ordinal=0; ordinal<families; ordinal++){
			final int rank=parent.familyRankAt(ordinal);
			require(rank>=0 && rank<families && ordinalByRank[rank]<0,
				"Baseline requires a complete unique physical/live rank mapping");
			ordinalByRank[rank]=ordinal;
			baselineMembers=Math.addExact(baselineMembers,parent.familyMembersAt(ordinal));
		}
		require(baselineMembers==members,"Baseline and preparation must retain the same complete assigned-member population");
		final DDLKeyRarity rarity=strategy ? null : DDLKeyRarity.fromProfiles(width,baselineRows);
		final MemberInputs source=new MemberInputs(inputPath("families"),required("familiessha80"),
			inputPath("chunks"),required("chunkssha80"),inputPath("consensus"),required("consensussha80"),
			families,members,buckets,batches,maxReference,memberRange,verificationMode);
		final int workers=Math.min(threads,families);
		final ProfileRebuildWorker[] builders=new ProfileRebuildWorker[workers];
		for(int i=0;i<workers;i++){
			builders[i]=strategy ? new SketchStrategyWorker(config,targetMode,consensusBonus) :
				new KeyRarityWorker(config,policy==DDLKeyRarity.Policy.RARITY ? rarity : null);
		}
		final ExecutorService pool=Executors.newFixedThreadPool(workers);
		Files.createDirectory(output); final Path profileDir=output.resolve("profiles"); Files.createDirectory(profileDir);
		final ByteStreamWriter writer=new ByteStreamWriter(profileDir.resolve("profiles.tsv.gz").toString(),
			false,false,false);
		writer.start();
		final boolean[] seen=new boolean[families];
		long totalMembers=0, totalResidues=0, totalWords=0, changedKeys=0, changedFamilies=0;
		try{
			emit(writer,new ByteBuilder(DDLFamilyProfiles.HEADER16).nl());
			for(int bucket=0; bucket<buckets; bucket++){
				final MemberInput[] rows=source.loadBucket(bucket,maxBucketMembers,maxBucketBytes);
				final KeyRarityRow[] rebuilt=new KeyRarityRow[rows.length];
				for(MemberInput row:rows){
					require(row.rank>=0 && row.rank<families && !seen[row.rank],"Each physical family must be rebuilt exactly once");
					seen[row.rank]=true;
				}
				parallel(pool,workers,rows.length,new Work(){
					@Override public void run(int worker,int first,int end){
						for(int f=first; f<end; f++){
							final MemberInput row=rows[f]; final int ordinal=ordinalByRank[row.rank];
							rebuilt[f]=builders[worker].build(row,config.families[ordinal],baselineRows[ordinal]);
							resultRows[ordinal]=rebuilt[f].values; rows[f]=null;
						}
					}
				});
				for(KeyRarityRow row:rebuilt){
					emit(writer,row.line);
					totalMembers=Math.addExact(totalMembers,row.members); totalResidues=Math.addExact(totalResidues,row.residues);
					totalWords=Math.addExact(totalWords,row.words); changedKeys=Math.addExact(changedKeys,row.changed);
					if(row.changed>0){changedFamilies++;}
				}
				log.println("Key tie-break rebuilt member bucket "+(bucket+1)+" of "+buckets+".");
			}
		}finally{shutdown(pool); require(!writer.poisonAndWait(),"Profile-only output writer failed; no success marker");}
		source.verifyComplete(); verify(); baseline.verify(); rowInputs.verify();
		for(boolean complete:seen){require(complete,"Profile-only output must cover the complete baseline roster");}
		require(totalMembers==members &&
				((strategy ? targetMode!=positionMode || consensusBonus>0 :
				policy!=DDLKeyRarity.Policy.SMALLEST) || changedKeys==0),
			"Profile rebuilding must preserve all members and reproduce baseline keys under smallest");
		final DDLPentamerBlacklist blacklist=config.blacklist;
		final boolean legacyK5=geometry.isExact(5);
		final ByteBuilder summary=new ByteBuilder(
			"namespace\tmanifest_sha80\tk\tseed\texponent\tmanifest_families\tbuilt_families"+
			"\tmembers\tresidues\tvalid_words\tsketch")
			.append(DDLFamilyProfiles.positionMetadataColumns(outputMode));
		if(!legacyK5){summary.append("\tgeometry");}
		if(blacklist!=null){
			summary.append("\tblacklist_sha80\tblacklist_words");
			if(FamilyBlacklistBinding.hasPolicy(options)){summary.tab().append(blacklist.metadataColumns());}
		}
		if(strategy){
			summary.append("\tkeytiebreak\tkeyencoding");
			if(consensusBonus>0){summary.append("\tconsensusbonus\tconsensusbonusrule");}
		}
		summary.nl().append("tail\t").append(required("familiessha80")).tab().append(geometry.k()).append("\t0\t6\t")
			.append(families).tab().append(families).tab().append(totalMembers).tab().append(totalResidues)
			.tab().append(totalWords).tab().append(outputMode.sketch)
			.append(DDLFamilyProfiles.positionMetadataValues(outputMode,partitions,overlap,geometry));
		if(!legacyK5){summary.tab().append(geometry.name());}
		if(blacklist!=null){
			summary.tab().append(blacklist.sha80()).tab().append(blacklist.size());
			if(FamilyBlacklistBinding.hasPolicy(options)){summary.tab().append(blacklist.metadataValues());}
		}
		if(strategy){
			summary.tab().append(DDLKeyTieBreak.Policy.CONSENSUS_DIVERSITY.label)
				.tab().append(DDLKeyTieBreak.encodingVersion(DDLKeyTieBreak.Policy.CONSENSUS_DIVERSITY));
			if(consensusBonus>0){summary.tab().append(consensusBonus).tab().append("all-distinct-consensus-bin-key-v1");}
		}
		write(profileDir.resolve("summary.tsv"),summary.nl());
		write(output.resolve("ddlprofiles.tsv"),new ByteBuilder("profiles/profiles.tsv.gz\tunmeasured\tunmeasured\n"));
		selected.put("ddlprofiles",output.resolve("ddlprofiles.tsv").toAbsolutePath().normalize().toString());
		selected.put("ddlprofilessha80","unmeasured");
		if(strategy){selected.put("positionmode",outputMode.flag);}
		final ByteBuilder recipe=new ByteBuilder("# keytiebreak=").append(policyLabel)
			.append("\n# baseline=").append(baseline.path.toString()).append("\n# member_recipe=").append(path.toString())
			.append("\nverification=structural\n");
		if(strategy){
			recipe.append("# families not rebuilt; lookup-only\n# keyencoding=")
				.append(DDLKeyTieBreak.encodingVersion(DDLKeyTieBreak.Policy.CONSENSUS_DIVERSITY)).nl();
			if(consensusBonus>0){recipe.append("# consensusbonus=").append(consensusBonus)
				.append("\n# consensusbonusrule=all-distinct-consensus-bin-key-v1\n");}
			recipe.append("# each-arm negative calibration required before evaluation\n");
		}
		for(String key:new TreeSet<String>(selected.keySet())){recipe.append(key).append('=').append(selected.get(key)).nl();}
		write(output.resolve("classifier.config"),recipe);
		final FamilyScoringPackage exported=new FamilyScoringPackage(output,verificationMode);
		final FamilyTailCascade.LoadedClassifier replay=exported.load(threads);
		require(replay.matchesProfiles(resultRows),"Serialized replacement profiles must match the computed keys");
		for(int ordinal=0; ordinal<families; ordinal++){
			require(replay.familyRankAt(ordinal)==parent.familyRankAt(ordinal) &&
				replay.familyIdAt(ordinal).equals(parent.familyIdAt(ordinal)) &&
				replay.familyMembersAt(ordinal)==parent.familyMembersAt(ordinal) &&
				replay.familyAlignmentEnabledAt(ordinal)==parent.familyAlignmentEnabledAt(ordinal) &&
				replay.familyIdentityGateAt(ordinal)==parent.familyIdentityGateAt(ordinal) &&
				replay.familyBlosumGateAt(ordinal)==parent.familyBlosumGateAt(ordinal) &&
				replay.familyEffectiveIdentityLowerAt(ordinal)==parent.familyEffectiveIdentityLowerAt(ordinal),
				"Profile-only reload must preserve exact baseline identities, counts and numeric alignment gates");
		}
		exported.verify(); baseline.verify(); rowInputs.verify(); verify();
		if(strategy){
			write(output.resolve("sketch_strategy.tsv"),new ByteBuilder(
				"scope\tposition_mode\tkeypolicy\tconsensusbonus\tfamilies\tmembers"+
				"\tchanged_families\tchanged_keys\n")
				.append("families not rebuilt; lookup-only\t").append(outputMode.flag).tab()
				.append(policyLabel).tab().append(consensusBonus)
				.tab().append(families).tab().append(totalMembers).tab().append(changedFamilies).tab().append(changedKeys).nl());
		}else{
			write(output.resolve("key_tiebreak.tsv"),new ByteBuilder(
				"policy\tfamilies\tmembers\tbaseline_postings\tchanged_families\tchanged_keys\n")
				.append(policy.label).tab().append(families).tab().append(totalMembers).tab().append(rarity.postings())
				.tab().append(changedFamilies).tab().append(changedKeys).nl());
		}
		write(output.resolve("preparation_inputs.tsv"),inputTable());
		write(output.resolve("PASS"),new ByteBuilder(strategy ? "FAMILY_SKETCH_STRATEGY_PROFILES_PASS\n" :
			"FAMILY_KEY_RARITY_PROFILES_PASS\n"));
	}

	/** Rebuilds the stored search rows as plain DDLs while retaining all reference/gate rows. */
	void writePlainProfiles(FamilyScoringPackage baseline,Path output,int threads,PrintStream log) throws Exception{
		require(!built && baseline!=null && output!=null && threads>0 && !Files.exists(output),
			"Plain profile rebuilding requires a fresh output, baseline and positive workers");
		require(verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL &&
			baseline.verificationMode()==verificationMode,"Plain profile rebuilding requires explicit structural verification");
		built=true; verify(); baseline.verify();
		final FamilyTailCascade.LoadedClassifier parent=baseline.load(threads);
		final FamilyTailCascade.SearchConfig config=parent.config();
		final HashMap<String,String> selected=baseline.resolvedOptions();
		require(selected.containsKey("ddlprofiles"),"Plain rebuilding requires an explicit stored-profile baseline");
		require(parent.familyCount()==families && memberRange.same(MemberLengthRange.from(selected)) &&
			config.width==width && config.seed==0 && config.geometry.equals(geometry) &&
			config.alignmentMetric.same(alignmentMetric),
			"Baseline and member recipe must name the same complete live roster, geometry, member range and alignment contract");
		for(String key:new String[]{"blacklist","blacklistsha80","blacklistk","blacklistmode",
				"blacklistwords","blacklistprovenance"}){
			require(java.util.Objects.equals(options.get(key),selected.get(key)),
				"Baseline/member blacklist binding differs: "+key);
		}
		final int[] ordinalByRank=new int[families]; Arrays.fill(ordinalByRank,-1);
		long parentMembers=0;
		for(int ordinal=0; ordinal<families; ordinal++){
			final int rank=parent.familyRankAt(ordinal);
			require(rank>=0 && rank<families && ordinalByRank[rank]<0,
				"Baseline requires a complete unique physical/live rank mapping");
			ordinalByRank[rank]=ordinal;
			parentMembers=Math.addExact(parentMembers,parent.familyMembersAt(ordinal));
		}
		require(parentMembers==members,"Baseline and preparation must retain the same complete assigned-member population");
		final MemberInputs source=new MemberInputs(inputPath("families"),required("familiessha80"),
			inputPath("chunks"),required("chunkssha80"),inputPath("consensus"),required("consensussha80"),
			families,members,buckets,batches,maxReference,memberRange,verificationMode);
		final int workers=Math.min(threads,families);
		final PlainProfileWorker[] builders=new PlainProfileWorker[workers];
		for(int i=0; i<workers; i++){builders[i]=new PlainProfileWorker(config);}
		final ExecutorService pool=Executors.newFixedThreadPool(workers);
		Files.createDirectory(output); final Path profileDir=output.resolve("profiles"); Files.createDirectory(profileDir);
		final ByteStreamWriter writer=new ByteStreamWriter(profileDir.resolve("profiles.tsv.gz").toString(),
			false,false,false);
		writer.start();
		final char[][] resultRows=new char[families][];
		final boolean[] seen=new boolean[families];
		long totalMembers=0, totalResidues=0, totalWords=0;
		try{
			emit(writer,new ByteBuilder(DDLFamilyProfiles.HEADER16).nl());
			for(int bucket=0; bucket<buckets; bucket++){
				final MemberInput[] rows=source.loadBucket(bucket,maxBucketMembers,maxBucketBytes);
				final KeyRarityRow[] rebuilt=new KeyRarityRow[rows.length];
				for(MemberInput row:rows){
					require(row.rank>=0 && row.rank<families && !seen[row.rank],"Each physical family must be rebuilt exactly once");
					seen[row.rank]=true;
				}
				parallel(pool,workers,rows.length,new Work(){
					@Override public void run(int worker,int first,int end){
						for(int f=first; f<end; f++){
							final MemberInput row=rows[f]; final int ordinal=ordinalByRank[row.rank];
							rebuilt[f]=builders[worker].build(row,config.families[ordinal]);
							resultRows[ordinal]=rebuilt[f].values; rows[f]=null;
						}
					}
				});
				for(KeyRarityRow row:rebuilt){
					emit(writer,row.line);
					totalMembers=Math.addExact(totalMembers,row.members);
					totalResidues=Math.addExact(totalResidues,row.residues);
					totalWords=Math.addExact(totalWords,row.words);
				}
				log.println("Plain profile rebuilt member bucket "+(bucket+1)+" of "+buckets+".");
			}
		}finally{shutdown(pool); require(!writer.poisonAndWait(),"Plain profile writer failed; no success marker");}
		source.verifyComplete(); verify(); baseline.verify();
		for(boolean complete:seen){require(complete,"Plain profile output must cover the complete baseline roster");}
		require(totalMembers==members,"Plain profile rebuilding must preserve all assigned members");
		final DDLPentamerBlacklist blacklist=config.blacklist;
		final boolean legacyK5=geometry.isExact(5);
		final ByteBuilder summary=new ByteBuilder(
			"namespace\tmanifest_sha80\tk\tseed\texponent\tmanifest_families\tbuilt_families"+
			"\tmembers\tresidues\tvalid_words\tsketch");
		if(!legacyK5){summary.append("\tgeometry");}
		if(blacklist!=null){
			summary.append("\tblacklist_sha80\tblacklist_words");
			if(FamilyBlacklistBinding.hasPolicy(options)){summary.tab().append(blacklist.metadataColumns());}
		}
		summary.nl().append("tail\t").append(required("familiessha80")).tab().append(geometry.k()).append("\t0\t6\t")
			.append(families).tab().append(families).tab().append(totalMembers).tab()
			.append(totalResidues).tab().append(totalWords)
			.append("\t").append(ProteinPositionMode.PLAIN.sketch);
		if(!legacyK5){summary.tab().append(geometry.name());}
		if(blacklist!=null){
			summary.tab().append(blacklist.sha80()).tab().append(blacklist.size());
			if(FamilyBlacklistBinding.hasPolicy(options)){summary.tab().append(blacklist.metadataValues());}
		}
		write(profileDir.resolve("summary.tsv"),summary.nl());
		write(output.resolve("ddlprofiles.tsv"),new ByteBuilder("profiles/profiles.tsv.gz\tunmeasured\tunmeasured\n"));
		selected.put("ddlprofiles",output.resolve("ddlprofiles.tsv").toAbsolutePath().normalize().toString());
		selected.put("ddlprofilessha80","unmeasured");
		selected.put("positional","false");
		selected.put("positionmode",ProteinPositionMode.PLAIN.flag);
		selected.remove("partitions");
		selected.remove("overlap");
		final ByteBuilder recipe=new ByteBuilder("# profilemode=plain\n# baseline=").append(baseline.path.toString())
			.append("\n# member_recipe=").append(path.toString()).append("\nverification=structural\n");
		for(String key:new TreeSet<String>(selected.keySet())){recipe.append(key).append('=').append(selected.get(key)).nl();}
		write(output.resolve("classifier.config"),recipe);
		final FamilyScoringPackage exported=new FamilyScoringPackage(output,verificationMode);
		final FamilyTailCascade.LoadedClassifier replay=exported.load(threads);
		require(replay.matchesProfiles(resultRows),"Serialized plain profiles must match the computed keys");
		for(int ordinal=0; ordinal<families; ordinal++){
			require(replay.familyRankAt(ordinal)==parent.familyRankAt(ordinal) &&
				replay.familyIdAt(ordinal).equals(parent.familyIdAt(ordinal)) &&
				replay.familyMembersAt(ordinal)==parent.familyMembersAt(ordinal) &&
				replay.familyAlignmentEnabledAt(ordinal)==parent.familyAlignmentEnabledAt(ordinal) &&
				replay.familyIdentityGateAt(ordinal)==parent.familyIdentityGateAt(ordinal) &&
				replay.familyBlosumGateAt(ordinal)==parent.familyBlosumGateAt(ordinal) &&
				replay.familyEffectiveIdentityLowerAt(ordinal)==parent.familyEffectiveIdentityLowerAt(ordinal),
				"Plain profile reload must preserve exact baseline identities, counts and numeric alignment gates");
		}
		exported.verify(); baseline.verify(); verify();
		write(output.resolve("plain_profiles.tsv"),new ByteBuilder("mode\tfamilies\tmembers\nplain\t")
			.append(families).tab().append(totalMembers).nl());
		write(output.resolve("preparation_inputs.tsv"),inputTable());
		write(output.resolve("PASS"),new ByteBuilder("FAMILY_PLAIN_PROFILES_PASS\n"));
	}

	/** Four-arm profile-only comparison. Count every original input tail once, then
	 * share that immutable resident table while traversing assigned members once.
	 */
	void writeKeyTieBreakComparison(FamilyScoringPackage baseline,Path output,Path rawManifest,
			long expectedRawTails,int threads,PrintStream log) throws Exception{
		require(!built && baseline!=null && output!=null && rawManifest!=null && threads>0 && !Files.exists(output),
			"Packed comparison requires fresh output, full raw manifest and positive workers");
		require(verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL &&
			baseline.verificationMode()==verificationMode,
			"Packed comparison requires explicit structural verification");
		built=true;verify();baseline.verify();
		final FamilyTailCascade.LoadedClassifier parent=baseline.load(threads);
		final FamilyTailCascade.SearchConfig config=parent.config();
		final HashMap<String,String> selected=baseline.resolvedOptions();
		require(selected.containsKey("ddlprofiles") && parent.familyCount()==families &&
			memberRange.same(MemberLengthRange.from(selected)) && config.width==width && config.seed==0 &&
			config.geometry.equals(geometry) && config.positionMode==positionMode && config.partitions==partitions &&
			config.overlap==overlap && config.alignmentMetric.same(alignmentMetric),
			"Packed baseline/member recipe geometry, roster, range and alignment contract must agree");
		for(String key:new String[]{"blacklist","blacklistsha80","blacklistk",
				"blacklistmode","blacklistwords","blacklistprovenance"}){
			require(java.util.Objects.equals(options.get(key),selected.get(key)),
				"Packed baseline/member blacklist differs: "+key);
		}
		final DDLFamilySearch.StructuralInputs tracked=new DDLFamilySearch.StructuralInputs();
		final char[][] baselineRows=loadKeyRarityBaseline(parent,
			java.nio.file.Paths.get(selected.get("ddlprofiles")),tracked);
		final int[] ordinalByRank=new int[families];Arrays.fill(ordinalByRank,-1);long parentMembers=0;
		for(int ordinal=0;ordinal<families;ordinal++){
			final int rank=parent.familyRankAt(ordinal);
			require(rank>=0 && rank<families && ordinalByRank[rank]<0,"Packed comparison needs a complete unique live roster");
			ordinalByRank[rank]=ordinal;parentMembers=Math.addExact(parentMembers,parent.familyMembersAt(ordinal));
		}
		require(parentMembers==members,"Packed comparison must retain the entire assigned-member population");
		// Parent geometry/filter are the very objects subsequently used on assigned members.
		final long retainedPayload=Math.addExact(Math.multiplyExact((long)width,65536L*6),
			Math.multiplyExact(Math.multiplyExact((long)families,width),2L*(DDLKeyTieBreak.Policy.values().length+1)));
		require(retainedPayload<=maxPreparedBytes,"Packed retained count/profile payload exceeds maxpreparedbytes");
		final DDLKeyTieBreak.RawTailCounts raw=countRawTails(rawManifest,expectedRawTails,config,tracked,log);
		verify();baseline.verify();tracked.verify();
		final DDLKeyTieBreak.Policy[] policies=DDLKeyTieBreak.Policy.values();
		final char[][][] resultRows=new char[policies.length][families][];
		final long[] changedFamilies=new long[policies.length],changedKeys=new long[policies.length];
		final Path[] outputs=new Path[policies.length];final ByteStreamWriter[] writers=new ByteStreamWriter[policies.length];
		final MemberInputs source=new MemberInputs(inputPath("families"),required("familiessha80"),inputPath("chunks"),
			required("chunkssha80"),inputPath("consensus"),required("consensussha80"),families,members,buckets,batches,
			maxReference,memberRange,verificationMode);
		final int workers=Math.min(threads,families);final PackedKeyWorker[] builders=new PackedKeyWorker[workers];
		for(int i=0;i<workers;i++){builders[i]=new PackedKeyWorker(config,raw);}
		Files.createDirectory(output);
		final ExecutorService pool=Executors.newFixedThreadPool(workers);
		final boolean[] seen=new boolean[families];long totalMembers=0,totalResidues=0,totalWords=0;
		try{
			for(int arm=0;arm<policies.length;arm++){
				outputs[arm]=output.resolve(policies[arm].label);
				Files.createDirectory(outputs[arm]);
				Files.createDirectory(outputs[arm].resolve("profiles"));
				writers[arm]=new ByteStreamWriter(outputs[arm].resolve("profiles/profiles.tsv.gz").toString(),false,false,false);
				writers[arm].start();emit(writers[arm],new ByteBuilder(DDLFamilyProfiles.HEADER16).nl());
			}
			for(int bucket=0;bucket<buckets;bucket++){
				final MemberInput[] rows=source.loadBucket(bucket,maxBucketMembers,maxBucketBytes);
				final KeyRarityRow[][] rebuilt=new KeyRarityRow[rows.length][];
				for(MemberInput row:rows){
					require(row.rank>=0 && row.rank<families && !seen[row.rank],
						"Each family must be rebuilt exactly once");
					seen[row.rank]=true;
				}
				parallel(pool,workers,rows.length,new Work(){
					@Override public void run(int worker,int first,int end){
						for(int f=first;f<end;f++){
							final MemberInput row=rows[f];final int ordinal=ordinalByRank[row.rank];
							rebuilt[f]=builders[worker].build(row,config.families[ordinal],baselineRows[ordinal]);
							for(int arm=0;arm<policies.length;arm++){resultRows[arm][ordinal]=rebuilt[f][arm].values;}
							rows[f]=null;
						}
					}
				});
				for(KeyRarityRow[] row:rebuilt){
					totalMembers=Math.addExact(totalMembers,row[0].members);totalResidues=Math.addExact(totalResidues,row[0].residues);
					totalWords=Math.addExact(totalWords,row[0].words);
					for(int arm=0;arm<policies.length;arm++){
						emit(writers[arm],row[arm].line);
						changedKeys[arm]+=row[arm].changed;
						if(row[arm].changed>0){changedFamilies[arm]++;}
					}
				}
				log.println("Packed four-arm member bucket "+(bucket+1)+" of "+buckets+" complete.");
			}
		}finally{
			try{shutdown(pool);}finally{
				boolean error=false;
				for(ByteStreamWriter writer:writers){if(writer!=null){error|=writer.poisonAndWait();}}
				require(!error,"Packed profile writer failed; no completed comparison");
			}
		}
		source.verifyComplete();verify();baseline.verify();tracked.verify();
		for(boolean present:seen){require(present,"Packed comparison omitted a live family");}
		require(totalMembers==members && changedKeys[DDLKeyTieBreak.Policy.SMALLEST.ordinal()]==0,
			"Packed default must reconstruct every baseline key and preserve all members");
		for(int arm=0;arm<policies.length;arm++){
			final DDLKeyTieBreak.Diagnostics diagnostics=new DDLKeyTieBreak.Diagnostics();
			for(PackedKeyWorker builder:builders){diagnostics.add(builder.diagnostics[arm]);}
			publishPackedArm(baseline,parent,outputs[arm],policies[arm],resultRows[arm],
				raw,totalMembers,totalResidues,totalWords,
				changedFamilies[arm],changedKeys[arm],threads,diagnostics);
			resultRows[arm]=null;
		}
		final ByteBuilder receipt=new ByteBuilder("key\tvalue\nraw_manifest\t")
			.append(rawManifest.toAbsolutePath().normalize().toString())
			.append("\nraw_tails\t").append(raw.tails()).append("\nheld_registers\t").append(raw.postings())
			.append("\nmember_recipe\t").append(path.toString()).append("\nbaseline\t").append(baseline.path.toString())
			.append("\ngeometry\t").append(config.geometry.name()).append("\nwidth\t").append(config.width)
			.append("\nsketchseed\t").append(config.seed).append("\nposition_mode\t").append(config.positionMode.flag)
			.append("\npartitions\t").append(config.partitions).append("\noverlap\t").append(config.overlap)
			.append("\nmember_min\t").append(memberRange.min).append("\nmember_max\t").append(memberRange.max)
			.append("\nhash_recipe\t").append(DDLKeyTieBreak.HASH_RECIPE)
			.append("\npivot\tunsigned40-votes-24-signal-descending-stable-canonical-uint16"
				+"\nraw_count_storage\tresident-shared-once\n");
		for(String key:new String[]{"aabits","blacklist","blacklistsha80",
				"blacklistmode","blacklistk","blacklistwords","blacklistprovenance"}){
			receipt.append(key).tab().append(selected.containsKey(key) ? selected.get(key) : "none").nl();
		}
		write(output.resolve("raw_count_receipt.tsv"),receipt);
		Files.copy(rawManifest,output.resolve("raw_inputs.tsv"));
		verify();baseline.verify();tracked.verify();
		write(output.resolve("PASS"),new ByteBuilder("FAMILY_KEY_TIEBREAK_COMPARISON_PASS\n"));
	}

	private DDLKeyTieBreak.RawTailCounts countRawTails(Path manifest,long expected,FamilyTailCascade.SearchConfig config,
			DDLFamilySearch.StructuralInputs tracked,PrintStream log) throws Exception{
		manifest=tracked.track(manifest);
		final ArrayList<Path> inputs=new ArrayList<Path>();
		final HashSet<Path> unique=new HashSet<Path>();
		final ByteFile reader=ByteFile.makeByteFile(manifest.toString(),false);
		try{
			final byte[] header=reader.nextLine();
			require(header!=null &&
				new String(header,StandardCharsets.US_ASCII).equals("order\tpath\tbytes\tmtime_seconds"),
				"Exact raw input manifest schema required");
			for(byte[] line=reader.nextLine();line!=null;line=reader.nextLine()){
				final String[] fields=new String(line,StandardCharsets.UTF_8).split("\t",-1);
				require(fields.length==4 && Integer.parseInt(fields[0])==inputs.size(),"Raw inputs require dense explicit order");
				final Path input=java.nio.file.Paths.get(fields[1]);
				require(input.isAbsolute() && input.normalize().equals(input) && unique.add(input),
					"Raw inputs must be distinct absolute normalized paths");
				tracked.track(input);
				require(Files.size(input)==Long.parseLong(fields[2]) &&
					Files.getLastModifiedTime(input).toMillis()/1000==Long.parseLong(fields[3]),
					"Raw input metadata differs from its declared population manifest: "+input);
				inputs.add(input);
			}
		}finally{require(!reader.close(),"Raw input manifest read failed");}
		require(!inputs.isEmpty(),"Complete raw input stream is required");
		final DDLKeyTieBreak.RawTailCounter counter=new DDLKeyTieBreak.RawTailCounter(config.width,expected);
		final char[] sketch=new char[config.width];long processed=0;
		for(Path input:inputs){
			final FileFormat format=FileFormat.testInput(input.toString(),FileFormat.FASTA,null,true,true);
			require(format.fasta() || format.fastq(),"Raw protein inputs must be FASTA/FASTQ");
			final Streamer stream=StreamerFactory.makeStreamer(format,null,true,-1,false,true,1);stream.start();
			try{
				for(ListNum<Read> batch=stream.nextList();batch!=null;batch=stream.nextList()){
					for(Read read:batch.list){
						require(read!=null && read.bases!=null && read.mate==null && memberRange.contains(read.bases.length),
							"Full raw population requires unpaired tails within the bound member range; no silent filtering");
						if(config.positionMode.positional()){
							DDLFamilyProfiles.sketchPositional(read.bases,config.geometry,config.seed,
								config.partitions,config.overlap,config.positionMode,sketch,config.blacklist);
						}else{DDLFamilyProfiles.sketch(read.bases,config.geometry,config.seed,sketch,config.blacklist);}
						counter.add(sketch);processed++;
					}
				}
			}finally{require(!ReadWrite.closeStream(stream),"Raw-tail stream failed: "+input);}
			tracked.verify();log.println("Raw-tail DDL count stream complete: "+input+"; cumulative tails="+processed);
		}
		tracked.verify();return counter.freeze();
	}

	private void publishPackedArm(FamilyScoringPackage baseline,FamilyTailCascade.LoadedClassifier parent,
			Path output, DDLKeyTieBreak.Policy policy, char[][] rows, DDLKeyTieBreak.RawTailCounts raw,
			long totalMembers, long residues, long words, long changedFamilies, long changedKeys,
			int threads, DDLKeyTieBreak.Diagnostics diagnostics) throws Exception{
		final FamilyTailCascade.SearchConfig config=parent.config();final DDLPentamerBlacklist blacklist=config.blacklist;
		final boolean legacyK5=geometry.isExact(5);
		final ByteBuilder summary=new ByteBuilder("namespace\tmanifest_sha80\tk\tseed\texponent"
			+"\tmanifest_families\tbuilt_families\tmembers\tresidues\tvalid_words\tsketch")
			.append(DDLFamilyProfiles.positionMetadataColumns(positionMode));
		if(!legacyK5){summary.append("\tgeometry");}
		if(blacklist!=null){
			summary.append("\tblacklist_sha80\tblacklist_words");
			if(FamilyBlacklistBinding.hasPolicy(options)){summary.tab().append(blacklist.metadataColumns());}
		}
		summary.nl().append("tail\t").append(required("familiessha80")).tab().append(geometry.k()).append("\t0\t6\t")
			.append(families).tab().append(families).tab().append(totalMembers).tab().append(residues).tab().append(words)
			.tab().append(positionMode.sketch)
			.append(DDLFamilyProfiles.positionMetadataValues(positionMode,partitions,overlap,geometry));
		if(!legacyK5){summary.tab().append(geometry.name());}
		if(blacklist!=null){
			summary.tab().append(blacklist.sha80()).tab().append(blacklist.size());
			if(FamilyBlacklistBinding.hasPolicy(options)){summary.tab().append(blacklist.metadataValues());}
		}
		write(output.resolve("profiles/summary.tsv"),summary.nl());
		write(output.resolve("ddlprofiles.tsv"),new ByteBuilder("profiles/profiles.tsv.gz\tunmeasured\tunmeasured\n"));
		final HashMap<String,String> selected=baseline.resolvedOptions();
		selected.put("ddlprofiles",output.resolve("ddlprofiles.tsv").toAbsolutePath().normalize().toString());
		selected.put("ddlprofilessha80","unmeasured");
		final ByteBuilder recipe=new ByteBuilder("# packed_keytiebreak=").append(policy.label)
			.append("\nverification=structural\n");
		for(String key:new TreeSet<String>(selected.keySet())){recipe.append(key).append('=').append(selected.get(key)).nl();}
		write(output.resolve("classifier.config"),recipe);
		final FamilyScoringPackage exported=new FamilyScoringPackage(output,verificationMode);
		final FamilyTailCascade.LoadedClassifier replay=exported.load(threads);
		require(replay.matchesProfiles(rows),"Packed serialized profile keys differ from computed values");
		for(int ordinal=0;ordinal<families;ordinal++){
			require(replay.familyRankAt(ordinal)==parent.familyRankAt(ordinal) &&
				replay.familyIdAt(ordinal).equals(parent.familyIdAt(ordinal)) &&
				replay.familyMembersAt(ordinal)==parent.familyMembersAt(ordinal) &&
				replay.familyAlignmentEnabledAt(ordinal)==parent.familyAlignmentEnabledAt(ordinal) &&
				replay.familyIdentityGateAt(ordinal)==parent.familyIdentityGateAt(ordinal) &&
				replay.familyBlosumGateAt(ordinal)==parent.familyBlosumGateAt(ordinal) &&
				replay.familyEffectiveIdentityLowerAt(ordinal)==parent.familyEffectiveIdentityLowerAt(ordinal),
				"Packed output changed fixed reference/gate fields");
		}
		exported.verify();baseline.verify();verify();
		write(output.resolve("key_tiebreak.tsv"),new ByteBuilder("policy\tfamilies\tmembers\traw_tails"
			+"\theld_registers\tchanged_families\tchanged_keys\n")
			.append(policy.label).tab().append(families).tab().append(totalMembers)
			.tab().append(raw.tails()).tab().append(raw.postings())
			.tab().append(changedFamilies).tab().append(changedKeys).nl());
		write(output.resolve("preparation_inputs.tsv"),inputTable());
		write(output.resolve("key_diagnostics.tsv"),new ByteBuilder("policy\t").append(DDLKeyTieBreak.Diagnostics.HEADER).nl()
			.append(policy.label).tab().append(diagnostics.summary()).nl());
		final ByteBuilder histogram=new ByteBuilder("stored_uint16_key\ttied_family_bin_wins\n");
		for(int key=1;key<=65535;key++){
			final long count=diagnostics.tiedWins(key);if(count>0){histogram.append(key).tab().append(count).nl();}
		}
		write(output.resolve("tied_winner_keys.tsv"),histogram);
		write(output.resolve("PASS"),new ByteBuilder("FAMILY_KEY_TIEBREAK_ARM_PASS\n"));
	}

	private static final class PackedKeyWorker{
		PackedKeyWorker(FamilyTailCascade.SearchConfig config_,DDLKeyTieBreak.RawTailCounts raw_){
			config=config_;raw=raw_;sketch=new char[config.width];plain=new DDLFamilyProfiles.Profile(config.width);
			consensus=new DDLKeyTieBreak.KeySet(config.width);scratch=new DDLKeyTieBreak.Workspace(config.width);
			diagnostics=new DDLKeyTieBreak.Diagnostics[DDLKeyTieBreak.Policy.values().length];
			for(int arm=0;arm<diagnostics.length;arm++){diagnostics[arm]=new DDLKeyTieBreak.Diagnostics();}
		}
		KeyRarityRow[] build(MemberInput row,FamilyReference reference,char[] baseline){
			require(row.id.equals(reference.id) && row.rank==reference.rank &&
				row.members.length==reference.members &&
				Arrays.equals(new ProteinSequence(row.id,row.consensus).enc,reference.consensus.enc),
				"Packed member/consensus identity differs: "+row.id);
			plain.clear();long residues=0,words=0;
			for(byte[] member:row.members){
				residues+=member.length;
				words+=config.positionMode.positional() ?
					DDLFamilyProfiles.sketchPositional(member,config.geometry,config.seed,
						config.partitions,config.overlap,config.positionMode,sketch,config.blacklist) :
					DDLFamilyProfiles.sketch(member,config.geometry,config.seed,sketch,config.blacklist);
				plain.add(sketch);
			}
			for(int bucket=0;bucket<config.width;bucket++){
				require(plain.value(bucket,false)==baseline[bucket],
					"Packed default member reconstruction differs: "+row.id+" bucket="+bucket);
			}
			DDLFamilyProfiles.enumerateKeys(row.consensus,config.geometry,config.seed,
				config.positionMode,config.partitions,config.overlap,config.blacklist,consensus);
			final long seed=DDLKeyTieBreak.familySeed(row.id);
			final DDLKeyTieBreak.Policy[] policies=DDLKeyTieBreak.Policy.values();
			final KeyRarityRow[] result=new KeyRarityRow[policies.length];
			for(int arm=0;arm<policies.length;arm++){
				final char[] values=plain.packedValues(policies[arm],raw,seed,consensus,scratch,diagnostics[arm]);
				int changed=0;
				final ByteBuilder line=new ByteBuilder("tail\t").append(row.id).tab()
					.append(values.length).tab().append(row.members.length)
					.tab().append(residues).tab().append(words).append("\tnonempty\t");
				for(int b=0;b<values.length;b++){
					if(b>0){line.append(',');}
					line.append((int)values[b]);
					if(values[b]!=baseline[b]){changed++;}
				}
				line.tab();
				for(int b=0;b<values.length;b++){if(b>0){line.append(',');}line.append(plain.votes(b,false));}
				line.tab();
				for(int b=0;b<values.length;b++){if(b>0){line.append(',');}line.append(plain.occupied(b));}
				result[arm]=new KeyRarityRow(values,line.nl(),row.members.length,residues,words,changed);
			}
			return result;
		}
		final FamilyTailCascade.SearchConfig config;final DDLKeyTieBreak.RawTailCounts raw;
		final char[] sketch;final DDLFamilyProfiles.Profile plain;
		final DDLKeyTieBreak.KeySet consensus;final DDLKeyTieBreak.Workspace scratch;
		final DDLKeyTieBreak.Diagnostics[] diagnostics;
	}

	/** Ordinary search indexes do not retain attribution rows. Load the same tracked
	 * profile list explicitly, preserving uint16 width and complete live-family order.
	 */
	private static char[][] loadKeyRarityBaseline(FamilyTailCascade.LoadedClassifier parent,Path list,
			DDLFamilySearch.StructuralInputs inputs) throws Exception{
		final FamilyTailCascade.SearchConfig config=parent.config();
		final DDLFamilySearch.Index source=DDLFamilySearch.loadForAttributionStructural(list,config.width,
			config.geometry,config.seed,config.positionMode.sketch,config.positional ? config.partitions : 1,
			config.positional ? config.overlap : 0,config.blacklist,inputs);
		require(source.families()==parent.familyCount(),"Stored profile list must cover the complete baseline live roster");
		final HashMap<String,Integer> ordinals=new HashMap<String,Integer>(source.families()*2);
		for(int i=0; i<source.names.length; i++){
			require(ordinals.put(source.names[i],i)==null,"Duplicate stored profile identity in the rarity baseline");
		}
		final char[][] rows=new char[parent.familyCount()][];
		for(int i=0; i<rows.length; i++){
			final Integer ordinal=ordinals.remove("tail:"+parent.familyIdAt(i));
			require(ordinal!=null,"Stored profile missing from the baseline live roster");
			rows[i]=source.profileCopy(ordinal.intValue());
		}
		require(ordinals.isEmpty() && parent.matchesProfiles(rows),
			"Frozen rarity rows must reproduce the baseline search index exactly");
		inputs.verify();
		return rows;
	}

	/** Thread-local vote maps and member sketch; emitted rows contain no member sequences. */
	private interface ProfileRebuildWorker{
		KeyRarityRow build(MemberInput row,FamilyReference reference,char[] baseline);
	}

	private static final class KeyRarityWorker implements ProfileRebuildWorker{
		KeyRarityWorker(FamilyTailCascade.SearchConfig config_,DDLKeyRarity rarity){
			config=config_; sketch=new char[config.width];
			plain=new DDLFamilyProfiles.Profile(config.width); selected=new DDLFamilyProfiles.Profile(config.width,rarity);
		}
		public KeyRarityRow build(MemberInput row,FamilyReference reference,char[] baseline){
			require(row.id.equals(reference.id) && row.rank==reference.rank && row.members.length==reference.members &&
				Arrays.equals(new ProteinSequence(row.id,row.consensus).enc,reference.consensus.enc),
				"Member roster, counts and consensus must match the unchanged baseline family: "+row.id);
			plain.clear(); selected.clear(); long residues=0, words=0;
			for(byte[] member:row.members){
				residues=Math.addExact(residues,member.length);
				words=Math.addExact(words,config.positionMode.positional() ?
					DDLFamilyProfiles.sketchPositional(member,config.geometry,config.seed,
						config.partitions,config.overlap,config.positionMode,sketch,config.blacklist) :
					DDLFamilyProfiles.sketch(member,config.geometry,config.seed,sketch,config.blacklist));
				plain.add(sketch); selected.add(sketch);
			}
			final char[] values=new char[config.width]; int changed=0;
			for(int bucket=0; bucket<values.length; bucket++){
				if(plain.value(bucket,false)!=baseline[bucket]){
					throw new IllegalArgumentException("Default member reconstruction differs from baseline: "+
						row.id+" bucket="+bucket);
				}
				values[bucket]=selected.value(bucket,false);
				if(values[bucket]!=baseline[bucket]){changed++;}
				assert(plain.votes(bucket,false)==selected.votes(bucket,false) &&
					plain.occupied(bucket)==selected.occupied(bucket)) :
					"Tie-breaking changes only the representative, never maximal support or occupancy";
			}
			final ByteBuilder line=new ByteBuilder("tail\t").append(row.id).tab()
				.append(values.length).tab().append(row.members.length)
				.tab().append(residues).tab().append(words).append("\tnonempty\t");
			for(int i=0; i<values.length; i++){if(i>0){line.append(',');} line.append((int)values[i]);}
			line.tab();
			for(int i=0; i<values.length; i++){if(i>0){line.append(',');} line.append(selected.votes(i,false));}
			line.tab();
			for(int i=0; i<values.length; i++){if(i>0){line.append(',');} line.append(selected.occupied(i));}
			return new KeyRarityRow(values,line.nl(),row.members.length,residues,words,changed);
		}
		private final FamilyTailCascade.SearchConfig config;
		private final char[] sketch;
		private final DDLFamilyProfiles.Profile plain, selected;
	}
	/** Reuses primitive maps per worker; baseline and target counts are always real members. */
	private static final class SketchStrategyWorker implements ProfileRebuildWorker{
		SketchStrategyWorker(FamilyTailCascade.SearchConfig config,ProteinPositionMode targetMode,long bonus){
			final String pin=config.blacklist==null ? null : config.blacklist.sha80();
			final DDLKeyTieBreak.Policy policy=DDLKeyTieBreak.Policy.CONSENSUS_DIVERSITY;
			baselineContext=new DDLKeyTieBreak.PolicyContext(policy,null,config.width,config.geometry,config.seed,
				config.positionMode,config.partitions,config.overlap,pin,null,0);
			targetContext=new DDLKeyTieBreak.PolicyContext(policy,null,config.width,config.geometry,config.seed,
				targetMode,config.partitions,config.overlap,pin,null,bonus);
			baselineBuilder=new ddl.FamilySketchBuilder(config.width,config.geometry,config.seed,config.positionMode,
				config.partitions,config.overlap,config.blacklist,baselineContext);
			targetBuilder=new ddl.FamilySketchBuilder(config.width,config.geometry,config.seed,targetMode,
				config.partitions,config.overlap,config.blacklist,targetContext);
			baselineCounts=new DDLFamilyProfiles.Profile(config.width);targetCounts=new DDLFamilyProfiles.Profile(config.width);
			scratch=new DDLKeyTieBreak.Workspace(config.width);consensus=new DDLKeyTieBreak.KeySet(config.width);
			sketch=new char[config.width];reconstructed=new char[config.width];
			this.config=config;
		}
		public KeyRarityRow build(MemberInput row,FamilyReference reference,char[] baseline){
			require(row.id.equals(reference.id) && row.rank==reference.rank && row.members.length==reference.members &&
				Arrays.equals(new ProteinSequence(row.id,row.consensus).enc,reference.consensus.enc),
				"Strategy member roster/count/consensus differs from accepted baseline: "+row.id);
			baselineCounts.clear();targetCounts.clear();long words=0,residues=0;
			for(byte[] member:row.members){
				residues=Math.addExact(residues,member.length);
				final long baselineWords=baselineBuilder.sketchRaw(member,sketch);baselineCounts.add(sketch);
				final long targetWords=targetBuilder.sketchRaw(member,sketch);targetCounts.add(sketch);
				require(baselineWords==targetWords,"Strategy must preserve the exact valid input word count");
				words=Math.addExact(words,targetWords);
			}
			final long seed=DDLKeyTieBreak.familySeed(row.id);
			DDLFamilyProfiles.enumerateKeys(row.consensus,config.geometry,config.seed,config.positionMode,
				config.partitions,config.overlap,config.blacklist,consensus);
			baselineCounts.writePackedValues(baselineContext,seed,consensus,scratch,null,reconstructed);
			require(Arrays.equals(reconstructed,baseline),"Accepted consensus-policy W0 member reconstruction differs: "+row.id);
			DDLFamilyProfiles.enumerateKeys(row.consensus,config.geometry,config.seed,targetBuilder.positionMode(),
				config.partitions,config.overlap,config.blacklist,consensus);
			final char[] values=new char[config.width];
			targetCounts.writePackedValues(targetContext,seed,consensus,scratch,null,values);
			final ByteBuilder line=new ByteBuilder("tail\t").append(row.id).tab()
				.append(values.length).tab().append(row.members.length)
				.tab().append(residues).tab().append(words).append("\tnonempty\t");
			int changed=0;
			for(int bin=0;bin<values.length;bin++){
				if(bin>0){line.append(',');}
				line.append((int)values[bin]);
				if(values[bin]!=baseline[bin]){changed++;}
			}
			line.tab();
			for(int bin=0;bin<values.length;bin++){
				if(bin>0){line.append(',');}
				line.append(values[bin]==0 ? 0 : targetCounts.memberSupport(bin,values[bin]));
			}
			line.tab();
			for(int bin=0;bin<values.length;bin++){
				if(bin>0){line.append(',');}
				line.append(targetCounts.occupied(bin));
			}
			return new KeyRarityRow(values,line.nl(),row.members.length,residues,words,changed);
		}
		private final FamilyTailCascade.SearchConfig config;
		private final DDLKeyTieBreak.PolicyContext baselineContext,targetContext;
		private final ddl.FamilySketchBuilder baselineBuilder,targetBuilder;
		private final DDLFamilyProfiles.Profile baselineCounts,targetCounts;
		private final DDLKeyTieBreak.Workspace scratch;
		private final DDLKeyTieBreak.KeySet consensus;
		private final char[] sketch,reconstructed;
	}

	private static final class PlainProfileWorker{
		/** Binds one reusable plain sketch/profile accumulator to the baseline classifier recipe. */
		PlainProfileWorker(FamilyTailCascade.SearchConfig config_){
			config=config_; sketch=new char[config.width]; plain=new DDLFamilyProfiles.Profile(config.width);
		}
		/** Rebuilds one physical family's modal row using ordinary position-free sketching. */
		KeyRarityRow build(MemberInput row,FamilyReference reference){
			require(row.id.equals(reference.id) && row.rank==reference.rank && row.members.length==reference.members &&
				Arrays.equals(new ProteinSequence(row.id,row.consensus).enc,reference.consensus.enc),
				"Member roster, counts and consensus must match the unchanged baseline family: "+row.id);
			plain.clear(); long residues=0, words=0;
			for(byte[] member:row.members){
				residues=Math.addExact(residues,member.length);
				words=Math.addExact(words,DDLFamilyProfiles.sketch(member,config.geometry,config.seed,sketch,config.blacklist));
				plain.add(sketch);
			}
			final char[] values=new char[config.width];
			final ByteBuilder line=new ByteBuilder("tail\t").append(row.id).tab()
				.append(values.length).tab().append(row.members.length)
				.tab().append(residues).tab().append(words).append("\tnonempty\t");
			for(int i=0; i<values.length; i++){
				if(i>0){line.append(',');}
				values[i]=plain.value(i,false);
				line.append((int)values[i]);
			}
			line.tab();
			for(int i=0; i<values.length; i++){if(i>0){line.append(',');} line.append(plain.votes(i,false));}
			line.tab();
			for(int i=0; i<values.length; i++){if(i>0){line.append(',');} line.append(plain.occupied(i));}
			return new KeyRarityRow(values,line.nl(),row.members.length,residues,words,0);
		}
		private final FamilyTailCascade.SearchConfig config;
		private final char[] sketch;
		private final DDLFamilyProfiles.Profile plain;
	}
	private static final class KeyRarityRow{
		KeyRarityRow(char[] values_,ByteBuilder line_,int members_,long residues_,long words_,int changed_){
			values=values_; line=line_; members=members_; residues=residues_; words=words_; changed=changed_;
		}
		final char[] values;
		final ByteBuilder line;
		final int members, changed;
		final long residues, words;
	}

	/** Verifies source and emitted dependencies, including profile contents beneath the pinned list. */
	void verify(){
		try{keyInputs.verify();}
		catch(java.io.IOException e){throw new IllegalArgumentException("Raw key-policy input changed",e);}
		for(Map.Entry<Path,String> entry:pins.entrySet()){
			if(verificationMode==FamilyCalibration.VerificationMode.SHA80){
				if(structuralOutputs.contains(entry.getKey())){
					require("unmeasured".equals(entry.getValue()),
						"Structural output must remain explicitly unmeasured: "+entry.getKey());
					require(Files.isRegularFile(entry.getKey()) && !Files.isSymbolicLink(entry.getKey()),
						"Missing/nonregular structural preparation artifact: "+entry.getKey());
					rememberStructuralFile(entry.getKey());
				}else{
					require(entry.getValue().equals(DigestSuffix.file(entry.getKey().toString())),
						"Preparation input/output changed: "+entry.getKey());
				}
			}else{
				require(Files.isRegularFile(entry.getKey()) && !Files.isSymbolicLink(entry.getKey()),
					"Missing/nonregular structural preparation artifact: "+entry.getKey());
				rememberStructuralFile(entry.getKey());
			}
		}
	}

	/** Selects the pooled raw-DDL union once, then aligns only explicit panel-local ordinals.
	 * The fuller frontier raised identity gates for most families relative to old
	 * sparse +1, but its +1 replay also admitted extra4/100 negatives, so this
	 * remains an explicit raw_pair_union path rather than the default. See
	 * results/identity_padding_fuller_replays_20261116.
	 */
	private RawPairSelection scoreMetricPairs(ExecutorService pool, int workers, final Prepared[] live,
			final NegativePanel[] metricPanels, final PackagedCalibration[] calibrations,
			final MetricResult[][] evidence, final RawPairStats[][] rawEvidence,
			final WitnessPair[] witnesses){
		final RawPairSelection rawPairs=collectRawPairs(pool, live, metricPanels, workers);
		final int[][][] selected=rawPairs.selected;
		for(int p=0; p<metricPanels.length; p++){
			final int panelIndex=p;
			final String panelName=panels[p].name;
			final NegativePanel panel=metricPanels[p];
			final MetricWorker[] scorers=new MetricWorker[workers];
			for(int w=0; w<workers; w++){scorers[w]=calibrations[p].newMetricWorker(upperTop);}
			final int[][] panelSelected=selected[p];
			final MetricResult[] panelEvidence=evidence[p];
			final RawPairStats[] panelRawEvidence=rawEvidence[p];
			for(int f=0; f<live.length; f++){panelRawEvidence[live[f].rank]=rawPairs.stats[p][f];}
			parallel(pool, workers, live.length, new Work(){
				@Override public void run(int worker, int first, int end){
					for(int f=first; f<end; f++){
						final Prepared family=live[f]; final int rank=family.rank;
						final int[] ordinals=panelSelected[f]==null ? EMPTY_ORDINALS : panelSelected[f];
						final WitnessPair witness=witnesses[rank];
						panelEvidence[rank]=scorers[worker].scoreSelectedForPreparation(f, ordinals, ordinals.length,
							new MetricPairObserver(){
								@Override public void observe(int familyIndex, int queryOrdinal, int hits, double cosine,
										FamilyAlignmentMetric.Observation alignment){
									if(alignment.measured()){
										witness.observe(panelIndex, queryOrdinal, panelName, panel.queryId(queryOrdinal), alignment);
									}
								}
							});
					}
				}
			});
		}
		return rawPairs;
	}

	/** Builds a pooled positive-score frontier without retaining every query sketch at once. */
	private RawPairSelection collectRawPairs(final ExecutorService pool, final Prepared[] live,
			final NegativePanel[] metricPanels, final int workers){
		final int totalQueries=queryTotal(metricPanels);
		final int ownerWorkers=rawOwnerWorkers(workers, live.length);
		final int queryWorkers=workers-ownerWorkers;
		final int queryTop=Math.min(RAW_QUERY_TOP, live.length);
		requireArrayCapacity((long)totalQueries*queryTop, "raw query pair frontier");
		final DDLFamilySearch.Index index=rawPairIndex(live);
		final FamilyRawPairCollector.FamilyStripe[] familyHeaps=
			FamilyRawPairCollector.makeFamilyStripes(live.length, upperTop, Math.max(1, ownerWorkers));
		final long[] queryPairs=new long[(int)((long)totalQueries*queryTop)];
		final int[] pooledPanels=new int[totalQueries], pooledLocalOrdinals=new int[totalQueries];
		final RawPairStatsBuilder stats=new RawPairStatsBuilder(metricPanels.length, live.length, pooledPanels);
		RawPositiveRoute positiveRoute=null;
		int firstQuery=0, queryPairCount=0; long positivePairs=0;
		try{
			positiveRoute=ownerWorkers>0 ? new RawPositiveOwners(familyHeaps, live.length, queryWorkers, stats) :
				new DirectRawPositiveRoute(familyHeaps, live.length, stats);
			final RawTileWorker[] rawWorkers=new RawTileWorker[queryWorkers];
			for(int w=0; w<queryWorkers; w++){rawWorkers[w]=new RawTileWorker(live.length, positiveRoute.sink());}
			for(int p=0; p<metricPanels.length; p++){
				final NegativePanel panel=metricPanels[p];
				final int fullTile=Math.min(RAW_PAIR_TILE, panel.size());
				final FamilyRawPairCollector.Query[] fullRows=new FamilyRawPairCollector.Query[fullTile];
				final char[][] fullSketches=rawSketchRows(fullTile, width);
				for(int from=0; from<panel.size(); from+=RAW_PAIR_TILE){
					final int count=Math.min(RAW_PAIR_TILE, panel.size()-from);
					final FamilyRawPairCollector.Query[] rows=count==fullTile ? fullRows : new FamilyRawPairCollector.Query[count];
					final char[][] sketches=count==fullTile ? fullSketches : rawSketchRows(count, width);
					panel.rawPairQueries(p, from, rows, sketches);
					for(int i=0; i<count; i++){
						final int global=firstQuery+from+i;
						pooledPanels[global]=rows[i].panel;
						pooledLocalOrdinals[global]=rows[i].localOrdinal;
					}
					final RawTileResult result=collectRawTile(pool, index, rows, firstQuery+from, queryTop, rawWorkers, queryWorkers);
					System.arraycopy(result.queryPairs, 0, queryPairs, queryPairCount, result.queryPairCount);
					queryPairCount+=result.queryPairCount;
					positivePairs=Math.addExact(positivePairs, result.positivePairs);
					positiveRoute.rethrow();
				}
				firstQuery+=panel.size();
			}
			for(RawTileWorker worker:rawWorkers){worker.flushPositives();}
		}finally{
			if(positiveRoute!=null){positiveRoute.close();}
		}
		final FamilyRawPairCollector.Result union=FamilyRawPairCollector.finishStreaming(live.length, RAW_QUERY_TOP, upperTop,
			positivePairs, queryPairs, queryPairCount, pooledPanels, pooledLocalOrdinals, familyHeaps);
		stats.recordQueryFrontier(union);
		stats.recordFamilyFrontier(union);
		final long querySlots=(long)totalQueries*queryTop;
		final long familySlots=(long)live.length*Math.min(upperTop, totalQueries);
		return new RawPairSelection(selectedOrdinals(union, metricPanels.length, live.length, stats),
			stats.finish(), union.positivePairs, union.queryPairs, union.familyPairs, union.size(),
			totalQueries, querySlots, familySlots);
	}

	private static int rawOwnerWorkers(final int workers, final int families){
		if(workers<1 || families<1){
			throw new IllegalArgumentException("Raw pair selection requires live worker and family counts");
		}
		return workers<2 ? 0 : Math.max(1, Math.min(families, workers/4));
	}

	/** Collects one query tile, then merges worker-private frontiers after every worker is done. */
	private static RawTileResult collectRawTile(final ExecutorService pool, final DDLFamilySearch.Index index,
			final FamilyRawPairCollector.Query[] rows, final int firstQuery, final int queryTop,
			final RawTileWorker[] rawWorkers, final int workers){
		for(RawTileWorker worker:rawWorkers){worker.clearTile();}
		parallel(pool, workers, rows.length, new Work(){
			@Override public void run(int worker, int first, int end){
				rawWorkers[worker].collect(index, rows, first, end, firstQuery, queryTop);
			}
		});
		final RawTileResult merged=new RawTileResult();
		for(RawTileWorker worker:rawWorkers){merged.absorb(worker);}
		return merged;
	}

	/** The raw collector needs indexed uint16 rows in the same order as the live preparation roster. */
	private DDLFamilySearch.Index rawPairIndex(final Prepared[] live){
		final String[] names=new String[live.length];
		final char[][] profiles=new char[live.length][];
		final boolean positional=positionMode.positional();
		for(int f=0; f<live.length; f++){
			names[f]=live[f].id;
			profiles[f]=positional ? live[f].positional : live[f].plain;
		}
		final DDLFamilySearch.Index index=new DDLFamilySearch.Index(names, profiles, null, true);
		if(!index.supportsInversePostingRarity()){
			throw new IllegalStateException("Raw pair collection needs indexed uint16 DDL rows");
		}
		return index;
	}

	private static char[][] rawSketchRows(final int rows, final int width){
		final char[][] sketches=new char[rows][];
		for(int i=0; i<rows; i++){sketches[i]=new char[width];}
		return sketches;
	}

	private static int queryTotal(final NegativePanel[] panels){
		int total=0;
		for(NegativePanel panel:panels){total=Math.addExact(total, panel.size());}
		return total;
	}

	/** Converts the pooled union into sorted per-panel, per-family local ordinals for metric alignment. */
	private static int[][][] selectedOrdinals(final FamilyRawPairCollector.Result pairs,
			final int panelCount, final int families, final RawPairStatsBuilder stats){
		final int[][] counts=new int[panelCount][families];
		for(int i=0; i<pairs.size(); i++){
			final int panel=pairs.panelAt(i), family=pairs.familyAt(i);
			counts[panel][family]++;
			stats.recordUnion(panel, family);
		}
		final int[][][] selected=new int[panelCount][families][];
		for(int p=0; p<panelCount; p++){
			for(int f=0; f<families; f++){
				if(counts[p][f]>0){selected[p][f]=new int[counts[p][f]];}
			}
			Arrays.fill(counts[p], 0);
		}
		for(int i=0; i<pairs.size(); i++){
			final int p=pairs.panelAt(i), f=pairs.familyAt(i);
			selected[p][f][counts[p][f]++]=pairs.localOrdinalAt(i);
		}
		for(int p=0; p<panelCount; p++){
			for(int f=0; f<families; f++){
				final int[] ordinals=selected[p][f];
				if(ordinals==null){continue;}
				Arrays.sort(ordinals);
				for(int i=1; i<ordinals.length; i++){
					require(ordinals[i-1]!=ordinals[i], "Pooled raw pairs must select each family/query at most once");
				}
			}
		}
		return selected;
	}

	private static void requireArrayCapacity(final long cells, final String label){
		if(cells>Integer.MAX_VALUE){throw new IllegalArgumentException(label+" exceeds Java array bounds: "+cells);}
	}

	/** Pooled raw selection plus aggregate counts retained for the sidecar summary. */
	private static final class RawPairSelection{
		RawPairSelection(final int[][][] selected_, final RawPairStats[][] stats_,
				final long positivePairs_, final int queryPairs_, final int familyPairs_,
				final int unionPairs_, final int totalQueries_, final long querySlots_,
				final long familySlots_){
			selected=selected_; stats=stats_; positivePairs=positivePairs_;
			queryPairs=queryPairs_; familyPairs=familyPairs_; unionPairs=unionPairs_;
			totalQueries=totalQueries_; querySlots=querySlots_; familySlots=familySlots_;
		}
		ByteBuilder appendSummary(final ByteBuilder bb){
			final long overlap=queryPairs+(long)familyPairs-unionPairs;
			require(overlap>=0 && queryPairs<=querySlots && familyPairs<=familySlots,
				"Raw-pair aggregate summary must conserve the retained frontiers");
			return bb.append("positive_pairs\tquery_frontier\tfamily_frontier\tfrontier_overlap"+
				"\tunion_pairs\tpooled_queries\tquery_slots\tquery_shortfall"+
				"\tfamily_slots\tfamily_shortfall\n")
				.append(positivePairs).tab().append(queryPairs).tab().append(familyPairs).tab()
				.append(overlap).tab().append(unionPairs).tab().append(totalQueries).tab()
				.append(querySlots).tab().append(querySlots-queryPairs).tab()
				.append(familySlots).tab().append(familySlots-familyPairs).nl();
		}
		private final int[][][] selected;
		private final RawPairStats[][] stats;
		private final long positivePairs, querySlots, familySlots;
		private final int queryPairs, familyPairs, unionPairs;
		private final int totalQueries;
	}

	/** Mutable per-panel statistics, written before owner drains and read only after finishStreaming. */
	private static final class RawPairStatsBuilder{
		RawPairStatsBuilder(final int panels, final int families, final int[] queryPanels_){
			positives=new long[panels][families]; querySelected=new int[panels][families];
			familySelected=new int[panels][families]; unionSelected=new int[panels][families];
			queryPanels=queryPanels_;
		}
		void recordPositive(final int family, final int queryOrdinal){
			final int panel=queryPanels[queryOrdinal];
			positives[panel][family]=Math.addExact(positives[panel][family], 1L);
		}
		void recordQueryFrontier(final FamilyRawPairCollector.Result pairs){
			for(int i=0; i<pairs.queryFrontierSize(); i++){record(querySelected, pairs.queryFrontierAt(i));}
		}
		void recordFamilyFrontier(final FamilyRawPairCollector.Result pairs){
			for(int i=0; i<pairs.familyFrontierSize(); i++){record(familySelected, pairs.familyFrontierAt(i));}
		}
		void recordUnion(final int panel, final int family){
			unionSelected[panel][family]=Math.addExact(unionSelected[panel][family], 1);
		}
		RawPairStats[][] finish(){
			final RawPairStats[][] stats=new RawPairStats[positives.length][];
			for(int p=0; p<stats.length; p++){
				stats[p]=new RawPairStats[positives[p].length];
				for(int f=0; f<stats[p].length; f++){
					stats[p][f]=new RawPairStats(positives[p][f], querySelected[p][f], familySelected[p][f], unionSelected[p][f]);
				}
			}
			return stats;
		}
		private void record(final int[][] counts, final long pair){
			final int query=FamilyRawPairCollector.query(pair);
			final int panel=queryPanels[query], family=FamilyRawPairCollector.family(pair);
			counts[panel][family]=Math.addExact(counts[panel][family], 1);
		}
		private final long[][] positives;
		private final int[][] querySelected, familySelected, unionSelected;
		private final int[] queryPanels;
	}

	/** Immutable per-panel, per-family evidence copied into raw_calibration rows. */
	private static final class RawPairStats{
		RawPairStats(final long positive_, final int query_, final int family_, final int union_){
			final int overlap=Math.addExact(query_, family_)-union_;
			require(positive_>=0 && query_>=0 && family_>=0 && union_>=0 && overlap>=0,
				"Raw-pair evidence must conserve positive and selected counts");
			positive=positive_; query=query_; family=family_; union=union_; this.overlap=overlap;
		}
		String status(final MetricResult result){
			if(result.eligible>0){return "OK";}
			return positive>0 ? "ZERO_MEASURED" : "ZERO_CANDIDATES";
		}
		private final long positive;
		private final int query, family, overlap, union;
	}

	/** Query worker scratch; no other thread touches its workspace or query frontier buffer. */
	private static final class RawTileWorker{
		RawTileWorker(final int families, final RawPositiveSink positives_){
			workspace=new DDLFamilySearch.Workspace(families);
			positives=positives_;
		}
		void clearTile(){queryPairCount=0; positivePairs=0;}
		void collect(final DDLFamilySearch.Index index, final FamilyRawPairCollector.Query[] rows,
				final int first, final int end, final int firstQuery, final int queryTop){
			ensureQueryCapacity((long)(end-first)*queryTop);
			final FamilyRawPairCollector.QueryBatch batch=FamilyRawPairCollector.collectQueryFrontier(index, rows,
				first, end, firstQuery+first, queryTop, workspace, positives, queryPairs, queryPairCount);
			queryPairCount=batch.nextOffset;
			positivePairs=Math.addExact(positivePairs, batch.positivePairs);
		}
		void flushPositives(){positives.flush();}
		private void ensureQueryCapacity(final long additional){
			requireArrayCapacity((long)queryPairCount+additional, "raw tile query frontier");
			final int needed=queryPairCount+(int)additional;
			if(queryPairs.length>=needed){return;}
			int capacity=Math.max(64, queryPairs.length);
			while(capacity<needed){capacity=Math.max(needed, capacity*2);}
			queryPairs=Arrays.copyOf(queryPairs, capacity);
		}
		private final DDLFamilySearch.Workspace workspace;
		private final RawPositiveSink positives;
		private long[] queryPairs=new long[0];
		private int queryPairCount;
		private long positivePairs;
	}

	/** Single-threaded accumulator for worker-local query frontiers and positive-pair counts. */
	private static final class RawTileResult{
		void absorb(final RawTileWorker source){
			ensureQueryCapacity(source.queryPairCount);
			System.arraycopy(source.queryPairs, 0, queryPairs, queryPairCount, source.queryPairCount);
			queryPairCount+=source.queryPairCount;
			positivePairs=Math.addExact(positivePairs, source.positivePairs);
		}
		private void ensureQueryCapacity(final long additional){
			requireArrayCapacity((long)queryPairCount+additional, "raw tile query frontier");
			final int needed=queryPairCount+(int)additional;
			if(queryPairs.length>=needed){return;}
			int capacity=Math.max(64, queryPairs.length);
			while(capacity<needed){capacity=Math.max(needed, capacity*2);}
			queryPairs=Arrays.copyOf(queryPairs, capacity);
		}
		private long[] queryPairs=new long[0];
		private int queryPairCount;
		private long positivePairs;
	}

	/** Thread-local positive sink; flush publishes any buffered batches before final owner shutdown. */
	private interface RawPositiveSink extends FamilyRawPairCollector.RawPositiveObserver{
		void flush();
	}

	/** Routes raw positives either directly or through exclusive family-owner queues. */
	private interface RawPositiveRoute{
		RawPositiveSink sink();
		void close();
		void rethrow();
	}

	/** t=1 route; the query thread owns every final family stripe directly. */
	private static final class DirectRawPositiveRoute implements RawPositiveRoute{
		DirectRawPositiveRoute(final FamilyRawPairCollector.FamilyStripe[] stripes_,
				final int families_, final RawPairStatsBuilder stats_){
			stripes=stripes_; families=families_; stats=stats_;
		}
		@Override public RawPositiveSink sink(){
			return new RawPositiveSink(){
				@Override public void observe(final int family, final int queryOrdinal, final int hits){
					stats.recordPositive(family, queryOrdinal);
					stripes[FamilyRawPairCollector.owner(family, families, stripes.length)].offer(family, queryOrdinal, hits);
				}
				@Override public void flush(){}
			};
		}
		@Override public void close(){}
		@Override public void rethrow(){}
		private final FamilyRawPairCollector.FamilyStripe[] stripes;
		private final RawPairStatsBuilder stats;
		private final int families;
	}

	/** Owns final family stripes on dedicated threads so query workers never race on family heaps. */
	private static final class RawPositiveOwners implements RawPositiveRoute{
		RawPositiveOwners(final FamilyRawPairCollector.FamilyStripe[] stripes, final int families_, final int workers,
				final RawPairStatsBuilder stats){
			families=families_;
			queues=new RawPositiveQueue[stripes.length];
			threads=new Thread[stripes.length];
			int started=0;
			try{
				for(; started<stripes.length; started++){
					queues[started]=new RawPositiveQueue(stripes[started], failure, workers, stats);
					threads[started]=new Thread(queues[started], "raw-positive-owner-"+started);
					threads[started].setDaemon(true);
					threads[started].start();
				}
			}catch(Throwable error){
				stopStarted(started);
				throw error;
			}
		}
		@Override public RawPositiveSink sink(){return new QueuedRawPositiveSink(this);}
		void enqueue(final int owner, final RawPositiveBatch batch){
			while(true){
				rethrow();
				try{
					if(queues[owner].queue.offer(batch, 1, TimeUnit.SECONDS)){return;}
				}catch(InterruptedException e){
					Thread.currentThread().interrupt();
					throw new IllegalStateException("Interrupted while routing raw positives", e);
				}
			}
		}
		@Override public void close(){
			stopStarted(threads.length);
			rethrow();
		}
		private void stopStarted(final int count){
			boolean interrupted=false;
			for(int i=0; i<count; i++){
				if(threads[i]==null){continue;}
				boolean sent=false;
				while(!sent && threads[i].isAlive()){
					try{sent=queues[i].queue.offer(RawPositiveBatch.POISON, 1, TimeUnit.SECONDS);}
					catch(InterruptedException e){interrupted=true;}
				}
			}
			for(int i=0; i<count; i++){
				final Thread thread=threads[i];
				if(thread==null){continue;}
				boolean joined=false;
				while(!joined){
					try{thread.join(); joined=true;}
					catch(InterruptedException e){interrupted=true;}
				}
			}
			if(interrupted){Thread.currentThread().interrupt();}
		}
		@Override public void rethrow(){
			final Throwable error=failure.get();
			if(error instanceof Error){throw (Error)error;}
			if(error instanceof RuntimeException){throw (RuntimeException)error;}
			if(error!=null){throw new IllegalStateException("Raw positive owner failed", error);}
		}
		private final int families;
		private final RawPositiveQueue[] queues;
		private final Thread[] threads;
		private final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
	}

	/** Buffers one bounded batch per owner; flush drains each partial batch before owners stop. */
	private static final class QueuedRawPositiveSink implements RawPositiveSink{
		QueuedRawPositiveSink(final RawPositiveOwners owners_){
			owners=owners_;
			pending=new RawPositiveBatch[owners.queues.length];
		}
		@Override public void observe(final int family, final int queryOrdinal, final int hits){
			final int owner=FamilyRawPairCollector.owner(family, owners.families, pending.length);
			RawPositiveBatch batch=pending[owner];
			if(batch==null){batch=pending[owner]=new RawPositiveBatch();}
			if(batch.full()){
				owners.enqueue(owner, batch);
				batch=pending[owner]=new RawPositiveBatch();
			}
			batch.add(family, queryOrdinal, hits);
		}
		@Override public void flush(){
			for(int owner=0; owner<pending.length; owner++){
				final RawPositiveBatch batch=pending[owner];
				if(batch==null || batch.size<1){continue;}
				pending[owner]=null;
				owners.enqueue(owner, batch);
			}
		}
		private final RawPositiveOwners owners;
		private final RawPositiveBatch[] pending;
	}

	/** Single-owner drain loop; all writes to its stripe happen-before owner threads join before finishStreaming. */
	private static final class RawPositiveQueue implements Runnable{
		RawPositiveQueue(final FamilyRawPairCollector.FamilyStripe stripe_, final AtomicReference<Throwable> failure_,
				final int workers, final RawPairStatsBuilder stats_){
			stripe=stripe_;
			failure=failure_; stats=stats_;
			queue=new ArrayBlockingQueue<RawPositiveBatch>(Math.max(2, workers*2));
		}
		@Override public void run(){
			try{
				while(true){
					final RawPositiveBatch batch=queue.take();
					if(batch==RawPositiveBatch.POISON){return;}
					for(int i=0; i<batch.size; i++){
						final int family=FamilyRawPairCollector.family(batch.pairs[i]);
						final int query=FamilyRawPairCollector.query(batch.pairs[i]);
						stats.recordPositive(family, query);
						stripe.offer(family, query, batch.scores[i]);
					}
				}
			}catch(Throwable error){failure.compareAndSet(null, error);}
		}
		private final FamilyRawPairCollector.FamilyStripe stripe;
		private final RawPairStatsBuilder stats;
		private final AtomicReference<Throwable> failure;
		private final BlockingQueue<RawPositiveBatch> queue;
	}

	/** Fixed-size transfer packet from a query worker to one family owner. */
	private static final class RawPositiveBatch{
		void add(final int family, final int query, final int score){
			pairs[size]=FamilyRawPairCollector.pack(family, query);
			scores[size++]=score;
		}
		boolean full(){return size==pairs.length;}
		private static final RawPositiveBatch POISON=new RawPositiveBatch(0);
		private RawPositiveBatch(){this(RAW_POSITIVE_BATCH);}
		private RawPositiveBatch(final int capacity){
			pairs=new long[capacity];
			scores=new int[capacity];
		}
		private final long[] pairs;
		private final int[] scores;
		private int size;
	}

	private FamilyTailCascade.LoadedClassifier writePayloads(Path out, Prepared[] prepared, int live,
			long members_, long residues, long words, double[] cosine, double[] identity,
			double[] blosum, MetricResult[][] evidence, RawPairStats[][] rawEvidence,
			RawPairSelection rawPairs, WitnessPair[] witnesses, DDLPentamerBlacklist blacklist,
			int threads) throws Exception{
		final Path profileDir=out.resolve("profiles"); Files.createDirectory(profileDir);
		final Path profiles=profileDir.resolve("profiles.tsv.gz"),
			raw=out.resolve("raw_calibration.tsv.gz");
		final FamilyReferenceCodec.Header header=new FamilyReferenceCodec.Header(
			alignmentMetric.minimumPaired>0 ? FamilyReferenceCodec.PAIRED_MONOMER_SCHEMA :
			!alignmentMetric.matrix.same(ProteinSubstitutionMatrix.BLOSUM62) ?
				FamilyReferenceCodec.MATRIX_MONOMER_SCHEMA :
			metricAware ? FamilyReferenceCodec.ALIGNMENT_MONOMER_SCHEMA :
				FamilyReferenceCodec.MONOMER_SCHEMA,
			geometry.k(),width,0,families,live,families-live,
			FamilyReferenceCodec.FORMULA_VERSION,FamilyReferenceCodec.SPECTRUM_ENCODING,
			scalarZ,alignmentMetric);
		final FamilyReference[] references=new FamilyReference[families];
		final ByteStreamWriter profileWriter=new ByteStreamWriter(profiles.toString(), false, false, false);
		final ByteStreamWriter rawWriter=new ByteStreamWriter(raw.toString(), false, false, false);
		profileWriter.start(); rawWriter.start();
		final ByteBuilder line=new ByteBuilder();
		final MetricResult[] familyEvidence=metricAware ? new MetricResult[panels.length] : null;
		try{
			emit(profileWriter, line.clear().append(DDLFamilyProfiles.HEADER16).nl());
			line.clear().append("roster_rank\tfamily_id\tmembers\tcosine_lower\tcosine_retained"+
				"\tidentity_lower\tidentity_retained\tblosum_lower\tblosum_retained"+
				"\tcosine_upper\tidentity_upper\tblosum_upper");
			if(metricAware){
				line.append("\talignment_mode\talignment_min_span\tmember_alignment_measured"+
					"\tmember_alignment_failed\tidentity_quota_achievable"+
					"\tblosum_quota_achievable");
				for(Panel panel:panels){
					for(String field:NEGATIVE_EVIDENCE_FIELDS){line.tab().append(panel.name).append('_').append(field);}
					if(rawPairUnion){
						for(String field:RAW_PAIR_EVIDENCE_FIELDS){line.tab().append(panel.name).append('_').append(field);}
					}
				}
				if(rawPairUnion){for(String field:RAW_PAIR_POOLED_FIELDS){line.tab().append(field);}}
				if(rawPairUnion){
					line.append("\tidentity_negative_set\tidentity_negative_id\tidentity_alignment_columns\tidentity_paired_columns")
						.append("\tblosum_negative_set\tblosum_negative_id\tblosum_alignment_columns\tblosum_paired_columns");
				}
			}
			emit(rawWriter,line.nl());
			for(Prepared f:prepared){
				if(metricAware){for(int p=0;p<panels.length;p++){familyEvidence[p]=evidence[p][f.rank];}}
				final FamilyReference reference=f.members==0 ? FamilyReference.empty(f.id, f.rank) :
					FamilyReferenceExport.projectMonomers(metricAware ?
						FamilyReferenceExport.fromPrepared(f,cosine[f.rank],familyEvidence,
							blosumFloor,scalarZ,identityUpperMultiplier,identityUpperAdditive) :
						FamilyReferenceExport.fromPrepared(f,cosine[f.rank],identity[f.rank],blosum[f.rank],
							blosumFloor,scalarZ,identityUpperMultiplier,identityUpperAdditive));
				references[f.rank]=reference;
				line.clear();
				if(f.appendPositionalProfile(line)){emit(profileWriter, line.nl());}
				line.clear().append(f.rank).tab().append(f.id).tab().append(f.members);
				if(f.members==0){for(int i=0; i<9+(metricAware ? 6+panels.length*(NEGATIVE_EVIDENCE_FIELDS.length+
							(rawPairUnion ? RAW_PAIR_EVIDENCE_FIELDS.length : 0))+
							(rawPairUnion ? RAW_PAIR_POOLED_FIELDS.length : 0)+
							(rawPairUnion ? WITNESS_EVIDENCE_FIELDS : 0) : 0); i++){line.append("\tNA");}}
				else if(metricAware){
					final RawPairStats[] familyRawEvidence=rawPairUnion ? new RawPairStats[panels.length] : null;
					if(rawPairUnion){for(int p=0; p<panels.length; p++){familyRawEvidence[p]=rawEvidence[p][f.rank];}}
					appendMetricEvidence(line,f,cosine[f.rank],familyEvidence,familyRawEvidence,
						rawPairUnion ? witnesses[f.rank] : null, upperTop);
				}
				else{
					for(FamilyCalibration.Quantile lower:f.lower){line.tab().appendSlow(lower.cutoff).tab().append(lower.retained);}
					line.tab().appendSlow(cosine[f.rank]).tab().appendSlow(identity[f.rank]).tab().appendSlow(blosum[f.rank]);
				}
				emit(rawWriter, line.nl());
			}
		}finally{
			final boolean profileError=profileWriter.poisonAndWait(), rawError=rawWriter.poisonAndWait();
			require(!profileError && !rawError, "Preparation output writer failed; no final success marker");
		}
		final boolean legacyK5=geometry.isExact(5);
		final ByteBuilder summary=new ByteBuilder("namespace\tmanifest_sha80\tk\tseed\texponent"
			+"\tmanifest_families\tbuilt_families\tmembers\tresidues\tvalid_words\tsketch")
			.append(DDLFamilyProfiles.positionMetadataColumns(positionMode));
		if(!legacyK5){summary.append("\tgeometry");}
		if(blacklist!=null){
			summary.append("\tblacklist_sha80\tblacklist_words");
			if(FamilyBlacklistBinding.hasPolicy(options)){summary.tab().append(blacklist.metadataColumns());}
		}
		summary.nl().append("tail\t").append(required("familiessha80")).tab()
			.append(geometry.k()).append("\t0\t6\t").append(families).tab().append(families)
			.tab().append(members_).tab().append(residues).tab().append(words)
			.tab().append(positionMode.sketch)
			.append(DDLFamilyProfiles.positionMetadataValues(positionMode, partitions, overlap, geometry));
		if(!legacyK5){summary.tab().append(geometry.name());}
		if(blacklist!=null){
			summary.tab().append(blacklist.sha80()).tab().append(blacklist.size());
			if(FamilyBlacklistBinding.hasPolicy(options)){summary.tab().append(blacklist.metadataValues());}
		}
		final Path metadata=profileDir.resolve("summary.tsv"); write(metadata, summary.nl());
		final String profilePin=recordOutput(profiles), summaryPin=recordOutput(metadata);
		recordOutput(raw);
		final Path rawPairSummary=rawPairUnion ? out.resolve("raw_pair_summary.tsv") : null;
		if(rawPairUnion){
			require(rawPairs!=null, "Raw-pair union preparation must publish retained frontier metadata");
			write(rawPairSummary, rawPairs.appendSummary(new ByteBuilder()));
		}
		final Path profileList=out.resolve("ddlprofiles.tsv");
		write(profileList, new ByteBuilder("profiles/profiles.tsv.gz\t").append(profilePin).tab().append(summaryPin).nl());
		final HashMap<String,String> config=classifierOptions();
		if(!alignmentMetric.matrix.same(ProteinSubstitutionMatrix.BLOSUM62)){
			config.put("matrix",alignmentMetric.matrix.name);config.put("matrixbinding",alignmentMetric.matrix.binding());
		}
		memberRange.put(config);
		config.put("familyreferences", out.resolve("references.tsv.gz").toAbsolutePath().normalize().toString());
		config.put("ddlprofiles", profileList.toAbsolutePath().normalize().toString());
		config.put("ddlprofilessha80", recordOutput(profileList));
		if(rawPairUnion){recordStructuralOutput(rawPairSummary);}
		outputOptions=config;
		return FamilyScoringPackage.bindRanking(FamilyTailCascade.loadPrepared(
			FamilyScoringPackage.baseOptions(config,verificationMode),threads,header,references,prepared,
			verificationMode,blacklist),config,verificationMode);
	}

	/** Final audit retains all assigned members and every panel's selected alignment outcome. */
	private static void appendMetricEvidence(ByteBuilder line,Prepared family,double cosineUpper,MetricResult[] evidence,
			RawPairStats[] rawEvidence, WitnessPair witnesses, int upperTop){
		for(FamilyCalibration.MemberCutoff lower:family.measuredLower){
			line.tab();if(lower.measured()){line.appendSlow(lower.cutoff());}else{line.append("NA");}
			line.tab().append(lower.retained);
		}
		line.tab().appendSlow(cosineUpper);
		boolean measured=false;double identity=0,blosum=0;
		for(MetricResult result:evidence){
			require(result!=null,"Every declared panel must leave an explicit result for each live family");
			if(!result.measured()){continue;}
			if(!measured || result.identityUpper()>identity){identity=result.identityUpper();}
			if(!measured || result.blosumUpper()>blosum){blosum=result.blosumUpper();}
			measured=true;
		}
		if(measured){line.tab().appendSlow(identity).tab().appendSlow(blosum);}else{line.append("\tNA\tNA");}
		line.tab().append(family.alignmentMetric.name()).tab().append(family.alignmentMetric.minimumSpan)
			.tab().append(family.measuredLower[1].measured).tab().append(family.measuredLower[1].failed)
			.tab().append(family.measuredLower[1].quotaAchievable ? "true" : "false")
			.tab().append(family.measuredLower[2].quotaAchievable ? "true" : "false");
		for(int i=0; i<evidence.length; i++){
			final MetricResult result=evidence[i];
			line.tab().append(result.scanned).tab().append(result.admitted)
				.tab().append(result.selected).tab().append(result.eligible);
			for(FamilyAlignmentMetric.State state:FAILED_ALIGNMENT_STATES){line.tab().append(result.count(state));}
			if(rawEvidence!=null){appendRawPairEvidence(line, rawEvidence[i], result);}
		}
		if(rawEvidence!=null){appendPooledRawPairEvidence(line, rawEvidence, evidence, upperTop);}
		if(witnesses!=null){
			appendWitness(line, witnesses.identity);
			appendWitness(line, witnesses.blosum);
		}
	}

	private static void appendRawPairEvidence(ByteBuilder line, RawPairStats raw, MetricResult result){
		require(raw!=null, "Every raw-pair panel must publish explicit selection evidence");
		line.tab().append(raw.positive).tab().append(raw.query).tab().append(raw.family).tab().append(raw.overlap)
			.tab().append(raw.union).tab().append(raw.status(result));
	}

	private static void appendPooledRawPairEvidence(ByteBuilder line, RawPairStats[] rawEvidence,
			MetricResult[] evidence, int upperTop){
		long positive=0; int union=0, measured=0;
		for(int i=0; i<rawEvidence.length; i++){
			final RawPairStats raw=rawEvidence[i];
			require(raw!=null, "Every raw-pair panel must contribute pooled selection evidence");
			positive=Math.addExact(positive, raw.positive);
			union=Math.addExact(union, raw.union);
			measured=Math.addExact(measured, evidence[i].eligible);
		}
		final int unionShortfall=Math.max(0, upperTop-union), measuredShortfall=Math.max(0, upperTop-measured);
		final String status=measured>0 ? "OK" : positive>0 ? "ZERO_MEASURED" : "ZERO_CANDIDATES";
		line.tab().append(positive).tab().append(union).tab().append(unionShortfall)
			.tab().append(measured).tab().append(measuredShortfall).tab().append(status);
	}

	private static void appendWitness(ByteBuilder line, NegativeWitness witness){
		if(witness==null){line.append("\tNA\tNA\tNA\tNA");}
		else{line.tab().append(witness.negativeSet).tab().append(witness.negativeId)
			.tab().append(witness.alignmentColumns).tab().append(witness.pairedColumns);}
	}

	/** Independent best measured identity and BLOSUM witnesses with panel/local ordinal ties. */
	private static final class WitnessPair{
		void observe(int panel, int ordinal, String negativeSet, String negativeId,
				FamilyAlignmentMetric.Observation alignment){
			final double identityValue=alignment.identityPct(), blosumValue=alignment.blosum();
			if(identityBetter(identityValue, panel, ordinal)){identity=new NegativeWitness(negativeSet, negativeId,
				panel, ordinal, alignment.columns, alignment.pairedColumns, identityValue, blosumValue);}
			if(blosumBetter(blosumValue, panel, ordinal)){blosum=new NegativeWitness(negativeSet, negativeId,
				panel, ordinal, alignment.columns, alignment.pairedColumns, identityValue, blosumValue);}
		}
		private boolean identityBetter(double value, int panel, int ordinal){
			return identity==null || value>identity.identity || value==identity.identity && before(panel, ordinal, identity);
		}
		private boolean blosumBetter(double value, int panel, int ordinal){
			return blosum==null || value>blosum.blosum || value==blosum.blosum && before(panel, ordinal, blosum);
		}
		private static boolean before(int panel, int ordinal, NegativeWitness previous){
			return panel<previous.panel || panel==previous.panel && ordinal<previous.ordinal;
		}
		NegativeWitness identity, blosum;
	}

	/** Scalar copy of the borrowed alignment observation that supplied one raw-pair maximum. */
	private static final class NegativeWitness{
		NegativeWitness(String negativeSet_, String negativeId_, int panel_, int ordinal_,
				int alignmentColumns_, int pairedColumns_, double identity_, double blosum_){
			negativeSet=negativeSet_; negativeId=negativeId_; panel=panel_; ordinal=ordinal_;
			alignmentColumns=alignmentColumns_; pairedColumns=pairedColumns_;
			identity=identity_; blosum=blosum_;
		}
		final String negativeSet, negativeId;
		final int panel, ordinal, alignmentColumns, pairedColumns;
		final double identity, blosum;
	}

	private HashMap<String,String> classifierOptions(){
		final HashMap<String,String> result=new HashMap<String,String>();
		for(String key:CLASSIFIER){if(options.containsKey(key)){result.put(key, options.get(key));}}
		if(!positionMode.positional()){result.remove("partitions"); result.remove("overlap");}
		return result;
	}

	/** The shared loader preserves FASTA ordinals; the enclosing input contract rejects duplicate identities. */
	private static void validateIds(NegativePanel panel, String name){
		final HashSet<String> seen=new HashSet<String>();
		for(int q=0; q<panel.size(); q++){
			final String id=panel.queryId(q);
			require(id!=null && !id.isEmpty() && seen.add(id), "Duplicate or empty negative query identity in "+name+": "+id);
		}
	}

	private void readPanels(ArrayList<Panel> list, ArrayList<Integer> counts){
		final Path table=inputPath("negativepanels"); final ByteFile input=ByteFile.makeByteFile(table.toString(), false);
		final HashSet<String> names=new HashSet<String>(); names.add("positive");
		final HashSet<Path> paths=new HashSet<Path>();
		try{
			final byte[] header=input.nextLine();
			require(header!=null &&
				new String(header, StandardCharsets.US_ASCII).equals("panel\tpath\tsha80\tqueries"),
				"Negative table requires panel/path/sha80/queries header");
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				final String[] row=new String(line, StandardCharsets.UTF_8).split("\t", -1);
				require(row.length==4 && names.add(row[0]), "Malformed or duplicate negative panel");
				final Path file=table.getParent().resolve(row[1]).normalize();
				require(!row[1].isEmpty() && paths.add(realPath(file)), "Repeated or missing negative path");
				final long count=Parse.parseKMG(row[3]);
				require(count>=upperTop && count<=FamilyNegativeUpperCalibration.MAX_PANEL_QUERIES,
					"Negative counts must be at least uppertop and at most"+FamilyNegativeUpperCalibration.MAX_PANEL_QUERIES);
				pin(file, row[2]);
				list.add(new Panel(row[0], file.toString(), row[2], false, verificationMode, count));
				counts.add((int)count);
			}
		}finally{require(!input.close(), "Negative table read failed");}
		require(!list.isEmpty(), "Preparation requires at least one negative panel");
	}

	private void resolveInput(String key, boolean mandatory){
		if(!options.containsKey(key)){
			require(!mandatory && !options.containsKey(key+"sha80"), "Missing preparation path: "+key); return;
		}
		final Path file=path.getParent().resolve(required(key)).normalize();
		final String identity=required(key+"sha80");
		pin(file,identity);
		if(!key.equals("blacklist") && identity.equals("unmeasured")){unmeasuredInputs=true;}
		options.put(key,file.toString());
	}
	private void pin(Path file, String expected){
		require(Files.isRegularFile(file) && expected!=null && (expected.matches("[0-9a-f]{20}") ||
			verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL && expected.equals("unmeasured")),
			"Preparation requires a regular file and valid declared identity: "+file);
		final String previous=pins.put(file, expected);
		require(previous==null || previous.equals(expected), "Conflicting fingerprints for one preparation file: "+file);
		if(verificationMode==FamilyCalibration.VerificationMode.SHA80){
			require(expected.equals(DigestSuffix.file(file.toString())), "Preparation fingerprint differs: "+file);
		}else{
			require(file.isAbsolute() && !Files.isSymbolicLink(file),
				"Structural inputs require absolute non-symlink paths: "+file);
			rememberStructuralFile(file);
		}
	}
	private String recordOutput(Path file){
		require(Files.isRegularFile(file) && !Files.isSymbolicLink(file), "Final preparation output must exist: "+file);
		final String pin=verificationMode==FamilyCalibration.VerificationMode.SHA80 ?
			DigestSuffix.file(file.toString()) : "unmeasured";
		if(verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL){rememberStructuralFile(file);}
		pins.put(file, pin); return pin;
	}
	/** Tracks a new sidecar as structural/unmeasured even when legacy artifacts are SHA80-bound. */
	private String recordStructuralOutput(Path file){
		require(Files.isRegularFile(file) && !Files.isSymbolicLink(file),
			"Final structural preparation output must exist: "+file);
		rememberStructuralFile(file);
		structuralOutputs.add(file);
		pins.put(file, "unmeasured"); return "unmeasured";
	}
	Path recipePath(){return path;}
	String matrixName(){return alignmentMetric.matrix.name;}
	/** Validated ranking bounds govern held-out queries; legacy recipes retain their100-aa contract. */
	int minimumHeldoutQueryAa(){
		return options.containsKey("minqueryaa") ? Integer.parseInt(options.get("minqueryaa")) : 100;
	}
	int maximumHeldoutQueryAa(){
		return options.containsKey("maxqueryaa") ? Integer.parseInt(options.get("maxqueryaa")) : 100;
	}
	/** Null denotes no filter; an unmeasured new list must not be reported as inherited input identity. */
	String blacklistIdentity(){return options.get("blacklistsha80");}
	/** Input provenance excludes newly written products and the unmeasured recipe itself. */
	String inputIdentity(){
		return unmeasuredInputs ? "inherited_with_unmeasured_inputs" :
			"unmeasured".equals(blacklistIdentity()) ? "inherited_with_unmeasured_blacklist" : "inherited";
	}
	/** Metadata detects replacement/truncation during a run; it is not a claim of byte identity. */
	private void rememberStructuralFile(Path file){
		final java.nio.file.attribute.BasicFileAttributes current;
		try{
			current=Files.readAttributes(file,java.nio.file.attribute.BasicFileAttributes.class,
				java.nio.file.LinkOption.NOFOLLOW_LINKS);
		}
		catch(java.io.IOException e){throw new IllegalArgumentException("Cannot inspect structural artifact: "+file,e);}
		require(current.isRegularFile(),"Structural artifact must remain a regular file: "+file);
		final java.nio.file.attribute.BasicFileAttributes previous=structuralFiles.putIfAbsent(file,current);
		require(previous==null || previous.size()==current.size() &&
			previous.lastModifiedTime().equals(current.lastModifiedTime()) &&
			java.util.Objects.equals(previous.fileKey(),current.fileKey()),
			"Structural artifact metadata changed during preparation: "+file);
	}
	private ByteBuilder inputTable(){
		final ByteBuilder table=new ByteBuilder("path\tsha80\n");
		for(Map.Entry<Path,String> e:pins.entrySet()){table.append(e.getKey().toString()).tab().append(e.getValue()).nl();}
		for(Path input:keyRawInputs){
			if(!pins.containsKey(input)){table.append(input.toString()).tab().append("unmeasured").nl();}
		}
		return table;
	}
	private void emit(ByteStreamWriter writer, ByteBuilder line){
		outputBytes=Math.addExact(outputBytes, line.length());
		require(outputBytes<=maxOutputBytes, "Decoded preparation output exceeds maxoutputbytes"); writer.print(line);
	}
	private void write(Path file, ByteBuilder bytes){
		final ByteStreamWriter writer=new ByteStreamWriter(file.toString(), false, false, false); writer.start();
		try{emit(writer, bytes);}finally{require(!writer.poisonAndWait(), "Preparation output failed: "+file);}
	}
	private Path inputPath(String key){return java.nio.file.Paths.get(required(key));}
	private String required(String key){
		final String value=options.get(key);
		require(value!=null, "Preparation must explicitly specify "+key);
		return value;
	}
	private long number(String key, long minimum){
		final long value=Parse.parseKMG(required(key));
		require(value>=minimum, "Preparation "+key+" must be at least "+minimum);
		return value;
	}
	private int integer(String key, int minimum, int maximum){
		final long value=number(key, minimum);
		require(value<=maximum, "Preparation "+key+" exceeds "+maximum);
		return (int)value;
	}
	private double decimal(String key, double minimum, double maximum){
		final double value=Double.parseDouble(required(key));
		require(Double.isFinite(value) && value>=minimum && value<=maximum,
			"Preparation "+key+" must be finite in ["+minimum+","+maximum+"]");
		return value;
	}
	private double identityPadding(String key,double fallback,double minimum){
		if(!options.containsKey(key)){return fallback;}
		final double value=Double.parseDouble(required(key));
		require(Double.isFinite(value) && value>=minimum, "Identity padding is below its finite minimum: "+key);
		return value;
	}
	private boolean bool(String key){
		final String value=required(key).toLowerCase(Locale.ROOT);
		require(value.equals("t") || value.equals("true") || value.equals("f") || value.equals("false"),
			"Preparation "+key+" requires t/f");
		return Parse.parseBoolean(value);
	}
	private static Path realPath(Path file){
		try{return file.toRealPath();}
		catch(java.io.IOException e){throw new IllegalArgumentException("Cannot resolve negative input: "+file, e);}
	}

	private interface Work{void run(int worker, int first, int end);}
	/** Coarse tasks claim up to eight families at a time; small rosters still exercise all workers. */
	private static void parallel(ExecutorService pool, int workers, int count, Work work){
		assert(workers>0 && count>=0) : "Preparation phases require bounded worker and family counts";
		final AtomicInteger next=new AtomicInteger();
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final int block=Math.max(1, Math.min(8, count/(workers*4)));
		final ArrayList<Future<?>> pending=new ArrayList<Future<?>>(workers);
		for(int i=0; i<workers; i++){
			final int worker=i;
			pending.add(pool.submit(new Runnable(){
				@Override public void run(){
					try{
						for(int first=next.getAndAdd(block); first<count && failure.get()==null; first=next.getAndAdd(block)){
							work.run(worker, first, Math.min(count, first+block));
						}
					}catch(Throwable error){failure.compareAndSet(null, error);}
				}
			}));
		}
		boolean interrupted=false;
		for(Future<?> future:pending){
			boolean joined=false;
			while(!joined){
				try{future.get(); joined=true;}
				catch(InterruptedException e){interrupted=true; failure.compareAndSet(null, e);}
				catch(ExecutionException e){failure.compareAndSet(null, e.getCause()); joined=true;}
			}
		}
		if(interrupted){Thread.currentThread().interrupt();}
		final Throwable error=failure.get();
		if(error instanceof Error){throw (Error)error;}
		if(error instanceof RuntimeException){throw (RuntimeException)error;}
		if(error!=null){throw new IllegalStateException("Preparation worker phase failed", error);}
	}
	private static void shutdown(ExecutorService pool){
		pool.shutdown(); boolean interrupted=false;
		while(!pool.isTerminated()){
			try{pool.awaitTermination(1, TimeUnit.SECONDS);}catch(InterruptedException e){interrupted=true;}
		}
		if(interrupted){Thread.currentThread().interrupt();}
	}
	private static void require(boolean valid, String reason){if(!valid){throw new IllegalArgumentException(reason);}}

	final Panel[] panels;
	private final Path path;
	private final FamilyCalibration.VerificationMode verificationMode;
	private boolean unmeasuredInputs;
	private final HashMap<String,String> options=new HashMap<String,String>();
	/** Final ranking dependencies and policy; the reference fingerprint is added after calibration. */
	private HashMap<String,String> outputOptions;
	private final LinkedHashMap<Path,String> pins=new LinkedHashMap<Path,String>();
	private final HashSet<Path> structuralOutputs=new HashSet<Path>();
	private final HashMap<Path,java.nio.file.attribute.BasicFileAttributes> structuralFiles=
		new HashMap<Path,java.nio.file.attribute.BasicFileAttributes>();
	private final ProteinKmerGeometry geometry;
	private final ProteinPositionMode positionMode;
	private final MemberLengthRange memberRange;
	private final int[] queryCounts;
	private final int families, buckets, batches, maxReference, upperTop, minNegative, maxNegative,
		partitions, overlap, cosinePanel, width;
	private final long members, maxBucketMembers, maxBucketBytes, maxPreparedBytes, maxOutputBytes;
	private final double blosumFloor, scalarZ, identityUpperMultiplier, identityUpperAdditive;
	private final FamilyCalibration.MemberRetention retention;
	private final FamilyAlignmentMetric alignmentMetric;
	private final boolean metricAware, rawPairUnion;
	private long outputBytes;
	private boolean built, published;
	private DDLKeyTieBreak.PolicyContext keyPolicy;
	private final DDLFamilySearch.StructuralInputs keyInputs=new DDLFamilySearch.StructuralInputs();
	private Path[] keyRawInputs=new Path[0];
	private static final String[] INPUTS={"families", "chunks", "consensus", "negativepanels", "scalarbounds"};
	private static final String[] CLASSIFIER={"k", "seed", "width", "aabits", "aapartition",
		"positional", "positionmode", "partitions", "overlap", "scalarbounds",
		"scalarboundssha80", "scalarscale", "scalarz", "top", "ddlreject", "ddlaccept",
		"cosine", "rankingfeatures", "scalargates", "identityacceptfloor",
		"blosumacceptfloor", "indexearlyreject", "indexdirect", "indexorder", "blacklist",
		"blacklistsha80", "blacklistk", "blacklistmode", "blacklistwords",
		"blacklistprovenance", "alignmentmode", "alignmentminspan", "alignmentminpaired",
		"heuristic", "lengthmedians", "lengthmedianssha80", "familyhydro", "familyhydrosha80",
		"minqueryaa", "maxqueryaa"};
	private static final String[] REQUIRED={"schema", "families", "familiessha80", "chunks",
		"chunkssha80", "consensus", "consensussha80", "negativepanels", "negativepanelssha80",
		"cosinepanel", "expectedfamilies", "expectedmembers", "buckets", "batches",
		"maxreferenceaa", "maxbucketmembers", "maxbucketbytes", "maxpreparedbytes",
		"maxoutputbytes", "retention", "uppertop", "minnegativeaa", "maxnegativeaa", "k",
		"seed", "width", "positional", "partitions", "overlap", "scalarbounds",
		"scalarboundssha80", "scalarscale", "blosumfloor", "scalarz", "monomer", "top",
		"ddlreject", "ddlaccept", "cosine", "rankingfeatures", "scalargates",
		"identityacceptfloor", "blosumacceptfloor", "indexearlyreject"};
	private static final HashSet<String> ALLOWED=new HashSet<String>(Arrays.asList(REQUIRED));
	static{ALLOWED.addAll(Arrays.asList("aabits", "aapartition", "blacklist", "blacklistsha80",
		"indexdirect", "indexorder", "minmemberaa", "maxmemberaa", "positionmode",
		"alignmentmode", "alignmentminspan", "alignmentminpaired", "heuristic",
		"lengthmedians", "lengthmedianssha80", "familyhydro", "familyhydrosha80",
		"minqueryaa", "maxqueryaa", "blacklistk", "blacklistmode", "blacklistwords",
		"blacklistprovenance", "matrix", "matrixfile", "matrixfilesha80", "matrixbinding",
		"keytiebreak", "consensusbonus", "rawmanifest", "expectedrawtails", "rawminlen", "rawmaxlen", "rawtailcrop",
		"identityuppermultiplier", "identityupperadditive", "negativeselection"));}
	private static final String[] NEGATIVE_EVIDENCE_FIELDS={"scanned", "admitted", "selected",
		"eligible", "input_too_short", "no_alignment", "span_too_short"};
	private static final String[] RAW_PAIR_EVIDENCE_FIELDS={"raw_positive", "raw_query",
		"raw_family", "raw_overlap", "raw_union", "coverage_status"};
	private static final String[] RAW_PAIR_POOLED_FIELDS={"raw_pooled_positive",
		"raw_pooled_union", "raw_pooled_union_shortfall",
		"raw_pooled_measured","raw_pooled_measured_shortfall","raw_pooled_coverage_status"};
	private static final int RAW_PAIR_TILE=8192;
	private static final int RAW_POSITIVE_BATCH=4096;
	private static final int RAW_QUERY_TOP=100;
	private static final int WITNESS_EVIDENCE_FIELDS=8;
	private static final int[] EMPTY_ORDINALS=new int[0];
	private static final FamilyAlignmentMetric.State[] FAILED_ALIGNMENT_STATES={
		FamilyAlignmentMetric.State.INPUT_TOO_SHORT,
		FamilyAlignmentMetric.State.NO_ALIGNMENT,
		FamilyAlignmentMetric.State.SPAN_TOO_SHORT};
}
