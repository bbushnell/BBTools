package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import ddl.DDLPentamerBlacklist;
import ddl.DDLFamilyProfiles;
import ddl.ProteinBlacklistPolicy;
import ddl.ProteinKmerGeometry;
import parse.LineParser1;
import prok.CallGenes;
import prok.StopTailCascadeResources;
import prok.StopTailCascadeSupport;
import prot.FamilyCalibration.VerificationMode;
import prot.FamilyCompactCodec;
import prot.FamilyCompactLoader;
import prot.FamilyEval;
import prot.FamilyReferenceCodec;
import prot.FamilyScoringPackage;
import prot.FamilyTailCascade;
import prot.ProteinSubstitutionMatrix;
import structures.ByteBuilder;

/** Direct callable fixture for embedded compact references and sidecar ranking inputs.
 * @author Yelan
 */
public final class FamilyCompactCallableTest{

	public static void main(String[] args)throws Exception{
		if(args.length!=2 || !args[0].startsWith("fixture=") || !args[1].startsWith("alphabet=")){
			throw new IllegalArgumentException("Require fixture=NEW_DIR alphabet=ACCEPTED_Q70_ALPHABET on Dori");
		}
		final Path root=Paths.get(args[0].substring(8)).toAbsolutePath().normalize();
		final String aabits=readAcceptedAlphabet(Paths.get(args[1].substring(9)));
		Files.createDirectory(root);

		final Descriptor policy=new Descriptor(aabits);
		final Path blacklist=blacklist(root.resolve("blacklist.faa"), policy);
		final HashMap<String,String> options=options(blacklist,policy);
		final FamilyCompactLoader.Loaded loaded=FamilyCompactLoader.load(
			compact(root.resolve("compact.txt"),options).toString(),3,3);
		final DDLPentamerBlacklist filter=DDLPentamerBlacklist.loadReducedStructural(
			blacklist,"unmeasured",policy.words,policy.geometry,policy.provenance);
		reject(()->FamilyTailCascade.loadCompactCallable(replaced(options,"blacklistprovenance","wrong"),
			2,loaded,filter),"metadata differs");
		reject(()->FamilyTailCascade.loadCompactCallable(replaced(options,"aabits",ProteinBlacklistPolicy.ALPHABET),
			2,loaded,filter),"metadata differs");
		final Descriptor wrongShapePolicy=new Descriptor(mergedGroups(aabits));
		final Path wrongShapeBlacklist=blacklist(root.resolve("wrong_shape.faa"),wrongShapePolicy);
		final HashMap<String,String> wrongShapeOptions=options(wrongShapeBlacklist,wrongShapePolicy);
		final FamilyCompactLoader.Loaded wrongShapeLoaded=FamilyCompactLoader.load(
			compact(root.resolve("wrong_shape.txt"),wrongShapeOptions).toString(),3,3);
		final DDLPentamerBlacklist wrongShapeFilter=DDLPentamerBlacklist.loadReducedStructural(
			wrongShapeBlacklist,"unmeasured",wrongShapePolicy.words,
			wrongShapePolicy.geometry,wrongShapePolicy.provenance);
		reject(()->FamilyTailCascade.loadCompactCallable(wrongShapeOptions,2,
			wrongShapeLoaded,wrongShapeFilter),"metadata differs");

		final FamilyTailCascade.LoadedClassifier bare=
			FamilyTailCascade.loadCompactCallable(options,2,loaded,filter);
		check(bare.familyCount()==2,"Compact callable loader keeps only live families");
		check(!bare.usesHeuristicRanking() && !bare.usesFamilyGcPenalty(),
			"Compact sidecars are opt-in");

		options.put("heuristic","hitfourth_h020_c150");
		options.put("minqueryaa","40");
		options.put("maxqueryaa","99");
		final FamilyTailCascade.LoadedClassifier heuristic=
			FamilyTailCascade.loadCompactCallable(options,2,loaded,filter);
		check(heuristic.familyCount()==2 && heuristic.usesHeuristicRanking() &&
			!heuristic.usesFamilyGcPenalty(),"Compact length/KD sidecars bind heuristic ranking");

		options.put("familygclambda","8");
		options.put("gcfilter","t");
		options.put("gc_z","3");
		final FamilyTailCascade.LoadedClassifier gc=
			FamilyTailCascade.loadCompactCallable(options,1,loaded,filter);
		check(gc.usesHeuristicRanking() && gc.usesFamilyGcPenalty(),
			"Compact GC sidecars bind package GC ranking");
		final Path pkg=packageRoot(root.resolve("package"),options,policy);
		reject(()->new FamilyScoringPackage(pkg.resolve("short").resolve("families.tsv")),
			"structural loading");
		final FamilyTailCascade.LoadedClassifier packaged=
			new FamilyScoringPackage(pkg.resolve("short").resolve("families.tsv"),
				VerificationMode.STRUCTURAL).load(2);
		check(packaged.familyCount()==2 && packaged.usesHeuristicRanking() &&
			packaged.usesFamilyGcPenalty(),"Compact package root binds shared resources and embedded sidecars");
		rejectMutation(packageRoot(root.resolve("mutate_family"),options,policy),
			"short/families.tsv");
		rejectMutation(packageRoot(root.resolve("mutate_alphabet"),options,policy),
			"alphabet.config");
		rejectMutation(packageRoot(root.resolve("mutate_blacklist"),options,policy),
			"blacklist.faa.gz");
		rejectMutation(packageRoot(root.resolve("mutate_matrix"),options,policy),
			"POSSUM55.tsv");
		final HashMap<String,String> bareOptions=new HashMap<String,String>(options);
		bareOptions.remove("heuristic");
		bareOptions.remove("minqueryaa");
		bareOptions.remove("maxqueryaa");
		bareOptions.remove("familygclambda");
		bareOptions.remove("gcfilter");
		bareOptions.remove("gc_z");
		bareOptions.put("heuristic","hitfourth_h020_c150");
		bareOptions.put("minqueryaa","40");
		bareOptions.put("maxqueryaa","99");
		final Path evalSource=packageRoot(root.resolve("eval_package"),bareOptions,policy)
			.resolve("short").resolve("families.tsv");
		replace(evalSource,"maxqueryaa=99","maxqueryaa=99\tmarginscale=0.985");
		compactEvaluation(root.resolve("eval"),evalSource);
		marginEvaluation(root.resolve("margin_eval"),bareOptions,policy);
		final Path gcEvalSource=packageRoot(root.resolve("eval_gc_package"),options,policy)
			.resolve("short").resolve("families.tsv");
		replace(gcEvalSource,"maxqueryaa=99","maxqueryaa=99\tmarginscale=0.985");
		compactQueryGcEvaluation(root.resolve("eval_gc"),gcEvalSource);
		reject(()->new FamilyScoringPackage(pkg),"ambiguous");
		reject(()->{
			final Path bad=packageRoot(root.resolve("absolute"),options,policy);
			replace(bad.resolve("short").resolve("families.tsv"),"alphabet=alphabet.config","alphabet=/tmp/a");
			new FamilyScoringPackage(bad.resolve("short").resolve("families.tsv"),
				VerificationMode.STRUCTURAL);
		},"relative");
		reject(()->{
			final Path bad=packageRoot(root.resolve("question"),options,policy);
			replace(bad.resolve("short").resolve("families.tsv"),"blacklist=blacklist.faa.gz","blacklist=?blacklist.faa.gz");
			new FamilyScoringPackage(bad.resolve("short").resolve("families.tsv"),
				VerificationMode.STRUCTURAL);
		},"relative");
		reject(()->{
			final Path bad=packageRoot(root.resolve("traverse"),options,policy);
			replace(bad.resolve("short").resolve("families.tsv"),"matrixresource=POSSUM55.tsv",
				"matrixresource=../POSSUM55.tsv");
			new FamilyScoringPackage(bad.resolve("short").resolve("families.tsv"),
				VerificationMode.STRUCTURAL);
		},"escapes");
		reject(()->{
			final Path bad=packageRoot(root.resolve("view"),options,policy);
			replace(bad.resolve("short").resolve("families.tsv"),"view=short","view=tails");
			new FamilyScoringPackage(bad.resolve("short").resolve("families.tsv"),
				VerificationMode.STRUCTURAL);
		},"view");

		boolean rejected=false;
		try{FamilyTailCascade.loadCallable(options,1);}
		catch(IllegalArgumentException expected){
			rejected=expected.getMessage()!=null &&
				expected.getMessage().contains("compact family references");
		}
		check(rejected,"Legacy callable loading rejects compact-only embedded sidecars");
		callerBridge(root.resolve("caller"),options,policy);

		System.out.println("FAMILY_COMPACT_CALLABLE_TEST_PASS checks="+checks);
	}

	/** The compact caller consumes both views without changing the direct loader's gates. */
	private static void callerBridge(Path root,HashMap<String,String> options,Descriptor policy)throws Exception{
		Files.createDirectory(root);
		final Path pack=callerPackage(root.resolve("valid"),options,policy);
		callerOptions(pack);
		final StopTailCascadeResources resources=StopTailCascadeResources.readStructural(pack.toString(),"unmeasured");
		check(resources.hasShortBins() && resources.queryLengths().equals("40-99,100") &&
			resources.firstStart && resources.allFamilySets,"Compact root retains both approved query views and first-start routing");
		final FamilyTailCascade.LoadedClassifier tails=resources.load(1),shorts=resources.loadShort40(1);
		check(tails.evaluationMarginScale()==1 && shorts.evaluationMarginScale()==.985,
			"Each compact view applies its absolute header margin exactly once");
		check(tails.usesHeuristicRanking() && shorts.usesHeuristicRanking() &&
			tails.usesFamilyGcPenalty() && shorts.usesFamilyGcPenalty(),
			"Caller loading retains embedded length/KD and GC ranking in both views");
		final StopTailCascadeSupport support=resources.loadSupport(1);
		try{
			check(support.maximumTailLength()==100 && support.usesAllFamilySets() && support.usesFirstStart(),
				"Compact resources feed the existing capped100 all-set scanner");
		}finally{support.close();}
		resources.protectOutputs(new String[]{root.resolve("outside.gff").toString()});
		reject(()->resources.protectOutputs(new String[]{pack.resolve("new.gff").toString()}),"inside stop-tail resources");
		final Path alias=root.resolve("blacklist-alias.faa");
		Files.createLink(alias,pack.resolve("blacklist.faa.gz"));
		reject(()->resources.protectOutputs(new String[]{alias.toString()}),"aliases stop-tail resource");
		reject(()->StopTailCascadeResources.readStructural(pack.toString(),"11111111111111111111"),"identity unmeasured");

		final FamilyScoringPackage shortPack=new FamilyScoringPackage(pack.resolve("short/families.tsv"),VerificationMode.STRUCTURAL);
		final Path[] detached=shortPack.compactInputs();
		detached[0]=root.resolve("not-an-input");
		check(shortPack.compactInputs()[0].equals(pack.resolve("short/families.tsv")),
			"Output-protection input discovery cannot mutate the package's retained roster");
		final Path reference=compactReference(root.resolve("short.reference.tsv"),
			pack.resolve("short/families.tsv"),"short");
		new FamilyScoringPackage(reference,VerificationMode.STRUCTURAL).load(1);
		final FamilyScoringPackage stable=new FamilyScoringPackage(reference,VerificationMode.STRUCTURAL);
		Files.write(reference,"extra=field\n".getBytes(StandardCharsets.US_ASCII),StandardOpenOption.APPEND);
		reject(()->stable.load(1),"changed");
		reject(()->new FamilyScoringPackage(compactReference(root.resolve("mismatch.reference.tsv"),
			pack.resolve("short/families.tsv"),"tails"),VerificationMode.STRUCTURAL),"exact compact view");
		reject(()->new FamilyScoringPackage(referenceFile(root.resolve("relative.reference.tsv"),
			"short","short/families.tsv"),VerificationMode.STRUCTURAL),"normalized absolute");
		reject(()->new FamilyScoringPackage(referenceFile(root.resolve("version2.reference.tsv"),
			"#bbgap_compact_classifier_reference\t2","view=short","family=classifier.config"),
			VerificationMode.STRUCTURAL),"Unsupported compact reference version");
		reject(()->new FamilyScoringPackage(compactReference(root.resolve("chain.reference.tsv"),
			reference,"short"),VerificationMode.STRUCTURAL),"exact compact view");
		final Path badMargin=callerPackage(root.resolve("wrong_margin"),options,policy);
		replace(badMargin.resolve("short/families.tsv"),"marginscale=0.985","marginscale=1");
		reject(()->StopTailCascadeResources.readCompact(badMargin),"SHORT40..99/margin0.985");
		final Path badRange=callerPackage(root.resolve("wrong_range"),options,policy);
		replace(badRange.resolve("tails/families.tsv"),"maxqueryaa=100","maxqueryaa=101");
		reject(()->StopTailCascadeResources.readCompact(badRange),"SHORT40..99/margin0.985");
		final String[][] gcChanges={{"familygclambda=24","familygclambda=8"},
			{"gcfilter=true","gcfilter=false"},{"gc_z=8","gc_z=3"},
			{"gc_minimum_members=30","gc_minimum_members=100"},
			{"taillow=40","taillow=41"},{"tailhigh=99","tailhigh=98"}};
		for(int i=0;i<gcChanges.length;i++){
			final Path badGc=callerPackage(root.resolve("wrong_gc_"+i),options,policy);
			replace(badGc.resolve("short/families.tsv"),gcChanges[i][0],gcChanges[i][1]);
			reject(()->StopTailCascadeResources.readCompact(badGc),"Compact caller GC policy requires");
		}
		final Path noGc=callerPackage(root.resolve("missing_gc"),options,policy);
		for(String view:new String[]{"short","tails"}){
			final Path familyFile=noGc.resolve(view).resolve("families.tsv");
			for(String field:new String[]{"familygclambda=24","gcfilter=true","gc_z=8"}){
				replace(familyFile,"\t"+field,"");
			}
		}
		reject(()->StopTailCascadeResources.readCompact(noGc),"Compact caller header requires gcfilter");
		final Path separate=callerPackage(root.resolve("separate_resource"),options,policy);
		Files.copy(separate.resolve("blacklist.faa.gz"),separate.resolve("other_blacklist.faa.gz"));
		replace(separate.resolve("tails/families.tsv"),"blacklist=blacklist.faa.gz","blacklist=other_blacklist.faa.gz");
		reject(()->StopTailCascadeResources.readCompact(separate),"share their alphabet, blacklist and matrix");
		final Path provenance=callerPackage(root.resolve("wrong_provenance"),options,policy);
		replace(provenance.resolve("tails/families.tsv"),"packageprovenance=caller_fixture","packageprovenance=other_fixture");
		reject(()->StopTailCascadeResources.readCompact(provenance),"different package provenance");
		final Path duplicate=callerPackage(root.resolve("duplicate_view"),options,policy);
		Files.copy(duplicate.resolve("short/families.tsv"),duplicate.resolve("short/families.tsv.gz"));
		reject(()->StopTailCascadeResources.readCompact(duplicate),"exactly one families.tsv");
		final Path missing=callerPackage(root.resolve("missing_view"),options,policy);
		Files.delete(missing.resolve("tails/families.tsv"));
		reject(()->StopTailCascadeResources.readCompact(missing),"exactly one families.tsv");
		Files.write(pack.resolve("short/families.tsv"),"\n".getBytes(StandardCharsets.US_ASCII),StandardOpenOption.APPEND);
		reject(resources::verify,"changed");
	}
	/** Exercise the caller's format-dependent default without loading nucleotide models. */
	private static void callerOptions(Path pack)throws Exception{
		final CallGenes.TailCascadeOptions omitted=new CallGenes.TailCascadeOptions();
		final StopTailCascadeResources automatic=omitted.readResources(pack.toString());
		check("unmeasured".equals(automatic.sha80) && automatic.hasShortBins(),
			"Omitted caller verification binds both compact views with structural identity");
		check(!omitted.verificationSpecified() && !omitted.structural() && omitted.referencePin()==null,
			"Format-dependent loading must not rewrite the caller's parsed options");
		final CallGenes.TailCascadeOptions structural=new CallGenes.TailCascadeOptions();
		check(structural.parse("tailverification","structural") &&
			structural.parse("tailreferencesha80","unmeasured"),"Caller consumes explicit compact verification flags");
		final StopTailCascadeResources explicit=structural.readResources(pack.toString());
		check(explicit.sha80.equals(automatic.sha80) && explicit.queryLengths().equals(automatic.queryLengths()),
			"Explicit structural loading preserves the same compact identity and query views");
		final CallGenes.TailCascadeOptions sha80=new CallGenes.TailCascadeOptions();
		sha80.parse("tailverification","sha80");
		reject(()->sha80.readResources(pack.toString()),"do not support tailverification=sha80");
		final CallGenes.TailCascadeOptions pinned=new CallGenes.TailCascadeOptions();
		pinned.parse("tailreferencesha80","11111111111111111111");
		reject(()->pinned.readResources(pack.toString()),"structural identity unmeasured");
		pinned.parse("tailverification","structural");
		reject(()->pinned.readResources(pack.toString()),"structural identity unmeasured");
		final CallGenes.TailCascadeOptions missingPin=new CallGenes.TailCascadeOptions();
		missingPin.parse("tailverification","structural");
		reject(()->missingPin.readResources(pack.toString()),"requires explicit tailreferencesha80");
		final CallGenes.TailCascadeOptions unmeasured=new CallGenes.TailCascadeOptions();
		unmeasured.parse("tailreferencesha80","unmeasured");
		reject(()->unmeasured.readResources(pack.toString()),"requires tailverification=structural");
		reject(()->structural.validate("off"),"require annotation=tailcascade");
	}

	private static Path callerPackage(Path root,HashMap<String,String> options,Descriptor policy)throws Exception{
		final Path pack=packageRoot(root,options,policy),shortFile=pack.resolve("short/families.tsv");
		Files.createDirectory(pack.resolve("tails"));
		final List<String> rows=Files.readAllLines(shortFile,StandardCharsets.US_ASCII);
		check(rows.size()>=FamilyCompactCodec.LINES_PER_RECORD && rows.size()%FamilyCompactCodec.LINES_PER_RECORD==0,
			"Caller fixture preserves complete ten-line header and family blocks");
		final byte[][] headerRows=new byte[FamilyCompactCodec.LINES_PER_RECORD][];
		for(int i=0;i<headerRows.length;i++){headerRows[i]=rows.get(i).getBytes(StandardCharsets.US_ASCII);}
		final FamilyCompactCodec.Header header=FamilyCompactCodec.parseHeader(headerRows,0,new LineParser1('\t'));
		for(String view:new String[]{"short","tails"}){
			final LinkedHashMap<String,String> fields=new LinkedHashMap<String,String>(header.fields);
			fields.put("view",view);
			fields.put("packageprovenance","caller_fixture");
			fields.put("minqueryaa",view.equals("short") ? "40" : "100");
			fields.put("maxqueryaa",view.equals("short") ? "99" : "100");
			fields.put("marginscale",view.equals("short") ? "0.985" : "1");
			fields.put("familygclambda","24");
			fields.put("gcfilter","true");
			fields.put("gc_z","8");
			fields.put("gc_minimum_members","30");
			fields.put("taillow",view.equals("short") ? "40" : "100");
			fields.put("tailhigh",view.equals("short") ? "99" : "100");
			final ByteBuilder out=FamilyCompactCodec.appendHeader(new FamilyCompactCodec.Header(fields),new ByteBuilder());
			for(int i=headerRows.length;i<rows.size();i++){out.append(rows.get(i)).nl();}
			Files.write(pack.resolve(view).resolve("families.tsv"),out.toBytes());
		}
		return pack;
	}

	private static Path compact(Path path,HashMap<String,String> options)throws Exception{
		final LinkedHashMap<String,String> fields=new LinkedHashMap<String,String>();
		fields.put("families","3");
		fields.put("view","short");
		fields.put("alphabet","alphabet.config");
		fields.put("blacklist","blacklist.faa.gz");
		fields.put("matrixresource","POSSUM55.tsv");
		fields.put("live","2");
		fields.put("empty","1");
		fields.put("k","6");
		fields.put("width","512");
		fields.put("seed","0");
		fields.put("ddlreject","6");
		fields.put("ddlaccept","1025");
		fields.put("scalargates","true");
		fields.put("indexearlyreject","false");
		fields.put("namespace","tail");
		fields.put("mode","nonempty");
		fields.put("sketch","positional-v1");
		fields.put("geometry",ddl.ProteinKmerGeometry.make(6,options.get("aabits")).name());
		fields.put("partitions","4");
		fields.put("overlap","2");
		fields.put("tag_bits","2");
		fields.put("segment_words","49");
		fields.put("formulaVersion",FamilyReferenceCodec.FORMULA_VERSION);
		fields.put("rankingfeatures","combined25");
		fields.put("cosine","ranking");
		fields.put("matrix","blosum62");
		fields.put("matrixbinding",ProteinSubstitutionMatrix.BLOSUM62.binding());
		fields.put("top","50");
		fields.put("scalarscale","0.1");
		fields.put("scalarz","3.8906");
		fields.put("identityacceptfloor","30");
		fields.put("blosumacceptfloor","0.13");
		fields.put("blacklistk","6");
		fields.put("blacklistmode","reduced_literal");
		fields.put("blacklistwords",options.get("blacklistwords"));
		fields.put("blacklistprovenance",options.get("blacklistprovenance"));
		fields.put("alignmentmode","pin_right_trim_left_aligned_query");
		fields.put("alignmentminspan","40");
		copy(fields,options,"heuristic");
		copy(fields,options,"minqueryaa");
		copy(fields,options,"maxqueryaa");
		copy(fields,options,"familygclambda");
		copy(fields,options,"gcfilter");
		copy(fields,options,"gc_z");
		for(int i=0;i<SCALARS.length;i++){
			fields.put(SCALARS[i]+"_min","0");
			fields.put(SCALARS[i]+"_max","1");
		}
		final ByteBuilder out=FamilyCompactCodec.appendHeader(new FamilyCompactCodec.Header(fields),
			new ByteBuilder());
		FamilyCompactCodec.appendRecord(record("family_0",0,700,30),out);
		FamilyCompactCodec.appendRecord(record("family_1",1,701,31),out);
		FamilyCompactCodec.appendRecord(new FamilyCompactCodec.Record("family_2",2,702,
			0,0,0,0,0,0,0,0,0,0,FamilyCompactCodec.AlignmentCalibration.UNAVAILABLE,
			new double[10],0,0,0,0,0,new byte[0],
			new char[FamilyCompactCodec.DDL_WIDTH],new float[FamilyCompactCodec.MONOMER_WIDTH]),out);
		Files.write(path,out.toBytes());
		return path;
	}

	private static FamilyCompactCodec.Record record(String id,int rank,int stableId,long members){
		final char[] ddl=new char[FamilyCompactCodec.DDL_WIDTH];
		ddl[rank]=(char)(rank+3);
		final float[] monomers=new float[FamilyCompactCodec.MONOMER_WIDTH];
		monomers[rank]=1;
		final long residues=members*80;
		return new FamilyCompactCodec.Record(id,rank,stableId,members,residues,residues-5,
			0.1f,30f,0f,0.75f,90f,0.13f,0,
			FamilyCompactCodec.AlignmentCalibration.ENABLED,new double[10],80+rank,0.25+rank,
			members,0.5+rank*0.01,members==1 ? 0 : 0.125,
			"ACDEFGHIKLMNPQRSTVWY".getBytes(StandardCharsets.US_ASCII),ddl,monomers);
	}

	private static Path blacklist(Path path,Descriptor policy)throws Exception{
		return blacklist(path,policy,policy.provenance);
	}

	private static Path blacklist(Path path,Descriptor policy,String provenance)throws Exception{
		final String[] groups=groups(policy.aabits);
		final String[] pair=new String[groups.length];
		for(int i=0; i<groups.length; i++){pair[i]=groups[i].split(",")[0];}
		final ByteBuilder out=new ByteBuilder();
		out.append(">bbgap_blacklist_v1;k=6;words=").append(policy.words)
			.append(";mode=").append(policy.mode.label())
			.append(";provenance=").append(encodeHeader(provenance))
			.append(";aabits=").append(encodeHeader(policy.aabits)).nl();
		for(int i=0; i<policy.words; i++){
			if(i>0){out.append(">b_").append(i).nl();}
			final int a=i%groups.length,b=(i/groups.length)%groups.length,c=i/(groups.length*groups.length);
			out.append(pair[a]).append(pair[b]).append(pair[c]).nl();
		}
		if(path.getFileName().toString().endsWith(".gz")){
			try(GZIPOutputStream gz=new GZIPOutputStream(Files.newOutputStream(path))){
				gz.write(out.toBytes());
			}
		}else{Files.write(path,out.toBytes());}
		return path;
	}

	private static HashMap<String,String> options(Path blacklist,Descriptor policy){
		final HashMap<String,String> options=new HashMap<String,String>();
		options.put("top","50");
		options.put("k","6");
		options.put("width","512");
		options.put("seed","0");
		options.put("ddlreject","6");
		options.put("ddlaccept","1025");
		options.put("cosine","ranking");
		options.put("rankingfeatures","combined25");
		options.put("scalarscale","0.1");
		options.put("scalargates","t");
		options.put("scalarz","3.8906");
		options.put("identityacceptfloor","30");
		options.put("blosumacceptfloor","0.13");
		options.put("positionmode","fixed");
		options.put("partitions","4");
		options.put("overlap","2");
		options.put("aabits",policy.aabits);
		options.put("blacklist",blacklist.toString());
		options.put("blacklistk",Integer.toString(policy.filterK));
		options.put("blacklistmode",policy.mode.label());
		options.put("blacklistwords",Integer.toString(policy.words));
		options.put("blacklistprovenance",policy.provenance);
		options.put("indexearlyreject","f");
		options.put("alignmentmode","pin_right_trim_left_aligned_query");
		options.put("alignmentminspan","40");
		options.put("matrix","blosum62");
		options.put("matrixbinding",ProteinSubstitutionMatrix.BLOSUM62.binding());
		return options;
	}

	private static void copy(LinkedHashMap<String,String> fields,HashMap<String,String> options,String key){
		final String value=options.get(key);
		if(value!=null){fields.put(key,value);}
	}

	private static Path packageRoot(Path root,HashMap<String,String> options,Descriptor policy)throws Exception{
		Files.createDirectory(root);
		Files.createDirectory(root.resolve("short"));
		Files.write(root.resolve("alphabet.config"),new ByteBuilder("#bbgap_alphabet_config\t1\n")
			.append("aabits=").append(policy.aabits).nl().toBytes());
		blacklist(root.resolve("blacklist.faa.gz"),policy,"bbgap_blacklist_k6_2048_20261116_v1");
		final ProteinSubstitutionMatrix matrix=possum55(root.resolve("POSSUM55.tsv"));
		final Path compact=compact(root.resolve("short").resolve("families.tsv"),options);
		replace(compact,"matrix=blosum62","matrix=POSSUM55");
		replace(compact,ProteinSubstitutionMatrix.BLOSUM62.binding(),matrix.binding());
		replace(compact,policy.provenance,"bbgap_blacklist_k6_2048_20261116_v1");
		return root;
	}

	private static void rejectMutation(Path root,String relative)throws Exception{
		final FamilyScoringPackage pack=new FamilyScoringPackage(
			root.resolve("short").resolve("families.tsv"),VerificationMode.STRUCTURAL);
		Files.write(root.resolve(relative),"\n".getBytes(StandardCharsets.US_ASCII),
			StandardOpenOption.APPEND);
		reject(()->pack.load(2),"changed");
	}

	private static void compactEvaluation(Path root,Path source)throws Exception{
		Files.createDirectory(root);
		final Path pos=root.resolve("pos.faa"),out=root.resolve("out");
		Files.write(pos,(">q\nACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n")
			.getBytes(StandardCharsets.US_ASCII));
		FamilyEval.main(new String[]{"in="+source,"pos="+pos,"out="+out,
			"verification=structural","t=2"});
		check(Files.exists(out.resolve("classifier_reference.tsv")) &&
			!Files.exists(out.resolve("classifier.config")),
			"Compact evaluation publishes a reloadable reference, not legacy classifier.config");
		final List<String> policy=Files.readAllLines(out.resolve("evaluation_policy.tsv"),StandardCharsets.US_ASCII);
		check(policy.contains("lengthmedians\tembedded") && policy.contains("familyhydro\tembedded"),
			"Compact evaluation records embedded ranking metadata truthfully");
		check(policy.contains("marginscale\t0.985"),
			"Compact evaluation records the inherited header margin");
		final Path explicit=root.resolve("explicit");
		FamilyEval.main(new String[]{"in="+source,"pos="+pos,"out="+explicit,
			"verification=structural","t=1","marginscale=1"});
		check(Files.readAllLines(explicit.resolve("evaluation_policy.tsv"),StandardCharsets.US_ASCII)
			.contains("marginscale\t1.0"),"Explicit compact margins are absolute one-shot overrides");
		final Path reloaded=root.resolve("reloaded");
		FamilyEval.main(new String[]{"in="+out.resolve("classifier_reference.tsv"),
			"pos="+pos,"out="+reloaded,"verification=structural","t=1"});
		check(Files.exists(reloaded.resolve("PASS")),"Compact reference reloads the original family file");
		reject(()->FamilyEval.main(new String[]{"in="+source,"pos="+pos,
			"out="+source.getParent().getParent().resolve("nested/out"),
			"verification=structural","t=1"}),"inside the compact package root");
		reject(()->FamilyEval.main(new String[]{"in="+source,"pos="+pos,
			"out="+root.resolve("mismatch"),"verification=structural","t=1",
			"heuristic=hitfourth_h020_c150","lengthmedians="+pos,
			"lengthmedianssha80=unmeasured","familyhydro="+pos,
			"familyhydrosha80=unmeasured","minqueryaa=40","maxqueryaa=99"}),"embed member medians");
	}

	/** Synthetic gate isolation: 90 percent identity passes 89.85, but not a strict 90 threshold. */
	private static void marginEvaluation(Path root,HashMap<String,String> options,Descriptor policy)throws Exception{
		Files.createDirectory(root);
		final Path pack=packageRoot(root.resolve("package"),options,policy);
		final Path source=pack.resolve("short/families.tsv"),pos=root.resolve("pos.faa");
		final byte[] consensus=("ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY").getBytes(StandardCharsets.US_ASCII);
		final byte[] query=consensus.clone(),substitutions="CDEF".getBytes(StandardCharsets.US_ASCII);
		int identities=0;
		for(int i=0;i<substitutions.length;i++){query[query.length-substitutions.length+i]=substitutions[i];}
		for(int i=0;i<query.length;i++){if(query[i]==consensus[i]){identities++;}}
		check(query.length==40 && identities==36,"Margin oracle needs exactly 36 of 40 identities");
		final DDLPentamerBlacklist filter=DDLPentamerBlacklist.loadReducedStructural(pack.resolve("blacklist.faa.gz"),
			"unmeasured",policy.words,policy.geometry,"bbgap_blacklist_k6_2048_20261116_v1");
		final char[] ddl=new char[FamilyCompactCodec.DDL_WIDTH];
		DDLFamilyProfiles.sketchPositional(query,policy.geometry,0,4,2,ddl,filter);
		int occupied=0;
		for(char value:ddl){if(value!=0){occupied++;}}
		check(occupied>=6 && occupied<1025,"Query-derived synthetic rows admit alignment without a DDL shortcut");
		final List<String> rows=Files.readAllLines(source,StandardCharsets.US_ASCII);
		final byte[][] headerRows=new byte[FamilyCompactCodec.LINES_PER_RECORD][];
		for(int i=0;i<headerRows.length;i++){headerRows[i]=rows.get(i).getBytes(StandardCharsets.US_ASCII);}
		final LinkedHashMap<String,String> fields=new LinkedHashMap<String,String>(
			FamilyCompactCodec.parseHeader(headerRows,0,new LineParser1('\t')).fields);
		fields.put("families","1");
		fields.put("live","1");
		fields.put("empty","0");
		fields.put("marginscale","0.985");
		final float[] monomers=new float[FamilyCompactCodec.MONOMER_WIDTH];
		java.util.Arrays.fill(monomers,0.05f);
		final ByteBuilder output=FamilyCompactCodec.appendHeader(new FamilyCompactCodec.Header(fields),new ByteBuilder());
		FamilyCompactCodec.appendRecord(new FamilyCompactCodec.Record("margin",0,700,30,1200,1050,
			0,80,0,1.1f,90,Float.MAX_VALUE,0,FamilyCompactCodec.AlignmentCalibration.ENABLED,
			new double[10],40,0,30,0.5,0.125,consensus,ddl,monomers),output);
		Files.write(source,output.toBytes());
		Files.write(pos,new ByteBuilder(">margin_query\n").append(query).nl().toBytes());
		marginOracle(source,query);
		for(boolean explicit:new boolean[]{false,true}){
			final Path out=root.resolve(explicit ? "explicit" : "inherited");
			final String[] args={"in="+source,"pos="+pos,"poscount=1","out="+out,
				"verification=structural","t=1","admissionstats=t"};
			if(explicit){
				final String[] scaled=java.util.Arrays.copyOf(args,args.length+1);
				scaled[args.length]="marginscale=1";
				FamilyEval.main(scaled);
			}else{FamilyEval.main(args);}
			final List<String> summary=Files.readAllLines(out.resolve("summary.tsv"),StandardCharsets.US_ASCII);
			check(summary.size()==2 && summary.get(0).startsWith("panel\tclass\tqueries\taccepted\trejected\t"),
				"Margin oracle needs the complete one-panel evaluator summary");
			check(summary.get(1).startsWith(explicit ? "positive\tpositive\t1\t0\t1\t" : "positive\tpositive\t1\t1\t0\t"),
				"Actual compact evaluation must apply the inherited margin and absolute scale-one override");
		}
		evaluationWithoutGc(root,source,pos,options,policy);
	}

	/** The same accepting/rejecting alignment controls must survive an explicit evaluation-only GC override. */
	private static void evaluationWithoutGc(Path root,Path baseline,Path pos,
			HashMap<String,String> options,Descriptor policy)throws Exception{
		final Path source=packageRoot(root.resolve("gc_package"),options,policy).resolve("short/families.tsv");
		Files.write(source,Files.readAllBytes(baseline));
		replace(source,"maxqueryaa=99","maxqueryaa=99\tfamilygclambda=24\tgcfilter=true\tgc_z=8");
		final byte[] original=Files.readAllBytes(source);
		final Path missing=root.resolve("gc_missing");
		final String[] args={"in="+source,"pos="+pos,"poscount=1","out="+missing,
			"verification=structural","t=1","admissionstats=t"};
		reject(()->FamilyEval.main(args),"evaluationgc=off");
		check(!Files.exists(missing),"A missing query-GC binding must fail before creating evaluation output");
		final String[] explicitPackage=java.util.Arrays.copyOf(args,args.length+1);
		explicitPackage[args.length]="evaluationgc=package";
		reject(()->FamilyEval.main(explicitPackage),"GC-ranked evaluation requires querygc");
		for(boolean strict:new boolean[]{false,true}){
			final Path out=root.resolve(strict ? "gc_off_strict" : "gc_off_inherited");
			final String[] off=java.util.Arrays.copyOf(args,args.length+(strict ? 2 : 1));
			off[3]="out="+out;
			off[args.length]="evaluationgc=off";
			if(strict){off[args.length+1]="marginscale=1";}
			FamilyEval.main(off);
			final Path control=root.resolve(strict ? "explicit" : "inherited");
			for(String table:new String[]{"summary.tsv","alignment_counts.tsv","admission_summary.tsv"}){
				check(java.util.Arrays.equals(Files.readAllBytes(control.resolve(table)),Files.readAllBytes(out.resolve(table))),
					"Explicit GC-off must preserve every count in the non-GC alignment control: "+table);
			}
			final List<String> effective=Files.readAllLines(out.resolve("evaluation_policy.tsv"),StandardCharsets.US_ASCII);
			check(effective.contains("evaluation_gc\toff") && effective.contains(
				"family_gc_policy\texplicitly_disabled_for_this_evaluation;package_unchanged"),
				"The effective evaluator policy must record explicit GC removal separately from the package");
			check(!Files.exists(out.resolve("gc_coverage.tsv")),"GC-off evaluation must not invent query-GC coverage");
			check(new FamilyScoringPackage(out.resolve("classifier_reference.tsv"),VerificationMode.STRUCTURAL)
				.load(1).usesFamilyGcPenalty(),"The emitted compact reference must retain the original package GC policy");
		}
		check(java.util.Arrays.equals(original,Files.readAllBytes(source)),
			"Evaluation-only GC removal must not change the compact header or records");
		for(String[] invalid:new String[][]{{"evaluationgc="},{"evaluationgc=auto"},
			{"evaluationgc=off","evaluationgc=off"},{"evaluationgc=package","evaluationgc=off"},
			{"evaluationgc=off","querygc="+root.resolve("unused_querygc.tsv")},
			{"evaluationgc=off","paired=t"},{"evaluationgc=off","floors=unused"}}){
			final String[] bad=java.util.Arrays.copyOf(args,args.length+invalid.length);
			System.arraycopy(invalid,0,bad,args.length,invalid.length);
			reject(()->new FamilyEval(bad),"evaluationgc");
		}
	}

	/** Direct control checks the alignment itself before interpreting evaluator acceptance counts. */
	private static void marginOracle(Path source,byte[] query)throws Exception{
		final FamilyTailCascade.LoadedClassifier base=new FamilyScoringPackage(source,VerificationMode.STRUCTURAL).load(1);
		final FamilyTailCascade.Assignment strict=base.withEvaluationMarginScale(1).newWorker().assign("strict",query);
		final FamilyTailCascade.Assignment relaxed=base.withEvaluationMarginScale(.985).newWorker().assign("relaxed",query);
		check(relaxed.identityMeasured && relaxed.identity==90 &&
			relaxed.querySpan==40 && relaxed.referenceSpan==40 && strict.dpAttempts==1 && relaxed.dpAttempts==1,
			"Margin control must reach 90 percent identity and retain the accepted winner's full 40-aa span");
		check(!strict.identityMeasured && strict.identity==0 && strict.querySpan==0 && strict.referenceSpan==0,
			"Exhausted alignPrepared calls clearAssignment, so rejected queries expose no winner metrics");
		check(!strict.assigned() && relaxed.assigned(),"Strict identity thresholds must distinguish scale one from .985");
	}

	private static void compactQueryGcEvaluation(Path root,Path source)throws Exception{
		Files.createDirectory(root);
		final Path pos=root.resolve("pos.faa"),manifest=root.resolve("querygc.tsv");
		Files.write(pos,(">q\nACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n")
			.getBytes(StandardCharsets.US_ASCII));
		writeNeutralQueryGc(manifest,pos,1);
		FamilyEval.main(new String[]{"in="+source,"pos="+pos,"poscount=1",
			"querygc="+manifest,"out="+root.resolve("out"),"verification=structural","t=2"});
		check(Files.exists(root.resolve("out").resolve("gc_coverage.tsv")),
			"Compact GC evaluation binds the explicit QueryGC manifest");
		reject(()->FamilyEval.main(new String[]{"in="+source,"pos="+pos,"poscount=1",
			"out="+root.resolve("missing"),"verification=structural","t=1"}),
			"GC-ranked evaluation requires querygc");
		final Path mismatch=root.resolve("querygc_mismatch.tsv");
		writeNeutralQueryGc(mismatch,pos,2);
		reject(()->FamilyEval.main(new String[]{"in="+source,"pos="+pos,"poscount=1",
			"querygc="+mismatch,"out="+root.resolve("mismatch"),"verification=structural","t=1"}),
			"panel name/path/count");
	}

	private static void writeNeutralQueryGc(Path path,Path panel,long count)throws Exception{
		Files.write(path,new ByteBuilder("panel\tpath\tqueries\tsource_kind\tquery_gc\tannotation_rows\tsidecar\tcontigs\tgenome_manifest\tgenomes\npositive\t")
			.append(panel.toAbsolutePath().normalize().toString()).tab().append(count)
			.append("\tnone\t-\t0\t-\t-\t-\t0\n").toBytes());
	}

	private static Path compactReference(Path path,Path family,String view)throws Exception{
		return referenceFile(path,"#bbgap_compact_classifier_reference\t1",
			"view="+view,"family="+family.toAbsolutePath().normalize().toString());
	}

	private static Path referenceFile(Path path,String view,String family)throws Exception{
		return referenceFile(path,"#bbgap_compact_classifier_reference\t1","view="+view,"family="+family);
	}

	private static Path referenceFile(Path path,String first,String... rows)throws Exception{
		final ByteBuilder out=new ByteBuilder(first).nl();
		for(String row:rows){out.append(row).nl();}
		Files.write(path,out.toBytes());
		return path;
	}

	private static ProteinSubstitutionMatrix possum55(Path path)throws Exception{
		final String letters="ARNDCEQGHILKMFPSTWYV";
		final ByteBuilder out=new ByteBuilder("#bbgap_substitution_matrix\t1\n");
		out.append("name\t").append(ProteinSubstitutionMatrix.POSSUM55_NAME).nl();
		out.append("description\t").append(ProteinSubstitutionMatrix.POSSUM55_DESCRIPTION).nl();
		out.append("scale\t4\ngap\t16\n");
		for(int i=0;i<letters.length();i++){
			for(int j=0;j<letters.length();j++){
				out.append("score").tab().append(letters.charAt(i)).tab().append(letters.charAt(j))
					.tab().append(i==j ? 20 : -1).nl();
			}
		}
		Files.write(path,out.toBytes());
		return ProteinSubstitutionMatrix.load(ProteinSubstitutionMatrix.POSSUM55_NAME,path.toString());
	}

	private static void replace(Path path,String old,String replacement)throws Exception{
		final String text=new String(Files.readAllBytes(path),StandardCharsets.US_ASCII);
		final String changed=text.replace(old,replacement);
		check(!text.equals(changed),"Replacement must change fixture text");
		Files.write(path,changed.getBytes(StandardCharsets.US_ASCII));
	}

	private static HashMap<String,String> replaced(HashMap<String,String> options,String key,String value){
		final HashMap<String,String> copy=new HashMap<String,String>(options);
		copy.put(key,value);
		return copy;
	}

	private static String readAcceptedAlphabet(Path path)throws Exception{
		final List<String> lines=Files.readAllLines(path,StandardCharsets.US_ASCII);
		check(lines.size()==2 && lines.get(0).equals("#bbgap_alphabet_config\t1"),
			"Accepted compact fixture alphabet must be a two-line versioned resource");
		final String line=lines.get(1);
		check(line.startsWith("aabits="),"Accepted compact fixture alphabet must declare aabits");
		final String aabits=line.substring(7);
		check(aabits.equals(aabits.trim()),"Accepted compact fixture alphabet must be canonical");
		check(!aabits.equals(ProteinBlacklistPolicy.DIMER120_ALPHABET),
			"Accepted compact fixture alphabet must not silently fall back to historical D120");
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(6,aabits);
		check(geometry.name().startsWith("dimer6+1,3,5+"),
			"Accepted compact fixture alphabet must be a K6 dimer projection");
		check(groups(aabits).length==120,"Accepted compact fixture alphabet must contain 120 dimer groups");
		return aabits;
	}

	private static String[] groups(String aabits){
		return aabits.substring(aabits.indexOf(':',aabits.indexOf(':')+1)+1).split("/");
	}

	private static String mergedGroups(String aabits){
		final int groupStart=aabits.indexOf(':',aabits.indexOf(':')+1)+1;
		final String[] groups=groups(aabits);
		final ByteBuilder merged=new ByteBuilder(aabits.substring(0,groupStart));
		for(int i=0; i<groups.length-2; i++){
			if(i>0){merged.append('/');}
			merged.append(groups[i]);
		}
		return merged.append('/').append(groups[groups.length-2]).append(',')
			.append(groups[groups.length-1]).toString();
	}

	private static void reject(Thrower action,String contains)throws Exception{
		boolean rejected=false;
		try{action.run();}
		catch(RuntimeException expected){
			rejected=expected.getMessage()!=null && expected.getMessage().contains(contains);
		}
		check(rejected,"Expected rejection containing "+contains);
	}

	private static String encodeHeader(String value){
		final StringBuilder encoded=new StringBuilder(value.length());
		for(int i=0; i<value.length(); i++){
			final char c=value.charAt(i);
			if(c==';' || c=='%'){encoded.append('%').append(Integer.toHexString(c));}
			else{encoded.append(c);}
		}
		return encoded.toString();
	}

	private static void check(boolean condition,String message){
		checks++;
		if(!condition){throw new AssertionError(message);}
	}

	private static final String[] SCALARS={
		"hydrophobicity","net_charge","aromaticity","fraction_small","pI_tail"
	};
	private static final class Descriptor{
		Descriptor(String aabits_){
			aabits=aabits_;
			geometry=ProteinKmerGeometry.make(6,aabits);
		}

		final String aabits;
		final ProteinKmerGeometry geometry;
		final int filterK=6;
		final int words=2048;
		final DDLPentamerBlacklist.MatchMode mode=DDLPentamerBlacklist.MatchMode.REDUCED_LITERAL;
		final String provenance="q70_compact_fixture";
	}
	private interface Thrower{void run()throws Exception;}
	private static int checks;
}
