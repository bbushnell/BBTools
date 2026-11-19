package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import ddl.DDLFamilyHoldout;
import ddl.DDLFamilySearch;
import dna.AminoAcid;
import fileIO.ByteFile;
import prot.ConsensusScalars;
import prot.FamilyCalibration;
import prot.FamilyAlignmentMetric;
import prot.FamilyCalibrate;
import prot.FamilyEval;
import prot.FamilyReference;
import prot.FamilyReferenceCodec;
import prot.FamilyReferenceLoader;
import prot.FamilyScoringPackage;
import prot.FamilyTailCascade;
import prot.ProteinSubstitutionMatrix;
import structures.ByteBuilder;

/** End-to-end resident preparation and recipe replay with explicit unavailable alignment evidence.
 * @author Brian Bushnell, Keqing
 */
public final class FamilyPreparationMetricTest{

	public static void main(String[] args)throws Exception{
		check(args.length==1,"Supply a fresh metric-preparation fixture directory");
		final Path root=Paths.get(args[0]);check(!Files.exists(root),"Preserve prior preparation evidence");Files.createDirectory(root);
		runWidth(root,1024,null,FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY);
		for(int reject:new int[]{4,5}){
			final Path top100=root.resolve("top100_hit"+reject);Files.createDirectory(top100);
			runWidth(top100,1024,5,null,FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY,null,100,reject,DIVERSE_CORE);
		}
		final Path narrow=root.resolve("width512");Files.createDirectory(narrow);
		runWidth(narrow,512,null,FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY);
		final Path fractional=root.resolve("fractional");Files.createDirectory(fractional);
		runWidth(fractional,1024,"hydro",FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY);
		final Path head=root.resolve("head");Files.createDirectory(head);
		runWidth(head,512,null,FamilyAlignmentMetric.Mode.PIN_LEFT_TRIM_RIGHT_ALIGNED_QUERY);
		for(int k:new int[]{5,6}){
			final Path reduced=root.resolve("reduced"+k);Files.createDirectory(reduced);
			runWidth(reduced,512,k,"all:"+REDUCED_GROUPS,FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY);
			checkReducedBinding(reduced,k);
		}
		for(String matrix:new String[]{"custom40","custom60"}){
			final Path custom=root.resolve(matrix);Files.createDirectory(custom);
			runWidth(custom,512,5,null,FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY,matrix);
		}
		same(root.resolve("one/raw_calibration.tsv.gz"),narrow.resolve("one/raw_calibration.tsv.gz"));
		final Path crossed=narrow.resolve("crossed_profiles.config");
		String recipe=text(narrow.resolve("one/classifier.config"));
		recipe=replaceConfig(recipe,"ddlprofiles",root.resolve("one/ddlprofiles.tsv").toAbsolutePath().toString());
		recipe=replaceConfig(recipe,"ddlprofilessha80",identity(root.resolve("one/ddlprofiles.tsv")));
		write(crossed,recipe);
		reject(()->new FamilyScoringPackage(crossed,FamilyCalibration.VerificationMode.STRUCTURAL).load(1),"A valid512 reference cannot load a1024 replacement profile");
		System.out.println("FAMILY_PREPARATION_METRIC_WIDTH_TEST_PASS geometries=6 fractional=true reduced_k5_k6=true mixed_profiles_rejected=true raw_evidence_parity=true");
		System.out.println("FAMILY_PREPARATION_CUSTOM_MATRIX_TEST_PASS matrices=2 portable=true threads=true updated=true selector_rejection=true");
	}

	/** Every geometry must finish resident calibration, file replay, a genuine gate update and evaluation. */
	private static void runWidth(Path root,int width,String aaBits,FamilyAlignmentMetric.Mode alignmentMode)throws Exception{
		runWidth(root,width,aaBits==null ? 5 : 4,aaBits,alignmentMode);
	}

	private static void runWidth(Path root,int width,int k,String aaBits,FamilyAlignmentMetric.Mode alignmentMode)throws Exception{
		runWidth(root,width,k,aaBits,alignmentMode,null);
	}

	private static void runWidth(Path root,int width,int k,String aaBits,FamilyAlignmentMetric.Mode alignmentMode,String matrixName)throws Exception{
		runWidth(root,width,k,aaBits,alignmentMode,matrixName,50);
	}

	private static void runWidth(Path root,int width,int k,String aaBits,FamilyAlignmentMetric.Mode alignmentMode,String matrixName,int top)throws Exception{
		runWidth(root,width,k,aaBits,alignmentMode,matrixName,top,1,k==6 && aaBits!=null ? DIVERSE_CORE : repeat('A',40));
	}

	private static void runWidth(Path root,int width,int k,String aaBits,FamilyAlignmentMetric.Mode alignmentMode,String matrixName,int top,int reject,String core)throws Exception{
		assert(width==512 || width==1024) : "The fixture covers exactly the validated member/query widths";
		final Path matrixFile=root.resolve("matrix.tsv");
		final boolean customMatrix=matrixName!=null;
		final ProteinSubstitutionMatrix matrix=customMatrix ? writeMatrix(matrixFile,matrixName) : ProteinSubstitutionMatrix.BLOSUM62;
		final Path families=root.resolve("families.tsv"),members=root.resolve("members.tsv"),chunks=root.resolve("chunks.tsv"),
			consensus=root.resolve("consensus.faa"),negative=root.resolve("negative.faa"),admitted=root.resolve("admitted.faa"),panels=root.resolve("panels.tsv"),
			bounds=root.resolve("bounds.tsv"),medians=root.resolve("medians.tsv"),hydro=root.resolve("hydro.tsv"),manifest=root.resolve("prepare.config");
		// The official reduced12 K6 filter excludes AAAAAA; retain a measured-upper witness.
		final double[] scalars=new double[ConsensusScalars.SCALAR_COUNT];
		check(ConsensusScalars.computeScalars(core.getBytes(StandardCharsets.US_ASCII),core.length(),scalars,0)==40,
			"The synthetic member metadata must describe all forty residues");
		final double hydrophobicity=scalars[ConsensusScalars.HYDRO];
		write(families,DDLFamilyHoldout.FAMILY_HEADER+"\n0\tnormal\t1\n1\tshortref\t1\n2\tempty\t0\n");
		write(members,DDLFamilyHoldout.MEMBER_HEADER+"\n0\tnormal\tg|q0\tg\t0\t"+core+"\n1\tshortref\tg|q1\tg\t1\t"+core+"\n");
		write(chunks,DDLFamilyHoldout.CHUNK_HEADER+"\nb00\t0\t"+members+"\t"+identity(members)+"\t2\t80\n");
		write(consensus,">normal\n"+core+"\n>shortref\n"+core.substring(0,20)+"\n>empty\n"+core+"\n");
		write(negative,">negative0\n"+repeat('R',40)+"\n");
		write(admitted,">negative1\n"+core+"\n");
		write(panels,"panel\tpath\tsha80\tqueries\nk5\t"+negative+"\t"+identity(negative)+"\t1\nshadow\t"+admitted+"\t"+identity(admitted)+"\t1\n");
		write(medians,"family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200\n"+
			"normal\t1\t40\t40\t40\t1\nshortref\t1\t40\t40\t40\t1\n");
		write(hydro,"family_rank\tfamily_id\tN\thydrophobicity_mean\n0\tnormal\t1\t"+hydrophobicity+"\n1\tshortref\t1\t"+hydrophobicity+"\n2\tempty\t0\tNA\n");
		final StringBuilder scalar=new StringBuilder("scalar\trows\tmin\tmax\tmin_family\tmax_family\n");
		for(String name:new String[]{"hydrophobicity","net_charge","aromaticity","fraction_small","pI_tail"}){
			scalar.append(name).append("\t79565\t-20\t20\tlo\thi\n");
		}
		write(bounds,scalar.toString());
		final String recipe="schema=family_preparation_v1\nfamilies="+families+"\nfamiliessha80="+identity(families)+
			"\nchunks="+chunks+"\nchunkssha80="+identity(chunks)+"\nconsensus="+consensus+"\nconsensussha80="+identity(consensus)+
			"\nnegativepanels="+panels+"\nnegativepanelssha80="+identity(panels)+"\ncosinepanel=k5\nexpectedfamilies=3\nexpectedmembers=2\n"+
			"buckets=1\nbatches=1\nmaxreferenceaa=1000\nmaxbucketmembers=4\nmaxbucketbytes=1048576\nmaxpreparedbytes=8388608\n"+
			"maxoutputbytes=8388608\nretention=metric_specific\nuppertop=1\nminmemberaa=40\nmaxmemberaa=99\nminnegativeaa=40\nmaxnegativeaa=99\n"+
			"k="+k+"\n"+(aaBits==null ? "" : "aabits="+aaBits+"\n")+"seed=0\nwidth="+width+"\npositional=t\npartitions=4\noverlap=2\nscalarbounds="+bounds+"\nscalarboundssha80="+identity(bounds)+
			"\nscalarscale=0.1\nblosumfloor=0\nscalarz=0\nmonomer=t\ntop="+top+"\nddlreject="+reject+"\nddlaccept=1025\ncosine=ranking\n"+
			"rankingfeatures=combined25\nscalargates=f\nidentityacceptfloor=30\nblosumacceptfloor=0\nindexearlyreject=f\n"+
			"alignmentmode="+alignmentMode.label+"\nalignmentminspan=28\nheuristic=hitfourth_h020_c150\n"+
			"lengthmedians="+medians+"\nlengthmedianssha80="+identity(medians)+"\nfamilyhydro="+hydro+"\nfamilyhydrosha80="+identity(hydro)+
			"\nminqueryaa=40\nmaxqueryaa=99\n"+
			(customMatrix ? "matrix="+matrixName+"\nmatrixfile="+matrixFile+"\nmatrixfilesha80="+identity(matrixFile)+"\n" : "matrix=blosum62\n");
		write(manifest,recipe);
		final Path invalidTop=root.resolve("invalid_top.config");write(invalidTop,replaceConfig(recipe,"top","0"));
		final Path invalidOut=root.resolve("invalid_top");
		reject(()->calibrate(new String[]{"prepare="+invalidTop,"out="+invalidOut,"t=1","rounds=1","targetfpr=1"}),
			"Metric preparation must reject a nonpositive top before running calibration");
		check(!Files.exists(invalidOut.resolve("FROZEN")) && !Files.exists(invalidOut.resolve("references.tsv.gz")),
			"A rejected top must not publish a final calibration package");
		if(width==1024 && k==5 && aaBits==null && !customMatrix && top==50 && reject==1){
			checkRawPairUnion(root, recipe, alignmentMode.label, core);
		}
		final Path one=root.resolve("one"),two=root.resolve("two");
		for(int t:new int[]{1,2}){
			calibrate(new String[]{"prepare="+manifest,"out="+(t==1 ? one : two),"t="+t,"rounds=1","targetfpr=1"});
		}
		if(customMatrix){Files.move(matrixFile,root.resolve("matrix.retained.tsv"));}
		check(!Files.exists(one.resolve("initial")) && Files.isRegularFile(one.resolve("FROZEN")),"Resident preparation writes only its final evaluated package");
		same(one.resolve("references.tsv.gz"),two.resolve("references.tsv.gz"));
		same(one.resolve("raw_calibration.tsv.gz"),two.resolve("raw_calibration.tsv.gz"));
		checkAudit(one.resolve("raw_calibration.tsv.gz"),alignmentMode.label);
		final FamilyReferenceLoader.Loaded loaded=FamilyReferenceLoader.loadMetricBound(one.resolve("references.tsv.gz").toString(),2,3);
		check(loaded.header.k==k,"Serialized references retain nominal k even when no exact prefix remains");
		check(loaded.header.width==width && loaded.all[0].ddl.length==width && loaded.all[1].ddl.length==width,
			"Resident serialization and parallel reload must preserve the declared sketch width");
		check(loaded.header.schema.equals(customMatrix ? FamilyReferenceCodec.MATRIX_MONOMER_SCHEMA : FamilyReferenceCodec.ALIGNMENT_MONOMER_SCHEMA) && loaded.header.alignmentMetric.minimumSpan==28,
			"Final monomer package retains its explicit alignment geometry and span");
		check(loaded.header.alignmentMetric.matrix.same(matrix),"Reference header retains all scores, name, scale and gap");
		check(loaded.all[0].alignmentCalibration==FamilyReference.AlignmentCalibration.ENABLED && loaded.all[0].identityUpper==Math.nextUp(110f) &&
			loaded.all[0].blosumFinal==prot.FamilyGateCodec.roundUpper(1.1),
			"A measured second panel contributes real upper margins after an unavailable first panel");
		check(loaded.all[1].alignmentCalibration==FamilyReference.AlignmentCalibration.UNAVAILABLE && loaded.all[2].members==0,
			"Unavailable live alignment and an empty family remain distinct through serialization");
		final Path config=one.resolve("classifier.config");
		checkConfigTop(config,top);
		final FamilyScoringPackage pack=new FamilyScoringPackage(config,FamilyCalibration.VerificationMode.STRUCTURAL);
		final FamilyTailCascade.LoadedClassifier replay=pack.load(2);
		check(replay.alignmentMatrix().same(matrix),"Package reload retains the selected matrix after removal of its external resource");
		final HashMap<String,String> badScalar=pack.resolvedOptions();
		for(String key:CALLABLE_WRAPPER_OPTIONS){badScalar.remove(key);}
		badScalar.put("scalarboundssha80","not-a-scalar-identity");
		rejectReason(()->FamilyTailCascade.loadCallableStructural(badScalar,1,new DDLFamilySearch.StructuralInputs()),
			"inherited scalar identity","Structural scalar identity must still be explicit");
		if(customMatrix){checkCustomBinding(root,config,negative,matrix);}
		check(replay.hasAlignmentOnlyAcceptance(),"The unchanged1025 shortcut must remain disabled at both sketch widths");
		check(replay.usesHeuristicRanking(),"Ordinary package replay inherits its calibrated common pivot");
		checkMetricGateUpdate(root,config,loaded,admitted,core,top);
		final Path evaluation=root.resolve("evaluation");
		evaluate(new String[]{"in="+config,"pos="+negative,"out="+evaluation,"t=2"});
		final String policy=text(evaluation.resolve("evaluation_policy.tsv"));
		check(policy.contains("top\t"+top+"\n") && policy.contains("heuristic\thitfourth_h020_c150\n") &&
			policy.contains("cosine_weight\t1.5\n") && policy.contains("min_query_aa\t40\n"),
			"Evaluator reports inherited package ranking and top rather than its default flags");
		// Structural mode accepts the matching calibrated heuristic from the config or CLI; off still cannot replace it.
		final Path inherited=root.resolve("matching_config_heuristic");
		evaluate(new String[]{"in="+config,"pos="+negative,"out="+inherited,"t=1"});
		same(evaluation.resolve("summary.tsv"),inherited.resolve("summary.tsv"));
		final Path cli=root.resolve("cli_heuristic");
		evaluate(new String[]{"in="+config,"pos="+negative,"out="+cli,"heuristic=hitfourth_h020_c150","t=1"});
		same(inherited.resolve("summary.tsv"),cli.resolve("summary.tsv"));
		reject(()->evaluate(new String[]{"in="+config,"pos="+negative,"out="+root.resolve("off"),"heuristic=off","t=1"}),
			"Explicit off cannot silently replace the calibrated package's ranking");
		final Path mismatch=root.resolve("mismatch.config");write(mismatch,text(config).replace("alignmentminspan=28","alignmentminspan=70"));
		reject(()->new FamilyScoringPackage(mismatch,FamilyCalibration.VerificationMode.STRUCTURAL).load(1),"Header and recipe span must agree");
		final Path wrongWidth=root.resolve("wrong_width.config");write(wrongWidth,replaceConfig(text(config),"width",Integer.toString(1536-width)));
		reject(()->new FamilyScoringPackage(wrongWidth,FamilyCalibration.VerificationMode.STRUCTURAL).load(1),"Header and recipe width must agree");
		final Path partial=root.resolve("partial.config");write(partial,text(config).replace("heuristic=hitfourth_h020_c150\n",""));
		reject(()->new FamilyScoringPackage(partial,FamilyCalibration.VerificationMode.STRUCTURAL),"A partial heuristic binding cannot silently load legacy ordering");
		System.out.println("FAMILY_PREPARATION_METRIC_TEST_PASS resident=true threads=true availability=true recipe=true evaluation=true verification=structural heuristic_coverage=config_bound");
	}

	private static void checkRawPairUnion(Path root,String recipe,String alignmentMode,String core)throws Exception{
		final Path manifest=root.resolve("raw_pair_union.config"), one=root.resolve("raw_one"), four=root.resolve("raw_four"),
			negative=root.resolve("raw_negative.faa"), admitted=root.resolve("raw_admitted.faa"), panels=root.resolve("raw_panels.tsv");
		write(negative,">raw_k5_0\n"+repeat('R',40)+"\n>raw_k5_1\n"+repeat('R',40)+"\n");
		write(admitted,">raw_shadow_0\n"+core.substring(0,39)+"W\n>raw_shadow_1\n"+
			core.substring(0,39)+"W\n>raw_shadow_2\n"+core.substring(0,38)+"SS\n");
		write(panels,"panel\tpath\tsha80\tqueries\nk5\t"+negative+"\tunmeasured\t2\nshadow\t"+admitted+"\tunmeasured\t3\n");
		recipe=replaceConfig(replaceConfig(recipe,"negativepanels",panels.toString()),"negativepanelssha80","unmeasured");
		write(manifest, recipe+"negativeselection=raw_pair_union\n");
		calibrate(new String[]{"prepare="+manifest,"out="+one,"t=1","rounds=1","targetfpr=1","verification=structural"});
		calibrate(new String[]{"prepare="+manifest,"out="+four,"t=4","rounds=1","targetfpr=1","verification=structural"});
		same(one.resolve("raw_calibration.tsv.gz"), four.resolve("raw_calibration.tsv.gz"));
		same(one.resolve("raw_pair_summary.tsv"), four.resolve("raw_pair_summary.tsv"));
		checkRawPairUnionAudit(one.resolve("raw_calibration.tsv.gz"), alignmentMode);
		checkRawPairSummary(one.resolve("raw_pair_summary.tsv"));
	}

	/** A matching nominal k cannot substitute a different alphabet or silently restore exact words. */
	private static void checkReducedBinding(Path root,int k)throws Exception{
		final String recipe=text(root.resolve("one/classifier.config"));
		check(recipe.contains("k="+k+"\n") && recipe.contains("aabits=all:"+REDUCED_GROUPS+"\n"),
			"Recipe replay retains nominal k and the complete ordered reduced alphabet");
		final String summary=text(root.resolve("one/profiles/summary.tsv"));
		check(summary.contains("reduced"+k+"+"+REDUCED_GROUPS),"Profile metadata identifies the full-word reduced alphabet");
		final Path reordered=root.resolve("reordered_groups.config");
		write(reordered,replaceConfig(recipe,"aabits","all:C/A/DE/FWY/G/IV/LM/NH/P/Q/RK/ST"));
		reject(()->new FamilyScoringPackage(reordered,FamilyCalibration.VerificationMode.STRUCTURAL).load(1),"Changing group order must reject the existing profile binding");
		final Path exact=root.resolve("exact_words.config");write(exact,replaceConfig(recipe,"aabits","none"));
		reject(()->new FamilyScoringPackage(exact,FamilyCalibration.VerificationMode.STRUCTURAL).load(1),"Exact words cannot consume reduced profiles at the same nominal k");
		final Path wrongK=root.resolve("wrong_k.config");write(wrongK,replaceConfig(recipe,"k",Integer.toString(11-k)));
		reject(()->new FamilyScoringPackage(wrongK,FamilyCalibration.VerificationMode.STRUCTURAL).load(1),"Reduced K5 and K6 references cannot share a classifier binding");
	}

	private static void checkMetricGateUpdate(Path root,Path config,FamilyReferenceLoader.Loaded loaded,Path admitted,String core,int top)throws Exception{
		final Path lowered=root.resolve("lowered.references.tsv"),loweredConfig=root.resolve("lowered.config"),out=root.resolve("metric_update");
		final ByteBuilder refs=new ByteBuilder();
		FamilyReferenceCodec.appendHeader(loaded.header,refs);
		for(FamilyReference f:loaded.all){
			final FamilyReference record=f.rank==0 ? new FamilyReference(f.id,f.rank,f.members,f.residues,f.validWords,
				f.consensus,f.ddl,f.spectrum,f.cosineLower,30,0,f.cosineAccept,30,f.blosumFinal,
				f.consensusRaw,f.scalarMask,f.scalarBounds,f.alignmentCalibration) : f;
			FamilyReferenceCodec.appendFiveLines(record,refs,loaded.header);
		}
		write(lowered,refs.toString());
		String recipe=replaceConfig(text(config),"familyreferences",lowered.toAbsolutePath().normalize().toString());
		recipe=replaceConfig(recipe,"familyreferencessha80",identity(lowered));
		final Path external=root.resolve("update.matrix.tsv");
		final boolean custom=!loaded.header.alignmentMetric.matrix.same(ProteinSubstitutionMatrix.BLOSUM62);
		if(custom){
			writeMatrix(external,loaded.header.alignmentMetric.matrix.name);
			recipe=recipe.replace("matrixbinding="+loaded.header.alignmentMetric.matrix.binding()+"\n","")+
				"matrixfile="+external+"\nmatrixfilesha80="+identity(external)+"\n";
		}
		write(loweredConfig,recipe);
		final FamilyTailCascade.LoadedClassifier loweredReplay=new FamilyScoringPackage(loweredConfig,FamilyCalibration.VerificationMode.STRUCTURAL).load(1);
		check(loweredReplay.usesHeuristicRanking(),"Synthetic metric update package must retain the calibrated heuristic before fitting");
		check(loweredReplay.newWorker().assign("accepted",core.getBytes(StandardCharsets.US_ASCII)).assigned(),
			"The synthetic lower gate must admit the measured negative before calibration");
		calibrate(new String[]{"in="+loweredConfig,"neg="+admitted,"out="+out,"t=1","rounds=2","targetfpr=0"});
		checkConfigTop(out.resolve("classifier.config"),top);
		if(custom){
			Files.move(external,root.resolve("update.matrix.retained.tsv"));
			check(!text(out.resolve("classifier.config")).contains("matrixfile="),"Calibration serializes an external resource into a portable score binding");
		}
		final java.util.List<String> rounds=Files.readAllLines(out.resolve("rounds.tsv"));
		check(rounds.size()>=3,"Metric-bound calibration must evaluate the initial and updated gate packages");
		final String[] update=rounds.get(1).split("\t");
		check(Integer.parseInt(update[5])==1,"Metric-bound calibration must update the newly exposed false-accept family");
		final FamilyReferenceLoader.Loaded updated=FamilyReferenceLoader.loadMetricBound(out.resolve("references.tsv.gz").toString(),1,3);
		check(updated.header.alignmentMetric.matrix.same(loaded.header.alignmentMetric.matrix),"Gate updates retain the exact matrix binding");
		check(updated.all[0].identityUpper>100,"Metric-bound calibration must publish a real raised identity gate");
		final FamilyTailCascade.LoadedClassifier replay=new FamilyScoringPackage(out,FamilyCalibration.VerificationMode.STRUCTURAL).load(1);
		check(replay.usesHeuristicRanking(),"Metric gate updates must retain the calibrated common-pivot heuristic");
		check(!replay.newWorker().assign("blocked",core.getBytes(StandardCharsets.US_ASCII)).assigned(),
			"A real metric-bound update must publish enough gate evidence to block the measured false accept");
		check(status(out).equals("converged") && Files.isRegularFile(out.resolve("references.tsv.gz")),
			"Metric-bound calibration must converge by publishing updated references");
	}

	/** Synthetic positive diagonals permit the same normalized oracle with distinct off-diagonal scores. */
	private static ProteinSubstitutionMatrix writeMatrix(Path path,String name)throws Exception{
		final String alphabet="ARNDCEQGHILKMFPSTWYV";
		final int[][] scores=new int[20][20];
		final StringBuilder out=new StringBuilder("#bbgap_substitution_matrix\t1\nname\t"+name+"\nscale\t4\ngap\t16\n");
		for(int i=0;i<20;i++){for(int j=0;j<20;j++){
			final byte a=AminoAcid.acidToNumberExtended[alphabet.charAt(i)],b=AminoAcid.acidToNumberExtended[alphabet.charAt(j)];
			final int s=4*ProteinSubstitutionMatrix.BLOSUM62.score(a,b)-(i==j ? 0 : name.equals("custom40") ? 1 : 2);
			scores[a][b]=s;out.append("score\t").append(alphabet.charAt(i)).append('\t').append(alphabet.charAt(j)).append('\t').append(s).append('\n');
		}}
		write(path,out.toString());return ProteinSubstitutionMatrix.fromScores(name,4,16,scores);
	}

	/** Public selectors assert package identity; they never replace calibrated score units. */
	private static void checkCustomBinding(Path root,Path config,Path queries,ProteinSubstitutionMatrix matrix)throws Exception{
		final String recipe=text(config);
		check(recipe.contains("matrix="+matrix.name+"\n") && recipe.contains("matrixbinding="+matrix.binding()+"\n") &&
			!recipe.contains("matrixfile="),"Prepared package is self-contained with no external matrix dependency");
		final Path missing=root.resolve("missing_matrix.config");
		write(missing,recipe.replace("matrixbinding="+matrix.binding()+"\n",""));
		reject(()->new FamilyScoringPackage(missing,FamilyCalibration.VerificationMode.STRUCTURAL).load(1),"A custom selector cannot silently fall back when scores are missing");
		final Path builtin=root.resolve("builtin_matrix.config");
		write(builtin,recipe.replace("matrix="+matrix.name+"\n","matrix=blosum62\n").replace("matrixbinding="+matrix.binding()+"\n",""));
		reject(()->new FamilyScoringPackage(builtin,FamilyCalibration.VerificationMode.STRUCTURAL).load(1),"A BLOSUM recipe cannot load custom-calibrated reference gates");
		final String[] altered=matrix.binding().split(":");altered[3]=Integer.toString(Integer.parseInt(altered[3])+1);
		final Path changed=root.resolve("changed_matrix.config");write(changed,replaceConfig(recipe,"matrixbinding",String.join(":",altered)));
		reject(()->new FamilyScoringPackage(changed,FamilyCalibration.VerificationMode.STRUCTURAL).load(1),"The same name cannot hide changed matrix scores");
		evaluate(new String[]{"in="+config,"pos="+queries,"out="+root.resolve("matrix_replay"),"matrix="+matrix.name,"admissionstats=t","t=1"});
		reject(()->evaluate(new String[]{"in="+config,"pos="+queries,"out="+root.resolve("wrong_matrix_eval"),"matrix=blosum62","t=1"}),
			"Evaluation cannot override a calibrated package matrix");
		reject(()->calibrate(new String[]{"in="+config,"neg="+queries,"out="+root.resolve("wrong_matrix_calibrate"),"matrix=blosum62","t=1","rounds=1","targetfpr=1"}),
			"Recalibration cannot silently reuse bounds with a different matrix");
	}
	private static void checkConfigTop(Path config,int top)throws Exception{
		check(text(config).contains("\ntop="+top+"\n"),"Serialized classifier must retain the prepared query top");
	}

	private static void checkAudit(Path path,String alignmentMode){
		final String base="roster_rank\tfamily_id\tmembers\tcosine_lower\tcosine_retained\tidentity_lower\tidentity_retained\tblosum_lower\tblosum_retained\tcosine_upper\tidentity_upper\tblosum_upper"+
			"\talignment_mode\talignment_min_span\tmember_alignment_measured\tmember_alignment_failed\tidentity_quota_achievable\tblosum_quota_achievable";
		final String fields="scanned,admitted,selected,eligible,input_too_short,no_alignment,span_too_short";
		final StringBuilder header=new StringBuilder(base);
		for(String panel:new String[]{"k5","shadow"}){for(String field:fields.split(",")){header.append('\t').append(panel).append('_').append(field);}}
		final ByteFile input=ByteFile.makeByteFile1(path.toString(),false);
		try{
			check(new String(input.nextLine(),StandardCharsets.US_ASCII).equals(header.toString()),"Raw audit schema identifies each measured population and failure count");
			for(int rank=0;rank<3;rank++){
				final String[] row=new String(input.nextLine(),StandardCharsets.US_ASCII).split("\t",-1);
				check(row.length==32 && row[0].equals(Integer.toString(rank)),"Audit rows preserve complete physical order and field count");
				check(row[1].equals(rank==0 ? "normal" : rank==1 ? "shortref" : "empty"),"Audit family identity follows the physical roster");
				if(rank==2){for(int i=3;i<row.length;i++){check(row[i].equals("NA"),"Empty families must not acquire empirical audit values");}continue;}
				check(row[12].equals(alignmentMode) && row[13].equals("28"),"Audit metric binding matches fitted measurements");
				check(row[14].equals(rank==0 ? "1" : "0") && row[15].equals(rank==0 ? "0" : "1"),"Member alignment denominator is retained even after failures");
				check(row[16].equals(rank==0 ? "true" : "false") && row[17].equals(row[16]),"All-failed members explicitly make both retention quotas unattainable");
				check(Double.parseDouble(row[3])==1 && row[4].equals("1") && Double.parseDouble(row[9])==0,
					"Raw cosine lower and the separately selected all-query upper remain measured");
				for(int col:new int[]{5,7,10,11}){
					if(rank==1){check(row[col].equals("NA"),"Unavailable alignment bounds are not numeric zeros");}
					else{check(Double.parseDouble(row[col])==(col==5 || col==10 ? 100 : 1),"Measured raw bounds retain original identity and score values");}
				}
				final int[] k5={1,0,0,0,0,0,0},shadow=rank==0 ? new int[]{1,1,1,1,0,0,0} : new int[]{1,1,1,0,1,0,0};
				for(int i=0;i<7;i++){
					check(Integer.parseInt(row[18+i])==k5[i] && Integer.parseInt(row[25+i])==shadow[i],"Each panel retains its literal scanned/admitted/selected/eligible outcomes");
				}
			}
			check(input.nextLine()==null,"Audit has exactly one row per physical family");
		}finally{check(!input.close(),"Audit read must complete successfully");}
	}
	private static void checkRawPairUnionAudit(Path path,String alignmentMode){
		final ByteFile input=ByteFile.makeByteFile1(path.toString(),false);
		try{
			final String[] header=new String(input.nextLine(),StandardCharsets.US_ASCII).split("\t",-1);
			final String[][] rows=new String[3][];
			for(int rank=0;rank<3;rank++){
				rows[rank]=new String(input.nextLine(),StandardCharsets.US_ASCII).split("\t",-1);
				check(rows[rank].length==header.length,"Raw-pair rows preserve the extended header width");
			}
			check(input.nextLine()==null,"Raw-pair audit has exactly one row per physical family");
			final int k5Positive=column(header,"k5_raw_positive"), k5Union=column(header,"k5_raw_union"),
				k5Status=column(header,"k5_coverage_status"),
				shadowPositive=column(header,"shadow_raw_positive"), shadowQuery=column(header,"shadow_raw_query"),
				shadowFamily=column(header,"shadow_raw_family"), shadowOverlap=column(header,"shadow_raw_overlap"),
				shadowUnion=column(header,"shadow_raw_union"),
				shadowStatus=column(header,"shadow_coverage_status"), shadowEligible=column(header,"shadow_eligible"),
				shadowInput=column(header,"shadow_input_too_short"), shadowSpan=column(header,"shadow_span_too_short"),
				identityWitness=column(header,"identity_negative_set"),
				identityId=column(header,"identity_negative_id"), identityColumns=column(header,"identity_alignment_columns"),
				identityPaired=column(header,"identity_paired_columns"), blosumWitness=column(header,"blosum_negative_set"),
				blosumId=column(header,"blosum_negative_id"), blosumColumns=column(header,"blosum_alignment_columns"),
				blosumPaired=column(header,"blosum_paired_columns"), pooledPositive=column(header,"raw_pooled_positive"),
				pooledUnion=column(header,"raw_pooled_union"), pooledUnionShortfall=column(header,"raw_pooled_union_shortfall"),
				pooledMeasured=column(header,"raw_pooled_measured"), pooledMeasuredShortfall=column(header,"raw_pooled_measured_shortfall"),
				pooledStatus=column(header,"raw_pooled_coverage_status");
			check(rows[0][1].equals("normal") && rows[0][12].equals(alignmentMode) && rows[0][13].equals("28"),
				"Raw-pair audit retains the explicit metric binding");
			check(rows[0][k5Positive].equals("0") && rows[0][k5Union].equals("0") &&
				rows[0][k5Status].equals("ZERO_CANDIDATES"), "Raw-pair audit separates panels with no raw positives");
			check(rows[0][shadowPositive].equals("3") && rows[0][shadowQuery].equals("3") && rows[0][shadowFamily].equals("1") &&
				rows[0][shadowOverlap].equals("1") && rows[0][shadowUnion].equals("3") &&
				rows[0][shadowEligible].equals("3") && rows[0][shadowStatus].equals("OK"),
				"Measured raw-pair union rows publish positive, query, family, overlap and union counts");
			check(rows[1][shadowPositive].equals("3") && rows[1][shadowUnion].equals("3") &&
				rows[1][shadowEligible].equals("0") && rows[1][shadowInput].equals("3") && rows[1][shadowSpan].equals("0") &&
				rows[1][shadowStatus].equals("ZERO_MEASURED"),
				"Selected raw positives that fail alignment are distinct from zero-candidate panels");
			check(rows[0][pooledPositive].equals("3") && rows[0][pooledUnion].equals("3") &&
				rows[0][pooledUnionShortfall].equals("0") && rows[0][pooledMeasured].equals("3") &&
				rows[0][pooledMeasuredShortfall].equals("0") && rows[0][pooledStatus].equals("OK") &&
				rows[1][pooledPositive].equals("3") && rows[1][pooledUnion].equals("3") &&
				rows[1][pooledUnionShortfall].equals("0") && rows[1][pooledMeasured].equals("0") &&
				rows[1][pooledMeasuredShortfall].equals("1") && rows[1][pooledStatus].equals("ZERO_MEASURED"),
				"Pooled raw-pair coverage fields publish actual family frontier and measured shortfalls");
			for(int i=3;i<rows[2].length;i++){check(rows[2][i].equals("NA"),"Empty raw-pair rows retain explicit NA evidence");}
			check(rows[0][identityWitness].equals("shadow") && rows[0][identityId].equals("raw_shadow_0") &&
				rows[0][identityColumns].equals("40") && rows[0][identityPaired].equals("40") &&
				rows[0][blosumWitness].equals("shadow") && rows[0][blosumId].equals("raw_shadow_2") &&
				rows[0][blosumColumns].equals("40") && rows[0][blosumPaired].equals("40") &&
				rows[1][identityWitness].equals("NA") && rows[1][blosumWitness].equals("NA"),
				"Raw-pair witnesses retain independent deterministic IDs and paired alignment columns");
		}finally{check(!input.close(),"Raw-pair audit read must complete successfully");}
	}
	private static void checkRawPairSummary(Path path)throws Exception{
		final java.util.List<String> rows=Files.readAllLines(path,StandardCharsets.US_ASCII);
		check(rows.size()==2,"Raw-pair summary must have one typed aggregate row");
		final String[] header=rows.get(0).split("\t",-1), values=rows.get(1).split("\t",-1);
		check(header.length==values.length && integer(header,values,"positive_pairs")==6 &&
			integer(header,values,"query_frontier")==6 && integer(header,values,"family_frontier")==2 &&
			integer(header,values,"frontier_overlap")==2 && integer(header,values,"union_pairs")==6 &&
			integer(header,values,"pooled_queries")==5 && integer(header,values,"query_slots")==10 &&
			integer(header,values,"query_shortfall")==4 && integer(header,values,"family_slots")==2 &&
			integer(header,values,"family_shortfall")==0,
			"Raw-pair summary preserves directional aggregate counts and slot shortfalls");
	}
	private static int column(String[] header,String field){
		for(int i=0;i<header.length;i++){if(header[i].equals(field)){return i;}}
		throw new AssertionError("Missing raw-pair audit column: "+field);
	}
	private static int integer(String[] header,String[] values,String field){
		return Integer.parseInt(values[column(header,field)]);
	}
	private static void same(Path a,Path b){
		final ByteFile x=ByteFile.makeByteFile1(a.toString(),false),y=ByteFile.makeByteFile1(b.toString(),false);
		try{for(byte[] p=x.nextLine(),q=y.nextLine();p!=null || q!=null;p=x.nextLine(),q=y.nextLine()){
			check(Arrays.equals(p,q),"Decoded preparation outputs must match across worker counts");
		}}finally{final boolean xe=x.close(),ye=y.close();check(!xe && !ye,"Parity readers must close without error");}
	}
	private static String repeat(char acid,int n){final char[] out=new char[n];Arrays.fill(out,acid);return new String(out);}
	private static String identity(Path p){
		check(Files.isRegularFile(p) && !Files.isSymbolicLink(p),"Fixture input must be a regular file: "+p);
		return "unmeasured";
	}
	/** Literal populations are declared before the product independently validates its input streams. */
	private static void write(Path p,String text)throws Exception{
		Files.write(p,text.getBytes(StandardCharsets.US_ASCII));
		if(p.toString().endsWith(".faa")){
			long count=0;
			for(String line : text.split("\n")){if(line.startsWith(">")){count++;}}
			populations.put(p.toAbsolutePath().normalize().toString(),count);
		}
	}
	private static long population(String path){
		final Long count=populations.get(Paths.get(path).toAbsolutePath().normalize().toString());
		check(count!=null && count>0,"Fixture panel must have a positive literal population: "+path);
		return count;
	}
	private static void evaluate(String[] args)throws Exception{
		FamilyEval.main(structuralArguments(args,false));
	}
	private static void calibrate(String[] args)throws Exception{
		FamilyCalibrate.main(structuralArguments(args,true));
	}
	private static String[] structuralArguments(String[] args,boolean calibration){
		final ArrayList<String> resolved=new ArrayList<String>(Arrays.asList(args));
		final HashMap<String,String> fields=new HashMap<String,String>();
		for(String arg : args){
			final int equal=arg.indexOf('=');
			if(equal>0){fields.put(arg.substring(0,equal),arg.substring(equal+1));}
		}
		if(!fields.containsKey("verification")){resolved.add("verification=structural");}
		if(fields.containsKey("pos")){
			if(!fields.containsKey("possha80")){resolved.add("possha80=unmeasured");}
			if(!fields.containsKey("poscount")){resolved.add("poscount="+population(fields.get("pos")));}
		}
		if(fields.containsKey("neg")){
			final StringBuilder counts=new StringBuilder(),identities=new StringBuilder();
			for(String panel : fields.get("neg").split(",",-1)){
				if(counts.length()>0){counts.append(',');identities.append(',');}
				counts.append(population(panel));identities.append("unmeasured");
			}
			if(!fields.containsKey("negsha80")){resolved.add("negsha80="+identities);}
			final String key=calibration ? "expectedpanelqueries" : "negcounts";
			if(!fields.containsKey(key)){resolved.add(key+"="+counts);}
		}
		return resolved.toArray(new String[resolved.size()]);
	}
	private static void rejectReason(Action action,String reason,String why)throws Exception{
		try{action.run();}catch(IllegalArgumentException expected){
			check(expected.getMessage()!=null && expected.getMessage().contains(reason),
					why+"; wrong rejection: "+expected.getMessage());
			return;
		}
		throw new AssertionError(why);
	}
	private static final HashMap<String,Long> populations=new HashMap<String,Long>();
	private static final HashSet<String> CALLABLE_WRAPPER_OPTIONS=new HashSet<String>(Arrays.asList(
		"heuristic","lengthmedians","lengthmedianssha80","familyhydro","familyhydrosha80","minqueryaa","maxqueryaa",
		"familygc","familygcsha80","familygcmap","familygcmapsha80","familygcfamilies","familygcmembers","familygclambda"));
	private static String text(Path p)throws Exception{return new String(Files.readAllBytes(p),StandardCharsets.US_ASCII);}
	private static String status(Path out)throws Exception{
		for(String row:Files.readAllLines(out.resolve("calibration.tsv"))){
			if(row.startsWith("status\t")){return row.substring(7);}
		}
		throw new AssertionError("Missing explicit calibration status");
	}
	private static String replaceConfig(String text,String key,String value){
		final StringBuilder out=new StringBuilder();
		for(String line:text.split("\n",-1)){
			if(line.isEmpty()){continue;}
			out.append(line.startsWith(key+"=") ? key+"="+value : line).append('\n');
		}
		return out.toString();
	}
	private static void check(boolean valid,String why){if(!valid){throw new AssertionError(why);}}
	private static void reject(Action action,String why)throws Exception{
		try{action.run();}catch(IllegalArgumentException|IllegalStateException expected){return;}
		throw new AssertionError(why);
	}
	private interface Action{void run()throws Exception;}
	// No R: the R40 negative still has zero dimer cosine. The full core's spectrum norm is exactly9.
	private static final String DIVERSE_CORE="ACDEFGHIKLMNPQSTVWYACDEFGHIKLMNPQSTVWYAC";
	private static final String REDUCED_GROUPS="A/C/DE/FWY/G/IV/LM/NH/P/Q/RK/ST";
}
