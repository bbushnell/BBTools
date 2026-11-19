package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import fileIO.ByteFile;
import ddl.DDLFamilyProfiles;
import ddl.DDLFamilySearch;
import prot.DigestSuffix;
import prot.FamilyEval;
import prot.FamilyReference;
import prot.FamilyReferenceCodec;
import prot.FamilyReferenceExport;
import prot.FamilyTailCascade;
import prot.ProteinSequence;
import prot.ProteinSpectrum;
import structures.ByteBuilder;

/** End-to-end parity with the existing scorer, independent counts, and panel failures.
 * @author Brian Bushnell, Keqing, Fischl
 */
public final class FamilyEvalTest{

	public static void main(String[] args) throws Exception{
		check(args.length==1 && args[0].startsWith("fixture="),"Require fixture=new_absolute_directory");
		final Path root=Paths.get(args[0].substring(8));
		check(root.isAbsolute() && !Files.exists(root),"Fixture directory must be fresh");
		Files.createDirectories(root);
		final String a=repeat("ACDEFGHIKLMNPQRSTVWY",10),b=repeat("WYVTSRQPNMLKIHGFEDCA",10),miss=repeat("A",200);
		final Path refs=root.resolve("references.tsv");
		final ByteBuilder reference=new ByteBuilder();
		FamilyReferenceCodec.appendHeader(new FamilyReferenceCodec.Header(FamilyReferenceCodec.SCHEMA,
			5,1024,0,3,2,1,FamilyReferenceCodec.FORMULA_VERSION,FamilyReferenceCodec.SPECTRUM_ENCODING),reference);
		FamilyReferenceCodec.appendFiveLines(family("f0",0,a),reference);
		FamilyReferenceCodec.appendFiveLines(FamilyReference.empty("empty",1),reference);
		FamilyReferenceCodec.appendFiveLines(family("f2",2,b),reference);
		write(refs,reference.toString());
		final Path config=root.resolve("classifier.config");
		final String policy="top=2\nk=5\nseed=0\nwidth=1024\nddlreject=8\nddlaccept=1025\n"+
			"cosine=ranking\nrankingfeatures=monomer\nscalargates=f\npositional=f\n"+
			"identityacceptfloor=30\nblosumacceptfloor=0\nindexearlyreject=f\n";
		write(config,"familyreferences=references.tsv\nfamilyreferencessha80=unmeasured\n"+policy);
		final Path legacyConfig=root.resolve("legacy.config");
		write(legacyConfig,"familyreferences=references.tsv\n"+policy);
		final Path sha80Config=root.resolve("sha80.config");
		write(sha80Config,"familyreferences=references.tsv\nfamilyreferencessha80="+hash(refs)+"\n"+policy);
		checkObserver(refs,policy,a);
		final Path positive=root.resolve("positive.faa"),negative=root.resolve("negative.faa"),clean=root.resolve("clean.faa");
		final StringBuilder pos=new StringBuilder(),neg=new StringBuilder(),none=new StringBuilder();
		for(int i=0;i<200;i++){fasta(pos,"a"+i,a);fasta(pos,"b"+i,b);fasta(pos,"miss"+i,miss);}
		for(int i=0;i<25;i++){fasta(neg,"false"+i,a);fasta(neg,"miss"+i,miss);}
		for(int i=0;i<40;i++){fasta(none,"clean"+i,miss);}
		write(positive,pos.toString());write(negative,neg.toString());write(clean,none.toString());
		final Path[] panels={positive,negative,clean};
		final int[] queries={600,50,40},accepts={400,25,0};
		final HashMap<String,Long>[] baseline=baselines(root,refs,panels,queries,policy);
		for(int i=0;i<panels.length;i++){
			long total=0;for(long n:baseline[i].values()){total+=n;}
			check(total==accepts[i],"Existing scorer must match the independently constructed fixture counts");
		}
		for(int t:new int[]{1,4}){
			final Path out=root.resolve("threads"+t);
			evaluate(arguments(config,positive,negative,clean,out,t));
			final List<String> summary=Files.readAllLines(out.resolve("summary.tsv"),StandardCharsets.US_ASCII);
			check(summary.size()==4,"Exactly three panel rows are required");
			for(int i=0;i<panels.length;i++){
				final String[] row=summary.get(i+1).split("\t",-1);
				check(Long.parseLong(row[2])==queries[i] && Long.parseLong(row[3])==accepts[i] &&
					Long.parseLong(row[4])==queries[i]-accepts[i],"Every panel needs its own observed denominator and exact outcome counts");
				check(!summary.get(i+1).contains("NaN") && !summary.get(i+1).contains("Infinity"),"Metrics must remain finite or explicitly unavailable");
			}
			check(summary.get(1).split("\t")[7].equals("NA"),"Positive rows do not invent negative FPR");
			check(Double.parseDouble(summary.get(2).split("\t")[8])==500000,"False accepts per million uses the negative panel denominator");
			final List<String> families=lines(out.resolve("family_false_accepts.tsv.gz"));
			check(families.size()==5,"Two live families, including zero-hit families, appear for each negative panel");
			for(int p=1;p<panels.length;p++){
				for(int f=0;f<2;f++){
					final String[] row=families.get(1+(p-1)*2+f).split("\t");
					final String id=f==0 ? "f0" : "f2";
					check(row[2].equals(id) && Integer.parseInt(row[1])==(f==0 ? 0 : 2),"Live ordinals must not replace original roster ranks");
					check(Long.parseLong(row[4])==baseline[p].getOrDefault(id,0L),"Per-family false accepts must match the existing scorer");
				}
			}
			check(Files.exists(out.resolve("PASS")),"Only a complete evaluation publishes PASS");
		}
		check(Files.readAllLines(root.resolve("threads1/summary.tsv")).equals(Files.readAllLines(root.resolve("threads4/summary.tsv"))),
			"Thread count must not change any panel statistic");
		check(lines(root.resolve("threads1/family_false_accepts.tsv.gz")).equals(lines(root.resolve("threads4/family_false_accepts.tsv.gz"))),
			"Thread count must not change the per-family table");
		final Path rawDefault=root.resolve("raw_default");
		FamilyEval.main(new String[]{"in="+legacyConfig,"pos="+positive,"neg="+negative+","+clean,"negnames=false,clean",
			"out="+rawDefault,"t=1"});
		check(Files.readAllLines(rawDefault.resolve("inputs.tsv"),StandardCharsets.US_ASCII).equals(Arrays.asList(
			"panel\tpath\texpected_sha80","package\t"+legacyConfig+"\tunmeasured","positive\t"+positive+"\tNA",
			"false\t"+negative+"\tNA","clean\t"+clean+"\tNA")) &&
			Files.readAllLines(rawDefault.resolve("input_populations.tsv"),StandardCharsets.US_ASCII).equals(Arrays.asList(
				"panel\texpected_queries","positive\tNA","false\tNA","clean\tNA")),
			"Raw default structural evaluation must accept omitted panel pins/populations as exact NA provenance");
		final Path sha80=root.resolve("sha80");
		FamilyEval.main(new String[]{"in="+sha80Config,"pos="+positive,"possha80="+hash(positive),
			"neg="+negative+","+clean,"negnames=false,clean","negsha80="+hash(negative)+","+hash(clean),
			"out="+sha80,"t=1","verification=sha80"});
		check(lines(sha80.resolve("summary.tsv")).equals(lines(rawDefault.resolve("summary.tsv"))),
			"Explicit SHA80 evaluation must retain the same ordinary panel counts");
		final Path timed=root.resolve("stage_cpu");
		final ArrayList<String> timedArgs=new ArrayList<String>(Arrays.asList(arguments(config,positive,negative,clean,timed,4)));
		timedArgs.add("stagecpu=t");
		evaluate(timedArgs.toArray(new String[timedArgs.size()]));
		check(lines(timed.resolve("summary.tsv")).equals(lines(root.resolve("threads4/summary.tsv"))),
			"Stage CPU timing must preserve panel results");
		final List<String> stage=lines(timed.resolve("stage_cpu.tsv"));
		check(stage.size()==6 && stage.get(0).equals("stage\tnanos\tseconds\tpercent_worker_cpu\tpercent_search_cpu"),
			"Stage CPU report needs the expected header and four-part partition");
		long parts=0,search=0,total=-1;
		for(int i=1;i<stage.size();i++){
			final String[] row=stage.get(i).split("\t");
			check(row.length==5 && Long.parseLong(row[1])>=0 && Double.parseDouble(row[2])>=0,
				"Stage CPU rows must contain nonnegative times");
			if(row[0].equals("total")){total=Long.parseLong(row[1]);}
			else{
				parts+=Long.parseLong(row[1]);
				if(!row[0].equals("other_worker")){search+=Long.parseLong(row[1]);}
			}
		}
		check(total>=parts && total==parts,"Stage CPU rows must partition total worker time");
		check(search>0,"Stage CPU timing must record resident index, shortlist and alignment worker time");
		final Path explicit=root.resolve("scale_one");
		final String[] defaultArgs=arguments(config,positive,negative,clean,explicit,4);
		final String[] scaledArgs=Arrays.copyOf(defaultArgs,defaultArgs.length+1);scaledArgs[defaultArgs.length]="marginscale=1";
		evaluate(scaledArgs);
		for(String name:new String[]{"summary.tsv","family_false_accepts.tsv.gz","classifier.config"}){
			check(lines(root.resolve("threads1").resolve(name)).equals(lines(explicit.resolve(name))),"Explicit scale1 preserves default scientific bytes: "+name);
		}
		final Path oldZero=root.resolve("margin_zero"), memberOne=root.resolve("member_m1");
		final String[] zeroBase=arguments(config,positive,negative,clean,oldZero,4);
		final String[] memberBase=arguments(config,positive,negative,clean,memberOne,4);
		final String[] zeroArgs=Arrays.copyOf(zeroBase,zeroBase.length+1);zeroArgs[zeroBase.length]="marginscale=0";
		final String[] memberOneArgs=Arrays.copyOf(memberBase,memberBase.length+1);memberOneArgs[memberBase.length]="membermultiplier=1";
		evaluate(zeroArgs);evaluate(memberOneArgs);
		for(String name:new String[]{"summary.tsv","family_false_accepts.tsv.gz","classifier.config"}){
			check(lines(oldZero.resolve(name)).equals(lines(memberOne.resolve(name))),"M1 matches old margin0 for matched stored lowers: "+name);
		}
		for(String text:new String[]{"0.9","1","1.75","0"}){
			final Path out=root.resolve("member_"+text.replace('.','_'));
			evaluate(new String[]{"in="+config,"pos="+positive,"neg="+negative,"out="+out,"membermultiplier="+text,"t=4"});
			final double multiplier=Double.parseDouble(text);
			final List<String> memberPolicy=lines(out.resolve("evaluation_policy.tsv"));
			check(memberPolicy.contains("member_multiplier_mode\ton") &&
				memberPolicy.contains("membermultiplier\t"+Double.toString(multiplier)) &&
				memberPolicy.contains("formula\tmax(runtime_floor,membermultiplier*stored_member_lower)") &&
				memberPolicy.contains("identity_runtime_floor\t30.0") &&
				memberPolicy.contains("blosum_runtime_floor\t0.0") &&
				memberPolicy.contains("identity_lower_rejection\tunchanged_stored_member_lower") &&
				memberPolicy.contains("negative_derived_upper_contribution\tnone") &&
				memberPolicy.contains("tie_definition\tfloor_binding_when_runtime_floor_ge_scaled_member_lower") &&
				memberPolicy.contains("zero_member_families\t0") &&
				memberPolicy.contains("alignment_unavailable_families\t0"),"Member multiplier policy records the exact replay schema for M="+text);
			if(multiplier==0){
				check(memberPolicy.contains("identity_floor_binding_families\t2") &&
					memberPolicy.contains("blosum_floor_binding_families\t2"),"M0 ties both families to both runtime floors");
			}
		}
		reject(new String[]{"in="+config,"neg="+negative,"out="+root.resolve("mixed_margin"),
			"marginscale=0","membermultiplier=1"},"Old and new replay modes must be explicit alternatives");
		for(String bad:new String[]{"-1","NaN","Infinity"}){
			reject(new String[]{"in="+config,"neg="+negative,"out="+root.resolve("bad_member"),
				"membermultiplier="+bad},"Invalid member multiplier rejected");
		}
		for(boolean positiveOnly:new boolean[]{true,false}){
			final Path out=root.resolve(positiveOnly ? "positive_only" : "shadow_only");
			evaluate(new String[]{"in="+config,(positiveOnly ? "pos=" : "neg=")+(positiveOnly ? positive : negative),
				"out="+out,"marginscale=2","t=4"});
			final List<String> summary=lines(out.resolve("summary.tsv"));
			final String[] row=summary.get(1).split("\t");
			check(summary.size()==2 && Integer.parseInt(row[2])==(positiveOnly ? 600 : 50) && Long.parseLong(row[3])==0,
				"Single-panel evaluation uses the scaled gates without a fabricated other panel");
			check(lines(out.resolve("evaluation_policy.tsv")).contains("marginscale\t2.0"),"Actual evaluation policy is recorded separately from frozen references");
		}
		for(String bad:new String[]{"-1","NaN","Infinity"}){
			reject(new String[]{"in="+config,"neg="+negative,"out="+root.resolve("bad_scale"),"marginscale="+bad},"Invalid margin scale rejected");
		}
		reject(new String[]{"in="+config,"out="+root.resolve("no_panels")},"At least one actual panel is required");
		reject(new String[]{"in="+config,"neg="+negative,"possha80=00000000000000000000","out="+root.resolve("dangling_pos")},"A positive pin requires a positive panel");
		reject(new String[]{"in="+config,"pos="+positive,"negnames=shadow","out="+root.resolve("dangling_neg")},"Negative labels require negative panels");
		final Path empty=root.resolve("empty.faa");write(empty,"");
		reject(arguments(config,positive,empty,clean,root.resolve("empty_result"),4),"Empty panel must not produce a valid rate");
		check(!Files.exists(root.resolve("empty_result/PASS")),"Empty panel failure cannot publish completion");
		final Path invalid=root.resolve("invalid.faa");write(invalid,">bad\nACDE*FGHI\n");
		reject(arguments(config,positive,invalid,clean,root.resolve("invalid_result"),4),"Worker input failure must propagate without a hang");
		check(!Files.exists(root.resolve("invalid_result/PASS")),"Worker failure cannot publish completion");
		final String[] duplicate=arguments(config,positive,negative,clean,root.resolve("duplicate_result"),1);
		for(int i=0;i<duplicate.length;i++){if(duplicate[i].startsWith("negnames=")){duplicate[i]="negnames=same,same";}}
		reject(duplicate,"Duplicate panel labels must be rejected");
		final Path badConfig=root.resolve("bad.config");
		write(badConfig,new String(Files.readAllBytes(config),StandardCharsets.US_ASCII).replace("ddlreject=8","ddlreject=8\nbogus=1"));
		reject(arguments(badConfig,positive,negative,clean,root.resolve("bad_result"),1),"Unknown package options must not be silently ignored");
		checkLengthBins(root, config, a);
		checkCapAndQueryBounds(root,config,a);
		checkBenchmarkRepeats(root,config,positive,negative,clean);
		System.out.println("FAMILY_EVAL_TEST_PASS checks="+checks+" verification=structural baseline=direct_callable");
	}

	/** Repetition adds timing work, never observations to the scientific panel denominators. */
	private static void checkBenchmarkRepeats(Path root,Path config,Path positive,Path negative,Path clean) throws Exception{
		assert(Files.exists(root.resolve("threads4/PASS"))) : "Repeated-work fixtures compare against the completed independent ordinary evaluation";
		check(!Files.exists(root.resolve("threads4/benchmark_time.tsv")),"Default evaluation must not enable benchmarking");
		for(int threads:new int[]{1,4}){
			for(int repeats:new int[]{1,3}){
				final Path out=root.resolve("benchmark_t"+threads+"_r"+repeats);
				final ArrayList<String> args=new ArrayList<String>(Arrays.asList(arguments(config,positive,negative,clean,out,threads)));
				args.add("benchmarkrepeats="+repeats);
				evaluate(args.toArray(new String[args.size()]));
				for(String file:new String[]{"summary.tsv","family_false_accepts.tsv.gz","classifier.config",
						"evaluation_policy.tsv","inputs.tsv","interpretation.txt","PASS"}){
					check(lines(root.resolve("threads4").resolve(file)).equals(lines(out.resolve(file))),
						"Repeated work must preserve ordinary scientific bytes and original denominators: "+file);
				}
				final HashMap<String,String> scoring=timingTable(out.resolve("scoring_time.tsv"));
				check(scoring.get("queries").equals("690") && scoring.get("worker_threads").equals(Integer.toString(threads)),
					"Ordinary scoring timing still describes one complete690-query pass");
				if(repeats==1){
					check(!Files.exists(out.resolve("benchmark_time.tsv")),"Explicit one pass keeps repeated-work timing disabled");
					continue;
				}
				final HashMap<String,String> bench=timingTable(out.resolve("benchmark_time.tsv"));
				check(bench.get("passes").equals("3") && bench.get("unique_queries").equals("690") &&
					bench.get("query_evaluations").equals("2070"),"Three full passes mean2070 evaluations of690 unique queries");
				check(bench.get("classifier_loads").equals("1") && bench.get("scientific_counts").equals("first_pass_only") &&
					bench.get("per_pass_count_parity").equals("exact") &&
					bench.get("repeated_queries_are_independent_samples").equals("false"),
					"The benchmark must expose reuse/parity and cannot imply extra independent observations");
				check(Double.parseDouble(bench.get("wall_seconds"))>=Double.parseDouble(scoring.get("wall_seconds")) &&
					Double.parseDouble(bench.get("process_cpu_seconds"))>=Double.parseDouble(scoring.get("process_cpu_seconds")),
					"Repeated-work timing must include the entire first scoring interval");
			}
		}
		for(String bad:new String[]{"0","-1","NaN","2147483648","inf"}){
			reject(new String[]{"in="+config,"pos="+positive,"out="+root.resolve("bad_repeats"),"benchmarkrepeats="+bad},
				"Invalid repetition counts must fail before scoring");
		}
		reject(new String[]{"in="+config,"pos="+positive,"out="+root.resolve("paired_repeats"),"paired=t","benchmarkrepeats=3"},
			"Ordinary repeated-work timing must not silently enter the separate paired-ranking path");
	}

	/** Parse the small explicit timing record, rejecting duplicate keys in fixture output. */
	private static HashMap<String,String> timingTable(Path path){
		assert(path!=null) : "A completed evaluator output path is required for timing assertions";
		final HashMap<String,String> result=new HashMap<String,String>();
		for(String line:lines(path)){
			final String[] fields=line.split("\t",-1);
			check(fields.length==2 && result.put(fields[0],fields[1])==null,"Timing records need unique two-column rows");
		}
		return result;
	}

	/** Explicit cap agrees with a separately written recipe; bounds validate normalized input without filtering. */
	private static void checkCapAndQueryBounds(Path root,Path config,String sequence) throws Exception{
		final byte[] original=Files.readAllBytes(config);
		final Path topOne=root.resolve("top_one.config"),panel=root.resolve("bounded.faa");
		write(topOne,new String(original,StandardCharsets.US_ASCII).replace("top=2\n","top=1\n"));
		write(panel,">full\n"+sequence+"*\n>short\n"+sequence.substring(0,20)+"\n");
		final Path explicit=root.resolve("cap_explicit"),recipe=root.resolve("cap_recipe");
		evaluate(new String[]{"in="+config,"pos="+panel,"out="+explicit,"top=1","heuristic=off","minqueryaa=1","maxqueryaa=200","t=2"});
		evaluate(new String[]{"in="+topOne,"pos="+panel,"out="+recipe,"t=2"});
		check(lines(explicit.resolve("summary.tsv")).equals(lines(recipe.resolve("summary.tsv"))),
			"Cap override must agree with independently changing only the recipe cap");
		check(lines(explicit.resolve("classifier.config")).contains("top=1") && Arrays.equals(original,Files.readAllBytes(config)),
			"Resolved output records the effective cap while preserving the parent recipe bytes");
		check(lines(explicit.resolve("summary.tsv")).get(1).split("\t")[2].equals("2"),
			"Both short input and200 residues plus terminal stop must be counted, never filtered or cropped");
		final List<String> policy=lines(explicit.resolve("evaluation_policy.tsv"));
		check(policy.contains("top\t1") && policy.contains("min_query_aa\t1") && policy.contains("max_query_aa\t200") &&
			policy.contains("query_length_policy\tvalidate_only_no_filter_no_truncate"),"Effective cap and validation bounds must be explicit");
		final HashMap<String,String> timing=new HashMap<String,String>();
		for(String line:lines(explicit.resolve("scoring_time.tsv"))){final String[] parts=line.split("\t");timing.put(parts[0],parts[1]);}
		check(timing.get("queries").equals("2") && timing.get("worker_threads").equals("2") &&
			Double.parseDouble(timing.get("wall_seconds"))>=0 && Double.parseDouble(timing.get("process_cpu_seconds"))>=0,
			"Scoring interval records nonnegative wall/processCPU and the actual full query/worker counts");
		final Path tooLong=root.resolve("too_long.faa");write(tooLong,">long\n"+sequence+"A\n");
		reject(new String[]{"in="+config,"pos="+tooLong,"out="+root.resolve("too_long_result"),"minqueryaa=1","maxqueryaa=200","t=2"},
			"An over-cap protein must fail instead of being cropped or skipped");
		check(!Files.exists(root.resolve("too_long_result/PASS")),"Out-of-range input cannot publish PASS");
		reject(new String[]{"in="+config,"pos="+panel,"out="+root.resolve("too_short_result"),"minqueryaa=40","maxqueryaa=99","t=2"},
			"Bounded evaluation cannot discard records outside its requested interval");
		for(String[] extra:new String[][]{{"top=0"},{"minqueryaa=1"},{"maxqueryaa=200"},{"minqueryaa=0","maxqueryaa=200"},
			{"minqueryaa=100","maxqueryaa=99"},{"minqueryaa=1","maxqueryaa=201"},{"heuristic=hitfourth_h020_c150"},
			{"paired=t","top=100"}}){
			final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+config,"pos="+panel,"out="+root.resolve("invalid_override")));
			args.addAll(Arrays.asList(extra));reject(args.toArray(new String[args.size()]),"Incomplete or conflicting evaluation overrides must fail");
		}
	}

	/** Compare one-pass length reports with independent complete evaluations of literal length subsets. */
	private static void checkLengthBins(Path root, Path config, String sequence) throws Exception{
		assert(sequence.length()==200) : "The literal boundary fixture uses prefixes of the existing200-aa reference";
		final int[] lengths={40,59,60,79,80,99,100,199,200,200};
		final int[] bins={1,1,2,2,3,3,4,6,7,7};
		final long[] queries={0,2,2,2,1,0,1,2}, residues={0,99,139,179,100,0,199,400};
		final StringBuilder all=new StringBuilder();
		final StringBuilder[] subsets=new StringBuilder[queries.length];
		for(int b=0; b<subsets.length; b++){subsets[b]=new StringBuilder();}
		for(int i=0; i<lengths.length; i++){
			final String peptide=sequence.substring(0, lengths[i])+(i==1 || i==9 ? "*" : "");
			fasta(all, "length"+i, peptide);fasta(subsets[bins[i]], "length"+i, peptide);
		}
		final Path positive=root.resolve("length_positive.faa"), negative=root.resolve("length_negative.faa");
		write(positive, all.toString());write(negative, all.toString());
		final Path baseline=root.resolve("length_baseline");
		evaluate(new String[]{"in="+config, "pos="+positive, "neg="+negative, "out="+baseline, "t=1"});
		check(!Files.exists(baseline.resolve("length_summary.tsv")), "Absent lengthbins leaves the default output set unchanged");
		final String[][] independent=new String[queries.length][];
		for(int b=0; b<subsets.length; b++){
			if(queries[b]==0){continue;}
			final Path input=root.resolve("length_subset"+b+".faa"), output=root.resolve("length_subset"+b);
			write(input, subsets[b].toString());
			evaluate(new String[]{"in="+config, "pos="+input, "out="+output, "t=1"});
			independent[b]=lines(output.resolve("summary.tsv")).get(1).split("\t");
		}
		for(int threads:new int[]{1,4}){
			final Path output=root.resolve("length_threads"+threads);
			evaluate(new String[]{"in="+config, "pos="+positive, "neg="+negative, "out="+output,
				"lengthbins=40,60,80,100,120,150,200", "t="+threads});
			for(String file:new String[]{"summary.tsv","family_false_accepts.tsv.gz","classifier.config","evaluation_policy.tsv"}){
				check(lines(baseline.resolve(file)).equals(lines(output.resolve(file))), "Length reporting preserves ordinary scientific bytes: "+file);
			}
			final List<String> table=lines(output.resolve("length_summary.tsv"));
			check(table.size()==17, "Two panels retain all eight bins, including empty and outside bins");
			final String[] lower={"0","40","60","80","100","120","150","200"};
			final String[] upper={"40","60","80","100","120","150","200","NA"};
			for(int p=0; p<2; p++){
				for(int b=0; b<queries.length; b++){
					final String[] row=table.get(1+p*queries.length+b).split("\t");
					check(row.length==20 && row[2].equals(lower[b]) && row[3].equals(upper[b]), "Length interval boundaries remain explicit and ordered");
					check(Long.parseLong(row[4])==queries[b] && Long.parseLong(row[7])==residues[b],
						"Every boundary query and normalized terminal-stop length belongs to its literal bin");
					if(queries[b]==0){
						for(int c=8; c<12; c++){check(row[c].equals("NA"), "An empty bin has unavailable rates, never NaN or zero-error evidence");}
						check(row[5].equals("0") && row[6].equals("0"), "Empty bins cannot invent outcomes");
					}else{
						check(row[5].equals(independent[b][3]) && row[6].equals(independent[b][4]),
							"One-pass length outcomes equal independent complete subset evaluation");
						for(int s=0; s<8; s++){check(row[s+12].equals(independent[b][s+9]), "One-pass terminal stages equal independent subset evaluation");}
						check(p==0 ? row[10].equals("NA") : row[8].equals("NA"), "Length reporting preserves positive/negative metric semantics");
					}
				}
			}
		}
		check(lines(root.resolve("length_threads1/length_summary.tsv")).equals(lines(root.resolve("length_threads4/length_summary.tsv"))),
			"One and four workers must emit identical length-bin tables");
		for(String bad:new String[]{"","0,60","60,60","80,40","-1","NaN","2147483648","40,"}){
			reject(new String[]{"in="+config, "pos="+positive, "out="+root.resolve("bad_length_bins"), "lengthbins="+bad},
				"Malformed or nonincreasing length boundaries must fail before scoring");
		}
		reject(new String[]{"in="+config, "pos="+positive, "out="+root.resolve("paired_length_bins"), "paired=t", "lengthbins=40,60"},
			"An ordinary length report must not silently ignore paired ranking mode");
	}

	private static String[] arguments(Path config,Path pos,Path neg,Path clean,Path out,int threads){
		return new String[]{"in="+config,"pos="+pos,"neg="+neg+","+clean,"negnames=challenge,clean","out="+out,"t="+threads};
	}
	private static String hash(Path p){return DigestSuffix.file(p.toString());}
	private static void checkObserver(Path refs,String policy,String sequence) throws Exception{
		final DDLFamilySearch.StructuralInputs inputs=new DDLFamilySearch.StructuralInputs();
		final FamilyTailCascade.LoadedClassifier classifier=loadFixture(refs,policy,2,inputs);
		check(classifier.familyCount()==2 && classifier.familyRankAt(1)==2 && classifier.familyIdAt(1).equals("f2") &&
			classifier.familyMembersAt(1)==1,"Public roster metadata must expose live ordinals without losing original ranks");
		final FamilyTailCascade.Worker worker=classifier.newWorker();
		final byte[] raw=(sequence+"*").getBytes(StandardCharsets.US_ASCII);
		final FamilyTailCascade.Assignment expected=worker.assign("observer",raw);
		final int[] calls={0};
		worker.assignObserve("observer",raw,(length,terminal,family,ddl,identity,blosum)->{
			check(length==sequence.length(),"Observer length follows the existing terminal-stop normalization");
			check(terminal.equals(expected.terminal) && ddl==expected.ddlMatches &&
				Float.compare(identity,expected.identity)==0 && Float.compare(blosum,expected.blosum)==0,
				"Observed terminal metrics must exactly match the existing immutable assignment API");
			check(family>=0 && classifier.familyRankAt(family)==expected.rosterRank &&
				classifier.familyIdAt(family).equals(expected.familyId),"Observer ordinal must name the same winning family");
			calls[0]++;
		});
		check(calls[0]==1,"Each query must invoke its observer exactly once");
		checkMargins(classifier,sequence);
		inputs.verify();
	}
	private static void checkMargins(FamilyTailCascade.LoadedClassifier classifier,String sequence){
		check(!classifier.evaluationMemberMultiplierMode() &&
			!classifier.withEvaluationMarginScale(0).evaluationMemberMultiplierMode(),
			"Default and margin-only classifiers do not define a member multiplier");
		check(classifier.withEvaluationMarginScale(1)==classifier,"Scale1 must take the original classifier path exactly");
		check(FamilyTailCascade.evaluationGate(30,90,40,0)==40 && FamilyTailCascade.evaluationGate(30,90,40,.5)==65 &&
			FamilyTailCascade.evaluationGate(30,90,40,1)==90 && FamilyTailCascade.evaluationGate(30,90,40,2)==140,"Literal identity margin/floor oracle");
		check(FamilyTailCascade.evaluationGate(-.5f,.5f,0,.5)==.25 && FamilyTailCascade.evaluationGate(40,50,60,100)==60,
			"BLOSUM uses its own units and a floor-bound zero margin cannot move");
		final double tie=FamilyTailCascade.evaluationGate(.5f,.75f,0,2);
		check(tie==1 && !(1f>tie) && Math.nextUp(1f)>tie,"Strict acceptance rejects exact threshold ties");
		for(double bad:new double[]{-1,Double.NaN,Double.POSITIVE_INFINITY}){
			rejectScale(()->classifier.withEvaluationMarginScale(bad));
			rejectScale(()->classifier.withMemberCutoffMultiplier(bad));
			rejectScale(()->FamilyTailCascade.memberMultiplierGate(50,30,bad));
		}
		rejectScale(()->FamilyTailCascade.evaluationGate(0,Float.MAX_VALUE,0,Double.MAX_VALUE));
		final byte[] exact=sequence.getBytes(StandardCharsets.US_ASCII);
		final byte[] partial=(sequence.substring(0,160)+repeat("A",40)).getBytes(StandardCharsets.US_ASCII);
		final byte[] belowLowerIdentityAccept=(sequence.substring(0,80)+repeat("A",120)).getBytes(StandardCharsets.US_ASCII);
		final byte[] belowLowerIdentityReject=(sequence.substring(0,20)+repeat("A",180)).getBytes(StandardCharsets.US_ASCII);
		final FamilyTailCascade.LoadedClassifier member09=classifier.withMemberCutoffMultiplier(.9);
		final FamilyTailCascade.LoadedClassifier member0=classifier.withMemberCutoffMultiplier(0);
		final FamilyTailCascade.LoadedClassifier member1=classifier.withMemberCutoffMultiplier(1);
		final FamilyTailCascade.LoadedClassifier member175=classifier.withMemberCutoffMultiplier(1.75);
		check(member0.evaluationMemberMultiplierMode() && member0.evaluationMemberMultiplier()==0 &&
			member175.evaluationMemberMultiplierMode() && member175.evaluationMemberMultiplier()==1.75,
			"Active zero and positive member multipliers remain available in member-multiplier mode");
		checkClose(FamilyTailCascade.memberMultiplierGate(50,30,.9),45,
			"Member M.9 scales stored lowers above the floor");
		checkClose(FamilyTailCascade.memberMultiplierGate(50,30,1),50,
			"Member M1 scales stored lowers above the floor");
		checkClose(FamilyTailCascade.memberMultiplierGate(50,30,1.75),87.5,
			"Member M1.75 scales stored lowers above the floor");
		checkClose(FamilyTailCascade.memberMultiplierGate(50,60,.9),60,
			"Member multipliers bind the existing runtime floor");
		check(member09.memberMultiplierIdentityFloorBindings()==0 && member09.memberMultiplierBlosumFloorBindings()==0 &&
			member1.memberMultiplierIdentityFloorBindings()==0 && member1.memberMultiplierBlosumFloorBindings()==0 &&
			member175.memberMultiplierIdentityFloorBindings()==0 && member175.memberMultiplierBlosumFloorBindings()==0 &&
			member0.memberMultiplierIdentityFloorBindings()==2 &&
			member0.memberMultiplierBlosumFloorBindings()==2,
			"Floor-binding counters use runtime-floor ties as bindings");
		checkClose(member1.familyIdentityGateAt(0),classifier.withEvaluationMarginScale(0).familyIdentityGateAt(0),
			"M1 identity gate matches old margin0 on the same lower");
		checkClose(member1.familyBlosumGateAt(0),classifier.withEvaluationMarginScale(0).familyBlosumGateAt(0),
			"M1 BLOSUM gate matches old margin0 on the same lower");
		checkClose(member09.familyIdentityGateAt(0),45,"M.9 identity gate uses the stored member lower");
		checkClose(member175.familyIdentityGateAt(0),87.5,"M1.75 identity gate uses the stored member lower");
		check(classifier.newWorker().assign("partial",partial).assigned()==false,"Partial fixture must be rejected at frozen thresholds");
		check(member1.newWorker().assign("partial",partial).assigned(),"Member M1 lowers only acceptance gates while preserving the identity lower rejection");
		check(classifier.withEvaluationMarginScale(0).newWorker().assign("partial",partial).assigned(),"Loosening margins recovers the partial fixture without lowering member rejection");
		final FamilyTailCascade.Assignment acceptBelowLower=member0.newWorker().assign("below_lower_accept",belowLowerIdentityAccept);
		check(acceptBelowLower.identityMeasured && acceptBelowLower.identity>member0.familyIdentityGateAt(0) &&
			acceptBelowLower.identity<member0.familyEffectiveIdentityLowerAt(0) && acceptBelowLower.assigned() &&
			acceptBelowLower.terminal.equals("identity_accept"),
			"Member M0 preserves identity-first acceptance below the stored member lower");
		final IdentityRejectTrace lowerRejectTrace=new IdentityRejectTrace();
		final FamilyTailCascade.Assignment lowerRejected=member0.newWorker().assignTrace("below_lower_reject",
			belowLowerIdentityReject,lowerRejectTrace);
		check(lowerRejectTrace.alignedCandidates>0 && lowerRejectTrace.belowAcceptanceAndLower &&
			!lowerRejected.assigned() && lowerRejected.terminal.equals("identity_reject"),
			"Member M0 identity-rejects measured candidates below both floor and stored member lower");
		for(byte[] query:new byte[][]{exact,partial,repeat("A",200).getBytes(StandardCharsets.US_ASCII)}){
			boolean previous=true;
			for(double scale:new double[]{0,.25,.5,1,2,4}){
				final boolean accepted=classifier.withEvaluationMarginScale(scale).newWorker().assign("monotone",query).assigned();
				check(previous || !accepted,"Tightening only acceptance margins cannot gain accepted queries");previous=accepted;
			}
		}
		final FamilyTailCascade.LoadedClassifier strict=classifier.withEvaluationMarginScale(2);
		check(!strict.newWorker().assign("strict",exact).assigned() && classifier.newWorker().assign("original",exact).assigned(),
			"An evaluation view cannot change the original classifier");
		check(strict.withEvaluationMarginScale(0).newWorker().assign("absolute",partial).assigned() &&
			strict.withEvaluationMarginScale(1).newWorker().assign("restored",exact).assigned(),"Scales are absolute against frozen gates, never compounded");
		checkClose(strict.withMemberCutoffMultiplier(1).familyIdentityGateAt(0),member1.familyIdentityGateAt(0),
			"Member multipliers ignore old upper-margin state and replay from stored lowers");
	}
	private static void rejectScale(Runnable action){
		try{action.run();}catch(IllegalArgumentException expected){checks++;return;}
		throw new AssertionError("Invalid evaluation scale accepted");
	}
	/** Captures measured candidate scores before a rejected assignment clears its terminal scores. */
	private static final class IdentityRejectTrace implements FamilyTailCascade.CandidateObserver{
		@Override
		public void alignedCandidate(int familyOrdinal,int attempt,float identity,float blosum,
				double identityUpper,double blosumUpper,float identityLower){
			alignedCandidates++;
			if(identity<=identityUpper && identity<identityLower){belowAcceptanceAndLower=true;}
		}
		@Override
		public void accept(int queryLength,String terminal,int liveFamilyOrdinal,int ddlMatches,float identity,float blosum){}
		int alignedCandidates;
		boolean belowAcceptanceAndLower;
	}
	private static void checkClose(double observed,double expected,String why){
		check(Math.abs(observed-expected)<1e-9,why+" observed="+observed+" expected="+expected);
	}
	/** Direct structural scorer API is the independent count baseline; no hash-only CLI is invoked. */
	@SuppressWarnings("unchecked")
	private static HashMap<String,Long>[] baselines(Path root,Path refs,Path[] panels,int[] queries,String policy) throws Exception{
		final HashMap<String,Long>[] result=new HashMap[panels.length];
		final DDLFamilySearch.StructuralInputs inputs=new DDLFamilySearch.StructuralInputs();
		final FamilyTailCascade.LoadedClassifier classifier=loadFixture(refs,policy,1,inputs);
		for(int i=0; i<panels.length; i++){
			final List<String> records=lines(panels[i]);
			check(records.size()==2*queries[i] && population(panels[i].toString())==queries[i],
					"Literal two-line baseline keeps its independent complete panel population");
			final FamilyTailCascade.Worker worker=classifier.newWorker();
			result[i]=new HashMap<String,Long>();
			for(int q=0; q<queries[i]; q++){
				final String header=records.get(2*q);
				check(header.startsWith(">"),"Literal baseline fixture requires a FASTA header");
				final FamilyTailCascade.Assignment assignment=worker.assign(header.substring(1),
						records.get(2*q+1).getBytes(StandardCharsets.US_ASCII));
				if(assignment.assigned()){
					result[i].put(assignment.familyId,result[i].getOrDefault(assignment.familyId,0L)+1);
				}
			}
		}
		inputs.verify();
		return result;
	}
	private static FamilyTailCascade.LoadedClassifier loadFixture(Path refs,String policy,int threads,
			DDLFamilySearch.StructuralInputs inputs) throws Exception{
		final HashMap<String,String> options=new HashMap<String,String>();
		for(String line : policy.trim().split("\n")){
			final int equal=line.indexOf('=');
			options.put(line.substring(0,equal),line.substring(equal+1));
		}
		options.put("familyreferences",refs.toString());
		options.put("familyreferencessha80","unmeasured");
		return FamilyTailCascade.loadCallableStructural(options,threads,inputs);
	}

	private static FamilyReference family(String id,int rank,String sequence){
		final byte[] raw=sequence.getBytes(StandardCharsets.US_ASCII);
		final char[] ddl=new char[1024];DDLFamilyProfiles.sketch(raw,5,0,ddl);
		final ProteinSpectrum.Spectrum counts=new ProteinSpectrum.Spectrum(2);counts.add(raw);
		return new FamilyReference(id,rank,1,raw.length,raw.length-4,new ProteinSequence(id,raw),ddl,
			FamilyReferenceExport.frequencies(counts),0,50,.5f,1.1f,99,.9f);
	}
	private static List<String> lines(Path path){
		final ArrayList<String> out=new ArrayList<String>();final ByteFile in=ByteFile.makeByteFile(path.toString(),false);
		try{for(byte[] row=in.nextLine();row!=null;row=in.nextLine()){out.add(new String(row,StandardCharsets.US_ASCII));}}
		finally{check(!in.close(),"Fixture table read must complete");}
		return out;
	}
	private static void reject(String[] args,String why) throws Exception{
		boolean failed=false;try{evaluate(args);}catch(IllegalArgumentException|IllegalStateException expected){failed=true;}
		check(failed,why);
	}
	private static String repeat(String text,int times){final StringBuilder out=new StringBuilder();for(int i=0;i<times;i++){out.append(text);}return out.toString();}
	private static void fasta(StringBuilder out,String id,String sequence){out.append('>').append(id).append('\n').append(sequence).append('\n');}
	/** Declare populations from the literal fixture text; the product independently verifies its stream. */
	private static void write(Path path,String text) throws Exception{
		Files.write(path,text.getBytes(StandardCharsets.US_ASCII));
		if(path.toString().endsWith(".faa")){
			long count=0;
			for(String line : text.split("\n")){if(line.startsWith(">")){count++;}}
			populations.put(path.toAbsolutePath().normalize().toString(),count);
		}
	}
	private static long population(String path){
		final Long count=populations.get(Paths.get(path).toAbsolutePath().normalize().toString());
		check(count!=null,"Fixture panel must have a declared literal population: "+path);
		// The intentional empty-input case must reach the product reader with a positive declaration.
		return Math.max(1,count);
	}
	private static void evaluate(String[] args) throws Exception{
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
			if(!fields.containsKey("negcounts")){resolved.add("negcounts="+counts);}
		}
		FamilyEval.main(resolved.toArray(new String[resolved.size()]));
	}
	private static final HashMap<String,Long> populations=new HashMap<String,Long>();
	private static void check(boolean valid,String why){if(!valid){throw new AssertionError(why);}checks++;}
	private static int checks;
}
