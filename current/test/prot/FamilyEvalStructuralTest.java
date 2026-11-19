package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.attribute.FileTime;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ddl.DDLFamilyProfiles;
import fileIO.ByteFile;
import prot.FamilyEval;
import prot.FamilyPanelRunner.Panel;
import prot.FamilyCalibration.VerificationMode;
import prot.FamilyReference;
import prot.FamilyReferenceCodec;
import prot.FamilyReferenceExport;
import prot.ProteinSequence;
import prot.ProteinSpectrum;
import structures.ByteBuilder;

/** Ordinary structural screens, including inherited common-pivot ranking.
 * Synthetic files use explicit unmeasured or dummy inherited labels; no hashes.
 * Execute only in the combined Dori validation build.
 * @author Yaoyao
 */
public final class FamilyEvalStructuralTest{
	public static void main(String[] args)throws Exception{
		check(args.length==1,"Fresh absolute fixture directory required");
		final Path root=Paths.get(args[0]);check(root.isAbsolute(),"Absolute fixture path");Files.createDirectory(root);
		prepare(root);
		preparePlain(root.resolve("plain"));
		guards(root);
		for(int threads:new int[]{1,2}){
			final String name="screen"+threads;
			FamilyEval.main(arguments(root,name,"0.985",threads));
			checkSummary(root.resolve(name),0);
			final List<String> policy=lines(root.resolve(name+"/evaluation_policy.tsv"));
			for(String row:new String[]{"marginscale\t0.985","heuristic\thitfourth_h020_c150","cosine_weight\t1.5",
				"hydrophobicity_weight\t0.2","verification\tstructural","top\t50","ddl_reject\t5","min_query_aa\t40","max_query_aa\t99"}){
				check(policy.contains(row),"Effective policy retains "+row);
			}
			check(lines(root.resolve(name+"/input_populations.tsv")).equals(Arrays.asList("panel\texpected_queries","positive\t1","k5\t2","shadow\t1","reversed\t1","shuffled\t1")),"Declared populations retained in panel order");
			check(lines(root.resolve(name+"/inputs.tsv")).contains("package\t"+root.resolve("classifier.config")+"\tunmeasured"),"Recipe identity is honestly unmeasured");
			check(lines(root.resolve(name+"/classifier.config")).contains("verification=structural"),"Emitted recipe cannot masquerade as a hash-verified package");
		}
		for(String file:new String[]{"summary.tsv","family_false_accepts.tsv.gz","input_populations.tsv"}){
			check(lines(root.resolve("screen1/"+file)).equals(lines(root.resolve("screen2/"+file))),"Thread parity: "+file);
		}
		final String[] omitted=remove(remove(remove(remove(remove(arguments(root,"omitted_identities","0.985",1),
			"verification"),"possha80"),"negsha80"),"poscount"),"negcounts");
		FamilyEval.main(omitted);checkSummary(root.resolve("omitted_identities"),0);
		check(lines(root.resolve("omitted_identities/inputs.tsv")).equals(Arrays.asList("panel\tpath\texpected_sha80",
			"package\t"+root.resolve("classifier.config")+"\tunmeasured","positive\t"+root.resolve("positive.faa")+"\tNA",
			"k5\t"+root.resolve("k5.faa")+"\tNA","shadow\t"+root.resolve("shadow.faa")+"\tNA",
			"reversed\t"+root.resolve("reversed.faa")+"\tNA","shuffled\t"+root.resolve("shuffled.faa")+"\tNA")),
			"Omitted structural panel identities must remain exact NA rows");
		check(lines(root.resolve("omitted_identities/input_populations.tsv")).equals(Arrays.asList("panel\texpected_queries",
			"positive\tNA","k5\tNA","shadow\tNA","reversed\tNA","shuffled\tNA")),"Omitted structural populations must remain exact NA rows");
		final String[] negativesOnly=remove(remove(remove(arguments(root,"negatives_only","0.985",1),"pos"),"possha80"),"poscount");
		FamilyEval.main(negativesOnly);
		check(lines(root.resolve("negatives_only/summary.tsv")).size()==5 &&
			lines(root.resolve("negatives_only/PASS")).get(0).equals("FAMILY_EVAL_PASS panels=4 queries=5 live_families=1"),"Four negative-only screen populations are retained");
		FamilyEval.main(remove(remove(arguments(root,"no_overrides","0.985",1),"top"),"ddlreject"));checkSummary(root.resolve("no_overrides"),0);
		check(lines(root.resolve("screen1/summary.tsv")).equals(lines(root.resolve("no_overrides/summary.tsv"))),"Equivalent explicit overrides preserve inherited ranking results");
		FamilyEval.main(explicitHeuristic(root,"explicit_heuristic","unranked.config"));checkSummary(root.resolve("explicit_heuristic"),0);
		FamilyEval.main(explicitHeuristic(root,"redundant_heuristic","classifier.config"));checkSummary(root.resolve("redundant_heuristic"),0);
		reject(replace(explicitHeuristic(root,"mismatched_heuristic","classifier.config"),"familyhydrosha80","11111111111111111111"),
			"Evaluator member metadata must match");
		FamilyEval.main(append(replace(arguments(root,"family_mean","0.985",1),"ddlreject","6"),"admission=family_mean"));
		checkSummary(root.resolve("family_mean"),0);
		FamilyEval.main(append(replace(arguments(root,"short_scaled","0.985",1),"ddlreject","6"),"admission=short_scaled"));
		checkSummary(root.resolve("short_scaled"),0);
		final Path plain=root.resolve("plain");
		FamilyEval.main(plainArguments(plain,"fixed",1));
		check(Files.exists(plain.resolve("fixed/PASS")),
			"Plain exact K5 fixtures must admit fixed scoring before length-admission rejection");
		reject(append(plainArguments(plain,"family_mean",1),"admission=family_mean"),
			"Length-based admission requires the released");
		reject(append(arguments(root,"admission_policy","0.985",1),"admission=family_mean"),
			"Family-mean admission requires the released new pivot, top50 and fixed6 parent anchor");
		// With lower50, upper100.8, one .985 scale leaves the identity gate above100;
		// multiplying twice would lower it below100. BLOSUM's1.02 upper also rejects
		// the perfect query at .985. A .97 scale must accept the exact same sequences.
		FamilyEval.main(arguments(root,"relaxed","0.97",1));checkSummary(root.resolve("relaxed"),1);
		final String[] repeat=append(arguments(root,"repeat","0.985",2),"benchmarkrepeats=2");
		FamilyEval.main(repeat);checkSummary(root.resolve("repeat"),0);
		check(lines(root.resolve("repeat/benchmark_time.tsv")).contains("query_evaluations\t12"),"Repeated work does not enlarge the six-query scientific denominator");
		for(String count:new String[]{"1,1,1,1","3,1,1,1"}){
			reject(replace(arguments(root,"bad_population","0.985",1),"negcounts",count),"count differs from its declared population");
		}
		final Path panel=root.resolve("positive.faa");final String original=text(panel);
		final FamilyEval changed=new FamilyEval(arguments(root,"changed","0.985",1));
		Files.setLastModifiedTime(panel,FileTime.fromMillis(Files.getLastModifiedTime(panel).toMillis()+5000));
		rejectAction(changed::process,"Structural input changed");check(!Files.exists(root.resolve("changed/PASS")),"Changed input cannot publish PASS");
		write(panel,">bad\nACDE*FGHI\n");reject(arguments(root,"invalid_sequence","0.985",1),"Family panel execution failed");
		write(panel,">long\n"+sequence()+sequence()+"\n");reject(arguments(root,"invalid_length","0.985",1),"Family panel execution failed");
		write(panel,original);
		final Path medians=root.resolve("medians.tsv");write(medians,text(medians).replace("f0\t1\t60", "f0\t2\t60"));
		reject(arguments(root,"metadata_roster","0.985",1),"Median identity and count");
		System.out.println("FAMILY_EVAL_STRUCTURAL_TEST_PASS checks="+checks+" hashes=0");
	}

	private static void prepare(Path root)throws Exception{
		final byte[] raw=sequence().getBytes(StandardCharsets.US_ASCII);final char[] ddl=new char[1024];
		DDLFamilyProfiles.sketchPositional(raw,5,0,4,2,ddl);
		final ProteinSpectrum.Spectrum counts=new ProteinSpectrum.Spectrum(2);counts.add(raw);
		final FamilyReference live=new FamilyReference("f0",0,1,raw.length,raw.length-4,new ProteinSequence("f0",raw),ddl,
			FamilyReferenceExport.frequencies(counts),0,50,.5f,1.1f,100.8f,1.02f);
		final ByteBuilder refs=new ByteBuilder();
		FamilyReferenceCodec.appendHeader(new FamilyReferenceCodec.Header(FamilyReferenceCodec.SCHEMA,5,1024,0,2,1,1,
			FamilyReferenceCodec.FORMULA_VERSION,FamilyReferenceCodec.SPECTRUM_ENCODING),refs);
		FamilyReferenceCodec.appendFiveLines(live,refs);FamilyReferenceCodec.appendFiveLines(FamilyReference.empty("empty",1),refs);
		write(root.resolve("references.tsv"),refs.toString());
		final ByteBuilder profiles=new ByteBuilder().append(DDLFamilyProfiles.HEADER16).nl()
			.append("tail\tf0\t1024\t1\t60\t56\tnonempty\t");
		for(int field=0; field<3; field++){
			if(field>0){profiles.tab();}
			for(int bin=0; bin<ddl.length; bin++){
				if(bin>0){profiles.append(',');} profiles.append(field==0 ? (int)ddl[bin] : ddl[bin]==0 ? 0 : 1);
			}
		} profiles.nl();
		write(root.resolve("profiles.tsv"),profiles.toString());
		write(root.resolve("summary.tsv"),"namespace\tk\tseed\texponent\tregisterbits\tsketch\tpartitions\toverlap\ttag_bits\tsegment_words\n"+
			"tail\t5\t0\t6\t16\t"+DDLFamilyProfiles.SKETCH_POSITIONAL+"\t4\t2\t"+
			DDLFamilyProfiles.positionalTagBits(4)+"\t"+DDLFamilyProfiles.positionalSegmentWords(4)+"\n");
		write(root.resolve("profile_list.tsv"),"profiles.tsv\tunmeasured\tunmeasured\n");
		final StringBuilder bounds=new StringBuilder("scalar\trows\tmin\tmax\tmin_family\tmax_family\n");
		for(String name:new String[]{"hydrophobicity","net_charge","aromaticity","fraction_small","pI_tail"}){bounds.append(name).append("\t79565\t-20\t20\tlo\thi\n");}
		write(root.resolve("bounds.tsv"),bounds.toString());
		write(root.resolve("medians.tsv"),"family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200\nf0\t1\t60\t60\t60\t1\n");
		write(root.resolve("hydro.tsv"),"family_rank\tfamily_id\tN\thydrophobicity_mean\n0\tf0\t1\t0.0\n1\tempty\t0\tNA\n");
		final String ranked="verification=structural\nfamilyreferences=references.tsv\nfamilyreferencessha80=unmeasured\n"+
			"top=50\nk=5\nseed=0\nwidth=1024\nddlreject=5\nddlaccept=1025\ncosine=ranking\nrankingfeatures=combined25\n"+
			"scalargates=f\npositional=t\npartitions=4\noverlap=2\nddlprofiles=profile_list.tsv\nddlprofilessha80=unmeasured\n"+
			"identityacceptfloor=30\nblosumacceptfloor=0\nindexearlyreject=f\n"+
			"scalarbounds=bounds.tsv\nscalarboundssha80="+PIN+"\nscalarscale=0.1\nheuristic=hitfourth_h020_c150\n"+
			"lengthmedians=medians.tsv\nlengthmedianssha80="+PIN+"\nfamilyhydro=hydro.tsv\nfamilyhydrosha80="+PIN+"\nminqueryaa=40\nmaxqueryaa=99\n";
		write(root.resolve("classifier.config"),ranked);
		write(root.resolve("unranked.config"),ranked.replace("heuristic=hitfourth_h020_c150\n"+
			"lengthmedians=medians.tsv\nlengthmedianssha80="+PIN+"\nfamilyhydro=hydro.tsv\nfamilyhydrosha80="+PIN+"\n"+
			"minqueryaa=40\nmaxqueryaa=99\n",""));
		write(root.resolve("positive.faa"),">positive\n"+sequence()+"\n");
		write(root.resolve("k5.faa"),">match\n"+sequence()+"\n>miss\n"+repeat("C",60)+"\n");
		for(String name:new String[]{"shadow","reversed","shuffled"}){write(root.resolve(name+".faa"),">match\n"+sequence()+"\n");}
	}
	private static void preparePlain(Path root)throws Exception{
		Files.createDirectory(root);
		final byte[] raw=sequence().getBytes(StandardCharsets.US_ASCII);
		final char[] ddl=new char[1024];
		DDLFamilyProfiles.sketch(raw,5,0,ddl);
		final ProteinSpectrum.Spectrum counts=new ProteinSpectrum.Spectrum(2);counts.add(raw);
		final FamilyReference live=new FamilyReference("f0",0,1,raw.length,raw.length-4,
			new ProteinSequence("f0",raw),ddl,FamilyReferenceExport.frequencies(counts),0,50,.5f,1.1f,100.8f,1.02f);
		final ByteBuilder refs=new ByteBuilder();
		FamilyReferenceCodec.appendHeader(new FamilyReferenceCodec.Header(FamilyReferenceCodec.SCHEMA,
			5,1024,0,2,1,1,FamilyReferenceCodec.FORMULA_VERSION,FamilyReferenceCodec.SPECTRUM_ENCODING),refs);
		FamilyReferenceCodec.appendFiveLines(live,refs);
		FamilyReferenceCodec.appendFiveLines(FamilyReference.empty("empty",1),refs);
		write(root.resolve("references.tsv"),refs.toString());
		final StringBuilder bounds=new StringBuilder("scalar\trows\tmin\tmax\tmin_family\tmax_family\n");
		for(String name:new String[]{"hydrophobicity","net_charge","aromaticity","fraction_small","pI_tail"}){
			bounds.append(name).append("\t79565\t-20\t20\tlo\thi\n");
		}
		write(root.resolve("bounds.tsv"),bounds.toString());
		write(root.resolve("medians.tsv"),
			"family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200\n"+
			"f0\t1\t60\t60\t60\t1\n");
		write(root.resolve("hydro.tsv"),
			"family_rank\tfamily_id\tN\thydrophobicity_mean\n0\tf0\t1\t0.0\n1\tempty\t0\tNA\n");
		write(root.resolve("classifier.config"),
			"verification=structural\nfamilyreferences=references.tsv\nfamilyreferencessha80=unmeasured\n"+
			"top=50\nk=5\nseed=0\nwidth=1024\nddlreject=6\nddlaccept=1025\ncosine=ranking\n"+
			"rankingfeatures=combined25\nscalargates=f\npositional=f\nidentityacceptfloor=30\n"+
			"blosumacceptfloor=0\nindexearlyreject=f\nscalarbounds=bounds.tsv\nscalarboundssha80="+PIN+"\n"+
			"scalarscale=0.1\nheuristic=hitfourth_h020_c150\nlengthmedians=medians.tsv\n"+
			"lengthmedianssha80="+PIN+"\nfamilyhydro=hydro.tsv\nfamilyhydrosha80="+PIN+"\n"+
			"minqueryaa=40\nmaxqueryaa=99\n");
		write(root.resolve("positive.faa"),">positive\n"+sequence()+"\n");
		write(root.resolve("k5.faa"),">match\n"+sequence()+"\n>miss\n"+repeat("C",60)+"\n");
		for(String name:new String[]{"shadow","reversed","shuffled"}){
			write(root.resolve(name+".faa"),">match\n"+sequence()+"\n");
		}
	}
	private static void guards(Path root)throws Exception{
		for(String[] options:new String[][]{{"floors=does-not-exist"},{"fitting=does-not-exist"},{"insha80="+PIN}}){
			reject(append(arguments(root,"incomplete_floor","0.985",1),options),"Saved floors require");
		}
		final String[] paired=remove(remove(remove(remove(remove(remove(arguments(root,"incomplete_paired","0.985",1),
			"top"),"ddlreject"),"marginscale"),"neg"),"negnames"),"negsha80");
		reject(append(remove(paired,"negcounts"),"paired=t"),
			"Paired evaluation requires pinned lengthmodel/lengthmedians and three explicit pairedscales");
		reject(append(arguments(root,"paired_ordinary_overrides","0.985",1),"paired=t"),
			"Paired learned-ranking mode cannot receive ordinary cap, heuristic, KD or query-bound overrides");
		for(String[] options:new String[][]{{"lengthmodel=does-not-exist"},{"pairedscales=1,1,1"}}){
			reject(append(arguments(root,"unsupported_length","0.985",1),options),"Learned length-model inputs require explicit paired=t");
		}
		reject(append(arguments(root,"unsupported_metadata","0.985",1),"familyhydro=does-not-exist"),
			"Member medians/KD require an explicit heuristic");
		for(String value:new String[]{"", "none", "STRUCTURAL"}){reject(replace(arguments(root,"mode","0.985",1),"verification",value),"verification requires");}
		reject(append(arguments(root,"duplicate_mode","0.985",1),"verification=structural"),"verification requires");
		for(String value:new String[]{"0","-1","1.5","9223372036854775808"}){
			reject(replace(arguments(root,"population_syntax","0.985",1),"poscount",value),null);
		}
		reject(replace(arguments(root,"count_arity","0.985",1),"negcounts","2,1,1"),
			"negcounts must match the negative panel count");
		reject(replace(arguments(root,"duplicate_labels","0.985",1),"negnames","same,same,reversed,shuffled"),"Duplicate panel label");
		reject(replace(arguments(root,"duplicate_paths","0.985",1),"neg",root.resolve("k5.faa")+","+root.resolve("k5.faa")+","+root.resolve("reversed.faa")+","+root.resolve("shuffled.faa")),"A sequence file cannot be used as two");
		rejectAction(()->new Panel("default",root.resolve("positive.faa").toString(),"unmeasured",true),"Panel fingerprint requires20hex");
		new Panel("structural",root.resolve("positive.faa").toString(),"unmeasured",true,VerificationMode.STRUCTURAL,1);checks++;
		final Path link=root.resolve("panel-link.faa");Files.createSymbolicLink(link,root.resolve("positive.faa"));
		reject(replace(arguments(root,"symlink","0.985",1),"pos",link.toString()),"Require regular non-symlink");
	}
	private static String[] arguments(Path root,String name,String scale,int threads){
		return new String[]{"in="+root.resolve("classifier.config"),"out="+root.resolve(name),"verification=structural","marginscale="+scale,"t="+threads,"top=50","ddlreject=5",
			"pos="+root.resolve("positive.faa"),"possha80=unmeasured","poscount=1",
			"neg="+root.resolve("k5.faa")+","+root.resolve("shadow.faa")+","+root.resolve("reversed.faa")+","+root.resolve("shuffled.faa"),
			"negnames=k5,shadow,reversed,shuffled","negsha80="+PIN+",unmeasured,"+PIN+",unmeasured","negcounts=2,1,1,1"};
	}
	private static String[] plainArguments(Path root,String name,int threads){
		return replace(replace(arguments(root,name,"0.985",threads),
			"in",root.resolve("classifier.config").toString()),"ddlreject","6");
	}
	private static String[] explicitHeuristic(Path root,String name,String config){
		return append(replace(arguments(root,name,"0.985",1),"in",root.resolve(config).toString()),
			"heuristic=hitfourth_h020_c150","lengthmedians="+root.resolve("medians.tsv"),
			"lengthmedianssha80="+PIN,"familyhydro="+root.resolve("hydro.tsv"),
			"familyhydrosha80="+PIN,"minqueryaa=40","maxqueryaa=99");
	}
	private static void checkSummary(Path out,int matches){
		final List<String> rows=lines(out.resolve("summary.tsv"));check(rows.size()==6,"One positive and four negative rows");
		for(int i=0;i<5;i++){
			final String[] row=rows.get(i+1).split("\t",-1);final long total=i==1?2:1;
			check(Long.parseLong(row[2])==total && Long.parseLong(row[3])==matches && Long.parseLong(row[4])==total-matches,"Independent observed denominator and threshold outcome for panel"+i);
			long stages=0;for(int j=9;j<row.length;j++){stages+=Long.parseLong(row[j]);}check(stages==total,"Terminal stages conserve each panel denominator");
			if(i>0){check(Double.parseDouble(row[8])==matches*1000000.0/total,"Each negative FPR uses its own denominator");}
		}
		check(lines(out.resolve("family_false_accepts.tsv.gz")).size()==5,"Exactly one live-family winner row for each negative panel");
		check(lines(out.resolve("PASS")).get(0).equals("FAMILY_EVAL_PASS panels=5 queries=6 live_families=1"),"Scientific counts exclude repeats and empty reference slots");
	}
	private static String[] append(String[] args,String... extra){final ArrayList<String> result=new ArrayList<String>(Arrays.asList(args));result.addAll(Arrays.asList(extra));return result.toArray(new String[0]);}
	private static String[] replace(String[] args,String key,String value){final String[] copy=args.clone();for(int i=0;i<copy.length;i++){if(copy[i].startsWith(key+"=")){copy[i]=key+"="+value;return copy;}}throw new AssertionError(key);}
	private static String[] remove(String[] args,String key){final ArrayList<String> result=new ArrayList<String>();for(String a:args){if(!a.startsWith(key+"=")){result.add(a);}}return result.toArray(new String[0]);}
	private static void reject(String[] args,String reason)throws Exception{
		rejectAction(()->FamilyEval.main(args),reason);
		for(String a:args){if(a.startsWith("out=")){check(!Files.exists(Paths.get(a.substring(4)).resolve("PASS")),"Rejected screen cannot publish PASS");}}
	}
	private interface Action{void run()throws Exception;}
	private static void rejectAction(Action action,String reason)throws Exception{
		boolean rejected=false;try{action.run();}catch(IllegalArgumentException|IllegalStateException e){rejected=true;check(reason==null || e.getMessage().contains(reason),"Expected "+reason+", got "+e.getMessage());}
		check(rejected,"Invalid evaluation must fail");
	}
	private static String sequence(){return repeat("ACDEFGHIKLMNPQRSTVWY",3);}
	private static String repeat(String text,int count){final StringBuilder b=new StringBuilder();for(int i=0;i<count;i++){b.append(text);}return b.toString();}
	private static String text(Path path)throws Exception{return new String(Files.readAllBytes(path),StandardCharsets.US_ASCII);}
	private static void write(Path path,String text)throws Exception{Files.write(path,text.getBytes(StandardCharsets.US_ASCII));}
	private static List<String> lines(Path path){
		final ArrayList<String> rows=new ArrayList<String>();final ByteFile in=ByteFile.makeByteFile1(path.toString(),false);
		try{for(byte[] row=in.nextLine();row!=null;row=in.nextLine()){rows.add(new String(row,StandardCharsets.US_ASCII));}}
		finally{check(!in.close(),"Fixture output read completes");}return rows;
	}
	private static void check(boolean valid,String why){checks++;if(!valid){throw new AssertionError(why);}}
	private static final String PIN="00000000000000000000";
	private static int checks;
}
