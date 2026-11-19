package test.prot.cluster;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

import ddl.DDLPentamerBlacklist;
import ddl.ProteinBlacklistPolicy;
import ddl.ProteinBlacklistPolicy.Descriptor;
import ddl.ProteinKmerGeometry;
import dna.Data;
import fileIO.ByteFile;
import prot.cluster.ClusterSketchCodec;
import prot.cluster.DdlClust;
import prot.cluster.DdlClustStageB;

/** Real seed, Cake and disk Stage-B policy-routing regressions; no scientific corpus.
 * Synthetic outputs use structural verification and retain existing official identities.
 * @author Qiqi
 */
public final class SelectedBlacklistClusteringTest{

	public static void main(String[] args)throws Exception{
		if(args.length==2 && args[0].equals("custom")){
			custom(Paths.get(args[1]));System.out.println("CUSTOM_BLACKLIST_CLUSTERING_TEST_PASS checks="+checks);return;
		}
		check(args.length==3 || args.length==4, "Require fresh absolute fixture directory, both monomer resources and optional selected d120 resource");
		final Path root=Paths.get(args[0]), resource12=Paths.get(args[1]), resource16=Paths.get(args[2]);
		check(root.isAbsolute() && !Files.exists(root), "Fresh absolute fixture directory");
		check(resource12.isAbsolute() && resource16.isAbsolute() &&
			resource12.getParent().resolve("bbgap").equals(resource16.getParent()),
			"Install-only official resources live under resources/bbgap");
		Data.DORI_RESOURCES=resource12.getParent().toString()+"/";
		final Descriptor p12=ProteinBlacklistPolicy.policy(ProteinBlacklistPolicy.GEOMETRY);
		final Descriptor p16=ProteinBlacklistPolicy.policy(ProteinKmerGeometry.make(5, ALPHABET16));
		check(resource12.getFileName().toString().equals(p12.resourceName) &&
			resource16.getFileName().toString().equals(p16.resourceName), "Canonical selected resource names");
		// Checkout resources precede DORI_RESOURCES; retained copies must match bytes, not paths.
		check(Arrays.equals(Files.readAllBytes(resource12), Files.readAllBytes(p12.resource())) &&
			Arrays.equals(Files.readAllBytes(resource16), Files.readAllBytes(p16.resource())), "Resolved selected resource bytes");
		Files.createDirectory(root);
		final Path input=root.resolve("tails.faa"), map=root.resolve("partition.tsv");
		final String a=repeat("ACDEFGHIKLMNPQRSTVWY", 6), b=repeat("YWVTSRQPNMLKIHGFEDCA", 6);
		write(input, ">a0\n"+a+"\n>a1\n"+a+"\n>b0\n"+b+"\n>b1\n"+b+"\n");
		write(map, "parent_input_index\tshard\tshard_input_index\tinput_id\n0\t0\t0\ta0\n1\t0\t1\ta1\n2\t0\t2\tb0\n3\t0\t3\tb1\n");
		final Path seed12=policy(root.resolve("p12"), input, map, p12), seed16=policy(root.resolve("p16"), input, map, p16);
		if(args.length==4){
			final Path resource120=Paths.get(args[3]);
			final Descriptor p120=ProteinBlacklistPolicy.policy(ProteinKmerGeometry.make(6, ProteinBlacklistPolicy.DIMER120_ALPHABET));
			check(resource120.isAbsolute() && resource120.getParent().equals(resource16.getParent()) &&
				resource120.getFileName().toString().equals(p120.resourceName), "Canonical selected d120 resource beside monomer resources");
			check(Arrays.equals(Files.readAllBytes(resource120), Files.readAllBytes(p120.resource())), "Resolved d120 resource bytes");
			final Path seed120=policy(root.resolve("p120"), input, map, p120);
			reject(()->run(cake(root.resolve("wrong_dimer_cake"), input, seed120, p12.geometry, explicit(p12))),
				"Seed model geometry or blacklist differs");
		}
		final ArrayList<String> wrongGeometry=seed(root.resolve("wrong_geometry"), input, p16.geometry, explicit(p16));
		replace(wrongGeometry, "passes", "1"); wrongGeometry.add("seeds="+seed12.resolve("families.pddl"));
		reject(()->run(wrongGeometry), "Seed model geometry or blacklist differs");
		final ArrayList<String> wrongCake=cake(root.resolve("wrong_cake_geometry"), input, seed16, p12.geometry, explicit(p12));
		reject(()->run(wrongCake), "Seed model geometry or blacklist differs");
		nonofficial(root.resolve("nonofficial"), input);
		System.out.println("SELECTED_BLACKLIST_CLUSTERING_TEST_PASS checks="+checks);
	}

	/** Explicit custom policies pass through real seed, Cake, disk union and preparation loading. */
	private static void custom(Path root)throws Exception{
		check(root.isAbsolute() && !Files.exists(root),"Custom fixture requires a fresh absolute directory");Files.createDirectory(root);
		final Path input=root.resolve("tails.faa"),map=root.resolve("partition.tsv");
		final String a=repeat("ACDEFGHIKLMNPQRSTVWY",6),b=repeat("YWVTSRQPNMLKIHGFEDCA",6);
		write(input,">a0\n"+a+"\n>a1\n"+a+"\n>b0\n"+b+"\n>b1\n"+b+"\n");
		write(map,"parent_input_index\tshard\tshard_input_index\tinput_id\n0\t0\t0\ta0\n1\t0\t1\ta1\n2\t0\t2\tb0\n3\t0\t3\tb1\n");
		final String aa="ARNDCEQGHILKMFPSTWYV";final StringBuilder dimer=new StringBuilder("dimer:1,3,5:");
		for(int c=0;c<120;c++){
			if(c>0){dimer.append('/');}
			for(int pair=c;pair<400;pair+=120){if(pair>c){dimer.append(',');}dimer.append(aa.charAt(pair/20)).append(aa.charAt(pair%20));}
		}
		final String[] alphabets={"all:A/C/DE/FY/G/IV/LM/NH/P/RQK/ST/W","all:A/C/D/E/FY/G/H/IV/LM/N/P/Q/RK/S/T/W",dimer.toString()};
		for(int arm=0;arm<alphabets.length;arm++){
			final int k=arm==1 ? 5 : 6;final String mode=arm==1 ? "literal" : "reduced_literal";
			final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(k,alphabets[arm]);
			check(ProteinBlacklistPolicy.policy(geometry)==null,"Fixture geometry must exercise the custom route");
			final Path dir=root.resolve("arm"+arm);Files.createDirectory(dir);
			final Path filter=dir.resolve("filter.faa");write(filter,">excluded\n"+(k==5 ? "AAAAA" : "AAAAAA")+"\n");
			final String[] policy={"blacklist="+filter,"blacklistsha80=unmeasured","blacklistwords=1","blacklistk="+k,"blacklistmode="+mode,"blacklistprovenance=custom_fixture"};
			final DDLPentamerBlacklist loaded=ProteinBlacklistPolicy.loadExplicit(geometry,filter.toString(),"unmeasured",1,k,mode,"custom_fixture",true);
			check(loaded.wordLength()==k && loaded.matchMode().label().equals(mode),"Loaded filter retains full-word k and matching space");
			loaded.validate(geometry);
			final java.util.HashMap<String,String> options=new java.util.HashMap<String,String>();
			for(String option:policy){final int at=option.indexOf('=');options.put(option.substring(0,at),option.substring(at+1));}
			options.put("k",Integer.toString(k));options.put("aabits",alphabets[arm]);
			final DDLPentamerBlacklist prepared=prot.FamilyBlacklistBinding.load(options,prot.FamilyCalibration.VerificationMode.STRUCTURAL);
			check(prepared.wordLength()==loaded.wordLength() && prepared.matchMode()==loaded.matchMode(),"Preparation uses the same explicit policy as clustering");
			final Path seed=dir.resolve("seed");run(seed(seed,input,geometry,policy));verifyCodec(seed.resolve("families.pddl"),geometry,"unmeasured");
			final Path cake=dir.resolve("cake");run(cake(cake,input,seed,geometry,policy));verifyCodec(cake.resolve("final/families.pddl.gz"),geometry,"unmeasured");
			check(setting(cake.resolve("layer1/run.config"),"blacklistmode").equals(mode) && setting(cake.resolve("layer1/run.config"),"blacklistk").equals(Integer.toString(k)),"Cake forwards mode and word length to child stages");
			final Path collected=dir.resolve("collected.faa"),roster=dir.resolve("roster.tsv"),next=dir.resolve("next");
			new DdlClustStageB(new String[]{"mode=collect","codecs="+seed.resolve("families.pddl"),"out="+collected,"roster="+roster,"ow=f"}).process();
			run(seed(next,collected,geometry,policy));final Path union=dir.resolve("union.pddl");
			stage(union(input,map,roster,seed,next,union,policy));verifyCodec(union,geometry,"unmeasured");
			final ArrayList<String> missing=seed(dir.resolve("missing_mode"),input,geometry,policy);
			missing.remove("blacklistmode="+mode);reject(()->run(missing),"requires structural verification");
			reject(()->ProteinBlacklistPolicy.loadExplicit(geometry,filter.toString(),"unmeasured",1,k,mode,"custom_fixture",false),"requires structural verification");
			if(k==6){reject(()->ProteinBlacklistPolicy.loadExplicit(geometry,filter.toString(),"unmeasured",1,5,mode,"custom_fixture",true),"word length must match");}
		}
		reject(()->new ClusterSketchCodec.Header(ProteinBlacklistPolicy.GEOMETRY,512,0,ddl.ProteinPositionMode.FIXED,4,2,0,"unmeasured"),"blacklist is required");
	}

	/** Exercises selected defaults, explicit bindings and rejected substitutions through all three entry points. */
	private static Path policy(Path root, Path input, Path map, Descriptor policy)throws Exception{
		Files.createDirectory(root);
		final Path supplied=root.resolve("seed_explicit"), implicit=root.resolve("seed_implicit");
		run(seed(supplied, input, policy.geometry, explicit(policy)));
		run(seed(implicit, input, policy.geometry, new String[0]));
		verifySeed(supplied, policy); verifySeed(implicit, policy);
		check(Arrays.equals(Files.readAllBytes(supplied.resolve("families.pddl")), Files.readAllBytes(implicit.resolve("families.pddl"))),
			"Explicit and implicit selected policy produce identical seed codec bytes");
		for(boolean specified:new boolean[]{false, true}){
			final Path out=root.resolve("cake_"+specified);
			run(cake(out, input, supplied, policy.geometry, specified ? explicit(policy) : new String[0]));
			verifyCodec(out.resolve("final/families.pddl.gz"), policy.geometry, policy.sha80);
			check(setting(out.resolve("layer1/run.config"), "blacklist").equals(policy.resource().toString()) &&
				setting(out.resolve("layer1/run.config"), "blacklistsha80").equals(policy.sha80) &&
				setting(out.resolve("layer1/run.config"), "blacklistwords").equals(Integer.toString(policy.words)),
				"Cake forwards the selected resource, identity and count through a real layer");
		}
		// Cake uses the resident Stage-B API. Exercise the separate disk loader too.
		final Path collected=root.resolve("collected.faa"), roster=root.resolve("roster.tsv"), stageB=root.resolve("stage_b");
		new DdlClustStageB(new String[]{"mode=collect", "codecs="+supplied.resolve("families.pddl"),
			"out="+collected, "roster="+roster, "ow=f"}).process();
		run(seed(stageB, collected, policy.geometry, explicit(policy)));
		for(boolean specified:new boolean[]{false, true}){
			final Path out=root.resolve("union_"+specified+".pddl");
			stage(union(input, map, roster, supplied, stageB, out, specified ? explicit(policy) : new String[0]));
			verifyCodec(out, policy.geometry, policy.sha80);
		}
		for(String key:new String[]{"blacklistsha80", "blacklistwords"}){
			final String wrong=key.equals("blacklistsha80") ? INHERITED : Integer.toString(policy.words+1);
			final ArrayList<String> badSeed=seed(root.resolve("bad_seed_"+key), input, policy.geometry, explicit(policy));
			replace(badSeed, key, wrong); reject(()->run(badSeed), "Clustering requires");
			final ArrayList<String> badCake=cake(root.resolve("bad_cake_"+key), input, supplied, policy.geometry, explicit(policy));
			replace(badCake, key, wrong); reject(()->run(badCake), "Cake requires");
			final ArrayList<String> badUnion=union(input, map, roster, supplied, stageB, root.resolve("bad_union_"+key+".pddl"), explicit(policy));
			replace(badUnion, key, wrong); reject(()->stage(badUnion), key.equals("blacklistsha80") ? "does not match" : "Continuation requires");
		}
		final String original=text(policy.resource());
		final String changed=policy.geometry.name().startsWith("dimer") ? forgeDimer(original, policy) :
			policy.filterK==5 ? original.replaceFirst("ARLRR", "ARMRR") : original.replaceFirst("AAAAAA", "CCCCCC");
		check(!changed.equals(original), "Forged resource changes one known key without changing its cardinality");
		final Path forged=root.resolve("forged.faa"); write(forged, changed);
		final ArrayList<String> forgedSeed=seed(root.resolve("forged_seed"), input, policy.geometry, explicit(policy));
		replace(forgedSeed, "blacklist", forged.toString());
		final ArrayList<String> forgedCake=cake(root.resolve("forged_cake"), input, supplied, policy.geometry, explicit(policy));
		replace(forgedCake, "blacklist", forged.toString());
		final ArrayList<String> forgedUnion=union(input, map, roster, supplied, stageB, root.resolve("forged_union.pddl"), explicit(policy));
		replace(forgedUnion, "blacklist", forged.toString());
		if(policy.geometry.name().startsWith("dimer")){
			run(forgedSeed);run(forgedCake);stage(forgedUnion);
			verifyCodec(root.resolve("forged_seed/families.pddl"), policy.geometry, policy.sha80);
			verifyCodec(root.resolve("forged_cake/final/families.pddl.gz"), policy.geometry, policy.sha80);
			verifyCodec(root.resolve("forged_union.pddl"), policy.geometry, policy.sha80);
		}else{
			reject(()->run(forgedSeed), "keys differ");
			reject(()->run(forgedCake), "keys differ");
			reject(()->stage(forgedUnion), "keys differ");
		}
		final DDLPentamerBlacklist wrongMode=DDLPentamerBlacklist.loadStructural(
			headerless(root.resolve("wrong_mode.faa"), policy.resource()), policy.sha80, policy.words,
			policy.filterK, policy.filterK==5 ? DDLPentamerBlacklist.MatchMode.PREFIX : DDLPentamerBlacklist.MatchMode.LITERAL, "wrong_mode_fixture");
		reject(()->ProteinBlacklistPolicy.resolve(policy.geometry, wrongMode), policy.label);
		return supplied;
	}

	/** Replaces one representative by an absent reduced key, preserving the selected count and claimed pin. */
	private static String forgeDimer(String original, Descriptor policy)throws Exception{
		final DDLPentamerBlacklist official=policy.load(policy.resource(), true);
		final String aa="ARNDCEQGHILKMFPSTWYV";
		final byte[] word=new byte[6];
		for(int value=0; value<64000000; value++){
			int remaining=value;
			for(int i=5; i>=0; i--){word[i]=(byte)aa.charAt(remaining%20); remaining/=20;}
			if(!official.contains(official.rawWord(policy.geometry, DDLPentamerBlacklist.pack(word, 6)))){
				return original.replaceFirst("(?m)^[A-Z]{6}$", new String(word, StandardCharsets.US_ASCII));
			}
		}
		throw new AssertionError("No word outside the selected dimer blacklist");
	}

	/** Preserves the existing raw K5/prefix route for exact and nonofficial dimer geometries. */
	private static void nonofficial(Path root, Path input)throws Exception{
		Files.createDirectory(root);
		final Path filter=root.resolve("raw_prefix.faa"); write(filter, ">raw\nAAAAA\n");
		final String[] explicit={"blacklist="+filter, "blacklistsha80="+INHERITED, "blacklistwords=1"};
		final String aa="ARNDCEQGHILKMFPSTWYV"; final StringBuilder groups=new StringBuilder("dimer:1,3,5:");
		for(int c=0; c<120; c++){
			if(c>0){groups.append('/');}
			for(int pair=c; pair<400; pair+=120){if(pair>c){groups.append(',');}groups.append(aa.charAt(pair/20)).append(aa.charAt(pair%20));}
		}
		for(ProteinKmerGeometry geometry:new ProteinKmerGeometry[]{ProteinKmerGeometry.exact(5), ProteinKmerGeometry.make(6, groups.toString())}){
			check(ProteinBlacklistPolicy.policy(geometry)==null, "Nonofficial geometry has no selected monomer policy");
			final Path seed=root.resolve("seed_k"+geometry.k()); run(seed(seed, input, geometry, explicit));
			verifyCodec(seed.resolve("families.pddl"), geometry, INHERITED);
			final Path cake=root.resolve("cake_k"+geometry.k()); run(cake(cake, input, seed, geometry, explicit));
			verifyCodec(cake.resolve("final/families.pddl.gz"), geometry, INHERITED);
		}
		// This is the pre-existing raw K5/prefix route, not reduced-dimer file-backed clustering.
	}

	/** Builds a fresh greedy-seed invocation with explicit artifact paths and no refinement passes. */
	private static ArrayList<String> seed(Path out, Path input, ProteinKmerGeometry geometry, String[] filter)throws Exception{
		Files.createDirectory(out);
		final ArrayList<String> args=common(input, geometry, filter);
		args.addAll(Arrays.asList("out=null", "members="+out.resolve("members.tsv"), "stats="+out.resolve("stats.tsv"),
			"ddlout="+out.resolve("families.pddl"), "consensus=null", "refinement=null", "quality=null", "hist=null", "pairhist=null",
			"checkpoints=null", "selectedseeds=null", "passes=0", "mincluster=1", "topn=0", "merge=f", "prunefinal=f", "ow=f"));
		return args;
	}
	/** Builds one real Cake layer while disabling structural proposals to isolate policy routing. */
	private static ArrayList<String> cake(Path out, Path input, Path seed, ProteinKmerGeometry geometry, String[] filter){
		final ArrayList<String> args=common(input, geometry, filter);
		args.addAll(Arrays.asList("mode=cake", "codecs="+seed.resolve("families.pddl"), "outdir="+out, "familybudget=2", "initialfraction=1",
			"fractions=.25", "passes=1", "stopchangefraction=0", "fissiontopfraction=0", "fusionbottomfraction=0",
			"mergeminid=1", "mergeminratio=0", "frozenindexthreads=1"));
		return args;
	}
	/** Shares synthetic geometry and strict alignment settings across seed and Cake invocations. */
	private static ArrayList<String> common(Path input, ProteinKmerGeometry geometry, String[] filter){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+input, "k="+geometry.k(), "aabits="+alphabet(geometry),
			"width=512", "sketchseed=0", "positional=t", "partitions=4", "overlap=2", "t=1", "memberranking=f",
			"verification=structural", "minlen=100", "minspan=70", "maxalignlen=100", "minid=1", "minblosumratio=0",
			"passhits=0", "greedyhits=0", "maxhydrodelta=9", "greedymaxhydrodelta=9", "trimleft=f", "consensustrimleft=f"));
		args.addAll(Arrays.asList(filter)); return args;
	}
	/** Binds a disk Stage-B union to the actual seed memberships, collected roster and original inputs. */
	private static ArrayList<String> union(Path input, Path map, Path roster, Path seed, Path stageB, Path out, String[] filter){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("mode=union", "in="+input, "out="+out,
			"roster="+roster, "partitionmap="+map, "stageacodecs="+seed.resolve("families.pddl"), "stageamembers="+seed.resolve("members.tsv"),
			"stagebcodec="+stageB.resolve("families.pddl"), "stagebmembers="+stageB.resolve("members.tsv"), "report="+out+".tsv",
			"minspan=70", "maxalignlen=100", "minid=1", "minblosumratio=0", "verification=structural", "t=1", "ow=f"));
		args.addAll(Arrays.asList(filter)); return args;
	}
	/** Supplies the accepted resource, inherited identity and exact cardinality explicitly. */
	private static String[] explicit(Descriptor policy){
		return new String[]{"blacklist="+policy.resource(), "blacklistsha80="+policy.sha80, "blacklistwords="+policy.words};
	}
	/** Converts the stored geometry name back to its corresponding CLI alphabet syntax. */
	private static String alphabet(ProteinKmerGeometry geometry){
		if(geometry.isReduced()){
			return "all:"+geometry.name().substring(geometry.name().indexOf('+')+1);
		}
		if(geometry.name().startsWith("dimer")){
			return "dimer:"+geometry.name().substring(geometry.name().indexOf('+')+1).replace('+', ':');
		}
		return "none";
	}
	/** Checks both persisted model identity and the seed statistics' effective blacklist count. */
	private static void verifySeed(Path root, Descriptor policy)throws Exception{
		verifyCodec(root.resolve("families.pddl"), policy.geometry, policy.sha80);
		final String stats=new String(Files.readAllBytes(root.resolve("stats.tsv")), StandardCharsets.US_ASCII);
		check(stats.contains("blacklist_sha80\t"+policy.sha80+"\n") && stats.contains("blacklist_words\t"+policy.words+"\n"),
			"Seed statistics retain selected identity and cardinality");
	}
	/** Verifies persisted geometry, inherited identity and conservation of all four synthetic members. */
	private static void verifyCodec(Path path, ProteinKmerGeometry geometry, String identity){
		final ClusterSketchCodec.Loaded loaded=ClusterSketchCodec.read(path.toString());
		check(loaded.header.geometry.equals(geometry) && identity.equals(loaded.header.blacklistSha80), "Codec retains geometry and blacklist identity");
		int members=0; for(ClusterSketchCodec.Entry entry:loaded.families){members+=entry.members;}
		check(members==4, "All four original synthetic inputs survive the policy round trip");
	}
	/** Reads a required emitted configuration field; absence is a fixture failure. */
	private static String setting(Path path, String key)throws Exception{
		for(String line:Files.readAllLines(path, StandardCharsets.US_ASCII)){if(line.startsWith(key+"=")){return line.substring(key.length()+1);}}
		throw new AssertionError("Missing setting "+key);
	}
	/** Replaces an existing argument so negative cases cannot accidentally append a second binding. */
	private static void replace(ArrayList<String> args, String key, String value){
		for(int i=0; i<args.size(); i++){if(args.get(i).startsWith(key+"=")){args.set(i, key+"="+value); return;}}
		throw new AssertionError("Missing argument "+key);
	}
	private static void run(ArrayList<String> args){DdlClust.main(args.toArray(new String[args.size()]));}
	private static void stage(ArrayList<String> args)throws Exception{new DdlClustStageB(args.toArray(new String[args.size()])).process();}
	/** Constructs deterministic synthetic sequences without a biological input population. */
	private static String repeat(String text, int count){
		final StringBuilder b=new StringBuilder();
		for(int i=0; i<count; i++){b.append(text);}
		return b.toString();
	}
	private static void write(Path path, String text)throws Exception{Files.write(path, text.getBytes(StandardCharsets.US_ASCII));}
	private static Path headerless(Path path, Path resource)throws Exception{
		write(path, text(resource).replaceFirst("(?m)^>[^\\n]*", ">headerless"));
		return path;
	}
	private static String text(Path path){
		final StringBuilder sb=new StringBuilder();
		final ByteFile in=ByteFile.makeByteFile(path.toString(), false);
		for(byte[] line=in.nextLine(); line!=null; line=in.nextLine()){
			sb.append(new String(line, StandardCharsets.US_ASCII)).append('\n');
		}
		if(in.close()){throw new RuntimeException("Error reading or closing "+path);}
		return sb.toString();
	}
	private static void check(boolean valid, String why){checks++; if(!valid){throw new AssertionError(why);}}
	/** Requires the intended argument rejection, including causes wrapped by a clustering entry point. */
	private static void reject(Action action, String reason)throws Exception{
		try{action.run();}catch(IllegalArgumentException expected){
			for(Throwable cause=expected; cause!=null; cause=cause.getCause()){
				if(cause.getMessage()!=null && cause.getMessage().contains(reason)){checks++; return;}
			}
			throw new AssertionError("Expected rejection containing "+reason, expected);
		}
		throw new AssertionError("Expected rejection containing "+reason);
	}
	private interface Action{void run()throws Exception;}
	private static final String ALPHABET16="all:A/C/D/E/FWY/G/H/IV/K/LM/N/P/Q/R/S/T", INHERITED="00000000000000000000";
	private static int checks;
}
