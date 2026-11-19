package test.ddl;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import ddl.DDLFamilyProfiles;
import ddl.DDLFamilySearch;
import ddl.DDLPentamerBlacklist;
import ddl.DDLQueryOrigins;
import ddl.FamilySketchBuilder;
import ddl.ProteinBlacklistPolicy;
import ddl.ProteinBlacklistPolicy.Descriptor;
import ddl.ProteinKmerGeometry;
import ddl.ProteinPositionMode;
import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import prot.FamilyBlacklistBinding;
import prot.FamilyCalibration.VerificationMode;
import prot.cluster.ClusterSketchCodec;
import prot.cluster.DdlClust;
import shared.Shared;

/** Official resource, kernel, producer and persisted-policy regressions.
 * Official resource loads retain their existing verification; synthetic artifacts use structural loading.
 * @author Yaoyao
 */
public final class ProteinBlacklistPolicyTest{

	public static void main(String[] args)throws Exception{
		if(args.length<2 || args.length>3){throw new IllegalArgumentException("Require NEW_FIXTURE_DIR, OFFICIAL_12K6_RESOURCE and optional OFFICIAL_16K5_RESOURCE");}
		final Path root=Paths.get(args[0]).toAbsolutePath(),resource=Paths.get(args[1]).toAbsolutePath();
		check(resource.getFileName().toString().equals(ProteinBlacklistPolicy.RESOURCE),"Canonical resource name");
		Files.createDirectory(root);Data.DORI_RESOURCES=resource.getParent().toString()+"/";
		Shared.setThreads(1);Shared.AMINO_IN=true;
		final DDLPentamerBlacklist official=ProteinBlacklistPolicy.resolve(GEOMETRY,null);
		check(official.size()==4096 && official.sha80().equals(ProteinBlacklistPolicy.SHA80) &&
			official.matchMode()==DDLPentamerBlacklist.MatchMode.REDUCED_LITERAL,"Official resource binding");
		final char[] raw=new char[512],encoded=new char[512],positions=new char[512];
		int records=0;
		for(String line:readLines(resource)){
			if(line.startsWith(">")){continue;}
			final byte[] sequence=ascii(line),codes=sequence.clone();
			for(int i=0;i<codes.length;i++){codes[i]=AminoAcid.acidToNumber[codes[i]];}
			check(sequence.length==6,"One excluded hexamer per record");
			check(DDLFamilyProfiles.sketch(sequence,GEOMETRY,0,raw,null)==1 && empty(raw),"Every official word is excluded by default");
			DDLFamilyProfiles.sketchEncoded(codes,GEOMETRY,0,encoded,null);
			check(Arrays.equals(raw,encoded),"Encoded and ASCII exclusion agree");
			DDLFamilyProfiles.sketchPositional(sequence,GEOMETRY,0,4,2,positions,null);
			check(empty(positions),"Filtering precedes positional tags");records++;
		}
		check(records==4096,"Complete official resource exercised");
		DDLFamilyProfiles.sketch(ascii("AAAAAK"),GEOMETRY,0,raw,null);
		check(empty(raw),"Canonical residue aliases are filtered in reduced coordinates");
		DDLFamilyProfiles.sketch(ascii("CCCCCC"),GEOMETRY,0,raw,null);
		check(!empty(raw),"A word outside the official set survives");
		DDLFamilyProfiles.sketchPositional(ascii("AAAAAA"),GEOMETRY,0,4,2,positions,null,new DDLQueryOrigins(512),ProteinPositionMode.FIXED);
		check(empty(positions),"Origin tracing uses the same filter");
		final FamilySketchBuilder builder=new FamilySketchBuilder(512,GEOMETRY,0,ProteinPositionMode.FIXED,4,2,null);
		builder.addRaw(ascii("AAAAAA"));builder.writeProfile(raw);
		check(empty(raw) && builder.fork().matchesGeometry(512,GEOMETRY,0,ProteinPositionMode.FIXED,4,2,official.sha80()),
			"Construction and fork retain official filter identity");
		for(ProteinKmerGeometry geometry:new ProteinKmerGeometry[]{ProteinKmerGeometry.exact(6),
			ProteinKmerGeometry.make(5,"all:A/C/DE/F/W/Y/G/I/V/L/M/NH/P/Q/RK/ST")}){
			check(ProteinBlacklistPolicy.resolve(geometry,null)==null,"Other geometries retain unfiltered default");
			DDLFamilyProfiles.sketch(ascii("AAAAAA"),geometry,0,raw,null);check(!empty(raw),"Non12K6 kernels remain unaffected");
		}
		final DDLPentamerBlacklist wrongMode=DDLPentamerBlacklist.loadStructural(resource,official.sha80(),4096,6,
			DDLPentamerBlacklist.MatchMode.LITERAL,"fixture_wrong_mode");
		reject(()->ProteinBlacklistPolicy.resolve(GEOMETRY,wrongMode),"official reduced-shadow4096");
		final Path wrongWords=root.resolve("wrong_words.faa");
		final String original=new String(Files.readAllBytes(resource),StandardCharsets.US_ASCII);
		write(wrongWords,original.replaceFirst("AAAAAA","CCCCCC"));
		final DDLPentamerBlacklist forged=DDLPentamerBlacklist.loadReducedStructural(wrongWords,official.sha80(),4096,GEOMETRY,ProteinBlacklistPolicy.PROVENANCE);
		reject(()->ProteinBlacklistPolicy.resolve(GEOMETRY,forged),"keys differ");
		final Descriptor policy12=ProteinBlacklistPolicy.policy(GEOMETRY);
		configuration(official,policy12);profiles(root.resolve("profiles"),official,policy12,"CCCCCC");cluster(root.resolve("cluster"),policy12);
		raw16(root.resolve("raw16"),resource,args.length==3 ? Paths.get(args[2]).toAbsolutePath() :
			resource.getParent().resolve("bbgap").resolve(ProteinBlacklistPolicy.policy(GEOMETRY16).resourceName));
		test.prot.FamilyResidentCalibrationTest.official12(root.resolve("resident"));
		System.out.println("PROTEIN_BLACKLIST_POLICY_TEST_PASS checks="+checks);
	}

	/** Raw matching precedes the16-letter projection; aliases are not automatically excluded. */
	private static void raw16(Path root,Path resource12,Path resource16)throws Exception{
		Files.createDirectory(root);
		final Descriptor policy=ProteinBlacklistPolicy.policy(GEOMETRY16),policy12=ProteinBlacklistPolicy.policy(GEOMETRY);
		check(policy!=null && policy!=policy12 && policy.filterK==5 && policy.words==1024 &&
			policy.mode==DDLPentamerBlacklist.MatchMode.LITERAL && policy.sha80.equals("f731f9dfe82458b93234"),"Exact16K5 descriptor");
		check(resource16.getFileName().toString().equals(policy.resourceName) &&
			resource16.getParent().equals(resource12.getParent().resolve("bbgap")),
			"The raw16 official resource resolves under resources/bbgap");
		check(ProteinBlacklistPolicy.resource().equals(policy12.resource()) && ProteinBlacklistPolicy.SHA80.equals(policy12.sha80) &&
			ProteinBlacklistPolicy.WORDS==policy12.words,"Legacy constants and resource retain12K6 semantics");
		check(ProteinBlacklistPolicy.policy(null)==null && ProteinBlacklistPolicy.policy(ProteinKmerGeometry.exact(5))==null &&
			ProteinBlacklistPolicy.policy(ProteinKmerGeometry.make(6,alphabet(GEOMETRY16)))==null,"No default for null/raw/different-K geometry");
		final DDLPentamerBlacklist official=ProteinBlacklistPolicy.resolve(GEOMETRY16,null);
		check(official.size()==1024 && official.matchMode()==policy.mode && official.sha80().equals(policy.sha80),"Selected raw-shadow1024 resource");
		final char[] raw=new char[512],encoded=new char[512],positional=new char[512];int records=0;
		for(String line:readLines(resource16)){
			if(line.startsWith(">")){continue;}
			final byte[] sequence=ascii(line),codes=sequence.clone();
			check(sequence.length==5,"One excluded pentamer per record");
			for(int i=0;i<codes.length;i++){codes[i]=AminoAcid.acidToNumber[codes[i]];}
			check(DDLFamilyProfiles.sketch(sequence,GEOMETRY16,0,raw,null)==1 && empty(raw),"Every selected raw word excluded by default");
			DDLFamilyProfiles.sketchEncoded(codes,GEOMETRY16,0,encoded,null);check(Arrays.equals(raw,encoded),"Encoded raw exclusion parity");
			DDLFamilyProfiles.sketchPositional(sequence,GEOMETRY16,0,4,2,positional,null);check(empty(positional),"Raw exclusion precedes positional tags");records++;
		}
		check(records==1024,"Complete raw resource exercised");
		final long excluded=DDLPentamerBlacklist.pack(ascii("ARLRR"),5),allowed=DDLPentamerBlacklist.pack(ascii("ARMRR"),5);
		check(official.contains(excluded) && !official.contains(allowed) && GEOMETRY16.word(excluded)==GEOMETRY16.word(allowed),"Literal witnesses share a reduced word but differ in raw membership");
		DDLFamilyProfiles.sketch(ascii("ARLRR"),GEOMETRY16,0,raw,null);check(empty(raw),"Selected raw witness excluded");
		DDLFamilyProfiles.sketch(ascii("ARMRR"),GEOMETRY16,0,raw,null);check(!empty(raw),"Unlisted reduced alias survives");
		DDLFamilyProfiles.sketchPositional(ascii("ARLRR"),GEOMETRY16,0,4,2,positional,null,new DDLQueryOrigins(512),ProteinPositionMode.FIXED);
		check(empty(positional),"Origin tracing honors raw policy");
		final FamilySketchBuilder builder=new FamilySketchBuilder(512,GEOMETRY16,0,ProteinPositionMode.FIXED,4,2,null);
		builder.addRaw(ascii("ARLRR"));builder.writeProfile(raw);
		check(empty(raw) && builder.fork().matchesGeometry(512,GEOMETRY16,0,ProteinPositionMode.FIXED,4,2,official.sha80()),"Construction and fork retain raw filter");
		reject(()->DDLPentamerBlacklist.loadStructural(resource16,policy.sha80,1024,5,
			DDLPentamerBlacklist.MatchMode.PREFIX,policy.provenance),"requested policy");
		final Path headerless16=headerless(root,"headerless16.faa",resource16);
		final DDLPentamerBlacklist wrongMode=DDLPentamerBlacklist.loadStructural(headerless16,policy.sha80,1024,5,
			DDLPentamerBlacklist.MatchMode.PREFIX,policy.provenance);
		reject(()->ProteinBlacklistPolicy.resolve(GEOMETRY16,wrongMode),policy.label);
		final Path wrongWords=root.resolve("wrong_words.faa");
		write(wrongWords,text(resource16).replaceFirst("ARLRR","ARMRR"));
		final DDLPentamerBlacklist forged=policy.load(wrongWords,true);
		reject(()->ProteinBlacklistPolicy.resolve(GEOMETRY16,forged),"keys differ");
		final DDLPentamerBlacklist copy12=policy12.load(resource12,true),copy16=policy.load(resource16,true);
		final DDLPentamerBlacklist structural12=ProteinBlacklistPolicy.resolve(GEOMETRY,copy12);
		check(structural12==ProteinBlacklistPolicy.resolve(GEOMETRY,copy12) &&
			structural12!=ProteinBlacklistPolicy.resolve(GEOMETRY,null),"Prime12K6 structural validation cache");
		final HashMap<String,String> residentOptions=new HashMap<String,String>();
		residentOptions.put("k","5");residentOptions.put("aabits",alphabet(GEOMETRY16));
		FamilyBlacklistBinding.applyDefaults(residentOptions);
		FamilyBlacklistBinding.verifyLoaded(residentOptions,copy16,VerificationMode.STRUCTURAL);checks++;
		reject(()->FamilyBlacklistBinding.verifyLoaded(residentOptions,forged,VerificationMode.STRUCTURAL),"keys differ");
		final DDLPentamerBlacklist structural16=ProteinBlacklistPolicy.resolve(GEOMETRY16,copy16);
		check(structural16==ProteinBlacklistPolicy.resolve(GEOMETRY16,copy16) && structural16!=official &&
			!structural16.identityMeasured() && structural16.sha80().equals(policy.sha80),"Prime16K5 structural validation cache");
		reject(()->ProteinBlacklistPolicy.resolve(GEOMETRY16,copy12),policy.label);
		reject(()->ProteinBlacklistPolicy.resolve(GEOMETRY,copy16),policy12.label);
		check(ProteinBlacklistPolicy.resolve(GEOMETRY16,copy16)==structural16,"Rejected cross-policy use does not poison valid structural cache");
		reject(()->ProteinBlacklistPolicy.requireStored(GEOMETRY16,policy12.sha80),policy.label);
		configuration(official,policy);profiles(root.resolve("profiles"),official,policy,"ARMRR");cluster(root.resolve("cluster"),policy);
	}

	/** Default materialization is explicit; bad historical or partial policy cannot be silently accepted. */
	private static void configuration(DDLPentamerBlacklist official,Descriptor policy)throws Exception{
		final HashMap<String,String> options=new HashMap<String,String>();
		options.put("k",Integer.toString(policy.filterK));options.put("aabits",alphabet(policy.geometry));
		reject(()->FamilyBlacklistBinding.requireOfficialOptions(options),policy.label);
		FamilyBlacklistBinding.applyDefaults(options);
		check(options.get("blacklistmode").equals(policy.mode.label()) && options.get("blacklistk").equals(Integer.toString(policy.filterK)) &&
			options.get("blacklistwords").equals(Integer.toString(policy.words)) && options.get("blacklistsha80").equals(policy.sha80) &&
			options.get("blacklistprovenance").equals(policy.provenance),"Defaults exactly materialize selected descriptor");
		final DDLPentamerBlacklist structural=FamilyBlacklistBinding.load(options,VerificationMode.STRUCTURAL);
		check(!structural.identityMeasured() && structural.sha80().equals(policy.sha80) &&
			ProteinBlacklistPolicy.resolve(policy.geometry,structural)==structural,"File-backed preparation retains structural official singleton");
		final HashMap<String,String> alias=new HashMap<String,String>();alias.put("k",Integer.toString(policy.filterK));
		alias.put("aapartition",alphabet(policy.geometry));FamilyBlacklistBinding.applyDefaults(alias);
		check(alias.get("blacklistsha80").equals(policy.sha80),"Alphabet argument alias selects identical policy");
		for(final String[] edit:new String[][]{{"blacklistwords","8192"},{"blacklistsha80","unmeasured"},{"blacklistk","4"},
			{"blacklistmode",policy.mode==DDLPentamerBlacklist.MatchMode.LITERAL ? "reduced_literal" : "literal"},
			{"blacklistprovenance","different"}}){
			final HashMap<String,String> bad=new HashMap<String,String>(options);bad.put(edit[0],edit[1]);
			reject(()->FamilyBlacklistBinding.applyDefaults(bad),policy.label);
		}
		final HashMap<String,String> partial=new HashMap<String,String>();partial.put("k",Integer.toString(policy.filterK));
		partial.put("aabits",alphabet(policy.geometry));partial.put("blacklistmode","prefix");
		reject(()->FamilyBlacklistBinding.applyDefaults(partial),"Partial blacklist policy");
		reject(()->new ClusterSketchCodec.Header(policy.geometry,512,0,ProteinPositionMode.FIXED,4,2,1,null),policy.label);
		new ClusterSketchCodec.Header(policy.geometry,512,0,ProteinPositionMode.FIXED,4,2,1,official.sha80());checks++;
		ProteinBlacklistPolicy.requireStored(policy.geometry,policy.sha80);checks++;
		reject(()->ProteinBlacklistPolicy.requireStored(policy.geometry,"unmeasured"),policy.label);
	}

	/** Reload a real profile table through structural attribution, then reject a historical unfiltered summary. */
	private static void profiles(Path root,DDLPentamerBlacklist official,Descriptor policy,String survivor)throws Exception{
		Files.createDirectory(root);
		final char[] row=new char[512];DDLFamilyProfiles.sketch(ascii(survivor),policy.geometry,0,row,null);
		final StringBuilder values=new StringBuilder(),votes=new StringBuilder();
		for(int i=0;i<row.length;i++){
			if(i>0){values.append(',');votes.append(',');}
			values.append((int)row[i]);votes.append(row[i]==0 ? 0 : 1);
		}
		write(root.resolve("profiles.tsv"),DDLFamilyProfiles.HEADER16+"\ntail\tf0\t512\t1\t"+policy.filterK+"\t1\tnonempty\t"+values+"\t"+votes+"\t"+votes+"\n");
		final String header="namespace\tk\tseed\texponent\tsketch\tgeometry";
		final String data="tail\t"+policy.filterK+"\t0\t6\t"+ProteinPositionMode.PLAIN.sketch+"\t"+policy.geometry.name();
		write(root.resolve("summary.tsv"),header+"\tblacklist_sha80\tblacklist_words\t"+official.metadataColumns()+"\n"+
			data+"\t"+official.sha80()+"\t"+policy.words+"\t"+official.metadataValues()+"\n");
		write(root.resolve("list.tsv"),root.resolve("profiles.tsv")+"\tunmeasured\tunmeasured\n");
		final DDLFamilySearch.Index loaded=DDLFamilySearch.loadForAttributionStructural(root.resolve("list.tsv"),512,policy.geometry,0,
			ProteinPositionMode.PLAIN.sketch,1,0,official,new DDLFamilySearch.StructuralInputs());
		check(loaded.families()==1 && Arrays.equals(row,loaded.profileCopy(0)),"Persisted profile/query parity");
		write(root.resolve("summary.tsv"),header+"\n"+data+"\n");
		reject(()->DDLFamilySearch.loadForAttributionStructural(root.resolve("list.tsv"),512,policy.geometry,0,
			ProteinPositionMode.PLAIN.sketch,1,0,new DDLFamilySearch.StructuralInputs()),"identity or size mismatch");
	}

	/** The real clustering entry point must emit the official filter without an explicit blacklist flag. */
	private static void cluster(Path root,Descriptor policy)throws Exception{
		Files.createDirectory(root);final Path input=root.resolve("input.faa");
		final StringBuilder sequence=new StringBuilder();for(int i=0;i<6;i++){sequence.append("ACDEFGHIKLMNPQRSTVWY");}
		write(input,">a\n"+sequence+"\n>b\n"+sequence+"\n");
		DdlClust.main(new String[]{"in="+input,"out=null","members="+root.resolve("members.tsv"),"stats="+root.resolve("stats.tsv"),
			"ddlout="+root.resolve("families.pddl"),"consensus=null","refinement=null","quality=null","hist=null","pairhist=null",
			"checkpoints=null","selectedseeds=null","k="+policy.filterK,"aabits="+alphabet(policy.geometry),"width=512","partitions=4","overlap=2",
			"t=1","memberranking=f","passes=1","consensuspasses=1","mincluster=1","topn=0","merge=f","prunefinal=f","verification=structural"});
		final ClusterSketchCodec.Loaded stored=ClusterSketchCodec.read(root.resolve("families.pddl").toString());
		check(policy.sha80.equals(stored.header.blacklistSha80),"Clustering emits official filter identity");
		final String stats=new String(Files.readAllBytes(root.resolve("stats.tsv")),StandardCharsets.US_ASCII);
		check(stats.contains("blacklist_sha80\t"+policy.sha80) && stats.contains("blacklist_words\t"+policy.words),"Clustering records effective filter");
	}

	private static String alphabet(ProteinKmerGeometry geometry){
		check(geometry.isReduced(),"Fixture uses a full-word alphabet");
		return "all:"+geometry.name().substring(geometry.name().indexOf('+')+1);
	}
	private static boolean empty(char[] row){for(char c:row){if(c!=0){return false;}}return true;}
	private static byte[] ascii(String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static void write(Path path,String text)throws Exception{Files.write(path,ascii(text));}
	private static Path headerless(Path root,String name,Path resource)throws Exception{
		final Path path=root.resolve(name);
		write(path,text(resource).replaceFirst("(?m)^>[^\\n]*",">headerless"));
		return path;
	}
	private static String text(Path path){
		final StringBuilder sb=new StringBuilder();
		for(String line:readLines(path)){sb.append(line).append('\n');}
		return sb.toString();
	}
	private static ArrayList<String> readLines(Path path){
		final ByteFile in=ByteFile.makeByteFile(path.toString(),false);
		final ArrayList<String> lines=new ArrayList<String>();
		for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
			lines.add(new String(line,StandardCharsets.US_ASCII));
		}
		if(in.close()){throw new RuntimeException("Error reading or closing "+path);}
		return lines;
	}
	private static void check(boolean valid,String why){checks++;if(!valid){throw new AssertionError(why);}}
	private static void reject(Action action,String reason)throws Exception{
		try{action.run();}catch(IllegalArgumentException expected){
			check(expected.getMessage()!=null && expected.getMessage().contains(reason),"Expected rejection: "+reason+"; got "+expected);return;
		}
		throw new AssertionError("Expected rejection: "+reason);
	}
	private interface Action{void run()throws Exception;}
	private static final ProteinKmerGeometry GEOMETRY=ProteinBlacklistPolicy.GEOMETRY;
	private static final ProteinKmerGeometry GEOMETRY16=ProteinKmerGeometry.make(5,"all:A/C/D/E/FWY/G/H/IV/K/LM/N/P/Q/R/S/T");
	private static int checks;
}
