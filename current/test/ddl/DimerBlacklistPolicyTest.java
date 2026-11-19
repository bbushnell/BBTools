package test.ddl;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import ddl.DDLFamilyProfiles;
import ddl.DDLFamilySearch;
import ddl.DDLPentamerBlacklist;
import ddl.DDLQueryOrigins;
import ddl.FamilySketchBuilder;
import ddl.ProteinBlacklistPolicy;
import ddl.ProteinBlacklistPolicy.Descriptor;
import ddl.ProteinKmerGeometry;
import ddl.ProteinPositionMode;
import dna.Data;
import fileIO.ByteFile;
import prot.FamilyBlacklistBinding;
import prot.FamilyCalibration.VerificationMode;
import prot.cluster.ClusterSketchCodec;
import shared.Shared;

/** Selected d120 resource, independent pair-coordinate oracle and file-backed policy regressions.
 * Accepted resource bytes retain measured verification; synthetic artifacts use structural verification.
 * @author Qiqi
 */
public final class DimerBlacklistPolicyTest{

	public static void main(String[] args)throws Exception{
		check(args.length==3, "Require fresh fixture directory, selected d120 resource and frozen candidate config");
		final Path root=Paths.get(args[0]), resource=Paths.get(args[1]), config=Paths.get(args[2]);
		check(root.isAbsolute() && !Files.exists(root) && resource.isAbsolute() && config.isAbsolute(), "Explicit fresh native fixture paths");
		Files.createDirectory(root); Shared.setThreads(1); Shared.AMINO_IN=true;
		Data.DORI_RESOURCES=resource.getParent().getParent().toString()+"/";
		final String alphabet=alphabet(config);
		check(alphabet.equals(ProteinBlacklistPolicy.DIMER120_ALPHABET), "Production mapping equals the independently frozen candidate mapping");
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(6, alphabet);
		final Descriptor policy=ProteinBlacklistPolicy.policy(geometry);
		check(policy!=null && policy.geometry.equals(geometry) && policy.filterK==6 && policy.words==2048 &&
			policy.mode==DDLPentamerBlacklist.MatchMode.REDUCED_LITERAL && policy.sha80.equals(IDENTITY) &&
			policy.provenance.equals(PROVENANCE), "Exact selected d120 descriptor and provenance");
		check(resource.getFileName().toString().equals(policy.resourceName) &&
			Arrays.equals(Files.readAllBytes(resource), Files.readAllBytes(policy.resource())), "Canonical resolved resource bytes");
		final DDLPentamerBlacklist official=ProteinBlacklistPolicy.resolve(geometry, null);
		check(official.size()==2048 && official.wordLength()==6 && official.matchMode()==policy.mode &&
			official.sha80().equals(IDENTITY) && official.provenance().equals(PROVENANCE), "Dedicated d120 official holder");
		final DDLPentamerBlacklist pinned=policy.load(resource, false), structural=policy.load(resource, true);
		check(ProteinBlacklistPolicy.resolve(geometry, pinned)==pinned &&
			ProteinBlacklistPolicy.resolve(geometry, structural)==structural, "Supplied selected d120 filters validate without singleton replacement");
		final String[][] groups=groups(alphabet);
		final int[] pairClasses=pairClasses(groups);
		final List<String> lines=readLines(resource);
		final boolean[] excluded=new boolean[1<<21];
		final char[] raw=new char[512], encoded=new char[512], positional=new char[512];
		int records=0, aliases=0; String first=null;
		for(String word:lines){
			if(word.startsWith(">")){continue;}
			check(word.length()==6, "Each selected representative spans three disjoint pairs");
			final int key=key(word, pairClasses);
			check(!excluded[key], "Independent pair-coordinate oracle has no duplicate selected key"); excluded[key]=true;
			check(official.rawWord(geometry, pack(word))==key && official.contains(key), "Native projection agrees with literal ordered-pair oracle");
			check(DDLFamilyProfiles.sketch(ascii(word), geometry, 0, raw, null)==1 && empty(raw), "Every selected word is excluded by the default ASCII kernel");
			final byte[] codes=new byte[6];
			for(int i=0; i<6; i++){codes[i]=(byte)ACIDS.indexOf(word.charAt(i));}
			DDLFamilyProfiles.sketchEncoded(codes, geometry, 0, encoded, null);
			check(Arrays.equals(raw, encoded), "Encoded and ASCII selected filtering agree");
			DDLFamilyProfiles.sketchPositional(ascii(word), geometry, 0, 4, 2, positional, null);
			check(empty(positional), "Selected filtering occurs before positional tags");
			final String alias=alias(word, groups, pairClasses);
			if(alias!=null){
				check(key(alias, pairClasses)==key && !alias.equals(word), "Different raw representative retains the independently computed key");
				DDLFamilyProfiles.sketch(ascii(alias), geometry, 0, raw, null);
				check(empty(raw), "Selected dimer aliases are excluded in reduced coordinates"); aliases++;
			}
			if(first==null){first=word;} records++;
		}
		check(records==2048 && aliases>0, "Complete selected list and real pair aliases exercised");
		final String survivor=survivor(groups, excluded);
		DDLFamilyProfiles.sketch(ascii(survivor), geometry, 0, raw, null);
		check(!empty(raw), "An independently unlisted dimer key survives");
		DDLFamilyProfiles.sketchPositional(ascii(first), geometry, 0, 4, 2, positional, null, new DDLQueryOrigins(512), ProteinPositionMode.FIXED);
		check(empty(positional), "Origin tracing uses the same selected exclusion");
		final FamilySketchBuilder builder=new FamilySketchBuilder(512, geometry, 0, ProteinPositionMode.FIXED, 4, 2, null);
		builder.addRaw(ascii(first)); builder.writeProfile(raw);
		check(empty(raw) && builder.fork().matchesGeometry(512, geometry, 0, ProteinPositionMode.FIXED, 4, 2, IDENTITY), "Member construction and forks retain selected identity");
		configuration(root, policy, official, alphabet);
		negatives(root, resource, policy, official, alphabet, first, survivor);
		profiles(root.resolve("profiles"), policy, official, survivor);
		System.out.println("DIMER_BLACKLIST_POLICY_TEST_PASS checks="+checks);
	}

	/** Reads the frozen scientific input rather than deriving the expected map from production constants. */
	private static String alphabet(Path config)throws Exception{
		String answer=null;
		for(String line:Files.readAllLines(config, StandardCharsets.US_ASCII)){
			if(line.startsWith("alphabet=")){check(answer==null, "Exactly one candidate alphabet"); answer=line.substring(9);}
		}
		check(answer!=null && answer.startsWith("dimer:1,3,5:"), "Frozen nonoverlapping K6 candidate alphabet");
		return answer;
	}

	/** Builds an independent literal pair-to-class table from all400 canonical pair names. */
	private static String[][] groups(String alphabet){
		final String[] text=alphabet.substring("dimer:1,3,5:".length()).split("/", -1);
		check(text.length==120, "Exactly120 ordered pair classes");
		final String[][] answer=new String[text.length][];
		for(int i=0; i<text.length; i++){answer[i]=text[i].split(",", -1);}
		return answer;
	}
	/** Maps every canonical ordered pair to its unique class in the frozen partition. */
	private static int[] pairClasses(String[][] groups){
		final int[] answer=new int[400]; Arrays.fill(answer, -1); int count=0;
		for(int group=0; group<groups.length; group++){
			for(String pair:groups[group]){
				check(pair.length()==2 && ACIDS.indexOf(pair.charAt(0))>=0 && ACIDS.indexOf(pair.charAt(1))>=0, "Canonical ordered pair");
				final int index=20*ACIDS.indexOf(pair.charAt(0))+ACIDS.indexOf(pair.charAt(1));
				check(answer[index]<0, "Each canonical pair occurs once"); answer[index]=group; count++;
			}
		}
		check(count==400, "Complete canonical pair partition"); return answer;
	}
	/** Packs the three nonoverlapping pair classes into independent seven-bit coordinates. */
	private static int key(String word, int[] classes){
		int key=0;
		for(int i=0; i<6; i+=2){key=(key<<7)|classes[20*ACIDS.indexOf(word.charAt(i))+ACIDS.indexOf(word.charAt(i+1))];}
		return key;
	}
	/** Substitutes a distinct pair from the same class without changing the reduced key. */
	private static String alias(String word, String[][] groups, int[] classes){
		for(int i=0; i<6; i+=2){
			final int group=classes[20*ACIDS.indexOf(word.charAt(i))+ACIDS.indexOf(word.charAt(i+1))];
			for(String pair:groups[group]){
				if(!pair.equals(word.substring(i, i+2))){return word.substring(0, i)+pair+word.substring(i+2);}
			}
		}
		return null;
	}
	/** Finds a valid unselected key without consulting the production filter or projection. */
	private static String survivor(String[][] groups, boolean[] excluded){
		for(int a=0; a<120; a++){
			for(int b=0; b<120; b++){
				for(int c=0; c<120; c++){if(!excluded[(a<<14)|(b<<7)|c]){return groups[a][0]+groups[b][0]+groups[c][0];}}
			}
		}
		throw new AssertionError("No unselected dimer key");
	}

	/** Checks actual file-backed defaults, explicit provenance, stale metadata and independent official holders. */
	private static void configuration(Path root, Descriptor policy, DDLPentamerBlacklist official, String alphabet)throws Exception{
		final HashMap<String,String> options=new HashMap<String,String>(); options.put("k", "6"); options.put("aabits", alphabet);
		reject(()->FamilyBlacklistBinding.requireOfficialOptions(options), policy.label);
		FamilyBlacklistBinding.applyDefaults(options);
		check(options.get("blacklistsha80").equals(IDENTITY) && options.get("blacklistwords").equals("2048") &&
			options.get("blacklistmode").equals("reduced_literal") && options.get("blacklistprovenance").equals(PROVENANCE), "Exact selected file-backed defaults");
		final DDLPentamerBlacklist pinned=FamilyBlacklistBinding.load(options, VerificationMode.SHA80);
		final DDLPentamerBlacklist structural=FamilyBlacklistBinding.load(options, VerificationMode.STRUCTURAL);
		check(pinned!=official && pinned.sha80().equals(IDENTITY) && pinned.size()==policy.words && pinned.matchMode()==policy.mode &&
			ProteinBlacklistPolicy.resolve(policy.geometry, pinned)==pinned && ProteinBlacklistPolicy.resolve(policy.geometry, structural)==structural,
			"File-backed SHA80 verifies shipped d120 bytes while structural d120 validates supplied filters");
		final HashMap<String,String> alias=new HashMap<String,String>(); alias.put("k", "6"); alias.put("aapartition", alphabet);
		FamilyBlacklistBinding.applyDefaults(alias); check(alias.equals(aliasOptions(options, alphabet)), "Alphabet argument alias has identical effective policy");
		for(String[] edit:new String[][]{{"blacklistsha80", ZERO}, {"blacklistwords", "2047"}, {"blacklistk", "5"},
			{"blacklistmode", "literal"}, {"blacklistprovenance", "different"}}){
			final HashMap<String,String> bad=new HashMap<String,String>(options); bad.put(edit[0], edit[1]);
			reject(()->FamilyBlacklistBinding.load(bad, VerificationMode.STRUCTURAL), policy.label);
		}
		final DDLPentamerBlacklist wrongProvenance=DDLPentamerBlacklist.loadReducedStructural(
			headerless(root, "wrong_provenance.faa", policy.resource()), IDENTITY, 2048, policy.geometry, "different");
		reject(()->FamilyBlacklistBinding.verifyLoaded(options, wrongProvenance, VerificationMode.STRUCTURAL), "provenance");
		for(String identity:new String[]{null, ZERO, ProteinBlacklistPolicy.SHA80, "f731f9dfe82458b93234"}){
			reject(()->ProteinBlacklistPolicy.requireStored(policy.geometry, identity), policy.label);
			reject(()->new ClusterSketchCodec.Header(policy.geometry, 512, 0, ProteinPositionMode.FIXED, 4, 2, 1, identity), policy.label);
		}
		new ClusterSketchCodec.Header(policy.geometry, 512, 0, ProteinPositionMode.FIXED, 4, 2, 1, IDENTITY); checks++;
		final DDLPentamerBlacklist monomer=ProteinBlacklistPolicy.resolve(ProteinBlacklistPolicy.GEOMETRY, null);
		check(monomer!=official && monomer.sha80().equals(ProteinBlacklistPolicy.SHA80), "Dimer holder does not replace legacy12K6 holder");
		final ProteinKmerGeometry g16=ProteinKmerGeometry.make(5, "all:A/C/D/E/FWY/G/H/IV/K/LM/N/P/Q/R/S/T");
		check(ProteinBlacklistPolicy.resolve(g16, null).sha80().equals("f731f9dfe82458b93234") &&
			ProteinBlacklistPolicy.resolve(policy.geometry, null)==official, "Interleaved16K5 and d120 lookups retain independent holders");
	}
	private static HashMap<String,String> aliasOptions(HashMap<String,String> options, String alphabet){
		final HashMap<String,String> answer=new HashMap<String,String>(options); answer.remove("aabits"); answer.put("aapartition", alphabet); return answer;
	}

	/** Rejects wrong coordinates, cardinality and matching mode while preserving supplied same-count content. */
	private static void negatives(Path root, Path resource, Descriptor policy, DDLPentamerBlacklist official,
			String alphabet, String first, String survivor)throws Exception{
		final String groups=alphabet.substring("dimer:1,3,5:".length());
		final int slash=groups.indexOf('/'), next=groups.indexOf('/', slash+1);
		final String swapped="dimer:1,3,5:"+groups.substring(slash+1, next)+"/"+groups.substring(0, slash)+groups.substring(next);
		final ProteinKmerGeometry alternate=ProteinKmerGeometry.make(6, swapped);
		check(ProteinBlacklistPolicy.policy(alternate)==null && ProteinBlacklistPolicy.resolve(alternate, null)==null, "Different ordered120-map gets no selected default");
		reject(()->official.validate(alternate), "alphabet and geometry");
		final ProteinKmerGeometry overlap=ProteinKmerGeometry.make(6, "dimer:1,2,3,4,5:"+groups);
		final ProteinKmerGeometry k5=ProteinKmerGeometry.make(5, "dimer:1,3,4:"+groups);
		for(ProteinKmerGeometry other:new ProteinKmerGeometry[]{overlap, k5, ProteinKmerGeometry.exact(6)}){
			check(ProteinBlacklistPolicy.policy(other)==null, "Different layout or K gets no selected default");
			reject(()->DDLPentamerBlacklist.loadReduced(resource, IDENTITY, 2048, other, PROVENANCE), "Pinned reduced loading");
			reject(()->DDLPentamerBlacklist.loadReducedStructural(resource, IDENTITY, 2048, other, PROVENANCE), "Reduced blacklist loading");
		}
		reject(()->DDLPentamerBlacklist.loadReducedStructural(resource, IDENTITY, 2047, policy.geometry, PROVENANCE), "requested policy");
		final Path headerless=headerless(root, "headerless.faa", resource);
		reject(()->DDLPentamerBlacklist.loadReducedStructural(headerless, IDENTITY, 2047, policy.geometry, PROVENANCE), "word count");
		final DDLPentamerBlacklist wrongIdentity=DDLPentamerBlacklist.loadReducedStructural(resource, ZERO, 2048, policy.geometry, PROVENANCE);
		reject(()->ProteinBlacklistPolicy.resolve(policy.geometry, wrongIdentity), policy.label);
		reject(()->DDLPentamerBlacklist.loadStructural(resource, IDENTITY, 2048, 6, DDLPentamerBlacklist.MatchMode.LITERAL, PROVENANCE),
			"requested policy");
		final DDLPentamerBlacklist raw=DDLPentamerBlacklist.loadStructural(headerless, IDENTITY, 2048, 6, DDLPentamerBlacklist.MatchMode.LITERAL, PROVENANCE);
		reject(()->ProteinBlacklistPolicy.resolve(policy.geometry, raw), policy.label);
		final String original=text(resource);
		final Path forged=root.resolve("forged.faa"); write(forged, original.replaceFirst("(?m)^"+first+"$", survivor));
		final DDLPentamerBlacklist changed=policy.load(forged, true);
		check(ProteinBlacklistPolicy.resolve(policy.geometry, changed)==changed && changed!=official,
			"Supplied selected d120 filters retain metadata-only structural validation");
		final Path duplicate=root.resolve("duplicate.faa"); write(duplicate, original.replaceFirst("words=2048","words=2049")+">duplicate\n"+first+"\n");
		reject(()->DDLPentamerBlacklist.loadReducedStructural(duplicate, IDENTITY, 2049, policy.geometry, PROVENANCE), "Duplicate blacklist word");
		check(ProteinBlacklistPolicy.resolve(policy.geometry, official)==official, "Rejected malformed inputs do not alter valid supplied d120 filters");
	}

	/** Loads a stored profile with complete d120 metadata, then rejects its stale unfiltered summary. */
	private static void profiles(Path root, Descriptor policy, DDLPentamerBlacklist official, String survivor)throws Exception{
		Files.createDirectory(root); final char[] row=new char[512];
		DDLFamilyProfiles.sketch(ascii(survivor), policy.geometry, 0, row, null);
		final StringBuilder values=new StringBuilder(), votes=new StringBuilder();
		for(int i=0; i<row.length; i++){
			if(i>0){values.append(','); votes.append(',');}
			values.append((int)row[i]); votes.append(row[i]==0 ? 0 : 1);
		}
		write(root.resolve("profiles.tsv"), DDLFamilyProfiles.HEADER16+"\ntail\tf0\t512\t1\t6\t1\tnonempty\t"+values+"\t"+votes+"\t"+votes+"\n");
		final String header="namespace\tk\tseed\texponent\tsketch\tgeometry";
		final String data="tail\t6\t0\t6\t"+ProteinPositionMode.PLAIN.sketch+"\t"+policy.geometry.name();
		write(root.resolve("summary.tsv"), header+"\tblacklist_sha80\tblacklist_words\t"+official.metadataColumns()+"\n"+
			data+"\t"+IDENTITY+"\t2048\t"+official.metadataValues()+"\n");
		write(root.resolve("list.tsv"), root.resolve("profiles.tsv")+"\tunmeasured\tunmeasured\n");
		final DDLFamilySearch.Index loaded=DDLFamilySearch.loadForAttributionStructural(root.resolve("list.tsv"), 512, policy.geometry, 0,
			ProteinPositionMode.PLAIN.sketch, 1, 0, official, new DDLFamilySearch.StructuralInputs());
		check(loaded.families()==1 && Arrays.equals(row, loaded.profileCopy(0)), "Stored d120 profile and query filtering agree");
		write(root.resolve("summary.tsv"), header+"\n"+data+"\n");
		reject(()->DDLFamilySearch.loadForAttributionStructural(root.resolve("list.tsv"), 512, policy.geometry, 0,
			ProteinPositionMode.PLAIN.sketch, 1, 0, new DDLFamilySearch.StructuralInputs()), "identity or size mismatch");
	}

	/** Packs raw canonical residues independently of the production blacklist helper. */
	private static long pack(String word){
		long answer=0;
		for(int i=0; i<word.length(); i++){
			answer=(answer<<5)|ACIDS.indexOf(word.charAt(i));
		}
		return answer;
	}
	private static boolean empty(char[] row){for(char c:row){if(c!=0){return false;}}return true;}
	private static byte[] ascii(String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static void write(Path file, String text)throws Exception{Files.write(file, ascii(text));}
	private static Path headerless(Path root, String name, Path resource)throws Exception{
		final Path path=root.resolve(name);
		write(path, text(resource).replaceFirst("(?m)^>[^\\n]*", ">headerless"));
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
	private static void check(boolean valid, String why){checks++; if(!valid){throw new AssertionError(why);}}
	private static void reject(Action action, String reason)throws Exception{
		try{action.run();}catch(IllegalArgumentException expected){
			check(expected.getMessage()!=null && expected.getMessage().contains(reason), "Expected "+reason+"; got "+expected); return;
		}
		throw new AssertionError("Expected "+reason);
	}
	private interface Action{void run()throws Exception;}
	private static final String ACIDS="ARNDCEQGHILKMFPSTWYV", IDENTITY="c0fb8ec0a91bb50535c5", ZERO="00000000000000000000";
	private static final String PROVENANCE="official_d120_k6_nonoverlap_shadow2048;selector=8126cccd270e71d3ead11e2d7b2f6416702bad0f";
	private static int checks;
}
