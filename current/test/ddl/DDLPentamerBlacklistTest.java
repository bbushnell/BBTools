package test.ddl;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.Arrays;
import java.util.HashMap;

import ddl.DDLPentamerBlacklist;
import ddl.DDLFamilyProfiles;
import ddl.DDLFamilySearch;
import ddl.DDLFamilyHoldout;
import ddl.DDLKeyFrequency;
import ddl.ProteinKmerGeometry;
import prot.FamilyCalibration;
import shared.KillSwitch;
import shared.Shared;
import stream.FastaStreamerST;

/** Independent raw-AA packing and native FASTA input checks for DDL exclusions.
 * Creates and retains tiny temporary fixture directories. Run once in an isolated
 * JVM with assertions enabled: the static check count accumulates, ordinary mode
 * sets Shared threads to one, and lifecycle assertions halt on stranded owned readers.
 * These are software fixtures, not biological calibration or performance panels.
 * @author Fischl
 */
public final class DDLPentamerBlacklistTest{

	/** Runs all checks by default, or selects structural/lifecycle-only validation.
	 * Structural mode avoids digest checks; both default and structural modes include lifecycle.
	 * @param args Empty, one literal structural, or one literal lifecycle argument */
	public static void main(String[] args) throws Exception{
		if(args.length==1 && args[0].equals("lifecycle")){
			lifecycle();
			System.out.println("DDL_PENTAMER_BLACKLIST_LIFECYCLE_TEST_PASS checks="+checks);
			return;
		}
		if(args.length==1 && args[0].equals("structural")){
			structuralOnly();
			return;
		}
		if(args.length!=0){throw new IllegalArgumentException("Expected no arguments, structural or lifecycle");}
		Shared.setThreads(1);
		final String alphabet="ARNDCEQGHILKMFPSTWYV";
		for(int i=0; i<alphabet.length(); i++){
			check(pack("AAAA"+alphabet.charAt(i))==i, "Last residue uses the native AA code");
			check(pack("aaaa"+Character.toLowerCase(alphabet.charAt(i)))==i, "Case does not change the word");
		}
		check(pack("ACDEF")==((4L<<15)|(3L<<10)|(5L<<5)|13), "Forward five-bit word, not nucleotide canonicalization");
		expectFailure(()->pack("AAAA"));
		expectFailure(()->pack("AAAAAA"));
		expectFailure(()->pack("AAAAX"));
		expectFailure(()->pack("AAAA*"));
		final Path dir=Files.createTempDirectory("ddl-pentamer-blacklist-");
		final Path good=write(dir, "good.faa", ">zero\nAA\nAAA\n>one\naaaAR\n>acid\nACDEF\n");
		final String original=pin(good);
		final DDLPentamerBlacklist list=DDLPentamerBlacklist.load(good, original);
		check(list.size()==3 && list.contains(0) && list.contains(1) && list.contains(pack("ACDEF")),
			"All three literal words, including raw zero");
		check(!list.contains(2) && original.equals(list.sha80()) && list.wordLength()==5 &&
			list.matchMode()==DDLPentamerBlacklist.MatchMode.PREFIX, "No unlisted word; recipe identity retained");
		sketches(list);
		modes(dir, list);
		dimerStructural(dir);
		metadata(list);
		integration(dir, good, original, list);
		final Path duplicate=write(dir, "duplicate.faa", ">a\nAAAAA\n>b\naaaaa\n");
		expectFailure(()->DDLPentamerBlacklist.load(duplicate, pin(duplicate)));
		final Path shortWord=write(dir, "short.faa", ">a\nAAAA\n");
		expectFailure(()->DDLPentamerBlacklist.load(shortWord, pin(shortWord)));
		final Path ambiguous=write(dir, "ambiguous.faa", ">a\nAAAAX\n");
		expectFailure(()->DDLPentamerBlacklist.load(ambiguous, pin(ambiguous)));
		Files.write(good, ">changed\nAAAAN\n".getBytes(StandardCharsets.US_ASCII));
		expectFailure(()->DDLPentamerBlacklist.load(good, original));
		check(list.size()==3 && list.contains(0) && !list.contains(2),
			"An already loaded set is immutable despite later file changes");
		lifecycle();
		System.out.println("DDL_PENTAMER_BLACKLIST_TEST_PASS checks="+checks);
	}

	/*--------------------------------------------------------------*/
	/*----------------      Structural loading      ----------------*/
	/*--------------------------------------------------------------*/

	/** Structural fixture path: no SHA-256 reads are needed for unmeasured filter validation. */
	private static void structuralOnly() throws Exception{
		final Path dir=Files.createTempDirectory("ddl-pentamer-blacklist-structural-");
		final Path six=write(dir, "six.faa", ">six\nACDEFG\n");
		final DDLPentamerBlacklist literal6=DDLPentamerBlacklist.loadStructural(six, "unmeasured", 1, 6,
			DDLPentamerBlacklist.MatchMode.LITERAL, "shadow_k6");
		check(literal6.size()==1 && literal6.contains(pack("ACDEFG", 6)) && literal6.provenance().equals("shadow_k6"),
			"Structural unmeasured loading parses the literal K6 list without computing a digest");
		expectFailure(()->DDLPentamerBlacklist.loadStructural(dir, "unmeasured", 1, 6,
			DDLPentamerBlacklist.MatchMode.LITERAL, "shadow_k6"));
		final Path link=dir.resolve("link.faa");
		Files.createSymbolicLink(link, six);
		expectFailure(()->DDLPentamerBlacklist.loadStructural(link, "unmeasured", 1, 6,
			DDLPentamerBlacklist.MatchMode.LITERAL, "shadow_k6"));
		expectFailure(()->DDLPentamerBlacklist.loadStructural(six, "unmeasured", 1, 6,
			DDLPentamerBlacklist.MatchMode.LITERAL, ""));
		expectFailure(()->DDLPentamerBlacklist.loadStructural(six, "unmeasured", 2, 6,
			DDLPentamerBlacklist.MatchMode.LITERAL, "shadow_k6"));
		expectFailure(()->DDLPentamerBlacklist.loadStructural(six, "unmeasured", 1, 7,
			DDLPentamerBlacklist.MatchMode.LITERAL, "shadow_k6"));
		check(DDLPentamerBlacklist.MatchMode.parse("literal")==DDLPentamerBlacklist.MatchMode.LITERAL &&
			DDLPentamerBlacklist.MatchMode.parse("prefix")==DDLPentamerBlacklist.MatchMode.PREFIX,
			"Mode parser accepts the explicit length-free literal and prefix spellings");
		expectFailure(()->DDLPentamerBlacklist.MatchMode.parse("literalK6"));
		expectFailure(()->DDLPentamerBlacklist.MatchMode.parse("prefixK5"));
		final Path good=write(dir, "good.faa", ">zero\nAA\nAAA\n>one\naaaAR\n>acid\nACDEF\n");
		final DDLPentamerBlacklist prefix5=DDLPentamerBlacklist.loadStructural(good, "unmeasured", 3, 5,
			DDLPentamerBlacklist.MatchMode.PREFIX, "real_k5");
		final DDLPentamerBlacklist literal5=DDLPentamerBlacklist.loadStructural(good, "unmeasured", 3, 5,
			DDLPentamerBlacklist.MatchMode.LITERAL, "literal5");
		final byte[] sequence="ACDEFGACDEFH".getBytes(StandardCharsets.US_ASCII);
		final char[] actual=new char[1024], expectedPrefix=new char[1024], expectedLiteral=new char[1024];
		addK6Window(sequence, 1, expectedPrefix);
		addK6Window(sequence, 2, expectedPrefix);
		addK6Window(sequence, 3, expectedPrefix);
		addK6Window(sequence, 4, expectedPrefix);
		addK6Window(sequence, 5, expectedPrefix);
		for(int start=1; start+6<=sequence.length; start++){addK6Window(sequence, start, expectedLiteral);}
		DDLFamilyProfiles.sketch(sequence, ProteinKmerGeometry.exact(6), 0, actual, prefix5);
		check(Arrays.equals(actual, expectedPrefix), "Structural K5 prefix mode filters K6 words by the first five residues");
		DDLFamilyProfiles.sketch(sequence, ProteinKmerGeometry.exact(6), 0, actual, literal6);
		check(Arrays.equals(actual, expectedLiteral), "Structural literal K6 mode filters only an exact six-residue word");
		expectFailure(()->DDLFamilyProfiles.sketch(sequence, ProteinKmerGeometry.exact(6), 0, actual, literal5));
		expectFailure(()->new ddl.FamilySketchBuilder(1024, ProteinKmerGeometry.exact(6), 0, false, 1, 0, literal5));
		final HashMap<String, String> legacyFields=new HashMap<String, String>();
		legacyFields.put("blacklist_sha80", literal5.sha80());
		legacyFields.put("blacklist_words", Integer.toString(literal5.size()));
		expectFailure(()->DDLFamilySearch.validateBlacklistMetadata(legacyFields, literal5));
		dimerStructural(dir);
		lifecycle();
		System.out.println("DDL_PENTAMER_BLACKLIST_STRUCTURAL_TEST_PASS checks="+checks);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Reader lifecycle       ----------------*/
	/*--------------------------------------------------------------*/

	/** Early rejection must release all readers owned by the load, not merely return its diagnostic. */
	private static void lifecycle() throws Exception{
		final int threads=Shared.threads(), buffers=Shared.READ_BUFFER_NUM_BUFFERS();
		final int batchSize=FastaStreamerST.TARGET_LIST_SIZE;
		final boolean simd=Shared.SIMD, lowMemory=Shared.LOW_MEMORY;
		final Path dir=Files.createTempDirectory("ddl-blacklist-lifecycle-");
		final String groups=dimerGroups(), provenance="bbgap_blacklist_k6_2048_20261116_v1";
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(6, "dimer:1,3,5:"+groups);
		final Path good=lifecycleInput(dir, "good.faa", groups, provenance, -1, false);
		final Path duplicate=lifecycleInput(dir, "duplicate.faa", groups, provenance, 1, false);
		final Path shortWord=lifecycleInput(dir, "short.faa", groups, provenance, 17, true);
		try{
			// Select the same bounded-queue scalar reader on both actual JVMs, with many pending batches.
			Shared.setThreads(4);
			Shared.SIMD=false;
			Shared.LOW_MEMORY=false;
			FastaStreamerST.TARGET_LIST_SIZE=8;
			check(2048>3*FastaStreamerST.TARGET_LIST_SIZE,
				"The rejection fixture must exceed the consumer batch plus both producer queue slots");
			for(int repeat=0; repeat<4; repeat++){
				ownedLoad("provenance-"+repeat, ()->expectFailureMessage(
					()->DDLPentamerBlacklist.loadReducedStructural(good, "unmeasured", 2048, geometry,
						"original_component_v1"), "BBGAP blacklist header does not match the requested policy"));
				ownedLoad("duplicate-"+repeat, ()->expectFailureMessage(
					()->DDLPentamerBlacklist.loadReducedStructural(duplicate, "unmeasured", 2048, geometry, provenance),
					"Duplicate blacklist word: b_1"));
				ownedLoad("short-word-"+repeat, ()->expectFailureMessage(
					()->DDLPentamerBlacklist.loadReducedStructural(shortWord, "unmeasured", 2048, geometry, provenance),
					"Each blacklist record must contain exactly 6 amino acids"));
			}
			ownedLoad("valid", ()->{
				final DDLPentamerBlacklist loaded=DDLPentamerBlacklist.loadReducedStructural(
					good, "unmeasured", 2048, geometry, provenance);
				boolean exact=loaded.size()==2048 && loaded.provenance().equals(provenance);
				for(int i=0; i<2048; i++){
					final long key=((long)(i%120)<<14)|((long)((i/120)%120)<<7)|(i/(120*120));
					exact&=loaded.contains(key);
				}
				check(exact, "Complete valid loading preserves all 2048 literal class-triplet keys and provenance");
			});
		}finally{
			FastaStreamerST.TARGET_LIST_SIZE=batchSize;
			Shared.SIMD=simd;
			Shared.LOW_MEMORY=lowMemory;
			Shared.setThreads(threads);
			Shared.setBuffers(buffers);
		}
	}

	/** Keep the real loader and its descendants in a private group so unrelated suite workers are excluded. */
	private static void ownedLoad(String name, Action action) throws Exception{
		final ThreadGroup group=new ThreadGroup("blacklist-lifecycle-"+name);
		final Throwable[] failure=new Throwable[1];
		final Thread load=new Thread(group, ()->{
			try{action.run();}catch(Throwable e){failure[0]=e;}
		}, "blacklist-load-"+name);
		load.start();
		load.join(10000);
		final boolean loadTerminated=!load.isAlive();
		assert(loadTerminated) : KillSwitch.assertDie("Blacklist lifecycle test cannot finish while owned load survives: "+
			load.getName()+" state="+load.getState());
		check(loadTerminated, "Blacklist load did not finish: "+name);
		final Thread[] readers=new Thread[64];
		final int count=group.enumerate(readers, true);
		assert(count<readers.length) : KillSwitch.assertDie("Blacklist lifecycle test cannot account for all owned readers: "+
			name+" group="+group.getName()+" enumerated="+count+" capacity="+readers.length);
		if(count==readers.length){throw new AssertionError("Too many owned blacklist reader threads: "+name);}
		boolean terminated=true;
		final StringBuilder stranded=new StringBuilder();
		for(int i=0; i<count; i++){
			readers[i].join(10000);
			if(readers[i].isAlive()){
				terminated=false;
				stranded.append(readers[i].getName()).append(" state=").append(readers[i].getState()).append(';');
			}
		}
		assert(terminated) : KillSwitch.assertDie("Blacklist lifecycle test cannot finish while owned readers survive: "+
			name+" readers="+stranded);
		check(terminated, "Blacklist rejection stranded an owned reader thread: "+name);
		if(failure[0] instanceof Exception){throw (Exception)failure[0];}
		if(failure[0] instanceof Error){throw (Error)failure[0];}
		if(failure[0]!=null){throw new AssertionError(failure[0]);}
	}

	/** A 2048-word native FASTA exceeds queue capacity; one optional mutation fails before EOF. */
	private static Path lifecycleInput(Path dir, String name, String groups, String provenance,
			int mutation, boolean shortWord) throws Exception{
		final String[] representatives=groups.split("/");
		assert(representatives.length==120) : "The literal lifecycle key oracle assumes 120 ordered dimer classes";
		for(int i=0; i<representatives.length; i++){
			representatives[i]=representatives[i].split(",")[0];
		}
		final StringBuilder text=new StringBuilder(">bbgap_blacklist_v1;k=6;words=2048;mode=reduced_literal;provenance=");
		text.append(provenance).append(";aabits=dimer:1,3,5:").append(groups).append('\n');
		for(int i=0; i<2048; i++){
			if(i>0){text.append(">b_").append(i).append('\n');}
			String word=representatives[i%120]+representatives[(i/120)%120]+representatives[i/(120*120)];
			if(i==mutation){word=shortWord ? "AAAAA" : "AAAAAA";}
			if(i==0){
				text.append(word.substring(0, 3)).append('\n').append(word.substring(3));
			}else{
				text.append(i==2047 ? word.toLowerCase(java.util.Locale.ROOT) : word);
			}
			text.append('\n');
		}
		return write(dir, name, text.toString());
	}

	/** The original validation reason must survive reader cleanup, not become a generic I/O error. */
	private static void expectFailureMessage(Action action, String message) throws Exception{
		try{action.run();}catch(IllegalArgumentException expected){
			check(message.equals(expected.getMessage()), "Original blacklist diagnostic changed: "+expected);
			return;
		}
		throw new AssertionError("Invalid blacklist accepted: "+message);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Dimer and mode checks     ----------------*/
	/*--------------------------------------------------------------*/

	/** Literal pair-class oracle exercises zero and the high end of a120-class packed keyspace. */
	private static void dimerStructural(Path dir)throws Exception{
		final String groups=dimerGroups();
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(6, "dimer:1,3,5:"+groups);
		final Path words=write(dir, "dimer.faa", ">zero\nAAAAAA\n>high\nEVEVEV\n");
		final DDLPentamerBlacklist filter=DDLPentamerBlacklist.loadReducedStructural(words, "unmeasured", 2, geometry,
			"dimer_fixture");
		final long high=(119L<<14)|(119L<<7)|119;
		check(filter.size()==2 && filter.contains(0) && filter.contains(high),
			"Three seven-bit class symbols retain zero and high packed keys");
		check(filter.rawWord(geometry, pack("QAQAQA", 6))==0 && filter.rawWord(geometry, pack("KVKVKV", 6))==high,
			"Distinct literal pairs at raw ordinals0/120 and119/239 share their declared classes");
		check(!filter.contains(filter.rawWord(geometry, pack("ARAAAA", 6))), "An unlisted class triplet remains available");
		final String header="bbgap_blacklist_v1;k=6;words=2;mode=reduced_literal;provenance=dimer_fixture;aabits=dimer:1,3,5:"+groups;
		final Path versioned=write(dir, "dimer_versioned.faa", ">"+header+"\nAAAAAA\n>2\nEVEVEV\n");
		final DDLPentamerBlacklist versionedFilter=DDLPentamerBlacklist.loadReducedStructural(versioned, "unmeasured",
			2, geometry, "dimer_fixture");
		check(versionedFilter.size()==filter.size() && versionedFilter.contains(0) && versionedFilter.contains(high),
			"Versioned FASTA header preserves the exact loaded dimer blacklist keys");
		final Path unknown=write(dir, "dimer_unknown.faa", ">"+header.replace("_v1", "_v2")+"\nAAAAAA\n>2\nEVEVEV\n");
		expectFailure(()->DDLPentamerBlacklist.loadReducedStructural(unknown, "unmeasured", 2, geometry, "dimer_fixture"));
		final Path badCount=write(dir, "dimer_bad_count.faa", ">"+header.replace("words=2", "words=3")+"\nAAAAAA\n>2\nEVEVEV\n");
		expectFailure(()->DDLPentamerBlacklist.loadReducedStructural(badCount, "unmeasured", 2, geometry, "dimer_fixture"));
		final Path badAlphabet=write(dir, "dimer_bad_alphabet.faa",
			">"+header.replace("aabits=dimer:1,3,5:"+groups, "aabits=none")+"\nAAAAAA\n>2\nEVEVEV\n");
		expectFailure(()->DDLPentamerBlacklist.loadReducedStructural(badAlphabet, "unmeasured", 2, geometry, "dimer_fixture"));
		check(filter.matchMode()==DDLPentamerBlacklist.MatchMode.REDUCED_LITERAL && filter.wordLength()==6,
			"Dimer provenance is reduced literal K6, not a raw representative exclusion");
		final char[] actual=new char[512], plain=new char[512], empty=new char[512];
		final byte[] sequence=new byte[100];
		Arrays.fill(sequence, (byte)'A');
		check(DDLFamilyProfiles.sketchPositional(sequence, geometry, 0, 4, 2, actual, filter)==95,
			"Excluded dimer windows still conserve the original window count");
		DDLFamilyProfiles.sketchPositional(sequence, geometry, 0, 4, 2, plain, null);
		check(Arrays.equals(actual, empty) && !Arrays.equals(plain, empty), "Every positional tag excludes the selected class");
		DDLFamilyProfiles.sketch("QAQAQA".getBytes(StandardCharsets.US_ASCII), geometry, 0, actual, filter);
		check(Arrays.equals(actual, empty), "An alternate raw representative is excluded during actual sketching");
		final Path duplicate=write(dir, "dimer_duplicate.faa", ">a\nAAAAAA\n>b\nQAQAQA\n");
		expectFailure(()->DDLPentamerBlacklist.loadReducedStructural(duplicate, "unmeasured", 2, geometry, "fixture"));
		expectFailure(()->DDLPentamerBlacklist.loadReducedStructural(words, "unmeasured", 3, geometry, "fixture"));
		final String[] split=groups.split("/", -1);
		final String first=split[0];
		split[0]=split[1];
		split[1]=first;
		final ProteinKmerGeometry reordered=ProteinKmerGeometry.make(6, "dimer:1,3,5:"+String.join("/", split));
		expectFailure(()->filter.validate(reordered));
		for(String layout:new String[]{"dimer:1,2,3,4,5:", "dimer:1,2,4,5:"}){
			final ProteinKmerGeometry unsupported=ProteinKmerGeometry.make(6, layout+groups);
			expectFailure(()->DDLPentamerBlacklist.loadReducedStructural(words, "unmeasured", 2, unsupported, "fixture"));
		}
		final ProteinKmerGeometry mixed=ProteinKmerGeometry.make(5, "dimer:1,3,4:"+groups);
		expectFailure(()->DDLPentamerBlacklist.loadReducedStructural(words, "unmeasured", 2, mixed, "fixture"));
		expectFailure(()->DDLPentamerBlacklist.loadReducedStructural(words, "unmeasured", 2,
			ProteinKmerGeometry.exact(6), "fixture"));
		final Path monomer=write(dir, "monomer_reduced.faa", ">a\nAAAAA\n");
		final DDLPentamerBlacklist old=DDLPentamerBlacklist.loadReducedStructural(monomer, "unmeasured", 1,
			ProteinKmerGeometry.make(5, "all:ARNDCEQGHI/LKMFPSTWYV"), "monomer_fixture");
		check(old.contains(0), "Existing full-word monomer structural loading remains supported");
	}

	/** Partitions the 400 ordered amino pairs by ordinal modulo120, preserving literal class order. */
	private static String dimerGroups(){
		final String aa="ARNDCEQGHILKMFPSTWYV";
		final StringBuilder result=new StringBuilder();
		for(int group=0; group<120; group++){
			if(group>0){result.append('/');}
			for(int pair=group; pair<400; pair+=120){
				if(pair>group){result.append(',');}
				result.append(aa.charAt(pair/20)).append(aa.charAt(pair%20));
			}
		}
		return result.toString();
	}

	/** Distinguish K5 prefixes on K6 from literal six-residue filtering and unmeasured structural identity. */
	private static void modes(Path dir, DDLPentamerBlacklist prefix) throws Exception{
		final Path six=write(dir, "six.faa", ">six\nACDEFG\n");
		final DDLPentamerBlacklist literal6=DDLPentamerBlacklist.load(six, pin(six), 6,
			DDLPentamerBlacklist.MatchMode.LITERAL, "shadow_k6");
		check(literal6.size()==1 && literal6.wordLength()==6 && literal6.provenance().equals("shadow_k6") &&
			literal6.contains(pack("ACDEFG", 6)) && !literal6.contains(pack("ACDEF", 5)),
				"Literal K6 loading keeps six-residue keys");
		final Path badDeclaredWords=write(dir, "literal_bad_declared_words.faa",
			">bbgap_blacklist_v1;k=6;words=2;mode=literal;provenance=shadow_k6;aabits=none\nACDEFG\n");
		expectFailure(()->DDLPentamerBlacklist.load(badDeclaredWords, pin(badDeclaredWords), 6,
			DDLPentamerBlacklist.MatchMode.LITERAL, "shadow_k6"));
		final DDLPentamerBlacklist structural=DDLPentamerBlacklist.loadStructural(six, "unmeasured", 1, 6,
			DDLPentamerBlacklist.MatchMode.LITERAL, "real_k6");
		check(structural.sha80().equals("unmeasured") && structural.provenance().equals("real_k6"),
			"Structural loading accepts explicit unmeasured identity with exact count and provenance");
		expectFailure(()->DDLPentamerBlacklist.loadStructural(six, "unmeasured", 2, 6,
			DDLPentamerBlacklist.MatchMode.LITERAL, "real_k6"));
		expectFailure(()->DDLPentamerBlacklist.loadStructural(six, "unmeasured", 1, 6,
			DDLPentamerBlacklist.MatchMode.LITERAL, null));
		expectFailure(()->DDLPentamerBlacklist.load(six, pin(six), 5, DDLPentamerBlacklist.MatchMode.PREFIX, "bad"));
		expectFailure(()->DDLPentamerBlacklist.load(six, pin(six), 6, DDLPentamerBlacklist.MatchMode.PREFIX, "bad"));
		check(DDLPentamerBlacklist.MatchMode.parse("literal")==DDLPentamerBlacklist.MatchMode.LITERAL &&
			DDLPentamerBlacklist.MatchMode.parse("prefix")==DDLPentamerBlacklist.MatchMode.PREFIX,
			"Mode parser accepts the explicit length-free literal and prefix spellings");
		expectFailure(()->DDLPentamerBlacklist.MatchMode.parse("literalK6"));
		expectFailure(()->DDLPentamerBlacklist.MatchMode.parse("prefixK5"));
		expectFailure(()->DDLPentamerBlacklist.loadStructural(six, "unmeasured", 1, 7,
			DDLPentamerBlacklist.MatchMode.LITERAL, "bad"));
		final byte[] sequence="ACDEFGACDEFH".getBytes(StandardCharsets.US_ASCII);
		final char[] actual=new char[1024], expectedPrefix=new char[1024], expectedLiteral=new char[1024];
		addK6Window(sequence, 1, expectedPrefix);
		addK6Window(sequence, 2, expectedPrefix);
		addK6Window(sequence, 3, expectedPrefix);
		addK6Window(sequence, 4, expectedPrefix);
		addK6Window(sequence, 5, expectedPrefix);
		for(int start=1; start+6<=sequence.length; start++){addK6Window(sequence, start, expectedLiteral);}
		DDLFamilyProfiles.sketch(sequence, ProteinKmerGeometry.exact(6), 0, actual, prefix);
		check(Arrays.equals(actual, expectedPrefix), "K5 prefix mode filters both six-mers starting with an excluded pentamer");
		DDLFamilyProfiles.sketch(sequence, ProteinKmerGeometry.exact(6), 0, actual, literal6);
		check(Arrays.equals(actual, expectedLiteral), "Literal K6 mode filters only the exact six-residue word");
		expectFailure(()->DDLFamilyProfiles.sketch(sequence, ProteinKmerGeometry.exact(6), 0, actual,
			DDLPentamerBlacklist.loadStructural(dir.resolve("good.faa"), "unmeasured", 3, 5,
				DDLPentamerBlacklist.MatchMode.LITERAL, "literal5")));
		expectFailure(()->new ddl.FamilySketchBuilder(1024, ProteinKmerGeometry.exact(6), 0,
			false, 1, 0, DDLPentamerBlacklist.loadStructural(dir.resolve("good.faa"), "unmeasured", 3, 5,
				DDLPentamerBlacklist.MatchMode.LITERAL, "literal5")));
	}

	/*--------------------------------------------------------------*/
	/*----------------   Integration and metadata   ----------------*/
	/*--------------------------------------------------------------*/

	/** Actual reference builder, metadata-aware loader, search CLI and frequency CLI share one filter. */
	private static void integration(Path dir, Path blacklistPath, String blacklistPin,
		DDLPentamerBlacklist blacklist) throws Exception{
		assert(blacklist.contains(0)) : "All-A members exercise complete exclusion without changing the family roster";
		final char[] residues=new char[100];
		Arrays.fill(residues, 'A');
		final String sequence=new String(residues);
		final Path families=write(dir, "families.tsv", DDLFamilyHoldout.FAMILY_HEADER+"\n0\tf0\t1\n");
		final Path part=write(dir, "part.tsv", DDLFamilyHoldout.MEMBER_HEADER+"\n0\tf0\tg0|q0\tg0\t0\t"+sequence+"\n");
		final Path chunks=write(dir, "chunks.tsv", DDLFamilyHoldout.CHUNK_HEADER+"\nb00\t0\t"+part+"\t"+pin(part)+"\t1\t100\n");
		final Path reference=dir.resolve("filtered_reference");
		FamilyCalibration.run(new String[]{"mode=references", "families="+families, "familiessha80="+pin(families),
			"chunks="+chunks, "chunkssha80="+pin(chunks), "bucket=0", "buckets=1", "expectedfamilies=1", "expectedbatches=1",
			"expectedmembers=1", "maxbucketmembers=10", "maxbucketbytes=1048576", "maxoutputbytes=1048576",
			"positional=t", "partitions=4", "overlap=8", "blacklist="+blacklistPath, "blacklistsha80="+blacklistPin,
			"out="+reference});
		final Path profile=reference.resolve("profiles.tsv"), summary=reference.resolve("summary.tsv");
		final Path list=write(dir, "filtered.list.tsv", profile+"\t"+pin(profile)+"\t"+pin(summary)+"\n");
		final DDLFamilySearch.Index index=DDLFamilySearch.loadForAttribution(list, 1024, 5, 0,
			DDLFamilyProfiles.SKETCH_POSITIONAL, 4, 8, blacklist);
		check(index.families()==1 && Arrays.equals(index.profileCopy(0), new char[1024]),
			"Regenerated modal family preserves roster but excludes every AAAAA word");
		expectFailure(()->DDLFamilySearch.loadForAttribution(list, 1024, 5, 0, DDLFamilyProfiles.SKETCH_POSITIONAL, 4, 8));
		final Path query=write(dir, "query.faa", ">q0\n"+sequence+"\n"), search=dir.resolve("filtered_search");
		DDLFamilySearch.main(new String[]{"profileslist="+list, "listsha80="+pin(list), "in="+query, "insha80="+pin(query),
			"buckets=1024", "mode=nonempty", "top=1", "positional=t", "partitions=4", "overlap=8",
			"blacklist="+blacklistPath, "blacklistsha80="+blacklistPin, "out="+search});
		final java.util.List<String> searchSummary=Files.readAllLines(search.resolve("summary.tsv"), StandardCharsets.US_ASCII);
		final String[] names=searchSummary.get(0).split("\t"), values=searchSummary.get(1).split("\t");
		check(names.length==values.length && values[Arrays.asList(names).indexOf("no_hits")].equals("1"),
			"Matching query filter leaves the all-excluded query as a no-hit");
		check(values[Arrays.asList(names).indexOf("blacklist_k")].equals("5") &&
			values[Arrays.asList(names).indexOf("blacklist_mode")].equals("prefix") &&
			values[Arrays.asList(names).indexOf("blacklist_provenance")].equals("unrecorded"),
			"Search output records explicit blacklist length, mode and provenance");
		final Path counts=dir.resolve("filtered_counts");
		DDLKeyFrequency.main(new String[]{"profileslist="+list, "listsha80="+pin(list), "expected=1", "sketch=positional-v1",
			"partitions=4", "overlap=8", "blacklist="+blacklistPath, "blacklistsha80="+blacklistPin, "out="+counts});
		final java.util.List<String> buckets=Files.readAllLines(counts.resolve("buckets.tsv"), StandardCharsets.US_ASCII);
		check(buckets.size()==1025 && buckets.get(1).equals("0\t1\t1\t0\t0\t0"),
			"Family collision counter retains empty slots without false zero-key matches");
	}

	/** Refuse both directions of filtered/unfiltered mismatch and incomplete identity metadata. */
	private static void metadata(DDLPentamerBlacklist blacklist) throws Exception{
		final HashMap<String, String> fields=new HashMap<String, String>();
		DDLFamilySearch.validateBlacklistMetadata(fields, null);
		expectFailure(()->DDLFamilySearch.validateBlacklistMetadata(fields, blacklist));
		fields.put("blacklist_sha80", blacklist.sha80());
		expectFailure(()->DDLFamilySearch.validateBlacklistMetadata(fields, blacklist));
		fields.put("blacklist_words", Integer.toString(blacklist.size()));
		DDLFamilySearch.validateBlacklistMetadata(fields, blacklist);
		fields.put("blacklist_k", Integer.toString(blacklist.wordLength()));
		fields.put("blacklist_mode", blacklist.matchMode().label());
		fields.put("blacklist_provenance", blacklist.provenance());
		DDLFamilySearch.validateBlacklistMetadata(fields, blacklist);
		fields.put("blacklist_mode", "literal");
		expectFailure(()->DDLFamilySearch.validateBlacklistMetadata(fields, blacklist));
		fields.put("blacklist_mode", blacklist.matchMode().label());
		fields.remove("blacklist_provenance");
		expectFailure(()->DDLFamilySearch.validateBlacklistMetadata(fields, blacklist));
		fields.put("blacklist_provenance", blacklist.provenance());
		final DDLPentamerBlacklist literal5=DDLPentamerBlacklist.loadStructural(write(Files.createTempDirectory("ddl-pentamer-blacklist-meta-"),
			"literal5.faa", ">a\nAAAAA\n"), "unmeasured", 1, 5, DDLPentamerBlacklist.MatchMode.LITERAL, "literal5");
		final HashMap<String, String> legacyFields=new HashMap<String, String>();
		legacyFields.put("blacklist_sha80", literal5.sha80());
		legacyFields.put("blacklist_words", Integer.toString(literal5.size()));
		expectFailure(()->DDLFamilySearch.validateBlacklistMetadata(legacyFields, literal5));
		expectFailure(()->DDLFamilySearch.validateBlacklistMetadata(fields, null));
		fields.put("blacklist_words", Integer.toString(blacklist.size()+1));
		expectFailure(()->DDLFamilySearch.validateBlacklistMetadata(fields, blacklist));
		fields.put("blacklist_words", Integer.toString(blacklist.size()));
		fields.put("blacklist_sha80", "00000000000000000000");
		expectFailure(()->DDLFamilySearch.validateBlacklistMetadata(fields, blacklist));
	}

	/** Check null parity, raw-zero exclusion under all tags, and independent per-window maxima. */
	private static void sketches(DDLPentamerBlacklist blacklist){
		final byte[] allA=new byte[200];
		Arrays.fill(allA, (byte)'A');
		final char[] old=new char[1024], actual=new char[1024], expected=new char[1024], one=new char[1024];
		for(int partitions : new int[]{2, 4, 8, 16}){
			DDLFamilyProfiles.sketchPositional(allA, 5, 0, partitions, 8, old);
			check(DDLFamilyProfiles.sketchPositional(allA, 5, 0, partitions, 8, actual, blacklist)==196,
				"Excluded words still count toward input-window conservation");
			check(Arrays.equals(actual, expected) && !Arrays.equals(old, expected),
				"Raw-zero exclusion applies to every positional and overlap tag");
			DDLFamilyProfiles.sketchPositional(allA, 5, 0, partitions, 8, actual, null);
			check(Arrays.equals(actual, old), "Null blacklist preserves positional bytes");
		}
		final byte[] sequence="AAAAAACDEFGHIKLMNPQRSTVWYAAAAARACDEF".getBytes(StandardCharsets.US_ASCII);
		for(int start=0; start+5<=sequence.length; start++){
			final byte[] word=Arrays.copyOfRange(sequence, start, start+5);
			if(blacklist.contains(DDLPentamerBlacklist.pack(word))){continue;}
			DDLFamilyProfiles.sketch(word, 5, 0, one);
			for(int b=0; b<one.length; b++){expected[b]=(char)Math.max(expected[b], one[b]);}
		}
		check(DDLFamilyProfiles.sketch(sequence, 5, 0, actual, blacklist)==sequence.length-4,
			"Plain filter preserves valid input-window count");
		check(Arrays.equals(actual, expected), "Unexcluded windows compete for register maxima normally");
		DDLFamilyProfiles.sketch(sequence, 5, 0, old);
		DDLFamilyProfiles.sketch(sequence, 5, 0, actual, null);
		check(Arrays.equals(actual, old), "Null blacklist preserves plain bytes");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Fixture helpers        ----------------*/
	/*--------------------------------------------------------------*/

	/** Encode test literals through the API under test; expected codes are specified independently. */
	private static long pack(String word){return DDLPentamerBlacklist.pack(word.getBytes(StandardCharsets.US_ASCII));}

	/** Encodes a literal of the explicitly requested length through the API under test. */
	private static long pack(String word, int length){
		return DDLPentamerBlacklist.pack(word.getBytes(StandardCharsets.US_ASCII), length);
	}

	/** Sketches one six-residue window and merges register maxima into caller-owned expected scratch. */
	private static void addK6Window(byte[] sequence, int start, char[] expected){
		final byte[] word=Arrays.copyOfRange(sequence, start, start+6);
		final char[] one=new char[expected.length];
		DDLFamilyProfiles.sketch(word, ProteinKmerGeometry.exact(6), 0, one, null);
		for(int b=0; b<one.length; b++){expected[b]=(char)Math.max(expected[b], one[b]);}
	}

	/** Write a tiny native-parser fixture, never a scientific corpus. */
	private static Path write(Path dir, String name, String text) throws Exception{
		final Path path=dir.resolve(name);
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII));
		return path;
	}

	/** Bind literal fixture bytes using the same final80 convention as recorded inputs. */
	private static String pin(Path path) throws Exception{
		final byte[] bytes=MessageDigest.getInstance("SHA-256").digest(Files.readAllBytes(path));
		final StringBuilder result=new StringBuilder(20);
		for(int i=bytes.length-10; i<bytes.length; i++){
			result.append(Character.forDigit((bytes[i]>>>4)&15, 16)).append(Character.forDigit(bytes[i]&15, 16));
		}
		return result.toString();
	}

	/** Reject bad inputs explicitly; unrelated I/O failures are not a passing negative fixture. */
	private static void expectFailure(Action action) throws Exception{
		boolean rejected=false;
		try{action.run();}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Invalid input must fail explicitly");
	}

	/** Count independent checks and fail even if the caller disabled assertions. */
	private static void check(boolean condition, String message){
		checks++;
		if(!condition){throw new AssertionError(message);}
	}

	/** One checked fixture operation; unexpected I/O failures propagate to the caller. */
	private interface Action{void run() throws Exception;}
	/** Attempted checks in this JVM; not reset by main or its submodes. */
	private static long checks;
}
