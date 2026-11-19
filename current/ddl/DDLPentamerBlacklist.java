package ddl;

import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.attribute.BasicFileAttributes;
import java.security.MessageDigest;
import java.util.Objects;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.LongHashSet;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;

/** Immutable exclusion set of forward raw or explicitly reduced amino-acid words.
 * Test membership before adding positional tags or hashing; never use nucleotide canonicalization.
 * @author Fischl
 */
public final class DDLPentamerBlacklist{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Keep the mutable primitive set private after loading so workers can share it read-only. */
	private DDLPentamerBlacklist(LongHashSet words_, String sha80_, int wordLength_,
			MatchMode mode_, String provenance_, ProteinKmerGeometry projection_,boolean measured_){
		assert(words_!=null && !words_.isEmpty()) : "A configured blacklist must exclude at least one word";
		words=words_; sha80=sha80_; wordLength=wordLength_; mode=mode_; provenance=provenance_;
		projection=projection_;measured=measured_;
	}

	/** Read a pinned protein FASTA containing one distinct five-residue word per record.
	 * Use a protein extension such as .faa, or enable the caller's normal amino input mode.
	 * File bytes are checked both before and after parsing; no global sketch state changes.
	 */
	public static DDLPentamerBlacklist load(Path path, String pin) throws Exception{
		return load(path, pin, 5, MatchMode.PREFIX, null);
	}

	/** Read a pinned protein FASTA using an explicit raw-word selection policy. */
	public static DDLPentamerBlacklist load(Path path, String pin, int wordLength, MatchMode mode,
			String provenance) throws Exception{
		return load(path, pin, -1, wordLength, mode, provenance, true, null);
	}

	/** Read structurally while retaining accepted upstream identity; no digest is computed. */
	public static DDLPentamerBlacklist loadStructural(Path path, String inheritedSha80, int expectedWords) throws Exception{
		return loadStructural(path, inheritedSha80, expectedWords, 5, MatchMode.PREFIX, null);
	}

	/** Read structurally while retaining explicit identity, count, word length, mode and provenance. */
	public static DDLPentamerBlacklist loadStructural(Path path, String inheritedIdentity, int expectedWords,
			int wordLength, MatchMode mode, String provenance) throws Exception{
		if(expectedWords<1){throw new IllegalArgumentException("Structural blacklist loading needs an exact positive word count");}
		return load(path, inheritedIdentity, expectedWords, wordLength, mode, provenance, false, null);
	}

	/** Load representative protein words in an explicit full-word reduced alphabet, without hashing.
	 * Any canonical member of a class can represent it; records collapsing to the same key are duplicates.
	 * Nonoverlapping K6 dimers use three independent ordered-pair classes.
	 * The immutable projection must also be used by both member and query sketchers.
	 */
	public static DDLPentamerBlacklist loadReducedStructural(Path path, String inheritedIdentity, long expectedWords,
			ProteinKmerGeometry projection, String provenance) throws Exception{
		if(expectedWords<1 || expectedWords>Integer.MAX_VALUE){
			throw new IllegalArgumentException("Reduced blacklist loading needs an exact positive word count within integer capacity");
		}
		if(!reducedGeometry(projection)){
			throw new IllegalArgumentException("Reduced blacklist loading requires a full-word monomer alphabet or nonoverlapping K6 dimers at1,3,5");
		}
		return load(path, inheritedIdentity, (int)expectedWords, projection.k(), MatchMode.REDUCED_LITERAL,
			provenance, false, projection);
	}

	/** Read the pinned reduced resource with explicit projection and exact cardinality. */
	public static DDLPentamerBlacklist loadReduced(Path path,String identity,int expectedWords,
			ProteinKmerGeometry projection,String provenance)throws Exception{
		if(expectedWords<1 || !reducedGeometry(projection)){
			throw new IllegalArgumentException("Pinned reduced loading requires a positive count and full-word reduced geometry or nonoverlapping K6 dimers at1,3,5");
		}
		return load(path,identity,expectedWords,projection.k(),MatchMode.REDUCED_LITERAL,provenance,true,projection);
	}

	/** Read published bytes while retaining the pre-header policy identity. */
	static DDLPentamerBlacklist loadPublished(Path path, String filePin, String inheritedIdentity,
			int expectedWords, int wordLength, MatchMode mode, String provenance,
			ProteinKmerGeometry projection) throws Exception{
		checkPinFormat(path, filePin, false);
		checkPinFormat(path, inheritedIdentity, false);
		return load(path, inheritedIdentity, expectedWords, wordLength, mode, provenance, filePin, projection);
	}

	/** Measured and structural resources use the same full-word projection contract. */
	private static boolean reducedGeometry(ProteinKmerGeometry projection){
		return projection!=null && (projection.isReduced() ||
			projection.k()==6 && projection.name().startsWith("dimer6+1,3,5+"));
	}

	/** Cache equivalence to one exact immutable official object, never just a claimed pin. */
	void requireOfficial(DDLPentamerBlacklist official,ProteinKmerGeometry geometry,String policyLabel){
		validate(geometry);
		if(officialValidatedAgainst==official){return;}
		if(mode!=official.mode || wordLength!=official.wordLength || words.size()!=official.words.size() ||
			!sha80.equals(official.sha80)){
			throw new IllegalArgumentException(policyLabel+" blacklist is required, not an omitted or different policy");
		}
		for(long word:words.toArray()){
			if(!official.words.contains(word)){throw new IllegalArgumentException(policyLabel+" blacklist keys differ from the official resource");}
		}
		officialValidatedAgainst=official;
	}

	private static DDLPentamerBlacklist load(Path path, String pin, int expectedWords, int wordLength,
			MatchMode mode, String provenance, boolean verifyPin, ProteinKmerGeometry projection) throws Exception{
		return load(path, pin, expectedWords, wordLength, mode, provenance, verifyPin ? pin : null, projection);
	}

	private static DDLPentamerBlacklist load(Path path, String pin, int expectedWords, int wordLength,
			MatchMode mode, String provenance, String filePin, ProteinKmerGeometry projection) throws Exception{
		checkConfiguration(wordLength, mode);
		if((mode==MatchMode.REDUCED_LITERAL)!=(projection!=null)){
			throw new IllegalArgumentException("Reduced literal blacklists require an explicit projection through loadReducedStructural");
		}
		if(filePin==null){checkPinFormat(path, pin, true);}
		else{checkPin(path, filePin); checkPinFormat(path, pin, true);}
		final BasicFileAttributes before=filePin==null ? regularFileAttributes(path) : null;
		if(filePin==null && pin.equals("unmeasured") && (provenance==null || provenance.length()==0)){
			throw new IllegalArgumentException("Unmeasured structural blacklists require explicit provenance");
		}
		final FileFormat format=FileFormat.testInput(path.toString(), FileFormat.FASTA, null, false, true);
		if(!format.fasta() || !(format.amino() || Shared.AMINO_IN)){
			throw new IllegalArgumentException("Protein blacklist requires protein FASTA (.faa or amino input mode)");
		}
		final LongHashSet words=new LongHashSet(8192);
		final Streamer input=StreamerFactory.makeStreamer(format, null, true, -1, false, true, 1);
		boolean checkedHeader=false;
		int declaredWords=-1;
		Throwable failure=null;
		input.start();
		try{
			for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
				for(Read read : batch.list){
					if(read==null || read.mate!=null){throw new IllegalArgumentException("Blacklist records must be unpaired protein words");}
					if(!checkedHeader){
						declaredWords=validateHeader(read.id, expectedWords, wordLength, mode, normalizeProvenance(provenance), projection);
						checkedHeader=true;
					}
					final long raw=pack(read.bases, wordLength), word=projection==null ? raw : projection.word(raw);
					if(!words.add(word)){throw new IllegalArgumentException("Duplicate blacklist word: "+read.id);}
				}
			}
		}catch(RuntimeException|Error e){
			failure=e;
			throw e;
		}finally{finishInput(input, failure);}
		if(words.isEmpty()){throw new IllegalArgumentException("Empty protein blacklist: "+path);}
		if(expectedWords>=0 && words.size()!=expectedWords){
			throw new IllegalArgumentException("Protein blacklist word count differs from expected count: "+words.size()+" != "+expectedWords);
		}
		if(declaredWords>=0 && words.size()!=declaredWords){
			throw new IllegalArgumentException("Protein blacklist word count differs from versioned header: "+words.size()+" != "+declaredWords);
		}
		if(filePin!=null){checkPin(path, filePin);}
		else{checkSameFile(path, before);}
		return new DDLPentamerBlacklist(words, pin, wordLength, mode, normalizeProvenance(provenance), projection,filePin!=null && filePin.equals(pin));
	}

	/** Consume the selected FASTA reader's terminal before closing, even after a rejected word.
	 * FastaStreamerST/2ST.close do not cancel their bounded-queue producers. All four
	 * FASTA implementations selected here keep hasMore true until terminal consumption;
	 * ST/2ST can also return null on interruption without consuming that terminal.
	 * Keep the original validation failure primary if draining or closure also fails.
	 */
	private static void finishInput(Streamer input, Throwable failure){
		assert(input!=null) : "Blacklist cleanup owns the FASTA reader started in load";
		Throwable primary=failure;
		boolean interrupted=Thread.interrupted();
		try{
			while(input.hasMore()){
				if(input.nextList()==null && input.hasMore()){interrupted=true;}
			}
		}catch(RuntimeException|Error e){
			if(primary==null){primary=e; throw e;}
			if(primary!=e){primary.addSuppressed(e);}
		}finally{
			try{
				if(ReadWrite.closeStream(input)){
					throw new IllegalStateException("Blacklist FASTA read failed");
				}
			}catch(RuntimeException|Error e){
				if(primary==null){throw e;}
				if(primary!=e){primary.addSuppressed(e);}
			}finally{if(interrupted){Thread.currentThread().interrupt();}}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Membership in this list's raw or reduced keyspace; use rawWord to project a rolling word first.
	 * A positional tag or hashed key must never enter this lookup.
	 */
	public boolean contains(long rawWord){
		assert(rawWord>=0 && rawWord<(1L<<(5*wordLength))) : "Blacklist lookup precedes positional tagging and hashing";
		return words.contains(rawWord);
	}

	/** Validate this list against a sketch geometry before a hot sketching loop starts. */
	public void validate(ProteinKmerGeometry geometry){
		if(geometry==null){throw new IllegalArgumentException("A blacklist requires an explicit protein geometry");}
		if(projection!=null && !projection.equals(geometry)){
			throw new IllegalArgumentException("Reduced blacklist alphabet and geometry must match the sketch geometry: "+
				projection.name()+" != "+geometry.name());
		}
		if(mode==MatchMode.LITERAL && geometry.span()!=wordLength){
			throw new IllegalArgumentException("Literal protein blacklist length "+wordLength+
				" does not match raw geometry span "+geometry.span());
		}
		if(mode==MatchMode.PREFIX && geometry.span()<wordLength){
			throw new IllegalArgumentException("Prefix protein blacklist length "+wordLength+
				" exceeds raw geometry span "+geometry.span());
		}
	}

	/** Derive the lookup coordinate from a raw rolling word, before positional tags or hashing.
	 * Raw modes retain their historical prefix behavior; reduced mode projects the entire word.
	 * Call validate once before the sketch loop to bind the supplied geometry.
	 */
	public long rawWord(ProteinKmerGeometry geometry, long rolling){
		return projection==null ? geometry.rawPrefix(rolling, wordLength) : projection.word(rolling);
	}

	/** Number of unique excluded words, not occurrences or DDL keys. */
	public int size(){return words.size();}

	/** Stable input identity for binding query and reference recipes. */
	public String sha80(){return sha80;}
	/** Whether this load measured the declared digest; inherited identities remain structural. */
	public boolean identityMeasured(){return measured;}

	public int wordLength(){return wordLength;}

	public MatchMode matchMode(){return mode;}

	public String provenance(){return provenance;}

	/** Extended metadata for new filtered profiles; old summaries carry only sha80 and size. */
	public String metadataColumns(){
		return "blacklist_k\tblacklist_mode\tblacklist_provenance";
	}

	public String metadataValues(){
		return wordLength+"\t"+mode.label()+"\t"+provenance;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Match DDLFamilyProfiles' forward five-bit encoding; ambiguity and stops are forbidden. */
	public static long pack(byte[] sequence){
		return pack(sequence, 5);
	}

	/** Match DDLFamilyProfiles' forward five-bit encoding; ambiguity and stops are forbidden. */
	public static long pack(byte[] sequence, int wordLength){
		checkWordLength(wordLength);
		if(sequence==null || sequence.length!=wordLength){
			throw new IllegalArgumentException("Each blacklist record must contain exactly "+wordLength+" amino acids");
		}
		long word=0;
		for(byte residue : sequence){
			final int aa=residue>=0 && residue<AminoAcid.acidToNumberNoStops.length ? AminoAcid.acidToNumberNoStops[residue] : -1;
			if(aa<0 || aa>=20){throw new IllegalArgumentException("Protein blacklist words require standard amino acids");}
			word=(word<<5)|aa;
		}
		assert(word>=0 && word<(1L<<(5*wordLength))) : "Five-bit AA codes must fit the raw blacklist keyspace";
		return word;
	}

	/** Compare the final80 digest bits without exposing a full digest or trusting a filename. */
	private static void checkPin(Path path, String expected) throws Exception{
		checkPinFormat(path, expected, false);
		final MessageDigest digest=MessageDigest.getInstance("SHA-256");
		final byte[] buffer=new byte[65536];
		try(InputStream input=Files.newInputStream(path)){
			for(int n=input.read(buffer); n>=0; n=input.read(buffer)){if(n>0){digest.update(buffer, 0, n);}}
		}
		final byte[] bytes=digest.digest();
		final StringBuilder actual=new StringBuilder(20);
		for(int i=bytes.length-10; i<bytes.length; i++){
			actual.append(Character.forDigit((bytes[i]>>>4)&15, 16)).append(Character.forDigit(bytes[i]&15, 16));
		}
		if(!expected.equals(actual.toString())){throw new IllegalArgumentException("Blacklist bytes differ from pinned input: "+path);}
	}

	/** Structural mode still needs a path and inherited identity. */
	private static void checkPinFormat(Path path, String expected, boolean allowUnmeasured){
		if(path==null || expected==null || !(expected.matches("[0-9a-f]{20}") || allowUnmeasured && expected.equals("unmeasured"))){
			throw new IllegalArgumentException("Require blacklist path and final80 SHA256 fingerprint or explicit unmeasured identity");
		}
	}

	private static BasicFileAttributes regularFileAttributes(Path path) throws Exception{
		final BasicFileAttributes attrs=Files.readAttributes(path, BasicFileAttributes.class, LinkOption.NOFOLLOW_LINKS);
		if(!attrs.isRegularFile()){
			throw new IllegalArgumentException("Structural blacklist path must be a regular file, not a symlink or directory: "+path);
		}
		return attrs;
	}

	private static void checkSameFile(Path path, BasicFileAttributes before) throws Exception{
		final BasicFileAttributes after=regularFileAttributes(path);
		if(before.size()!=after.size() || !before.lastModifiedTime().equals(after.lastModifiedTime()) ||
				!Objects.equals(before.fileKey(), after.fileKey())){
			throw new IllegalArgumentException("Structural blacklist changed while being read: "+path);
		}
	}

	private static void checkConfiguration(int wordLength, MatchMode mode){
		checkWordLength(wordLength);
		if(mode==null){throw new IllegalArgumentException("A blacklist requires a matching mode");}
		if(mode==MatchMode.PREFIX && wordLength!=5){
			throw new IllegalArgumentException("Prefix blacklists are supported only for raw K5");
		}
	}

	private static void checkWordLength(int wordLength){
		if(wordLength<5 || wordLength>6){throw new IllegalArgumentException("Protein blacklists support only K5 or K6 words");}
	}

	private static String normalizeProvenance(String provenance){
		if(provenance==null || provenance.length()==0){return "unrecorded";}
		if(provenance.indexOf('\t')>=0 || provenance.indexOf('\n')>=0 || provenance.indexOf('\r')>=0){
			throw new IllegalArgumentException("Blacklist provenance must fit one TSV field");
		}
		return provenance;
	}

	/** Validate an optional versioned header embedded in the first FASTA record ID. */
	private static int validateHeader(String id, int expectedWords, int wordLength, MatchMode mode,
			String provenance, ProteinKmerGeometry projection){
		if(id==null || !id.startsWith("bbgap_blacklist_v")){return -1;}
		final String[] fields=id.split(";",-1);
		if(fields.length<2 || !fields[0].equals("bbgap_blacklist_v1")){
			throw new IllegalArgumentException("Unknown BBGAP blacklist header version: "+id);
		}
		int k=-1,words=-1;MatchMode declaredMode=null;String declaredProvenance=null,aabits=null;
		for(int i=1; i<fields.length; i++){
			final String field=fields[i];
			final int split=field.indexOf('=');
			if(split<1 || split==field.length()-1){throw new IllegalArgumentException("Malformed BBGAP blacklist header field: "+field);}
			final String key=field.substring(0, split), value=decodeHeaderValue(field.substring(split+1));
			if(key.equals("k")){k=Integer.parseInt(value);}
			else if(key.equals("words")){words=Integer.parseInt(value);}
			else if(key.equals("mode")){declaredMode=MatchMode.parse(value);}
			else if(key.equals("provenance")){declaredProvenance=value;}
			else if(key.equals("aabits")){aabits=value;}
			else{throw new IllegalArgumentException("Unknown BBGAP blacklist header field: "+key);}
		}
		if(k!=wordLength || words<1 || declaredMode!=mode || aabits==null ||
				expectedWords>=0 && words!=expectedWords || !provenance.equals(declaredProvenance)){
			throw new IllegalArgumentException("BBGAP blacklist header does not match the requested policy");
		}
		final ProteinKmerGeometry headerProjection=aabits==null || aabits.equals("none") ? null : ProteinKmerGeometry.make(k, aabits);
		if(projection==null ? headerProjection!=null : !projection.equals(headerProjection)){
			throw new IllegalArgumentException("BBGAP blacklist header alphabet does not match the requested projection");
		}
		return words;
	}

	/** Decode a semicolon-delimited FASTA header value without allowing raw whitespace. */
	private static String decodeHeaderValue(String value){
		final StringBuilder decoded=new StringBuilder(value.length());
		for(int i=0; i<value.length(); i++){
			final char c=value.charAt(i);
			if(c=='%'){
				if(i+2>=value.length()){throw new IllegalArgumentException("Truncated BBGAP blacklist header escape: "+value);}
				final int hi=Character.digit(value.charAt(++i), 16), lo=Character.digit(value.charAt(++i), 16);
				if(hi<0 || lo<0){throw new IllegalArgumentException("Malformed BBGAP blacklist header escape: "+value);}
				final char decodedChar=(char)((hi<<4)|lo);
				if(Character.isWhitespace(decodedChar)){
					throw new IllegalArgumentException("Escaped whitespace is not allowed in BBGAP blacklist header fields");
				}
				decoded.append(decodedChar);
			}else{
				if(Character.isWhitespace(c)){
					throw new IllegalArgumentException("Raw whitespace is not allowed in BBGAP blacklist header fields");
				}
				decoded.append(c);
			}
		}
		return decoded.toString();
	}

	public enum MatchMode{
		LITERAL("literal"),
		PREFIX("prefix"),
		REDUCED_LITERAL("reduced_literal");

		MatchMode(String label_){label=label_;}

		public String label(){return label;}

		public static MatchMode parse(String text){
			if(text==null){return PREFIX;}
			final String clean=text.trim();
			if(clean.equalsIgnoreCase("literal")){
				return LITERAL;
			}
			if(clean.equalsIgnoreCase("prefix")){
				return PREFIX;
			}
			if(clean.equalsIgnoreCase("reduced_literal")){return REDUCED_LITERAL;}
			throw new IllegalArgumentException("Unknown blacklist mode: "+text);
		}

		private final String label;
	}

	private final LongHashSet words;
	private final String sha80;
	private final boolean measured;
	private final int wordLength;
	private final MatchMode mode;
	private final String provenance;
	/** Null for raw policies; otherwise the exact full-word alphabet used for stored keys. */
	private final ProteinKmerGeometry projection;
	/** Memoized validation only; the set and its scientific meaning never change. */
	private volatile DDLPentamerBlacklist officialValidatedAgainst;
}
