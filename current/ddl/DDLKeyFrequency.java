package ddl;

import java.io.InputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;
import template.ThreadWaiter;

/** Exact background frequency of unsigned16 keys in family or raw-member sketches.
 * Matching-register profile pairs are a collision proxy, not false positives.
 * Family mode loads a declared plain or positional recipe; member mode is plain,
 * unfiltered K5/seed0. This is a descriptive census, not an accuracy measurement.
 * A chunk limit selects a manifest prefix and reports it as partial; it is not
 * random sampling. Outputs are created only after the selected record count matches.
 * Counters are mutable and not thread-safe: chunk workers own separate counters,
 * which are merged after joining. Member counters also retain a dense 25-bit word table.
 * The CLI changes shared amino-input, thread and compression settings without restoring them;
 * use a separate process. Legacy file pins remain part of this tool's input contract.
 * @author Fischl
 */
public final class DDLKeyFrequency{

	/** Run the descriptive analysis; never alter reference profiles or calibration. */
	public static void main(String[] args) throws Exception{
		final DDLKeyFrequency tool=new DDLKeyFrequency(args);
		try{
			tool.process();
		}finally{
			Shared.closeStream(tool.outstream);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses help/config preprocessing output and validates the mutually exclusive input modes.
	 * Family mode requires profileslist/listsha80; member mode requires explicit expected
	 * and either in/insha80 or chunks/chunkssha80. Parallel counting requires chunks.
	 * This constructor validates settings; process performs I/O and count verification. */
	public DDLKeyFrequency(String[] args){
		assert(args!=null) : "CLI arguments bind the input recipe for reproducible counts";
		final PreParser pp=new PreParser(args, getClass(), false);
		outstream=pp.outstream;
		for(String arg : pp.args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String a=pair[0].toLowerCase(Locale.ROOT), b=pair[1];
			if(a.equals("mode")){
				mode=b.toLowerCase(Locale.ROOT);
			}else if(a.equals("profileslist")){
				profilesList=Paths.get(b).toAbsolutePath();
			}else if(a.equals("listsha80")){
				listSha80=b;
			}else if(a.equals("sketch")){
				sketchMode=b.toLowerCase(Locale.ROOT);
			}else if(a.equals("partitions")){
				partitions=Parse.parseIntKMG(b);
			}else if(a.equals("overlap")){
				overlap=Parse.parseIntKMG(b);
			}else if(a.equals("blacklist")){
				blacklistPath=Paths.get(b).toAbsolutePath();
			}else if(a.equals("blacklistsha80")){
				blacklistPin=b;
			}else if(a.equals("in")){
				input=Paths.get(b).toAbsolutePath();
			}else if(a.equals("insha80")){
				inputSha80=b;
			}else if(a.equals("chunks")){
				chunks=Paths.get(b).toAbsolutePath();
			}else if(a.equals("chunkssha80")){
				chunksSha80=b;
			}else if(a.equals("chunklimit")){
				chunkLimit=Parse.parseIntKMG(b);
			}else if(a.equals("out")){
				output=Paths.get(b).toAbsolutePath();
			}else if(a.equals("buckets")){
				width=Parse.parseIntKMG(b);
			}else if(a.equals("top")){
				top=Parse.parseIntKMG(b);
			}else if(a.equals("threads") || a.equals("t")){
				threads=Parse.parseIntKMG(b);
			}else if(a.equals("expected")){
				expected=Parse.parseIntKMG(b);
				expectedSet=true;
			}else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		if(output==null || (!mode.equals("families") && !mode.equals("members"))){
			throw new IllegalArgumentException("Require mode=families/members and a new out directory");
		}
		if(mode.equals("families")){
			if(profilesList==null || !validPin(listSha80) || input!=null || inputSha80!=null || chunks!=null ||
				chunksSha80!=null || chunkLimit!=Integer.MAX_VALUE){
				throw new IllegalArgumentException("Family mode requires only profileslist/listsha80 inputs");
			}
		}else{
			if(!expectedSet || profilesList!=null || listSha80!=null || (input==null)==(chunks==null)){
				throw new IllegalArgumentException("Member mode requires explicit expected and either in or chunks, not family profiles");
			}
			if(input!=null ? !validPin(inputSha80) || chunksSha80!=null ||
				chunkLimit!=Integer.MAX_VALUE : !validPin(chunksSha80) || inputSha80!=null){
				throw new IllegalArgumentException("Require matching input sha80; chunklimit applies only to a chunk manifest");
			}
		}
		if(width<128 || width>16384 || Integer.bitCount(width)!=1 || top<1 || top>65535 || expected<1 || chunkLimit<1 ||
			threads<1){
			throw new IllegalArgumentException("Require power-of-two buckets=128..16384, top=1..65535, expected>0, threads>0");
		}
		if(threads>1 && chunks==null){
			throw new IllegalArgumentException("Parallel counting requires the independently readable chunks input");
		}
		if(!sketchMode.equals(DDLFamilyProfiles.SKETCH_PLAIN) && !sketchMode.equals(DDLFamilyProfiles.SKETCH_POSITIONAL)){
			throw new IllegalArgumentException("Unknown family sketch mode");
		}
		if((blacklistPath==null)!=(blacklistPin==null) || (blacklistPin!=null && !validPin(blacklistPin))){
			throw new IllegalArgumentException("Blacklist path and fingerprint must be supplied together");
		}
		if(mode.equals("members") && (!sketchMode.equals(DDLFamilyProfiles.SKETCH_PLAIN) || blacklistPath!=null)){
			throw new IllegalArgumentException("Member census remains plain/unfiltered; regenerate family profiles separately");
		}
		DDLFamilyProfiles.validatePositionalLayout(partitions, overlap);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Counts one observation per occupied register and publishes tables after population validation.
	 * Requires a nonexistent output directory with an existing parent. Output failures can
	 * leave partial files; PASS is written only after all tables and metadata succeed.
	 * This instance is single-use: chunk-selection counters are not reset here. */
	private void process() throws Exception{
		assert(expected>0) : "Expected population size catches missing inputs before publishing counts";
		if(Files.exists(output)){throw new IllegalArgumentException("Output directory already exists: "+output);}
		// Parallelism belongs to the chunk pool, not another decoder pool per open file.
		Shared.setThreads(1);
		Shared.AMINO_IN=true;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		final long started=System.nanoTime();
		final Counts counts=new Counts(width);
		final MemberCounts members=mode.equals("members") ? new MemberCounts(counts, width) : null;
		if(members==null){
			readFamilies(counts);
		}else if(chunks!=null){
			readChunks(members);
		}else{
			readFasta(members);
		}
		if(counts.records!=selectedExpected){
			throw new IllegalStateException("Expected "+selectedExpected+" selected profiles, counted "+counts.records);
		}
		Files.createDirectory(output);
		writeTables(counts);
		if(members!=null){writeWords(members);}
		final String recipe="mode="+mode+"\nprofileslist="+profilesList+"\nlistsha80="+listSha80+"\nexpected="+expected+
			"\nin="+input+"\ninsha80="+inputSha80+"\nchunks="+chunks+"\nchunkssha80="+chunksSha80+"\nchunklimit="+chunkLimit+
			"\nbuckets="+width+"\ntop="+top+"\nthreads="+threads+"\nk=5\nseed=0\nfamilymode=nonempty\nsketch="+sketchMode+"\n"+
			(sketchMode.equals(DDLFamilyProfiles.SKETCH_POSITIONAL) ? "partitions="+partitions+"\noverlap="+overlap+"\n" : "")+
			(blacklistPath==null ? "" : "blacklist="+blacklistPath+"\nblacklistsha80="+blacklistPin+"\n");
		Files.write(output.resolve("recipe.txt"), recipe.getBytes(StandardCharsets.UTF_8));
		final String summary="population\trecords\tresidues\tvalid_words\tselected_chunks\tmanifest_chunks\tpartial\n"+
			mode+"\t"+counts.records+"\t"+(members==null ? "." : members.residues)+"\t"+(members==null ? "." : members.validWords)+
			"\t"+selectedChunks+"\t"+manifestChunks+"\t"+partial+"\n";
		Files.write(output.resolve("summary.tsv"), summary.getBytes(StandardCharsets.US_ASCII));
		Files.write(output.resolve("PASS"),
			(partial ? "DDL_KEY_FREQUENCY_PILOT_PASS\n" : "DDL_KEY_FREQUENCY_PASS\n").getBytes(StandardCharsets.US_ASCII));
		outstream.println("Counted "+counts.records+" "+mode+" profiles in "+((System.nanoTime()-started)*1e-9)+" seconds.");
	}

	/** Load matched modal profiles without inverted postings; validate the complete sketch recipe. */
	private void readFamilies(Counts counts) throws Exception{
		assert(profilesList!=null) : "Family mode requires the pinned profile list";
		checkPin(profilesList, listSha80);
		final DDLPentamerBlacklist blacklist=blacklistPath==null ? null : DDLPentamerBlacklist.load(blacklistPath, blacklistPin);
		final DDLFamilySearch.Index index=DDLFamilySearch.loadForAttribution(profilesList, width, 5, 0, sketchMode,
			partitions, overlap, blacklist);
		if(index.families()!=expected){
			throw new IllegalArgumentException("Expected "+expected+" profiles, loaded "+index.families());
		}
		for(int family=0; family<index.families(); family++){counts.add(index.profileCopy(family));}
		selectedExpected=expected;
		checkPin(profilesList, listSha80);
	}

	/** Stream an explicitly pinned protein FASTA without retaining the member corpus. */
	private void readFasta(MemberCounts members) throws Exception{
		assert(input!=null) : "FASTA member input must be selected explicitly";
		checkPin(input, inputSha80);
		final FileFormat format=FileFormat.testInput(input.toString(), FileFormat.FASTA, null, false, true);
		final Streamer reader=StreamerFactory.makeStreamer(format, null, true, -1, false, true, 1);
		reader.start();
		try{
			for(ListNum<Read> batch=reader.nextList(); batch!=null; batch=reader.nextList()){
				for(Read read : batch.list){
					if(read==null || read.bases==null || read.mate!=null){
						throw new IllegalArgumentException("Require unpaired protein sequences");
					}
					members.add(read.bases);
				}
			}
		}finally{
			if(ReadWrite.closeStream(reader)){throw new IllegalStateException("Protein FASTA read failed");}
		}
		selectedExpected=expected;
		checkPin(input, inputSha80);
	}

	/** Validates the entire joined-member manifest, then consumes its selected prefix.
	 * Declared totals cover all rows; selectedExpected covers only the first chunkLimit
	 * rows. Member lengths must be 100..200 aa, and bucket IDs use rank modulo256. */
	private void readChunks(MemberCounts members) throws Exception{
		assert(chunks!=null) : "Joined member input requires a pinned manifest";
		checkPin(chunks, chunksSha80);
		final ByteFile file=ByteFile.makeByteFile1(chunks.toString(), false);
		final LineParser1 lp=new LineParser1('\t');
		final HashSet<String> paths=new HashSet<String>();
		final ArrayList<Chunk> selected=new ArrayList<Chunk>();
		long total=0;
		try{
			requireHeader(file, DDLFamilyHoldout.CHUNK_HEADER);
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				lp.set(line);
				if(lp.terms()!=6){throw new IllegalArgumentException("Joined chunk manifest requires six fields");}
				final String path=lp.parseString(2), pin=lp.parseString(3);
				final long n=lp.parseLong(4), residues=lp.parseLong(5);
				final int bucket=lp.parseInt(1);
				if(!Paths.get(path).isAbsolute() || !paths.add(path) || !validPin(pin) || n<0 || bucket<0 || bucket>=256 ||
					residues<Math.multiplyExact(100L, n) || residues>Math.multiplyExact(200L, n)){
					throw new IllegalArgumentException("Invalid/duplicate chunk manifest entry");
				}
				total=Math.addExact(total, n);
				if(manifestChunks<chunkLimit){
					selected.add(new Chunk(Paths.get(path), pin, n, residues, bucket));
					selectedExpected=Math.addExact(selectedExpected, n);
					selectedChunks++;
				}
				manifestChunks++;
			}
		}finally{
			if(file.close()){throw new IllegalStateException("Chunk manifest read failed");}
		}
		if(total!=expected || manifestChunks==0){
			throw new IllegalArgumentException("Manifest member population mismatch: "+total+" != "+expected);
		}
		partial=selectedChunks<manifestChunks;
		if(threads==1){
			for(Chunk chunk : selected){
				readChunk(chunk.path, chunk.pin, chunk.records, chunk.residues, chunk.bucket, members);
			}
		}else{readParallel(selected, members);}
		checkPin(chunks, chunksSha80);
	}

	/** Give each worker private counters; merge once after every worker successfully joins. */
	private void readParallel(ArrayList<Chunk> selected, MemberCounts result){
		assert(!selected.isEmpty()) : "A validated nonempty chunk manifest supplies parallel work";
		final int workers=Math.min(threads, selected.size());
		final AtomicInteger next=new AtomicInteger(), completed=new AtomicInteger();
		final AtomicBoolean failed=new AtomicBoolean();
		final ArrayList<ChunkWorker> pool=new ArrayList<ChunkWorker>(workers);
		for(int i=0; i<workers; i++){
			final MemberCounts local=i==0 ? result : new MemberCounts(new Counts(width), width);
			pool.add(new ChunkWorker(selected, next, completed, failed, local, outstream));
		}
		ThreadWaiter.startAndWait(pool);
		for(ChunkWorker worker : pool){
			if(worker.failure!=null){
				throw new IllegalStateException("Chunk worker failed; no output published", worker.failure);
			}
		}
		if(completed.get()!=selected.size()){throw new IllegalStateException("Not all selected chunks completed");}
		for(int i=1; i<pool.size(); i++){result.merge(pool.get(i).members);}
	}

	/** Read a pinned chunk, checking its declared population and sequence byte totals. */
	private static void readChunk(Path path, String pin, long expected, long residues, int bucket,
		MemberCounts members) throws Exception{
		assert(bucket>=0 && bucket<256) : "Joined member layout assigns family rank modulo256 to its bucket";
		checkPin(path, pin);
		final ByteFile file=ByteFile.makeByteFile1(path.toString(), false);
		final LineParser1 lp=new LineParser1('\t');
		final byte[][] buffers=new byte[201][];
		final int before=members.counts.records;
		final long basesBefore=members.residues;
		try{
			requireHeader(file, DDLFamilyHoldout.MEMBER_HEADER);
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				lp.set(line);
				if(lp.terms()!=6){throw new IllegalArgumentException("Joined member requires six fields");}
				final int length=lp.length(5), rank=lp.parseInt(0);
				if(length<100 || length>200 || rank<0 || rank%256!=bucket){
					throw new IllegalArgumentException("Joined member length/rank mismatch");
				}
				if(buffers[length]==null){buffers[length]=new byte[length];}
				System.arraycopy(line, line.length-length, buffers[length], 0, length);
				members.add(buffers[length]);
			}
		}finally{
			if(file.close()){throw new IllegalStateException("Joined member read failed: "+path);}
		}
		if(members.counts.records-before!=expected || members.residues-basesBefore!=residues){
			throw new IllegalArgumentException("Joined member/residue count mismatch: "+path);
		}
		checkPin(path, pin);
	}

	/** Reject another TSV schema before interpreting its sequence column. */
	private static void requireHeader(ByteFile file, String header){
		if(!Arrays.equals(file.nextLine(), header.getBytes(StandardCharsets.US_ASCII))){
			throw new IllegalArgumentException("Unexpected header; required "+header);
		}
	}

	/** Emit small top-key tables and cumulative pair-occurrence coverage curves. */
	private void writeTables(Counts counts){
		assert(counts.records==selectedExpected) : "Coverage denominators must use the validated selected population";
		final ByteStreamWriter keys=writer("keys.tsv"), buckets=writer("buckets.tsv"), curve=writer("coverage.tsv"),
			all=writer("counts.tsv.gz");
		final ByteBuilder line=new ByteBuilder();
		final long[] rankCells=new long[top+1], rankPairs=new long[top+1], rankKeys=new long[top+1];
		long allCells=0, allPairs=0, allKeys=0;
		final long[] ranked=new long[65535];
		try{
			keys.println("position\trank\tkey\tpacked_key\trecords\tmatching_register_pairs");
			all.println("packed_key\trecords");
			buckets.println("position\trecords\tempty\toccupied\tdistinct_keys\tmatching_register_pairs");
			curve.println("tracked_per_position\tkey_capacity_bytes\ttracked_keys\tall_keys\tcovered_cells\tall_cells\tcell_fraction\tcovered_pair_occurrences\tall_pair_occurrences\tpair_occurrence_fraction");
			for(int bucket=0; bucket<width; bucket++){
				final int distinct=counts.rank(bucket, ranked);
				long cells=0, pairs=0;
				for(int rank=1; rank<=distinct; rank++){
					final long packed=ranked[distinct-rank];
					final int key=key(packed), frequency=frequency(packed);
					final long matches=pairCount(frequency);
					cells+=frequency;
					pairs=Math.addExact(pairs, matches);
					line.clear().append((bucket<<16)|key).tab().append(frequency).nl();
					all.print(line);
					if(rank<=top){
						rankCells[rank]=Math.addExact(rankCells[rank], frequency);
						rankPairs[rank]=Math.addExact(rankPairs[rank], matches);
						rankKeys[rank]++;
						line.clear().append(bucket).tab().append(rank).tab().append(key).tab().append((bucket<<16)|key)
							.tab().append(frequency).tab().append(matches).nl();
						keys.print(line);
					}
				}
				assert(cells<=counts.records) : "Each profile has at most one nonempty key at a position";
				line.clear().append(bucket).tab().append(counts.records).tab().append(counts.records-cells)
					.tab().append(cells).tab().append(distinct).tab().append(pairs).nl();
				buckets.print(line);
				allCells=Math.addExact(allCells, cells);
				allPairs=Math.addExact(allPairs, pairs);
				allKeys+=distinct;
			}
			long coveredCells=0, coveredPairs=0, coveredKeys=0;
			for(int n=0; n<=top; n++){
				coveredCells+=rankCells[n];
				coveredPairs+=rankPairs[n];
				coveredKeys+=rankKeys[n];
				assert(coveredCells<=allCells && coveredPairs<=allPairs) : "Top-N covers a subset of the full register histogram";
				line.clear().append(n).tab().append(2L*n*width).tab().append(coveredKeys).tab().append(allKeys)
					.tab().append(coveredCells).tab().append(allCells).tab().append(fraction(coveredCells, allCells), 9)
					.tab().append(coveredPairs).tab().append(allPairs).tab().append(fraction(coveredPairs, allPairs), 9).nl();
				curve.print(line);
			}
		}finally{
			if(keys.poisonAndWait() | buckets.poisonAndWait() | curve.poisonAndWait() | all.poisonAndWait()){
				throw new IllegalStateException("Frequency table output failed; results are incomplete");
			}
		}
	}

	/** Retain all observed raw pentamer counts, sorted by packed word, without a giant sort. */
	private void writeWords(MemberCounts members){
		assert(members.validWords>=0) : "Raw word counts must conserve accepted windows from the sketch kernel";
		final ByteStreamWriter writer=writer("pentamers.tsv.gz");
		final ByteBuilder line=new ByteBuilder();
		long total=0;
		try{
			writer.println("packed_word\tpentamer\toccurrences");
			for(int word=0; word<members.words.length; word++){
				final long count=members.words[word];
				if(count==0){continue;}
				total=Math.addExact(total, count);
				line.clear().append(word).tab();
				for(int shift=20; shift>=0; shift-=5){line.append((char)AminoAcid.numberToAcid[(word>>>shift)&31]);}
				line.tab().append(count).nl();
				writer.print(line);
			}
			if(total!=members.validWords){
				throw new IllegalStateException("Raw word table failed occurrence conservation");
			}
		}finally{
			if(writer.poisonAndWait()){throw new IllegalStateException("Pentamer output failed");}
		}
	}

	/** Start a non-overwriting native writer inside the new output directory. */
	private ByteStreamWriter writer(String name){
		assert(output!=null) : "Analysis writers require the validated output directory";
		final ByteStreamWriter writer=new ByteStreamWriter(output.resolve(name).toString(), false, false, false);
		writer.start();
		return writer;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Counting Core         ----------------*/
	/*--------------------------------------------------------------*/

	/** Mutable dense histogram indexed by position in the high bits and unsigned16 key below.
	 * Values are int record counts; key zero denotes an empty register and is never counted.
	 * Borrowed input profiles are read but not retained. Callers must serialize updates. */
	public static final class Counts{
		/** Allocates one int count per unsigned16 key, including an unused zero slot per position. */
		public Counts(int width){
			if(width<1 || width>16384){throw new IllegalArgumentException("Invalid histogram width: "+width);}
			this.width=width;
			histogram=new int[width<<16];
		}

		/** Add one profile observation, excluding every empty register. */
		public void add(char[] profile){
			if(profile==null || profile.length!=width){throw new IllegalArgumentException("Mixed profile widths");}
			if(records==Integer.MAX_VALUE){throw new IllegalStateException("Profile counts exceed int range");}
			for(int b=0; b<profile.length; b++){
				final int key=profile[b];
				if(key!=0){histogram[(b<<16)|key]++;}
			}
			records++;
		}

		/** Fills and sorts the occupied-key prefix of caller-owned scratch, returning its size.
		 * Read that prefix backwards for count-descending/key-ascending order; the remaining
		 * scratch suffix is unspecified. Histogram contents are unchanged. */
		public int rank(int position, long[] scratch){
			if(position<0 || position>=width || scratch==null || scratch.length<65535){
				throw new IllegalArgumentException("Ranking needs a valid position and 65535 scratch slots");
			}
			int size=0;
			final int offset=position<<16;
			for(int key=1; key<65536; key++){
				final int count=histogram[offset|key];
				if(count>0){scratch[size++]=((long)count<<16)|(65535-key);}
			}
			Arrays.sort(scratch, 0, size);
			assert(histogram[offset]==0) : "Empty registers must never contribute matches or frequent-key candidates";
			return size;
		}

		/** Number of positions, fixed at construction. */
		private final int width;
		/** Dense position/key counts; each position reserves 65536 entries. */
		private final int[] histogram;
		/** Profiles accumulated; callers must preserve consistency with the histogram. */
		public int records=0;
	}

	/** Raw occurrence counts and unmerged member-sketch prevalence from the same sequence stream. */
	public static final class MemberCounts{
		/** Allocate one rolling-word table and one reusable sketch, not a sketch per protein. */
		public MemberCounts(Counts counts, int width){
			if(counts==null || counts.width!=width){
				throw new IllegalArgumentException("Member sketch/count widths must match");
			}
			this.counts=counts;
			sketch=new char[width];
		}

		/** Independent rolling pentamers must count the same accepted windows as production DDL. */
		public void add(byte[] sequence){
			if(sequence==null){throw new IllegalArgumentException("Missing protein sequence");}
			final long expected=DDLFamilyProfiles.sketch(sequence, 5, 0, sketch);
			int word=0, run=0;
			long added=0;
			for(byte residue : sequence){
				final int aa=residue>=0 && residue<AminoAcid.acidToNumberNoStops.length ? AminoAcid.acidToNumberNoStops[residue] : -1;
				if(aa<0){
					word=0;
					run=0;
					continue;
				}
				word=((word<<5)|aa)&WORD_MASK;
				if(++run>=5){
					words[word]=Math.addExact(words[word], 1);
					added++;
				}
			}
			if(added!=expected){throw new IllegalStateException("Raw word/sketch valid-window mismatch");}
			counts.add(sketch);
			residues=Math.addExact(residues, sequence.length);
			validWords=Math.addExact(validWords, added);
		}

		/** Adds disjoint equal-width counters after worker joins, leaving other unchanged.
		 * Requires distinct, exclusively owned counters. Overflow throws; partially merged
		 * state is not rolled back and must not be published or reused after failure. */
		public void merge(MemberCounts other){
			if(other==null || other==this || counts.width!=other.counts.width){
				throw new IllegalArgumentException("Independent equal-width member counters required");
			}
			counts.records=Math.addExact(counts.records, other.counts.records);
			for(int i=0; i<counts.histogram.length; i++){
				counts.histogram[i]=Math.addExact(counts.histogram[i], other.counts.histogram[i]);
			}
			for(int i=0; i<words.length; i++){words[i]=Math.addExact(words[i], other.words[i]);}
			residues=Math.addExact(residues, other.residues);
			validWords=Math.addExact(validWords, other.validWords);
		}

		/** Exact raw occurrence count for fixtures and future blacklist selection. */
		public long occurrences(int word){
			if(word<0 || word>=words.length){throw new IllegalArgumentException("Pentamer must fit25 bits");}
			return words[word];
		}

		/** Register counts updated together with raw word totals; exclusively owned by this counter. */
		private final Counts counts;
		/** Reusable plain K5 sketch; no input sequence is retained. */
		private final char[] sketch;
		/** Dense raw pentamer occurrence table indexed by five-bit amino codes. */
		private final long[] words=new long[1<<25];
		/** Total input residue bytes and accepted K5 windows, respectively. */
		public long residues=0, validWords=0;
	}

	/** Recover the full unsigned16 key; zero is never ranked. */
	public static int key(long packed){return 65535-(int)(packed&65535);}

	/** Recover the exact number of profile registers containing a ranked key. */
	public static int frequency(long packed){return (int)(packed>>>16);}

	/** Unordered profile-pair occurrences contributed by one occupied register value. */
	public static long pairCount(int count){
		if(count<0){throw new IllegalArgumentException("Negative frequency");}
		return (long)count*(count-1)/2;
	}

	/** Empty populations have explicit zero coverage, with their denominator retained. */
	public static double fraction(long numerator, long denominator){
		if(numerator<0 || numerator>denominator){
			throw new IllegalArgumentException("Coverage must be a nonnegative subset");
		}
		return denominator==0 ? 0 : (double)numerator/denominator;
	}

	/** Streams and verifies the supplied path against its required legacy SHA80 value.
	 * Used for manifests and sequence/chunk files, not only small metadata; referenced
	 * chunks are checked separately. The local stream is closed on success or failure. */
	private static void checkPin(Path path, String expected) throws Exception{
		assert(expected.length()==20) : "Manifest pins use the project's final-80-bit SHA256 convention";
		final MessageDigest digest=MessageDigest.getInstance("SHA-256");
		final byte[] buffer=new byte[65536];
		try(InputStream input=Files.newInputStream(path)){
			for(int n=input.read(buffer); n!=-1; n=input.read(buffer)){digest.update(buffer, 0, n);}
		}
		final byte[] hash=digest.digest();
		final StringBuilder actual=new StringBuilder(20);
		for(int i=hash.length-10; i<hash.length; i++){
			actual.append(Character.forDigit((hash[i]>>>4)&15, 16)).append(Character.forDigit(hash[i]&15, 16));
		}
		if(!actual.toString().equals(expected)){throw new IllegalArgumentException("Manifest sha80 mismatch: "+path);}
	}

	/** Require exact project digest spelling before any expensive input pass. */
	private static boolean validPin(String value){return value!=null && value.matches("[0-9a-f]{20}");}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Family-list input and fresh output directory. */
	private Path profilesList, output;
	/** Mutually exclusive member FASTA and joined-member manifest inputs. */
	private Path input, chunks;
	/** Required family-list fingerprint in legacy family mode. */
	private String listSha80;
	/** Required fingerprint for the selected member input form. */
	private String inputSha80, chunksSha80;
	/** Population to count: family profiles or individual members. */
	private String mode="families";
	/** Family profile recipe and optional blacklist binding; members remain plain/unfiltered. */
	private String sketchMode=DDLFamilyProfiles.SKETCH_PLAIN, blacklistPin;
	/** Optional family-mode blacklist file. */
	private Path blacklistPath;
	/** Positional layout bound to selected family profiles. */
	private int partitions=4, overlap=8;
	/** Register width, per-position report depth and full expected record population. */
	private int width=1024, top=64, expected=79563;
	/** Maximum independently reading chunk workers; decoder threads remain one. */
	private int threads=1;
	/** Manifest-prefix limit and observed selected/full chunk counts. */
	private int chunkLimit=Integer.MAX_VALUE, selectedChunks=0, manifestChunks=0;
	/** Expected records in the selected input, distinct from full-manifest expected. */
	private long selectedExpected=0;
	/** Explicit member population supplied; selected prefix omits manifest rows. */
	private boolean expectedSet=false, partial=false;
	/** Diagnostic stream from PreParser, closed by main through Shared.closeStream. */
	private final PrintStream outstream;
	/** Five amino-acid codes at five bits each. */
	private static final int WORD_MASK=(1<<25)-1;

	/*--------------------------------------------------------------*/
	/*----------------        Chunk Workers         ----------------*/
	/*--------------------------------------------------------------*/
	/** Immutable chunk metadata is shared; sequence buffers and histograms are never shared. */
	private static final class Chunk{
		/** Retains immutable values from one validated manifest row. */
		Chunk(Path path_, String pin_, long records_, long residues_, int bucket_){
			assert(path_!=null && records_>=0) : "Chunks originate in the validated manifest";
			path=path_;
			pin=pin_;
			records=records_;
			residues=residues_;
			bucket=bucket_;
		}
		/** Absolute chunk input path. */
		final Path path;
		/** Legacy content binding for this chunk. */
		final String pin;
		/** Expected member rows and total amino-acid bytes. */
		final long records, residues;
		/** Required family-rank modulo256 bucket. */
		final int bucket;
	}

	/** An atomic claim occurs once per file, never per residue or counted key. */
	private static final class ChunkWorker extends Thread{
		/** Borrows fixed chunk metadata, shared claim/completion state and one exclusive counter. */
		ChunkWorker(ArrayList<Chunk> chunks_, AtomicInteger next_, AtomicInteger completed_,
				AtomicBoolean failed_, MemberCounts members_, PrintStream out_){
			assert(chunks_!=null && members_!=null) : "Workers require immutable work metadata and private counters";
			chunks=chunks_;
			next=next_;
			completed=completed_;
			failed=failed_;
			members=members_;
			out=out_;
		}

		/** Capture all failures for the joining thread; no consumer waits on a dead worker. */
		@Override
		public void run(){
			try{
				while(!failed.get()){
					final int i=next.getAndIncrement();
					if(i>=chunks.size()){break;}
					final Chunk chunk=chunks.get(i);
					readChunk(chunk.path, chunk.pin, chunk.records, chunk.residues, chunk.bucket, members);
					final int done=completed.incrementAndGet();
					if((done&1023)==0 || done==chunks.size()){out.println("Completed chunks: "+done+"/"+chunks.size());}
				}
			}catch(Throwable problem){
				failure=problem;
				failed.set(true);
			}
		}
		/** Fixed chunk metadata shared by workers; no worker modifies this list. */
		final ArrayList<Chunk> chunks;
		/** Shared next-row claim and successfully completed chunk counts. */
		final AtomicInteger next, completed;
		/** Shared failure flag checked before each new chunk claim. */
		final AtomicBoolean failed;
		/** Counter exclusively updated by this worker, merged by its joining caller. */
		final MemberCounts members;
		/** Shared diagnostic stream for coarse completion messages. */
		final PrintStream out;
		/** First escaping failure for this worker, read only after joining. */
		Throwable failure;
	}
}
