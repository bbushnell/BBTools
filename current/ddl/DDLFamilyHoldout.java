package ddl;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Random;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import prot.DigestSuffix;
import prot.MemberLengthRange;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;

/** Exact per-family member holdout and training-only legacy16 profile kernels.
 * Input order is the pinned winner-assignment/source join order, not old retained.faa.
 * Reads one bounded bucket of pinned winner-member TSV chunks in canonical batch order.
 * Each family contributes exactly floor(memberCount/10) held-out records; all
 * remaining records vote in K5, seed-zero, 1024-register profiles. Families with
 * fewer than ten members are training-only; zero-member families retain split
 * metadata but have no profile. Both literal and nonempty modal rows are written.
 * The command retains one bucket's TSV rows and builds one family at a time.
 * Its byte cap counts retained line bytes, excluding array/object overhead and
 * profile storage; it is not a total-heap bound.
 * Split and FamilyBuilder are mutable, single-owner objects. The command also
 * changes process-wide thread, amino-input and compression settings without
 * restoring them, so run() is not an isolated concurrent library operation.
 * @author Sayu
 */
public final class DDLFamilyHoldout{

	/** Static command and nested per-family helpers only. */
	private DDLFamilyHoldout(){}

	/** Applies standard preprocessing; reports any failure and exits nonzero. */
	public static void main(String[] args){
		try{
			final PreParser pp=new PreParser(args, DDLFamilyHoldout.class, false);
			try{run(pp.args);}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	/** Runs the same validated command for the launcher and integration fixtures.
	 * Requires a dense ordered family roster, complete batch-then-bucket manifest,
	 * matching declared counts/pins and a fresh output directory under an existing
	 * parent. Selected members must have canonical amino residues and lie within
	 * minmemberaa..maxmemberaa (100..200 when both options are omitted).
	 * Writes profiles, split membership, held-out FASTA, recipe and count summaries;
	 * PASS is written last. A failure can leave partial output without PASS.
	 * Sets threads=1 and amino input, and disables pigz/unpigz process-wide without restoration.
	 * @param args Unique supported flag=value arguments; see the launcher for required keys.
	 * @throws Exception If arguments, input bindings, conservation or I/O checks fail.
	 */
	public static void run(String[] args) throws Exception{
		final HashMap<String,String> options=new HashMap<String,String>();
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value");
			final String key=pair[0].toLowerCase(Locale.ROOT);
			require(Arrays.asList("families","familiessha80","chunks","chunkssha80","bucket","buckets",
				"expectedfamilies","expectedmembers","expectedbatches","maxbucketmembers","maxbucketbytes","out","splitseed","minmemberaa","maxmemberaa").contains(key), "Unknown option: "+key);
			require(options.put(key, pair[1])==null, "Duplicate option: "+key);
		}
		final MemberLengthRange memberRange=MemberLengthRange.from(options);
		final int bucket=integer(options, "bucket"), buckets=integer(options, "buckets");
		final int nf=integer(options, "expectedfamilies"), nb=integer(options, "expectedbatches");
		final long total=number(options, "expectedmembers"), maxMembers=number(options, "maxbucketmembers"), maxBytes=number(options, "maxbucketbytes");
		final long seed=Long.parseLong(required(options, "splitseed"));
		require(buckets>0 && bucket<buckets && nf>0 && nb>0 && (long)buckets*nb<=Integer.MAX_VALUE && maxMembers>0 && maxBytes>0, "Invalid bucket/count/cap dimensions");
		final Path output=Paths.get(required(options, "out")).toAbsolutePath();
		require(!Files.exists(output) && Files.isDirectory(output.getParent()), "Output must be fresh under an existing parent");
		final Path families=Paths.get(required(options, "families")), chunks=Paths.get(required(options, "chunks"));
		pin(families, required(options, "familiessha80")); pin(chunks, required(options, "chunkssha80"));
		Shared.setThreads(1); Shared.AMINO_IN=true; ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		final String[] ids=new String[nf]; final int[] counts=new int[nf];
		final LineParser1 lp=new LineParser1('\t'); final HashSet<String> seen=new HashSet<String>();
		final ByteFile metadata=ByteFile.makeByteFile1(families.toString(), false);
		long countSum=0, selected=0; int rank=0;
		try{
			require(Arrays.equals(metadata.nextLine(), FAMILY_HEADER.getBytes(StandardCharsets.US_ASCII)), "Wrong family count header");
			for(byte[] line=metadata.nextLine(); line!=null; line=metadata.nextLine()){
				lp.set(line);
				require(lp.terms()==3 && rank<nf && lp.parseInt(0)==rank, "Family roster must be dense and ordered");
				final String id=lp.parseString(1); final int n=lp.parseInt(2);
				require(id.matches("[A-Za-z0-9_.-]+") && seen.add(id) && n>=0, "Unsafe/duplicate family or negative count");
				ids[rank]=id; counts[rank]=n; countSum=Math.addExact(countSum, n);
				if(rank%buckets==bucket){selected=Math.addExact(selected, n);}
				rank++;
			}
		}finally{require(!metadata.close(), "Family count input failed");}
		require(rank==nf && countSum==total && selected<=maxMembers, "Family totals or bucket member cap failed");
		final byte[][][] rows=new byte[nf][][]; final int[] loaded=new int[nf];
		for(int i=bucket; i<nf; i+=buckets){rows[i]=new byte[counts[i]][];}
		final ByteFile parts=ByteFile.makeByteFile1(chunks.toString(), false);
		final HashSet<String> paths=new HashSet<String>();
		final HashSet<String> genomes=new HashSet<String>();
		long selectedRows=0, selectedBases=0, retainedBytes=0, declared=0; int part=0;
		try{
			require(Arrays.equals(parts.nextLine(), CHUNK_HEADER.getBytes(StandardCharsets.US_ASCII)), "Wrong chunk manifest header");
			for(byte[] line=parts.nextLine(); line!=null; line=parts.nextLine()){
				lp.set(line);
				require(lp.terms()==6 && part<(long)nb*buckets, "Wrong or extra chunk manifest row");
				final int b=part%buckets, batch=part/buckets;
				require(lp.termEquals(String.format(Locale.ROOT, "b%02d", batch), 0) && lp.parseInt(1)==b, "Chunks must be ordered batch then bucket");
				final String path=lp.parseString(2), hash=lp.parseString(3);
				final long n=lp.parseLong(4), bases=lp.parseLong(5);
				require(Paths.get(path).isAbsolute() && paths.add(path) && hash.matches("[0-9a-f]{20}") && n>=0
					&& memberRange.containsTotal(n, bases), "Invalid chunk binding/counts or member range "+memberRange);
				declared=Math.addExact(declared, n);
				if(b==bucket){
					final long[] got=loadChunk(Paths.get(path), hash, rows, loaded, ids, bucket, buckets, maxBytes-retainedBytes, genomes, memberRange);
					require(got[0]==n && got[1]==bases, "Chunk member/residue counts differ from pinned manifest");
					selectedRows+=n; selectedBases+=bases; retainedBytes+=got[2];
					System.err.println("holdout_chunk_loaded\t"+batch+"\t"+bucket+"\t"+selectedRows);
				}
				part++;
			}
		}finally{require(!parts.close(), "Chunk manifest input failed");}
		require(part==(long)nb*buckets && declared==total && selectedRows==selected, "Chunk/global member conservation failed");
		for(int i=bucket; i<nf; i+=buckets){require(loaded[i]==counts[i], "Family join count differs: "+ids[i]);}
		Files.createDirectory(output);
		final ByteStreamWriter profiles=writer(output.resolve("profiles.tsv")), split=writer(output.resolve("split.tsv")), held=writer(output.resolve("heldout.faa"));
		long train=0, test=0, trainBases=0, testBases=0, trainWords=0; int nonempty=0, empty=0, small=0;
		final ByteBuilder lineOut=new ByteBuilder(); final byte[][] buffers=new byte[201][];
		try{
			profiles.println(DDLFamilyProfiles.HEADER16); split.println(SPLIT_HEADER);
			for(int i=bucket; i<nf; i+=buckets){
				final FamilyBuilder builder=new FamilyBuilder(counts[i], seed, "tail", ids[i]);
				for(int ordinal=0; ordinal<rows[i].length; ordinal++){
					final byte[] row=rows[i][ordinal]; lp.set(row); final int len=lp.length(5);
					if(buffers[len]==null){buffers[len]=new byte[len];}
					final byte[] sequence=buffers[len];
					// Sequence is the last field; avoid a fresh Read/String/sequence allocation per member.
					System.arraycopy(row, row.length-len, sequence, 0, len);
					if(builder.add(sequence)){
						lineOut.clear().append(">ddlholdout|").append(i).append('|').append(ids[i]).append('|').append(ordinal).append(" query_id=");
						lp.appendTerm(lineOut, 2).append(" genome_id="); lp.appendTerm(lineOut, 3).append(" source_ordinal="); lp.appendTerm(lineOut, 4);
						lineOut.nl().append(sequence).nl(); held.print(lineOut);
					}
					rows[i][ordinal]=null;
				}
				builder.finish(); rows[i]=null;
				if(builder.profile==null){empty++;}else{writeProfile(profiles, ids[i], builder, lineOut); nonempty++;}
				if(counts[i]>0 && counts[i]<10){small++;}
				lineOut.clear().append(i).tab().append(ids[i]).tab().append(counts[i]).tab().append(builder.split.train).tab().append(builder.split.test)
					.tab().append(builder.trainResidues()).tab().append(builder.testResidues()).tab().append(builder.trainWords()).tab()
					.append(counts[i]==0 ? "NO_ASSIGNED_MEMBERS" : counts[i]<10 ? "TRAIN_ONLY" : "SPLIT").tab();
				for(byte value : builder.split.membership()){lineOut.append(HEX[(value&255)>>>4]).append(HEX[value&15]);}
				split.println(lineOut);
				train+=builder.split.train; test+=builder.split.test; trainBases+=builder.trainResidues(); testBases+=builder.testResidues();
				trainWords=Math.addExact(trainWords, builder.trainWords());
				System.err.println("holdout_family_complete\t"+i+"\t"+ids[i]+"\t"+builder.split.train+"\t"+builder.split.test);
			}
		}finally{
			final boolean failed=profiles.poisonAndWait() | split.poisonAndWait() | held.poisonAndWait();
			require(!failed, "Held-out profile/split/query output failed; no PASS");
		}
		require(train+test==selected && trainBases+testBases==selectedBases, "Train/test member/residue conservation failed");
		pin(families, options.get("familiessha80")); pin(chunks, options.get("chunkssha80"));
		Files.write(output.resolve("recipe.txt"), (String.join("\n", args)+"\nsplit=uniform_exact_floor_n_over_10\nfamily_seed=SHA256_DDLFamilyHoldout-v1_last8bytes_JavaRandom\norder=canonical_batch_genome_source\nprofile=k5_width1024_legacy16_seed0\n").getBytes(StandardCharsets.US_ASCII));
		Files.write(output.resolve("summary.tsv"), ("namespace\tfamilies_sha80\tchunks_sha80\tk\tseed\texponent\tregisterbits\tbuckets\tshard\tshards\tmanifest_families\tbuilt_families\tmembers\tresidues\tvalid_words\tsplitseed\n"+
			"tail\t"+options.get("familiessha80")+"\t"+options.get("chunkssha80")+"\t5\t0\t6\t16\t1024\t"+bucket+"\t"+buckets+"\t"+nf+"\t"+nonempty+"\t"+train+"\t"+trainBases+"\t"+trainWords+"\t"+seed+"\n").getBytes(StandardCharsets.US_ASCII));
		Files.write(output.resolve("split_summary.tsv"), ("bucket\tassigned\ttrain\ttest\ttrain_residues\ttest_residues\tprofiles\tzero_families\ttrain_only_families\tretained_row_bytes\n"+
			bucket+"\t"+selected+"\t"+train+"\t"+test+"\t"+trainBases+"\t"+testBases+"\t"+nonempty+"\t"+empty+"\t"+small+"\t"+retainedBytes+"\n").getBytes(StandardCharsets.US_ASCII));
		Files.write(output.resolve("PASS"), "DDL_FAMILY_HOLDOUT_PASS\n".getBytes(StandardCharsets.US_ASCII));
	}

	/** Validates one selected chunk and retains its original row arrays by family.
	 * The shared genome set rejects a repeated genome block, including across chunks;
	 * source ordinals must strictly increase within each genome. Loaded counters and
	 * row slots are updated as records pass, so a failure does not roll state back.
	 * @param byteCap Remaining retained-line-byte allowance for the bucket.
	 * @return Member count, residue count and retained line bytes, in that order.
	 */
	private static long[] loadChunk(Path path, String hash, byte[][][] rows, int[] loaded, String[] ids, int bucket, int buckets, long byteCap, HashSet<String> genomes, MemberLengthRange memberRange){
		pin(path, hash);
		final ByteFile input=ByteFile.makeByteFile1(path.toString(), false); final LineParser1 lp=new LineParser1('\t');
		long count=0, bases=0, bytes=0, lastOrdinal=-1; String genome=null;
		try{
			require(Arrays.equals(input.nextLine(), MEMBER_HEADER.getBytes(StandardCharsets.US_ASCII)), "Wrong joined member header");
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line); require(lp.terms()==6, "Joined member must have six fields");
				final int rank=lp.parseInt(0), len=lp.length(5);
				require(rank>=0 && rank<ids.length && rank%buckets==bucket && lp.termEquals(ids[rank], 1), "Joined family rank/identity/bucket mismatch");
				if(loaded[rank]>=rows[rank].length || !memberRange.contains(len) || lp.parseLong(4)<0){
					throw new IllegalArgumentException("Extra member, sequence outside "+memberRange+" or invalid source ordinal; length="+len);
				}
				for(int term=2; term<=3; term++){
					require(lp.length(term)>0, "Empty query/genome identity");
					for(int j=0; j<lp.length(term); j++){require(lp.parseByte(term, j)>32 && lp.parseByte(term, j)<127, "Whitespace/non-ASCII query/genome identity");}
				}
				if(genome==null || !lp.termEquals(genome, 3)){
					genome=lp.parseString(3); require(genomes.add(genome), "Reordered/repeated genome across joined chunks"); lastOrdinal=-1;
				}
				final long sourceOrdinal=lp.parseLong(4);
				require(sourceOrdinal>lastOrdinal, "Duplicate or reordered original member ordinal"); lastOrdinal=sourceOrdinal;
				require(lp.length(2)>genome.length()+1 && lp.parseByte(2, genome.length())=='|', "Query ID lacks genome namespace");
				for(int j=0; j<genome.length(); j++){require(lp.parseByte(2, j)==genome.charAt(j), "Query/genome namespace mismatch");}
				for(int j=line.length-len; j<line.length; j++){
					final int aa=line[j]&255;
					require(aa<128 && AminoAcid.acidToNumberNoStops[aa]>=0, "Ambiguous/stop residue in canonical assigned tail");
				}
				bytes=Math.addExact(bytes, line.length); require(bytes<=byteCap, "Bucket retained-byte cap exceeded");
				rows[rank][loaded[rank]++]=line; count++; bases+=len;
			}
		}finally{require(!input.close(), "Joined member read failed");}
		pin(path, hash);
		return new long[]{count,bases,bytes};
	}

	/** Emits literal/nonempty modal values, support and occupancy using reusable scratch. */
	private static void writeProfile(ByteStreamWriter writer, String id, FamilyBuilder builder, ByteBuilder line){
		assert(builder.profile!=null && builder.profile.members()>0) : "Only nonempty training families can be searchable profiles";
		for(boolean literal : new boolean[]{true,false}){
			line.clear().append("tail\t").append(id).append("\t1024\t").append(builder.split.train).tab().append(builder.trainResidues()).tab().append(builder.trainWords())
				.tab().append(literal ? "literal" : "nonempty");
			for(int kind=0; kind<3; kind++){
				line.tab();
				for(int b=0; b<1024; b++){
					if(b>0){line.append(',');}
					line.append(kind==0 ? (int)builder.profile.value(b, literal) : kind==1 ? builder.profile.votes(b, literal) : builder.profile.occupied(b));
				}
			}
			writer.println(line);
		}
	}

	/** Starts a fresh-file writer; the caller must poison and join it. */
	private static ByteStreamWriter writer(Path path){
		final ByteStreamWriter result=new ByteStreamWriter(path.toString(), false, false, false); result.start(); return result;
	}
	/** Verifies the legacy command's explicit absolute-path/content binding. */
	private static void pin(Path path, String expected){
		require(path.isAbsolute() && expected.matches("[0-9a-f]{20}") && Files.isRegularFile(path), "Invalid input path/pin");
		require(DigestSuffix.file(path.toString()).equals(expected), "Input sha80 mismatch: "+path);
	}
	/** Gets a mandatory nonempty option without changing the option map. */
	private static String required(HashMap<String,String> options, String key){
		final String value=options.get(key); require(value!=null && !value.isEmpty(), "Missing option: "+key); return value;
	}
	/** Parses a mandatory nonnegative decimal long. */
	private static long number(HashMap<String,String> options, String key){
		final long value=Long.parseLong(required(options, key)); require(value>=0, "Negative option: "+key); return value;
	}
	/** Parses a nonnegative count and rejects values outside the int range. */
	private static int integer(HashMap<String,String> options, String key){return Math.toIntExact(number(options, key));}

	/** Exact schemas for the roster, chunk manifest, joined members and split output. */
	public static final String FAMILY_HEADER="roster_rank\tfamily_id\tassigned";
	public static final String CHUNK_HEADER="batch_id\tbucket\tfile\tsha80\tmembers\tresidues";
	public static final String MEMBER_HEADER="family_rank\tfamily_id\tquery_id\tgenome_id\tsource_ordinal\tsequence";
	public static final String SPLIT_HEADER="roster_rank\tfamily_id\ttotal\ttrain\ttest\ttrain_residues\ttest_residues\ttrain_words\tstatus\ttest_bits_hex";
	private static final char[] HEX="0123456789abcdef".toCharArray();

	/** Domain-separated deterministic seed; Java Random supplies its specified48-bit stream.
	 * Family identity, not worker/shard number, controls the stream. Input order remains pinned.
	 * Uses the final eight SHA-256 bytes in big-endian order; Random subsequently
	 * initializes its 48-bit state from that long. This derives a seed, not a file pin.
	 * @param seed Master split seed; all long values are supported.
	 * @param namespace Nonempty identifier containing only letters, digits, underscore, dot or hyphen.
	 * @param family Family identifier with the same character restrictions.
	 * @return Deterministic signed long supplied to the per-family Random.
	 */
	public static long familySeed(long seed, String namespace, String family){
		require(namespace!=null && namespace.matches("[A-Za-z0-9_.-]+"), "Unsafe split namespace");
		require(family!=null && family.matches("[A-Za-z0-9_.-]+"), "Unsafe split family ID");
		try{
			final MessageDigest digest=MessageDigest.getInstance("SHA-256");
			final String binding="DDLFamilyHoldout-v1\n"+seed+"\n"+namespace+"\n"+family+"\n";
			final byte[] bytes=digest.digest(binding.getBytes(StandardCharsets.US_ASCII));
			long value=0;
			for(int i=bytes.length-8; i<bytes.length; i++){value=(value<<8)|(bytes[i]&255L);}
			return value;
		}catch(NoSuchAlgorithmException e){throw new IllegalStateException("Required SHA-256 unavailable", e);}
	}

	/** Uniform subset of exactly floor(n/10) ordinals, without sorting/shuffling records.
	 * At every step P(test)=remainingTest/remainingMembers; selection is independent of bases.
	 * Bits are little-bit-order within bytes; unused high bits stay zero.
	 */
	public static final class Split{
		/** Allocates membership bits and a deterministic stream for a nonnegative quota. */
		public Split(int members, long seed, String namespace, String family){
			require(members>=0, "Assigned member count cannot be negative");
			total=members; test=members/10; train=members-test;
			remaining=members; remainingTest=test;
			bits=new byte[(int)(((long)members+7)/8)];
			random=new Random(familySeed(seed, namespace, family));
		}

		/** Consumes the next ordinal and returns true for held-out, false for training.
		 * Calls beyond the declared member count fail; records themselves are not inspected.
		 */
		public boolean next(){
			require(remaining>0, "Extra assigned member beyond pinned family count");
			assert(remainingTest>=0 && remainingTest<=remaining) :
				"Remaining-count sampling must preserve feasibility of the exact test quota";
			final boolean held=remainingTest>0 && (remainingTest==remaining || random.nextInt(remaining)<remainingTest);
			if(held){bits[position>>>3]|=(byte)(1<<(position&7)); remainingTest--;}
			position++; remaining--;
			return held;
		}

		/** Checks complete consumption and exact quota; does not reset or mutate the split. */
		public void finish(){
			require(remaining==0 && remainingTest==0 && position==total,
				"Missing assigned members or unfulfilled exact held-out quota");
		}

		/** Returns a previously consumed ordinal's bit; future/out-of-range ordinals fail. */
		public boolean heldOut(int ordinal){
			require(ordinal>=0 && ordinal<position, "Split membership requested before assignment or outside family");
			return (bits[ordinal>>>3]&(1<<(ordinal&7)))!=0;
		}

		/** Returns a detached bit-array copy only after all members have been consumed. */
		public byte[] membership(){finish(); return bits.clone();}
		/** Number of ordinals assigned so far. */
		public int consumed(){return position;}
		/** Declared total and fixed training/test quotas, not incremental counters. */
		public final int total, train, test;
		private final byte[] bits;
		private final Random random;
		private int remaining, remainingTest, position;
	}

	/** One active family keeps vote-map memory bounded. Held-out reads never reach sketch/add.
	 * The caller owns query serialization, using the returned split bit and original identity.
	 * Sequence arrays are consumed synchronously and never retained; callers may reuse
	 * them after add returns. This helper does not apply the command's length/residue
	 * validation. Its exposed split/profile are live mutable objects, not snapshots;
	 * callers must not mutate them independently while building the family.
	 */
	public static final class FamilyBuilder{
		/** Creates the exact split and an empty profile, or a null profile for zero members. */
		public FamilyBuilder(int members, long splitSeed, String namespace, String family){
			split=new Split(members, splitSeed, namespace, family);
			profile=members==0 ? null : new DDLFamilyProfiles.Profile(1024);
		}

		/** Borrows an identified unpaired read's bases; the ID does not affect selection. */
		public boolean add(Read read){
			require(read!=null && read.bases!=null && read.mate==null && read.id!=null,
				"Assigned member must be an identified unpaired sequence");
			return add(read.bases);
		}

		/** Consumes one member; true means held-out, false means added to training.
		 * Byte-native joined TSV consumers reuse length-specific sequence buffers.
		 * Training words count valid forward K5 occurrences, not distinct words.
		 * A failed add may already have advanced the split; discard failed builders.
		 */
		public boolean add(byte[] sequence){
			require(sequence!=null, "Missing assigned sequence");
			final boolean held=split.next();
			if(held){testResidues=Math.addExact(testResidues, sequence.length);}
			else{
				trainResidues=Math.addExact(trainResidues, sequence.length);
				trainWords=Math.addExact(trainWords, DDLFamilyProfiles.sketch(sequence, 5, 0, sketch));
				profile.add(sketch);
			}
			return held;
		}

		/** Checks exact split completion and training-vote conservation without resetting state. */
		public void finish(){
			split.finish();
			require(profile==null ? split.total==0 : profile.members()==split.train,
				"Training votes must equal assigned members minus held-out members");
			assert(split.total==split.train+split.test) : "Every assigned record belongs to exactly one partition";
		}

		/** Sum of lengths of training sequences consumed so far. */
		public long trainResidues(){return trainResidues;}
		/** Sum of lengths of held-out sequences consumed so far. */
		public long testResidues(){return testResidues;}
		/** Valid forward K5 occurrence count from training sequences consumed so far. */
		public long trainWords(){return trainWords;}
		/** Live ordinal assignment state owned by this builder. */
		public final Split split;
		/** Null for a zero-assigned family: identity remains in metadata, no fabricated profile. */
		public final DDLFamilyProfiles.Profile profile;
		private final char[] sketch=new char[1024];
		private long trainResidues, testResidues, trainWords;
	}

	/** Unconditional command/helper contract check, independent of assertion settings. */
	private static void require(boolean condition, String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}
}
