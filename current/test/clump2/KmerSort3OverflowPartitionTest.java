package test.clump2;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.HashMap;
import java.util.HashSet;

import clump2.Clump;
import clump2.KmerComparator;
import clump2.KmerSort;
import clump2.KmerSort3;
import shared.Shared;
import stream.Read;

/**
 * KmerSort3 skew-streaming regression for the next-pass partition key.
 * @author Jean
 */
public final class KmerSort3OverflowPartitionTest {

	/** Runs the normal split path and the overflow-streaming path. */
	public static void main(String[] args) throws Exception{
		check(args.length==2, "Expected baseline|fixed and a fresh fixture directory");
		final boolean fixed=args[0].equals("fixed");
		check(fixed || args[0].equals("baseline"), "Unknown expectation: "+args[0]);

		final Path root=Paths.get(args[1]);
		Files.createDirectory(root);

		final StaticState state=new StaticState();
		try{
			state.setDeterministic();

			System.out.println("case\treads\tmisplaced");
			final HashMap<String, Integer> expected=verifyNormalSplitOracle();
			final HashMap<String, String> bases=expectedBases();

			final int normalMisplaced=runKmerSort3Case(root.resolve("normal"),
					expected, bases, true);
			check(normalMisplaced==0,
					"Expected 0 misplaced forceprocess reads, got "+normalMisplaced);
			System.out.println("normal-control\t"+expected.size()+'\t'+normalMisplaced);

			final int misplaced=runKmerSort3Case(root.resolve("overflow"),
					expected, bases, false);
			check((misplaced==0)==fixed,
					"Expected "+(fixed ? "0" : ">0")+" misplaced overflow reads, got "+misplaced);
			System.out.println("overflow-partition\t"+expected.size()+'\t'+misplaced);
		}finally{
			state.restore();
		}
	}

	/** Verifies the selected literals differ between raw and normal next-pass buckets. */
	private static HashMap<String, Integer> verifyNormalSplitOracle() throws Exception{
		final KmerComparator current=new KmerComparator(K, 1, 1, 4, false, false);
		final KmerComparator next=new KmerComparator(K, 2, 0, 4, false, false);
		final HashMap<String, Integer> expected=new HashMap<String, Integer>();
		int rawMismatches=0;

		for(int i=0; i<READS.length; i++){
			final String id="r"+i;
			final Read rawRead=read(id, READS[i], i);
			final long rawPivot=current.hash(rawRead, null, 0, false);
			final Read normalRead=read(id, READS[i], i);
			next.hash(normalRead, null, 0, true);
			final long keyKmer=KEY_KMER.getLong(normalRead.obj());
			final int rawGroup=(int)(rawPivot%GROUPS);
			final int normalGroup=(int)(next.hash(keyKmer)%GROUPS);
			check(normalGroup==EXPECTED_GROUPS[i],
					"Unexpected literal bucket for "+id+": expected "+
					EXPECTED_GROUPS[i]+", got "+normalGroup);
			if(rawGroup!=normalGroup){rawMismatches++;}
			expected.put(id, EXPECTED_GROUPS[i]);
		}

		check(rawMismatches==READS.length,
				"All selected literals must exercise the raw-vs-normal partition distinction");
		System.out.println("literal-oracle\t"+READS.length+"\t0");
		return expected;
	}

	/** Runs KmerSort3 with the selected literals and returns mispartitioned reads. */
	private static int runKmerSort3Case(Path root, HashMap<String, Integer> expected,
			HashMap<String, String> expectedBases, boolean forceProcess) throws Exception{
		Files.createDirectory(root);
		final Path inPattern=root.resolve("in_%.fq");
		final Path outPattern=root.resolve("out_%.fq");

		writeInputGroup(root, 0, fastq());
		writeInputGroup(root, 1, "");
		writeInputGroup(root, 2, "");

		KmerSort3.main(Long.MAX_VALUE/4, 1, 2, new String[] {
				"in="+inPattern,
				"out="+outPattern,
				"groups="+GROUPS,
				"k="+K,
				"seed=1",
				"border=1",
				"hashes=4",
				"ow=t",
				"qin=33",
				"qout=33",
				"interleaved=f",
				"ecco=f",
				"dedupe=f",
				"condense=f",
				"correct=f",
				"repair=f",
				"namesort=f",
				"forceprocess="+forceProcess,
				"fetchthreads=1",
				"verbose=f"
		});

		final HashSet<String> seen=new HashSet<String>();
		int misplaced=0;
		for(int group=0; group<GROUPS; group++){
			final Path output=root.resolve("out_"+group+".fq");
			for(String id : readIds(output, expectedBases)){
				check(seen.add(id), "Duplicate output read: "+id);
				final Integer expectedGroup=expected.get(id);
				check(expectedGroup!=null, "Unexpected output read: "+id);
				if(expectedGroup.intValue()!=group){misplaced++;}
			}
		}
		check(seen.size()==expected.size(),
				"Expected "+expected.size()+" output reads, observed "+seen.size());
		return misplaced;
	}

	/** Writes one KmerSort3 input shard. */
	private static void writeInputGroup(Path root, int group, String contents)
			throws Exception{
		Files.write(root.resolve("in_"+group+".fq"),
				contents.getBytes(StandardCharsets.US_ASCII),
				StandardOpenOption.CREATE_NEW);
	}

	/** Returns the read identifiers from a single FASTQ file. */
	private static HashSet<String> readIds(final Path path,
			final HashMap<String, String> expectedBases) throws Exception{
		check(Files.exists(path), "Missing output file: "+path);
		final java.util.List<String> lines=Files.readAllLines(path, StandardCharsets.US_ASCII);
		check(lines.size()%4==0, "Malformed FASTQ output line count: "+path);

		final HashSet<String> ids=new HashSet<String>();
		for(int i=0; i<lines.size(); i+=4){
			final String line=lines.get(i);
			check(line.startsWith("@"), "Malformed FASTQ header in "+path+": "+line);
			final String id=line.substring(1).split("\\s+", 2)[0];
			final String bases=lines.get(i+1);
			final String separator=lines.get(i+2);
			final String quals=lines.get(i+3);
			check(separator.startsWith("+"), "Malformed FASTQ separator in "+path+": "+separator);
			check(bases.equals(expectedBases.get(id)), "Unexpected bases for "+id+": "+bases);
			check(quals.length()==bases.length(), "Quality length mismatch for "+id);
			for(int j=0; j<quals.length(); j++){
				check(quals.charAt(j)=='I', "Unexpected quality for "+id+": "+quals);
			}
			check(ids.add(id), "Duplicate output read in "+path+": "+id);
		}
		return ids;
	}

	/** Builds the group-0 FASTQ fixture. */
	private static String fastq(){
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<READS.length; i++){
			sb.append("@r").append(i).append('\n');
			sb.append(READS[i]).append('\n');
			sb.append("+\n");
			for(int j=0; j<READS[i].length(); j++){sb.append('I');}
			sb.append('\n');
		}
		return sb.toString();
	}

	/** Maps fixture IDs to their selected literal bases. */
	private static HashMap<String, String> expectedBases(){
		final HashMap<String, String> map=new HashMap<String, String>();
		for(int i=0; i<READS.length; i++){
			map.put("r"+i, READS[i]);
		}
		return map;
	}

	/** Builds one in-memory read for the normal partition oracle. */
	private static Read read(final String id, final String bases, final long numericID){
		return new Read(bases.getBytes(StandardCharsets.US_ASCII), quals(bases.length()),
				id, numericID);
	}

	/** Creates Q40 quality values matching the fixture FASTQ. */
	private static byte[] quals(final int len){
		final byte[] quals=new byte[len];
		for(int i=0; i<len; i++){quals[i]=40;}
		return quals;
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final class StaticState {

		StaticState() throws Exception{
			threads=Shared.threads();
			doHashAndSplit=KmerSort.doHashAndSplit;
			kcVerbose=KmerComparator.verbose;
			forceProcess=CLUMP_FORCE_PROCESS.getBoolean(null);
			fetchThreads=KMER_SORT3_FETCH_THREADS.getInt(null);
			defaultSeed=KMER_COMPARATOR_DEFAULT_SEED.getLong(null);
			defaultBorder=KMER_COMPARATOR_DEFAULT_BORDER.getInt(null);
			defaultHashes=KMER_COMPARATOR_DEFAULT_HASHES.getInt(null);
		}

		void setDeterministic() throws Exception{
			Shared.setThreads(1);
			KmerSort.doHashAndSplit=true;
			KmerComparator.verbose=false;
			CLUMP_FORCE_PROCESS.setBoolean(null, false);
			KMER_SORT3_FETCH_THREADS.setInt(null, 1);
			KMER_COMPARATOR_DEFAULT_SEED.setLong(null, 1);
			KMER_COMPARATOR_DEFAULT_BORDER.setInt(null, 1);
			KMER_COMPARATOR_DEFAULT_HASHES.setInt(null, 4);
		}

		void restore() throws Exception{
			Shared.setThreads(threads);
			KmerSort.doHashAndSplit=doHashAndSplit;
			KmerComparator.verbose=kcVerbose;
			CLUMP_FORCE_PROCESS.setBoolean(null, forceProcess);
			KMER_SORT3_FETCH_THREADS.setInt(null, fetchThreads);
			KMER_COMPARATOR_DEFAULT_SEED.setLong(null, defaultSeed);
			KMER_COMPARATOR_DEFAULT_BORDER.setInt(null, defaultBorder);
			KMER_COMPARATOR_DEFAULT_HASHES.setInt(null, defaultHashes);
		}

		private final int threads;
		private final boolean doHashAndSplit;
		private final boolean kcVerbose;
		private final boolean forceProcess;
		private final int fetchThreads;
		private final long defaultSeed;
		private final int defaultBorder;
		private final int defaultHashes;
	}

	private static Field field(final Class<?> type, final String name) throws Exception{
		final Field field=type.getDeclaredField(name);
		field.setAccessible(true);
		return field;
	}

	private static final int K=7;
	private static final int GROUPS=3;
	private static final String[] READS=new String[] {
			"AAAAAAAAAAAA",
			"AGCCAAAAAAAA",
			"AACAAAAAAAAA",
			"AGAAAAAAAAAA",
			"GTGTAAAAAAAA",
			"ACCAAAAAAAAA"
	};
	private static final int[] EXPECTED_GROUPS=new int[] {
			//The expected groups were generated from an independent enumerator:
			//repivot the literal with KmerComparator(k=7, seed=2, border=0,
			//hashes=4), read ReadKey.kmer, then bucket by hash(kmer)%3.
			2,
			0,
			1,
			2,
			0,
			1
	};

	private static final Field KEY_KMER;
	private static final Field CLUMP_FORCE_PROCESS;
	private static final Field KMER_SORT3_FETCH_THREADS;
	private static final Field KMER_COMPARATOR_DEFAULT_SEED;
	private static final Field KMER_COMPARATOR_DEFAULT_BORDER;
	private static final Field KMER_COMPARATOR_DEFAULT_HASHES;

	static{
		try{
			KEY_KMER=field(Class.forName("clump2.ReadKey"), "kmer");
			CLUMP_FORCE_PROCESS=field(Clump.class, "forceProcess");
			KMER_SORT3_FETCH_THREADS=field(KmerSort3.class, "fetchThreads");
			KMER_COMPARATOR_DEFAULT_SEED=field(KmerComparator.class, "defaultSeed");
			KMER_COMPARATOR_DEFAULT_BORDER=field(KmerComparator.class, "defaultBorder");
			KMER_COMPARATOR_DEFAULT_HASHES=field(KmerComparator.class, "defaultHashes");
		}catch(Exception e){
			throw new RuntimeException(e);
		}
	}

	private KmerSort3OverflowPartitionTest(){}
}
