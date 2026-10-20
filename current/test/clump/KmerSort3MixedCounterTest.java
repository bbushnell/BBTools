package test.clump;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import clump.Clump;
import clump.KmerComparator;
import clump.KmerSort;
import clump.KmerSort3;
import shared.Shared;

/** Mixed normal/overflow KmerSort3 counter oracle with exact payload checks.
 * @author Jean
 */
public final class KmerSort3MixedCounterTest{

	/** Runs all-normal and mixed normal/overflow counter cases. */
	public static void main(String[] args) throws Exception{
		check(args.length==1, "Expected a fresh fixture directory");
		final Path root=Paths.get(args[0]);
		Files.createDirectory(root);

		final StaticState state=new StaticState();
		try{
			state.setDeterministic();
			final CaseData data=new CaseData();
			System.out.println("case\treads_in\tbases_in\treads_out\tbases_out\tpayloads");
			run(root.resolve("normal"), data, true, true);
			run(root.resolve("mixed"), data, false, true);
			run(root.resolve("mixed-null"), data, false, false);
		}finally{
			state.restore();
		}
	}

	/** Runs KmerSort3 on one fresh directory. */
	private static void run(Path root, CaseData data, boolean forceProcess, boolean write)
			throws Exception{
		Files.createDirectory(root);
		for(int group=0; group<GROUPS; group++){
			write(root.resolve("in_"+group+".fq"), data.input[group]);
		}

		final long fileMem=Shared.memAvailable(1);
		KmerSort3.main(fileMem, 1, 2, new String[] {
				"in="+root.resolve("in_%.fq"),
				write ? "out="+root.resolve("out_%.fq") : "out=null",
				"groups="+GROUPS,
				"k=7",
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
				"fetchthreads=2",
				"verbose=f"
		});

		check(KmerSort.lastReadsIn==data.reads,
				"Wrong readsIn for "+root+": "+KmerSort.lastReadsIn);
		check(KmerSort.lastBasesIn==data.bases,
				"Wrong basesIn for "+root+": "+KmerSort.lastBasesIn);
		check(KmerSort.lastReadsOut==data.reads,
				"Wrong output-eligible reads for "+root+": "+KmerSort.lastReadsOut);
		check(KmerSort.lastBasesOut==data.bases,
				"Wrong output-eligible bases for "+root+": "+KmerSort.lastBasesOut);

		final String payloadState;
		if(write){
			final LinkedHashMap<String, String> observed=new LinkedHashMap<String, String>();
			for(int group=0; group<GROUPS; group++){
				for(Map.Entry<String, String> e : read(root.resolve("out_"+group+".fq")).entrySet()){
					final String prior=observed.put(e.getKey(), e.getValue());
					check(prior==null, "Duplicate output read across groups: "+e.getKey());
				}
			}
			check(observed.equals(data.expected),
					"Output payload mismatch for "+root+": "+observed);
			payloadState="exact";
		}else{
			payloadState="none";
		}
		System.out.println(root.getFileName().toString()+'\t'+KmerSort.lastReadsIn+'\t'+
				KmerSort.lastBasesIn+'\t'+KmerSort.lastReadsOut+'\t'+
				KmerSort.lastBasesOut+'\t'+payloadState);
	}

	/** Writes one complete input shard. */
	private static void write(Path path, String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII),
				StandardOpenOption.CREATE_NEW);
	}

	/** Strict four-line FASTQ parser keyed by read ID. */
	private static LinkedHashMap<String, String> read(final Path path) throws Exception{
		check(Files.exists(path), "Missing output: "+path);
		final List<String> lines=Files.readAllLines(path, StandardCharsets.US_ASCII);
		check(lines.size()%4==0, "Malformed FASTQ line count in "+path);
		final LinkedHashMap<String, String> map=new LinkedHashMap<String, String>();
		for(int i=0; i<lines.size(); i+=4){
			check(lines.get(i).startsWith("@"), "Malformed header in "+path+": "+lines.get(i));
			check(lines.get(i+2).equals("+"), "Malformed separator in "+path+": "+lines.get(i+2));
			check(lines.get(i+1).length()==lines.get(i+3).length(),
					"Sequence/quality length mismatch in "+path);
			final String id=lines.get(i).substring(1).split("\\s+", 2)[0];
			final String prior=map.put(id, lines.get(i+1)+'\n'+lines.get(i+3));
			check(prior==null, "Duplicate output read "+id+" in "+path);
		}
		return map;
	}

	/** Lightweight assertion independent of JVM assertion flags. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Generates one bounded skewed four-shard fixture. */
	private static final class CaseData{

		CaseData(){
			int id=0;
			for(int group=0; group<GROUPS; group++){
				final int limit=(group==0 ? BIG_GROUP_READS : SMALL_GROUP_READS);
				final StringBuilder sb=new StringBuilder(limit*120);
				for(int i=0; i<limit; i++, id++){
					final String name="r"+id;
					final String bases=bases(group, i);
					final String quals=quals(bases.length());
					sb.append('@').append(name).append('\n');
					sb.append(bases).append('\n');
					sb.append("+\n");
					sb.append(quals).append('\n');
					expected.put(name, bases+'\n'+quals);
					reads++;
					this.bases+=bases.length();
				}
				input[group]=sb.toString();
			}
		}

		private String bases(final int group, final int read){
			final StringBuilder sb=new StringBuilder(READ_LEN);
			for(int i=0; i<READ_LEN; i++){
				sb.append(ALPHABET[(group*3+read*5+i*7)&3]);
			}
			return sb.toString();
		}

		private String quals(final int len){
			final StringBuilder sb=new StringBuilder(len);
			for(int i=0; i<len; i++){sb.append((char)('I'-(i&3)));}
			return sb.toString();
		}

		final String[] input=new String[GROUPS];
		final LinkedHashMap<String, String> expected=new LinkedHashMap<String, String>();
		long reads, bases;
	}

	/** Saves and restores the static state touched by in-JVM KmerSort3 calls. */
	private static final class StaticState{

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
			Shared.setThreads(2);
			KmerSort.doHashAndSplit=true;
			KmerComparator.verbose=false;
			CLUMP_FORCE_PROCESS.setBoolean(null, false);
			KMER_SORT3_FETCH_THREADS.setInt(null, 2);
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

	private static final int GROUPS=4;
	private static final int BIG_GROUP_READS=48;
	private static final int SMALL_GROUP_READS=12;
	private static final int READ_LEN=40;
	private static final char[] ALPHABET=new char[] {'A', 'C', 'G', 'T'};

	private static final Field CLUMP_FORCE_PROCESS;
	private static final Field KMER_SORT3_FETCH_THREADS;
	private static final Field KMER_COMPARATOR_DEFAULT_SEED;
	private static final Field KMER_COMPARATOR_DEFAULT_BORDER;
	private static final Field KMER_COMPARATOR_DEFAULT_HASHES;

	static{
		try{
			CLUMP_FORCE_PROCESS=field(Clump.class, "forceProcess");
			KMER_SORT3_FETCH_THREADS=field(KmerSort3.class, "fetchThreads");
			KMER_COMPARATOR_DEFAULT_SEED=field(KmerComparator.class, "defaultSeed");
			KMER_COMPARATOR_DEFAULT_BORDER=field(KmerComparator.class, "defaultBorder");
			KMER_COMPARATOR_DEFAULT_HASHES=field(KmerComparator.class, "defaultHashes");
		}catch(Exception e){
			throw new RuntimeException(e);
		}
	}

	private KmerSort3MixedCounterTest(){}
}
