package test.jgi;

import java.io.IOException;
import java.lang.invoke.MethodHandle;
import java.lang.invoke.MethodHandles;
import java.lang.invoke.MethodType;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import jgi.BBDuk;
import kmer.AbstractKmerTable;
import kmer.HashArray;
import stream.Read;

/** Synthetic timing of the unchanged matched primary/short edit-mode control.
 * One ordinary owner and one NEW shard2 loader reuse ACGA/Q40 with both edit
 * budgets1. An untimed actual addToMap(Read,0) call seeds seven literal keys.
 * Warmup and measurement then use the same primitive bound-handle loop; duplicate
 * insertions return zero. Complete length3/4 word vectors across all seven maps,
 * payloads, NEW state and cumulative counters are checked outside measurement.
 * No counters/tables are reset, no final/configuration/counter fields are seeded,
 * and no queues, loader start/run/fetch, processing, streams or reports are used.
 * Only normal constructor table allocators may start and join. This measures
 * cache-hot recursive generation and duplicate lookup on one representative
 * shard, not first insertion, reference-build scaling or whole-tool performance.
 * Bash owns fresh calibration and three matched JVM samples per source phase.
 * @author Sucrose
 */
public final class BBDukShortEditModeTimingTest{

	/** phase=baseline|fixed mode=pilot|measure warmup=POS calls=POS outdir=NEW. */
	public static void main(final String[] args) throws Throwable{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final int split=arg.indexOf('=');
			require(split>0 && split<arg.length()-1, "Use nonempty flag=value arguments: "+arg);
			require(options.put(arg.substring(0, split), arg.substring(split+1))==null, "Duplicate option: "+arg);
		}
		require(options.size()==5 && options.keySet().containsAll(Arrays.asList("phase", "mode", "warmup", "calls", "outdir")), "Use exactly phase, mode, warmup, calls and outdir.");
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbdukshorteditmodetimingtest.sh with assertions enabled.");
		final String phase=options.get("phase"), mode=options.get("mode");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require("pilot".equals(mode) || "measure".equals(mode), "Set mode=pilot|measure.");
		final long warmup=positive(options.get("warmup")), calls=positive(options.get("calls"));
		final long total=Math.addExact(1L, Math.addExact(warmup, calls));
		require(total<=MAX_TOTAL_CALLS, "Seed plus warmup and measured calls must total at most "+MAX_TOTAL_CALLS+" for bounded cumulative counters.");
		Math.multiplyExact(total, 4L);
		final ThreadMXBean clock=ManagementFactory.getThreadMXBean();
		require(clock.isCurrentThreadCpuTimeSupported() && clock.isThreadCpuTimeEnabled() && clock.getCurrentThreadCpuTime()>=0, "Require supported and already enabled main-thread CPU timing; do not substitute another metric.");
		final BBDukShortEditModeTimingTest test=new BBDukShortEditModeTimingTest(path(options.get("outdir")));
		test.run(clock, phase, mode, warmup, calls, total);
	}

	private BBDukShortEditModeTimingTest(final Path out_) throws Exception{
		out=out_; Files.createDirectory(out);
		write(out.resolve("settings.tsv"), String.join("\n", SETTINGS)+"\n");
		write(out.resolve("constructor-input.fq"), "@constructor_input\nACGA\n+\nIIII\n");
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(SETTINGS));
		args.add("literal=ACGA"); args.add("in="+out.resolve("constructor-input.fq"));
		write(out.resolve("constructor-args.txt"), String.join("\n", args)+"\n");
		tool=new BBDuk(args.toArray(new String[args.size()]));
		sets=(AbstractKmerTable[])field(BBDuk.class, "keySets").get(tool);
		require(sets.length==7, "Normal constructor must retain seven tables; one loader owns shard2.");
		require(field(BBDuk.class, "k").getInt(tool)==4 && field(BBDuk.class, "mink").getInt(tool)==3 && field(BBDuk.class, "editDistance").getInt(tool)==1 && field(BBDuk.class, "editDistance2").getInt(tool)==1 && field(BBDuk.class, "hammingDistance").getInt(tool)==1 && field(BBDuk.class, "hammingDistance2").getInt(tool)==1, "Normal constructor must select the unchanged matched edit1 control and promote each Hamming budget to1.");
		final Class<?> type=Class.forName("jgi.BBDuk$LoadThread");
		final Constructor<?> constructor=type.getDeclaredConstructor(BBDuk.class, int.class); constructor.setAccessible(true);
		final Method method=type.getDeclaredMethod("addToMap", Read.class, int.class); method.setAccessible(true);
		worker=(Thread)constructor.newInstance(tool, Integer.valueOf(2));
		bound=MethodHandles.lookup().unreflect(method).bindTo(worker);
		require(bound.type().equals(MethodType.methodType(long.class, Read.class, int.class)), "Bind the exact primitive stable addToMap(Read,int) entry; no per-call reflection arrays or boxing.");
		refReads=field(type, "refReadsT"); refBases=field(type, "refBasesT"); refKmers=field(type, "refKmersT"); added=field(type, "addedT"); stored=field(BBDuk.class, "storedKmers");
		expectedBases=bytes("ACGA"); expectedQuality=new byte[]{40, 40, 40, 40};
		read=new Read(expectedBases.clone(), expectedQuality.clone(), "measured_reference", 0); read.obj=Integer.valueOf(1);
		write(out.resolve("inputs.tsv"), "id\tnumeric_id\tscaffold_id\tbases\tqualities\nmeasured_reference\t0\t1\tACGA\tIIII\n");
		final StringBuilder keys=new StringBuilder("word\tkey\tshard\tscaffold_id\n");
		for(int i=0; i<KEYS.length; i++){
			require(key(WORDS[i])==KEYS[i] && KEYS[i]%7==2, "Independently packed literal word/key pair must belong to shard2: "+WORDS[i]);
			keys.append(WORDS[i]).append('\t').append(KEYS[i]).append("\t2\t1\n");
		}
		write(out.resolve("literal-keys.tsv"), keys.toString());
		write(out.resolve("reflection.tsv"), "operation\tmember\nREAD_ONLY\tBBDuk.keySets,storedKmers,k,mink,editDistance,editDistance2,hammingDistance,hammingDistance2\nCONSTRUCT_ONLY\tBBDuk.LoadThread(shard2)\nUNREFLECT_BIND_ONCE\tLoadThread.addToMap(Read,int)\nOBSERVE\tLoadThread.refReadsT,refBasesT,refKmersT,addedT\nFIELD_WRITES\tNONE\n");
		write(out.resolve("controls.tsv"), "configuration\tprimary_edit1_short_edit1\nreference\tACGA\nshard\t2\nseed_return\t7\nduplicate_return\t0\nexpected_keys\t7\nloader_starts\t0\n");
		write(out.resolve("observations.tsv"), "stage\tcompleted_calls\tref_reads\tref_bases\tref_kmers\tadded_field\towner_stored_kmers\tloader_state\tmap_key_counts\n");
	}

	/** The same primitive loop performs warmup and measurement; all observations surround it. */
	private void run(final ThreadMXBean clock, final String phase, final String mode, final long warmup, final long calls, final long total) throws Throwable{
		final long main=Thread.currentThread().getId(), warmupTotal=Math.addExact(1L, warmup);
		write(out.resolve("config.tsv"), "version\t1\nphase\t"+phase+"\nmode\t"+mode+"\nwarmup\t"+warmup+"\ncalls\t"+calls+"\nreference_length\t4\nk\t4\nmink\t3\nfull_windows_per_call\t1\nshort_seeds_per_call\t2\nshard\t2\nexpected_keys\t7\nmax_total_calls\t"+MAX_TOTAL_CALLS+"\n");
		write(out.resolve("clock.tsv"), "metric\tcurrent_thread_cpu_nanoseconds\nmain_thread_id\t"+main+"\nsupported\ttrue\nenabled\ttrue\n");
		observe("before_seed", 0, false);
		final long seed=(long)bound.invokeExact(read, 0);
		require(seed==7, "Untimed real addToMap must create exactly the seven independently listed shard2 keys; return="+seed);
		observe("after_seed", 1, true);
		write(out.resolve("seed.tsv"), "returned_new_keys\t"+seed+"\nexpected_keys\t7\nuntimed\ttrue\n");
		final long warmupChecksum=loop(warmup);
		require(warmupChecksum==0, "After the untimed seed, warmup may insert no additional keys; checksum="+warmupChecksum);
		observe("after_warmup", warmupTotal, true);
		final long cpu0=clock.getCurrentThreadCpuTime(), wall0=System.nanoTime();
		final long checksum=loop(calls);
		final long wall=System.nanoTime()-wall0, cpu1=clock.getCurrentThreadCpuTime(), cpu=cpu1-cpu0;
		require(Thread.currentThread().getId()==main && cpu0>=0 && cpu1>cpu0 && wall>0, "Complete the loop on the original main thread with positive CPU/wall deltas.");
		require(checksum==0, "Every measured call must revisit the already stored literal keys; checksum="+checksum);
		observe("after_measure", total, true);
		write(out.resolve("totals.tsv"), "seed_new_keys\t"+seed+"\nwarmup_checksum\t"+warmupChecksum+"\nmeasured_checksum\t"+checksum+"\ntotal_reference_calls\t"+total+"\nexpected_ref_reads\t"+total+"\nexpected_ref_bases\t"+Math.multiplyExact(total, 4L)+"\nexpected_ref_kmers\t"+total+"\nexpected_map_keys\t7\n");
		require("pilot".equals(mode) || cpu>=10000000000L, "Measured loop below ten main-thread CPU seconds; enlarge the common count and redo BOTH matched source cohorts.");
		final String summary="BBDUK_SHORT_EDIT_MODE_TIMING_PASS version=1 phase="+phase+" mode="+mode+" warmup="+warmup+" calls="+calls+" checksum="+checksum+" cpuNanos="+cpu+" wallNanos="+wall+" seedKeys=7 owners=1 loaders=1 loaderStarts=0 failures=0\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
	}

	/** No per-call allocation, counter/table resets, reflection argument arrays, checks or output. */
	private long loop(final long calls) throws Throwable{
		assert(calls>0 && calls<=MAX_TOTAL_CALLS) : "Bounded positive calls protect cumulative reference-counter arithmetic; calls="+calls;
		long checksum=0;
		for(long i=0; i<calls; i++){checksum+=(long)bound.invokeExact(read, 0);}
		return checksum;
	}
	private void observe(final String stage, final long completed, final boolean seeded) throws Exception{
		require(worker.getState()==Thread.State.NEW && !worker.isAlive(), "Never start/run/fetch or queue the single shard2 loader.");
		require(refReads.getLong(worker)==completed && refBases.getLong(worker)==Math.multiplyExact(completed, 4L) && refKmers.getLong(worker)==completed && added.getLong(worker)==0, "Direct addToMap must retain cumulative N/4N/N counters and leave run()-owned addedT zero at "+stage);
		require(stored.getLong(tool)==0, "No owner aggregation may run during direct-loader timing.");
		require("measured_reference".equals(read.id) && read.numericID==0 && Integer.valueOf(1).equals(read.obj) && read.mate==null && Arrays.equals(read.bases, expectedBases) && Arrays.equals(read.quality, expectedQuality), "Preserve independent ACGA/Q40 payload,identity,assigned scaffold1 andnullmate at "+stage);
		final long[] totals=new long[7];
		for(int shard=0; shard<7; shard++){
			require(sets[shard] instanceof HashArray, "Normal array=t must expose primary and victim key counts.");
			final HashArray table=(HashArray)sets[shard]; totals[shard]=table.size()+table.victims().size();
			require(totals[shard]==(seeded && shard==2 ? 7 : 0), "Only shard2 may contain the seven literal keys; stage="+stage+" shard="+shard+" keys="+totals[shard]);
		}
		final StringBuilder vectors=new StringBuilder("word\tkey\texpected_present\tshard_values\n");
		for(int length=3; length<=4; length++){
			for(int value=0; value<(1<<(2*length)); value++){
				final long key=(1L<<(2*length))|value; final boolean present=seeded && contains(key);
				final int[] actual=new int[7], expected=new int[7]; Arrays.fill(expected, -1); if(present){expected[2]=1;}
				for(int shard=0; shard<7; shard++){actual[shard]=sets[shard].getValue(key);}
				require(Arrays.equals(actual, expected), "Complete literal key/shard vector mismatch at "+stage+" key="+key+" actual="+Arrays.toString(actual));
				vectors.append(word(value, length)).append('\t').append(key).append('\t').append(present).append('\t').append(Arrays.toString(actual)).append('\n');
			}
		}
		write(out.resolve(stage+".keys.tsv"), vectors.toString());
		append(out.resolve("observations.tsv"), stage+"\t"+completed+"\t"+refReads.getLong(worker)+"\t"+refBases.getLong(worker)+"\t"+refKmers.getLong(worker)+"\t"+added.getLong(worker)+"\t"+stored.getLong(tool)+"\t"+worker.getState()+"\t"+Arrays.toString(totals)+"\n");
	}
	private static boolean contains(final long key){for(final long literal : KEYS){if(key==literal){return true;}}return false;}
	/** Independent literal DNA packing; never call production key/mask/recursion helpers for an oracle. */
	private static long key(final String word){long value=0; for(int i=0; i<word.length(); i++){final int base="ACGT".indexOf(word.charAt(i)); require(base>=0, "Use literal DNA oracle words."); value=(value<<2)|base;} return value|(1L<<(2*word.length()));}
	private static String word(int value, final int length){final char[] result=new char[length]; for(int i=length-1; i>=0; i--){result[i]="ACGT".charAt(value&3); value>>>=2;} return new String(result);}
	private static Field field(final Class<?> type, final String name) throws NoSuchFieldException{final Field result=type.getDeclaredField(name); result.setAccessible(true); return result;}
	private static long positive(final String text){require(text!=null && text.matches("[0-9]+"), "Use positive integer call counts: "+text); final long value=Long.parseLong(text); require(value>0, "Call counts must be positive."); return value;}
	private static Path path(final String text){require(text!=null && !text.isEmpty() && text.indexOf('=')<0 && text.indexOf('\t')<0 && text.indexOf('\n')<0, "Use a nonempty path without equals,tabs,newlines for the legacy constructor."); return Paths.get(text).toAbsolutePath().normalize();}
	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static void write(final Path file, final String text) throws IOException{Files.write(file, bytes(text), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path file, final String text) throws IOException{Files.write(file, bytes(text), StandardOpenOption.APPEND);}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}

	private final Path out;
	private final BBDuk tool;
	private final AbstractKmerTable[] sets;
	private final Thread worker;
	private final MethodHandle bound;
	private final Read read;
	private final byte[] expectedBases, expectedQuality;
	private final Field refReads, refBases, refKmers, added, stored;
	private static final long MAX_TOTAL_CALLS=1000000000L;
	private static final String[] WORDS={"AAC", "AGA", "CCG", "ACGG", "AGGA", "ATCG", "GCGA"};
	private static final long[] KEYS={65L, 72L, 86L, 282L, 296L, 310L, 408L};
	private static final String[] SETTINGS={"k=4", "mink=3", "edist=1", "edist2=1", "hdist=0", "hdist2=0", "ktrim=r", "maskmiddle=f", "rcomp=f", "speed=0", "rskip=1", "qskip=1", "qhdist=0", "qhdist2=0", "findbestmatch=f", "rename=f", "skipr1=f", "skipr2=f", "restrictleft=0", "restrictright=0", "minlen=1", "t=1", "simd=f", "qin=33", "qout=33", "array=t", "prealloc=f", "ow=f", "app=f"};
}
