package test.bloom;

import java.io.ObjectOutputStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

import bloom.BloomFilter;
import bloom.BloomFilterWrapper;
import bloom.KCountArray7MTA;
import shared.Shared;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;

/** Checks that loaded filter geometry also governs worker scratch and count diagnostics.
 * Uses an actual serialized counting filter, then invokes the real worker helpers
 * without starting a worker or output stream. Does not test pipeline scheduling.
 * Run in a fresh JVM because wrapper construction changes process-wide settings.
 * @author Brian Bushnell, Shinobu */
public final class BloomSerializedFilterTest {

	/** Takes one new scratch-directory path and checks four CLI/snapshot combinations. */
	public static void main(final String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Expected one new scratch-directory path.");}
		final Path dir=Files.createDirectory(Paths.get(args[0]));
		final Path snapshot=dir.resolve("filter.ser"), input=dir.resolve("input.fq");
		Files.write(input, "@query\nCCCCCC\n+\nIIIIII\n".getBytes(StandardCharsets.US_ASCII));
		Shared.setThreads(1);
		BloomFilter.OVERRIDE_CELLS=10007;
		KCountArray7MTA.setSeed(1);
		final BloomFilter filter=new BloomFilter(null, null, null, 3, 3, 4, 1, 2, true, false, false, 0.01f, false);
		filter.filter.increment(42L, 5);//Canonical CCC is GGG, literal packed key42.
		check(filter.getCount(42L)==5 && filter.getCount(63L)==0,
				"The literal table must support GGG and reject canonical AAA/TTT before serialization.");
		try(ObjectOutputStream out=new ObjectOutputStream(Files.newOutputStream(snapshot))){out.writeObject(filter);}
		final byte[] saved=Files.readAllBytes(snapshot);
		int failures=0;
		failures+=exercise(input, snapshot, "matching-options", 3, true, 2);
		failures+=exercise(input, snapshot, "different-minhits", 3, true, 1);
		failures+=exercise(input, snapshot, "different-k", 4, true, 2);
		failures+=exercise(input, snapshot, "different-strand", 3, false, 2);
		check(Arrays.equals(saved, Files.readAllBytes(snapshot)), "Loading/querying must leave snapshot bytes unchanged.");
		if(failures!=0){throw new AssertionError("Loaded filter contracts failed "+failures+" of 4 cases.");}
		System.out.println("BLOOM_SERIALIZED_FILTER_TEST_PASS cases=4");
	}

	/** Compares worker matches and every diagnostic record against literal snapshot expectations. */
	private static int exercise(final Path input, final Path snapshot, final String name,
			final int k, final boolean rcomp, final int minhits) throws Exception{
		final BloomFilterWrapper wrapper=new BloomFilterWrapper(new String[]{"in="+input, "serialin="+snapshot,
				"k="+k, "rcomp="+rcomp, "minhits="+minhits, "mincount=1", "t=1", "interleaved=f", "pigz=f", "unpigz=f"});
		final Class<?> workerClass=Class.forName("bloom.BloomFilterWrapper$ProcessThread");
		final Constructor<?> constructor=workerClass.getDeclaredConstructor(BloomFilterWrapper.class,
				ConcurrentReadInputStream.class, ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class,
				ConcurrentReadOutputStream.class, int.class);
		constructor.setAccessible(true);
		final Object worker=constructor.newInstance(wrapper, null, null, null, null, 0);
		final Method matches=workerClass.getDeclaredMethod("matchesFilter", Read.class, Read.class);
		final Method addCounts=workerClass.getDeclaredMethod("addCounts", Read.class, ArrayList.class);
		matches.setAccessible(true);
		addCounts.setAccessible(true);
		final Read supported=new Read("CCCCCC".getBytes(StandardCharsets.US_ASCII), null, "query", 17, false);
		final Read absent=new Read("AAAAAA".getBytes(StandardCharsets.US_ASCII), null, "absent", 18, false);
		final ArrayList<Read> counts=new ArrayList<Read>();
		boolean matched=false, rejected=false;
		String error="none";
		try{
			matched=(Boolean)matches.invoke(worker, supported, null);
			rejected=!(Boolean)matches.invoke(worker, absent, null);
			addCounts.invoke(worker, supported, counts);
		}catch(InvocationTargetException e){error=e.getCause().getClass().getSimpleName();}
		boolean pass=error.equals("none") && matched && rejected && counts.size()==4;
		for(Read read : counts){
			pass&=Arrays.equals(read.bases, GGG) && read.quality==null && "5".equals(read.id) && read.numericID==17;
		}
		System.out.println("case="+name+" match="+matched+" rejectAbsent="+rejected+
				" expectedWindows=4 actualWindows="+counts.size()+" exception="+error+" pass="+pass);
		return pass ? 0 : 1;
	}

	/** Rejects an invalid fixture rather than interpreting its output as a product failure. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final byte[] GGG={'G', 'G', 'G'};
}
