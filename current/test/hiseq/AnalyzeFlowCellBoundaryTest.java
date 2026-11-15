package test.hiseq;

import java.io.File;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import bloom.BloomFilter;
import hiseq.AnalyzeFlowCell;
import hiseq.FlowCell;
import shared.Timer;
import stream.Read;
import stream.SamLine;
import stream.Streamer;
import structures.AtomicStringNum;
import structures.IntList;
import structures.ListNum;

/**
 * Finite observations of reused depth buffers, projection bounds, output aliases
 * and first-observation barcode insertion. Reflection exercises existing private
 * paths without adding a production-only testing interface.
 * @author Brian Bushnell, Shinobu
 */
public class AnalyzeFlowCellBoundaryTest{

	/** Runs baseline observations or fixed-contract assertions in one dedicated JVM. */
	public static void main(String[] args) throws Exception{
		File out=null;
		boolean fixed=false;
		for(String arg : args){
			if(arg.startsWith("outdir=")){out=new File(arg.substring(7));}
			else if(arg.equals("phase=fixed")){fixed=true;}
			else if(!arg.equals("phase=baseline")){throw new IllegalArgumentException(arg);}
		}
		require(out!=null && out.mkdir(), "Use a fresh fixture directory");
		final File input=new File(out, "input.fq");
		final byte[] original=("@INST:1:FC:1:1101:100:100 1:N:0:ACGT+TGCA\nACGT\n+\nIIII\n")
			.getBytes(StandardCharsets.US_ASCII);
		Files.write(input.toPath(), original);
		System.out.println("case\tobserved\texpected\tresult");
		final AnalyzeFlowCell tool=create(input, "idmaskread=0");
		final Class<?>[] depthTypes={Read.class, IntList.class};
		for(String name : new String[] {"short", "null"}){
			final IntList list=new IntList(4);
			list.add(7); list.add(8);
			call(tool, "fillDepthList", depthTypes, name.equals("null") ? null : read(1, 0), list);
			check(name+"-buffer", list.size, fixed ? 0 : 2);
		}
		final long previousCells=BloomFilter.OVERRIDE_CELLS;
		final BloomFilter filter;
		try{
			BloomFilter.OVERRIDE_CELLS=4096;
			filter=new BloomFilter(31, 31, 2, 1, 1, true, false, false, 0.001f);
		}finally{BloomFilter.OVERRIDE_CELLS=previousCells;}
		set(tool, "bloomFilter", filter);
		try{
			for(int length : new int[] {31,32}){
				final IntList list=new IntList(4);
				call(tool, "fillDepthList", depthTypes, read(length, 0), list);
				check("all-kmers-"+length, list.size, length-31+(fixed ? 1 : 0));
				for(int i=0; i<list.size; i++){require(list.get(i)==0, "Empty Bloom filter must have zero depth");}
			}
			set(tool, "idmask_read", 7);
			for(long id : new long[] {0,6}){
				final IntList list=new IntList(4);
				call(tool, "fillDepthList", depthTypes, read(31, id), list);
				check("sampled-first-id-"+id, list.size, fixed && id==6 ? 1 : 0);
			}
		}finally{filter.filter.shutdown();}
		final Class<?>[] projectionTypes={IntList.class, IntList.class};
		final IntList kmers=new IntList(2), bases=new IntList(40);
		call(tool, "deblur", projectionTypes, kmers, bases);
		check("empty-projection", bases.size, fixed ? 0 : 31);
		kmers.add(2); kmers.add(1);
		call(tool, "deblur", projectionTypes, kmers, bases);
		check("two-kmer-projection", bases.size, fixed ? 32 : 33);
		for(int i=0; i<31; i++){require(bases.get(i)==2, "First kmer covers positions0..30");}
		require(bases.get(31)==1, "Second kmer alone covers position31");
		for(String flag : new String[] {"coords", "alignout", "out"}){
			boolean rejected=false;
			try{create(input, flag+"="+input.getPath());}
			catch(RuntimeException e){rejected=true;}
			check(flag+"-input-alias", rejected ? 1 : 0, fixed || flag.equals("out") ? 1 : 0);
			require(Arrays.equals(original, Files.readAllBytes(input.toPath())), "Constructor changed its input file");
		}
		final File equalsPath=new File(out, "input=literal.fq");
		Files.write(equalsPath.toPath(), original);
		boolean accepted=false;
		try{create(equalsPath, "align=f"); accepted=true;}
		catch(RuntimeException e){accepted=false;}
		check("equals-in-input-path", accepted ? 1 : 0, fixed ? 1 : 0);
		barcodeRace(input, fixed);
		barcodeWriteFailure(fixed);
		emptyBatch(input, fixed);
		sampledProjection(input, fixed);
		extraReferences(out, fixed);
		emptySummary(out, fixed);
	}

	/** Empty filtering input must report finite discard percentages and emit no records. */
	private static void emptySummary(File out, boolean fixed) throws Exception{
		final File input=new File(out, "empty-report.fq"), output=new File(out, "empty-output.fq");
		Files.write(input.toPath(), new byte[0]);
		final AnalyzeFlowCell tool=create(input, "out="+output, "align=f");
		set(tool, "flowcell", new FlowCell(31));
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		final PrintStream report=new PrintStream(bytes, true, "UTF-8");
		set(tool, "outstream", report);
		try{call(tool, "processReads", new Class<?>[] {Timer.class}, new Timer());}
		finally{report.close();}
		final String text=new String(bytes.toByteArray(), StandardCharsets.UTF_8);
		System.out.println("OBSERVED-empty-summary\n"+text);
		require(!tool.errorState && tool.readsProcessed==0 && tool.basesProcessed==0
			&& tool.readsDiscarded==0 && tool.basesDiscarded==0,
			"Empty filtering input must complete without processing records");
		require(output.isFile() && output.length()==0, "Empty input must emit an empty output file");
		int lines=0, nonfinite=0;
		for(String line : text.split("\n")){
			if(line.startsWith("Reads Discarded:") || line.startsWith("Bases Discarded:")){
				lines++;
				if(line.contains("NaN%")){nonfinite++;}
				if(fixed){require(line.trim().endsWith("\t0.000%"), "Zero observations must report zero percent: "+line);}
			}
		}
		require(lines==2, "Both discard summary lines must be observed");
		check("empty-summary-nonfinite-percentages", nonfinite, fixed ? 0 : 2);
	}

	/** Four sampled windows cover bases2..56 of a64-base read, not adjacent starts0..3. */
	private static void sampledProjection(File input, boolean fixed) throws Exception{
		final AnalyzeFlowCell tool=create(input, "idmaskread=7");
		set(tool, "deblurDepths", true);
		set(tool, "smoothDepths", 0);
		final long previousCells=BloomFilter.OVERRIDE_CELLS;
		final BloomFilter filter;
		try{
			BloomFilter.OVERRIDE_CELLS=4096;
			filter=new BloomFilter(31, 31, 2, 1, 1, true, false, false, 0.001f);
		}finally{BloomFilter.OVERRIDE_CELLS=previousCells;}
		try{
			//A31 and its reverse-complement T31 canonicalize to the larger62-bit all-ones key.
			final long key=(1L<<62)-1;
			filter.filter.increment(key);
			require(filter.getCount(key)==1, "One inserted homopolymer must have depth1");
			set(tool, "bloomFilter", filter);
			final IntList kmers=new IntList(8), bases=new IntList(64);
			call(tool, "processTileKmers", new Class<?>[] {Read.class,IntList.class,IntList.class},
				read(64, 0), kmers, bases);
			check("sampled-window-count", kmers.size, 4);
			check("sampled-projection-length", bases.size, fixed ? 64 : 35);
			for(int i=0; i<kmers.size; i++){require(kmers.get(i)==1, "All four sampled windows are A31");}
			for(int i=0; i<bases.size; i++){
				final int expected=fixed ? (i>=2 && i<=56 ? 1 : 0) : (i<34 ? 1 : 0);
				require(bases.get(i)==expected, "Projection at base"+i+": observed="+bases.get(i)+", expected="+expected);
			}
		}finally{filter.filter.shutdown();}
	}

	/** Two extra references must both load even when the primary reads limit is zero. */
	private static void extraReferences(File out, boolean fixed) throws Exception{
		final File primary=new File(out, "empty-primary.fq");
		Files.write(primary.toPath(), new byte[0]);
		final File extra1=new File(out, "extra-one.fq"), extra2=new File(out, "extra-two.fq");
		final char[] bases=new char[31], qualities=new char[31];
		Arrays.fill(bases, 'C'); Arrays.fill(qualities, 'I');
		final byte[] record=("@extra\n"+new String(bases)+"\n+\n"+new String(qualities)+"\n")
			.getBytes(StandardCharsets.US_ASCII);
		Files.write(extra1.toPath(), record); Files.write(extra2.toPath(), record);
		final AnalyzeFlowCell tool=create(primary, "loadthreads=1", "idmaskwrite=0", "reads=0",
			"extra="+extra1+","+extra2);
		final long previousCells=BloomFilter.OVERRIDE_CELLS;
		try{
			BloomFilter.OVERRIDE_CELLS=4096;
			call(tool, "loadKmers", new Class<?>[0]);
			final Field field=AnalyzeFlowCell.class.getDeclaredField("bloomFilter");
			field.setAccessible(true);
			final BloomFilter filter=(BloomFilter)field.get(tool);
			//C31 reverse-complements to G31, whose two-bit code is repeated10.
			check("extra-reference-count", filter.getCount(0x2aaaaaaaaaaaaaaaL), fixed ? 2 : 0);
			require(!tool.errorState, "Valid empty primary and two literal extras must not report an input error");
		}finally{BloomFilter.OVERRIDE_CELLS=previousCells;}
	}

	/** An empty nonterminal batch must not hide the subsequent valid read. */
	private static void emptyBatch(File input, boolean fixed) throws Exception{
		final AnalyzeFlowCell tool=create(input, "align=f");
		set(tool, "flowcell", new FlowCell(31));
		set(tool, "fillThreads", 1);
		final EmptyThenRead stream=new EmptyThenRead();
		stream.start();
		try{tool.fillTilesInner(stream, null);}
		finally{stream.close();}
		require(!tool.errorState, "The finite stream and valid read must not report failure");
		check("empty-batch-read-count", tool.readsProcessed, fixed ? 1 : 0);
		check("empty-batch-base-count", tool.basesProcessed, fixed ? 40 : 0);
	}

	/** Single-consumer fixture implementing the documented empty-batch contract. */
	private static class EmptyThenRead implements Streamer{
		@Override
		public String fname(){return "finite-empty-then-read";}
		@Override
		public void start(){}
		@Override
		public void close(){closed=true;}
		@Override
		public boolean paired(){return false;}
		@Override
		public int pairnum(){return 0;}
		@Override
		public long readsProcessed(){return next>1 ? 1 : 0;}
		@Override
		public long basesProcessed(){return readsProcessed()*40;}
		@Override
		public void setSampleRate(float rate, long seed){throw new UnsupportedOperationException();}
		@Override
		public ListNum<Read> nextList(){
			if(!hasMore()){return null;}
			final ArrayList<Read> reads=new ArrayList<Read>();
			if(next==1){reads.add(read(40, 0));}
			return new ListNum<Read>(reads, next++);
		}
		@Override
		public ListNum<SamLine> nextLines(){throw new UnsupportedOperationException();}
		@Override
		public boolean hasMore(){return !closed && next<2;}
		@Override
		public boolean errorState(){return false;}
		private int next;
		private boolean closed;
	}

	/** Creates the tool without starting input, output, alignment or worker streams. */
	private static AnalyzeFlowCell create(File input, String... extra){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+input,
			"loadkmers=f", "threads=1", "simd=f", "overwrite=t", "qin=33"));
		args.addAll(Arrays.asList(extra));
		return new AnalyzeFlowCell(args.toArray(new String[args.size()]));
	}

	/** Makes a read with literal valid bases, qualities and coordinates. */
	private static Read read(int length, long id){
		final byte[] bases=new byte[length], qualities=new byte[length];
		Arrays.fill(bases, (byte)'A'); Arrays.fill(qualities, (byte)40);
		return new Read(bases, qualities, "INST:1:FC:1:1101:100:100 1:N:0:ACGT+TGCA", id);
	}

	/** Forces both workers to observe absence before either inserts the shared code. */
	private static void barcodeRace(File input, boolean fixed) throws Exception{
		final AnalyzeFlowCell tool=create(input, "align=f");
		final FlowCell flowcell=new FlowCell(31);
		set(tool, "flowcell", flowcell);
		set(tool, "barcodeStats", call(tool, "loadBarcodes", new Class<?>[] {String.class}, (Object)null));
		final FirstLookupMap counts=new FirstLookupMap();
		set(tool, "barcodeMap", counts);
		final Class<?> workerClass=Class.forName("hiseq.AnalyzeFlowCell$ProcessThread");
		final Constructor<?> constructor=workerClass.getDeclaredConstructor(AnalyzeFlowCell.class,
			Streamer.class, Streamer.class, int.class, FlowCell.class);
		constructor.setAccessible(true);
		final Throwable[] failures=new Throwable[2];
		final Thread[] runners=new Thread[2];
		for(int i=0; i<2; i++){
			final int index=i;
			final Object worker=constructor.newInstance(tool, null, null, i, flowcell);
			runners[i]=new Thread(new Runnable(){
				@Override
				public void run(){
					try{call(worker, "processReadPair", new Class<?>[] {Read.class,Read.class}, read(40,index), null);}
					catch(Throwable t){failures[index]=t;}
				}
			}, "barcode-first-lookup-"+i);
			runners[i].setDaemon(true);
		}
		for(Thread runner : runners){runner.start();}
		for(Thread runner : runners){runner.join(15000); require(!runner.isAlive(), "Barcode fixture did not complete");}
		for(Throwable failure : failures){if(failure!=null){throw new AssertionError("Barcode worker failed", failure);}}
		require(counts.size()==1, "Both literal headers carry the same barcode");
		check("barcode-first-observation", counts.values().iterator().next().n.get(), fixed ? 2 : 1);
	}

	/** Verifies that the real barcode writer's failure reaches its caller. */
	private static void barcodeWriteFailure(boolean fixed) throws Exception{
		require(new File("/dev/full").exists(), "This fixture requires Linux /dev/full");
		final ArrayList<AtomicStringNum> counts=new ArrayList<AtomicStringNum>();
		counts.add(new AtomicStringNum("ACGT+TGCA", 2));
		boolean threw=false;
		try{call(null, "dumpBarcodes", new Class<?>[] {Collection.class,String.class,boolean.class}, counts, "/dev/full", true);}
		catch(RuntimeException e){threw=true;}
		check("barcode-write-error", threw ? 1 : 0, fixed ? 1 : 0);
	}

	/** Returns the pre-barrier lookup value to force the valid first-insertion race. */
	private static class FirstLookupMap extends ConcurrentHashMap<String,AtomicStringNum>{
		@Override
		public AtomicStringNum get(Object key){
			final AtomicStringNum value=super.get(key);
			if(calls.getAndIncrement()<2){
				try{barrier.await(10, TimeUnit.SECONDS);}
				catch(Exception e){throw new RuntimeException("First-lookup barrier failed", e);}
			}
			return value;
		}
		private final CyclicBarrier barrier=new CyclicBarrier(2);
		private final AtomicInteger calls=new AtomicInteger();
		private static final long serialVersionUID=1L;
	}

	/** Invokes an existing private path and preserves the underlying failure type. */
	private static Object call(Object target, String name, Class<?>[] types, Object...args) throws Exception{
		final Method method=(target==null ? AnalyzeFlowCell.class : target.getClass()).getDeclaredMethod(name, types);
		method.setAccessible(true);
		try{return method.invoke(target, args);}
		catch(InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(cause instanceof Exception){throw (Exception)cause;}
			if(cause instanceof Error){throw (Error)cause;}
			throw new RuntimeException(cause);
		}
	}

	/** Installs fixture-owned state without modifying production APIs. */
	private static void set(Object target, String name, Object value) throws Exception{
		final Field field=target.getClass().getDeclaredField(name);
		field.setAccessible(true); field.set(target, value);
	}

	/** Checks exact finite results, independently of Java assertion-enable flags. */
	private static void check(String name, long observed, long expected){
		require(observed==expected, name+": observed="+observed+", expected="+expected);
		System.out.println(name+"\t"+observed+"\t"+expected+"\tPASS");
	}

	/** Fails loudly if a fixture's expected result is not observed. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
