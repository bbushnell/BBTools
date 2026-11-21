package test.jgi;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import parse.Parse;
import processor.ReformatProcessor;
import shared.Shared;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.SamLine;
import stream.Writer;
import structures.ListNum;

/** Finite driver boundaries, including a deterministic out-of-order worker schedule.
 * Recording output checks batch submission without deliberately hanging an ordered writer.
 * Each invocation runs in its own JVM; reflection stays confined to this test package.
 * @author Shinobu
 */
public final class ReformatDriverBoundaryTest{

	/** Requires phase, variant, mode and a new output directory with an existing parent. */
	public static void main(final String[] args) throws Exception{
		String phase=null, variant=null, mode=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException(arg);}
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("variant")){variant=pair[1];}
			else if(pair[0].equals("mode")){mode=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException(arg);}
		}
		if(!"baseline".equals(phase) && !"fixed".equals(phase)){throw new IllegalArgumentException("Set phase=baseline|fixed");}
		if(!"streamer".equals(variant) && !"cris".equals(variant)){throw new IllegalArgumentException("Set variant=streamer|cris");}
		if(!"methods".equals(mode) && !"sam".equals(mode)){throw new IllegalArgumentException("Set mode=methods|sam");}
		if(directory==null){throw new IllegalArgumentException("Set outdir to a new directory.");}
		final boolean baseline=phase.equals("baseline"), streamer=variant.equals("streamer");
		final Class<?> type=Class.forName(streamer ? "jgi.ReformatStreamer" : "jgi.Reformat2");
		final Path root=Files.createDirectory(Paths.get(directory)).toAbsolutePath();
		final Path input=root.resolve("input.fq");
		Files.write(input, ("@r0\nACGT\n+\nIIII\n@r1\nTGCA\n+\nIIII\n"+
			"@r2\nAAAA\n+\nIIII\n@r3\nCCCC\n+\nIIII\n").getBytes(StandardCharsets.US_ASCII));
		if(mode.equals("sam")){samSampling(type, root, baseline, streamer);}
		else{
			argumentValue(type, root, input, baseline);
			workerPolicy(type, root, input);
			skipOrder(type, root, input, baseline, streamer);
			pairedSampling(type, root, input);
			precount(type, root, input, baseline, streamer);
			if(streamer){singletonBatches(type, root, input, baseline);}
		}
		System.out.println("REFORMAT_DRIVER_BOUNDARY_PASS phase="+phase+" variant="+variant+" mode="+mode+" checks="+checks);
	}

	/** A value after the first equals sign remains one filename. */
	private static void argumentValue(final Class<?> type, final Path root, final Path input, final boolean baseline) throws Exception{
		final Path output=root.resolve("out=literal.fq");
		final Object driver=driver(type, input, output);
		check(get(driver, "out1").equals(baseline ? root.resolve("out").toString() : output.toString()),
			"An output filename must preserve embedded equals signs.");
	}

	/** Name history and exact-sampling state retain their existing single-worker policy. */
	private static void workerPolicy(final Class<?> type, final Path root, final Path input) throws Exception{
		check((Integer)get(driver(type, input, root.resolve("unique.fq"), "wt=4", "uniquenames=t"), "workers")==1,
			"Unique-name maps require serial processing even with an explicit larger worker request.");
		check((Integer)get(driver(type, input, root.resolve("sample.fq"), "wt=4", "srt=1"), "workers")==1,
			"Exact-sampling counters require serial processing.");
	}

	/** A later batch may execute before an earlier one; skip is an immutable input boundary. */
	private static void skipOrder(final Class<?> type, final Path root, final Path input,
		final boolean baseline, final boolean streamer) throws Exception{
		final Object driver=driver(type, input, root.resolve("skip.fq"), "wt=2", "skipreads=1");
		final ReformatProcessor proc=(ReformatProcessor)get(driver, "processor");
		final Object worker=streamer ? null : worker(driver, proc);
		final ListNum<Read> later=batch(1, read("later", "ACGT", 3));
		final ListNum<Read> earlier=batch(0, read("earlier", "ACGT", 0));
		final ListNum<Read> boundary=batch(2, read("boundary", "ACGT", 1));
		process(driver, worker, proc, later, null, null);
		process(driver, worker, proc, earlier, null, null);
		process(driver, worker, proc, boundary, null, null);
		check(nonnull(later)==1, "Input ordinal three is beyond skipreads=1.");
		check(nonnull(earlier)==(baseline ? 1 : 0), "Input ordinal zero must be skipped regardless of batch arrival order.");
		check(nonnull(boundary)==1, "The first ordinal after the skipped prefix must remain.");
	}

	/** One selected fragment retains both mates and contributes both base counts. */
	private static void pairedSampling(final Class<?> type, final Path root, final Path input) throws Exception{
		final Object driver=driver(type, input, root.resolve("paired.fq"), "srt=1");
		set(driver, "readsRemaining", 1L);
		set(driver, "randy", Shared.random(1));
		final Read left=read("pair/1", "ACGT", 0), right=read("pair/2", "TGCA", 0);
		pair(left, right);
		final ArrayList<Read> in=new ArrayList<Read>(), out=new ArrayList<Read>();
		in.add(left);
		invoke(driver, "sampleExact", new Class<?>[]{ArrayList.class, ArrayList.class, ArrayList.class}, in, out, null);
		check(out.size()==1 && out.get(0)==left && left.mate==right, "One selected fragment retains its two linked mates.");
		check((Long)get(driver, "sampledBasesOut")==8L, "Both four-base mates contribute to the base count.");
		check((Long)get(driver, "sampleReadsTarget")==0L && (Long)get(driver, "readsRemaining")==0L,
			"Selecting one fragment exhausts a one-fragment target and population.");
	}

	/** A precount must use the same input quota as the processing pass. */
	private static void precount(final Class<?> type, final Path root, final Path input,
		final boolean baseline, final boolean streamer) throws Exception{
		final Object driver=driver(type, input, root.resolve("limited.fq"), "reads=1");
		final long[] counts=(long[])invoke(driver, "countReads", new Class<?>[]{long.class}, 1L);
		final long expected=baseline && streamer ? 4 : 1;
		check(counts.length>=3 && counts[0]==expected && counts[1]==expected && counts[2]==4*expected,
			"The precount must stop at one unpaired record; observed "+Arrays.toString(counts));
		final long[] unlimited=(long[])invoke(driver, "countReads", new Class<?>[]{long.class}, -1L);
		check(unlimited.length>=3 && unlimited[0]==4 && unlimited[1]==4 && unlimited[2]==16,
			"An unlimited precount retains all four records and sixteen bases.");
	}

	/** Empty singleton batches carry ordering IDs even when no reads pass to that output. */
	private static void singletonBatches(final Class<?> type, final Path root, final Path input,
		final boolean baseline) throws Exception{
		final Object driver=driver(type, input, root.resolve("pairs.fq"), "wt=2", "minlength=3",
			"requirebothbad=t", "outsingle="+root.resolve("single.fq"));
		final ReformatProcessor proc=(ReformatProcessor)get(driver, "processor");
		final RecordingWriter pairs=new RecordingWriter(), singles=new RecordingWriter();
		for(int i=0; i<4; i++){
			final Read left=read("p"+i+"/1", i<2 ? "ACGT" : "AC", i);
			final Read right=read("p"+i+"/2", i==0 || i==2 ? "TGCA" : "TG", i);
			pair(left, right);
			process(driver, null, proc, batch(i, left), pairs, singles);
		}
		check(pairs.ids.equals(Arrays.asList(0L, 1L, 2L, 3L)), "Primary output receives every input batch ID.");
		check(singles.ids.equals(baseline ? Arrays.asList(1L, 2L) : Arrays.asList(0L, 1L, 2L, 3L)),
			"Ordered singleton output needs dense IDs including empty batches: "+singles.ids);
		check(singles.names.equals(Arrays.asList("p1/1", "p2/2")), "The surviving mate of each mixed pair reaches singleton output.");
	}

	/** Valid unmapped SAM records have the same one-record sampling target as reads. */
	private static void samSampling(final Class<?> type, final Path root, final boolean baseline,
		final boolean streamer) throws Exception{
		final Path input=root.resolve("input.sam"), output=root.resolve("sample.sam");
		Files.write(input, ("@HD\tVN:1.6\tSO:unknown\n"+
			"r0\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"+
			"r1\t4\t*\t0\t0\t*\t*\t0\t0\tTGCA\tIIII\n"+
			"r2\t4\t*\t0\t0\t*\t*\t0\t0\tAAAA\tIIII\n").getBytes(StandardCharsets.US_ASCII));
		final String[] args=options(input, output, "srt=1", "sampleseed=1");
		try{type.getMethod("main", String[].class).invoke(null, (Object)args);}
		catch(InvocationTargetException e){throw new AssertionError("SAM processing failed", e.getCause());}
		int records=0;
		for(final String line : Files.readAllLines(output, StandardCharsets.US_ASCII)){
			if(!line.isEmpty() && line.charAt(0)!='@'){records++;}
		}
		check(records==(baseline && streamer ? 3 : 1), "A one-record target applies to SAM output too; observed "+records);
		pairedOutput(type, root, baseline);
	}

	/** Checks public output statistics, without prescribing the sampler's internal count representation. */
	private static void pairedOutput(final Class<?> type, final Path root, final boolean baseline) throws Exception{
		final Path input=root.resolve("pair-input.fq"), output=root.resolve("pair-output.fq");
		final String expected="@pair/1\nACGT\n+\nIIII\n@pair/2\nTGCA\n+\nIIII\n";
		Files.write(input, expected.getBytes(StandardCharsets.US_ASCII));
		final String[] args=options(input, output, "srt=1", "int=t", "qin=33", "qout=33", "sampleseed=1");
		final PrintStream previous=System.err;
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		try(PrintStream capture=new PrintStream(bytes, true, "UTF-8")){
			System.setErr(capture);
			try{type.getMethod("main", String[].class).invoke(null, (Object)args);}
			catch(InvocationTargetException e){throw new AssertionError("Paired processing failed", e.getCause());}
		}finally{System.setErr(previous);}
		final String report=new String(bytes.toByteArray(), StandardCharsets.UTF_8);
		Files.write(root.resolve("paired-report.txt"), bytes.toByteArray());
		System.err.print(report);
		check(expected.equals(new String(Files.readAllBytes(output), StandardCharsets.US_ASCII)),
			"Sampling one available pair must preserve both complete FASTQ records.");
		final Matcher matcher=Pattern.compile("(?m)^Output:\\s+(\\d+) reads?\\b").matcher(report);
		check(matcher.find() && Long.parseLong(matcher.group(1))==(baseline ? 1L : 2L),
			"The output report must count both mates as individual reads: "+report);
	}

	private static Object driver(final Class<?> type, final Path input, final Path output, final String... extra) throws Exception{
		try{return type.getConstructor(String[].class).newInstance((Object)options(input, output, extra));}
		catch(InvocationTargetException e){throw new AssertionError("Driver construction failed", e.getCause());}
	}

	private static String[] options(final Path input, final Path output, final String... extra){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+input, "out="+output,
			"t=2", "wt=1", "tin=1", "tout=1", "int=f", "ow=f", "simd=f", "minlength=0"));
		args.addAll(Arrays.asList(extra));
		return args.toArray(new String[args.size()]);
	}

	private static Object worker(final Object driver, final ReformatProcessor proc) throws Exception{
		final Class<?> type=Class.forName("jgi.Reformat2$ProcessThread");
		final Constructor<?> constructor=type.getDeclaredConstructor(driver.getClass(), ConcurrentReadInputStream.class,
			ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class, ReformatProcessor.class, int.class);
		constructor.setAccessible(true);
		return constructor.newInstance(driver, null, null, null, proc.clone(), 0);
	}

	private static void process(final Object driver, final Object worker, final ReformatProcessor proc,
		final ListNum<Read> batch, final Writer main, final Writer singles) throws Exception{
		if(worker!=null){invoke(worker, "processList", new Class<?>[]{ListNum.class}, batch);}
		else{invoke(driver, "processReadList", new Class<?>[]{ListNum.class, ReformatProcessor.class, Writer.class, Writer.class}, batch, proc, main, singles);}
	}

	private static Object invoke(final Object owner, final String name, final Class<?>[] types, final Object... args) throws Exception{
		final Method method=owner.getClass().getDeclaredMethod(name, types);
		method.setAccessible(true);
		try{return method.invoke(owner, args);}
		catch(InvocationTargetException e){throw new AssertionError("Production call failed: "+name, e.getCause());}
	}
	private static Object get(final Object owner, final String name) throws Exception{
		final Field field=owner.getClass().getDeclaredField(name);
		field.setAccessible(true);
		return field.get(owner);
	}
	private static void set(final Object owner, final String name, final Object value) throws Exception{
		final Field field=owner.getClass().getDeclaredField(name);
		field.setAccessible(true);
		field.set(owner, value);
	}
	private static ListNum<Read> batch(final long id, final Read read){
		final ArrayList<Read> list=new ArrayList<Read>(); list.add(read); return new ListNum<Read>(list, id);
	}
	private static Read read(final String name, final String bases, final long id){
		return new Read(bases.getBytes(StandardCharsets.US_ASCII), null, name, id);
	}
	private static void pair(final Read left, final Read right){left.mate=right; right.mate=left; left.setPairnum(0); right.setPairnum(1);}
	private static int nonnull(final ListNum<Read> batch){int count=0; for(final Read read : batch){if(read!=null){count++;}} return count;}
	private static void check(final boolean condition, final String message){if(!condition){throw new AssertionError(message);} checks++;}

	/** Records submissions immediately; it does not emulate an ordered queue or launch threads. */
	private static final class RecordingWriter implements Writer{
		final ArrayList<Long> ids=new ArrayList<Long>();
		final ArrayList<String> names=new ArrayList<String>();
		public void add(final ArrayList<Read> reads, final long id){addReads(new ListNum<Read>(reads, id));}
		public void addReads(final ListNum<Read> batch){ids.add(batch.id); for(final Read read : batch){if(read!=null){names.add(read.id);}}}
		public void addLines(final ListNum<SamLine> lines){throw new UnsupportedOperationException("Read-batch fixture only.");}
		public void start(){}
		public long readsWritten(){return names.size();}
		public long basesWritten(){return 0;}
		public void poison(){}
		public boolean waitForFinish(){return false;}
		public boolean poisonAndWait(){return false;}
		public void finishError(){throw new AssertionError("Unexpected writer abandonment.");}
		public String fname(){return "recording";}
		public boolean errorState(){return false;}
		public boolean finishedSuccessfully(){return true;}
	}
	private static int checks;
}
