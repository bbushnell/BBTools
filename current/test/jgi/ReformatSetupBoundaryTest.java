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
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.ReadStreamByteWriter;
import stream.SamLine;
import structures.ListNum;

/** Finite constructor checks for SAM setup and QUAL path collisions.
 * Invalid configurations are never processed or allowed to open output streams.
 * The paired mode reuses the original output oracle independently of its SAM test.
 * Additional modes check finite ordered-singleton completion and skip-aware input counts.
 * Sampling modes check skipped-prefix length thresholds and filter-generated singletons.
 * Each invocation requires a fresh JVM and new output directory.
 * @author Shinobu, Sucrose
 */
public final class ReformatSetupBoundaryTest{

	/** Requires phase=baseline|fixed, variant=streamer|cris, a test mode and a fresh outdir. */
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
		if(!"setup".equals(mode) && !"paired".equals(mode) && !"singletons".equals(mode) && !"skipcount".equals(mode) &&
			!"histogramskip".equals(mode) && !"sampledsingletons".equals(mode)){
			throw new IllegalArgumentException("Set mode=setup|paired|singletons|skipcount|histogramskip|sampledsingletons");
		}
		if("singletons".equals(mode) && !"cris".equals(variant)){
			throw new IllegalArgumentException("Singleton control requires variant=cris");
		}
		if(directory==null){throw new IllegalArgumentException("Set outdir to a new directory.");}
		final boolean baseline=phase.equals("baseline"), streamer=variant.equals("streamer");
		final Class<?> type=Class.forName(streamer ? "jgi.ReformatStreamer" : "jgi.Reformat2");
		final Path root=Files.createDirectory(Paths.get(directory)).toAbsolutePath();
		if(mode.equals("setup")){
			samSetup(type, root, baseline && !streamer);
			qualPaths(type, root, baseline);
		}else if(mode.equals("singletons")){
			crisSingletons(type, root, baseline);
		}else if(mode.equals("skipcount")){
			skipCounts(type, root, baseline, streamer);
		}else if(mode.equals("histogramskip")){
			histogramSkip(type, root, baseline);
		}else if(mode.equals("sampledsingletons")){
			sampledSingletons(type, root, baseline);
		}else{
			final Method method=ReformatDriverBoundaryTest.class.getDeclaredMethod("pairedOutput", Class.class, Path.class, boolean.class);
			method.setAccessible(true);
			try{method.invoke(null, type, root, baseline);}
			catch(InvocationTargetException e){throw new AssertionError("Paired output control failed", e.getCause());}
			checks+=2; //The frozen helper checks exact FASTQ bytes and the public individual-read count.
		}
		System.out.println("REFORMAT_SETUP_BOUNDARY_PASS phase="+phase+" variant="+variant+" mode="+mode+" checks="+checks);
	}

	/** Four ordered batches leave only two pending lists on baseline, well below capacity.
	 * Closing the stream diagnoses the missing IDs; no blocked-add/timeout oracle is used.
	 */
	private static void crisSingletons(final Class<?> type, final Path root, final boolean baseline) throws Exception{
		final Path input=root.resolve("singleton-input.fq"), output=root.resolve("singletons.fq");
		Files.write(input, "@seed\nACGT\n+\nIIII\n".getBytes(StandardCharsets.US_ASCII));
		final String[] args={"in="+input, "out=null", "outsingle="+output, "t=2", "wt=2", "tin=1", "tout=1",
			"ordered=t", "int=f", "minlength=3", "requirebothbad=t", "qin=33", "qout=33", "simd=f", "ow=f"};
		final Object driver=type.getConstructor(String[].class).newInstance((Object)args);
		final Method make=type.getDeclaredMethod("makeRosb", boolean.class);
		make.setAccessible(true);
		final ConcurrentReadOutputStream singles=(ConcurrentReadOutputStream)make.invoke(driver, true);
		try{
			check(singles.ordered && singles.ff1.ordered(), "The singleton control must exercise actual ordered CROS output.");
			final Field field=type.getDeclaredField("processor");
			field.setAccessible(true);
			final ReformatProcessor proc=((ReformatProcessor)field.get(driver)).clone();
			final Class<?> workerType=Class.forName("jgi.Reformat2$ProcessThread");
			final Constructor<?> constructor=workerType.getDeclaredConstructor(type, ConcurrentReadInputStream.class,
				ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class, ReformatProcessor.class, int.class);
			constructor.setAccessible(true);
			final Object worker=constructor.newInstance(driver, null, null, singles, proc, 0);
			final Method process=workerType.getDeclaredMethod("processList", ListNum.class);
			process.setAccessible(true);
			for(int i=0; i<4; i++){
				final String left=i<2 ? "ACGT" : "AC", right=i==0 || i==2 ? "TGCA" : "TG";
				final byte[] q1=new byte[left.length()], q2=new byte[right.length()];
				Arrays.fill(q1, (byte)40);
				Arrays.fill(q2, (byte)40);
				final Read r1=new Read(left.getBytes(StandardCharsets.US_ASCII), q1, "p"+i+"/1", i);
				final Read r2=new Read(right.getBytes(StandardCharsets.US_ASCII), q2, "p"+i+"/2", i);
				r1.mate=r2;
				r2.mate=r1;
				r1.setPairnum(0);
				r2.setPairnum(1);
				final ArrayList<Read> reads=new ArrayList<Read>();
				reads.add(r1);
				try{process.invoke(worker, new ListNum<Read>(reads, i));}
				catch(InvocationTargetException e){throw new AssertionError("Singleton processing failed", e.getCause());}
			}
		}finally{
			singles.close();
			try{singles.join();}
			catch(AssertionError e){
				final StackTraceElement first=e.getStackTrace()[0];
				if(!baseline || !singles.errorState() || !first.getClassName().equals("stream.ConcurrentGenericReadOutputStream") ||
					!first.getMethodName().equals("join")){throw e;}
				//The real join has already joined both backend threads before checking its pending table.
				System.out.println("CRIS_SINGLETON_BASELINE_PENDING_IDS");
			}
		}
		check(singles.errorState()==baseline && (baseline || singles.finishedSuccessfully()),
			"Ordered singleton close must diagnose baseline ID gaps and complete successfully after repair.");
		final byte[] actual=Files.exists(output) ? Files.readAllBytes(output) : new byte[0];
		final String expected=baseline ? "" : "@p1/1\nACGT\n+\nIIII\n@p2/2\nTGCA\n+\nIIII\n";
		check(expected.equals(new String(actual, StandardCharsets.US_ASCII)),
			"Ordered singleton output must contain exactly the two surviving mates after repair.");
	}

	/** A precount excludes the skipped prefix within the same raw-input quota as processing. */
	private static void skipCounts(final Class<?> type, final Path root, final boolean baseline,
		final boolean streamer) throws Exception{
		final Path input=root.resolve("skipcount-input.fq");
		Files.write(input, ("@r0\nACGT\n+\nIIII\n@r1\nTGCA\n+\nIIII\n"+
			"@r2\nAAAA\n+\nIIII\n@r3\nCCCC\n+\nIIII\n").getBytes(StandardCharsets.US_ASCII));
		final String[] args={"in="+input, "out=null", "skipreads=1", "t=2", "wt=1", "tin=1", "tout=1",
			"int=f", "minlength=0", "qin=33", "qout=33", "simd=f"};
		final Object driver=type.getConstructor(String[].class).newInstance((Object)args);
		final Method count=type.getDeclaredMethod("countReads", long.class);
		count.setAccessible(true);
		final long[] unlimited, limited;
		try{
			unlimited=(long[])count.invoke(driver, -1L);
			limited=(long[])count.invoke(driver, 2L);
		}catch(InvocationTargetException e){throw new AssertionError("Skip precount failed", e.getCause());}
		final long[] expectedAll=baseline ? new long[]{4, 4, 16} : new long[]{3, 3, 12};
		final long[] expectedLimit=baseline ? (streamer ? new long[]{4, 4, 16} : new long[]{2, 2, 8}) : new long[]{1, 1, 4};
		check(unlimited.length>=3 && unlimited[0]==expectedAll[0] && unlimited[1]==expectedAll[1] && unlimited[2]==expectedAll[2],
			"Unlimited skip precount must match literal fragment/read/base counts: "+Arrays.toString(unlimited));
		check(limited.length>=3 && limited[0]==expectedLimit[0] && limited[1]==expectedLimit[1] && limited[2]==expectedLimit[2],
			"A raw two-record quota with one skipped record must leave one fragment: "+Arrays.toString(limited));
	}

	/** A skipped longest read must not set the threshold for the retained suffix. */
	private static void histogramSkip(final Class<?> type, final Path root, final boolean baseline) throws Exception{
		final Path input=root.resolve("histogram-input.fq");
		Files.write(input, ("@r0\nAAAAAAAA\n+\nIIIIIIII\n@r1\nACGT\n+\nIIII\n"+
			"@r2\nAC\n+\nII\n@r3\nACGTAC\n+\nIIIIII\n").getBytes(StandardCharsets.US_ASCII));
		for(int i=0; i<2; i++){
			final Path output=root.resolve("histogram-"+i+".fq");
			final String[] args={"in="+input, "out="+output, "skipreads=1", "srt=1", "prioritizelength=t",
				"reads="+(i==0 ? -1 : 2), "t=2", "wt=1", "tin=1", "tout=1", "int=f", "qin=33", "qout=33", "simd=f", "ow=f"};
			runDriver(type, args);
			final String expected=baseline ? "" : (i==0 ? "@r3\nACGTAC\n+\nIIIIII\n" : "@r1\nACGT\n+\nIIII\n");
			check(expected.equals(new String(Files.readAllBytes(output), StandardCharsets.US_ASCII)),
				"Length selection must use the unskipped suffix within the raw quota: "+output);
		}
	}

	/** N filtering produces two singletons without using excluded trim/length/quality options. */
	private static void sampledSingletons(final Class<?> type, final Path root, final boolean baseline) throws Exception{
		final Path input=root.resolve("sampled-singletons-input.fq");
		Files.write(input, ("@p0/1\nACGT\n+\nIIII\n@p0/2\nTGCA\n+\nIIII\n"+
			"@p1/1\nACGT\n+\nIIII\n@p1/2\nNNNN\n+\nIIII\n"+
			"@p2/1\nNNNN\n+\nIIII\n@p2/2\nTGCA\n+\nIIII\n").getBytes(StandardCharsets.US_ASCII));
		final String pair="@p0/1\nACGT\n+\nIIII\n@p0/2\nTGCA\n+\nIIII\n";
		final String singles="@p1/1\nACGT\n+\nIIII\n@p2/2\nTGCA\n+\nIIII\n";
		for(int i=0; i<2; i++){
			final boolean exact=i==1;
			final Path output=root.resolve("sampled-pairs-"+i+".fq"), single=root.resolve("sampled-singles-"+i+".fq");
			final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+input, "out="+output, "outsingle="+single,
				"t=2", "wt=1", "tin=1", "tout=1", "int=t", "qin=33", "qout=33", "simd=f", "ow=f", "maxns=0", "requirebothbad=t"));
			if(exact){args.add("srt=3"); args.add("sampleseed=1");}
			final String report=runDriverWithReport(type, args.toArray(new String[args.size()]), root.resolve("sampled-report-"+i+".txt"));
			check(pair.equals(new String(Files.readAllBytes(output), StandardCharsets.US_ASCII)),
				"Selecting all available fragments must retain the complete passing pair.");
			check((baseline && exact ? "" : singles).equals(new String(Files.readAllBytes(single), StandardCharsets.US_ASCII)),
				"The two N-filter survivors must reach singleton output in ordinary and exact modes.");
			final Matcher outputCount=Pattern.compile("(?m)^Output:\\s+(\\d+) reads?\\b[^\\r\\n]*?\\s+(\\d+) bases?\\b").matcher(report);
			final long expectedReads=baseline && exact ? 2 : 4, expectedBases=baseline && exact ? 8 : 16;
			check(outputCount.find() && Long.parseLong(outputCount.group(1))==expectedReads &&
				Long.parseLong(outputCount.group(2))==expectedBases, "Public output totals must match both files: "+report);
		}
	}

	/** Saves and replays actual diagnostics while leaving public output-count representation unconstrained. */
	private static String runDriverWithReport(final Class<?> type, final String[] args, final Path diagnostics) throws Exception{
		final PrintStream previous=System.err;
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		try(PrintStream capture=new PrintStream(bytes, true, "UTF-8")){
			System.setErr(capture);
			runDriver(type, args);
		}finally{
			System.setErr(previous);
			Files.write(diagnostics, bytes.toByteArray());
			previous.print(new String(bytes.toByteArray(), StandardCharsets.UTF_8));
		}
		return new String(bytes.toByteArray(), StandardCharsets.UTF_8);
	}

	/** Invokes the real driver lifecycle; processing exceptions cannot become passing observations. */
	private static void runDriver(final Class<?> type, final String[] args) throws Exception{
		try{type.getMethod("main", String[].class).invoke(null, (Object)args);}
		catch(InvocationTargetException e){throw new AssertionError("Finite driver control failed", e.getCause());}
	}

	/** SAM forwarding must be initialized without requiring a reference genome. */
	private static void samSetup(final Class<?> type, final Path root, final boolean missing) throws Exception{
		final Path input=root.resolve("input.sam"), output=root.resolve("output.sam");
		Files.write(input, ("@HD\tVN:1.6\tSO:unknown\n"+
			"r0\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n").getBytes(StandardCharsets.US_ASCII));
		SamLine.SET_FROM_OK=false;
		ReadStreamByteWriter.USE_ATTACHED_SAMLINE=false;
		construct(type, input, output, null, null);
		check(SamLine.SET_FROM_OK==!missing, "SAM input must enable forwarding without a loaded reference.");
		check(ReadStreamByteWriter.USE_ATTACHED_SAMLINE==!missing, "SAM output must retain the input alignment record.");
		check(!Files.exists(output), "A constructor-only SAM check must not open output.");
	}

	/** QUAL aliases must be rejected before processing, even when overwrite is enabled. */
	private static void qualPaths(final Class<?> type, final Path root, final boolean baseline) throws Exception{
		final Path input=root.resolve("input.fa"), qual=root.resolve("input.qual");
		final byte[] bases=">r0\nACGT\n".getBytes(StandardCharsets.US_ASCII);
		final byte[] scores=">r0\n40 40 40 40\n".getBytes(StandardCharsets.US_ASCII);
		Files.write(input, bases);
		Files.write(qual, scores);
		final Path output=root.resolve("output.fa"), outputQual=root.resolve("output.qual");
		construct(type, input, output, qual, outputQual);
		check(!Files.exists(output) && !Files.exists(outputQual), "Valid preflight must not create output files.");
		alias(type, input, output, qual, qual, baseline);
		alias(type, input, qual, qual, outputQual, baseline);
		check(Arrays.equals(bases, Files.readAllBytes(input)), "Preflight must preserve the source FASTA bytes.");
		check(Arrays.equals(scores, Files.readAllBytes(qual)), "Preflight must preserve the source QUAL bytes.");
		check(!Files.exists(output) && !Files.exists(outputQual), "Rejected or baseline-accepted aliases must not open outputs.");
	}

	/** Only the explicit duplicate-path diagnostic establishes rejection of this alias. */
	private static void alias(final Class<?> type, final Path input, final Path output,
		final Path qual, final Path outputQual, final boolean baseline) throws Exception{
		boolean rejected=false;
		try{construct(type, input, output, qual, outputQual);}
		catch(InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(!(cause instanceof RuntimeException) || cause.getMessage()==null ||
				!cause.getMessage().contains("specified multiple times")){throw e;}
			rejected=true;
		}
		check(rejected==!baseline, "Input QUAL must not also be output: "+output+", "+outputQual);
	}

	/** Constructs only; no process/main method is invoked for any path-collision case. */
	private static void construct(final Class<?> type, final Path input, final Path output,
		final Path qual, final Path outputQual) throws Exception{
		final String[] args={"in="+input, "out="+output, "qfin="+(qual==null ? "null" : qual),
			"qfout="+(outputQual==null ? "null" : outputQual), "ow=t", "t=1", "wt=1", "tin=1", "tout=1", "int=f", "simd=f"};
		type.getConstructor(String[].class).newInstance((Object)args);
	}

	/** Counts independent expectations and fails regardless of assertion settings. */
	private static void check(final boolean condition, final String message){if(!condition){throw new AssertionError(message);} checks++;}
	private static int checks;
}
