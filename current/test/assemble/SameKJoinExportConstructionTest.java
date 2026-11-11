package test.assemble;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.lang.management.GarbageCollectorMXBean;
import java.lang.management.ManagementFactory;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.List;

import fileIO.ByteStreamWriter;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Tests ownership of real plain-file outputs when same-K exporter construction fails.
 * Run one case per fresh JVM on Linux. Reflection accesses package-private production
 * code without widening public visibility. Exact-path descriptors, not process-wide
 * descriptor totals, establish closure. A GC during a failure observation invalidates
 * setup because collection could hide an unreachable, unclosed stream.
 * Startup cases inject named synchronous failures into the real start/header helper
 * using real plain-file writers. They are not ordinary-input reproductions. Their
 * observations precede fixture rescue cleanup, which must finish both workers.
 * This is a constructor/resource check, not an assembly or throughput test.
 * @author Brian Bushnell, Sucrose
 */
public final class SameKJoinExportConstructionTest{

	/** Accepts the four original constructor cases or a named inject-* startup case, plus outdir=NEW. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, SameKJoinExportConstructionTest.class, false);
		try{
			String name=null;
			Path root=null;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(pair.length==2 && !pair[1].isEmpty(), "Expected nonempty flag=value: "+arg);
				if(pair[0].equalsIgnoreCase("case")){name=pair[1];}
				else if(pair[0].equalsIgnoreCase("outdir")){root=Paths.get(pair[1]).toAbsolutePath().normalize();}
				else{throw new IllegalArgumentException("Unknown fixture option: "+arg);}
			}
			check(root!=null && ("feature".equals(name) || "pair".equals(name) ||
					"first-failure".equals(name) || "second-failure".equals(name) || injectionCase(name)),
					"Required a supported constructor or inject-* case and outdir=NEW.");
			Files.createDirectory(root);
			new SameKJoinExportConstructionTest().run(name, root);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Resolves existing constructors and GC counters before measuring resource ownership. */
	private SameKJoinExportConstructionTest() throws Exception{
		check(Files.isDirectory(FD_ROOT), "Fixture requires Linux /proc/self/fd.");
		final Class<?> momentsClass=Class.forName("assemble.SameKJoinFeatures$Moments");
		momentsConstructor=momentsClass.getDeclaredConstructor();
		momentsConstructor.setAccessible(true);
		exportConstructor=Class.forName("assemble.SameKJoinExport").getDeclaredConstructor(
				String.class, String.class, int.class, int.class, momentsClass);
		exportConstructor.setAccessible(true);
		startup=Class.forName("assemble.SameKJoinExport").getDeclaredMethod("startAndWriteHeaders",
				ByteStreamWriter.class, ByteStreamWriter.class, momentsClass, ByteBuilder.class);
		startup.setAccessible(true);
		collectors=ManagementFactory.getGarbageCollectorMXBeans();
		check(!collectors.isEmpty(), "GC counters are required to exclude collection-masked closure.");
		gcCounts();
	}

	/** Warms the actual constructor/FD observation path, then records one independent case. */
	private void run(final String name, final Path root) throws Exception{
		final Path warmup=root.resolve("warmup");
		Files.createDirectory(warmup);
		healthy(warmup, true);
		if(injectionCase(name)){injectedStartup(name, root);return;}
		final boolean failure=name.endsWith("failure");
		final int remaining;
		if(failure){remaining=failedConstruction(root, name.equals("first-failure"));}
		else{remaining=healthy(root, name.equals("pair"));}
		final boolean passed=remaining==0;
		final String report="case\texpected_writer_fds\tobserved_writer_fds\tpassed\n"+
				name+"\t0\t"+remaining+"\t"+passed+"\n";
		write(root.resolve("checks.tsv"), report);
		System.out.print(report);
		System.out.println("EXPORT_CONSTRUCTION case="+name+" assertions="+
				SameKJoinExportConstructionTest.class.desiredAssertionStatus()+" setup=true passed="+passed);
		if(!passed){throw new AssertionError("Same-K exporter left a writer open after construction failure: "+remaining);}
	}

	/** Selects four baseline lifecycle points and one additional fixed-only cleanup-failure case. */
	private static boolean injectionCase(final String name){
		return "inject-before-first-start".equals(name) || "inject-before-second-start".equals(name) ||
				"inject-after-second-start".equals(name) || "inject-header".equals(name) ||
				"inject-cleanup-failure".equals(name);
	}

	/** Records the real helper's failure ownership before fixture cleanup, then checks the desired contract. */
	private void injectedStartup(final String name, final Path root) throws Exception{
		final Path feature=root.resolve("features.tsv"), control=root.resolve("controls.tsv");
		final RuntimeException marker=new RuntimeException("Injected exporter startup failure: "+name);
		final RuntimeException cleanupMarker=name.equals("inject-cleanup-failure") ?
				new RuntimeException("Injected exporter cleanup failure after feature writer closed.") : null;
		FaultWriter featureWriter=null, controlWriter=null;
		int featureFds=-1, controlFds=-1;
		boolean workersStopped=false;
		try{
			featureWriter=new FaultWriter(feature, name.equals("inject-before-first-start") ? "before-start" :
					name.equals("inject-header") || cleanupMarker!=null ? "header" : null, marker, cleanupMarker);
			controlWriter=new FaultWriter(control, name.equals("inject-before-second-start") ? "before-start" :
					name.equals("inject-after-second-start") ? "after-start" : null, marker, null);
			check(countWriterFds(feature)==1 && countWriterFds(control)==1, "Injection requires two real open plain-file outputs.");
			final long[] before=gcCounts();
			Throwable observed=null;
			try{startup.invoke(null, featureWriter, controlWriter, momentsConstructor.newInstance(), new ByteBuilder());}
			catch(InvocationTargetException e){observed=e.getCause();}
			featureFds=countWriterFds(feature);controlFds=countWriterFds(control);
			final long[] after=gcCounts();
			workersStopped=!featureWriter.isAlive() && !controlWriter.isAlive();
			write(root.resolve("observation.tsv"),
					"feature_fds_before_rescue\tcontrol_fds_before_rescue\tfeature_state\tcontrol_state\tfeature_alive\tcontrol_alive\texception_identity\n"+
					featureFds+"\t"+controlFds+"\t"+featureWriter.getState()+"\t"+controlWriter.getState()+"\t"+
					featureWriter.isAlive()+"\t"+controlWriter.isAlive()+"\t"+(observed==marker)+"\n");
			final StringBuilder gc=new StringBuilder("collector\tbefore\tafter\n");
			for(int i=0; i<before.length; i++){
				gc.append(collectors.get(i).getName()).append('\t').append(before[i]).append('\t').append(after[i]).append('\n');
			}
			write(root.resolve("gc.tsv"), gc.toString());
			if(observed!=null){
				final StringWriter stack=new StringWriter();
				final PrintWriter writer=new PrintWriter(stack);
				observed.printStackTrace(writer);writer.flush();
				write(root.resolve("construction-error.txt"), stack.toString());
			}
			check(Arrays.equals(before, after), "GC changed during injected-failure observation.");
			check(observed==marker, "Expected the original injected exception object, observed: "+observed);
			if(cleanupMarker!=null){
				final Throwable[] suppressed=observed.getSuppressed();
				final boolean exact=suppressed.length==1 && suppressed[0]==cleanupMarker;
				write(root.resolve("cleanup.tsv"), "cleanup_thrown\tsuppressed_count\tsuppressed_identity\n"+
						featureWriter.cleanupThrown+"\t"+suppressed.length+"\t"+exact+"\n");
				check(featureWriter.cleanupThrown && exact,
						"Original startup failure must retain the exact post-close cleanup failure as its only suppressed exception.");
			}
		}finally{
			// Always attempt both: a failure rescuing the first must not skip the second.
			try{rescue(featureWriter);}finally{rescue(controlWriter);}
		}
		final int rescuedFeature=countWriterFds(feature), rescuedControl=countWriterFds(control);
		write(root.resolve("rescue.tsv"), "feature_fds_after_rescue\tcontrol_fds_after_rescue\tfeature_alive\tcontrol_alive\n"+
				rescuedFeature+"\t"+rescuedControl+"\t"+featureWriter.isAlive()+"\t"+controlWriter.isAlive()+"\n");
		check(rescuedFeature==0 && rescuedControl==0 && !featureWriter.isAlive() && !controlWriter.isAlive(),
				"Fixture rescue must release both outputs and finish both owned workers.");
		final boolean passed=featureFds==0 && controlFds==0 && workersStopped;
		final String report="case\texpected_writer_fds\tfeature_fds_before_rescue\tcontrol_fds_before_rescue\tworkers_stopped\tpassed\n"+
				name+"\t0\t"+featureFds+"\t"+controlFds+"\t"+workersStopped+"\t"+passed+"\n";
		write(root.resolve("checks.tsv"), report);System.out.print(report);
		System.out.println("EXPORT_STARTUP_INJECTION case="+name+" assertions="+
				SameKJoinExportConstructionTest.class.desiredAssertionStatus()+" identity=true rescued=true passed="+passed);
		if(!passed){throw new AssertionError("Same-K exporter startup cleanup contract failed: "+name);}
	}

	/** Fixture-only recovery uses actual thread state, never whether start() returned. */
	private static void rescue(final ByteStreamWriter writer){
		if(writer==null){return;}
		final boolean failed=writer.getState()==Thread.State.NEW ? writer.closeBeforeStart() : writer.poisonAndWait();
		check(!failed && !writer.isAlive(), "Fixture rescue failed to finish an owned plain-file writer.");
	}

	/** Real output/worker implementation with named faults; optional cleanup fault follows actual closure. */
	private static final class FaultWriter extends ByteStreamWriter{
		/** Opens a fresh plain output; point may be null for an ordinary writer. */
		FaultWriter(final Path path, final String point_, final RuntimeException marker_, final RuntimeException cleanupMarker_){
			super(path.toString(), false, false, false);
			point=point_;marker=marker_;cleanupMarker=cleanupMarker_;
		}
		/** The after-start fault occurs only after the real worker's startup handshake. */
		@Override
		public synchronized void start(){
			if("before-start".equals(point)){throw marker;}
			super.start();
			if("after-start".equals(point)){throw marker;}
		}
		/** Fails before the first feature-header append; the control header path stays real. */
		@Override
		public ByteStreamWriter print(final ByteBuilder row){
			if("header".equals(point)){throw marker;}
			return super.print(row);
		}
		/** Throws once after the real shutdown/join, so a later fixture rescue remains safe. */
		@Override
		public boolean poisonAndWait(){
			final boolean failed=super.poisonAndWait();
			if(cleanupMarker!=null && !cleanupThrown){
				check(!failed && !isAlive() && getState()==Thread.State.TERMINATED,
						"Cleanup injection must follow successful real writer shutdown and join.");
				cleanupThrown=true;
				throw cleanupMarker;
			}
			return failed;
		}
		/** Selected failure point and the exact exception object the helper must preserve. */
		private final String point;
		private final RuntimeException marker;
		/** Null for baseline cases; otherwise thrown once after real feature cleanup. */
		private final RuntimeException cleanupMarker;
		/** Producer-owned injection latch; repeated fixture rescue delegates without another fault. */
		private boolean cleanupThrown;
	}

	/** Confirms live descriptors are visible, then exact headers and both closures. */
	private int healthy(final Path root, final boolean pair) throws Exception{
		final Path feature=root.resolve("features.tsv"), control=root.resolve("controls.tsv");
		check(countWriterFds(feature)==0 && countWriterFds(control)==0, "New control paths already have writers.");
		try(AutoCloseable export=create(feature, pair ? control : null)){
			check(countWriterFds(feature)==1, "Live exporter must own exactly one plain feature writer.");
			check(countWriterFds(control)==(pair ? 1 : 0), "Unexpected live control-writer ownership.");
		}
		final int featureFds=countWriterFds(feature), controlFds=countWriterFds(control);
		check(featureFds==0 && controlFds==0, "Healthy exporter close failed to release its writers.");
		checkBytes(feature, featureHeader());
		if(pair){checkBytes(control, CONTROL_HEADER);}
		else{check(!Files.exists(control), "Feature-only construction created a control output.");}
		write(root.resolve("observation.tsv"), "feature_fds_after\tcontrol_fds_after\n"+
				featureFds+"\t"+controlFds+"\n");
		return featureFds+controlFds;
	}

	/** Captures an actual second/first open failure before a GC can hide a leaked writer. */
	private int failedConstruction(final Path root, final boolean first) throws Exception{
		final Path blocker=root.resolve("blocker");
		write(blocker, BLOCKER);
		final Path feature=first ? blocker.resolve("features.tsv") : root.resolve("features.tsv");
		final Path control=first ? root.resolve("controls.tsv") : blocker.resolve("controls.tsv");
		final Path rejected=first ? feature : control;
		check(countWriterFds(feature)==0 && countWriterFds(control)==0, "Failure paths already have writers.");
		// Resolve paths before the observed interval; both may legitimately be absent.
		feature.toFile().getCanonicalPath();
		control.toFile().getCanonicalPath();
		final long[] before=gcCounts();
		Throwable failure=null;
		AutoCloseable unexpected=null;
		try{unexpected=create(feature, control);}
		catch(RuntimeException e){failure=e;}
		finally{if(unexpected!=null){unexpected.close();}}
		final int featureFds=countWriterFds(feature), controlFds=countWriterFds(control);
		final long[] after=gcCounts();
		final StringBuilder gc=new StringBuilder("collector\tbefore\tafter\n");
		for(int i=0; i<before.length; i++){
			gc.append(collectors.get(i).getName()).append('\t').append(before[i]).append('\t').append(after[i]).append('\n');
		}
		write(root.resolve("gc.tsv"), gc.toString());
		write(root.resolve("observation.tsv"), "feature_fds_after\tcontrol_fds_after\n"+
				featureFds+"\t"+controlFds+"\n");
		if(failure!=null){
			final StringWriter stack=new StringWriter();
			final PrintWriter writer=new PrintWriter(stack);
			failure.printStackTrace(writer);
			writer.flush();
			write(root.resolve("construction-error.txt"), stack.toString());
		}
		check(Arrays.equals(before, after), "GC changed during the failure observation; cleanup result is invalid.");
		check(failure!=null && failure.getClass()==RuntimeException.class &&
				failure.getCause() instanceof FileNotFoundException &&
				failure.getCause().getMessage()!=null && failure.getCause().getMessage().contains(rejected.toString()),
				"Expected the actual blocked-output FileNotFoundException, observed: "+failure);
		checkBytes(blocker, BLOCKER);
		check(controlFds==0 && !Files.exists(control), "Failed setup opened or created the control output.");
		if(first){check(!Files.exists(feature), "First-open failure created a feature output.");}
		else{check(Files.isRegularFile(feature) && Files.size(feature)==0, "Unstarted feature output must exist and remain empty.");}
		return featureFds;
	}

	/** Invokes the real exporter and rethrows its original unchecked constructor failure. */
	private AutoCloseable create(final Path feature, final Path control) throws Exception{
		try{
			return (AutoCloseable)exportConstructor.newInstance(feature.toString(),
					control==null ? null : control.toString(), 31, 4, momentsConstructor.newInstance());
		}catch(InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(cause instanceof RuntimeException){throw (RuntimeException)cause;}
			if(cause instanceof Error){throw (Error)cause;}
			throw e;
		}
	}

	/** Records every collector; unavailable counters cannot establish an unchanged-GC interval. */
	private long[] gcCounts(){
		final long[] counts=new long[collectors.size()];
		for(int i=0; i<counts.length; i++){
			counts[i]=collectors.get(i).getCollectionCount();
			check(counts[i]>=0, "Unavailable GC collection count: "+collectors.get(i).getName());
		}
		return counts;
	}

	/** Counts writable descriptors for exactly one canonical pathname, including absent paths.
	 * Only disappearance between directory enumeration and readlink is ignored. */
	private static int countWriterFds(final Path file) throws Exception{
		final Path real=file.toFile().getCanonicalFile().toPath();
		int count=0;
		try(DirectoryStream<Path> fds=Files.newDirectoryStream(FD_ROOT)){
			for(Path fd : fds){
				final Path target;
				try{target=Files.readSymbolicLink(fd);}
				catch(NoSuchFileException disappeared){continue;}
				if(target.equals(real)){
					final int access=fdFlags(fd)&3;
					check(access!=3, "Invalid descriptor access mode: "+fd);
					if(access==1 || access==2){count++;}
				}
			}
		}
		return count;
	}

	/** Reads one actual procfs octal flags field; missing or malformed metadata is a setup error. */
	private static int fdFlags(final Path fd) throws Exception{
		final Path info=Paths.get("/proc/self/fdinfo", fd.getFileName().toString());
		int flags=-1;
		for(String line : Files.readAllLines(info, StandardCharsets.US_ASCII)){
			if(line.startsWith("flags:")){
				check(flags==-1, "Duplicate descriptor flags: "+info);
				flags=Integer.parseInt(line.substring(6).trim(), 8);
			}
		}
		check(flags>=0, "Missing descriptor flags: "+info);
		return flags;
	}

	/** Independent literal schema expectation for empty, valid depth moments. */
	private static String featureHeader(){
		final StringBuilder text=new StringBuilder("#schema=samek27-v1 unique=0 observations=0.0 squares=0.0\nkind\tevent");
		for(int i=0; i<27; i++){text.append("\tf").append(i);}
		return text.append('\n').toString();
	}

	/** Checks complete tiny ASCII fixtures, not merely file existence or line counts. */
	private static void checkBytes(final Path path, final String expected) throws Exception{
		check(Arrays.equals(Files.readAllBytes(path), expected.getBytes(StandardCharsets.US_ASCII)),
				"Unexpected bytes: "+path);
	}

	/** Saves a fresh tiny artifact without overwriting prior evidence. */
	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Setup/oracle checks remain active under both assertion modes. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Reflective entry points into unchanged package-private production code. */
	private final Constructor<?> exportConstructor, momentsConstructor;
	/** Real constructor-used helper; reflection leaves production API visibility unchanged. */
	private final Method startup;
	/** Fixed collector roster sampled around each failure observation. */
	private final List<GarbageCollectorMXBean> collectors;
	/** Linux process descriptor directory; never use whole-process descriptor totals. */
	private static final Path FD_ROOT=Paths.get("/proc/self/fd");
	/** Immutable regular-file parent used to make child output creation impossible. */
	private static final String BLOCKER="unchanged blocked parent\n";
	/** Independent literal control-context header from the documented output schema. */
	private static final String CONTROL_HEADER="event\tk\tunitig\tcut\tsource_len\tdest_len\tsource_word\tdest_word\tleft_context\tright_context\n";
}
