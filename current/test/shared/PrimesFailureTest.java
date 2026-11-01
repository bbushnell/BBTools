package test.shared;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.management.ManagementFactory;
import java.lang.reflect.InvocationTargetException;
import java.net.URL;
import java.net.URLClassLoader;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPOutputStream;

/** Isolated resource-read and makePrimes error/cleanup checks.
 * Every invocation runs one bounded child JVM. Only tiny input/output fixtures
 * are created; a class loader reads the existing runtime without copying it.
 * A synthetic Data.class resource URL redirects resource discovery into the
 * fixture directory while all actual class bytes come from the original runtime.
 * Linux /proc/self/fd is used only to inspect the child's own input descriptors.
 * @author Sucrose
 */
public final class PrimesFailureTest{

	/** Requires case=name work=directory; child=true is reserved for the harness. */
	public static void main(String[] args) throws Exception{
		String name=null, directory=null;
		boolean child=false;
		for(String arg : args){
			if(arg.startsWith("case=")){name=arg.substring(5);}
			else if(arg.startsWith("work=")){directory=arg.substring(5);}
			else if(arg.equals("child=true")){child=true;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(name==null || directory==null || !Arrays.asList(CASES).contains(name)){
			throw new IllegalArgumentException("Expected case="+Arrays.toString(CASES)+" work=directory");
		}
		final Path work=Paths.get(directory).toAbsolutePath();
		Files.createDirectories(work);
		if(child){
			try{runChild(name, work);System.out.println("PrimesFailureTest "+name+" PASS");System.exit(0);}
			catch(Throwable failure){failure.printStackTrace(System.err);System.exit(1);}
			return;
		}
		boolean assertions=false;
		assert(assertions=true);
		final ArrayList<String> command=new ArrayList<String>();
		command.add(Paths.get(System.getProperty("java.home"), "bin", "java").toString());
		command.add(assertions ? "-ea" : "-da");
		command.add("-XX:ActiveProcessorCount=1");
		command.add("-XX:+ExitOnOutOfMemoryError");
		command.add("-Xmx256m");
		for(String arg : ManagementFactory.getRuntimeMXBean().getInputArguments()){
			if(arg.startsWith("--add-modules=")){command.add(arg);}
		}
		command.add("-cp");command.add(System.getProperty("java.class.path"));
		command.add(PrimesFailureTest.class.getName());
		command.add("case="+name);command.add("work="+work);command.add("child=true");
		final Path log=work.resolve(name+"-child.log");
		require(!Files.exists(log), "Refusing to replace prior child log: "+log);
		final Process process=new ProcessBuilder(command).redirectErrorStream(true).redirectOutput(log.toFile()).start();
		process.getOutputStream().close();
		if(!process.waitFor(20, TimeUnit.SECONDS)){
			process.destroyForcibly();
			process.waitFor(5, TimeUnit.SECONDS);
			Files.copy(log, System.out);
			throw new AssertionError("Child exceeded finite test deadline: "+name);
		}
		Files.copy(log, System.out);
		require(process.exitValue()==0, "Child contract check failed: "+name+", exit="+process.exitValue());
	}

	/** Creates a private resource root, then checks one real Primes initialization or invocation. */
	private static void runChild(String name, Path work) throws Exception{
		final Path root=Files.createDirectory(work.resolve("isolated-"+name));
		final Path resources=Files.createDirectory(root.resolve("resources"));
		final Path resource=resources.resolve("primes.txt.gz");
		final byte[] compressed=gzip(PRIMES);
		if(name.equals("resource-corrupt")){compressed[compressed.length-8]^=1;}
		Files.write(resource, compressed);
		final URL runtime=PrimesFailureTest.class.getProtectionDomain().getCodeSource().getLocation();
		try(ResourceLoader loader=new ResourceLoader(runtime, root)){
			setBoolean(loader, "shared.Shared", "SIMD", false);
			setBoolean(loader, "fileIO.ReadWrite", "ALLOW_NATIVE_BGZF", false);
			setBoolean(loader, "fileIO.ReadWrite", "USE_UNBGZIP", false);
			setBoolean(loader, "fileIO.ReadWrite", "USE_UNPIGZ", false);
			setBoolean(loader, "fileIO.ReadWrite", "USE_GUNZIP", false);
			final Class<?> data=Class.forName("dna.Data", true, loader);
			final String selected=(String)data.getMethod("findPath", String.class).invoke(null, "?primes.txt.gz");
			require(Files.isSameFile(resource, Paths.get(selected)), "Resource lookup escaped the isolated fixture: "+selected);
			Class<?> primes=null;
			Throwable initializationFailure=null;
			try{primes=Class.forName("shared.Primes", true, loader);}
			catch(ExceptionInInitializerError failure){initializationFailure=failure.getCause();}
			if(name.equals("resource-corrupt")){
				System.out.println("resource initialization failure="+initializationFailure);
				require(initializationFailure instanceof RuntimeException, "Corrupt prime resource was accepted");
				return;
			}
			require(initializationFailure==null && primes!=null, "Valid resource initialization failed: "+initializationFailure);
			require(Arrays.equals((long[])primes.getField("primes").get(null), new long[] {2, 3, 5, 7, 11, 13}), "Valid resource contents differ");
			if(name.equals("resource-control")){return;}
			checkMakePrimes(name, work, primes);
		}
	}

	/** Tests literal output, exceptional reader cleanup, or asynchronous writer-flush failure.
	 * The injected failure is on borrowed stdout flush, not a regular-file close.
	 * ReadWrite.finishWriting owns flushing/checkError but preserves borrowed stdout.
	 * Return-time observations are reported separately; failure judgment waits for
	 * writer termination and proof that the injected flush actually ran. */
	private static void checkMakePrimes(String name, Path work, Class<?> primes) throws Exception{
		final boolean readError=name.equals("make-read-error"), parseError=name.equals("make-parse-error");
		final boolean writeError=name.equals("make-write-error");
		final Path input=work.resolve(readError ? "input.txt.gz" : "input.txt");
		if(readError){
			final byte[] compressed=gzip(PRIMES);compressed[compressed.length-8]^=1;Files.write(input, compressed);
		}else{
			Files.write(input, (parseError ? "2\nnot-a-number\n" : PRIMES).getBytes(StandardCharsets.US_ASCII));
		}
		final Path output=work.resolve("output.txt");
		final PrintStream original=System.out;
		final FlushFailure injected=new FlushFailure();
		final PrintStream failing=new PrintStream(injected);
		Throwable failure=null;
		int liveWriters, writersAtReturn;
		boolean flushAtReturn;
		try{
			if(writeError){System.setOut(failing);}
			try{
				primes.getMethod("makePrimes", String[].class).invoke(null,
						(Object)new String[] {input.toString(), writeError ? "stdout" : output.toString(), "2"});
			}catch(InvocationTargetException invocation){failure=invocation.getCause();}
			writersAtReturn=countWriters();
			flushAtReturn=injected.attempted;
			liveWriters=waitForWriters();
		}finally{System.setOut(original);}
		final int openInputs=countOpenInput(input);
		System.out.println("make result failure="+failure+" liveWriters="+liveWriters+" openInputs="+openInputs+
				" injectedFlush="+injected.attempted+" writersAtReturn="+writersAtReturn+" flushAtReturn="+flushAtReturn);
		require(liveWriters==0 && openInputs==0, "makePrimes leaked a writer thread or owned reader");
		if(writeError){
			require(injected.attempted, "Writer did not exercise injected flush failure");
			require(failure instanceof RuntimeException, "Writer/close failure was not propagated to makePrimes caller");
		}else if(parseError){
			require(failure instanceof NumberFormatException, "Original numeric parse failure was not preserved: "+failure);
		}else if(readError){
			require(failure instanceof RuntimeException, "Corrupt input was not rejected");
		}else{
			require(failure==null, "Valid makePrimes call failed: "+failure);
			require(Arrays.equals(Files.readAllBytes(output), "2\n5\n11\n".getBytes(StandardCharsets.US_ASCII)), "Spaced-prime output differs");
		}
	}

	/** Observes return-time state without making a scheduling-sensitive correctness claim. */
	private static int countWriters(){
		int count=0;
		for(Thread thread : Thread.getAllStackTraces().keySet()){
			if(thread.getClass().getName().equals("fileIO.TextStreamWriter") && thread.isAlive()){count++;}
		}
		return count;
	}

	/** Gives the sole tiny writer a bounded chance to terminate, then reports leaked threads. */
	private static int waitForWriters() throws InterruptedException{
		int live=0;
		for(Thread thread : Thread.getAllStackTraces().keySet()){
			if(thread.getClass().getName().equals("fileIO.TextStreamWriter")){
				thread.join(1500);
				if(thread.isAlive()){live++;}
			}
		}
		return live;
	}

	/** Counts descriptors pointing at this child's exact input path, without invoking GC. */
	private static int countOpenInput(Path input) throws IOException{
		int count=0;
		final Path exact=input.toRealPath();
		try(DirectoryStream<Path> descriptors=Files.newDirectoryStream(Paths.get("/proc/self/fd"))){
			for(Path descriptor : descriptors){
				try{if(Files.readSymbolicLink(descriptor).equals(exact)){count++;}}
				catch(NoSuchFileException closedDuringScan){/* A descriptor may close during enumeration. */}
			}
		}
		return count;
	}

	/** Configures the isolated runtime before resource loading; parent static state is untouched. */
	private static void setBoolean(ClassLoader loader, String owner, String name, boolean value) throws Exception{
		Class.forName(owner, true, loader).getField(name).setBoolean(null, value);
	}

	/** Produces a tiny gzip member; callers may flip its CRC while preserving payload bytes. */
	private static byte[] gzip(String text) throws IOException{
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		try(GZIPOutputStream gzip=new GZIPOutputStream(bytes)){gzip.write(text.getBytes(StandardCharsets.US_ASCII));}
		return bytes.toByteArray();
	}

	/** Enforces checks under either JVM assertion setting. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Uses genuine runtime classes with only Data's root-discovery URL substituted. */
	private static final class ResourceLoader extends URLClassLoader{
		/** Runtime classes load privately; the parent supplies JDK and optional module classes. */
		ResourceLoader(URL runtime, Path root) throws IOException{
			super(new URL[] {runtime}, ClassLoader.getSystemClassLoader());
			dataLocation=root.resolve("current/dna/Data.class").toUri().toURL();
			boolean assertions=false;assert(assertions=true);setDefaultAssertionStatus(assertions);
		}
		/** Loads BBTools classes from the real runtime first, preserving access to JDK modules. */
		@Override
		protected Class<?> loadClass(String name, boolean resolve) throws ClassNotFoundException{
			if(name.startsWith("java.") || name.startsWith("javax.") || name.startsWith("jdk.")
					|| name.startsWith("sun.") || name.startsWith("com.sun.")){return super.loadClass(name, resolve);}
			synchronized(getClassLoadingLock(name)){
				Class<?> found=findLoadedClass(name);
				if(found==null){
					try{found=findClass(name);}
					catch(ClassNotFoundException absent){found=super.loadClass(name, false);}
				}
				if(resolve){resolveClass(found);}
				return found;
			}
		}
		/** Returns the fixture location for root discovery; all other resources are genuine. */
		@Override
		public URL getResource(String name){return name.equals("dna/Data.class") ? dataLocation : super.getResource(name);}
		/** Synthetic location only; no class file is copied here. */
		private final URL dataLocation;
	}

	/** Discards writes but fails on flush, exposing errors reported during writer shutdown. */
	private static final class FlushFailure extends OutputStream{
		/** Data writes succeed; the test isolates the closure/flush path. */
		@Override
		public void write(int value){}
		/** Signals the failure before throwing so the test proves this path ran. */
		@Override
		public void flush() throws IOException{attempted=true;throw new IOException("Injected prime-output flush failure");}
		/** Volatile because the writer thread sets it and the main thread reads it. */
		volatile boolean attempted;
	}

	/** Literal sorted prime inputs, shared by valid and CRC-corrupt members. */
	private static final String PRIMES="2\n3\n5\n7\n11\n13\n";
	/** Each case is run in a fresh child process and private resource directory. */
	private static final String[] CASES={"resource-control", "resource-corrupt", "make-control", "make-read-error", "make-parse-error", "make-write-error"};
}
