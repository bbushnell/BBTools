package test.tax;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.TimeUnit;

import fileIO.ReadWrite;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import tax.SplitByTaxa;
import tax.TaxTree;

/** Checks baseline and fixed literal taxon-name substitution in SplitByTaxa outputs.
 * @author Sucrose */
public final class SplitByTaxaOutputNameLiteralTest{

	/** Runs parent or child checks for output-name literal replacement.
	 * @param args Optional expect=baseline|fixed or child-mode arguments
	 * @throws Exception if fixture setup, child execution, or cleanup fails */
	public static void main(String[] args) throws Exception{
		if(args.length>0 && "child".equals(args[0])){
			runChildCase(args[1], Boolean.parseBoolean(args[2]));
			return;
		}

		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		boolean success=false;
		final Path root=Files.createTempDirectory("split-by-taxa-literal-parent-");
		try{
			runParentCase(root, "normal", baseline, false);
			runParentCase(root, "normal", baseline, true);
			runParentCase(root, "backslash", baseline, false);
			runParentCase(root, "backslash", baseline, true);
			runParentCase(root, "dollar", baseline, false);
			runParentCase(root, "dollar", baseline, true);
			success=true;
		}finally{
			if(success){deleteTempTree(root);}
			else{System.err.println("Retained SplitByTaxa child logs in "+root);}
		}
		System.out.println("SPLIT_BY_TAXA_OUTPUT_NAME_LITERAL_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private SplitByTaxaOutputNameLiteralTest(){}

	/** Runs one bounded child JVM and validates its oracle.
	 * @param root Parent-owned temporary directory for child logs
	 * @param caseName Child case selector
	 * @param baseline Whether to assert the baseline or fixed expectation
	 * @param assertions Whether to enable child assertions
	 * @throws Exception if the child cannot be launched, read, or verified */
	private static void runParentCase(final Path root, final String caseName,
			final boolean baseline, final boolean assertions) throws Exception{
		final String label=caseName+(assertions ? "-ea" : "-da");
		final Path log=root.resolve(label+".log");
		Process child=null;
		final long started=System.nanoTime();
		try{
			final ProcessBuilder pb=javaProcess(assertions,
					SplitByTaxaOutputNameLiteralTest.class.getName(), "child", caseName,
					Boolean.toString(baseline));
			pb.redirectErrorStream(true);
			pb.redirectOutput(log.toFile());
			child=pb.start();
			child.getOutputStream().close();
			if(!child.waitFor(CHILD_TIMEOUT_MS, TimeUnit.MILLISECONDS)){
				terminateChild(child, label, log);
				throw new AssertionError("Timed out waiting for "+label+" child: "+
						readAscii(log));
			}
			final String output=readAscii(log);
			final long elapsedMs=TimeUnit.NANOSECONDS.toMillis(System.nanoTime()-started);
			System.out.println(label+" child exit="+child.exitValue()+
					" elapsed_ms="+elapsedMs+" output='"+compact(output)+"'");
			check(child.exitValue()==0,
					"Child "+label+" failed with exit "+child.exitValue()+": "+output);
			check(output.contains("case="+caseName),
					"Child "+label+" did not report an oracle marker: "+output);
		}finally{
			if(child!=null && child.isAlive()){terminateChild(child, label, log);}
		}
	}

	/** Runs one literal-name case inside a fresh child JVM.
	 * @param caseName Case selector
	 * @param baseline Whether to assert the baseline or fixed expectation
	 * @throws Exception if the case diverges from the oracle */
	private static void runChildCase(final String caseName, final boolean baseline)
			throws Exception{
		if("normal".equals(caseName)){
			checkResolvedPath("normal", 3, "split-gamma_delta-tail-%.fq");
		}else if("backslash".equals(caseName)){
			checkResolvedPath("backslash", 4, "split-backslash-tail-%.fq");
		}else if("dollar".equals(caseName)){
			checkDollarName(baseline);
		}else{
			throw new IllegalArgumentException("Unknown child case: "+caseName);
		}
	}

	/** Requires a normal sanitized taxon name to produce one exact output path.
	 * @param caseName Case name for child diagnostics
	 * @param taxID TaxID whose tree row should be read
	 * @param outputName Expected output filename relative to the temp directory
	 * @throws Exception if SplitByTaxa does not create exactly the expected writer */
	private static void checkResolvedPath(final String caseName, final int taxID,
			final String outputName) throws Exception{
		final CaseResult result=runCase(taxID);
		check(result.failure==null, caseName+" should process: "+result.failure);
		check(result.writerCount==1, caseName+" should open one writer");
		check(result.pathCount==1, caseName+" should open one output path");
		check(result.returnedLists==2, caseName+" should return read and tail lists");
		final String expected=result.root.resolve(outputName).toString();
		check(expected.equals(result.onlyPath),
				caseName+" should substitute the first % only; expected "+expected+
				" but found "+result.onlyPath);
		System.out.println("case="+caseName+" path="+result.onlyPath);
	}

	/** Requires the baseline dollar-name failure or fixed literal-dollar path.
	 * @param baseline Whether to assert the baseline or fixed expectation
	 * @throws Exception if the dollar-sign taxon oracle diverges */
	private static void checkDollarName(final boolean baseline) throws Exception{
		final CaseResult result=runCase(5);
		if(baseline){
			check(result.failure instanceof IllegalArgumentException,
					"Baseline dollar replacement should throw IllegalArgumentException: "+
					result.failure);
			checkStackContains(result.failure, "java.lang.String", "replaceFirst");
			checkStackContains(result.failure, "tax.SplitByTaxa", "makeOutputFile");
			check(result.writerCount==0,
					"Baseline dollar replacement should fail before opening a writer");
			check(result.pathCount==0,
					"Baseline dollar replacement should fail before creating a path");
			check(result.returnedLists==0,
					"Baseline dollar replacement should fail before returning a list");
			System.out.println("case=dollar baseline_exception="+
					result.failure.getClass().getSimpleName());
		}else{
			check(result.failure==null, "Fixed dollar replacement should process: "+
					result.failure);
			check(result.writerCount==1, "Fixed dollar replacement should open one writer");
			check(result.pathCount==1, "Fixed dollar replacement should open one path");
			check(result.returnedLists==2,
					"Fixed dollar replacement should return read and tail lists");
			final String expected=result.root.resolve("split-cash$taxon-tail-%.fq")
					.toString();
			check(expected.equals(result.onlyPath),
					"Fixed dollar replacement should preserve literal $ and second %; "+
					"expected "+expected+" but found "+result.onlyPath);
			System.out.println("case=dollar path="+result.onlyPath);
		}
	}

	/** Runs processInner over one explicit taxID and captures the writer it constructs.
	 * @param taxID TaxID embedded in the fixture read's header
	 * @return Writer counts, output path, and any processInner failure
	 * @throws Exception if setup or cleanup fails */
	private static CaseResult runCase(final int taxID) throws Exception{
		final Path root=Files.createTempDirectory("split-by-taxa-literal-");
		final HashMap<String, ConcurrentReadOutputStream> map=
				new HashMap<String, ConcurrentReadOutputStream>();
		final FakeCris cris=new FakeCris(taxID);
		Throwable failure=null;
		try{
			final Path tree=root.resolve("tree.tsv");
			final Path input=root.resolve("input.fa");
			Files.write(tree, treeText().toBytes());
			Files.write(input, new byte[0]);
			final SplitByTaxa split=new SplitByTaxa(new String[]{
					"in="+input,
					"out="+root.resolve("split-%-tail-%.fq"),
					"tree="+tree,
					"level=species",
					"ow=t"
			});
			final TaxTree loaded=(TaxTree)TREE_FIELD.get(split);
			check(loaded.getNode(taxID)!=null,
					"Loaded tree did not preserve taxID "+taxID);
			try{
				PROCESS_INNER.invoke(split, cris, map);
			}catch(InvocationTargetException e){
				failure=e.getCause();
			}
			final HashSet<String> paths=new HashSet<String>();
			String onlyPath=null;
			for(ConcurrentReadOutputStream ros : map.values()){
				onlyPath=ros.ff1.name();
				paths.add(onlyPath);
			}
			return new CaseResult(root, failure, map.size(), paths.size(),
					onlyPath, cris.returnedLists);
		}finally{
			boolean closeError=false;
			try{
				for(ConcurrentReadOutputStream ros : map.values()){
					closeError|=ReadWrite.closeStream(ros);
				}
			}finally{
				deleteTempTree(root);
			}
			check(!closeError, "A SplitByTaxa test writer failed to close cleanly");
		}
	}

	/** Builds a TaxTreeText tree containing literal output-name test cases.
	 * @return Serialized TaxTreeText bytes */
	private static ByteBuilder treeText(){
		final ByteBuilder bb=new ByteBuilder();
		bb.append("##TaxTree2\t1\n##Encoding\tA48\n##Tree\t6\t5\n");
		bb.append("##Config\t0\t0\t0\t0\t0\n##MergedMap\t0\t0\n");
		bb.append("##LevelExtended\tTID\tPID\tNumChildren\t");
		bb.append("MinParentLevelExtended\tMaxChildLevelExtended\tFlag\tName\n");
		row(bb, 1, 1, TaxTree.LIFE_E, 1, "Life");
		row(bb, 2, 1, TaxTree.DOMAIN_E, 3, "Bacteria");
		row(bb, 3, 2, TaxTree.SPECIES_E, 0, "gamma delta");
		row(bb, 4, 2, TaxTree.SPECIES_E, 0, "back\\slash");
		row(bb, 5, 2, TaxTree.SPECIES_E, 0, "cash$taxon");
		return bb;
	}

	/** Appends one A48 TaxTreeText row for the fixed test tree.
	 * @param bb Destination tree text
	 * @param id Node TaxID
	 * @param parent Parent TaxID
	 * @param rank Extended taxonomic rank
	 * @param children Number of direct children
	 * @param name Literal taxon name, including any filename-sensitive characters */
	private static void row(final ByteBuilder bb, final int id, final int parent,
			final int rank, final int children, final String name){
		bb.appendA48(rank).tab().appendA48(id).tab().appendA48(parent)
				.tab().appendA48(children).tab().appendA48(TaxTree.LIFE_E)
				.tab().appendA48(0).tab().appendA48(rank).tab().append(name).nl();
	}

	/** Creates a bounded child JVM command.
	 * @param assertions Whether to enable child assertions
	 * @param className Fully qualified class name to run
	 * @param args Arguments passed after the class name
	 * @return configured process builder */
	private static ProcessBuilder javaProcess(final boolean assertions,
			final String className, final String... args){
		final ArrayList<String> command=new ArrayList<String>();
		command.add(System.getProperty("java.home")+"/bin/java");
		command.add(assertions ? "-ea" : "-da");
		command.add("-Xmx128m");
		command.add("-XX:ActiveProcessorCount=1");
		command.add("-XX:+ExitOnOutOfMemoryError");
		addVectorModule(command);
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(className);
		for(String arg : args){command.add(arg);}
		return new ProcessBuilder(command);
	}

	/** Adds the Vector incubator module when this JDK exposes it.
	 * @param command Child JVM command to extend */
	private static void addVectorModule(final ArrayList<String> command){
		try{
			final Class<?> moduleLayer=Class.forName("java.lang.ModuleLayer");
			final Object boot=moduleLayer.getMethod("boot").invoke(null);
			final Object optional=moduleLayer.getMethod("findModule", String.class)
					.invoke(boot, "jdk.incubator.vector");
			final boolean present=((Boolean)optional.getClass()
					.getMethod("isPresent").invoke(optional)).booleanValue();
			if(present){
				command.add("--add-modules");
				command.add("jdk.incubator.vector");
			}
		}catch(ClassNotFoundException e){
			return;
		}catch(Exception e){
			throw new RuntimeException("Could not probe jdk.incubator.vector", e);
		}
	}

	/** Terminates a child JVM with bounded graceful and forced waits.
	 * @param child Child process to terminate
	 * @param label Child case label for diagnostics
	 * @param log Child log retained for diagnostics
	 * @throws Exception if the child cannot be terminated or its log cannot be read */
	private static void terminateChild(final Process child, final String label,
			final Path log) throws Exception{
		child.destroy();
		if(child.waitFor(DESTROY_TIMEOUT_MS, TimeUnit.MILLISECONDS)){return;}

		child.destroyForcibly();
		if(child.waitFor(DESTROY_TIMEOUT_MS, TimeUnit.MILLISECONDS)){return;}

		throw new AssertionError("Could not terminate "+label+" child: "+readAscii(log));
	}

	/** Deletes an owned temporary tree after the tested writers have closed.
	 * @param root Temporary directory whose descendants are owned by this fixture
	 * @throws Exception if directory traversal or deletion fails */
	private static void deleteTempTree(final Path root) throws Exception{
		if(!Files.exists(root)){return;}
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths.sorted(java.util.Comparator.reverseOrder())::iterator){
				Files.deleteIfExists(path);
			}
		}
	}

	/** Reads an ASCII text file.
	 * @param path Path to read
	 * @return File text decoded as ASCII
	 * @throws Exception if reading fails */
	private static String readAscii(final Path path) throws Exception{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}

	/** Compacts child output into one parent log line.
	 * @param text Text to compact
	 * @return Single-line representation of the text */
	private static String compact(final String text){
		return text.replace('\n', '|').replace('\r', '|');
	}

	/** Checks that a stack contains a specific frame.
	 * @param thrown Throwable to inspect
	 * @param className Required frame class name
	 * @param methodName Required frame method name */
	private static void checkStackContains(final Throwable thrown, final String className,
			final String methodName){
		for(StackTraceElement frame : thrown.getStackTrace()){
			if(className.equals(frame.getClassName()) &&
					methodName.equals(frame.getMethodName())){return;}
		}
		throw new AssertionError("Expected stack to contain "+className+"."+
				methodName+" but found "+java.util.Arrays.toString(thrown.getStackTrace()));
	}

	/** Keeps test checks active under -da.
	 * @param condition Required condition
	 * @param message Failure message */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Reflected package-private SplitByTaxa.processInner method. */
	private static final Method PROCESS_INNER;
	/** Reflected private SplitByTaxa.tree field. */
	private static final Field TREE_FIELD;

	static{
		try{
			PROCESS_INNER=SplitByTaxa.class.getDeclaredMethod("processInner",
					ConcurrentReadInputStream.class, HashMap.class);
			PROCESS_INNER.setAccessible(true);
			TREE_FIELD=SplitByTaxa.class.getDeclaredField("tree");
			TREE_FIELD.setAccessible(true);
		}catch(Exception e){throw new RuntimeException(e);}
	}

	/** Captured behavior from one processInner invocation. */
	private static final class CaseResult{
		/** Captures the outcome without changing the tested streams.
		 * @param root_ Temporary directory used to resolve exact expected paths
		 * @param failure_ Cause from reflective invocation, or null on normal return
		 * @param writerCount_ Number of raw-name entries in the output map
		 * @param pathCount_ Number of distinct first-mate output path strings
		 * @param onlyPath_ Only first-mate output path, or null when no writer opened
		 * @param returnedLists_ Number of returnList calls observed by the fake input */
		CaseResult(final Path root_, final Throwable failure_, final int writerCount_,
				final int pathCount_, final String onlyPath_, final int returnedLists_){
			root=root_;
			failure=failure_;
			writerCount=writerCount_;
			pathCount=pathCount_;
			onlyPath=onlyPath_;
			returnedLists=returnedLists_;
		}
		/** Temporary directory used to resolve exact expected paths. */
		final Path root;
		/** Invocation failure, or null if processInner returned normally. */
		final Throwable failure;
		/** Output-map entries, counted before stream cleanup. */
		final int writerCount;
		/** Distinct first-mate path strings among the constructed output streams. */
		final int pathCount;
		/** Only first-mate output path, or null when no writer opened. */
		final String onlyPath;
		/** Input-list returns observed during processing. */
		final int returnedLists;
	}

	/** Finite read stream that yields exactly one nonempty list and one empty tail. */
	private static final class FakeCris extends ConcurrentReadInputStream{
		/** Creates one unpaired read with an explicit taxonomy header.
		 * @param taxID TaxID embedded in the read header */
		FakeCris(final int taxID){
			super("SplitByTaxaOutputNameLiteralTest");
			reads.add(new Read("ACGT".getBytes(StandardCharsets.US_ASCII),
					null, "tid|"+taxID+"|", 1));
		}
		@Override
		public ListNum<Read> nextList(){
			if(next==0){next++; return new ListNum<Read>(reads, 0);}
			if(next==1){next++; return new ListNum<Read>(new ArrayList<Read>(), 1);}
			return null;
		}
		@Override
		public void returnList(final long listNum, final boolean poison){returnedLists++;}
		@Override
		public void run(){}
		@Override
		public void shutdown(){}
		@Override
		public void restart(){next=0;}
		@Override
		public void close(){}
		@Override
		public boolean paired(){return false;}
		@Override
		public Object[] producers(){return new Object[0];}
		@Override
		public boolean errorState(){return false;}
		@Override
		public void setSampleRate(final float rate, final long seed){}
		@Override
		public long basesIn(){return 0;}
		@Override
		public long readsIn(){return 0;}
		@Override
		public boolean verbose(){return false;}

		/** Owned one-read payload returned as the first list. */
		final ArrayList<Read> reads=new ArrayList<Read>();
		/** Delivery cursor and cumulative returnList call count. */
		int next=0, returnedLists=0;
	}

	/** Maximum time to wait for one child JVM. */
	private static final long CHILD_TIMEOUT_MS=30000;
	/** Time to wait for graceful child termination after a timeout or parent failure. */
	private static final long DESTROY_TIMEOUT_MS=5000;
}
