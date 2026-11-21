package test.jgi;

import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jgi.BBMask;
import parse.Parse;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * Finite BBMask SAM worker lifecycle checks for borrowed-batch return and fatal
 * worker failure.  Direct rows call the private SAM helper with a synthetic CRIS;
 * the fixed-only child exercises the actual MaskSamThread without real producer
 * queues or owner close semantics.
 * @author Jean
 */
public final class BBMaskWorkerLifecycleTest{

	/** Requires phase=baseline|fixed, a new outdir and the BBTools root. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null, rootPath=null, mode="parent", receipt=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("bbtools")){rootPath=pair[1];}
			else if(pair[0].equals("mode")){mode=pair[1];}
			else if(pair[0].equals("receipt")){receipt=pair[1];}
			else{throw new IllegalArgumentException("Unknown BBMask worker lifecycle option: "+arg);}
		}
		if(mode.equals("worker-child")){
			require(receipt!=null && !receipt.isEmpty(), "Set receipt=<path> for worker-child.");
			actualWorkerChild(Paths.get(receipt).toAbsolutePath().normalize());
			return;
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		require(rootPath!=null && !rootPath.isEmpty(), "Set bbtools=<repo root>.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/bbmaskworkerlifecycletest.sh with assertions enabled.");
		final Path out=Paths.get(directory).toAbsolutePath().normalize();
		final Path repo=Paths.get(rootPath).toAbsolutePath().normalize();
		Files.createDirectory(out);
		final BBMaskWorkerLifecycleTest test=new BBMaskWorkerLifecycleTest(out, repo,
			phase.equals("baseline"));
		test.runDirectRows();
		test.runActualWorkerRow();
		test.check(test.directCases==7, "Expected seven direct helper cases: "+test.directCases);
		test.check(test.childCases==(test.baseline ? 0 : 1), "Unexpected fixed-only child count: "+test.childCases);
		final String summary="BBMASK_WORKER_LIFECYCLE_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" direct="+test.directCases+" children="+test.childCases+
			" skipped="+test.skipped+" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(out.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+out);}
		System.out.print(summary);
	}

	/** Creates the evidence table and binds the repository under test. */
	private BBMaskWorkerLifecycleTest(final Path out_, final Path repo_) throws IOException{
		out=out_; repo=repo_; ref=out.resolve("tiny-ref.fa");
		require(Files.isDirectory(repo.resolve("current")), "Missing BBTools current/ at "+repo);
		Files.write(ref, REF.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		Files.write(out.resolve("cases.tsv"), ("case\tphase\tprimaryType\texpectedPrimaryType"+
			"\tprimary\texpectedPrimary\tsuppressed\texpectedSuppressed"+
			"\tsuppressedCount\texpectedSuppressedCount\treturns\texpectedReturns\tmatch\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	private BBMaskWorkerLifecycleTest(final Path out_, final Path repo_, final boolean baseline_)
			throws IOException{
		this(out_, repo_);
		baseline=baseline_;
	}

	/** Runs all finite direct-helper ownership rows. */
	private void runDirectRows() throws Exception{
		direct("normal-batch-then-terminal", cris(batch(1, unmapped()), terminal(2)),
			null, "", "1:false,2:true", 0, null, "");
		direct("null-eof", cris(), null, "", "", 0, null, "");
		direct("empty-terminal", cris(terminal(3)), null, "", "3:true", 0, null, "");
		direct("mapped-null-samline", cris(batch(4, mappedWithoutSamLine())),
			AssertionError.class, "No sam line", baseline ? "" : "4:false", 0, null, "");
		direct("nextlist-fails-after-return", cris(batch(5, unmapped()),
			new Sentinel("nextList sentinel")), Sentinel.class, "nextList sentinel",
			"5:false", 0, null, "");
		direct("return-fails-no-retry", cris(new Sentinel("return sentinel"),
			batch(6, unmapped())), Sentinel.class, "return sentinel", "6:false", 0, null, "");
		direct("processing-and-return-fail", cris(new Sentinel("cleanup sentinel"),
			batch(7, mappedWithoutSamLine())), AssertionError.class, "No sam line",
			baseline ? "" : "7:false", baseline ? 0 : 1, baseline ? null : Sentinel.class,
			baseline ? "" : "cleanup sentinel");
	}

	/** Invokes maskSam and checks primary, suppressed and return ownership. */
	private void direct(final String name, final CountingCris cris,
			final Class<? extends Throwable> expectedPrimaryType, final String expectedPrimary,
			final String expectedReturns, final int expectedSuppressedCount,
			final Class<? extends Throwable> expectedSuppressedType,
			final String expectedSuppressed) throws Exception{
		directCases++;
		final Throwable primary=invokeMaskSam(cris);
		final String primaryText=primary==null ? "" : String.valueOf(primary.getMessage());
		final String primaryType=primary==null ? "" : primary.getClass().getName();
		final String expectedPrimaryTypeName=expectedPrimaryType==null ? "" :
			expectedPrimaryType.getName();
		final String suppressedText=suppressedMessages(primary);
		final int suppressedCount=primary==null ? 0 : primary.getSuppressed().length;
		final boolean suppressedTypes=suppressedTypeMatches(primary, expectedSuppressedType);
		final String returns=join(cris.returns, ",");
		final boolean match=primaryType.equals(expectedPrimaryTypeName) &&
			containsOrEmpty(primaryText, expectedPrimary) &&
			containsOrEmpty(suppressedText, expectedSuppressed) &&
			suppressedCount==expectedSuppressedCount && suppressedTypes &&
			returns.equals(expectedReturns);
		Files.write(out.resolve("cases.tsv"), (name+"\t"+(baseline ? "baseline" : "fixed")+
			"\t"+primaryType+"\t"+expectedPrimaryTypeName+"\t"+
			clean(primaryText)+"\t"+clean(expectedPrimary)+"\t"+
			clean(suppressedText)+"\t"+clean(expectedSuppressed)+"\t"+
			suppressedCount+"\t"+expectedSuppressedCount+"\t"+returns+"\t"+
			expectedReturns+"\t"+match+"\n").getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.APPEND);
		check(match, name+" did not match expected primary/suppressed/returns: primary="+
			primaryType+"/"+primaryText+", suppressed="+suppressedText+", returns="+returns);
	}

	/** Runs the fixed-only actual MaskSamThread fatal child. */
	private void runActualWorkerRow() throws Exception{
		final Path caseDir=out.resolve("actual-worker-fatal");
		Files.createDirectory(caseDir);
		Files.write(caseDir.resolve("expected.tsv"), ("runInBaseline\texit\treceipt\tstderr\n"+
			"false\t1\ttrue\tNo sam line + BBMask SAM worker failed\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		if(baseline){
			skipped++;
			Files.write(caseDir.resolve("SKIPPED.txt"),
				"Fixed-only MaskSamThread fatal path is not executed in baseline.\n"
				.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
			return;
		}
		childCases++;
		final Path receipt=caseDir.resolve("return-receipt.tsv");
		final Path stdout=caseDir.resolve("stdout.txt");
		final Path stderr=caseDir.resolve("stderr.txt");
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(repo.resolve("dev/bbmaskworkerlifecycletest.sh").toString());
		command.add("mode=worker-child");
		command.add("receipt="+receipt);
		Files.write(caseDir.resolve("argv.txt"), join(command, "\n").getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		final ProcessBuilder builder=new ProcessBuilder(command).directory(repo.toFile())
			.redirectOutput(stdout.toFile()).redirectError(stderr.toFile());
		final int exit=builder.start().waitFor();
		Files.write(caseDir.resolve("exit.txt"), (Integer.toString(exit)+"\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		final byte[] outBytes=Files.readAllBytes(stdout);
		final String errText=new String(Files.readAllBytes(stderr), StandardCharsets.UTF_8);
		final boolean receiptPresent=Files.isRegularFile(receipt) &&
			Arrays.equals(Files.readAllBytes(receipt),
				"8\tfalse\n".getBytes(StandardCharsets.UTF_8));
		final boolean stderrMatches=errText.contains("java.lang.AssertionError: No sam line") &&
			errText.contains("java.lang.Exception: "+FATAL_DIAGNOSTIC) &&
			!errText.contains("BBMASK_WORKER_LIFECYCLE_PASS") &&
			!errText.contains("MaskSamThread returned without fatal KillSwitch");
		check(exit==1, "actual-worker-fatal child exit was "+exit);
		check(receiptPresent, "actual-worker-fatal did not return its borrowed batch before halt");
		check(stderrMatches, "actual-worker-fatal stderr missing original cause or KillSwitch message: "+errText);
		check(outBytes.length==0, "actual-worker-fatal wrote stdout");
		Files.write(out.resolve("cases.tsv"), ("actual-worker-fatal\tfixed\texit="+exit+
			"\texit=1\tstderr="+stderrMatches+"\tstderr=true\treceipt="+
			receiptPresent+"\treceipt=true\t"+(exit==1 && receiptPresent && stderrMatches &&
			outBytes.length==0)+"\n").getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.APPEND);
	}

	/** The fixed worker fatal path exits here through KillSwitch.kill(String). */
	private static void actualWorkerChild(final Path receipt) throws Exception{
		final CountingCris cris=cris(receipt, batch(8, mappedWithoutSamLine()));
		final Class<?> workerType=Class.forName("jgi.BBMask$MaskSamThread");
		final Constructor<?> constructor=workerType.getDeclaredConstructor(BBMask.class,
			ConcurrentReadInputStream.class);
		constructor.setAccessible(true);
		final Thread worker=(Thread)constructor.newInstance(owner(receipt.getParent()), cris);
		worker.start();
		worker.join();
		throw new AssertionError("MaskSamThread returned without fatal KillSwitch");
	}

	/** Calls BBMask.maskSam(ConcurrentReadInputStream) by reflection. */
	private Throwable invokeMaskSam(final CountingCris cris) throws Exception{
		final Method method=BBMask.class.getDeclaredMethod("maskSam", ConcurrentReadInputStream.class);
		method.setAccessible(true);
		try{
			method.invoke(owner(out), cris);
			return null;
		}catch(InvocationTargetException e){
			return e.getCause();
		}
	}

	/** Constructs BBMask without processing the reference or opening producer queues. */
	private static BBMask owner(final Path directory) throws Exception{
		final Path ref=directory.resolve("tiny-ref.fa");
		if(!Files.isRegularFile(ref)){
			Files.write(ref, REF.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		}
		return new BBMask(new String[]{"in="+ref, "out=null", "mr=f", "me=f",
			"minscaf=1", "ow=t", "t=1", "simd=f"});
	}

	private static ListNum<Read> batch(final long id, final Read... reads){
		return new ListNum<Read>(new ArrayList<Read>(Arrays.asList(reads)), id);
	}

	private static ListNum<Read> terminal(final long id){
		return new ListNum<Read>(new ArrayList<Read>(), id);
	}

	private static Read unmapped(){
		final Read r=read();
		r.setMapped(false);
		return r;
	}

	private static Read mappedWithoutSamLine(){
		final Read r=read();
		r.setMapped(true);
		return r;
	}

	private static Read read(){
		return new Read("ACGT".getBytes(StandardCharsets.US_ASCII), null, "read", 1, false);
	}

	private static CountingCris cris(final Object... events){return new CountingCris(null, null, events);}

	private static CountingCris cris(final Sentinel returnFailure, final Object... events){
		return new CountingCris(returnFailure, null, events);
	}

	private static CountingCris cris(final Path receipt, final Object... events){
		return new CountingCris(null, receipt, events);
	}

	private static String suppressedMessages(final Throwable t){
		if(t==null || t.getSuppressed().length<1){return "";}
		final ArrayList<String> list=new ArrayList<String>();
		for(final Throwable s : t.getSuppressed()){list.add(String.valueOf(s.getMessage()));}
		return join(list, ",");
	}

	private static boolean suppressedTypeMatches(final Throwable primary,
			final Class<? extends Throwable> expectedType){
		if(primary==null || primary.getSuppressed().length<1){return expectedType==null;}
		if(expectedType==null){return false;}
		for(final Throwable t : primary.getSuppressed()){
			if(t.getClass()!=expectedType){return false;}
		}
		return true;
	}

	private static boolean containsOrEmpty(final String text, final String token){
		return token==null || token.isEmpty() ? text.isEmpty() : text.contains(token);
	}

	private static String join(final List<String> list, final String delimiter){
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<list.size(); i++){
			if(i>0){sb.append(delimiter);}
			sb.append(list.get(i));
		}
		return sb.toString();
	}

	private static String clean(final String text){
		return (text==null ? "" : text).replace('\t', ' ').replace('\n', ';').replace('\r', ' ');
	}

	private void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println("FAIL: "+message);}
	}

	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Synthetic CRIS with finite nextList events and exact return accounting. */
	private static final class CountingCris extends ConcurrentReadInputStream{

		CountingCris(final Sentinel returnFailure_, final Path receipt_, final Object... events_){
			super("counting");
			returnFailure=returnFailure_;
			receipt=receipt_;
			events=events_;
		}

		@Override
		public ListNum<Read> nextList(){
			if(next>=events.length){return null;}
			final Object event=events[next++];
			if(event instanceof RuntimeException){throw (RuntimeException)event;}
			return (ListNum<Read>)event;
		}

		@Override
		public void returnList(final long listNum, final boolean poison){
			returns.add(listNum+":"+poison);
			if(receipt!=null){
				try{
					Files.write(receipt, (listNum+"\t"+poison+"\n").getBytes(StandardCharsets.UTF_8),
						StandardOpenOption.CREATE_NEW);
				}catch(IOException e){throw new RuntimeException(e);}
			}
			if(returnFailure!=null){throw returnFailure;}
		}

		@Override public void run(){}
		@Override public void shutdown(){}
		@Override public void restart(){next=0;}
		@Override public void close(){closed=true;}
		@Override public boolean paired(){return false;}
		@Override public Object[] producers(){return new Object[0];}
		@Override public boolean errorState(){return false;}
		@Override public void setSampleRate(final float rate, final long seed){}
		@Override public long basesIn(){return 0;}
		@Override public long readsIn(){return 0;}
		@Override public boolean verbose(){return false;}

		final Object[] events;
		final Sentinel returnFailure;
		final Path receipt;
		final ArrayList<String> returns=new ArrayList<String>();
		int next=0;
		boolean closed=false;
	}

	/** Distinct fixture sentinel exception. */
	private static final class Sentinel extends RuntimeException{
		Sentinel(final String message){super(message);}
		private static final long serialVersionUID=1L;
	}

	private static final String REF=">ref\nACGT\n";
	private static final String FATAL_DIAGNOSTIC=
		"BBMask SAM worker failed; aborting because input cleanup cannot safely recover from worker failure.";

	private final Path out, repo, ref;
	private boolean baseline;
	private int directCases, childCases, skipped, checks, failures;
}
