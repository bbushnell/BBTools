package test.jgi;

import java.io.IOException;
import java.io.InputStream;
import java.io.InterruptedIOException;
import java.lang.management.LockInfo;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadInfo;
import java.lang.management.ThreadMXBean;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.LockSupport;

import fileIO.ByteFile1;
import fileIO.ReadWrite;
import jgi.BBMerge;
import jgi.BBMergeOverlapper;
import parse.Parse;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.Writer;

/** Real BBMerge owner-join interruption cases for JGI089, in fresh Linux children.
 * Three phase children receive zero/one/two signals; two main children receive zero/one.
 * Two output-disabled phase failures use a real unpaired FASTA remnant. No production
 * hooks or extracted wait helpers are used: finite standard input holds EOF while the
 * owner joins its actual MateThreads. Only that fixture-owned owner is interrupted.
 * Explicit t=2/simd=f/bf1=t selects synchronous FASTA readers and unthreaded plain
 * FASTA outputs. Quiescence assertions apply only to those observed backends.
 * Owner status is captured immediately at return/throw, before fixture observations.
 * Main flag checks observe final values; redirected descriptor closure is separate.
 * Main's asynchronous histogram is drained later by an uninterrupted observer; that
 * counted drain is not a join or completion at main return. No arbitrary codec,
 * early-close, missing-output-batch, merger-accuracy or throughput claim is made.
 * @author Sucrose
 */
public final class BBMergeOwnerInterruptionTest{

	/** phase, root and new outdir; fixed needs {@code reference=<baseline fixture directory>}. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && Arrays.asList("mode", "phase", "root", "outdir", "reference", "case").contains(pair[0]), "Expected documented flag=value: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate fixture flag: "+pair[0]);
		}
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbmergeownerinterruptiontest.sh with assertions enabled.");
		final String phase=options.get("phase"), mode=options.containsKey("mode") ? options.get("mode") : "suite";
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		final Path root=path(options.get("root")), out=path(options.get("outdir"));
		require(Files.isDirectory(root.resolve("current")), "Missing runtime root.");
		if(mode.equals("child")){
			require(!options.containsKey("reference") && Files.isDirectory(out), "Child needs existing case directory and no reference.");
			child(out, spec(options.get("case")), phase.equals("baseline"));
			return;
		}
		require(mode.equals("suite") && !options.containsKey("case"), "Use suite or internal child mode.");
		final Path reference=options.containsKey("reference") ? path(options.get("reference")) : null;
		require(phase.equals("baseline") ? reference==null : reference!=null, "Only fixed needs a baseline fixture reference.");
		if(reference!=null){require(read(reference.resolve("summary.txt")).startsWith("BBMERGE_OWNER_INTERRUPTION_PASS phase=baseline cases=7 children=7 "), "Reference must be a successful seven-child baseline fixture directory.");}
		Files.createDirectory(out);
		write(out.resolve("cases.tsv"), "case\tentry\tsignals\tchildExit\townerInterrupted\texpectedInterrupted\tworkerFailure\n");
		int children=0, signals=0, phaseCases=0, mainCases=0, failureCases=0, childChecks=0;
		for(final Spec s : SPECS){
			final Path work=Files.createDirectory(out.resolve(s.name));
			write(work.resolve("input-r1.fa"), R1);
			write(work.resolve("input-r2.fa"), s.failure ? "" : R2);
			final List<String> command=Arrays.asList("bash", root.resolve("dev/bbmergeownerinterruptiontest.sh").toString(),
				"mode=child", "phase="+phase, "root="+root, "outdir="+work, "case="+s.name);
			write(work.resolve("child.command.tsv"), indexed(command));
			final ProcessBuilder builder=new ProcessBuilder(command).directory(work.toFile());
			builder.redirectOutput(work.resolve("child.stdout").toFile()); builder.redirectError(work.resolve("child.stderr").toFile());
			builder.environment().remove("_JAVA_OPTIONS"); builder.environment().remove("JDK_JAVA_OPTIONS");
			builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=2 -Dfile.encoding=UTF-8");
			final Process process=builder.start(); children++; process.getOutputStream().close();
			boolean ended=false;
			try{ended=process.waitFor(90, TimeUnit.SECONDS);}
			finally{if(process.isAlive()){stopOwnedChild(process, s, work);}}
			require(ended, s.name+": child exceeded its finite failure budget; inspect saved stacks/identity.");
			final int exit=process.exitValue(); write(work.resolve("child.exit"), exit+"\n");
			require(exit==0, s.name+": child exit="+exit+"; inspect child.stderr.");
			final Map<String, String> result=table(work.resolve("result.tsv"));
			final boolean expected=phase.equals("fixed") && s.signals>0;
			require(s.name.equals(result.get("case")) && phase.equals(result.get("phase")) && Boolean.toString(expected).equals(result.get("ownerInterrupted")), "Wrong child identity/phase/return status: "+s.name);
			require(read(work.resolve("child.stdout")).equals("BBMERGE_OWNER_INTERRUPTION_CHILD_PASS case="+s.name+"\n"), "Missing exact child completion marker: "+s.name);
			require(read(work.resolve("input-r1.fa")).equals(R1) && read(work.resolve("input-r2.fa")).equals(s.failure ? "" : R2), "Input file changed: "+s.name);
			if(reference!=null){
				final Path old=reference.resolve(s.name);
				for(final String file : outputFiles(s)){same(work.resolve(file), old.resolve(file), "Matched output changed: "+s.name+"/"+file);}
				for(final String file : new String[]{"input-r1.fa", "input-r2.fa"}){same(work.resolve(file), old.resolve(file), "Matched input changed: "+s.name+"/"+file);}
			}
			childChecks+=Integer.parseInt(result.get("checks")); signals+=s.signals;
			if(s.failure){failureCases++;}else if(s.main){mainCases++;}else{phaseCases++;}
			append(out.resolve("cases.tsv"), s.name+"\t"+(s.main ? "main" : "runPhase")+"\t"+s.signals+"\t"+exit+"\t"+result.get("ownerInterrupted")+"\t"+expected+"\t"+s.failure+"\n");
		}
		require(children==7 && phaseCases==3 && mainCases==2 && failureCases==2 && signals==5, "Keep seven children and five actual owner signals.");
		final String summary="BBMERGE_OWNER_INTERRUPTION_PASS phase="+phase+" cases=7 children="+children+" phaseCases="+phaseCases+
			" mainCases="+mainCases+" failureCases="+failureCases+" signals="+signals+" checks="+(checks+childChecks)+" failures=0\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
	}

	/** Uses only normal construction/main and real runPhase worker/stream ownership. */
	private static void child(final Path work, final Spec s, final boolean baseline) throws Exception{
		require(!Files.exists(Paths.get("stdin.fa")), "The stdin.fa alias must not name a real file in the child cwd.");
		require(mates(null).isEmpty(), "Fresh child already has BBMerge workers.");
		final byte[] first=Files.readAllBytes(work.resolve("input-r1.fa")), second=Files.readAllBytes(work.resolve("input-r2.fa"));
		require(Arrays.equals(first, R1.getBytes(StandardCharsets.US_ASCII)) && Arrays.equals(second, (s.failure ? "" : R2).getBytes(StandardCharsets.US_ASCII)), "Use the literal finite inputs.");
		final GateInput gate=new GateInput(first);
		final InputStream original=System.in;
		System.setIn(gate);
		final List<String> arguments=new ArrayList<String>(Arrays.asList("in1=stdin.fa", "in2="+work.resolve("input-r2.fa"),
			"t=2", "simd=f", "bf1=t", "nn=f", "hollow=f", "qtrim=f", "minlength=17", "ordered=t", "ow=t", "app=f", "parsecustom=f"));
		if(!s.failure){arguments.addAll(Arrays.asList("out="+work.resolve("merged.fa"), "outu1="+work.resolve("unmerged1.fa"), "outu2="+work.resolve("unmerged2.fa")));}
		if(s.main){arguments.addAll(Arrays.asList("ihist="+work.resolve("histogram.tsv"), "showhiststats=f", "nzo=t",
			"outa="+work.resolve("adapters.fa"), "outc="+work.resolve("cardinality.txt"), "loglogtype=LogLog", "outstream="+work.resolve("diagnostics.txt")));}
		write(work.resolve("constructor.command.tsv"), indexed(arguments));
		final Owner owner=new Owner(s, arguments.toArray(new String[0]));
		owner.setDaemon(true);// Real workers inherit this fixture owner's daemon status.
		Throwable problem=null;
		Capture captured=null;
		int sent=0;
		try{
			if(!s.main){owner.merge=new BBMerge(owner.arguments);}
			owner.start();
			await(gate.atEof, "real input reader reached held EOF");
			captured=capture(owner, s);
			owner.captured=captured;
			write(work.resolve("identities.tsv"), captured.identities);
			require(!owner.done() && gate.release.getCount()==1, "Owner returned before controlled EOF release.");
			final Thread target=waitingOn(owner, captured.workers);
			write(work.resolve("join-before.txt"), info(owner));
			for(int i=0; i<s.signals; i++){
				require(target.isAlive() && !owner.done() && gate.release.getCount()==1, "A real worker must remain gated before each signal.");
				owner.interrupt(); sent++;
				consumed(owner, target);
				write(work.resolve("join-consumed-"+sent+".txt"), info(owner));
				require(!owner.done() && target.isAlive() && gate.release.getCount()==1, "Owner escaped or worker terminated while input was held.");
			}
		}catch(Throwable t){problem=t;}
		finally{
			gate.release.countDown();
			try{join(owner, "owner");}catch(Throwable t){if(problem==null){problem=t;}else{problem.addSuppressed(t);}}
			try{
				final List<Thread> owned=captured==null ? mates(owner.merge) : captured.workers;
				for(final Thread worker : owned){try{join(worker, "owned MateThread");}catch(Throwable t){if(problem==null){problem=t;}else{problem.addSuppressed(t);}}}
			}catch(Throwable t){if(problem==null){problem=t;}else{problem.addSuppressed(t);}}
			System.setIn(original);
		}
		if(problem!=null){throw new AssertionError("Finite child setup/gate/cleanup failed: "+s.name, problem);}
		require(captured!=null && owner.done() && sent==s.signals, "Missing complete owned lifecycle.");
		require(owner.snapshotFailure==null, "Return-boundary snapshot failed: "+owner.snapshotFailure);
		write(work.resolve("return-states.tsv"), owner.returnStates);
		final Throwable cause=owner.thrown==null ? null : owner.thrown.getCause();
		final StringBuilder failureEvidence=new StringBuilder("role\ttype\tmessage\tidentity\tisOwnerCause\n");
		failureEvidence.append(failureRow("owner", owner.thrown, cause)).append(failureRow("owner-cause", cause, cause));
		int workerFailures=0; boolean originalCause=false;
		for(final Thread worker : captured.workers){final Throwable failure=(Throwable)value(worker, "failure"); if(failure!=null){workerFailures++;} originalCause|=cause!=null && failure==cause; failureEvidence.append(failureRow("worker-"+worker.getId(), failure, cause));}
		write(work.resolve("failures.tsv"), failureEvidence.toString());
		require(owner.interrupted==(!baseline && sent>0), "Owner interruption was not preserved at its exact return/throw boundary.");
		require(owner.workersTerminated && owner.inputsClosed && owner.writersFinished, "Selected actual workers/readers/writers were not complete at owner return.");
		require(owner.pairs==(s.failure ? 0 : 1) && owner.bases==(s.failure ? 0 : 32), "Literal processed pair/base counters differ: "+owner.pairs+"/"+owner.bases);
		require(owner.correct==0 && owner.incorrect==0 && owner.ambiguous==0 && owner.shortCount==0 && owner.noSolution==0, "The minlength guard must bypass overlap work without assigning another outcome counter.");
		require(owner.error==s.failure && gate.release.getCount()==0 && !gate.readInterrupted, "Unexpected input/phase error or reader interruption.");
		if(s.failure){
			require(owner.thrown!=null && owner.thrown.getClass()==RuntimeException.class && "BBMerge worker thread failed.".equals(owner.thrown.getMessage()), "Preserve the existing worker-failure wrapper.");
			require(cause instanceof AssertionError && cause.getMessage().startsWith("Paired files have different read counts (unpaired remnant):"), "Expected the actual paired-count assertion as cause: "+cause);
			require(workerFailures==2 && originalCause, "Cause must be an original captured worker failure, with both workers terminated.");
		}else{
			require(owner.thrown==null, "Unexpected owner failure: "+owner.thrown);
			for(final Thread worker : captured.workers){require(value(worker, "failure")==null, "Successful control has a worker failure.");}
			require(read(work.resolve("merged.fa")).isEmpty(), "The minlength control must not merge a read.");
			require(read(work.resolve("unmerged1.fa")).equals(R1) && read(work.resolve("unmerged2.fa")).equals(R2), "Expected the exact original unmerged FASTAs.");
			for(int i=0; i<captured.writers.size(); i++){final Writer writer=captured.writers.get(i); require(writer.readsWritten()==(i==0 ? 0 : 1) && writer.basesWritten()==(i==0 ? 0 : 16), "Selected writer counters differ from the literal empty-merged/one-read-per-mate outputs.");}
		}
		if(s.main){
			require(owner.mainFlagsExpected && owner.diagnosticsClosed, "Main flag postconditions or selected plain diagnostic closure differ at return.");
			require(!Thread.currentThread().isInterrupted(), "Histogram observation must run on the un-interrupted fixture observer.");
			ReadWrite.waitForWritingToFinish();// Counted write completion, deliberately AFTER owner return; not a join.
			require(ReadWrite.countActiveThreads()==0, "Counted histogram writes remain active after observer drain.");
			require(read(work.resolve("histogram.tsv")).equals("#InsertSize\tCount\n"), "Expected the literal empty histogram header.");
			require(read(work.resolve("adapters.fa")).equals(">Read1_adapter\nN\n>Read2_adapter\nN\n"), "No accepted overlap means two unsupported adapter consensuses.");
			require(read(work.resolve("cardinality.txt")).equals("0\n"), "Sixteen-base reads supply no31-mers to the selected LogLog tracker.");
			final String diagnostics=read(work.resolve("diagnostics.txt"));
			require(diagnostics.matches("(?s).*\\nPairs:[ \\t]+1\\R.*") && diagnostics.contains("Adapters counted: 0"), "Main diagnostics omit the literal pair/adapter counters.");
		}
		require(Arrays.equals(first, gate.data) && Arrays.equals(first, Files.readAllBytes(work.resolve("input-r1.fa"))) && Arrays.equals(second, Files.readAllBytes(work.resolve("input-r2.fa"))), "Input payload or files changed.");
		final String result="case\t"+s.name+"\nphase\t"+(baseline ? "baseline" : "fixed")+"\nsignals\t"+sent+
			"\nownerInterrupted\t"+owner.interrupted+"\nworkersTerminatedAtReturn\t"+owner.workersTerminated+"\ninputsClosedAtReturn\t"+owner.inputsClosed+
			"\nwritersFinishedAtReturn\t"+owner.writersFinished+"\nworkerFailure\t"+(owner.thrown!=null)+"\nworkerCount\t"+captured.workers.size()+"\nreaderCount\t"+captured.readers.length+"\nwriterCount\t"+captured.writers.size()+"\npairs\t"+owner.pairs+"\nbases\t"+owner.bases+
			"\nmainFlagsExpected\t"+owner.mainFlagsExpected+"\ndiagnosticsClosed\t"+owner.diagnosticsClosed+"\nhistogramDrainAfterOwner\t"+s.main+"\nchecks\t"+checks+"\n";
		write(work.resolve("result.tsv"), result);
		System.out.print("BBMERGE_OWNER_INTERRUPTION_CHILD_PASS case="+s.name+"\n");
	}

	/** Captures live production identities without changing a production field. */
	private static Capture capture(final Owner owner, final Spec s) throws Exception{
		final long end=System.nanoTime()+BUDGET;
		List<Thread> workers=null;
		while(System.nanoTime()<end){
			workers=mates(owner.merge);
			if(workers.size()==2){break;}
			if(owner.done()){throw new AssertionError("Owner returned before its two actual workers were captured.");}
			LockSupport.parkNanos(PARK);
		}
		require(workers!=null && workers.size()==2, "Expected exactly two actual MateThreads.");
		final BBMerge merge=(BBMerge)value(workers.get(0), "this$0");
		if(owner.merge==null){owner.merge=merge;}
		require(owner.merge==merge && !Shared.SIMD && Shared.threads()==2 && !(Boolean)value(merge, "hollow"), "Unexpected merger/global selection.");
		final Object input=value(workers.get(0), "cris");
		require(input.getClass().getName().equals("stream.PairStreamer"), "Require actual paired input.");
		final StringBuilder identities=new StringBuilder("role\tclass\tidentity\n");
		identities.append("input\t").append(input.getClass().getName()).append('\t').append(System.identityHashCode(input)).append('\n');
		final ByteFile1[] readers=new ByteFile1[2];
		for(int i=0; i<2; i++){
			final Object child=value(input, i==0 ? "s1" : "s2");
			require(child.getClass().getName().equals("stream.FastaStreamerZT"), "Require selected synchronous FASTA reader.");
			final Object bf=value(child, "bf");
			require(bf!=null && bf.getClass()==ByteFile1.class && ((ByteFile1)bf).isOpen(), "Require a captured open ByteFile1.");
			readers[i]=(ByteFile1)bf;
			identities.append("fasta-"+i).append('\t').append(child.getClass().getName()).append('\t').append(System.identityHashCode(child)).append('\n');
		}
		final List<Writer> writers=new ArrayList<Writer>();
		for(final Thread worker : workers){
			require(worker.isAlive() && worker.isDaemon() && value(worker, "this$0")==merge && value(worker, "cris")==input, "Wrong worker identity or ownership.");
			identities.append("worker-").append(worker.getId()).append('\t').append(worker.getClass().getName()).append('\t').append(System.identityHashCode(worker)).append('\n');
		}
		final Writer good=(Writer)value(workers.get(0), "rosgood"), bad=(Writer)value(workers.get(0), "rosbad");
		require(value(workers.get(0), "rosi")==null && value(workers.get(1), "rosgood")==good && value(workers.get(1), "rosbad")==bad && value(workers.get(1), "rosi")==null, "Unexpected output ownership.");
		if(s.failure){require(good==null && bad==null, "Failure controls must have no production output queues.");}
		else{
			require(good!=null && bad!=null && bad.getClass().getName().equals("stream.PairedWriter"), "Require merged and paired-unmerged outputs.");
			writers.add(good); writers.add((Writer)value(bad, "w1")); writers.add((Writer)value(bad, "w2"));
			for(final Writer writer : writers){require(writer.getClass().getName().equals("stream.FastqWriterST2") && !(Boolean)value(writer, "threaded") && value(writer, "writerThread")==null, "Require actual unthreaded FastqWriterST2.");}
		}
		for(final ByteFile1 bf : readers){identities.append("reader\t").append(bf.getClass().getName()).append('\t').append(System.identityHashCode(bf)).append('\n');}
		for(final Writer writer : writers){identities.append("writer\t").append(writer.getClass().getName()).append('\t').append(System.identityHashCode(writer)).append('\n');}
		if(s.main){require(openDescriptors(owner.diagnostics())>0, "Redirected diagnostics must be open before main cleanup.");}
		return new Capture(workers, readers, writers, (Streamer)input, identities.toString());
	}

	/** Measures owner status before any observation that might perform fixture I/O. */
	private static final class Owner extends Thread{
		Owner(final Spec spec_, final String[] args_){super("fixture-owner-"+spec_.name); spec=spec_; arguments=args_;}
		@Override
		public void run(){
			try{if(spec.main){BBMerge.main(arguments);}else{merge.runPhase(true, -1, false);}}
			catch(Throwable t){thrown=t;}
			finally{
				interrupted=isInterrupted();
				try{
					final Capture c=captured;
					require(c!=null, "Missing live captures at owner boundary.");
					final StringBuilder states=new StringBuilder("role\tstate\treads\tbases\n");
					workersTerminated=true; for(final Thread worker : c.workers){final State state=worker.getState(); workersTerminated&=state==State.TERMINATED; states.append("worker-").append(worker.getId()).append('\t').append(state).append("\tNA\tNA\n");}
					inputsClosed=true; for(final ByteFile1 reader : c.readers){final boolean open=reader.isOpen(); inputsClosed&=!open; states.append("reader-open\t").append(open).append("\tNA\tNA\n");}
					writersFinished=true; for(final Writer writer : c.writers){final boolean finished=writer.finishedSuccessfully() && !writer.errorState(); writersFinished&=finished; states.append("writer-finished\t").append(finished).append('\t').append(writer.readsWritten()).append('\t').append(writer.basesWritten()).append('\n');}
					returnStates=states.toString();
					pairs=number(merge, "readsProcessedTotal"); bases=number(merge, "basesProcessedTotal"); correct=number(merge, "correctCountTotal");
					incorrect=number(merge, "incorrectCountTotal"); ambiguous=number(merge, "ambiguousCountTotal"); shortCount=number(merge, "tooShortCountTotal"); noSolution=number(merge, "noSolutionCountTotal");
					error=(Boolean)value(merge, "errorState");
					mainFlagsExpected=!(Boolean)value(BBMerge.class, "MAKE_VECTOR") && !(Boolean)value(BBMergeOverlapper.class, "MAKE_VECTOR") && Read.VALIDATE_IN_CONSTRUCTOR;
					diagnosticsClosed=!spec.main || openDescriptors(diagnostics())==0;
				}catch(Throwable t){snapshotFailure=t;}
				done.countDown();
			}
		}
		Path diagnostics(){for(final String arg : arguments){if(arg.startsWith("outstream=")){return Paths.get(arg.substring(10));}} throw new AssertionError("No redirected diagnostics.");}
		boolean done(){return done.getCount()==0;}
		final Spec spec; final String[] arguments; final CountDownLatch done=new CountDownLatch(1);
		volatile BBMerge merge; volatile Capture captured; volatile Throwable thrown, snapshotFailure;
		volatile String returnStates;
		volatile boolean interrupted, workersTerminated, inputsClosed, writersFinished, error, mainFlagsExpected, diagnosticsClosed;
		volatile long pairs, bases, correct, incorrect, ambiguous, shortCount, noSolution;
	}

	/** Ordinary finite standard input with a latch-controlled terminal EOF. */
	private static final class GateInput extends InputStream{
		GateInput(final byte[] source){data=source.clone();}
		@Override
		public int read() throws IOException{final byte[] one=new byte[1]; return read(one, 0, 1)<0 ? -1 : one[0]&255;}
		@Override
		public int read(final byte[] b, final int off, final int len) throws IOException{
			if(len==0){return 0;}
			if(position<data.length){final int n=Math.min(len, data.length-position); System.arraycopy(data, position, b, off, n); position+=n; return n;}
			atEof.countDown();
			try{release.await();}catch(InterruptedException e){readInterrupted=true; Thread.currentThread().interrupt(); throw new InterruptedIOException("Fixture reader was interrupted instead of its owner.");}
			return -1;
		}
		final byte[] data; int position; volatile boolean readInterrupted;
		final CountDownLatch atEof=new CountDownLatch(1), release=new CountDownLatch(1);
	}

	private static final class Capture{
		Capture(final List<Thread> w, final ByteFile1[] r, final List<Writer> o, final Streamer i, final String text){workers=w; readers=r; writers=o; input=i; identities=text;}
		final List<Thread> workers; final ByteFile1[] readers; final List<Writer> writers; final Streamer input; final String identities;
	}
	private static final class Spec{
		Spec(final String n, final boolean m, final boolean f, final int s){name=n; main=m; failure=f; signals=s;}
		final String name; final boolean main, failure; final int signals;
	}

	/** Finds only this child's actual MateThreads, then pins their synthetic outer owner. */
	private static List<Thread> mates(final BBMerge merge) throws Exception{
		final List<Thread> list=new ArrayList<Thread>();
		for(final Thread thread : Thread.getAllStackTraces().keySet()){
			if(thread.getClass().getName().equals("jgi.BBMerge$MateThread") && (merge==null || value(thread, "this$0")==merge)){list.add(thread);}
		}
		return list;
	}
	private static Thread waitingOn(final Thread owner, final List<Thread> workers) throws Exception{
		final long end=System.nanoTime()+BUDGET;
		while(System.nanoTime()<end){for(final Thread worker : workers){if(joinWait(owner, worker)){return worker;}} if(!owner.isAlive()){throw new AssertionError("Owner terminated before its live join.");} LockSupport.parkNanos(PARK);}
		throw new AssertionError("Owner did not enter a join on a captured worker: "+info(owner));
	}
	private static boolean joinWait(final Thread owner, final Thread worker){
		final ThreadInfo info=BEAN.getThreadInfo(owner.getId(), 30); if(info==null){return false;}
		final LockInfo lock=info.getLockInfo(); boolean join=false, phase=false;
		for(final StackTraceElement frame : info.getStackTrace()){join|=frame.getClassName().equals("java.lang.Thread") && frame.getMethodName().equals("join"); phase|=frame.getClassName().equals("jgi.BBMerge") && frame.getMethodName().equals("runPhase");}
		return info.getThreadState()==Thread.State.TIMED_WAITING && lock!=null && lock.getIdentityHashCode()==System.identityHashCode(worker) && join && phase;
	}
	private static void consumed(final Thread owner, final Thread worker) throws Exception{
		final long end=System.nanoTime()+BUDGET;
		while(System.nanoTime()<end){if(!owner.isAlive()){throw new AssertionError("Owner exited before interrupt consumption.");} if(!owner.isInterrupted() && joinWait(owner, worker)){return;} LockSupport.parkNanos(PARK);}
		throw new AssertionError("Owner did not consume the signal and rejoin: "+info(owner));
	}
	private static String info(final Thread owner){final ThreadInfo info=BEAN.getThreadInfo(owner.getId(), 40); return info==null ? "terminated\n" : info.toString();}
	private static void await(final CountDownLatch latch, final String what) throws Exception{require(latch.await(BUDGET, TimeUnit.NANOSECONDS), "Timed out: "+what);}
	private static void join(final Thread thread, final String what) throws Exception{final long end=System.nanoTime()+BUDGET; while(thread.isAlive() && System.nanoTime()<end){thread.join(20);} require(!thread.isAlive(), "Timed out joining "+what);}

	/** Uses an exact own Process handle, verified current command and parent/session lineage. */
	private static void stopOwnedChild(final Process process, final Spec s, final Path work) throws Exception{
		if(!process.isAlive()){return;}
		long pid;
		try{pid=((Number)Process.class.getMethod("pid").invoke(process)).longValue();}
		catch(NoSuchMethodException e){final Field f=process.getClass().getDeclaredField("pid"); f.setAccessible(true); pid=((Number)f.get(process)).longValue();}
		final String[] childStat=procStat(Paths.get("/proc/"+pid+"/stat")), selfStat=procStat(Paths.get("/proc/self/stat"));
		final String self=ManagementFactory.getRuntimeMXBean().getName().split("@")[0];
		final String command=new String(Files.readAllBytes(Paths.get("/proc/"+pid+"/cmdline")), StandardCharsets.UTF_8).replace('\0', ' ');
		require(childStat[1].equals(self) && childStat[3].equals(selfStat[3]) && command.contains(BBMergeOwnerInterruptionTest.class.getName()) && command.contains("mode=child") && command.contains("case="+s.name), "Cannot establish current child command/parent/session; no signal sent.");
		write(work.resolve("timeout-identity.txt"), "pid="+pid+"\nparent="+childStat[1]+"\nsession="+childStat[3]+"\ncommand="+command+"\n");
		if(process.isAlive()){process.destroyForcibly(); require(process.waitFor(10, TimeUnit.SECONDS), "Owned timed-out child did not terminate.");}
	}
	private static String[] procStat(final Path path) throws IOException{final String text=read(path); return text.substring(text.lastIndexOf(')')+2).trim().split(" +");}
	/** Reads only this process's descriptor links; does not open the target file. */
	private static int openDescriptors(final Path target) throws IOException{
		int count=0;
		final Path canonical=target.toRealPath();
		try(DirectoryStream<Path> files=Files.newDirectoryStream(Paths.get("/proc/self/fd"))){for(final Path file : files){try{final Path link=Files.readSymbolicLink(file); if(link.isAbsolute() && link.toRealPath().equals(canonical)){count++;}}catch(java.nio.file.NoSuchFileException e){/* descriptor closed between list and link read */}}}
		return count;
	}
	private static Spec spec(final String name){for(final Spec s : SPECS){if(s.name.equals(name)){return s;}} throw new AssertionError("Unknown child case: "+name);}
	private static String[] outputFiles(final Spec s){return s.failure ? new String[0] : s.main ? new String[]{"merged.fa", "unmerged1.fa", "unmerged2.fa", "histogram.tsv", "adapters.fa", "cardinality.txt"} : new String[]{"merged.fa", "unmerged1.fa", "unmerged2.fa"};}
	private static String failureRow(final String role, final Throwable failure, final Throwable cause){return role+"\t"+(failure==null ? "NONE\t" : failure.getClass().getName()+"\t"+String.valueOf(failure.getMessage()).replace("\t", "\\t").replace("\n", "\\n"))+"\t"+System.identityHashCode(failure)+"\t"+(cause!=null && failure==cause)+"\n";}
	private static Map<String, String> table(final Path path) throws IOException{final Map<String, String> map=new LinkedHashMap<String, String>(); for(final String row : read(path).split("\n")){final String[] pair=row.split("\t", 2); require(pair.length==2 && map.put(pair[0], pair[1])==null, "Invalid child result row: "+row);} return map;}
	private static Object value(final Object owner, final String name) throws Exception{final Class<?> type=owner instanceof Class<?> ? (Class<?>)owner : owner.getClass(); final Field f=type.getDeclaredField(name); f.setAccessible(true); return f.get(owner instanceof Class<?> ? null : owner);}
	private static long number(final Object owner, final String name) throws Exception{return ((Number)value(owner, name)).longValue();}
	private static String indexed(final List<String> args){final StringBuilder b=new StringBuilder(); for(int i=0; i<args.size(); i++){b.append(i).append('\t').append(args.get(i)).append('\n');} return b.toString();}
	private static void same(final Path a, final Path b, final String message) throws IOException{require(Arrays.equals(Files.readAllBytes(a), Files.readAllBytes(b)), message);}
	private static Path path(final String s){require(s!=null && !s.isEmpty(), "Missing required path."); return Paths.get(s).toAbsolutePath().normalize();}
	private static void require(final boolean condition, final String message){checks++; if(!condition){throw new AssertionError(message);}}
	private static void write(final Path path, final String s) throws IOException{Files.write(path, s.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path path, final String s) throws IOException{Files.write(path, s.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);}
	private static String read(final Path path) throws IOException{return new String(Files.readAllBytes(path), StandardCharsets.UTF_8);}

	private static int checks;
	private static final ThreadMXBean BEAN=ManagementFactory.getThreadMXBean();
	private static final long BUDGET=TimeUnit.SECONDS.toNanos(10), PARK=TimeUnit.MILLISECONDS.toNanos(1);
	private static final String R1=">pair/1\nACGTTGCAACGTTGCA\n", R2=">pair/2\nTGCAACGTTGCAACGT\n";
	private static final Spec[] SPECS={new Spec("phase-zero", false, false, 0), new Spec("phase-one", false, false, 1),
		new Spec("phase-two", false, false, 2), new Spec("main-zero", true, false, 0), new Spec("main-one", true, false, 1),
		new Spec("failure-zero", false, true, 0), new Spec("failure-one", false, true, 1)};
}
