package test.prok;

import java.lang.management.ManagementFactory;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.util.Collections;
import java.util.Locale;
import java.util.concurrent.CountDownLatch;

import fileIO.ByteStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.CallGenes;
import prok.GeneModel;
import prok.NcrnaStageDiagSink;
import prok.RefinementAttemptSink;
import shared.KillSwitch;
import stream.ConcurrentCollectionReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ListNum;

/** Subprocess-only test of real heap exhaustion on the ordinary caller worker.
 * Reuses ProteinTrainingFailureTest's native collection-stream/worker fixture.
 * No production fault flag or replacement implementation of the handler.
 * Expected result: original worker OOM,memory guidance,whole-JVM exit1.
 * @author Brian Bushnell, Raiden
 */
public final class CallGenesWorkerOomTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, CallGenesWorkerOomTest.class, false);
		final Parser parser=new Parser();
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			check(parser.parse(arg, a, b), "Unknown worker OOM fixture option: "+arg);
		}
		check(parser.in1!=null && parser.in2==null && parser.out1==null && parser.out2==null,
			"Use in=one_readable_tiny_FASTA without output or paired input");
		boolean assertions=false;
		assert(assertions=true);
		check(assertions && Runtime.getRuntime().maxMemory()<=32L*1024*1024,
			"Use assertions and -Xmx32m so the worker's valid64MiB array must exhaust the heap");
		for(final String arg:ManagementFactory.getRuntimeMXBean().getInputArguments()){
			check(!arg.equals("-XX:+ExitOnOutOfMemoryError") && !arg.equals("-XX:+CrashOnOutOfMemoryError"),
				"The worker catch must handle the real allocation failure; omit JVM OOM-exit flags");
		}
		final CallGenes owner=new CallGenes(new String[]{"in="+parser.in1, "out=null",
			"cdsonly=t", "trnaalign=f", "trnaboundarynet=f", "tax=f", "t=1"});
		final GeneModel model=new GeneModel(false);
		final Class<?> type=Class.forName("prok.CallGenes$ProcessThread");
		// The native constructor is private; reflect only in this subprocess test,
		// rather than exposing a production fault-injection API or starting file workers.
		final Constructor<?> constructor=type.getDeclaredConstructor(CallGenes.class,
			ConcurrentReadInputStream.class, ByteStreamWriter.class, ConcurrentReadOutputStream.class,
			ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class,
			GeneModel.class, int.class, int.class, RefinementAttemptSink.class, NcrnaStageDiagSink.class,
			Class.forName("prok.NcrnaVoteDiagSink"), long.class);
		constructor.setAccessible(true);
		final Thread worker=(Thread)constructor.newInstance(owner, new OomInput(), null,
			null, null, null, null, model, 60, 0, null, null, null, 0L);
		check(unset(type, worker, "cost") && unset(type, worker, "fileOutput"),
			"CallGenes.runInner must take the ordinary cost-null/fileOutput-null path under test");
		final CountDownLatch ready=new CountDownLatch(1), blocked=new CountDownLatch(1);
		final Thread consumer=new Thread(()->{
			ready.countDown();
			while(true){
				try{blocked.await();}
				catch(InterruptedException ignored){/* An escaping worker OOM must not release this consumer. */}
			}
		}, "ordinary-worker-oom-live-consumer");
		consumer.setDaemon(false);
		KillSwitch.addBallast();
		consumer.start();
		ready.await();
		System.err.println("CALLGENES_ORDINARY_WORKER_OOM_ARMED_WITH_LIVE_CONSUMER");
		worker.start();
		worker.join();
		throw new AssertionError("Ordinary worker returned without halting the live-consumer subprocess");
	}

	/** Asserts the actual native worker route,not a helper-only or file-worker test. */
	private static boolean unset(final Class<?> type, final Object worker, final String name) throws Exception{
		assert(type!=null && worker!=null && name!=null) : "The native worker constructor supplies both lifecycle fields";
		final Field field=type.getDeclaredField(name);
		field.setAccessible(true);
		return field.get(worker)==null;
	}

	private static void check(final boolean condition, final String why){
		if(!condition){throw new AssertionError(why);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	private static final class OomInput extends ConcurrentCollectionReadInputStream{
		OomInput(){super(Collections.<Read>emptyList(), null, -1);}
		@Override
		public ListNum<Read> nextList(){
			System.err.println("CALLGENES_ORDINARY_WORKER_REQUESTED_64_MIB");
			allocation=new byte[64*1024*1024];
			throw new AssertionError("The64MiB allocation unexpectedly fit the at-most32MiB heap");
		}
		private byte[] allocation;
	}
}
