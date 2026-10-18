package test.prok;

import prok.CallGenes;
import prok.GeneModel;

import java.lang.reflect.Constructor;
import java.util.Collections;
import java.util.concurrent.CountDownLatch;

import stream.ConcurrentCollectionReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ListNum;

/** Subprocess-only regression: the actual training worker must halt on an escaping failure.
 * A live non-daemon consumer makes a swallowed worker exception fail by wrapper timeout.
 * No production fault flag, HBM input, or modification of the biological runtime is needed.
 * @author Sayu
 */
public final class ProteinTrainingFailureTest {

	public static void main(final String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Expected one readable tiny FASTA fixture");}
		final CallGenes owner=new CallGenes(new String[]{"in="+args[0], "out=null", "annotation=off",
			"cdsonly=t", "tax=f", "t=1"});
		final GeneModel base=new GeneModel(false);
		base.statsCDS.setInner(2, 3);
		final ConcurrentReadInputStream input=new ConcurrentCollectionReadInputStream(Collections.<Read>emptyList(), null, -1){
			@Override
			public ListNum<Read> nextList(){throw new IllegalStateException("PROTEIN_TRAINING_WORKER_SENTINEL");}
		};
		final Class<?> type=Class.forName("prok.CallGenes$ProteinTrainingThread");
		final Constructor<?> constructor=type.getDeclaredConstructor(CallGenes.class,
			ConcurrentReadInputStream.class, GeneModel.class, String.class);
		constructor.setAccessible(true);
		final Thread worker=(Thread)constructor.newInstance(owner, input, base, args[0]);
		final CountDownLatch ready=new CountDownLatch(1), blocked=new CountDownLatch(1);
		final Thread consumer=new Thread(()->{
			ready.countDown();
			while(true){
				try{blocked.await();}catch(InterruptedException ignored){/* A swallowed failure must not release this consumer. */}
			}
		}, "training-failure-live-consumer");
		consumer.setDaemon(false); consumer.start(); ready.await();
		System.err.println("PROTEIN_TRAINING_FAILURE_ARMED consumer_daemon="+consumer.isDaemon());
		worker.start(); worker.join();
		throw new AssertionError("Training worker returned without halting the live-consumer subprocess");
	}
}
