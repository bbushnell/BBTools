package test.consensus;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.charset.StandardCharsets;
import java.util.LinkedHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicReference;

import consensus.Lilypad;
import stream.SamLine;

/** Exact counter checks for Lilypad mixed-pair edge updates. @author Jean */
public final class LilypadMixedPairEdgeTest{

	/** Builds real Lilypad contigs and checks reciprocal mixed-pair edge counters. */
	public static void main(String[] args) throws Exception{
		runScenario(false);
		runScenario(true);
		System.out.println("LILYPAD_MIXED_PAIR_EDGE_PASS");
	}

	/** Runs either the immutable sequential control or a four-thread shared-contig accumulation. */
	private static void runScenario(final boolean concurrent) throws Exception{
		final Object lilypad=makeLilypad();
		final LinkedHashMap<?, ?> refMap=(LinkedHashMap<?, ?>)field(lilypad, "refMap");
		final Object left=refMap.get("left");
		final Object right=refMap.get("right");
		set(lilypad, "maxPairDist", 3000);
		set(lilypad, "sameStrandPairs", false);
		if(concurrent){accumulateConcurrent(left, right);}
		else{accumulateSequential(left, right, 0, 1024);}
		check(edge(left, "right", false), right, false, 1024, 61440, 33280);
		check(edge(right, "left", true), left, true, 1024, 61440, 33280);
	}

	private static void accumulateSequential(final Object left, final Object right,
			final int start, final int count){
		for(int i=start, lim=start+count; i<lim; i++){
			handle(left, sam(99, "left", 20, "right", (i&7)+1));
			handle(right, sam(83, "right", 8, "left", 20+(i&7)));
		}
	}

	/** Starts four workers together against the same two real Contigs and propagates failures. */
	private static void accumulateConcurrent(final Object left, final Object right) throws Exception{
		final CountDownLatch start=new CountDownLatch(1);
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final Thread[] threads=new Thread[4];
		for(int tid=0; tid<threads.length; tid++){
			final int offset=tid*256;
			threads[tid]=new Thread(new Runnable(){
				@Override
				public void run(){
					try{
						start.await();
						accumulateSequential(left, right, offset, 256);
					}catch(Throwable t){failure.compareAndSet(null, t);}
				}
			}, "LilypadMixedPairEdgeTest-"+tid);
			threads[tid].start();
		}
		start.countDown();
		for(final Thread thread : threads){thread.join();}
		final Throwable thrown=failure.get();
		if(thrown!=null){throw failure(thrown);}
	}

	/** Invokes the real mixed-pair method. */
	private static void handle(final Object contig, final SamLine line){
		call(contig, "handleMixedPair", new Class<?>[]{SamLine.class}, line);
	}

	/** Reads one actual Edge out of the selected Contig map. */
	private static Object edge(final Object contig, final String name, final boolean left){
		final LinkedHashMap<?, ?> map=(LinkedHashMap<?, ?>)field(contig, left ? "leftEdgeMap" : "rightEdgeMap");
		final Object edge=map.get(name);
		if(edge==null){throw new AssertionError("Missing "+(left ? "left" : "right")+" edge to "+name);}
		return edge;
	}

	/** Checks exact Edge counter fields after deterministic mixed-pair accumulation. */
	private static void check(final Object edge, final Object dest, final boolean left,
			final long count, final long weight, final long distance){
		final long same=(Long)field(edge, "sameStrandCount");
		final long dif=(Long)field(edge, "difStrandCount");
		if(field(edge, "dest")!=dest || (Boolean)field(edge, "leftEdge")!=left || same!=count || dif!=0 ||
				(Long)field(edge, "weight")!=weight || (Long)field(edge, "distanceSum")!=distance ||
				(Long)field(edge, "badCount")!=0){
			throw new AssertionError("Unexpected edge counters: dest="+field(edge, "dest")+
					" left="+field(edge, "leftEdge")+" same="+same+" dif="+dif+
					" weight="+field(edge, "weight")+" distance="+field(edge, "distanceSum")+
					" bad="+field(edge, "badCount"));
		}
	}

	/** Creates the SAM fields used by Lilypad and Edge.add. */
	private static SamLine sam(final int flag, final String rname, final int pos, final String rnext, final int pnext){
		final SamLine line=new SamLine();
		line.flag=flag;
		line.pos=pos;
		line.mapq=60;
		line.setCigar("10M");
		line.pnext=pnext;
		line.setRnameS(rname);
		line.setRnextS(rnext);
		line.setSeq("AAAAAAAAAA".getBytes(StandardCharsets.US_ASCII));
		return line;
	}

	private static Object call(final Object target, final String name, final Class<?>[] parameters, final Object... args){
		try{
			final Method method=target.getClass().getDeclaredMethod(name, parameters);
			method.setAccessible(true);
			return method.invoke(target, args);
		}catch(InvocationTargetException e){throw failure(e.getCause());}
		catch(ReflectiveOperationException e){throw new AssertionError("Cannot access "+name, e);}
	}

	private static Object field(final Object target, final String name){
		try{
			final Field field=target.getClass().getDeclaredField(name);
			field.setAccessible(true);
			return field.get(target);
		}catch(ReflectiveOperationException e){throw new AssertionError("Cannot read "+name, e);}
	}

	private static void set(final Object target, final String name, final Object value){
		try{
			final Field field=target.getClass().getDeclaredField(name);
			field.setAccessible(true);
			field.set(target, value);
		}catch(ReflectiveOperationException e){throw new AssertionError("Cannot set "+name, e);}
	}

	private static RuntimeException failure(final Throwable cause){
		if(cause instanceof Error){throw (Error)cause;}
		if(cause instanceof RuntimeException){return (RuntimeException)cause;}
		throw new AssertionError("Unexpected checked failure from Lilypad mixed-pair helper", cause);
	}

	private static Object makeLilypad() throws Exception{
		final Path dir=Files.createTempDirectory("lilypad-edge");
		final Path ref=dir.resolve("ref.fa");
		final Path sam=dir.resolve("empty.sam");
		final Path out=dir.resolve("out.fa");
		write(ref, ">left\n"+repeat('A', 40)+"\n>right\n"+repeat('C', 40)+"\n");
		write(sam, "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:left\tLN:40\n@SQ\tSN:right\tLN:40\n");
		final Lilypad lilypad=new Lilypad(new String[]{"in="+sam, "ref="+ref, "out="+out, "ow=t", "append=f"});
		call(lilypad, "loadReferenceCustom", new Class<?>[0]);
		return lilypad;
	}

	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII));
	}

	private static String repeat(final char c, final int length){
		final char[] chars=new char[length];
		for(int i=0; i<length; i++){chars[i]=c;}
		return new String(chars);
	}
}
