package test.structures;

import structures.LongHeap;
import structures.LongHeapMap;
import structures.LongHeapSet;

/** Baseline and fixed oracles for zero-capacity LongHeap retainers.
 * @author Jean
 */
public final class LongHeapZeroCapacityCase{

	/** Runs finite zero-capacity checks. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		checkAssertionsEnabled();

		longHeapZeroPositive(expect);
		longHeapZeroZero(expect);
		longHeapOnePositive(expect);
		longHeapSetZeroPositive(expect);
		longHeapSetOnePositive(expect);
		longHeapMapZeroPositive(expect);
		longHeapMapOnePositive(expect);
	}

	/** Exposes zero-capacity heap construction or insertion assertion. */
	private static void longHeapZeroPositive(final String expect){
		if("fixed".equals(expect)){
			expectAssertion("longheap-zero-positive", new Runnable(){
				@Override
				public void run(){new LongHeap(0);}
			});
			pass("longheap-zero-positive", "ConstructorAssertionError", expect);
			return;
		}
		final LongHeap heap=new LongHeap(0);
		expectAssertion("longheap-zero-positive", new Runnable(){
			@Override
			public void run(){heap.add(1);}
		});
		checkEquals("longheap-zero-positive-size", 0, heap.size());
		checkEquals("longheap-zero-positive-capacity", 0, heap.capacity());
		checkEquals("longheap-zero-positive-peek", LongHeap.EMPTY, heap.peek());
		pass("longheap-zero-positive", "AssertionError", expect);
	}

	/** Checks zero against the zero-capacity heap constructor or default root value. */
	private static void longHeapZeroZero(final String expect){
		if("fixed".equals(expect)){
			expectAssertion("longheap-zero-zero", new Runnable(){
				@Override
				public void run(){new LongHeap(0);}
			});
			pass("longheap-zero-zero", "ConstructorAssertionError", expect);
			return;
		}
		final LongHeap heap=new LongHeap(0);
		checkEquals("longheap-zero-zero-add", false, heap.add(0));
		checkEquals("longheap-zero-zero-size", 0, heap.size());
		checkEquals("longheap-zero-zero-capacity", 0, heap.capacity());
		checkEquals("longheap-zero-zero-peek", LongHeap.EMPTY, heap.peek());
		pass("longheap-zero-zero", "false", expect);
	}

	/** Checks ordinary one-slot heap insertion. */
	private static void longHeapOnePositive(final String expect){
		final LongHeap heap=new LongHeap(1);
		checkEquals("longheap-one-positive-add", true, heap.add(1));
		checkEquals("longheap-one-positive-size", 1, heap.size());
		checkEquals("longheap-one-positive-capacity", 1, heap.capacity());
		checkEquals("longheap-one-positive-peek", 1, heap.peek());
		pass("longheap-one-positive", "true", expect);
	}

	/** Exposes zero-capacity set construction assertion. */
	private static void longHeapSetZeroPositive(final String expect){
		expectAssertion("longheapset-zero-positive", new Runnable(){
			@Override
			public void run(){new LongHeapSet(0);}
		});
		pass("longheapset-zero-positive", "ConstructorAssertionError", expect);
	}

	/** Checks ordinary one-slot heap-set insertion. */
	private static void longHeapSetOnePositive(final String expect){
		final LongHeapSet set=new LongHeapSet(1);
		checkEquals("longheapset-one-positive-add", true, set.add(1));
		checkEquals("longheapset-one-positive-contains", true, set.contains(1));
		checkEquals("longheapset-one-positive-size", 1, set.size());
		checkEquals("longheapset-one-positive-set-size", 1, set.set.size());
		checkEquals("longheapset-one-positive-heap-size", 1, set.heap.size());
		checkEquals("longheapset-one-positive-capacity", 1, set.capacity());
		checkEquals("longheapset-one-positive-peek", 1, set.peek());
		pass("longheapset-one-positive", "true", expect);
	}

	/** Exposes zero-capacity map construction assertion. */
	private static void longHeapMapZeroPositive(final String expect){
		expectAssertion("longheapmap-zero-positive", new Runnable(){
			@Override
			public void run(){new LongHeapMap(0);}
		});
		pass("longheapmap-zero-positive", "ConstructorAssertionError", expect);
	}

	/** Checks ordinary one-slot heap-map insertion. */
	private static void longHeapMapOnePositive(final String expect){
		final LongHeapMap map=new LongHeapMap(1);
		checkEquals("longheapmap-one-positive-add", true, map.add(1));
		checkEquals("longheapmap-one-positive-contains", true, map.contains(1));
		checkEquals("longheapmap-one-positive-count", 1, map.map.get(1));
		checkEquals("longheapmap-one-positive-size", 1, map.size());
		checkEquals("longheapmap-one-positive-map-size", 1, map.map.size());
		checkEquals("longheapmap-one-positive-heap-size", 1, map.heap.size());
		checkEquals("longheapmap-one-positive-capacity", 1, map.capacity());
		checkEquals("longheapmap-one-positive-peek", 1, map.peek());
		pass("longheapmap-one-positive", "true", expect);
	}

	/** Checks that an action throws AssertionError. */
	private static void expectAssertion(final String name, final Runnable action){
		try{
			action.run();
		}catch(final AssertionError expected){
			return;
		}
		throw new AssertionError(name+" expected AssertionError");
	}

	/** Checks an int oracle. */
	private static void checkEquals(final String name, final int expected, final int observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}

	/** Checks a long oracle. */
	private static void checkEquals(final String name, final long expected, final long observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}

	/** Checks a boolean oracle. */
	private static void checkEquals(final String name, final boolean expected,
			final boolean observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}

	/** Checks that JVM assertions are enabled. */
	private static void checkAssertionsEnabled(){
		check(LongHeapZeroCapacityCase.class.desiredAssertionStatus(),
				"Zero-capacity baseline checks require -ea");
	}

	/** Prints one successful case row. */
	private static void pass(final String name, final String observed,
			final String expect){
		System.out.println("PASS\t"+expect+'\t'+name+'\t'+observed);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

}
