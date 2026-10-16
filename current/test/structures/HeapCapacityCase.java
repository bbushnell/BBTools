package test.structures;

import structures.Heap;
import structures.HeapLoc;
import structures.SetLoc;

/** Finite capacity-boundary cases for Heap and HeapLoc.
 * @author Jean
 */
public final class HeapCapacityCase{

	/** Runs finite constructor and HeapLoc resize-capacity checks. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		checkAssertionsEnabled();

		heapZeroRollover(expect);
		heapNegativeRollover(expect);
		heapLocZeroRollover(expect);
		heapLocNegativeRollover(expect);
		heapOneNonRollover(expect);
		heapOneRollover(expect);
		heapLocOneNonRollover(expect);
		heapLocResizeExact(expect);
		heapLocResizeShrinkSpare(expect);
		heapLocResizeUndersized(false, expect);
		heapLocResizeUndersized(true, expect);
	}

	/** Exposes zero-capacity Heap rollover construction or insertion failure. */
	private static void heapZeroRollover(final String expect){
		if("fixed".equals(expect)){
			expectAssertion("heap-zero-rollover", new Runnable(){
				@Override
				public void run(){new Heap<Integer>(0, true);}
			});
			pass("heap-zero-rollover", "ConstructorAssertionError", expect);
			return;
		}
		final Heap<Integer> heap=new Heap<Integer>(0, true);
		checkEquals("heap-zero-rollover-add", "NullPointerException", observe(new Runnable(){
			@Override
			public void run(){heap.add(Integer.valueOf(1));}
		}));
		checkEquals("heap-zero-rollover-size", 0, heap.size());
		pass("heap-zero-rollover", "NullPointerException", expect);
	}

	/** Exposes negative-capacity Heap rollover construction or insertion failure. */
	private static void heapNegativeRollover(final String expect){
		if("fixed".equals(expect)){
			expectAssertion("heap-negative-rollover", new Runnable(){
				@Override
				public void run(){new Heap<Integer>(-1, true);}
			});
			pass("heap-negative-rollover", "ConstructorAssertionError", expect);
			return;
		}
		final Heap<Integer> heap=new Heap<Integer>(-1, true);
		checkEquals("heap-negative-rollover-add", "ArrayIndexOutOfBoundsException",
				observe(new Runnable(){
					@Override
					public void run(){heap.add(Integer.valueOf(1));}
				}));
		checkEquals("heap-negative-rollover-size", 0, heap.size());
		pass("heap-negative-rollover", "ArrayIndexOutOfBoundsException", expect);
	}

	/** Exposes zero-capacity HeapLoc rollover construction or insertion failure. */
	private static void heapLocZeroRollover(final String expect){
		if("fixed".equals(expect)){
			expectAssertion("heaploc-zero-rollover", new Runnable(){
				@Override
				public void run(){new HeapLoc<Item>(0, true);}
			});
			pass("heaploc-zero-rollover", "ConstructorAssertionError", expect);
			return;
		}
		final HeapLoc<Item> heap=new HeapLoc<Item>(0, true);
		final Item item=new Item(1);
		checkEquals("heaploc-zero-rollover-add", "NullPointerException", observe(new Runnable(){
			@Override
			public void run(){heap.add(item);}
		}));
		checkEquals("heaploc-zero-rollover-size", 0, heap.size());
		checkEquals("heaploc-zero-rollover-loc", -1, item.loc());
		pass("heaploc-zero-rollover", "NullPointerException", expect);
	}

	/** Exposes negative-capacity HeapLoc rollover construction or insertion failure. */
	private static void heapLocNegativeRollover(final String expect){
		if("fixed".equals(expect)){
			expectAssertion("heaploc-negative-rollover", new Runnable(){
				@Override
				public void run(){new HeapLoc<Item>(-1, true);}
			});
			pass("heaploc-negative-rollover", "ConstructorAssertionError", expect);
			return;
		}
		final HeapLoc<Item> heap=new HeapLoc<Item>(-1, true);
		final Item item=new Item(1);
		checkEquals("heaploc-negative-rollover-add", "ArrayIndexOutOfBoundsException",
				observe(new Runnable(){
					@Override
					public void run(){heap.add(item);}
				}));
		checkEquals("heaploc-negative-rollover-size", 0, heap.size());
		checkEquals("heaploc-negative-rollover-loc", -1, item.loc());
		pass("heaploc-negative-rollover", "ArrayIndexOutOfBoundsException", expect);
	}

	/** Checks ordinary one-slot Heap insertion and removal without rollover. */
	private static void heapOneNonRollover(final String expect){
		final Heap<Integer> heap=new Heap<Integer>(1, false);
		checkEquals("heap-one-nonrollover-add", true, heap.add(Integer.valueOf(1)));
		checkEquals("heap-one-nonrollover-size", 1, heap.size());
		checkEquals("heap-one-nonrollover-peek", Integer.valueOf(1), heap.peek());
		checkEquals("heap-one-nonrollover-poll", Integer.valueOf(1), heap.poll());
		checkEquals("heap-one-nonrollover-empty", true, heap.isEmpty());
		pass("heap-one-nonrollover", "true", expect);
	}

	/** Checks one-slot Heap rollover rejection and eviction. */
	private static void heapOneRollover(final String expect){
		final Heap<Integer> heap=new Heap<Integer>(1, true);
		checkEquals("heap-one-rollover-add-one", true, heap.add(Integer.valueOf(1)));
		checkEquals("heap-one-rollover-add-zero", false, heap.add(Integer.valueOf(0)));
		checkEquals("heap-one-rollover-peek-one", Integer.valueOf(1), heap.peek());
		checkEquals("heap-one-rollover-add-two", true, heap.add(Integer.valueOf(2)));
		checkEquals("heap-one-rollover-size", 1, heap.size());
		checkEquals("heap-one-rollover-peek-two", Integer.valueOf(2), heap.peek());
		pass("heap-one-rollover", "true", expect);
	}

	/** Checks ordinary one-slot HeapLoc insertion and removal without rollover. */
	private static void heapLocOneNonRollover(final String expect){
		final HeapLoc<Item> heap=new HeapLoc<Item>(1, false);
		final Item item=new Item(7);
		checkEquals("heaploc-one-nonrollover-add", true, heap.add(item));
		checkEquals("heaploc-one-nonrollover-loc", 1, item.loc());
		checkEquals("heaploc-one-nonrollover-poll", item, heap.poll());
		checkEquals("heaploc-one-nonrollover-polled-loc", -1, item.loc());
		checkEquals("heaploc-one-nonrollover-empty", true, heap.isEmpty());
		pass("heaploc-one-nonrollover", "true", expect);
	}

	/** Checks exact-capacity HeapLoc resize transfer and source clearing. */
	private static void heapLocResizeExact(final String expect){
		final HeapLoc<Item> heap=new HeapLoc<Item>(2, true);
		final Item one=new Item(1);
		final Item two=new Item(2);
		checkEquals("heaploc-resize-exact-add-one", true, heap.add(one));
		checkEquals("heaploc-resize-exact-add-two", true, heap.add(two));
		final HeapLoc<Item> resized=heap.resizeNew(2);
		checkEquals("heaploc-resize-exact-old-empty", true, heap.isEmpty());
		checkEquals("heaploc-resize-exact-size", 2, resized.size());
		checkEquals("heaploc-resize-exact-one-loc", 1, one.loc());
		checkEquals("heaploc-resize-exact-two-loc", 2, two.loc());
		checkEquals("heaploc-resize-exact-poll-one", one, resized.poll());
		checkEquals("heaploc-resize-exact-poll-two", two, resized.poll());
		checkEquals("heaploc-resize-exact-one-polled-loc", -1, one.loc());
		checkEquals("heaploc-resize-exact-two-polled-loc", -1, two.loc());
		pass("heaploc-resize-exact", "true", expect);
	}

	/** Checks shrink-to-fit HeapLoc resize transfer and source clearing. */
	private static void heapLocResizeShrinkSpare(final String expect){
		final HeapLoc<Item> heap=new HeapLoc<Item>(3, true);
		final Item one=new Item(1);
		final Item two=new Item(2);
		checkEquals("heaploc-resize-shrink-add-one", true, heap.add(one));
		checkEquals("heaploc-resize-shrink-add-two", true, heap.add(two));
		final HeapLoc<Item> resized=heap.resizeNew(2);
		checkEquals("heaploc-resize-shrink-old-empty", true, heap.isEmpty());
		checkEquals("heaploc-resize-shrink-capacity", 2, resized.CAPACITY);
		checkEquals("heaploc-resize-shrink-size", 2, resized.size());
		checkEquals("heaploc-resize-shrink-one-loc", 1, one.loc());
		checkEquals("heaploc-resize-shrink-two-loc", 2, two.loc());
		checkEquals("heaploc-resize-shrink-poll-one", one, resized.poll());
		checkEquals("heaploc-resize-shrink-poll-two", two, resized.poll());
		checkEquals("heaploc-resize-shrink-one-polled-loc", -1, one.loc());
		checkEquals("heaploc-resize-shrink-two-polled-loc", -1, two.loc());
		pass("heaploc-resize-shrink", "true", expect);
	}

	/** Exposes resizeNew's undersized replacement mutation path. */
	private static void heapLocResizeUndersized(final boolean rollover, final String expect){
		final HeapLoc<Item> heap=new HeapLoc<Item>(2, rollover);
		final Item one=new Item(1);
		final Item two=new Item(2);
		checkEquals("heaploc-resize-undersized-add-one", true, heap.add(one));
		checkEquals("heaploc-resize-undersized-add-two", true, heap.add(two));
		final String observed=observe(new Runnable(){
			@Override
			public void run(){heap.resizeNew(1);}
		});
		checkEquals("heaploc-resize-undersized-throw", "AssertionError", observed);
		checkEquals("heaploc-resize-undersized-size", 2, heap.size());
		final String name=rollover ? "heaploc-resize-undersized-rollover" :
				"heaploc-resize-undersized-nonrollover";
		if("fixed".equals(expect)){
			checkEquals(name+"-one-loc", 1, one.loc());
			checkEquals(name+"-two-loc", 2, two.loc());
		}else if(rollover){
			checkEquals(name+"-one-loc", -1, one.loc());
			checkEquals(name+"-two-loc", 1, two.loc());
		}else{
			checkEquals(name+"-one-loc", 1, one.loc());
			checkEquals(name+"-two-loc", -1, two.loc());
		}
		pass(name, observed+",size="+heap.size()+",oneLoc="+one.loc()+",twoLoc="+two.loc(), expect);
	}

	/** Checks that an action throws AssertionError. */
	private static void expectAssertion(final String name, final Runnable action){
		checkEquals(name, "AssertionError", observe(action));
	}

	/** Returns the thrown class simple name, or no-throw when an action completes. */
	private static String observe(final Runnable action){
		try{
			action.run();
		}catch(final Throwable observed){
			return observed.getClass().getSimpleName();
		}
		return "no-throw";
	}

	/** Checks an int oracle. */
	private static void checkEquals(final String name, final int expected, final int observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}

	/** Checks a boolean oracle. */
	private static void checkEquals(final String name, final boolean expected,
			final boolean observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}

	/** Checks an object oracle. */
	private static void checkEquals(final String name, final Object expected,
			final Object observed){
		check(expected==observed || (expected!=null && expected.equals(observed)),
				name+" expected "+expected+", got "+observed);
	}

	/** Checks that JVM assertions are enabled for every relevant class. */
	private static void checkAssertionsEnabled(){
		check(HeapCapacityCase.class.desiredAssertionStatus(),
				"HeapCapacityCase assertions must be enabled");
		check(Heap.class.desiredAssertionStatus(), "Heap assertions must be enabled");
		check(HeapLoc.class.desiredAssertionStatus(), "HeapLoc assertions must be enabled");
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

	/** Small HeapLoc element with a mutable one-based backing location. */
	private static final class Item implements SetLoc<Item>{

		/** Creates an absent item with the given ordering value. */
		Item(final int value_){
			value=value_;
		}

		@Override
		public int compareTo(final Item other){
			return value-other.value;
		}

		@Override
		public void setLoc(final int newLoc){
			loc=newLoc;
		}

		@Override
		public int loc(){
			return loc;
		}

		private final int value;
		private int loc=-1;
	}

}
