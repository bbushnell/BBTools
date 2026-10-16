package test.structures;

import structures.IntRingBuffer;
import structures.RingBuffer;
import structures.RingBufferMask;

/** Baseline and fixed oracles for ring constructor bounds.
 * @author Jean
 */
public final class RingBufferBoundsCase{

	/** Runs tiny valid wrap controls and bounded unsupported-constructor checks. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		final boolean fixed="fixed".equals(expect);

		if(fixed){checkProductionAssertions();}

		intLogical1(expect);
		intLogical3(expect);
		longLogical1(expect);
		longLogical3(expect);
		maskBits0(expect);
		maskBits2(expect);

		intSize0(fixed, expect);
		intSizeNegative1(fixed, expect);
		longSize0(fixed, expect);
		longSizeNegative1(fixed, expect);
		maskBits32(fixed, expect);
		maskBitsNegative32(fixed, expect);

		if(fixed){
			intSizeTooLarge(expect);
			longSizeTooLarge(expect);
			maskBits31(expect);
		}
	}

	/** Checks that logical capacity one keeps only the second int value. */
	private static void intLogical1(final String expect){
		final IntRingBuffer ring=new IntRingBuffer(1);
		ring.add(11);
		ring.add(22);
		checkObservation(expect, "int_logical1_wrap", "1,22,22",
				ring.size()+","+ring.getOldest()+","+ring.get(0));
	}

	/** Checks that logical capacity three drops 11 after four int adds. */
	private static void intLogical3(final String expect){
		final IntRingBuffer ring=new IntRingBuffer(3);
		ring.add(11);
		ring.add(12);
		ring.add(13);
		ring.add(14);
		checkObservation(expect, "int_logical3_wrap", "3,12,14,13,12",
				ring.size()+","+ring.getOldest()+","+ring.get(0)+","+
				ring.get(1)+","+ring.get(2));
	}

	/** Checks that logical capacity one keeps only the second long value. */
	private static void longLogical1(final String expect){
		final RingBuffer ring=new RingBuffer(1);
		ring.add(11);
		ring.add(22);
		checkObservation(expect, "long_logical1_wrap", "1,22,22",
				ring.size()+","+ring.getOldest()+","+ring.get(0));
	}

	/** Checks that logical capacity three drops 11 after four long adds. */
	private static void longLogical3(final String expect){
		final RingBuffer ring=new RingBuffer(3);
		ring.add(11);
		ring.add(12);
		ring.add(13);
		ring.add(14);
		checkObservation(expect, "long_logical3_wrap", "3,12,14,13,12",
				ring.size()+","+ring.getOldest()+","+ring.get(0)+","+
				ring.get(1)+","+ring.get(2));
	}

	/** Checks that mask bit count zero creates one slot and keeps only 22. */
	private static void maskBits0(final String expect){
		final RingBufferMask ring=new RingBufferMask(0);
		ring.add(11);
		ring.add(22);
		checkObservation(expect, "mask_bits0_wrap", "1,22,22",
				ring.size()+","+ring.getOldest()+","+ring.get(0));
	}

	/** Checks that mask bit count two creates four slots and drops 11 after five adds. */
	private static void maskBits2(final String expect){
		final RingBufferMask ring=new RingBufferMask(2);
		ring.add(11);
		ring.add(12);
		ring.add(13);
		ring.add(14);
		ring.add(15);
		checkObservation(expect, "mask_bits2_wrap", "4,12,15,14,13,12",
				ring.size()+","+ring.getOldest()+","+ring.get(0)+","+
				ring.get(1)+","+ring.get(2)+","+ring.get(3));
	}

	/** Records that unsupported int size zero currently uses a two-slot backing. */
	private static void intSize0(final boolean fixed, final String expect){
		final String name="int_size0_backing2";
		if(fixed){expectConstructorAssertion(expect, name, new Runnable(){
			@Override
			public void run(){new IntRingBuffer(0);}
		}); return;}
		final IntRingBuffer ring=new IntRingBuffer(0);
		ring.add(11);
		ring.add(22);
		ring.add(33);
		checkObservation(expect, name, "0,33,22,22",
				ring.size()+","+ring.get(0)+","+ring.get(1)+","+ring.getCurrent());
	}

	/** Records that unsupported int size -1 currently uses a two-slot backing. */
	private static void intSizeNegative1(final boolean fixed, final String expect){
		final String name="int_size_negative1_backing2";
		if(fixed){expectConstructorAssertion(expect, name, new Runnable(){
			@Override
			public void run(){new IntRingBuffer(-1);}
		}); return;}
		final IntRingBuffer ring=new IntRingBuffer(-1);
		ring.add(11);
		ring.add(22);
		ring.add(33);
		checkObservation(expect, name, "-1,33,22,22",
				ring.size()+","+ring.get(0)+","+ring.get(1)+","+ring.getCurrent());
	}

	/** Records that unsupported long size zero currently uses a two-slot backing. */
	private static void longSize0(final boolean fixed, final String expect){
		final String name="long_size0_backing2";
		if(fixed){expectConstructorAssertion(expect, name, new Runnable(){
			@Override
			public void run(){new RingBuffer(0);}
		}); return;}
		final RingBuffer ring=new RingBuffer(0);
		ring.add(11);
		ring.add(22);
		ring.add(33);
		checkObservation(expect, name, "0,33,22,22",
				ring.size()+","+ring.get(0)+","+ring.get(1)+","+ring.getCurrent());
	}

	/** Records that unsupported long size -1 currently uses a two-slot backing. */
	private static void longSizeNegative1(final boolean fixed, final String expect){
		final String name="long_size_negative1_backing2";
		if(fixed){expectConstructorAssertion(expect, name, new Runnable(){
			@Override
			public void run(){new RingBuffer(-1);}
		}); return;}
		final RingBuffer ring=new RingBuffer(-1);
		ring.add(11);
		ring.add(22);
		ring.add(33);
		checkObservation(expect, name, "-1,33,22,22",
				ring.size()+","+ring.get(0)+","+ring.get(1)+","+ring.getCurrent());
	}

	/** Records that unsupported bit count 32 currently wraps to one slot. */
	private static void maskBits32(final boolean fixed, final String expect){
		final String name="mask_bits32_backing1";
		if(fixed){expectConstructorAssertion(expect, name, new Runnable(){
			@Override
			public void run(){new RingBufferMask(32);}
		}); return;}
		final RingBufferMask ring=new RingBufferMask(32);
		ring.add(11);
		ring.add(22);
		ring.add(33);
		checkObservation(expect, name, "1,33,33,33",
				ring.size()+","+ring.getOldest()+","+ring.get(0)+","+ring.getCurrent());
	}

	/** Records that unsupported bit count -32 currently wraps to one slot. */
	private static void maskBitsNegative32(final boolean fixed, final String expect){
		final String name="mask_bits_negative32_backing1";
		if(fixed){expectConstructorAssertion(expect, name, new Runnable(){
			@Override
			public void run(){new RingBufferMask(-32);}
		}); return;}
		final RingBufferMask ring=new RingBufferMask(-32);
		ring.add(11);
		ring.add(22);
		ring.add(33);
		checkObservation(expect, name, "1,33,33,33",
				ring.size()+","+ring.getOldest()+","+ring.get(0)+","+ring.getCurrent());
	}

	/** Checks fixed-only rejection above the maximum int-ring logical capacity. */
	private static void intSizeTooLarge(final String expect){
		expectConstructorAssertion(expect, "int_size_too_large", new Runnable(){
			@Override
			public void run(){new IntRingBuffer((1<<30)+1);}
		});
	}

	/** Checks fixed-only rejection above the maximum long-ring logical capacity. */
	private static void longSizeTooLarge(final String expect){
		expectConstructorAssertion(expect, "long_size_too_large", new Runnable(){
			@Override
			public void run(){new RingBuffer((1<<30)+1);}
		});
	}

	/** Checks fixed-only rejection of the first mask bit count above the range. */
	private static void maskBits31(final String expect){
		expectConstructorAssertion(expect, "mask_bits31", new Runnable(){
			@Override
			public void run(){new RingBufferMask(31);}
		});
	}

	/** Checks that an unsupported constructor throws AssertionError. */
	private static void expectConstructorAssertion(final String expect,
			final String name, final Runnable action){
		try{
			action.run();
		}catch(final AssertionError expected){
			pass(expect, name, "ConstructorAssertionError");
			return;
		}
		throw new AssertionError(name+" expected AssertionError");
	}

	/** Verifies production assertions before fixed invalid-constructor checks. */
	private static void checkProductionAssertions(){
		check(IntRingBuffer.class.desiredAssertionStatus(),
				"IntRingBuffer fixed bounds checks require -ea");
		check(RingBuffer.class.desiredAssertionStatus(),
				"RingBuffer fixed bounds checks require -ea");
		check(RingBufferMask.class.desiredAssertionStatus(),
				"RingBufferMask fixed bounds checks require -ea");
	}

	/** Checks and prints one observed value. */
	private static void checkObservation(final String expect, final String name,
			final String expected, final String observed){
		check(expected.equals(observed), name+" expected "+expected+", got "+observed);
		pass(expect, name, observed);
	}

	/** Prints one successful case row. */
	private static void pass(final String expect, final String name,
			final String observed){
		System.out.println("PASS\t"+expect+'\t'+name+'\t'+observed);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

}
