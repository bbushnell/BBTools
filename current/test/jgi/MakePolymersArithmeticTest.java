package test.jgi;

import java.math.BigInteger;
import java.util.TreeSet;

import jgi.MakePolymers;
import structures.ByteBuilder;

/**
 * MakePolymers arithmetic oracle for k limits, repeat counts and record caps.
 * Directly calls only static production helpers; no generator is constructed.
 * @author Jean
 */
public final class MakePolymersArithmeticTest{

	/** Runs deterministic BigInteger-backed boundary checks. */
	public static void main(final String[] args){
		require(args.length==0, "MakePolymersArithmeticTest takes no arguments.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/makepolymersarithmetictest.sh with assertions enabled.");
		final MakePolymersArithmeticTest test=new MakePolymersArithmeticTest();
		test.checkMaxKeys();
		test.checkLengths();
		if(test.failures>0){
			throw new AssertionError("MakePolymers arithmetic failed "+test.failures+
				" of "+test.checks+" checks");
		}
		System.out.println("MAKE_POLYMERS_ARITHMETIC_PASS maxKeyChecks="+
			test.maxKeyChecks+" lengthCases="+test.lengthCases+" repeatChecks="+
			test.repeatChecks+" recordChecks="+test.recordChecks+" checks="+test.checks);
	}

	/** Checks every accepted two-bit key width against BigInteger shift math. */
	private void checkMaxKeys(){
		for(int k=1; k<=MakePolymers.MAX_K; k++){
			maxKeyChecks++;
			check(MakePolymers.maxKey(k)==oracleMaxKey(k), "maxKey mismatch for k="+k);
		}
		check(maxKeyChecks==31, "Expected k=1..31 max-key checks; saw "+maxKeyChecks);
	}

	/** Checks small, INT_MAX-adjacent and ByteBuilder-cap-adjacent lengths. */
	private void checkLengths(){
		for(int k=1; k<=MakePolymers.MAX_K; k++){
			final TreeSet<Integer> lengths=new TreeSet<Integer>();
			for(int len=1; len<=128; len++){lengths.add(len);}
			addIntMaxNeighbors(lengths, k);
			addRecordSizeBoundaries(lengths, k);
			for(final int minLen : lengths){checkLength(k, minLen);}
		}
		check(lengthCases>31*128, "Boundary lengths were not added; cases="+lengthCases);
	}

	/** Adds positive lengths near the signed int ceiling. */
	private static void addIntMaxNeighbors(final TreeSet<Integer> lengths, final int k){
		for(long delta=-2L*k-4; delta<=4; delta++){
			addIfInt(lengths, Integer.MAX_VALUE+delta);
		}
	}

	/** Adds positive lengths straddling the maximum ByteBuilder-backed record. */
	private static void addRecordSizeBoundaries(final TreeSet<Integer> lengths, final int k){
		final long overhead=Integer.toString(k).length()+
			BigInteger.ONE.shiftLeft(2*k).subtract(BigInteger.ONE).toString().length()+4;
		final long maxRepeat=(ByteBuilder.MAXLEN-overhead)/k;
		for(long repeat=maxRepeat-2; repeat<=maxRepeat+2; repeat++){
			final long firstForRepeat=(repeat-2)*k+2;
			final long lastForRepeat=(repeat-1)*k+1;
			for(long delta=-2; delta<=2; delta++){
				addIfInt(lengths, firstForRepeat+delta);
				addIfInt(lengths, lastForRepeat+delta);
			}
		}
		for(long delta=-2; delta<=2; delta++){
			addIfInt(lengths, ByteBuilder.MAXLEN-overhead+delta);
		}
	}

	/** Adds x if it can be supplied through MakePolymers' int minLen parser. */
	private static void addIfInt(final TreeSet<Integer> lengths, final long x){
		if(x>0 && x<=Integer.MAX_VALUE){lengths.add((int)x);}
	}

	/** Compares one production helper row to independent arithmetic oracles. */
	private void checkLength(final int k, final int minLen){
		lengthCases++;
		repeatChecks++;
		final int repeat=MakePolymers.repeatCount(k, minLen);
		final int expectedRepeat=oracleRepeatCount(k, minLen);
		check(repeat==expectedRepeat,
			"repeatCount mismatch for k="+k+", minLen="+minLen+
			": expected "+expectedRepeat+", saw "+repeat);
		final long target=(long)minLen+k-1;
		check((long)(repeat-1)*k<target && (long)repeat*k>=target,
			"repeatCount is not the ceiling repeat for k="+k+", minLen="+minLen+
			", repeat="+repeat);

		recordChecks++;
		final long recordBytes=MakePolymers.maxRecordBytes(k, minLen);
		final long expectedBytes=oracleRecordBytes(k, minLen);
		check(recordBytes==expectedBytes,
			"maxRecordBytes mismatch for k="+k+", minLen="+minLen+
			": expected "+expectedBytes+", saw "+recordBytes);
	}

	/** Returns ceil((minLen+k-1)/k) using BigInteger quotient/remainder. */
	private static int oracleRepeatCount(final int k, final int minLen){
		final BigInteger divisor=BigInteger.valueOf(k);
		final BigInteger target=BigInteger.valueOf(minLen).add(divisor).subtract(BigInteger.ONE);
		final BigInteger[] qr=target.divideAndRemainder(divisor);
		final BigInteger repeat=qr[1].signum()==0 ? qr[0] : qr[0].add(BigInteger.ONE);
		return repeat.intValueExact();
	}

	/** Returns 4^k-1 using BigInteger shift math. */
	private static long oracleMaxKey(final int k){
		return BigInteger.ONE.shiftLeft(2*k).subtract(BigInteger.ONE).longValueExact();
	}

	/** Returns the largest FASTA record size using decimal string lengths. */
	private static long oracleRecordBytes(final int k, final int minLen){
		final BigInteger repeat=BigInteger.valueOf(oracleRepeatCount(k, minLen));
		final BigInteger seqBases=repeat.multiply(BigInteger.valueOf(k));
		final int overhead=Integer.toString(k).length()+
			BigInteger.ONE.shiftLeft(2*k).subtract(BigInteger.ONE).toString().length()+4;
		return seqBases.add(BigInteger.valueOf(overhead)).longValueExact();
	}

	/** Accumulates errors so every boundary row is checked before failing. */
	private void check(final boolean correct, final String message){
		checks++;
		if(!correct){failures++; System.err.println("FAIL: "+message);}
	}

	/** Rejects invalid fixture or launcher state. */
	private static void require(final boolean correct, final String message){
		if(!correct){throw new IllegalArgumentException(message);}
	}

	private int maxKeyChecks, lengthCases, repeatChecks, recordChecks, checks, failures;
}
