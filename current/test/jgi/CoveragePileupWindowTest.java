package test.jgi;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;

import jgi.CoveragePileup;
import parse.Parse;

/**
 * Direct positive-window fixture and timing carrier for CoveragePileup.
 *
 * @author Jean
 */
public final class CoveragePileupWindowTest{

	/** oracle=old|fixed [warmup=1000] [iterations=10000] [trials=3] */
	public static void main(final String[] args){
		String oracle=null;
		long warmup=DEFAULT_WARMUP;
		long iterations=DEFAULT_ITERATIONS;
		int trials=DEFAULT_TRIALS;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected key=value argument: "+arg);
			final String key=pair[0].toLowerCase();
			final String value=pair[1];
			if("oracle".equals(key) || "expected".equals(key)){
				oracle=value;
			}else if("warmup".equals(key)){
				warmup=Long.parseLong(value);
			}else if("iterations".equals(key) || "iters".equals(key)){
				iterations=Long.parseLong(value);
			}else if("trials".equals(key)){
				trials=Integer.parseInt(value);
			}else{
				throw new IllegalArgumentException("Unknown CoveragePileupWindowTest option: "+arg);
			}
		}
		final boolean fixed;
		if("old".equals(oracle) || "baseline".equals(oracle)){
			fixed=false;
		}else if("fixed".equals(oracle) || "future".equals(oracle)){
			fixed=true;
		}else{
			throw new IllegalArgumentException("Set oracle=old or oracle=fixed.");
		}
		require(warmup>0, "warmup must be positive.");
		require(iterations>0, "iterations must be positive.");
		require(trials>0, "trials must be positive.");
		final ThreadMXBean clock=ManagementFactory.getThreadMXBean();
		require(clock.isCurrentThreadCpuTimeSupported() && clock.isThreadCpuTimeEnabled() &&
				clock.getCurrentThreadCpuTime()>=0,
			"Current-thread CPU time must be supported and enabled.");

		final CoveragePileup pileup=new CoveragePileup(new String[] {"in=stdin"});
		checkTinyCases(pileup, fixed);
		System.out.println("#stage\ttrial\toracle\tarray\tshape\tlength\twidth\tthreshold\titerations\texpected\tchecksum\tcpuNanos\twallNanos");
		for(final Shape shape : TIMING_SHAPES){
			measure(pileup, shape, false, fixed, "warmup", 0, warmup, clock);
			measure(pileup, shape, true, fixed, "warmup", 0, warmup, clock);
			for(int trial=1; trial<=trials; trial++){
				measure(pileup, shape, false, fixed, "trial", trial, iterations, clock);
				measure(pileup, shape, true, fixed, "trial", trial, iterations, clock);
			}
		}
		System.out.println("COVERAGE_PILEUP_WINDOW_TEST_PASS oracle="+(fixed ? "fixed" : "old")+
				" shapes="+TIMING_SHAPES.length+" trials="+trials+" warmup="+warmup+
				" iterations="+iterations);
	}

	/** Asserts the six public 1598 vectors outside the measured loops. */
	private static void checkTinyCases(final CoveragePileup pileup, final boolean fixed){
		for(final TinyCase tiny : TINY_CASES){
			final int charActual=pileup.basesUnderAverageCoverage(tiny.chars(), TINY_THRESHOLD, TINY_WIDTH);
			final int intActual=pileup.basesUnderAverageCoverage(tiny.ints(), TINY_THRESHOLD, TINY_WIDTH);
			final int charExpected=fixed ? tiny.fixed : tiny.oldChar;
			final int intExpected=fixed ? tiny.fixed : tiny.oldInt;
			require(charActual==charExpected,
				"char tiny case "+tiny.name+" expected "+charExpected+" but got "+charActual);
			require(intActual==intExpected,
				"int tiny case "+tiny.name+" expected "+intExpected+" but got "+intActual);
		}
	}

	/** Times one fixed array shape with the requested public overload. */
	private static long measure(final CoveragePileup pileup, final Shape shape,
			final boolean useInt, final boolean fixed, final String stage, final int trial,
			final long iterations, final ThreadMXBean clock){
		long checksum=0;
		shape.resetFirst(useInt);
		final long cpuStart=clock.getCurrentThreadCpuTime();
		final long wallStart=System.nanoTime();
		if(useInt){
			for(long i=0; i<iterations; i++){
				shape.alternateFirst(useInt, i);
				checksum+=pileup.basesUnderAverageCoverage(shape.ints, TIMING_THRESHOLD, TIMING_WIDTH);
			}
		}else{
			for(long i=0; i<iterations; i++){
				shape.alternateFirst(useInt, i);
				checksum+=pileup.basesUnderAverageCoverage(shape.chars, TIMING_THRESHOLD, TIMING_WIDTH);
			}
		}
		final long wallNanos=System.nanoTime()-wallStart;
		final long cpuNanos=clock.getCurrentThreadCpuTime()-cpuStart;
		final long expected=(fixed ? shape.fixed : (useInt ? shape.oldInt : shape.oldChar));
		require(checksum==expected*iterations,
				"checksum mismatch for "+shape.name+" "+(useInt ? "int" : "char")+
				": expected "+(expected*iterations)+" but got "+checksum);
		System.out.println(stage+'\t'+trial+'\t'+(fixed ? "fixed" : "old")+'\t'+
				(useInt ? "int" : "char")+'\t'+shape.name+'\t'+shape.length+'\t'+
				TIMING_WIDTH+'\t'+TIMING_THRESHOLD+'\t'+iterations+'\t'+expected+'\t'+
				checksum+'\t'+cpuNanos+'\t'+wallNanos);
		return checksum;
	}

	/** Returns an all-high coverage shape with no low window. */
	private static Shape allHigh(){
		final int[] array=new int[2000];
		for(int i=0; i<array.length; i++){array[i]=1;}
		return new Shape("all_high", array, 0, 0, 0, 1, 2);
	}

	/** Returns an all-low coverage shape whose low-window union spans the scaffold. */
	private static Shape allLow(){
		return new Shape("all_low", new int[2000], 2000, 2000, 2000, 0, 1);
	}

	/** Returns two separated single-zero low-window regions. */
	private static Shape separated(){
		final int[] array=new int[3000];
		for(int i=0; i<array.length; i++){array[i]=1;}
		array[600]=0;
		array[1800]=0;
		// Zeroes at 600 and 1800 independently make starts 101..600 and 1301..1800 low.
		// Fixed unions are [101,1100) and [1301,2300), for 999+999 bases.
		// Old char closes each run at 1000 bases; old int overwrites the first run.
		return new Shape("separated_low", array, 2000, 1000, 1998, 1, 2);
	}

	private static void require(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final class TinyCase{
		TinyCase(final String name_, final int[] coverage_, final int oldChar_,
				final int oldInt_, final int fixed_){
			name=name_;
			coverage=coverage_;
			oldChar=oldChar_;
			oldInt=oldInt_;
			fixed=fixed_;
		}

		char[] chars(){
			final char[] array=new char[coverage.length];
			for(int i=0; i<array.length; i++){array[i]=(char)coverage[i];}
			return array;
		}

		int[] ints(){return coverage.clone();}

		final String name;
		final int[] coverage;
		final int oldChar;
		final int oldInt;
		final int fixed;
	}

	private static final class Shape{
		Shape(final String name_, final int[] ints_, final int oldChar_,
				final int oldInt_, final int fixed_, final int firstA_, final int firstB_){
			name=name_;
			ints=ints_;
			chars=new char[ints.length];
			for(int i=0; i<ints.length; i++){chars[i]=(char)ints[i];}
			length=ints.length;
			oldChar=oldChar_;
			oldInt=oldInt_;
			fixed=fixed_;
			firstA=firstA_;
			firstB=firstB_;
			resetFirst(false);
			resetFirst(true);
		}

		/** Reset element zero before every timing segment. */
		void resetFirst(final boolean useInt){setFirst(useInt, firstA);}

		/**
		 * Alternates only element zero so repeated scans have an input dependency.
		 * The identical deterministic sequence is used for char/int and old/fixed
		 * timings, and values are chosen to preserve every expected low-window count.
		 */
		void alternateFirst(final boolean useInt, final long iteration){
			setFirst(useInt, ((iteration&1)==0) ? firstA : firstB);
		}

		private void setFirst(final boolean useInt, final int value){
			if(useInt){
				ints[0]=value;
			}else{
				chars[0]=(char)value;
			}
		}

		final String name;
		final int[] ints;
		final char[] chars;
		final int length;
		final int oldChar;
		final int oldInt;
		final int fixed;
		final int firstA;
		final int firstB;
	}

	private static final double TINY_THRESHOLD=1;
	private static final int TINY_WIDTH=3;
	private static final double TIMING_THRESHOLD=1;
	private static final int TIMING_WIDTH=500;
	private static final long DEFAULT_WARMUP=1000;
	private static final long DEFAULT_ITERATIONS=10000;
	private static final int DEFAULT_TRIALS=3;

	private static final TinyCase[] TINY_CASES={
		new TinyCase("exact", new int[] {2,0,0}, 0, 0, 3),
		new TinyCase("final", new int[] {3,0,0,0}, 0, 0, 3),
		new TinyCase("separated", new int[] {0,0,0,3,0,0,0,3,3,3}, 8, 4, 6),
		new TinyCase("overlap", new int[] {0,0,0,0,3,3}, 5, 5, 4),
		new TinyCase("nolow", new int[] {1,1,1,1,1}, 0, 0, 0),
		new TinyCase("opentail", new int[] {0,0,2,1,0,0,0}, 7, 4, 7),
	};
	private static final Shape[] TIMING_SHAPES={allHigh(), allLow(), separated()};
}
