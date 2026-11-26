package test.jgi;

import java.util.Arrays;

import jgi.Info;

/**
 * Literal count-width checks for the Info multinomial helpers.
 */
public final class InfoCountTest{

	public static void main(String[] args){
		check("all-A length=0", repeat('A', 0), 1, 2, -1);
		check("all-A length=1", repeat('A', 1), 1, 2, -1);
		check("all-A length=32767", repeat('A', 32767), 1, 2, -1);
		check("all-A length=32768", repeat('A', 32768), 1, 2, -1);
		check("all-A length=32769", repeat('A', 32769), 1, 2, -1);
		check("all-A length=65536", repeat('A', 65536), 1, 2, -1);
		check("AACC", "AACC".getBytes(), 6, 6, 4);
		check("AaNCC", "AaNCC".getBytes(), 6, 6, 5);

		System.out.println("InfoCountTest cases="+cases+" failures="+failures);
		if(failures>0){throw new AssertionError(failures+" failed cases");}
	}

	private static void check(String name, byte[] bases, long expectedInfo,
			long threshold, int expectedPrefix){
		assert(name!=null) : "test name is required for stable fixture output";
		assert(bases!=null) : "test bases are required for "+name;
		cases++;
		final LongResult actualInfo=info(bases);
		final IntResult actualPrefix=prefix(bases, threshold);
		final boolean pass=actualInfo.matches(expectedInfo) &&
				actualPrefix.matches(expectedPrefix);
		if(!pass){failures++;}

		System.out.println((pass ? "PASS " : "FAIL ")+name+
				" info="+actualInfo+" expected="+expectedInfo+
				" prefix="+actualPrefix+" expected="+expectedPrefix);
	}

	private static LongResult info(byte[] bases){
		try{
			return new LongResult(Info.info(bases));
		}catch(RuntimeException e){
			return new LongResult(e);
		}
	}

	private static IntResult prefix(byte[] bases, long threshold){
		try{
			return new IntResult(Info.prefixForInfo(bases, threshold));
		}catch(RuntimeException e){
			return new IntResult(e);
		}
	}

	private static byte[] repeat(char base, int length){
		final byte[] bases=new byte[length];
		Arrays.fill(bases, (byte)base);
		return bases;
	}

	private static class LongResult{

		LongResult(long value_){
			value=value_;
			exception=null;
		}

		LongResult(RuntimeException exception_){
			value=0;
			exception=exception_;
		}

		boolean matches(long expected){
			return exception==null && value==expected;
		}

		@Override
		public String toString(){
			return exception==null ? Long.toString(value) :
				exception.getClass().getSimpleName();
		}

		final long value;
		final RuntimeException exception;
	}

	private static class IntResult{

		IntResult(int value_){
			value=value_;
			exception=null;
		}

		IntResult(RuntimeException exception_){
			value=0;
			exception=exception_;
		}

		boolean matches(int expected){
			return exception==null && value==expected;
		}

		@Override
		public String toString(){
			return exception==null ? Integer.toString(value) :
				exception.getClass().getSimpleName();
		}

		final int value;
		final RuntimeException exception;
	}

	private static int cases=0;
	private static int failures=0;

}
