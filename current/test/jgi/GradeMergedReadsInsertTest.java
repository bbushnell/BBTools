package test.jgi;

import jgi.GradeMergedReads;

/** Literal positive-integer range fixture for {@link GradeMergedReads#parseInsert(String)}. */
public final class GradeMergedReadsInsertTest{

	public static void main(String[] args){
		if(args.length!=0){
			throw new IllegalArgumentException("No arguments expected");
		}

		final Case[] cases={
			new Case("ordinary f_ insert1", randomReadsHeader("1"), 1L),
			new Case("f_ Integer.MAX_VALUE", randomReadsHeader("2147483647"), 2147483647L),
			new Case("tab-prefixed f_ insert300", "prefix\t"+randomReadsHeader("300"), 300L),
			new Case("leading-zero f_ insert0000000300", randomReadsHeader("0000000300"), 300L),
			new Case("insert=300_suffix control", "insert=300_suffix", 300L),
			new Case("f_2147483648", randomReadsHeader("2147483648"), 2147483648L),
			new Case("f_4294967296", randomReadsHeader("4294967296"), 4294967296L),
			new Case("f_4294967297", randomReadsHeader("4294967297"), 4294967297L),
			new Case("f_4294967298", randomReadsHeader("4294967298"), 4294967298L)
		};

		int failures=0;
		for(int i=0; i<cases.length; i++){
			final Case c=cases[i];
			final Outcome actual=Outcome.parse(c.header);
			final boolean expectException=c.expected>Integer.MAX_VALUE;
			final boolean pass=expectException ? actual.exception!=null : actual.value==c.expected;
			if(pass){
				System.out.println("PASS "+c.name+": expected "+expected(c)+", observed "+actual);
			}else{
				failures++;
				System.out.println("FAIL "+c.name+": expected "+expected(c)+", observed "+actual+
						", header="+c.header);
			}
		}

		System.out.println("GradeMergedReadsInsertTest cases="+cases.length+" failures="+failures);
		if(failures!=0){
			throw new AssertionError("GradeMergedReadsInsertTest failures="+failures);
		}
	}

	private static String expected(final Case c){
		return c.expected>Integer.MAX_VALUE ? "RuntimeException" : Long.toString(c.expected);
	}

	private static String randomReadsHeader(final String insert){
		return "f_0_c_0_s_0_p_0_i_"+insert+"_r_100";
	}

	private static final class Case{

		Case(final String name_, final String header_, final long expected_){
			// The exact long oracle treats values above Integer.MAX_VALUE as
			// positive overflow, so every row must stay on that side of zero.
			assert(expected_>0) : "Positive-only overflow oracle requires a positive expected value: "+expected_;
			name=name_;
			header=header_;
			expected=expected_;
		}

		final String name;
		final String header;
		final long expected;
	}

	private static final class Outcome{

		private static Outcome parse(final String header){
			try{
				return new Outcome(GradeMergedReads.parseInsert(header), null);
			}catch(RuntimeException e){
				return new Outcome(Long.MIN_VALUE, e);
			}
		}

		private Outcome(final long value_, final RuntimeException exception_){
			value=value_;
			exception=exception_;
		}

		@Override
		public String toString(){
			return exception==null ? Long.toString(value) : exception.getClass().getSimpleName()+": "+exception.getMessage();
		}

		final long value;
		final RuntimeException exception;
	}
}
