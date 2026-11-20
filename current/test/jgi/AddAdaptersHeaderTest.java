package test.jgi;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;

import jgi.AddAdapters;
import stream.FASTQ;

/**
 * Fixed-only scalar header-decoder cases with independent literal length components.
 * Reflects the real private static decodeExpectedLengths(String,int); no Read,
 * AddAdapters instance, stream or generic header parser is constructed. Each case
 * runs with custom parsing false and true and must preserve that setting.
 * Reflection boxes the returned long; this is functional evidence, not an allocation benchmark.
 * @author Sucrose
 */
public final class AddAdaptersHeaderTest{

	/** Accepts no arguments or {@code phase=fixed}; the decoder must already be implemented. */
	public static void main(final String[] args) throws Exception{
		if(args.length!=0 && !(args.length==1 && args[0].equals("phase=fixed"))){
			throw new IllegalArgumentException("This scalar decoder fixture accepts only phase=fixed.");
		}
		final Method decoder=AddAdapters.class.getDeclaredMethod("decodeExpectedLengths", String.class, int.class);
		if(decoder.getReturnType()!=long.class || !Modifier.isStatic(decoder.getModifiers()) || !Modifier.isPrivate(decoder.getModifiers())){
			throw new AssertionError("Expected the settled private static long String/int decoder interface.");
		}
		decoder.setAccessible(true);
		final boolean original=FASTQ.PARSE_CUSTOM;
		try{
			for(final boolean flag : new boolean[]{false, true}){
				for(final HeaderCase row : VALID){run(decoder, row, flag, true);}
				for(final HeaderCase row : INVALID){run(decoder, row, flag, false);}
			}
		}finally{
			FASTQ.PARSE_CUSTOM=original;
		}
		check(FASTQ.PARSE_CUSTOM==original, "Fixture must restore its inherited custom-parser flag.");
		check(VALID.length==15 && INVALID.length==27 && calls==84,
			"Frozen literal matrix requires15 valid/27 invalid cases in both parser states: "+VALID.length+"/"+INVALID.length+"/"+calls);
		if(failures>0){throw new AssertionError("Header decoder failed "+failures+" of "+checks+" checks.");}
		System.out.println("ADD_ADAPTERS_HEADER_PASS phase=fixed validCases=15 invalidCases=27 calls="+calls+" parserStates=2 checks="+checks);
	}

	/** Unexpected errors cannot masquerade as malformed-header rejection. */
	private static void run(final Method decoder, final HeaderCase row, final boolean flag, final boolean valid) throws Exception{
		FASTQ.PARSE_CUSTOM=flag;
		long original=-1, remaining=-1;
		boolean rejected=false;
		String error="<none>";
		try{
			final long packed=((Long)decoder.invoke(null, row.id, row.mate)).longValue();
			original=packed>>>32;
			remaining=packed&0xffffffffL;
		}catch(InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(valid || !(cause instanceof IllegalArgumentException)){
				throw new AssertionError("Unexpected decoder failure for "+row.name+", header="+display(row.id), cause);
			}
			rejected=true;
			error=cause.getClass().getSimpleName()+":"+display(cause.getMessage());
		}
		calls++;
		System.out.println("header_case\t"+row.name+"\tparserFlag="+flag+"\tparserAfter="+FASTQ.PARSE_CUSTOM+
			"\tmate="+row.mate+"\tid="+display(row.id)+"\toutcome="+(rejected ? "rejected" : "value")+
			"\toriginal="+original+"\tremaining="+remaining+"\terror="+error);
		check(FASTQ.PARSE_CUSTOM==flag, row.name+": scalar decoding must not change FASTQ.PARSE_CUSTOM="+flag);
		if(valid){
			check(original==row.original, row.name+": upper32 must equal literal original length "+row.original+", observed="+original);
			check(remaining==row.remaining, row.name+": lower32 must equal literal remaining length "+row.remaining+", observed="+remaining);
		}else{
			check(rejected, row.name+": malformed or ambiguous input must throw IllegalArgumentException, not return "+original+"/"+remaining);
		}
	}

	/** Keeps whitespace and non-ASCII invalid inputs visible in one ASCII observation row. */
	private static String display(final String text){
		if(text==null){return "<null>";}
		final StringBuilder out=new StringBuilder("[");
		for(int i=0; i<text.length(); i++){
			final char c=text.charAt(i);
			if(c>=32 && c<127 && c!='\\'){out.append(c);}
			else{
				out.append("\\u");
				for(int shift=12; shift>=0; shift-=4){out.append(HEX[(c>>>shift)&15]);}
			}
		}
		return out.append(']').toString();
	}
	private static void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println("CHECK_FAILED "+message);}
	}

	/** No oracle uses the decoder or reconstructs a packed expected value. */
	private static final class HeaderCase{
		HeaderCase(final String name_, final String id_, final int mate_, final int original_, final int remaining_){
			name=name_; id=id_; mate=mate_; original=original_; remaining=remaining_;
		}
		HeaderCase(final String name_, final String id_, final int mate_){this(name_, id_, mate_, -1, -1);}
		final String name, id;
		final int mate, original, remaining;
	}
	private static final String FIRST="7_0_5_3_+_100_1_._ref", SECOND="7_0_7_3_-_200_1_._ref";
	private static final String SINGLE="SYN_"+FIRST, PAIRED=SINGLE+"&"+SECOND;
	private static final HeaderCase[] VALID={
		new HeaderCase("legacy-zero", "7_6_0", -1, 6, 0),
		new HeaderCase("legacy-empty", "7_0_0", -1, 0, 0),
		new HeaderCase("legacy-untrimmed", "7_6_6", -1, 6, 6),
		new HeaderCase("rename-first", "7_6_2 1:", -1, 6, 2),
		new HeaderCase("rename-second", "7_8_5 2:", -1, 8, 5),
		new HeaderCase("shared-suffix-first", "7_6_2_8_5 /1", 1, 6, 2),
		new HeaderCase("shared-suffix-second", "7_6_2_8_5 /2", 0, 8, 5),
		new HeaderCase("shared-context-first", "7_6_2_8_5", 0, 6, 2),
		new HeaderCase("shared-context-second", "7_6_2_8_5", 1, 8, 5),
		new HeaderCase("legacy-maximum", "7_2147483647_2147483647", -1, Integer.MAX_VALUE, Integer.MAX_VALUE),
		new HeaderCase("synthetic-single", SINGLE, -1, 6, 3),
		new HeaderCase("synthetic-clamped", "SYN_7_0_5_9_+_100_1_._ref", -1, 6, 6),
		new HeaderCase("synthetic-suffix-second", PAIRED+" 2:", 0, 8, 3),
		new HeaderCase("synthetic-context-first", PAIRED, 0, 6, 3),
		new HeaderCase("synthetic-context-second", PAIRED, 1, 8, 3)
	};
	private static final HeaderCase[] INVALID={
		new HeaderCase("null", null, -1),
		new HeaderCase("empty", "", -1),
		new HeaderCase("two-fields", "7_6", -1),
		new HeaderCase("four-fields", "7_6_2_8", 0),
		new HeaderCase("six-fields", "7_6_2_8_5_9", 0),
		new HeaderCase("empty-original", "7__2", -1),
		new HeaderCase("empty-remaining", "7_6_", -1),
		new HeaderCase("nondigit", "7_6_x", -1),
		new HeaderCase("nonascii-digit", "7_6_\u0662", -1),
		new HeaderCase("negative-original", "7_-1_0", -1),
		new HeaderCase("negative-remaining", "7_6_-1", -1),
		new HeaderCase("remaining-over-original", "7_6_7", -1),
		new HeaderCase("original-overflow", "7_2147483648_0", -1),
		new HeaderCase("remaining-overflow", "7_2147483647_2147483648", -1),
		new HeaderCase("long-overflow", "7_999999999999999999999999999999_0", -1),
		new HeaderCase("shared-ambiguous", "7_6_2_8_5", -1),
		new HeaderCase("synthetic-ambiguous", PAIRED, -1),
		new HeaderCase("synthetic-empty-insert", "SYN_7_0_5__+_100_1_._ref", -1),
		new HeaderCase("synthetic-zero-insert", "SYN_7_0_5_0_+_100_1_._ref", -1),
		new HeaderCase("synthetic-negative-insert", "SYN_7_0_5_-1_+_100_1_._ref", -1),
		new HeaderCase("synthetic-reversed-span", "SYN_7_5_0_3_+_100_1_._ref", -1),
		new HeaderCase("synthetic-span-overflow", "SYN_7_-2147483648_2147483647_3_+_0_1_._ref", -1),
		new HeaderCase("synthetic-eight-fields", "SYN_7_0_5_3_+_100_1_.", -1),
		new HeaderCase("synthetic-number-overflow", "SYN_7_0_5_999999999999999999999999_+_100_1_._ref", -1),
		new HeaderCase("mate-below-range", "7_6_2", -2),
		new HeaderCase("mate-above-range", "7_6_2", 2),
		new HeaderCase("synthetic-three-bodies", PAIRED+"&"+FIRST, 0)
	};
	private static final char[] HEX="0123456789abcdef".toCharArray();
	private static int checks, failures, calls;
}
