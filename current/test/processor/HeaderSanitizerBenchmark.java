package test.processor;

import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;

import jgi.ReformatReads;
import jgi.ReformatReadsOld;
import parse.Parse;
import processor.ReformatProcessor;

/**
 * Single-thread timing of actual legacy String sanitizers on two ASCII profiles.
 * Allowed headers retain identity; dirty headers require literal mapped values.
 * Each profile has64 distinct128-char headers. Value, identity, modern and Unicode
 * controls run before timing. Every timed result contributes its hash and a ring
 * reference. Use fresh matched JVMs and counts giving at least ten CPU seconds;
 * short calibration does not support performance claims.
 * @author Sucrose
 */
public final class HeaderSanitizerBenchmark{

	/** Requires phase, target and profile; iterations default to one million, capped at one billion. */
	public static void main(final String[] args){
		String phase=null, target=null, profile=null;
		long iterations=1000000L;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException(arg);}
			if(pair[0].equalsIgnoreCase("phase")){phase=pair[1];}
			else if(pair[0].equalsIgnoreCase("target")){target=pair[1];}
			else if(pair[0].equalsIgnoreCase("profile")){profile=pair[1];}
			else if(pair[0].equalsIgnoreCase("iterations")){iterations=Parse.parseKMG(pair[1]);}
			else{throw new IllegalArgumentException(arg);}
		}
		if(!"baseline".equals(phase) && !"fixed".equals(phase)){
			throw new IllegalArgumentException("Set phase=baseline|fixed.");
		}
		if(!"legacy".equals(target) && !"legacyold".equals(target)){
			throw new IllegalArgumentException("Set target=legacy|legacyold.");
		}
		if(!"allowed".equals(profile) && !"dirty".equals(profile)){
			throw new IllegalArgumentException("Set profile=allowed|dirty.");
		}
		if(iterations<1 || iterations>1000000000L){
			throw new IllegalArgumentException("iterations must be in [1,1000000000]: "+iterations);
		}
		final Charset charset=Charset.defaultCharset();
		require(charset.equals(StandardCharsets.UTF_8), "The fixed ASCII-byte String controls require UTF-8: "+charset);
		require(Runtime.getRuntime().availableProcessors()==1,
			"Matched single-CPU JVMs require -XX:ActiveProcessorCount=1, including compiler and GC work.");
		final String unicode=controls(target, phase.equals("baseline"));
		System.out.println("sanitizer_controls\tencoding="+charset.name()+"\tunicode="+unicode+
			"\tmodern_unicode=A_B\tnull=true\tempty_identity=true\tallowed_identity=true");

		final boolean dirty=profile.equals("dirty");
		final String inputBody=dirty ? body(DIRTY_BLOCK, DIRTY_TAIL) : body(ALLOWED_BLOCK, ALLOWED_TAIL);
		final String expectedBody=dirty ? body(EXPECTED_BLOCK, EXPECTED_TAIL) : inputBody;
		final String[] headers=new String[HEADER_COUNT], results=new String[HEADER_COUNT];
		final int[] expectedHashes=new int[HEADER_COUNT];
		long cycleChecksum=0;
		for(int i=0; i<HEADER_COUNT; i++){
			final String prefix=(i<10 ? "0" : "")+i;
			final String header=prefix+inputBody;
			final String expected=prefix+expectedBody;
			require(header.length()==128 && expected.length()==128 && ascii(header) && ascii(expected),
				"Both profiles must time exactly 128 ASCII chars per conversion: "+i);
			final String actual=convert(target, header);
			final String modern=ReformatProcessor.fixHeader(header);
			require(expected.equals(actual) && (actual==header)==!dirty,
				"Selected sanitizer must match independent literal value/identity controls: "+target+"/"+profile+"/"+i);
			require(expected.equals(modern) && (modern==header)==!dirty,
				"Modern counterpart must match the same independent profile controls: "+profile+"/"+i);
			headers[i]=header;
			expectedHashes[i]=expected.hashCode();
			cycleChecksum+=expectedHashes[i];
		}
		long expectedChecksum=(iterations/HEADER_COUNT)*cycleChecksum;
		for(int i=0; i<iterations%HEADER_COUNT; i++){expectedChecksum+=expectedHashes[i];}
		System.out.println("sanitizer_pool\tprofile="+profile+"\theaders=64\tlength=128\tidentity="+!dirty+
			"\texpected_checksum="+expectedChecksum);
		//Select a direct static loop before its timer starts; no per-conversion dispatch.
		final long checksum=target.equals("legacy") ? timeLegacy(headers, results, iterations) : timeLegacyOld(headers, results, iterations);
		resultSink=results;
		checksumSink=checksum;
		require(checksum==expectedChecksum, "Every timed conversion must match the independent expected checksum: "+
			checksum+" != "+expectedChecksum);
		System.out.println("HEADER_SANITIZER_BENCHMARK_PASS phase="+phase+" target="+target+" profile="+profile+
			" iterations="+iterations+" checksum="+checksum+" loop_elapsed_ns="+elapsedNanos);
	}

	/** Checks phase, null and original-reference behavior without timing any Unicode exception. */
	private static String controls(final String target, final boolean baseline){
		final String empty=new String(""), allowed=new String("AZaz09_.#-()~");
		require(convert(target, null)==null, "Null sanitizer input must return null.");
		require(convert(target, empty)==empty, "Empty sanitizer input must retain its original reference.");
		require(convert(target, allowed)==allowed, "Allowed ASCII must retain its original reference.");
		final String input="A\u0080B";
		require("A_B".equals(ReformatProcessor.fixHeader(input)), "The modern non-ASCII control must map char128 to underscore.");
		String actual=null;
		boolean failed=false;
		try{actual=convert(target, input);}
		catch(ArrayIndexOutOfBoundsException e){
			require(baseline, "Fixed sanitizer must not index outside the ASCII table.");
			final StackTraceElement top=e.getStackTrace()[0];
			final String owner=target.equals("legacy") ? "jgi.ReformatReads" : "jgi.ReformatReadsOld";
			require(top.getClassName().equals(owner) && top.getMethodName().equals("fixHeader"),
				"Expected array-index failure must originate in the selected legacy sanitizer: "+top);
			failed=true;
		}
		require(failed==baseline, "The selected source must match phase=baseline|fixed before timing.");
		if(!baseline){require("A_B".equals(actual), "Fixed sanitizer must replace only the non-ASCII char in the setup control.");}
		return failed ? "array-index-failure" : actual;
	}

	/** Dispatches only setup controls; both timed methods call their actual helper directly. */
	private static String convert(final String target, final String header){
		assert(target.equals("legacy") || target.equals("legacyold")) :
			"Setup and timed calls must use the same validated legacy target: "+target;
		return target.equals("legacy") ? ReformatReads.fixHeader(header) : ReformatReadsOld.fixHeader(header);
	}

	/** Two prefix digits plus seven16-char blocks and a14-char tail make128 chars. */
	private static String body(final String block, final String tail){
		assert(block.length()==16 && tail.length()==14) :
			"The fixed 128-char profiles require 2+7*16+14 chars before timing.";
		final StringBuilder out=new StringBuilder(126);
		for(int i=0; i<7; i++){out.append(block);}
		return out.append(tail).toString();
	}

	/** Refuses accidental non-ASCII workload data; Unicode belongs only in setup controls. */
	private static boolean ascii(final String value){
		assert(value!=null) : "Each setup header and independent expected value must be present.";
		for(int i=0; i<value.length(); i++){if(value.charAt(i)>=128){return false;}}
		return true;
	}

	/** Times conversion, every-result hashing and ring retention for the first legacy helper. */
	private static long timeLegacy(final String[] headers, final String[] results, final long iterations){
		assert(headers.length==HEADER_COUNT && results.length==HEADER_COUNT) :
			"The checksum cycle and masked ring need 64 input/result slots.";
		long checksum=0;
		final long start=System.nanoTime();
		for(long i=0; i<iterations; i++){
			final int slot=(int)(i&HEADER_MASK);
			final String result=ReformatReads.fixHeader(headers[slot]);
			results[slot]=result;
			checksum+=result.hashCode();
		}
		elapsedNanos=System.nanoTime()-start;
		return checksum;
	}

	/** Uses exactly the same timed work for the older legacy helper. */
	private static long timeLegacyOld(final String[] headers, final String[] results, final long iterations){
		assert(headers.length==HEADER_COUNT && results.length==HEADER_COUNT) :
			"The checksum cycle and masked ring need 64 input/result slots.";
		long checksum=0;
		final long start=System.nanoTime();
		for(long i=0; i<iterations; i++){
			final int slot=(int)(i&HEADER_MASK);
			final String result=ReformatReadsOld.fixHeader(headers[slot]);
			results[slot]=result;
			checksum+=result.hashCode();
		}
		elapsedNanos=System.nanoTime()-start;
		return checksum;
	}

	/** Keeps correctness checks active independently of the JVM assertion setting. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static volatile String[] resultSink;
	private static volatile long checksumSink;
	private static long elapsedNanos;
	private static final int HEADER_COUNT=64, HEADER_MASK=HEADER_COUNT-1;
	private static final String ALLOWED_BLOCK="read_ID.09#-()~_", ALLOWED_TAIL="read_ID.09#-()";
	private static final String DIRTY_BLOCK="read ID/09:\t?\u007f_X", DIRTY_TAIL="read ID/09:\t?\u007f";
	private static final String EXPECTED_BLOCK="read_ID_09_____X", EXPECTED_TAIL="read_ID_09____";
}
