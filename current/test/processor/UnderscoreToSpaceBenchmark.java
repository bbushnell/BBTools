package test.processor;

import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;

import jgi.ReformatReads;
import jgi.ReformatReadsOld;
import parse.Parse;
import processor.ReformatProcessor;

/** Single-thread timing of actual last-underscore conversion on a fixed ASCII header pool.
 * target=processor is the unchanged default; legacy and legacyold select the two old drivers.
 * Includes result hashing and ring retention in the measured loop on both versions.
 * Unicode correctness and no-underscore identity controls run before timing.
 * Use fresh matched JVMs; choose iterations giving at least ten seconds of CPU first.
 * @author Sucrose
 */
public final class UnderscoreToSpaceBenchmark{

	/** Requires phase=baseline|fixed; optional target defaults to processor, iterations to ten million (cap one billion). */
	public static void main(final String[] args){
		String phase=null;
		String target="processor";
		long iterations=10000000L;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException(arg);}
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("target")){target=pair[1];}
			else if(pair[0].equals("iterations")){iterations=Parse.parseKMG(pair[1]);}
			else{throw new IllegalArgumentException(arg);}
		}
		if(!"baseline".equals(phase) && !"fixed".equals(phase)){
			throw new IllegalArgumentException("Set phase=baseline|fixed.");
		}
		if(!"processor".equals(target) && !"legacy".equals(target) && !"legacyold".equals(target)){
			throw new IllegalArgumentException("Set target=processor|legacy|legacyold.");
		}
		if(iterations<1 || iterations>1000000000L){
			throw new IllegalArgumentException("iterations must be in [1,1000000000]: "+iterations);
		}
		final Charset charset=Charset.defaultCharset();
		if(!charset.equals(StandardCharsets.UTF_8)){
			throw new IllegalStateException("The baseline Unicode control requires UTF-8, observed "+charset);
		}
		final String unicode="\u00e9_a_b";
		final String unicodeResult=convert(target, unicode);
		final String unicodeExpected=phase.equals("baseline") ? "\u00e9_ _b" : "\u00e9_a b";
		require(unicodeExpected.equals(unicodeResult), "Unexpected Unicode control: "+unicodeResult);
		final String plain="plain\u00e9";
		require(convert(target, plain)==plain,
			"The no-underscore path must return its input reference without conversion.");
		if(!target.equals("processor")){System.out.println("header_target\ttarget="+target);}
		System.out.println("header_controls\tencoding="+charset.name()+"\tunicode="+unicodeResult+"\tno_underscore_identity=true");

		final String[] headers=new String[HEADER_COUNT], results=new String[HEADER_COUNT];
		final int[] expectedHashes=new int[HEADER_COUNT];
		long cycleChecksum=0;
		for(int i=0; i<HEADER_COUNT; i++){
			final String prefix="instrument:42:flowcell:"+i+"_read";
			headers[i]=prefix+"_"+i;
			final String expected=prefix+" "+i;
			require(expected.equals(convert(target, headers[i])),
				"ASCII control must replace only the final underscore: "+headers[i]);
			expectedHashes[i]=expected.hashCode();
			cycleChecksum+=expectedHashes[i];
		}
		long expectedChecksum=(iterations/HEADER_COUNT)*cycleChecksum;
		for(int i=0; i<iterations%HEADER_COUNT; i++){expectedChecksum+=expectedHashes[i];}

		// Dispatch before the selected method starts its timer; timed loops call one static helper directly.
		final long checksum;
		if(target.equals("legacy")){checksum=timeLegacy(headers, results, iterations);}
		else if(target.equals("legacyold")){checksum=timeLegacyOld(headers, results, iterations);}
		else{checksum=timeProcessor(headers, results, iterations);}
		resultSink=results;
		checksumSink=checksum;
		require(checksum==expectedChecksum, "Every conversion must contribute to the literal-oracle checksum: "+
			checksum+" != "+expectedChecksum);
		System.out.println("UNDERSCORE_BENCHMARK_PASS phase="+phase+" iterations="+iterations+
			" checksum="+checksum+" loop_elapsed_ns="+elapsedNanos);
	}

	/** Selects the same helper for untimed Unicode, identity and literal ASCII controls. */
	private static String convert(final String target, final String header){
		assert(target.equals("processor") || target.equals("legacy") || target.equals("legacyold")) :
			"Setup controls must use the same validated target as the selected timed loop: "+target;
		if(target.equals("legacy")){return ReformatReads.underscoreToSpace(header);}
		if(target.equals("legacyold")){return ReformatReadsOld.underscoreToSpace(header);}
		return ReformatProcessor.underscoreToSpace(header);
	}

	/** Times only direct processor conversion, hashing and retention, preserving the default workload. */
	private static long timeProcessor(final String[] headers, final String[] results, final long iterations){
		assert(headers.length==HEADER_COUNT && results.length==HEADER_COUNT) :
			"The shared checksum pool and masked ring indexing require 64 input/result slots.";
		long checksum=0;
		final long start=System.nanoTime();
		for(long i=0; i<iterations; i++){
			final int slot=(int)(i&HEADER_MASK);
			final String result=ReformatProcessor.underscoreToSpace(headers[slot]);
			results[slot]=result;
			checksum+=result.hashCode();
		}
		elapsedNanos=System.nanoTime()-start;
		return checksum;
	}

	/** Uses the identical ASCII loop with the first legacy driver's actual static helper. */
	private static long timeLegacy(final String[] headers, final String[] results, final long iterations){
		assert(headers.length==HEADER_COUNT && results.length==HEADER_COUNT) :
			"The shared checksum pool and masked ring indexing require 64 input/result slots.";
		long checksum=0;
		final long start=System.nanoTime();
		for(long i=0; i<iterations; i++){
			final int slot=(int)(i&HEADER_MASK);
			final String result=ReformatReads.underscoreToSpace(headers[slot]);
			results[slot]=result;
			checksum+=result.hashCode();
		}
		elapsedNanos=System.nanoTime()-start;
		return checksum;
	}

	/** Uses the identical ASCII loop with the older legacy driver's actual static helper. */
	private static long timeLegacyOld(final String[] headers, final String[] results, final long iterations){
		assert(headers.length==HEADER_COUNT && results.length==HEADER_COUNT) :
			"The shared checksum pool and masked ring indexing require 64 input/result slots.";
		long checksum=0;
		final long start=System.nanoTime();
		for(long i=0; i<iterations; i++){
			final int slot=(int)(i&HEADER_MASK);
			final String result=ReformatReadsOld.underscoreToSpace(headers[slot]);
			results[slot]=result;
			checksum+=result.hashCode();
		}
		elapsedNanos=System.nanoTime()-start;
		return checksum;
	}

	/** Fails loudly independently of assertion settings. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Retains produced objects and the checksum beyond the measured loop. */
	private static volatile String[] resultSink;
	private static volatile long checksumSink;
	/** Written after each selected timer stops, then reported by the same main thread. */
	private static long elapsedNanos;
	private static final int HEADER_COUNT=64;
	private static final int HEADER_MASK=HEADER_COUNT-1;
}
