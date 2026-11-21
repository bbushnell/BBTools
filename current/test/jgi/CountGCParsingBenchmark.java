package test.jgi;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;

import jgi.CountGC;
import parse.Parse;
import structures.ByteBuilder;

/** Matched in-memory timing of actual CountGC readers with two equal-base ASCII shapes.
 * Exact formatted output is checked before timing; summary parsing, reader allocations
 * and six-bucket checks are timed together. Generated inputs are small and reused.
 * Use fresh matched JVMs and calibrate each format/shape to at least ten CPU seconds.
 * @author Sucrose
 */
public final class CountGCParsingBenchmark{

	/** Requires format=fasta|fastq and shape=short|long; passes defaults to ten thousand. */
	public static void main(final String[] args) throws Exception{
		String format=null, shape=null;
		long passes=10000;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException(arg);}
			final String key=pair[0].toLowerCase(Locale.ROOT);
			if(key.equals("format")){format=pair[1];}
			else if(key.equals("shape")){shape=pair[1];}
			else if(key.equals("passes")){passes=Parse.parseKMG(pair[1]);}
			else{throw new IllegalArgumentException(arg);}
		}
		if(!"fasta".equals(format) && !"fastq".equals(format)){
			throw new IllegalArgumentException("Set format=fasta|fastq.");
		}
		if(!"short".equals(shape) && !"long".equals(shape)){
			throw new IllegalArgumentException("Set shape=short|long.");
		}
		if(passes<1 || passes>1000000000L){
			throw new IllegalArgumentException("passes must be in [1,1000000000]: "+passes);
		}
		final boolean fasta=format.equals("fasta"), longShape=shape.equals("long");
		final int records=longShape ? 4 : 256, length=longShape ? 9600 : 150;
		final byte[] input=makeInput(fasta, longShape);
		final int expectedBytes=fasta ? (longShape ? 38428 : 40192) : (longShape ? 76840 : 79360);
		require(input.length==expectedBytes, "Generated input size differs: "+input.length+" != "+expectedBytes);
		final long[] expected={9472, 9472, 9472, 9472, 512, fasta ? records : 0};
		final long expectedPerPass=fasta ? (longShape ? 97304 : 98816) : 97280;
		final String row="read\t"+length+"\t0.25000\t0.25000\t0.25000\t0.25000\t0.01333\t0.50000\n";
		final StringBuilder expectedOutput=new StringBuilder();
		for(int i=0; i<records; i++){expectedOutput.append(row);}
		if(fasta){expectedOutput.append('\n');}

		final int oldFormat=CountGC.FORMAT;
		final boolean oldSummary=CountGC.SUMMARY_ONLY;
		final PrintStream oldOut=System.out;
		try{
			CountGC.FORMAT=1;
			CountGC.SUMMARY_ONLY=false;
			final ByteArrayOutputStream capturedBytes=new ByteArrayOutputStream();
			final PrintStream captured=new PrintStream(capturedBytes, true, StandardCharsets.UTF_8.name());
			try{
				System.setOut(captured);
				requireCounts(count(fasta, new ByteArrayInputStream(input)), expected);
			}finally{
				System.setOut(oldOut);
				captured.close();
			}
			require(expectedOutput.toString().equals(capturedBytes.toString(StandardCharsets.UTF_8.name())),
				"Formatted output differs for "+format+"/"+shape);
			oldOut.println("countgc_controls format="+format+" shape="+shape+" input_bytes="+input.length+
				" records="+records+" bases=38400 output_exact=true counts="+Arrays.toString(expected));

			CountGC.SUMMARY_ONLY=true;
			final PrintStream silent=new PrintStream(new NoOutput());
			final ByteArrayInputStream stream=new ByteArrayInputStream(input);
			long checksum=0, elapsed;
			try{
				System.setOut(silent);
				final long start=System.nanoTime();
				for(long pass=0; pass<passes; pass++){
					stream.reset();
					final long[] counts=count(fasta, stream);
					requireCounts(counts, expected);
					for(int i=0; i<counts.length; i++){checksum+=(i+1)*counts[i];}
				}
				elapsed=System.nanoTime()-start;
			}finally{
				System.setOut(oldOut);
				silent.close();
			}
			checksumSink=checksum;
			require(checksum==passes*expectedPerPass,
				"Every returned bucket must contribute to checksum: "+checksum+" != "+passes*expectedPerPass);
			oldOut.println("COUNTGC_PARSING_BENCHMARK_PASS format="+format+" shape="+shape+
				" passes="+passes+" checksum="+checksum+" nanos="+elapsed);
		}finally{
			System.setOut(oldOut);
			CountGC.FORMAT=oldFormat;
			CountGC.SUMMARY_ONLY=oldSummary;
		}
	}

	/** Builds 38400 bases in either 256 short records or four long records, outside timing. */
	private static byte[] makeInput(final boolean fasta, final boolean longShape){
		final ByteBuilder unit=new ByteBuilder();
		for(int i=0; i<37; i++){unit.append("ACGT");}
		unit.append("NN");
		final byte[] bases=unit.toBytes();
		final int records=longShape ? 4 : 256, units=longShape ? 64 : 1;
		final ByteBuilder input=new ByteBuilder();
		for(int record=0; record<records; record++){
			input.append(fasta ? '>' : '@').append("read").nl();
			for(int i=0; i<units; i++){input.append(bases);}
			input.nl();
			if(!fasta){
				input.append('+').nl();
				for(int i=0; i<units*150; i++){input.append('I');}
				input.nl();
			}
		}
		return input.toBytes();
	}

	/** Calls the production reader; the supplied stream and its immutable backing input are borrowed. */
	private static long[] count(final boolean fasta, final ByteArrayInputStream input) throws Exception{
		return fasta ? CountGC.countFasta(input, null) : CountGC.countFastq(input, null);
	}

	/** Checks all six independently predicted buckets on every measured pass. */
	private static void requireCounts(final long[] actual, final long[] expected){
		if(!Arrays.equals(actual, expected)){
			throw new AssertionError("Unexpected base counts: "+Arrays.toString(actual));
		}
	}

	/** Fails loudly even if assertions were disabled. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Summary-mode parsing must emit no output, including on the measured path. */
	private static final class NoOutput extends OutputStream{
		@Override
		public void write(final int value){throw new AssertionError("Summary parsing wrote output.");}
	}

	private static volatile long checksumSink;
}
