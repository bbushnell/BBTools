package test.icecream;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import icecream.ZMW;
import icecream.ZMWStreamer;
import parse.Parse;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;

/** Matched same-movie grouping workload through the real FASTA/CRIS or SAM backend.
 * Preparation writes reusable inputs once; timing mode only reads an existing input.
 * Fixed passes permit external process-CPU measurement across all backend threads.
 * Counts and ordered header/base/group checksums are verified after every pass.
 * Checksums intentionally wrap modulo 2^64. External process CPU includes startup,
 * oracle setup and warmup; elapsed_ns covers only the fixed measured passes.
 * A separate preparetrim mode writes compressed input for end-trimming CLI timings.
 * The grouping benchmark measures the complete input/grouping/consumption path, including its oracle,
 * rather than isolated movie-comparison instructions or a complete analysis pipeline.
 * @author Brian Bushnell, Sucrose, Shinobu
 */
public final class ZMWStreamerTimingTest{

	/** Creates deterministic fixtures or consumes one fixture a fixed number of times.
	 * @param args mode=prepare outdir=new-directory, or mode=bench in=file backend=cris|streamer;
	 * optional wells=8192 subreads=4 passes=8 warmup=1; preparetrim requires subreads=1
	 * @throws Exception If fixture creation or producer joining fails
	 */
	public static void main(String[] args) throws Exception{
		if(args.length==0){
			System.err.println("mode=prepare outdir=new-directory [wells=8192 subreads=4]\n"+
					"mode=preparetrim outdir=new-directory wells=count subreads=1\n"+
					"mode=bench in=file backend=cris|streamer [wells=8192 subreads=4 passes=8 warmup=1]");
			return;
		}
		String mode="bench", backend="cris";
		Path input=null, outdir=null;
		int wells=8192, subreads=4, passes=8, warmup=1;
		for(String arg : args){
			final int equals=arg.indexOf('=');
			require(equals>0, "Expected flag=value: "+arg);
			final String key=arg.substring(0, equals).toLowerCase(Locale.ROOT);
			final String value=arg.substring(equals+1);
			if(key.equals("mode")){mode=value.toLowerCase(Locale.ROOT);}
			else if(key.equals("backend")){backend=value.toLowerCase(Locale.ROOT);}
			else if(key.equals("in")){input=Paths.get(value);}
			else if(key.equals("outdir")){outdir=Paths.get(value);}
			else if(key.equals("wells")){wells=Parse.parseIntKMG(value);}
			else if(key.equals("subreads")){subreads=Parse.parseIntKMG(value);}
			else if(key.equals("passes")){passes=Parse.parseIntKMG(value);}
			else if(key.equals("warmup")){warmup=Parse.parseIntKMG(value);}
			else{throw new IllegalArgumentException("Unknown fixture argument: "+arg);}
		}
		final boolean trimInput=mode.equals("preparetrim");
		require(wells>0 && subreads>0 && subreads<=16 &&
				(trimInput ? subreads==1 && wells<=20000000 : (long)wells*subreads<=262144),
				"Grouping input allows262144 records; trimming input requires subreads=1 and at most20000000 wells");
		require(passes>0 && warmup>=0, "passes must be positive and warmup nonnegative");
		require(backend.equals("cris") || backend.equals("streamer"), "backend must be cris or streamer");
		require(ZMWStreamer.class.desiredAssertionStatus(), "Production assertions must be enabled");
		Shared.setThreads(1);
		if(trimInput){
			require(outdir!=null && input==null, "Preparation requires outdir=new-directory and no in argument");
			Files.createDirectory(outdir);
			writeTrimInput(outdir.resolve("heavy-trim.fa.gz"), wells);
			System.out.println("wells\treads\tinput_bases\ttrimmed_output_bases");
			System.out.println(wells+"\t"+wells+"\t"+(200L*wells)+"\t"+(80L*wells));
			return;
		}
		if(mode.equals("prepare")){
			require(outdir!=null && input==null, "Preparation requires outdir=new-directory and no in argument");
			Files.createDirectory(outdir);
			writeInput(outdir.resolve("same-movie.fa"), false, wells, subreads);
			writeInput(outdir.resolve("same-movie.sam"), true, wells, subreads);
			System.out.println("wells\treads\tbases\tchecksum");
			final Observation expected=expected(wells, subreads);
			System.out.println(expected.wells+"\t"+expected.reads+"\t"+expected.bases+"\t"+expected.checksum);
			return;
		}
		require(mode.equals("bench"), "mode must be prepare or bench");
		require(input!=null && outdir==null && Files.isRegularFile(input),
				"Timing requires in=existing-regular-file and no outdir argument");
		final boolean sam=backend.equals("streamer");
		final String name=input.toString().toLowerCase(Locale.ROOT);
		require(name.endsWith(sam ? ".sam" : ".fa"), "Use the prepared uncompressed .fa or .sam input for this backend");
		final Observation expected=expected(wells, subreads);
		final boolean originalStreamer=ZMWStreamer.useStreamer;
		long checksum=0;
		long elapsed;
		try{
			ZMWStreamer.useStreamer=sam;
			for(int i=0; i<warmup; i++){verify(consume(input, sam), expected, "warmup", i);}
			final long start=System.nanoTime();
			for(int i=0; i<passes; i++){
				final Observation observed=consume(input, sam);
				verify(observed, expected, "measured", i);
				checksum+=observed.checksum;
			}
			elapsed=System.nanoTime()-start;
		}finally{ZMWStreamer.useStreamer=originalStreamer;}
		System.out.println("backend\twells_per_pass\treads_per_pass\tbases_per_pass\tpasses\twarmup\texpected_checksum\tchecksum\telapsed_ns");
		System.out.println(backend+"\t"+expected.wells+"\t"+expected.reads+"\t"+expected.bases+"\t"+
				passes+"\t"+warmup+"\t"+expected.checksum+"\t"+checksum+"\t"+elapsed);
	}

	/** Writes identical names/bases in both formats; movie prefix never changes.
	 * The writer copies the reused builder before it is cleared for the next record.
	 */
	private static void writeInput(Path path, boolean sam, int wells, int subreads){
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(), false, false, false);
		writer.start();
		boolean error;
		try{
			if(sam){writer.print("@HD\tVN:1.6\tSO:unknown\n");}
			final ByteBuilder row=new ByteBuilder(256);
			for(int well=1; well<=wells; well++){
				for(int sub=0; sub<subreads; sub++){
					row.clear();
					if(!sam){row.append('>');}
					row.append(MOVIE).append('/').append(well).append('/').append(sub*BASES.length).
							append('_').append((sub+1)*BASES.length);
					if(sam){row.append("\t4\t*\t0\t0\t*\t*\t0\t0\t").append(BASES).append('\t').append(QUALITIES).append('\n');}
					else{row.append('\n').append(BASES).append('\n');}
					writer.print(row);
				}
			}
		}finally{error=writer.poisonAndWait();}
		require(!error, "Could not finish deterministic input: "+path);
	}

	/** Writes one200-base read per well:80C followed by120N, with coordinates100_300.
	 * Gzip keeps the reusable timing input small; no uncompressed copy is retained.
	 */
	private static void writeTrimInput(Path path, int wells){
		final byte[] bases=new byte[200];
		Arrays.fill(bases, 0, 80, (byte)'C');
		Arrays.fill(bases, 80, 200, (byte)'N');
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(), false, false, false);
		writer.start();
		boolean error;
		try{
			final ByteBuilder row=new ByteBuilder(256);
			for(int well=1; well<=wells; well++){
				row.clear();
				row.append('>').append(MOVIE).append('/').append(well).append("/100_300\n");
				row.append(bases).append('\n');
				writer.print(row);
			}
		}finally{error=writer.poisonAndWait();}
		require(!error, "Could not finish trimming input: "+path);
	}

	/** Drains through the normal terminal and joins this pass's producer before checking counters.
	 * No per-read strings or retained well lists are allocated by the consumer oracle.
	 */
	private static Observation consume(Path input, boolean sam) throws InterruptedException{
		final FileFormat ff=FileFormat.testInput(input.toString(), sam ? FileFormat.SAM : FileFormat.FASTA, null, false, true);
		final ZMWStreamer streamer=new ZMWStreamer(ff, 16, -1, -1);
		final Thread producer=streamer.runStreamer(true);
		final Observation observed=new Observation();
		for(ZMW well=streamer.nextZMW(); well!=null; well=streamer.nextZMW()){
			observed.wells++;
			observed.checksum=observed.checksum*MIX+well.size();
			for(Read read : well){
				observed.reads++;
				observed.bases+=read.length();
				observed.checksum=observed.checksum*MIX+(read.id==null ? 0 : read.id.hashCode());
				observed.checksum=observed.checksum*MIX+baseHash(read.bases);
			}
		}
		producer.join();
		require(!streamer.errorState && streamer.ZMWs==observed.wells,
				"Input error or published/consumed well count mismatch: "+input+", "+streamer.ZMWs+"/"+observed.wells);
		return observed;
	}

	/** Computes the layout oracle without parsing headers or invoking production grouping.
	 * Decimal digits are hashed arithmetically, independently of ByteBuilder formatting.
	 */
	private static Observation expected(int wells, int subreads){
		final Observation expected=new Observation();
		expected.wells=wells;
		expected.reads=(long)wells*subreads;
		expected.bases=expected.reads*BASES.length;
		final int prefixHash=31*MOVIE.hashCode()+'/';
		final int sequenceHash=baseHash(BASES);
		for(int well=1; well<=wells; well++){
			expected.checksum=expected.checksum*MIX+subreads;
			final int wellHash=31*decimalHash(prefixHash, well)+'/';
			for(int sub=0; sub<subreads; sub++){
				final int headerHash=decimalHash(31*decimalHash(wellHash, sub*BASES.length)+'_', (sub+1)*BASES.length);
				expected.checksum=expected.checksum*MIX+headerHash;
				expected.checksum=expected.checksum*MIX+sequenceHash;
			}
		}
		return expected;
	}

	/** Extends Java's String hash with the nonnegative value's decimal digits, without a String. */
	private static int decimalHash(int hash, int value){
		assert(value>=0) : "Fixture coordinates and positive well IDs must be nonnegative for decimal formatting: "+value;
		int divisor=1;
		while(value/divisor>=10){divisor*=10;}
		for(; divisor>0; divisor/=10){hash=31*hash+'0'+value/divisor%10;}
		return hash;
	}

	/** Hashes every base; a missing sequence differs from the nonempty literal sequence. */
	private static int baseHash(byte[] bases){
		int hash=1;
		if(bases==null){return 0;}
		for(byte base : bases){hash=31*hash+(base&255);}
		return hash;
	}

	/** Checks each pass independently, preventing compensating errors across repetitions. */
	private static void verify(Observation actual, Observation expected, String phase, int pass){
		require(actual.wells==expected.wells && actual.reads==expected.reads && actual.bases==expected.bases &&
				actual.checksum==expected.checksum, "Unexpected "+phase+" pass "+pass+": wells="+actual.wells+
				", reads="+actual.reads+", bases="+actual.bases+", checksum="+actual.checksum+
				"; expected "+expected.wells+", "+expected.reads+", "+expected.bases+", "+expected.checksum);
	}

	/** Reports fixture/oracle failures even when a caller accidentally disables assertions. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Primitive counters; one instance per complete pass, never per read. */
	private static final class Observation{
		long wells, reads, bases, checksum;
	}

	private static final long MIX=1000003;
	private static final String MOVIE="m64000_200101_000000_fixture_same_movie";
	private static final byte[] BASES="ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] QUALITIES=new byte[BASES.length];
	static{Arrays.fill(QUALITIES, (byte)'I');}
	private ZMWStreamerTimingTest(){}
}
