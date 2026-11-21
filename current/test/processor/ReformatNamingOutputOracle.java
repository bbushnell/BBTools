package test.processor;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.security.MessageDigest;

import structures.ByteBuilder;

/**
 * Independent byte/hash oracle for the frozen sixteen-read naming benchmark.
 * The first output block is literal; later blocks use closed ordinal formulas
 * for this exact input order, rather than invoking or reproducing the allocator.
 * Only a small specification file is written. No FASTQ, output-name registry,
 * production processor, driver or timing work is created or invoked here.
 * The runner supplies trd=f, default identifier scope and unchanged normalization.
 * @author Sucrose
 */
public final class ReformatNamingOutputOracle{

	/** Requires mode=spec blocks=1..1048576 out=NEW_SPEC; refuses to replace an existing file. */
	public static void main(final String[] args) throws Exception{
		String mode=null, output=null;
		int blocks=-1;
		for(final String arg : args){
			final int equals=arg.indexOf('=');
			if(equals<1){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String key=arg.substring(0, equals), value=arg.substring(equals+1);
			if(key.equalsIgnoreCase("mode")){mode=value;}
			else if(key.equalsIgnoreCase("blocks")){blocks=Integer.parseInt(value);}
			else if(key.equalsIgnoreCase("out")){output=value;}
			else{throw new IllegalArgumentException("Unknown oracle option: "+arg);}
		}
		if(!"spec".equals(mode) || blocks<1 || blocks>1048576 || output==null || output.isEmpty()){
			throw new IllegalArgumentException("Set mode=spec blocks=1..1048576 out=NEW_SPEC.");
		}
		final Path destination=Paths.get(output);
		if(Files.exists(destination)){throw new IllegalArgumentException("Oracle output already exists: "+destination);}
		final ByteBuilder input=new ByteBuilder(512), enabled=new ByteBuilder(512);
		literalBlock(input, INPUT);
		literalBlock(enabled, FIRST);
		check(input.length()==271, "The frozen sixteen-read input block must contain271 bytes.");
		check(enabled.length()==289, "The literal first fixed block must contain289 bytes.");
		final MessageDigest inputDigest=MessageDigest.getInstance("SHA-256");
		final MessageDigest enabledDigest=MessageDigest.getInstance("SHA-256");
		long enabledBytes=0;
		for(int block=0; block<blocks; block++){
			inputDigest.update(input.array, 0, input.length());
			if(block>0){laterBlock(enabled, block);}
			enabledDigest.update(enabled.array, 0, enabled.length());
			enabledBytes+=enabled.length();
		}
		final long reads=16L*blocks, bases=4L*reads, inputBytes=271L*blocks;
		final String inputHash=sha80(inputDigest.digest()), enabledHash=sha80(enabledDigest.digest());
		final String spec="phase\tfixed\nblocks\t"+blocks+"\nreads\t"+reads+"\nbases\t"+bases+
			"\ninput_bytes\t"+inputBytes+"\ninput_sha80\t"+inputHash+
			"\nexpected_off_bytes\t"+inputBytes+"\nexpected_off_sha80\t"+inputHash+
			"\nexpected_on_bytes\t"+enabledBytes+"\nexpected_on_sha80\t"+enabledHash+
			"\ntrd\tfalse\nuniquenamescope\tidentifier\n";
		Files.write(destination, spec.getBytes(StandardCharsets.US_ASCII),
			StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE);
		System.out.println("REFORMAT_NAMING_OUTPUT_ORACLE_PASS blocks="+blocks+" reads="+reads+
			" bases="+bases+" input_bytes="+inputBytes+" enabled_bytes="+enabledBytes);
	}

	/** Appends one literal ASCII block; every record has the same ACGT/IIII payload. */
	private static void literalBlock(final ByteBuilder bytes, final String[] headers){
		check(headers.length==16, "Each frozen benchmark block must contain sixteen reads.");
		bytes.clear();
		for(final String header : headers){bytes.append('@').append(header).append(PAYLOAD);}
	}

	/** Formats block b>=1 using its sixteen independent slot formulas.
	 * Each completed block advances x and x_2 by three ordinals, p and n by two,
	 * and the remaining natural input keys by one. The first block's natural x_2,
	 * x_2_2, n_2, p_2/1 and q_2 inputs account for the offsets below.
	 */
	private static void laterBlock(final ByteBuilder bytes, final int b){
		assert(b>0 && b<1048576) : "Closed later-block formulas require zero-based blocks1..1048575; b="+b;
		bytes.clear();
		record(bytes, "x_", 3*b+1, "");
		record(bytes, "x_", 3*b+2, "");
		record(bytes, "x_2_", 3*b+2, "");
		record(bytes, "x_", 3*b+3, "");
		record(bytes, "x_02_", b+1, "");
		record(bytes, "x_2_", 3*b+3, "");
		record(bytes, "x_2_2_", b+2, "");
		record(bytes, "x_2_", 3*b+4, "");
		record(bytes, "n_", 2*b+2, " a");
		record(bytes, "n_2_", b+1, " b");
		record(bytes, "n_", 2*b+3, " c");
		record(bytes, "p_", 2*b+1, "/1");
		record(bytes, "p_", 2*b+2, "/1");
		record(bytes, "p_2_", b+2, "/1");
		record(bytes, "q_", b+2, "");
		record(bytes, "q_2_", b+1, "");
	}

	/** Appends primitive ordinals directly to one reused byte buffer, preserving literal descriptions/mate tags. */
	private static void record(final ByteBuilder bytes, final String prefix, final int ordinal, final String tail){
		assert(ordinal>=2) : "Repeated-block generated ordinals start at2; ordinal="+ordinal;
		bytes.append('@').append(prefix).append(ordinal).append(tail).append(PAYLOAD);
	}

	/** Formats only the final80 digest bits, matching the retained benchmark specification. */
	private static String sha80(final byte[] digest){
		assert(digest.length==32) : "The oracle requests SHA-256 before retaining its final80 bits.";
		final StringBuilder text=new StringBuilder(20);
		for(int i=digest.length-10; i<digest.length; i++){
			text.append(Character.forDigit((digest[i]>>>4)&15, 16)).append(Character.forDigit(digest[i]&15, 16));
		}
		return text.toString();
	}

	/** Literal shape checks remain active even when JVM assertions are disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final String PAYLOAD="\nACGT\n+\nIIII\n";
	/** Exact same input order as the unchanged ReformatNameCollisionTest writer. */
	private static final String[] INPUT={"x", "x", "x_2", "x", "x_02", "x_2", "x_2_2", "x_2",
		"n a", "n_2 b", "n c", "p/1", "p/1", "p_2/1", "q", "q_2"};
	/** Literal first output block; the original ID is preserved until occupied. */
	private static final String[] FIRST={"x", "x_2", "x_2_2", "x_3", "x_02", "x_2_3", "x_2_2_2", "x_2_4",
		"n a", "n_2 b", "n_3 c", "p/1", "p_2/1", "p_2_2/1", "q", "q_2"};
}
