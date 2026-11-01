package test.shared;

import java.io.ByteArrayOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.zip.GZIPOutputStream;

import fileIO.ByteFile;
import fileIO.ReadWrite;
import shared.Tools;
import stream.Read;
import structures.IntList;

/** Finite utility regressions with literal independent expectations.
 * Generated files remain in the supplied work directory for inspection.
 * @author Shinobu
 */
public final class ToolsBoundaryTest{

	/** Runs every case with explicit checks under both assertion modes.
	 * @param args Exactly one work=directory argument */
	public static void main(String[] args) throws Exception{
		if(args.length!=1 || !args[0].startsWith("work=")){
			throw new IllegalArgumentException("Expected work=directory");
		}
		final Path work=Paths.get(args[0].substring(5));
		Files.createDirectories(work);
		final ToolsBoundaryTest test=new ToolsBoundaryTest();
		test.check("low throughput", () -> {
			require(Tools.thingsProcessed(1000000000L, 1, 1, "items").endsWith("1.00 items/sec"), "one item per second");
			require(Tools.thingsProcessed(1000000000L, 999, 1, "items").endsWith("999.00 items/sec"), "999 items per second");
		});
		test.check("scaled throughput controls", () -> {
			require(Tools.thingsProcessed(1000000000L, 1000, 1, "items").endsWith("1.00k items/sec"), "1000 items per second");
			require(Tools.thingsProcessed(1000000000L, 1000000, 1, "items").endsWith("1.00m items/sec"), "one million items per second");
		});
		test.check("minimum-size final fragment", () -> fragments(4, new String[] {"ACGTAC", "GTAA"}));
		test.check("fragment controls", () -> {
			fragments(3, new String[] {"ACGTAC", "GTAA"});
			fragments(5, new String[] {"ACGTAC"});
		});
		test.check("all byte uppercase values", () -> {
			int failures=0;
			for(int value=0; value<256; value++){
				final byte input=(byte)value;
				final byte expected=(byte)(value>='a' && value<='z' ? value-32 : value);
				try{if(Tools.toUpperCase(input)!=expected){failures++;}}
				catch(RuntimeException e){failures++;}
			}
			require(failures==0, "uppercase mismatches or exceptions="+failures);
		});
		ByteFile.FORCE_MODE_BF1=true;
		ByteFile.FORCE_MODE_BF2=ByteFile.FORCE_MODE_BF3=ByteFile.FORCE_MODE_BF4=false;
		ReadWrite.USE_UNPIGZ=false;
		ReadWrite.ALLOW_NATIVE_BGZF=false;
		final byte[] compressed=gzip("1\n2\n");
		final Path valid=work.resolve("valid-ids.txt.gz"), corrupt=work.resolve("corrupt-ids.txt.gz");
		Files.write(valid, compressed);
		compressed[compressed.length-8]^=1;// Deliberately corrupt the CRC, keeping the payload intact.
		Files.write(corrupt, compressed);
		test.check("integer file control", () -> {
			final IntList list=Tools.loadIntegers(valid.toString());
			require(list!=null && list.size==2 && list.get(0)==1 && list.get(1)==2, "expected IDs 1,2");
		});
		test.check("corrupt integer file rejection", () -> {
			boolean rejected=false;
			try{Tools.loadIntegers(corrupt.toString());}
			catch(RuntimeException e){rejected=true;}
			require(rejected, "corrupt ID file must not be accepted as a complete list");
		});
		System.out.println("ToolsBoundaryTest cases="+test.cases+" failures="+test.failures);
		if(test.failures!=0){throw new AssertionError("Utility contract failures="+test.failures);}
	}

	/** Checks fragment bases, qualities and identifiers from a literal ten-base read. */
	private static void fragments(int minimum, String[] expected){
		final byte[] bases="ACGTACGTAA".getBytes(StandardCharsets.US_ASCII);
		final byte[] quality={20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
		final ArrayList<Read> reads=new ArrayList<Read>();
		reads.add(new Read(bases, quality, "read", 0, false));
		Tools.breakReads(reads, 6, minimum, null);
		require(reads.size()==expected.length, "fragment count="+reads.size()+", expected="+expected.length);
		int offset=0;
		for(int i=0; i<expected.length; i++){
			final Read read=reads.get(i);
			require(Arrays.equals(read.bases, expected[i].getBytes(StandardCharsets.US_ASCII)), "fragment bases "+i);
			require(Arrays.equals(read.quality, Arrays.copyOfRange(quality, offset, offset+expected[i].length())), "fragment qualities "+i);
			require(read.id.equals("read_"+(i+1)), "fragment ID "+i);
			offset+=expected[i].length();
		}
	}

	/** Creates a tiny gzip fixture without external compressors. */
	private static byte[] gzip(String text) throws Exception{
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		try(GZIPOutputStream gzip=new GZIPOutputStream(bytes)){
			gzip.write(text.getBytes(StandardCharsets.US_ASCII));
		}
		return bytes.toByteArray();
	}

	/** Runs an independent case and reports all failures before exiting. */
	private void check(String name, Action action){
		cases++;
		try{action.run();}
		catch(Exception | AssertionError e){
			failures++;
			System.out.println("FAIL "+name+": "+e);
		}
	}

	/** Checks expectations even when JVM assertions are disabled. */
	private static void require(boolean condition, String reason){
		if(!condition){throw new AssertionError(reason);}
	}

	/** A finite case that may create or read a fixture file. */
	private interface Action{void run() throws Exception;}
	/** Number of attempted cases. */
	private int cases;
	/** Number of failed cases. */
	private int failures;
}
