package test.bin;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipException;

import bin.DataLoader;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ReadWrite;
import map.ObjectDoubleMap;
import shared.Shared;
import shared.Tools;

/** Checked-close coverage for bin.DataLoader.loadCovFile2's owned ByteFile input.
 * Run each case in a fresh JVM because global ByteFile/gzip state is pinned to
 * synchronous Java gzip:
 *
 * case=valid validates literal two-sample depth sums.
 * case=crc expects full-read CRC corruption to propagate instead of returning
 * a partial map.
 * case=parse expects the owned input to close when line parsing throws.
 *
 * @author Jean
 */
public final class DataLoaderCovFile2Test{

	/** Entrypoint for one isolated valid, CRC, or parse-failure case. */
	public static void main(String[] args) throws Exception{
		String mode=null;
		Path work=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("work=")){work=Paths.get(arg.substring(5));}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		require("valid".equals(mode) || "crc".equals(mode) || "parse".equals(mode),
				"Choose case=valid, case=crc, or case=parse");
		require(work!=null, "work=fresh-directory is required");

		Shared.setThreads(1);
		ByteFile.FORCE_MODE_BF1=true;
		ByteFile.FORCE_MODE_BF2=false;
		ReadWrite.RAWMODE=false;
		ReadWrite.USE_UNPIGZ=ReadWrite.USE_GUNZIP=ReadWrite.USE_UNBGZIP=false;
		ReadWrite.ALLOW_NATIVE_BGZF=false;
		Tools.FORCE_JAVA_PARSE_DOUBLE=true;

		Files.createDirectory(work);
		if("valid".equals(mode)){checkValid(work);}
		else if("crc".equals(mode)){checkCrc(work);}
		else{checkParseFailure(work);}
		System.out.println("DataLoaderCovFile2Test "+mode+" PASS");
	}

	/** Verifies that loadCovFile2 sums all depth columns for each contig name. */
	private static void checkValid(Path work) throws Exception{
		final Path input=work.resolve("depth.cov");
		Files.write(input, validPayload().getBytes(StandardCharsets.US_ASCII));
		final ObjectDoubleMap<String> depths=DataLoader.loadCovFile2(input.toString());
		require(depths.size()==2, "Expected two contig depths, found "+depths.size());
		require(depths.get("alpha")==4.0, "alpha depth should be 1.25 + 2.75, actual="+depths.get("alpha"));
		require(depths.get("beta")==4.125, "beta depth should be 0.5 + 3.625, actual="+depths.get("beta"));
	}

	/** Requires a CRC failure latched at end-of-stream to be reported to the caller. */
	private static void checkCrc(Path work) throws Exception{
		final byte[] raw=largePayload();
		final byte[] validCompressed=gzip(raw);
		require(raw.length>ByteFile1.bufferlen,
				"Fixture must exceed the configured ByteFile1 buffer: "+ByteFile1.bufferlen);
		verifyGzip(validCompressed, raw, false);
		final Path validInput=work.resolve("valid-depth.cov.gz");
		Files.write(validInput, validCompressed);
		checkLargeDepths(DataLoader.loadCovFile2(validInput.toString()));

		final byte[] compressed=Arrays.copyOf(validCompressed, validCompressed.length);
		compressed[compressed.length-8]^=1;//CRC32's first trailer byte; payload and ISIZE stay unchanged.
		verifyGzip(compressed, raw, true);
		final Path input=work.resolve("depth.cov.gz");
		Files.write(input, compressed);

		RuntimeException failure=null;
		try{DataLoader.loadCovFile2(input.toString());}
		catch(RuntimeException expected){failure=expected;}
		require(failure!=null, "CRC-corrupt coverage returned normally instead of propagating the close failure");
	}

	/** Requires the reader to close even if a row throws before EOF. */
	private static void checkParseFailure(Path work) throws Exception{
		final Path input=work.resolve("bad-depth.cov");
		Files.write(input, badPayload().getBytes(StandardCharsets.US_ASCII));

		NumberFormatException failure=null;
		try{DataLoader.loadCovFile2(input.toString());}
		catch(NumberFormatException expected){failure=expected;}
		final int open=openInputs(input);
		System.out.println("case=parse forceJavaParseDouble="+Tools.FORCE_JAVA_PARSE_DOUBLE+
				" failure="+failure+" openInputs="+open);
		require(failure!=null, "Malformed coverage returned normally instead of failing with NumberFormatException");
		require(isParseFailure(failure),
				"Malformed coverage failed outside the expected input parser path: "+failure);
		require(open==0, "loadCovFile2 must close its owned input when parsing fails");
	}

	/** Small literal coverage file with two samples and two contigs. */
	private static String validPayload(){
		return "#Contigs\t2\n"
				+"#Depths\t2\n"
				+"#ShortName\tID\tSize\tCov_0\tCov_1\n"
				+"alpha\t0\t100\t1.25\t2.75\n"
				+"beta\t1\t80\t0.5\t3.625\n";
	}

	/** Verifies literal values from the large CRC fixture's unmutated gzip path. */
	private static void checkLargeDepths(ObjectDoubleMap<String> depths){
		require(depths.size()==20000, "Expected 20000 contig depths, found "+depths.size());
		require(depths.get("contig_0")==3.0, "contig_0 depth should be 1.0 + 2.0, actual="+depths.get("contig_0"));
		require(depths.get("contig_19999")==3.0,
				"contig_19999 depth should be 1.0 + 2.0, actual="+depths.get("contig_19999"));
	}

	/** The same logical coverage layout repeated past ByteFile1's buffer length. */
	private static byte[] largePayload(){
		final StringBuilder text=new StringBuilder();
		text.append("#Contigs\t20000\n");
		text.append("#Depths\t2\n");
		text.append("#ShortName\tID\tSize\tCov_0\tCov_1\n");
		for(int i=0; i<20000; i++){
			text.append("contig_").append(i).append('\t').append(i).append("\t100\t1.0\t2.0\n");
		}
		return text.toString().getBytes(StandardCharsets.US_ASCII);
	}

	/** A large coverage file with an early malformed Java-parser depth token. */
	private static String badPayload(){
		final StringBuilder text=new StringBuilder();
		text.append("#Contigs\t20000\n");
		text.append("#Depths\t2\n");
		text.append("#ShortName\tID\tSize\tCov_0\tCov_1\n");
		text.append("bad\t0\t10\tbad-depth\t1.0\n");
		for(int i=1; i<20000; i++){
			text.append("contig_").append(i).append('\t').append(i).append("\t100\t1.0\t2.0\n");
		}
		return text.toString();
	}

	/** Produces one complete gzip member before the optional trailer mutation. */
	private static byte[] gzip(byte[] raw) throws Exception{
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		try(GZIPOutputStream stream=new GZIPOutputStream(bytes)){stream.write(raw);}
		return bytes.toByteArray();
	}

	/** Proves the corrupt trailer is independently rejected by Java gzip. */
	private static void verifyGzip(byte[] compressed, byte[] raw, boolean corrupt) throws Exception{
		final ByteArrayOutputStream decoded=new ByteArrayOutputStream();
		ZipException failure=null;
		try(GZIPInputStream stream=new GZIPInputStream(new ByteArrayInputStream(compressed))){
			final byte[] buffer=new byte[8192];
			for(int n; (n=stream.read(buffer))>=0;){decoded.write(buffer, 0, n);}
		}catch(ZipException expected){failure=expected;}
		if(corrupt){require(failure!=null, "The independent Java gzip reader accepted the corrupted CRC");}
		else{require(failure==null && Arrays.equals(raw, decoded.toByteArray()),
				"Valid gzip payload differs from its literal input");}
	}

	/** Confines case=parse acceptance to LineParser1/Parse failures from loadCovFile2. */
	private static boolean isParseFailure(Throwable failure){
		if(failure==null){return false;}
		boolean sawDataLoader=false, sawParser=false;
		for(StackTraceElement frame : failure.getStackTrace()){
			if("bin.DataLoader".equals(frame.getClassName()) && "loadCovFile2".equals(frame.getMethodName())){
				sawDataLoader=true;
			}else if(frame.getClassName().startsWith("parse.")){
				sawParser=true;
			}
		}
		return sawDataLoader && sawParser;
	}

	/** Counts descriptors still naming this fixture on Linux. */
	private static int openInputs(Path input) throws Exception{
		final Path directory=Paths.get("/proc/self/fd"), target=input.toAbsolutePath().normalize();
		require(Files.isDirectory(directory), "Linux /proc/self/fd is required for the closure check");
		int count=0;
		try(DirectoryStream<Path> entries=Files.newDirectoryStream(directory)){
			for(Path entry : entries){
				try{if(Files.readSymbolicLink(entry).equals(target)){count++;}}
				catch(java.nio.file.NoSuchFileException closed){/* Descriptor closed while inspecting it. */}
			}
		}
		return count;
	}

	/** Throws AssertionError with a focused message when a fixture invariant fails. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

}
