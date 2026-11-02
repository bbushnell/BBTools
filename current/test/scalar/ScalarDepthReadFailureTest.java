package test.scalar;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipException;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ReadWrite;
import map.ObjectDoubleMap;
import scalar.ScalarData;
import shared.Shared;

/** Finite full-read error propagation and literal header-depth checks for ScalarData.
 * Run each case in a fresh JVM: global reader/decompression settings are configured
 * for synchronous ByteFile1 and Java gzip. Inputs exceed one reader buffer but
 * compress to small fixtures. Only the CRC trailer differs in each corrupt case.
 * No external decompressor, sequence-processing worker or production file is used.
 * @author Sucrose
 */
public final class ScalarDepthReadFailureTest{

	/** Requires work=fresh-directory and case=pileup-valid|pileup-crc|sam-valid|sam-crc,
	 * or header-tadpole|header-spades|header-generic|header-missing. */
	public static void main(String[] args) throws Exception{
		String mode=null;
		Path work=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("work=")){work=Paths.get(arg.substring(5));}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		final boolean header="header-tadpole".equals(mode) || "header-spades".equals(mode)
				|| "header-generic".equals(mode) || "header-missing".equals(mode);
		require(header || "pileup-valid".equals(mode) || "pileup-crc".equals(mode)
				|| "sam-valid".equals(mode) || "sam-crc".equals(mode), "Choose a supported depth loader/header case");
		require(work!=null, "work=fresh-directory is required");
		if(header){
			Files.createDirectory(work);
			checkHeader(mode);
			return;
		}
		Shared.setThreads(1);
		ByteFile.FORCE_MODE_BF1=true;
		ByteFile.FORCE_MODE_BF2=false;
		ReadWrite.RAWMODE=false;
		ReadWrite.USE_UNPIGZ=ReadWrite.USE_GUNZIP=ReadWrite.USE_UNBGZIP=false;
		ReadWrite.ALLOW_NATIVE_BGZF=false;
		final boolean pileup=mode.startsWith("pileup"), corrupt=mode.endsWith("crc");
		final byte[] raw=payload(pileup), compressed=gzip(raw);
		require(raw.length>ByteFile1.bufferlen, "Fixture must exceed the configured ByteFile1 buffer: "+ByteFile1.bufferlen);
		if(corrupt){compressed[compressed.length-8]^=1;}//CRC32's first trailer byte; payload and ISIZE are unchanged.
		verifyGzip(compressed, raw, corrupt);
		Files.createDirectory(work);
		final Path input=work.resolve(pileup ? "depth.cov.gz" : "depth.sam.gz");
		Files.write(input, compressed);
		final Constructor<ScalarData> constructor=ScalarData.class.getDeclaredConstructor(boolean.class, long.class);
		constructor.setAccessible(true);
		final ScalarData data=constructor.newInstance(false, -1L);
		final Method loader=ScalarData.class.getDeclaredMethod(pileup ? "loadPileupCoverage" : "loadDepthFromAlignment", String.class);
		loader.setAccessible(true);
		Object result=null;
		RuntimeException failure=null;
		try{result=invoke(loader, data, input.toString());}
		catch(RuntimeException expected){failure=expected;}
		final int open=openInputs(input);
		System.out.println("case="+mode+" rawBytes="+raw.length+" gzipBytes="+compressed.length
				+" bufferBytes="+ByteFile1.bufferlen+" failure="+failure+" openInputs="+open);
		require(open==0, "The full-read loader must close its owned input before returning/throwing");
		if(corrupt){
			require(failure!=null, "CRC-corrupt full input returned normally instead of propagating its latched read failure");
		}else{
			require(failure==null, "Valid full input failed: "+failure);
			require(result instanceof ObjectDoubleMap<?>, "Expected a depth map");
			@SuppressWarnings("unchecked")
			final ObjectDoubleMap<String> depths=(ObjectDoubleMap<String>)result;
			final double expected=pileup ? 7 : 100;
			require(depths.size()==1 && depths.get("ctg")==expected,
					"Expected only ctg at depth "+expected+", actual="+depths.get("ctg")+" entries="+depths.size());
			System.out.println("ctgDepth="+depths.get("ctg")+" expected="+expected);
		}
		System.out.println("ScalarDepthReadFailureTest "+mode+" PASS");
	}

	/** Uses literal values independent of the production field parser.
	 * Contig.toHeader emits contig_ID,len=...,cov=...,gc=... before optional suffixes;
	 * the ordered =,=, separators put coverage in zero-based field3. */
	private static void checkHeader(String mode) throws Exception{
		final String name;
		final float expected;
		if(mode.equals("header-tadpole")){
			name="contig_1,len=100,cov=2.5,gc=0.400";
			expected=2.5f;
		}else if(mode.equals("header-spades")){
			name="NODE_1_length_100_cov_2.5";
			expected=2.5f;
		}else if(mode.equals("header-generic")){
			name="sample_cov_3.5_tail";
			expected=3.5f;
		}else{
			name="contig_1,len=100,gc=0.400";
			expected=-1;
		}
		final Constructor<ScalarData> constructor=ScalarData.class.getDeclaredConstructor(boolean.class, long.class);
		constructor.setAccessible(true);
		final ScalarData data=constructor.newInstance(false, -1L);
		final float actual=data.parseDepthFromHeader(name);
		require(actual==expected, "Header depth must come from the literal coverage value: "+name+" expected="+expected+" actual="+actual);
		System.out.println("ScalarDepthReadFailureTest "+mode+" PASS expected="+expected+" actual="+actual);
	}

	/** Pileup repeatedly assigns literal depth7; SAM has10000 ten-base mappings to a1000-base reference. */
	private static byte[] payload(boolean pileup){
		final StringBuilder text=new StringBuilder();
		if(pileup){
			text.append("#ID\tAvg_fold\n");
			for(int i=0; i<20000; i++){text.append("ctg\t7\n");}
		}else{
			text.append("@HD\tVN:1.6\n@SQ\tSN:ctg\tLN:1000\n");
			for(int i=0; i<10000; i++){
				text.append("read\t0\tctg\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n");
			}
		}
		return text.toString().getBytes(StandardCharsets.US_ASCII);
	}

	/** Produces the same complete gzip member before the optional one-byte trailer mutation. */
	private static byte[] gzip(byte[] raw) throws Exception{
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		try(GZIPOutputStream stream=new GZIPOutputStream(bytes)){stream.write(raw);}
		return bytes.toByteArray();
	}

	/** Independently validates exact good payload bytes and proves the corrupt trailer is rejected by Java gzip. */
	private static void verifyGzip(byte[] compressed, byte[] raw, boolean corrupt) throws Exception{
		final ByteArrayOutputStream decoded=new ByteArrayOutputStream();
		ZipException failure=null;
		try(GZIPInputStream stream=new GZIPInputStream(new ByteArrayInputStream(compressed))){
			final byte[] buffer=new byte[8192];
			for(int n; (n=stream.read(buffer))>=0;){decoded.write(buffer, 0, n);}
		}catch(ZipException expected){failure=expected;}
		if(corrupt){require(failure!=null, "The independent Java gzip reader accepted the corrupted CRC");}
		else{require(failure==null && Arrays.equals(raw, decoded.toByteArray()), "Valid gzip payload differs from its literal input");}
	}

	/** Preserves the production loader's underlying exception rather than accepting a reflection wrapper. */
	private static Object invoke(Method method, Object target, String filename) throws Exception{
		try{return method.invoke(target, filename);}
		catch(InvocationTargetException wrapper){
			final Throwable cause=wrapper.getCause();
			if(cause instanceof Exception){throw (Exception)cause;}
			if(cause instanceof Error){throw (Error)cause;}
			throw new AssertionError("Unexpected reflected failure", cause);
		}
	}

	/** Counts descriptors still naming this fixture on Linux; no arbitrary process state is inspected. */
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

	/** Keeps independent oracles active under both JVM assertion modes. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
