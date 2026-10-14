package fileIO;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

/**
 * Finite Java ZIP ownership and retained-error checks. Run one case per JVM with -ea.
 * Arguments are mode=baseline|fixed, api=factory|writer, name=stdout|stderr|file and
 * sink=normal|write|flush. File controls require normal. Global standard streams are
 * temporarily replaced and restored; concurrent users of those streams are unsupported.
 * Creates output.zip for file controls and captured.zip for healthy borrowed output
 * in the current directory. ZIP timestamps are not compared; one entry, its name and
 * literal ABCD payload are checked through ZipFile, including the central directory.
 * mode selects the expected pre/post-repair observations, not the implementation.
 * No other codec, performance, asynchronous or exceptional-writer-cleanup claim.
 * @author Shinobu
 */
public final class ZipOutputOwnershipTest{

	/*--------------------------------------------------------------*/
	/*----------------            Startup           ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs one isolated case and prints api, name, sink, threw and backend closes. */
	public static void main(String[] args) throws Exception{
		boolean assertions=false;
		assert(assertions=true);
		if(!assertions){throw new IllegalStateException("Run ZIP ownership checks with -ea.");}
		String mode=null, api=null, name=null, sink=null;
		for(String arg : args){
			final int eq=arg.indexOf('=');
			if(eq<1){throw new IllegalArgumentException("Expected key=value: "+arg);}
			final String key=arg.substring(0, eq), value=arg.substring(eq+1);
			if(key.equals("mode")){
				mode=value;
			}else if(key.equals("api")){
				api=value;
			}else if(key.equals("name")){
				name=value;
			}else if(key.equals("sink")){
				sink=value;
			}else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(!"baseline".equals(mode) && !"fixed".equals(mode)){
			throw new IllegalArgumentException("mode=baseline|fixed required");
		}
		if(!"factory".equals(api) && !"writer".equals(api)){
			throw new IllegalArgumentException("api=factory|writer required");
		}
		if(!"stdout".equals(name) && !"stderr".equals(name) && !"file".equals(name)){
			throw new IllegalArgumentException("name=stdout|stderr|file required");
		}
		if(!"normal".equals(sink) && !"write".equals(sink) && !"flush".equals(sink)){
			throw new IllegalArgumentException("sink=normal|write|flush required");
		}
		final boolean file=name.equals("file"), fixed=mode.equals("fixed"), normal=sink.equals("normal");
		if(file && !normal){throw new IllegalArgumentException("File controls require normal");}
		final String filename=file ? "output.zip" : name+".zip";
		final ControlledOutput backend=new ControlledOutput(sink);
		final PrintStream controlled=new PrintStream(backend, false, "US-ASCII");
		final PrintStream savedOut=System.out, savedErr=System.err;
		Exception failure=null;
		try{
			if(!file){
				if(name.equals("stdout")){System.setOut(controlled);}else{System.setErr(controlled);}
			}
			try{write(api, filename);}
			catch(IOException | RuntimeException e){failure=e;}
			assert((failure!=null)==(!normal && fixed)) :
				"ZIP retained errors must become visible after repair; api="+api+" sink="+sink+" failure="+failure;
			if(failure!=null){
				assert(recognized(failure, filename)) : "Unexpected exception instead of borrowed-output failure: "+failure;
			}
			if(!file){
				assert(backend.closes==(fixed ? 0 : 1)) :
					"ZIP must preserve the borrowed backend; mode="+mode+" closes="+backend.closes;
			}
			if(normal){
				final Path archive=Paths.get(file ? filename : "captured.zip");
				final byte[] bytes=file ? Files.readAllBytes(archive) : backend.bytes.toByteArray();
				if(!file){Files.write(archive, bytes);}
				verifyZip(archive, file ? "output" : name);
				if(!file){
					controlled.write(AFTER);
					controlled.flush();
					final byte[] expected=Arrays.copyOf(bytes, bytes.length+(fixed ? AFTER.length : 0));
					if(fixed){System.arraycopy(AFTER, 0, expected, bytes.length, AFTER.length);}
					assert(Arrays.equals(expected, backend.bytes.toByteArray())) :
						"Subsequent caller marker must reach the still-open borrowed stream";
				}
			}
		}finally{
			System.setOut(savedOut);
			System.setErr(savedErr);
		}
		System.out.println(api+"\t"+name+"\t"+sink+"\t"+(failure!=null)+"\t"+(file ? "NA" : backend.closes));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Helpers        ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes the same literal payload through the factory or synchronous character helper. */
	private static void write(final String api, final String filename) throws IOException{
		if(api.equals("writer")){ReadWrite.writeString("ABCD", filename, true, false); return;}
		try(OutputStream out=ReadWrite.getZipOutputStream(filename, true, false)){
			assert(out.getClass()==ZipOutputStream.class) : "Existing callers require exact ZipOutputStream.class";
			out.write(PAYLOAD);
			out.flush();
		}
	}

	/** Accepts only the intended retained-error message or its explicit IOException cause. */
	private static boolean recognized(final Exception failure, final String filename){
		if(failure instanceof RuntimeException && ("Error writing to "+filename).equals(failure.getMessage())){return true;}
		final Throwable cause=failure instanceof IOException ? failure : failure.getCause();
		return cause instanceof IOException && "Error writing to borrowed standard output".equals(cause.getMessage());
	}

	/** Requires one complete archive entry and the independent literal payload, without timestamp assumptions. */
	private static void verifyZip(final Path path, final String name) throws IOException{
		try(ZipFile zip=new ZipFile(path.toFile())){
			final Enumeration<? extends ZipEntry> entries=zip.entries();
			assert(entries.hasMoreElements()) : "Expected one ZIP entry in "+path;
			final ZipEntry entry=entries.nextElement();
			assert(!entries.hasMoreElements() && entry.getName().equals(name)) : "Expected sole ZIP entry "+name;
			assert(entry.getMethod()==ZipEntry.DEFLATED) : "Java ZIP control must remain deflated";
			try(InputStream in=zip.getInputStream(entry)){
				for(byte expected : PAYLOAD){assert(in.read()==expected) : "ZIP payload must be literal ABCD";}
				assert(in.read()==-1) : "ZIP payload must end after ABCD";
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Controlled Sink       ----------------*/
	/*--------------------------------------------------------------*/

	/** Counts raw closure and independently injects write or flush IOExceptions under PrintStream. */
	private static final class ControlledOutput extends OutputStream{
		/** Selects normal operation, write failure or flush failure. */
		ControlledOutput(final String sink_){sink=sink_;}
		/** Retains one byte or throws before retaining it. */
		@Override
		public void write(final int value) throws IOException{
			if(sink.equals("write")){throw new IOException("Controlled ZIP write failure");}
			bytes.write(value);
		}
		/** Retains one complete block or throws before retaining any of it. */
		@Override
		public void write(final byte[] data, final int offset, final int length) throws IOException{
			if(sink.equals("write")){throw new IOException("Controlled ZIP write failure");}
			bytes.write(data, offset, length);
		}
		/** Injects only the selected flush failure. */
		@Override
		public void flush() throws IOException{
			if(sink.equals("flush")){throw new IOException("Controlled ZIP flush failure");}
		}
		/** Counts closure without discarding evidence. */
		@Override
		public void close(){closes++;}
		/** Selected failure site; independent of production behavior. */
		final String sink;
		/** Bytes retained before any subsequent caller marker. */
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		/** Number of raw close delegations, including on failure. */
		int closes;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Independent payload bytes, unaffected by encoder or compression implementation. */
	private static final byte[] PAYLOAD={65, 66, 67, 68};
	/** Bytes written by a later owner to prove continued usability. */
	private static final byte[] AFTER={124, 97, 102, 116, 101, 114};
}
