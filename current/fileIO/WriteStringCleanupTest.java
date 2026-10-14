package fileIO;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

/**
 * Finite synchronous writer cleanup checks. Run each case in a fresh JVM with -ea.
 * Required arguments: mode=baseline|fixed, name=file|filezip|stdout|stderr|stdoutzip|stderrzip,
 * fault=normal|body|status|close|same. Owned files accept normal or body only.
 * Body cases throw one sentinel from CharSequence.toString after the writer is opened.
 * status injects an IOException during cleanup; close throws a distinct unchecked
 * cleanup sentinel; same rethrows the original sentinel during cleanup.
 * Standard streams are temporarily replaced and restored. Borrowed error state is
 * explicitly cleared before the later-marker check, which tests ownership, not recovery.
 * Owned descriptor checks require Linux /proc/self/fd and an Epsilon-GC JVM so garbage
 * collection cannot hide an unclosed descriptor. The saved runner supplies that JVM.
 * Creates output.txt/output.zip/captured.zip in cwd; no concurrent use of global streams.
 * No asynchronous, arbitrary backend, permanent-leak, performance or other-codec claim.
 * @author Shinobu
 */
public final class WriteStringCleanupTest{

	/*--------------------------------------------------------------*/
	/*----------------            Startup           ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs one case; prints name, fault, primary failure, suppressed count and target descriptors. */
	public static void main(String[] args) throws Exception{
		boolean assertions=false;
		assert(assertions=true);
		if(!assertions){throw new IllegalStateException("Run cleanup checks with -ea.");}
		String mode=null, name=null, fault=null;
		for(String arg : args){
			final int eq=arg.indexOf('=');
			if(eq<1){throw new IllegalArgumentException("Expected key=value: "+arg);}
			final String key=arg.substring(0, eq), value=arg.substring(eq+1);
			if(key.equals("mode")){
				mode=value;
			}else if(key.equals("name")){
				name=value;
			}else if(key.equals("fault")){
				fault=value;
			}else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(!"baseline".equals(mode) && !"fixed".equals(mode)){
			throw new IllegalArgumentException("mode=baseline|fixed required");
		}
		if(!Arrays.asList("file", "filezip", "stdout", "stderr", "stdoutzip", "stderrzip").contains(name)){
			throw new IllegalArgumentException("Unsupported name: "+name);
		}
		if(!Arrays.asList("normal", "body", "status", "close", "same").contains(fault)){
			throw new IllegalArgumentException("Unsupported fault: "+fault);
		}
		final boolean fixed=mode.equals("fixed"), owned=name.startsWith("file");
		final boolean zipped=name.endsWith("zip"), normal=fault.equals("normal");
		if(owned && !normal && !fault.equals("body")){
			throw new IllegalArgumentException("Owned files accept normal/body only");
		}
		final String stem=owned ? "output" : name.startsWith("stdout") ? "stdout" : "stderr";
		final String filename=stem+(zipped ? ".zip" : owned ? ".txt" : "");
		final Path path=Paths.get(filename).toAbsolutePath().normalize();
		final RuntimeException primary=new RuntimeException("Controlled body failure");
		final RuntimeException closing=new RuntimeException("Controlled unchecked cleanup failure");
		final Sink sink=new Sink();
		final FaultPrintStream stream=new FaultPrintStream(sink, primary, closing, fault);
		final PrintStream savedOut=System.out, savedErr=System.err;
		RuntimeException failure=null;
		int descriptors=-1;
		try{
			if(!owned){
				if(stem.equals("stdout")){System.setOut(stream);}else{System.setErr(stream);}
			}else{assert(countOpen(path)==0) : "Fresh owned output must have no preexisting descriptor";}
			final CharSequence text=normal ? "ABCD" : new FailingText(primary, stream);
			try{ReadWrite.writeString(text, filename, true, false);}
			catch(RuntimeException e){failure=e;}
			assert(normal ? failure==null : failure==primary) :
				"Cleanup must preserve the exact original exception; fault="+fault+" observed="+failure;
			final int expectedSuppressed=fixed && (fault.equals("status") || fault.equals("close")) ? 1 : 0;
			assert(primary.getSuppressed().length==expectedSuppressed) :
				"Only a distinct cleanup failure may be suppressed; mode="+mode+" fault="+fault;
			if(expectedSuppressed>0){
				final Throwable actual=primary.getSuppressed()[0];
				assert(fault.equals("close") ? isCloseFailure(actual, closing) :
					actual instanceof RuntimeException && ("Error writing to "+filename).equals(actual.getMessage())) :
					"Suppressed failure must be the known close sentinel/JDK wrapper or retained-error report: "+actual;
				assert(actual.getSuppressed().length==0) : "Cleanup must not suppress the primary back beneath itself";
			}
			if(owned){
				descriptors=countOpen(path);
				assert(descriptors==(!normal && !fixed ? 1 : 0)) :
					"Owned writer must close its descriptor on body failure; mode="+mode+" descriptors="+descriptors;
				if(normal){verify(path, stem, zipped, true);}else if(fixed){verify(path, stem, zipped, false);}
			}else{
				assert(sink.closes==0) : "Temporary writer cleanup must never close caller-owned standard streams";
				final byte[] before=sink.bytes.toByteArray();
				if(normal || (fixed && fault.equals("body"))){
					final Path captured=Paths.get(zipped ? "captured.zip" : "captured.txt");
					Files.write(captured, before);
					verify(captured, stem, zipped, normal);
				}
				stream.disarm();
				stream.write(AFTER);
				stream.flush();
				final byte[] expected=Arrays.copyOf(before, before.length+AFTER.length);
				System.arraycopy(AFTER, 0, expected, before.length, AFTER.length);
				assert(Arrays.equals(expected, sink.bytes.toByteArray())) : "Borrowed stream must accept a later owner's marker";
			}
		}finally{
			System.setOut(savedOut);
			System.setErr(savedErr);
		}
		System.out.println(name+"\t"+fault+"\t"+(failure==primary)+"\t"+primary.getSuppressed().length+"\t"+descriptors);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Helpers        ----------------*/
	/*--------------------------------------------------------------*/

	/** JDK25 StreamEncoder can wrap a repeated flush sentinel when try-with-resources suppresses itself. */
	private static boolean isCloseFailure(final Throwable actual, final Throwable closing){
		return actual==closing || (actual instanceof IllegalArgumentException && actual.getCause()==closing
			&& "Self-suppression not permitted".equals(actual.getMessage()));
	}

	/** Counts only descriptors naming the exact output path, closing the directory iterator itself. */
	private static int countOpen(final Path target) throws IOException{
		int count=0;
		try(DirectoryStream<Path> entries=Files.newDirectoryStream(Paths.get("/proc/self/fd"))){
			for(Path entry : entries){
				try{if(Files.readSymbolicLink(entry).equals(target)){count++;}}
				catch(NoSuchFileException e){/* A transient JVM descriptor may have closed during enumeration. */}
			}
		}
		return count;
	}

	/** Checks literal output or an empty, finalized container after a pre-print failure. */
	private static void verify(final Path path, final String entryName, final boolean zipped,
			final boolean normal) throws IOException{
		final byte[] expected=normal ? PAYLOAD : new byte[0];
		if(!zipped){
			assert(Arrays.equals(Files.readAllBytes(path), expected)) : "Raw output must match the literal body expectation";
			return;
		}
		try(ZipFile zip=new ZipFile(path.toFile())){
			final Enumeration<? extends ZipEntry> entries=zip.entries();
			assert(entries.hasMoreElements()) : "Closing the writer must finalize its ZIP central directory";
			final ZipEntry entry=entries.nextElement();
			assert(entry.getName().equals(entryName) && !entries.hasMoreElements()) : "Expected exactly the named ZIP entry";
			try(InputStream in=zip.getInputStream(entry)){
				for(byte value : expected){assert(in.read()==value) : "ZIP payload must match literal expected bytes";}
				assert(in.read()==-1) : "ZIP entry must end after the expected payload";
			}
		}
	}

	/** Arms cleanup faults only after acquisition, then throws the exact body sentinel. */
	private static final class FailingText implements CharSequence{
		/** Retains the sentinels without arming them during factory construction. */
		FailingText(final RuntimeException failure_, final FaultPrintStream stream_){failure=failure_; stream=stream_;}
		@Override
		public String toString(){stream.armed=true; throw failure;}
		@Override
		public int length(){return 4;}
		@Override
		public char charAt(final int index){return "ABCD".charAt(index);}
		@Override
		public CharSequence subSequence(final int start, final int end){return "ABCD".subSequence(start, end);}
		/** The sentinel whose identity must survive cleanup. */
		final RuntimeException failure;
		/** Test-owned standard stream; arms its fault after output acquisition. */
		final FaultPrintStream stream;
	}

	/** Injects unchecked or retained status failures only after the body sentinel is armed. */
	private static final class FaultPrintStream extends PrintStream{
		/** Builds a test-owned ASCII stream with independent body and cleanup failures. */
		FaultPrintStream(final Sink sink_, final RuntimeException primary_, final RuntimeException closing_,
				final String fault_) throws IOException{
			super(sink_, false, "US-ASCII");
			sink=sink_; primary=primary_; closing=closing_; fault=fault_;
		}
		@Override
		public void flush(){
			if(armed && fault.equals("close")){throw closing;}
			if(armed && fault.equals("same")){throw primary;}
			sink.fail=armed && fault.equals("status");
			super.flush();
		}
		/** Clears only test faults/status before testing ownership with a subsequent write. */
		void disarm(){armed=false; sink.fail=false; clearError();}
		/** Retained byte sink, also used to inject checked flush failures. */
		final Sink sink;
		/** Original body and distinct cleanup identities. */
		final RuntimeException primary, closing;
		/** Selected fault site and whether the body has reached it. */
		final String fault;
		boolean armed;
	}

	/** Retains bytes, counts raw closure, and optionally rejects flush with an IOException. */
	private static final class Sink extends OutputStream{
		@Override
		public void write(final int value){bytes.write(value);}
		@Override
		public void write(final byte[] data, final int off, final int len){bytes.write(data, off, len);}
		@Override
		public void flush() throws IOException{if(fail){throw new IOException("Controlled cleanup status failure");}}
		@Override
		public void close(){closes++;}
		/** Actual bytes, retained even when closure is counted. */
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		/** Whether flush throws the controlled IOException. */
		boolean fail;
		/** Number of delegated raw close calls. */
		int closes;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Independent literal successful body and subsequent-owner marker. */
	private static final byte[] PAYLOAD="ABCD".getBytes(StandardCharsets.US_ASCII),
		AFTER="|after".getBytes(StandardCharsets.US_ASCII);
}
