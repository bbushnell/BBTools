package fileIO;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FilterInputStream;
import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.Reader;
import java.lang.ProcessBuilder.Redirect;
import java.net.MalformedURLException;
import java.net.URL;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;

import dna.Data;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import stream.ConcurrentReadOutputStream;
import stream.ConcurrentReadStreamInterface;
import stream.MultiCros;
import stream.Streamer;
import stream.Writer;
import stream.bam.BamInputStream;
import stream.bam.BamOutputStream;
import stream.bam.BgzfInputStream;
import stream.bam.BgzfInputStreamMT2;
import stream.bam.BgzfOutputStream;
import stream.bam.BgzfOutputStreamMT;
import stream.bam.BgzfOutputStreamMT2;
import stream.bam.BgzfSettings;
import structures.ByteBuilder;

/**
 * Stream factories, whole-file helpers, filename utilities and legacy writer accounting.
 * Compression selection combines case-sensitive factory suffix checks, mutable global
 * settings and external-tool availability; filename metadata helpers have different
 * alias/case rules. Factories may open files or start processes before returning.
 * Callers own returned streams and must complete the applicable close/process protocol.
 * Standard streams are borrowed, but only some close helpers protect their direct identities.
 *
 * The methods named writeStringAsync and writeAsync execute in the calling thread;
 * only the InThread wrappers start workers. Global settings are startup configuration,
 * not per-call snapshots. Synchronization protects selected counters/maps and write/copy
 * sections, not all I/O. Failure cleanup and propagation differ among legacy helpers;
 * the source TODOs identify unresolved cases rather than promising universal cleanup.
 * @author Brian Bushnell, Shinobu
 */
public class ReadWrite{
	
	
	/**
	 * Program entry point for file copying functionality.
	 * Asserts that the destination does not exist, then delegates to copyFile.
	 * @param args Command-line arguments [source, destination]
	 */
	public static void main(String[] args){
		File f=new File(args[1]);
		assert(!f.exists()) : "Destination file already exists.";
		copyFile(args[0], args[1]);
	}
	
	/**
	 * Starts a counted background character writer and returns without joining it.
	 * @param x Character sequence to write
	 * @param fname Output filename
	 */
	public static void writeStringInThread(CharSequence x, String fname){
		writeStringInThread(x, fname, false);
	}
	
	/**
	 * Starts a counted background character writer and returns without joining it.
	 * @param x Character sequence to write
	 * @param fname Output filename
	 * @param append Whether to append to existing file or overwrite
	 */
	public static void writeStringInThread(CharSequence x, String fname, boolean append){
		startWritingThread(new WriteStringThread(x, fname, append));
	}
	
	/**
	 * Starts a counted background object writer and returns without a completion handle.
	 * @param x Object to serialize and write
	 * @param fname Output filename
	 * @param allowSubprocess Whether to allow subprocess compression
	 */
	public static void writeObjectInThread(Object x, String fname, boolean allowSubprocess){
		startWritingThread(new WriteObjectThread(x, fname, allowSubprocess));
	}
	
	/** Reserves a counted worker only after construction; rolls back a failed start.
	 * @param writer Worker that transfers its reservation to running on entry */
	private static void startWritingThread(final Runnable writer){
		final Thread thread=new Thread(writer);
		addThread(1);
		boolean started=false;
		try{
			thread.start();
			started=true;
		}finally{
			if(!started){
				synchronized(activeThreads){
					//No worker ran, so this reservation is waiting, not running.
					activeThreads[0]--;
					activeThreads[1]--;
					assert(activeThreads[0]>=0 && activeThreads[1]>=0) :
						"Failed writer start must release exactly its waiting reservation: "+Arrays.toString(activeThreads);
					activeThreads.notifyAll();
				}
			}
		}
	}

	/** Counted worker retaining the supplied CharSequence reference without a defensive copy. */
	private static class WriteStringThread implements Runnable{
		
		private final CharSequence x;
		private final String fname;
		private final boolean append;
		WriteStringThread(CharSequence x_, String fname_, boolean append_){
			x=x_;
			fname=fname_;
			append=append_;
		}
		
		@Override
		public void run(){
			if(verbose){System.err.println("WriteStringThread.run() started for fname "+fname);}
			addRunningThread(1);
			try{writeStringAsync(x, fname, append);}
			finally{addThread(-1);}
			if(verbose){System.err.println("WriteStringThread.run() finished for fname "+fname);}
		}
		
	}
	
	/** Counted worker retaining the supplied object until serialization completes. */
	private static class WriteObjectThread implements Runnable{
		
		private final Object x;
		private final String fname;
		private final boolean allowSubprocess;
		WriteObjectThread(Object x_, String fname_, boolean allowSubprocess_){
			x=x_;
			fname=fname_;
			allowSubprocess=allowSubprocess_;
		}
		
		@Override
		public void run(){
			if(verbose){System.err.println("WriteObjectThread.run() started for fname "+fname);}
			addRunningThread(1);
			try{writeAsync(x, fname, allowSubprocess);}
			finally{addThread(-1);}
			if(verbose){System.err.println("WriteObjectThread.run() finished for fname "+fname);}
		}
		
	}
	
	/**
	 * Attempts to set read, write and execute permissions on an existing path.
	 * Setter boolean results are currently ignored; true does not establish success (FIO-045).
	 *
	 * @param fname File path
	 * @param read Whether to grant read permission
	 * @param write Whether to grant write permission
	 * @param execute Whether to grant execute permission
	 * @param ownerOnly Whether permissions apply to owner only
	 * @return False for a missing path or caught Exception; otherwise true even if a setter returned false
	 */
	public static boolean setPermissions(String fname, boolean read, boolean write, boolean execute, boolean ownerOnly){
		File f=new File(fname);
		if(!f.exists()){return false;}
		try{
			//TODO: Probable bug - FIO-045: ignored boolean setter results can report success after failure.
			f.setReadable(read, ownerOnly);
			f.setWritable(write, ownerOnly);
			f.setExecutable(execute, ownerOnly);
		}catch(Exception e){
			return false;
		}
		return true;
	}

	/** Writes character sequence to file, overwriting any existing content */
	public static void writeString(CharSequence x, String fname){writeString(x, fname, false);}
	/**
	 * Writes character sequence to file with append option.
	 * @param x Character sequence to write
	 * @param fname Output filename
	 * @param append Whether to append to existing file or overwrite
	 */
	public static void writeString(CharSequence x, String fname, boolean append){
		writeString(x, fname, !append, append);
	}
	/**
	 * Writes character sequence to file with explicit overwrite and append control.
	 * Uses PrintWriter charset rules under diskSync. Overwrite is assertion-checked.
	 * Cleanup is attempted on normal and exceptional paths; retained writer errors are reported.
	 * An original failure stays primary, with distinct cleanup failures suppressed beneath it.
	 * Direct System.out/System.err backends remain open while the temporary writer is finalized.
	 *
	 * @param x Character sequence to write
	 * @param fname Output filename
	 * @param overwrite Whether to allow overwriting existing files
	 * @param append Whether to append to existing file
	 * @throws RuntimeException If output creation or writing fails
	 */
	public static void writeString(CharSequence x, String fname, boolean overwrite, boolean append){
		if(verbose){System.err.println("writeString(x, "+fname+", "+append+")");}
		File f=new File(fname);
		assert(overwrite || append || !f.exists()) : "File "+fname+" exists and overwrite=f";
		OutputStream os=getOutputStream(fname, append, true, false);
		final PrintStream borrowed=(os==System.out || os==System.err ? (PrintStream)os : null);
		
		try{

			synchronized(diskSync){
				PrintWriter out=null;
				Throwable failure=null;
				try{
					out=(borrowed==null ? new PrintWriter(os) : borrowedWriter(borrowed));
					out.print(x);
					out.flush();

					if(os.getClass()==ZipOutputStream.class){
						ZipOutputStream zos=(ZipOutputStream)os;
						zos.closeEntry();
						zos.finish();
					}
//				else if(PROCESS_XZ && os.getClass()==org.tukaani.xz.XZOutputStream.class){
//					org.tukaani.xz.XZOutputStream zos=(org.tukaani.xz.XZOutputStream)os;
//					zos.finish();
//				}
				}catch(Throwable t){
					failure=t;
					throw t;
				}finally{
					try{closeStringWriter(out, os, borrowed, fname, failure);}
					catch(Throwable closing){
						if(failure==null){throw closing;}
						//JDK encoder cleanup may wrap a repeated primary in a self-suppression exception.
						if(closing!=failure && closing.getCause()!=failure){failure.addSuppressed(closing);}
					}
				}
			}
//			System.out.println("Wrote to "+fname);
			
//			String read=readString(fname);
//			assert(x.equals(read)) : x.length()+", "+read.length();
			
		}catch(FileNotFoundException e){
			throw new RuntimeException(e);
		}catch(IOException e){
			throw new RuntimeException(e);
		}catch(OutOfMemoryError e){
			KillSwitch.memKill(e);
		}
	}

	/** Attempts writer cleanup while retaining caller ownership and reporting stored errors.
	 * If construction or unchecked writer closure fails, owned/wrapped backend closure
	 * is attempted. Borrowed encoders still close normally to finalize their characters.
	 * @param out Constructed writer, or null if construction failed
	 * @param backend Output returned by the factory, possibly null on a faulty factory path
	 * @param borrowed Captured direct standard stream, or null for owned/wrapped backends
	 * @param fname Filename used in the retained-error diagnostic
	 * @param primary Original body/construction failure, excluded from nested suppression
	 * @throws IOException If fallback backend closure fails
	 */
	private static void closeStringWriter(final PrintWriter out, final OutputStream backend,
			final PrintStream borrowed, final String fname, final Throwable primary) throws IOException{
		assert(borrowed==null || backend==borrowed) : "writeString marks only direct standard backends as borrowed; wrappers require closure.";
		if(out==null){
			if(backend!=null && borrowed==null){backend.close();}
		}else{
			try{out.close();}
			catch(Throwable failure){
				try{if(backend!=null && borrowed==null){backend.close();}}
				catch(Throwable cleanup){
					//Stored exceptions can recur; do not build a suppression cycle through the primary.
					if(cleanup!=failure && cleanup!=primary && cleanup.getCause()!=failure
							&& (primary==null || cleanup.getCause()!=primary)){
						failure.addSuppressed(cleanup);
					}
				}
				throw failure;
			}
		}
		boolean error=out!=null && out.checkError();
		if(borrowed!=null){error|=borrowed.checkError();}
		if(error){throw new RuntimeException("Error writing to "+fname);}
	}

	/** Builds an owned encoder whose closure cannot close the borrowed PrintStream. */
	private static PrintWriter borrowedWriter(final PrintStream stream){
		assert(stream!=null) : "writeString passes a recognized borrowed standard stream.";
		Charset charset=Charset.defaultCharset();
		try{
			// Java18+ PrintWriter adopts PrintStream.charset(); older constructors use the default.
			charset=(Charset)PrintStream.class.getMethod("charset").invoke(stream);
		}catch(NoSuchMethodException e){
			// Keep the older PrintWriter constructor's default-charset behavior on Java8.
		}catch(ReflectiveOperationException e){
			throw new RuntimeException("Could not obtain the borrowed output stream charset", e);
		}
		return new PrintWriter(new BufferedWriter(new OutputStreamWriter(new NonClosingOutputStream(stream), charset)));
	}

	/** Borrows an output backend while allowing its owning encoder to finish normally. */
	private static final class NonClosingOutputStream extends FilterOutputStream{
		NonClosingOutputStream(final OutputStream out){super(out);}
		/** Delegates a block without FilterOutputStream's per-byte loop. */
		@Override
		public void write(byte[] bytes, int offset, int length) throws IOException{
			out.write(bytes, offset, length);
		}
		/** Flushes the borrowed backend; the caller retains responsibility for closing it. */
		@Override
		public void close() throws IOException{flush();}
	}

	/** Writes characters in the caller with append=false; failure policy follows the overload. */
	public static void writeStringAsync(CharSequence x, String fname){writeStringAsync(x, fname, false);}
	/**
	 * Writes in the calling thread using writeString's locking, cleanup and borrowed-stream rules.
	 * The InThread wrapper supplies a worker. A write failure halts the JVM through
	 * KillSwitch when assertions are enabled; with -da the original unchecked failure
	 * propagates without explicitly terminating the JVM. Counted workers release their
	 * reservation in finally when the JVM continues. No completion handle is returned.
	 * @param x Character sequence to write
	 * @param fname Output filename
	 * @param append Whether to append to existing file or overwrite
	 */
	public static void writeStringAsync(CharSequence x, String fname, boolean append){
		try{writeString(x, fname, append);}
		catch(RuntimeException | Error failure){
			//IO-042: worker failure must stop scientific runs; -da servers may survive log failure.
			assert(false) : KillSwitch.assertDie("Error writing to "+fname+": "+failure);
			throw failure;
		}
	}
	
	/**
	 * Writes an object to file using serialization with disk synchronization.
	 * Owns the returned backend and object wrapper. Attempts wrapper closure once and
	 * falls back to backend closure if wrapper construction or closure fails.
	 * Only distinct cleanup failures are suppressed onto the original failure.
	 *
	 * @param <X> Type of object to write
	 * @param x Object to serialize and write
	 * @param fname Output filename
	 * @param allowSubprocess Whether to allow subprocess compression
	 * @throws RuntimeException If output creation, serialization or closure fails
	 */
	public static <X> void write(X x, String fname, boolean allowSubprocess){
		if(verbose){System.err.println("write(x, "+fname+", "+allowSubprocess+")");}
		
		OutputStream os=getOutputStream(fname, false, true, allowSubprocess);
		
		try{

			synchronized(diskSync){
				ObjectOutputStream out=null;
				Throwable failure=null;
				try{
					out=new ObjectOutputStream(os);
					out.writeObject(x);
				}catch(Throwable t){
					failure=t;
					throw t;
				}finally{
					try{closeObjectOutput(out, os, failure);}
					catch(Throwable closing){
						if(failure==null){throw closing;}
						//Stored backend errors can be rethrown by close; never mask them with self-suppression.
						if(closing!=failure){failure.addSuppressed(closing);}
					}
				}
			}
			
		}catch(FileNotFoundException e){
			throw new RuntimeException(e);
		}catch(IOException e){
			throw new RuntimeException(e);
		}catch(OutOfMemoryError e){
			KillSwitch.memKill(e);
		}
	}
	
	/** Closes the wrapper, using the backend only if construction or wrapper closure failed.
	 * Avoids adding the original body failure back into a cleanup exception's suppressed list.
	 * @param out Object wrapper, or null when its construction failed
	 * @param backend Returned output backend, or null if the factory returned null
	 * @param primary Original construction or serialization failure, or null
	 * @throws IOException If closure fails */
	private static void closeObjectOutput(final ObjectOutputStream out, final OutputStream backend,
			final Throwable primary) throws IOException{
		if(out==null){
			if(backend!=null){backend.close();}
			return;
		}
		try{out.close();}
		catch(Throwable failure){
			try{if(backend!=null){backend.close();}}
			catch(Throwable cleanup){
				//BgzfOutputStreamMT/MT2 can rethrow one stored IOException at several boundaries.
				if(cleanup!=failure && cleanup!=primary){failure.addSuppressed(cleanup);}
			}
			throw failure;
		}
	}
	
	/**
	 * Serializes in the calling thread without diskSync; this method creates no worker.
	 * Owns the object wrapper/backend and attempts cleanup before reporting failure.
	 * Write failures halt the JVM only when assertions are enabled; with -da an
	 * unchecked failure propagates without explicitly terminating the JVM. The
	 * InThread wrapper releases its counted reservation when execution continues.
	 * Existing OutOfMemoryError handling remains unconditional.
	 * @param <X> Type of object to write
	 * @param x Object to serialize and write
	 * @param fname Output filename
	 * @param allowSubprocess Whether to allow subprocess compression
	 */
	public static <X> void writeAsync(X x, String fname, boolean allowSubprocess){
		if(verbose){System.err.println("writeAsync(x, "+fname+", "+allowSubprocess+")");}
		try{
			final OutputStream os=getOutputStream(fname, false, true, allowSubprocess);
			ObjectOutputStream out=null;
			Throwable failure=null;
			try{
				out=new ObjectOutputStream(os);
				out.writeObject(x);
			}catch(Throwable t){
				failure=t;
				throw t;
			}finally{
				try{closeObjectOutput(out, os, failure);}
				catch(Throwable closing){
					if(failure==null){throw closing;}
					if(closing!=failure){failure.addSuppressed(closing);}
				}
			}
		}catch(IOException failure){
			assert(false) : KillSwitch.assertDie("Error writing to "+fname+": "+failure);
			throw new RuntimeException(failure);
		}catch(OutOfMemoryError failure){
			KillSwitch.memKill(failure);
		}catch(RuntimeException | Error failure){
			//IO-042: the assertion message is deliberately unevaluated with -da.
			assert(false) : KillSwitch.assertDie("Error writing to "+fname+": "+failure);
			throw failure;
		}
	}
	
	/**
	 * Closes supplied readers first, then delegates backend/process completion.
	 * Reader IOExceptions are printed and accumulated; unchecked failures propagate.
	 * Reader entries must be nonnull. A wrapper can close System.in before the backend's
	 * direct-identity protection is reached; callers must account for borrowed ownership.
	 *
	 * @param is Input stream to close
	 * @param fname Filename for process identification
	 * @param killProcess Whether to kill associated subprocess
	 * @param ra Readers to close in order; null array is allowed except with verbose logging
	 * @return true if any errors occurred during closure
	 */
	public static final boolean finishReading(InputStream is, String fname, boolean killProcess, Reader...ra){
		if(verbose){System.err.println("finishReading("+is+", "+fname+", "+killProcess+", "+ra.length+")");}
//		assert(!killProcess);
		boolean error=false;
		if(ra!=null){
			for(Reader r : ra){
				try{
					r.close();
				}catch(IOException e){
					error=true;
					e.printStackTrace();
				}
			}
		}
		error|=finishReading(is, fname, killProcess);
		if(verbose){System.err.println("finishReading("+is+", "+fname+", "+killProcess+", "+ra.length+") returned "+error);}
		return error;
	}
	
	/**
	 * Closes a nonnull backend unless it is exactly System.in, then optionally waits
	 * for the filename-associated process. Caught close IOExceptions are printed and
	 * reported; unchecked failures propagate. A wrapper around System.in is not exempt.
	 *
	 * @param is Input stream to close
	 * @param fname Filename for process identification
	 * @param killProcess Whether to kill associated subprocess
	 * @return true if any errors occurred during closure
	 */
	public static final boolean finishReading(InputStream is, String fname, boolean killProcess){
		if(verbose){System.err.println("finishReading("+is+", "+fname+", "+killProcess+")");}
//		assert(!killProcess);
		boolean error=false;
		//FIXED [fileIO/ReadWrite#006]: null-guard 'is'. A streamer interrupted mid-shutdown (e.g. bbnorm ecc's
		//FastqScanStreamer, whose non-volatile 'is' may not yet be published when close() races in) can pass null
		//here; without the guard is.close() NPEs. A null stream simply means "nothing to close."
		if(is!=null && is!=System.in){
			try{
				is.close();
			}catch(IOException e){
				error=true;
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		if(killProcess && fname!=null && is!=null && is!=System.in){error|=ReadWrite.killProcess(fname);}
		if(verbose){System.err.println("finishReading("+is+", "+fname+", "+killProcess+") returned "+error);}
		return error;
	}
	
//	public static final boolean finishWriting(PrintWriter writer, OutputStream outStream, String fname){
//		return finishWriting(writer, outStream, fname, fname!=null);
//	}
	
	/**
	 * Checks PrintWriter status, closes the backend, then closes/checks the writer.
	 * Direct System.out/System.err backends are flushed but kept open; wrappers are
	 * not identified as borrowed. Finally requests filename-associated process completion.
	 * Combines boolean failures without short-circuiting; unchecked exceptions can abort cleanup.
	 *
	 * @param writer Print writer to close (may be null)
	 * @param outStream Nonnull backend to flush/close
	 * @param fname Filename for process identification
	 * @param killProcess Whether to kill associated subprocess
	 * @return True if the writer, stream closure or subprocess reports an error
	 */
	public static final boolean finishWriting(PrintWriter writer, OutputStream outStream, String fname, boolean killProcess){
		if(verbose){System.err.println("finishWriting("+writer+", "+outStream+" , "+fname+", "+killProcess+")");}
		boolean error=false;
		//FIO-033: checkError flushes and reports errors swallowed by PrintWriter.
		if(writer!=null){error|=writer.checkError();}
		error|=close(outStream);
		if(writer!=null && outStream!=System.out && outStream!=System.err){
			writer.close();
			error|=writer.checkError();
		}
		if(killProcess && fname!=null && outStream!=System.err && outStream!=System.out){error|=ReadWrite.killProcess(fname);}
		if(verbose){System.err.println("finishWriting("+writer+", "+outStream+" , "+fname+", "+killProcess+") returned "+error);}
		return error;
	}
	
//	/**
//	 * Closes output stream and kills associated subprocess if specified.
//	 * @param os Output stream to close
//	 * @param fname Filename for process identification
//	 * @return true if any errors occurred during closure
//	 */
//	public static final boolean close(OutputStream os, String fname){
//		if(verbose){System.err.println("close("+os+", "+fname+")");}
//		boolean error=false;
//		if(os!=null){error|=close(os);}
//		if(fname!=null && os!=System.err && os!=System.out){error|=killProcess(fname);}
//		if(verbose){System.err.println("close("+os+", "+fname+") returned "+error);}
//		return error;
//	}
	
	/**
	 * Flushes a nonnull backend and closes it unless it is exactly System.out/System.err.
	 * An exact ZipOutputStream also has its entry and archive finished before closure.
	 * IOExceptions are printed and accumulated; unchecked failures propagate. This does
	 * not detect PrintStream retained errors or wait for registered processes.
	 * @param os Nonnull output backend
	 * @return true if any errors occurred during closure
	 */
	public static final boolean close(OutputStream os){
		if(verbose){System.err.println("close("+os+")");}
		boolean error=false;
		try{
			os.flush();
		}catch(IOException e1){
			// TODO Auto-generated catch block
			e1.printStackTrace();
			error=true;
		}
		if(os.getClass()==ZipOutputStream.class){
			ZipOutputStream zos=(ZipOutputStream)os;
			try{
				zos.closeEntry();
				zos.finish();
			}catch(IOException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
				error=true;
			}
		}
//		else if(PROCESS_XZ && os.getClass()==org.tukaani.xz.XZOutputStream.class){
//			org.tukaani.xz.XZOutputStream zos=(org.tukaani.xz.XZOutputStream)os;
//			try {
//				zos.finish();
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		}
		if(os!=System.out && os!=System.err){
			try{
				os.close();
			}catch(IOException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
				error=true;
			}
		}
		if(verbose){System.err.println("close("+os+") returned "+error);}
		return error;
	}
	
	/**
	 * Creates output stream from FileFormat specification.
	 * @param ff FileFormat containing filename and options
	 * @param buffered Whether to use buffered output
	 * @return Configured output stream
	 */
	public static OutputStream getOutputStream(FileFormat ff, boolean buffered){
		return getOutputStream(ff.name(), ff.append(), buffered, ff.allowSubprocess());
	}

	/**
	 * Selects an output backend using case-sensitive suffixes and global codec settings.
	 * Normalizes backslashes and attempts parent creation before selection. The buffered
	 * flag applies only to paths that forward it. allowSubprocess is reduced when fewer
	 * than two threads are configured, but several codec/force branches ignore that hint.
	 * ZIP/BZ2/XZ/DSRC/FQZ/Alapy append restrictions here are assertions, not universal guards.
	 * BAM output has its own factory and is not selected by this method.
	 *
	 * @param fname Output filename with extension indicating compression type
	 * @param append Whether to append to existing file
	 * @param buffered Whether to use buffered output
	 * @param allowSubprocess Whether to allow external compression tools
	 * @return Configured output stream for the specified format
	 */
	public static OutputStream getOutputStream(String fname, boolean append, boolean buffered, boolean allowSubprocess){
		
		if(verbose){
			System.err.println("getOutputStream("+fname+", "+append+", "+buffered+", "+allowSubprocess+")");
			new Exception().printStackTrace(System.err);
		}
		
//		assert(false) : fname; //TODO: for testing
//		fname=fname.replaceAll("\\\\", "/");
		fname=fname.replace('\\', '/');
		assert(fname.indexOf('\\')<0);
//		assert(!fname.contains("//"));
		
		{//Create directories if needed.
			final int index=fname.lastIndexOf('/');
			if(index>0){
				File f=new File(fname.substring(0, index+1));
				if(!f.exists()){f.mkdirs();}
			}
		}
		
		boolean gzipped=gzipSuffix(fname);
		boolean zipped=fname.endsWith(".zip");
		boolean bzipped=PROCESS_BZ2 && fname.endsWith(".bz2");
		boolean xz=PROCESS_XZ && fname.endsWith(".xz");
		boolean dsrced=fname.endsWith(".dsrc");
		boolean fqz=USE_FQZ && fname.endsWith(".fqz");
		boolean alapy=USE_ALAPY && fname.endsWith(".ac");
		boolean zst=USE_ZSTD && fname.endsWith(".zst");//FIXED [fileIO/ReadWrite#002]: was USE_ALAPY (conflated zstd with Alapy) — now its own flag
		
//		assert(false) : fname;
		
		allowSubprocess=(allowSubprocess && Shared.threads()>1);
		
		if(gzipped){
//			assert(!append);
			return getGZipOutputStream(fname, append, allowSubprocess);
		}else if(zipped){
			assert(!append) : "Append is not allowed for zip archives.";
			return getZipOutputStream(fname, buffered, allowSubprocess);
		}else if(bzipped){
			assert(!append) : "Append is not allowed for bz2 archives.";//TODO: This might be OK; try it.
			return getBZipOutputStream(fname, buffered, append, allowSubprocess);
		}else if(zst){
//			assert(!append) : "Append is not allowed for zstd archives.";
			return getZstdStream(fname, append);
		}else if(xz){
			assert(!append) : "Append is not allowed for xz archives.";
			return getXZOutputStream(fname, buffered, allowSubprocess);
		}else if(dsrced){
			assert(!append) : "Append is not allowed for dsrc archives.";
			return getDsrcOutputStream(fname, buffered, allowSubprocess);
		}else if(fqz){
			assert(!append) : "Append is not allowed for fqz archives.";
			return getFqzStream(fname);
		}else if(alapy){
			assert(!append) : "Append is not allowed for alapy archives.";
			return getAlapyStream(fname);
		}
		return getRawOutputStream(fname, append, buffered);
	}
	
	/** Chooses native BAM output or samtools, forwarding append to the selected backend.
	 * Native construction IOExceptions are printed, then the method throws unsupported-output
	 * rather than falling back to samtools. Ignores MAX_SAMTOOLS_THREADS in this path. */
	public static OutputStream getBamOutputStream(String fname, boolean append){
		int zl=Tools.min(ZIPLEVEL, 6);
		if(nativeBamOut()){
			int threads=Tools.mid(1, Shared.threads(), zl>6 ? 16 : zl>4 ? 16 : zl>3 ? 8 : 4);
			try{
				//append must be forwarded: this path previously ignored it and truncated (2026-09-05)
				return new BamOutputStream(fname, zl, threads, append);
			}catch(IOException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}else if(Data.SAMTOOLS()){
			int threads=Tools.mid(1, Shared.threads(), zl>4 ? 16 : zl>3 ? 8 : 4);
			return getOutputStreamFromProcess(fname, "samtools view -@ "+threads+" -S -b -h - ", true, append, true, true);
		}else if(false && Data.SAMBAMBA()){
			return ReadWrite.getOutputStreamFromProcess(fname, "sambamba view -S -f bam -h ", true, append, true, true); //Sambamba does not support stdin
		}
		throw new RuntimeException("No bam output support available");
	}
	
	/**
	 * Creates uncompressed file output stream with special handling for stdout/stderr.
	 * Creates parent directories as needed and handles special filenames.
	 *
	 * @param fname Output filename ("stdout", "stderr", or regular path)
	 * @param append Whether to append to existing file
	 * @param buffered Whether to use buffered output
	 * @return Raw file output stream
	 */
	public static OutputStream getRawOutputStream(String fname, boolean append, boolean buffered){
		
		if(verbose){System.err.println("getRawOutputStream("+fname+", "+append+", "+buffered+")");}
		
		if(fname.equals("stdout") || fname.startsWith("stdout.")){
			return System.out;
		}else if(fname.equals("stderr") || fname.startsWith("stderr.")){
			return System.err;
		}else if(fname.startsWith("/dev/null/")){
			fname="/dev/null/";
		}
		
		if(fname.indexOf('|')>=0){fname=fname.replace('|', '_');}
		
		FileOutputStream fos=null;
		try{
			fos=new FileOutputStream(fname, append);
		}catch(FileNotFoundException e){
			synchronized(ReadWrite.class){
				try{
					File f=new File(fname);
					String parent=f.getParent();
					
					if(parent!=null){
						f=new File(parent);
						if(!f.exists()){
							boolean b=f.mkdirs();
							if(!b){System.err.println("Warning - could not create directory "+f.getAbsolutePath());}
						}
					}
					fos=new FileOutputStream(fname, append);
				}catch(Exception e2){
					throw new RuntimeException(e2);
				}
			}
		}
		assert(fos!=null);
		if(buffered){return new BufferedOutputStream(fos);}
		return fos;
	}
	
	/**
	 * Opens raw output, returning it in RAWMODE; otherwise throws for unsupported XZ.
	 * The unsupported path currently leaves that opened backend unclosed (FIO-046).
	 *
	 * @param fname Output filename
	 * @param buffered Whether to use buffered output
	 * @param allowSubprocess Unused
	 * @return Raw output only when RAWMODE is true
	 * @throws RuntimeException XZ format currently unsupported
	 */
	public static OutputStream getXZOutputStream(String fname, boolean buffered, boolean allowSubprocess){
		//TODO: Probable bug - FIO-046: the unsupported branch throws after opening raw output without closing it.
		final OutputStream raw=getRawOutputStream(fname, false, buffered);
		if(RAWMODE){return raw;}
		throw new RuntimeException("Unsupported format: XZ");
//		try {
//			org.tukaani.xz.LZMA2Options options = new org.tukaani.xz.LZMA2Options();
//			options.setPreset(ZIPLEVEL);
//			org.tukaani.xz.XZOutputStream out=new org.tukaani.xz.XZOutputStream(raw, options);
//			return out;
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		assert(false);
//		return null;
	}
	
	/**
	 * Creates bzip2 compressed output stream using external tools.
	 * Tries lbzip2, pbzip2, or bzip2 in order of preference.
	 *
	 * @param fname Output filename
	 * @param buffered Whether to use buffered output
	 * @param append Forwarded to external output; RAWMODE ignores it and opens with append=false
	 * @param allowSubprocess Unused; available external tools are selected regardless
	 * @return Bzip2 compressed output stream
	 * @throws RuntimeException if no bzip2 tools are available
	 */
	public static OutputStream getBZipOutputStream(String fname, boolean buffered, boolean append, boolean allowSubprocess){
		if(verbose){System.err.println("getBZipOutputStream("+fname+", "+buffered+", "+append+", "+allowSubprocess+")");}
//		assert(false) : ReadWrite.ZIPLEVEL+", "+Shared.threads()+", "+MAX_ZIP_THREADS+", "+ZIP_THREAD_MULT+", "+allowSubprocess+", "+USE_PIGZ+", "+Data.PIGZ();
		
		if(RAWMODE){
			final OutputStream raw=getRawOutputStream(fname, false, buffered);
			return raw;
		}
		
		if(USE_LBZIP2 && Data.LBZIP2()){return getLbzip2Stream(fname, append);}
		if(USE_PBZIP2 && Data.PBZIP2()){return getPbzip2Stream(fname, append);}
		if(USE_BZIP2 && Data.BZIP2()){return getBzip2Stream(fname, append);}
		
		throw new RuntimeException("bz2 compression not supported in this version, unless lbzip2, pbzip2 or bzip2 is installed.");
		
		
//		getBzip2Stream
		
//		{//comment to disable BZip2
//			try {
//				raw.write('B');
//				raw.write('Z');
//				CBZip2OutputStream out=new CBZip2OutputStream(raw, 8192);
//				return out;
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//			assert(false);
//			return null;
//		}
	}
	
	/**
	 * Creates DSRC compressed output stream for FASTQ files.
	 *
	 * @param fname Output filename
	 * @param buffered Whether to use buffered output
	 * @param append Forwarded to DSRC output; RAWMODE ignores it and opens with append=false
	 * @return DSRC compressed output stream
	 * @throws RuntimeException if DSRC tool not available
	 */
	public static OutputStream getDsrcOutputStream(String fname, boolean buffered, boolean append){
		if(verbose){System.err.println("getDsrcOutputStream("+fname+", "+buffered+", "+append+")");}
		if(RAWMODE){
			final OutputStream raw=getRawOutputStream(fname, false, buffered);
			return raw;
		}
		
		if(USE_DSRC && Data.DSRC() /*&& (Data.SH() || fname.equals("stdout") || fname.startsWith("stdout."))*/){return getDsrcOutputStream2(fname, append);}
		
		throw new RuntimeException("dsrc compression requires dsrc in the path.");
	}
	
	/**
	 * Creates ZIP compressed output stream using Java's built-in ZIP support.
	 * Borrowed stdout/stderr remain open; their retained errors become IOExceptions.
	 * Error checks flush the borrowed PrintStream after each write. Outside RAWMODE,
	 * the wrapper remains exactly ZipOutputStream for concrete-class checks.
	 *
	 * @param fname Output filename
	 * @param buffered Whether to use buffered output
	 * @param allowSubprocess Unused; this factory uses Java ZIP
	 * @return ZIP compressed output stream
	 * @throws RuntimeException If creation or entry setup fails; owned output cleanup is attempted
	 */
	public static OutputStream getZipOutputStream(String fname, boolean buffered, boolean allowSubprocess){
		if(verbose){System.err.println("getZipOutputStream("+fname+", "+buffered+", "+allowSubprocess+")");}
		final OutputStream raw=getRawOutputStream(fname, false, buffered);
		if(RAWMODE){return raw;}
		final OutputStream backend=(raw==System.out || raw==System.err ?
			new CheckedBorrowedOutputStream((PrintStream)raw) : raw);
		ZipOutputStream out=null;
		try{
			out=new ZipOutputStream(backend);
			out.setLevel(Tools.min(ZIPLEVEL, 9));
			final String basename=basename(fname);
			out.putNextEntry(new ZipEntry(basename));
			return out;
		}catch(IOException e){
			closeFailedCompression(out, backend, e);
			throw new RuntimeException(e);
		}catch(RuntimeException | Error e){
			closeFailedCompression(out, backend, e);
			throw e;
		}
	}

	/** Preserves borrowed ZIP backend ownership and exposes PrintStream's retained errors. */
	private static final class CheckedBorrowedOutputStream extends OutputStream{
		/** Captures the borrowed standard stream selected by the ZIP factory. */
		CheckedBorrowedOutputStream(final PrintStream stream_){
			assert(stream_!=null) : "ZIP borrowed backend must exist for byte delegation and error checks.";
			stream=stream_;
		}
		/** Delegates one byte and surfaces any retained write/flush error. */
		@Override
		public void write(final int value) throws IOException{
			stream.write(value);
			checkError();
		}
		/** Delegates a complete block without introducing a per-byte adapter loop. */
		@Override
		public void write(final byte[] bytes, final int offset, final int length) throws IOException{
			stream.write(bytes, offset, length);
			checkError();
		}
		/** checkError itself flushes the PrintStream before reporting its retained status. */
		@Override
		public void flush() throws IOException{checkError();}
		/** Completes this adapter without closing the caller-owned standard stream. */
		@Override
		public void close() throws IOException{flush();}
		/** Converts the borrowed PrintStream's silent error flag to the OutputStream contract. */
		private void checkError() throws IOException{
			if(stream.checkError()){throw new IOException("Error writing to borrowed standard output");}
		}
		/** Caller-owned byte sink; never closed by this adapter. */
		private final PrintStream stream;
	}
	
	/**
	 * Selects native/external BGZF, pigz or gzip using global preferences, then falls
	 * back to Java GZIP. Forced/preferred BGZF is considered before allowSubprocess;
	 * FORCE_PIGZ also enters the external-selection branch despite that hint.
	 * RAWMODE is checked only after earlier selection paths, not as a universal bypass.
	 *
	 * @param fname Output filename
	 * @param append Whether to append to existing file
	 * @param allowSubprocess Whether to allow external compression tools
	 * @return Gzip compressed output stream
	 * @throws RuntimeException If Java gzip construction fails; owned output cleanup is attempted
	 */
	public static OutputStream getGZipOutputStream(String fname, boolean append, boolean allowSubprocess){
		if(verbose){System.err.println("getGZipOutputStream("+fname+", "+append+", "+allowSubprocess+"); "+FORCE_BGZIP+", "+USE_BGZIP+", "+Data.BGZIP()+", "+USE_PIGZ+", "+USE_GZIP+", "+RAWMODE);}
		final boolean bgzip=(USE_BGZF && (ALLOW_NATIVE_BGZF || (USE_BGZIP && Data.BGZIP())));
//		assert(false) : fname+", "+allowSubprocess+", "+bgzip+", "+ZIPLEVEL;
		if(bgzip && (FORCE_BGZIP || (PREFER_BGZIP && ZIPLEVEL<10))){return getBgzipStream(fname, append);}
		if(FORCE_PIGZ || (allowSubprocess && Shared.threads()>=2)){
			if((fname.endsWith(".vcf.gz") || fname.endsWith(".sam.gz") || (PREFER_BGZIP && ZIPLEVEL<10)) && bgzip){return getBgzipStream(fname, append);}
			if(USE_PIGZ && Data.PIGZ()){return getPigzStream(fname, append);}
			if(bgzip){return getBgzipStream(fname, append);}
			if(USE_GZIP && Data.GZIP()/* && (Data.SH() /*|| fname.equals("stdout") || fname.startsWith("stdout."))*/){return getGzipStream(fname, append);}
		}
		final OutputStream raw=getRawOutputStream(fname, append, false);
		if(RAWMODE){return raw;}
		try{
			final GZIPOutputStream out=new GZIPOutputStream(raw, 8192){
				{
					//					        def.setLevel(Deflater.DEFAULT_COMPRESSION);
					def.setLevel(Tools.min(ZIPLEVEL, 9));
				}
			};
			return out;
		}catch(IOException e){
			closeFailedCompression(null, raw, e);
			throw new RuntimeException(e);
		}catch(RuntimeException | Error e){
			closeFailedCompression(null, raw, e);
			throw e;
		}
	}

	/** Closes owned output after factory failure without replacing the primary cause.
	 * An assigned wrapper is tried once; raw closure is the fallback if it cannot close.
	 * Borrowed stdout/stderr are left open, including when hidden beneath a wrapper.
	 * @param wrapper Constructed wrapper, or null if construction did not return
	 * @param raw Underlying output opened by the factory
	 * @param primary Construction or setup failure that will be propagated */
	private static void closeFailedCompression(final OutputStream wrapper, final OutputStream raw,
			final Throwable primary){
		if(raw==System.out || raw==System.err){return;}
		if(wrapper!=null){
			try{
				wrapper.close();
				return;
			}catch(Throwable closing){
				if(closing!=primary){primary.addSuppressed(closing);}
			}
		}
		try{raw.close();}
		catch(Throwable closing){
			if(closing!=primary){primary.addSuppressed(closing);}
		}
	}
	
	/**
	 * Creates pigz (parallel gzip) output stream with configured thread and
	 * compression-level adjustments based on configured thread count.
	 *
	 * @param fname Output filename
	 * @param append Whether to append to existing file
	 * @return Pigz compressed output stream
	 */
	public static OutputStream getPigzStream(String fname, boolean append){
		if(verbose){System.err.println("getPigzStream("+fname+")");}
//		System.err.println(MAX_ZIP_THREADS); //123
		int threads=Tools.min(MAX_ZIP_THREADS, Tools.max((int)((Shared.threads()+1)*ZIP_THREAD_MULT), 1));
//		System.err.println(threads); //123
		threads=Tools.max(1, Tools.min(Shared.threads(), threads));
//		System.err.println(threads); //123
		int zl=(ZIPLEVEL<10 ? ZIPLEVEL : Data.PIGZ_VERSION_23plus ? 11 : 9);
		if(ALLOW_ZIPLEVEL_CHANGE && threads>=4 && zl>0 && zl<4){zl=4;}
		if(zl<3){threads=Tools.min(threads, 12);}
		else if(zl<5){threads=Tools.min(threads, 24);}
		else if(zl<7){threads=Tools.min(threads, 40);}
//		System.err.println(threads); //123
		OutputStream out;
		String command="pigz -c -p "+threads+" -"+zl;
		if(PIGZ_BLOCKSIZE!=128){
			command=command+" -b "+PIGZ_BLOCKSIZE;
		}
		if(PIGZ_ITERATIONS>0 && Data.PIGZ_VERSION_231plus){
			command=command+" -I "+PIGZ_ITERATIONS;
		}
		
//		System.err.println("*** "+command);
		
		//Sample command on command line, without piping: pigz -11 -k -f -b 256 -I 25000 file.fa
//		assert(false) : MAX_ZIP_THREADS+", "+Shared.threads()+", "+ZIP_THREAD_MULT+", "+ZIPLEVEL+", "+command;
		out=getOutputStreamFromProcess(fname, command, true, append, true, true);
		
		return out;
	}
	
	/**
	 * Creates FQZ compressed output stream for FASTQ files.
	 * @param fname Output filename
	 * @return FQZ compressed output stream
	 */
	public static OutputStream getFqzStream(String fname){
		if(verbose){System.err.println("getFqzStream("+fname+")");}
		String command="fqz_comp -s"+Tools.mid(1, ZIPLEVEL, 8)+"+"; //9 gives bad compression
		if(ZIPLEVEL>5){command=command+" -q3";}
		//if(ZIPLEVEL>5){command=command+" -b -q3";} //b does not seem to work
		OutputStream out=getOutputStreamFromProcess(fname, command, true, false, true, true);
		return out;
	}
	
	/**
	 * Creates Alapy compressed output stream with compression level selection.
	 * @param fname Output filename
	 * @return Alapy compressed output stream
	 */
	public static OutputStream getAlapyStream(String fname){
		if(verbose){System.err.println("getAlapyStream("+fname+")");}
		String compression=(ZIPLEVEL>6 ? "-l best" : ZIPLEVEL<4 ? "-l fast" : "-l medium");
//		String compression="";
		String command="alapy_arc "+compression+" -n "+fname+" -q -c - ";
		OutputStream out=getOutputStreamFromProcess(fname, command, true, false, true, false);
		return out;
	}
	
	/**
	 * Creates standard gzip output stream using external gzip command.
	 * @param fname Output filename
	 * @param append Whether to append to existing file
	 * @return Gzip compressed output stream
	 */
	public static OutputStream getGzipStream(String fname, boolean append){
		if(verbose){System.err.println("getGzipStream("+fname+")");}
		OutputStream out=getOutputStreamFromProcess(fname, "gzip -c -"+Tools.min(ZIPLEVEL, 9), true, append, true, true);
		return out;
	}
	
	/**
	 * Creates bgzip output stream with configured thread allocation for
	 * block-gzip format compatible with tabix indexing.
	 *
	 * @param fname Output filename
	 * @param append Whether to append to existing file
	 * @return Bgzip compressed output stream
	 */
	public static OutputStream getBgzipStream(String fname, boolean append){
		if(verbose){System.err.println("getBgzipStream("+fname+")");}
		
		int threads=Tools.min(MAX_ZIP_THREADS, Tools.max((int)((Shared.threads()+1)*ZIP_THREAD_MULT), 1));
		threads=Tools.max(1, Tools.min(Shared.threads(), threads));
		int zl=Tools.mid(ZIPLEVEL, 1, 9);
//		System.err.println("A: ZIPLEVEL="+ZIPLEVEL+", ALLOW_CHANGE="+ALLOW_ZIPLEVEL_CHANGE+", zl="+zl);
		if(!USE_BGZIP || nativeBgzfOut() || !Data.BGZIP()){
			if(zl>5){zl=5;}//Required for native bgzip.
			if(ALLOW_ZIPLEVEL_CHANGE){
				if(zl<4 && zl>0 && threads>=16){zl=4;}
			}
			if(zl<3){threads=Tools.min(threads, 16);}
			else if(zl<5){threads=Tools.min(threads, 24);}
			else if(zl<6){threads=Tools.min(threads, 64);}
//			System.err.println("B: ZIPLEVEL="+ZIPLEVEL+", ALLOW_CHANGE="+ALLOW_ZIPLEVEL_CHANGE+", zl="+zl);
			final OutputStream raw=getRawOutputStream(fname, append, false);//TODO - should it be true or false?
			if(RAWMODE){return raw;}
			final OutputStream out;
			if(!BgzfSettings.USE_MULTITHREADED_BGZF){out=new BgzfOutputStream(raw);}
			else if(BgzfSettings.USE_BGZFOS_MT2){
				out=new BgzfOutputStreamMT2(raw, Tools.mid(1, 64, threads), zl);
			}else{
				out=new BgzfOutputStreamMT(raw, Tools.mid(1, 64, threads), zl);
			}
			return out;
		}
		
		if(ALLOW_ZIPLEVEL_CHANGE && threads>=4 && zl>0 && zl<4){zl=4;}
		if(zl<3){threads=Tools.min(threads, 24);}
		else if(zl<5){threads=Tools.min(threads, 48);}//RandomReadsMG can go way over 16
		else if(zl<7){threads=Tools.min(threads, 64);}

//		System.err.println("C: ZIPLEVEL="+ZIPLEVEL+", ALLOW_CHANGE="+ALLOW_ZIPLEVEL_CHANGE+", zl="+zl);
//		if(threads<32 && zl>6 || true) {
//			System.err.println("***BGZIP: "+threads+", "+zl+", "+Shared.threads()+", "+MAX_ZIP_THREADS+", "+ZIP_THREAD_MULT);
//			new Exception().printStackTrace();
//		}
//		assert(false) : Data.BGZIP()+", "+Data.PIGZ();
		String command="bgzip -c "+(append ? "" : "-f ")+(Data.BGZIP_VERSION_levelFlag ? "-l "+zl+" " : "")+(Data.BGZIP_VERSION_threadsFlag ? "-@ "+threads+" " : "");
		if(verbose){System.err.println(command);}
		OutputStream out=getOutputStreamFromProcess(fname, command, true, append, true, true);
		if(verbose){System.err.println("fetched bgzip stream.");}
		return out;
	}
	
//	Usage :
//	      zstd [args] [FILE(s)] [-o file]
//
//	FILE    : a filename
//	          with no FILE, or when FILE is - , read standard input
//	Arguments :
//	 -#     : # compression level (1-19, default: 3)
//	 -d     : decompression
//	 -D DICT: use DICT as Dictionary for compression or decompression
//	 -o file: result stored into `file` (only 1 output file)
//	 -f     : disable input and output checks. Allows overwriting existing files,
//	          input from console, output to stdout, operating on links,
//	          block devices, etc.
//	--rm    : remove source file(s) after successful de/compression
//	 -k     : preserve source file(s) (default)
//	 -h/-H  : display help/long help and exit
//
	
	//This works correctly, but decompression doesn't
	/**
	 * Creates Zstandard compressed output stream with configured thread selection.
	 * @param fname Output filename
	 * @param append Whether to append to existing file
	 * @return Zstd compressed output stream
	 */
	public static OutputStream getZstdStream(String fname, boolean append){
		if(verbose){System.err.println("getZstdStream("+fname+")");}
		
		int threads=Tools.min(MAX_ZIP_THREADS, Tools.max((int)((Shared.threads()+1)*ZIP_THREAD_MULT), 1));
		threads=Tools.max(1, Tools.min(Shared.threads(), threads));
		int zl=Tools.mid(ZIPLEVEL, 1, 19);
		if(ALLOW_ZIPLEVEL_CHANGE && threads>=4 && zl>0 && zl<4){zl=4;}
		if(zl<3){threads=Tools.min(threads, 12);}
		else if(zl<5){threads=Tools.min(threads, 16);}
		else if(zl<7){threads=Tools.min(threads, 32);}
		
		String command="zstd -"+zl+" -T"+threads;
		if(verbose){System.err.println(command);}
		OutputStream out=getOutputStreamFromProcess(fname, command, true, append, true, true);
		if(verbose){System.err.println("fetched zstd stream.");}
		return out;
	}
	
	/**
	 * Creates bzip2 output stream using external bzip2 command.
	 * @param fname Output filename
	 * @param append Whether to append to existing file
	 * @return Bzip2 compressed output stream
	 */
	public static OutputStream getBzip2Stream(String fname, boolean append){
		if(verbose){System.err.println("getBzip2Stream("+fname+")");}
		String command="bzip2 -c -"+Tools.min(BZIPLEVEL, 9);
		OutputStream out=getOutputStreamFromProcess(fname, command, true, append, true, true);
		return out;
	}
	
	/**
	 * Creates parallel bzip2 output stream with configured thread selection.
	 * @param fname Output filename
	 * @param append Whether to append to existing file
	 * @return Pbzip2 compressed output stream
	 */
	public static OutputStream getPbzip2Stream(String fname, boolean append){
		if(verbose){System.err.println("getPbzip2Stream("+fname+")");}
		int threads=Tools.min(MAX_ZIP_THREADS, Tools.max((int)((Shared.threads()+1)*ZIP_THREAD_MULT), 1));
		threads=Tools.max(1, Tools.min(Shared.threads(), threads));
		String command="pbzip2 -c -p"+threads+" -"+Tools.min(BZIPLEVEL, 9);
		OutputStream out=getOutputStreamFromProcess(fname, command, true, append, true, true);
		return out;
	}
	
	/**
	 * Creates lbzip2 output stream for bzip2 compression.
	 * @param fname Output filename
	 * @param append Whether to append to existing file
	 * @return Lbzip2 compressed output stream
	 */
	public static OutputStream getLbzip2Stream(String fname, boolean append){
		if(verbose){System.err.println("getLbzip2Stream("+fname+")");}
		String command="lbzip2 -"+Tools.min(BZIPLEVEL, 9);
		OutputStream out=getOutputStreamFromProcess(fname, command, true, append, true, true);
		return out;
	}
	
	/**
	 * Creates DSRC output stream with compression parameter selection
	 * based on compression level. DSRC is specialized for FASTQ files.
	 *
	 * @param fname Output filename
	 * @param append Whether to append to existing file
	 * @return DSRC compressed output stream
	 */
	public static OutputStream getDsrcOutputStream2(String fname, boolean append){
		if(verbose){System.err.println("getDsrcOutpustream2("+fname+")");}
		int threads=Tools.min(MAX_ZIP_THREADS, Tools.max((int)((Shared.threads()+1)*ZIP_THREAD_MULT), 1));
		threads=Tools.max(1, Tools.min(Shared.threads()-1, threads));
		String params=null;
		if(ZIPLEVEL<=2){
			params="-d0 -q0 -b8";
		}else if(ZIPLEVEL<=4){
			params="-d1 -q1 -b16";
		}else if(ZIPLEVEL<=8){
			params="-d2 -q2 -b32";
		}else{
			params="-d3 -q2 -b64";
		}
		//TODO [fileIO/ReadWrite#005] DSRC is abandoned (not installed, interface forgotten, effectively unused) ->
		//safe_to_delete candidate. KNOWN ISSUE, intentionally NOT fixed (dsrc CLI unverifiable, dsrc dead): fname is
		//appended to 'command' here AND again at the call below (command+" "+fname), so the command carries fname twice.
		String command="dsrc c -t"+threads+" "+params+" -s";
		if(fname.equals("stdout") || fname.startsWith("stdout.")){
			//???
			assert(false) : "Undefined dsrc option.";
		}else{
			command+=" "+fname;
		}
		//[fileIO/ReadWrite#005a FIXED]: removed a leftover 'System.err.println(command);//123' debug print here.
//		OutputStream out=getOutputStreamFromProcess(fname, command, true, append, true);
		OutputStream out=getOutputStreamFromProcess(fname, command+" "+fname, true, append, true, false);
		return out;
	}
	
	/** Starts a named output process and returns its stdin pipe; callers must close it
	 * before completing the registered process with killProcess. ProcessBuilder splits
	 * command on literal spaces, inherits stderr and redirects stdout when useFname is true;
	 * stdout aliases inherit stdout regardless of useFname. The legacy Runtime path may
	 * use sh redirection. Startup IOExceptions are printed, followed by an assertion or
	 * null dereference if startup failed. Duplicate filename registration is fatal. */
	public static OutputStream getOutputStreamFromProcess(final String fname, final String command, boolean sh, boolean append, boolean useProcessBuilder, boolean useFname){
		if(verbose){System.err.println("getOutputStreamFromProcess("+fname+", "+command+", "+sh+", "+useProcessBuilder+")");}
		
		OutputStream out=null;
		Process p=null;
		if(useProcessBuilder){
			ProcessBuilder pb=new ProcessBuilder();
			pb.redirectError(Redirect.INHERIT);
			
			if(fname.equals("stdout") || fname.startsWith("stdout.")){
				pb.redirectOutput(Redirect.INHERIT);
				pb.command(command.split(" "));
			}else{
				
				if(useFname){
					if(append){
						pb.redirectOutput(ProcessBuilder.Redirect.appendTo(new File(fname)));
					}else{
						pb.redirectOutput(new File(fname));
					}
				}
				
				pb.command(command.split(" "));
			}
			try{
				p=pb.start();
			}catch(IOException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			assert(p!=null) : "Could not execute "+command;
			addProcess(fname, p);
			out=p.getOutputStream();
			{
				out=p.getOutputStream();
				InputStream es=p.getErrorStream();
				assert(es!=null);
				PipeThread et=new PipeThread(es, System.err);
				addPipeThread(fname, et);
				et.start();
			}
			return out;
		}
		
		if(fname.equals("stdout") || fname.startsWith("stdout.")){
			try{
				p=Runtime.getRuntime().exec(command);
			}catch(IOException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			assert(p!=null) : "Could not execute "+command;
			InputStream is=p.getInputStream();
			PipeThread it=new PipeThread(is, System.out);
			addPipeThread(fname, it);
			it.start();
//		}else if(fname.equals("stderr") || fname.startsWith("stderr.")){
//			try {
//				p = Runtime.getRuntime().exec(command);
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//			InputStream is=p.getErrorStream();
//			PipeThread it=new PipeThread(is, System.err);
//			it.start();
		}else{
			try{
				if(sh){
					String[] cmd={
							"sh",
							"-c",
							command+(useFname ? " 1"+(append ? ">>" : ">")+fname : "")
					};
					p=Runtime.getRuntime().exec(cmd);
				}else{
					//TODO: append won't work here...
					assert(false) : command;
					p=Runtime.getRuntime().exec(command);
				}
			}catch(IOException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		assert(p!=null) : "Could not execute "+command;
		addProcess(fname, p);
		out=p.getOutputStream();
		InputStream es=p.getErrorStream();
		assert(es!=null);
		PipeThread et=new PipeThread(es, System.err);
		addPipeThread(fname, et);
		et.start();

		return out;
	}
	
	/**
	 * Reads decoded lines using the default charset and appends one LF after every line.
	 * Original line endings and the absence of a final terminator are not preserved.
	 * Raw, Java ZIP/GZIP, native BGZF, bzip2 and Zstandard inputs receive normal and exceptional cleanup.
	 * Java/native wrappers protect borrowed stdin; bzip2 and Zstandard also collect decoder exit status.
	 * Other codecs retain legacy normal-only closure without process completion.
	 * Closing wrappers does not preserve bytes already consumed by their read-ahead.
	 * @param fname Input filename
	 * @return Complete file content as string
	 */
	public static String readString(String fname){
		//TODO: FIO-046: other codec cleanup and process completion remain separate.
		if(verbose){System.err.println("readString("+fname+")");}
		String x=null;
		final HelperInput input=new HelperInput();
		//Failed acquisition remains the factory's responsibility.
		input.stream=getInputStream(fname, false, false, true, input);
		try{
			StringBuilder sb=new StringBuilder();
			InputStreamReader decoder=new InputStreamReader(input.stream);
			input.outer=decoder;
			BufferedReader in=new BufferedReader(decoder, INBUF);
			input.outer=in;
			String temp=in.readLine();
			while(temp!=null){
				sb.append(temp).append('\n');
				temp=in.readLine();
			}
			input.closeNormally();
			x=sb.toString();
		}catch(FileNotFoundException e){
			input.closeAfterFailure(e);
			throw new RuntimeException(e);
		}catch(IOException e){
			input.closeAfterFailure(e);
			throw new RuntimeException(e);
		}catch(OutOfMemoryError e){
			input.closeAfterFailure(e);
			KillSwitch.memKill(e);
		}catch(RuntimeException | Error e){
			input.closeAfterFailure(e);
			throw e;
		}
		return x;
	}
	
	/**
	 * Reads one serialized object. Raw, Java ZIP/GZIP, native BGZF, bzip2 and external
	 * pigz/bgzip and Zstandard inputs receive normal and exceptional cleanup. Java/native wrappers
	 * protect borrowed stdin; bzip2, pigz/bgzip and Zstandard also collect decoder status. Other
	 * codecs retain normal-only closure without process completion. Reader prefetch
	 * may consume bytes beyond the object.
	 * @param fname Input filename
	 * @param allowSubprocess Whether to allow subprocess decompression
	 * @return Deserialized object, which may be null
	 */
	public static Object readObject(String fname, boolean allowSubprocess){
		//TODO: FIO-046: other codec cleanup and process completion remain separate.
		if(verbose){System.err.println("readObject("+fname+")");}
		Object x=null;
		final HelperInput input=new HelperInput();
		//Failed acquisition remains the factory's responsibility.
		input.stream=getInputStream(fname, true, allowSubprocess, true, input);
		try{
			ObjectInputStream in=new ObjectInputStream(input.stream);
			input.outer=in;
			x=in.readObject();
			input.closeNormally();
		}catch(IOException e){
			input.closeAfterFailure(e);
			throw new RuntimeException(e);
		}catch(ClassNotFoundException e){
			input.closeAfterFailure(e);
			throw new RuntimeException(e);
		}catch(OutOfMemoryError e){
			input.closeAfterFailure(e);
			KillSwitch.memKill(e);
		}catch(RuntimeException | Error e){
			input.closeAfterFailure(e);
			throw e;
		}
		return x;
	}
	
	/**
	 * Selects an input backend with ordered BAM decoding requested by default.
	 * See the four-argument overload for codec selection and ownership.
	 *
	 * @param fname Input filename with extension indicating compression type
	 * @param buffer Whether to use buffered input
	 * @param allowSubprocess Whether to allow external decompression tools
	 * @return Configured input stream for the specified format
	 */
	public static InputStream getInputStream(String fname, boolean buffer, boolean allowSubprocess){
		return getInputStream(fname, buffer, allowSubprocess, true);
	}
	
	/** Selects decompression using case-sensitive suffixes and global flags, unless RAWMODE.
	 * Only raw fallback honors buffer directly; codec branches choose their own buffering.
	 * ordered is forwarded only to native BAM. allowSubprocess is a hint honored by some
	 * paths, not a guarantee that no external tool starts. Caller owns completion/closure. */
	public static InputStream getInputStream(String fname, boolean buffer, boolean allowSubprocess, 
			boolean ordered){
		return getInputStream(fname, buffer, allowSubprocess, ordered, null);
	}

	/** Carries helper ownership through raw, Java/native gzip, pigz/bgzip, bzip2 and Zstandard acquisition. */
	private static InputStream getInputStream(String fname, boolean buffer, boolean allowSubprocess,
			boolean ordered, final HelperInput helper){
		if(verbose){System.err.println("getInputStream("+fname+", "+buffer+", "+allowSubprocess+")");}
		//[fileIO/ReadWrite#003 FIXED]: removed a dead 'boolean xz' local here — XZ input is unsupported (PROCESS_XZ=false,
		//getXZInputStream returns null) and the local was never referenced below. No order-sensitive change (it was inert).
		boolean gzipped=gzipSuffix(fname);
		boolean zipped=fname.endsWith(".zip");
		boolean bzipped=PROCESS_BZ2 && fname.endsWith(".bz2");
		boolean dsrced=fname.endsWith(".dsrc");
		boolean bam=fname.endsWith(".bam") && Data.BAM_SUPPORT_IN();
		boolean fqz=fname.endsWith(".fqz");
		boolean alapy=fname.endsWith(".ac");
		boolean zst=USE_ZSTD && fname.endsWith(".zst");//[fileIO/ReadWrite#002 FIXED]: gate on USE_ZSTD (was ungated; output was on USE_ALAPY) — now symmetric
		
		allowSubprocess=(allowSubprocess && Shared.threads()>1);
		
		if(!RAWMODE){
			if(zipped){return getZipInputStream(fname, true, helper);}
			if(gzipped){return getGZipInputStream(fname, allowSubprocess, false, helper);}
			if(bzipped){return getBZipInputStream(fname, allowSubprocess, helper);}
			if(dsrced){return getDsrcInputStream(fname);}
			if(zst){return getUnzstdStream(fname, helper);}
			if(bam){return getBamInputStream(fname, ordered);}

			if(fqz){return getInputStreamFromProcess(fname, "fqz_comp -d ", false, true, true);}
			if(alapy){
				return getInputStreamFromProcess(fname, "alapy_arc -q -d "+fname+" -", false, false, true);
			}
		}
		
		return getRawInputStream(fname, buffer, helper);
	}
	
	/** Chooses native BAM decoding, then enabled sambamba, otherwise samtools.
	 * ordered affects only native decoding; external tools receive SAM flag filtering.
	 * Returned external streams require registered-process completion by the caller. */
	public static InputStream getBamInputStream(String fname, boolean ordered){
		if(nativeBamIn()){
			return new BamInputStream(fname, ordered);
		}else if(SAMBAMBA()){
			String command="sambamba -q view -h";//Adding -t 8 did not change speed
//			new Exception().printStackTrace(); //123
			if(SAMTOOLS_IGNORE_FLAG!=0){
				command=command+" --num-filter=0/"+SAMTOOLS_IGNORE_FLAG;
			}
			return getInputStreamFromProcess(fname, command, false, true, true);
		}else{
			String command="samtools view -h";//Adding -@ 4 or 16 did not change speed
//			new Exception().printStackTrace(); //123
			if(SAMTOOLS_IGNORE_FLAG!=0){
				//					command=command+" -F 4";
				command=command+" -F 0x"+Integer.toHexString(SAMTOOLS_IGNORE_FLAG);
			}
			String version=Data.SAMTOOLS_VERSION;
			if(Shared.threads()>1 && version!=null && version.startsWith("1.") && version.length()>2){
				try{
					String[] split=version.split("\\.");
					int number=-1;
					try{
						number=Integer.parseInt(split[1]);
					}catch(Exception e){}
					if(number<0){
						try{
							number=Integer.parseInt(split[1].substring(0, 1));
						}catch(Exception e1){}
					}
					if(number>3){
						command=command+" -@ 2";
					}
				}catch(NumberFormatException e){
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			//				System.err.println(command);
			return getInputStreamFromProcess(fname, command, false, true, true);
		}
	}
	
	/**
	 * Creates uncompressed file input stream with support for stdin, JAR resources,
	 * and automatic file path resolution. Handles special filenames and missing files.
	 *
	 * @param fname Input filename ("stdin" or regular path, JAR URLs supported)
	 * @param buffer Whether to use buffered input
	 * @return Raw file input stream
	 */
	public static InputStream getRawInputStream(String fname, boolean buffer){
		return getRawInputStream(fname, buffer, null);
	}

	/** Marks actual raw acquisition and optionally shields the borrowed standard input. */
	private static InputStream getRawInputStream(String fname, boolean buffer, final HelperInput helper){
		if(verbose){System.err.println("getRawInputStream("+fname+", "+buffer+")");}
		
		assert(fname!=null);
		fname=fname.replace('\\', '/');
		assert(fname.indexOf('\\')<0);
		assert(!fname.contains("\\\\"));
//		assert(!fname.contains("//")) : fname;
		
		final boolean jar=fname.startsWith("jar:");
		
		if(!jar){
			boolean failed=false;
			File f=new File(fname);
			if(!f.exists()){
				String f2=fname.toLowerCase();
				if(f2.equals("stdin") || f2.startsWith("stdin.")){
					//				System.err.println("Returning stdin: A");
					if(helper==null){return System.in;}
					helper.managed=true;
					return new NonClosingInputStream(System.in);
				}
				
				if(fname.indexOf('/')<0){
					f2=Data.ROOT_CURRENT+"/"+fname;
					if(!new File(f2).exists()){
						failed=true;
					}else{
						fname=f2;
					}
				}else{
					failed=true;
				}
			}
//			if(failed){throw new RuntimeException("Can't find file "+fname);}
			if(failed){
				shared.KillSwitch.exceptionKill(new RuntimeException("Can't find file "+fname));
			}
		}
		
//		System.err.println("Getting input stream for "+fname);
//		assert(!fname.contains("\\"));
//		assert(!loadedFiles.contains(fname)) : "Already loaded "+fname;
//		loadedFiles.add(fname);
		
		InputStream in=null;
		if(jar){
			try{
				
				URL url=new URL(fname);
				
				InputStream is=url.openStream();

				if(buffer){
					BufferedInputStream bis=new BufferedInputStream(is, INBUF);
					in=bis;
				}else{
					in=is;
				}

			}catch(FileNotFoundException e){
				System.err.println("Error when attempting to read "+fname);
				throw new RuntimeException(e);
			}catch(MalformedURLException e){
				System.err.println("Error when attempting to read "+fname);
				throw new RuntimeException(e);
			}catch(IOException e){
				System.err.println("Error when attempting to read "+fname);
				throw new RuntimeException(e);
			}
		}else{
			try{

				FileInputStream fis=new FileInputStream(fname);

				if(buffer){
					BufferedInputStream bis=new BufferedInputStream(fis, INBUF);
					in=bis;
				}else{
					in=fis;
				}

			}catch(FileNotFoundException e){
				throw new RuntimeException(e);
			}
		}
		
		if(helper!=null){helper.managed=true;}
		return in;
	}
	
	/** Creates ZIP input stream with default buffering */
	public static InputStream getZipInputStream(String fname){return getZipInputStream(fname, true);}
	/**
	 * Opens Java ZIP input at its first entry. Assertions require an entry whose name
	 * equals basename(fname); later entries are not traversed here. Setup failures
	 * attempt owned wrapper/backend cleanup while preserving the original failure.
	 * Borrowed System.in is left open, including beneath a partially initialized wrapper.
	 *
	 * @param fname Input ZIP filename
	 * @param buffer Whether to use buffered input
	 * @return ZIP input stream positioned at first entry
	 */
	public static InputStream getZipInputStream(String fname, boolean buffer){
		return getZipInputStream(fname, buffer, null);
	}

	/** Acquires the Java ZIP backend using the helper's optional borrowed-input shield. */
	private static InputStream getZipInputStream(String fname, boolean buffer, final HelperInput helper){
		if(verbose){System.err.println("getZipInputStream("+fname+", "+buffer+")");}
		InputStream raw=getRawInputStream(fname, buffer, helper);
		InputStream in=null;

		try{
			final String basename=basename(fname);
			ZipInputStream zis=new ZipInputStream(raw);
			in=zis;//FIO-046: retain ownership before entry parsing or assertions can fail.
			ZipEntry ze=zis.getNextEntry();
			assert(ze!=null);
			assert(basename.equals(ze.getName())) : basename+" != "+ze.getName();

		}catch(FileNotFoundException e){
			closeFailedInput(in, raw, e);
			System.err.println("Error when attempting to read "+fname);
			throw new RuntimeException(e);
		}catch(IOException e){
			closeFailedInput(in, raw, e);
			System.err.println("Error when attempting to read "+fname);
			throw new RuntimeException(e);
		}catch(RuntimeException | Error e){
			closeFailedInput(in, raw, e);
			throw e;
		}

		return in;
	}
	
	/**
	 * Selects native/external decompression from configured threads, low-memory mode,
	 * JAR status and codec preferences, then falls back to Java GZIP. Native BGZF may
	 * be selected without subprocess permission. Java-wrapper construction failures
	 * attempt owned backend cleanup; borrowed System.in and the original failure are preserved.
	 *
	 * @param fname Input filename
	 * @param allowSubprocess Whether to allow external decompression tools
	 * @param buffer Whether to use buffered input
	 * @return Gzip decompressed input stream
	 */
	public static InputStream getGZipInputStream(String fname, boolean allowSubprocess, boolean buffer){
		return getGZipInputStream(fname, allowSubprocess, buffer, null);
	}

	/** Preserves codec selection; Java/native BGZF and external pigz/bgzip carry helper ownership. */
	private static InputStream getGZipInputStream(String fname, boolean allowSubprocess, boolean buffer,
			final HelperInput helper){
		if(verbose){System.err.println("getGZipInputStream("+fname+", "+allowSubprocess+")");}
		
		if(Shared.threads()<2 || Shared.LOW_MEMORY || fname.startsWith("jar:") || 
				(!allowSubprocess && !ALLOW_NATIVE_BGZF)){
			//Use plain gzip stream
		}else if(USE_UNBGZIP && ALLOW_NATIVE_BGZF && PREFER_NATIVE_BGZF_IN){
			return getUnbgzipStream(fname, helper);
		}else if(allowSubprocess){
			if(USE_UNBGZIP && (ALLOW_NATIVE_BGZF || Data.BGZIP())){return getUnbgzipStream(fname, helper);}
			if(USE_UNPIGZ && Data.PIGZ()){return getUnpigzStream(fname, helper);}
			if(USE_GUNZIP && Data.GUNZIP()){return getGunzipStream(fname);}
			//Fallthrough to plain
		}
		
		InputStream raw=getRawInputStream(fname, buffer, helper);
		InputStream in=null;
		try{
			in=new GZIPInputStream(raw, INBUF);
		}catch(Exception e){
			closeFailedInput(in, raw, e);
			System.err.println("Error when attempting to read "+fname);
			throw new RuntimeException(e);
		}catch(Error e){
			closeFailedInput(in, raw, e);
			throw e;
		}

		return in;
	}

	/** Attempts cleanup after Java input-wrapper setup fails, preserving the primary failure.
	 * An assigned wrapper is closed first, with raw closure as fallback. Direct borrowed
	 * System.in is never closed, including through its wrapper. This does not manage
	 * registered processes or promise immediate reclamation of an unassigned inflater.
	 * @param wrapper Constructed wrapper, or null when construction failed
	 * @param raw Backend returned by getRawInputStream
	 * @param primary Original setup failure; distinct cleanup failures are suppressed onto it */
	private static void closeFailedInput(final InputStream wrapper, final InputStream raw,
			final Throwable primary){
		if(raw==System.in){return;}
		Throwable wrapperFailure=null;
		if(wrapper!=null){
			try{
				wrapper.close();
				return;
			}catch(Throwable closing){
				wrapperFailure=closing;
				if(closing!=primary){primary.addSuppressed(closing);}
			}
		}
		try{raw.close();}
		catch(Throwable closing){
			if(closing!=primary && closing!=wrapperFailure){primary.addSuppressed(closing);}
		}
	}
	
	/** One whole-file helper acquisition; no shared or thread-local ownership state. */
	private static final class HelperInput{

		/** Preserves normal close-before-conversion ordering for every selected backend. */
		void closeNormally() throws IOException{
			assert(outer!=null) : "Whole-file helpers assign their outer wrapper before normal closure";
			closeAttempted=true;
			outer.close();
			closed=true;
			finishProcess();
		}

		/** Cleans only classified inputs and preserves the primary failure and its type. */
		void closeAfterFailure(final Throwable primary){
			if(!managed || closed){return;}
			try{
				Throwable wrapperFailure=null;
				if(!closeAttempted && outer!=null){
					closeAttempted=true;
					try{
						outer.close();
						closed=true;
						return;
					}catch(Throwable closing){
						wrapperFailure=closing;
						if(closing!=primary){primary.addSuppressed(closing);}
					}
				}
				if(stream==outer){return;}//The same object must not be closed twice after a failed attempt.
				try{stream.close(); closed=true;}
				catch(Throwable closing){
					if(closing!=primary && closing!=wrapperFailure){primary.addSuppressed(closing);}
				}
			}finally{
				try{finishProcess();}
				catch(Throwable completion){if(completion!=primary){primary.addSuppressed(completion);}}
			}
		}

		/** Records one successfully acquired external decoder; acquisition failures remain with the factory. */
		void acquiredProcess(final String name, final String codec){
			assert(processName==null) : "Each HelperInput owns one acquisition; replacing its registered process would lose completion: "+processName;
			managed=true;
			processName=name;
			processCodec=codec;
		}

		/** Collects the acquired decoder once, after its stdout has been closed.
		 * Uses the existing registry policy: exit0/141 accepted, unbounded wait.
		 * @throws IOException If the registered decoder reports failure */
		private void finishProcess() throws IOException{
			final String name=processName;
			if(name==null){return;}
			processName=null;
			assert(processCodec!=null && (name.endsWith(".bz2") || gzipSuffix(name) || name.endsWith(".zst"))) :
				"HelperInput tracks bzip2/pigz/bgzip/Zstandard dispatch; killProcess must recognize the suffix and the error must identify its codec: "+name;
			if(ReadWrite.killProcess(name)){throw new IOException(processCodec+" decoder failed for "+name);}
		}

		/** Successfully acquired backend, including any Java decoder. */
		InputStream stream;
		/** Most recently completed reader construction, for exceptional cleanup. */
		Closeable outer;
		/** Set at raw/native or successful bzip2/pigz/bgzip/Zstandard acquisition; other codecs leave it false. */
		boolean managed;
		/** Exact registered decoder filename; cleared before the one completion attempt. */
		String processName;
		/** Decoder name used in a completion error; preserves the existing Bzip2 wording. */
		String processCodec;
		/** Prevents retrying an outer close that already failed on the normal path. */
		boolean closeAttempted;
		/** Prevents duplicate closure if result conversion fails after a successful close. */
		boolean closed;
	}

	/** Lets helper-owned readers/decoders close without closing their borrowed raw input. */
	private static final class NonClosingInputStream extends FilterInputStream{
		NonClosingInputStream(final InputStream in){super(in);}
		@Override
		public int read(final byte[] b, final int off, final int len) throws IOException{
			return in.read(b, off, len);
		}
		@Override
		public void close(){}//System.in remains owned by the caller, below any decoder.
	}
	
	/**
	 * Creates gunzip input stream using external gzip command.
	 * @param fname Input filename
	 * @return Gunzip decompressed input stream
	 */
	public static InputStream getGunzipStream(String fname){
		if(verbose){System.err.println("getGunzipStream("+fname+")");}
		return getInputStreamFromProcess(fname, "gzip -c -d", false, true, true);
	}
	
	/**
	 * Creates unpigz input stream for parallel gzip decompression.
	 * @param fname Input filename
	 * @return Unpigz decompressed input stream
	 */
	public static InputStream getUnpigzStream(String fname){return getUnpigzStream(fname, null);}

	/** Records helper ownership only after the registered pigz stdout has been acquired. */
	private static InputStream getUnpigzStream(String fname, final HelperInput helper){
		if(verbose){System.err.println("getUnpigzStream("+fname+")");}
		final InputStream in=getInputStreamFromProcess(fname, "pigz -c -d", false, true, true);
		if(helper!=null){helper.acquiredProcess(fname, "pigz");}
		return in;
	}
	
	/**
	 * Creates unbgzip input stream with configured thread selection for
	 * block-gzip decompression.
	 * @param fname Input filename
	 * @return Unbgzip decompressed input stream
	 */
	public static InputStream getUnbgzipStream(String fname){return getUnbgzipStream(fname, null);}

	/** Threads helper ownership through native raw acquisition or successful external stdout acquisition.
	 * Native constructor failures remain separate.
	 * Decoder close behavior, including MT close-error suppression, is unchanged. */
	private static InputStream getUnbgzipStream(String fname, final HelperInput helper){
		if(verbose){System.err.println("getUnbgzipStream("+fname+")");}
		boolean stdin=FileFormat.isStdin(fname);
		boolean bgz=stdin ? true : isBGZip(fname);

		int threads=Tools.mid(BgzfSettings.READ_THREADS, 1, Shared.threads());
		if(!bgz && !stdin){threads=2;}
		
		if(nativeBgzfIn()){
//			System.err.println("Native BGZF");
			final InputStream in, raw=getRawInputStream(fname, true, helper);
			if(!BgzfSettings.USE_MULTITHREADED_BGZF){in=new BgzfInputStream(raw);}
			else{in=new BgzfInputStreamMT2(raw, threads);}
			return in;
		}
		Data.BGZIP();//Ensure that threads capability was detected
		final InputStream in=getInputStreamFromProcess(fname, "bgzip -c -d"+
			(Data.BGZIP_VERSION_threadsFlag ? " -@ "+threads : ""), false, true, true);
		if(helper!=null){helper.acquiredProcess(fname, "bgzip");}
		return in;
	}
	
	//Does not seem to work; just makes a big file somewhere (?).
	//Very slow, too.
	/**
	 * Creates strict Zstandard stdout; file redirection also accepts ordinary input symlinks.
	 * Historical warnings above describe the old named-file command; no speed claim is made.
	 * @param fname Input filename
	 * @return Caller-owned stream; finishReading with process collection reports decoder status
	 */
	public static InputStream getUnzstdStream(String fname){return getUnzstdStream(fname, null);}

	/** Records whole-file helper ownership after successful Zstandard stdout acquisition. */
	private static InputStream getUnzstdStream(final String fname, final HelperInput helper){
		if(verbose){System.err.println("getUnzstdStream("+fname+")");}
		final InputStream in=getInputStreamFromProcess(fname, "zstd -d -c", false, false, true);//FIO-123: no derived file or force-mode raw pass-through.
		if(helper!=null){helper.acquiredProcess(fname, "Zstandard");}
		return in;
	}
	
	/**
	 * Creates unpbzip2 input stream for parallel bzip2 decompression.
	 * @param fname Input filename
	 * @return Unpbzip2 decompressed input stream
	 */
	public static InputStream getUnpbzip2Stream(String fname){
		if(verbose){System.err.println("getUnpbzip2Stream("+fname+")");}
		return getInputStreamFromProcess(fname, "pbzip2 -c -d", false, true, true);
	}
	
	/**
	 * Creates unlbzip2 input stream for fast bzip2 decompression.
	 * @param fname Input filename
	 * @return Unlbzip2 decompressed input stream
	 */
	public static InputStream getUnlbzip2Stream(String fname){
		if(verbose){System.err.println("getUnlbzip2Stream("+fname+")");}
		return getInputStreamFromProcess(fname, "lbzip2 -c -d", false, true, true);
	}
	
	/**
	 * Creates unbzip2 input stream using external bzip2 command.
	 * @param fname Input filename
	 * @return Unbzip2 decompressed input stream
	 */
	public static InputStream getUnbzip2Stream(String fname){
		if(verbose){System.err.println("getUnbzip2Stream("+fname+")");}
		return getInputStreamFromProcess(fname, "bzip2 -c -d", false, true, true);
	}
	
	/**
	 * Creates DSRC decompression input stream with configured thread selection.
	 * @param fname Input filename
	 * @return DSRC decompressed input stream
	 */
	public static InputStream getUnDsrcStream(String fname){
		if(verbose){System.err.println("getUnDsrcStream("+fname+")");}
		int threads=Tools.min(MAX_ZIP_THREADS, Tools.max((int)((Shared.threads()+1)*ZIP_THREAD_MULT), 1));
		threads=Tools.max(1, Tools.min(Shared.threads()-1, threads));
		return getInputStreamFromProcess(fname, "dsrc d -s -t"+threads, false, true, true);
	}
	
	
	/** Starts a named input process and returns its stdout pipe. ProcessBuilder splits
	 * command on literal spaces; stdin aliases inherit stdin, otherwise appendFname adds
	 * the path or redirects process input from it. cat is used only in the legacy Runtime
	 * branch, whose cat mode is assertion-disabled/untested. The process is registered
	 * separately from its stream; callers must perform applicable killProcess completion. */
	public static InputStream getInputStreamFromProcess(final String fname, String command, boolean cat, final boolean appendFname, final boolean useProcessBuilder){
		if(verbose){System.err.println("getInputStreamFromProcess("+fname+", "+command+", "+cat+")");}

		//InputStream raw=getRawInputStream(fname, false);
		InputStream in=null;
		
		Process p=null;
		
		if(useProcessBuilder){
			ProcessBuilder pb=new ProcessBuilder();
			pb.redirectError(Redirect.INHERIT);
			
			if(fname.equals("stdin") || fname.startsWith("stdin.")){
				pb.redirectInput(Redirect.INHERIT);
				pb.command(command.split(" "));
			}else{
				if(appendFname){
					command=command+" "+fname;
				}else{
					pb.redirectInput(new File(fname));
				}
//				System.err.println(command+", "+appendFname);
				pb.command(command.split(" "));
			}
			try{
				p=pb.start();
			}catch(IOException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			assert(p!=null) : "Could not execute "+command;
			addProcess(fname, p);
			in=p.getInputStream();
//			{
//				out=p.getOutputStream();
//				InputStream es=p.getErrorStream();
//				assert(es!=null);
//				PipeThread et=new PipeThread(es, System.err);
//				addPipeThread(fname, et);
//				et.start();
//			}
			return in;
		}
		
		if(!appendFname){
			try{
				p=Runtime.getRuntime().exec(command);
			}catch(IOException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}else if(fname.equals("stdin") || fname.startsWith("stdin.")){
			try{
				if(cat){
					throw new RuntimeException();
				}else{
					p=Runtime.getRuntime().exec(command);
				}
			}catch(IOException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			assert(p!=null) : "Could not execute "+command;
			OutputStream os=p.getOutputStream();
			PipeThread it=new PipeThread(System.in, os);
			addPipeThread(fname, it);
			it.start();
		}else{
			try{
				if(cat){
					assert(false) : "This mode is untested.";
					String[] cmd={
							"sh","cat "+fname,
							" | "+command
					};
					p=Runtime.getRuntime().exec(cmd);
				}else{
					p=Runtime.getRuntime().exec(command+" "+fname);
				}
			}catch(IOException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		assert(p!=null) : "Could not execute "+command;
		
		addProcess(fname, p);
		in=p.getInputStream();
		InputStream es=p.getErrorStream();
		assert(es!=null);
		PipeThread et=new PipeThread(es, System.err);
		addPipeThread(fname, et);
		et.start();
		
		return in;
	}
	
	
	/**
	 * Creates bzip2 input stream with error handling wrapper.
	 * Tries available bzip2 decompression tools in preference order.
	 *
	 * @param fname Input filename
	 * @param allowSubprocess Unused; available external tools are selected regardless
	 * @return Bzip2 decompressed input stream
	 */
	public static InputStream getBZipInputStream(String fname, boolean allowSubprocess){return getBZipInputStream(fname, allowSubprocess, null);}

	/** Records helper ownership only after a registered bzip2 stdout stream is returned. */
	private static InputStream getBZipInputStream(String fname, boolean allowSubprocess, final HelperInput helper){
		if(verbose){System.err.println("getBZipInputStream("+fname+")");}
		InputStream in=null;
		
		try{in=getBZipInputStream2(fname, allowSubprocess);}
		catch(IOException e){
			System.err.println("Error when attempting to read "+fname);
			throw new RuntimeException(e);
		}catch(NullPointerException e){
			System.err.println("Error when attempting to read "+fname);
			throw new RuntimeException(e);
		}
		
		assert(in!=null);
		if(helper!=null){helper.acquiredProcess(fname, "Bzip2");}
		return in;
	}
	
	/**
	 * Creates bzip2 input stream using external tools in preference order:
	 * lbzip2, pbzip2, or bzip2.
	 *
	 * @param fname Input filename
	 * @param allowSubprocess Unused; available external tools are selected regardless
	 * @return Bzip2 decompressed input stream
	 * @throws IOException if no bzip2 tools are available
	 */
	private static InputStream getBZipInputStream2(String fname, boolean allowSubprocess) throws IOException{
		if(verbose){
			if(verbose){System.err.println("getBZipInputStream("+fname+")");}
		}
		
		if(!fname.startsWith("jar:")){
			if(verbose){System.err.println("Fetching bz2 input stream: "+fname+", "+USE_PBZIP2+", "+USE_BZIP2+", "+Data.PBZIP2()+Data.BZIP2());}
			if(USE_LBZIP2 && Data.LBZIP2()){return getUnlbzip2Stream(fname);}
			if(USE_PBZIP2 && Data.PBZIP2()){return getUnpbzip2Stream(fname);}
			if(USE_BZIP2 && Data.BZIP2()){return getUnbzip2Stream(fname);}
		}
		
		throw new IOException("\nlbzip2, pbzip2, or bzip2 must be in the path to read bz2 files:\n"+fname+"\n");
	}
	
	/**
	 * Creates DSRC input stream with error handling wrapper.
	 * @param fname Input filename
	 * @return DSRC decompressed input stream
	 */
	public static InputStream getDsrcInputStream(String fname){
		if(verbose){System.err.println("getDsrcInputStream("+fname+")");}
		InputStream in=null;
		
		try{in=getDsrcInputStream2(fname);}
		catch(IOException e){
			System.err.println("Error when attempting to read "+fname);
			throw new RuntimeException(e);
		}catch(NullPointerException e){
			System.err.println("Error when attempting to read "+fname);
			throw new RuntimeException(e);
		}
		
		assert(in!=null);
		return in;
	}
	
	/**
	 * Creates DSRC input stream using external DSRC tool.
	 * @param fname Input filename
	 * @return DSRC decompressed input stream
	 * @throws IOException if DSRC tool not available
	 */
	private static InputStream getDsrcInputStream2(String fname) throws IOException{
		if(verbose){
			if(verbose){System.err.println("getDsrcInputStream2("+fname+")");}
		}
		
		if(USE_DSRC && Data.DSRC()){return getUnDsrcStream(fname);}
		
		throw new IOException("\nDsrc must be in the path to read Dsrc files:\n"+fname+"\n");
	}
	
	/**
	 * Creates XZ input stream (currently disabled).
	 * @param fname Input filename
	 * @return null (XZ support disabled)
	 */
	public static InputStream getXZInputStream(String fname){
		
		InputStream in=null;
		
//		if(PROCESS_XZ){
//			InputStream raw=getRawInputStream(fname, true);
//			try {
//				in=new org.tukaani.xz.XZInputStream(raw);
//			} catch (FileNotFoundException e) {
//				throw new RuntimeException(e);
//			} catch (IOException e) {
//				throw new RuntimeException(e);
//			}
//		}

		return in;
	}

	/**
	 * Collects raw bytes until a read returns zero or EOF; no decompression is selected.
	 * Attempts closure on normal return or a propagated failure, preserving borrowed
	 * stdin. The result must fit a single byte array; VM termination cannot be cleaned up.
	 * @param fname Input filename
	 * @return File bytes through the first zero read or EOF
	 * @throws IOException if reading or closing the acquired stream fails
	 */
	public static byte[] readRaw(String fname) throws IOException{
		final HelperInput input=new HelperInput();
		//Failed acquisition remains the factory's responsibility.
		input.stream=getRawInputStream(fname, false, input);
		input.outer=input.stream;
		try{
			ByteBuilder bb=new ByteBuilder();
			byte[] buffer=new byte[16384];
			int x=input.stream.read(buffer);
			while(x>0){
				bb.append(buffer, x);
				x=input.stream.read(buffer);
			}
			input.closeNormally();
			return bb.toBytes();
		}catch(IOException e){
			input.closeAfterFailure(e);
			throw e;
		}catch(RuntimeException | Error e){
			input.closeAfterFailure(e);
			throw e;
		}
	}
	
	/**
	 * Reads an object and performs an unchecked generic cast; cx is not consulted.
	 *
	 * @param <X> Expected object type
	 * @param cx Unused type token
	 * @param fname Input filename
	 * @param allowSubprocess Whether to allow subprocess decompression
	 * @return Deserialized object of specified type
	 */
	public static <X> X read(Class<X> cx, String fname, boolean allowSubprocess){
		X x=(X)readObject(fname, allowSubprocess);
		return x;
	}
	
	/**
	 * Reads an object and casts to an object array; element class cx is not consulted.
	 *
	 * @param <X> Expected array element type
	 * @param cx Unused element type token
	 * @param fname Input filename
	 * @param allowSubprocess Whether to allow subprocess decompression
	 * @return Deserialized array of specified type
	 */
	public static <X> X[] readArray(Class<X> cx, String fname, boolean allowSubprocess){
		X[] x=(X[])readObject(fname, allowSubprocess);
		return x;
	}
	
	/**
	 * Reads an object and casts to a two-dimensional object array; cx is not consulted.
	 *
	 * @param <X> Expected array element type
	 * @param cx Unused element type token
	 * @param fname Input filename
	 * @param allowSubprocess Whether to allow subprocess decompression
	 * @return Deserialized 2D array of specified type
	 */
	public static <X> X[][] readArray2(Class<X> cx, String fname, boolean allowSubprocess){
		X[][] x=(X[][])readObject(fname, allowSubprocess);
		return x;
	}
	
	/**
	 * Reads an object and casts to a three-dimensional object array; cx is not consulted.
	 *
	 * @param <X> Expected array element type
	 * @param cx Unused element type token
	 * @param fname Input filename
	 * @param allowSubprocess Whether to allow subprocess decompression
	 * @return Deserialized 3D array of specified type
	 */
	public static <X> X[][][] readArray3(Class<X> cx, String fname, boolean allowSubprocess){
		X[][][] x=(X[][][])readObject(fname, allowSubprocess);
		return x;
	}
	
	
	/**
	 * Normalizes backslashes, removes directories and strips one recognized suffix.
	 * Recognizes lowercase .gz, .zip, .dsrc and, when enabled, .bz2. The xz local is
	 * unused: .xz and other compression aliases remain in the returned name.
	 * @param fname Full file path
	 * @return Last path component after the supported single-suffix removal
	 */
	public static String basename(String fname){
		fname=fname.replace('\\', '/');
		boolean xz=fname.endsWith(".xz");
		boolean gzipped=fname.endsWith(".gz");
		boolean zipped=fname.endsWith(".zip");
		boolean bzipped=PROCESS_BZ2 && fname.endsWith(".bz2");
		boolean dsrced=fname.endsWith(".dsrc");
		String basename=fname;
//		if(basename.contains("\\")){basename=basename.substring(basename.lastIndexOf("\\")+1);}
		if(basename.contains("/")){basename=basename.substring(basename.lastIndexOf('/')+1);}
		if(zipped || bzipped){basename=basename.substring(0, basename.length()-4);}
		else if(gzipped){basename=basename.substring(0, basename.length()-3);}
		else if(dsrced){basename=basename.substring(0, basename.length()-5);}
		return basename;
	}
	
	/**
	 * Visits compressedExtensions once in table order, stripping repeated matches for
	 * each suffix. Matching is case-sensitive; mixed suffix layers can remain when
	 * stripping a later table entry exposes an earlier one. Does not remove directories.
	 * @param fname Filename with potential compression extensions
	 * @return Name after this ordered suffix pass
	 */
	public static String rawName(String fname){
		for(String s : compressedExtensions){
			while(fname.endsWith(s)){fname=fname.substring(0, fname.length()-s.length());}
		}
		return fname;
	}
	
	/** 
	 * Returns the path without the file extension.
	 * Only strips known extensions. */
	public static String stripExtension(String fname){
		if(fname==null){return null;}
		for(String ext : FileFormat.EXTENSION_LIST){
			String s="."+ext;
			if(fname.endsWith(s)){return stripExtension(fname.substring(0, fname.length()-s.length()));}
		}
		return fname;
	}
	
	/** Returns the whole extension, include compression and raw type */
	public static String getExtension(String fname){
		if(fname==null){return null;}
		String stripped=stripExtension(fname);
		if(stripped==null){return fname;}
		if(stripped.length()==fname.length()){return "";}
		return fname.substring(stripped.length());
	}
	
	/**
	 * Strips both path and extension, returning just the core filename.
	 * @param fname Full file path
	 * @return Core filename without path or extension
	 */
	public static String stripToCore(String fname){
		fname=stripPath(fname);
		return stripExtension(fname);
	}
	
	/**
	 * Strips the directories, leaving only a filename
	 * @param fname Path, or null
	 * @return File name without directories
	 */
	public static String stripPath(String fname){
		if(fname==null){return null;}
		fname=fname.replace('\\', '/');
		int idx=fname.lastIndexOf('/');
		if(idx>=0){fname=fname.substring(idx+1);}
		return fname;
	}
	
	/**
	 * Extracts directory path from full filename.
	 * @param fname Full file path
	 * @return Normalized directory prefix including separator, empty string if none, or null for null input
	 */
	public static String getPath(String fname){
		if(fname==null){return null;}
		fname=fname.replace('\\', '/');
		int idx=fname.lastIndexOf('/');
		if(idx>=0){return fname.substring(0, idx+1);}
		return "";
	}
	
	/**
	 * Determines compression type from filename extension.
	 * @param fname Filename to analyze
	 * @return Compression type string or null if not compressed
	 */
	public static String compressionType(String fname){
		if(fname==null || fname.indexOf('.')<0){return null;}
		fname=fname.toLowerCase(Locale.ENGLISH);
		for(int i=0; i<compressedExtensions.length; i++){
			if(fname.endsWith(compressedExtensions[i])){return compressedExtensionMap[i];}
		}
		return null;
	}
	
	/**
	 * Checks if filename indicates a compressed file format.
	 * @param fname Filename to check
	 * @return true if filename has compression extension
	 */
	public static boolean isCompressed(String fname){
		return compressionType(fname)!=null;
	}
	
	/**
	 * Checks if filename indicates SAM format, including compressed SAM.
	 * @param fname Filename to check
	 * @return true if filename indicates SAM format
	 */
	public static boolean isSam(String fname){
		fname=fname.toLowerCase(Locale.ENGLISH);
		if(fname.endsWith(".sam")){return true;}
		String s=compressionType(fname);
		if(s==null){return false;}
		return fname.substring(0, fname.lastIndexOf('.')).endsWith(".sam");
	}
	
	/** Returns text after the last dot following rawName, lowercased with Locale.ENGLISH.
	 * A dotless name returns the whole remaining name; directory dots are not excluded. */
	public static String rawExtension(String fname){
		fname=rawName(fname);
		int x=fname.lastIndexOf('.');
		//if(x<0){return "";}
		return fname.substring(x+1).toLowerCase(Locale.ENGLISH);
	}
	
	/**
	 * Extracts root directory from file path, ensuring trailing separator.
	 * Throws exception if path doesn't exist.
	 *
	 * @param path File or directory path
	 * @return Root directory path with trailing separator
	 * @throws RuntimeException if path not found
	 */
	public static String parseRoot(String path){
		File f=new File(path);
		if(f.isDirectory()){
			if(!path.endsWith(FILESEP)){
				path=path+FILESEP;
			}
			return path;
		}else if(f.isFile()){
			int slash=path.lastIndexOf(FILESEP);
			if(slash<0){
				return "";
			}else{
				return path.substring(0, slash+1);
			}
		}else{
			throw new RuntimeException("Can't find "+path); //Try using parseRoot2 instead.
		}
	}
	
	/** This one does not throw an exception for non-existing paths */
	public static String parseRoot2(String path){
		File f=new File(path);
		
		if(!f.exists()){
			if(path.endsWith(FILESEP)){return path;}
			int slash=path.lastIndexOf(FILESEP);
			if(slash<0){
				return "";
			}else{
				return path.substring(0, slash+1);
			}
		}
		
		if(f.isDirectory()){
			if(!path.endsWith(FILESEP)){
				path=path+FILESEP;
			}
			return path;
		}else if(f.isFile()){
			int slash=path.lastIndexOf(FILESEP);
			if(slash<0){
				return "";
			}else{
				return path.substring(0, slash+1);
			}
		}else{
			throw new RuntimeException("Can't find "+path);
		}
	}
	
	/**
	 * Finds existing file by trying various compression extensions.
	 * Searches for uncompressed, gzip, zip, bzip2, and xz variants.
	 * @param fname Base filename to search for
	 * @return Path to existing file or original name if none found
	 */
	public static String findFileExtension(final String fname){

		File file=new File(fname);
		if(file.exists()){return fname;}

		String basename=fname, temp;
		if(fname.endsWith(".zip") || fname.endsWith(".gz") || (PROCESS_BZ2 && fname.endsWith(".bz2")) || (PROCESS_XZ && fname.endsWith(".xz"))){
			basename=fname.substring(0, fname.lastIndexOf('.'));
		}
		temp=basename;
		file=new File(temp);
		if(!file.exists()){
			temp=basename+".gz";
			file=new File(temp);
		}
//		System.err.println(temp+" "+(file.exists() ? " exists" : " does not exist"));
		if(!file.exists()){
			temp=basename+".zip";
			file=new File(temp);
		}
//		System.err.println(temp+" "+(file.exists() ? " exists" : " does not exist"));
		if(!file.exists() && PROCESS_BZ2){
			temp=basename+".bz2";
			file=new File(temp);
		}
//		System.err.println(temp+" "+(file.exists() ? " exists" : " does not exist"));
		if(!file.exists() && PROCESS_XZ){
			temp=basename+".xz";
			file=new File(temp);
		}
//		System.err.println(temp+" "+(file.exists() ? " exists" : " does not exist"));
		if(!file.exists()){temp=fname;}
		
		return temp;
	}
	
	/**
	 * Attempts nonrecursive deletion of an existing path and reports success.
	 * A nonempty directory remains intact and returns false.
	 * @param path Path, or null
	 * @param verbose Whether to print the attempted path
	 * @return True only when File.delete() succeeds; false for null/missing paths,
	 * failed deletion, or a caught Exception
	 */
	public static boolean delete(String path, boolean verbose){
		if(path==null){return false;}
		if(verbose){System.err.println("Trying to delete "+path);}
		File f=new File(path);
		if(f.exists()){
			try{
				return f.delete();
			}catch(Exception e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return false;
	}
	
	/**
	 * Copies file from source to destination without creating intermediate paths
	 */
	public static synchronized void copyFile(String source, String dest){copyFile(source, dest, false);}
	/**
	 * Copies file from source to destination with optional path creation.
	 * Uses global RAWMODE for same-suffix ZIP/GZIP/BZ2/XZ copies; otherwise factories
	 * may transform compression. Destination nonexistence is assertion-checked.
	 * The class lock serializes copies, not all I/O. Normal and caught-IOException
	 * paths restore RAWMODE; unchecked failures can leave it changed. Exceptional
	 * stream cleanup is not guaranteed. The global-state repair remains deferred (#004).
	 *
	 * @param source Source file path
	 * @param dest Destination file path (must not exist)
	 * @param createPathIfNeeded Whether to create parent directories
	 */
	public static synchronized void copyFile(String source, String dest, boolean createPathIfNeeded){
		
		assert(!new File(dest).exists()) : "Destination file already exists: "+dest;
		if(createPathIfNeeded){
			File parent=new File(dest).getParentFile();
			if(parent!=null && !parent.exists()){
				parent.mkdirs();
			}
		}
		
		//TODO [fileIO/ReadWrite#004] copyFile toggles the GLOBAL static RAWMODE during a synchronized copy; a concurrent
		//compressed-I/O thread would transiently see it flipped. DEFERRED per Brian (never bitten; a clean fix is a small
		//refactor — thread-local rawmode or an explicit raw param — that needs its own testing). Used only for byte-copying
		//already-compressed files (same extension in+out).
		final boolean oldRawmode=RAWMODE;
		if((source.endsWith(".zip") && dest.endsWith(".zip"))
				 || (source.endsWith(".gz") && dest.endsWith(".gz")
						 || (source.endsWith(".bz2") && dest.endsWith(".bz2"))
						 || (source.endsWith(".xz") && dest.endsWith(".xz")))){
			RAWMODE=true;
		}
		
		try{
			InputStream in=getInputStream(source, false, false);
			OutputStream out=getOutputStream(dest, false, false, true);

			byte[] buffer=new byte[INBUF];
			int len;
			
			while((len=in.read(buffer))>0){
				out.write(buffer, 0, len);
			}
			
			in.close();
			out.flush();
			if(out.getClass()==ZipOutputStream.class){
				ZipOutputStream zos=(ZipOutputStream)out;
				zos.closeEntry();
				zos.finish();
			}
//			else if(PROCESS_XZ && out.getClass()==org.tukaani.xz.XZOutputStream.class){
//				org.tukaani.xz.XZOutputStream zos=(org.tukaani.xz.XZOutputStream)out;
//				zos.finish();
//			}
			out.close();
			
		}catch(FileNotFoundException e){
			RAWMODE=oldRawmode;
			throw new RuntimeException(e);
		}catch(IOException e){
			RAWMODE=oldRawmode;
			throw new RuntimeException(e);
		}
		
		RAWMODE=oldRawmode;
	}
	
	/**
	 * Recursively copies all contents from source directory to destination.
	 * Creates destination directories as needed and handles nested structures.
	 * @param from Source directory path
	 * @param to Destination directory path
	 */
	public static void copyDirectoryContents(String from, String to){
		assert(!from.equalsIgnoreCase(to));
		
		if(to.indexOf('\\')>0){to=to.replace('\\', '/');}
		
		File d1=new File(from);
		assert(d1.exists());
		assert(d1.isDirectory());
		
		File d2=new File(to);
		assert(!d1.equals(d2));
		if(d2.exists()){
			assert(d2.isDirectory());
		}else{
			d2.mkdirs();
		}
		if(!to.endsWith("/")){to=to+"/";}
		
		File[] array=d1.listFiles();
		
		for(File f : array){
			String name=f.getName();
			String dest=to+name;
			if(f.isFile()){
				copyFile(f.getAbsolutePath(), dest);
			}else{
				assert(f.isDirectory());
				File f2=new File(dest);
				if(!f2.exists()){
					f2.mkdir();
				}else{
					assert(f2.isDirectory());
				}
				copyDirectoryContents(f.getAbsolutePath(), f2.getAbsolutePath());
			}
		}
		
	}
	
	
	/** Reserves positive active/waiting counts, or delegates negative completion to
	 * addRunningThread. Uses the activeThreads array as both mutable counters and monitor. */
	static final int addThread(int x){
		if(verbose){System.err.println("addThread("+x+")");}
		synchronized(activeThreads){
			assert(x!=0);
			if(x>0){
				activeThreads[0]+=x;
				activeThreads[1]+=x;
			}else{
				addRunningThread(x);
			}
			assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 &&
					activeThreads[2]>=0 && activeThreads[2]<=maxWriteThreads) : Arrays.toString(activeThreads);
					
			return activeThreads[0];
		}
	}
	
	/** Moves positive counts from waiting to running after capacity becomes available;
	 * negative counts leave both running and active. LOW_MEMORY caps running at one.
	 * Interrupted waits print and retry without restoring interruption. */
	static final int addRunningThread(int x){
		if(verbose){System.err.println("addRunningThread("+x+")");}
		final int max=(Shared.LOW_MEMORY ? 1 : maxWriteThreads);
		synchronized(activeThreads){
			assert(x!=0);
			if(x>0){
				assert(activeThreads[1]>=x);
				while(activeThreads[2]>=max){
					try{
						activeThreads.wait();
					}catch(InterruptedException e){
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				activeThreads[1]-=x; //Remove from waiting
			}else{
				activeThreads[0]+=x; //Remove from active
			}
			activeThreads[2]+=x; //Change number running
			
			assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 &&
					activeThreads[2]>=0 && activeThreads[2]<=max) : Arrays.toString(activeThreads);
			
			if(activeThreads[2]==0 || (activeThreads[2]<max && activeThreads[1]>0)){activeThreads.notifyAll();}//FIXED [fileIO/ReadWrite#007]: notify()->notifyAll(); mixed worker(no-timeout)+main(waitForWritingToFinish,8s) waiters, a wasted single-notify to the main thread stalled a ready worker until the 8s timeout. Latency fix (self-healed, output was correct). Twin of LoadThread#002.
			return activeThreads[2];
		}
	}
	
	/** Returns current count of active write threads.
	 * @return Number of active threads (running + waiting) */
	public static final int countActiveThreads(){
		if(verbose){System.err.println("countActiveThreads()");}
		synchronized(activeThreads){
			assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 &&
					activeThreads[2]>=0 && activeThreads[2]<=maxWriteThreads) : Arrays.toString(activeThreads);
			return activeThreads[0];
		}
	}
	
	/** Waits for the counted InThread operations to decrement the active count to zero.
	 * This is not a worker join or failure report. Interrupted waits print and retry.
	 * Worker finally blocks release their counts on write failure when the JVM continues;
	 * assertions-enabled write failures instead halt the JVM under IO-042. */
	public static final void waitForWritingToFinish(){
		if(verbose){System.err.println("waitForWritingToFinish()");}
		synchronized(activeThreads){
			while(activeThreads[0]>0){
				assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 &&
						activeThreads[2]>=0 && activeThreads[2]<=maxWriteThreads) : Arrays.toString(activeThreads);
				try{
					activeThreads.wait(8000);
				}catch(InterruptedException e){
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(activeThreads[2]==0 || (activeThreads[2]<maxWriteThreads && activeThreads[1]>0)){activeThreads.notifyAll();}//FIXED [fileIO/ReadWrite#007]: notify()->notifyAll(), see addRunningThread.
			}
		}
	}


	/** Closes the input if present and returns its recorded error state; no explicit join here. */
	public static final boolean closeStream(Streamer st){return closeStreams(st, (Writer[])null);}
	/** Calls poisonAndWait on each nonnull Writer and combines failure states. */
	public static final boolean closeOutputStreams(Writer...ross){return closeStreams(null, ross);}

	
	/** Closes single concurrent read stream and reports error state */
	public static final boolean closeStream(ConcurrentReadStreamInterface cris){return closeStreams(cris, (ConcurrentReadOutputStream[])null);}
	/** Closes single concurrent read output stream and reports error state */
	public static final boolean closeStream(ConcurrentReadOutputStream ross){return closeStreams((ConcurrentReadStreamInterface)null, ross);}
	/** Closes multiple concurrent read output streams and reports error state */
	public static final boolean closeOutputStreams(ConcurrentReadOutputStream...ross){return closeStreams(null, ross);}

	/**
	 * Closes all streams in MultiCros collection.
	 * @param mc MultiCros containing streams to close
	 * @return true if any errors occurred during closure
	 */
	public static final boolean closeStreams(MultiCros mc){
		if(mc==null){return false;}
		return closeStreams(null, mc.streamList.toArray(new ConcurrentReadOutputStream[0]));
	}
	
	/** 
	 * Closes the input and reads its error state, then closes and joins each output.
	 * Null entries are skipped. Unchecked exceptions propagate and can skip later streams.
	 * @param cris An input stream.  May be null.
	 * @param ross Zero or more output streams.
	 * @return True if an error was encountered.
	 */
	public static final boolean closeStreams(ConcurrentReadStreamInterface cris, ConcurrentReadOutputStream...ross){
		if(verbose){
			System.err.println("closeStreams("+cris+", "+(ross==null ? "null" : ross.length)+")");
			new Exception().printStackTrace(System.err);
		}
		boolean errorState=false;
		if(cris!=null){
			if(verbose){System.err.println("Closing cris; error="+errorState+"; c.error="+cris.errorState());}
			cris.close();
			errorState|=cris.errorState();
//			Object[] prods=cris.producers();
//			for(Object o : prods){
//				if(o!=null && o.getClass()==ReadInputStream.class){
//					ReadInputStream ris=(ReadInputStream)o;
//					ris.
//				}
//			}
			if(verbose){System.err.println("Closed cris; error="+errorState);}
		}
		if(ross!=null){
			for(ConcurrentReadOutputStream ros : ross){
				if(ros!=null){
					if(verbose){System.err.println("Closing ros "+ros+"; error="+errorState);}
					ros.close();
					ros.join();
					errorState|=(ros.errorState() || !ros.finishedSuccessfully());
					if(verbose){System.err.println("Closed ros; error="+errorState);}
				}
			}
		}
		return errorState;
	}
	
	/** 
	 * Closes the input and reads its error state, then calls poisonAndWait on each writer.
	 * No separate input join occurs here. Null entries are skipped; unchecked exceptions
	 * propagate and can skip later writers.
	 * @param st An input stream.  May be null.
	 * @param writers Zero or more output streams.
	 * @return True if an error was encountered.
	 */
	public static final boolean closeStreams(Streamer st, Writer...writers){
		if(verbose){
			System.err.println("closeStreams("+st+", "+(writers==null ? "null" : writers.length)+")");
			new Exception().printStackTrace(System.err);
		}
		boolean errorState=false;
		if(st!=null){
			if(verbose){System.err.println("Closing streamer "+st.fname()+
				"; error="+errorState+"; c.error="+st.errorState());}
			st.close();
			errorState|=st.errorState();
			if(verbose){System.err.println("Closed streamer"+st.fname()+"; error="+errorState);}
		}
		if(writers!=null){
			for(Writer w : writers){
				if(w!=null){
					if(verbose){System.err.println("Closing writer "+w.fname()+"; error="+errorState);}
					w.poisonAndWait();
					boolean a=w.errorState(), b=w.finishedSuccessfully();
					errorState|=(a || !b);
					if(verbose){System.err.println("Closed writer "+w.fname()+"; error="+errorState
						+", a="+a+", b="+b+", "+w);}
				}
			}
		}
		return errorState;
	}
	
	/**
	 * Removes the registered process, waits for it and terminates registered pipe threads.
	 * The name must indicate recognized compression or .bam unless FORCE_KILL is set.
	 * waitFor holds the process-map lock and has no timeout; only twenty interruptions
	 * trigger destroy. Exit 0 or 141 is accepted. Callers must close process input first
	 * where required; this method is not an unconditional immediate kill.
	 * @param fname Filename used to identify associated process
	 * @return true if errors occurred during process termination
	 */
	public static boolean killProcess(String fname){
		if(verbose){
			System.err.println("killProcess("+fname+")");
			new Exception().printStackTrace(System.err);
			System.err.println("processMap before: "+processMap.keySet());
		}
		if(fname==null || (!isCompressed(fname) && !fname.endsWith(".bam") && !FORCE_KILL)){return false;}
		
		boolean error=false;
		synchronized(processMap){
			Process p=processMap.remove(fname);
			if(p!=null){
				if(verbose){System.err.println("Found Process for "+fname);}
				int x=-1, tries=0;
				for(; tries<20; tries++){
					if(verbose){System.err.println("Trying p.waitFor()");}
					try{
//						long t=System.nanoTime();
//						Thread.sleep(4000);
						if(verbose){System.err.println("p.isAlive()="+p.isAlive());}
						x=p.waitFor();
//						if(verbose){System.err.println(System.nanoTime()-t+" ns");}
						if(verbose){System.err.println("success; return="+x);}
						break;
					}catch(InterruptedException e){
						if(verbose){System.err.println("Failed.");}
						e.printStackTrace();
					}
				}
				error|=(tries>=20 || (x!=0 && x!=141));//141 is sigpipe and appears to be OK when forcibly closing a pipe.
				if(verbose){System.err.println("killProcess("+fname+") returned "+error+"; tries="+tries+", code="+x);}
				if(tries>=20){
					if(verbose){System.err.println("Calling p.destroy because tries=="+tries+"; error="+error);}
					p.destroy();
					if(verbose){System.err.println("destroyed");}
				}
			}else{
				if(verbose){System.err.println("WARNING: Could not find process for "+fname);}
			}
			if(verbose){
				System.err.println("processMap after: "+processMap.keySet());
			}
		}
		synchronized(pipeThreadMap){
			if(verbose){System.err.println("pipeMap before: "+processMap.keySet());}
			ArrayList<PipeThread> atp=pipeThreadMap.remove(fname);
			if(atp!=null){
				for(PipeThread p : atp){
					if(p!=null){
						if(verbose){System.err.println("Found PipeThread for "+fname);}
						p.terminate();
						if(verbose){System.err.println("Terminated PipeThread");}
					}else{
						if(verbose){System.err.println("WARNING: Could not find process for "+fname);}
					}
				}
			}
			if(verbose){System.err.println("pipeMap after: "+processMap.keySet());}
		}
		if(verbose){System.err.println("killProcess("+fname+") returned "+error);}
		return error;
	}
	
	/** Registers one process per exact filename under the map lock. Replacing an entry
	 * destroys the previous process and invokes KillSwitch; it is not normal replacement. */
	private static void addProcess(String fname, Process p){
		if(verbose){
			System.err.println("addProcess("+fname+", "+p+")");
			new Exception().printStackTrace();
		}
		synchronized(processMap){
			Process old=processMap.put(fname, p);
			if(old!=null){
				old.destroy();
//				throw new RuntimeException("Duplicate process for file "+fname);
				KillSwitch.kill("Duplicate process for file "+fname);
			}
		}
	}
	
	/** Adds a pipe thread to the list for an exact filename under the pipe-map lock. */
	private static void addPipeThread(String fname, PipeThread pt){
		if(verbose){System.err.println("addPipeThread("+fname+", "+pt+")");}
		synchronized(pipeThreadMap){
//			System.err.println("Adding PipeThread for "+fname);
			ArrayList<PipeThread> atp=pipeThreadMap.get(fname);
			if(atp==null){
				atp=new ArrayList<PipeThread>(2);
				pipeThreadMap.put(fname, atp);
			}
			atp.add(pt);
		}
	}
	
	/** Returns false for case-sensitive stdin-prefixed names; otherwise compares the
	 * raw four-byte signature with 529205252. This is a heuristic, not complete BGZF validation. */
	public static boolean isBGZip(String fname){
		if(fname.startsWith("stdin")){return false;}
		int magic=getMagicNumber(fname);
		return magic==529205252;
//		byte[] header=getFirst16Bytes(fname);
//		if(header[10]==0x42 && header[11]==0x43){ //BC subfield
//			return true; //Is BGZF
//		}
	}
	
	/** Legacy private helper that reads n raw bytes and closes on the normal path.
	 * Returns null for a short read, except the empty-file negative-offset defect below.
	 * Caught Exceptions invoke KillSwitch; this method has no current call sites here. */
	private static byte[] getFirstNBytes(String fname, int n){
		byte[] array=new byte[n];
		FileInputStream fis;
		try{
			fis=new FileInputStream(fname);
			//TODO: Probable bug - FIO-047: immediate EOF leaves read=-1, then the next read uses a negative offset.
			int read=fis.read(array, 0, n);
			while(read<n){
				int r=fis.read(array, read, n-read);
				read+=r;
				if(r<1){break;}
			}
			fis.close();
			return read<n ? null : array;
		}catch(Exception e){
			KillSwitch.exceptionKill(e);
			return null;
		}
	}

	/** Reads four raw bytes in big-endian order and attempts closure.
	 * Missing files print an exception and return zero. EOF bytes become 255; caught
	 * read/close IOExceptions are printed, so this is not strict header validation.
	 * Magic number commonly used for bgzip files (first 4 bytes):
	 * 1f 8b 08 04
	 * 31 139 8 4
	 *  = 529205252
	 * 
	 * gzip/pigz:
	 * 1f 8b 08 00
	 * 31 139 8 0
	 *  = 529205248
	 * 
	 * od --format=x1 --read-bytes=16 names.txt_gzip.gz
	 */
	//[fileIO/ReadWrite#001 FIXED 2026-06-15]: 'is' was opened but never closed -> file-descriptor leak on
	//every call, on the common gzip-read path (getUnbgzipStream->isBGZip->getMagicNumber) and FileFormat.bgzip().
	//Now closes 'is' after reading and guards the null case (missing file) instead of NPEing at is.read().
	//Mirrors getFirstNBytes() above, which already closes its fis.
	public static int getMagicNumber(String fname){
		InputStream is=null;
		try{
			FileInputStream fis=new FileInputStream(fname);
			is=new BufferedInputStream(fis);
		}catch(FileNotFoundException e){
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		if(is==null){return 0;}

		int x=0;
		for(int i=0; i<4; i++){
			try{
				x=(x<<8)|(is.read()&255);
			}catch(IOException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		try{is.close();}catch(IOException e){e.printStackTrace();}
//		System.err.println("Returning "+x);
//		new Exception().printStackTrace();
		return x;
	}
	
	/** Selects native BAM input when allowed and either preferred or both external tools are absent. */
	public static boolean nativeBamIn(){
		return ALLOW_NATIVE_BAM_IN && (PREFER_NATIVE_BAM_IN || (!Data.SAMBAMBA() && !Data.SAMTOOLS()));
	}
	
	/** Selects native BAM output when allowed and either preferred or both external tools are absent. */
	public static boolean nativeBamOut(){
		return ALLOW_NATIVE_BAM_OUT && (PREFER_NATIVE_BAM_OUT || (!Data.SAMBAMBA() && !Data.SAMTOOLS()));
	}
	
	/** Selects native BGZF input when allowed and either preferred or threaded bgzip is unavailable. */
	public static boolean nativeBgzfIn(){
		return ALLOW_NATIVE_BGZF && (PREFER_NATIVE_BGZF_IN || !Data.BGZIP_THREADED());
	}
	
	/** Selects native BGZF output when allowed and either preferred at {@code ZIPLEVEL<=5} or
	 * threaded bgzip is unavailable; does not itself open output. */
	public static boolean nativeBgzfOut(){
		return ALLOW_NATIVE_BGZF && ((PREFER_NATIVE_BGZF_OUT && ZIPLEVEL<=5) || !Data.BGZIP_THREADED());
	}
	
	/** Counters {active, waiting, running} for InThread wrappers; active=waiting+running.
	 * This array is also the synchronization monitor. External replacement/mutation is
	 * unsupported by the accounting invariants, despite its public mutable declaration.
	 */
	public static int[] activeThreads={0, 0, 0};
	/** Maximum number of concurrent write threads allowed */
	public static int maxWriteThreads=Shared.threads();
	
	/** Whether to enable verbose debugging output for I/O operations */
	public static boolean verbose=false;
	
	/** Global raw-mode request; honored by selected factory paths, not all earlier codec branches. */
	public static boolean RAWMODE=false; //Honored by selected factory paths; see the selection contracts.

	//For killing subprocesses that are neither compression nor samtools
	/** Permits killProcess lookup for names without a recognized compression/.bam suffix. */
	public static boolean FORCE_KILL=false;

	/** Whether to use external gzip command for compression */
	public static boolean USE_GZIP=false;
	/** Whether to use bgzip for block-gzip compression */
	public static boolean USE_BGZIP=true;
	/** Whether to compress output in bgzf */
	public static boolean USE_BGZF=true;
	/** Whether to use pigz for parallel gzip compression */
	public static boolean USE_PIGZ=true;
	/** Whether to use external gunzip command for decompression */
	public static boolean USE_GUNZIP=false;
	/** Whether to use unbgzip for block-gzip decompression */
	public static boolean USE_UNBGZIP=true;
	/** Whether to use unpigz for parallel gzip decompression */
	public static boolean USE_UNPIGZ=true;
	
	public static boolean ALLOW_NATIVE_BGZF=true;
	public static boolean PREFER_NATIVE_BGZF_IN=true;
	public static boolean PREFER_NATIVE_BGZF_OUT=false;

	public static boolean USE_READ_STREAM_SAM_WRITER=true;
	public static boolean ALLOW_NATIVE_BAM_IN=true;
	public static boolean ALLOW_NATIVE_BAM_OUT=true;
	public static boolean PREFER_NATIVE_BAM_IN=true;
	public static boolean PREFER_NATIVE_BAM_OUT=true;
	
	/** Enters gzip external-selection logic even without permission; earlier BGZF preference can still win. */
	public static boolean FORCE_PIGZ=false;
	/** Requests the early BGZF path when that backend is enabled/available. */
	public static boolean FORCE_BGZIP=false;
	
	/** Prefers bgzip over other compression tools when available */
	public static boolean PREFER_BGZIP=true;
	/** Prefers unbgzip over other decompression tools when available */
	public static boolean PREFER_UNBGZIP=true;
	
	/** Whether to use external bzip2 command */
	public static boolean USE_BZIP2=true;
	/** Whether to use pbzip2 for parallel bzip2 compression */
	public static boolean USE_PBZIP2=true;
	/** Whether to use lbzip2 for bzip2 compression */
	public static boolean USE_LBZIP2=true;
	/** Whether to use DSRC compression for FASTQ files */
	public static boolean USE_DSRC=true;
	/** Whether to use FQZ compression for FASTQ files */
	public static boolean USE_FQZ=true;
	/** Whether to use Alapy compression */
	public static boolean USE_ALAPY=true;
	/** Whether to use zstd compression for .zst files. [fileIO/ReadWrite#002]: deconflated from USE_ALAPY (they were
	 * accidentally sharing a flag). Gates both .zst input (getInputStream) and output (getOutputStream). */
	public static boolean USE_ZSTD=true;
	/** Whether to use sambamba for BAM file processing */
	public static boolean USE_SAMBAMBA=true;
	/** Returns true if both USE_SAMBAMBA is enabled and sambamba is available */
	public static boolean SAMBAMBA(){return USE_SAMBAMBA && Data.SAMBAMBA();}
	
//	public static boolean SAMTOOLS_IGNORE_UNMAPPED_INPUT=false;
	/** SAM flags to ignore when processing BAM files with samtools */
	public static int SAMTOOLS_IGNORE_FLAG=0;
	/** SAM flag constant for unmapped reads */
	public static final int SAM_UNMAPPED=0x4;
	/** SAM flag constant for duplicate reads */
	public static final int SAM_DUPLICATE=0x400;
	public static final int SAM_SUPPLEMENTARY=0x800;
	/** SAM flag constant for secondary alignments */
	public static final int SAM_SECONDARY=0x100;
	/** SAM flag constant for reads that failed quality checks */
	public static final int SAM_QFAIL=0x200;
	
	/** Whether bzip2 format processing is enabled */
	public static boolean PROCESS_BZ2=true;
	/** Whether XZ format processing is enabled (currently disabled) */
	public static final boolean PROCESS_XZ=false;
	
	/** Default input buffer size in bytes */
	public static final int INBUF=65536;
	/** Default output buffer size in bytes */
	public static final int OUTBUF=65536;

	/** Gzip compression level */
	public static int ZIPLEVEL=4;
	/** Bzip2 compression level */
	public static int BZIPLEVEL=9;
	private static int MAX_ZIP_THREADS=96;
	/** Maximum threads for samtools operations */
	public static int MAX_SAMTOOLS_THREADS=64;
	/** Block size in KB for pigz compression */
	public static int PIGZ_BLOCKSIZE=128;
	/** Number of iterations for pigz optimization (-1 = default) */
	public static int PIGZ_ITERATIONS=-1;

	/** Whether zip thread multiplier has been explicitly set */
	public static boolean SET_ZIP_THREAD_MULT=false;
	/** Whether zip threads have been explicitly set */
	public static boolean SET_ZIP_THREADS=false;
	
	/** Returns maximum number of threads allowed for compression operations */
	public static int MAX_ZIP_THREADS(){return MAX_ZIP_THREADS;}
	/** Sets maximum number of compression threads with bounds checking.
	 * @param x Number of threads (clamped to 1-96 range) */
	public static void setZipThreads(int x){
//		System.err.println("setZipThreads("+x+")");new Exception().printStackTrace();
		MAX_ZIP_THREADS=Tools.mid(1, x, 96);
		SET_ZIP_THREADS=true;
	}
	/** Returns compression thread multiplier factor */
	public static float ZIP_THREAD_MULT(){return ZIP_THREAD_MULT;}
	/** Sets compression thread multiplier with bounds checking.
	 * @param x Multiplier factor (clamped to 0.125-1.0 range) */
	public static void setZipThreadMult(float x){
//		System.err.println("setZipThreadMult("+x+")");new Exception().printStackTrace();
		ZIP_THREAD_MULT=Tools.mid(0.125f, x, 1f);
		SET_ZIP_THREAD_MULT=true;
	}
	private static float ZIP_THREAD_MULT=1f;
	/**
	 * Whether compression level can be automatically adjusted based on thread count
	 */
	public static boolean ALLOW_ZIPLEVEL_CHANGE=true;
	
	/** System-specific file separator character */
	public static final String FILESEP=System.getProperty("file.separator");

	/** Monitor shared by character writes (including Async) and synchronous object writes. */
	private static final String diskSync=new String("DISKSYNC");
	
	/** Legacy debugging set; current load-tracking calls in this class are commented out. */
	public static final HashSet<String> loadedFiles=new HashSet<String>();

	/** BGZF aliases need the same gzip-compatible stream dispatch as the suffix inventory below. */
	private static boolean gzipSuffix(String name){
		return name.endsWith(".gz") || name.endsWith(".gzip") || name.endsWith(".bgz") || name.endsWith(".bgzip");
	}
	private static final String[] compressedExtensions=new String[]{".gz", ".gzip", ".bgz", ".bgzip", ".zip", ".bz2", ".xz", ".dsrc", ".fqz", ".ac", ".7z", ".zst"};
	private static final String[] compressedExtensionMap=new String[]{"gz", "gz", "gz", "gz", "zip", "bz2", "xz", "dsrc", "fqz", "ac", "7z", "zst"};

//	private static HashMap<String, Process> inputProcesses=new HashMap<String, Process>(8);
//	private static HashMap<String, Process> outputProcesses=new HashMap<String, Process>(8);
	/** Filename-keyed processes; normal access is synchronized on this map. */
	private static HashMap<String, Process> processMap=new HashMap<String, Process>(8);
	/** Filename-keyed pipe helpers; normal access is synchronized on this map. */
	private static HashMap<String, ArrayList<PipeThread>> pipeThreadMap=new HashMap<String, ArrayList<PipeThread>>(8);
	
}
