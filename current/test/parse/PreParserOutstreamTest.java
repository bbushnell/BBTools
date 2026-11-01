package test.parse;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import parse.PreParser;

/** Linux regression tests for repeated PreParser outstream ownership.
 * Observes /proc/self/fd; run each case in a separate JVM with a non-collecting
 * GC so collection cannot hide an abandoned stream before observation. */
public final class PreParserOutstreamTest{

	/** Runs one replacement case per JVM so baseline failures are independent.
	 * @param args Required case=file|stderr|null|borrowed */
	public static void main(String[] args) throws IOException{
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Expected exactly one case= argument");
		}
		final String name=args[0].substring("case=".length());
		final File dir=makeTempDir();
		try{
			if(name.equals("file")){checkFileReplacementClosesPrevious(dir);}
			else if(name.equals("stderr")){checkSystemReplacementClosesPrevious(dir);}
			else if(name.equals("null")){checkNullReplacementClosesPrevious(dir);}
			else if(name.equals("borrowed")){checkBorrowedDefaultIsPreserved(dir);}
			else{throw new IllegalArgumentException("Unknown case="+name);}
		}finally{
			deleteContents(dir);
			dir.delete();
		}
		System.out.println("PASS PreParserOutstreamTest "+name);
	}

	/** File replacement should close the first file while retaining the second. */
	private static void checkFileReplacementClosesPrevious(File dir) throws IOException{
		final File first=new File(dir, "first.txt");
		final File second=new File(dir, "second.txt");
		PreParser pp=null;
		try{
			pp=new PreParser(new String[] {"outstream="+first, "outstream="+second},
					System.err, PreParserOutstreamTest.class, false, true, false);
			require(!isOpen(first), "first file outstream was not closed on replacement");
			require(isOpen(second), "final file outstream was not retained");
		}finally{
			if(pp!=null){pp.outstream.close();}
		}
		requireOpenCount(dir, 0, "file-to-file replacement cleanup");
	}

	/** Borrowed System streams should displace but not inherit an owned file. */
	private static void checkSystemReplacementClosesPrevious(File dir) throws IOException{
		final File file=new File(dir, "stderr.txt");
		new PreParser(new String[] {"outstream="+file, "outstream=stderr"},
				System.err, PreParserOutstreamTest.class, false, true, false);
		require(!isOpen(file), "file outstream was not closed when stderr replaced it");
		requireOpenCount(dir, 0, "file-to-stderr replacement cleanup");
	}

	/** Null output streams are owned but should still displace an earlier owned file. */
	private static void checkNullReplacementClosesPrevious(File dir) throws IOException{
		final File file=new File(dir, "null.txt");
		final PreParser pp=new PreParser(new String[] {"outstream="+file, "outstream=null"},
				System.err, PreParserOutstreamTest.class, false, true, false);
		try{
			require(!isOpen(file), "file outstream was not closed when null replaced it");
			requireOpenCount(dir, 0, "file-to-null replacement cleanup");
		}finally{
			pp.outstream.close();
		}
	}

	/** Replacement may not close caller-owned streams. */
	private static void checkBorrowedDefaultIsPreserved(File dir) throws IOException{
		checkBorrowedDefault(dir);
		checkBorrowedSystem(dir, true);
		checkBorrowedSystem(dir, false);
	}

	/** Replacement may not close the caller-owned default stream. */
	private static void checkBorrowedDefault(File dir) throws IOException{
		final TrackingPrintStream defaultStream=new TrackingPrintStream();
		final File file=new File(dir, "default.txt");
		final PreParser pp=new PreParser(new String[] {"outstream="+file},
				defaultStream, PreParserOutstreamTest.class, false, true, false);
		try{
			require(!defaultStream.closed, "borrowed default stream was closed");
			require(isOpen(file), "file outstream was not retained");
		}finally{
			pp.outstream.close();
			defaultStream.close();
		}
		requireOpenCount(dir, 0, "borrowed-default replacement cleanup");
	}

	/** Explicit System stream selection must remain borrowed. */
	private static void checkBorrowedSystem(File dir, boolean stdout) throws IOException{
		final TrackingPrintStream stream=new TrackingPrintStream();
		final PrintStream old=stdout ? System.out : System.err;
		final File file=new File(dir, stdout ? "stdout.txt" : "stderr.txt");
		PreParser pp=null;
		try{
			setSystemStream(stdout, stream);
			pp=new PreParser(new String[] {"outstream="+(stdout ? "stdout" : "stderr"),
					"outstream="+file}, System.err, PreParserOutstreamTest.class, false, true,
					false);
			require(!stream.closed, "borrowed System stream was closed");
			require(isOpen(file), "file outstream was not retained");
		}finally{
			setSystemStream(stdout, old);
			if(pp!=null){pp.outstream.close();}
			stream.close();
		}
		requireOpenCount(dir, 0, "borrowed-System replacement cleanup");
	}

	/** Counts temp-dir descriptors currently open through /proc/self/fd. */
	private static int countOpen(File dir) throws IOException{
		final String prefix=dir.getCanonicalPath()+File.separator;
		int count=0;
		final Path fdRoot=Paths.get("/proc/self/fd");
		try(final DirectoryStream<Path> stream=Files.newDirectoryStream(fdRoot)){
			for(final Path fd : stream){
				try{
					final String target=Files.readSymbolicLink(fd).toString();
					if(target.startsWith(prefix)){count++;}
				}catch(IOException | UnsupportedOperationException ignored){
					// Descriptors can disappear while /proc/self/fd is being scanned.
				}
			}
		}
		return count;
	}

	/** True iff this JVM still has the requested file open. */
	private static boolean isOpen(File file) throws IOException{
		final String target=file.getCanonicalPath();
		final Path fdRoot=Paths.get("/proc/self/fd");
		try(final DirectoryStream<Path> stream=Files.newDirectoryStream(fdRoot)){
			for(final Path fd : stream){
				try{
					if(target.equals(Files.readSymbolicLink(fd).toString())){return true;}
				}catch(IOException | UnsupportedOperationException ignored){
					// Descriptors can disappear while /proc/self/fd is being scanned.
				}
			}
		}
		return false;
	}

	/** Requires the current temp-dir descriptor count to equal an expected value. */
	private static void requireOpenCount(File dir, int expected, String label) throws IOException{
		final int observed=countOpen(dir);
		if(observed!=expected){
			throw new AssertionError(label+" left "+observed+
					" open descriptors, expected "+expected);
		}
	}

	/** Replaces a System stream without closing the previous value. */
	private static void setSystemStream(boolean stdout, PrintStream stream){
		if(stdout){System.setOut(stream);}
		else{System.setErr(stream);}
	}

	/** Requires a condition. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Makes a temporary directory for descriptor-target filtering. */
	private static File makeTempDir() throws IOException{
		final File dir=File.createTempFile("preparser-outstream-", "");
		if(!dir.delete() || !dir.mkdirs()){
			throw new IOException("Could not create temp directory "+dir);
		}
		return dir.getCanonicalFile();
	}

	/** Deletes ordinary files in a temporary test directory. */
	private static void deleteContents(File dir){
		final File[] files=dir.listFiles();
		if(files==null){return;}
		for(final File f : files){f.delete();}
	}

	/** PrintStream that records unexpected closure by PreParser. */
	private static final class TrackingPrintStream extends PrintStream{
		TrackingPrintStream(){
			super(new ByteArrayOutputStream());
		}

		@Override
		public void close(){
			closed=true;
			super.close();
		}

		/** Whether close has been invoked on this borrowed-stream control. */
		boolean closed=false;
	}

}
