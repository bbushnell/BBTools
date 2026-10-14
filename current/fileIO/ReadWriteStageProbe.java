package fileIO;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

/**
 * Linux diagnostic for finishWriting and buffered-writer error reporting.
 * Checks the normal in-memory output, then prints normal/close/buffered error flags.
 * Expected flags are false/true/true; the two failure flags are reported, not asserted.
 * Intentionally writes a short record to /dev/full; no ordinary files are created.
 */
public final class ReadWriteStageProbe{
	/**
	 * Exercises normal closure, a controlled close failure, and the full-device sink.
	 * Starts and joins the buffered writer before printing the flags to stdout.
	 * @param args Ignored
	 * @throws Exception If an operation fails before the flags can be reported
	 */
	public static void main(String[] args) throws Exception{
		if(!new File("/dev/full").exists()){
			throw new IllegalStateException("This diagnostic requires Linux /dev/full");
		}
		ByteArrayOutputStream good=new ByteArrayOutputStream();
		good.write("valid".getBytes("UTF-8"));
		boolean normal=ReadWrite.finishWriting(null, good, null, false);
		if(normal || !"valid".equals(good.toString("UTF-8"))){throw new AssertionError("normal output changed");}
		OutputStream broken=new OutputStream(){
			@Override
			public void write(int b){}
			@Override
			public void close() throws IOException{throw new IOException("intentional close failure");}
		};
		boolean closeError=ReadWrite.finishWriting(null, broken, null, false);
		ByteStreamWriter writer=new ByteStreamWriter("/dev/full", true, false, true);
		writer.start();
		writer.print("intentional buffered write failure\n");
		boolean writerError=writer.poisonAndWait();
		System.out.println("normal_error="+normal+" close_error="+closeError+" buffered_writer_error="+writerError);
	}
}
