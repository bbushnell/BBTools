package fileIO;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/**
 * Listens to an output stream and copies it to an input stream.
 * For example, redirects the error stream of some process to stderr.
 * @author Brian Bushnell
 * @date Jan 22, 2013
 *
 */
public class PipeThread extends Thread {
	
//	public PipeThread(InputStream is_){this(is_, System.err);}
	
	/**
	 * Constructs a PipeThread with specified input and output streams.
	 * Validates that both streams are non-null and initializes the thread for data transfer.
	 *
	 * @param is_ The input stream to read data from
	 * @param os_ The output stream to write data to
	 * @throws RuntimeException If either stream is null
	 */
	public PipeThread(InputStream is_, OutputStream os_){
		is=is_;
		os=os_;
		if(is==null){throw new RuntimeException("Null input stream.");}
		if(os==null){throw new RuntimeException("Null output stream.");}
//		synchronized(list){list.add(this);}
	}
	
	@Override
	public void run(){
		final byte[] buf=new byte[32768];
		try {
			for(int len=is.read(buf); !finished && len>0; len=is.read(buf)){
				os.write(buf, 0, len);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		if(is!=System.in){
			try {
				is.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		if(os!=System.out && os!=System.err){
			ReadWrite.close(os);
		}
		
		synchronized(this){
			finished=true;
			this.notify();
		}
	}
	
	/**
	 * Checks if the data transfer operation has completed.
	 * Thread-safe method using synchronization to ensure accurate status reporting.
	 * @return true if the pipe operation has finished, false otherwise
	 */
	public boolean finished(){
		synchronized(this){
			return finished;
		}
	}
	
	/**
	 * Terminates the pipe operation by setting the finished flag and interrupting the thread.
	 * Thread-safe method that ensures clean shutdown of the data transfer operation.
	 * Only interrupts if the thread has not already finished.
	 */
	public void terminate(){
		synchronized(this){
			if(!finished){
				finished=true;
				interrupt();
			}
		}
	}
	
	/** The input stream from which data is read during the pipe operation */
	public final InputStream is;
	/** The output stream to which data is written during the pipe operation */
	public final OutputStream os;
	/**
	 * Flag indicating whether the pipe operation has completed, marked volatile for thread safety
	 */
	private volatile boolean finished=false;
	
//	private static ArrayList<PipeThread> list=new ArrayList<PipeThread>(8);
	
}
