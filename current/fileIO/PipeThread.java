package fileIO;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/**
 * A specialized thread for efficiently copying data between input and output streams.
 * Facilitates stream redirection and pipeline management by running data transfer operations
 * in the background without blocking the calling thread.
 *
 * @author Brian Bushnell
 * @date Jan 22, 2013
 */
public class PipeThread extends Thread {
	
//	public PipeThread(InputStream is_){this(is_, System.err);}
	
	/** Constructs a PipeThread that copies all bytes from is_ to os_ in the background; both streams must be non-null. */
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
			errorState=true;
			e.printStackTrace();
		}
		
		if(is!=System.in){
			try {
				is.close();
			} catch (IOException e) {
				errorState=true;
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
	
	/** @return true once the copy loop has completed or been terminated. */
	public boolean finished(){
		synchronized(this){
			return finished;
		}
	}
	
	/** Signals the copy loop to stop at the next read boundary and interrupts this thread. */
	public void terminate(){
		synchronized(this){
			if(!finished){
				finished=true;
				interrupt();
			}
		}
	}
	
	/** The source stream this thread reads from. */
	public final InputStream is;
	/** The destination stream this thread writes to. */
	public final OutputStream os;
	private volatile boolean finished=false;
	/** Set true if an IOException occurs during copying or stream close; lets callers detect a failed/partial pipe. */
	public volatile boolean errorState=false;
	
//	private static ArrayList<PipeThread> list=new ArrayList<PipeThread>(8);
	
}
