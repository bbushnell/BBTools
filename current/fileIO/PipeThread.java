package fileIO;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/**
 * Copies bytes between streams using a 32 KiB buffer when the caller starts this thread.
 * Also permits synchronous execution through run(). The copy stops on a nonpositive
 * read or a termination request observed after a read. Reads and writes may block.
 * Normal completion and caught IOExceptions lead to stream cleanup; unchecked failures
 * can skip cleanup and completion publication. Exact System streams remain borrowed.
 *
 * @author Brian Bushnell
 * @date Jan 22, 2013
 */
public class PipeThread extends Thread{
	
//	public PipeThread(InputStream is_){this(is_, System.err);}
	
	/** Stores nonnull streams without starting the thread.
	 * @param is_ Source; closed after copying unless it is exactly System.in
	 * @param os_ Destination; flushed/closed after copying unless exactly System.out or System.err
	 * @throws RuntimeException If either stream is null
	 */
	public PipeThread(InputStream is_, OutputStream os_){
		is=is_;
		os=os_;
		if(is==null){throw new RuntimeException("Null input stream.");}
		if(os==null){throw new RuntimeException("Null output stream.");}
//		synchronized(list){list.add(this);}
	}
	
	/** Copies bytes, attempts cleanup and publishes completion after handled IOExceptions.
	 * Exact System.out/System.err receive no explicit final flush here. Output-cleanup
	 * errors reported by ReadWrite.close contribute to errorState.
	 */
	@Override
	public void run(){
		//TODO: Probable bug - FIO-068: unchecked failures skip remaining cleanup and finished publication; failure policy remains unresolved.
		final byte[] buf=new byte[32768];
		try{
			for(int len=is.read(buf); !finished && len>0; len=is.read(buf)){
				os.write(buf, 0, len);
			}
		}catch(IOException e){
			errorState=true;
			e.printStackTrace();
		}
		
		if(is!=System.in){
			try{
				is.close();
			}catch(IOException e){
				errorState=true;
				e.printStackTrace();
			}
		}
		
		if(os!=System.out && os!=System.err){
			errorState|=ReadWrite.close(os);
		}
		
		synchronized(this){
			finished=true;
			this.notify();
		}
	}
	
	/** Returns the stop/completion flag, not a thread-termination or cleanup guarantee.
	 * @return true after terminate() requests a stop or run() reaches its final block
	 */
	public boolean finished(){
		synchronized(this){
			return finished;
		}
	}
	
	/** Sets the stop flag and interrupts this thread without joining or closing streams.
	 * The loop tests the flag only after reading; an interrupt does not guarantee that
	 * the particular stream unblocks. An in-progress write may finish before another read.
	 */
	public void terminate(){
		synchronized(this){
			if(!finished){
				finished=true;
				interrupt();
			}
		}
	}
	
	/** Source; exact System.in is borrowed, other inputs are closed after handled copy failures or EOF. */
	public final InputStream is;
	/** Destination; exact System.out/System.err are borrowed, other outputs receive ReadWrite.close. */
	public final OutputStream os;
	/** Stop requested or final run block reached; neither state proves thread termination. */
	private volatile boolean finished=false;
	/** Sticky flag for caught copy/input-close IOExceptions and reported output-cleanup errors; PrintStream retained errors are not included. */
	public volatile boolean errorState=false;
	
//	private static ArrayList<PipeThread> list=new ArrayList<PipeThread>(8);
	
}
