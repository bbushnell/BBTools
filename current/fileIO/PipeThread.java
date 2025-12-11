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
	
	public boolean finished(){
		synchronized(this){
			return finished;
		}
	}
	
	public void terminate(){
		synchronized(this){
			if(!finished){
				finished=true;
				interrupt();
			}
		}
	}
	
	public final InputStream is;
	public final OutputStream os;
	private volatile boolean finished=false;
	
//	private static ArrayList<PipeThread> list=new ArrayList<PipeThread>(8);
	
}
