package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import shared.Shared;

/**
 * Abstract superclass for ConcurrentReadOutputStream implementations.
 * These manage ReadStreamWriters, which write reads to a file in their own thread.
 * ConcurrentReadOutputStreams allow paired reads output to twin files to be treated as a single stream.
 * @author Brian Bushnell
 * @date Jan 26, 2015
 *
 */
public abstract class ConcurrentReadOutputStream {
	
	/*--------------------------------------------------------------*/
	/*----------------           Factory            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** @See primary method */
	public static ConcurrentReadOutputStream getStream(FileFormat ff1, int rswBuffers, CharSequence header, boolean useSharedHeader){
		return getStream(ff1, null, null, null, rswBuffers, header, useSharedHeader, Shared.USE_MPI, Shared.MPI_KEEP_ALL);
	}
	
	/** @See primary method */
	public static ConcurrentReadOutputStream getStream(FileFormat ff1, FileFormat ff2, int rswBuffers, CharSequence header, boolean useSharedHeader){
		return getStream(ff1, ff2, null, null, rswBuffers, header, useSharedHeader, Shared.USE_MPI, Shared.MPI_KEEP_ALL);
	}
	
	/** @See primary method */
	public static ConcurrentReadOutputStream getStream(FileFormat ff1, FileFormat ff2, String qf1, String qf2,
			int rswBuffers, CharSequence header, boolean useSharedHeader){
		return getStream(ff1, ff2, qf1, qf2, rswBuffers, header, useSharedHeader, Shared.USE_MPI, Shared.MPI_KEEP_ALL);
	}
	
	/**
	 * Create a ConcurrentReadOutputStream.
	 * @param ff1 Read 1 file (required)
	 * @param ff2 Read 2 file (optional)
	 * @param qf1 Qual file 1 (optional)
	 * @param qf2 Qual file 2 (optional)
	 * @param rswBuffers Maximum number of lists to buffer for each ReadStreamWriter
	 * @param header A header to write to each output file before anything else
	 * @param useSharedHeader Write the shared header to each output file (mainly for sam output)
	 * @param mpi True if MPI will be used
	 * @param keepAll In MPI mode, tells this stream to keep all reads instead of just a fraction
	 * @return
	 */
	public static ConcurrentReadOutputStream getStream(FileFormat ff1, FileFormat ff2, String qf1, String qf2,
			int rswBuffers, CharSequence header, boolean useSharedHeader, final boolean mpi, final boolean keepAll){
		if(mpi){
			final int rank=Shared.MPI_RANK;
			final ConcurrentReadOutputStream cros0;
			if(rank==0){
				cros0=new ConcurrentGenericReadOutputStream(ff1, ff2, qf1, qf2, rswBuffers, header, useSharedHeader);
			}else{
				cros0=null;
			}
			final ConcurrentReadOutputStream crosD;
			if(Shared.USE_CRISMPI){
				assert(false) : "To support MPI, uncomment this.";
				crosD=null;
//				crosD=new ConcurrentReadOutputStreamMPI(cros0, rank==0);
			}else{
				crosD=new ConcurrentReadOutputStreamD(cros0, rank==0);
			}
			return crosD;
		}else{
			return new ConcurrentGenericReadOutputStream(ff1, ff2, qf1, qf2, rswBuffers, header, useSharedHeader);
		}
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	ConcurrentReadOutputStream(FileFormat ff1_, FileFormat ff2_){
		ff1=ff1_;
		ff2=ff2_;
		ordered=(ff1==null ? true : ff1.ordered());
	}
	
	/** Must be called before writing to the stream */
	public abstract void start();
	
	public final boolean started(){return started;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Enqueue this list to be written.
	 * @param list List of reads
	 * @param listnum A number, starting at 0.  In ordered mode, lists will only be written in numeric order, regardless of adding order.
	 */
	public abstract void add(ArrayList<Read> list, long listnum);
	
	public abstract void close();
	
	public abstract void join();
	
	public abstract void resetNextListID();
	
	public abstract String fname();
	
	/** Return true if this stream has detected an error */
	public abstract boolean errorState();

	public abstract boolean finishedSuccessfully();
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	public long basesWritten(){
		long x=0;
		ReadStreamWriter rsw1=getRS1();
		ReadStreamWriter rsw2=getRS2();
		if(rsw1!=null){x+=rsw1.basesWritten();}
		if(rsw2!=null){x+=rsw2.basesWritten();}
		return x;
	}
	
	public long readsWritten(){
		long x=0;
		ReadStreamWriter rsw1=getRS1();
		ReadStreamWriter rsw2=getRS2();
		if(rsw1!=null){x+=rsw1.readsWritten();}
		if(rsw2!=null){x+=rsw2.readsWritten();}
		return x;
	}
	
	public abstract ReadStreamWriter getRS1();
	public abstract ReadStreamWriter getRS2();
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	public final FileFormat ff1, ff2;
	public final boolean ordered;
	
	boolean errorState=false;
	boolean finishedSuccessfully=false;
	boolean started=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean verbose=false;
	
}
