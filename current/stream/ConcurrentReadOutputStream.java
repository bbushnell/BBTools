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
	
	/**
	 * Protected constructor for subclasses.
	 * Initializes file formats and determines if output should be ordered.
	 * @param ff1_ Primary file format
	 * @param ff2_ Secondary file format (may be null)
	 */
	ConcurrentReadOutputStream(FileFormat ff1_, FileFormat ff2_){
		ff1=ff1_;
		ff2=ff2_;
		ordered=(ff1==null ? true : ff1.ordered());
	}
	
	/** Must be called before writing to the stream */
	public abstract void start();
	
	/** Returns whether the stream has been started */
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
	
	/** Closes the output stream and releases resources */
	public abstract void close();
	
	/** Waits for all writing threads to complete */
	public abstract void join();
	
	/** Resets the next list ID counter for ordered output */
	public abstract void resetNextListID();
	
	/** Returns the filename of the primary output file */
	public abstract String fname();
	
	/** Return true if this stream has detected an error */
	public abstract boolean errorState();

	/** Returns true if the stream finished successfully without errors */
	public abstract boolean finishedSuccessfully();
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Returns the total number of bases written across all output streams.
	 * Sums bases from both primary and secondary ReadStreamWriters.
	 * @return Total bases written
	 */
	public long basesWritten(){
		long x=0;
		ReadStreamWriter rsw1=getRS1();
		ReadStreamWriter rsw2=getRS2();
		if(rsw1!=null){x+=rsw1.basesWritten();}
		if(rsw2!=null){x+=rsw2.basesWritten();}
		return x;
	}
	
	/**
	 * Returns the total number of reads written across all output streams.
	 * Sums reads from both primary and secondary ReadStreamWriters.
	 * @return Total reads written
	 */
	public long readsWritten(){
		long x=0;
		ReadStreamWriter rsw1=getRS1();
		ReadStreamWriter rsw2=getRS2();
		if(rsw1!=null){x+=rsw1.readsWritten();}
		if(rsw2!=null){x+=rsw2.readsWritten();}
		return x;
	}
	
	/** Returns the primary ReadStreamWriter */
	public abstract ReadStreamWriter getRS1();
	/** Returns the secondary ReadStreamWriter (may be null) */
	public abstract ReadStreamWriter getRS2();
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	public final FileFormat ff1, ff2;
	/** Whether output should maintain input order */
	public final boolean ordered;
	
	/** Tracks whether an error has occurred in the stream */
	boolean errorState=false;
	/** Tracks whether the stream finished successfully */
	boolean finishedSuccessfully=false;
	/** Tracks whether the stream has been started */
	boolean started=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Global flag for verbose output during stream operations */
	public static boolean verbose=false;
	
}
