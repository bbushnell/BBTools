package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import shared.Shared;

/**
 * Abstract base for concurrent read output streams that wrap ReadStreamWriters.
 * Supports single/paired outputs, shared headers, ordered writing, and MPI passthrough.
 * @author Brian Bushnell
 * @date Jan 26, 2015
 */
public abstract class ConcurrentReadOutputStream {
	
	/*--------------------------------------------------------------*/
	/*----------------           Factory            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates a single-file output stream with shared-header option.
	 *
	 * @param ff1 Primary output format
	 * @param rswBuffers Max buffered lists per writer
	 * @param header Header text to prepend
	 * @param useSharedHeader Whether to write the shared header
	 * @return ConcurrentReadOutputStream instance
	 */
	public static ConcurrentReadOutputStream getStream(FileFormat ff1, int rswBuffers, CharSequence header, boolean useSharedHeader){
		return getStream(ff1, null, null, null, rswBuffers, header, useSharedHeader, Shared.USE_MPI, Shared.MPI_KEEP_ALL);
	}
	
	/**
	 * Creates a paired-file output stream.
	 *
	 * @param ff1 Read 1 format
	 * @param ff2 Read 2 format (optional)
	 * @param rswBuffers Max buffered lists per writer
	 * @param header Header text to prepend
	 * @param useSharedHeader Whether to write the shared header
	 * @return ConcurrentReadOutputStream instance
	 */
	public static ConcurrentReadOutputStream getStream(FileFormat ff1, FileFormat ff2, int rswBuffers, CharSequence header, boolean useSharedHeader){
		return getStream(ff1, ff2, null, null, rswBuffers, header, useSharedHeader, Shared.USE_MPI, Shared.MPI_KEEP_ALL);
	}
	
	/**
	 * Creates an output stream with optional quality files for paired output.
	 *
	 * @param ff1 Read 1 format
	 * @param ff2 Read 2 format (optional)
	 * @param qf1 Quality file 1 (optional)
	 * @param qf2 Quality file 2 (optional)
	 * @param rswBuffers Max buffered lists per writer
	 * @param header Header text to prepend
	 * @param useSharedHeader Whether to write the shared header
	 * @return ConcurrentReadOutputStream instance
	 */
	public static ConcurrentReadOutputStream getStream(FileFormat ff1, FileFormat ff2, String qf1, String qf2,
			int rswBuffers, CharSequence header, boolean useSharedHeader){
		return getStream(ff1, ff2, qf1, qf2, rswBuffers, header, useSharedHeader, Shared.USE_MPI, Shared.MPI_KEEP_ALL);
	}
	
	/**
	 * Primary factory creating standard or MPI-backed concurrent output streams.
	 *
	 * @param ff1 Read 1 format (required)
	 * @param ff2 Read 2 format (optional)
	 * @param qf1 Quality file 1 (optional)
	 * @param qf2 Quality file 2 (optional)
	 * @param rswBuffers Max buffered lists per writer
	 * @param header Header text to prepend
	 * @param useSharedHeader Write shared header (SAM)
	 * @param mpi Use MPI-backed stream
	 * @param keepAll In MPI mode, keep all reads instead of a fraction
	 * @return ConcurrentReadOutputStream implementation
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
	 * Protected base constructor setting formats and ordered flag.
	 * @param ff1_ Primary output format
	 * @param ff2_ Secondary output format (may be null)
	 */
	ConcurrentReadOutputStream(FileFormat ff1_, FileFormat ff2_){
		ff1=ff1_;
		ff2=ff2_;
		ordered=(ff1==null ? true : ff1.ordered());
	}
	
	/** Starts underlying writers/threads; must be called before adding reads. */
	public abstract void start();
	
	/** Indicates whether the stream has been started.
	 * @return true if start() was called */
	public final boolean started(){return started;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Enqueues a list of reads to be written; ordered streams enforce listnum order.
	 * @param list Reads to write
	 * @param listnum Sequential list number starting at 0
	 */
	public abstract void add(ArrayList<Read> list, long listnum);
	
	/** Closes the output stream and releases resources. */
	public abstract void close();
	
	/** Waits for all writer threads to finish. */
	public abstract void join();
	
	/** Resets ordered output list numbering back to zero. */
	public abstract void resetNextListID();
	
	/** Returns the primary output filename.
	 * @return Output filename */
	public abstract String fname();
	
	/** Indicates whether an error has been detected.
	 * @return true if error occurred */
	public abstract boolean errorState();

	/** Indicates whether the stream completed without errors.
	 * @return true if finished cleanly */
	public abstract boolean finishedSuccessfully();
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns total bases written across all writers.
	 * @return Base count written */
	public long basesWritten(){
		long x=0;
		ReadStreamWriter rsw1=getRS1();
		ReadStreamWriter rsw2=getRS2();
		if(rsw1!=null){x+=rsw1.basesWritten();}
		if(rsw2!=null){x+=rsw2.basesWritten();}
		return x;
	}
	
	/** Returns total reads written across all writers.
	 * @return Read count written */
	public long readsWritten(){
		long x=0;
		ReadStreamWriter rsw1=getRS1();
		ReadStreamWriter rsw2=getRS2();
		if(rsw1!=null){x+=rsw1.readsWritten();}
		if(rsw2!=null){x+=rsw2.readsWritten();}
		return x;
	}
	
	/** Returns the primary ReadStreamWriter.
	 * @return Primary writer or null */
	public abstract ReadStreamWriter getRS1();
	/** Returns the secondary ReadStreamWriter, if any.
	 * @return Secondary writer or null */
	public abstract ReadStreamWriter getRS2();
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Secondary output file format (may be null). */
	/** Primary output file format. */
	public final FileFormat ff1, ff2;
	/** Whether output must preserve list order. */
	public final boolean ordered;
	
	/** Tracks whether an error was encountered. */
	boolean errorState=false;
	/** Tracks whether writing completed successfully. */
	boolean finishedSuccessfully=false;
	/** Tracks whether start() has been called. */
	boolean started=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Enables verbose logging for stream operations. */
	public static boolean verbose=false;
	
}
