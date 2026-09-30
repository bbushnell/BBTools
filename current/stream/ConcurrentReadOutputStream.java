package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import shared.Shared;

/**
 * Abstract base for concurrent read output streams that wrap ReadStreamWriters.
 * Provides factories and lifecycle contracts for local single/paired output,
 * shared headers and ordered lists. Legacy distributed routing is retained;
 * the true-MPI implementation remains fenced off in the factory.
 * @author Brian Bushnell
 * @date Jan 26, 2015
 */
public abstract class ConcurrentReadOutputStream{
	
	/*--------------------------------------------------------------*/
	/*----------------           Factory            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates an unstarted single-file output stream with shared-header option.
	 * Uses the current Shared MPI settings.
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
	 * Creates an unstarted output stream with an optional second file.
	 * Uses the current Shared MPI settings; writer selection depends on ff1.
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
	 * Creates an unstarted output stream with optional second and quality files.
	 * Uses the current Shared MPI settings; the selected writer determines which
	 * secondary and quality-file options apply.
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
	 * Primary factory creating an unstarted local or legacy distributed stream.
	 * Local output delegates writer selection to ConcurrentGenericReadOutputStream.
	 * The true-MPI route is deliberately disabled by the assertion below.
	 *
	 * @param ff1 Read 1 format (required)
	 * @param ff2 Read 2 format (optional)
	 * @param qf1 Quality file 1 (optional)
	 * @param qf2 Quality file 2 (optional)
	 * @param rswBuffers Max buffered lists per writer
	 * @param header Header text to prepend
	 * @param useSharedHeader Write shared header (SAM)
	 * @param mpi Use MPI-backed stream
	 * @param keepAll Retained compatibility argument; unused by this factory
	 * @return ConcurrentReadOutputStream implementation
	 */
	public static ConcurrentReadOutputStream getStream(FileFormat ff1, FileFormat ff2, String qf1, String qf2,
			int rswBuffers, CharSequence header, boolean useSharedHeader, final boolean mpi, final boolean keepAll){
		//MPI path: rank 0 builds the real local ConcurrentGenericReadOutputStream; other ranks get null. Both wrap in ConcurrentReadOutputStreamD (the distributed cros). Non-MPI (else branch) returns the generic cros directly.
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
				//Deliberate fence (NOT a bug): the true-MPI cros is unfinished/disabled. The message IS the acceptance test - validate before uncommenting; do not delete as "obsolete".
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
	 * Package-private base constructor storing formats and the ordered flag.
	 * A null primary format defaults ordered to true; otherwise ff1_.ordered()
	 * is captured once. Subclasses may impose stronger format requirements.
	 * @param ff1_ Primary output format (may be null in this base constructor)
	 * @param ff2_ Secondary output format (may be null)
	 */
	ConcurrentReadOutputStream(FileFormat ff1_, FileFormat ff2_){
		ff1=ff1_;
		ff2=ff2_;
		ordered=(ff1==null ? true : ff1.ordered());
	}
	
	/** Starts underlying writers/threads; must be called before adding reads. */
	public abstract void start();
	
	/** Returns the started flag maintained by the implementation.
	 * @return Current stored start state; this getter does not wait */
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
	
	/** Requests output shutdown after all producer adds; call join() to await completion. */
	public abstract void close();

	/**
	 * Abandons buffered output after an upstream failure and wakes blocked add calls.
	 * Implementations must mark an error and tolerate repeated calls.
	 */
	public abstract void abort();
	
	/** Waits for writer completion after shutdown has been requested. */
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
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Sums the currently reported base counts of nonnull primary/secondary writers.
	 * Does not wait for pending output; inspect totals after shutdown and join.
	 * @return Sum of observed writer base counts, or zero if both writers are null */
	public long basesWritten(){
		long x=0;
		ReadStreamWriter rsw1=getRS1();
		ReadStreamWriter rsw2=getRS2();
		if(rsw1!=null){x+=rsw1.basesWritten();}
		if(rsw2!=null){x+=rsw2.basesWritten();}
		return x;
	}
	
	/** Sums the currently reported read counts of nonnull primary/secondary writers.
	 * Does not wait for pending output; inspect totals after shutdown and join.
	 * @return Sum of observed writer read counts, or zero if both writers are null */
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
	
	/**
	 * Primary and secondary output formats; either may be null in the base class.
	 * The generic local implementation requires ff1.
	 * Historical #001: combined declaration replaces swapped, stacked field Javadocs.
	 */
	public final FileFormat ff1, ff2;
	/** Ordering captured from ff1 at construction; true when ff1 is null. */
	public final boolean ordered;
	
	/** Tracks whether an error was encountered. */
	boolean errorState=false;
	/** Tracks whether writing completed successfully. */
	boolean finishedSuccessfully=false;
	/** Start state maintained by concrete implementations. */
	boolean started=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Enables verbose logging for stream operations. */
	public static boolean verbose=false;
	
}
