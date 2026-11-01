package shared;

import java.util.concurrent.locks.ReentrantReadWriteLock;

/**
 * Synchronized accessors for selected mutable configuration values in BBTools.
 * The private ReadWriteLock allows multiple concurrent readers and exclusive
 * writers among callers that use this class.  Direct access to the underlying
 * Shared fields is outside this lock.
 * 
 * @author Brian Bushnell
 * @author Isla
 * @date June 18, 2025
 */
public class SyncHeart{
	
	/** ReadWriteLock for thread-safe access to configuration values. */
	private static final ReentrantReadWriteLock configLock=new ReentrantReadWriteLock();
	
	/** Acquires read lock for concurrent access through SyncHeart. */
	private static final void readLock(){
		configLock.readLock().lock();
	}
	
	/** Releases read lock. */
	private static final void readUnlock(){
		configLock.readLock().unlock();
	}
	
	/** Acquires write lock for access through SyncHeart. */
	private static final void writeLock(){
		configLock.writeLock().lock();
	}
	
	/** Releases write lock. */
	private static final void writeUnlock(){
		configLock.writeLock().unlock();
	}
	
	// Threading and Performance Configuration
	
	/** Locked setter for amino acid input mode. */
	public static void setAminoIn(boolean b){
		writeLock();
		try{Shared.AMINO_IN=b;}
		finally{writeUnlock();}
	}
	
	/** Locked getter for amino acid input mode. */
	public static boolean aminoIn(){
		readLock();
		try{return Shared.AMINO_IN;}
		finally{readUnlock();}
	}
	
	/** Locked setter for low memory mode. */
	public static void setLowMemory(boolean b){
		writeLock();
		try{Shared.LOW_MEMORY=b;}
		finally{writeUnlock();}
	}
	
	/** Locked getter for low memory mode. */
	public static boolean lowMemory(){
		readLock();
		try{return Shared.LOW_MEMORY;}
		finally{readUnlock();}
	}
	
	/** Locked setter for garbage collection before memory printing. */
	public static void setGcBeforePrintMemory(boolean b){
		writeLock();
		try{Shared.GC_BEFORE_PRINT_MEMORY=b;}
		finally{writeUnlock();}
	}
	
	/** Locked getter for garbage collection before memory printing. */
	public static boolean gcBeforePrintMemory(){
		readLock();
		try{return Shared.GC_BEFORE_PRINT_MEMORY;}
		finally{readUnlock();}
	}
	
	/** Locked setter for parallel sort usage. */
	public static void setParallelSort(boolean b){
		writeLock();
		try{Shared.parallelSort=b;}
		finally{writeUnlock();}
	}
	
	/** Locked getter for parallel sort usage. */
	public static boolean parallelSort(){
		readLock();
		try{return Shared.parallelSort;}
		finally{readUnlock();}
	}
	
	/** Locked setter for SIMD optimization usage. */
	public static void setSimd(boolean b){
		writeLock();
		try{Shared.SIMD=b;}
		finally{writeUnlock();}
	}
	
	/** Locked getter for SIMD optimization usage. */
	public static boolean simd(){
		readLock();
		try{return Shared.SIMD;}
		finally{readUnlock();}
	}
	
	// MPI Configuration
	
	/** Locked setter for MPI usage. */
	public static void setUseMpi(boolean b){
		writeLock();
		try{Shared.USE_MPI=b;}
		finally{writeUnlock();}
	}
	
	/** Locked getter for MPI usage. */
	public static boolean useMpi(){
		readLock();
		try{return Shared.USE_MPI;}
		finally{readUnlock();}
	}
	
	/** Locked setter for MPI keep all mode. */
	public static void setMpiKeepAll(boolean b){
		writeLock();
		try{Shared.MPI_KEEP_ALL=b;}
		finally{writeUnlock();}
	}
	
	/** Locked getter for MPI keep all mode. */
	public static boolean mpiKeepAll(){
		readLock();
		try{return Shared.MPI_KEEP_ALL;}
		finally{readUnlock();}
	}
	
	/** Locked setter for CRISMPI usage. */
	public static void setUseCrismpi(boolean b){
		writeLock();
		try{Shared.USE_CRISMPI=b;}
		finally{writeUnlock();}
	}
	
	/** Locked getter for CRISMPI usage. */
	public static boolean useCrismpi(){
		readLock();
		try{return Shared.USE_CRISMPI;}
		finally{readUnlock();}
	}
	
	/** Locked setter for MPI rank. */
	public static void setMpiRank(int rank){
		writeLock();
		try{Shared.MPI_RANK=rank;}
		finally{writeUnlock();}
	}
	
	/** Locked getter for MPI rank. */
	public static int mpiRank(){
		readLock();
		try{return Shared.MPI_RANK;}
		finally{readUnlock();}
	}
	
	/** Locked setter for number of MPI ranks. */
	public static void setMpiNumRanks(int ranks){
		writeLock();
		try{Shared.MPI_NUM_RANKS=ranks;}
		finally{writeUnlock();}
	}
	
	/** Locked getter for number of MPI ranks. */
	public static int mpiNumRanks(){
		readLock();
		try{return Shared.MPI_NUM_RANKS;}
		finally{readUnlock();}
	}
	
	// File and I/O Configuration
	
	/** Locked setter for FASTA line wrap length. */
	public static void setFastaWrap(int wrap){
		writeLock();
		try{Shared.FASTA_WRAP=wrap;}
		finally{writeUnlock();}
	}
	
	/** Locked getter for FASTA line wrap length. */
	public static int fastaWrap(){
		readLock();
		try{return Shared.FASTA_WRAP;}
		finally{readUnlock();}
	}
	
	/** Locked setter for fake quality score. */
	public static void setFakeQual(byte qual){
		writeLock();
		try{Shared.FAKE_QUAL=qual;}
		finally{writeUnlock();}
	}
	
	/** Locked getter for fake quality score. */
	public static byte fakeQual(){
		readLock();
		try{return Shared.FAKE_QUAL;}
		finally{readUnlock();}
	}
	
	/** Locked setter for file extension fixing. */
	public static void setFixExtensions(boolean b){
		writeLock();
		try{Shared.FIX_EXTENSIONS=b;}
		finally{writeUnlock();}
	}
	
	/** Locked getter for file extension fixing. */
	public static boolean fixExtensions(){
		readLock();
		try{return Shared.FIX_EXTENSIONS;}
		finally{readUnlock();}
	}
	
	/** Locked setter for read comment trimming. */
	public static void setTrimReadComments(boolean b){
		writeLock();
		try{Shared.TRIM_READ_DESCRIPTION=b;}
		finally{writeUnlock();}
	}
	
	/** Locked getter for read comment trimming. */
	public static boolean trimReadComments(){
		readLock();
		try{return Shared.TRIM_READ_DESCRIPTION;}
		finally{readUnlock();}
	}
	
	/** Locked setter for RNAME trimming in SAM reads. */
	public static void setTrimRname(boolean b){
		writeLock();
		try{Shared.TRIM_RNAME=b;}
		finally{writeUnlock();}
	}
	
	/** Locked getter for RNAME trimming in SAM reads. */
	public static boolean trimRname(){
		readLock();
		try{return Shared.TRIM_RNAME;}
		finally{readUnlock();}
	}
	
	/** Locked setter for KMG output formatting. */
	public static void setOutputKmg(boolean b){
		writeLock();
		try{Shared.OUTPUT_KMG=b;}
		finally{writeUnlock();}
	}
	
	/** Locked getter for KMG output formatting. */
	public static boolean outputKmg(){
		readLock();
		try{return Shared.OUTPUT_KMG;}
		finally{readUnlock();}
	}
	
	// Buffer Configuration
	
	/** Locked setter for read buffer length. */
	public static void setBufferLen(int len){
		writeLock();
		try{Shared.setBufferLen(len);}
		finally{writeUnlock();}
	}
	
	/** Locked getter for read buffer length. */
	public static int bufferLen(){
		readLock();
		try{return Shared.bufferLen();}
		finally{readUnlock();}
	}
	
	/** Locked setter for maximum buffer data. */
	public static void setBufferData(int data){
		writeLock();
		try{Shared.setBufferData(data);}
		finally{writeUnlock();}
	}
	
	/** Locked getter for maximum buffer data. */
	public static long bufferData(){
		readLock();
		try{return Shared.bufferData();}
		finally{readUnlock();}
	}
}
