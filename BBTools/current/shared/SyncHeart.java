package shared;

import java.util.concurrent.locks.ReentrantReadWriteLock;

/**
 * Thread-safe synchronization utilities for mutable configuration values in BBTools.
 * Provides safe concurrent access to shared state with optimized read performance.
 * Uses ReadWriteLock to allow multiple concurrent readers while ensuring exclusive writes.
 * 
 * @author Brian Bushnell
 * @author Isla Winglet  
 * @date June 18, 2025
 */
public class SyncHeart {
	
	/** ReadWriteLock for thread-safe access to configuration values */
	private static final ReentrantReadWriteLock configLock = new ReentrantReadWriteLock();
	
	/** Acquires read lock for safe concurrent read access */
	private static final void readLock() {
		configLock.readLock().lock();
	}
	
	/** Releases read lock */
	private static final void readUnlock() {
		configLock.readLock().unlock();
	}
	
	/** Acquires write lock for exclusive write access */
	private static final void writeLock() {
		configLock.writeLock().lock();
	}
	
	/** Releases write lock */
	private static final void writeUnlock() {
		configLock.writeLock().unlock();
	}
	
	// Threading and Performance Configuration
	
	/** Thread-safe setter for amino acid input mode */
	public static void setAminoIn(boolean b) {
		writeLock();
		Shared.AMINO_IN = b;
		writeUnlock();
	}
	
	/** Thread-safe getter for amino acid input mode */
	public static boolean aminoIn() {
		readLock();
		boolean result = Shared.AMINO_IN;
		readUnlock();
		return result;
	}
	
	/** Thread-safe setter for low memory mode */
	public static void setLowMemory(boolean b) {
		writeLock();
		Shared.LOW_MEMORY = b;
		writeUnlock();
	}
	
	/** Thread-safe getter for low memory mode */
	public static boolean lowMemory() {
		readLock();
		boolean result = Shared.LOW_MEMORY;
		readUnlock();
		return result;
	}
	
	/** Thread-safe setter for garbage collection before memory printing */
	public static void setGcBeforePrintMemory(boolean b) {
		writeLock();
		Shared.GC_BEFORE_PRINT_MEMORY = b;
		writeUnlock();
	}
	
	/** Thread-safe getter for garbage collection before memory printing */
	public static boolean gcBeforePrintMemory() {
		readLock();
		boolean result = Shared.GC_BEFORE_PRINT_MEMORY;
		readUnlock();
		return result;
	}
	
	/** Thread-safe setter for parallel sort usage */
	public static void setParallelSort(boolean b) {
		writeLock();
		Shared.parallelSort = b;
		writeUnlock();
	}
	
	/** Thread-safe getter for parallel sort usage */
	public static boolean parallelSort() {
		readLock();
		boolean result = Shared.parallelSort;
		readUnlock();
		return result;
	}
	
	/** Thread-safe setter for SIMD optimization usage */
	public static void setSimd(boolean b) {
		writeLock();
		Shared.SIMD = b;
		writeUnlock();
	}
	
	/** Thread-safe getter for SIMD optimization usage */
	public static boolean simd() {
		readLock();
		boolean result = Shared.SIMD;
		readUnlock();
		return result;
	}
	
	// MPI Configuration
	
	/** Thread-safe setter for MPI usage */
	public static void setUseMpi(boolean b) {
		writeLock();
		Shared.USE_MPI = b;
		writeUnlock();
	}
	
	/** Thread-safe getter for MPI usage */
	public static boolean useMpi() {
		readLock();
		boolean result = Shared.USE_MPI;
		readUnlock();
		return result;
	}
	
	/** Thread-safe setter for MPI keep all mode */
	public static void setMpiKeepAll(boolean b) {
		writeLock();
		Shared.MPI_KEEP_ALL = b;
		writeUnlock();
	}
	
	/** Thread-safe getter for MPI keep all mode */
	public static boolean mpiKeepAll() {
		readLock();
		boolean result = Shared.MPI_KEEP_ALL;
		readUnlock();
		return result;
	}
	
	/** Thread-safe setter for CRISMPI usage */
	public static void setUseCrismpi(boolean b) {
		writeLock();
		Shared.USE_CRISMPI = b;
		writeUnlock();
	}
	
	/** Thread-safe getter for CRISMPI usage */
	public static boolean useCrismpi() {
		readLock();
		boolean result = Shared.USE_CRISMPI;
		readUnlock();
		return result;
	}
	
	/** Thread-safe setter for MPI rank */
	public static void setMpiRank(int rank) {
		writeLock();
		Shared.MPI_RANK = rank;
		writeUnlock();
	}
	
	/** Thread-safe getter for MPI rank */
	public static int mpiRank() {
		readLock();
		int result = Shared.MPI_RANK;
		readUnlock();
		return result;
	}
	
	/** Thread-safe setter for number of MPI ranks */
	public static void setMpiNumRanks(int ranks) {
		writeLock();
		Shared.MPI_NUM_RANKS = ranks;
		writeUnlock();
	}
	
	/** Thread-safe getter for number of MPI ranks */
	public static int mpiNumRanks() {
		readLock();
		int result = Shared.MPI_NUM_RANKS;
		readUnlock();
		return result;
	}
	
	// File and I/O Configuration
	
	/** Thread-safe setter for FASTA line wrap length */
	public static void setFastaWrap(int wrap) {
		writeLock();
		Shared.FASTA_WRAP = wrap;
		writeUnlock();
	}
	
	/** Thread-safe getter for FASTA line wrap length */
	public static int fastaWrap() {
		readLock();
		int result = Shared.FASTA_WRAP;
		readUnlock();
		return result;
	}
	
	/** Thread-safe setter for fake quality score */
	public static void setFakeQual(byte qual) {
		writeLock();
		Shared.FAKE_QUAL = qual;
		writeUnlock();
	}
	
	/** Thread-safe getter for fake quality score */
	public static byte fakeQual() {
		readLock();
		byte result = Shared.FAKE_QUAL;
		readUnlock();
		return result;
	}
	
	/** Thread-safe setter for file extension fixing */
	public static void setFixExtensions(boolean b) {
		writeLock();
		Shared.FIX_EXTENSIONS = b;
		writeUnlock();
	}
	
	/** Thread-safe getter for file extension fixing */
	public static boolean fixExtensions() {
		readLock();
		boolean result = Shared.FIX_EXTENSIONS;
		readUnlock();
		return result;
	}
	
	/** Thread-safe setter for read comment trimming */
	public static void setTrimReadComments(boolean b) {
		writeLock();
		Shared.TRIM_READ_COMMENTS = b;
		writeUnlock();
	}
	
	/** Thread-safe getter for read comment trimming */
	public static boolean trimReadComments() {
		readLock();
		boolean result = Shared.TRIM_READ_COMMENTS;
		readUnlock();
		return result;
	}
	
	/** Thread-safe setter for RNAME trimming in SAM reads */
	public static void setTrimRname(boolean b) {
		writeLock();
		Shared.TRIM_RNAME = b;
		writeUnlock();
	}
	
	/** Thread-safe getter for RNAME trimming in SAM reads */
	public static boolean trimRname() {
		readLock();
		boolean result = Shared.TRIM_RNAME;
		readUnlock();
		return result;
	}
	
	/** Thread-safe setter for KMG output formatting */
	public static void setOutputKmg(boolean b) {
		writeLock();
		Shared.OUTPUT_KMG = b;
		writeUnlock();
	}
	
	/** Thread-safe getter for KMG output formatting */
	public static boolean outputKmg() {
		readLock();
		boolean result = Shared.OUTPUT_KMG;
		readUnlock();
		return result;
	}
	
	// Buffer Configuration
	
	/** Thread-safe setter for read buffer length */
	public static void setBufferLen(int len) {
		writeLock();
		Shared.setBufferLen(len);
		writeUnlock();
	}
	
	/** Thread-safe getter for read buffer length */
	public static int bufferLen() {
		readLock();
		int result = Shared.bufferLen();
		readUnlock();
		return result;
	}
	
	/** Thread-safe setter for maximum buffer data */
	public static void setBufferData(long data) {
		writeLock();
		Shared.setBufferData(data);
		writeUnlock();
	}
	
	/** Thread-safe getter for maximum buffer data */
	public static long bufferData() {
		readLock();
		long result = Shared.bufferData();
		readUnlock();
		return result;
	}
}