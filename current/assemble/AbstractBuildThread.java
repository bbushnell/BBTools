package assemble;

import java.util.ArrayList;

import shared.KillSwitch;
import stream.ConcurrentReadInputStream;
import structures.ByteBuilder;
import structures.LongList;

/**
 * @author Brian Bushnell
 * @date Jul 18, 2015
 *
 */
abstract class AbstractBuildThread extends Thread {
	
	/**
	 * Constructs a new abstract build thread with specified identifier, mode,
	 * and input streams. Initializes thread state for concurrent k-mer processing
	 * and contig assembly operations.
	 *
	 * @param id_ Unique thread identifier for coordination and debugging
	 * @param mode_ Assembly mode configuration encoding processing flags
	 * @param crisa_ Array of concurrent read input streams for thread-safe access
	 */
	public AbstractBuildThread(int id_, int mode_, ConcurrentReadInputStream[] crisa_){
		id=id_;
		crisa=crisa_;
		mode=mode_;
	}
	
	/** Input read stream */
	final ConcurrentReadInputStream[] crisa;
	
	/**
	 * Assembly mode configuration encoding processing flags and behavior settings
	 */
	final int mode;
	/**
	 * Current minimum k-mer coverage threshold for seed selection during assembly
	 */
	int minCountSeedCurrent;

	/**
	 * Memory-tracked array for left-side nucleotide frequency counting (A,C,G,T)
	 */
	final int[] leftCounts=KillSwitch.allocInt1D(4);
	/**
	 * Memory-tracked array for right-side nucleotide frequency counting (A,C,G,T)
	 */
	final int[] rightCounts=KillSwitch.allocInt1D(4);
	/**
	 * Thread-local sequence builder for dynamic contig construction with automatic capacity management
	 */
	final ByteBuilder builderT=new ByteBuilder();
//	final Contig tempContig=new Contig(null);
	
	/**
	 * Dynamic collection tracking insert size statistics during paired-read processing
	 */
	final LongList insertSizes=new LongList();
	
	/** Thread-local contig storage for assembled sequences before aggregation */
	ArrayList<Contig> contigs=new ArrayList<Contig>();
	
	/** Thread-local counter tracking total reads processed by this worker */
	long readsInT=0;
	/** Thread-local counter for total nucleotides processed */
	long basesInT=0;
	/** Counter for reads failing quality thresholds */
	long lowqReadsT=0;
	/** Counter for individual bases below quality cutoffs */
	long lowqBasesT=0;
	/**
	 * Unique thread identifier for coordination, debugging, and performance tracking
	 */
	final int id;
	
}
