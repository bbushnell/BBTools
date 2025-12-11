package assemble;

import java.util.ArrayList;

import shared.KillSwitch;
import stream.ConcurrentReadInputStream;
import structures.ByteBuilder;
import structures.LongList;

/**
 * Abstract base class providing essential infrastructure for parallel De Bruijn
 * graph construction in Tadpole assembler. Implements thread-safe k-mer
 * processing and contig building operations with shared data structures,
 * performance monitoring, and resource management for multithreaded assembly
 * workers.
 *
 * @author Brian Bushnell
 * @date Jul 18, 2015
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
	
	/**
	 * Array of concurrent input streams enabling thread-safe read access across workers
	 */
	final ConcurrentReadInputStream[] crisa;
	
	final int mode;
	int minCountSeedCurrent;

	final int[] leftCounts=KillSwitch.allocInt1D(4);
	final int[] rightCounts=KillSwitch.allocInt1D(4);
	final ByteBuilder builderT=new ByteBuilder();
//	final Contig tempContig=new Contig(null);
	
	final LongList insertSizes=new LongList();
	
	ArrayList<Contig> contigs=new ArrayList<Contig>();
	
	long readsInT=0;
	long basesInT=0;
	long lowqReadsT=0;
	long lowqBasesT=0;
	final int id;
	
}
