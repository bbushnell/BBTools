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
	
	public AbstractBuildThread(int id_, int mode_, ConcurrentReadInputStream[] crisa_){
		id=id_;
		crisa=crisa_;
		mode=mode_;
	}
	
	/** Input read stream */
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
