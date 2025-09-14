package assemble;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import shared.KillSwitch;
import structures.ByteBuilder;

/**
 * Thread for exploring connectivity graph between contigs.
 * @author Brian Bushnell
 * @date July 12, 2018
 *
 */
public abstract class AbstractProcessContigThread extends Thread {

	AbstractProcessContigThread(ArrayList<Contig> contigs_, AtomicInteger next_){
		contigs=contigs_;
		next=next_;
	}
	
	@Override
	public void run(){
		processContigs(contigs);
	}

	public final void processContigs(ArrayList<Contig> contigs){
		for(int cnum=next.getAndIncrement(); cnum<contigs.size(); cnum=next.getAndIncrement()){
			Contig c=contigs.get(cnum);
			processContigLeft(c, leftCounts, rightCounts, extraCounts, bb);
			processContigRight(c, leftCounts, rightCounts, extraCounts, bb);
		}
	}

	abstract void processContigLeft(Contig c, int[] leftCounts, int[] rightCounts, int[] extraCounts, ByteBuilder bb);

	abstract void processContigRight(Contig c, int[] leftCounts, int[] rightCounts, int[] extraCounts, ByteBuilder bb);

	final int[] leftCounts=KillSwitch.allocInt1D(4);
	final int[] rightCounts=KillSwitch.allocInt1D(4);
	final int[] extraCounts=KillSwitch.allocInt1D(4);

	final ArrayList<Contig> contigs;
	final AtomicInteger next;

	int lastLength=-1;
	int lastTarget=-1;
	int lastExitCondition=-1;
	int lastOrientation=-1;
	ByteBuilder bb=new ByteBuilder();
	long edgesMadeT=0;

}
