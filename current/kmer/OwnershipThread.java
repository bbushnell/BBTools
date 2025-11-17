package kmer;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import shared.KillSwitch;
import shared.Shared;
import shared.Tools;

/**
 * Thread for parallel processing of ownership operations on k-mer tables.
 * Coordinates initialization and clearing of ownership data structures
 * across multiple AbstractKmerTable instances using worker threads.
 *
 * @author Brian Bushnell
 * @date December 5, 2024
 */
public class OwnershipThread extends Thread {
	
	/** Clears ownership data structures from all tables using parallel processing.
	 * @param tables Array of k-mer tables to clear ownership data from */
	public static void clear(AbstractKmerTable[] tables){
		process(tables, CLEAR);
	}
	
	/** Initializes ownership data structures in all tables using parallel processing.
	 * @param tables Array of k-mer tables to initialize ownership data for */
	public static void initialize(AbstractKmerTable[] tables){
		process(tables, INITIALIZE);
	}
	
	/**
	 * Core processing method that coordinates parallel ownership operations.
	 * Creates worker threads and manages their execution across table array.
	 * For single table arrays, processes directly without threading overhead.
	 *
	 * @param tables Array of k-mer tables to process
	 * @param mode Operation mode: INITIALIZE or CLEAR
	 */
	private static void process(AbstractKmerTable[] tables, int mode){
		if(tables.length<2){
			if(mode==INITIALIZE){
				for(AbstractKmerTable akt : tables){akt.initializeOwnership();}
			}else if(mode==CLEAR){
				for(AbstractKmerTable akt : tables){akt.clearOwnership();}
			}else{
				KillSwitch.kill("Bad mode: "+mode);
			}
			return;
		}
		final int threads=Tools.min(Shared.threads(), tables.length);
		final AtomicInteger next=new AtomicInteger(0);
		ArrayList<OwnershipThread> alpt=new ArrayList<OwnershipThread>(threads);
		for(int i=0; i<threads; i++){alpt.add(new OwnershipThread(tables, mode, next));}
		for(OwnershipThread pt : alpt){pt.start();}
		
		for(OwnershipThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	
	/**
	 * Constructs a worker thread for processing ownership operations.
	 * @param tables_ Array of k-mer tables to process
	 * @param mode_ Operation mode: INITIALIZE or CLEAR
	 * @param next_ Atomic counter for work distribution among threads
	 */
	public OwnershipThread(AbstractKmerTable[] tables_, int mode_, AtomicInteger next_){
		tables=tables_;
		mode=mode_;
		next=next_;
	}
	
	@Override
	public void run(){
		for(int i=next.getAndIncrement(); i<tables.length; i=next.getAndIncrement()){
			if(mode==INITIALIZE){
				tables[i].initializeOwnership();
			}else if(mode==CLEAR){
				tables[i].clearOwnership();
			}else{
				KillSwitch.kill("Bad mode: "+mode);
			}
		}
	}
	
	/** Array of k-mer tables this thread will process */
	private final AbstractKmerTable[] tables;
	/** Atomic counter for distributing work among threads */
	private final AtomicInteger next;
	/** Operation mode: INITIALIZE or CLEAR */
	private final int mode;

	/** Mode constant for initializing ownership data structures */
	public static final int INITIALIZE=0;
	/** Mode constant for clearing ownership data structures */
	public static final int CLEAR=1;
	
}
