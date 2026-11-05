package ukmer;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import shared.KillSwitch;
import shared.Shared;
import shared.Tools;

/**
 * Initializes or clears ownership storage across distinct, stable long-k-mer tables.
 * Static entry points join their workers and propagate failures to the caller;
 * completed table mutations are not rolled back. The caller must exclude other
 * table mutations and ownership users for the duration of the operation.
 * @author Brian Bushnell
 */
public class OwnershipThread extends Thread{

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Releases ownership storage and waits for all helper workers before returning.
	 * @param tables Array of k-mer tables to clear ownership data from
	 */
	public static void clear(AbstractKmerTableU[] tables){
		process(tables, CLEAR);
	}

	/**
	 * Initializes ownership storage; each table must satisfy its initialization preconditions.
	 * @param tables Array of k-mer tables to initialize ownership data for
	 */
	public static void initialize(AbstractKmerTableU[] tables){
		process(tables, INITIALIZE);
	}

	/**
	 * Processes zero or one table directly, otherwise starts up to Shared.threads()
	 * workers. Joins all workers before rethrowing the first failure found in worker
	 * order. Interrupted joins are logged and retried; successful mutations remain.
	 *
	 * @param tables Array of k-mer tables to process
	 * @param mode Operation mode (INITIALIZE or CLEAR)
	 */
	private static void process(AbstractKmerTableU[] tables, int mode){
		if(tables.length<2){
			if(mode==INITIALIZE){
				for(AbstractKmerTableU akt : tables){akt.initializeOwnership();}
			}else if(mode==CLEAR){
				for(AbstractKmerTableU akt : tables){akt.clearOwnership();}
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

		Throwable failure=null;
		for(OwnershipThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try{
					pt.join();
				}catch(InterruptedException e){
					//Preserve the existing policy of waiting for all owned workers.
					e.printStackTrace();
				}
			}
			if(failure==null){failure=pt.failure;}
		}
		rethrowFailure(failure);
	}

	/** Reports failure only after all workers finish; earlier table mutations are retained. */
	private static void rethrowFailure(final Throwable failure){
		if(failure==null){return;}
		if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
		if(failure instanceof Error){throw (Error)failure;}
		throw new RuntimeException("Long-k-mer ownership worker failed.", failure);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Borrows the tables and shared work counter without starting this worker.
	 * Use initialize/clear for coordinated failure reporting; directly joining a
	 * constructed worker does not itself rethrow its captured failure.
	 * @param tables_ Array of k-mer tables to process
	 * @param mode_ Operation mode (INITIALIZE or CLEAR)
	 * @param next_ Atomic counter for work distribution across threads
	 */
	public OwnershipThread(AbstractKmerTableU[] tables_, int mode_, AtomicInteger next_){
		tables=tables_;
		mode=mode_;
		next=next_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Claims table indices until exhausted or failed; publishes failure to the coordinator. */
	@Override
	public void run(){
		try{
			for(int i=next.getAndIncrement(); i<tables.length; i=next.getAndIncrement()){
				if(mode==INITIALIZE){
					tables[i].initializeOwnership();
				}else if(mode==CLEAR){
					tables[i].clearOwnership();
				}else{
					KillSwitch.kill("Bad mode: "+mode);
				}
			}
		}catch(Throwable t){failure=t;}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Array of k-mer tables to process */
	private final AbstractKmerTableU[] tables;
	/** Atomic counter for thread-safe work distribution across tables */
	private final AtomicInteger next;
	/** Operation mode: INITIALIZE (0) or CLEAR (1) */
	private final int mode;
	/** Exception/error made visible to the coordinator after this worker finishes. */
	private volatile Throwable failure;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Mode constant for initializing k-mer table ownership data */
	public static final int INITIALIZE=0;
	/** Mode constant for clearing k-mer table ownership data */
	public static final int CLEAR=1;

}
