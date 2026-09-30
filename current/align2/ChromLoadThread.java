package align2;

import dna.ChromosomeArray;
import shared.Shared;

/** Loads chromosome arrays with a shared concurrency limit.
 * loadAll joins every loader it starts, including after a failed read or launch,
 * and propagates failure instead of waiting for a result slot that cannot fill.
 * Callers must own the target slots; Data.loadChromosomes serializes access using
 * CHROMLOCKS. A failed batch can leave successfully loaded slots populated.
 * @author Brian Bushnell, Collei
 * @date Dec 31, 2012
 */
public class ChromLoadThread extends Thread{

	public static void main(String[] args){}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructs an unstarted loader. load reserves capacity before starting;
	 * a directly started instance acquires its slot at the beginning of run. */
	public ChromLoadThread(final String fname_, final int id_, final ChromosomeArray[] r_){
		fname=fname_; id=id_; array=r_;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Reserves capacity and starts a loader for an empty target slot.
	 * Under assertions, an occupied slot is a caller error. With assertions
	 * disabled, the historical occupied-slot behavior is to return null. */
	public static ChromLoadThread load(final String fname, final int id, final ChromosomeArray[] r){
		assert(r[id]==null) : "Chromosome loader requires an empty target slot: "+id;
		if(r[id]!=null){return null;}
		final ChromLoadThread loader=new ChromLoadThread(fname, id, r);
		increment(1);
		loader.reserved=true;
		try{loader.start();}
		catch(RuntimeException|Error problem){
			loader.reserved=false;
			increment(-1);
			throw problem;
		}
		return loader;
	}

	/** Loads the inclusive range using the last '#' in pattern as a placeholder.
	 * The last chromosome is read synchronously. Every started worker is joined
	 * before returning or throwing; callers never receive a silently incomplete batch.
	 * @param r Destination, or null to allocate max+1 slots
	 * @return Destination containing all requested chromosomes
	 * @throws RuntimeException If a read, launch, or interruption prevents completion */
	public static ChromosomeArray[] loadAll(final String pattern, final int min, final int max, ChromosomeArray[] r){
		if(r==null){r=new ChromosomeArray[max+1];}
		assert(r.length>max) : "Chromosome destination must include maximum index "+max+"; length="+r.length;
		final int pound=pattern.lastIndexOf('#');
		if(pound<0){throw new IllegalArgumentException("Chromosome filename pattern requires '#': "+pattern);}
		final String a=pattern.substring(0, pound), b=pattern.substring(pound+1);
		final ChromLoadThread[] loaders=new ChromLoadThread[max];
		Throwable failure=null;
		try{
			for(int i=min; i<max; i++){loaders[i]=load(a+i+b, i, r);}
			if(max>=min){
				increment(1);
				try{r[max]=ChromosomeArray.read(a+max+b);}
				finally{increment(-1);}
			}
		}catch(RuntimeException|Error problem){failure=problem;}

		// A failed synchronous read or launch must not abandon already started loaders.
		boolean interrupted=false;
		for(int i=min; i<max; i++){
			final ChromLoadThread loader=loaders[i];
			if(loader==null){continue;}
			while(loader.isAlive()){
				try{loader.join();}
				catch(InterruptedException problem){
					interrupted=true;
					if(failure==null){failure=new RuntimeException("Interrupted while joining chromosome loader "+i, problem);}
				}
			}
			if(loader.failure!=null){
				if(failure==null){failure=loader.failure;}
				else if(failure!=loader.failure){failure.addSuppressed(loader.failure);}
			}
		}
		if(interrupted){Thread.currentThread().interrupt();}
		if(failure!=null){rethrow(failure);}
		for(int i=min; i<=max; i++){
			if(r[i]==null){throw new IllegalStateException("Chromosome loader completed without a result for "+a+i+b);}
		}
		return r;
	}

	/** Preserves unchecked failure type when passing a worker failure to its caller. */
	private static void rethrow(final Throwable problem){
		if(problem instanceof Error){throw (Error)problem;}
		if(problem instanceof RuntimeException){throw (RuntimeException)problem;}
		throw new RuntimeException(problem);
	}

	/** Acquires/releases one live-load slot. Waiters all recheck the same capacity.
	 * Interruption before acquisition leaves the counter unchanged. */
	private static int increment(final int delta){
		synchronized(lock){
			if(delta>0){
				if(MAX_CONCURRENT<1){throw new IllegalStateException("Chromosome load concurrency must be positive: "+MAX_CONCURRENT);}
				while(lock[0]>=MAX_CONCURRENT){
					try{lock.wait();}
					catch(InterruptedException problem){
						Thread.currentThread().interrupt();
						throw new RuntimeException("Interrupted while acquiring a chromosome load slot", problem);
					}
				}
			}
			lock[0]+=delta;
			assert(lock[0]>=0) : "Each chromosome loader releases exactly one acquired slot; count="+lock[0];
			if(delta<0){lock.notifyAll();}
			return lock[0];
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Publishes a chromosome or records/rethrows its failure; always releases its slot. */
	@Override
	public void run(){
		if(!reserved){increment(1); reserved=true;}
		try{array[id]=ChromosomeArray.read(fname);}
		catch(RuntimeException|Error problem){failure=problem; throw problem;}
		finally{reserved=false; increment(-1);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int id;
	private final String fname;
	private final ChromosomeArray[] array;
	/** Published by the worker and inspected only after join. */
	private Throwable failure;
	/** Set before factory start, or acquired by run for direct construction. */
	private boolean reserved;
	/** Shared live-load count; reads and writes require synchronization on this array. */
	public static final int[] lock=new int[1];
	/** Set before loading; changing it during a batch is unsupported. */
	public static int MAX_CONCURRENT=Shared.threads();
}
