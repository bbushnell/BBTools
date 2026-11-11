package assemble;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import kmer.AbstractKmerTableSet;
import kmer.KmerTableSet;
import shared.Timer;
import ukmer.KmerTableSetU;

/**
 * Coordinates exploration and removal passes that shave dead ends and bubbles
 * from a k-mer graph. Concrete subclasses bind this scheduler to a specific
 * table implementation.
 *
 * @author Brian Bushnell
 * @date Jun 26, 2015
 */
public abstract class Shaver extends ShaveObject{


	/*--------------------------------------------------------------*/
	/*----------------          Factories           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a Shaver with the default pass thresholds.
	 *
	 * @param tables K-mer table set containing the assembly graph
	 * @param threads Number of threads to use for each phase
	 * @return Shaver1 for exact {@link KmerTableSet} instances or Shaver2 for
	 *         exact {@link KmerTableSetU} instances
	 */
	public static final Shaver makeShaver(AbstractKmerTableSet tables, int threads){
		return makeShaver(tables, threads, 1, 1, 1, 1, 3, 100, 100, true, true);
	}

	/**
	 * Creates a Shaver bound to the supplied table implementation. The exact
	 * runtime class selects the implementation, and subclasses are rejected.
	 *
	 * @param tables K-mer table set containing the assembly graph
	 * @param threads Number of threads to use for each phase
	 * @param minCount Minimum k-mer count to consider for removal
	 * @param maxCount Maximum k-mer count to consider for removal
	 * @param minSeed Minimum k-mer count to start exploration from
	 * @param minCountExtend Minimum count for extending paths during exploration
	 * @param branchMult2 Branch multiplier for bubble detection
	 * @param maxLengthToDiscard Maximum length of paths to remove
	 * @param maxDistanceToExplore Maximum distance to explore from seed k-mers
	 * @param removeHair Whether to remove dead-end paths (hairs)
	 * @param removeBubbles Whether to remove bubble structures
	 * @return Shaver1 for exact {@link KmerTableSet} instances or Shaver2 for
	 *         exact {@link KmerTableSetU} instances
	 * @throws RuntimeException If the exact table type is not supported
	 */
	public static final Shaver makeShaver(AbstractKmerTableSet tables, int threads,
			int minCount, int maxCount, int minSeed, int minCountExtend, float branchMult2, int maxLengthToDiscard, int maxDistanceToExplore,
			boolean removeHair, boolean removeBubbles){
		final Class<?> c=tables.getClass();
		if(c==KmerTableSet.class){
			return new Shaver1((KmerTableSet)tables, threads, minCount, maxCount, minSeed, minCountExtend, branchMult2,
					maxLengthToDiscard, maxDistanceToExplore, removeHair, removeBubbles);
		}else if(c==KmerTableSetU.class){
			return new Shaver2((KmerTableSetU)tables, threads, minCount, maxCount, minSeed, minCountExtend, branchMult2,
					maxLengthToDiscard, maxDistanceToExplore, removeHair, removeBubbles);
		}
		throw new RuntimeException("Unhandled class "+c);
	}


	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Stores pass thresholds shared by table-specific subclasses.
	 *
	 * @param tables_ K-mer table set containing the assembly graph
	 * @param threads_ Number of threads to use for each phase
	 * @param minCount_ Minimum k-mer count to consider for removal
	 * @param maxCount_ Maximum k-mer count to consider for removal
	 * @param minSeed_ Minimum k-mer count to start exploration from
	 * @param minCountExtend_ Minimum count for extending paths during exploration
	 * @param branchMult2_ Branch multiplier for bubble detection
	 * @param maxLengthToDiscard_ Maximum length of paths to remove
	 * @param maxDistanceToExplore_ Maximum distance to explore from seed k-mers
	 * @param removeHair_ Whether to remove dead-end paths (hairs)
	 * @param removeBubbles_ Whether to remove bubble structures
	 */
	public Shaver(AbstractKmerTableSet tables_, int threads_,
			int minCount_, int maxCount_, int minSeed_, int minCountExtend_, float branchMult2_, int maxLengthToDiscard_, int maxDistanceToExplore_,
			boolean removeHair_, boolean removeBubbles_){
		threads=threads_;
		minCount=minCount_;
		maxCount=maxCount_;
		minSeed=minSeed_;
		minCountExtend=minCountExtend_;
		branchMult2=branchMult2_;
		maxLengthToDiscard=maxLengthToDiscard_;
		maxDistanceToExplore=maxDistanceToExplore_;
		removeHair=removeHair_;
		removeBubbles=removeBubbles_;

		kbig=tables_.kbig();
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Allocates an exploration worker for this table implementation.
	 * @param id_ Worker identifier for the current exploration phase
	 * @return Configured exploration thread instance
	 */
	protected abstract AbstractExploreThread makeExploreThread(int id_);
	/**
	 * Allocates a removal worker for this table implementation.
	 * @param id_ Worker identifier for the current removal phase
	 * @return Configured shaving thread instance
	 */
	protected abstract AbstractShaveThread makeShaveThread(int id_);

	/**
	 * Performs shaving with updated count parameters.
	 * Convenience method that updates thresholds and calls main shave method.
	 *
	 * @param minCount_ New minimum k-mer count threshold
	 * @param maxCount_ New maximum k-mer count threshold
	 * @return Number of k-mers removed from the graph
	 */
	public final long shave(int minCount_, int maxCount_){
		minCount=minCount_;
		maxCount=maxCount_;
		return shave();
	}

	/**
	 * Removes selected dead-end and bubble paths from the graph in two phases.
	 * The exploration phase initializes fresh ownership storage, marks candidate
	 * kmers and records selection statistics; the removal phase then deletes
	 * claimed kmers. Callers must run shave invocations on the same coordinator
	 * and table serially and must satisfy the table's ownership-initialization
	 * preconditions before this method allocates new ownership storage. A worker
	 * failure is rethrown after every worker in that phase has joined and before
	 * that phase's counters are aggregated. Kmers already selected or removed by
	 * other workers are not rolled back.
	 * @return Number of k-mers removed from the graph
	 * @throws RuntimeException If an exploration or removal worker fails
	 */
	public final long shave(){
		assert(minSeed>=minCount) : "Required: minCount <= minSeed <= maxCount. "+minCount+", "+minSeed+", "+maxCount;//msg was backwards ('mincount >= minSeed >= maxCount'); code asserts minCount<=minSeed<=maxCount [assemble/Shaver#001 FIXED]
		assert(minSeed<=maxCount) : "Required: minCount <= minSeed <= maxCount. "+minCount+", "+minSeed+", "+maxCount;
		assert(removeHair || removeBubbles);

		Timer t=new Timer();

		long kmersTestedTemp=0;
		long deadEndsFoundTemp=0;
		long kmersRemovedTemp=0;
		long bubblesFoundTemp=0;

		tables().initializeOwnership();


		// Aggregate every valid event code, including branch codes 17-19.
		countMatrix=new long[MAX_CODE][MAX_CODE];
		removeMatrix=new long[MAX_CODE][MAX_CODE];

		{
			nextTable.set(0);
			nextVictims.set(0);

			/* Create exploration workers. */
			ArrayList<AbstractExploreThread> alpt=new ArrayList<AbstractExploreThread>(threads);
			for(int i=0; i<threads; i++){alpt.add(makeExploreThread(i));}
			for(AbstractExploreThread pt : alpt){pt.start();}

			/* Join every exploration worker before rethrowing the first failure. */
			Throwable workerFailure=null;
			for(AbstractExploreThread pt : alpt){
				while(pt.getState()!=Thread.State.TERMINATED){
					try{
						pt.join();
					}catch(InterruptedException e){
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				if(workerFailure==null && pt.failure!=null){workerFailure=pt.failure;}
			}
			if(workerFailure!=null){
				throw new RuntimeException(getClass().getSimpleName()+" exploration worker thread failed.",
						workerFailure);
			}

			/* Aggregate only after a clean exploration phase. */
			for(AbstractExploreThread pt : alpt){
				kmersTestedTemp+=pt.kmersTestedT;
				deadEndsFoundTemp+=pt.deadEndsFoundT;
				bubblesFoundTemp+=pt.bubblesFoundT;

				for(int i=0; i<countMatrix.length; i++){
					for(int j=0; j<countMatrix[i].length; j++){
						countMatrix[i][j]+=pt.countMatrixT[i][j];
						removeMatrix[i][j]+=pt.removeMatrixT[i][j];
					}
				}
			}
			kmersTested+=kmersTestedTemp;
			deadEndsFound+=deadEndsFoundTemp;
			bubblesFound+=bubblesFoundTemp;

			t.stop();

			outstream.println("Tested "+kmersTestedTemp+" kmers.");
			outstream.println("Found "+deadEndsFoundTemp+" dead ends.");
			outstream.println("Found "+bubblesFoundTemp+" bubbles.");

			outstream.println("Search time: "+t);
		}

		{
			t.start();

			nextTable.set(0);
			nextVictims.set(0);

			/* Create removal workers. */
			ArrayList<AbstractShaveThread> alpt=new ArrayList<AbstractShaveThread>(threads);
			for(int i=0; i<threads; i++){alpt.add(makeShaveThread(i));}
			for(AbstractShaveThread pt : alpt){pt.start();}

			/* Join every removal worker before rethrowing the first failure. */
			Throwable workerFailure=null;
			for(AbstractShaveThread pt : alpt){
				while(pt.getState()!=Thread.State.TERMINATED){
					try{
						pt.join();
					}catch(InterruptedException e){
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				if(workerFailure==null && pt.failure!=null){workerFailure=pt.failure;}
			}
			if(workerFailure!=null){
				throw new RuntimeException(getClass().getSimpleName()+" removal worker thread failed.",
						workerFailure);
			}

			/* Aggregate only after a clean removal phase. */
			for(AbstractShaveThread pt : alpt){
				kmersRemovedTemp+=pt.kmersRemovedT;
			}

			kmersRemoved+=kmersRemovedTemp;

			outstream.println("Removed "+kmersRemovedTemp+" kmers.");
			t.stop();
			outstream.println("Shave time: "+t);
		}

		if(printEventCounts){
			outstream.println("\nEvent counts:");
			for(int i=0; i<countMatrix.length; i++){
				for(int j=0; j<countMatrix[i].length; j++){
					outstream.print(countMatrix[i][j]+" ");
				}
				outstream.println();
			}
			outstream.println("\nRemoval counts:");
			for(int i=0; i<removeMatrix.length; i++){
				for(int j=0; j<removeMatrix[i].length; j++){
					outstream.print(removeMatrix[i][j]+" ");
				}
				outstream.println();
			}
		}

		return kmersRemovedTemp;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Kmers tested by all completed exploration phases. */
	public long kmersTested=0;
	/** Paths selected as dead ends; concrete searches may include bubble paths. */
	public long deadEndsFound=0;
	/** Bubble paths selected by all completed exploration phases. */
	public long bubblesFound=0;
	/** Kmers deleted by all completed removal phases. */
	public long kmersRemoved=0;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Provides access to the k-mer table set containing the assembly graph.
	 * Implementation varies by concrete subclass.
	 * @return The k-mer table set being processed
	 */
	protected abstract AbstractKmerTableSet tables();
	/** K-mer length supplied by the concrete table set. */
	final int kbig;
	/** Number of exploration and shaving workers to create per phase. */
	final int threads;
	/** Inclusive minimum count eligible for exploration/removal. */
	int minCount;
	/** Inclusive maximum count eligible for exploration/removal. */
	int maxCount;
	/** Minimum count for starting an explored candidate path. */
	final int minSeed;
	/** Minimum count for extending along a candidate path. */
	final int minCountExtend;
	/** Count-ratio cutoff used by bubble classifiers. */
	final float branchMult2;
	/** Maximum removable hair or bubble path length. */
	final int maxLengthToDiscard;
	/** Maximum graph distance explored before a candidate is abandoned. */
	final int maxDistanceToExplore;
	/** Enables dead-end path removal during this pass. */
	final boolean removeHair;
	/** Enables bubble path removal during this pass. */
	final boolean removeBubbles;
	/** Starts exploration from high-count kmers when enabled. */
	static boolean startFromHighCounts=true; //True is much faster, but decreases contiguity.
	/** Enables the fast shaving search mode. */
	static boolean shaveFast=true;
	/** Enables the fastest shaving search mode when compiled true. */
	static final boolean shaveVFast=false; //True is faster, but slightly decreases contiguity.

	/** Explored path counts by sorted minimum and maximum termination codes. */
	private long[][] countMatrix;
	/** Paths selected for removal, counted before final ownership validation. */
	private long[][] removeMatrix;

	/** Atomic counter for coordinating thread access to k-mer tables. */
	final AtomicInteger nextTable=new AtomicInteger(0);

	/** Atomic counter for coordinating access to overflow victim arrays. */
	final AtomicInteger nextVictims=new AtomicInteger(0);

}
