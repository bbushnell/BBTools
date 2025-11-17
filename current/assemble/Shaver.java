package assemble;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import kmer.AbstractKmerTableSet;
import kmer.KmerTableSet;
import shared.Timer;
import ukmer.KmerTableSetU;

/**
 * Designed for removal of dead ends (aka hairs).
 * @author Brian Bushnell
 * @date Jun 26, 2015
 *
 */
public abstract class Shaver extends ShaveObject {
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Factory          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Factory method to create a Shaver with default parameters.
	 * Uses conservative settings suitable for most assembly scenarios.
	 *
	 * @param tables K-mer table set containing the assembly graph
	 * @param threads Number of threads to use for parallel processing
	 * @return Configured Shaver instance (Shaver1 or Shaver2 based on table type)
	 */
	public static final Shaver makeShaver(AbstractKmerTableSet tables, int threads){
		return makeShaver(tables, threads, 1, 1, 1, 1, 3, 100, 100, true, true);
	}
	
	/**
	 * Factory method to create a Shaver with custom parameters.
	 * Returns appropriate concrete implementation based on table type.
	 *
	 * @param tables K-mer table set containing the assembly graph
	 * @param threads Number of threads for parallel processing
	 * @param minCount Minimum k-mer count to consider for removal
	 * @param maxCount Maximum k-mer count to consider for removal
	 * @param minSeed Minimum k-mer count to start exploration from
	 * @param minCountExtend Minimum count for extending paths during exploration
	 * @param branchMult2 Branch multiplier for bubble detection
	 * @param maxLengthToDiscard Maximum length of paths to remove
	 * @param maxDistanceToExplore Maximum distance to explore from seed k-mers
	 * @param removeHair Whether to remove dead-end paths (hairs)
	 * @param removeBubbles Whether to remove bubble structures
	 * @return Configured Shaver instance matching the table type
	 * @throws RuntimeException If table type is not supported
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
	 * Constructor for Shaver instances.
	 * Initializes all configuration parameters and extracts k-mer size from tables.
	 *
	 * @param tables_ K-mer table set containing the assembly graph
	 * @param threads_ Number of threads for parallel processing
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
	 * Creates exploration threads for finding dead ends and bubbles.
	 * Implementation varies by concrete subclass based on k-mer table type.
	 * @param id_ Thread identifier for the exploration worker
	 * @return Configured exploration thread instance
	 */
	abstract AbstractExploreThread makeExploreThread(int id_);
	/**
	 * Creates shaving threads for removing identified graph artifacts.
	 * Implementation varies by concrete subclass based on k-mer table type.
	 * @param id_ Thread identifier for the shaving worker
	 * @return Configured shaving thread instance
	 */
	abstract AbstractShaveThread makeShaveThread(int id_);
	
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
	 * Main shaving method that removes dead ends and bubbles from k-mer graph.
	 * Uses two-phase approach: exploration to identify artifacts, then removal.
	 * Creates worker threads for parallel processing and collects statistics.
	 * @return Number of k-mers removed from the graph
	 */
	public final long shave(){
		assert(minSeed>=minCount) : "Required: mincount >= minSeed >= maxCount. "+minCount+", "+minSeed+", "+maxCount;
		assert(minSeed<=maxCount) : "Required: mincount >= minSeed >= maxCount. "+minCount+", "+minSeed+", "+maxCount;
		assert(removeHair || removeBubbles);
		
		Timer t=new Timer();
		
		long kmersTestedTemp=0;
		long deadEndsFoundTemp=0;
		long kmersRemovedTemp=0;
		long bubblesFoundTemp=0;
		
		tables().initializeOwnership();
		

		countMatrix=new long[8][8];
		removeMatrix=new long[8][8];
		
		{
			nextTable.set(0);
			nextVictims.set(0);
			
			/* Create Explorethreads */
			ArrayList<AbstractExploreThread> alpt=new ArrayList<AbstractExploreThread>(threads);
			for(int i=0; i<threads; i++){alpt.add(makeExploreThread(i));}
			for(AbstractExploreThread pt : alpt){pt.start();}

			/* Wait for threads to die, and gather statistics */
			for(AbstractExploreThread pt : alpt){
				while(pt.getState()!=Thread.State.TERMINATED){
					try {
						pt.join();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}

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
			
			/* Create Shavethreads */
			ArrayList<AbstractShaveThread> alpt=new ArrayList<AbstractShaveThread>(threads);
			for(int i=0; i<threads; i++){alpt.add(makeShaveThread(i));}
			for(AbstractShaveThread pt : alpt){pt.start();}
			
			for(AbstractShaveThread pt : alpt){
				while(pt.getState()!=Thread.State.TERMINATED){
					try {
						pt.join();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				
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
	
	/** Cumulative count of k-mers tested across all shaving operations */
	public long kmersTested=0;
	/** Cumulative count of dead ends (hairs) found across all operations */
	public long deadEndsFound=0;
	/** Cumulative count of bubble structures found across all operations */
	public long bubblesFound=0;
	/** Cumulative count of k-mers removed across all shaving operations */
	public long kmersRemoved=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Provides access to the k-mer table set containing the assembly graph.
	 * Implementation varies by concrete subclass.
	 * @return The k-mer table set being processed
	 */
	abstract AbstractKmerTableSet tables();
	/** K-mer size used in the assembly graph */
	final int kbig;
	/** Number of threads to use for parallel processing */
	final int threads;
	/** Minimum k-mer count to consider for removal operations */
	int minCount;
	/** Maximum k-mer count to consider for removal operations */
	int maxCount;
	/** Minimum k-mer count required to start path exploration */
	final int minSeed;
	/** Minimum k-mer count required to extend paths during exploration */
	final int minCountExtend;
	/** Branch multiplier threshold for bubble detection algorithms */
	final float branchMult2;
	/** Maximum length in bases of paths that can be removed */
	final int maxLengthToDiscard;
	/** Maximum distance to explore outward from seed k-mers */
	final int maxDistanceToExplore;
	/** Whether to remove dead-end paths (hairs) from the graph */
	final boolean removeHair;
	/** Whether to remove bubble structures from the graph */
	final boolean removeBubbles;
	/**
	 * Whether to start exploration from high-count k-mers for better performance
	 */
	static boolean startFromHighCounts=true; //True is much faster, but decreases contiguity.
	/** Whether to use fast shaving algorithm optimizations */
	static boolean shaveFast=true;
	/** Whether to use very fast shaving with potential contiguity trade-offs */
	static final boolean shaveVFast=false; //True is faster, but slightly decreases contiguity.

	/**
	 * Matrix tracking counts of different types of graph events during exploration
	 */
	private long[][] countMatrix;
	/** Matrix tracking removal counts for different types of graph artifacts */
	private long[][] removeMatrix;
	
	/** For controlling access to tables */
	final AtomicInteger nextTable=new AtomicInteger(0);
	
	/** For controlling access to victim buffers */
	final AtomicInteger nextVictims=new AtomicInteger(0);
	
}
