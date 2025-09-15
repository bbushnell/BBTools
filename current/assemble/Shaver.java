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

	public static final Shaver makeShaver(AbstractKmerTableSet tables, int threads){
		return makeShaver(tables, threads, 1, 1, 1, 1, 3, 100, 100, true, true);
	}
	
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

	abstract AbstractExploreThread makeExploreThread(int id_);
	abstract AbstractShaveThread makeShaveThread(int id_);
	
	public final long shave(int minCount_, int maxCount_){
		minCount=minCount_;
		maxCount=maxCount_;
		return shave();
	}
	
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
	
	public long kmersTested=0;
	public long deadEndsFound=0;
	public long bubblesFound=0;
	public long kmersRemoved=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	abstract AbstractKmerTableSet tables();
	final int kbig;
	final int threads;
	int minCount;
	int maxCount;
	final int minSeed;
	final int minCountExtend;
	final float branchMult2;
	final int maxLengthToDiscard;
	final int maxDistanceToExplore;
	final boolean removeHair;
	final boolean removeBubbles;
	static boolean startFromHighCounts=true; //True is much faster, but decreases contiguity.
	static boolean shaveFast=true;
	static final boolean shaveVFast=false; //True is faster, but slightly decreases contiguity.

	private long[][] countMatrix;
	private long[][] removeMatrix;
	
	/** For controlling access to tables */
	final AtomicInteger nextTable=new AtomicInteger(0);
	
	/** For controlling access to victim buffers */
	final AtomicInteger nextVictims=new AtomicInteger(0);
	
}
