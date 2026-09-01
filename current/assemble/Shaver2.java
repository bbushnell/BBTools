package assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicIntegerArray;

import dna.AminoAcid;
import kmer.AbstractKmerTableSet;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import ukmer.AbstractKmerTableU;
import ukmer.HashArrayU1D;
import ukmer.HashForestU;
import ukmer.Kmer;
import ukmer.KmerNodeU;
import ukmer.KmerTableSetU;

/**
 * Designed for removal of dead ends (aka hairs).
 * @author Brian Bushnell
 * @date Jun 26, 2015
 *
 */
public class Shaver2 extends Shaver {
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Creates a Shaver2 with default conservative trimming parameters.
	 * Uses minimal thresholds with both hair and bubble removal enabled.
	 * @param tables_ K-mer hash table set for unlimited k-mer sizes
	 * @param threads_ Number of worker threads for parallel processing
	 */
	public Shaver2(KmerTableSetU tables_, int threads_){
		this(tables_, threads_, 1, 1, 1, 1, 3, 100, 100, true, true);
	}
	
	/**
	 * Creates a Shaver2 with comprehensive parameter configuration.
	 * Enables fine-tuned control over trimming behavior and thresholds.
	 *
	 * @param tables_ K-mer hash table set for unlimited k-mer sizes
	 * @param threads_ Number of worker threads for parallel processing
	 * @param minCount_ Minimum k-mer count for inclusion in operations
	 * @param maxCount_ Maximum k-mer count threshold
	 * @param minSeed_ Minimum coverage for seed k-mer selection
	 * @param minCountExtend_ Coverage threshold for extension operations
	 * @param branchMult2_ Coverage ratio threshold for branch removal
	 * @param maxLengthToDiscard_ Maximum path length eligible for removal
	 * @param maxDistanceToExplore_ Search radius limit during exploration
	 * @param removeHair_ Enable removal of dead-end paths
	 * @param removeBubbles_ Enable removal of bubble structures
	 */
	public Shaver2(KmerTableSetU tables_, int threads_,
			int minCount_, int maxCount_, int minSeed_, int minCountExtend_, float branchMult2_, int maxLengthToDiscard_, int maxDistanceToExplore_,
			boolean removeHair_, boolean removeBubbles_){
		super(tables_, threads_, minCount_, maxCount_, minSeed_, minCountExtend_, branchMult2_, maxLengthToDiscard_, maxDistanceToExplore_, removeHair_, removeBubbles_);
		tables=tables_;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	final AbstractExploreThread makeExploreThread(int id_){return new ExploreThread(id_);}
	@Override
	final AbstractShaveThread makeShaveThread(int id_){return new ShaveThread(id_);}

	/** Washes a compact table from contig-tip seeds without enumerating stored keys. */
	public long shaveFromContigs(final ArrayList<Contig> contigs, final boolean crossK){
		if(!tables.hashOnly()){throw new IllegalStateException("Tip-seeded washing is specific to compact kmer tables.");}
		assert(minSeed>=minCount && minSeed<=maxCount);
		final Timer timer=new Timer();
		tables.initializeOwnership();
		final AtomicInteger nextTip=new AtomicInteger(0);
		final ArrayList<CompactExploreThread> explorers=new ArrayList<CompactExploreThread>(threads);
		for(int i=0; i<threads; i++){explorers.add(new CompactExploreThread(contigs, nextTip, crossK));}
		for(CompactExploreThread thread : explorers){thread.start();}
		long tested=0, paths=0;
		Throwable failure=null;
		for(CompactExploreThread thread : explorers){
			while(thread.getState()!=Thread.State.TERMINATED){
				try{thread.join();}
				catch(InterruptedException e){Thread.currentThread().interrupt(); throw new RuntimeException(e);}
			}
			tested+=thread.tested;
			paths+=thread.paths;
			if(failure==null){failure=thread.failure;}
		}
		if(failure!=null){
			tables.clearOwnership();
			throw new RuntimeException("Compact wash exploration failed.", failure);
		}
		timer.stop();
		outstream.println("Tested "+tested+" tip-seeded paths.");
		outstream.println("Found "+paths+" removable paths.");
		outstream.println("Search time: "+timer);

		timer.start();
		final AtomicInteger nextWay=new AtomicInteger(0);
		final ArrayList<CompactRemoveThread> removers=new ArrayList<CompactRemoveThread>(threads);
		for(int i=0; i<threads; i++){removers.add(new CompactRemoveThread(nextWay));}
		for(CompactRemoveThread thread : removers){thread.start();}
		long removed=0;
		failure=null;
		for(CompactRemoveThread thread : removers){
			while(thread.getState()!=Thread.State.TERMINATED){
				try{thread.join();}
				catch(InterruptedException e){Thread.currentThread().interrupt(); throw new RuntimeException(e);}
			}
			removed+=thread.removed;
			if(failure==null){failure=thread.failure;}
		}
		tables.clearOwnership();
		if(failure!=null){throw new RuntimeException("Compact wash removal failed.", failure);}
		timer.stop();
		outstream.println("Kmers removed:              \t"+removed);
		outstream.println("Removal time:               \t"+timer);
		return removed;
	}

	/** Explores all eligible tips assigned through a shared atomic index. */
	private class CompactExploreThread extends Thread {
		CompactExploreThread(ArrayList<Contig> contigs_, AtomicInteger nextTip_, boolean crossK_){
			contigs=contigs_;
			nextTip=nextTip_;
			crossK=crossK_;
		}

		@Override
		public void run(){
			try{
				for(int x=nextTip.getAndIncrement(), limit=contigs.size()*2; x<limit; x=nextTip.getAndIncrement()){
					final Contig c=contigs.get(x>>1);
					final boolean right=(x&1)==1;
					if(c.length()<kbig || (crossK && !(right ? c.rightBridgeEndpoint : c.leftBridgeEndpoint))){continue;}
					if(right){c.rightKmer(main);}
					else{c.leftKmer(main); main.rcomp();}
					scanTip();
				}
			}catch(Throwable t){failure=t;}
		}

		/** Scans a bounded unique high-depth path and tests all adjacent low-depth paths. */
		private void scanTip(){
			for(int distance=0; distance<=maxDistanceToExplore; distance++){
				if(getCount(main)<=maxCount){return;}
				final int maxPos=fillRightCounts(main, scanCounts);
				inspectLowNeighbors(main, scanCounts, candidate);
				main.rcomp();
				fillRightCounts(main, reverseCounts);
				inspectLowNeighbors(main, reverseCounts, candidate);
				main.rcomp();
				int high=0;
				for(int count : scanCounts){if(count>maxCount){high++;}}
				if(high!=1 || distance==maxDistanceToExplore){return;}
				main.addRightNumeric(maxPos);
				fillLeftCounts(main, reverseCounts);
				high=0;
				for(int count : reverseCounts){if(count>maxCount){high++;}}
				if(high!=1){return;}
			}
		}

		private void inspectLowNeighbors(final Kmer source, final int[] counts, final Kmer dest){
			for(int base=0; base<4; base++){
				final int count=counts[base];
				if(count<minCount || count>maxCount){continue;}
				dest.setFrom(source);
				dest.addRightNumeric(base);
				if(findOwner(dest)==STATUS_KEEP){continue;}
				tested++;
				if(exploreAndMark(dest, builder, leftCounts, rightCounts, minCount, maxCount,
						maxLengthToDiscard, maxDistanceToExplore, true, true, countMatrix, removeMatrix)){paths++;}
			}
		}

		private final ArrayList<Contig> contigs;
		private final AtomicInteger nextTip;
		private final boolean crossK;
		private final Kmer main=new Kmer(kbig), candidate=new Kmer(kbig);
		private final int[] scanCounts=new int[4], reverseCounts=new int[4], leftCounts=new int[4], rightCounts=new int[4];
		private final ByteBuilder builder=new ByteBuilder();
		private final long[][] countMatrix=new long[MAX_CODE+1][MAX_CODE+1];
		private final long[][] removeMatrix=new long[MAX_CODE+1][MAX_CODE+1];
		private long tested=0, paths=0;
		private Throwable failure;
	}

	/** Converts final compact-table REMOVE ownership states to deletion tombstones. */
	private class CompactRemoveThread extends Thread {
		CompactRemoveThread(AtomicInteger nextWay_){nextWay=nextWay_;}
		@Override
		public void run(){
			try{
				for(int way=nextWay.getAndIncrement(); way<tables.ways; way=nextWay.getAndIncrement()){
					removed+=tables.removeHashEntriesByOwner(way, STATUS_REMOVE);
				}
			}catch(Throwable t){failure=t;}
		}
		private final AtomicInteger nextWay;
		private long removed=0;
		private Throwable failure;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------       Dead-End Removal       ----------------*/
	/*--------------------------------------------------------------*/
	
//	private boolean valid(ByteBuilder bb, boolean doAssertion){
//		Kmer kmer=new Kmer(kbig);
//		if(bb.length()<kbig){return false;}
//		kmer.clear();
//		for(int i=0; i<bb.length; i++){
//			byte b=bb.array[i];
//			kmer.addRight(b);
//			if(kmer.len()>=kbig){
//				int count=getCount(kmer);
//				if(count<1){
//					assert(!doAssertion || false) : "count="+count+", minCount="+minCount+", maxCount="+maxCount+", kbig="+kbig+", kmer.kbig="+kmer.kbig+"\n"
//							+"kmer="+kmer.toString()+"\nbb=  "+bb.toString()+"\n";
//					kmer.clear();
//					return false;
//				}
//			}
//		}
//		kmer.clear();
//		return true;
//	}
	
	/**
	 * Explores paths bidirectionally and marks k-mers for removal based on topology.
	 * Performs complete path characterization by exploring forward from seed k-mer,
	 * then reverse complement exploration for comprehensive structure analysis.
	 *
	 * @param kmer Starting k-mer seed for exploration
	 * @param bb Sequence builder for path construction
	 * @param leftCounts Array for left-side nucleotide frequency analysis
	 * @param rightCounts Array for right-side nucleotide frequency analysis
	 * @param minCount Minimum k-mer coverage threshold
	 * @param maxCount Maximum k-mer coverage threshold
	 * @param maxLengthToDiscard Maximum path length eligible for removal
	 * @param maxDistanceToExplore Search radius limit during exploration
	 * @param prune Whether to enable pruning during exploration
	 * @param countMatrixT Matrix for tracking exploration event combinations
	 * @param removeMatrixT Matrix for tracking removal event combinations
	 * @return true if path was marked for removal, false otherwise
	 */
	public boolean exploreAndMark(Kmer kmer, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int minCount, int maxCount,
			int maxLengthToDiscard, int maxDistanceToExplore, boolean prune,
			long[][] countMatrixT, long[][] removeMatrixT){
		return exploreAndMark(kmer, bb, leftCounts, rightCounts, minCount, maxCount,
				maxLengthToDiscard, maxDistanceToExplore, prune, false, countMatrixT, removeMatrixT);
	}

	/** Compact-table variant can conservatively protect paths rejected by another seed. */
	private boolean exploreAndMark(Kmer kmer, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int minCount, int maxCount,
			int maxLengthToDiscard, int maxDistanceToExplore, boolean prune, boolean protectRejected,
			long[][] countMatrixT, long[][] removeMatrixT){
		bb.clear();
		assert(kmer.len>=kmer.kbig);
		final int initialOwner=findOwner(kmer);
		if(initialOwner>STATUS_UNEXPLORED && (!protectRejected || initialOwner==STATUS_KEEP)){return false;}
		
		assert(countWithinLimits(kmer)) : "count="+getCount(kmer)+", minCount="+minCount+", maxCount="+maxCount+"\n"+kmer.toString();
		
		bb.appendKmer(kmer);
//		assert(kmer.toString().equals(bb.toString())) : "\n"+kmer+"\n"+bb+"\n";//123
		//assert(valid(bb, true));
//		assert(tables.getCount(kmer)==1) : tables.getCount(kmer);//count>0 && count<=maxCount
		final int rightCode=explore(kmer, bb, leftCounts, rightCounts, minCount, maxCount, maxDistanceToExplore);
		//assert(valid(bb, true));
		
		bb.reverseComplementInPlace();
		//assert(valid(bb, true));
		kmer=tables.rightmostKmer(bb, kmer);
		assert(getCount(kmer)>0) : "count="+getCount(kmer)+", minCount="+minCount+", maxCount="+maxCount+", rightCode="+rightCode+"\n"+kmer.toString();//123
//		assert(tables.getCount(kmer)==1) : tables.getCount(kmer);//count>0 && count<=maxCount
		final int leftCode=explore(kmer, bb, leftCounts, rightCounts, minCount, maxCount, maxDistanceToExplore);
		//assert(valid(bb, true));
		
		kmer=tables.rightmostKmer(bb, kmer);//123
		assert(getCount(kmer)>0) : "count="+getCount(kmer)+", minCount="+minCount+", maxCount="+maxCount+", leftCode="+leftCode+", rightCode="+rightCode+"\n"+kmer.toString();//123

		final int min=Tools.min(rightCode, leftCode);
		final int max=Tools.max(rightCode, leftCode);
		
		countMatrixT[min][max]++;
		
		if(rightCode==TOO_LONG || rightCode==TOO_DEEP || rightCode==LOOP || rightCode==F_BRANCH){
			claim(bb, protectRejected ? STATUS_KEEP : STATUS_EXPLORED, false, kmer);
			return false;
		}
		
		if(leftCode==TOO_LONG || leftCode==TOO_DEEP || leftCode==LOOP || leftCode==F_BRANCH){
			claim(bb, protectRejected ? STATUS_KEEP : STATUS_EXPLORED, false, kmer);
			return false;
		}
		
		if(bb.length()-kbig>maxLengthToDiscard){
			claim(bb, protectRejected ? STATUS_KEEP : STATUS_EXPLORED, false, kmer);
			return false;
		}
		
		if(removeHair && min==DEAD_END){
			if(max==DEAD_END || max==B_BRANCH){
				removeMatrixT[min][max]++;
				claim(bb, STATUS_REMOVE, false, kmer);
				final boolean success=!protectRejected || findOwner(bb, STATUS_REMOVE, kmer)<=STATUS_REMOVE;
				if(verbose || verbose2){System.err.println("Claiming ("+rightCode+","+leftCode+") length "+bb.length()+": "+bb);}
				if(!success && protectRejected){claim(bb, STATUS_KEEP, false, kmer);}
				else{assert(success);}
				return success;
			}
		}
		
		if(removeBubbles){
			if(rightCode==B_BRANCH && leftCode==B_BRANCH){
				removeMatrixT[min][max]++;
				claim(bb, STATUS_REMOVE, false, kmer);
				final boolean success=!protectRejected || findOwner(bb, STATUS_REMOVE, kmer)<=STATUS_REMOVE;
				if(verbose || verbose2){System.err.println("Claiming ("+rightCode+","+leftCode+") length "+bb.length()+": "+bb);}
				if(!success && protectRejected){claim(bb, STATUS_KEEP, false, kmer);}
				else{assert(success);}
				return success;
			}
		}
		
		claim(bb, protectRejected ? STATUS_KEEP : STATUS_EXPLORED, false, kmer);
		return false;
	}

	/** Classifies one low-depth path without changing shared counts or ownership. */
	public int classifyPath(Kmer kmer, ByteBuilder bb, int[] leftCounts, int[] rightCounts){
		bb.clear();
		if(!countWithinLimits(kmer)){return 0;}
		bb.appendKmer(kmer);
		final int rightCode=explore(kmer, bb, leftCounts, rightCounts, minCount, maxCount, maxDistanceToExplore);
		bb.reverseComplementInPlace();
		kmer=tables.rightmostKmer(bb, kmer);
		final int leftCode=explore(kmer, bb, leftCounts, rightCounts, minCount, maxCount, maxDistanceToExplore);
		if(rightCode==TOO_LONG || rightCode==TOO_DEEP || rightCode==LOOP || rightCode==F_BRANCH
				|| leftCode==TOO_LONG || leftCode==TOO_DEEP || leftCode==LOOP || leftCode==F_BRANCH
				|| bb.length()-kbig>maxLengthToDiscard){return 0;}
		final int min=Tools.min(rightCode, leftCode), max=Tools.max(rightCode, leftCode);
		if(removeHair && min==DEAD_END && (max==DEAD_END || max==B_BRANCH)){return 1;}
		if(removeBubbles && rightCode==B_BRANCH && leftCode==B_BRANCH){return 2;}
		return 0;
	}
	
	/** Explores a single unbranching path in the forward (right) direction.
	 * @param kmer
	 * @param bb
	 * @param leftCounts
	 * @param rightCounts
	 * @param minCount
	 * @param maxCount
	 * @param maxLength0
	 * @return A termination code such as DEAD_END
	 */
	public int explore(Kmer kmer, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int minCount, int maxCount, int maxLength0){
		if(verbose){outstream.println("Entering explore with bb.length()="+bb.length());}
		if(bb.length()==0){bb.appendKmer(kmer);}
		assert(bb.length()>=kmer.kbig && kmer.len>=kmer.kbig) : bb.length()+", "+kmer.len+", "+kmer.kbig;
		//assert(valid(bb, true));
		
		final int initialLength=bb.length();
		final int maxLength=maxLength0+kbig;
		
		final long firstKey=kmer.xor(), firstKey2=kmer.xor2();
		AbstractKmerTableU table=tables.getTableAbstract(kmer);
		int count=table.getValue(kmer);
//		kmer.verify(false);
//		kmer.verify(true);
		assert(count>=minCount && count<=maxCount) : count+", "+Kmer.MASK_CORE+", "+kmer.verify(false)+", "+kmer.verify(true);
		
		int nextRightMaxPos=fillRightCounts(kmer, rightCounts);
		int nextRightMax=rightCounts[nextRightMaxPos];
		if(nextRightMax<minCount){
			if(verbose){outstream.println("Returning DEAD_END: rightMax="+nextRightMax);}
			return DEAD_END;
		}
		
		while(bb.length()<=maxLength){
			
			final int rightMaxPos=nextRightMaxPos;
			final int rightMax=rightCounts[rightMaxPos];
			final int rightSecondPos=Tools.secondHighestPosition(rightCounts);
			final int rightSecond=rightCounts[rightSecondPos];
			
			if(verbose){
				outstream.println("kmer: "+toText(kmer));
				outstream.println("Right counts: "+count+", "+Arrays.toString(rightCounts));
				outstream.println("rightMaxPos="+rightMaxPos);
				outstream.println("rightMax="+rightMax);
				outstream.println("rightSecondPos="+rightSecondPos);
				outstream.println("rightSecond="+rightSecond);
			}
			assert(count>0);
//			assert(getCount(kmer)==count); //123
			
			final int prevCount=count;
			
			//Generate the new base
			final byte b=AminoAcid.numberToBase[rightMaxPos];
			final long x=rightMaxPos;
			long evicted=kmer.addRightNumeric(x);
			
//			assert(getCount(kmer)==rightMax); //123
			
			//Now consider the next kmer
			if(kmer.xor()==firstKey && kmer.xor2()==firstKey2){
				if(verbose){outstream.println("Returning LOOP");}
				return LOOP;
			}
			table=tables.getTableAbstract(kmer);
			
			assert(table.getValue(kmer)==rightMax || rightMax==0);
			count=rightMax;
			assert(count>0);
			
			
			{//Fill right and look for dead end
				nextRightMaxPos=fillRightCounts(kmer, rightCounts);
				nextRightMax=rightCounts[nextRightMaxPos];
				if(nextRightMax<minCount){
					if(verbose){outstream.println("Returning DEAD_END: rightMax="+rightMax);}
					return DEAD_END;
				}
			}
			
			
			{//Look left
				final int leftMaxPos=fillLeftCounts(kmer, leftCounts);
				final int leftMax=leftCounts[leftMaxPos];
				final int leftSecondPos=Tools.secondHighestPosition(leftCounts);
				final int leftSecond=leftCounts[leftSecondPos];
				
//				assert(leftMax==1 || leftMax==0) : prevCount+" -> "+Arrays.toString(leftCounts)+", "+count+", "+Arrays.toString(rightCounts);
				
				if(verbose){
					outstream.println("Left counts: "+count+", "+Arrays.toString(leftCounts));
					outstream.println("leftMaxPos="+leftMaxPos);
					outstream.println("leftMax="+leftMax);
					outstream.println("leftSecondPos="+leftSecondPos);
					outstream.println("leftSecond="+leftSecond);
				}
				
				if(leftSecond>=minCount || leftMax>prevCount){//Backward branch
//					assert(leftSecond==1) : prevCount+" -> "+Arrays.toString(leftCounts)+", "+count+", "+Arrays.toString(rightCounts);
					if(leftMax>prevCount){
						if(verbose){outstream.println("Returning B_BRANCH_LOWER: " +
								"count="+count+", prevCount="+prevCount+", leftMax="+leftMax+", leftSecond="+leftSecond);}
						return B_BRANCH;
					}else{
						assert(leftMax==prevCount);
						if(leftMax>=2*leftSecond){//This constant is adjustable
//							assert(false) : prevCount+" -> "+Arrays.toString(leftCounts)+", "+count+", "+Arrays.toString(rightCounts);
							//keep going
						}else{
							if(verbose){outstream.println("Returning B_BRANCH_SIMILAR: " +
								"count="+count+", prevCount="+prevCount+", leftMax="+leftMax+", leftSecond="+leftSecond);}
							return B_BRANCH;
						}
					}
				}
				
			}
			
			//Look right
			if(rightSecond>=minCount){
				if(verbose){outstream.println("Returning F_BRANCH: rightSecond="+rightSecond);}
				return F_BRANCH;
			}
			
			if(count>maxCount){
				if(verbose){outstream.println("Returning TOO_DEEP: rightMax="+rightMax);}
				return TOO_DEEP;
			}
			
			bb.append(b);
//			assert(valid(bb, false)); //123
//			if(!valid(bb, false)){
//				System.err.println("kbig="+kbig+", "+kmer.kbig+"\nb="+
//						(char)b+"\n"+Arrays.toString(rightCounts)+"\ncount="+count+", "+getCount(kmer)+"\nrightMax="+rightMax+"\nrightMax="+rightMax);
//				System.err.println("\nkmer="+kmer+"\nbb=  "+bb);
//				kmer.addLeftNumeric(evicted);
//				System.err.println("prev="+kmer+"\nprevCount="+prevCount+", "+getCount(kmer));
//				valid(bb, true);
//			}
			if(verbose){outstream.println("Added base "+(char)b);}
		}
		
		assert(bb.length()>maxLength);
		if(verbose){outstream.println("Returning TOO_LONG: length="+bb.length());}
		return TOO_LONG;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         ExploreThread        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Searches for dead ends.
	 */
	class ExploreThread extends AbstractExploreThread{
		
		/**
		 * Constructor
		 */
		public ExploreThread(int id_){
			super(id_, kbig);
		}
		
		@Override
		boolean processNextTable(final Kmer kmer, Kmer temp){
			final int tnum=nextTable.getAndAdd(1);
			if(tnum>=tables.ways){return false;}
			final HashArrayU1D table=tables.getTable(tnum);
			final int[] counts=table.values();
			final int max=counts.length;
//			final int max=table.arrayLength();
			if(startFromHighCounts){
//				for(int cell=0; cell<max; cell++){
//					int x=processCell_high(table, cell, kmer, temp);
//					deadEndsFoundT+=x;
//				}
				for(int cell=0; cell<max; cell++){//Not noticeably faster
					final int count=counts[cell];
					if(count>maxCount){
						int x=processCell_high(table, cell, kmer, temp, count);
						deadEndsFoundT+=x;
					}
				}
			}else{
				for(int cell=0; cell<max; cell++){
					int x=processCell_low(table, cell, kmer, temp);
					deadEndsFoundT+=x;
				}
			}
			return true;
		}
		
		@Override
		boolean processNextVictims(final Kmer kmer, Kmer temp){
			final int tnum=nextVictims.getAndAdd(1);
			if(tnum>=tables.ways){return false;}
			final HashArrayU1D table=tables.getTable(tnum);
			final HashForestU forest=table.victims();
			final int max=forest.arrayLength();
			if(startFromHighCounts){
				for(int cell=0; cell<max; cell++){
					KmerNodeU kn=forest.getNode(cell);
					int x=traverseKmerNodeU_high(kn, kmer, temp);
					deadEndsFoundT+=x;
				}
			}else{
				for(int cell=0; cell<max; cell++){
					KmerNodeU kn=forest.getNode(cell);
					int x=traverseKmerNodeU_low(kn, kmer, temp);
					deadEndsFoundT+=x;
				}
			}
			return true;
		}
		
		/**
		 * Recursively traverses k-mer node tree for high-coverage processing.
		 * Processes current node then recursively handles left and right subtrees.
		 *
		 * @param kn K-mer node to process (may be null)
		 * @param kmer Primary k-mer object for operations
		 * @param temp Temporary k-mer object for calculations
		 * @return Number of dead ends found in this subtree
		 */
		private int traverseKmerNodeU_high(KmerNodeU kn, Kmer kmer, Kmer temp){
			int sum=0;
			if(kn!=null){
				sum+=processKmerNodeU_high(kn, kmer, temp);
				if(kn.left()!=null){
					sum+=traverseKmerNodeU_high(kn.left(), kmer, temp);
				}
				if(kn.right()!=null){
					sum+=traverseKmerNodeU_high(kn.right(), kmer, temp);
				}
			}
			return sum;
		}
		
		/**
		 * Recursively traverses k-mer node tree for low-coverage processing.
		 * Processes current node then recursively handles left and right subtrees.
		 *
		 * @param kn K-mer node to process (may be null)
		 * @param kmer Primary k-mer object for operations
		 * @param temp Temporary k-mer object for calculations
		 * @return Number of dead ends found in this subtree
		 */
		private int traverseKmerNodeU_low(KmerNodeU kn, Kmer kmer, Kmer temp){
			int sum=0;
			if(kn!=null){
				sum+=processKmerNodeU_low(kn, kmer, temp);
				if(kn.left()!=null){
					sum+=traverseKmerNodeU_low(kn.left(), kmer, temp);
				}
				if(kn.right()!=null){
					sum+=traverseKmerNodeU_low(kn.right(), kmer, temp);
				}
			}
			return sum;
		}
		
		/*--------------------------------------------------------------*/
		
		//old
		/**
		 * Processes individual hash table cell with low-coverage algorithm.
		 * Validates coverage thresholds and ownership before exploration.
		 *
		 * @param table Hash table containing the cell
		 * @param cell Cell index to process
		 * @param kmer0 K-mer object to populate with cell data
		 * @param temp Temporary k-mer for calculations
		 * @return 1 if dead end found and processed, 0 otherwise
		 */
		private int processCell_low(HashArrayU1D table, int cell, Kmer kmer0, Kmer temp){
			int count=table.readCellValue(cell);
			if(count<minSeed || count>maxCount){return 0;}
			int owner=table.getCellOwner(cell);
			if(owner>STATUS_UNEXPLORED){return 0;}
			Kmer kmer=table.fillKmer(cell, kmer0);
			if(kmer==null){return 0;}
			if(verbose){outstream.println("id="+id+" processing cell "+cell+"; \tkmer="+kmer);}
			
			return processKmer_low(kmer, temp);
		}
		
		//old
		/**
		 * Processes k-mer node with low-coverage algorithm.
		 * Validates coverage and ownership before delegating to k-mer processing.
		 *
		 * @param kn K-mer node containing data to process
		 * @param kmer K-mer object to populate with node data
		 * @param temp Temporary k-mer for calculations
		 * @return 1 if dead end found and processed, 0 otherwise
		 */
		private int processKmerNodeU_low(KmerNodeU kn, Kmer kmer, Kmer temp){
			kmer.setFrom(kn.pivot());
			final int count=kn.getValue(kmer);
			if(count<minSeed || count>maxCount){return 0;}
			int owner=kn.getOwner(kmer);
			if(owner>STATUS_UNEXPLORED){return 0;}
			
			return processKmer_low(kmer, temp);
		}
		
		//old
		/**
		 * Core k-mer processing with exploration and marking logic.
		 * Increments test counter and delegates to exploreAndMark for analysis.
		 *
		 * @param original K-mer to analyze for dead-end characteristics
		 * @param temp Temporary k-mer for calculations
		 * @return 1 if k-mer was marked for removal, 0 otherwise
		 */
		private int processKmer_low(Kmer original, Kmer temp){
			kmersTestedT++;
			boolean b=exploreAndMark(original, builderT, leftCounts, rightCounts, minCount, maxCount, maxLengthToDiscard, maxDistanceToExplore, true
					, countMatrixT, removeMatrixT
				);
			return b ? 1 : 0;
		}
		

		
		/*--------------------------------------------------------------*/
		
		//new
		/**
		 * Processes hash table cell using high-coverage optimization strategy.
		 * Handles k-mers exceeding maxCount threshold through neighbor exploration.
		 *
		 * @param table Hash table containing the cell
		 * @param cell Cell index to process
		 * @param kmer0 K-mer object to populate with cell data
		 * @param temp Temporary k-mer for calculations
		 * @param count Pre-calculated k-mer count for efficiency
		 * @return Number of dead ends found through neighbor processing
		 */
		private int processCell_high(HashArrayU1D table, int cell, Kmer kmer0, Kmer temp, int count){
			
//		private int processCell_high(HashArrayU1D table, int cell, Kmer kmer0, Kmer temp){
//			int count=table.readCellValue(cell);
			
			if(count<=maxCount){return 0;}
//			if(shaveFast && maxCount==1 && count>=6){return 0;}
//			int owner=table.getCellOwner(cell);
//			if(owner>STATUS_UNEXPLORED){return 0;}
			Kmer kmer=table.fillKmer(cell, kmer0);
			if(kmer==null){return 0;}
			if(verbose){outstream.println("id="+id+" processing cell "+cell+"; \tkmer="+kmer);}
//			assert(false) : shaveFast;
			
			return processKmer_high(kmer, temp, count);
		}
		
		//new
		/**
		 * Processes k-mer node using high-coverage optimization strategy.
		 * Handles nodes with counts exceeding maxCount threshold.
		 *
		 * @param kn K-mer node containing data to process
		 * @param kmer K-mer object to populate with node data
		 * @param temp Temporary k-mer for calculations
		 * @return Number of dead ends found through neighbor processing
		 */
		private int processKmerNodeU_high(KmerNodeU kn, Kmer kmer, Kmer temp){
			kmer.setFrom(kn.pivot());
			final int count=kn.getValue(kmer);
			if(count<=maxCount){return 0;}
//			if(shaveFast && maxCount==1 && count>=6){return 0;}
//			int owner=kn.getOwner(kmer);
//			if(owner>STATUS_UNEXPLORED){return 0;}
			
			return processKmer_high(kmer, temp, count);
		}
		
		//new
		/**
		 * Processes high-coverage k-mers by exploring their neighbors.
		 * Uses left-side extension in both orientations to find removal candidates.
		 *
		 * @param original High-coverage k-mer to analyze
		 * @param kmer Working k-mer for neighbor exploration
		 * @param count0 Original k-mer count for reference
		 * @return Total number of dead ends found through neighbor analysis
		 */
		private int processKmer_high(final Kmer original, Kmer kmer, final int count0){
			int sum=0;
			
//			assert(false) : "This loop is broken and removes kmers with counts outside the range, including non-existing kmers.";
//			for(long i=0; i<4; i++){
//				kmer.setFrom(original);
//				long old=kmer.addRightNumeric(i);
//				HashArrayU1D table=(HashArrayU1D) tables.getTable(kmer);
//				int count=table.getCount(kmer);
//				if(count>0 && count<=maxCount && tables.getOwner(kmer)<=STATUS_UNEXPLORED){
//					kmersTestedT++;
//					boolean b=exploreAndMark(kmer, builderT, leftCounts, rightCounts,
//							minCount, maxCount, maxLengthToDiscard, maxDistanceToExplore, true,
//							countMatrixT, removeMatrixT);
//					if(b){sum++;}
//				}
//			}
			
			assert(original.len()==kbig);
			
			//This loop appears to work correctly.
			sum+=processKmer_high_leftLoop(original, kmer, count0);
			
			original.rcomp();
			assert(original.len()==kbig);
			sum+=processKmer_high_leftLoop(original, kmer, count0);
			return sum;
		}
		
		/**
		 * Explores left-side neighbors of high-coverage k-mers for removal candidates.
		 * Implements optional shaveVFast acceleration through branch pre-filtering.
		 *
		 * @param original High-coverage k-mer to analyze neighbors of
		 * @param kmer Working k-mer for neighbor exploration
		 * @param count0 Original k-mer count for comparison
		 * @return Number of removal candidates found among neighbors
		 */
		private int processKmer_high_leftLoop(final Kmer original, Kmer kmer, final int count0){
			if(shaveVFast){//Optional accelerator; only examines kmers that actually branch
				int inf=Integer.MAX_VALUE;
				int maxHighBranch=-1, minHighBranch=inf, highBranches=0;
				int maxLowBranch=-1, minLowBranch=inf, lowBranches=0;
				for(long i=0; i<4; i++){
					kmer.setFrom(original);
					long old=kmer.addLeftNumeric(i);
					assert(kmer.len()>=kbig) : kmer.len()+", "+kbig;
					HashArrayU1D table=(HashArrayU1D) tables.getTable(kmer);
					int count=table.getValue(kmer);
					if(count>0){
						if(count<=maxCount){
							minLowBranch=Tools.min(count, minLowBranch);
							maxLowBranch=Tools.max(count, maxLowBranch);
							lowBranches++;
						}else{
							minHighBranch=Tools.min(count, minHighBranch);
							maxHighBranch=Tools.max(count, maxHighBranch);
							highBranches++;
						}
					}
				}

				if(maxLowBranch<0){return 0;} //This is fine
//				if(maxHighBranch<0){return 0;} //Speed increase but quality decrease

				if(highBranches+lowBranches<2){return 0;} //Speed increase (25% for shave) but contiguity decrease (~1%).
//				if(maxLowBranch>0 && maxLowBranch*16<minHighBranch){return 0;} //Speed increase but quality decrease
			}
			
			int sum=0;
			for(long i=0; i<4; i++){
				kmer.setFrom(original);
				long old=kmer.addLeftNumeric(i);
				assert(kmer.len()>=kbig) : kmer.len()+", "+kbig;
				HashArrayU1D table=(HashArrayU1D) tables.getTable(kmer);
				int count=table.getValue(kmer);
				if(count>0 && count<=maxCount && table.getOwner(kmer)<=STATUS_UNEXPLORED){
					kmersTestedT++;
					boolean b=exploreAndMark(kmer, builderT, leftCounts, rightCounts,
							minCount, maxCount, maxLengthToDiscard, maxDistanceToExplore, true,
							countMatrixT, removeMatrixT);
					if(b){sum++;}
				}
			}
			return sum;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          ShaveThread         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Removes dead-end kmers.
	 */
	private class ShaveThread extends AbstractShaveThread{

		/**
		 * Constructor
		 */
		public ShaveThread(int id_){
			super(id_);
		}
		
		@Override
		boolean processNextTable(){
			final int tnum=nextTable.getAndAdd(1);
			if(tnum>=tables.ways){return false;}
//			long x=0;
			final HashArrayU1D table=tables.getTable(tnum);
			final AtomicIntegerArray owners=table.owners();
			final int[] values=table.values();
			final int max=table.arrayLength();
			for(int cell=0; cell<max; cell++){
				if(owners.get(cell)==STATUS_REMOVE){
//					x++;
					values[cell]=0;
				}
			}
			for(KmerNodeU kn : table.victims().array()){
				if(kn!=null){traverseKmerNodeU(kn);}
			}
			
			table.clearOwnership();
			kmersRemovedT+=table.regenerate(0);
//			outstream.println(x);
			return true;
		}
		
		/**
		 * Recursively removes k-mer nodes marked for deletion.
		 * Zeros counts for nodes with STATUS_REMOVE ownership.
		 * @param kn K-mer node to process (null-safe)
		 */
		private void traverseKmerNodeU(KmerNodeU kn){
			if(kn==null){return;}
			if(kn.owner()==STATUS_REMOVE){kn.set(0);}
			traverseKmerNodeU(kn.left());
			traverseKmerNodeU(kn.right());
		}
		
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Recall Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Claims ownership of all k-mers in a sequence path.
	 * Delegates to array-based claiming with sequence length.
	 *
	 * @param bb Sequence builder containing path to claim
	 * @param id Thread identifier for ownership assignment
	 * @param exitEarly Whether to stop on first ownership conflict
	 * @param kmer Working k-mer object for operations
	 * @return true if all k-mers successfully claimed, false otherwise
	 */
	public boolean claim(final ByteBuilder bb, final int id, final boolean exitEarly, Kmer kmer){
		return claim(bb.array, bb.length(), id, exitEarly, kmer);
	}
	
	/**
	 * Claims ownership of all k-mers in a sequence array.
	 * Iterates through sequence maintaining rolling k-mer state and
	 * attempts ownership claiming for each valid k-mer position.
	 *
	 * @param bases Sequence array containing k-mers to claim
	 * @param blen Length of valid sequence data in array
	 * @param id Thread identifier for ownership assignment
	 * @param exitEarly Whether to stop on first ownership conflict
	 * @param kmer Working k-mer object for rolling state maintenance
	 * @return true if all k-mers successfully claimed, false on any conflict
	 */
	public boolean claim(final byte[] bases, final int blen, final int id, boolean exitEarly, final Kmer kmer){
		if(blen<kbig){return false;}
		if(verbose){outstream.println("Thread "+id+" claim start.");}
		int len=0;
		kmer.clear();
		boolean success=true;
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=0; i<blen && success; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			kmer.addRight(b);

			if(x<0){len=0;}
			else{len++;}
			assert(len==kmer.len);

			if(len>=kbig){
				assert(countWithinLimits(kmer)) : "count="+getCount(kmer)+", minCount="+minCount+", maxCount="+maxCount+"\n"
						+ "len="+len+", i="+i+", blen="+blen+"\n"
								+ ""+kmer.toString();
				success=claim(kmer, id/*, rid, i*/);
				success=(success || !exitEarly);
			}
		}
		return success;
	}
	
	/**
	 * Validates k-mer count falls within acceptable range.
	 * @param kmer K-mer to validate count for
	 * @return true if count between minCount and maxCount inclusive
	 */
	final boolean countWithinLimits(Kmer kmer){
		int count=getCount(kmer);
		return count>=minCount && count<=maxCount;
	}
	
	/** Retrieves k-mer frequency count from hash tables */
	int getCount(Kmer kmer){return tables.getCount(kmer);}
	/** Claims ownership of single k-mer for thread coordination */
	boolean claim(Kmer kmer, int id){return tables.claim(kmer, id);}
	/** Claims ownership of sequence in both orientations */
	boolean doubleClaim(ByteBuilder bb, int id/*, long rid*/, Kmer kmer){return tables.doubleClaim(bb, id, kmer/*, rid*/);}
//	boolean claim(ByteBuilder bb, int id, /*long rid, */boolean earlyExit, Kmer kmer){return tables.claim(bb, id/*, rid*/, earlyExit, kmer);}
//	boolean claim(byte[] array, int len, int id, /*long rid, */boolean earlyExit, Kmer kmer){return tables.claim(array, len, id/*, rid*/, earlyExit, kmer);}
	/** Retrieves current ownership status of k-mer */
	int findOwner(Kmer kmer){return tables.findOwner(kmer);}
	/** Finds ownership status along sequence path */
	int findOwner(ByteBuilder bb, int id, Kmer kmer){return tables.findOwner(bb, id, kmer);}
	/** Finds ownership status along sequence array */
	int findOwner(byte[] array, int len, int id, Kmer kmer){return tables.findOwner(array, len, id, kmer);}
	/** Releases ownership of k-mers in sequence */
	void release(ByteBuilder bb, int id, Kmer kmer){tables.release(bb, id, kmer);}
	/** Releases ownership of k-mers in sequence array */
	void release(byte[] array, int len, int id, Kmer kmer){tables.release(array, len, id, kmer);}
	/** Populates array with counts for right-side k-mer extensions */
	int fillRightCounts(Kmer kmer, int[] counts){return tables.fillRightCounts(kmer, counts);}
	/** Populates array with counts for left-side k-mer extensions */
	int fillLeftCounts(Kmer kmer, int[] counts){return tables.fillLeftCounts(kmer, counts);}
	/** Converts k-mer to human-readable text representation */
	static StringBuilder toText(Kmer kmer){return AbstractKmerTableU.toText(kmer);}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	AbstractKmerTableSet tables(){return tables;}
	
	/** K-mer hash table set supporting unlimited k-mer lengths */
	final KmerTableSetU tables;
	
}
