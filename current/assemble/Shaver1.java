package assemble;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicIntegerArray;

import dna.AminoAcid;
import kmer.AbstractKmerTable;
import kmer.AbstractKmerTableSet;
import kmer.HashArray1D;
import kmer.HashForest;
import kmer.KmerNode;
import kmer.KmerTableSet;
import shared.Tools;
import structures.ByteBuilder;
import ukmer.Kmer;

/**
 * Designed for removal of dead ends (aka hairs).
 * @author Brian Bushnell
 * @date Jun 26, 2015
 *
 */
public class Shaver1 extends Shaver {
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/

	
	/**
	 * Constructor with default parameters for dead-end removal.
	 * Uses default values for count thresholds, branch multiplier, and length limits.
	 * @param tables_ The k-mer table set to operate on
	 * @param threads_ Number of threads for parallel processing
	 */
	public Shaver1(KmerTableSet tables_, int threads_){
		this(tables_, threads_, 1, 1, 1, 1, 3, 100, 100, true, true);
	}
	
	/**
	 * Full constructor with all parameters for customized dead-end removal.
	 *
	 * @param tables_ The k-mer table set to operate on
	 * @param threads_ Number of threads for parallel processing
	 * @param minCount_ Minimum k-mer count to consider for extension
	 * @param maxCount_ Maximum k-mer count to consider (filters high-coverage repeats)
	 * @param minSeed_ Minimum count for seed k-mers to start exploration
	 * @param minCountExtend_ Minimum count required to extend a path
	 * @param branchMult2_ Branch multiplier threshold for branch detection
	 * @param maxLengthToDiscard_ Maximum length of paths that can be removed
	 * @param maxDistanceToExplore_ Maximum distance to explore in each direction
	 * @param removeHair_ Whether to remove hair/dead-end structures
	 * @param removeBubbles_ Whether to remove bubble structures between branches
	 */
	public Shaver1(KmerTableSet tables_, int threads_,
			int minCount_, int maxCount_, int minSeed_, int minCountExtend_, float branchMult2_, int maxLengthToDiscard_, int maxDistanceToExplore_,
			boolean removeHair_, boolean removeBubbles_){
		super(tables_, threads_, minCount_, maxCount_, minSeed_, minCountExtend_, branchMult2_, maxLengthToDiscard_, maxDistanceToExplore_, removeHair_, removeBubbles_);
		tables=tables_;
		k=tables.k;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	final AbstractExploreThread makeExploreThread(int id_){return new ExploreThread(id_);}
	@Override
	final AbstractShaveThread makeShaveThread(int id_){return new ShaveThread(id_);}
	
	
	/*--------------------------------------------------------------*/
	/*----------------       Dead-End Removal       ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Explores bidirectional paths from a k-mer and determines if they should be removed.
	 * Explores in both directions from the starting k-mer to determine the full extent
	 * of the unbranching path. Classifies termination conditions (dead end, branch, etc.)
	 * and decides whether to mark the path for removal based on the criteria.
	 *
	 * @param kmer Starting k-mer for exploration
	 * @param bb ByteBuilder for sequence construction
	 * @param leftCounts Array to store left extension counts
	 * @param rightCounts Array to store right extension counts
	 * @param minCount Minimum count threshold for valid extensions
	 * @param maxCount Maximum count threshold to avoid repeats
	 * @param maxLengthToDiscard Maximum length of paths that can be discarded
	 * @param maxDistanceToExplore Maximum distance to explore in each direction
	 * @param prune Whether to mark paths for removal
	 * @param countMatrixT Matrix to count termination condition combinations
	 * @param removeMatrixT Matrix to count actual removals by termination type
	 * @return true if the path was marked for removal, false otherwise
	 */
	public boolean exploreAndMark(long kmer, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int minCount, int maxCount,
			int maxLengthToDiscard, int maxDistanceToExplore, boolean prune,
			long[][] countMatrixT, long[][] removeMatrixT){
		bb.clear();
		if(findOwner(kmer)>STATUS_UNEXPLORED){return false;}
		
		bb.appendKmer(kmer, k);
		final int rightCode=explore(kmer, bb, leftCounts, rightCounts, minCount, maxCount, maxDistanceToExplore);
		
		bb.reverseComplementInPlace();
		kmer=tables.rightmostKmer(bb);
		final int leftCode=explore(kmer, bb, leftCounts, rightCounts, minCount, maxCount, maxDistanceToExplore);

		final int min=Tools.min(rightCode, leftCode);
		final int max=Tools.max(rightCode, leftCode);
		
		countMatrixT[min][max]++;
		
		if(rightCode==TOO_LONG || rightCode==TOO_DEEP || rightCode==LOOP || rightCode==F_BRANCH){
			claim(bb, STATUS_EXPLORED, false);
			return false;
		}
		
		if(leftCode==TOO_LONG || leftCode==TOO_DEEP || leftCode==LOOP || leftCode==F_BRANCH){
			claim(bb, STATUS_EXPLORED, false);
			return false;
		}
		
		if(bb.length()-k>maxLengthToDiscard){
			claim(bb, STATUS_EXPLORED, false);
			return false;
		}
		
		if(removeHair && min==DEAD_END){
			if(max==DEAD_END || max==B_BRANCH){
				removeMatrixT[min][max]++;
				boolean success=claim(bb, STATUS_REMOVE, false);
				if(verbose || verbose2){System.err.println("Claiming ("+rightCode+","+leftCode+") length "+bb.length()+": "+bb);}
				assert(success);
				return true;
			}
		}
		
		if(removeBubbles){
			if(rightCode==B_BRANCH && leftCode==B_BRANCH){
				removeMatrixT[min][max]++;
				boolean success=claim(bb, STATUS_REMOVE, false);
				if(verbose || verbose2){System.err.println("Claiming ("+rightCode+","+leftCode+") length "+bb.length()+": "+bb);}
				assert(success);
				return true;
			}
		}
		
		claim(bb, STATUS_EXPLORED, false);
		return false;
	}
	
	/** Explores a single unbranching path in the forward direction.
	 * @param kmer
	 * @param bb
	 * @param leftCounts
	 * @param rightCounts
	 * @param minCount
	 * @param maxCount
	 * @param maxLength0
	 * @return A termination code such as DEAD_END
	 */
	public int explore(long kmer, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int minCount, int maxCount, int maxLength0){
		if(verbose){outstream.println("Entering explore with bb.length()="+bb.length());}
		assert(bb.length()==0 || tables.rightmostKmer(bb)==kmer);
		if(bb.length()==0){bb.appendKmer(kmer, k);}
		
		final int initialLength=bb.length();
		final int maxLength=maxLength0+k;
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		
		long rkmer=AminoAcid.reverseComplementBinaryFast(kmer, k);
		
		long key=toValue(kmer, rkmer);
		final long firstKey=key;
		HashArray1D table=tables.getTableForKey(key);
		int count=table.getValue(key);
		assert(count>=minCount && count<=maxCount);
		
		int nextRightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
		int nextRightMax=rightCounts[nextRightMaxPos];
		if(nextRightMax<minCount){
			if(verbose){outstream.println("Returning DEAD_END: rightMax="+nextRightMax);}
			return DEAD_END;
		}
		
		while(bb.length()<=maxLength){
			
			final int rightMaxPos=nextRightMaxPos;
			final int rightMax=nextRightMax;
			final int rightSecondPos=Tools.secondHighestPosition(rightCounts);
			final int rightSecond=rightCounts[rightSecondPos];
			final int prevCount=count;
			count=rightMax;
			
			if(verbose){
				outstream.println("kmer: "+toText(kmer, k)+", "+toText(rkmer, k));
				outstream.println("Right counts: "+prevCount+", "+Arrays.toString(rightCounts));
				outstream.println("rightMaxPos="+rightMaxPos);
				outstream.println("rightMax="+rightMax);
				outstream.println("rightSecondPos="+rightSecondPos);
				outstream.println("rightSecond="+rightSecond);
			}
			
			//Generate the new base
			final byte b=AminoAcid.numberToBase[rightMaxPos];
			final long x=rightMaxPos;
			final long x2=AminoAcid.numberToComplement[rightMaxPos];
			kmer=((kmer<<2)|(long)x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			
			//Now consider the next kmer
			key=toValue(kmer, rkmer);
			if(key==firstKey){
				if(verbose){outstream.println("Returning LOOP");}
				return LOOP;
			}
			table=tables.getTableForKey(key);
			
			assert(table.getValue(key)==rightMax || rightMax==0);
			
			{//Fill right and look for dead end
				nextRightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
				nextRightMax=rightCounts[nextRightMaxPos];
				if(nextRightMax<minCount){
					if(verbose){outstream.println("Returning DEAD_END: rightMax="+rightMax);}
					return DEAD_END;
				}
			}
			
			
			{//Look left
				final int leftMaxPos=fillLeftCounts(kmer, rkmer, leftCounts, mask, shift2);
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
			if(verbose){outstream.println("Added base "+(char)b);}
		}
		
		assert(bb.length()>maxLength);
		if(verbose){outstream.println("Returning TOO_LONG: length="+bb.length());}
		return TOO_LONG;
	}
	
	/** Explores a single unbranching path in the forward direction.
	 * Returns reason for ending in this direction:
	 *  DEAD_END, TOO_LONG, TOO_DEEP, F_BRANCH, B_BRANCH */
	public int explore2(long kmer, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int minCount, int maxCount, int maxLength0){
		if(verbose){outstream.println("Entering explore with bb.length()="+bb.length());}
		assert(bb.length()==0 || tables.rightmostKmer(bb)==kmer);
		if(bb.length()==0){bb.appendKmer(kmer, k);}
		
		final int maxLength=maxLength0+k;
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		
		long rkmer=AminoAcid.reverseComplementBinaryFast(kmer, k);
		
		long key=toValue(kmer, rkmer);
		final long firstKey=key;
		HashArray1D table=tables.getTableForKey(key);
		int count=table.getValue(key);
		assert(count>=minCount && count<=maxCount);

		int rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
		int rightMax=rightCounts[rightMaxPos];
		int rightSecondPos=Tools.secondHighestPosition(rightCounts);
		int rightSecond=rightCounts[rightSecondPos];

		if(verbose){
			outstream.println("kmer: "+toText(kmer, k)+", "+toText(rkmer, k));
			outstream.println("Right counts: "+count+", "+Arrays.toString(rightCounts));
			outstream.println("rightMaxPos="+rightMaxPos);
			outstream.println("rightMax="+rightMax);
			outstream.println("rightSecondPos="+rightSecondPos);
			outstream.println("rightSecond="+rightSecond);
		}

		if(rightMax<minCount){
			if(verbose){outstream.println("Returning DEAD_END: rightMax="+rightMax);}
			return DEAD_END;
		}else if(rightSecond>=minCount){
			if(verbose){outstream.println("Returning F_BRANCH: rightSecond="+rightSecond);}
			return F_BRANCH;
		}else if(rightMax>maxCount){
			if(verbose){outstream.println("Returning TOO_DEEP: rightMax="+rightMax);}
			return TOO_DEEP;
		}
		
		while(bb.length()<=maxLength){
			final int prevCount=count;
			
			//Generate the new base
			final byte b=AminoAcid.numberToBase[rightMaxPos];
			final long x=rightMaxPos;
			final long x2=AminoAcid.numberToComplement[rightMaxPos];
			kmer=((kmer<<2)|(long)x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			
			//Now consider the next kmer
			key=toValue(kmer, rkmer);
			if(key==firstKey){
				if(verbose){outstream.println("Returning LOOP");}
				return LOOP;
			}
			table=tables.getTableForKey(key);
			
			assert(table.getValue(key)==rightMax || rightMax==0);
			count=rightMax;
			
			{//Look right
				rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
				rightMax=rightCounts[rightMaxPos];
				rightSecondPos=Tools.secondHighestPosition(rightCounts);
				rightSecond=rightCounts[rightSecondPos];

				if(verbose){
					outstream.println("kmer: "+toText(kmer, k)+", "+toText(rkmer, k));
					outstream.println("Right counts: "+count+", "+Arrays.toString(rightCounts));
					outstream.println("rightMaxPos="+rightMaxPos);
					outstream.println("rightMax="+rightMax);
					outstream.println("rightSecondPos="+rightSecondPos);
					outstream.println("rightSecond="+rightSecond);
				}

				if(rightMax<minCount){
					if(verbose){outstream.println("Returning DEAD_END: rightMax="+rightMax);}
					return DEAD_END;
				}else if(rightSecond>=minCount){
					if(verbose){outstream.println("Returning F_BRANCH: rightSecond="+rightSecond);}
					return F_BRANCH;
				}
			}
			
			{//Look left
				int leftMaxPos=fillLeftCounts(kmer, rkmer, leftCounts, mask, shift2);
				int leftMax=leftCounts[leftMaxPos];
				int leftSecondPos=Tools.secondHighestPosition(leftCounts);
				int leftSecond=leftCounts[leftSecondPos];
				
//				assert(leftMax==1 || leftMax==0) : prevCount+" -> "+Arrays.toString(leftCounts)+", "+count+", "+Arrays.toString(rightCounts);
				
				if(verbose){
					outstream.println("Left counts: "+count+", "+Arrays.toString(leftCounts));
					outstream.println("leftMaxPos="+leftMaxPos);
					outstream.println("leftMax="+leftMax);
					outstream.println("leftSecondPos="+leftSecondPos);
					outstream.println("leftSecond="+leftSecond);
				}
				
				if(leftSecond>=minCount){//Backward branch
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
			
			if(count>maxCount){
				if(verbose){outstream.println("Returning TOO_DEEP: rightMax="+rightMax);}
				return TOO_DEEP;
			}
			
			bb.append(b);
			if(verbose){outstream.println("Added base "+(char)b);}
		}
		
		assert(bb.length()>maxLength);
		if(verbose){outstream.println("Returning TOO_LONG: length="+bb.length()+", rightMax="+rightMax);}
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
			super(id_, k);		
			shift=2*k;
			shift2=shift-2;
			mask=(shift>63 ? -1L : ~((-1L)<<shift));
		}
		
		@Override
		boolean processNextTable(Kmer kmer, Kmer temp){
			final int tnum=nextTable.getAndAdd(1);
			if(tnum>=tables.ways){return false;}
			final HashArray1D table=tables.getTable(tnum);
			final int max=table.arrayLength();
			if(startFromHighCounts){
				for(int cell=0; cell<max; cell++){
					int x=processCell_high(table, cell);
					deadEndsFoundT+=x;
				}
			}else{
				for(int cell=0; cell<max; cell++){
					int x=processCell_low(table, cell);
					deadEndsFoundT+=x;
				}
			}
			return true;
		}
		
		@Override
		boolean processNextVictims(Kmer kmer, Kmer temp){
			final int tnum=nextVictims.getAndAdd(1);
			if(tnum>=tables.ways){return false;}
			final HashArray1D table=tables.getTable(tnum);
			final HashForest forest=table.victims();
			final int max=forest.arrayLength();
			if(startFromHighCounts){
				for(int cell=0; cell<max; cell++){
					KmerNode kn=forest.getNode(cell);
					int x=traverseKmerNode_high(kn);
					deadEndsFoundT+=x;
				}
			}else{
				for(int cell=0; cell<max; cell++){
					KmerNode kn=forest.getNode(cell);
					int x=traverseKmerNode_low(kn);
					deadEndsFoundT+=x;
				}
			}
			return true;
		}

		
		/**
		 * Recursively traverses k-mer nodes in high-count processing mode.
		 * Visits the current node and recursively processes left and right children.
		 * @param kn K-mer node to traverse (may be null)
		 * @return Total number of dead ends found in this subtree
		 */
		private int traverseKmerNode_high(KmerNode kn){
			int sum=0;
			if(kn!=null){
				sum+=processKmerNode_high(kn);
				if(kn.left()!=null){
					sum+=traverseKmerNode_high(kn.left());
				}
				if(kn.right()!=null){
					sum+=traverseKmerNode_high(kn.right());
				}
			}
			return sum;
		}

		
		/**
		 * Recursively traverses k-mer nodes in low-count processing mode.
		 * Visits the current node and recursively processes left and right children.
		 * @param kn K-mer node to traverse (may be null)
		 * @return Total number of dead ends found in this subtree
		 */
		private int traverseKmerNode_low(KmerNode kn){
			int sum=0;
			if(kn!=null){
				sum+=processKmerNode_low(kn);
				if(kn.left()!=null){
					sum+=traverseKmerNode_low(kn.left());
				}
				if(kn.right()!=null){
					sum+=traverseKmerNode_low(kn.right());
				}
			}
			return sum;
		}
		
		/*--------------------------------------------------------------*/
		
		//old
		/**
		 * Processes a single hash table cell in low-count mode.
		 * Checks if the k-mer meets count criteria and ownership requirements
		 * before attempting dead-end exploration.
		 *
		 * @param table Hash table containing the cell
		 * @param cell Cell index within the table
		 * @return 1 if a dead end was found and marked, 0 otherwise
		 */
		private int processCell_low(HashArray1D table, int cell){
			int count=table.readCellValue(cell);
			if(count<minSeed || count>maxCount){return 0;}
			int owner=table.getCellOwner(cell);
			if(owner>STATUS_UNEXPLORED){return 0;}
			long key=table.getKmer(cell);
			if(verbose){outstream.println("id="+id+" processing cell "+cell+"; \tkmer="+key+"\t"+AminoAcid.kmerToString(key, k));}
			
			return processKmer_low(key);
		}
		
		//old
		/**
		 * Processes a k-mer node in low-count mode.
		 * Extracts the k-mer and count from the node, validates criteria,
		 * then attempts dead-end exploration.
		 *
		 * @param kn K-mer node to process
		 * @return 1 if a dead end was found and marked, 0 otherwise
		 */
		private int processKmerNode_low(KmerNode kn){
			final long key=kn.pivot();
			final int count=kn.getValue(key);
			if(count<minSeed || count>maxCount){return 0;}
			int owner=kn.getOwner(key);
			if(owner>STATUS_UNEXPLORED){return 0;}
			
			return processKmer_low(key);
		}
		
		//old
		/**
		 * Processes a single k-mer in low-count mode.
		 * Calls exploreAndMark to determine if the k-mer is part of a removable dead end.
		 * @param key K-mer to process
		 * @return 1 if a dead end was found and marked, 0 otherwise
		 */
		private int processKmer_low(long key){
			kmersTestedT++;
			boolean b=exploreAndMark(key, builderT, leftCounts, rightCounts, minCount, maxCount, maxLengthToDiscard, maxDistanceToExplore, true
					, countMatrixT, removeMatrixT
				);
			return b ? 1 : 0;
		}
		
		/*--------------------------------------------------------------*/
		
		//new
		/**
		 * Processes a single hash table cell in high-count mode.
		 * Focuses on high-coverage k-mers that may have low-coverage neighbors
		 * representing dead ends extending from repetitive sequences.
		 *
		 * @param table Hash table containing the cell
		 * @param cell Cell index within the table
		 * @return Number of dead ends found from this high-coverage k-mer
		 */
		private int processCell_high(HashArray1D table, int cell){
			int count=table.readCellValue(cell);
			if(count<=maxCount){return 0;}
			if(shaveFast && maxCount==1 && count>=6){return 0;}
//			int owner=table.getCellOwner(cell);
//			if(owner>STATUS_UNEXPLORED){return 0;}
			long key=table.getKmer(cell);
			if(verbose){outstream.println("id="+id+" processing cell "+cell+"; \tkmer="+key+"\t"+AminoAcid.kmerToString(key, k));}
			
			return processKmer_high(key);
		}
		
		//new
		/**
		 * Processes a k-mer node in high-count mode.
		 * Looks for high-coverage k-mers that may have low-coverage dead-end neighbors.
		 * @param kn K-mer node to process
		 * @return Number of dead ends found from this high-coverage k-mer
		 */
		private int processKmerNode_high(KmerNode kn){
			final long key=kn.pivot();
			final int count=kn.getValue(key);
			if(count<=maxCount){return 0;}
			if(shaveFast && maxCount==1 && count>=6){return 0;}
//			int owner=kn.getOwner(key);
//			if(owner>STATUS_UNEXPLORED){return 0;}
			
			return processKmer_high(key);
		}
		
//		//new
//		private int processKmer_high(final long key0){
//			if(!tables.MASK_CORE){return processKmer_high_safe(key0);}
//			final int shift=2*k;
//			final int shift2=shift-2;
//			final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
//			int sum=0;
//			final long kmer0=key0;
//			final long rkmer0=AminoAcid.reverseComplementBinaryFast(key0, k);
////			fillRightCounts(kmer0, kmer0, rightCounts, mask, shift2);
//			{
//				final long kmer1=(kmer0<<2)&mask;
//				final long rkmer1=(rkmer0>>2);
//				if((kmer1&tables.coreMask)==(rkmer1&tables.coreMask)){return processKmer_high_safe(key0);}
//				final long key1=toValue(kmer1|3, rkmer1);
//				final int way=tables.kmerToWay(key1);
//				final HashArray1D table=(HashArray1D) tables.tables()[way];
//				final int cell=table.kmerToCell(key1);
//				for(long i=0, j=3; i<4; i++, j--){
//					long kmer=kmer1|i;
//					long rkmer=rkmer1|(j<<shift2);
////					long kmer=((kmer0<<2)|i)&mask;
////					long rkmer=(rkmer0>>2)|(j<<shift2);
//					long key=toValue(kmer,  rkmer);
////					int way=tables.kmerToWay(key);
////					HashArray1D table=(HashArray1D) tables.tables()[way];
//					int count=table.getValue(key, cell);
//					if(count>0 && count<=maxCount && table.getOwner(key)<=STATUS_UNEXPLORED){
//						kmersTestedT++;
//						boolean b=exploreAndMark(key, builderT, leftCounts, rightCounts,
//								minCount, maxCount, maxLengthToDiscard, maxDistanceToExplore, true,
//								countMatrixT, removeMatrixT);
//						if(b){sum++;}
//					}
//				}
//			}
//			{
//				final long rkmer1=(rkmer0<<2)&mask;
//				final long kmer1=(kmer0>>2);
//				if((kmer1&tables.coreMask)==(rkmer1&tables.coreMask)){return processKmer_high_safe(key0);}
//				final long key1=toValue(kmer1, rkmer1|3);
//				final int way=tables.kmerToWay(key1);
//				final HashArray1D table=(HashArray1D) tables.tables()[way];
//				final int cell=table.kmerToCell(key1);
//				for(long i=0, j=3; i<4; i++, j--){
//					long rkmer=rkmer1|i;
//					long kmer=kmer1|(j<<shift2);
////					long rkmer=((rkmer0<<2)|i)&mask;
////					long kmer=(kmer0>>2)|(j<<shift2);
//					long key=toValue(kmer,  rkmer);
////					int way=tables.kmerToWay(key);
////					HashArray1D table=(HashArray1D) tables.tables()[way];
//					int count=table.getValue(key, cell);
//					if(count>0 && count<=maxCount && table.getOwner(key)<=STATUS_UNEXPLORED){
//						kmersTestedT++;
//						boolean b=exploreAndMark(key, builderT, leftCounts, rightCounts,
//								minCount, maxCount, maxLengthToDiscard, maxDistanceToExplore, true,
//								countMatrixT, removeMatrixT);
//						if(b){sum++;}
//					}
//				}
//			}
//			return sum;
//		}
		
		//new
		private int processKmer_high(final long key0){
			if(!tables.MASK_CORE){return processKmer_high_safe(key0);}
			final long kmer0=key0;
			final long rkmer0=AminoAcid.reverseComplementBinaryFast(key0, k);
			int sum=0;
			sum+=processKmer_high_leftLoop(kmer0, rkmer0);
			sum+=processKmer_high_leftLoop(rkmer0, kmer0);
			return sum;
		}
		
		private int processKmer_high_leftLoop(final long kmer0, final long rkmer0){
			int sum=0;
			final long rkmer1=(rkmer0<<2)&mask;
			final long kmer1=(kmer0>>2);

			if((kmer1&tables.coreMask)==(rkmer1&tables.coreMask)){return processKmer_high_safe(kmer0);}
			final long key1=toValue(kmer1, rkmer1|3);
			final int way=tables.kmerToWay(key1);
			final HashArray1D table=(HashArray1D) tables.tables()[way];
			final int cell=table.kmerToCell(key1);
			
			for(long i=0, j=3; i<4; i++, j--){
				long rkmer=rkmer1|i;
				long kmer=kmer1|(j<<shift2);
				//				long rkmer=((rkmer0<<2)|i)&mask;
				//				long kmer=(kmer0>>2)|(j<<shift2);
				long key=toValue(kmer, rkmer);
				//				int way=tables.kmerToWay(key);
				//				HashArray1D table=(HashArray1D) tables.tables()[way];
				int count=table.getValue(key, cell);
				if(count>0 && count<=maxCount && table.getOwner(key)<=STATUS_UNEXPLORED){
					kmersTestedT++;
					boolean b=exploreAndMark(key, builderT, leftCounts, rightCounts,
							minCount, maxCount, maxLengthToDiscard, maxDistanceToExplore, true,
							countMatrixT, removeMatrixT);
					if(b){sum++;}
				}
			}
			return sum;
		}
		
//		//new
//		private int processKmer_high_safe(final long key0){
//			final int shift=2*k;
//			final int shift2=shift-2;
//			final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
//			int sum=0;
//			final long kmer0=key0;
//			final long rkmer0=AminoAcid.reverseComplementBinaryFast(key0, k);
//			{
//				final long kmer1=(kmer0<<2)&mask;
//				final long rkmer1=(rkmer0>>2);
//				for(long i=0, j=3; i<4; i++, j--){
//					long kmer=kmer1|i;
//					long rkmer=rkmer1|(j<<shift2);
//					long key=toValue(kmer, rkmer);
//					int way=tables.kmerToWay(key);
//					HashArray1D table=(HashArray1D) tables.tables()[way];
//					int count=table.getValue(key);
//					if(count>0 && count<=maxCount && table.getOwner(key)<=STATUS_UNEXPLORED){
//						kmersTestedT++;
//						boolean b=exploreAndMark(key, builderT, leftCounts, rightCounts,
//								minCount, maxCount, maxLengthToDiscard, maxDistanceToExplore, true,
//								countMatrixT, removeMatrixT);
//						if(b){sum++;}
//					}
//				}
//			}
//			{
//				final long rkmer1=(rkmer0<<2)&mask;
//				final long kmer1=(kmer0>>2);
//				for(long i=0, j=3; i<4; i++, j--){
//					long rkmer=rkmer1|i;
//					long kmer=kmer1|(j<<shift2);
//					long key=toValue(kmer, rkmer);
//					int way=tables.kmerToWay(key);
//					HashArray1D table=(HashArray1D) tables.tables()[way];
//					int count=table.getValue(key);
//					if(count>0 && count<=maxCount && table.getOwner(key)<=STATUS_UNEXPLORED){
//						kmersTestedT++;
//						boolean b=exploreAndMark(key, builderT, leftCounts, rightCounts,
//								minCount, maxCount, maxLengthToDiscard, maxDistanceToExplore, true,
//								countMatrixT, removeMatrixT);
//						if(b){sum++;}
//					}
//				}
//			}
//			return sum;
//		}
		
		//new
		private int processKmer_high_safe(final long key0){
			int sum=0;
			final long kmer0=key0;
			final long rkmer0=AminoAcid.reverseComplementBinaryFast(key0, k);
			sum+=processKmer_high_safe_leftLoop(kmer0, rkmer0);
			sum+=processKmer_high_safe_leftLoop(rkmer0, kmer0);
			return sum;
		}
		
		private int processKmer_high_safe_leftLoop(final long kmer0, final long rkmer0){
			int sum=0;
			final long kmer1=(kmer0<<2)&mask;
			final long rkmer1=(rkmer0>>2);
			for(long i=0, j=3; i<4; i++, j--){
				final long kmer=kmer1|i;
				final long rkmer=rkmer1|(j<<shift2);
				final long key=toValue(kmer, rkmer);
				final int way=tables.kmerToWay(key);
				HashArray1D table=(HashArray1D) tables.tables()[way];
				final int count=table.getValue(key);
				if(count>0 && count<=maxCount && table.getOwner(key)<=STATUS_UNEXPLORED){
					kmersTestedT++;
					boolean b=exploreAndMark(key, builderT, leftCounts, rightCounts,
							minCount, maxCount, maxLengthToDiscard, maxDistanceToExplore, true,
							countMatrixT, removeMatrixT);
					if(b){sum++;}
				}
			}
			return sum;
		}
		
		/** Bit shift amount for k-mer rolling (2*k bits) */
		final int shift;
		/** Secondary bit shift amount for complement rolling (shift-2 bits) */
		final int shift2;
		/** Bit mask for k-mer length limitation */
		final long mask;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          ShaveThread         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Removes dead-end kmers.
	 */
	class ShaveThread extends AbstractShaveThread{

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
			final HashArray1D table=tables.getTable(tnum);
			final AtomicIntegerArray owners=table.owners();
			final int[] values=table.values();
			final int max=table.arrayLength();
			for(int cell=0; cell<max; cell++){
				if(owners.get(cell)==STATUS_REMOVE){
//					x++;
					values[cell]=0;
				}
			}
			for(KmerNode kn : table.victims().array()){
				if(kn!=null){traverseKmerNode(kn);}
			}
			
			table.clearOwnership();
			kmersRemovedT+=table.regenerate(0);
//			outstream.println(x);
			return true;
		}
		
		/**
		 * Recursively traverses k-mer nodes to remove marked entries.
		 * Sets the count to 0 for nodes marked with STATUS_REMOVE
		 * and recursively processes child nodes.
		 * @param kn K-mer node to process (may be null)
		 */
		private void traverseKmerNode(KmerNode kn){
			if(kn==null){return;}
			if(kn.owner()==STATUS_REMOVE){kn.set(0);}
			traverseKmerNode(kn.left());
			traverseKmerNode(kn.right());
		}
		
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Recall Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Converts k-mer and its reverse complement to a canonical value.
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @return Canonical representation for table lookup
	 */
	private final long toValue(long kmer, long rkmer){return tables.toValue(kmer, rkmer);}
	/**
	 * Gets the count for a k-mer pair.
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @return Count value from the k-mer tables
	 */
	int getCount(long kmer, long rkmer){return tables.getCount(kmer, rkmer);}
	/**
	 * Claims ownership of a k-mer.
	 * @param kmer K-mer to claim
	 * @param id Claiming thread/process ID
	 * @return true if successfully claimed, false if already owned
	 */
	boolean claim(long kmer, int id){return claim(kmer, AminoAcid.reverseComplementBinaryFast(kmer, k), id);}
	/**
	 * Claims ownership of a k-mer pair.
	 *
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @param id Claiming thread/process ID
	 * @return true if successfully claimed, false if already owned
	 */
	boolean claim(long kmer, long rkmer, int id){return tables.claim(kmer, rkmer, id);}
	boolean doubleClaim(ByteBuilder bb, int id/*, long rid*/){return tables.doubleClaim(bb, id/*, rid*/);}
	boolean claim(ByteBuilder bb, int id, /*long rid, */boolean earlyExit){return tables.claim(bb, id/*, rid*/, earlyExit);}
	boolean claim(byte[] array, int len, int id, /*long rid, */boolean earlyExit){return tables.claim(array, len, id/*, rid*/, earlyExit);}
	/**
	 * Finds the current owner of a k-mer.
	 * @param kmer K-mer to check
	 * @return Owner ID or status code
	 */
	int findOwner(long kmer){return tables.findOwner(kmer);}
	int findOwner(ByteBuilder bb, int id){return tables.findOwner(bb, id);}
	int findOwner(byte[] array, int len, int id){return tables.findOwner(array, len, id);}
	void release(ByteBuilder bb, int id){tables.release(bb, id);}
	void release(byte[] array, int len, int id){tables.release(array, len, id);}
	/**
	 * Fills array with counts of right-extension k-mers.
	 *
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @param counts Array to fill with extension counts
	 * @param mask Bit mask for k-mer operations
	 * @param shift2 Bit shift amount for reverse complement
	 * @return Position of maximum count
	 */
	int fillRightCounts(long kmer, long rkmer, int[] counts, long mask, int shift2){return tables.fillRightCounts(kmer, rkmer, counts, mask, shift2);}
	/**
	 * Fills array with counts of left-extension k-mers.
	 *
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @param counts Array to fill with extension counts
	 * @param mask Bit mask for k-mer operations
	 * @param shift2 Bit shift amount for reverse complement
	 * @return Position of maximum count
	 */
	int fillLeftCounts(long kmer, long rkmer, int[] counts, long mask, int shift2){return tables.fillLeftCounts(kmer, rkmer, counts, mask, shift2);}
	/**
	 * Converts a k-mer to text representation.
	 * @param kmer K-mer to convert
	 * @param k Length of k-mer
	 * @return Text representation of the k-mer sequence
	 */
	static StringBuilder toText(long kmer, int k){return AbstractKmerTable.toText(kmer, k);}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	AbstractKmerTableSet tables(){return tables;}
	
	/** The k-mer table set containing sequence data for shaving operations */
	final KmerTableSet tables;
	/** Length of k-mers used in this shaver instance */
	final int k;
	
}
