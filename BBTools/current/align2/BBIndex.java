package align2;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import dna.AminoAcid;
import dna.Data;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import stream.SiteScore;
import structures.LongM;


/**
 * Based on Index11f
 * Index stored in single array per block.
 * 
 * 
 * 
 * @author Brian Bushnell
 * @date Dec 22, 2012
 *
 */
public final class BBIndex extends AbstractIndex {
	
	
	public static void main(String[] args){
		
		int k=13;
		
		for(int i=0; i<args.length; i++){
			String s=args[i].toLowerCase();
			if(s.contains("=")){
				String[] split=s.split("=");
				String a=split[0];
				String b=split[1];
				if(a.equals("build") || a.equals("b")){
					Data.setGenome(Integer.parseInt(b));
				}else if(a.equals("minchrom")){
					MINCHROM=Integer.parseInt(b);
				}else if(a.equals("maxchrom")){
					MAXCHROM=Integer.parseInt(b);
				}else if(a.equals("keylen") || a.equals("k")){
					k=Integer.parseInt(b);
				}
			}
		}

		if(MINCHROM==-1){MINCHROM=1;}
		if(MAXCHROM==-1){
			assert(Data.numChroms<=Byte.MAX_VALUE) : "TODO";
			MAXCHROM=Data.numChroms;
		}
		
		
		System.err.println("Writing build "+Data.GENOME_BUILD+" "+
				"BASESPACE index, keylen="+k+", chrom bits="+NUM_CHROM_BITS);
		
		
		int first=(NUM_CHROM_BITS==0 ? 1 : 0);
		
		
		Data.sysout.println("Loading index for chunk "+first+"-"+MAXCHROM+", build "+Data.GENOME_BUILD);
		index=IndexMaker4.makeIndex(Data.GENOME_BUILD, first, MAXCHROM,
				k, NUM_CHROM_BITS, MAX_ALLOWED_CHROM_INDEX, CHROM_MASK_LOW, CHROM_MASK_HIGH, SITE_MASK, SHIFT_LENGTH, true, false, index);
		
		
		System.err.println("Finished all chroms, may still be writing.");
	}
	
	
	public BBIndex(int k_, int minChrom_, int maxChrom_, int kfilter_, MSA msa_){
		super(k_, kfilter_, BASE_HIT_SCORE, minChrom_, maxChrom_, msa_);
		INV_BASE_KEY_HIT_SCORE=1f/BASE_KEY_HIT_SCORE;
		INDEL_PENALTY=(BASE_KEY_HIT_SCORE/2)-1; //default (HIT_SCORE/2)-1
		INDEL_PENALTY_MULT=20; //default 20; penalty for indel length
		MAX_PENALTY_FOR_MISALIGNED_HIT=BASE_KEY_HIT_SCORE-(1+BASE_KEY_HIT_SCORE/8);
		SCOREZ_1KEY=Z_SCORE_MULT*KEYLEN;
		{
			int cyc=0;
			for(int chrom=minChrom; chrom<=maxChrom; chrom=((chrom&CHROM_MASK_HIGH)+CHROMS_PER_BLOCK)){cyc+=2;}
			cycles=cyc;
		}
		prescoreArray=new int[cycles];
		precountArray=new int[cycles];
	}
	
	/** Load or generate index from minChrom to maxChrom, inclusive, with keylength k.
	 * This range can encompass multiple blocks.
	 * Should only be called once in a process. */
	public static final synchronized void loadIndex(int minChrom, int maxChrom, int k, boolean writeToDisk, boolean diskInvalid){
		if(minChrom<1){minChrom=1;}
		if(maxChrom>Data.numChroms){maxChrom=Data.numChroms;}
		assert(minChrom<=maxChrom) : minChrom+", "+maxChrom+"\nThe reference file appears to be empty.";
		Data.sysout.println("Loading index for chunk "+minChrom+"-"+maxChrom+", build "+Data.GENOME_BUILD);
		index=IndexMaker4.makeIndex(Data.GENOME_BUILD, minChrom, maxChrom,
				k, NUM_CHROM_BITS, MAX_ALLOWED_CHROM_INDEX, CHROM_MASK_LOW, CHROM_MASK_HIGH, SITE_MASK, SHIFT_LENGTH, writeToDisk, diskInvalid, index);
	}
	
	/** Calculate statistics of index, such as list lengths, and find clumpy keys */
	public static final synchronized void analyzeIndex(int minChrom, int maxChrom, float fractionToExclude, int k){
		assert(lengthHistogram==null);
		assert(COUNTS==null);
		
		int KEYSPACE=1<<(2*k);
		COUNTS=new int[KEYSPACE];
		maxChrom=maxChrom(maxChrom);
		
		HashMap<Integer, LongM> cmap=new HashMap<Integer, LongM>();
		
		for(int chrom=minChrom; chrom<=maxChrom; chrom=((chrom&CHROM_MASK_HIGH)+CHROMS_PER_BLOCK)){
			Block b=index[chrom];
			final int[] sites=b.sites;
			final int[] starts=b.starts;

			for(int key=0; key<KEYSPACE; key++){

				long clumps=0;

				final int start1=starts[key];
				final int stop1=starts[key+1];
				final int len1=stop1-start1;
				COUNTS[key]=Tools.min(Integer.MAX_VALUE, COUNTS[key]+len1);

				if(REMOVE_CLUMPY){
					for(int i=start1+1; i<stop1; i++){
						int dif=sites[i]-sites[i-1];
						assert(dif!=0);
						if(dif>0 && dif<=CLUMPY_MAX_DIST){
							clumps++;
						}
					}
					if(clumps>0){
						final int x=Tools.min(key, AminoAcid.reverseComplementBinaryFast(key, k));
						final Integer ko=x;
						LongM lm=cmap.get(ko);
						if(lm==null){
							lm=new LongM(0);
							cmap.put(ko, lm);
						}
						lm.increment(clumps);
					}
				}
			}
		}
		
		for(int key=0; key<COUNTS.length; key++){
			int rkey=AminoAcid.reverseComplementBinaryFast(key, k);
			if(key<rkey){
				int x=(int)Tools.min(Integer.MAX_VALUE, COUNTS[key]+(long)COUNTS[rkey]);
				COUNTS[key]=COUNTS[rkey]=x;
			}
		}
		
		if(REMOVE_CLUMPY){
			Integer[] keys=cmap.keySet().toArray(new Integer[cmap.size()]);
			Arrays.sort(keys);
			
			for(Integer key : keys){
				long clumps=cmap.get(key).value();
				long len=COUNTS[key];
				if((len>CLUMPY_MIN_LENGTH_INDEX && clumps>CLUMPY_FRACTION*len)/* || (len>8*CLUMPY_MIN_LENGTH_INDEX && clumps>.75f*CLUMPY_FRACTION*len)*/){
					int rkey=AminoAcid.reverseComplementBinaryFast(key, k);
					assert(key<=rkey);
					assert(key==KeyRing.reverseComplementKey(rkey, k));
					COUNTS[key]=0;
					COUNTS[rkey]=0;
				}
			}
		}
		
		lengthHistogram=Tools.makeLengthHistogram3(COUNTS, 1000, verbose2);
		
		//if(verbose2){System.err.println("lengthHistogram: "+Arrays.toString(lengthHistogram));}
		
		if(REMOVE_FREQUENT_GENOME_FRACTION){

			int lengthLimitIndex=(int)((1-fractionToExclude)*(lengthHistogram.length-1));
			int lengthLimitIndex2=(int)((1-fractionToExclude*DOUBLE_SEARCH_THRESH_MULT)*(lengthHistogram.length-1));
			
			MAX_USABLE_LENGTH=Tools.max(2*SMALL_GENOME_LIST, lengthHistogram[lengthLimitIndex]);
			MAX_USABLE_LENGTH2=Tools.max(6*SMALL_GENOME_LIST, lengthHistogram[lengthLimitIndex2]);
			
			if(verbose2){System.err.println("MAX_USABLE_LENGTH:  "+MAX_USABLE_LENGTH+"\nMAX_USABLE_LENGTH2: "+MAX_USABLE_LENGTH2);}
		}
		
		Solver.POINTS_PER_SITE=(int)Math.floor((Solver.BASE_POINTS_PER_SITE*4000f)/Tools.max(2*SMALL_GENOME_LIST, lengthHistogram[MAX_AVERAGE_LIST_TO_SEARCH]));
		if(Solver.POINTS_PER_SITE==0){Solver.POINTS_PER_SITE=-1;}
		if(verbose2){System.err.println("POINTS_PER_SITE:  "+Solver.POINTS_PER_SITE);}
		assert(Solver.POINTS_PER_SITE<0) : Solver.POINTS_PER_SITE;
	}
	
	
	/** Returns the filename for the block holding this chrom */
	public static final String fname(int chrom, int k){
		return IndexMaker4.fname(minChrom(chrom), maxChrom(chrom), k, NUM_CHROM_BITS);
	}
	
	/** Ensure key offsets are strictly ascending. */
	private static boolean checkOffsets(int[] offsets){
		for(int i=1; i<offsets.length; i++){
			if(offsets[i]<=offsets[i-1]){return false;}
		}
		return true;
	}
	
	@Deprecated
	private final int trimExcessHitLists(int[] keys, int[][] hits){
		
		assert(false) : "Needs to be redone because hits are no longer sorted by length.";
		
		assert(hits.length==keys.length);
//		assert(false) : "modify this function so that it gives more weight to trimming lists over highly covered baits";
		//And also, incorporate the "remove the longest list" function
		
		final int limit=Tools.max(SMALL_GENOME_LIST, lengthHistogram[MAX_AVERAGE_LIST_TO_SEARCH])*keys.length;
		final int limit2=Tools.max(SMALL_GENOME_LIST, lengthHistogram[MAX_AVERAGE_LIST_TO_SEARCH2]);
		final int limit3=Tools.max(SMALL_GENOME_LIST, lengthHistogram[MAX_SHORTEST_LIST_TO_SEARCH]);
		
		int sum=0;
		int initialHitCount=0;

		int shortest=Integer.MAX_VALUE-1;
		int shortest2=Integer.MAX_VALUE;
		
		for(int i=0; i<keys.length; i++){
			int key=keys[i];
			int x=count(key);
			sum+=x;
			initialHitCount+=(x==0 ? 0 : 1);
			if(x>0){
				if(x<shortest2){
					shortest2=x;
					if(shortest2<shortest){
						shortest2=shortest;
						shortest=x;
					}
				}
			}
		}
		assert(shortest2>=shortest);
		if(initialHitCount<MIN_APPROX_HITS_TO_KEEP){return initialHitCount;}
		if(shortest>limit3 && !SLOW){
			for(int i=0; i<hits.length; i++){hits[i]=null;}
			return 0;
		}
		if(sum<=limit && sum/initialHitCount<=limit2){return initialHitCount;}
		
		Pointer[] ptrs=Pointer.loadMatrix(hits);
//		ptrs[0].value/=2;
//		ptrs[ptrs.length-1].value/=2;
		Arrays.sort(ptrs);
		
		int finalHitCount=initialHitCount;
		for(int i=ptrs.length-1; sum>limit || sum/finalHitCount>limit2; i--){
			Pointer p=ptrs[i];
			sum-=hits[p.key].length;
			hits[p.key]=null;
			finalHitCount--;
		}
		
		return finalHitCount;
	}
	
	/** Remove least useful keys to accelerate search */
	private final int trimExcessHitListsByGreedy(int[] offsets, int[] keyScores, int maxHitLists, int[] keys){
		
		float[] keyWeights=getKeyWeightArray(keyScores.length);
		for(int i=0; i<keyScores.length; i++){
			keyWeights[i]=keyScores[i]*INV_BASE_KEY_HIT_SCORE;
		}
		
//		assert(false) : "modify this function so that it gives more weight to trimming lists over highly covered baits";
		//And also, incorporate the "remove the longest list" function
		
		//Remove lists until the length sum is below this value
		final int limit=Tools.max(SMALL_GENOME_LIST, lengthHistogram[MAX_AVERAGE_LIST_TO_SEARCH])*keys.length;
		final int limit2=Tools.max(SMALL_GENOME_LIST, lengthHistogram[MAX_AVERAGE_LIST_TO_SEARCH2]);
		
		//Do not search if the shortest list is longer than this
		final int limit3=Tools.max(SMALL_GENOME_LIST, lengthHistogram[MAX_SHORTEST_LIST_TO_SEARCH]);
//		final int limitS=lengthHistogram[chrom][MAX_SINGLE_LIST_TO_SEARCH];
		
		int sum=0;
		int initialHitCount=0;
		
		//Shortest list found
		int shortest=Integer.MAX_VALUE-1;
		//Second shortest list found
		int shortest2=Integer.MAX_VALUE;
		
//		for(int i=0; i<hits.length; i++){
//			if(hits[i]!=null && hits[i].length>limitS){hits[i]=null;}
//		}
		
		final int[] lengths=getGenericArray(keys.length);
		
		//Find shortest lists
		for(int i=0; i<keys.length; i++){
			int key=keys[i];
			int x=count(key);
			lengths[i]=x;
			sum+=x;
			initialHitCount+=(x==0 ? 0 : 1);
			if(x>0){
				if(x<shortest2){
					shortest2=x;
					if(shortest2<shortest){
						shortest2=shortest;
						shortest=x;
					}
				}
			}
		}
//		assert(false) : limit+", "+limit2+", "+limit3+", "+shortest2+", "+shortest+", "+initialHitCount+", "+MIN_APPROX_HITS_TO_KEEP+"\n"+Arrays.toString(lengths);
		assert(shortest2>=shortest);
		if(initialHitCount<MIN_APPROX_HITS_TO_KEEP){return initialHitCount;}
		
		//Don't search if the shortest list is too long
		if(shortest>limit3 && !SLOW){
			for(int i=0; i<keys.length; i++){keys[i]=-1;}
			return 0;
		}
		
		int hitsCount=initialHitCount;
		int worstValue=Integer.MIN_VALUE;
		
		while(hitsCount>=MIN_APPROX_HITS_TO_KEEP && (sum>limit || sum/initialHitCount>limit2 || hitsCount>maxHitLists/* || worstValue<0*/)){
			final int[] lists=getGreedyListArray(hitsCount);
			for(int i=0, j=0; j<lists.length; i++){
				if(lengths[i]>0){
					lists[j]=i;
					j++;
				}
			}
			
			Solver.findWorstGreedy(offsets, lengths, keyWeights, KEYLEN, lists, greedyReturn);
			int worstIndex=greedyReturn[0];
			int worst=lists[worstIndex];
			worstValue=greedyReturn[1];
			sum-=lengths[worst];

//			if(worstValue>0 && (hitsCount<=maxHitLists || lengths[worst]<excessListLimit)){return hitsCount;}
			if(worstValue>0 || lengths[worst]<SMALL_GENOME_LIST){return hitsCount;} //This line increases accuracy at expense of speed.  Lower constant = more accurate, default 0.
			hitsCount--;
			lengths[worst]=0;
			keys[worst]=-1;
		}
		return hitsCount;
	}
	
	
	private final int getHits(final int[] keys, final int chrom, final int maxLen, final int[] starts, final int[] stops){
		int numHits=0;
		final Block b=index[chrom];
		for(int i=0; i<keys.length; i++){
			final int key=keys[i];
			starts[i]=-1;
			stops[i]=-1;
			if(key>=0){
				final int len=count(key);
				if(len>0 && len<maxLen){
					final int len2=b.length(key);
					if(len2>0){
						starts[i]=b.starts[key];
						stops[i]=starts[i]+len2;
						numHits++;
					}
				}
			}
		}
		return numHits;
	}
	
	
	private final int countHits(final int[] keys, final int maxLen, boolean clearBadKeys){
		int numHits=0;
		for(int i=0; i<keys.length; i++){
			final int key=keys[i];
			if(key>=0){
				final int len=count(key);
//				System.err.println(len);
				if(len>0 && len<maxLen){
					numHits++;
				}else if(clearBadKeys){
					keys[i]=-1;
				}
			}
		}
		return numHits;
	}
	
	
	@Override
	public final ArrayList<SiteScore> findAdvanced(byte[] basesP, byte[] basesM, byte[] qual, byte[] baseScoresP, int[] keyScoresP, int[] offsets, long id){
		assert(minChrom<=maxChrom && minChrom>=0);
		ArrayList<SiteScore> result=find(basesP, basesM, qual, baseScoresP, keyScoresP, offsets, true, id);
		if(DOUBLE_SEARCH_NO_HIT && (result==null || result.isEmpty())){result=find(basesP, basesM, qual, baseScoresP, keyScoresP, offsets, false, id);}
		
		return result;
	}
	
	
	public final ArrayList<SiteScore> find(byte[] basesP, byte[] basesM, byte[] qual,  byte[] baseScoresP, int[] keyScoresP, int[] offsetsP, boolean obeyLimits, long id){
		
		assert(checkOffsets(offsetsP)) : Arrays.toString(offsetsP);
		final int[] keysOriginal=KeyRing.makeKeys(basesP, offsetsP, KEYLEN);
		int[] keysP=Arrays.copyOf(keysOriginal, keysOriginal.length);

		initialKeys+=offsetsP.length;
		initialKeyIterations++;
		
		final int maxLen=(obeyLimits ? MAX_USABLE_LENGTH : MAX_USABLE_LENGTH2);
		
		int numHits=0;
		numHits=countHits(keysP, maxLen, true);
		
		if(verbose){
			System.err.println("initial hits: "+numHits);
			System.err.println("initial keys: "+keysP.length+"\n"+Arrays.toString(keysP));
		}
		
		if(numHits>0){ //TODO: Change these to higher numbers
			int trigger=(3*keysP.length)/4;
			if(numHits<4 && numHits<trigger){
				for(int i=0; i<keysP.length; i++){keysP[i]=keysOriginal[i];}
				numHits=countHits(keysP, (maxLen*3)/2, true);
			}
			if(numHits<3 && numHits<trigger){
				for(int i=0; i<keysP.length; i++){keysP[i]=keysOriginal[i];}
				numHits=countHits(keysP, maxLen*2, true);
			}
			if(numHits<3 && numHits<trigger){
				for(int i=0; i<keysP.length; i++){keysP[i]=keysOriginal[i];}
				numHits=countHits(keysP, maxLen*3, true);
			}
			if(numHits<2 && numHits<trigger){
				for(int i=0; i<keysP.length; i++){keysP[i]=keysOriginal[i];}
				numHits=countHits(keysP, maxLen*5, true);
			}
		}
//		assert(checkOffsets(offsetsP)) : Arrays.toString(offsetsP);
		if(numHits<keysP.length){
			int[][] r=shrink2(offsetsP, keysP, keyScoresP);
			assert(r!=null);
			if(r!=null){
				offsetsP=r[0];
				keysP=r[1];
				keyScoresP=r[2];
			}
		}else{
			assert(shrink2(offsetsP, keysP, keyScoresP)==null);
		}
		initialKeys2+=numHits;
		//assert(checkOffsets(offsetsP)) : Arrays.toString(offsetsP);
		
//		assert(checkOffsets(offsets)) : Arrays.toString(offsets);
		if(TRIM_BY_GREEDY && obeyLimits){
			int maxLists=Tools.max((int)(HIT_FRACTION_TO_RETAIN*keysP.length), MIN_HIT_LISTS_TO_RETAIN);
			numHits=trimExcessHitListsByGreedy(offsetsP, keyScoresP, maxLists, keysP);
		}
//		System.out.println("After greedy: numHits = "+numHits);
		
		if(verbose){
			System.err.println("final hits: "+numHits);
			System.err.println("final keys: "+keysP.length+"\n"+Arrays.toString(keysP));
		}
		
		if(TRIM_BY_TOTAL_SITE_COUNT && obeyLimits){
			throw new RuntimeException("Needs to be redone.");
//			numHits=trimExcessHitLists(keys, hits);
		}
		
		if(TRIM_LONG_HIT_LISTS && obeyLimits && numHits>MIN_APPROX_HITS_TO_KEEP){
			int cutoffIndex=((int) (HIT_FRACTION_TO_RETAIN*(keysP.length)-0.01f))+(keysP.length-numHits);

			int zeroes=keysP.length-numHits;
			int altMinIndex=(zeroes+(MIN_HIT_LISTS_TO_RETAIN-1));
			cutoffIndex=Tools.max(cutoffIndex, altMinIndex);
			
			assert(cutoffIndex>0) : cutoffIndex+"\n"+numHits;

			if(cutoffIndex<(keysP.length-1)){
				int[] lens=getGenericArray(keysP.length);
				for(int i=0; i<keysP.length; i++){lens[i]=count(keysP[i]);}
				Arrays.sort(lens);
				int cutoff=lens[cutoffIndex];
				
				cutoff=Tools.max(lengthHistogram[MIN_INDEX_TO_DROP_LONG_HIT_LIST], cutoff);
				
				int removed=0;
				
				for(int i=0; i<keysP.length; i++){
					int key=keysP[i];
					if(count(key)>cutoff){
						keysP[i]=-1;
						removed++;
						numHits--;
					}
				}
			}
		}
//		assert(checkOffsets(offsets)) : Arrays.toString(offsets);
		
		final ArrayList<SiteScore> result=new ArrayList<SiteScore>(8);
		if(numHits<MIN_APPROX_HITS_TO_KEEP){return result;}
		//assert(checkOffsets(offsetsP)) : Arrays.toString(offsetsP);
		if(numHits<keysP.length){
			int[][] r=shrink2(offsetsP, keysP, keyScoresP);
			assert(r!=null);
			if(r!=null){
				offsetsP=r[0];
				keysP=r[1];
				keyScoresP=r[2];
			}
		}else{
			assert(shrink2(offsetsP, keysP, keyScoresP)==null);
		}
		assert(keysP.length==numHits);
		//assert(checkOffsets(offsetsP)) : Arrays.toString(offsetsP);
		//Reverse the offsets for minus-strand mapping, since they are generated based on quality
		int[] offsetsM=KeyRing.reverseOffsets(offsetsP, KEYLEN, basesP.length);
		if(verbose){
			System.err.println("Reversed offsets: \n"+Arrays.toString(offsetsP)+"  ->\n"+Arrays.toString(offsetsM));
		}
		final int[] keysM=KeyRing.reverseComplementKeys(keysP, KEYLEN);
		
//		assert(checkOffsets(offsetsP)) : Arrays.toString(offsetsP);
//		assert(checkOffsets(offsetsM)) : Arrays.toString(offsetsM);

		assert(!USE_EXTENDED_SCORE || (baseScoresP!=null && (qual==null || baseScoresP.length==qual.length)));
		assert(keyScoresP!=null);
		assert(keyScoresP.length==offsetsP.length) : keyScoresP.length+", "+offsetsP.length+", "+Arrays.toString(keyScoresP);
		final byte[] baseScoresM=Tools.reverseAndCopy(baseScoresP, getBaseScoreArray(baseScoresP.length, 1));
		final int[] keyScoresM=Tools.reverseAndCopy(keyScoresP, getKeyScoreArray(keyScoresP.length, 1));
		final int maxQuickScore=maxQuickScore(offsetsP, keyScoresP);
		
		assert(offsetsM.length==offsetsP.length);
		assert(maxQuickScore==maxQuickScore(offsetsM, keyScoresM));
		
		/*
		 * bestScores:
		 * 
		 * bestScores[0]	currentTopScore
		 * bestScores[1]	maxHits
		 * bestScores[2]	qcutoff
		 * bestScores[3]	bestqscore
		 * bestScores[4]	maxQuickScore
		 * bestScores[5]	perfectsFound
		 */
		final int[] bestScores=KillSwitch.allocInt1D(6);
		
		//This prevents filtering by qscore when a low-quality read only uses a few keys.
		//In that case, extending is more important.
		final boolean prescan_qscore=(PRESCAN_QSCORE && numHits>=5);
		
		int[][] prescanResults=null;
		int[] precounts=null;
		int[] prescores=null;

		int hitsCutoff=0;
		int qscoreCutoff=(int)(MIN_QSCORE_MULT*maxQuickScore);
		
		boolean allBasesCovered=true;
		{
			if(offsetsP[0]!=0){allBasesCovered=false;}
			else if(offsetsP[offsetsP.length-1]!=(basesP.length-KEYLEN)){allBasesCovered=false;}
			else{
				for(int i=1; i<offsetsP.length; i++){
					if(offsetsP[i]>offsetsP[i-1]+KEYLEN){
						allBasesCovered=false;
						break;
					}
				}
			}
		}
		
		//TODO I don't understand this logic
		final boolean pretendAllBasesAreCovered=//false;
					(allBasesCovered ||
					keysP.length>=keysOriginal.length-4 ||
					(keysP.length>=9 && (offsetsP[offsetsP.length-1]-offsetsP[0]+KEYLEN)>Tools.max(40, (int)(basesP.length*.75f))));
		
//		System.err.println(allBasesCovered+"\t"+Arrays.toString(offsetsP));
//		assert(allBasesCovered);
		
		if(prescan_qscore){
			prescanResults=prescanAllBlocks(bestScores,
					keysP, keyScoresP, offsetsP,
					keysM, keyScoresM, offsetsM,
					pretendAllBasesAreCovered);
			
			if(prescanResults!=null){
				precounts=prescanResults[0];
				prescores=prescanResults[1];
			}
			
			if(bestScores[1]<MIN_APPROX_HITS_TO_KEEP){return result;}
			if(bestScores[3]<maxQuickScore*MIN_QSCORE_MULT2){return result;} //if(bestScores[3]<maxQuickScore(offsetsP, keyScoresP)*.10f){return result;}
			
			if(bestScores[3]>=maxQuickScore && pretendAllBasesAreCovered){
				assert(bestScores[3]==maxQuickScore);
				assert(bestScores[1]==numHits);

				hitsCutoff=calcApproxHitsCutoff(keysP.length, bestScores[1], MIN_APPROX_HITS_TO_KEEP, true);
				qscoreCutoff=Tools.max(qscoreCutoff, (int)(bestScores[3]*DYNAMIC_QSCORE_THRESH_PERFECT));
			}else{
				hitsCutoff=calcApproxHitsCutoff(keysP.length, bestScores[1], MIN_APPROX_HITS_TO_KEEP, false);
				qscoreCutoff=Tools.max(qscoreCutoff, (int)(bestScores[3]*PRESCAN_QSCORE_THRESH));
			}
		}
		
		final int maxScore=maxScore(offsetsP, baseScoresP, keyScoresP, basesP.length, true);
		final boolean fullyDefined=AminoAcid.isFullyDefined(basesP);
		assert(bestScores[2]<=0) : Arrays.toString(bestScores);
		
		int cycle=0;
		for(int chrom=minChrom; chrom<=maxChrom; chrom=((chrom&CHROM_MASK_HIGH)+CHROMS_PER_BLOCK)){
			if(precounts==null || precounts[cycle]>=hitsCutoff || prescores[cycle]>=qscoreCutoff){
				find(keysP, basesP, baseScoresP, keyScoresP, chrom, Shared.PLUS,
						offsetsP, obeyLimits, result, bestScores, allBasesCovered, maxScore, fullyDefined);
			}
			if(QUIT_AFTER_TWO_PERFECTS){
				if(bestScores[5]>=2){break;} //if(bestScores[5]>=3 || (bestScores[5]>=2 && chrom<24)){break;} //for human
			}
			cycle++;
			if(precounts==null || precounts[cycle]>=hitsCutoff || prescores[cycle]>=qscoreCutoff){
				find(keysM, basesM, baseScoresM, keyScoresM, chrom, Shared.MINUS,
						offsetsM, obeyLimits, result, bestScores, allBasesCovered, maxScore, fullyDefined);
			}
			if(QUIT_AFTER_TWO_PERFECTS){
				if(bestScores[5]>=2){break;} //if(bestScores[5]>=3 || (bestScores[5]>=2 && chrom<24)){break;} //for human
			}
			cycle++;
		}

//		assert(Read.CHECKSITES(result, basesP, basesM, id, false)); //TODO: Comment out once checked
		
		return result;
	}
	
	/** Search blocks rapidly to find max hits, and perfect sites.  May indicate some blocks can be skipped. */
	private final int[][] prescanAllBlocks(int[] bestScores,
			int[] keysP, int[] keyScoresP, int[] offsetsP,
			int[] keysM, int[] keyScoresM, int[] offsetsM,
			final boolean allBasesCovered){
		
		int[][][] pm=new int[][][] {{keysP, keyScoresP, offsetsP}, {keysM, keyScoresM, offsetsM}};
		
		int bestqscore=0;
		int maxHits=0;
		int minHitsToScore=MIN_APPROX_HITS_TO_KEEP;
		
		final int maxQuickScore=maxQuickScore(offsetsP, keyScoresP);
		
		final int[] counts=precountArray;
		final int[] scores=prescoreArray;
		final int[][] ret=prescanReturn;
		Arrays.fill(counts, keysP.length);
		Arrays.fill(scores, maxQuickScore);
		ret[0]=counts;
		ret[1]=scores;
		
		int cycle=0;
		for(int chrom=minChrom; chrom<=maxChrom; chrom=((chrom&CHROM_MASK_HIGH)+CHROMS_PER_BLOCK)){
			final int baseChrom=baseChrom(chrom);
			for(int pmi=0; pmi<2; pmi++, cycle++){
				
				int[] keys=pm[pmi][0];
				int[] keyScores=pm[pmi][1];
				int[] offsets=pm[pmi][2];
//				int[][] hits=getHitArray(offsets.length);
				
				int[] starts=startArray;
				int[] stops=stopArray;
				int numHits=getHits(keys, chrom, Integer.MAX_VALUE, starts, stops);
				
				if(numHits<minHitsToScore){
					scores[cycle]=-9999;
					counts[cycle]=0;
				}else{
					
//					final int maxQuickScore=maxQuickScore(offsets, keyScores);
					//				System.err.println("maxScore = "+maxScore);
					
					if(numHits<keys.length){
						int[][] r=shrink(starts, stops, offsets, keyScores, offsets.length);
						if(r!=null){
							starts=r[0];
							stops=r[1];
							offsets=r[2];
							keyScores=r[4];
						}
					}
					
					assert(numHits==offsets.length);
					assert(numHits==keyScores.length);
					heap.clear();
					final Quad[] triples=tripleStorage;
					final int[] values=valueArray;
					
					int[] temp=findMaxQscore2(starts, stops, offsets, keyScores, baseChrom, triples, values, minHitsToScore, true,
							bestqscore>=maxQuickScore && allBasesCovered);

					scores[cycle]=temp[0];
					counts[cycle]=temp[1];
					
					bestqscore=Tools.max(temp[0], bestqscore);
					maxHits=Tools.max(maxHits, temp[1]);
					if(bestqscore>=maxQuickScore && allBasesCovered){
						assert(bestqscore==maxQuickScore);
						assert(maxHits==keysP.length) :
							"\nTemp: \t"+Arrays.toString(temp)+", cycle="+cycle+"\n" +
							"Scores: \t"+Arrays.toString(scores)+
							"Counts: \t"+Arrays.toString(counts)+
							"bestqscore: \t"+bestqscore+
							"maxHits: \t"+maxHits+
							"maxQuickScore: \t"+maxQuickScore+
							"numHits: \t"+numHits+
							"minHitsToScore: \t"+minHitsToScore+
							"keys.length: \t"+keys.length;
						
						minHitsToScore=Tools.max(minHitsToScore, maxHits);
						
						{
							//This early exit is optional.  Does not seem to impact speed much either way.
							bestScores[1]=Tools.max(bestScores[1], maxHits);
							bestScores[3]=Tools.max(bestScores[3], bestqscore);
							return ret;
						}
					}
				}
			}
		}
		
		bestScores[1]=Tools.max(bestScores[1], maxHits);
		bestScores[3]=Tools.max(bestScores[3], bestqscore);
		
		if(!RETAIN_BEST_QCUTOFF){bestScores[2]=-9999;}
		
		return ret;
	}
	
	
	/** Search a single block and strand */
	public final ArrayList<SiteScore> find(int[] keys, final byte[] bases, final byte[] baseScores, int[] keyScores,
			final int chrom, final byte strand,
			int[] offsets, final boolean obeyLimits, ArrayList<SiteScore> ssl, int[] bestScores,
			final boolean allBasesCovered, final int maxScore, final boolean fullyDefined){
		
		assert(checkOffsets(offsets)) : Arrays.toString(offsets);
		
		//Index of first location of each key
		int[] starts=startArray;
		//Index of first location of next key (i.e., (last location of key)+1)
		int[] stops=stopArray;
		
		int numHits=getHits(keys, chrom, Integer.MAX_VALUE, starts, stops);
		if(numHits<MIN_APPROX_HITS_TO_KEEP){return ssl;}
//		assert(false) : "\n"+Data.getChromosome(1).minIndex+"\n"+Data.getChromosome(1).maxIndex+"\n"+Data.getChromosome(1).array.length;
//		assert(false) : "\n"+Data.getChromosome(1).minIndex+"\n"+Data.getChromosome(1).maxIndex+"\n"+Data.getChromosome(1).array.length+"\n"+(char)Data.getChromosome(1).array[0];
//		assert(false) : numHits+"\n"+Arrays.toString(starts)+"\n"+Arrays.toString(stops)+"\n"+Arrays.toString(index[0].getHitList(starts[0], stops[0]))+"\n";
		if(USE_CAMELWALK){
			if(USE_SLOWALK3){
				if(!RETAIN_BEST_SCORES){Arrays.fill(bestScores, 0);}
				ssl=camelWalk3(starts, stops, bases, baseScores, keyScores, offsets, chrom, strand, obeyLimits, ssl, bestScores, allBasesCovered, maxScore, fullyDefined);
			}else{
				assert(false) : "TODO";
				ssl=slowWalk2(starts, stops, bases, baseScores, keyScores, offsets, chrom, strand, obeyLimits, ssl, fullyDefined);
			}
		}else{
			if(USE_SLOWALK3){
				if(!RETAIN_BEST_SCORES){Arrays.fill(bestScores, 0);}
				ssl=slowWalk3(starts, stops, bases, baseScores, keyScores, offsets, chrom, strand, obeyLimits, ssl, bestScores, allBasesCovered, maxScore, fullyDefined);
			}else{
				ssl=slowWalk2(starts, stops, bases, baseScores, keyScores, offsets, chrom, strand, obeyLimits, ssl, fullyDefined);
			}
		}
		
		return ssl;
	}
	
	/** Compress arrays by removing null/empty lists */
	private final int[][] shrink(int[] starts, int[] stops, int[] offsets, int[] keyScores, final int len){
		int numHits=0;
		for(int i=0; i<len; i++){
			if(starts[i]>=0){numHits++;}
		}

		if(numHits==offsets.length){
			return null;
		}
		int[][] r=shrinkReturn3;
		int[] starts2=startArray;
		int[] stops2=stopArray;
		int[] offsets2=getOffsetArray(numHits);
		int[] keyScores2=new int[numHits];

		for(int i=0, j=0; i<len; i++){
			if(starts[i]>=0){
				starts2[j]=starts[i];
				stops2[j]=stops[i];
				offsets2[j]=offsets[i];
				keyScores2[j]=keyScores[i];
				j++;
			}
		}
		r[0]=starts2;
		r[1]=stops2;
		r[2]=offsets2;
		r[4]=keyScores2;
		return r;
	}
	
	/** Removes "-1" keys. */
	private final int[][] shrink2(int[] offsets, int[] keys, int[] keyScores){


		int numHits=0;
		for(int i=0; i<keys.length; i++){
			if(keys[i]>=0){numHits++;}
		}


		assert(checkOffsets(offsets)) : Arrays.toString(offsets);
		if(numHits==keys.length){
			return null;
		}else{
			int[][] r=shrinkReturn2;
			int[] offsets2=getOffsetArray(numHits);
			assert(offsets2!=offsets);
			assert(offsets2.length<offsets.length);
			int[] keys2=new int[numHits];
			int[] keyScores2=new int[numHits];

			for(int i=0, j=0; i<keys.length; i++){
				if(keys[i]>=0){
					offsets2[j]=offsets[i];
					keys2[j]=keys[i];
					keyScores2[j]=keyScores[i];
					j++;
				}
			}
			assert(checkOffsets(offsets2)) : "\nnumHits="+numHits+"\n"+Arrays.toString(offsets)+" -> \n"+Arrays.toString(offsets2)+"\n"+
				"\n"+Arrays.toString(keys)+" -> \n"+Arrays.toString(keys2)+"\n";
			r[0]=offsets2;
			r[1]=keys2;
			r[2]=keyScores2;
			return r;
		}
	}
	
	
	/** This uses a heap to track next column to increment */
	private final ArrayList<SiteScore> slowWalk2(int[] starts, int[] stops, final byte[] bases,
			final byte[] baseScores, int[] keyScores, int[] offsets,
			final int baseChrom_, final byte strand, final boolean obeyLimits, ArrayList<SiteScore> ssl, final boolean fullyDefined){
		
		if(SHRINK_BEFORE_WALK){
			int[][] r=shrink(starts, stops, offsets, keyScores, offsets.length);
			if(r!=null){
				starts=r[0];
				stops=r[1];
				offsets=r[2];
				keyScores=r[4];
			}
		}
		
		final int numHits=offsets.length; //After shrink
		
		assert(numHits==offsets.length);
		assert(numHits==keyScores.length);
		
		assert(!(!SHRINK_BEFORE_WALK && ADD_SCORE_Z));
		
		//This can be done before or after shrinking, but the results will change depending on MIN_SCORE_MULT and etc.
		final int maxScore=maxScore(offsets, baseScores, keyScores, bases.length, true);
//		final int maxQuickScore=(!USE_EXTENDED_SCORE ? maxScore : maxQuickScore(offsets));
//		System.err.println("maxScore = "+maxScore);
		
//		final int minScore=(obeyLimits ? (int)(MIN_SCORE_MULT*maxScore) : (int)(MIN_SCORE_MULT*0.85f*maxScore));
		final int minScore=(obeyLimits ? (int)(MIN_SCORE_MULT*maxScore) : (int)(MIN_SCORE_MULT*1.25f*maxScore));
//		final int minQuickScore=(!USE_EXTENDED_SCORE ? minScore : (int)(maxQuickScore*0.15f));
//		final int minScore=(int)(MIN_SCORE_MULT*maxScore);
//		System.err.println("minScore = "+minScore);
		
		final int baseChrom=baseChrom(baseChrom_);
		
		
		heap.clear();
		final Quad[] triples=tripleStorage;
		
		final Block b=index[baseChrom];
		final int[] values=valueArray;
		final int[] sizes=sizeArray;
		final int[] locArray=(USE_EXTENDED_SCORE ? getLocArray(bases.length) : null);
		
		for(int i=0; i<numHits; i++){
			final int[] sites=b.sites;
			final int start=starts[i];
			sizes[i]=b.length(start, stops[i]);
			assert(sizes[i]>0);
			
			int a=sites[start];
			int a2;
			if((a&SITE_MASK)>=offsets[i]){
				a2=a-offsets[i];
			}else{
				int ch=numberToChrom(a, baseChrom);
				int st=numberToSite(a);
				int st2=Tools.max(st-offsets[i], 0);
				a2=toNumber(st2, ch);
			}
			assert(numberToChrom(a, baseChrom) == numberToChrom(a2, baseChrom));

			Quad t=triples[i];
			assert(t!=null) : "Should be using tripleStorage";
			assert(i==t.column);
			t.row=start;
			t.site=a2;
			t.list=sites;
			values[i]=a2;

			heap.add(t);
		}
		
		if(ssl==null){ssl=new ArrayList<SiteScore>(8);}
		
		int currentTopScore=-999999999;
		
		int cutoff=minScore;
		
		int maxHits=0;
		int approxHitsCutoff=MIN_APPROX_HITS_TO_KEEP;
		
//		System.out.println("\nEntering SS loop:");
//		System.out.println("maxScore="+maxScore+"\tminScore="+minScore+"\tcurrentTopScore="+currentTopScore+"\n" +
//				"cutoff="+cutoff+"\tmaxHits="+maxHits+"\tapproxHitsCutoff="+approxHitsCutoff);
//		System.out.println();
		
		SiteScore prevSS=null;
		while(!heap.isEmpty()){
			Quad t=heap.peek();
			final int site=t.site;
			final int centerIndex=t.column;
			
			int maxNearbySite=site;
			int approxHits=0;
			
			{//Inner loop
				final int minsite=site-MAX_INDEL, maxsite=site+MAX_INDEL2;
				for(int column=0, chances=numHits-approxHitsCutoff; column<numHits && chances>=0; column++){
					final int x=values[column];
					assert(x==triples[column].site);
					if(x>=minsite && x<=maxsite){
						maxNearbySite=(x>maxNearbySite ? x : maxNearbySite);
						approxHits++;
					}else{chances--;}
				}
			}
			hist_hits[Tools.min(HIT_HIST_LEN, approxHits)]++;

			assert(centerIndex>=0) : centerIndex;
			
			//I don't remember what this assertion was for or why, but it's causing trouble.
			//assert(approxHits>=1 || approxHitsCutoff>1) : approxHits+", "+approxHitsCutoff+", "+numHits+", "+t.column;
			if(approxHits>=approxHitsCutoff){
				
				int score;
				
				int mapStart=site, mapStop=maxNearbySite;
				
				if(USE_EXTENDED_SCORE){
					final int chrom=numberToChrom(site, baseChrom);
					score=extendScore(bases, baseScores, offsets, values, chrom, centerIndex, locArray, numHits, approxHits);
					if(true/*USE_AFFINE_SCORE*/){
						//Correct begin and end positions if they changed.
						int min=Integer.MAX_VALUE;
						int max=Integer.MIN_VALUE;
						for(int i=0; i<locArray.length; i++){
							int x=locArray[i];
							if(x>-1){
								if(x<min){min=x;}
								if(x>max){max=x;}
							}
						}
						
//						assert(min>-1 && max>-1) : Arrays.toString(locArray); //TODO: How did this assertion trigger?
						if(min<0 || max<0){
							//Note: This error can trigger if minChrom and maxChrom do not align to block boundaries
							if(!Shared.anomaly){
								Shared.anomaly=true;
								System.err.println("Anomaly in "+getClass().getName()+".slowWalk: "+
										chrom+", "+mapStart+", "+mapStop+", "+centerIndex+", "+
										Arrays.toString(locArray)+"\n"+
										Arrays.toString(values)+"\n"+
										new String(bases)+"\nstrand="+strand+"\n");
								System.err.println();
							}
							score=-99999;
						}
						
						
						mapStart=toNumber(min, chrom);
						mapStop=toNumber(max, chrom);
						
						
						
//						System.err.println("site="+site+", maxNearbySite="+maxNearbySite+", min="+min+", max="+max+
//								", mapStart="+mapStart+", mapStop="+mapStop);
						
//						if(chrom==17 && absdif(min, 30354420)<2000){
//							System.err.println("\n*****\n");
//							System.err.println("site="+site+" ("+numberToSite(site)+"), maxNearbySite="+maxNearbySite+
//									" ("+numberToSite(maxNearbySite)+"), min="+min+", max="+max+
//									", mapStart="+mapStart+", mapStop="+mapStop);
//							System.err.println();
//							System.err.println(Arrays.toString(locArray));
//							System.err.println();
//							System.err.println("chrom="+chrom);
//							System.err.println("score="+score);
//						}
					}
				}else{
					score=quickScore(values, keyScores, centerIndex, offsets, sizes, true, approxHits, numHits);
					if(ADD_SCORE_Z){
						int scoreZ=scoreZ2(values, centerIndex, offsets, approxHits, numHits);
						score+=scoreZ;
					}
				}
				
				if(score>=cutoff){

					if(score>currentTopScore){
//						System.err.println("New top score!");

						if(DYNAMICALLY_TRIM_LOW_SCORES){

							maxHits=Tools.max(approxHits, maxHits);
//							approxHitsCutoff=Tools.max(approxHitsCutoff, maxHits);
							approxHitsCutoff=Tools.max(approxHitsCutoff, maxHits-1); //More sensitive, but slower

							cutoff=Tools.max(cutoff, (int)(score*DYNAMIC_SCORE_THRESH));

							//									cutoff=Tools.max(cutoff, minScore+(int)((score-minScore)*DYNAMIC_SCORE_THRESH));
							if(USE_EXTENDED_SCORE && score>=maxScore){
								cutoff=Tools.max(cutoff, (int)(score*0.95f));
							}
						}

						currentTopScore=score;
						
//						System.out.println("New top score: "+currentTopScore+" \t("+cutoff+")");
					}

					final int chrom=numberToChrom(mapStart, baseChrom);
					final int site2=numberToSite(mapStart);
					final int site3=numberToSite(mapStop)+bases.length-1;
					
					assert(NUM_CHROM_BITS==0 || site2<SITE_MASK-1000) : "chrom="+chrom+", strand="+strand+", site="+site+
						", maxNearbySite="+maxNearbySite+", site2="+site2+", site3="+site3+", read.length="+bases.length+
						"\n\n"+Arrays.toString(b.getHitList(centerIndex));
					assert(site2<site3) : "chrom="+chrom+", strand="+strand+", site="+site+
						", maxNearbySite="+maxNearbySite+", site2="+site2+", site3="+site3+", read.length="+bases.length;
					
					
					//Note:  I could also do this as soon as score is calculated.
//					if(ADD_SCORE_Z){
//						int scoreZ=scoreZ2(values, centerIndex, offsets);
//						score+=scoreZ;
//					}
					
					//This block is optional, but tries to eliminate multiple identical alignments
//					SiteScore prevSS=(ssl.size()<1 ? null : ssl.get(ssl.size()-1));
					
					SiteScore ss=null;
					final boolean perfect1=USE_EXTENDED_SCORE && score==maxScore && fullyDefined;
					
					int[] gapArray=null;
					if(site3-site2>=MINGAP+bases.length){
						gapArray=makeGapArray(locArray, mapStart, MINGAP);
						if(gapArray!=null){
							int sub=site2-mapStart;//thus site2=mapStart+sub
							for(int i=0; i<gapArray.length; i++){
								gapArray[i]+=sub;
							}
							assert(gapArray[0]==mapStart) : gapArray[0]+", "+mapStart;
							assert(gapArray[gapArray.length-1]==mapStop);
						}
						assert(false) : Arrays.toString(locArray);
					}
					
					if(gapArray==null && prevSS!=null && prevSS.gaps==null &&
							prevSS.chrom==chrom && prevSS.strand==strand && overlap(prevSS.start, prevSS.stop, site2, site3)){

						int betterScore=Tools.max(score, prevSS.score);
						int minStart=Tools.min(prevSS.start, site2);
						int maxStop=Tools.max(prevSS.stop, site3);
						final boolean perfect2=USE_EXTENDED_SCORE && prevSS.score==maxScore && fullyDefined;
						assert(!USE_EXTENDED_SCORE || perfect2==prevSS.perfect);

						boolean shortEnough=(!LIMIT_SUBSUMPTION_LENGTH_TO_2X  || (maxStop-minStart<2*bases.length));

						if(prevSS.start==site2 && prevSS.stop==site3){
							prevSS.score=prevSS.quickScore=betterScore;
						}else if(SUBSUME_SAME_START_SITES && shortEnough && prevSS.start==site2){
							if(perfect2){
								//do nothing
							}else if(perfect1){
								prevSS.setStop(site3);
								prevSS.setPerfect();
							}else{
								prevSS.setStop(maxStop);
							}
							prevSS.score=prevSS.quickScore=betterScore;
						}else if(SUBSUME_SAME_STOP_SITES && shortEnough && prevSS.stop==site3){
							if(perfect2){
								//do nothing
							}else if(perfect1){
								prevSS.setStart(site2);
								prevSS.setPerfect();
							}else{
								prevSS.setStart(minStart);
							}
							prevSS.score=prevSS.quickScore=betterScore;
						}else if(SUBSUME_OVERLAPPING_SITES && shortEnough && (maxStop-minStart<=bases.length+MAX_SUBSUMPTION_LENGTH)
								&& !perfect1 && !perfect2){
							prevSS.setLimits(minStart, maxStop);
							prevSS.score=prevSS.quickScore=betterScore;
						}else{
							ss=new SiteScore(chrom, strand, site2, site3, approxHits, score, false, perfect1);
							if(!perfect1){ss.setPerfect(bases);}
							assert(!perfect1 || ss.stop-ss.start==bases.length-1);
						}
						assert(!perfect2 || prevSS.stop-prevSS.start==bases.length-1);
					}else{
						ss=new SiteScore(chrom, strand, site2, site3, approxHits, score, false, perfect1);
						if(!perfect1){ss.setPerfect(bases);}
						ss.gaps=gapArray;
						if(verbose && gapArray!=null){
							System.err.println(ss.toText()+"\t"+Arrays.toString(gapArray)+"\n"+Arrays.toString(locArray)+"\n");
						}
					}
					
					if(ss!=null){
//						System.out.println("Added site "+ss.toText());
						ssl.add(ss);
						prevSS=ss;
					}else{
//						System.out.println("Subsumed site "+new SiteScore(chrom, strand, site2, site3, score).toText());
					}
					
//					if(prevSS!=null && prevSS.chrom==chrom && prevSS.strand==strand && overlap(prevSS.start, prevSS.stop, site2, site3)){
//						int betterScore=Tools.max(score, prevSS.score);
//						if(prevSS.start==site2 && prevSS.stop==site3){
//							prevSS.score=prevSS.quickScore=betterScore;
//						}else if(prevSS.start==site2
//								/*isWithin(prevSS.start, prevSS.stop, site2, site3) ||
//								isWithin(site2, site3, prevSS.start, prevSS.stop)*/){
//							prevSS.score=prevSS.quickScore=betterScore;
//							assert(prevSS.start<prevSS.stop);
////							prevSS.start=Tools.min(prevSS.start, site2);
//							prevSS.stop=Tools.max(prevSS.stop, site3);
//							assert(prevSS.start<prevSS.stop);
//						}else{
//							SiteScore ss=new SiteScore(chrom, strand, site2, site3, score);
//							ssl.add(ss);
//							prevSS=ss;
//						}
//					}else{
//						SiteScore ss=new SiteScore(chrom, strand, site2, site3, score);
//						ssl.add(ss);
//						prevSS=ss;
//					}

				}
			}

			while(heap.peek().site==site){ //Remove all identical elements, and add subsequent elements
				final Quad t2=heap.poll();
				final int row=t2.row+1, col=t2.column;
				if(row<stops[col]){
					t2.row=row;
					
					int a=t2.list[row];
					int a2;
					if((a&SITE_MASK)>=offsets[col]){
						a2=a-offsets[col];
						
						assert(numberToChrom(a, baseChrom) == numberToChrom(a2, baseChrom)) :
							"baseChrom="+baseChrom+", chrom="+numberToChrom(a, baseChrom)+", strand="+strand+", site="+site+
							", maxNearbySite="+maxNearbySite+", a="+a+", a2="+a2+", offsets["+col+"]="+offsets[col];
					}else{
						int ch=numberToChrom(a, baseChrom);
						int st=numberToSite(a);
						int st2=Tools.max(st-offsets[col], 0);
						a2=toNumber(st2, ch);
						
						assert(numberToChrom(a, baseChrom) == numberToChrom(a2, baseChrom)) :
							"baseChrom="+baseChrom+", chrom="+numberToChrom(a, baseChrom)+", strand="+strand+", site="+site+
							", maxNearbySite="+maxNearbySite+", a="+a+", a2="+a2+", offsets["+col+"]="+offsets[col];
					}
					
					assert(numberToChrom(a, baseChrom) == numberToChrom(a2, baseChrom)) :
						"baseChrom="+baseChrom+", chrom="+numberToChrom(a, baseChrom)+", strand="+strand+", site="+site+
						", maxNearbySite="+maxNearbySite+", a="+a+", a2="+a2+", offsets["+col+"]="+offsets[col];

					t2.site=a2;
					values[col]=a2;
					heap.add(t2);
				}
				if(heap.isEmpty()){break;}
			}

		}
		
		return ssl;
	}
	
	/** This uses a heap to track next column to increment */
	private final ArrayList<SiteScore> slowWalk3(int[] starts, int[] stops, final byte[] bases,
			final byte[] baseScores, int[] keyScores, int[] offsets,
			final int baseChrom_, final byte strand, final boolean obeyLimits, ArrayList<SiteScore> ssl,
			int[] bestScores, final boolean allBasesCovered, final int maxScore, final boolean fullyDefined){
		assert(USE_EXTENDED_SCORE);
		
		final int numKeys=offsets.length; //Before shrink
		
		//This can be done before or after shrinking, but the results will change depending on MIN_SCORE_MULT and etc.
		final int maxQuickScore=maxQuickScore(offsets, keyScores);
		
		if(SHRINK_BEFORE_WALK){
			int[][] r=shrink(starts, stops, offsets, keyScores, offsets.length);
			if(r!=null){
				starts=r[0];
				stops=r[1];
				offsets=r[2];
				keyScores=r[4];
			}
		}
		
		final int numHits=offsets.length; //After shrink
		
		
		assert(numHits==offsets.length);
		assert(numHits==keyScores.length);
		
		usedKeys+=numHits;
		usedKeyIterations++;
		
		final boolean filter_by_qscore=(FILTER_BY_QSCORE && numKeys>=5);
		
		assert(!(!SHRINK_BEFORE_WALK && ADD_SCORE_Z));
		
		
//		final int minScore=(obeyLimits ? (int)(MIN_SCORE_MULT*maxScore) : (int)(MIN_SCORE_MULT*0.85f*maxScore));
		final int minScore=(obeyLimits ? (int)(MIN_SCORE_MULT*maxScore) : (int)(MIN_SCORE_MULT*1.25f*maxScore));
		final int minQuickScore=(int)(MIN_QSCORE_MULT*maxQuickScore);
		
		final int baseChrom=baseChrom(baseChrom_);
		
		heap.clear();
		
		final Quad[] triples=tripleStorage;

		final int[] values=valueArray;
		final int[] sizes=sizeArray;
		final int[] locArray=(USE_EXTENDED_SCORE ? getLocArray(bases.length) : null);
		final Block b=index[baseChrom];
		
		if(ssl==null){ssl=new ArrayList<SiteScore>(8);}
		
		int currentTopScore=bestScores[0];
		int cutoff=Tools.max(minScore, (int)(currentTopScore*DYNAMIC_SCORE_THRESH));
		
		int qcutoff=Tools.max(bestScores[2], minQuickScore);
		int bestqscore=bestScores[3];
		int maxHits=bestScores[1];
		int perfectsFound=bestScores[5];
		assert((currentTopScore>=maxScore) == (perfectsFound>0)) : currentTopScore+", "+maxScore+", "+perfectsFound+", "+maxHits+", "+numHits+", "+new String(bases);
		int approxHitsCutoff=calcApproxHitsCutoff(numKeys, maxHits, MIN_APPROX_HITS_TO_KEEP, currentTopScore>=maxScore);
		if(approxHitsCutoff>numHits){return ssl;}
		
		final boolean shortCircuit=(allBasesCovered && numKeys==numHits && filter_by_qscore);
		
		
		assert(USE_EXTENDED_SCORE);
		
		if(currentTopScore>=maxScore){
			assert(currentTopScore==maxScore);
			qcutoff=Tools.max(qcutoff, (int)(maxQuickScore*DYNAMIC_QSCORE_THRESH_PERFECT));
		}
			
		
		for(int i=0; i<numHits; i++){
			final int[] sites=b.sites;
			final int start=starts[i];
			sizes[i]=b.length(start, stops[i]);
			assert(sizes[i]>0);
			
			int a=sites[start];
			int a2;
			if((a&SITE_MASK)>=offsets[i]){
				a2=a-offsets[i];
			}else{
				int ch=numberToChrom(a, baseChrom);
				int st=numberToSite(a);
				int st2=Tools.max(st-offsets[i], 0);
				a2=toNumber(st2, ch);
			}
			assert(numberToChrom(a, baseChrom) == numberToChrom(a2, baseChrom));

			Quad t=triples[i];
			assert(t!=null) : "Should be using tripleStorage";
			assert(i==t.column);
			t.row=start;
			t.site=a2;
			t.list=sites;
			values[i]=a2;

			heap.add(t);
		}

//		System.out.println("\nEntering SS loop:");
//		System.out.println("maxScore="+maxScore+"\tminScore="+minScore+"\tcurrentTopScore="+currentTopScore+"\n" +
//				"cutoff="+cutoff+"\tmaxHits="+maxHits+"\tapproxHitsCutoff="+approxHitsCutoff);
//		System.out.println("maxQuickScore="+maxQuickScore+"\tminQuickScore="+minQuickScore+"\tqcutoff="+qcutoff);
		
		
		SiteScore prevSS=null;
		while(!heap.isEmpty()){
			Quad t=heap.peek();
			final int site=t.site;
			final int centerIndex=t.column;
			
			int maxNearbySite=site;


			int approxHits=0;
			
			{//Inner loop
				final int minsite=site-MAX_INDEL, maxsite=site+MAX_INDEL2;
				for(int column=0, chances=numHits-approxHitsCutoff; column<numHits && chances>=0; column++){
					final int x=values[column];
					assert(x==triples[column].site);
					if(x>=minsite && x<=maxsite){
						maxNearbySite=(x>maxNearbySite ? x : maxNearbySite);
						approxHits++;
					}else{chances--;}
				}
			}
			hist_hits[Tools.min(HIT_HIST_LEN, approxHits)]++;

			assert(centerIndex>=0) : centerIndex;
			
			//I don't remember what this assertion was for or why, but it's causing trouble.
			//assert(approxHits>=1 || approxHitsCutoff>1) : approxHits+", "+approxHitsCutoff+", "+numHits+", "+t.column;
			if(approxHits>=approxHitsCutoff){
				
				int score;
				int qscore=(filter_by_qscore ? quickScore(values, keyScores, centerIndex, offsets, sizes, true, approxHits, numHits) : qcutoff);
				if(ADD_SCORE_Z){
					int scoreZ=scoreZ2(values, centerIndex, offsets, approxHits, numHits);
					qscore+=scoreZ;
				}
				
				int mapStart=site, mapStop=maxNearbySite;
				
				assert(USE_EXTENDED_SCORE);
				
				boolean locArrayValid=false;
				if(qscore<qcutoff){
					score=-1;
				}else{

					final int chrom=numberToChrom(site, baseChrom);
					
					//TODO Note that disabling the shortCircuit code seems to make things run 2% faster (with identical results).
					//However, theoretically, shortCircuit should be far more efficient.  Test both ways on cluster and on a larger run.
					//May have something to do with compiler loop optimizations.
					if(shortCircuit && qscore==maxQuickScore){
						assert(approxHits==numKeys);
						score=maxScore;
					}else{
						if(verbose){
							System.err.println("numHits="+numHits+", approxHits="+approxHits+", keys="+numKeys+", centerIndex="+centerIndex+
									", qscore="+qscore+", qcutoff="+qcutoff+", filter_by_qscore="+filter_by_qscore);
							System.err.println("Extending "+Arrays.toString(values));
						}
						score=extendScore(bases, baseScores, offsets, values, chrom, centerIndex, locArray, numHits, approxHits);
						locArrayValid=true;

						if(verbose){
							System.err.println("score: "+score);
							System.err.println("locArray: "+Arrays.toString(locArray));
						}
					
						//Correct begin and end positions if they changed.
						int min=Integer.MAX_VALUE;
						int max=Integer.MIN_VALUE;
						for(int i=0; i<locArray.length; i++){
							int x=locArray[i];
							if(x>-1){
								if(x<min){min=x;}
								if(x>max){max=x;}
							}
						}

						if(score>=maxScore){
							assert(min==max && min>-1) : "\n"+score+", "+maxScore+", "+min+", "+max+
							", "+(max-min)+", "+bases.length+"\n"+Arrays.toString(locArray)+"\n";
						}

						//							assert(min>-1 && max>-1) : Arrays.toString(locArray); //TODO: How did this assertion trigger?
						if(min<0 || max<0){
							if(!Shared.anomaly){
								Shared.anomaly=true;
								System.err.println("\nAnomaly in "+getClass().getName()+".slowWalk:\n"+
										"chrom="+chrom+", mapStart="+mapStart+", mapStop="+mapStop+", centerIndex="+centerIndex+", strand="+strand+"\n"+
										"score="+score+", maxScore="+maxScore+", qscore="+qscore+", filter_by_qscore="+filter_by_qscore+"\n"+
										"numHits="+approxHits+", approxHits="+approxHits+"\n"+
										"min="+min+", max="+max+", (max-min)="+(max-min)+"\n"+
										"bases.length="+bases.length+"\n"+Arrays.toString(locArray)+"\n"+
										"locArray:\t"+Arrays.toString(locArray)+"\n"+
										"values:\t"+Arrays.toString(values)+"\n"+
										"bases:\t"+new String(bases));
								System.err.println();
								assert(false);
							}
							score=-99999;
						}

						//mapStart and mapStop are indices
						mapStart=toNumber(min, chrom);
						mapStop=toNumber(max, chrom);

						if(score>=maxScore){
							assert(mapStop-mapStart==0) : "\n"+score+", "+maxScore+", "+min+", "+max+
							", "+(max-min)+", "+(mapStop-mapStart)+", "+bases.length+"\n"+Arrays.toString(locArray)+"\n";
						}
					}

					if(score==maxScore){
						qcutoff=Tools.max(qcutoff, (int)(maxQuickScore*DYNAMIC_QSCORE_THRESH_PERFECT));
						approxHitsCutoff=calcApproxHitsCutoff(numKeys, maxHits, MIN_APPROX_HITS_TO_KEEP, true);
					}
					
					if(score>=cutoff){
						qcutoff=Tools.max(qcutoff, (int)(qscore*DYNAMIC_QSCORE_THRESH));
						bestqscore=Tools.max(qscore, bestqscore);
					}
				}
				
				if(score>=cutoff){

					if(score>currentTopScore){
//						System.err.println("New top score!");

						if(DYNAMICALLY_TRIM_LOW_SCORES){

							maxHits=Tools.max(approxHits, maxHits);
							approxHitsCutoff=calcApproxHitsCutoff(numKeys, maxHits, approxHitsCutoff, currentTopScore>=maxScore);

							cutoff=Tools.max(cutoff, (int)(score*DYNAMIC_SCORE_THRESH));
							if(score>=maxScore){
								assert(USE_EXTENDED_SCORE);
								cutoff=Tools.max(cutoff, (int)(score*0.95f));
							}
						}

						currentTopScore=score;

//						System.out.println("New top score: "+currentTopScore+" \t("+cutoff+")");
					}

					final int chrom=numberToChrom(mapStart, baseChrom);
					final int site2=numberToSite(mapStart);
					final int site3=numberToSite(mapStop)+bases.length-1;
					
					assert(site2!=site3) : site2+", "+site3+", "+mapStart+", "+mapStop;
					
					assert(NUM_CHROM_BITS==0 || site2<SITE_MASK-1000) : "chrom="+chrom+", strand="+strand+", site="+site+
						", maxNearbySite="+maxNearbySite+", site2="+site2+", site3="+site3+", read.length="+bases.length+
						"\n\n"+Arrays.toString(b.getHitList(centerIndex));
					assert(site2<site3) : "chrom="+chrom+", strand="+strand+", site="+site+
						", maxNearbySite="+maxNearbySite+", site2="+site2+", site3="+site3+", read.length="+bases.length;
					
					
					int[] gapArray=null;
					if(site3-site2>=MINGAP+bases.length){
						assert(locArrayValid) : "Loc array was not filled.";
//						System.err.println("****\n"+Arrays.toString(locArray)+"\n");
//						int[] clone=locArray.clone();
						gapArray=makeGapArray(locArray, site2, MINGAP);
						if(gapArray!=null){
//							System.err.println(Arrays.toString(locArray)+"\n");
//							System.err.println(Arrays.toString(gapArray));
//
////							int sub=site2-mapStart;//thus site2=mapStart+sub
////							for(int i=0; i<gapArray.length; i++){
////								gapArray[i]+=sub;
////							}
////							System.err.println(Arrays.toString(gapArray));
//
//							System.err.println(mapStart+" -> "+site2);
//							System.err.println(mapStop+" -> "+site3);
							
							assert(gapArray[0]>=site2 && gapArray[0]-site2<bases.length);
							assert(gapArray[gapArray.length-1]<=site3 && site3-gapArray[gapArray.length-1]<bases.length) : "\n"+
								mapStart+" -> "+site2+"\n"+
								mapStop+" -> "+site3+"\n\n"+
								Arrays.toString(gapArray)+"\n\n"+
//								Arrays.toString(clone)+"\n\n"+
								Arrays.toString(locArray)+"\n"+
								"numHits="+numHits+", "+
								"heap.size="+heap.size()+", "+
								"numHits="+numHits+", "+
								"approxHits="+approxHits+"\n";
							gapArray[0]=Tools.min(gapArray[0], site2);
							gapArray[gapArray.length-1]=Tools.max(gapArray[gapArray.length-1], site3);
						}
						if(verbose){System.err.println("@ site "+site2+", made gap array: "+Arrays.toString(gapArray));}
//						assert(false) : Arrays.toString(locArray);
					}
					
					
					//This block is optional, but tries to eliminate multiple identical alignments
					
					SiteScore ss=null;
					final boolean perfect1=USE_EXTENDED_SCORE && score==maxScore && fullyDefined;
					final boolean inbounds=(site2>=0 && site3<Data.chromLengths[chrom]);
//					if(!inbounds){System.err.println("Index tossed out-of-bounds site chr"+chrom+", "+site2+"-"+site3);}
					
					if(verbose){
						System.err.println("Considering site chr"+chrom+", site2="+site2+", site3="+site3+", mapStart="+mapStart+", mapStop="+mapStop);
						if(!inbounds){System.err.println("Index tossed out-of-bounds site.");}
						assert(site2!=site3) : site2+", "+site3+", "+mapStart+", "+mapStop;
					}
					//assert((ss==null || !ss.semiperfect) && (prevSS==null || !prevSS.semiperfect)) : (ss==null ? false : ss.semiperfect)+", "+(prevSS==null ? false : prevSS.semiperfect); //***
					
					if(inbounds && !SEMIPERFECTMODE && !PERFECTMODE && gapArray==null && prevSS!=null &&
							prevSS.chrom==chrom && prevSS.strand==strand && overlap(prevSS.start, prevSS.stop, site2, site3)){
						
						if(verbose){System.err.println("Considering overlapping site chr"+chrom+", "+site2+"-"+site3);}

						final int betterScore=Tools.max(score, prevSS.score);
						final int minStart=Tools.min(prevSS.start, site2);
						final int maxStop=Tools.max(prevSS.stop, site3);
						final boolean perfect2=USE_EXTENDED_SCORE && prevSS.score==maxScore && fullyDefined;
						assert(!USE_EXTENDED_SCORE || perfect2==prevSS.perfect);

						final boolean shortEnough=(!LIMIT_SUBSUMPTION_LENGTH_TO_2X  || (maxStop-minStart<2*bases.length));

						if(prevSS.start==site2 && prevSS.stop==site3){
							if(verbose){System.err.println("Class 1: Same bounds as last site.");}
							prevSS.score=prevSS.quickScore=betterScore;
							prevSS.perfect=(prevSS.perfect || perfect1 || perfect2);
							if(prevSS.perfect){prevSS.semiperfect=true;}
							//assert((ss==null || !ss.semiperfect) && (prevSS==null || !prevSS.semiperfect)) : (ss==null ? false : ss.semiperfect)+", "+(prevSS==null ? false : prevSS.semiperfect); //***
						}else if(SUBSUME_SAME_START_SITES && shortEnough && prevSS.start==site2 && !prevSS.semiperfect){
							if(verbose){System.err.println("Class 2: Same start as last site.");}
							if(perfect2){
								//do nothing
							}else if(perfect1){
								prevSS.setStop(site3);
								if(!prevSS.perfect){perfectsFound++;}//***$
								prevSS.perfect=prevSS.semiperfect=true;
							}else{
								prevSS.setStop(maxStop);
								prevSS.setPerfect(bases);
							}
							//assert((ss==null || !ss.semiperfect) && (prevSS==null || !prevSS.semiperfect)) : (ss==null ? false : ss.semiperfect)+", "+(prevSS==null ? false : prevSS.semiperfect); //***
							prevSS.score=prevSS.quickScore=betterScore;
						}else if(SUBSUME_SAME_STOP_SITES && shortEnough && prevSS.stop==site3 && !prevSS.semiperfect){
							if(verbose){System.err.println("Class 3: Same stop as last site.");}
							if(perfect2){
								//do nothing
							}else if(perfect1){
								prevSS.setStart(site2);
								if(!prevSS.perfect){perfectsFound++;}//***$
								prevSS.perfect=prevSS.semiperfect=true;
							}else{
								prevSS.setStart(minStart);
								prevSS.setPerfect(bases);
							}
							//assert((ss==null || !ss.semiperfect) && (prevSS==null || !prevSS.semiperfect)) : (ss==null ? false : ss.semiperfect)+", "+(prevSS==null ? false : prevSS.semiperfect); //***
							prevSS.score=prevSS.quickScore=betterScore;
						}else if(SUBSUME_OVERLAPPING_SITES && shortEnough && (maxStop-minStart<=bases.length+MAX_SUBSUMPTION_LENGTH)
								&& !perfect1 && !perfect2 && !prevSS.semiperfect){
							if(verbose){System.err.println("Class 4.");}
							prevSS.setLimits(minStart, maxStop);
							prevSS.score=prevSS.quickScore=betterScore;
							prevSS.setPerfect(bases);
							//assert((ss==null || !ss.semiperfect) && (prevSS==null || !prevSS.semiperfect)) : (ss==null ? false : ss.semiperfect)+", "+(prevSS==null ? false : prevSS.semiperfect); //***
						}else{
							if(verbose){System.err.println("Class 5: Making new site");}
							ss=new SiteScore(chrom, strand, site2, site3, approxHits, score, false, perfect1);
							if(!perfect1){ss.setPerfect(bases);}
							//assert((ss==null || !ss.semiperfect) && (prevSS==null || !prevSS.semiperfect)) : (ss==null ? false : ss.semiperfect)+", "+(prevSS==null ? false : prevSS.semiperfect); //***
//							assert(Read.CHECKSITE(ss, bases));
							if(verbose){System.err.println("A) Index made SiteScore "+ss.toText()+", "+Arrays.toString(ss.gaps));}
							assert(!perfect1 || ss.stop-ss.start==bases.length-1);
						}
						assert(!perfect2 || prevSS.stop-prevSS.start==bases.length-1);
					}else if(inbounds){
						if(verbose){System.err.println("Considering new site chr"+chrom+", "+site2+"-"+site3);}
						ss=new SiteScore(chrom, strand, site2, site3, approxHits, score, false, perfect1);
						if(!perfect1){ss.setPerfect(bases);}
						//assert((ss==null || !ss.semiperfect) && (prevSS==null || !prevSS.semiperfect)) : (ss==null ? false : ss.semiperfect)+", "+(prevSS==null ? false : prevSS.semiperfect); //***
//						assert(Read.CHECKSITE(ss, bases));
						ss.gaps=gapArray;
						if(verbose){System.err.println("B) Index made SiteScore "+ss.toText()+", "+Arrays.toString(ss.gaps));}
					}
					//assert((ss==null || !ss.semiperfect) && (prevSS==null || !prevSS.semiperfect)) : (ss==null ? false : ss.semiperfect)+", "+(prevSS==null ? false : prevSS.semiperfect); //***
					assert(ss==null || !ss.perfect || ss.semiperfect) : ss;
					assert(prevSS==null || !prevSS.perfect || prevSS.semiperfect) : "\n"+SiteScore.header()+"\n"+ss+"\n"+prevSS;
					if(ss!=null && ((SEMIPERFECTMODE && !ss.semiperfect) || (PERFECTMODE && !ss.perfect))){ss=null;}
					
					
					if(ss!=null){
//						System.out.println("Added site "+ss.toText()+", qscore="+qscore);
						ssl.add(ss);
						if(ss.perfect){
							
							if(prevSS==null || !prevSS.perfect || !ss.overlaps(prevSS)){
								if(prevSS==null){assert ssl.size()<2 || !ss.overlaps(ssl.get(ssl.size()-2));}
								perfectsFound++;

								//Human-specific code
//								if(QUIT_AFTER_TWO_PERFECTS){
//									if(perfectsFound>=3 || (perfectsFound>=2 && chrom<24)){break;}
//								}

								if(QUIT_AFTER_TWO_PERFECTS && perfectsFound>=2){break;}
							}
						}
						
						prevSS=ss;
					}else{
//						System.out.println("Subsumed site "+new SiteScore(chrom, strand, site2, site3, score).toText());
					}
					//assert((ss==null || !ss.semiperfect) && (prevSS==null || !prevSS.semiperfect)) : (ss==null ? false : ss.semiperfect)+", "+(prevSS==null ? false : prevSS.semiperfect); //***
				}
			}

			while(heap.peek().site==site){ //Remove all identical elements, and add subsequent elements
				final Quad t2=heap.poll();
				final int row=t2.row+1, col=t2.column;
				if(row<stops[col]){
					t2.row=row;
					
					int a=t2.list[row];
					int a2;
					if((a&SITE_MASK)>=offsets[col]){
						a2=a-offsets[col];
						
						assert(numberToChrom(a, baseChrom) == numberToChrom(a2, baseChrom)) :
							"baseChrom="+baseChrom+", chrom="+numberToChrom(a, baseChrom)+", strand="+strand+", site="+site+
							", maxNearbySite="+maxNearbySite+", a="+a+", a2="+a2+", offsets["+col+"]="+offsets[col];
					}else{
						int ch=numberToChrom(a, baseChrom);
						int st=numberToSite(a);
						int st2=Tools.max(st-offsets[col], 0);
						a2=toNumber(st2, ch);
						
						assert(numberToChrom(a, baseChrom) == numberToChrom(a2, baseChrom)) :
							"baseChrom="+baseChrom+", chrom="+numberToChrom(a, baseChrom)+", strand="+strand+", site="+site+
							", maxNearbySite="+maxNearbySite+", a="+a+", a2="+a2+", offsets["+col+"]="+offsets[col];
					}
					
					assert(numberToChrom(a, baseChrom) == numberToChrom(a2, baseChrom)) :
						"baseChrom="+baseChrom+", chrom="+numberToChrom(a, baseChrom)+", strand="+strand+", site="+site+
						", maxNearbySite="+maxNearbySite+", a="+a+", a2="+a2+", offsets["+col+"]="+offsets[col];
					
					t2.site=a2;
					values[col]=a2;
					heap.add(t2);
				}else if(heap.size()<approxHitsCutoff || PERFECTMODE){
					assert(USE_EXTENDED_SCORE);
					bestScores[0]=Tools.max(bestScores[0], currentTopScore);
					bestScores[1]=Tools.max(bestScores[1], maxHits);
					bestScores[2]=Tools.max(bestScores[2], qcutoff);
					bestScores[3]=Tools.max(bestScores[3], bestqscore);
					
					bestScores[4]=maxQuickScore;
					bestScores[5]=perfectsFound; //***$ fixed by adding this line
					if(!RETAIN_BEST_QCUTOFF){bestScores[2]=-9999;}
					
					return ssl;
				}
				if(heap.isEmpty()){
					assert(false) : heap.size()+", "+approxHitsCutoff;
					break;
				}
			}

		}
		
		assert(USE_EXTENDED_SCORE);
		bestScores[0]=Tools.max(bestScores[0], currentTopScore);
		bestScores[1]=Tools.max(bestScores[1], maxHits);
		bestScores[2]=Tools.max(bestScores[2], qcutoff);
		bestScores[3]=Tools.max(bestScores[3], bestqscore);

		bestScores[4]=maxQuickScore;
		bestScores[5]=perfectsFound;
		if(!RETAIN_BEST_QCUTOFF){bestScores[2]=-9999;}
		
		return ssl;
	}
	
	
	/** Uses dual heaps - reserve (big) and active (small). */
	private final ArrayList<SiteScore> camelWalk3(int[] starts, int[] stops, final byte[] bases,
			final byte[] baseScores, int[] keyScores, int[] offsets,
			final int baseChrom_, final byte strand, final boolean obeyLimits, ArrayList<SiteScore> ssl,
			int[] bestScores, final boolean allBasesCovered, final int maxScore, final boolean fullyDefined){
		assert(USE_EXTENDED_SCORE);
		
		final int numKeys=offsets.length; //Before shrink
		
		//This can be done before or after shrinking, but the results will change depending on MIN_SCORE_MULT and etc.
		final int maxQuickScore=maxQuickScore(offsets, keyScores);
		
		if(SHRINK_BEFORE_WALK){
			int[][] r=shrink(starts, stops, offsets, keyScores, offsets.length);
			if(r!=null){
				starts=r[0];
				stops=r[1];
				offsets=r[2];
				keyScores=r[4];
			}
		}
		
		final int numHits=offsets.length; //After shrink
		
		
		assert(numHits==offsets.length);
		assert(numHits==keyScores.length);
		
		usedKeys+=numHits;
		usedKeyIterations++;
		
		final boolean filter_by_qscore=(FILTER_BY_QSCORE && numKeys>=5);
		
		assert(!(!SHRINK_BEFORE_WALK && ADD_SCORE_Z));
		
		
//		final int minScore=(obeyLimits ? (int)(MIN_SCORE_MULT*maxScore) : (int)(MIN_SCORE_MULT*0.85f*maxScore));
		final int minScore=(obeyLimits ? (int)(MIN_SCORE_MULT*maxScore) : (int)(MIN_SCORE_MULT*1.25f*maxScore));
		final int minQuickScore=(int)(MIN_QSCORE_MULT*maxQuickScore);
		
		final int baseChrom=baseChrom(baseChrom_);
		
		heap.clear();
		active.clear();
		
		final Quad[] triples=tripleStorage;

		final int[] values=valueArray;
		final int[] sizes=sizeArray;
		final int[] locArray=(USE_EXTENDED_SCORE ? getLocArray(bases.length) : null);
		final Block b=index[baseChrom];
		
		if(ssl==null){ssl=new ArrayList<SiteScore>(8);}
		
		int currentTopScore=bestScores[0];
		int cutoff=Tools.max(minScore, (int)(currentTopScore*DYNAMIC_SCORE_THRESH));
		
		int qcutoff=Tools.max(bestScores[2], minQuickScore);
		int bestqscore=bestScores[3];
		int maxHits=bestScores[1];
		int perfectsFound=bestScores[5];
		assert((currentTopScore>=maxScore) == (perfectsFound>0)) : currentTopScore+", "+maxScore+", "+perfectsFound+", "+maxHits+", "+numHits;
		int approxHitsCutoff=calcApproxHitsCutoff(numKeys, maxHits, MIN_APPROX_HITS_TO_KEEP, currentTopScore>=maxScore);
		if(approxHitsCutoff>numHits){return ssl;}
		
		final boolean shortCircuit=(allBasesCovered && numKeys==numHits && filter_by_qscore);
		
		if(currentTopScore>=maxScore){
			assert(currentTopScore==maxScore);
			qcutoff=Tools.max(qcutoff, (int)(maxQuickScore*DYNAMIC_QSCORE_THRESH_PERFECT));
		}
			
		
		for(int i=0; i<numHits; i++){
			final int[] sites=b.sites;
			final int start=starts[i];
			sizes[i]=b.length(start, stops[i]);
			assert(sizes[i]>0);
			
			int a=sites[start];
			int a2;
			if((a&SITE_MASK)>=offsets[i]){
				a2=a-offsets[i];
			}else{
				int ch=numberToChrom(a, baseChrom);
				int st=numberToSite(a);
				int st2=Tools.max(st-offsets[i], 0);
				a2=toNumber(st2, ch);
			}
			assert(numberToChrom(a, baseChrom) == numberToChrom(a2, baseChrom));

			Quad t=triples[i];
			assert(t!=null) : "Should be using tripleStorage";
			assert(i==t.column);
			t.row=start;
			t.site=a2;
			t.list=sites;
			values[i]=a2;

			heap.add(t);
		}
		
		assert(numHits>0);
		assert(heap.size()==numHits);
		
		/* Tracks largest element allowed in 'active' */
		
//		System.err.println("\nEntering SS loop:");
//		System.err.println("maxScore="+maxScore+"\tminScore="+minScore+"\tcurrentTopScore="+currentTopScore+"\n" +
//				"cutoff="+cutoff+"\tmaxHits="+maxHits+"\tapproxHitsCutoff="+approxHitsCutoff);
//		System.err.println("maxQuickScore="+maxQuickScore+"\tminQuickScore="+minQuickScore+"\tqcutoff="+qcutoff);
		
//		int iter=0;
		SiteScore prevSS=null;
		int maxNearbySite=0;
		int site=0;
		int horizon=0;
		assert(active.isEmpty());
		while(!heap.isEmpty() || !active.isEmpty()){
//			iter++;
			
			do{
				while(!active.isEmpty() && active.peek().site==site){ //Remove all identical elements, and add subsequent elements
					final Quad t2=active.poll();
					final int row=t2.row+1, col=t2.column;

					//This is called the "increment" operation.  Very messy and slow due to rare cases at beginning of a chrom.
					if(row<stops[col]){//then increment and return to the heap
						t2.row=row;

						int a=t2.list[row];
						int a2;
						if((a&SITE_MASK)>=offsets[col]){
							a2=a-offsets[col];

							assert(numberToChrom(a, baseChrom) == numberToChrom(a2, baseChrom)) :
								"baseChrom="+baseChrom+", chrom="+numberToChrom(a, baseChrom)+", strand="+strand+", site="+site+
								", maxNearbySite="+maxNearbySite+", a="+a+", a2="+a2+", offsets["+col+"]="+offsets[col];
						}else{
							int ch=numberToChrom(a, baseChrom);
							int st=numberToSite(a);
							int st2=Tools.max(st-offsets[col], 0);
							a2=toNumber(st2, ch);

							assert(numberToChrom(a, baseChrom) == numberToChrom(a2, baseChrom)) :
								"baseChrom="+baseChrom+", chrom="+numberToChrom(a, baseChrom)+", strand="+strand+", site="+site+
								", maxNearbySite="+maxNearbySite+", a="+a+", a2="+a2+", offsets["+col+"]="+offsets[col];
						}

						assert(numberToChrom(a, baseChrom) == numberToChrom(a2, baseChrom)) :
							"baseChrom="+baseChrom+", chrom="+numberToChrom(a, baseChrom)+", strand="+strand+", site="+site+
							", maxNearbySite="+maxNearbySite+", a="+a+", a2="+a2+", offsets["+col+"]="+offsets[col];

						t2.site=a2;
						values[col]=a2;
						if(a2<=horizon){
							active.add(t2);
							maxNearbySite=Tools.max(t2.site, maxNearbySite);
						}else{heap.add(t2);}
					}else if((heap.size()+active.size())<approxHitsCutoff || PERFECTMODE){ //Then there are not enough keys remaining for a site
						assert(USE_EXTENDED_SCORE);
						bestScores[0]=Tools.max(bestScores[0], currentTopScore);
						bestScores[1]=Tools.max(bestScores[1], maxHits);
						bestScores[2]=Tools.max(bestScores[2], qcutoff);
						bestScores[3]=Tools.max(bestScores[3], bestqscore);

						bestScores[4]=maxQuickScore;
						bestScores[5]=perfectsFound; //***$ fixed by adding this line
						if(!RETAIN_BEST_QCUTOFF){bestScores[2]=-9999;}

						return ssl;
					}
//					if(heap.isEmpty() && active.isEmpty()){
//						assert(false) : heap.size()+", "+active.size()+", "+approxHitsCutoff;
//						break;
//					}
				}
				
				final Quad t;
				if(active.isEmpty()){
					t=heap.poll();
					active.add(t);
					maxNearbySite=t.site;
				}else{
					t=active.peek();
				}
				site=t.site;
				horizon=(int)Tools.min(Integer.MAX_VALUE, site+(long)MAX_INDEL2);
				while(!heap.isEmpty() && heap.peek().site<=horizon){
					Quad t2=heap.poll();
					active.add(t2);
					maxNearbySite=t2.site;
				}
				
				if(verbose){System.err.println("\nFinished loop iteration. active="+active+"\nheap="+heap+"\nfloor="+site+", horizon="+horizon+", maxNearbySite="+maxNearbySite);}
				
			}while(active.size()<approxHitsCutoff);
			
			final int approxHits=active.size();
			final Quad t=active.peek();
			final int centerIndex=t.column;
			
			if(verbose){
				System.err.println("\nLeft loop. active="+active+"\nheap="+heap+"\n" +
						"floor="+site+", maxNearbySite="+maxNearbySite+", approxHits="+approxHits+", approxHitsCutoff="+approxHitsCutoff+"\n");
			}
			
			
			
//			approxHits=0;
//			{//Inner loop
//				final int minsite=site, maxsite=site+MAX_INDEL2;
//				for(int column=0, chances=numHits-approxHitsCutoff; column<numHits && chances>=0; column++){
//					final int x=values[column];
//					assert(x==triples[column].site);
//					if(x>=minsite && x<=maxsite){
//						maxNearbySite=(x>maxNearbySite ? x : maxNearbySite);
//						approxHits++;
//					}else{chances--;}
////					if(verbose){
////						System.err.println("column="+column+", numHits="+numHits+", approxHits="+approxHits+
////								", approxHitsCutoff="+approxHitsCutoff+", chances="+chances);
////					}
//				}
//				//Invalid assertion due to loop early exit
////				assert(approxHits>0) : "\niter="+iter+", maxHits="+maxHits+", numHits="+numHits+", approxHitsCutoff="+approxHitsCutoff+
////					"\nheap.size()="+heap.size()+", minsite="+minsite+", maxsite="+maxsite+", values[center]="+values[centerIndex]+", t="+t;
//			}
//			assert(approxHits<=active.size()) : "approxHits="+approxHits+", active.size()="+active.size()+", maxNearbySite="+maxNearbySite+"\nvalues="+Arrays.toString(values);
			
			hist_hits[Tools.min(HIT_HIST_LEN, approxHits)]++;

			assert(centerIndex>=0) : centerIndex;
			
			//I don't remember what this assertion was for or why, but it's causing trouble.
			//assert(approxHits>=1 || approxHitsCutoff>1) : approxHits+", "+approxHitsCutoff+", "+numHits+", "+t.column;
			if(approxHits>=approxHitsCutoff){
				
//				if(verbose){System.err.println("A");}
				
				int score;
				int qscore=(filter_by_qscore ? quickScore(values, keyScores, centerIndex, offsets, sizes, true, approxHits, numHits) : qcutoff);
				if(ADD_SCORE_Z){
					int scoreZ=scoreZ2(values, centerIndex, offsets, approxHits, numHits);
					qscore+=scoreZ;
				}
				
				int mapStart=site, mapStop=maxNearbySite;
				assert(mapStart<=mapStop);
				
				assert(USE_EXTENDED_SCORE);
				
				boolean locArrayValid=false;
				if(qscore<qcutoff){
					score=-1;
//					if(verbose){System.err.println("B");}
				}else{
//					if(verbose){System.err.println("C");}
					final int chrom=numberToChrom(site, baseChrom);
					
					//TODO Note that disabling the shortCircuit code seems to make things run 2% faster (with identical results).
					//However, theoretically, shortCircuit should be far more efficient.  Test both ways on cluster and on a larger run.
					//May have something to do with compiler loop optimizations.
					if(shortCircuit && qscore==maxQuickScore){
//						if(verbose){System.err.println("D");}
						assert(approxHits==numKeys);
						score=maxScore;
					}else{
//						if(verbose){System.err.println("E");}
						if(verbose){
							System.err.println("numHits="+numHits+", approxHits="+approxHits+", qscore="+qscore+", qcutoff="+qcutoff+", filter_by_qscore="+filter_by_qscore);
							System.err.println("Extending "+Arrays.toString(values));
						}
						score=extendScore(bases, baseScores, offsets, values, chrom, centerIndex, locArray, numHits, approxHits);
						locArrayValid=true;

						if(verbose){
							System.err.println("score: "+score);
							System.err.println("locArray: "+Arrays.toString(locArray));
						}
					
						//Correct begin and end positions if they changed.
						int min=Integer.MAX_VALUE;
						int max=Integer.MIN_VALUE;
						for(int i=0; i<locArray.length; i++){
							int x=locArray[i];
							if(x>-1){
								if(x<min){min=x;}
								if(x>max){max=x;}
							}
						}

						if(score>=maxScore){
							assert(min==max && min>-1) : "\n"+score+", "+maxScore+", "+min+", "+max+
							", "+(max-min)+", "+bases.length+"\n"+Arrays.toString(locArray)+"\n";
						}
						
						//							assert(min>-1 && max>-1) : Arrays.toString(locArray); //TODO: How did this assertion trigger?
						if(min<0 || max<0){
							if(!Shared.anomaly){
								Shared.anomaly=true;
								System.err.println("\nAnomaly in "+getClass().getName()+".slowWalk:\n"+
										"chrom="+chrom+", mapStart="+mapStart+", mapStop="+mapStop+", centerIndex="+centerIndex+", strand="+strand+"\n"+
										"score="+score+", maxScore="+maxScore+", qscore="+qscore+", filter_by_qscore="+filter_by_qscore+"\n"+
										"numHits="+approxHits+", approxHits="+approxHits+"\n"+
										"min="+min+", max="+max+", (max-min)="+(max-min)+"\n"+
										"bases.length="+bases.length+"\n"+Arrays.toString(locArray)+"\n"+
										"locArray:\t"+Arrays.toString(locArray)+"\n"+
										"values:\t"+Arrays.toString(values)+"\n"+
										"bases:\t"+new String(bases));
								System.err.println();
								assert(false);
							}
							score=-99999;
						}

						//mapStart and mapStop are indices
						mapStart=toNumber(min, chrom);
						mapStop=toNumber(max, chrom);
						assert(mapStart<=mapStop);

						if(score>=maxScore){
							assert(mapStop-mapStart==0) : "\n"+score+", "+maxScore+", "+min+", "+max+
							", "+(max-min)+", "+(mapStop-mapStart)+", "+bases.length+"\n"+Arrays.toString(locArray)+"\n";
						}
					}
//					if(verbose){System.err.println("F");}

					if(score==maxScore){
						qcutoff=Tools.max(qcutoff, (int)(maxQuickScore*DYNAMIC_QSCORE_THRESH_PERFECT));
						approxHitsCutoff=calcApproxHitsCutoff(numKeys, maxHits, MIN_APPROX_HITS_TO_KEEP, true);
					}
					
					if(score>=cutoff){
						qcutoff=Tools.max(qcutoff, (int)(qscore*DYNAMIC_QSCORE_THRESH));
						bestqscore=Tools.max(qscore, bestqscore);
					}
				}
//				if(verbose){System.err.println("G");}
				
				if(score>=cutoff){
//					if(verbose){System.err.println("H");}

					if(score>currentTopScore){
//						if(verbose){System.err.println("I");}
//						System.err.println("New top score!");

						if(DYNAMICALLY_TRIM_LOW_SCORES){

							maxHits=Tools.max(approxHits, maxHits);
							approxHitsCutoff=calcApproxHitsCutoff(numKeys, maxHits, approxHitsCutoff, currentTopScore>=maxScore);

							cutoff=Tools.max(cutoff, (int)(score*DYNAMIC_SCORE_THRESH));
							if(score>=maxScore){
								assert(USE_EXTENDED_SCORE);
								cutoff=Tools.max(cutoff, (int)(score*0.95f));
							}
						}

						currentTopScore=score;

//						System.err.println("New top score: "+currentTopScore+" \t("+cutoff+")");
					}
//					if(verbose){System.err.println("J");}

					final int chrom=numberToChrom(mapStart, baseChrom);
					final int site2=numberToSite(mapStart);
					final int site3=numberToSite(mapStop)+bases.length-1;
					
					assert(NUM_CHROM_BITS==0 || site2<SITE_MASK-1000) : "chrom="+chrom+", strand="+strand+", site="+site+
						", maxNearbySite="+maxNearbySite+", site2="+site2+", site3="+site3+", read.length="+bases.length+
						"\n\n"+Arrays.toString(b.getHitList(centerIndex));
					assert(site2<site3) : "\nchrom="+chrom+", strand="+strand+", site="+site+", maxNearbySite="+maxNearbySite+"\n"+
						"mapStart="+mapStart+", mapStop="+mapStop+", site2="+site2+", site3="+site3+", read.length="+bases.length+"\n"+
						"numberToSite("+mapStart+")="+numberToSite(mapStart)+", numberToSite("+mapStop+")="+numberToSite(mapStop)+"\n"+
						"\n"+new String(bases)+"\n";
					
					if(verbose){
						System.err.println("chrom="+chrom+", strand="+strand+", site="+site+", maxNearbySite="+maxNearbySite+"\n"+
								"mapStart="+mapStart+", mapStop="+mapStop+", site2="+site2+", site3="+site3+", read.length="+bases.length);
					}
					
					int[] gapArray=null;
					if(site3-site2>=MINGAP+bases.length){
						assert(locArrayValid) : "Loc array was not filled.";
//						System.err.println("****\n"+Arrays.toString(locArray)+"\n");
//						int[] clone=locArray.clone();
						gapArray=makeGapArray(locArray, site2, MINGAP);
						if(gapArray!=null){
//							System.err.println(Arrays.toString(locArray)+"\n");
//							System.err.println(Arrays.toString(gapArray));
//
////							int sub=site2-mapStart;//thus site2=mapStart+sub
////							for(int i=0; i<gapArray.length; i++){
////								gapArray[i]+=sub;
////							}
////							System.err.println(Arrays.toString(gapArray));
//
//							System.err.println(mapStart+" -> "+site2);
//							System.err.println(mapStop+" -> "+site3);
							
							assert(gapArray[0]>=site2 && gapArray[0]-site2<bases.length);
							assert(gapArray[gapArray.length-1]<=site3 && site3-gapArray[gapArray.length-1]<bases.length) : "\n"+
								mapStart+" -> "+site2+"\n"+
								mapStop+" -> "+site3+"\n\n"+
								Arrays.toString(gapArray)+"\n\n"+
//								Arrays.toString(clone)+"\n\n"+
								Arrays.toString(locArray)+"\n"+
								"numHits="+numHits+", "+
								"heap.size="+heap.size()+", "+
								"numHits="+numHits+", "+
								"approxHits="+approxHits+"\n";
							gapArray[0]=Tools.min(gapArray[0], site2);
							gapArray[gapArray.length-1]=Tools.max(gapArray[gapArray.length-1], site3);
						}
						if(verbose){System.err.println("@ site "+site2+", made gap array: "+Arrays.toString(gapArray));}
//						assert(false) : Arrays.toString(locArray);
					}
					
					
					//This block is optional, but tries to eliminate multiple identical alignments
					
					SiteScore ss=null;
					final boolean perfect1=USE_EXTENDED_SCORE && score==maxScore && fullyDefined;
					final boolean inbounds=(site2>=0 && site3<Data.chromLengths[chrom]);
//					if(!inbounds){System.err.println("Index tossed out-of-bounds site chr"+chrom+", "+site2+"-"+site3);}
					
					if(inbounds && !SEMIPERFECTMODE && !PERFECTMODE && gapArray==null && prevSS!=null &&
							prevSS.chrom==chrom && prevSS.strand==strand && overlap(prevSS.start, prevSS.stop, site2, site3)){

						final int betterScore=Tools.max(score, prevSS.score);
						final int minStart=Tools.min(prevSS.start, site2);
						final int maxStop=Tools.max(prevSS.stop, site3);
						final boolean perfect2=USE_EXTENDED_SCORE && prevSS.score==maxScore && fullyDefined;
						assert(!USE_EXTENDED_SCORE || perfect2==prevSS.perfect);

						final boolean shortEnough=(!LIMIT_SUBSUMPTION_LENGTH_TO_2X  || (maxStop-minStart<2*bases.length));

						if(prevSS.start==site2 && prevSS.stop==site3){
							prevSS.score=prevSS.quickScore=betterScore;
							prevSS.perfect=(prevSS.perfect || perfect1 || perfect2);
							if(prevSS.perfect){prevSS.semiperfect=true;}
						}else if(SUBSUME_SAME_START_SITES && shortEnough && prevSS.start==site2 && !prevSS.semiperfect){
							if(perfect2){
								//do nothing
							}else if(perfect1){
								prevSS.setStop(site3);
								if(!prevSS.perfect){perfectsFound++;}//***$
								prevSS.perfect=prevSS.semiperfect=true;
							}else{
								prevSS.setStop(maxStop);
								prevSS.setPerfect(bases);
							}
							prevSS.score=prevSS.quickScore=betterScore;
						}else if(SUBSUME_SAME_STOP_SITES && shortEnough && prevSS.stop==site3 && !prevSS.semiperfect){
							if(perfect2){
								//do nothing
							}else if(perfect1){
								prevSS.setStart(site2);
								if(!prevSS.perfect){perfectsFound++;}//***$
								prevSS.perfect=prevSS.semiperfect=true;
							}else{
								prevSS.setStart(minStart);
								prevSS.setPerfect(bases);
							}
							prevSS.score=prevSS.quickScore=betterScore;
						}else if(SUBSUME_OVERLAPPING_SITES && shortEnough && (maxStop-minStart<=bases.length+MAX_SUBSUMPTION_LENGTH)
								&& !perfect1 && !perfect2 && !prevSS.semiperfect){
							prevSS.setLimits(minStart, maxStop);
							prevSS.score=prevSS.quickScore=betterScore;
							prevSS.setPerfect(bases);
						}else{
							ss=new SiteScore(chrom, strand, site2, site3, approxHits, score, false, perfect1);
							if(!perfect1){ss.setPerfect(bases);}
//							assert(Read.CHECKSITE(ss, bases));
							if(verbose){System.err.println("A) Index made SiteScore "+ss.toText()+", "+Arrays.toString(ss.gaps));}
							assert(!perfect1 || ss.stop-ss.start==bases.length-1);
						}
						assert(!perfect2 || prevSS.stop-prevSS.start==bases.length-1);
					}else if(inbounds){
						ss=new SiteScore(chrom, strand, site2, site3, approxHits, score, false, perfect1);
						if(!perfect1){ss.setPerfect(bases);}
//						assert(Read.CHECKSITE(ss, bases));
						ss.gaps=gapArray;
						if(verbose){System.err.println("B) Index made SiteScore "+ss.toText()+", "+Arrays.toString(ss.gaps));}
					}
					
					assert(ss==null || !ss.perfect || ss.semiperfect) : ss;
					assert(prevSS==null || !prevSS.perfect || prevSS.semiperfect) : "\n"+SiteScore.header()+"\n"+ss+"\n"+prevSS;
					if(ss!=null && ((SEMIPERFECTMODE && !ss.semiperfect) || (PERFECTMODE && !ss.perfect))){ss=null;}
					
					
					if(ss!=null){
//						System.err.println("Added site "+ss.toText()+", qscore="+qscore);
						ssl.add(ss);
						if(ss.perfect){
							
							if(prevSS==null || !prevSS.perfect || !ss.overlaps(prevSS)){
								if(prevSS==null){assert ssl.size()<2 || !ss.overlaps(ssl.get(ssl.size()-2));}
								perfectsFound++;

								//Human-specific code
//								if(QUIT_AFTER_TWO_PERFECTS){
//									if(perfectsFound>=3 || (perfectsFound>=2 && chrom<24)){break;}
//								}

								if(QUIT_AFTER_TWO_PERFECTS && perfectsFound>=2){break;}
							}
						}
						
						prevSS=ss;
					}else{
//						System.err.println("Subsumed site "+new SiteScore(chrom, strand, site2, site3, score).toText());
					}
				}
			}
			
//			while(!active.isEmpty() && active.peek().site==site){ //Remove all identical elements, and add subsequent elements
//				final Quad t2=active.poll();
//				final int row=t2.row+1, col=t2.column;
//
//				//This is called the "increment" operation.  Very messy and slow due to rare cases at beginning of a chrom.
//				if(row<stops[col]){//then increment and return to the heap
//					t2.row=row;
//
//					int a=t2.list[row];
//					int a2;
//					if((a&SITE_MASK)>=offsets[col]){
//						a2=a-offsets[col];
//
//						assert(numberToChrom(a, baseChrom) == numberToChrom(a2, baseChrom)) :
//							"baseChrom="+baseChrom+", chrom="+numberToChrom(a, baseChrom)+", strand="+strand+", site="+site+
//							", maxNearbySite="+maxNearbySite+", a="+a+", a2="+a2+", offsets["+col+"]="+offsets[col];
//					}else{
//						int ch=numberToChrom(a, baseChrom);
//						int st=numberToSite(a);
//						int st2=Tools.max(st-offsets[col], 0);
//						a2=toNumber(st2, ch);
//
//						assert(numberToChrom(a, baseChrom) == numberToChrom(a2, baseChrom)) :
//							"baseChrom="+baseChrom+", chrom="+numberToChrom(a, baseChrom)+", strand="+strand+", site="+site+
//							", maxNearbySite="+maxNearbySite+", a="+a+", a2="+a2+", offsets["+col+"]="+offsets[col];
//					}
//
//					assert(numberToChrom(a, baseChrom) == numberToChrom(a2, baseChrom)) :
//						"baseChrom="+baseChrom+", chrom="+numberToChrom(a, baseChrom)+", strand="+strand+", site="+site+
//						", maxNearbySite="+maxNearbySite+", a="+a+", a2="+a2+", offsets["+col+"]="+offsets[col];
//
//					t2.site=a2;
//					values[col]=a2;
//					if()
//					heap.add(t2);
//				}else if((heap.size()+active.size())<approxHitsCutoff || PERFECTMODE){ //Then there are not enough keys remaining for a site
//					assert(USE_EXTENDED_SCORE);
//					bestScores[0]=Tools.max(bestScores[0], currentTopScore);
//					bestScores[1]=Tools.max(bestScores[1], maxHits);
//					bestScores[2]=Tools.max(bestScores[2], qcutoff);
//					bestScores[3]=Tools.max(bestScores[3], bestqscore);
//
//					bestScores[4]=maxQuickScore;
//					bestScores[5]=perfectsFound; //***$ fixed by adding this line
//					if(!RETAIN_BEST_QCUTOFF){bestScores[2]=-9999;}
//
//					return ssl;
//				}
//				if(heap.isEmpty() && active.isEmpty()){
//					assert(false) : heap.size()+", "+active.size()+", "+approxHitsCutoff;
//					break;
//				}
//			}

		}
		
		assert(USE_EXTENDED_SCORE);
		bestScores[0]=Tools.max(bestScores[0], currentTopScore);
		bestScores[1]=Tools.max(bestScores[1], maxHits);
		bestScores[2]=Tools.max(bestScores[2], qcutoff);
		bestScores[3]=Tools.max(bestScores[3], bestqscore);

		bestScores[4]=maxQuickScore;
		bestScores[5]=perfectsFound;
		if(!RETAIN_BEST_QCUTOFF){bestScores[2]=-9999;}
		
		return ssl;
	}
	
	
	private final int[] findMaxQscore2(final int[] starts, final int[] stops, final int[] offsets, final int[] keyScores,
			final int baseChrom_, final Quad[] triples, final int[] values, final int prevMaxHits,
			boolean earlyExit, boolean perfectOnly){
		
		final int numHits=offsets.length;
		assert(numHits>=prevMaxHits);
		
		final int baseChrom=baseChrom(baseChrom_);
		final Block b=index[baseChrom];
		final int[] sizes=sizeArray;
		
		heap.clear();
		for(int i=0; i<numHits; i++){
			final int[] sites=b.sites;
			final int start=starts[i];
			sizes[i]=b.length(start, stops[i]);
			assert(sizes[i]>0);
			
			int a=sites[start];
			int a2;
			if((a&SITE_MASK)>=offsets[i]){
				a2=a-offsets[i];
			}else{
				int ch=numberToChrom(a, baseChrom);
				int st=numberToSite(a);
				int st2=Tools.max(st-offsets[i], 0);
				a2=toNumber(st2, ch);
			}
			assert(numberToChrom(a, baseChrom) == numberToChrom(a2, baseChrom));

			Quad t=triples[i];
			assert(t!=null) : "Should be using tripleStorage";
			assert(i==t.column);
			t.row=start;
			t.site=a2;
			t.list=sites;
			values[i]=a2;

			heap.add(t);
		}
		
		final int maxQuickScore=maxQuickScore(offsets, keyScores);
		
		int topQscore=-999999999;
		
		int maxHits=0;
//		int approxHitsCutoff=MIN_APPROX_HITS_TO_KEEP;
		
		
		int approxHitsCutoff;
		final int indelCutoff;
		if(perfectOnly){
			approxHitsCutoff=numHits;
			indelCutoff=0;
		}else{
			approxHitsCutoff=Tools.max(prevMaxHits, Tools.min(MIN_APPROX_HITS_TO_KEEP, numHits-1)); //Faster, same accuracy
			indelCutoff=MAX_INDEL2;
		}
		
		
		while(!heap.isEmpty()){
			Quad t=heap.peek();
			final int site=t.site;
			final int centerIndex=t.column;
			
			int maxNearbySite=site;


			int approxHits=0;
			{//Inner loop
				final int minsite=site-Tools.min(MAX_INDEL, indelCutoff), maxsite=site+MAX_INDEL2;
				for(int column=0, chances=numHits-approxHitsCutoff; column<numHits && chances>=0; column++){
					final int x=values[column];
					assert(x==triples[column].site);
					if(x>=minsite && x<=maxsite){
						maxNearbySite=(x>maxNearbySite ? x : maxNearbySite);
						approxHits++;
					}else{chances--;}
				}
			}
			hist_hits[Tools.min(HIT_HIST_LEN, approxHits)]++;

			assert(centerIndex>=0) : centerIndex;
			
			//I don't remember what this assertion was for or why, but it's causing trouble.
			//assert(approxHits>=1 || approxHitsCutoff>1) : approxHits+", "+approxHitsCutoff+", "+numHits+", "+t.column;
			if(approxHits>=approxHitsCutoff){
				
				int qscore=quickScore(values, keyScores, centerIndex, offsets, sizes, true, approxHits, numHits);
				
				if(ADD_SCORE_Z){
					int scoreZ=scoreZ2(values, centerIndex, offsets, approxHits, numHits);
					qscore+=scoreZ;
				}
				
				if(qscore>topQscore){
					
//					maxHits=Tools.max(approxHits, maxHits);
//					approxHitsCutoff=Tools.max(approxHitsCutoff, maxHits); //Best setting for pre-scan
					
					maxHits=Tools.max(approxHits, maxHits);
					approxHitsCutoff=Tools.max(approxHitsCutoff, approxHits-1); //Best setting for pre-scan
					
					topQscore=qscore;
					
					if(qscore>=maxQuickScore){
						assert(qscore==maxQuickScore);
						assert(approxHits==numHits);
						if(earlyExit){
							return new int[] {topQscore, maxHits};
						}
					}
				}
			}

			while(heap.peek().site==site){ //Remove all identical elements, and add subsequent elements
				final Quad t2=heap.poll();
				final int row=t2.row+1, col=t2.column;
				if(row<stops[col]){
					t2.row=row;
					
					int a=t2.list[row];
					int a2;
					if((a&SITE_MASK)>=offsets[col]){
						a2=a-offsets[col];
						
						assert(numberToChrom(a, baseChrom) == numberToChrom(a2, baseChrom)) :
							"baseChrom="+baseChrom+", chrom="+numberToChrom(a, baseChrom)+", site="+site+
							", maxNearbySite="+maxNearbySite+", a="+a+", a2="+a2+", offsets["+col+"]="+offsets[col];
					}else{
						int ch=numberToChrom(a, baseChrom);
						int st=numberToSite(a);
						int st2=Tools.max(st-offsets[col], 0);
						a2=toNumber(st2, ch);
						
						assert(numberToChrom(a, baseChrom) == numberToChrom(a2, baseChrom)) :
							"baseChrom="+baseChrom+", chrom="+numberToChrom(a, baseChrom)+", site="+site+
							", maxNearbySite="+maxNearbySite+", a="+a+", a2="+a2+", offsets["+col+"]="+offsets[col];
					}
					
					assert(numberToChrom(a, baseChrom) == numberToChrom(a2, baseChrom)) :
						"baseChrom="+baseChrom+", chrom="+numberToChrom(a, baseChrom)+", site="+site+
						", maxNearbySite="+maxNearbySite+", a="+a+", a2="+a2+", offsets["+col+"]="+offsets[col];

					t2.site=a2;
					values[col]=a2;
					heap.add(t2);
				}else if(earlyExit && (perfectOnly || heap.size()<approxHitsCutoff)){
					return new int[] {topQscore, maxHits};
				}
				if(heap.isEmpty()){break;}
			}

		}
		
		
		
		return new int[] {topQscore, maxHits};
	}
	
	
	private static final int absdif(int a, int b){
		return a>b ? a-b : b-a;
	}
	
	
	@Override
	final int maxScore(int[] offsets, byte[] baseScores, int[] keyScores, int readlen, boolean useQuality){
		
		if(useQuality){
			//These lines apparently MUST be used if quality is used later on for slow align.
			if(USE_AFFINE_SCORE){return msa.maxQuality(baseScores);}
			if(USE_EXTENDED_SCORE){return readlen*(BASE_HIT_SCORE+BASE_HIT_SCORE/5)+Tools.sumInt(baseScores);}
		}else{
			if(USE_AFFINE_SCORE){return msa.maxQuality(readlen);}
			if(USE_EXTENDED_SCORE){return readlen*(BASE_HIT_SCORE+BASE_HIT_SCORE/5);}
		}
		
		return maxQuickScore(offsets, keyScores);
	}
	
	
	public final int maxQuickScore(int[] offsets, int[] keyScores){

//		int x=offsets.length*BASE_KEY_HIT_SCORE;
		int x=Tools.intSum(keyScores);
		int y=Y_SCORE_MULT*(offsets[offsets.length-1]-offsets[0]);
//		if(ADD_LIST_SIZE_BONUS){x+=(LIST_SIZE_BONUS[1]*offsets.length);}
//		assert(!ADD_SCORE_Z) : "Need to make sure this is correct...";

//		if(ADD_SCORE_Z){x+=((offsets[offsets.length-1]+CHUNKSIZE)*Z_SCORE_MULT);}
		if(ADD_SCORE_Z){x+=maxScoreZ(offsets);}
		
		return x+y;
//		int bonus=(2*(HIT_SCORE/2)); //For matching both ends
//		return x+y+bonus;
	}
	
	
	private final int quickScore(final int[] locs, final int[] keyScores, final int centerIndex, final int offsets[],
			int[] sizes, final boolean penalizeIndels, final int numApproxHits, final int numHits){
		
		hist_hits_score[Tools.min(HIT_HIST_LEN, numApproxHits)]++;
		if(numApproxHits==1){return keyScores[centerIndex];}
		
		//Done!
		//Correct way to calculate score:
		//Find the first chunk that exactly hits the center.
		//Then, align leftward of it, and align rightward of it, and sum the scores.
		
		//"-centerIndex" is a disambiguating term that, given otherwise identical match patterns
		//(for example, a small indel will generate two valid site candidates), choose the lower site.

		int x=keyScores[centerIndex]+scoreLeft(locs, keyScores, centerIndex, sizes, penalizeIndels)+
			scoreRight(locs, keyScores, centerIndex, sizes, penalizeIndels, numHits)-centerIndex;
		
		int y=Y_SCORE_MULT*scoreY(locs, centerIndex, offsets);
		if(ADD_LIST_SIZE_BONUS){x+=calcListSizeBonus(sizes[centerIndex]);}
//		int z=scoreZ(locs, hits);
		return x+y;
	}
	
	
//	/** Generates a term that increases score with how many bases in the read match the ref. */
//	private static final int scoreZ(int[] locs, int centerIndex, int offsets[]){
//		final int center=locs[centerIndex];
//
//		final int[] refLoc=new int[offsets[offsets.length-1]+CHUNKSIZE];
//
//		final int maxLoc=center+MAX_INDEL2;
//		final int minLoc=Tools.max(0, center-MAX_INDEL);
//
//		int score=0;
//
//		for(int i=0; i<locs.length; i++){
//			int loc=locs[i];
////			int dif=absdif(loc, center);
//			if(loc>=minLoc && loc<=maxLoc){
////				assert(loc>=center) : "loc="+loc+"\ni="+i+"\ncenterIndex="+centerIndex+
////					"\nmaxLoc="+maxLoc+"\nlocs:\t"+Arrays.toString(locs)+"\noffsets:\t"+Arrays.toString(offsets);
//
//				int offset=offsets[i];
//				int max=CHUNKSIZE+offset;
//
//				for(int j=offset; j<max; j++){
//					int old=refLoc[j];
//					if(old==0){
//						refLoc[j]=loc;
//						score+=4;
//					}else if(old>loc){
//						refLoc[j]=loc;
//						score-=2;
//					}else if(old==loc){
//						score-=1;
//						//do nothing, perhaps, or add 1?
//					}else{
//						score-=2;
//						assert(old<loc);
//					}
//				}
//			}
//		}
//		return score;
//	}
	
	
	
	private final int extendScore(final byte[] bases, final byte[] baseScores, final int[] offsets, final int[] values,
			final int chrom, final int centerIndex, final int[] locArray, final int numHits, final int numApproxHits){
		callsToExtendScore++;
		hist_hits_extend[Tools.min(HIT_HIST_LEN, numApproxHits)]++;
		
		final int centerVal=values[centerIndex];
		final int centerLoc=numberToSite(centerVal);
		
		final int minLoc=Tools.max(0, centerLoc-MAX_INDEL); //Legacy, for assertions
		final int maxLoc=centerLoc+MAX_INDEL2; //Legacy, for assertions

		final int minVal=centerVal-MAX_INDEL;
		final int maxVal=centerVal+MAX_INDEL2;
		
		final byte[] ref=Data.getChromosome(chrom).array;
		
		if(verbose){
			System.err.println("\n");
			System.err.println("minLoc="+minLoc+", maxLoc="+ maxLoc+", centerIndex="+centerIndex+", centerVal="+centerVal+", centerLoc="+centerLoc);
			System.err.println("minVal="+minVal+", maxVal="+ maxVal+", numHits="+numHits+", numApproxHits="+numApproxHits);
			System.err.println("offsets:\t"+Arrays.toString(offsets));
			System.err.println("values:\t"+Arrays.toString(values));
			System.err.println();
			int centerOffset=offsets[centerIndex];

			for(int i=0; i<centerOffset; i++){System.err.print(" ");}
			for(int i=centerOffset; i<centerOffset+KEYLEN; i++){System.err.print((char)bases[i]);}
			System.err.println();

			System.err.println(new String(bases));
			System.err.println(new String(KillSwitch.copyOfRange(ref, centerLoc, centerLoc+bases.length)));
			System.err.println();
		}
		
//		int[] locArray=new int[bases.length];
		Arrays.fill(locArray, -1);
		

//		if(verbose){
//			System.err.println("Reverse fill:");
//		}
		
		//First fill in reverse
		for(int i=0, keynum=0; i<numHits; i++){
			final int value=values[i];
//			if(verbose){System.err.println("values["+i+"]="+value);}
			
			if(value>=minVal && value<=maxVal){
				final int refbase=numberToSite(value);
//				if(verbose){System.err.println("refbase="+refbase);}
				
//				Exception in thread "Thread-23" java.lang.AssertionError: 71543, 536356470, 536956470
//		        at align2.BBIndex.extendScore(BBIndex.java:2620)
//		        at align2.BBIndex.slowWalk3(BBIndex.java:1393)
//		        at align2.BBIndex.find(BBIndex.java:777)
//		        at align2.BBIndex.find(BBIndex.java:623)
//		        at align2.BBIndex.findAdvanced(BBIndex.java:400)
//		        at align2.AbstractMapThread.quickMap(AbstractMapThread.java:761)
//		        at align2.BBMapThread.processRead(BBMapThread.java:408)
//		        at align2.AbstractMapThread.run(AbstractMapThread.java:519)
				
				//TODO - figure out why this assertion fires
//				assert(refbase>=minLoc && refbase<=maxLoc) : refbase+", "+minLoc+", "+maxLoc; //Apparently not a correct assumption

				//			System.err.println("Reverse: Trying key "+refbase+" @ "+offsets[i]);
				//				System.err.println("Passed!");
				keynum++;
				final int callbase=offsets[i];

				int misses=0;
				for(int cloc=callbase+KEYLEN-1, rloc=refbase+cloc; cloc>=0 && rloc>=0 && rloc<ref.length; cloc--, rloc--){
					int old=locArray[cloc];
					if(old==refbase){
//						if(verbose){System.err.println("Broke because old="+old+", refbase="+refbase);}
						break;
					} //Already filled with present value
					if(misses>0 && old>=0){
//						if(verbose){System.err.println("Broke because old="+old+", misses="+misses);}
						break;
					} //Already filled with something that has no errors
					byte c=bases[cloc];
					byte r=ref[rloc];

					if(c==r){
						if(old<0 || refbase==centerLoc){ //If the cell is empty or this key corresponds to center
							locArray[cloc]=refbase;
						}
					}else{
						misses++;
						//Only extends first key all the way back.  Others stop at the first error.
						if(old>=0 || keynum>1){
//							if(verbose){System.err.println("Broke because old="+old+", keynum="+keynum);}
							break;
						}
					}
				}
			}
		}
		
//		if(verbose){
//			System.err.println("locArray:\t"+Arrays.toString(locArray));
//			System.err.println("Forward fill:");
//		}
		
		//Then fill forward
		for(int i=0; i<numHits; i++){
			final int value=values[i];
//			if(verbose){System.err.println("values["+i+"]="+value);}
			
			if(value>=minVal && value<=maxVal){
				final int refbase=numberToSite(value);
//				if(verbose){System.err.println("refbase="+refbase);}
				//TODO - figure out why this assertion fires
//				assert(refbase>=minLoc && refbase<=maxLoc) : refbase+", "+minLoc+", "+maxLoc; //Apparently not a correct assumption
				final int callbase=offsets[i];
				
				int misses=0;
				for(int cloc=callbase+KEYLEN, rloc=refbase+cloc; cloc<bases.length && rloc<ref.length; cloc++, rloc++){
					int old=locArray[cloc];
					if(old==refbase){break;} //Already filled with present value
					if(misses>0 && old>=0){break;} //Already filled with something that has no errors
					byte c=bases[cloc];
					byte r=ref[rloc];
					
					if(c==r){
						if(old<0 || refbase==centerLoc){ //If the cell is empty or this key corresponds to center
							locArray[cloc]=refbase;
						}
					}else{
						misses++;
						if(old>=0){break;} //Already filled with something that has no errors
					}
				}
			}
		}
		
		//Try to subsume out-of-order locs where higher numbers come before lower numbers. Made things worse.
//		for(int i=1; i<locArray.length; i++){
//			final int loc=locArray[i];
//			final int last=locArray[i-1];
//			if(loc<last && last>-1){
//				final byte c=bases[i];
////				final int rloc1=loc+i;
//				final int rloc2=last+i;
//				byte r2=ref[rloc2];
//				if(r2==c){
//					locArray[i]=last;
//				}
//			}
//		}
		
//		//Change 'N' to -2.  A bit slow.
//		{
//			int firstMatch=0;
//			while(firstMatch<locArray.length && locArray[firstMatch]<0){firstMatch++;}
//			assert(firstMatch<locArray.length) : new String(bases);
//			int last=locArray[firstMatch];
//			for(int i=firstMatch-1; i>=0; i--){
//				final byte c=bases[i];
//				if(c=='N'){locArray[i]=-2;}
//				else{
//					assert(locArray[i]==-1);
//					final int rloc=last+i;
//					byte r=ref[rloc];
//					if(r=='N'){locArray[i]=-2;}
//				}
//			}
//			for(int i=firstMatch; i<locArray.length; i++){
//				final int loc=locArray[i];
//				if(last<1){last=loc;}
//				final byte c=bases[i];
//				if(c=='N'){locArray[i]=-2;}
//				else if(loc==-1 && last>0){
//					final int rloc=last+i;
//					byte r=ref[rloc];
//					if(r=='N'){locArray[i]=-2;}
//				}
//			}
//		}
		
		//Change 'N' to -2, but only for nocalls, not norefs.  Much faster.
		{
			final byte nb=(byte)'N';
			for(int i=0; i<bases.length; i++){
				if(bases[i]==nb){locArray[i]=-2;}
			}
		}
		
//
//		{
//			int last=locArray[0];
//			for(int i=1; i<locArray.length; i++){
//				final int loc=locArray[i];
//				if(loc>0){
//					if(last<1){last=loc;}
//				}else{
//					if(last>0){
//						final byte c=bases[i];
//						final int rloc2=last+i;
//						byte r2=ref[rloc2];
//						if(c=='N'){
//
//						}else if(c==r2){
//							locArray[i]=last;
//						}
//					}
//					last=-1;
//				}
//			}
//		}
//		{
//			int last=locArray[locArray.length-1];
//			for(int i=locArray.length-2; i>=0; i--){
//				final int loc=locArray[i];
//				if(loc>0){
//					if(last<1){last=loc;}
//				}else{
//					if(last>0){
//						final byte c=bases[i];
//						final int rloc2=last+i;
//						byte r2=ref[rloc2];
//						if(c=='N'){
//
//						}else if(c==r2){
//							locArray[i]=last;
//						}
//					}
//					last=-1;
//				}
//			}
//		}
		
//		for(int i=locArray.length-2; i>=0; i--){
//			final int loc=locArray[i];
//			final int last=locArray[i+1];
//			if(loc>last && last>-1){
//				final byte c=bases[i];
////				final int rloc1=loc+i;
//				final int rloc2=last+i;
//				byte r2=ref[rloc2];
//				if(r2==c){
//					locArray[i]=last;
//				}
//			}
//		}
		
		if(verbose){
//			System.err.println("locArray:\t"+Arrays.toString(locArray));
			
			int centerOffset=offsets[centerIndex];
			int lim=centerOffset+KEYLEN;
			for(int i=centerOffset; i<lim; i++){
				assert(locArray[i]!=-1) : " ( "+centerOffset+" < "+i+" < "+lim+" ) "+
						"\nlocArray: "+Arrays.toString(locArray)+
						"\nvalues: "+Arrays.toString(values)+
						"\noffsets: "+Arrays.toString(offsets);
			}
		}
		
		if(USE_AFFINE_SCORE){
			/* TODO - sometimes returns a higher score than actual alignment.  This should never happen. */
			int score=(KFILTER<2 ? msa.calcAffineScore(locArray, baseScores, bases) :
				msa.calcAffineScore(locArray, baseScores, bases, KFILTER));
			return score;
		}
		
		int score=0;
		int lastLoc=-1;
		int centerBonus=BASE_HIT_SCORE/5;
		for(int i=0; i<locArray.length; i++){
			int loc=locArray[i];
			if(loc>=0){
				score+=BASE_HIT_SCORE+baseScores[i];
				if(loc==centerLoc){score+=centerBonus;}
				if(loc!=lastLoc && lastLoc>=0){
					int dif=absdif(loc, lastLoc);
					int penalty=Tools.min(INDEL_PENALTY+INDEL_PENALTY_MULT*dif, MAX_PENALTY_FOR_MISALIGNED_HIT);
					score-=penalty;
				}
				lastLoc=loc;
			}
		}

//		System.err.println("Extended score: "+score);
//		System.err.println(Arrays.toString(locArray));
		
		
		return score;
	}
	
	
	/** NOTE!  This destroys the locArray, so use a copy if needed. */
	private static final int[] makeGapArray(int[] locArray, int minLoc, int minGap){
		int gaps=0;
		boolean doSort=false;
		
		if(locArray[0]<0){locArray[0]=minLoc;}
		for(int i=1; i<locArray.length; i++){
			if(locArray[i]<0){locArray[i]=locArray[i-1]+1;}
			else{locArray[i]+=i;}
			if(locArray[i]<locArray[i-1]){doSort=true;}
		}
		
//		System.err.println(Arrays.toString(locArray)+"\n");
		
		if(doSort){
//			System.err.println("*");
			Arrays.sort(locArray);
		}
//		System.err.println(Arrays.toString(locArray)+"\n");
		
		for(int i=1; i<locArray.length; i++){
			int dif=locArray[i]-locArray[i-1];
			assert(dif>=0);
			if(dif>minGap){
				gaps++;
			}
		}
		if(gaps<1){return null;}
		int[] out=new int[2+gaps*2];
		out[0]=locArray[0];
		out[out.length-1]=locArray[locArray.length-1];
		
		for(int i=1, j=1; i<locArray.length; i++){
			int dif=locArray[i]-locArray[i-1];
			assert(dif>=0);
			if(dif>minGap){
				out[j]=locArray[i-1];
				out[j+1]=locArray[i];
				j+=2;
			}
		}
		return out;
	}
	
	
	/** Generates a term that increases score with how many bases in the read match the ref. */
	private final int scoreZ2(int[] locs, int centerIndex, int offsets[], int numApproxHits, int numHits){
		
		if(numApproxHits==1){return SCOREZ_1KEY;}
		
		final int center=locs[centerIndex];

		final int maxLoc=center+MAX_INDEL2;
		final int minLoc=Tools.max(0, center-MAX_INDEL);
		
		int score=0;
		
		int a0=-1, b0=-1;
		
		for(int i=0; i<numHits; i++){
			int loc=locs[i];
//			int dif=absdif(loc, center);
			if(loc>=minLoc && loc<=maxLoc){
//				assert(loc>=center) : "loc="+loc+"\ni="+i+"\ncenterIndex="+centerIndex+
//					"\nmaxLoc="+maxLoc+"\nlocs:\t"+Arrays.toString(locs)+"\noffsets:\t"+Arrays.toString(offsets);
				int a=offsets[i];
				
				if(b0<a){
					score+=b0-a0;
					a0=a;
				}
				b0=a+KEYLEN;
			}
		}
		score+=b0-a0;
		score=score*Z_SCORE_MULT;
//		assert(score==scoreZslow(locs, centerIndex, offsets, false)) : scoreZslow(locs, centerIndex, offsets, true)+" != "+score;
		return score;
	}
	
	@Deprecated
	/** This was just to verify scoreZ2. */
	private final int scoreZslow(int[] locs, int centerIndex, int offsets[], boolean display, int numHits){
		final int center=locs[centerIndex];

		final int maxLoc=center+MAX_INDEL2;
		final int minLoc=Tools.max(0, center-MAX_INDEL);
		
		byte[] array=new byte[offsets[offsets.length-1]+KEYLEN];
		int score=0;
		
		for(int i=0; i<numHits; i++){
			int loc=locs[i];
//			int dif=absdif(loc, center);
			if(loc>=minLoc && loc<=maxLoc){
				int pos=offsets[i];
//				if(true){
//					System.err.println("\ni="+i+", pos="+pos+", array=["+array.length+"], limit="+(pos+CHUNKSIZE-1));
//				}
				for(int j=pos; j<pos+KEYLEN; j++){
					if(array[j]==0){score++;}
					array[j]=1;
				}
			}
		}
		
		if(display){System.err.println("\n"+Arrays.toString(array)+"\n");}
		
		return score*Z_SCORE_MULT;
	}
	
	/** Generates a term that increases score with how many bases in the read match the ref. */
	private final int maxScoreZ(int offsets[]){
		int score=0;
		int a0=-1, b0=-1;

		for(int i=0; i<offsets.length; i++){
			int a=offsets[i];

			if(b0<a){
				score+=b0-a0;
				a0=a;
			}
			b0=a+KEYLEN;

		}
		score+=b0-a0;
		return score*Z_SCORE_MULT;
	}
	

	private final int scoreRight(int[] locs, int[] keyScores, int centerIndex, int[] sizes, boolean penalizeIndels, int numHits){
		
		int score=0;
		
		int prev, loc=locs[centerIndex];
		
		for(int i=centerIndex+1; i<numHits; i++){
			
			if(locs[i]>=0){
				prev=loc;
				loc=locs[i];
				
				int offset=absdif(loc, prev);
				
				if(offset<=MAX_INDEL){
					score+=keyScores[i];
					if(ADD_LIST_SIZE_BONUS){score+=calcListSizeBonus(sizes[i]);}
					
//					if(i==locs.length-1){score+=HIT_SCORE/2;} //Adds a bonus for matching the first or last key
					
					if(penalizeIndels && offset!=0){
						int penalty=Tools.min(INDEL_PENALTY+INDEL_PENALTY_MULT*offset, MAX_PENALTY_FOR_MISALIGNED_HIT);
						score-=penalty;
//						score-=(INDEL_PENALTY+Tools.min(INDEL_PENALTY_MULT*offset, 1+HIT_SCORE/4));
					}
				}else{
					loc=prev;
				}
			}
			
		}
		return score;
		
	}
	
	private final int scoreLeft(int[] locs, int[] keyScores, int centerIndex, int[] sizes, boolean penalizeIndels){
		
		callsToScore++;
		
		int score=0;
		
		int prev, loc=locs[centerIndex];
		
		for(int i=centerIndex-1; i>=0; i--){
			
			if(locs[i]>=0){
				prev=loc;
				loc=locs[i];
				
				int offset=absdif(loc, prev);
				
				if(offset<=MAX_INDEL){
					score+=keyScores[i];
					if(ADD_LIST_SIZE_BONUS){score+=calcListSizeBonus(sizes[i]);}
					
//					if(i==0){score+=HIT_SCORE/2;} //Adds a bonus for matching the first or last key
					if(penalizeIndels && offset!=0){
						int penalty=Tools.min(INDEL_PENALTY+INDEL_PENALTY_MULT*offset, MAX_PENALTY_FOR_MISALIGNED_HIT);
						score-=penalty;
					}
				}else{
					loc=prev;
				}
			}
			
		}
		return score;
		
	}
	
	/** Encode a (location, chrom) pair to an index */
	private static final int toNumber(int site, int chrom){
		int out=(chrom&CHROM_MASK_LOW);
		out=out<<SHIFT_LENGTH;
		out=(out|site);
		return out;
	}
	
	/** Decode an (index, baseChrom) pair to a chromosome */
	private static final int numberToChrom(int number, int baseChrom){
		assert((baseChrom&CHROM_MASK_LOW)==0) : Integer.toHexString(number)+", baseChrom="+baseChrom;
		assert(baseChrom>=0) : Integer.toHexString(number)+", baseChrom="+baseChrom;
		int out=(number>>>SHIFT_LENGTH);
		out=out+(baseChrom&CHROM_MASK_HIGH);
		return out;
	}
	
	/** Decode an index to a location */
	private static final int numberToSite(int number){
		return (number&SITE_MASK);
	}

	private static final int minChrom(int chrom){return Tools.max(MINCHROM, chrom&CHROM_MASK_HIGH);}
	private static final int baseChrom(int chrom){return Tools.max(0, chrom&CHROM_MASK_HIGH);}
	private static final int maxChrom(int chrom){return Tools.max(MINCHROM, Tools.min(MAXCHROM, chrom|CHROM_MASK_LOW));}
	
	
	private final int[] getOffsetArray(int len){
		if(offsetArrays[len]==null){offsetArrays[len]=new int[len];}
		return offsetArrays[len];
	}
	private final int[] getLocArray(int len){
		if(len>=locArrays.length){return new int[len];}
		if(locArrays[len]==null){locArrays[len]=new int[len];}
		return locArrays[len];
	}
	private final int[] getGreedyListArray(int len){
		if(greedyListArrays[len]==null){greedyListArrays[len]=new int[len];}
		return greedyListArrays[len];
	}
	private final int[] getGenericArray(int len){
		if(genericArrays[len]==null){genericArrays[len]=new int[len];}
		return genericArrays[len];
	}

	@Override
	final byte[] getBaseScoreArray(int len, int strand){
		if(len>=baseScoreArrays[0].length){return new byte[len];}
		if(baseScoreArrays[strand][len]==null){baseScoreArrays[strand][len]=new byte[len];}
		return baseScoreArrays[strand][len];
	}
	@Override
	final int[] getKeyScoreArray(int len, int strand){
		if(len>=keyScoreArrays.length){return new int[len];}
		if(keyScoreArrays[strand][len]==null){keyScoreArrays[strand][len]=new int[len];}
		return keyScoreArrays[strand][len];
	}
	private final float[] getKeyWeightArray(int len){
		if(len>=keyWeightArrays.length){return new float[len];}
		if(keyWeightArrays[len]==null){keyWeightArrays[len]=new float[len];}
		return keyWeightArrays[len];
	}
	@Override
	float[] keyProbArray() {
		return keyProbArray;
	}
	
	private static final int KEY_BUFFER_LENGTH=256; //Boosted from 128 to 256 to deal with crashes from short kmers, long reads, and vslow.
	
	private final int[][] locArrays=new int[601][];
	private final int[] valueArray=new int[KEY_BUFFER_LENGTH];
	private final int[] sizeArray=new int[KEY_BUFFER_LENGTH];
	private final int[][] offsetArrays=new int[KEY_BUFFER_LENGTH][];
	private final int[][] greedyListArrays=new int[KEY_BUFFER_LENGTH][];
	private final int[][] genericArrays=new int[KEY_BUFFER_LENGTH][];
	private final int[] startArray=new int[KEY_BUFFER_LENGTH];
	private final int[] stopArray=new int[KEY_BUFFER_LENGTH];
	private final Quad[] tripleStorage=makeQuadStorage(KEY_BUFFER_LENGTH);
	private final int[] greedyReturn=KillSwitch.allocInt1D(2);
	private final int[][] shrinkReturn2=new int[3][];
	private final int[][] shrinkReturn3=new int[5][];
	private final int[][] prescanReturn=new int[2][];
	private final int[] prescoreArray;
	private final int[] precountArray;

	private final byte[][][] baseScoreArrays=new byte[2][601][];
	private final int[][][] keyScoreArrays=new int[2][KEY_BUFFER_LENGTH][];
	final float[] keyProbArray=new float[601];
	private final float[][] keyWeightArrays=new float[KEY_BUFFER_LENGTH][];
	
	
	private final static Quad[] makeQuadStorage(int number){
		Quad[] r=new Quad[number];
		for(int i=0; i<number; i++){r[i]=new Quad(i, 0, 0);}
		return r;
	}
	

	private final QuadHeap heap=new QuadHeap(KEY_BUFFER_LENGTH-1);
	private final QuadHeap active=new QuadHeap(KEY_BUFFER_LENGTH-1);
	
	static int SHIFT_LENGTH=(32-1-NUM_CHROM_BITS);
	static int MAX_ALLOWED_CHROM_INDEX=~((-1)<<SHIFT_LENGTH);
	
	/** Mask the number to get the site, which is in the lower bits */
	static int SITE_MASK=((-1)>>>(NUM_CHROM_BITS+1));
	
	/** Mask the chromosome's high bits to get the low bits */
	static int CHROM_MASK_LOW=CHROMS_PER_BLOCK-1;
	
	/** Mask the chromosome's lower bits to get the high bits */
	static int CHROM_MASK_HIGH=~CHROM_MASK_LOW;
	
	static void setChromBits(int x){
		
		NUM_CHROM_BITS=x;
		CHROMS_PER_BLOCK=(1<<(NUM_CHROM_BITS));
		SHIFT_LENGTH=(32-1-NUM_CHROM_BITS);
		MAX_ALLOWED_CHROM_INDEX=~((-1)<<SHIFT_LENGTH);
		SITE_MASK=((-1)>>>(NUM_CHROM_BITS+1));
		CHROM_MASK_LOW=CHROMS_PER_BLOCK-1;
		CHROM_MASK_HIGH=~CHROM_MASK_LOW;
		
//		assert(NUM_CHROM_BITS<30);
		assert(NUM_CHROM_BITS>=0); //max is 3 for human; perhaps more for other organisms
//		assert((1<<(NUM_CHROM_BITS))>=CHROMSPERBLOCK) : (1<<(NUM_CHROM_BITS))+" < "+CHROMSPERBLOCK;
		assert((1<<(NUM_CHROM_BITS))==CHROMS_PER_BLOCK) : (1<<(NUM_CHROM_BITS))+" < "+CHROMS_PER_BLOCK;
		assert(Integer.bitCount(CHROMS_PER_BLOCK)==1);
		assert(Integer.numberOfLeadingZeros(SITE_MASK)==(NUM_CHROM_BITS+1)) : Integer.toHexString(SITE_MASK);
	}
	
	private final int cycles;

	public static final int BASE_HIT_SCORE=100;
	public static final int ALIGN_COLUMNS=3000;
	public static int MAX_INDEL=16000; //Max indel length, min 0, default 400; longer is more accurate
	public static int MAX_INDEL2=2*MAX_INDEL;
	
	private final float INV_BASE_KEY_HIT_SCORE;
	private final int INDEL_PENALTY; //default (HIT_SCORE/2)-1
	private final int INDEL_PENALTY_MULT; //default 20; penalty for indel length
	private final int MAX_PENALTY_FOR_MISALIGNED_HIT;
	private final int SCOREZ_1KEY;

	public static final boolean ADD_SCORE_Z=true; //Increases quality, decreases speed
	public static final int Z_SCORE_MULT=20;
	public static final int Y_SCORE_MULT=10;
	
	
	/**
	 * Return only sites that match completely or with partial no-reference
	 */
	public static void setSemiperfectMode() {
		assert(!PERFECTMODE);
		SEMIPERFECTMODE=true;
		PRESCAN_QSCORE=false;
//		MIN_APPROX_HITS_TO_KEEP++;
		
		
		
		MAX_INDEL=0;
		MAX_INDEL2=0;
	}
	
	/**
	 * Return only sites that match completely
	 */
	public static void setPerfectMode() {
		assert(!SEMIPERFECTMODE);
		PERFECTMODE=true;
		PRESCAN_QSCORE=false;
//		MIN_APPROX_HITS_TO_KEEP++;
		
		
		
		MAX_INDEL=0;
		MAX_INDEL2=0;
	}
	
	static float FRACTION_GENOME_TO_EXCLUDE=0.03f; //Default .03; lower is slower and more accurate.  For perfect reads and small genomes, lower is FASTER.
	
	public static final void setFractionToExclude(float f){
		assert(f>=0 && f<1);
		FRACTION_GENOME_TO_EXCLUDE=f;
		MIN_INDEX_TO_DROP_LONG_HIT_LIST=(int)(1000*(1-3.5*FRACTION_GENOME_TO_EXCLUDE)); //default 810
		MAX_AVERAGE_LIST_TO_SEARCH=(int)(1000*(1-2.3*FRACTION_GENOME_TO_EXCLUDE)); //lower is faster, default 840
		MAX_AVERAGE_LIST_TO_SEARCH2=(int)(1000*(1-1.4*FRACTION_GENOME_TO_EXCLUDE)); //default 910
		MAX_SINGLE_LIST_TO_SEARCH=(int)(1000*(1-1.0*FRACTION_GENOME_TO_EXCLUDE)); //default 935
		MAX_SHORTEST_LIST_TO_SEARCH=(int)(1000*(1-2.8*FRACTION_GENOME_TO_EXCLUDE)); //Default 860
	}
	
	
	/** Default .75.  Range: 0 to 1 (but 0 will break everything).  Lower is faster and less accurate. */
	static final float HIT_FRACTION_TO_RETAIN=0.85f; //default: .8
	/** Range: 0 to 1000.  Lower should be faster and less accurate. */
	static int MIN_INDEX_TO_DROP_LONG_HIT_LIST=(int)(1000*(1-3.5*FRACTION_GENOME_TO_EXCLUDE)); //default 810
	/** Range: 2 to infinity.  Lower should be faster and less accurate. */
	static final int MIN_HIT_LISTS_TO_RETAIN=6;
	
	static int MAX_AVERAGE_LIST_TO_SEARCH=(int)(1000*(1-2.3*FRACTION_GENOME_TO_EXCLUDE)); //lower is faster, default 840
	//lower is faster
	static int MAX_AVERAGE_LIST_TO_SEARCH2=(int)(1000*(1-1.4*FRACTION_GENOME_TO_EXCLUDE)); //default 910
	//lower is faster
	static int MAX_SINGLE_LIST_TO_SEARCH=(int)(1000*(1-1.0*FRACTION_GENOME_TO_EXCLUDE)); //default 935
	//lower is faster
	static int MAX_SHORTEST_LIST_TO_SEARCH=(int)(1000*(1-2.8*FRACTION_GENOME_TO_EXCLUDE)); //Default 860
	
	/** To increase accuracy on small genomes, override greedy list dismissal when the list is at most this long. */
	public static final int SMALL_GENOME_LIST=20;
	
	static{assert(!(TRIM_BY_GREEDY && TRIM_BY_TOTAL_SITE_COUNT)) : "Pick one.";}
	
	static final int CLUMPY_MAX_DIST=5; //Keys repeating over intervals of this or less are clumpy.
	
	/** Minimum length of list before clumpiness is considered. This is an index in the length histogram, from 0 to 1000. */
	static final int CLUMPY_MIN_LENGTH_INDEX=2000;
	static final float CLUMPY_FRACTION=0.75f; //0 to 1; higher is slower but more stringent. 0.5 means the median distance is clumpy.
	
	static final int MAX_SUBSUMPTION_LENGTH=MAX_INDEL2;
	
	/** approxHitsCutoff=maxHits-MAX_HITS_REDUCTION1 when slowWalk3 is first entered */
	public static final int MAX_HITS_REDUCTION1=0;
	
	/** approxHitsCutoff=maxHits-MAX_HITS_REDUCTION2 dynamically when best score is exceeded */
	public static int MAX_HITS_REDUCTION2=2; //default 1; higher is more accurate (more mapping and less FP) but slower
	
	/** approxHitsCutoff=maxHits-MAX_HITS_REDUCTION_PERFECT when perfect score is found */
	public static final int MAX_HITS_REDUCTION_PERFECT=0;

	public static int MAXIMUM_MAX_HITS_REDUCTION=3;
	public static int HIT_REDUCTION_DIV=5;
	
	private static final int calcApproxHitsCutoff(final int keys, final int hits, int currentCutoff, final boolean perfect){ //***$
		assert(keys>=hits) : keys+", "+hits;
		assert(hits>=0);
		
		int mahtk=MIN_APPROX_HITS_TO_KEEP;
		if(SEMIPERFECTMODE || PERFECTMODE){
			if(keys==1){return 1;}
			else if(MIN_APPROX_HITS_TO_KEEP<keys){
				mahtk++;
				if(currentCutoff==MIN_APPROX_HITS_TO_KEEP){currentCutoff++;}
			}
		}
		
		int reduction=Tools.min(Tools.max((hits)/HIT_REDUCTION_DIV, MAX_HITS_REDUCTION2), Tools.max(MAXIMUM_MAX_HITS_REDUCTION, keys/8));
//		if(hits<3){reduction=0;}
//		else if(hits<4){reduction=Tools.min(1, reduction);}
//		else if(hits<5){reduction=Tools.min(2, reduction);}
//		else if(hits<6){reduction=Tools.min(3, reduction);}
		assert(reduction>=0);
		int r=hits-reduction;
		
		r=Tools.max(mahtk, currentCutoff, r);
		
		if(perfect){
			r=Tools.max(r, keys-MAX_HITS_REDUCTION_PERFECT);
		}
		return r;
	}
	
	public static final boolean USE_SLOWALK3=true && USE_EXTENDED_SCORE;
	public static boolean PRESCAN_QSCORE=true && USE_EXTENDED_SCORE; //Decrease quality and increase speed
	public static final boolean FILTER_BY_QSCORE=true; //Slightly lower quality, but very fast.
	public static final float MIN_SCORE_MULT=(USE_AFFINE_SCORE ? 0.15f : USE_EXTENDED_SCORE ? .3f : 0.10f);  //Fraction of max score to use as cutoff.  Default 0.15, max is 1; lower is more accurate
	public static final float MIN_QSCORE_MULT=0.025f;  //Fraction of max score to use as cutoff.  Default 0.15, max is 1; lower is more accurate
	public static final float MIN_QSCORE_MULT2=0.1f;
	static final float DYNAMIC_SCORE_THRESH=(USE_AFFINE_SCORE ? 0.84f : USE_EXTENDED_SCORE ? .84f : 0.6f); //Default .85f; lower is more accurate
	static final float DYNAMIC_QSCORE_THRESH=0.6f; //default .58f
	static final float DYNAMIC_QSCORE_THRESH_PERFECT=0.8f; //***$
	static final float PRESCAN_QSCORE_THRESH=DYNAMIC_QSCORE_THRESH*.95f; //default 1.0f; lower is more accurate and 0 essentially sets PRESCAN_QSCORE=false
	static{
		assert(MIN_SCORE_MULT>=0 && MIN_SCORE_MULT<1);
		assert(DYNAMIC_SCORE_THRESH>=0 && DYNAMIC_SCORE_THRESH<1);
	}
	
	
}
