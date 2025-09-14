package clump;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import dna.AminoAcid;
import shared.KillSwitch;
import shared.Tools;
import stream.Read;
import structures.IntList;
import structures.LongList;

/**
 * A tool for splitting clumps by allele.
 * @author Brian Bushnell
 * @date September 26, 2016
 *
 */
class Splitter {
	
	static ArrayList<Clump> splitOnPivot(Clump c){
		ArrayList<Clump> list=new ArrayList<Clump>(3);
		list.add(c);
		if(c.size()<minSizeSplit){
//			assert(findBestPivot(c)<0);
			return list;
		}
		return splitOnPivot(list);
	}
	
	static ArrayList<Clump> splitOnPivot(ArrayList<Clump> list){
		ArrayList<Clump> out=new ArrayList<Clump>();
		
		final IntList pivots=new IntList(2);
		for(int i=0; i<list.size(); i++){
			Clump clump=list.get(i);
			list.set(i, null);
			int pivot=findBestPivots(clump, FIND_CORRELATIONS, pivots);
			if(pivot<0){
				assert(pivots.size==0);
				out.add(clump);
			}else{
				assert(pivots.size==1 || pivots.size==2) : pivot+", "+pivots.size+", "+pivots;
				assert(pivots.get(0)==pivot);
				int added=splitAndAdd(clump, pivots.get(0), (pivots.size>1 ? pivots.get(1) : -1), list);
				if(added<2) {//Rare case that caused an assertion once
					assert(added==1);
					out.add(clump);
					i++;
				}
			}
		}
		return out;
	}
	
	static int splitAndAdd(Clump c, final int var1, final int var2, ArrayList<Clump> list) {
		final int maxLeft=c.maxLeft();
		
		final Clump major=Clump.makeClump(c.kmer), minor=Clump.makeClump(c.kmer);
		
		for(Read r : c){
			if(containsVar(var1, r, maxLeft) || (var2>=0 && containsVar(var2, r, maxLeft))){
				minor.add(r);
			}else{
				major.add(r);
			}
		}
		
//		assert(major.size()>0) : c.size()+", "+major.size()+", "+minor.size()+
//			", "+(var1>>2)+", "+(var1&alleleMask)+", "+(var2>>2)+", "+(var2&alleleMask)+
//			", "+list.size()+", "+c.get(0);
//		assert(minor.size()>0) : c.size()+", "+major.size()+", "+minor.size()+
//			", "+(var1>>2)+", "+(var1&alleleMask)+", "+(var2>>2)+", "+(var2&alleleMask)+
//			", "+list.size()+", "+c.get(0);
//		assert(major.size()+minor.size()==c.size()) : c.size()+", "+major.size()+", "+minor.size()+
//			", "+(var1>>2)+", "+(var1&alleleMask)+", "+(var2>>2)+", "+(var2&alleleMask)+
//			", "+list.size()+", "+c.get(0);
		
//		Exception in thread "Thread-110" java.lang.AssertionError: 21, 0, 21, 163, 1, 104, 2, 71,      
//		TTTTGGCAGCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCGCCACGCCGCCCCCCCCCCCC <-90
//		CCCCCCCCCCGGGGGGGGGGGGCAGGGGGGGGGGGGGGGGGGGGGGGCGCCCCCCCCCC
//		        at clump.Splitter.splitAndAdd(Splitter.java:72)
//		        at clump.Splitter.splitOnPivot(Splitter.java:46)
//		        at clump.Splitter.splitOnPivot(Splitter.java:29)
//		        at clump.Clump.splitAndErrorCorrect(Clump.java:612)
//		        at clump.ClumpList$ProcessThread.run(ClumpList.java:342)
		
		
		int added=(major.size()>0 ? 1 : 0)+(minor.size()>0 ? 1 : 0);
		if(major.size()>0) {list.add(major);}
		if(minor.size()>0) {list.add(minor);}
		assert(added>0) : c.size()+", "+major.size()+", "+minor.size();
		return added;
	}
	
	//Returns the c
	static int countVariants(Clump c, LongList varCounts){
		varCounts.clear();
		final int[][] bcounts=c.baseCounts();
		final int[][] qcounts=c.qualityCounts();
		final int len=bcounts[0].length;
		for(int i=0; i<len; i++){
			final int major=c.getConsensusAtPosition(qcounts, i);
			for(int x=0; x<4; x++){
				final int bcount=bcounts[x][i];
				if(bcount>1 && x!=major){
					final long var=(((long)bcount)<<32)|((i<<shift)|x);
					varCounts.add(var);
				}
			}
		}
		if(varCounts.size()<1){return 0;}
		varCounts.sort();
		varCounts.reverse();
		return (int)(varCounts.get(0)>>>32);
	}
	
	static LinkedHashMap<Integer, ArrayList<Read>> findReadVariants(Clump c, boolean makeMap){
		if(c.size()<5){return null;} //Not really needed with tiny clumps
//		assert(c.size()>3); //TODO
		LinkedHashMap<Integer, ArrayList<Read>> map=null;
		map=(makeMap ? new LinkedHashMap<Integer, ArrayList<Read>>() : null);
//		if(makeMap){
//			map=localMap.get();
//			if(map==null){
//				map=new LinkedHashMap<Integer, ArrayList<Read>>();
//				localMap.set(map);
//			}
//			map.clear();
//		}
		
		final int[][] bcounts=c.baseCounts();
		final Read ref=c.consensusRead();
		final byte[] rbases=ref.bases;
		for(Read r : c){
			final byte[] bases=r.bases;
			final ReadKey key=(ReadKey) r.obj;
			IntList list=key.vars;
			if(list!=null){list.clear();}
			
			
			final int cStart=0, rStart=c.maxLeft()-key.position;
			
			for(int i=cStart, j=rStart; i<bases.length; i++, j++){
				final byte cb=bases[i], rb=rbases[j];
				if(cb!=rb){
					byte x=AminoAcid.baseToNumber[cb];
					if(x>=0){
						int count=bcounts[x][j];
						if(count>1){
							int var=(j<<2)|x;
							if(list==null){list=key.vars=new IntList(4);}
							list.add(var);
							
							if(map!=null){
								Integer mapkey=var;//mapKeys[var];
								ArrayList<Read> alr=map.get(mapkey);
								if(alr==null){
									alr=new ArrayList<Read>(4);
									map.put(mapkey, alr);
								}
								alr.add(r);
							}
							
						}
					}
				}
			}
		}
		return map==null || map.isEmpty() ? null : map;
	}
	
	static int findBestPivot_Correlated(Clump c, IntList pivots){
		assert(pivots.size==0);
		LinkedHashMap<Integer, ArrayList<Read>> map=findReadVariants(c, true);
		if(map==null){return -1;}
		
		IntList collection=new IntList(32);
		int[] rvector=KillSwitch.allocInt1D(5);
		
		int bestVar=-1;
		int bestVarCount=-1;
		int bestVar2=-1;
		int bestVar2Count=-1;
		int bestDifferent=-1;
		float bestCorrelation=-1;
		float bestScore=-1;
		
		final float minCorrelation=0.75f;
		final int minVar2Count=2;
		
		int max=0;
		for(Entry<Integer, ArrayList<Read>> entry : map.entrySet()){
			max=Tools.max(max, entry.getValue().size());
		}
		if(max<2){return -1;}
		final int thresh=Tools.max(2, max/2);
		final int thresh2=max/4;
		int numOverThresh=0;
		ArrayList<Integer> remove=new ArrayList<Integer>();
		for(Entry<Integer, ArrayList<Read>> entry : map.entrySet()){
			int x=entry.getValue().size();
			if(x>=thresh){numOverThresh++;}
			else if(x<thresh2){remove.add(entry.getKey());}
		}
		for(Integer key : remove){map.remove(key);}
		if(numOverThresh>MAX_CORRELATIONS){return -1;}
		
		for(Entry<Integer, ArrayList<Read>> entry : map.entrySet()){
			final ArrayList<Read> rlist=entry.getValue();
			final Integer key=entry.getKey();
			final int var=key;
			if(rlist.size()>=thresh){
				
				final int var2=examineVar(var, rlist, collection, rvector, map);

				if(var2>=0){

					final int varCount=rvector[1];
					final int var2Count=rvector[3];
					final int different=rvector[4];

					final int var2reads;
					final int var2ReadsWithoutVar;
					{
						ArrayList<Read> temp=map.get(var2);
						var2reads=(temp==null ? 0 : temp.size());

//						if(var2reads==var2Count){
//							var2ReadsWithoutVar=0;
//						}else{
//							var2ReadsWithoutVar=countDifferentAlleles(var, temp);
//						}

						var2ReadsWithoutVar=var2reads-varCount;

					}

//					final float correlation=(var2Count-0.05f)/(float)varCount;
//					final float correlation=(varCount-different)/(float)varCount;
					final float correlation=(Tools.max(varCount, var2reads)-different)/(float)Tools.max(varCount, var2reads);
					final int distance=Tools.absdif(var>>2, var2>>2);

//					final float score=correlation*((var2Count-1)-0.5f*var2ReadsWithoutVar-different+1.0f*(varCount));

					final float score=correlation*(var2reads/*var2Count*/+varCount-different+0.5f*var2ReadsWithoutVar)*(distance+250);

//					final float score=correlation*((var2Count-1)-1.0f*var2ReadsWithoutVar+1.0f*(varCount));
					if(correlation>=minCorrelation && var2Count>=minVar2Count){
						if(score>bestScore || (score==bestScore && varCount>bestVarCount)){
							bestVar=var;
							bestVarCount=varCount;
							bestVar2=var2;
							bestVar2Count=var2Count;
							bestCorrelation=correlation;
							bestScore=score;
							bestDifferent=different;
						}
					}
				}
			}
		}
		
		if(bestVar2Count>=minVar2Count && bestCorrelation>=minCorrelation){
			pivots.add(bestVar);
			pivots.add(bestVar2);
			return bestVar;
		}
		return -1;
	}
	
	static boolean containsVar(final int var, final Read r, final int maxLeft){
		final int varPos=var>>2;
		final int varAllele=var&alleleMask;
		final ReadKey rk=(ReadKey) r.obj;
		final int rloc=toReadLocation(varPos, maxLeft, rk.position);
		if(rloc<0 || rloc>=r.length()){
			return false;
		}
		final int readAllele=AminoAcid.baseToNumber[r.bases[rloc]];
		return readAllele==varAllele;
	}
	
	static boolean hasDifferentAllele(final int var, final Read r/*, final Clump c*/){
		final int varPos=var>>2;
		final int varAllele=var&alleleMask;
		final ReadKey rk=(ReadKey) r.obj;
		final IntList vars=rk.vars;
		final Clump c=rk.clump;
		assert(c==rk.clump);
		final int maxLeft=c.maxLeft();
		final int rloc=toReadLocation(varPos, maxLeft, rk.position);
		if(rloc<0 || rloc>=r.length()){
			assert(!vars.contains(var));
			return false;
		}
		final int readAllele=AminoAcid.baseToNumber[r.bases[rloc]];
		assert((readAllele==varAllele)==vars.contains(var));
		
		return readAllele!=varAllele;
	}
	
	static int countDifferentAlleles(final int var, ArrayList<Read> list){
		if(list==null){return 0;}
		int sum=0;
		for(Read r : list){
			if(hasDifferentAllele(var, r)){sum++;}
		}
		return sum;
	}
	
	static int examineVar(final int var, final ArrayList<Read> list, final IntList collection, final int[] rvector, LinkedHashMap<Integer, ArrayList<Read>> map){
		collection.clear();
		
		for(Read r : list){
			final ReadKey rk=(ReadKey) r.obj;
			final IntList vars=rk.vars;
			
			for(int i=0; i<vars.size; i++){
				final int v2=vars.get(i);
				if(v2!=var){
					collection.add(v2);
				}
			}
		}
		collection.sort();
		
		final int varCount=list.size();
		
		int lastVar2=-1, bestVar2=-1;
		int sharedCount=0, bestSharedCount=0, bestDifferent=999;
		for(int i=0; i<collection.size; i++){//TODO: Note that not all reads actually cover a given var
			int currentVar2=collection.get(i);
			if(currentVar2==lastVar2){sharedCount++;}
			else{
				if(sharedCount>bestSharedCount){
					final int different1=(sharedCount==varCount ? 0 : countDifferentAlleles(lastVar2, list));
					if(different1*8<varCount){
						ArrayList<Read> list2=map.get(lastVar2);
						final int varCount2=(list2==null ? 0 : list2.size());
						final int different2=(sharedCount==varCount2 ? 0 : countDifferentAlleles(var, list2));
						if(different2*8<varCount2){
							bestVar2=lastVar2;
							bestSharedCount=sharedCount;
							bestDifferent=Tools.max(different1, different2);
						}
					}
				}
				sharedCount=1;
			}
			lastVar2=currentVar2;
		}
		if(sharedCount>bestSharedCount){
			final int different1=(sharedCount==varCount ? 0 : countDifferentAlleles(lastVar2, list));
			if(different1*8<varCount){
				ArrayList<Read> list2=map.get(lastVar2);
				final int varCount2=(list2==null ? 0 : list2.size());
				final int different2=(sharedCount==varCount2 ? 0 : countDifferentAlleles(var, list2));
				if(different2*8<varCount2){
					bestVar2=lastVar2;
					bestSharedCount=sharedCount;
					bestDifferent=Tools.max(different1, different2);
				}
			}
		}
		rvector[0]=var;
		rvector[1]=list.size();
		rvector[2]=bestVar2;
		rvector[3]=sharedCount;
		rvector[4]=bestDifferent;
		
		return bestVar2;
	}
	
	static final int toReadLocation(final int clumpLocation, final int maxLeft, final int readPos){
		final int readLocation=clumpLocation+readPos-maxLeft;
		return readLocation;
	}
	
	static final int toClumpLocation(final int readLocation, final int maxLeft, final int readPos){
		final int clumpLocation=readLocation-readPos+maxLeft;
		assert(readLocation==toReadLocation(clumpLocation, maxLeft, readPos));
		return clumpLocation;
	}
	
	static int findBestPivots(Clump c, boolean findCorrelations, IntList pivots){
		pivots.clear();
		final int size=c.size();
		if(size<minSizeSplit){return -1;}
		
		if(findCorrelations){
			int x=findBestPivot_Correlated(c, pivots);
			if(x>-1){return x;}
		}
		
		final int[][] bcounts=c.baseCounts();
		final int[][] qcounts=c.qualityCounts();
		final float[][] qAvgs=c.qualityAverages();
		final int width=c.width();
		
		int bestPosition=-1;
		int bestVar=-1;
		long bestBsecond=0;
		long bestQsecond=0;

//		final float bmult=8f, bmult2=15f, qmult=(c.useQuality() ? 1.5f : 100f);
		final boolean useQuality=c.useQuality();
		final float bmult=20f, bmult2=20f;
		final float qmult=20f, qamult=1.5f;
		
		final int minPivotDepth=Tools.max(4, (int)(minSizeFractionSplit*size));
//		final int minMinorAllele=Tools.max((conservative ? 1 : 2), (int)(0.25f+size/bmult2));
		final int minMinorAllele=Tools.max(2, (int)(size/bmult2));
		final int minMinorAlleleQ=minMinorAllele*10;
		
//		assert(false) : size+", "+minSizeSplit+", "+minPivotDepth+", "+minMinorAllele;
		
		for(int i=0; i<width; i++){
			final long sum=c.getSumAtPosition(bcounts, i);
			if(sum>=minPivotDepth){
				final int pmajor=c.getConsensusAtPosition(bcounts, i);
				final int pminor=c.getSecondHighest(bcounts, i);
				if(pmajor>=0){
					final long bmajor=bcounts[pmajor][i];
					final long bsecond=bcounts[pminor][i];

					final long qmajor=qcounts[pmajor][i];
					final long qsecond=qcounts[pminor][i];

					final float qamajor=qAvgs[pmajor][i];
					final float qasecond=qAvgs[pminor][i];

					if(bsecond*bmult>=bmajor && bsecond>=minMinorAllele){
						if(!useQuality || (qsecond*qmult>=qmajor && qasecond*qamult>=qamajor && qsecond>=minMinorAlleleQ)){//candidate
							if(bsecond>bestBsecond || (bsecond==bestBsecond && qsecond>bestQsecond)){//Perhaps Qsecond should be more important...?

//								assert(false) : size+", "+minSizeSplit+", "+minPivotDepth+", "+minMinorAllele+", "+bmajor+", "+bsecond+", "+qmajor+", "+qsecond+", "+minMinorAlleleQ;
								
								bestBsecond=bsecond;
								bestQsecond=qsecond;
								bestPosition=i;
								bestVar=(bestPosition<<2)|pminor;
							}
						}
					}
				}
			}
		}
		
//		if(bestVar<0 && findCorrelations){
//			bestVar=findBestPivot_Correlated(c);
//		}
		
		if(bestVar>=0){pivots.add(bestVar);}
		return bestVar;
	}

	static int minSizeSplit=4; //5 is actually better than 4 in allele separation tests...?
	static float minSizeFractionSplit=0.17f; //0.2 is substantially worse, 0.14 is a tiny bit better than 0.17
	static boolean conservative=false;
	
	private static final int alleleMask=0x3;
	private static final int posMask=~alleleMask;
	private static final int shift=2;
	
	static boolean FIND_CORRELATIONS=true;
	static int MAX_CORRELATIONS=12;
	
//	private static final ThreadLocal<LinkedHashMap<Integer, ArrayList<Read>>> localMap=new ThreadLocal<LinkedHashMap<Integer, ArrayList<Read>>>();
//	private static final Integer[] mapKeys=new Integer[2000];
//	static{
//		for(int i=0; i<mapKeys.length; i++){mapKeys[i]=i;}
//	}
	
}
