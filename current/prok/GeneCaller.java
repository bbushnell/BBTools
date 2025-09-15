package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import aligner.IDAligner;
import aligner.SingleStateAlignerFlat2;
import aligner.SingleStateAlignerFlat3;
import aligner.SingleStateAlignerFlatFloat;
import dna.AminoAcid;
import shared.KillSwitch;
import shared.Tools;
import stream.Read;
import structures.FloatList;
import structures.IntList;
import structures.LongHashSet;


/**
 * This class calls genes within a single thread.
 * @author Brian Bushnell
 * @date Sep 24, 2018
 *
 */
public class GeneCaller extends ProkObject {
	
	/*--------------------------------------------------------------*/
	/*----------------             Init             ----------------*/
	/*--------------------------------------------------------------*/
	
	GeneCaller(int minLen_, int maxOverlapSameStrand_, int maxOverlapOppositeStrand_, 
			float minStartScore_, float minStopScore_, float minInnerScore_,
			float minOrfScore_, float minAvgScore_, GeneModel pgm_){
		minLen=minLen_;
		maxOverlapSameStrand=maxOverlapSameStrand_;
		maxOverlapOppositeStrand=maxOverlapOppositeStrand_;
		pgm=pgm_;
		
		minStartScore=minStartScore_;
		minStopScore=minStopScore_;
		minInnerScore=minInnerScore_;
		minOrfScore=minOrfScore_;
		minAvgScore=minAvgScore_;
		assert(pgm!=null);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public ArrayList<Orf> callGenes(ArrayList<Read> reads){
		ArrayList<Orf> orfs=new ArrayList<Orf>();
		for(Read r : reads) {
			ArrayList<Orf> list=callGenes(r);
			if(list!=null) {orfs.addAll(list);}
		}
		return orfs;
	}
	
	public ArrayList<Orf> callGenes(Read r){
		return callGenes(r, pgm, true);
	}
	
	public ArrayList<Orf> callGenes(Read r, boolean breakOrfs){
		return callGenes(r, pgm, breakOrfs);
	}
	
	public ArrayList<Orf> callGenes(Read r, GeneModel pgm_, boolean breakOrfs){
		pgm=pgm_;
		assert(pgm!=null);
		
		final String name=r.id;
		final byte[] bases=r.bases;

		//Lists of all longest orfs per frame
		ArrayList<Orf>[] frameLists=makeOrfs(name, bases, minLen);
		//Lists of all high-scoring orfs per frame, with potentially multiple orfs sharing stops.
		ArrayList<Orf>[] brokenLists=(breakOrfs ? breakOrfs(frameLists, bases) : frameLists);
		if(!breakOrfs) {
			for(int i=0; i<brokenLists.length; i++) {
				if(brokenLists[i]==null) {brokenLists[i]=new ArrayList<Orf>(1);}
			}
		}
		
		ArrayList<Orf>[] rnaLists=null;
		final int rlen=r.length();
		if(calltRNA || (call16S && rlen>800) || (call23S && rlen>1500) || call5S || (call18S && rlen>1000)){
			rnaLists=makeRnas(name, bases);
			if(rnaLists[0]!=null && !rnaLists[0].isEmpty()) {
				brokenLists[0].addAll(rnaLists[0]);
				Collections.sort(brokenLists[0]);
			}
			if(rnaLists[1]!=null && !rnaLists[1].isEmpty()) {
				brokenLists[3].addAll(rnaLists[1]);
				Collections.sort(brokenLists[3]);
			}
		}
		
		boolean printAllOrfs=false;
		boolean printRnas=false;
		if(printAllOrfs){
			ArrayList<Orf> temp=new ArrayList<Orf>();
			for(ArrayList<Orf> broken : brokenLists){
				temp.addAll(broken);
			}
			Collections.sort(temp);
//			System.err.print(temp.size()+"; ");
			return temp;
		}
		
		if(printRnas && rnaLists!=null){
			ArrayList<Orf> temp=new ArrayList<Orf>();
			for(ArrayList<Orf> list : rnaLists){
				temp.addAll(list);
			}
			Collections.sort(temp);
			return temp;
		}
		
		stCds2.add(brokenLists);
		
		//Find the optimal path through Orfs
		ArrayList<Orf> path=findPath(brokenLists, bases);
//		geneStartsOut+=path.size();

		if(callCDS){stCdsPass.add(path);}
		if(calltRNA){sttRNA.add(path);}
		if(call16S){st16s.add(path);}
		if(call23S){st23s.add(path);}
		if(call5S){st5s.add(path);}
		if(call18S){st18s.add(path);}
		
		return path;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Generates lists of all max-length non-overlapping Orfs per frame.
	 * There IS overlap between frames.
	 * All Orfs come out flipped to + orientation. 
	 * */
	ArrayList<Orf>[] makeRnas(String name, byte[] bases){
		@SuppressWarnings("unchecked")
		ArrayList<Orf>[] array=new ArrayList[2];
		array[0]=new ArrayList<Orf>();
		array[1]=new ArrayList<Orf>();
		final float[] scores=new float[bases.length];
		final int[] kmersSeen=(lsuKmers==null && ssuKmers==null && trnaKmers==null && r5SKmers==null) ? null : new int[bases.length];
		for(int strand=0; strand<2; strand++){
			for(StatsContainer sc : pgm.rnaContainers){
				if(ProkObject.callType(sc.type)){
					ArrayList<Orf> list=makeRnasForStrand(name, bases, strand, sc, scores, (sc.kmerSet()==null ? null : kmersSeen), false, -1);//TODO: Make this loop through all RNA types
					if(strand==1 && list!=null){
						for(Orf orf : list){
							assert(orf.strand==strand);
							orf.flip();
						}
					}
					if(list!=null){array[strand].addAll(list);}
				}
			}
			Collections.sort(array[strand]);
			AminoAcid.reverseComplementBasesInPlace(bases);
		}
		return array;
	}
	
	/** Designed for quickly calling a single SSU */
	public Orf makeRna(String name, byte[] bases, int type){
		final float[] scores=new float[bases.length];//TODO: Big and slow; make a FloatList?
		StatsContainer sc=pgm.allContainers[type];
		final int[] kmersSeen=(sc.kmerSet()==null ? null : new int[bases.length]);//TODO: IntList?
		
		int strand=0;
		ArrayList<Orf> list=makeRnasForStrand(name, bases, strand, sc, scores, kmersSeen, true, -1);
		final Orf best1=pickBest(list);
		assert(best1==null || best1.start>=0 && best1.stop<bases.length) : bases.length+"\n"+best1;
		if(best1!=null && best1.orfScore>-999){return best1;}
		
		strand++;
		AminoAcid.reverseComplementBasesInPlace(bases);
		list=makeRnasForStrand(name, bases, strand, sc, scores, kmersSeen, true, -1);
		AminoAcid.reverseComplementBasesInPlace(bases);
		if(strand==1 && list!=null){
			for(Orf orf : list){
				assert(orf.strand==strand);
				orf.flip();
			}
		}
		final Orf best2=pickBest(list);
		assert(best2==null || best2.start>=0 && best2.stop<bases.length) : bases.length+"\n"+best2;
		if(best2!=null && best2.orfScore>-999){return best2;}
		return best1!=null ? best1 : best2;
	}
	
	final Orf pickBest(ArrayList<Orf> list){
		if(list==null){return null;}
		Orf best=null;
		for(Orf orf : list){
			if(best==null || orf.orfScore>best.orfScore){
				best=orf;
			}
		}
		return best;
	}
	
	/** 
	 * Generates lists of all max-length non-overlapping Orfs per frame.
	 * There IS overlap between frames.
	 * All Orfs come out flipped to + orientation. 
	 * */
	static ArrayList<Orf>[] makeOrfs(String name, byte[] bases, int minlen){
		@SuppressWarnings("unchecked")
		int orfs=0;
		ArrayList<Orf>[] array=new ArrayList[6];
		for(int strand=0; strand<2; strand++){
			for(int frame=0; frame<3; frame++){
				ArrayList<Orf> list=makeOrfsForFrame(name, bases, frame, strand, minlen);
				orfs+=(list==null ? 0 : list.size());
				array[frame+3*strand]=list;
				if(strand==1 && list!=null){
					for(Orf orf : list){
						assert(orf.frame==frame);
						assert(orf.strand==strand);
						orf.flip();
					}
				}
			}
			AminoAcid.reverseComplementBasesInPlace(bases);
		}
//		System.err.print("orfs="+orfs+", ");
		return array;
	}
	
	/**
	 * Dynamic programming phase.
	 * @param frameLists
	 * @param bases
	 * @return
	 */
	private ArrayList<Orf> findPath(ArrayList<Orf>[] frameLists, byte[] bases){
		ArrayList<Orf> all=new ArrayList<Orf>();
		for(ArrayList<Orf> list : frameLists){all.addAll(list);}
		if(all.isEmpty()){return all;}
		Collections.sort(all);
		
		for(Orf orf : all){
			orf.pathScorePlus=-999999;
			orf.pathScoreMinus=-999999;
		}
		
		int[] lastPositionScored=KillSwitch.allocInt1D(6);
		Arrays.fill(lastPositionScored, -1);

		//Index of highest-scoring ORF in this frame, with prev on the plus strand
		int[] bestIndexPlus=KillSwitch.allocInt1D(6);
		//Index of highest-scoring ORF in this frame, with prev on the minus strand
		int[] bestIndexMinus=KillSwitch.allocInt1D(6);
		//Highest-scoring ORF in this frame, with prev on the plus strand
		Orf[] bestOrfPlus=new Orf[6];
		//Highest-scoring ORF in this frame, with prev on the minus strand
		Orf[] bestOrfMinus=new Orf[6];

		int[][] bestIndex=new int[][] {bestIndexPlus, bestIndexMinus};
		Orf[][] bestOrf=new Orf[][] {bestOrfPlus, bestOrfMinus};
		
		for(Orf orf : all){
			final int myListNum=3*orf.strand+orf.frame;
			calcPathScore(orf, frameLists, lastPositionScored, bestIndex);
			if(bestOrfPlus[myListNum]==null || orf.pathScorePlus>=bestOrfPlus[myListNum].pathScorePlus){
				bestOrfPlus[myListNum]=orf;
				bestIndexPlus[myListNum]=lastPositionScored[myListNum];
				assert(frameLists[myListNum].get(lastPositionScored[myListNum])==orf);
			}
			if(bestOrfMinus[myListNum]==null || orf.pathScoreMinus>=bestOrfMinus[myListNum].pathScoreMinus){
				bestOrfMinus[myListNum]=orf;
				bestIndexMinus[myListNum]=lastPositionScored[myListNum];
				assert(frameLists[myListNum].get(lastPositionScored[myListNum])==orf);
			}
		}
		
		Orf best=bestOrf[0][0];
		for(Orf[] array : bestOrf){
			for(Orf orf : array){
				if(best==null || (orf!=null && orf.pathScore()>best.pathScore())){
					best=orf;
				}
			}
		}
		ArrayList<Orf> bestPath=new ArrayList<Orf>();
		for(Orf orf=best; orf!=null; orf=orf.prev()){
			bestPath.add(orf);
			if(orf.type==CDS){geneStartsOut++;}
			else if(orf.type==tRNA){tRNAOut++;}
			else if(orf.type==r16S){r16SOut++;}
			else if(orf.type==r23S){r23SOut++;}
			else if(orf.type==r5S){r5SOut++;}
			else if(orf.type==r18S){r18SOut++;}
		}
		Collections.sort(bestPath);
		return bestPath;
	}
	
	/**
	 * Calculate the best path to this ORF.
	 * @param orf
	 * @param frameLists
	 * @param lastPositionScored
	 * @param bestIndex
	 */
	private void calcPathScore(Orf orf, ArrayList<Orf>[] frameLists, int[] lastPositionScored, int[][] bestIndex){
		final int myListNum=3*orf.strand+orf.frame;

//		System.err.println("* "+orf);
//		System.err.println("* "+Arrays.toString(lastPositionScored));
//		System.err.println();
		
		for(int listStrand=0; listStrand<2; listStrand++){
			for(int listFrame=0; listFrame<3; listFrame++){
				int listNum=listFrame+3*listStrand;
				ArrayList<Orf> list=frameLists[listNum];
				int lastPos=lastPositionScored[listNum];
				int bestPos=bestIndex[listStrand][listNum];
				if(listStrand==0){
					calcPathScorePlus(orf, list, listStrand, lastPos, bestPos);
				}else{
					calcPathScoreMinus(orf, list, listStrand, lastPos, bestPos);
				}
			}
		}
		
//		System.err.println(myListNum+", "+Arrays.toString(lastPositionScored)+", "+frameLists[myListNum].size());
		
		lastPositionScored[myListNum]++;
		assert(frameLists[myListNum].get(lastPositionScored[myListNum])==orf) : myListNum+"\n"+orf+"\n"+frameLists[myListNum].get(lastPositionScored[myListNum])+"\n"
			+Arrays.toString(lastPositionScored)+"\n"+frameLists[myListNum].get(lastPositionScored[myListNum]+1);
		
		//These are sanity checks to make sure that the path did not break in the middle.
		//Safe to disable.
//		assert(orf.prevPlus!=null || orf.stop<100000);
//		assert(orf.prevMinus!=null || orf.stop<100000);
//		assert(orf.pathScore>-10) : orf.pathScore+"\n"+orf+"\n"+orf.prev+"\n";
	}
	
	/**
	 * Calculate the best path to this ORF from a plus-strand previous ORF.
	 * @param orf
	 * @param list
	 * @param listStrand
	 * @param lastPos
	 * @param bestPos
	 */
	private void calcPathScorePlus(final Orf orf, final ArrayList<Orf> list, final int listStrand, final int lastPos, final int bestPos){
		assert(listStrand==0);
		if(lastPos<0){
			if(orf.prevPlus==null){
				orf.pathScorePlus=orf.orfScore;
				orf.pathLengthPlus=1;
			}
			return;
		}
		if(list.isEmpty()){return;}
		
//		System.err.println("\nExamining   \t"+orf+"\nlastPos="+lastPos+", bestPos="+bestPos+", sameFrame="+sameFrame);
		boolean found=false;
		final boolean sameStrand=(orf.strand==listStrand);
		final int maxOverlap=(sameStrand ? maxOverlapSameStrand : maxOverlapOppositeStrand);
		for(int i=lastPos, min=Tools.max(0, bestPos-lookbackPlus); i>=min || (i>0 && !found); i--){
			Orf prev=list.get(i);
			assert(prev!=orf) : prev;
//			System.err.println("Comparing to \t"+prev);
			if(orf.isValidPrev(prev, maxOverlap)){
				int overlap=Tools.max(0, prev.stop-orf.start+1);
				float orfScore=overlap==0 ? orf.orfScore : orf.calcOrfScore(overlap);
				
				final float prevScore=prev.pathScore();
				final int prevLength=prev.pathLength();
				
				float pathScore;
				final int pathLength;
				if(sameStrand){
					pathLength=prevLength+1;
					pathScore=prevScore+orfScore;
					pathScore+=p0+p1*(Tools.mid(p5*(p2+pathLength), p6*(p3-pathLength), p4));
				}else{
					pathLength=1;
					pathScore=prev.pathScore()+orfScore;
					pathScore+=q1+Tools.mid(q2*prevLength, q3+q4*prevLength, q5);
				}
				
				if(overlap<1 && prevScore>0){found=true;}
				if(pathScore>=orf.pathScorePlus){
					orf.pathScorePlus=pathScore;
					orf.prevPlus=prev;
					orf.pathLengthPlus=pathLength;
//					System.err.println("Set as best");
				}
			}
//			if(found && prev.stop<maxOverlap-2000 && orf.prevPlus!=null){
//				System.err.println("Breaking");
//				break;
//			}
		}
	}
	
	/**
	 * Calculate the best path to this ORF from a minus-strand previous ORF.
	 * @param orf
	 * @param list
	 * @param listStrand
	 * @param lastPos
	 * @param bestPos
	 */
	private void calcPathScoreMinus(final Orf orf, final ArrayList<Orf> list, final int listStrand, final int lastPos, final int bestPos){
		assert(listStrand==1);
		if(lastPos<0){
			if(orf.prevMinus==null){
				orf.pathScoreMinus=orf.orfScore;
				orf.pathLengthMinus=1;
			}
			return;
		}
		if(list.isEmpty()){return;}
		
//		System.err.println("\nExamining   \t"+orf+"\nlastPos="+lastPos+", bestPos="+bestPos+", sameFrame="+sameFrame);
		boolean found=false;
		final boolean sameStrand=(orf.strand==listStrand);
		final int maxOverlap=(sameStrand ? maxOverlapSameStrand : maxOverlapOppositeStrand);
		for(int i=lastPos, min=Tools.max(0, bestPos-lookbackMinus); i>=min || (i>0 && !found); i--){
			Orf prev=list.get(i);
			assert(prev!=orf) : prev;
//			System.err.println("Comparing to \t"+prev);
			if(orf.isValidPrev(prev, maxOverlap)){
				int overlap=Tools.max(0, prev.stop-orf.start+1);
				float orfScore=overlap==0 ? orf.orfScore : orf.calcOrfScore(overlap);
				
				final float prevScore=prev.pathScore();
				final int prevLength=prev.pathLength();
				
				float pathScore;
				final int pathLength;
				if(sameStrand){
					pathLength=prevLength+1;
					pathScore=prevScore+orfScore;
					pathScore+=p0+p1*(Tools.mid(p5*(p2+pathLength), p6*(p3-pathLength), p4));
				}else{
					pathLength=1;
					pathScore=prev.pathScore()+orfScore;
					pathScore+=q1+Tools.mid(q2*prevLength, q3+q4*prevLength, q5);
				}
				if(overlap<1 && prevScore>0){found=true;}
				if(pathScore>=orf.pathScoreMinus){
					orf.pathScoreMinus=pathScore;
					orf.prevMinus=prev;
					orf.pathLengthMinus=pathLength;
//					System.err.println("Set as best");
				}
			}
//			if(found && prev.stop<maxOverlap-2000 && orf.prevMinus!=null){
//				System.err.println("Breaking");
//				break;
//			}
		}
	}
	
	/** 
	 * Generates a list of maximal-length Orfs only (non-overlapping).
	 * All Orfs come out in native orientation (unflipped). 
	 * */
	static ArrayList<Orf> makeOrfsForFrame(String name, byte[] bases, int startFrame, int strand, int minlen){
//		assert(false) : "TODO";
		assert(minlen>=3);
		if(bases==null || bases.length<minlen){return null;}
		ArrayList<Orf> orfs=new ArrayList<Orf>();
		if(!ProkObject.callCDS){return orfs;}
//		int mask=63;
		int code=0;
		int start=-2;
		int frame=0;
		int pos=startFrame;
		
		
		for(; pos<bases.length; pos++){
			byte b=bases[pos];
			int x=AminoAcid.baseToNumber[b];
//			code=((code<<2)|x)&mask;
			code=((code<<2)|x);
			frame++;
			if(frame==3){
				frame=0;
				if(start>=0){
					if(GeneModel.isStopCodon(code) || code<0){//NOTE: This adds a stop codon wherever there are Ns.
						int len=pos-start+1;
						if(len>=minlen){
							Orf f=new Orf(name, start, pos, strand, startFrame, bases, true, CDS);
							orfs.add(f);
						}
						start=-1;
					}
				}else{
					if(start==-2 || (start<0 && GeneModel.isStartCodon(code))){
						start=pos-2;
					}
				}
				code=0;
			}
		}

		//Add a stop codon at the sequence end.
		if(start>=0){
			pos--;
			while(frame!=3 && frame!=-1){
				pos--;
				frame--;
			}
			int len=pos-start+1;
			if(len>=minlen){
				assert(pos<bases.length) : start+", "+pos+", "+bases.length;
				Orf f=new Orf(name, start, pos, strand, startFrame, bases, true, CDS);
				orfs.add(f);
			}
		}

//		System.err.print(" *"+orfs.size()+"* ");
		return orfs;
	}
	
	/** 
	 * Generates a list of maximal-length RNAs (non-overlapping).
	 * All RNAs come out in native orientation (unflipped). 
	 * */
	ArrayList<Orf> makeRnasForStrand(String name, byte[] bases, int strand, StatsContainer sc, float[] scores, int[] kmersSeen, boolean quitEarly, float bias){
		final int window=sc.lengthAvg;
		if(bases==null || bases.length*2<window){return null;}
		ArrayList<Orf> orfs=new ArrayList<Orf>(sc.type==tRNA ? 32 : 8);
		
		final FrameStats inner=sc.inner;
//		final FrameStats start=sc.start;
//		final FrameStats stop=sc.stop;
		
		final int k=inner.k;
		final int mask=inner.mask;
//		final float invLen=sc.invLengthAvg;
		final int halfWindow=window/2;
		final int maxWindow=(int)(window*1.5f);
		final int maxWindow2=(int)(window*2.5f);
//		final int slop=Tools.max(50, window/8);
		int len=0;
		int kmer=0;
		float currentScore=0;
		float currentScoreAbs=0;
		bias=(bias>-1 ? bias : biases[sc.type]);
		final float maxBias=biases[sc.type]*1.45f;
		
		float thresh=cutoff1[sc.type];
		float prevScore=0;
		int lastStart=0;
		
		float max=0;
		int maxPos=0;
		
		for(int pos=0; pos<bases.length; pos++){
			final byte b=bases[pos];
			assert(b>=0 && b<128) : "Invalid base b="+((int)b)+"; pos="+pos+"\n"+new String(bases)+"\n";
			final int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			
			if(x>=0){
				len++;
				if(len>=k){
					float prob=inner.probs[0][kmer];
					float dif=prob-bias;//Prob above 1 is more likely than average
					currentScoreAbs+=prob;
					currentScore=Tools.max(0, currentScore+dif);
				}
				
				if(currentScore>0){
					if(currentScore>max){
						max=currentScore;
						maxPos=pos;
					}
					if(prevScore<=0){
						lastStart=pos;
					}
				}else{
					int rnaLen=maxPos-lastStart;
					if(max>thresh && rnaLen>=halfWindow){
						if(rnaLen>maxWindow){
							if(bias<=maxBias){
								orfs=null;
								float biasMult=(rnaLen>8*window ? 1.2f : rnaLen>4*window ? 1.1f : 1.05f);
								return makeRnasForStrand(name, bases, strand, sc, scores, kmersSeen, quitEarly, bias*biasMult);
							}
						}
						if(rnaLen<=maxWindow2){
							Orf orf=new Orf(name, lastStart, maxPos, strand, 0, bases, false, sc.type);
							orfs.add(orf);
							orf.orfScore=max;
							if(verbose){
								final int start2=(strand==0 ? lastStart : bases.length-maxPos-1);
								final int stop2=(strand==0 ? maxPos : bases.length-lastStart-1);
								System.err.println("Made Orf "+start2+"\t"+stop2+"\t"+max);
							}
						}
					}
					max=0;
					lastStart=pos;
				}
//				System.err.println("i="+i+"\tscore="+score+"\tmax="+max+"\tmaxPos="+maxPos+"\tprevScore="+prevScore+"\tlastStart="+lastStart);
				prevScore=currentScore;
				
//				if(pos>=223000 && pos<232000){
////					System.err.println("i="+i+"\tscore="+score+"\tmax="+max+"\tmaxPos="+maxPos+"\tprevScore="+prevScore+"\tlastStart="+lastStart);
//					System.out.println(pos+"\t"+currentScore);
//				}
				
			}else{
				len=0;
				kmer=0;
			}
			scores[pos]=currentScoreAbs;
		}
		
//		System.err.println("size="+orfs.size()+", type="+Orf.typeStrings[sc.type]);
		
		
		{
			int rnaLen=maxPos-lastStart;
			if(max>thresh && rnaLen>=halfWindow){
				if(rnaLen>maxWindow){
					if(bias<=maxBias){
						orfs=null;
						float biasMult=(rnaLen>8*window ? 1.2f : rnaLen>4*window ? 1.1f : 1.05f);
						return makeRnasForStrand(name, bases, strand, sc, scores, kmersSeen, quitEarly, bias*biasMult);
					}
				}
				if(rnaLen<=maxWindow2){
					Orf orf=new Orf(name, lastStart, maxPos, strand, 0, bases, false, sc.type);
					orfs.add(orf);
					orf.orfScore=max;
					if(verbose){
						final int start2=(strand==0 ? lastStart : bases.length-maxPos-1);
						final int stop2=(strand==0 ? maxPos : bases.length-lastStart-1);
						System.err.println(start2+"\t"+stop2+"\t"+max);
					}
				}
			}
		}
		
		if(kmersSeen!=null && orfs.size()>0 && sc.kmerSet()!=null){
			fillKmersSeen(bases, kmersSeen, sc.kmerSet(), sc.kLongLen());
		}
		
		float cutoff=cutoff2[sc.type];
		
		for(int i=0; i<orfs.size(); i++){
			Orf orf=orfs.get(i);
//			System.err.println(orf.orfScore);
			boolean good=refineRna(orf, bases, strand, sc, scores, kmersSeen);
			if(orf.orfScore<cutoff || !good){
				if(verbose){System.err.println("REJECT: "+orf.toStringFlipped());}
				orfs.set(i, null);
			}else{
				if(verbose){System.err.println("ACCEPT: "+orf.toStringFlipped());}
				if(quitEarly){
					orfs.clear();
					orfs.add(orf);
					return orfs;
				}
			}
		}
		Tools.condenseStrict(orfs);
		
//		assert(false);
		
//		for(int pos=0; pos<bases.length; pos++){
//			final byte b=bases[pos];
//			final int x=AminoAcid.baseToNumber[b];
//			kmer=((kmer<<2)|x)&mask;
//			
//			if(x>=0){
//				len++;
//				if(len>=k){
//					float prob=inner.probs[0][kmer];
//					float dif=prob-1.2f;//Prob above 1 is more likely than average
//					currentScore=Tools.max(0, currentScore+dif);
//					if(currentScore>0){
//						currentStreak++;
//					}else{
//						currentStreak=0;
//					}
//					if(currentScore>200 && currentStreak>1500 && currentStreak<1700){
//						Orf orf=new Orf(name, pos-currentStreak-1, pos, strand, 0, bases, false);
//						orfs.add(orf);
//						orf.orfScore=currentScore;
//						orf.startScore=start.scorePoint(orf.start, bases);
//						orf.stopScore=stop.scorePoint(orf.stop, bases);
//						currentStreak=0;
//						currentScore=0;
//					}
//				}
//			}else{
//				len=0;
//				kmer=0;
//			}
//		}
		
		return orfs;
	}
	
	void fillKmersSeen(byte[] bases, int[] kmersSeen, LongHashSet set, int k){
		final long mask=~((-1L)<<(2*k));
		long kmer=0;
		int len=0;
		int seen=0;
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			final int num=AminoAcid.baseToNumber[b];
			if(num>=0){
				len++;
				kmer=((kmer<<2)|num)&mask;
				if(len>=k && set.contains(kmer)){seen++;}
			}else{
				len=0;
			}
			kmersSeen[i]=seen;
		}
	}
	
	boolean refineRna(Orf orf, byte[] bases, int strand, StatsContainer sc, float[] scores, int[] kmersSeen){
		if(orf==null){return false;}
		if(verbose){System.err.println("REFINE: "+orf.toStringFlipped());}
		final int window=sc.lengthAvg;
//		final int halfWindow=window/2;
//		final int maxWindow=(int)(window*1.5f);
		final int slop=Tools.max(60, 10+window/8);
		final float invWindow=sc.invLengthAvg;
		final float oldScore=orf.orfScore;
		IntList starts=new IntList();
		IntList stops=new IntList();
		FloatList startScores=new FloatList();
		FloatList stopScores=new FloatList();

		final int leftmost=Tools.max(0, orf.start-slop);
		final int rightmost=Tools.min(bases.length-1, orf.stop+slop);
		if(kmersSeen!=null){
			if(kmersSeen[leftmost]>=kmersSeen[rightmost]){
//				System.err.println("Bad: "+oldScore);
				orf.orfScore=-999;
				return false;
			}else{
//				System.err.println("Good: "+oldScore);
			}
		}
		if(verbose){System.err.println("REFINE2");}
		
		{//starts
			final int left=leftmost;
			final int right=Tools.min(bases.length-1, orf.stop+slop-window);
			final float thresh=cutoff3[sc.type];
			fillPoints(left, right, bases, sc.start, thresh, starts, startScores);
		}
		if(verbose){System.err.println("starts: "+starts.size);}
//		if((orf.start+"").startsWith("146") || true){System.err.println(starts);}
		if(verbose){System.err.println(startScores);}

		{//stops
			final int left=Tools.max(0, orf.start-slop+window);
			final int right=rightmost;
			final float thresh=cutoff4[sc.type];
			fillPoints(left, right, bases, sc.stop, thresh, stops, stopScores);
		}
		if(verbose){System.err.println("stops: "+stops.size);}
//		if((orf.start+"").startsWith("146") || true){System.err.println(stops);}
		if(verbose){System.err.println(stopScores);}

		final float innerThresh=cutoff5[sc.type];
		
		final int minlen=Tools.max(window/2, window-slop);
		final int maxlen=window+slop;
		
		orf.orfScore=0;
		int lastStop=0;
		for(int i=0; i<starts.size; i++){
			final int start=starts.get(i);
			final int startSeen=kmersSeen==null ? 0 : kmersSeen[start];
			final float startScore=startScores.get(i);
//			System.err.println("start="+start);
			for(int j=lastStop; j<stops.size; j++){
				final int stop=stops.get(j);
				final int rnalen=stop-start+1;
//				System.err.println("stop="+stop);
				if(rnalen<minlen){
					lastStop=j+1;
				}else if(rnalen>maxlen){
//					System.err.println("broke");
					break;
				}else if(kmersSeen!=null && kmersSeen[stop]<=startSeen){//TODO: Test this
//					//skip
				}else{
					final int len=stop-start+1;
					final float stopScore=stopScores.get(j);
					final float innerScore=(scores[stop]-scores[start])/len;
//					System.err.println("innerScore="+innerScore);
					assert(rnalen<=maxlen) : "start="+start+", stop="+stop+", rnalen="+rnalen+", minlen="+minlen+", maxlen="+maxlen;
					if(innerScore>=innerThresh){
						final float a=Tools.max(startScore+0.25f, 0.25f);
						final float b=Tools.max(stopScore+0.25f, 0.25f);
						final float c=innerScore*innerScore;
						final float d=(window-(2.4f*Tools.absdif(len, window)))*invWindow;
						final float score=c*d*(float)Math.sqrt(a*b);
						if(verbose && score>=0.2f*orf.orfScore){
							final int start2=(strand==0 ? start : bases.length-stop-1);
							final int stop2=(strand==0 ? stop : bases.length-start-1);
							System.err.println(start2+"-"+stop2+", "+startScore+", "+stopScore+", "+innerScore+", "+(score*scoreMult[sc.type])+", "+oldScore);
						}
						if(score>orf.orfScore){
							orf.start=start;
							orf.stop=stop;
							orf.kmerScore=innerScore*len;
//							if(verbose){
//								final int start2=(strand==0 ? start : bases.length-stop-1);
//								final int stop2=(strand==0 ? stop : bases.length-start-1);
//								System.err.println(start2+"-"+stop2+", "+startScore+", "+stopScore+", "+innerScore+", "+(score*scoreMult[sc.type])+", "+oldScore);
//							}
							orf.startScore=startScore;
							orf.stopScore=stopScore;
							orf.orfScore=score;
						}
					}else{
						assert(true);
					}
				}
			}
		}
		orf.orfScore*=scoreMult[sc.type];
		
		if(starts.isEmpty() || stops.isEmpty()){
			if(verbose){System.err.println("No starts or stops.");}
			orf.orfScore=Tools.min(-999, orf.orfScore-9999);
			return false;
		}
		return refineByAlignment(orf, bases, strand, sc);//Returns true if no consensus is present
	}
	
	boolean refineByAlignment(Orf orf, byte[] bases, int strand, StatsContainer sc){
		if(verbose){System.err.println("ALIGN");}
		Read[] consensus=ProkObject.consensusReads(sc.type);
		if(consensus==null || consensus.length==0){return true;}
		boolean refined=false;
//		System.err.println("Initial: "+orf.start+", "+orf.stop);
		for(Read r : consensus){
//			refined=refineByAlignment(orf, bases, strand, sc, r.bases, 15, 15, 2);
			refined=refineByAlignment(orf, bases, strand, sc, r.bases, sc.startSlop(), sc.stopSlop(), 2);
			if(refined || sc.type==r18S || true){break;}
		}
		if(refined){
			if(verbose){System.err.println("Aligned to: "+orf.start+", "+orf.stop);}
		}else{
			if(verbose){System.err.println("Alignment failed.");}
			orf.orfScore=Tools.min(-999, orf.orfScore-9999);
		}
		return refined;
	}
	
	boolean refineByAlignment(Orf orf, byte[] bases, int strand, StatsContainer sc, byte[] consensus, 
			final int startSlop, final int stopSlop, int recurLimit){
		if(useIDAligner) {
			return refineByAlignment_IDA(orf, bases, strand, 
					sc, consensus, startSlop, stopSlop, recurLimit);
		}else {
			return refineByAlignment_SSA(orf, bases, strand, 
					sc, consensus, startSlop, stopSlop, recurLimit);
		}
	}
	
	boolean refineByAlignment_IDA(Orf orf, byte[] bases, int strand, StatsContainer sc, byte[] consensus, 
			final int startSlop, final int stopSlop, int recurLimit){
		final int start0=orf.start;
		final int stop0=orf.stop;
		
		assert(start0>=0 && start0<bases.length) : start0+", "+stop0;
		assert(stop0>=0 && stop0<bases.length) : start0+", "+stop0;
		
		final float minID=sc.minIdentity();
		final int padding=Tools.min(alignmentPaddingIDA, 30+sc.lengthAvg/4);
		final int a=Tools.max(0, orf.start-padding);
		final int b=Tools.min(bases.length-1, orf.stop+padding);
		final int reflen=b-a+1;
		assert(reflen>0) : reflen;
		if(reflen>10*sc.lengthAvg && reflen>20000){
			System.err.println("Skipped reflen "+reflen+"/"+sc.lengthAvg+" for "
					+ "seqlen="+bases.length+", orf="+orf.toString());
			assert(false);
			//TODO: Possibly change return to -1, 0, 1 ("can't align")
			//Should be a limit on window size...
			//Also consider shrinking matrix after jumbo alignments
			return false;
		}
		assert(a>=0 && b<bases.length) : a+", "+b;
		IDAligner ida=aligner.Factory.makeIDAligner();
		int[] pos=new int[2];
		float id=ida.align(consensus, bases, pos, a, b);
		if(id<minID){return false;}
		
		final int rstart=Tools.max(pos[0], 0);
		final int rstop=Tools.min(pos[1], bases.length-1);
		
		if(Tools.absdif(rstart, start0)>startSlop){orf.start=rstart;}
		if(Tools.absdif(rstop, stop0)>stopSlop){orf.stop=rstop;}
		assert(orf.length()>0) : "\nreflen="+reflen+", orflen="+orf.length()+
			", id"+id+", minid="+minID+", strand="+orf.strand+"\n"+
			"qlen="+consensus.length+", rlen="+bases.length+", a="+a+", b="+b+"\n"+
			"pos="+Arrays.toString(pos)+", start0="+start0+", stop0="+stop0+"\n"+orf;
			//+new String(consensus)+"\n"+new String(bases, a, b-a+1)+"\n";
		return orf.length()>0;
	}
	
	boolean refineByAlignment_SSA(Orf orf, byte[] bases, int strand, StatsContainer sc, byte[] consensus, 
			final int startSlop, final int stopSlop, int recurLimit){
		final int start0=orf.start;
		final int stop0=orf.stop;
		
		assert(start0>=0 && start0<bases.length) : start0+", "+stop0;
		assert(stop0>=0 && stop0<bases.length) : start0+", "+stop0;
		
		final float minID=sc.minIdentity();
		final int padding=Tools.min(alignmentPadding, 30+sc.lengthAvg/4);
		final int a=Tools.max(0, orf.start-padding);
		final int b=Tools.min(bases.length-1, orf.stop+padding);
		final int reflen=b-a+1;
		if(reflen>10*sc.lengthAvg && reflen>20000){
			System.err.println("Skipped reflen "+reflen+"/"+sc.lengthAvg+" for "
					+ "seqlen="+bases.length+", orf="+orf.toString());
			assert(false);
			//TODO: Possibly change return to -1, 0, 1 ("can't align")
			//Should be a limit on window size...
			//Also consider shrinking matrix after jumbo alignments
			return false;
		}
		assert(a>=0 && b<bases.length) : a+", "+b;
		SingleStateAlignerFlat2 ssa=getSSA();
		final int minScore=ssa.minScoreByIdentity(consensus.length, minID);
		int[] max=ssa.fillUnlimited(consensus, bases, a, b, minScore);
		if(max==null){return false;}
		
		final int rows=max[0];
		final int maxCol=max[1];
		final int maxState=max[2];
//		final int maxScore=max[3];
//		final int maxStart=max[4];
		
		//returns {score, bestRefStart, bestRefStop} 
		//padded: {score, bestRefStart, bestRefStop, padLeft, padRight};
		final int[] score=ssa.score(consensus, bases, a, b, rows, maxCol, maxState);
		final int rstart=Tools.max(score[1], 0);
		final int rstop=Tools.min(score[2], bases.length-1);
		final float id=ssa.tracebackIdentity(consensus, bases, a, b, rows, maxCol, maxState, null);
		
		assert(rstart>=0 && rstart<bases.length) : "bases="+bases.length+
			", rstart="+rstart+", rstop="+rstop+", a="+a+", b="+b+"\n"+"score="+Arrays.toString(score)+", id="+id;
		assert(rstop>=0 && rstop<bases.length) : rstart+", "+rstop;
		
		if(score.length>3 && recurLimit>0){
			final int padLeft=score[3];
			final int padRight=score[4];
			//TODO:  This takes extra memory.  It may be better to cap the width or ignore start0/stop0 here.
			int rstart2=Tools.max(0, Tools.min(start0, rstart)-10);
			int rstop2=Tools.min(bases.length-1, Tools.max(stop0, rstop)+10);
			assert(rstart2>=0 && rstart2<bases.length) : rstart2+", "+rstop2;
			assert(rstop2>=0 && rstop2<bases.length) : rstart2+", "+rstop2;
			orf.start=rstart2;
			orf.stop=rstop2;
			if(orf.start<a || orf.stop>b){
				boolean ret=refineByAlignment(orf, bases, strand, sc, consensus, 0, 0, recurLimit-1);
				if(ret){return true;}
			}
			orf.start=start0;
			orf.stop=stop0;
//			assert(false) : "Don't traceback after recur.";
		}
		if(verbose){
			System.err.println("Identity: "+id);
		}
//		assert(score.length==3) : "TODO: Handle padding requests.";
		
//		System.err.println("Identity: "+String.format("%.2f", 100*id)+"; location: "+rstart+"-"+rstop);
		if(id<minID){return false;}
		
		
		if(Tools.absdif(rstart, start0)>startSlop){orf.start=rstart;}
		if(Tools.absdif(rstop, stop0)>stopSlop){orf.stop=rstop;}
		return true;
	}
	
	void fillPoints(final int left, final int right, final byte[] bases, final FrameStats fs, float thresh, final IntList points, final FloatList scores){
		points.clear();
		scores.clear();
		final float minThresh=thresh;//thresh*0.05f;
//		System.err.println(left+", "+right+", "+thresh+", "+minThresh);
		while(points.size()<8 && thresh>=minThresh){
//			System.err.println("Running with thresh="+thresh);
			points.clear();
			scores.clear();
			for(int i=left; i<right; i++){
				float score=fs.scorePoint(i, bases);
//				System.err.println(i+", "+score);
				if(score>=thresh){
					points.add(i);
					scores.add(score);
				}
			}
			thresh=thresh*0.75f;
		}
	}
	
	/**
	 * Generate all possible genes from each Orf, and return them in a new set of lists.
	 * @param frameLists
	 * @param bases
	 * @return Lists of orfs.
	 */
	private ArrayList<Orf>[] breakOrfs(ArrayList<Orf>[] frameLists, byte[] bases){

		@SuppressWarnings("unchecked")
		ArrayList<Orf>[] brokenLists=new ArrayList[6];
		for(int strand=0; strand<2; strand++){
			for(int frame=0; frame<3; frame++){
				int fnum=frame+3*strand;
				ArrayList<Orf> longest=frameLists[fnum]; //Longest Orf per stop
				ArrayList<Orf> broken=new ArrayList<Orf>(); //All high scoring Orfs, including multiple per stop, with different starts
				if(longest!=null) {
					for(Orf orf : longest){
						assert(orf.frame==frame);
						assert(orf.strand==strand);
						ArrayList<Orf> temp=breakOrf(orf, bases);
						if(temp!=null){
							broken.addAll(temp);
						}
					}
				}
				Collections.sort(broken);
				brokenLists[fnum]=broken;
			}
			//Reverse-complement bases after processing each strand
			AminoAcid.reverseComplementBasesInPlace(bases);
		}
		return brokenLists;
	}
	
	/**
	 * Generate an Orf for each possible start codon.
	 * Retain only the high-scoring ones.
	 * @param longest Longest open reading frame for a given stop.
	 * @param bases Bases, oriented for this Orf.
	 * @return List of Orfs.
	 */
	private ArrayList<Orf> breakOrf(Orf longest, byte[] bases){
		assert(longest.start<longest.stop);
		final int flipped=longest.flipped();
		if(flipped==1){longest.flip();}//Now the orf is aligned to its native strand
		
		geneStopsMade++;
		
		final FrameStats innerStats=pgm.statsCDS.inner;
		final FrameStats startStats=pgm.statsCDS.start;
		final FrameStats stopStats=pgm.statsCDS.stop;
		
		final String name=longest.scafName;
		final int start=longest.start;
		final int stop=longest.stop;
		final int strand=longest.strand;
		final int max=Tools.min(longest.stop-2, longest.stop-minLen+4);

		assert(pgm!=null) : pgm;
		assert(pgm.statsCDS!=null) : pgm;
		assert(pgm.statsCDS.inner!=null) : pgm.statsCDS;
		assert(pgm.statsCDS.inner.k>0) : pgm.statsCDS.inner;
		
		final int k=innerStats.k;
		final int mask=~((-1)<<(2*k));

		final float stopScore=stopStats.scorePoint(longest.stop, bases);
		stCds.geneStopScoreSum+=stopScore;
		stCds.geneStopScoreCount++;
		
		ArrayList<Orf> broken=new ArrayList<Orf>();
		int created=0;
		
		int codon=0;
		int kmer=0;
		int len=0;
		int numKmers=0;
		float currentScore=0;
		for(int pos=start, currentFrame=0; pos<=stop; pos++){
			final byte b=bases[pos];
			final int x=AminoAcid.baseToNumber[b];
			codon=((codon<<2)|x);
			kmer=((kmer<<2)|x)&mask;
			
			if(x>=0){
				len++;
				if(len>=k){
					float prob=innerStats.probs[currentFrame][kmer];
					float dif=prob-0.99f;//Prob above 1 is more likely than average
					currentScore+=dif;
//					outstream.println("pos="+pos+"\tdif="+Tools.format("%.4f", dif)+",\tscore="+Tools.format("%.4f", currentScore)+
//							"\tasStart="+Tools.format("%.4f", pgm.calcStartScore(pos-2, bases))+"\tasStop="+Tools.format("%.4f", stopStats.scorePoint(pos, bases))+
//							"\tcodon="+AminoAcid.kmerToString(kmer, 3)+" frame="+(currentFrame));
				}else{
//					outstream.println("pos="+pos+"\tdif="+Tools.format("%.4f", 0.0)+",\tscore="+Tools.format("%.4f", 0.0)+
//							"\tasStart="+Tools.format("%.4f", pgm.calcStartScore(pos-2, bases))+"\tasStop="+Tools.format("%.4f", stopStats.scorePoint(pos, bases))+
//							"\tcodon="+AminoAcid.kmerToString(kmer, 3)+" frame="+(currentFrame));
				}
			}else{
				len=0;
				kmer=0;
			}
			
			currentFrame++;
//			outstream.println("pos="+pos+", codon="+AminoAcid.kmerToString(kmer, 3)+", frame="+currentFrame+", start="+start+", isStartCodon="+pgm.isStartCodon(codon));
			if(currentFrame>2){
				currentFrame=0;
				if(pos<max && created<breakLimit && (pos==start+2 || pgm.isStartCodon(codon))){
//					outstream.println(x);
					int glen=stop-pos+3;
					assert(glen>=minLen) : "glen="+glen+", minLen="+minLen+", pos="+pos+", max="+max+", start="+start;
					
					int oStart=pos-2;
					float startScore=startStats.scorePoint(oStart, bases);

					stCds.geneStartScoreSum+=startScore;
					stCds.geneStartScoreCount++;
					
					stCds.lengthSum+=(stop-(pos-2)+1);
					stCds.lengthCount++;
					
					if((startScore>=minStartScore || pos<6) /* && stopScore>=minStopScore /*|| broken.isEmpty()*/){
						Orf orf=new Orf(name, pos-2, stop, strand, longest.frame, bases, false, longest.type);
						
						geneStartsMade++;
						orf.kmerScore=currentScore;
						orf.startScore=startScore;
						orf.stopScore=stopScore;

						assert(orf.frame==longest.frame);
						assert(orf.strand==strand);

						if(strand==1){orf.flip();}
						broken.add(orf);
						created++;
					}
				}
				codon=0;
			}
		}
		
		final int size=broken.size();
		final int sizeCutoff=Tools.max(5, size/2);
		if(size<1){return broken;}
		Orf best=broken.get(0);
		Orf bestStart=broken.get(0);
		for(int i=0; i<size; i++){
			Orf orf=broken.get(i);
			
			//This fixes scores because they were generated together, from start to stop, to make this O(N) instead of O(N^2).
			orf.kmerScore=currentScore-orf.kmerScore;
			orf.orfScore=orf.calcOrfScore();
			if(orf.orfScore>=best.orfScore){best=orf;}
			if(orf.startScore>=bestStart.startScore){bestStart=orf;}
			
			stCds.geneInnerScoreSum+=orf.averageKmerScore();
			stCds.geneInnerScoreCount++;
		}
		
		//Sort to by score descending to eliminate low-scoring copies. 
		if(broken.size()>1){Collections.sort(broken, PFeature.featureComparatorScore);}
		
		int removed=0;
		for(int i=0; i<size; i++){//Fix scores because they were generated together, from start to stop, to make this O(N) instead of O(N^2).
			Orf orf=broken.get(i);
			if(keepAtLeastOneOrf && i==0) {
				//skip
			}else if(orf.averageKmerScore()<minInnerScore || orf.orfScore<minOrfScore || 
					orf.orfScore/orf.length()<minAvgScore || 
					orf.orfScore<0.5f*best.orfScore-10 || (orf.startScore<bestStart.startScore-0.55f && orf.kmerScore<best.kmerScore*1.1f && orf!=best)){
				broken.set(i, null);
				removed++;
			}else if(i>sizeCutoff){
				broken.set(i, null);
				removed++;
			}
		}
		if(removed>0){
			Tools.condenseStrict(broken);
		}
		
		geneStartsRetained+=broken.size();
		geneStopsRetained+=(broken.size()>0 ? 1 : 0);
		
		if(flipped==1){longest.flip();}
		return broken;
	}
	
	/**
	 * Return a score of how well this sequence matches a feature based on inner kmer frequencies.
	 * Intended to see whether a short read matches a feature, or which side matches better.
	 * Generally this would be called twice, with bases being reverse-complemented in between.
	 * @param bases Bases, oriented on the intended strand.
	 * @param type Type of feature, e.g. ProkObject.CDS.
	 * @return Average score of best-scoring frame for this strand.
	 */
	public float scoreFeature(byte[] bases, int type){
		assert(pgm!=null) : pgm;
		StatsContainer sc=pgm.allContainers[type];
		assert(sc!=null) : type+", "+pgm;
		final FrameStats innerStats=pgm.allContainers[type].inner;
		assert(innerStats!=null) : sc;
		assert(innerStats.k>0) : innerStats;
		final int frames=(type==CDS ? 3 : 1);
		
		final int k=innerStats.k;
		final int mask=~((-1)<<(2*k));
		
		int codon=0;
		int kmer=0;
		int len=0;
		int numKmers=0;
		float[] frameScores=new float[frames];
		for(int pos=0, currentFrame=0; pos<bases.length; pos++){
			final byte b=bases[pos];
			final int x=AminoAcid.baseToNumber[b];
			codon=((codon<<2)|x);
			kmer=((kmer<<2)|x)&mask;
			
			if(x>=0){
				len++;
				if(len>=k){
					for(int frame2=0; frame2<frames; frame2++) {
						int frame=(currentFrame+frame2)%frames;
						float prob=innerStats.probs[frame][kmer];
//						float dif=prob-0.99f;//Prob above 1 is more likely than average
						frameScores[frame]+=prob;
						numKmers++;
					}
				}
			}else{
				len=0;
				kmer=0;
			}
			
			currentFrame++;
			if(currentFrame>2){currentFrame=0;}
		}
//		System.err.println(Arrays.toString(frameScores));
		return (Tools.max(frameScores)*frames)/Tools.max(1, numKmers);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** 
	 * Current gene model.
	 * TODO: Dynamically swap this as needed for contigs with varying GC.
	 */
	GeneModel pgm;
	
	//Gene-calling cutoffs
	final int minLen;
	final int maxOverlapSameStrand;
	final int maxOverlapOppositeStrand;
	final float minStartScore;
	final float minStopScore;
	final float minInnerScore;
	final float minOrfScore;
	final float minAvgScore;

//	static float[] cutoff1=new float[] {0, 40f, 200f, 150f, 45f};
//	static float[] cutoff2=new float[] {0, 44f, 300f, 150f, 60f};
//	static float[] cutoff3=new float[] {0, 1.7f, 1.5f, 1.4f, 1.8f};
//	static float[] cutoff4=new float[] {0, 1.6f, 1.5f, 0.6f, 1.1f};
//	static float[] cutoff5=new float[] {0, 1.0f, 1.0f, 1.0f, 1.55f};//inner score
//	static float[] scoreMult=new float[] {1f, 8f, 200f, 175f, 12f};
//	static float[] biases=new float[] {1f, 1.17f, 1.2f, 1.2f, 1.15f};

//	static float[] cutoff1=new float[] {0, 40f, 140f, 150f, 45f};
//	static float[] cutoff2=new float[] {0, 44f, 300f, 150f, 45f};
//	static float[] cutoff3=new float[] {0, 1.7f, 1.1f, 1.3f, 1.8f};
//	static float[] cutoff4=new float[] {0, 1.6f, 1.3f, 0.5f, 1.1f};
//	static float[] cutoff5=new float[] {0, 1.0f, 1.0f, 1.0f, 1.3f};//inner score
//	static float[] scoreMult=new float[] {1f, 8f, 200f, 175f, 12f};
//	static float[] biases=new float[] {1f, 1.17f, 1.2f, 1.2f, 1.15f};

	//for k=4,2,2 
//	static float[] cutoff1=new float[] {0, 40f, 140f, 150f, 40f};
//	static float[] cutoff2=new float[] {0, 44f, 300f, 150f, 45f};
//	static float[] cutoff3=new float[] {0, 1.7f, 1.1f, 1.1f, 1.8f};
//	static float[] cutoff4=new float[] {0, 1.6f, 1.3f, 0.45f, 1.1f};
//	static float[] cutoff5=new float[] {0, 1.0f, 1.0f, 1.0f, 1.3f};//inner score
//	static float[] scoreMult=new float[] {1f, 8f, 200f, 175f, 12f};
//	static float[] biases=new float[] {1f, 1.17f, 1.2f, 1.2f, 1.22f};

//	//for k=5,3,3
//	static float[] cutoff1=new float[] {0, 40f, 300f, 300f, 135f};//sum score
//	static float[] cutoff2=new float[] {0, 44f, 300f, 400f, 100f};//orf score
//	static float[] cutoff3=new float[] {0, 1.7f, 1.5f, 1.5f, 3.5f};//start score
//	static float[] cutoff4=new float[] {0, 1.6f, 1.5f, 1.4f, 1.8f};//stop score
//	static float[] cutoff5=new float[] {0, 1.0f, 1.1f, 1.15f, 1.4f};//inner score
//	static float[] scoreMult=new float[] {1f, 8f, 80f, 120f, 5f};//score mult
//	static float[] biases=new float[] {1f, 1.17f, 1.2f, 1.2f, 1.22f};//bias for sum score

	//for k=6,3,3 - this has low scores for correct 23s stops; may try k=2 or k=4 for that end 
	//Also, 5S seems to score low very for archaea
//	public static int CDS=0, tRNA=1, /*r16S=2,*/ r23S=3, r5S=4, r18S=5, r28S=6, RNA=7;
//	static float[] cutoff1=new float[] {0, 100f, 600f, 400f, 300f};//sum score
//	static float[] cutoff2=new float[] {0, 48f, 300f, 300f, 32f};//orf score
//	static float[] cutoff3=new float[] {0, 3.0f, 1.8f, 1.6f, 4.0f};//start score
//	static float[] cutoff4=new float[] {0, 2.8f, 2.0f, 0.90f, 1.9f};//stop score
//	static float[] cutoff5=new float[] {0, 2.8f, 1.1f, 1.55f, 1.4f};//inner score
//	static float[] scoreMult=new float[] {1f, 1f, 35f, 80f, 1.0f};//score mult
//	static float[] biases=new float[] {1f, 1.25f, 1.30f, 1.30f, 1.22f};//16S bias for sum score: 1.25 best for archs, 1.4 best for bacts

////	{"CDS", "tRNA", "16S", "23S", "5S", "18S", "28S", "RNA"};
//	static float[] cutoff1=new float[] {0, 90f, 300f, 400f, 270f, 300f};//sum score
//	static float[] cutoff2=new float[] {0, 48f, 300f, 300f, 32f, 300f};//orf score
//	static float[] cutoff3=new float[] {0, 2.8f, 1.65f, 1.6f, 3.8f, 1.65f};//start score
//	static float[] cutoff4=new float[] {0, 2.6f, 1.70f, 0.90f, 1.8f, 1.70f};//stop score
//	static float[] cutoff5=new float[] {0, 2.8f, 1.70f, 1.55f, 1.4f, 1.70f};//inner score
//	static float[] scoreMult=new float[] {1f, 1f, 35f, 80f, 1.0f, 35f};//score mult
//	static float[] biases=new float[] {1f, 1.50f, 1.30f, 1.30f, 1.30f, 1.30f};

////	{"CDS", "tRNA", "16S", "23S", "5S", "18S", "28S", "RNA"};
//	static float[] cutoff1=new float[] {0, 90f, 300f, 400f, 90f, 300f};//sum score
//	static float[] cutoff2=new float[] {0, 48f, 300f, 300f, 24f, 300f};//orf score
//	static float[] cutoff3=new float[] {0, 2.8f, 1.65f, 1.6f, 2.0f, 1.65f};//start score
//	static float[] cutoff4=new float[] {0, 2.6f, 1.70f, 0.90f, 1.0f, 1.70f};//stop score
//	static float[] cutoff5=new float[] {0, 2.8f, 1.70f, 1.55f, 2.6f, 1.70f};//inner score
//	static float[] scoreMult=new float[] {1f, 1f, 35f, 80f, 1.25f, 35f};//score mult
//	static float[] biases=new float[] {1f, 1.50f, 1.30f, 1.30f, 1.50f, 1.30f};

////	{"CDS", "tRNA", "16S", "23S", "5S", "18S", "28S", "RNA"};
//	static float[] cutoff1=new float[] {0, 90f, 300f, 400f, 100f, 300f};//sum score
//	static float[] cutoff2=new float[] {0, 48f, 300f, 300f, 32f, 300f};//orf score
//	static float[] cutoff3=new float[] {0, 2.8f, 1.65f, 1.6f, 1.8f, 1.65f};//start score
//	static float[] cutoff4=new float[] {0, 2.6f, 1.70f, 0.90f, 1.0f, 1.70f};//stop score
//	static float[] cutoff5=new float[] {0, 2.8f, 1.70f, 1.55f, 2.6f, 1.70f};//inner score
//	static float[] scoreMult=new float[] {1f, 1f, 35f, 80f, 1.25f, 35f};//score mult
//	static float[] biases=new float[] {1f, 1.50f, 1.30f, 1.30f, 1.55f, 1.30f};

//	{"CDS", "tRNA", "16S", "23S", "5S", "18S", "28S", "RNA"};
	static float[] cutoff1=new float[] {0, 20f, 300f, 400f, 100f, 400f};//sum score
	static float[] cutoff2=new float[] {0, 36f, 300f, 300f, 32f, 300f};//orf score //tRNA is very sensitive here
	static float[] cutoff3=new float[] {0, 2.4f, 1.65f, 1.6f, 1.8f, 1.5f};//start score
	static float[] cutoff4=new float[] {0, 1.5f, 1.70f, 0.90f, 1.0f, 1.1f};//stop score
	static float[] cutoff5=new float[] {0, 2.2f, 1.70f, 1.55f, 2.6f, 1.5f};//inner score
	static float[] scoreMult=new float[] {1f, 1.0f, 35f, 80f, 1.25f, 35f};//score mult
	static float[] biases=new float[] {1f, 1.45f, 1.30f, 1.30f, 1.55f, 1.50f};
	
	long geneStopsMade=0;
	long geneStartsMade=0;
	long geneStartsRetained=0;
	long geneStopsRetained=0;
	long geneStartsOut=0;

	long tRNAOut=0;
	long r16SOut=0;
	long r23SOut=0;
	long r5SOut=0;
	long r18SOut=0;	
	
	ScoreTracker stCds=new ScoreTracker(CDS);
	ScoreTracker stCds2=new ScoreTracker(CDS);
	ScoreTracker stCdsPass=new ScoreTracker(CDS);
	ScoreTracker sttRNA=new ScoreTracker(tRNA);
	ScoreTracker st16s=new ScoreTracker(r16S);
	ScoreTracker st23s=new ScoreTracker(r23S);
	ScoreTracker st5s=new ScoreTracker(r5S);
	ScoreTracker st18s=new ScoreTracker(r18S);
	
	ScoreTracker[] trackers=new ScoreTracker[] {stCdsPass, sttRNA, st16s, st23s, st5s, st18s};
	
	public static SingleStateAlignerFlat2 getSSA(){
		SingleStateAlignerFlat2 ssa=localSSA.get();
		if(ssa==null){
			synchronized(localSSA){
				ssa=localSSA.get();
				if(ssa==null){
					ssa=new SingleStateAlignerFlat2();
					localSSA.set(ssa);
				}
			}
		}
		return ssa;
	}
	
	public static SingleStateAlignerFlat3 getSSA3(){
		SingleStateAlignerFlat3 ssa=localSSA3.get();
		if(ssa==null){
			synchronized(localSSA3){
				ssa=localSSA3.get();
				if(ssa==null){
					ssa=new SingleStateAlignerFlat3();
					localSSA3.set(ssa);
				}
			}
		}
		return ssa;
	}
	
	public static SingleStateAlignerFlatFloat getSSAF(){
		SingleStateAlignerFlatFloat ssa=localSSAF.get();
		if(ssa==null){
			synchronized(localSSAF){
				ssa=localSSAF.get();
				if(ssa==null){
					ssa=new SingleStateAlignerFlatFloat();
					localSSAF.set(ssa);
				}
			}
		}
		return ssa;
	}
	
	public boolean keepAtLeastOneOrf=false;

	private static ThreadLocal<SingleStateAlignerFlat2> localSSA=new ThreadLocal<SingleStateAlignerFlat2>();
	private static ThreadLocal<SingleStateAlignerFlat3> localSSA3=new ThreadLocal<SingleStateAlignerFlat3>();
	private static ThreadLocal<SingleStateAlignerFlatFloat> localSSAF=new ThreadLocal<SingleStateAlignerFlatFloat>();
//	public static int maxAlignmentEndpointDifference=15;
	public static int alignmentPadding=300;
	public static int alignmentPaddingIDA=80;//TODO: Recur for edge alignments, and add start hint
	public static boolean useIDAligner=false;//Quantum/Drifting are fastest
	
	public static int breakLimit=12;
	public static int lookbackPlus=70;
	public static int lookbackMinus=25;
	
//	pathScore+=p0+p1*(Tools.mid(p5*(p2+pathLength), p6*(p3-pathLength), p4));
	
	//Operon length modifiers for same strand
	public static float p0=-30f;
	public static float p1=-0.35f; //Sensitive - higher decreases FP and TP
	public static float p2=4.0f;//Insensitive - higher positive decreases FP and TP
	public static float p3=12f; //Higher decreases FP (substantially) and TP (slightly)
	public static float p4=-10f; //Lower decreases FP (weakly) and TP (greatly)
	public static float p5=2.0f; //Insensitive - lower increases FP and TP
	public static float p6=2f; //Greater increases FP and TP
	
//	pathScore+=q1+Tools.mid(q2*prevLength, q3+q4*prevLength, q5);
	
	//Operon length modifiers for opposite strand
	public static float q1=-36f;
	public static float q2=-1.6f; //q2 and q4 must have opposite signs
	public static float q3=-12f;
	public static float q4=3.0f; //(Lower [even negative] decreases FP and TP)
	public static float q5=-40f;
	
	private static PrintStream outstream=System.err;
	static boolean verbose;
	
}
