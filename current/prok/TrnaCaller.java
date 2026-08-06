package prok;

import java.util.ArrayList;
import java.util.Arrays;

import consensus.BaseGraph;
import dna.AminoAcid;
import idaligner.ScrabbleAligner;
import shared.Tools;
import structures.FloatList;
import structures.IntList;

/**
 * Dedicated tRNA gene finder optimized for clustered tRNAs.
 * Unlike the general RNA caller, this handles arbitrarily large
 * score-positive regions and extracts multiple non-overlapping
 * tRNAs from each cluster via greedy selection.
 * @author Neptune, Brian Bushnell
 */
public class TrnaCaller extends ProkObject {

	public TrnaCaller(GeneModel pgm_){
		this(pgm_, null);
	}

	public TrnaCaller(GeneModel pgm_, byte[][] trnaLibrary_){
		this(pgm_, trnaLibrary_, null);
	}

	public TrnaCaller(GeneModel pgm_, byte[][] trnaLibrary_, BaseGraph[] trnaModels_){
		this(pgm_, trnaLibrary_, trnaModels_, null);
	}

	public TrnaCaller(GeneModel pgm_, byte[][] trnaLibrary_, BaseGraph[] trnaModels_, String[] modelNames_){
		pgm=pgm_;
		sc=pgm.statstRNA;
		inner=sc.inner;
		k=inner.k;
		mask=inner.mask;
		trnaLibrary=trnaLibrary_;
		trnaModels=trnaModels_;
		modelNames=modelNames_;
		annotate=(modelNames!=null);
		kmerIndex=(trnaLibrary!=null ? buildKmerIndex(trnaLibrary) : null);
	}

	public static long alignmentCount(){return alignmentCount;}

	public ArrayList<Orf> callTrnas(String name, byte[] bases, int strand){
		ArrayList<Orf> results=new ArrayList<>();
		if(bases==null || bases.length<MIN_TRNA){return results;}

		float[] scores=scanInner(bases);
		ArrayList<int[]> regions=findRegions(scores, bases);

		for(int[] region : regions){
			ArrayList<Orf> found=extractTrnas(name, bases, strand, scores, region[0], region[1]);
			results.addAll(found);
		}

		if(trnaLibrary!=null && trnaLibrary.length>0){
			results=verifyByAlignment(results, bases);
		}

		return results;
	}

	/**
	 * Scans the sequence with the tRNA inner k-mer model.
	 * Returns cumulative scores for use in inner-score computation.
	 */
	private float[] scanInner(byte[] bases){
		float[] cumulative=new float[bases.length];
		int kmer=0, len=0;
		float sum=0;
		for(int i=0; i<bases.length; i++){
			int x=AminoAcid.baseToNumber[bases[i]];
			kmer=((kmer<<2)|x)&mask;
			if(x>=0){
				len++;
				if(len>=k){
					sum+=inner.probs[0][kmer];
				}
			}else{
				len=0; kmer=0;
			}
			cumulative[i]=sum;
		}
		return cumulative;
	}

	/**
	 * Finds score-positive regions using a local accumulator.
	 * No max-window filter — clusters of any size are kept.
	 */
	private ArrayList<int[]> findRegions(float[] cumulative, byte[] bases){
		ArrayList<int[]> regions=new ArrayList<>();
		int kmer=0, len=0;
		float currentScore=0, prevScore=0, max=0;
		int lastStart=0, maxPos=0;
		float bias=GeneCaller.biases[tRNA];

		for(int pos=0; pos<bases.length; pos++){
			int x=AminoAcid.baseToNumber[bases[pos]];
			kmer=((kmer<<2)|x)&mask;
			if(x>=0){
				len++;
				if(len>=k){
					float prob=inner.probs[0][kmer];
					currentScore=Tools.max(0, currentScore+(prob-bias));
				}
				if(currentScore>0){
					if(currentScore>max){max=currentScore; maxPos=pos;}
					if(prevScore<=0){lastStart=pos;}
				}else{
					int rnaLen=maxPos-lastStart;
					if(max>REGION_THRESH && rnaLen>=MIN_TRNA){
						regions.add(new int[]{lastStart, maxPos});
					}
					max=0; lastStart=pos;
				}
				prevScore=currentScore;
			}else{
				len=0; kmer=0;
			}
		}
		if(max>REGION_THRESH && (maxPos-lastStart)>=MIN_TRNA){
			regions.add(new int[]{lastStart, maxPos});
		}
		return regions;
	}

	/**
	 * Extracts multiple non-overlapping tRNAs from a score-positive region.
	 * Enumerates all valid (start,stop) pairs using the start/stop point models,
	 * scores each, then greedily selects non-overlapping ones by score.
	 */
	private ArrayList<Orf> extractTrnas(String name, byte[] bases, int strand,
			float[] cumulative, int regionStart, int regionStop){
		ArrayList<Orf> results=new ArrayList<>();

		final int slop=Tools.max(30, sc.lengthAvg/4);
		final int left=Tools.max(0, regionStart-slop);
		final int right=Tools.min(bases.length-1, regionStop+slop);

		IntList starts=new IntList();
		IntList stops=new IntList();
		FloatList startScores=new FloatList();
		FloatList stopScores=new FloatList();

		fillPoints(left, right, bases, sc.start, GeneCaller.cutoff3[tRNA], starts, startScores);
		fillPoints(left, right, bases, sc.stop, GeneCaller.cutoff4[tRNA], stops, stopScores);

		if(starts.isEmpty() || stops.isEmpty()){return results;}

		final int window=sc.lengthAvg;
		final float invWindow=sc.invLengthAvg;
		final int minlen=Tools.max(MIN_TRNA, window/2);
		final int maxlen=MAX_TRNA;

		ArrayList<float[]> candidates=new ArrayList<>();

		for(int i=0; i<starts.size; i++){
			final int start=starts.get(i);
			final float sScore=startScores.get(i);
			for(int j=0; j<stops.size; j++){
				final int stop=stops.get(j);
				final int len=stop-start+1;
				if(len<minlen){continue;}
				if(len>maxlen){break;}
				final float pScore=stopScores.get(j);
				final float innerScore=(cumulative[stop]-cumulative[start])/len;
				if(innerScore>=INNER_THRESH){
					final float a=Tools.max(sScore+0.25f, 0.25f);
					final float b=Tools.max(pScore+0.25f, 0.25f);
					final float c=innerScore*innerScore;
					final float d=(window-(2.4f*Tools.absdif(len, window)))*invWindow;
					final float score=c*d*(float)Math.sqrt(a*b);
					if(score>CANDIDATE_THRESH){
						candidates.add(new float[]{start, stop, score, sScore, pScore, innerScore*len});
					}
				}
			}
		}

		if(candidates.isEmpty()){return results;}
		candidates.sort((x, y)->Float.compare(y[2], x[2]));

		ArrayList<int[]> accepted=new ArrayList<>();
		for(float[] cand : candidates){
			final int cStart=(int)cand[0], cStop=(int)cand[1];
			boolean overlaps=false;
			for(int[] prev : accepted){
				if(cStart<=prev[1] && cStop>=prev[0]){overlaps=true; break;}
			}
			if(overlaps){continue;}

			Orf orf=new Orf(name, cStart, cStop, strand, 0, bases, false, tRNA);
			orf.orfScore=cand[2]*GeneCaller.scoreMult[tRNA];
			orf.startScore=cand[3];
			orf.stopScore=cand[4];
			orf.kmerScore=cand[5];

			accepted.add(new int[]{cStart, cStop});
			results.add(orf);
		}
		return results;
	}

	private void fillPoints(int left, int right, byte[] bases, FrameStats fs,
			float thresh, IntList points, FloatList scores){
		points.clear();
		scores.clear();
		final float minThresh=thresh*0.5f;
		while(points.size()<16 && thresh>=minThresh){
			points.clear();
			scores.clear();
			for(int i=left; i<=right; i++){
				float score=fs.scorePoint(i, bases);
				if(score>=thresh){
					points.add(i);
					scores.add(score);
				}
			}
			thresh=thresh*0.75f;
		}
	}

	private ArrayList<Orf> verifyByAlignment(ArrayList<Orf> candidates, byte[] bases){
		ArrayList<Orf> verified=new ArrayList<>();
		final int topN=INDEX_TOP_N_OVERRIDE>0 ? INDEX_TOP_N_OVERRIDE : INDEX_TOP_N_DEFAULT;
		for(Orf orf : candidates){
			byte[] seq=Arrays.copyOfRange(bases, orf.start, orf.stop+1);

			// Get shortlist via k-mer index, sorted by hit count descending
			int[] shortlist=shortlistByKmer(seq, topN);

			// Fast pass: align in order, with patience-based early exit
			float bestId=0;
			int bestModel=-1;
			int[] borderlineModels=new int[shortlist.length];
			int borderlineCount=0;
			boolean passed=false;
			int sinceImproved=0;

			for(int j=0; j<shortlist.length; j++){
				int i=shortlist[j];
				float id=ScrabbleAligner.alignStatic(seq, trnaLibrary[i], null);
				alignmentCount++;
				if(id>bestId){bestId=id; bestModel=i; sinceImproved=0;}
				else{sinceImproved++;}
				if(id>=ID_PASS){
					passed=true;
					if(earlyExit && sinceImproved>=earlyExitPatience){break;}
				}
				if(trnaModels!=null && id>=ID_BORDERLINE){
					borderlineModels[borderlineCount++]=i;
				}
			}

			if(passed){
				if(annotate && bestModel>=0 && modelNames!=null && bestModel<modelNames.length){
					orf.trnaModel=modelNames[bestModel];
				}
				verified.add(orf);
			}else if(bestId>=ID_BORDERLINE && trnaModels!=null && borderlineCount>0){
				float bestHbm=-999;
				int bestHbmModel=-1;
				for(int j=0; j<borderlineCount; j++){
					int idx=borderlineModels[j];
					if(idx<trnaModels.length){
						float score=TrnaConsensusBuilder.scoreAgainstModel(seq, trnaModels[idx]);
						if(score>bestHbm){bestHbm=score; bestHbmModel=idx;}
					}
				}
				if(bestHbm>=HBM_PASS){
					if(annotate && bestHbmModel>=0 && modelNames!=null && bestHbmModel<modelNames.length){
						orf.trnaModel=modelNames[bestHbmModel];
					}
					verified.add(orf);
				}
			}
		}
		return verified;
	}

	private static int[][] buildKmerIndex(byte[][] library){
		final int numKmers=1<<(2*INDEX_K);
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] lists=new ArrayList[numKmers];
		for(int i=0; i<numKmers; i++){lists[i]=new ArrayList<>();}

		for(int refId=0; refId<library.length; refId++){
			byte[] seq=library[refId];
			int kmer=0, len=0;
			for(int j=0; j<seq.length; j++){
				int x=AminoAcid.baseToNumber[seq[j]];
				if(x>=0){
					kmer=((kmer<<2)|x)&(numKmers-1);
					len++;
					if(len>=INDEX_K){
						lists[kmer].add(refId);
					}
				}else{
					len=0; kmer=0;
				}
			}
		}

		int[][] index=new int[numKmers][];
		for(int i=0; i<numKmers; i++){
			ArrayList<Integer> list=lists[i];
			if(list.isEmpty()){index[i]=EMPTY; continue;}
			int[] arr=new int[list.size()];
			for(int j=0; j<list.size(); j++){arr[j]=list.get(j);}
			index[i]=arr;
		}
		return index;
	}

	private int[] shortlistByKmer(byte[] seq, int topN){
		if(kmerIndex==null || trnaLibrary==null){
			int[] all=new int[trnaLibrary.length];
			for(int i=0; i<all.length; i++){all[i]=i;}
			return all;
		}
		final int numKmers=1<<(2*INDEX_K);
		int[] hits=new int[trnaLibrary.length];
		int kmer=0, len=0;
		for(int j=0; j<seq.length; j++){
			int x=AminoAcid.baseToNumber[seq[j]];
			if(x>=0){
				kmer=((kmer<<2)|x)&(numKmers-1);
				len++;
				if(len>=INDEX_K){
					int[] refIds=kmerIndex[kmer];
					for(int i=0; i<refIds.length; i++){
						hits[refIds[i]]++;
					}
				}
			}else{
				len=0; kmer=0;
			}
		}

		// Find top-N by hit count, cut off at minhits
		final int minHits=INDEX_MINHITS_OVERRIDE>0 ? INDEX_MINHITS_OVERRIDE : INDEX_MINHITS_DEFAULT;
		int[][] scored=new int[trnaLibrary.length][2];
		for(int i=0; i<trnaLibrary.length; i++){
			scored[i][0]=i;
			scored[i][1]=hits[i];
		}
		Arrays.sort(scored, (a, b)->b[1]-a[1]);

		int count=0;
		int limit=Tools.min(topN, trnaLibrary.length);
		while(count<limit && scored[count][1]>=minHits){count++;}
		if(count<1 && scored[0][1]>0){count=1;}
		int[] result=new int[Tools.max(0, count)];
		for(int i=0; i<result.length; i++){
			result[i]=scored[i][0];
		}
		return result;
	}

	/*--------------------------------------------------------------*/

	private final GeneModel pgm;
	private final StatsContainer sc;
	private final FrameStats inner;
	private final int k;
	private final int mask;
	private final byte[][] trnaLibrary;
	private final BaseGraph[] trnaModels;
	private final String[] modelNames;
	private final boolean annotate;
	private final int[][] kmerIndex;

	private static final int MIN_TRNA=40;
	private static final int MAX_TRNA=120;
	private static final int INDEX_K=5;
	static int INDEX_TOP_N_OVERRIDE=-1;
	static int INDEX_MINHITS_OVERRIDE=-1;
	private static final int INDEX_TOP_N_DEFAULT=20;
	private static final int INDEX_MINHITS_DEFAULT=3;
	private static final int[] EMPTY=new int[0];
	private static final float REGION_THRESH=GeneCaller.cutoff1[tRNA];
	private static final float INNER_THRESH=GeneCaller.cutoff5[tRNA];
	private static final float CANDIDATE_THRESH=GeneCaller.cutoff2[tRNA];
	static float ID_PASS=0.75f;
	static float ID_BORDERLINE=0.65f;
	static float HBM_PASS=0.70f;
	static boolean earlyExit=true;
	static int earlyExitPatience=10;
	static long alignmentCount=0;
}
