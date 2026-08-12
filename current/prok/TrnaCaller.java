package prok;

import java.util.ArrayList;
import java.util.Arrays;

import consensus.BaseGraph;
import dna.AminoAcid;
import idaligner.AlignmentStats;
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
		acPositions=(annotate && trnaLibrary!=null ? findAnticodonPositions(trnaLibrary, modelNames) : null);
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
		final int maxlen=(MAX_TRNA_OVERRIDE>0 ? MAX_TRNA_OVERRIDE : MAX_TRNA);
		//The length-deviation term d goes negative near ~1.4x window, which is the
		//REAL length cap (MAX_TRNA never binds).  When the cap is raised, floor d
		//so long candidates (e.g. intron-containing archaeal tRNAs) survive scoring
		//and let alignment verification decide their fate.
		final boolean longMode=(maxlen>MAX_TRNA);

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
					float d=(window-(2.4f*Tools.absdif(len, window)))*invWindow;
					if(longMode && len>window){d=Tools.max(d, LONG_D_FLOOR);}
					final float score=c*d*(float)Math.sqrt(a*b);
					if(score>CANDIDATE_THRESH){
						candidates.add(new float[]{start, stop, score, sScore, pScore, innerScore*len});
					}
				}
			}
		}

		if(candidates.isEmpty()){return results;}
		candidates.sort((x, y)->Float.compare(y[2], x[2]));

		//Verify DURING greedy selection: a candidate that fails alignment verification
		//does not block the overlapping candidates it outscored, so a valid suppressed
		//tRNA still gets its chance.  Verification cost only rises on the failure path.
		final boolean verify=(trnaLibrary!=null && trnaLibrary.length>0);
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

			if(verify && !verifyOrf(orf, bases)){continue;}

			//Record the POST-TRIM span: bases trimmed off a verified tRNA no longer
			//block adjacent candidates, so tightly packed neighbors (tRNA arrays)
			//that were overlapped only by scanner slop can still be called.
			accepted.add(new int[]{orf.start, orf.stop});
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

	/**
	 * Verifies a single tRNA candidate by alignment against the library,
	 * annotating it (model name, extracted anticodon) on success.
	 * @return True if the candidate passed verification.
	 */
	private boolean verifyOrf(Orf orf, byte[] bases){
		final int topN=INDEX_TOP_N_OVERRIDE>0 ? INDEX_TOP_N_OVERRIDE : INDEX_TOP_N_DEFAULT;
		byte[] seq=Arrays.copyOfRange(bases, orf.start, orf.stop+1);

		// Get shortlist via k-mer index, sorted by hit count descending
		int[] shortlist=shortlistByKmer(seq, topN);

		// Fast pass: align in order, with patience-based early exit.
		// Model choice for annotation/trimming ranks by identity*lengthSimilarity:
		// plain identity favors TRUNCATED models (a short model matching perfectly
		// inside the candidate scores ~1.0 via the swap path, beating the correct
		// full-length model), which then poisons boundary trimming.
		float bestId=0, bestCombo=0;
		int bestModel=-1, bestTrimModel=-1;
		int[] borderlineModels=new int[shortlist.length];
		int borderlineCount=0;
		boolean passed=false;
		int sinceImproved=0;

		for(int j=0; j<shortlist.length; j++){
			int i=shortlist[j];
			float id=ScrabbleAligner.alignStatic(seq, trnaLibrary[i], null);
			alignmentCount++;
			final float lenSim=Tools.min(seq.length, trnaLibrary[i].length)/(float)Tools.max(seq.length, trnaLibrary[i].length);
			final float combo=id*lenSim;
			if(id>bestId){bestId=id; bestModel=i; sinceImproved=0;}
			else{sinceImproved++;}
			if(combo>bestCombo){bestCombo=combo; bestTrimModel=i;}
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
				annotateAndTrim(orf, bases, bestModel, bestTrimModel);
			}
			return true;
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
					annotateAndTrim(orf, bases, bestHbmModel, bestHbmModel);
				}
				return true;
			}
		}
		return false;
	}

	/**
	 * Precomputes the anticodon position of each library consensus via structural scan,
	 * constrained to the triplet parsed from the model name when available.
	 * @return Per-model 0-based anticodon start positions; -1 where none was confidently found.
	 */
	public static int[] findAnticodonPositions(byte[][] library, String[] names){
		int[] positions=new int[library.length];
		for(int i=0; i<library.length; i++){
			byte[] ac=(names!=null && i<names.length ? parseLibraryAnticodon(names[i]) : null);
			positions[i]=findAnticodonPosition(library[i], ac);
		}
		return positions;
	}

	/**
	 * Parses the anticodon triplet from a consensus name like "tRNA_consensus_UCC_c2 n=5",
	 * falling back to RefSeq-style parsing for custom libraries.
	 * @return Uppercase DNA-alphabet triplet, or null if unavailable (e.g. UNK groups).
	 */
	public static byte[] parseLibraryAnticodon(String name){
		if(name==null){return null;}
		String token;
		int idx=name.indexOf("tRNA_consensus_");
		if(idx>=0){
			int start=idx+15;
			int end=start;
			while(end<name.length() && name.charAt(end)!='_' && name.charAt(end)!=' '){end++;}
			token=name.substring(start, end);
		}else{
			token=TrnaConsensusBuilder.parseAnticodon(name);
		}
		if(token==null || token.length()!=3){return null;}
		byte[] ac=new byte[3];
		for(int i=0; i<3; i++){
			int x=AminoAcid.baseToNumber[token.charAt(i)];//handles lowercase and U->T
			if(x<0){return null;}//amino-acid names and UNK land here
			ac[i]=AminoAcid.numberToBase[x];
		}
		return ac;
	}

	/**
	 * Structural plausibility score for an anticodon loop with the anticodon at position t.
	 * Scores the 5bp anticodon stem (WC pair=2, GU wobble=1), U33 immediately 5' (+3),
	 * and purine 37 immediately 3' of the triplet (+2); max 15.
	 * @return Score, or -1 if t is too close to the sequence ends to have a full arm.
	 */
	public static int scoreAnticodonLoop(byte[] seq, int t){
		if(t<7 || t+9>=seq.length){return -1;}
		final byte[] bton=AminoAcid.baseToNumber;
		int score=0;
		for(int d=3; d<=7; d++){//stem pairs (t-3,t+5)..(t-7,t+9) around the 7nt loop
			int x=bton[seq[t-d]], y=bton[seq[t+2+d]];
			if(x<0 || y<0){continue;}
			if(x+y==3){score+=2;}
			else if((x==2 && y==3) || (x==3 && y==2)){score+=1;}
		}
		if(bton[seq[t-1]]==3){score+=3;}//U33
		int p=bton[seq[t+3]];
		if(p==0 || p==2){score+=2;}//purine 37
		return score;
	}

	/**
	 * Locates the anticodon in a consensus sequence.  When the triplet is known from
	 * the model name, only matching positions are considered; otherwise pure structure
	 * is used with a stricter threshold.
	 * @return 0-based anticodon start position, or -1 if no confident position exists.
	 */
	public static int findAnticodonPosition(byte[] seq, byte[] ac){
		final byte[] bton=AminoAcid.baseToNumber;
		final int thresh=(ac!=null ? AC_FIND_MATCH : AC_FIND_BLIND);
		final int minT=15, maxT=Tools.min(60, seq.length-10);//anticodon sits ~27-40 from the 5' end
		int best=-1, bestScore=-1;
		for(int t=minT; t<=maxT; t++){
			if(ac!=null && (bton[seq[t]]!=bton[ac[0]] || bton[seq[t+1]]!=bton[ac[1]] || bton[seq[t+2]]!=bton[ac[2]])){continue;}
			int s=scoreAnticodonLoop(seq, t);
			if(s<thresh){continue;}
			if(s>bestScore || (s==bestScore && Tools.absdif(t, 34)<Tools.absdif(best, 34))){bestScore=s; best=t;}
		}
		return best;
	}

	/**
	 * Aligns a verified candidate to a library model with traceback, using whichever
	 * sequence is shorter as the query (the glocal aligner requires query<=ref).
	 */
	private AlignmentStats alignToModel(byte[] seq, int model){
		final byte[] cons=trnaLibrary[model];
		AlignmentStats stats=new AlignmentStats(true);
		stats.doTrace=true;
		if(seq.length>=cons.length){ScrabbleAligner.alignAndTraceStatic(cons, seq, stats);}
		else{ScrabbleAligner.alignAndTraceStatic(seq, cons, stats);}
		alignmentCount++;
		return stats;
	}

	/**
	 * Trims a verified tRNA's boundaries to the extent of its consensus alignment,
	 * cutting scanner slop and overhang beyond the family's depth-trimmed termini.
	 * Only applies when the candidate is the reference side of the alignment
	 * (candidate at least as long as the consensus), where rStart/rStop are
	 * candidate coordinates.
	 */
	/**
	 * Aligns against an EXTENDED window around the call (so boundary refinement can
	 * extend a late-started call upstream, not only shrink slop), then trims and
	 * extracts.  Custody is split: nameModel (identity-best) annotates and hosts
	 * anticodon extraction; trimModel (identity*lengthSimilarity-best) drives the
	 * boundary trim, since plain identity favors truncated models via the swap path.
	 */
	private void annotateAndTrim(Orf orf, byte[] bases, int nameModel, int trimModel){
		final int xFrom=Tools.max(0, orf.start-TRIM_EXT);
		final int xTo=Tools.min(bases.length-1, orf.stop+TRIM_EXT);
		byte[] seqX=Arrays.copyOfRange(bases, xFrom, xTo+1);
		AlignmentStats stats=alignToModel(seqX, nameModel);
		if(trimToAlignment){
			if(trimModel<0 || trimModel==nameModel){trimOrf(orf, seqX, xFrom, nameModel, stats);}
			else{trimOrf(orf, seqX, xFrom, trimModel, alignToModel(seqX, trimModel));}
		}
		if(extractAnticodons){orf.trnaAnticodon=extractAnticodon(seqX, nameModel, stats);}
	}

	private void trimOrf(Orf orf, byte[] seqX, int xFrom, int model, AlignmentStats stats){
		if(stats.matchString==null){return;}
		if(seqX.length<trnaLibrary[model].length){return;}//window within consensus span
		final int rStart=stats.rStart, rStop=stats.rStop;
		if(rStart<0 || rStop<=rStart || rStop>=seqX.length){return;}
		//Acceptor-stem snap: alignment boundaries carry a systematic late bias from
		//score-tie resolution in the DP, so refine both cuts by 7bp acceptor-stem
		//pairing — the tRNA's own structure defines its termini.  The scanner's
		//boundaries are RefSeq-trained and usually exact, so without a CONFIDENT
		//stem the scanner boundaries are kept (never the biased raw alignment).
		int bestS=-1, bestE=-1, bestScore=-1, bestDist=999;
		for(int s=Tools.max(0, rStart-5); s<=rStart+5; s++){
			for(int e=Tools.max(s+20, rStop-4); e<=Tools.min(seqX.length-1, rStop+4); e++){
				int sc=acceptorStemScore(seqX, s, e);
				int dist=Tools.absdif(s, rStart)+Tools.absdif(e, rStop);
				if(sc>bestScore || (sc==bestScore && dist<bestDist)){
					bestScore=sc; bestS=s; bestE=e; bestDist=dist;
				}
			}
		}
		if(bestScore<AC_STEM_MIN){return;}//no confident stem: keep scanner boundaries
		orf.start=xFrom+bestS;
		orf.stop=xFrom+bestE;
	}

	/**
	 * Best acceptor-stem pairing score for a candidate spanning [s,e]: pairs the
	 * first 7 bases against the 3' stem side, testing both the CCA-encoded register
	 * (discriminator+CCA after the stem) and the CCA-absent register.  Max 14.
	 */
	public static int acceptorStemScore(byte[] seq, int s, int e){
		final byte[] bton=AminoAcid.baseToNumber;
		//A literal CCA tail anchors the CCA register: without this, the two registers
		//can tie and shuffle the 3' cut by +-3 in CCA-encoding genomes.
		final boolean cca=(e>=2 && bton[seq[e-2]]==1 && bton[seq[e-1]]==1 && bton[seq[e]]==0);
		int best=0;
		for(int mode=0; mode<2; mode++){
			final int off=(mode==0 ? 4 : 1);//3' partner of 5' base 0: e-4 with CCA, e-1 without
			if(e-off-6<=s+6){continue;}
			int sc=(mode==0 && cca ? 2 : 0);
			for(int k=0; k<7; k++){
				int x=bton[seq[s+k]], y=bton[seq[e-off-k]];
				if(x<0 || y<0){continue;}
				if(x+y==3){sc+=2;}
				else if((x==2 && y==3) || (x==3 && y==2)){sc+=1;}
			}
			if(sc>best){best=sc;}
		}
		return best;
	}

	/**
	 * Extracts the anticodon triplet from a verified candidate by projecting the model's
	 * anticodon position through a traceback alignment, then validating anticodon-loop
	 * structure on the candidate itself.  Any failure returns null, which leaves the
	 * model-name annotation as the sole source (identical to pre-extraction behavior).
	 * @param stats Traceback alignment from alignToModel (shared with trimOrf).
	 * @return Uppercase DNA triplet, or null if extraction failed validation.
	 */
	private String extractAnticodon(byte[] seq, int model, AlignmentStats stats){
		if(acPositions==null || model<0 || model>=acPositions.length){return null;}
		final int acPos=acPositions[model];
		if(acPos<0){return null;}
		final byte[] cons=trnaLibrary[model];
		final boolean consIsQuery=(seq.length>=cons.length);
		if(stats==null || stats.matchString==null){return null;}
		//Walk the match string (aligned region only; ref-consuming ops start at rStart),
		//projecting consensus positions acPos..acPos+2 onto the candidate.
		int consPos=(consIsQuery ? 0 : stats.rStart), candPos=(consIsQuery ? stats.rStart : 0);
		int qa=-1, qc=-1;
		for(byte b : stats.matchString){
			switch(b){
				case 'm': case 'S': case 'N':
					if(consPos==acPos){qa=candPos;}
					else if(consPos==acPos+2){qc=candPos;}
					consPos++; candPos++;
					break;
				case 'D': if(consIsQuery){candPos++;}else{consPos++;} break;
				case 'I': if(consIsQuery){consPos++;}else{candPos++;} break;
			}
			if(consPos>acPos+2){break;}
		}
		if(qa<0 || qc!=qa+2){return null;}//anticodon not cleanly aligned (indel or clipped)
		//Local refinement: projection locates the anticodon to ~+-1 (banded-path wobble);
		//the loop structure resolves the exact register.  Without this, anticodons starting
		//with U produce shifted extractions whose stolen U34 passes the U33 check.
		int bestT=qa, bestScore=scoreAnticodonLoop(seq, qa), secondScore=-1;
		for(int t=qa-1; t<=qa+1; t+=2){
			int s=scoreAnticodonLoop(seq, t);
			if(s>bestScore){secondScore=bestScore; bestScore=s; bestT=t;}
			else if(s>secondScore){secondScore=s;}
		}
		if(bestScore<AC_VALIDATE){return null;}
		if(bestScore-secondScore<AC_MARGIN){return null;}//ambiguous register: decline, keep model annotation
		final byte[] bton=AminoAcid.baseToNumber;
		byte[] triplet=new byte[3];
		for(int i=0; i<3; i++){
			int x=bton[seq[bestT+i]];
			if(x<0){return null;}
			triplet[i]=AminoAcid.numberToBase[x];
		}
		return new String(triplet);
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
	/** Per-model anticodon start position in the consensus, or -1; null when not annotating */
	private final int[] acPositions;

	private static final int MIN_TRNA=40;
	private static final int MAX_TRNA=120;
	/** Raises the candidate length cap (flag maxtrna=); enables the d-floor for
	 * over-window candidates so intron-containing tRNAs can reach verification */
	static int MAX_TRNA_OVERRIDE=-1;
	/** Floor for the length-deviation score term in long mode */
	private static final float LONG_D_FLOOR=0.05f;
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
	static float HBM_PASS=0.75f;
	static boolean earlyExit=true;
	static int earlyExitPatience=10;
	static long alignmentCount=0;
	/** Master switch for structural anticodon extraction */
	static boolean extractAnticodons=true;
	/** Trim verified tRNA boundaries to the consensus alignment extent */
	static boolean trimToAlignment=true;
	/** Min structural score to accept a projected anticodon on a candidate */
	static int AC_VALIDATE=12;
	/** Min score margin between best and runner-up register; smaller = ambiguous, decline */
	static int AC_MARGIN=3;
	/** Min acceptor-stem pairing score (max 14) to snap trim boundaries to the stem */
	static int AC_STEM_MIN=10;
	/** Extension of the alignment window beyond the call for boundary refinement */
	static int TRIM_EXT=10;
	/** Min structural score for a name-matched anticodon position in a consensus */
	static int AC_FIND_MATCH=5;
	/** Min structural score for a structure-only anticodon position (UNK groups) */
	static int AC_FIND_BLIND=12;
}
