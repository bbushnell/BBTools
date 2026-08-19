package prok;

import java.util.ArrayList;
import java.util.Arrays;

import consensus.BaseGraph;
import dna.AminoAcid;
import idaligner.AlignmentStats;
import idaligner.QuantumAligner;
import idaligner.ScrabbleAligner;
import map.LongHashSet;
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
		kmerIndex=(trnaLibrary!=null ? new TrnaKmerIndex(trnaLibrary, INDEX_K, ADAPTIVE_MINHITS,
			ADAPT_FLOOR, ADAPT_TOPFRAC, ADAPT_QFRAC,
			(INDEX_MINHITS_OVERRIDE>0 ? INDEX_MINHITS_OVERRIDE : INDEX_MINHITS_DEFAULT)) : null);
		acPositions=(annotate && trnaLibrary!=null ? findAnticodonPositions(trnaLibrary, modelNames) : null);
	}

	public static long alignmentCount(){return alignmentCount;}

	public ArrayList<Orf> callTrnas(String name, byte[] bases, int strand){
		ArrayList<Orf> results=new ArrayList<>();
		if(bases==null || bases.length<MIN_TRNA){return results;}

		float[] scores=scanInner(bases);
		ArrayList<int[]> regions=findRegions(scores, bases);

		if(DEBUG){for(int[] rg : regions){System.err.println("REGION\t"+name+"\t"+strand+"\t"+rg[0]+"\t"+rg[1]);}}

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

		//Track the top-3 long (intron-suspect) start/stop pairs by inner score, for a BOUNDED intron
		//pass below.  Splice+verify only these 3 -- NEVER every long pair -- so the added alignment
		//cost stays negligible and tRNA calling stays <0.5s/genome.
		//DECOUPLE the intron-span cap from the normal candidate cap: enumerate stops up to intronMax so
		//the intron pass can SEE long spans, but form/verify NORMAL candidates only up to maxlen (the
		//default 120).  This recovers long introns WITHOUT the cost of verifying every long unspliced
		//candidate, so no global maxtrna= is needed and per-genome cost stays at baseline.
		final int minLong=window+MIN_INTRON_SPLICE;
		final int intronMax=Tools.max(maxlen, INTRON_MAX_SPAN);
		final int[] ilStart={-1,-1,-1}, ilStop={-1,-1,-1};
		final float[] ilInner={-1,-1,-1}, ilS=new float[3], ilP=new float[3];

		for(int i=0; i<starts.size; i++){
			final int start=starts.get(i);
			final float sScore=startScores.get(i);
			for(int j=0; j<stops.size; j++){
				final int stop=stops.get(j);
				final int len=stop-start+1;
				if(len<minlen){continue;}
				if(len>intronMax){break;}
				final float pScore=stopScores.get(j);
				final float innerScore=(cumulative[stop]-cumulative[start])/len;
				if(DEBUG && len>100){System.err.println("LONGPAIR\t"+name+"\t"+strand+"\t"+start+"\t"+stop
					+"\t"+len+"\t"+String.format("%.4f", innerScore)+"\t"+(innerScore>=INNER_THRESH?"passInner":"failInner")
					+"\tinnerThr="+String.format("%.4f", INNER_THRESH));}
				if(len>=minLong && innerScore>ilInner[2]){//maintain top-3 by inner score (O(3), no alloc); cheap, no alignment
					final int p=(innerScore>ilInner[0] ? 0 : innerScore>ilInner[1] ? 1 : 2);
					for(int q=2; q>p; q--){ilInner[q]=ilInner[q-1]; ilStart[q]=ilStart[q-1]; ilStop[q]=ilStop[q-1]; ilS[q]=ilS[q-1]; ilP[q]=ilP[q-1];}
					ilInner[p]=innerScore; ilStart[p]=start; ilStop[p]=stop; ilS[p]=sScore; ilP[p]=pScore;
				}
				if(len>maxlen){continue;}//normal candidates only up to the normal cap; long spans go to the intron pass
				if(innerScore>=INNER_THRESH){
					final float a=Tools.max(sScore+0.25f, 0.25f);
					final float b=Tools.max(pScore+0.25f, 0.25f);
					final float c=innerScore*innerScore;
					float d=(window-(2.4f*Tools.absdif(len, window)))*invWindow;
					if(longMode && len>window){d=Tools.max(d, LONG_D_FLOOR);}
					final float score=c*d*(float)Math.sqrt(a*b);
					if(DEBUG && len>100){System.err.println("LONGCAND\t"+name+"\t"+strand+"\t"+start+"\t"+stop
						+"\t"+len+"\t"+String.format("%.3f", score)+"\t"+(score>CANDIDATE_THRESH?"passCand":"failCand")
						+"\tcandThr="+String.format("%.1f", CANDIDATE_THRESH)+"\td="+String.format("%.4f", d));}
					if(score>CANDIDATE_THRESH){
						candidates.add(new float[]{start, stop, score, sScore, pScore, innerScore*len});
					}
				}
			}
		}

		final boolean verify=(trnaLibrary!=null && trnaLibrary.length>0);
		ArrayList<int[]> accepted=new ArrayList<>();

		if(!candidates.isEmpty()){
			candidates.sort((x, y)->Float.compare(y[2], x[2]));

			//Verify DURING greedy selection: a candidate that fails alignment verification
			//does not block the overlapping candidates it outscored, so a valid suppressed
			//tRNA still gets its chance.  Verification cost only rises on the failure path.
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
		}

		//INTRON PASS (evidence recall_gap_subgate) -- BOUNDED: splice+verify only the top-3 intron-suspect
		//long spans tracked above, never every long pair.  Recovers intron-bearing tRNAs whose unspliced
		//span fails INNER_THRESH and/or CANDIDATE_THRESH at generation (the intron depresses the average
		//inner score AND inflates length past the window).  Keeps tRNA calling <0.5s/genome.
		if(verify && INTRON_PASS){
			for(int t=0; t<ilStart.length; t++){
				if(ilStart[t]<0){continue;}
				intronAttempt(name, bases, strand, cumulative, ilStart[t], ilStop[t], ilS[t], ilP[t],
					window, invWindow, accepted, results);
			}
		}

		return results;
	}

	/**
	 * INTRON PASS (single span): given one long intron-suspect (start,stop), locate the intron as the
	 * low-inner segment, splice it out, and if the spliced product clears the generation gates
	 * (INNER_THRESH on the all-exon average; CANDIDATE_THRESH once the length ~= window), verify the
	 * SPLICED product against the mature library.  Accepted calls span the FULL locus (intron included,
	 * matching the reference annotation).  Called only for a bounded number of top spans per region.
	 */
	private void intronAttempt(String name, byte[] bases, int strand, float[] cumulative,
			int start, int stop, float sScore, float pScore, int window, float invWindow,
			ArrayList<int[]> accepted, ArrayList<Orf> results){
		for(int[] p : accepted){if(start<=p[1] && stop>=p[0]){return;}}//already called here
		final int[] intron=findIntronByInner(cumulative, start, stop, window);
		if(intron==null){return;}
		final int gis=intron[0], gie=intron[1];
		final int splicedLen=(gis-start)+(stop-gie);
		if(splicedLen<MIN_TRNA){return;}
		final float exon1=cumulative[gis-1]-cumulative[start];
		final float exon2=cumulative[stop]-cumulative[gie];
		final float splicedInner=(exon1+exon2)/splicedLen;
		if(splicedInner<INNER_THRESH){return;}//not tRNA-like even after splicing
		final float a=Tools.max(sScore+0.25f, 0.25f);
		final float b=Tools.max(pScore+0.25f, 0.25f);
		final float c=splicedInner*splicedInner;
		final float d=(window-(2.4f*Tools.absdif(splicedLen, window)))*invWindow;
		final float score=c*d*(float)Math.sqrt(a*b);
		if(score<=CANDIDATE_THRESH){return;}
		final byte[] product=spliceOut(bases, start, stop, gis, gie);
		if(product==null || product.length<MIN_TRNA){return;}
		if(DEBUG){System.err.println("INTRONTRY\t"+name+"\t"+strand+"\t"+start+"\t"+stop
			+"\tintron="+gis+"-"+gie+"\tsplicedLen="+splicedLen+"\tsplicedInner="+String.format("%.3f", splicedInner)
			+"\tscore="+String.format("%.2f", score));}
		final Orf orf=new Orf(name, start, stop, strand, 0, bases, false, tRNA);
		orf.orfScore=score*GeneCaller.scoreMult[tRNA];
		orf.startScore=sScore;
		orf.stopScore=pScore;
		orf.kmerScore=splicedInner*splicedLen;
		if(verifyIntronOrf(orf, product)){
			if(DEBUG){System.err.println("INTRONHIT\t"+name+"\t"+strand+"\t"+start+"\t"+stop+"\t"+orf.trnaModel);}
			accepted.add(new int[]{orf.start, orf.stop});
			results.add(orf);
		}
	}

	/**
	 * Locate the intron within genomic span [start,stop] as the contiguous segment whose removal
	 * maximizes the spliced product's average inner k-mer score, using the precomputed inner prefix
	 * sums (cumulative).  Intron length is searched around (span - window); each side keeps at least
	 * MIN_EXON_SPLICE of exon.  Returns {intronStart,intronEnd} (genomic, inclusive), or null.
	 */
	private int[] findIntronByInner(float[] cumulative, int start, int stop, int window){
		final int span=stop-start+1;
		final int estIntron=span-window;//mature length ~= window
		if(estIntron<MIN_INTRON_SPLICE){return null;}
		final int lo=start+MIN_EXON_SPLICE, hi=stop-MIN_EXON_SPLICE;
		final int ilMin=Tools.max(MIN_INTRON_SPLICE, estIntron-INTRON_LEN_SLACK);
		final int ilMax=Tools.min(span-2*MIN_EXON_SPLICE, estIntron+INTRON_LEN_SLACK);
		final float total=cumulative[stop]-cumulative[start];
		int bestGis=-1, bestGie=-1; float bestAvg=-1;
		for(int il=ilMin; il<=ilMax; il++){
			final int splicedLen=span-il;
			for(int gis=lo; gis+il-1<=hi; gis++){
				final int gie=gis+il-1;
				final float removed=cumulative[gie]-cumulative[gis-1];
				final float avg=(total-removed)/splicedLen;
				if(avg>bestAvg){bestAvg=avg; bestGis=gis; bestGie=gie;}
			}
		}
		return bestGis<0 ? null : new int[]{bestGis, bestGie};
	}

	/** Extract genomic span [start,stop] with the intron [gis,gie] removed (exons joined), or null. */
	private static byte[] spliceOut(byte[] bases, int start, int stop, int gis, int gie){
		final int e1=gis-start, e2=stop-gie;
		if(e1<1 || e2<1){return null;}
		final byte[] product=new byte[e1+e2];
		System.arraycopy(bases, start, product, 0, e1);
		System.arraycopy(bases, gie+1, product, e1, e2);
		return product;
	}

	/**
	 * Verify an intron-bearing candidate by aligning the SPLICED product to the library.  Uses a
	 * reduced borderline (ID_BORDERLINE_LONG) since a few residual intron bases can still depress
	 * identity even after splicing (guard 3), and a length-similarity floor (guard 1: a short product
	 * cannot match a long model via free terminal deletions).  On pass, annotates the model name; the
	 * tRNA span keeps the full locus (intron included), so no boundary trim is applied.
	 */
	private boolean verifyIntronOrf(Orf orf, byte[] product){
		//Long-kmer pre-filter (Brian): the SPLICED product must carry conserved tRNA 15-mers before align.
		final int khits=trnaKmerHits(product);
		if(khits<MIN_TRNA_KHITS){
			if(DEBUG){System.err.println("KFILTER_REJECT_INTRON\t"+orf.start+"\t"+orf.stop+"\t"+product.length+"\thits="+khits);}
			return false;
		}
		//Cap the shortlist tightly (intron verify runs a few times per region): bounds added alignments.
		final int topN=Tools.min(INTRON_VERIFY_TOPN, INDEX_TOP_N_OVERRIDE>0 ? INDEX_TOP_N_OVERRIDE : INDEX_TOP_N_DEFAULT);
		final int[] shortlist=shortlistByKmer(product, topN);
		float bestId=0; int bestModel=-1, sinceImproved=0;
		final int[] borderlineModels=new int[shortlist.length];
		int borderlineCount=0;
		boolean passed=false;
		for(int idx=0; idx<shortlist.length; idx++){
			final int m=shortlist[idx];
			final float lenSim=Tools.min(product.length, trnaLibrary[m].length)/(float)Tools.max(product.length, trnaLibrary[m].length);
			if(lenSim<MIN_LEN_SIM){continue;}//coverage guard: no zero-cost-deletion match
			final float id=ScrabbleAligner.alignStatic(product, trnaLibrary[m], null);
			alignmentCount++;
			if(id>bestId){bestId=id; bestModel=m; sinceImproved=0;}else{sinceImproved++;}
			if(id>=ID_PASS){passed=true; if(earlyExit && sinceImproved>=earlyExitPatience){break;}}
			if(trnaModels!=null && id>=ID_BORDERLINE_LONG){borderlineModels[borderlineCount++]=m;}
		}
		if(passed){
			if(annotate && bestModel>=0 && modelNames!=null && bestModel<modelNames.length){orf.trnaModel=modelNames[bestModel];}
			return true;
		}else if(bestId>=ID_BORDERLINE_LONG && trnaModels!=null && borderlineCount>0){
			float bestHbm=-999; int bestHbmModel=-1;
			for(int b=0; b<borderlineCount; b++){
				final int idx=borderlineModels[b];
				if(idx<trnaModels.length){
					final float sc=TrnaConsensusBuilder.scoreAgainstModel(product, trnaModels[idx]);
					if(sc>bestHbm){bestHbm=sc; bestHbmModel=idx;}
				}
			}
			if(bestHbm>=HBM_PASS){
				if(annotate && bestHbmModel>=0 && modelNames!=null && bestHbmModel<modelNames.length){orf.trnaModel=modelNames[bestHbmModel];}
				return true;
			}
		}
		return false;
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

		//Long-kmer pre-filter (Brian): never align/trim a candidate lacking conserved tRNA 15-mers.
		final int khits=trnaKmerHits(seq);
		if(khits<MIN_TRNA_KHITS){
			if(DEBUG){System.err.println("KFILTER_REJECT\t"+orf.start+"\t"+orf.stop+"\t"+seq.length+"\thits="+khits);}
			return false;
		}

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

	/*--------------------------------------------------------------*/
	/*----------------       tRNA Scavenger Pass       ----------------*/
	/*--------------------------------------------------------------*/

	public ArrayList<Orf> scavengeTrnas(String name, byte[] bases, int strand, ArrayList<int[]> called){
		ArrayList<Orf> results=new ArrayList<>();
		if((!SCAVENGE && !SCAVENGE_ONLY) || bases==null || bases.length<MIN_TRNA || trnaLibrary==null || ProkObject.trnaKmers==null){return results;}
		int[] hitPositions=findKmerHitPositions(bases);
		if(hitPositions.length==0){return results;}
		ArrayList<int[]> windows=buildCandidateWindows(hitPositions, bases.length);
		windows=collapseByIntersection(windows);
		windows=subtractClaimed(windows, called);
		for(int[] w : windows){
			Orf orf=alignWindow(name, bases, strand, w[0], w[1]);
			if(orf!=null){
				called.add(new int[]{orf.start, orf.stop});
				results.add(orf);
			}
		}
		if(SCAVENGE_PASS2){
			int[] nearHits=findNearbyUnclaimed(hitPositions, called, bases.length);
			if(nearHits.length>0){
				ArrayList<int[]> pass2Windows=buildCandidateWindows(nearHits, bases.length);
				pass2Windows=collapseByIntersection(pass2Windows);
				pass2Windows=subtractClaimed(pass2Windows, called);
				for(int[] w : pass2Windows){
					Orf orf=alignWindow(name, bases, strand, w[0], w[1]);
					if(orf!=null){
						called.add(new int[]{orf.start, orf.stop});
						results.add(orf);
					}
				}
			}
		}
		return results;
	}

	private int[] findKmerHitPositions(byte[] bases){
		final LongHashSet set=ProkObject.trnaKmers;
		final int kLong=ProkObject.kLongTRna;
		if(set==null || kLong<=0 || kLong>31 || bases.length<kLong){return EMPTY;}
		final long kmask=~((-1L)<<(2*kLong));
		final byte[] bton=AminoAcid.baseToNumber;
		IntList hits=new IntList();
		long kmer=0; int len=0;
		for(int i=0; i<bases.length; i++){
			final int x=bton[bases[i]];
			if(x>=0){
				kmer=((kmer<<2)|x)&kmask; len++;
				if(len>=kLong && set.contains(kmer)){hits.add(i-kLong/2);}
			}else{len=0; kmer=0;}
		}
		return hits.toArray();
	}

	private ArrayList<int[]> buildCandidateWindows(int[] hitPositions, int seqLen){
		ArrayList<int[]> windows=new ArrayList<>();
		for(int center : hitPositions){
			int start=Tools.max(0, center-SCAV_PAD);
			int stop=Tools.min(seqLen-1, center+SCAV_PAD+ProkObject.kLongTRna);
			windows.add(new int[]{start, stop});
		}
		return windows;
	}

	private static ArrayList<int[]> collapseByIntersection(ArrayList<int[]> windows){
		if(windows.size()<2){return windows;}
		windows.sort((a, b)->a[0]-b[0]);
		ArrayList<int[]> result=new ArrayList<>();
		int[] current=windows.get(0);
		for(int i=1; i<windows.size(); i++){
			int[] next=windows.get(i);
			int overlap=Tools.min(current[1], next[1])-Tools.max(current[0], next[0]);
			int shorter=Tools.min(current[1]-current[0], next[1]-next[0]);
			if(overlap>0 && overlap>=shorter*0.9f){
				current=new int[]{Tools.max(current[0], next[0]), Tools.min(current[1], next[1])};
			}else{
				if(current[1]-current[0]>=MIN_TRNA){result.add(current);}
				current=next;
			}
		}
		if(current[1]-current[0]>=MIN_TRNA){result.add(current);}
		return result;
	}

	private static ArrayList<int[]> subtractClaimed(ArrayList<int[]> windows, ArrayList<int[]> claimed){
		ArrayList<int[]> result=new ArrayList<>();
		for(int[] w : windows){
			int lo=w[0], hi=w[1];
			for(int[] c : claimed){
				if(lo<=c[1] && hi>=c[0]){
					if(lo>=c[0]){lo=c[1]+1;}
					else if(hi<=c[1]){hi=c[0]-1;}
					else{lo=c[1]+1;}
				}
			}
			if(hi-lo>=MIN_TRNA){result.add(new int[]{lo, hi});}
		}
		return result;
	}

	private static int[] findNearbyUnclaimed(int[] allHits, ArrayList<int[]> claimed, int seqLen){
		IntList result=new IntList();
		for(int pos : allHits){
			boolean inside=false, nearby=false;
			for(int[] c : claimed){
				if(pos>=c[0] && pos<=c[1]){inside=true; break;}
				if(pos>=c[0]-SCAV_NEARBY && pos<=c[1]+SCAV_NEARBY){nearby=true;}
			}
			if(!inside && nearby){result.add(pos);}
		}
		return result.toArray();
	}

	/** SHORTLIST_STATS: emit the accepted model's shared-5-mer count with the query (from the last shortlist),
	 * with qlen and the top hit, so the safe shortlist cutoff can be read off the distribution.
	 * Columns: SLSTAT <sharedHits> <qlen> <maxHitsShared>. */
	private void logShortlistStat(int model, int qlen){
		if(kmerIndex!=null && model>=0){
			System.err.println("SLSTAT\t"+kmerIndex.lastSharedCount(model)+"\t"+qlen+"\t"+kmerIndex.lastMaxShared());
		}
	}

	private Orf alignWindow(String name, byte[] bases, int strand, int wStart, int wStop){
		final int wLen=wStop-wStart+1;
		if(wLen<MIN_TRNA){return null;}
		byte[] seq=Arrays.copyOfRange(bases, wStart, wStop+1);
		final int khits=trnaKmerHits(seq);
		if(khits<MIN_TRNA_KHITS){return null;}
		final int topN=INDEX_TOP_N_OVERRIDE>0 ? INDEX_TOP_N_OVERRIDE : INDEX_TOP_N_DEFAULT;
		int[] shortlist=shortlistByKmer(seq, topN);
		float bestId=0; int bestModel=-1;
		int bestStart=0, bestStop=wLen-1;
		if(wLen>SCAV_QUANTUM_THRESH){
			int[] pos=new int[4];
			for(int j=0; j<shortlist.length; j++){
				int m=shortlist[j];
				float id=QuantumAligner.alignStatic(trnaLibrary[m], seq, pos);
				alignmentCount++;
				if(id>bestId){bestId=id; bestModel=m; bestStart=pos[0]; bestStop=pos[1];}
			}
		}else{
			for(int j=0; j<shortlist.length; j++){
				int m=shortlist[j];
				float id=ScrabbleAligner.alignStatic(seq, trnaLibrary[m], null);
				alignmentCount++;
				if(id>bestId){bestId=id; bestModel=m;}
			}
		}
		if(bestId<ID_BORDERLINE || bestModel<0){return null;}
		final int orfStart=wStart+bestStart;
		final int orfStop=wStart+bestStop;
		if(orfStop-orfStart<MIN_TRNA){return null;}
		Orf orf=new Orf(name, orfStart, orfStop, strand, 0, bases, false, tRNA);
		orf.orfScore=bestId*100;

		//Skip-verify-on-pass (Brian, 2026-08-16): a clear pass (bestId>=ID_PASS) trusts the window
		//alignment's identity+coordinates -- no second alignment.  Only borderline windows (ID_BORDERLINE
		//..ID_PASS) pay a ScrabbleAligner+HBM pass on the trimmed span.  Most tRNAs pass, so the common
		//path adds zero aligns beyond the window search (E.coli 4934->3167 aligns).
		if(bestId>=ID_PASS){
			if(annotate && modelNames!=null && bestModel<modelNames.length){
				orf.trnaModel=modelNames[bestModel];
				annotateAndTrim(orf, bases, bestModel, bestModel);
			}
			if(SHORTLIST_STATS){logShortlistStat(bestModel, seq.length);}
			return orf;
		}else{
			//Borderline (window bestId in [ID_BORDERLINE, ID_PASS)): re-align the trimmed span with
			//ScrabbleAligner.  Accept if any model hits ID_PASS on the trimmed span (restores verifyOrf's
			//recall -- a candidate whose padded window scored low but whose trimmed ORF is a clear match),
			//else fall back to HBM.  Early-exit on the first clear pass.  (Neptune fix, 2026-08-16.)
			byte[] orfSeq=Arrays.copyOfRange(bases, orfStart, orfStop+1);
			float bestReId=0; int bestReModel=-1;
			float bestHbm=-999; int bestHbmModel=-1;
			for(int j=0; j<shortlist.length; j++){
				final int m=shortlist[j];
				final float id=ScrabbleAligner.alignStatic(orfSeq, trnaLibrary[m], null);
				alignmentCount++;
				if(id>bestReId){bestReId=id; bestReModel=m;}
				if(id>=ID_PASS){break;}//clear pass -- stop searching
				if(trnaModels!=null && m<trnaModels.length && id>=ID_BORDERLINE){
					final float hbm=TrnaConsensusBuilder.scoreAgainstModel(orfSeq, trnaModels[m]);
					if(hbm>bestHbm){bestHbm=hbm; bestHbmModel=m;}
				}
			}
			if(bestReId>=ID_PASS && bestReModel>=0){
				if(annotate && modelNames!=null && bestReModel<modelNames.length){
					orf.trnaModel=modelNames[bestReModel];
					annotateAndTrim(orf, bases, bestReModel, bestReModel);
				}
				if(SHORTLIST_STATS){logShortlistStat(bestReModel, seq.length);}
				return orf;
			}else if(bestHbm>=HBM_PASS && bestHbmModel>=0){
				if(annotate && modelNames!=null && bestHbmModel<modelNames.length){
					orf.trnaModel=modelNames[bestHbmModel];
					annotateAndTrim(orf, bases, bestHbmModel, bestHbmModel);
				}
				if(SHORTLIST_STATS){logShortlistStat(bestHbmModel, seq.length);}
				return orf;
			}
		}
		return null;
	}

	/**
	 * Long-kmer pre-filter (Brian): counts how many of seq's forward k-mers (k=kLongTRna, default 15)
	 * are in the conserved tRNA long-kmer set (ProkObject.trnaKmers, from resources/tRNA_15mers.fa).
	 * Candidates below the threshold are never aligned -- this rejects the many non-tRNA candidates on
	 * repetitive/GC-rich genomes cheaply, before any alignment.  Returns MAX_VALUE (filter disabled) if
	 * the set was not loaded, so behavior is unchanged when the resource is absent.
	 */
	private static int trnaKmerHits(byte[] seq){
		if(ProkObject.trnaKmers==null){return Integer.MAX_VALUE;}
		final int k=ProkObject.kLongTRna;
		if(k<=0 || k>31 || seq.length<k){return 0;}
		final long kmask=~((-1L)<<(2*k));
		final byte[] bton=AminoAcid.baseToNumber;
		long kmer=0; int len=0, hits=0;
		for(int i=0; i<seq.length; i++){
			final int x=bton[seq[i]];
			if(x>=0){
				kmer=((kmer<<2)|x)&kmask; len++;
				//Forward-only membership: this REQUIRES resources/tRNA_15mers.fa to be stranded (forward,
				//rcomp=f).  The current file is canonicalized, which drops forward matching to ~87.8% (vs
				//99.7% canonical) and caused the k-filter's net-negative recall -- regenerate the set
				//stranded so the sense kmers match here directly (Brian, 2026-08-16).
				if(len>=k && ProkObject.trnaKmers.contains(kmer)){hits++;}
			}else{len=0; kmer=0;}
		}
		return hits;
	}

	private int[] shortlistByKmer(byte[] seq, int topN){
		if(kmerIndex==null || trnaLibrary==null){
			int[] all=new int[trnaLibrary.length];
			for(int i=0; i<all.length; i++){all[i]=i;}
			return all;
		}
		//Count + counting-top-N select live in TrnaKmerIndex (reused byte[] buffer, no per-query alloc,
		//no boxed sort).  Output order is byte-identical to the old scored-sort (count desc, id asc), so
		//the order-sensitive borderline early-exit in alignWindow is unchanged.
		final int[] result=kmerIndex.shortlist(seq, topN);
		if(REFHIST){System.err.println("REFHIST\t"+kmerIndex.lastDistinctRefs()+"\t"
			+Tools.max(0, seq.length-INDEX_K+1));}
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
	private final TrnaKmerIndex kmerIndex;
	/** Per-model anticodon start position in the consensus, or -1; null when not annotating */
	private final int[] acPositions;

	private static final int MIN_TRNA=40;
	private static final int MAX_TRNA=120;
	/** Raises the candidate length cap (flag maxtrna=); enables the d-floor for
	 * over-window candidates so intron-containing tRNAs can reach verification */
	static int MAX_TRNA_OVERRIDE=-1;
	/** Floor for the length-deviation score term in long mode */
	private static final float LONG_D_FLOOR=0.05f;
	/** Shortlist index k-mer length (flag indexk=).  Longer k => more selective shortlist => fewer aligns,
	 * at the risk of missing divergent tRNAs.  Must be set before the caller is constructed (buildKmerIndex
	 * reads it).  numKmers=1<<(2*INDEX_K): k=5->1024, k=6->4096, k=7->16384.  (Runtime param, Brian 2026-08-17.) */
	static int INDEX_K=7;
	static int INDEX_TOP_N_OVERRIDE=-1;
	static int INDEX_MINHITS_OVERRIDE=-1;
	//Library-search breadth defaults raised to the measured-best scavenger eval config (Brian, 2026-08-16):
	//align up to 60 shortlisted models per candidate (was 20), keyed on >=12 shared index k-mers (was 3).
	//Flags indextopn=/indexminhits= still override.  Provisional -- to be re-swept after the next library.
	private static final int INDEX_TOP_N_DEFAULT=60;
	private static final int INDEX_MINHITS_DEFAULT=12;
	/** Adaptive shortlist cutoff (flag adaptiveminhits=): use min(ADAPTIVE_MINHITS_CAP, qlen/2,
	 * 0.75*maxHitsShared) instead of the flat indexminhits.  Experimental (Brian, 2026-08-17). */
	//SHIPPED default (Brian, 2026-08-17): adaptive shortlist cutoff ON at (floor=11, topfrac=0.48, qfrac=0.072)
	//with indexk=7.  Measured (203-bench, taxonomy=f): recall 94.6% / prec 98.1% / 334K aligns (~4.7x fewer
	//than the pre-adaptive 1.57M; +5.7pp recall over the 88.9% baseline).  Chosen over k=8 (Pareto-dominated:
	//k=8 needs ~37% more aligns for equal recall; its lighter index-ops don't offset that).
	static boolean ADAPTIVE_MINHITS=true;
	/** Adaptive-cutoff params (Brian, 2026-08-17): minHits = ceil(max(ADAPT_FLOOR, ADAPT_TOPFRAC*maxShared,
	 * ADAPT_QFRAC*queryKmers)).  Flags: adaptfloor=/adapttopfrac=/adaptqfrac=. */
	static float ADAPT_FLOOR=11;
	static float ADAPT_TOPFRAC=0.48f;
	static float ADAPT_QFRAC=0.072f;
	/** Evidence gathering (flag shortliststats=): for each ACCEPTED tRNA, print the winning model's shared
	 * 5-mer count with the query -- absolute, and enough to compute /qlen and /maxHitsShared -- so the safe
	 * shortlist cutoff can be set from the distribution instead of a guessed formula (Brian, 2026-08-17). */
	static boolean SHORTLIST_STATS=false;
	/** Evidence (flag refhist=): per shortlist query, log how many library refs share >=1 index-kmer
	 * (REFHIST <distinctRefs> <queryKmers>) -- direct-array vs sparse IntHashMap counter decision. */
	static boolean REFHIST=false;
	private static final int[] EMPTY=new int[0];
	//Per-INSTANCE (not static) so CallGenes flags that set GeneCaller.cutoffN[tRNA] before this caller is
	//constructed take effect.  cutoff1=region open, cutoff5=avg inner, cutoff2=composite candidate.
	private final float REGION_THRESH=GeneCaller.cutoff1[tRNA];
	private final float INNER_THRESH=GeneCaller.cutoff5[tRNA];
	private final float CANDIDATE_THRESH=GeneCaller.cutoff2[tRNA];
	static float ID_PASS=0.75f;
	static float ID_BORDERLINE=0.65f;
	static float HBM_PASS=0.75f;
	/** Reduced identity borderline for SPLICED intron-bearing products: a few residual intron bases
	 * depress identity even after splicing (guard 3). Flag: idborderlinelong= */
	static float ID_BORDERLINE_LONG=0.55f;
	/** Min length-similarity (spliced product vs model) to accept an intron verification -- guard 1:
	 * a short product cannot match a long model via free terminal deletions. */
	static final float MIN_LEN_SIM=0.70f;
	/** Intron-splice geometry: min intron length, min exon on each side, and +/- search slack. */
	static final int MIN_INTRON_SPLICE=10;
	static final int MIN_EXON_SPLICE=15;
	static final int INTRON_LEN_SLACK=20;
	/** Max models aligned per intron verification -- caps the bounded intron pass's alignment cost. */
	static final int INTRON_VERIFY_TOPN=12;
	/** Master on/off for the two-half intron pass (flag trnaintron=). */
	static boolean INTRON_PASS=true;
	/** Intron-span tracking cap: the intron pass considers spans up to this even when the normal
	 * candidate cap (MAX_TRNA) is lower, so long introns are found without a global maxtrna=. */
	static final int INTRON_MAX_SPAN=260;
	/** Long-kmer pre-filter threshold (flag mintrnakhits=): min conserved tRNA 15-mers a candidate must
	 * carry before it is aligned.  DEFAULT 1 = ON.  Requires the STRANDED (forward, rcomp=f) set
	 * resources/tRNA_15mers.fa -- forward-only trnaKmerHits then matches the sense kmers directly at
	 * 99.74% coverage.  (An earlier CANONICAL set gave only 87.8% forward -> -10.7pp recall; the fix was
	 * to regenerate the set stranded, not to canonicalize the lookup.)  Measured on the 203-genome
	 * benchmark: recall 88.8% (baseline 88.9, neutral), precision +0.2pp, and 365x fewer tRNA alignments
	 * on outlier genomes (Pyrodictium 9.2M -> 25.3K) -- keeps tRNA calling <0.5s/genome.  See
	 * results/recall_gap_kfilter_intron_20260816.md. */
	static int MIN_TRNA_KHITS=1;
	/** tRNA scavenger pass (Neptune): recover tRNAs at conserved-kmer-hit positions the PGM-based candidate
	 * generation missed.  SCAVENGE=augment the normal call; SCAVENGE_ONLY=replace it entirely.  Flags:
	 * trnascavenge=/scavenge=, trnascavengeonly=/scavengeonly=.  DEFAULT: SCAVENGE_ONLY=true -- the scavenger
	 * alone measured 98% precision / 95% recall / >80% both-ends-exact on the benchmark (Neptune), matching or
	 * beating the PGM candidate generator, so it replaces the normal tRNA path by default (Brian, 2026-08-16).
	 * NOTE: under SCAVENGE_ONLY the normal callTrnas path (incl. the INTRON_PASS) is bypassed; the measured
	 * recall already reflects that.  Provisional -- re-grade after the next (2.5M) library before finalizing. */
	static boolean SCAVENGE=false;
	static boolean SCAVENGE_ONLY=true;
	static boolean SCAVENGE_PASS2=true;
	static final int SCAV_PAD=83;
	static final int SCAV_QUANTUM_THRESH=120;
	static final int SCAV_NEARBY=200;
	static boolean earlyExit=true;
	//Raised 10->20 to the measured-best scavenger eval config (Brian, 2026-08-16); flag patience= overrides.
	static int earlyExitPatience=20;
	static long alignmentCount=0;
	/** Debug trace (flag trnadebug=): dumps findRegions regions + long-candidate INNER_THRESH
	 * decisions to stderr to pin which sub-gate rejects intron-bearing candidates. Zero-cost off. */
	public static boolean DEBUG=false;
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
