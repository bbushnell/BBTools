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
import structures.IntList;

/**
 * Generic conserved-ncRNA scavenger: the family-agnostic core of TrnaCaller's
 * scavenger pass (conserved-kmer seed -> candidate window -> kmer-index
 * shortlist -> ScrabbleAligner/QuantumAligner verification -> HBM fallback),
 * factored out per Noire's read of TrnaCaller (2026-08-23) so any conserved
 * ncRNA family with a consensus library, HBM model set, and conserved
 * long-kmer set can be called the same way -- without tRNA's PGM-based
 * candidate generation (callTrnas/scanInner/findRegions/extractTrnas),
 * anticodon logic, acceptor-stem trimming, or intron handling, all of which
 * stay tRNA-only in TrnaCaller.
 *
 * <p>Boundary trimming here is deliberately simplified from TrnaCaller's:
 * an extended window is aligned to the winning model with traceback and the
 * ORF is snapped to the alignment's own extent (rStart/rStop) -- no
 * acceptor-stem search, no anticodon extraction. That structural refinement
 * is tRNA-specific and stays in TrnaCaller.
 *
 * <p>Per-family tunables (window pad, min length, kmer set/k) are
 * constructor parameters since they genuinely differ per family (tRNA
 * windowPad=83; RNase P/SRP need their own). The remaining knobs (index k,
 * identity thresholds, etc.) are mutable instance fields defaulted to
 * TrnaCaller's measured-best tRNA values, per Noire's "start with tRNA's
 * defaults, tune per family later" -- not yet wired to per-family values.
 * @author Neptune, Noire, Brian Bushnell
 */
public class NcrnaScavenger {

	public NcrnaScavenger(byte[][] library_, BaseGraph[] models_, String[] modelNames_,
			LongHashSet kmerSet_, int kLong_, int minLen_, int windowPad_){
		this(library_, models_, modelNames_, kmerSet_, kLong_, minLen_, windowPad_,
				7, 60, true, 11f, 0.48f, 0.072f, 12);
	}

	public NcrnaScavenger(byte[][] library_, BaseGraph[] models_, String[] modelNames_,
			LongHashSet kmerSet_, int kLong_, int minLen_, int windowPad_,
			int indexK_, int indexTopN_, boolean adaptive_,
			float adaptFloor_, float adaptTopFrac_, float adaptQFrac_, int fixedMinHits_){
		library=library_;
		models=models_;
		modelNames=modelNames_;
		annotate=(modelNames!=null);
		kmerSet=kmerSet_;
		kLong=kLong_;
		minLen=minLen_;
		windowPad=windowPad_;
		indexK=indexK_;
		indexTopN=indexTopN_;
		adaptiveMinHits=adaptive_;
		adaptFloor=adaptFloor_;
		adaptTopFrac=adaptTopFrac_;
		adaptQFrac=adaptQFrac_;
		indexMinHitsDefault=fixedMinHits_;
		kmerIndex=(library!=null ? new TrnaKmerIndex(library, indexK, adaptiveMinHits,
			adaptFloor, adaptTopFrac, adaptQFrac, indexMinHitsDefault) : null);
	}

	public long alignmentCount(){return alignmentCount;}

	/*--------------------------------------------------------------*/
	/*----------------       Scavenger Pass       ----------------*/
	/*--------------------------------------------------------------*/

	public ArrayList<Orf> scavenge(String name, byte[] bases, int strand, ArrayList<int[]> called){
		ArrayList<Orf> results=new ArrayList<>();
		if(bases==null || bases.length<minLen || library==null || kmerSet==null){return results;}
		int[] hitPositions=findKmerHitPositions(bases);
		if(DEBUG){System.err.println("DEBUG scavenge name="+name+" strand="+strand+" bases.length="+bases.length
			+" hitPositions="+Arrays.toString(hitPositions));}
		if(hitPositions.length==0){return results;}
		ArrayList<int[]> windows=buildCandidateWindows(hitPositions, bases.length);
		if(DEBUG){System.err.println("DEBUG windows before collapse: "+dumpWindows(windows));}
		windows=collapseByIntersection(windows);
		if(DEBUG){System.err.println("DEBUG windows after collapse: "+dumpWindows(windows));}
		windows=subtractClaimed(windows, called);
		if(DEBUG){System.err.println("DEBUG windows after subtractClaimed: "+dumpWindows(windows));}
		for(int[] w : windows){
			Orf orf=alignWindow(name, bases, strand, w[0], w[1]);
			if(DEBUG){System.err.println("DEBUG alignWindow("+w[0]+","+w[1]+") -> "+(orf==null ? "null" : (orf.start+"-"+orf.stop+" score="+orf.orfScore)));}
			if(orf!=null){
				called.add(new int[]{orf.start, orf.stop});
				results.add(orf);
			}
		}
		if(scavengePass2){
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

	/** Temporary diagnostic flag (Neptune, 2026-08-23): traces kmer hits, candidate window
	 * construction, and per-model alignment scores for FN root-cause diagnosis. Default false,
	 * zero cost when off. Not wired to a command-line flag -- toggle+recompile for one-off use. */
	static boolean DEBUG=Boolean.getBoolean("ncrna.debug");

	private static String dumpWindows(ArrayList<int[]> windows){
		StringBuilder sb=new StringBuilder("[");
		for(int[] w : windows){sb.append("(").append(w[0]).append(",").append(w[1]).append(") ");}
		sb.append("]");
		return sb.toString();
	}

	private int[] findKmerHitPositions(byte[] bases){
		if(kmerSet==null || kLong<=0 || kLong>31 || bases.length<kLong){return EMPTY;}
		final long kmask=~((-1L)<<(2*kLong));
		final byte[] bton=AminoAcid.baseToNumber;
		IntList hits=new IntList();
		long kmer=0; int len=0;
		for(int i=0; i<bases.length; i++){
			final int x=bton[bases[i]];
			if(x>=0){
				kmer=((kmer<<2)|x)&kmask; len++;
				if(len>=kLong && kmerSet.contains(kmer)){hits.add(i-kLong/2);}
			}else{len=0; kmer=0;}
		}
		return hits.toArray();
	}

	private ArrayList<int[]> buildCandidateWindows(int[] hitPositions, int seqLen){
		ArrayList<int[]> windows=new ArrayList<>();
		for(int center : hitPositions){
			int start=Tools.max(0, center-windowPad);
			int stop=Tools.min(seqLen-1, center+windowPad+kLong);
			windows.add(new int[]{start, stop});
		}
		return windows;
	}

	/** No explicit upper bound on the merged window size -- matches TrnaCaller's own
	 * collapseByIntersection, which has never capped this either (SCAV_QUANTUM_THRESH
	 * only switches the aligner for large windows, it doesn't reject them). */
	private ArrayList<int[]> collapseByIntersection(ArrayList<int[]> windows){
		if(windows.size()<2){return windows;}
		windows.sort((a, b)->a[0]-b[0]);
		ArrayList<int[]> result=new ArrayList<>();
		int[] current=windows.get(0);
		for(int i=1; i<windows.size(); i++){
			int[] next=windows.get(i);
			int overlap=Tools.min(current[1], next[1])-Tools.max(current[0], next[0]);
			int shorter=Tools.min(current[1]-current[0], next[1]-next[0]);
			if(overlap>0 && overlap>=shorter*collapseFrac){
				current=new int[]{Tools.max(current[0], next[0]), Tools.min(current[1], next[1])};
			}else{
				if(current[1]-current[0]>=minLen){result.add(current);}
				current=next;
			}
		}
		if(current[1]-current[0]>=minLen){result.add(current);}
		return result;
	}

	/** Subtracts every claimed interval from every window, keeping BOTH surviving remainders when a
	 * claim lands strictly inside a window (2026-08-27 fix, mirrors the identical fix in
	 * TrnaCaller.subtractClaimed -- see that javadoc for the root-cause trace; found via the tRNA path
	 * first, applied here because this class shares the same buggy pattern). Carries a list of
	 * surviving segments through every claim so a claim strictly inside a window correctly produces
	 * TWO output windows (left + right), not zero or one. */
	private ArrayList<int[]> subtractClaimed(ArrayList<int[]> windows, ArrayList<int[]> claimed){
		ArrayList<int[]> result=new ArrayList<>();
		for(int[] w : windows){
			ArrayList<int[]> segments=new ArrayList<>();
			segments.add(new int[]{w[0], w[1]});
			for(int[] c : claimed){
				ArrayList<int[]> next=new ArrayList<>();
				for(int[] seg : segments){
					final int lo=seg[0], hi=seg[1];
					if(hi<c[0] || lo>c[1]){
						next.add(seg);//no overlap with this claim -- unchanged
						continue;
					}
					if(lo<c[0]){next.add(new int[]{lo, c[0]-1});}//left remainder survives
					if(hi>c[1]){next.add(new int[]{c[1]+1, hi});}//right remainder survives
					//neither branch fires -> claim fully covers this segment, it's consumed
				}
				segments=next;
			}
			//TODO: Probable bug (pre-existing, not fixed here -- see the matching TODO in
			//TrnaCaller.subtractClaimed) -- hi-lo>=minLen treats an inclusive [lo,hi] range as if it
			//had hi-lo bases, when it actually has hi-lo+1. Not changed here; flagging only.
			for(int[] seg : segments){
				if(seg[1]-seg[0]>=minLen){result.add(seg);}
			}
		}
		return result;
	}

	private int[] findNearbyUnclaimed(int[] allHits, ArrayList<int[]> claimed, int seqLen){
		IntList result=new IntList();
		for(int pos : allHits){
			boolean inside=false, nearby=false;
			for(int[] c : claimed){
				if(pos>=c[0] && pos<=c[1]){inside=true; break;}
				if(pos>=c[0]-nearbyPad && pos<=c[1]+nearbyPad){nearby=true;}
			}
			if(!inside && nearby){result.add(pos);}
		}
		return result.toArray();
	}

	private Orf alignWindow(String name, byte[] bases, int strand, int wStart, int wStop){
		final int wLen=wStop-wStart+1;
		if(wLen<minLen){return null;}
		byte[] seq=Arrays.copyOfRange(bases, wStart, wStop+1);
		final int khits=kmerHits(seq);
		if(DEBUG){System.err.println("DEBUG alignWindow wLen="+wLen+" khits="+khits+" minKmerHits="+minKmerHits
			+" quantumThresh="+quantumThresh+" usingQuantum="+(wLen>quantumThresh));}
		if(khits<minKmerHits){return null;}
		int[] shortlist=shortlistByKmer(seq, indexTopN);
		if(DEBUG){System.err.println("DEBUG shortlist.length="+shortlist.length+" shortlist="+Arrays.toString(shortlist));}
		float bestId=0; int bestModel=-1;
		int bestStart=0, bestStop=wLen-1;
		if(wLen>quantumThresh){
			int[] pos=new int[4];
			for(int j=0; j<shortlist.length; j++){
				int m=shortlist[j];
				float id=QuantumAligner.alignStatic(library[m], seq, pos);
				alignmentCount++;
				if(DEBUG){System.err.println("DEBUG   quantum model="+m+" modelLen="+library[m].length+" id="+id+" pos="+Arrays.toString(pos));}
				if(id>bestId){bestId=id; bestModel=m; bestStart=pos[0]; bestStop=pos[1];}
			}
		}else{
			for(int j=0; j<shortlist.length; j++){
				int m=shortlist[j];
				float id=ScrabbleAligner.alignStatic(seq, library[m], null);
				alignmentCount++;
				if(DEBUG){System.err.println("DEBUG   scrabble model="+m+" modelLen="+library[m].length+" id="+id);}
				if(id>bestId){bestId=id; bestModel=m;}
			}
		}
		if(DEBUG){System.err.println("DEBUG best: bestId="+bestId+" bestModel="+bestModel+" idBorderline="+idBorderline+" idPass="+idPass);}
		if(bestId<idBorderline || bestModel<0){return null;}
		final int orfStart=wStart+bestStart;
		final int orfStop=wStart+bestStop;
		if(orfStop-orfStart<minLen){return null;}
		Orf orf=new Orf(name, orfStart, orfStop, strand, 0, bases, false, ProkObject.RNA);
		orf.orfScore=bestId*100;

		if(bestId>=idPass){
			if(annotate && modelNames!=null && bestModel<modelNames.length){
				orf.trnaModel=modelNames[bestModel];
				trimToAlignmentExtent(orf, bases, bestModel);
			}
			return orf;
		}else{
			byte[] orfSeq=Arrays.copyOfRange(bases, orfStart, orfStop+1);
			if(DEBUG){System.err.println("DEBUG rescue orfStart="+orfStart+" orfStop="+orfStop
				+" orfSeq.length="+orfSeq.length+" (from Quantum pos bestStart="+bestStart+" bestStop="+bestStop+")");}
			float bestReId=0; int bestReModel=-1;
			float bestHbm=-999; int bestHbmModel=-1;
			for(int j=0; j<shortlist.length; j++){
				final int m=shortlist[j];
				final float id=ScrabbleAligner.alignStatic(orfSeq, library[m], null);
				alignmentCount++;
				if(DEBUG){System.err.println("DEBUG   rescue-scrabble model="+m+" modelLen="+library[m].length+" id="+id);}
				if(id>bestReId){bestReId=id; bestReModel=m;}
				if(id>=idPass){break;}
				if(models!=null && m<models.length && id>=idBorderline){
					final float hbm=TrnaConsensusBuilder.scoreAgainstModel(orfSeq, models[m]);
					if(DEBUG){System.err.println("DEBUG   rescue-hbm model="+m+" hbm="+hbm);}
					if(hbm>bestHbm){bestHbm=hbm; bestHbmModel=m;}
				}
			}
			if(DEBUG){System.err.println("DEBUG rescue result: bestReId="+bestReId+" bestReModel="+bestReModel
				+" bestHbm="+bestHbm+" bestHbmModel="+bestHbmModel+" idPass="+idPass+" hbmPass="+hbmPass);}
			if(bestReId>=idPass && bestReModel>=0){
				if(annotate && modelNames!=null && bestReModel<modelNames.length){
					orf.trnaModel=modelNames[bestReModel];
					trimToAlignmentExtent(orf, bases, bestReModel);
				}
				return orf;
			}else if(bestHbm>=hbmPass && bestHbmModel>=0){
				if(annotate && modelNames!=null && bestHbmModel<modelNames.length){
					orf.trnaModel=modelNames[bestHbmModel];
					trimToAlignmentExtent(orf, bases, bestHbmModel);
				}
				return orf;
			}
		}
		return null;
	}

	/**
	 * Generic boundary trim: aligns an extended window around the candidate to
	 * the winning model with traceback, and snaps the ORF boundaries to the
	 * alignment's own extent (rStart/rStop) directly -- no acceptor-stem
	 * search, no anticodon extraction (both tRNA-specific; stay in
	 * TrnaCaller). Mirrors TrnaCaller.trimOrf's "window within consensus
	 * span" guard.
	 */
	private void trimToAlignmentExtent(Orf orf, byte[] bases, int model){
		final int xFrom=Tools.max(0, orf.start-trimExt);
		final int xTo=Tools.min(bases.length-1, orf.stop+trimExt);
		byte[] seqX=Arrays.copyOfRange(bases, xFrom, xTo+1);
		final byte[] cons=library[model];
		AlignmentStats stats=new AlignmentStats(true);
		stats.doTrace=true;
		if(seqX.length>=cons.length){ScrabbleAligner.alignAndTraceStatic(cons, seqX, stats);}
		else{ScrabbleAligner.alignAndTraceStatic(seqX, cons, stats);}
		alignmentCount++;
		if(stats.matchString==null || seqX.length<cons.length){return;}
		final int rStart=stats.rStart, rStop=stats.rStop;
		if(rStart<0 || rStop<=rStart || rStop>=seqX.length){return;}
		orf.start=xFrom+rStart;
		orf.stop=xFrom+rStop;
	}

	private int kmerHits(byte[] seq){
		if(kmerSet==null){return Integer.MAX_VALUE;}
		if(kLong<=0 || kLong>31 || seq.length<kLong){return 0;}
		final long kmask=~((-1L)<<(2*kLong));
		final byte[] bton=AminoAcid.baseToNumber;
		long kmer=0; int len=0, hits=0;
		for(int i=0; i<seq.length; i++){
			final int x=bton[seq[i]];
			if(x>=0){
				kmer=((kmer<<2)|x)&kmask; len++;
				if(len>=kLong && kmerSet.contains(kmer)){hits++;}
			}else{len=0; kmer=0;}
		}
		return hits;
	}

	private int[] shortlistByKmer(byte[] seq, int topN){
		if(kmerIndex==null || library==null){
			int[] all=new int[library.length];
			for(int i=0; i<all.length; i++){all[i]=i;}
			return all;
		}
		return kmerIndex.shortlist(seq, topN);
	}

	/*--------------------------------------------------------------*/

	private final byte[][] library;
	private final BaseGraph[] models;
	private final String[] modelNames;
	private final boolean annotate;
	private final LongHashSet kmerSet;
	private final int kLong;
	private final int minLen;
	private final int windowPad;
	private final TrnaKmerIndex kmerIndex;

	private static final int[] EMPTY=new int[0];

	//Per-family tunables not yet wired to per-family values (Noire, 2026-08-23: "start
	//with tRNA's defaults, tune per family later") -- default to TrnaCaller's measured-best
	//tRNA config. Mutable instance fields (not static): unlike windowPad/minLen/kLong, these
	//don't yet have measured per-family values, but must stay per-instance since coexisting
	//scavengers for different families will eventually need different ones.
	int indexK=7;
	int indexTopN=60;
	int indexMinHitsDefault=12;
	boolean adaptiveMinHits=true;
	float adaptFloor=11;
	float adaptTopFrac=0.48f;
	float adaptQFrac=0.072f;
	int minKmerHits=1;
	float idPass=0.75f;
	float idBorderline=0.65f;
	float hbmPass=0.75f;
	int quantumThresh=120;
	int nearbyPad=200;
	float collapseFrac=0.9f;
	int trimExt=10;
	boolean scavengePass2=true;

	private long alignmentCount=0;
}
