package prok;

import java.util.ArrayList;
import java.util.Arrays;

import consensus.BaseGraph;
import dna.AminoAcid;
import idaligner.AlignmentStats;
import idaligner.QuantumAligner;
import idaligner.ScrabbleAligner;
import map.LongHashSet;
import ml.CellNet;
import shared.KillSwitch;
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
		this(library_, models_, modelNames_, kmerSet_, kLong_, minLen_, windowPad_,
				indexK_, indexTopN_, adaptive_, adaptFloor_, adaptTopFrac_, adaptQFrac_, fixedMinHits_,
				0f, 20f, 0.75f, 0.65f);
	}

	/** Forward-ported from Noire's ncRNA-family-loading tree (2026-08-28, C3 merge) -- see
	 * NcrnaFamily's matching constructor javadoc for the scoreA/scoreB/idPass/idBorderline
	 * rationale (Noire's #1 recall lever, +11.5pp rnasep). idPass/idBorderline are now
	 * FINAL, set here -- previously mutable instance fields defaulted to 0.75f/0.65f
	 * (unchanged as the fallback defaults for this overload's own delegation below). */
	public NcrnaScavenger(byte[][] library_, BaseGraph[] models_, String[] modelNames_,
			LongHashSet kmerSet_, int kLong_, int minLen_, int windowPad_,
			int indexK_, int indexTopN_, boolean adaptive_,
			float adaptFloor_, float adaptTopFrac_, float adaptQFrac_, int fixedMinHits_,
			float scoreA_, float scoreB_, float idPass_, float idBorderline_){
		this(library_, models_, modelNames_, kmerSet_, kLong_, minLen_, windowPad_,
				indexK_, indexTopN_, adaptive_, adaptFloor_, adaptTopFrac_, adaptQFrac_, fixedMinHits_,
				scoreA_, scoreB_, idPass_, idBorderline_,
				null, null, null, null, -1, -1, -1, -1, 0f);
	}

	/** Full constructor, adding C3's boundary-precision-NN resources (Noire's spec,
	 * plans/c3_ncrnaboundaryscorer_spec.md; G11, 2026-08-28). boundary5NetTemplate==null
	 * means OFF -- structural, matches NcrnaFamily's own default. Per-instance CLONES its
	 * own net copies from the shared read-only templates (mirrors TrnaCaller's
	 * boundary5Net/boundary3Net constructor-clone pattern exactly, same reason: CellNet.
	 * feedForward mutates per-Cell state, so concurrent use of ONE net object across
	 * per-thread NcrnaScavenger instances would corrupt it). meanLen is asserted >0
	 * whenever a template is given -- there is no safe shared default across families of
	 * very different length (rnasep ~380bp vs srp_small ~95bp), so a caller passing a real
	 * net without a real meanLen is a construction-time bug, not a runtime one -- fail
	 * loud immediately rather than silently score every candidate with lengthRatio=0. */
	public NcrnaScavenger(byte[][] library_, BaseGraph[] models_, String[] modelNames_,
			LongHashSet kmerSet_, int kLong_, int minLen_, int windowPad_,
			int indexK_, int indexTopN_, boolean adaptive_,
			float adaptFloor_, float adaptTopFrac_, float adaptQFrac_, int fixedMinHits_,
			float scoreA_, float scoreB_, float idPass_, float idBorderline_,
			CellNet boundary5NetTemplate, CellNet boundary3NetTemplate,
			TrnaBoundaryFeatures.NinemerTable boundaryStartTable_, TrnaBoundaryFeatures.NinemerTable boundaryStopTable_,
			int boundaryStartInside_, int boundaryStartOutside_, int boundaryStopInside_, int boundaryStopOutside_,
			float boundaryMeanLen_){
		this(library_, models_, modelNames_, kmerSet_, kLong_, minLen_, windowPad_,
				indexK_, indexTopN_, adaptive_, adaptFloor_, adaptTopFrac_, adaptQFrac_, fixedMinHits_,
				scoreA_, scoreB_, idPass_, idBorderline_, boundary5NetTemplate, boundary3NetTemplate,
				boundaryStartTable_, boundaryStopTable_, boundaryStartInside_, boundaryStartOutside_,
				boundaryStopInside_, boundaryStopOutside_, boundaryMeanLen_,
				NcrnaFamily.LEGACY_START_OFFSETS, NcrnaFamily.LEGACY_STOP_OFFSETS);
	}

	public NcrnaScavenger(byte[][] library_, BaseGraph[] models_, String[] modelNames_,
			LongHashSet kmerSet_, int kLong_, int minLen_, int windowPad_,
			int indexK_, int indexTopN_, boolean adaptive_,
			float adaptFloor_, float adaptTopFrac_, float adaptQFrac_, int fixedMinHits_,
			float scoreA_, float scoreB_, float idPass_, float idBorderline_,
			CellNet boundary5NetTemplate, CellNet boundary3NetTemplate,
			TrnaBoundaryFeatures.NinemerTable boundaryStartTable_, TrnaBoundaryFeatures.NinemerTable boundaryStopTable_,
			int boundaryStartInside_, int boundaryStartOutside_, int boundaryStopInside_, int boundaryStopOutside_,
			float boundaryMeanLen_, int[] boundaryStartOffsets_, int[] boundaryStopOffsets_){
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
		scoreA=scoreA_;
		scoreB=scoreB_;
		idPass=idPass_;
		idBorderline=idBorderline_;
		kmerIndex=(library!=null ? new TrnaKmerIndex(library, indexK, adaptiveMinHits,
			adaptFloor, adaptTopFrac, adaptQFrac, indexMinHitsDefault) : null);
		boundary5Net=(boundary5NetTemplate!=null ? boundary5NetTemplate.copy(false) : null);
		boundary3Net=(boundary3NetTemplate!=null ? boundary3NetTemplate.copy(false) : null);
		boundaryStartTable=boundaryStartTable_;
		boundaryStopTable=boundaryStopTable_;
		boundaryStartInside=boundaryStartInside_;
		boundaryStartOutside=boundaryStartOutside_;
		boundaryStopInside=boundaryStopInside_;
		boundaryStopOutside=boundaryStopOutside_;
		boundaryMeanLen=boundaryMeanLen_;
		boundaryStartOffsets=boundaryStartOffsets_.clone();
		boundaryStopOffsets=boundaryStopOffsets_.clone();
		assert(boundary5Net==null || boundaryMeanLen>0) : KillSwitch.assertDie(
			"NcrnaScavenger built with a boundary-precision net but boundaryMeanLen="+boundaryMeanLen_
			+" (must be >0) -- lengthRatio=(e-s+1)/meanLen would silently divide by a bogus value for "
			+"every scored candidate. Fix the caller (NcrnaFamily/CallGenes.loadNcrnaResources).");
	}

	public long alignmentCount(){return alignmentCount;}

	/*--------------------------------------------------------------*/
	/*----------------       Scavenger Pass       ----------------*/
	/*--------------------------------------------------------------*/

	public ArrayList<Orf> scavenge(String name, byte[] bases, int strand, ArrayList<int[]> called){
		ArrayList<Orf> results=new ArrayList<>();
		if(bases==null || bases.length<minLen || library==null || kmerSet==null){return results;}
		int[] hitPositions=findKmerHitPositions(bases);
		if(workloadSink!=null){workloadSink.seedHits(name, strand, Arrays.copyOf(hitPositions, hitPositions.length));}
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
			if(workloadSink!=null){workloadSink.scheduledWindow(name, strand, 1, w[0], w[1]);}
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
					if(workloadSink!=null){workloadSink.scheduledWindow(name, strand, 2, w[0], w[1]);}
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
		//Forward-ported from Noire's tree (2026-08-28, C3 merge): inclusive-length check (+1)
		//alongside the scoreA/scoreB formula below -- the OLD `orfStop-orfStart<minLen` undercounted
		//an inclusive [orfStart,orfStop] span by one base, the same off-by-one class flagged (but
		//deliberately NOT fixed) in subtractClaimed's TODO comment; this one Noire did fix, as part
		//of the same scoreA/scoreB commit that introduced orfLen.
		if(orfStop-orfStart+1<minLen){return null;}
		Orf orf=new Orf(name, orfStart, orfStop, strand, 0, bases, false, ProkObject.RNA);
		final int orfLen=orfStop-orfStart+1;
		//Forward-ported from Noire's tree (2026-08-28, C3 merge): scoreA/scoreB per-family score
		//formula replaces the flat bestId*100 -- Noire's #1 recall lever (+11.5pp rnasep).
		//NcrnaScavenger-only; TrnaCaller's tRNA orfScore (bestId*100) is untouched.
		orf.orfScore=scoreA+scoreB*orfLen*bestId*bestId;

		if(bestId>=idPass){
			if(annotate && modelNames!=null && bestModel<modelNames.length){
				orf.trnaModel=modelNames[bestModel];
				trimToAlignmentExtent(orf, bases, bestModel, wStart, wStop);
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
					trimToAlignmentExtent(orf, bases, bestReModel, wStart, wStop);
				}
				return orf;
			}else if(bestHbm>=hbmPass && bestHbmModel>=0){
				if(annotate && modelNames!=null && bestHbmModel<modelNames.length){
					orf.trnaModel=modelNames[bestHbmModel];
					trimToAlignmentExtent(orf, bases, bestHbmModel, wStart, wStop);
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
	/** Package-visible (was private) so NcrnaBoundaryInstrumentSinkTest can directly construct
	 * the trim-no-op case (extended window shorter than the model) without needing to coax it
	 * out of the full scavenge() pipeline's alignment dynamics -- deterministic and fast versus
	 * fragile fixture engineering for the same real code path. */
	void trimToAlignmentExtent(Orf orf, byte[] bases, int model, int wStart, int wStop){
		final int xFrom=Tools.max(0, orf.start-trimExt);
		final int xTo=Tools.min(bases.length-1, orf.stop+trimExt);
		byte[] seqX=Arrays.copyOfRange(bases, xFrom, xTo+1);
		final byte[] cons=library[model];
		AlignmentStats stats=new AlignmentStats(true);
		stats.doTrace=true;
		if(seqX.length>=cons.length){ScrabbleAligner.alignAndTraceStatic(cons, seqX, stats);}
		else{ScrabbleAligner.alignAndTraceStatic(seqX, cons, stats);}
		alignmentCount++;
		//Citan, 2026-08-28: the original code `return`ed here on either guard, which ALSO
		//skipped instrumentation capture for these loci -- silently undercounting the accepted
		//denominator for exactly the accepted-but-unrefined cases. Restructured to a boolean so
		//capture always fires once per accepted locus (see below). trimSucceeded is the
		//STRUCTURAL eligibility signal for downstream bootstrap work: true means orf.start/
		//orf.stop below are a real post-trim position (usable, pending the driver's own
		//sweep-reachability check against truth) -- false means they are the untouched raw
		//alignWindow span and must not be used as if they were a trim result.
		boolean trimSucceeded=false;
		if(stats.matchString!=null && seqX.length>=cons.length){
			final int rStart=stats.rStart, rStop=stats.rStop;
			if(rStart>=0 && rStop>rStart && rStop<seqX.length){
				orf.start=xFrom+rStart;
				orf.stop=xFrom+rStop;
				trimSucceeded=true;
			}
		}
		//C3 boundary-precision NN (Noire's spec, plans/c3_ncrnaboundaryscorer_spec.md; G11,
		//2026-08-28): a further small adjustment on top of the alignment-extent snap above,
		//mirroring TrnaCaller.refineBoundaryNN's placement immediately after its own trim
		//(TrnaCaller.java:644-658). nnInvoked is EXACTLY the original condition under
		//which refineBoundaryNN used to run (unreachable at all when either early-return guard
		//above had fired, AND boundary5Net/boundary3Net non-null) -- boundary5Net==null (the
		//structural off-by-default) still means this never executes and orf.start/orf.stop are
		//exactly what the snap above left them, byte-identical to pre-C3 behavior. nnInvoked is
		//a PRODUCTION-RUN FACT (did this exact call actually run the net), NOT a bootstrap-
		//suitability signal -- bootstrap capture runs are intentionally NN-off (no net staged
		//yet), so trimSucceeded=true with nnInvoked=false is the NORMAL, expected bootstrap
		//state, not a degraded one. Citan, 2026-08-28: renamed from refinementEligible, which
		//conflated "was the net literally invoked this run" with "is this locus structurally
		//usable for training" -- they are different questions with different answers during
		//bootstrap capture.
		final boolean nnInvoked=trimSucceeded && boundary5Net!=null && boundary3Net!=null;
		//Boundary-NN instrumentation (Citan/Brian, 2026-08-28): fires for EVERY accepted locus
		//this method is called on -- not gated on trimSucceeded -- so the accepted denominator
		//is never undercounted. Single null-check, off by default -- see instrumentSink's
		//javadoc.
		if(instrumentSink!=null){captureInstrumentation(orf, bases, model, wStart, wStop, trimSucceeded, nnInvoked);}
		if(nnInvoked){refineBoundaryNN(orf, bases, model);}
	}

	/** Builds the private window copy (never a live bases[] reference -- see
	 * NcrnaBoundaryInstrumentSink's thread-safety javadoc) and invokes the sink. Split out of
	 * trimToAlignmentExtent so the hot (instrumentation-off) path is a single null-check with
	 * no other cost. */
	private void captureInstrumentation(Orf orf, byte[] bases, int model, int wStart, int wStop,
			boolean trimSucceeded, boolean nnInvoked){
		final int copyFrom=Tools.max(0, orf.start-INSTRUMENT_CAPTURE_PAD);
		final int copyTo=Tools.min(bases.length-1, orf.stop+INSTRUMENT_CAPTURE_PAD);
		final byte[] windowCopy=Arrays.copyOfRange(bases, copyFrom, copyTo+1);
		instrumentSink.capture(orf.scafName, orf.strand, model, wStart, wStop,
			orf.start, orf.stop, windowCopy, copyFrom, trimSucceeded, nnInvoked);
	}

	/** Applies the boundary-precision NN's refinement to an already-trimmed, already-verified
	 * Orf. Builds a padded window around the current [orf.start,orf.stop] (PAD=10, matching
	 * TrnaCaller.refineBoundaryNN's convention -- comfortably covers the current per-family
	 * candidate ranges), then defers to NcrnaBoundaryScorer.refineBoundaries.
	 * No-op (leaves orf untouched) if the padded window can't hold a valid base candidate. */
	private void refineBoundaryNN(Orf orf, byte[] bases, int model){
		final int PAD=10;
		final int winStart=Tools.max(0, orf.start-PAD);
		final int winStop=Tools.min(bases.length-1, orf.stop+PAD);
		final byte[] window=Arrays.copyOfRange(bases, winStart, winStop+1);
		final int s=orf.start-winStart, e=orf.stop-winStart;
		if(s<0 || e>=window.length || e-s<15){return;}
		final float contigGC=contigGC(bases);
		final BaseGraph modelGraph=(models!=null && model<models.length ? models[model] : null);
		final int[] offsets=NcrnaBoundaryScorer.refineBoundaries(boundary5Net, boundary3Net, window, s, e,
			library[model], modelGraph, boundaryStartTable, boundaryStopTable,
			boundaryStartInside, boundaryStartOutside, boundaryStopInside, boundaryStopOutside,
			contigGC, boundaryMeanLen, boundaryStartOffsets, boundaryStopOffsets);
		orf.start+=offsets[0];
		orf.stop+=offsets[1];
	}

	/** Per-contig GC cache (identity-keyed on the bases[] reference), mirrors TrnaCaller's own
	 * contigGC cache exactly -- a NcrnaScavenger instance processes one contig/strand's bases[]
	 * across many calls within one scavenge() invocation, so recomputing GC from scratch per
	 * call would rescan the whole contig once per locus. */
	private float contigGC(byte[] bases){
		if(bases!=gcCacheBases){gcCacheValue=shared.Tools.calcGC(bases); gcCacheBases=bases;}
		return gcCacheValue;
	}
	private byte[] gcCacheBases=null;
	private float gcCacheValue=0;

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
	//Forward-ported from Noire's tree (2026-08-28, C3 merge): idPass/idBorderline are now
	//FINAL, set unconditionally by every constructor (0.75f/0.65f remain the fallback values
	//for the 8-arg and 15-arg delegating overloads) -- previously mutable with the same
	//defaults but never actually reassigned outside a constructor, so this is a safety
	//tightening, not a behavior change for any existing caller.
	final float idPass;
	final float idBorderline;
	//Forward-ported from Noire's tree (2026-08-28, C3 merge): per-family orfScore formula
	//constants -- see the 15-arg constructor's javadoc for the full rationale.
	final float scoreA;
	final float scoreB;
	float hbmPass=0.75f;
	int quantumThresh=120;
	int nearbyPad=200;
	float collapseFrac=0.9f;
	int trimExt=10;
	boolean scavengePass2=true;

	private long alignmentCount=0;

	//C3 boundary-precision-NN resources (G11, 2026-08-28) -- boundary5Net==null (this
	//instance's default unless the full constructor is used with real templates) means OFF,
	//structurally: refineBoundaryNN is never called (see trimToAlignmentExtent's pre-call
	//guard) and these fields are never read. Per-instance CLONES of the family's shared
	//read-only templates (thread safety -- see the full constructor's javadoc).
	private final CellNet boundary5Net, boundary3Net;
	private final TrnaBoundaryFeatures.NinemerTable boundaryStartTable, boundaryStopTable;
	private final int boundaryStartInside, boundaryStartOutside, boundaryStopInside, boundaryStopOutside;
	private final float boundaryMeanLen;
	private final int[] boundaryStartOffsets, boundaryStopOffsets;

	//Boundary-NN instrumentation (Citan/Brian, 2026-08-28) -- opt-in, off by default. A setter
	//rather than another constructor param: NcrnaScavenger already has 3 overloaded
	//constructors with 8-16 params each; threading one more opt-in field through all of them
	//would touch every call site (CallGenes/NcrnaFamily/GeneCaller) for a feature that's off in
	//every production run. Null-checked once per ACCEPTED locus only (see
	//trimToAlignmentExtent) -- not the per-candidate-window hot path, so this costs nothing on
	//the scan loop either way.
	private NcrnaBoundaryInstrumentSink instrumentSink=null;

	//B4 seed-trigger/workload instrumentation (Citan/G11, 2026-08-29) -- opt-in, off by default,
	//same cost model as instrumentSink above (single null-check per call site, zero cost when
	//off). Separate field from instrumentSink: different question (workload/gate-metric capture
	//on EVERY scanned strand and EVERY scheduled window, vs accepted-locus-only boundary state),
	//different consumer (a B4/B5 family evaluator driver, not boundary-NN training).
	private NcrnaWorkloadInstrumentSink workloadSink=null;

	/** Arms (or disarms, via null) B4 workload instrumentation capture. Package-visible: only an
	 * evaluation driver in this package should call this, mirroring setInstrumentSink's scoping
	 * rationale (never wired to a production CallGenes flag without an explicit opt-in gate). No
	 * arm-time validation needed here (unlike setInstrumentSink's modelNames-length assert) --
	 * seedHits/scheduledWindow never index into modelNames, so there is no equivalent silent-skip
	 * hazard to guard against at arm time. */
	void setWorkloadSink(NcrnaWorkloadInstrumentSink sink){workloadSink=sink;}

	/** Arms (or disarms, via null) boundary-NN instrumentation capture. Package-visible: only
	 * an instrumentation driver in this package should call this, never production CallGenes
	 * flag plumbing without an explicit opt-in flag gating it. */
	/** Fail-loud arming, per Citan (2026-08-28): all 3 accepted branches in alignWindow call
	 * trimToAlignmentExtent only inside `if(annotate && modelNames!=null && bestModel&lt;
	 * modelNames.length)` -- if modelNames were null or shorter than library, an accepted locus
	 * would silently skip trim AND capture entirely, undercounting the accepted denominator
	 * again, exactly the class of bug just fixed for the trim-guard case. Validating at ARM time
	 * (not silently at run time) guarantees that once instrumentation is successfully armed,
	 * every accepted bestModel is GUARANTEED to be a valid modelNames index, so this specific
	 * gate can never again be the reason a capture is skipped. A plain assert (not
	 * KillSwitch.assertDie): this runs on the single calling thread that builds and arms an
	 * instrumentation driver, never a producer/consumer worker thread, so an AssertionError here
	 * cannot leave anything else silently hung -- the textbook case for a plain assert per the
	 * assertions skill. Sink-off (sink==null) never runs this check -- zero change to production
	 * arming behavior (there is none) or cost. */
	void setInstrumentSink(NcrnaBoundaryInstrumentSink sink){
		if(sink!=null){
			assert(modelNames!=null && modelNames.length==library.length) : "Cannot arm boundary-NN "
				+"instrumentation: modelNames is "+(modelNames==null ? "null" : "length "+modelNames.length)
				+" but library has "+(library==null ? "null" : ""+library.length)+" models -- they must be "
				+"non-null and equal length, or alignWindow's annotate/modelNames guard would silently skip "
				+"trim+capture for some accepted loci (NcrnaScavenger.java, the 3 trimToAlignmentExtent call "
				+"sites), undercounting the accepted denominator just like the trim-guard bug this replaces.";
		}
		instrumentSink=sink;
	}

	/** Padding beyond the post-trim [start,stop] captured into the instrumentation window copy
	 * -- must cover both the family-configured boundary candidate arrays and the enrichment profile's own local
	 * +-2 radius plus the widest currently-staged k-mer window (k=11, srp_small) -- 4(sweep)+
	 * 2(local radius)+11(k)=17 is the true minimum reach past the boundary; 30 leaves real
	 * margin without meaningfully growing the copy. */
	static final int INSTRUMENT_CAPTURE_PAD=30;
}
