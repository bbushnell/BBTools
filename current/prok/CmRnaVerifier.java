package prok;

import java.io.File;
import java.util.ArrayList;
import aligner.CovarianceModel;
import aligner.CovarianceModelGate;
import aligner.CovarianceModelParser;
import aligner.CovarianceModelWindowGate;
import aligner.CovarianceModelConsensusMap;
import aligner.CovarianceModelConsensusMaps;
import aligner.CovarianceModelPlacementGate;
import idaligner.AlignmentStats;
import map.ObjectMap;
import dna.Data;
import fileIO.ByteStreamWriter;
import shared.KillSwitch;
import structures.ByteBuilder;

/** Per-worker CM verification; shared immutable model/scorer topology.
 * Caller endpoints/path scores remain unchanged. Each candidate gets its own
 * best half-shorter-overlapping interval; any model passing its provisional GA
 * admits it. The highest-scoring evaluated passing model is selected for diagnostics.
 * 18S uses the refined extent, not window search. Its scoring mode is independent
 * of small-family defaults and consistent across exact, placement and fallback paths.
 * @author Brian, Raiden
 */
final class CmRnaVerifier {
	CmRnaVerifier(Library library_, ByteStreamWriter out_){
		this(library_, out_, null);
	}
	CmRnaVerifier(Library library_, ByteStreamWriter out_, ByteStreamWriter placementOut_){
		if(library_==null){throw new IllegalArgumentException("CM worker requires a bound model library");}
		library=library_;out=out_;placementOut=placementOut_;
		greedyBuffer=library.trnaExtendedGreedy || library.smallGreedy ? new ByteBuilder() : null;
		int max=0;for(CovarianceModelWindowGate[] models : library.gates){if(models!=null){max=Math.max(max, models.length);}}
		decisions=new CovarianceModelWindowGate.Decision[max];
		wideDecisions=new CovarianceModelWindowGate.Decision[max];
		extendedDecisions=new CovarianceModelWindowGate.Decision[max];
	}
	boolean supports(int type){return type>=0 && type<library.gates.length
		&& (library.gates[type]!=null || type==ProkObject.r18S && library.exact18S!=null);}
	/** Explicit window on the already oriented chromosome; candidate may extend
	 * beyond its window after boundary refinement, without inventing new padding. */
	boolean verify(Orf orf, byte[] bases, int windowStart, int windowStop, String route){
		return verify(orf, bases, windowStart, windowStop, route, orf.trnaModel, null);
	}
	boolean verify(Orf orf, byte[] bases, int windowStart, int windowStop, String route, String modelName, AlignmentStats stats){
		if(orf.type==ProkObject.r18S && library.placement18S!=null){return verifyPlacement(orf, bases, windowStart, route, modelName, stats);}
		if(orf.type==ProkObject.r18S && library.exact18S!=null){return verify18S(orf, bases, route);}
		return verify(orf, bases, windowStart, windowStop, orf.start, orf.stop, 0, Integer.MAX_VALUE, 0, route, null, 0, 0);
	}
	/** Exact scoring consumes the full refined extent. No searched optimum
	 * or traceback is available; record this explicitly, leaving caller ends intact. */
	private boolean verify18S(Orf orf, byte[] bases, String route){
		assert(orf.flipped()==0 && orf.type==ProkObject.r18S && library.exact18S!=null):
			KillSwitch.assertDie("18S exact scoring requires an oriented refined candidate and its bound RF01960 gate");
		final CovarianceModelGate gate=library.exact18S;final CovarianceModelGate.Decision r;
		try{r=gate.evaluate(bases, orf.start, orf.stop);}
		catch(RuntimeException | AssertionError e){
			KillSwitch.assertDie("18S exact verification failed on "+orf.scafName+":"+orf.start+"-"+orf.stop+": "+e);throw e;
		}
		candidates[orf.type]++;if(r.accepted){accepted[orf.type]++;}
		if(out!=null){
			buffer.clear();safe(orf.scafName);buffer.tab().append(orf.strand==0 ? '+' : '-');coordinates(orf, orf.start, orf.stop);
			buffer.tab().append(route).append(gate.local ? "_local_exact_extent" : "_exact_extent").tab().append(ProkObject.typeStrings[orf.type]).tab().append(gate.model.accession)
				.tab().append(r.length).tab().appendSlow(r.bits).tab().appendSlow(gate.threshold).tab().append(r.accepted ? "PASS" : "REJECT")
				.tab().append(r.seconds, 9).tab().append(r.peakScoreCells).tab();safe(orf.ncrnaFamily);
			coordinates(orf, orf.start, orf.stop);
			// Both modes consume the supplied extent; no claim that EL is homologous.
			coordinates(orf, Float.isFinite(r.bits) ? orf.start : -1, Float.isFinite(r.bits) ? orf.stop : -1);
			buffer.append("\tNA\tNA\tNA\t0\ttrue\t").append(Math.multiplyExact(r.peakScoreCells, 4L)).append("\tfalse\n");
			synchronized(out){out.print(buffer);}
		}
		return r.accepted;
	}
	private boolean verifyPlacement(Orf orf, byte[] bases, int windowStart, String route, String modelName, AlignmentStats stats){
		assert(orf.flipped()==0 && orf.type==ProkObject.r18S):
			KillSwitch.assertDie("Placement must use oriented candidate coordinates before caller strand conversion");
		final CovarianceModelPlacementGate gate=library.placement18S;
		final CovarianceModelConsensusMap map=modelName==null ? null : library.placementMaps.get(modelName);
		final CovarianceModelPlacementGate.Decision r;
		try{r=gate.evaluate(bases, orf.start, orf.stop, map, stats, windowStart);}
		catch(RuntimeException | AssertionError e){KillSwitch.assertDie("18S placement failed on "+orf.scafName+":"+orf.start+"-"+orf.stop+": "+e);throw e;}
		candidates[orf.type]++;if(r.accepted){accepted[orf.type]++;}
		if(out!=null){
			final long cells=Math.max(r.bandScoreCells, r.exactScoreCells);
			final long bytes=Math.max(Math.addExact(Math.multiplyExact(r.bandScoreCells, 4L), r.traceBytes), Math.multiplyExact(r.exactScoreCells, 4L));
			buffer.clear();safe(orf.scafName);buffer.tab().append(orf.strand==0 ? '+' : '-');coordinates(orf, orf.start, orf.stop);
			buffer.tab().append(route).append(gate.local ? "_local_placement_extent" : "_placement_extent").tab().append(ProkObject.typeStrings[orf.type]).tab().append(gate.model.accession)
				.tab().append(r.length).tab().appendSlow(r.bits).tab().appendSlow(gate.threshold).tab().append(r.accepted ? "PASS" : "REJECT")
				.tab().append(r.seconds, 9).tab().append(cells).tab();safe(orf.ncrnaFamily);
			coordinates(orf, orf.start, orf.stop);coordinates(orf, Float.isFinite(r.bits) ? orf.start : -1, Float.isFinite(r.bits) ? orf.stop : -1);
			buffer.append("\tNA\tNA\tNA\t0\ttrue\t").append(bytes).append("\tfalse\n");
			synchronized(out){out.print(buffer);}
		}
		if(placementOut!=null){
			buffer.clear();safe(orf.scafName);buffer.tab().append(orf.strand==0 ? '+' : '-');coordinates(orf, orf.start, orf.stop);
			buffer.tab().append(route).append(gate.local ? "_local" : "").tab();safe(modelName);buffer.tab().append(gate.model.accession).tab().append(r.length).tab().append(r.finalRadius).tab().append(r.fallbackReason).tab();
			if(Float.isNaN(r.rawBits)){buffer.append("NA");}else{buffer.appendSlow(r.rawBits);}
			buffer.tab().appendSlow(r.bits).tab().append(r.accepted ? "PASS" : "REJECT")
				.tab().append(r.anchorSeconds, 9).tab().append(r.bandSeconds, 9).tab().append(r.exactSeconds, 9).tab().append(r.seconds, 9)
				.tab().append(r.bandScoreCells).tab().append(r.traceBytes).tab().append(r.exactScoreCells)
				.tab().append(r.contacts.steps).tab().append(r.contacts.startLow).tab().append(r.contacts.startHigh)
				.tab().append(r.contacts.endLow).tab().append(r.contacts.endHigh).tab().append(r.bandAttempts).nl();
			synchronized(placementOut){placementOut.print(buffer);}
		}
		return r.accepted;
	}
	/** Preserve the caller's splice and pad its outer genomic bounds. Association
	 * is measured in the mature-product coordinate system used by the CM. */
	boolean verifySpliced(Orf orf, byte[] genome, int intronStart, int intronStop, int pad, String route){
		assert(intronStart>orf.start && intronStart<=intronStop && intronStop<orf.stop && pad>=0):
			KillSwitch.assertDie("A spliced CM product requires the caller's internal genomic intron and nonnegative flank");
		final int start=Math.max(0, orf.start-pad), stop=Math.min(genome.length-1, orf.stop+pad);
		final int gap=intronStop-intronStart+1, left=intronStart-start;
		final byte[] product=new byte[stop-start+1-gap];
		System.arraycopy(genome, start, product, 0, left);
		System.arraycopy(genome, intronStop+1, product, left, stop-intronStop);
		return verify(orf, product, 0, product.length-1, orf.start-start, orf.stop-start-gap, start, left, gap, route, genome, start, stop);
	}
	private boolean verify(Orf orf, byte[] sequence, int from, int to, int candidateStart, int candidateStop,
			int offset, int gapAt, int gapLength, String route, byte[] unspliced, int genomicStart, int genomicStop){
		assert(orf!=null && orf.flipped()==0 && supports(orf.type)):
			KillSwitch.assertDie("CM must score an oriented candidate before strand-coordinate flipping");
		final CovarianceModelWindowGate[] models=library.gates[orf.type];
		final int fullFrom=from, fullTo=to;
		if(library.windowSlack>=0){
			from=(int)Math.max(from, (long)candidateStart-library.windowSlack);
			to=(int)Math.min(to, (long)candidateStop+library.windowSlack);
		}
		assert(from<=to):KillSwitch.assertDie("Candidate slack must intersect its original CM search window");
		int selected=-1, evaluated=0, wideSelected=-1, wideEvaluated=0, extendedSelected=-1, extendedEvaluated=0;
		final boolean spliced=gapLength>0;
		final boolean smallGreedy=library.smallGreedy && (orf.type==ProkObject.tRNA || orf.type==ProkObject.r5S);
		final int wideFrom=spliced ? genomicStart : fullFrom, wideTo=spliced ? genomicStop : fullTo;
		final byte[] genomic=spliced ? unspliced : sequence;
		final int extendedFrom=(int)Math.max(0L, (long)wideFrom-library.localExtend);
		final int extendedTo=(int)Math.min(genomic.length-1L, (long)wideTo+library.localExtend);
		try{
			evaluated=evaluateModels(orf.type, models, decisions, sequence, from, to, candidateStart, candidateStop,smallGreedy);
			selected=select(decisions, evaluated);
			if(library.localWide && (spliced || !decisions[selected].accepted && (from!=fullFrom || to!=fullTo))){
				// D35: retain original padding; spliced candidates additionally search
				// the unexcised genomic window even when their mature product passed.
				wideEvaluated=evaluateModels(orf.type, models, wideDecisions, spliced ? unspliced : sequence,
					wideFrom, wideTo, spliced ? orf.start : candidateStart, spliced ? orf.stop : candidateStop,smallGreedy);
				wideSelected=select(wideDecisions, wideEvaluated);
			}
			if(library.localExtend>0 && !decisions[selected].accepted && (wideSelected<0 || !wideDecisions[wideSelected].accepted)
					&& (extendedFrom!=wideFrom || extendedTo!=wideTo)){
				// Retry only rejected candidates, once, against the unexcised genomic sequence.
				extendedEvaluated=evaluateModels(orf.type, models, extendedDecisions, genomic,
					extendedFrom, extendedTo, spliced ? orf.start : candidateStart, spliced ? orf.stop : candidateStop,
					smallGreedy || library.trnaExtendedGreedy && orf.type==ProkObject.tRNA);
				extendedSelected=select(extendedDecisions, extendedEvaluated);
			}
		}catch(RuntimeException | AssertionError e){
			KillSwitch.assertDie("CM window verification failed on "+orf.scafName+":"+orf.start+"-"+orf.stop+": "+e);throw e;
		}
		assert(selected>=0):KillSwitch.assertDie("A supported RNA family must evaluate at least one bound CM");
		final boolean pass=decisions[selected].accepted || wideSelected>=0 && wideDecisions[wideSelected].accepted
			|| extendedSelected>=0 && extendedDecisions[extendedSelected].accepted;
		final boolean chooseWide=wideSelected>=0 && better(wideDecisions[wideSelected], decisions[selected]);
		final boolean chooseExtended=extendedSelected>=0 && better(extendedDecisions[extendedSelected], chooseWide ? wideDecisions[wideSelected] : decisions[selected]);
		candidates[orf.type]++;if(pass){accepted[orf.type]++;}
		if(out!=null){
			buffer.clear();
			appendDiagnostics(orf, models, decisions, evaluated, chooseWide || chooseExtended ? -1 : selected,
				from, to, offset, gapAt, gapLength, route);
			if(wideEvaluated>0){appendDiagnostics(orf, models, wideDecisions, wideEvaluated, chooseWide && !chooseExtended ? wideSelected : -1,
				wideFrom, wideTo, 0, Integer.MAX_VALUE, 0, route+"_wide");}
			if(extendedEvaluated>0){appendDiagnostics(orf, models, extendedDecisions, extendedEvaluated, chooseExtended ? extendedSelected : -1,
				extendedFrom, extendedTo, 0, Integer.MAX_VALUE, 0, route+"_extended");}
			// One optional write per candidate; all attempted model/window costs remain visible.
			synchronized(out){out.print(buffer);}
		}
		return pass;
	}
	private int evaluateModels(int type, CovarianceModelWindowGate[] models, CovarianceModelWindowGate.Decision[] results,
			byte[] sequence, int from, int to, int candidateStart, int candidateStop,boolean greedy){
		assert(models.length>0 && results.length>=models.length):"Each supported family needs space for every bound CM decision";
		int count=0;
		for(int i=0; i<models.length; i++){
			if(i>0 && type==ProkObject.tRNA && library.secFallback && results[0].accepted){break;}
			results[i]=greedy ? models[i].evaluateGreedy(sequence,from,to,candidateStart,candidateStop)
				: models[i].evaluate(sequence, from, to, candidateStart, candidateStop);count++;
		}
		return count;
	}
	private static int select(CovarianceModelWindowGate.Decision[] results, int count){
		assert(count>0):"A supported family must complete at least one model evaluation before selection";
		int selected=0;for(int i=1; i<count; i++){if(better(results[i], results[selected])){selected=i;}}return selected;
	}
	private static boolean better(CovarianceModelWindowGate.Decision a, CovarianceModelWindowGate.Decision b){
		assert(a!=null && b!=null):"Selection compares completed model evaluations, never stale array slots";
		return a.accepted && !b.accepted || a.accepted==b.accepted && a.bits>b.bits;
	}
	private void appendDiagnostics(Orf orf, CovarianceModelWindowGate[] models, CovarianceModelWindowGate.Decision[] results,
			int evaluated, int selected, int from, int to, int offset, int gapAt, int gapLength, String route){
		assert(evaluated>0 && selected<evaluated):"Only evaluated rows may be selected; -1 denotes the other window's winner";
		for(int i=0; i<evaluated; i++){
				final CovarianceModelWindowGate.Decision r=results[i];
				safe(orf.scafName);buffer.tab().append(orf.strand==0 ? '+' : '-');
				coordinates(orf, orf.start, orf.stop);
				buffer.tab().append(route).append(library.local ? "_local" : "").tab().append(ProkObject.typeStrings[orf.type]).tab().append(models[i].model.accession)
					.tab().append(r.windowLength).tab().appendSlow(r.bits).tab().appendSlow(models[i].threshold)
					.tab().append(r.accepted ? "PASS" : r.unrestrictedBits>=models[i].threshold ? "REJECT_OFF_CANDIDATE" : "REJECT")
					.tab().append(r.seconds, 9).tab().append(r.cells).tab();safe(orf.ncrnaFamily);
				coordinates(orf, genomic(from, offset, gapAt, gapLength), genomic(to, offset, gapAt, gapLength));
				coordinates(orf, genomic(r.from, offset, gapAt, gapLength), genomic(r.to, offset, gapAt, gapLength));
				buffer.tab().appendSlow(r.unrestrictedBits);
				coordinates(orf, genomic(r.unrestrictedFrom, offset, gapAt, gapLength), genomic(r.unrestrictedTo, offset, gapAt, gapLength));
				buffer.tab().append(r.ties).tab().append(i==selected).tab().append(r.matrixBytes()).tab().append(gapLength>0).nl();
				if(r.greedyHits!=null){appendGreedyDiagnostics(orf,models[i],r,from,to,offset,gapAt,gapLength,route);}
			}
	}
	/** Supplemental evidence preserves the stable25-column CM stream. One
	 * synchronized stderr write per model contains all extracted nonoverlapping hits. */
	private void appendGreedyDiagnostics(Orf orf,CovarianceModelWindowGate model,CovarianceModelWindowGate.Decision r,
			int from,int to,int offset,int gapAt,int gapLength,String route){
		assert(r.greedyHits!=null && greedyBuffer!=null):"Greedy diagnostic hits require the enabled per-worker buffer";
		// Association is measured where the CM searched. A mature product excludes
		// the caller's intron; expanding first would change the half-shorter test.
		final int candidateStart=orf.start-offset, candidateStop=orf.stop-offset-gapLength;
		final int windowStart=genomic(from,offset,gapAt,gapLength), windowStop=genomic(to,offset,gapAt,gapLength);
		greedyBuffer.clear();
		for(int h=0;h<r.greedyHits.starts.size;h++){
			final int hitStart=r.greedyOffset+r.greedyHits.starts.get(h), hitStop=r.greedyOffset+r.greedyHits.stops.get(h);
			final int a=genomic(hitStart,offset,gapAt,gapLength), z=genomic(hitStop,offset,gapAt,gapLength);
			final int overlap=Math.max(0,Math.min(hitStop,candidateStop)-Math.max(hitStart,candidateStart)+1);
			greedyBuffer.append("CM_GREEDY_HIT\t").append(RrnaResourceIO.first(orf.scafName)).tab().append(orf.strand==0?'+':'-')
				.tab().append(orf.strand==0?orf.start+1:orf.scaflen-orf.stop).tab().append(orf.strand==0?orf.stop+1:orf.scaflen-orf.start)
				.tab().append(route).tab().append(model.model.accession).tab().append(h+1)
				.tab().append(orf.strand==0?a+1:orf.scaflen-z).tab().append(orf.strand==0?z+1:orf.scaflen-a)
				.tab().appendSlow(r.greedyHits.bits.get(h)).tab().appendSlow(model.threshold)
				.tab().append(2L*overlap>=Math.min((long)hitStop-hitStart+1,candidateStop-candidateStart+1L)).tab().append(h==r.greedyHits.selected)
				.tab().append(orf.strand==0?windowStart+1:orf.scaflen-windowStop).tab().append(orf.strand==0?windowStop+1:orf.scaflen-windowStart).nl();
		}
		if(greedyBuffer.length>0){synchronized(System.err){System.err.print(greedyBuffer.toString());}}
	}
	private static int genomic(int position, int offset, int gapAt, int gapLength){
		return position<0 ? -1 : offset+position+(position>=gapAt ? gapLength : 0);
	}
	private void coordinates(Orf orf, int from, int to){
		if(from<0 || to<0){buffer.append("\tNA\tNA");return;}
		assert(from<=to && to<orf.scaflen):KillSwitch.assertDie("CM coordinate conversion must fit the original genomic scaffold");
		final int a=orf.strand==0 ? from : orf.scaflen-to-1, b=orf.strand==0 ? to : orf.scaflen-from-1;
		buffer.tab().append(a+1).tab().append(b+1);
	}
	private void safe(String text){
		if(text==null){buffer.append('.');return;}
		for(int i=0; i<text.length(); i++){
			final char c=text.charAt(i);buffer.append(c=='\t' || c=='\n' || c=='\r' ? ' ' : c);
		}
	}
	static final class Library {
		/** Called once during resource loading, before workers share this library. */
		void loadPlacement(String path, int radius, long maxCells, ArrayList<NcrnaFamily> families){
			loadPlacement(path, radius, maxCells, families, false);
		}
		void loadPlacement(String path, int radius, long maxCells, ArrayList<NcrnaFamily> families, boolean adaptive){
			if(exact18S==null || placement18S!=null){throw new IllegalArgumentException("Placement resources require enabled18S and a single initialization");}
			final ArrayList<String> names=new ArrayList<String>();final ArrayList<byte[]> sequences=new ArrayList<byte[]>();
			for(NcrnaFamily family:families){if(family.outputType==ProkObject.r18S){
				if(family.modelNames==null || family.modelNames.length!=family.library.length){throw new IllegalArgumentException("18S anchors require stable consensus names");}
				for(int i=0; i<family.library.length; i++){names.add(family.modelNames[i]);sequences.add(family.library[i]);}
			}}
			final CovarianceModelConsensusMap[] maps=CovarianceModelConsensusMaps.read(exact18S.model, path, names.toArray(new String[0]), sequences.toArray(new byte[0][]));
			for(int i=0; i<maps.length; i++){if(maps[i]!=null){placementMaps.put(names.get(i), maps[i]);}}
			placement18S=new CovarianceModelPlacementGate(exact18S.model, exact18S.threshold, radius, maxCells, maxCells, adaptive, exact18S.local);paths.add(path);
		}
		Library(String directory, boolean trna, boolean fiveS, long maxCells){
			this(directory, trna, fiveS, maxCells, -1, false);
		}
		Library(String directory, boolean trna, boolean fiveS, long maxCells, int slack, boolean fallback){
			this(directory, trna, fiveS, false, maxCells, slack, fallback);
		}
		Library(String directory, boolean trna, boolean fiveS, boolean eighteenS, long maxCells, int slack, boolean fallback){
			this(directory, trna, fiveS, eighteenS, maxCells, slack, fallback, false);
		}
		Library(String directory, boolean trna, boolean fiveS, boolean eighteenS, long maxCells, int slack, boolean fallback, boolean local_){
			this(directory, trna, fiveS, eighteenS, maxCells, slack, fallback, local_, false);
		}
		Library(String directory, boolean trna, boolean fiveS, boolean eighteenS, long maxCells, int slack, boolean fallback, boolean local_, boolean localWide_){
			this(directory, trna, fiveS, eighteenS, maxCells, slack, fallback, local_, localWide_, 0);
		}
		Library(String directory, boolean trna, boolean fiveS, boolean eighteenS, long maxCells, int slack, boolean fallback, boolean local_, boolean localWide_, int localExtend_){
			this(directory, trna, fiveS, eighteenS, maxCells, slack, fallback, local_, localWide_, localExtend_, local_);
		}
		/** CLI defaults may differ by family; legacy constructors retain one explicit mode. */
		Library(String directory, boolean trna, boolean fiveS, boolean eighteenS, long maxCells, int slack, boolean fallback, boolean local_, boolean localWide_, int localExtend_, boolean local18S){
			this(directory,trna,fiveS,eighteenS,maxCells,slack,fallback,local_,localWide_,localExtend_,local18S,false);
		}
		Library(String directory, boolean trna, boolean fiveS, boolean eighteenS, long maxCells, int slack, boolean fallback, boolean local_, boolean localWide_, int localExtend_, boolean local18S,boolean trnaExtendedGreedy_){
			this(directory,trna,fiveS,eighteenS,maxCells,slack,fallback,local_,localWide_,localExtend_,local18S,trnaExtendedGreedy_,false);
		}
		Library(String directory, boolean trna, boolean fiveS, boolean eighteenS, long maxCells, int slack, boolean fallback, boolean local_, boolean localWide_, int localExtend_, boolean local18S,boolean trnaExtendedGreedy_,boolean smallGreedy_){
			if(maxCells<=0){throw new IllegalArgumentException("cmmaxcells must be positive");}
			if(slack< -1){throw new IllegalArgumentException("cmwindowslack must be -1 (full window) or nonnegative");}
			if(localWide_ && !local_){throw new IllegalArgumentException("cmlocalwide requires cmlocal=t");}
			if(localExtend_<0 || localExtend_>0 && (!local_ || !localWide_)){throw new IllegalArgumentException("cmlocalextend requires nonnegative padding and cmlocal=t cmlocalwide=t");}
			windowSlack=slack;secFallback=fallback;local=local_;localWide=localWide_;
			localExtend=localExtend_;
			if(trnaExtendedGreedy_ && (!local_ || !localWide_ || localExtend_<=0)){throw new IllegalArgumentException("cmtrnaextendedgreedy requires local/wide scoring and a positive extended retry");}
			trnaExtendedGreedy=trnaExtendedGreedy_;
			if(smallGreedy_ && !local_){throw new IllegalArgumentException("cmsmallgreedy requires local small-family scoring");}
			smallGreedy=smallGreedy_;
			if(trna){gates[ProkObject.tRNA]=new CovarianceModelWindowGate[]{load(directory, "RF00005", maxCells), load(directory, "RF01852", maxCells)};}
			if(fiveS){gates[ProkObject.r5S]=new CovarianceModelWindowGate[]{load(directory, "RF00001", maxCells)};}
			final CovarianceModel model=eighteenS ? readModel(directory, "RF01960") : null;
			exact18S=model==null ? null : new CovarianceModelGate(model, CovarianceModelGate.gatheringThreshold(model), maxCells, local18S);
		}
		private CovarianceModelWindowGate load(String directory, String accession, long maxCells){
			final CovarianceModel model=readModel(directory, accession);
			return new CovarianceModelWindowGate(model, CovarianceModelGate.gatheringThreshold(model), maxCells, local);
		}
		private CovarianceModel readModel(String directory, String accession){
			final String path=directory==null ? Data.findPath("?"+accession+".cm") : new File(directory, accession+".cm").getPath();
			final CovarianceModel model=CovarianceModelParser.read(path);
			if(!accession.equals(model.accession)){throw new IllegalArgumentException("Wrong CM family: expected "+accession+", found "+model.accession);}
			paths.add(path);return model;
		}
		final CovarianceModelWindowGate[][] gates=new CovarianceModelWindowGate[ProkObject.typeStrings.length][];
		final ArrayList<String> paths=new ArrayList<String>();
		final int windowSlack;
		final boolean secFallback;
		final boolean local;
		final boolean localWide;
		final int localExtend;
		final boolean trnaExtendedGreedy;
		final boolean smallGreedy;
		final CovarianceModelGate exact18S;
		CovarianceModelPlacementGate placement18S;
		final ObjectMap<String,CovarianceModelConsensusMap> placementMaps=new ObjectMap<String,CovarianceModelConsensusMap>(String.class, CovarianceModelConsensusMap.class);
	}
	static final String PLACEMENT_HEADER="contig\tstrand\tstart1\tstop1\tcallerRoute\tconsensus\tcmAccession\tlength\tradius\tfallbackReason\tbandBits\teffectiveBits\tdecision\tanchorSeconds\tbandSeconds\texactSeconds\ttotalSeconds\tbandScoreCells\ttraceBytes\texactScoreCells"
		+"\tedgeContactSteps\tstartLowContacts\tstartHighContacts\tendLowContacts\tendHighContacts\tbandAttempts";
	static final String HEADER="contig\tstrand\tstart1\tstop1\tcallerRoute\tfamily\tcmAccession\tscoredLength\tbits\tprovisionalCutoff\tdecision\tseconds\tpeakScoreFloatCells\tncrnaFamily"
		+"\twindowStart1\twindowStop1\tcmStart1\tcmStop1\tunrestrictedBits\tunrestrictedStart1\tunrestrictedStop1\tbestIntervalTies\tselectedModel\tmatrixBytes\tspliced";
	final long[] candidates=new long[ProkObject.typeStrings.length], accepted=new long[ProkObject.typeStrings.length];
	private final Library library;
	private final ByteStreamWriter out;
	private final ByteStreamWriter placementOut;
	private final ByteBuilder buffer=new ByteBuilder();
	private final ByteBuilder greedyBuffer;
	private final CovarianceModelWindowGate.Decision[] decisions;
	private final CovarianceModelWindowGate.Decision[] wideDecisions;
	private final CovarianceModelWindowGate.Decision[] extendedDecisions;
}
