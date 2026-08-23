package prok;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import consensus.BaseGraph;
import dna.AminoAcid;
import fileIO.FileFormat;
import idaligner.AlignmentStats;
import idaligner.ScrabbleAligner;
import shared.Shared;
import stream.Read;
import stream.ReadInputStream;

/**
 * Training-vector generator for the tRNA boundary-precision NN (Brian's idea #3,
 * trna repo plans/boundary_precision_ideas.md). Reads a genome FASTA + a reference
 * GFF of true tRNA loci, and for each known start and stop boundary emits a labeled
 * feature vector for the true position (label=1) and for shifted
 * candidates (label=0), in train.sh's '#dims' TSV format (see nntrain skill).
 * Candidate ranges are asymmetric per boundary (START_OFFSETS/STOP_OFFSETS,
 * Aug 22 2026 -- see their javadoc), not a uniform +-2bp.
 *
 * This is offline machinery assembled while Dori is down (Aug 19, 2026): it is
 * Dori-BLIND by design (loads the shipped library/kmer-set from local BBTools
 * resources only), but 2 of the ~10 planned features need corpus-scale data that
 * doesn't exist yet -- see TrnaBoundaryFeatures for exactly which and why. Model
 * identification per locus reuses the real production path (TrnaKmerIndex
 * shortlist -> ScrabbleAligner identity, same as TrnaCaller), not a shortcut.
 *
 * Two input modes:
 * <p>1. fasta=<genome.fna[.gz]> gff=<ref.gff[.gz]> out=<vectors.tsv> [tablestart=<t.tsv> tablestop=<t.tsv>]
 *   Single-genome mode (original design): one genome + its reference GFF.
 * <p>2. fasta=<flanked.fa[.gz]> out=<vectors.tsv> tablestart=<t.tsv> tablestop=<t.tsv>
 *   Flanked-batch mode (Noire, Aug 21): reads a cutgff flank= corpus directly (the
 *   SAME artifact TrnaNinemerTableBuilder's fasta= table-build mode consumes) -- one
 *   pass, no genome+gff re-read, boundaries from lflank=/rflank= headers exactly like
 *   the table builder. Selected by OMITTING gff=. tablestart=/tablestop= are REQUIRED
 *   in this mode (no smoke fallback makes sense for a full-corpus training run).
 * <p>tablestart=/tablestop=: REAL corpus-derived tables (prok.TrnaNinemerTableBuilder's
 *   output). Both required together or omitted together (mode 1 only -- mode 2
 *   requires them). Each table is self-describing (boundary type + insideCount +
 *   outsideCount recorded in its own header, see TrnaNinemerTableBuilder.LoadedTable)
 *   -- no k= or offset= to pass or mismatch; the loaded table IS the window spec
 *   (Noire's hardening, Aug 21). Mode 1 without tablestart=/tablestop= falls back to
 *   the SMOKE-TEST table (self-derived from this one genome, see
 *   buildSmokeNinemerTable) -- never mistake this fallback for real training data.
 * @author Neptune
 */
public class TrnaBoundaryVectorGen {

	public static void main(String[] args) throws IOException {
		String fastaPath=null, gffPath=null, outPath=null, tableStartPath=null, tableStopPath=null;
		for(String a : args){
			String[] kv=a.split("=", 2);
			if(kv.length<2){continue;}
			if(kv[0].equalsIgnoreCase("fasta")){fastaPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("gff")){gffPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("out")){outPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("tablestart")){tableStartPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("tablestop")){tableStopPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("includecrossenrichment")){INCLUDE_CROSS_ENRICHMENT=parse.Parse.parseBoolean(kv[1]);}
			else if(kv[0].equalsIgnoreCase("noisyopposite")){NOISY_OPPOSITE_BOUNDARY=parse.Parse.parseBoolean(kv[1]);}
			else if(kv[0].equalsIgnoreCase("noisystartrate")){NOISY_START_RATE=Float.parseFloat(kv[1]);}
			else if(kv[0].equalsIgnoreCase("noisystoprate")){NOISY_STOP_RATE=Float.parseFloat(kv[1]);}
			else if(kv[0].equalsIgnoreCase("noisyseed")){NOISY_SEED=Long.parseLong(kv[1]);}
		}
		if(fastaPath==null || outPath==null){
			System.err.println("Usage (single-genome): fasta=<genome.fna[.gz]> gff=<ref.gff[.gz]> out=<vectors.tsv> [tablestart=<t.tsv> tablestop=<t.tsv>]");
			System.err.println("Usage (flanked-batch):  fasta=<flanked.fa[.gz]> out=<vectors.tsv> tablestart=<t.tsv> tablestop=<t.tsv>");
			System.exit(1);
		}
		if((tableStartPath==null)!=(tableStopPath==null)){
			System.err.println("tablestart= and tablestop= must be given together or not at all.");
			System.exit(1);
		}

		CallGenes.loadTrnaResources();
		CallGenes.getPhylumPGM(null);
		ProkObject.loadLongKmers();
		if(GeneCaller.trnaLibrary==null || GeneCaller.trnaLibrary.length==0){
			System.err.println("FATAL: trnaLibrary failed to load -- can't run this tool.");
			System.exit(2);
		}
		final byte[][] library=GeneCaller.trnaLibrary;
		final TrnaKmerIndex index=new TrnaKmerIndex(library, TrnaCaller.INDEX_K, TrnaCaller.ADAPTIVE_MINHITS,
			TrnaCaller.ADAPT_FLOOR, TrnaCaller.ADAPT_TOPFRAC, TrnaCaller.ADAPT_QFRAC,
			(TrnaCaller.INDEX_MINHITS_OVERRIDE>0 ? TrnaCaller.INDEX_MINHITS_OVERRIDE : 12));//12 mirrors
			//TrnaCaller's private INDEX_MINHITS_DEFAULT; kept in sync by eye, not by reference (private).

		if(gffPath==null){
			if(tableStartPath==null){
				System.err.println("Flanked-batch mode (no gff=) requires tablestart=/tablestop= -- no smoke fallback for a full-corpus run.");
				System.exit(1);
			}
			runFlankedBatchMode(fastaPath, outPath, tableStartPath, tableStopPath, library, index);
		}else{
			runSingleGenomeMode(fastaPath, gffPath, outPath, tableStartPath, tableStopPath, library, index);
		}
	}

	/** Flanked-batch mode (Noire, Aug 21): one pass over a cutgff flank= corpus, no
	 * genome+gff re-read. Mirrors TrnaNinemerTableBuilder's runFlankedFastaMode for
	 * header parsing/boundary extraction (same parseFlankValue, package-shared so the
	 * two readers can't drift on the header format). */
	private static void runFlankedBatchMode(String fastaPath, String outPath, String tableStartPath,
			String tableStopPath, byte[][] library, TrnaKmerIndex index) throws IOException {
		Shared.TRIM_READ_DESCRIPTION=false;//need the full header for lflank=/rflank=
		Shared.TRIM_RNAME=true;
		Read.TO_UPPER_CASE=true;

		final TrnaNinemerTableBuilder.LoadedTable lt1=TrnaNinemerTableBuilder.loadTable(tableStartPath);
		final TrnaNinemerTableBuilder.LoadedTable lt2=TrnaNinemerTableBuilder.loadTable(tableStopPath);
		assert(lt1.type==TrnaBoundaryFeatures.BoundaryType.START) : "tablestart= ("+tableStartPath
			+") is not labeled a start-boundary table -- wrong file, or tablestart=/tablestop= swapped.";
		assert(lt2.type==TrnaBoundaryFeatures.BoundaryType.STOP) : "tablestop= ("+tableStopPath
			+") is not labeled a stop-boundary table -- wrong file, or tablestart=/tablestop= swapped.";
		final TrnaBoundaryFeatures.NinemerTable startTable=lt1.table, stopTable=lt2.table;
		final int startInside=lt1.insideCount, startOutside=lt1.outsideCount;
		final int stopInside=lt2.insideCount, stopOutside=lt2.outsideCount;
		System.err.println("Using REAL corpus tables: "+tableStartPath+" (inside="+startInside+" outside="+startOutside
			+") / "+tableStopPath+" (inside="+stopInside+" outside="+stopOutside+")");

		final ArrayList<Read> reads=ReadInputStream.toReads(fastaPath, FileFormat.FA, -1);
		assert(!reads.isEmpty()) : "Empty or unreadable flanked fasta: "+fastaPath;
		System.err.println("Loaded "+reads.size()+" records from "+fastaPath);

		long noFlankHeader=0, malformedFlank=0, skippedNoModel=0, noContigGC=0;
		int written=0;
		try(PrintStream out=new PrintStream(new FileOutputStream(outPath))){
			//v3 feature set (Brian+Noire, Aug 21): stem(1) + ani(1) + enrichmentProfile(3) +
			//isStop(1) + tipFuzziness(3) + lengthRatio(1) + contigGC(1) = 11. Requires the
			//re-extracted gc-tagged corpus (contig_gc= in the header) -- the plain dedup/split
			//fastas built earlier this session do NOT carry it (see reextract_{dedup,split}_gc.sh).
			final int nf=11+(INCLUDE_CROSS_ENRICHMENT ? 3 : 0);
			out.println("#dims\t"+nf+"\t1");
			for(Read r : reads){
				final int lflank=TrnaNinemerTableBuilder.parseFlankValue(r.id, "lflank=");
				final int rflank=TrnaNinemerTableBuilder.parseFlankValue(r.id, "rflank=");
				if(lflank<0 || rflank<0){noFlankHeader++; continue;}
				final byte[] bases=r.bases;
				final int trueStart=lflank;
				final int trueStop=bases.length-rflank-1;
				if(trueStart<0 || trueStop>=bases.length || trueStop<=trueStart){malformedFlank++; continue;}
				final float contigGC=parseHeaderFloat(r.id, "contig_gc=");
				if(Float.isNaN(contigGC)){noContigGC++; continue;}

				final byte[] plusSeq=java.util.Arrays.copyOfRange(bases, trueStart, trueStop+1);
				final int[] shortlist=index.shortlist(plusSeq, 24);
				if(shortlist.length==0){skippedNoModel++; continue;}
				int bestModel=-1; float bestId=-1;
				for(int m : shortlist){
					float id=TrnaBoundaryFeatures.aniFeature(plusSeq, library[m]);
					if(id>bestId){bestId=id; bestModel=m;}
				}
				if(bestModel<0){skippedNoModel++; continue;}
				//Same model already picked for the ANI feature above -- library/GeneCaller.trnaModels
				//are index-aligned (both loaded from the SAME anticodon-cluster loop in
				//TrnaConsensusBuilder.process(), see CallGenes.loadTrnaResources:1562-1572).
				final BaseGraph model=(GeneCaller.trnaModels!=null && bestModel<GeneCaller.trnaModels.length
					? GeneCaller.trnaModels[bestModel] : null);

				written+=emitBoundaryVectors(out, bases, trueStart, trueStop, true, library[bestModel], startTable, startInside, startOutside, model, contigGC, stopTable, stopInside, stopOutside);
				written+=emitBoundaryVectors(out, bases, trueStart, trueStop, false, library[bestModel], stopTable, stopInside, stopOutside, model, contigGC, startTable, startInside, startOutside);
			}
		}
		System.err.println("Wrote "+written+" vectors. Skipped: "+noFlankHeader+" no-flank-header, "
			+malformedFlank+" malformed-flank, "+noContigGC+" no-contig-gc-header, "+skippedNoModel+" no-shortlist-model.");
	}

	/** Single-genome mode (original design): one genome + its reference GFF. */
	private static void runSingleGenomeMode(String fastaPath, String gffPath, String outPath,
			String tableStartPath, String tableStopPath, byte[][] library, TrnaKmerIndex index) throws IOException {
		Map<String,byte[]> contigs=readFasta(fastaPath);
		List<TrnaLocus> loci=readGff(gffPath);
		System.err.println("Loaded "+contigs.size()+" contigs, "+loci.size()+" reference tRNA loci.");

		final TrnaBoundaryFeatures.NinemerTable startTable, stopTable;
		final int startInside, startOutside, stopInside, stopOutside;
		if(tableStartPath!=null){
			final TrnaNinemerTableBuilder.LoadedTable lt1=TrnaNinemerTableBuilder.loadTable(tableStartPath);
			final TrnaNinemerTableBuilder.LoadedTable lt2=TrnaNinemerTableBuilder.loadTable(tableStopPath);
			assert(lt1.type==TrnaBoundaryFeatures.BoundaryType.START) : "tablestart= ("+tableStartPath
				+") is not labeled a start-boundary table -- wrong file, or tablestart=/tablestop= swapped.";
			assert(lt2.type==TrnaBoundaryFeatures.BoundaryType.STOP) : "tablestop= ("+tableStopPath
				+") is not labeled a stop-boundary table -- wrong file, or tablestart=/tablestop= swapped.";
			startTable=lt1.table; stopTable=lt2.table;
			startInside=lt1.insideCount; startOutside=lt1.outsideCount;
			stopInside=lt2.insideCount; stopOutside=lt2.outsideCount;
			System.err.println("Using REAL corpus tables: "+tableStartPath+" (inside="+startInside+" outside="+startOutside
				+") / "+tableStopPath+" (inside="+stopInside+" outside="+stopOutside+")");
		}else{
			//SMOKE-TEST 9-mer table: self-derived from this SAME genome's own true boundaries
			//only, ONE table shared for both boundary types. This exercises the
			//ninemerFeatures plumbing end-to-end -- it is NOT the real corpus-derived
			//enrichment signal and must never be mistaken for one. buildSmokeNinemerTable's
			//own build formula is boundaryPos-5 (fixed, doesn't go through kmerWindowStart at
			//all -- see its TODO-annotated known limitation). To keep the lookup side an
			//EXACT match for both boundary types, these (inside,outside) pairs are chosen so
			//TrnaBoundaryFeatures.kmerWindowStart reduces to that same boundaryPos-5 literal
			//for BOTH types -- historical-compatibility values, not a real design choice, and
			//not used anywhere but this smoke path (START: boundaryPos-outsideCount=
			//boundaryPos-5 -> outsideCount=5,insideCount=4; STOP: boundaryPos-(insideCount-1)=
			//boundaryPos-5 -> insideCount=6,outsideCount=3; both width 9).
			final TrnaBoundaryFeatures.NinemerTable smokeTable=buildSmokeNinemerTable(contigs, loci);
			startTable=smokeTable; stopTable=smokeTable;
			startInside=4; startOutside=5;
			stopInside=6; stopOutside=3;
			System.err.println("Using SMOKE-TEST table (NOT real training data).");
		}

		int written=0, skippedNoModel=0, skippedNoContig=0, skippedShort=0;
		final Map<String,Float> gcCache=new HashMap<>();
		try(PrintStream out=new PrintStream(new FileOutputStream(outPath))){
			//v3 feature set (Brian+Noire, Aug 21): stem(1) + ani(1) + enrichmentProfile(3) +
			//isStop(1) + tipFuzziness(3) + lengthRatio(1) + contigGC(1) = 11. Unlike flanked-batch
			//mode, contig-GC is computed directly from the full contig (already in hand here,
			//unlike the flanked-batch mode's small flanked snippet) via shared.Tools.calcGC --
			//same formula gff.CutGff's gccontig= flag uses -- cached per contig id since several
			//loci commonly share one contig.
			final int nf=11+(INCLUDE_CROSS_ENRICHMENT ? 3 : 0);
			out.println("#dims\t"+nf+"\t1");
			for(TrnaLocus locus : loci){
				byte[] contig=contigs.get(locus.contig);
				if(contig==null){skippedNoContig++; continue;}
				final float contigGC=gcCache.computeIfAbsent(locus.contig, k->shared.Tools.calcGC(contig));

				final int PAD=10;
				int winStart=Math.max(0, locus.start-PAD);
				int winStop=Math.min(contig.length-1, locus.stop+PAD);
				if(winStop-winStart<20){skippedShort++; continue;}
				byte[] window=java.util.Arrays.copyOfRange(contig, winStart, winStop+1);
				if(locus.strand=='-'){window=revcomp(window);}
				//In the reoriented (5'->3') frame: true start/stop local coordinates.
				final int trueStart=(locus.strand=='-' ? (winStop-locus.stop) : (locus.start-winStart));
				final int trueStop=(locus.strand=='-' ? (winStop-locus.start) : (locus.stop-winStart));
				if(trueStart<0 || trueStop>=window.length || trueStop<=trueStart){skippedShort++; continue;}

				byte[] plusSeq=java.util.Arrays.copyOfRange(window, trueStart, trueStop+1);
				int[] shortlist=index.shortlist(plusSeq, 24);
				if(shortlist.length==0){skippedNoModel++; continue;}
				int bestModel=-1; float bestId=-1;
				for(int m : shortlist){
					float id=TrnaBoundaryFeatures.aniFeature(plusSeq, library[m]);
					if(id>bestId){bestId=id; bestModel=m;}
				}
				if(bestModel<0){skippedNoModel++; continue;}
				final BaseGraph model=(GeneCaller.trnaModels!=null && bestModel<GeneCaller.trnaModels.length
					? GeneCaller.trnaModels[bestModel] : null);

				written+=emitBoundaryVectors(out, window, trueStart, trueStop, true, library[bestModel], startTable, startInside, startOutside, model, contigGC, stopTable, stopInside, stopOutside);
				written+=emitBoundaryVectors(out, window, trueStart, trueStop, false, library[bestModel], stopTable, stopInside, stopOutside, model, contigGC, startTable, startInside, startOutside);
			}
		}
		System.err.println("Wrote "+written+" vectors. Skipped: "+skippedNoContig+" no-contig, "
			+skippedShort+" too-short/malformed, "+skippedNoModel+" no-shortlist-model.");
	}

	/** Candidate search ranges (Brian, Aug 22 2026 -- boundary-offset histogram on the
	 * 203-genome bench AND independently on the full 2.01M-record flanked corpus, two
	 * measurements agreeing to within 0.1pp): START skews slightly negative (a real
	 * -3bp secondary peak), STOP has a large, sharp +3/+4bp secondary peak (~7.3% of
	 * loci on BOTH measurements) -- asymmetric per boundary, not the old symmetric
	 * +-2. Both are 6 candidates wide. Range chosen for ~99.3% coverage while keeping
	 * candidate-set size (and therefore false-positive opportunity) as small as the
	 * data allows -- see plans/boundary_offset_histogram_result_20260822.md for the
	 * full histograms and the coverage-vs-width tradeoff Brian weighed. */
	/** Package-private, not private -- TrnaBoundaryScorer.refineBoundaries' inference-time
	 * sweep MUST query the identical candidate set training scored, or the net gets
	 * evaluated on offsets it never saw a label for. Shared reference, not a duplicated
	 * literal, so the two sides structurally cannot drift apart (same principle as
	 * TrnaBoundaryFeatures.kmerWindowStart being factored out for the table builder).
	 *
	 * <p>STOP_OFFSETS DIRECTION CORRECTED (Neptune, Aug 22 2026, failure analysis --
	 * plans/failure_analysis_20260822.md): the refiner's own formula is
	 * `e=e0+offset`, applied relative to the CURRENT (already-imperfect) call, not
	 * relative to truth. STOP's dominant pre-refinement error is the base call
	 * sitting +3/+4 PAST true (the same bias that motivated widening the range in
	 * the first place) -- correcting that needs offset=-3/-4, which the original
	 * {-1,0,1,2,3,4} could never reach (only -1 negative). 86% of remaining STOP
	 * failures after the first with-refiner grading run were this exact,
	 * structurally-unreachable case, hand-verified on a real locus. Re-deriving the
	 * best contiguous 6-wide window from the CORRECTLY-signed (needed-correction)
	 * distribution gives {-4,-3,-2,-1,0,+1} at the same 99.30% aggregate coverage --
	 * same bins, opposite direction. START_OFFSETS was independently re-verified
	 * correct by the same method (its histogram convention already matches the
	 * formula's needed sign; no change). */
	static final int[] START_OFFSETS={-3,-2,-1,0,1,2};
	static final int[] STOP_OFFSETS={-4,-3,-2,-1,0,1};

	/** Cross-boundary enrichment feature (Brian, Aug 22 2026, queued directive #2):
	 * when scoring a START candidate, ALSO include the STOP boundary's enrichment
	 * profile at its OWN current position (and vice versa) -- 3 extra values, dims
	 * 11->14 when enabled. Gives the net direct cross-boundary context: "the other
	 * end is strongly/weakly enriched" calibrates how much to trust THIS end's
	 * palindrome score. DEFAULT false: dims stay 11, byte-identical to every vector
	 * file generated before this change -- toggleable so the STOP_OFFSETS-direction
	 * retrain in progress is untouched; flip on for a dedicated 14-dim regen. Package-
	 * private, shared with TrnaBoundaryScorer for the same drift-proofing reason as
	 * START_OFFSETS/STOP_OFFSETS above. */
	static boolean INCLUDE_CROSS_ENRICHMENT=false;

	/** Noisy-opposite-boundary training augmentation (Brian's key insight, Aug 22 2026,
	 * queued directive from the failure-analysis follow-up): training vectors are
	 * currently generated with the OPPOSITE boundary always exactly correct (from truth),
	 * but at real inference time CallGenes' own opposite-boundary call is wrong ~5.7% of
	 * the time for starts, ~10.6% for stops (203-genome-bench 5'/3'-non-exact rates,
	 * v1 no-refiner baseline) -- a genuine train/inference distribution mismatch on
	 * load-bearing features (stem pairing, length, ANI, and cross-boundary enrichment if
	 * that's also on). Fix: with the real error-rate probability, perturb the FIXED
	 * opposite boundary away from truth by a realistic magnitude/direction BEFORE
	 * generating this boundary's 6 candidate rows -- the VARIED boundary itself always
	 * stays anchored to its own true position (offset=0 must still mean "exactly
	 * correct," or the label is wrong).
	 *
	 * <p>Perturbation magnitude/direction is sampled from the REAL full-corpus offset
	 * histogram (Part 1, boundary_offset_histogram_result_20260822.md), non-zero bins
	 * only (the "given an error happens, how big/which way" distribution) -- NOT
	 * uniform, matching the actual observed bias shape (STOP overwhelmingly +3/+4,
	 * START skewed -3).
	 *
	 * <p>SIGN CONVENTION -- read this before touching either delta table (this exact
	 * class of mistake cost real rework earlier today, see the STOP_OFFSETS postmortem
	 * above): a perturbation delta is defined as calledPos-truePos (the thing you ADD
	 * to the true position to simulate a realistic called position). The full-corpus
	 * histogram's STOP convention (stopOff=calledStop-trueStop) IS this delta directly
	 * -- no sign flip. Its START convention (startOff=trueStart-calledStart) is the
	 * NEGATIVE of this delta -- MUST flip sign before using it as a perturbation
	 * (delta=calledStart-trueStart=-startOff). START_NOISE_DELTAS below is already the
	 * flipped (correct-for-this-use) array; do not re-derive from startOff without
	 * flipping again. */
	static boolean NOISY_OPPOSITE_BOUNDARY=false;
	/** Real 203-genome-bench, v1 no-refiner-baseline non-exact rates (1-5'exact, 1-3'exact).
	 * Overridable so the augmentation can be re-tuned without a recompile. */
	static float NOISY_START_RATE=0.057f;
	static float NOISY_STOP_RATE=0.106f;
	static long NOISY_SEED=1;
	private static java.util.Random noisyRandom=null;

	private static java.util.Random noisyRandom(){
		if(noisyRandom==null){noisyRandom=new java.util.Random(NOISY_SEED);}
		return noisyRandom;
	}

	//STOP delta=calledStop-trueStop directly (no flip -- see class comment above).
	//Bin counts are the REAL full-corpus histogram's non-zero bins (Aug 22 2026).
	private static final int[] STOP_NOISE_DELTAS={-4, -3, -2, -1, 1, 2, 3, 4};
	private static final long[] STOP_NOISE_WEIGHTS={1280, 3174, 1196, 5723, 20001, 6565, 125293, 20156};
	//START delta=calledStart-trueStart=-startOff (FLIPPED from the histogram's own
	//startOff=trueStart-calledStart convention -- see class comment above).
	private static final int[] START_NOISE_DELTAS={4, 3, 2, 1, -1, -2, -3, -4};
	private static final long[] START_NOISE_WEIGHTS={3183, 19498, 1890, 8944, 38643, 6845, 2683, 1092};

	/** Weighted-random pick from a parallel (delta,weight) table -- one draw per call,
	 * reused for both START and STOP noise sampling. */
	private static int sampleWeighted(int[] deltas, long[] weights, java.util.Random rnd){
		long total=0;
		for(long w : weights){total+=w;}
		long r=(long)(rnd.nextDouble()*total);
		long cum=0;
		for(int i=0; i<deltas.length; i++){
			cum+=weights[i];
			if(r<cum){return deltas[i];}
		}
		return deltas[deltas.length-1];//unreachable in practice (r<total by construction)
	}

	/** Emits 6 vectors (true + 5 shifted) for either the start (varyStart=true) or the
	 * stop (varyStart=false) boundary of one locus, holding the other boundary fixed.
	 * trueBoundaryPos (this locus's OWN unshifted true start/stop) is passed to
	 * enrichmentProfile's TRAINING overload for leave-one-out -- see its javadoc.
	 * otherTable/otherInside/otherOutside are the OPPOSITE boundary type's table
	 * (always available -- both start and stop tables are already loaded at both
	 * call sites), used only when INCLUDE_CROSS_ENRICHMENT is on. */
	private static int emitBoundaryVectors(PrintStream out, byte[] window, int trueStart, int trueStop,
			boolean varyStart, byte[] modelConsensus, TrnaBoundaryFeatures.NinemerTable table, int insideCount, int outsideCount,
			BaseGraph model, float contigGC, TrnaBoundaryFeatures.NinemerTable otherTable, int otherInside, int otherOutside){
		final TrnaBoundaryFeatures.BoundaryType type=(varyStart
			? TrnaBoundaryFeatures.BoundaryType.START : TrnaBoundaryFeatures.BoundaryType.STOP);
		final int trueBoundaryPos=(varyStart ? trueStart : trueStop);
		final int[] offsets=(varyStart ? START_OFFSETS : STOP_OFFSETS);

		//NOISY OPPOSITE BOUNDARY (see class-level comment for the full rationale and the
		//sign-convention warning): the FIXED opposite boundary used for the rest of this
		//call is, with realistic probability, perturbed away from truth. The VARIED
		//boundary (trueStart when varyStart, trueStop otherwise) is NEVER perturbed here --
		//it stays the exact anchor offset=0 is measured against, or the label breaks.
		int effectiveTrueStart=trueStart, effectiveTrueStop=trueStop;
		if(NOISY_OPPOSITE_BOUNDARY){
			if(varyStart){//opposite boundary is STOP
				if(noisyRandom().nextFloat()<NOISY_STOP_RATE){
					effectiveTrueStop=trueStop+sampleWeighted(STOP_NOISE_DELTAS, STOP_NOISE_WEIGHTS, noisyRandom());
				}
			}else{//opposite boundary is START
				if(noisyRandom().nextFloat()<NOISY_START_RATE){
					effectiveTrueStart=trueStart+sampleWeighted(START_NOISE_DELTAS, START_NOISE_WEIGHTS, noisyRandom());
				}
			}
		}

		//Cross-boundary enrichment (queued directive #2): the OTHER boundary stays fixed at
		//its EFFECTIVE (possibly noisy) position throughout this whole loop (only THIS
		//boundary varies), so its profile is the SAME for every offset -- compute it once
		//outside the loop, not per row. Reading the effective (not true) position here is
		//deliberate: if the opposite boundary is simulated-wrong, its enrichment profile
		//there should look however a genuinely-wrong position looks, giving the net real
		//signal that "the other end doesn't look right."
		final float[] otherProf;
		if(INCLUDE_CROSS_ENRICHMENT){
			final TrnaBoundaryFeatures.BoundaryType otherType=(varyStart
				? TrnaBoundaryFeatures.BoundaryType.STOP : TrnaBoundaryFeatures.BoundaryType.START);
			final int otherPos=(varyStart ? effectiveTrueStop : effectiveTrueStart);
			otherProf=TrnaBoundaryFeatures.enrichmentProfile(window, otherPos, otherPos,
				otherType, otherInside, otherOutside, otherTable);
		}else{
			otherProf=null;
		}
		int n=0;
		for(int offset : offsets){
			int s=(varyStart ? trueStart+offset : effectiveTrueStart);
			int e=(varyStart ? effectiveTrueStop : trueStop+offset);
			if(s<0 || e>=window.length || e-s<15){continue;}
			int label=(offset==0 ? 1 : 0);
			float stem=TrnaBoundaryFeatures.stemFeature(window, s, e);
			int boundaryPos=(varyStart ? s : e);
			float[] prof=TrnaBoundaryFeatures.enrichmentProfile(window, boundaryPos, trueBoundaryPos,
				type, insideCount, outsideCount, table);
			float isStop=(varyStart ? 0f : 1f);
			//LEAKAGE FIX (Brian, Aug 21 2026): ani and fuzziness must both be recomputed from THIS
			//candidate's own shifted span, never reused from the true span -- at real inference
			//time (CallGenes), the true boundary is unknown, so every feature has to be computable
			//from a candidate's own guessed span alone. Model CHOICE stays fixed per locus (matches
			//TrnaCaller's real behavior: pick a model once via the kmer shortlist, then refine
			//boundaries against that same model) -- only the ALIGNMENT is redone per candidate.
			//ani used to be computed ONCE by the caller from the true span and passed in unchanged
			//for all 10 rows -- a real leak, caught by Brian's identical-to-CallGenes-output check
			//before any training ran on it.
			final byte[] candSeq=java.util.Arrays.copyOfRange(window, s, e+1);
			float ani=TrnaBoundaryFeatures.aniFeature(candSeq, modelConsensus);
			float[] fuzz=TrnaBoundaryFeatures.tipFuzzinessFeature(candSeq, model, varyStart);
			float lengthRatio=TrnaBoundaryFeatures.lengthRatioFeature(s, e);
			if(INCLUDE_CROSS_ENRICHMENT){
				out.printf("%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d%n",
					stem, ani, prof[0], prof[1], prof[2], isStop, fuzz[0], fuzz[1], fuzz[2], lengthRatio, contigGC,
					otherProf[0], otherProf[1], otherProf[2], label);
			}else{
				out.printf("%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d%n",
					stem, ani, prof[0], prof[1], prof[2], isStop, fuzz[0], fuzz[1], fuzz[2], lengthRatio, contigGC, label);
			}
			n++;
		}
		return n;
	}

	/** Builds a 9-mer count table from ONLY this genome's own true boundaries -- a
	 * smoke-test stand-in, not the real corpus enrichment table (see class header).
	 * //TODO: Probable bug (found incidentally while refactoring the shared kmerWindowStart
	 * helper into TrnaBoundaryFeatures, Aug 21 2026, Noire's review): this reads RAW genomic
	 * contig bytes at locus.start/locus.stop directly, with no strand revcomp and no
	 * start/stop role-swap for minus-strand loci -- but the REAL lookup path
	 * (emitBoundaryVectors -> ninemerFeatures) queries the REORIENTED (revcomp'd-when-minus)
	 * window at trueStart/trueStop. For a minus-strand locus this table's key is the coding-
	 * strand sequence, the lookup's key is the reverse-complement -- different sequences, so
	 * this table under-counts/miskeys roughly the minus-strand half of its own loci. Not fixed
	 * here: this table is EXPLICITLY smoke-test-only (never trains anything real; the real
	 * corpus table is prok.TrnaNinemerTableBuilder, which DOES revcomp correctly -- verified,
	 * see its readGffLoci/window-reorientation code).
	 * UPDATE (Aug 21): this bug now CRASHES the smoke path loudly (AssertionError in
	 * enrichmentProfile's leave-one-out check, verified: reproduces on tid_1001994's
	 * minus-strand loci) instead of silently returning a wrong zero -- confirms
	 * crash-loud-never-wrong is working as intended, not a regression to fix. Smoke
	 * mode has no practical use now that real locked tables exist (see
	 * ~/Neptune/ninemer/locked_{start,stop}.tsv) -- not worth patching a path nothing
	 * real depends on. */
	private static TrnaBoundaryFeatures.NinemerTable buildSmokeNinemerTable(
			Map<String,byte[]> contigs, List<TrnaLocus> loci){
		final Map<String,Integer> counts=new HashMap<>();
		for(TrnaLocus locus : loci){
			byte[] contig=contigs.get(locus.contig);
			if(contig==null){continue;}
			for(int genomicBoundary : new int[]{locus.start, locus.stop}){
				int center=genomicBoundary-5;
				if(center<0 || center+9>contig.length){continue;}
				String kmer=new String(contig, center, 9);
				counts.merge(kmer, 1, Integer::sum);
			}
		}
		return (seq, pos)->{
			if(pos<0 || pos+9>seq.length){return 0;}
			String kmer=new String(seq, pos, 9);
			return counts.getOrDefault(kmer, 0);
		};
	}

	/** Parses a float-valued header field like "contig_gc=0.3267" (decimal point, unlike
	 * TrnaNinemerTableBuilder.parseFlankValue's integer-only digit scan). Returns NaN
	 * (not a negative sentinel -- the value is a fraction that could legitimately be
	 * near 0) if the key is absent or malformed. */
	private static float parseHeaderFloat(String header, String key){
		final int idx=header.indexOf(key);
		if(idx<0){return Float.NaN;}
		final int start=idx+key.length();
		int end=start;
		while(end<header.length() && (Character.isDigit(header.charAt(end)) || header.charAt(end)=='.')){end++;}
		if(end==start){return Float.NaN;}
		try{return Float.parseFloat(header.substring(start, end));}
		catch(NumberFormatException e){return Float.NaN;}
	}

	private static byte[] revcomp(byte[] seq){
		byte[] out=new byte[seq.length];
		for(int i=0; i<seq.length; i++){out[i]=AminoAcid.baseToComplementExtended[seq[seq.length-1-i]];}
		return out;
	}

	/*-- Minimal, self-contained FASTA/GFF readers (this is a standalone offline tool,
	     not part of the production I/O path) --*/

	static final class TrnaLocus {
		String contig; int start; int stop; char strand;//start/stop are 0-based, inclusive
	}

	private static Map<String,byte[]> readFasta(String path) throws IOException {
		Map<String,byte[]> out=new HashMap<>();
		StringBuilder seq=new StringBuilder();
		String name=null;
		try(BufferedReader br=openText(path)){
			String line;
			while((line=br.readLine())!=null){
				if(line.isEmpty()){continue;}
				if(line.charAt(0)=='>'){
					if(name!=null){out.put(name, seq.toString().getBytes());}
					int sp=line.indexOf(' ');
					name=line.substring(1, sp<0 ? line.length() : sp);
					seq.setLength(0);
				}else{
					seq.append(line.trim());
				}
			}
			if(name!=null){out.put(name, seq.toString().getBytes());}
		}
		return out;
	}

	private static List<TrnaLocus> readGff(String path) throws IOException {
		List<TrnaLocus> out=new ArrayList<>();
		try(BufferedReader br=openText(path)){
			String line;
			while((line=br.readLine())!=null){
				if(line.isEmpty() || line.charAt(0)=='#'){continue;}
				String[] f=line.split("\t");
				if(f.length<7 || !f[2].equals("tRNA")){continue;}
				TrnaLocus loc=new TrnaLocus();
				loc.contig=f[0];
				loc.start=Integer.parseInt(f[3])-1;//GFF is 1-based inclusive
				loc.stop=Integer.parseInt(f[4])-1;
				loc.strand=f[6].charAt(0);
				out.add(loc);
			}
		}
		return out;
	}

	private static BufferedReader openText(String path) throws IOException {
		InputStream is=new FileInputStream(path);
		if(path.endsWith(".gz")){is=new GZIPInputStream(is);}
		return new BufferedReader(new InputStreamReader(is));
	}
}
