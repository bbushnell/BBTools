package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import parse.LineParser1;
import shared.Shared;
import stream.Read;
import stream.ReadInputStream;
import structures.IntHashMap;

/**
 * Builds the real corpus-derived boundary 9-mer enrichment tables for the tRNA
 * boundary-precision NN (Brian's idea #2, trna repo plans/boundary_precision_ideas.md):
 * counts the 9-mer = 5bp OUTSIDE + 4bp INSIDE at every known tRNA START and every known
 * tRNA STOP, tabulated SEPARATELY (two tables, per the design doc), over a corpus of
 * genome+GFF pairs. TrnaBoundaryVectorGen's smoke-test table (self-derived from 1-2
 * genomes, one merged table for both boundary types) is explicitly NOT this.
 *
 * <p>Genomic-frame convention: each locus's window is reoriented to the tRNA's
 * 5'-&gt;3' reading frame (reverse-complemented as a whole when strand=='-'), matching
 * TrnaBoundaryVectorGen. The kmer window start is computed by
 * TrnaBoundaryFeatures.kmerWindowStart(boundaryPos, type, insideCount, outsideCount)
 * -- a SHARED, type-safe helper this class and TrnaBoundaryFeatures.ninemerFeatures
 * both call (hardened Aug 21, Noire's review: an earlier raw-signed-offset version
 * silently flipped its own inside/outside meaning between START and STOP; the
 * semantic form makes that structurally impossible -- see BoundaryType's javadoc).
 * LOCKED PRODUCTION VALUES (2.01M-genome dedup-corpus discrimination grid, Aug 21):
 * k=9, START insideCount=9/outsideCount=0 (fully inside the tRNA), STOP
 * insideCount=8/outsideCount=1 (one base into the downstream flank).
 *
 * <p>Fasta-&gt;gff pairing: same-directory, same-basename convention (verified live on
 * the bench5x corpus, archaea4/tid_X.fna.gz + archaea4/tid_X.gff.gz), matching
 * gff.AnnotateAnticodon.inferGff's rule exactly.
 *
 * <p>Three modes:
 * <p>1. list=&lt;fna_list.txt&gt; outstart=&lt;start_table.tsv&gt; outstop=&lt;stop_table.tsv&gt;
 *   startinside=N startoutside=N stopinside=N stopoutside=N
 *   Genome-walk mode: fna_list.txt is one .fna[.gz] path per line (e.g.
 *   results/bench5x_list.txt); pairs each with its GFF and walks every tRNA locus.
 *   Used for the 934-genome holdout preview (subset, not full-corpus distribution).
 * <p>2. fasta=&lt;flanked.fa[.gz]&gt; outstart=&lt;start_table.tsv&gt; outstop=&lt;stop_table.tsv&gt;
 *   startinside=N startoutside=N stopinside=N stopoutside=N
 *   Flanked-record mode: reads a cutgff flank=/pad= output fasta directly (Brian's pivot,
 *   Aug 21 -- one artifact serves both the tRNA junction DB and this table, full 2.7M-tRNA
 *   distribution, no genome re-read). Each header carries " lflank=N rflank=M" (see
 *   gff.CutGff.java line ~577) recording the ACTUAL applied flank (may be less than
 *   requested if clipped at a contig edge) on a sequence ALREADY reoriented to 5'-&gt;3' by
 *   cutgff (no revcomp needed here). START boundary local offset = lflank; STOP boundary
 *   local offset (last included base) = length-rflank-1.
 *   startinside=/startoutside=/stopinside=/stopoutside= are ALL required (no default) --
 *   an explicit window per boundary type, never a raw shared offset.
 * <p>3. fasta=&lt;flanked.fa[.gz]&gt; ks=&lt;comma-separated k values, e.g. 7,8,9,10&gt;
 *   outdiscrim=&lt;report.tsv&gt; [ms=1.5,2,3,5,10]
 *   Grid discrimination-sweep mode (Noire+Brian, Aug 21): sweeps k AND a RAW offset
 *   (out in 0..k, deliberately unsafe/exploratory -- this mode exists BECAUSE the
 *   safe composition wasn't known yet; see runGridDiscriminationMode's javadoc) and
 *   reports a TPR/FPR decisiveness curve at each threshold m, not plain argmax. Reads
 *   the fasta once, sweeps the whole grid against the same in-memory records. This
 *   mode already answered the design question (see LOCKED PRODUCTION VALUES above) --
 *   only rerun it if revisiting the window choice itself.
 *
 * @author Neptune
 */
public class TrnaNinemerTableBuilder {

	public static void main(String[] args){
		String listPath=null, flankedFastaPath=null, outStartPath=null, outStopPath=null, ksArg=null, msArg=null, outDiscrimPath=null;
		Integer startInside=null, startOutside=null, stopInside=null, stopOutside=null;
		for(String a : args){
			String[] kv=a.split("=", 2);
			if(kv.length<2){continue;}
			if(kv[0].equalsIgnoreCase("list")){listPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("fasta")){flankedFastaPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("outstart")){outStartPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("outstop")){outStopPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("startinside")){startInside=Integer.parseInt(kv[1]);}
			else if(kv[0].equalsIgnoreCase("startoutside")){startOutside=Integer.parseInt(kv[1]);}
			else if(kv[0].equalsIgnoreCase("stopinside")){stopInside=Integer.parseInt(kv[1]);}
			else if(kv[0].equalsIgnoreCase("stopoutside")){stopOutside=Integer.parseInt(kv[1]);}
			else if(kv[0].equalsIgnoreCase("ks")){ksArg=kv[1];}
			else if(kv[0].equalsIgnoreCase("ms")){msArg=kv[1];}
			else if(kv[0].equalsIgnoreCase("outdiscrim")){outDiscrimPath=kv[1];}
		}

		if(ksArg!=null){
			//Grid discrimination-sweep mode (Noire+Brian, Aug 21): flanked-fasta only.
			if(flankedFastaPath==null || outDiscrimPath==null){
				System.err.println("Usage (grid discrimination): fasta=<flanked.fa[.gz]> ks=<comma-separated k values, e.g. 7,8,9,10> outdiscrim=<report.tsv> [ms=1.5,2,3,5,10]");
				System.exit(1);
			}
			final String[] parts=ksArg.split(",");
			final int[] ks=new int[parts.length];
			for(int i=0; i<parts.length; i++){
				ks[i]=Integer.parseInt(parts[i].trim());
				assert(ks[i]>=1 && ks[i]<=15) : "k out of supported range: "+ks[i];
			}
			final double[] ms;
			if(msArg!=null){
				final String[] mparts=msArg.split(",");
				ms=new double[mparts.length];
				for(int i=0; i<mparts.length; i++){
					ms[i]=Double.parseDouble(mparts[i].trim());
					assert(ms[i]>1.0) : "Threshold m must be >1 (m<=1 makes TPR/FPR non-exclusive): "+ms[i];
				}
			}else{
				ms=new double[]{1.5, 2, 3, 5, 10};
			}
			runGridDiscriminationMode(flankedFastaPath, ks, ms, outDiscrimPath);
			return;
		}

		final boolean haveAllCounts=(startInside!=null && startOutside!=null && stopInside!=null && stopOutside!=null);
		final boolean haveAnyCounts=(startInside!=null || startOutside!=null || stopInside!=null || stopOutside!=null);
		if(outStartPath==null || outStopPath==null || (listPath==null)==(flankedFastaPath==null) || !haveAllCounts){
			System.err.println("Usage (genome-walk):     list=<fna_list.txt> outstart=<start_table.tsv> outstop=<stop_table.tsv> startinside=9 startoutside=0 stopinside=8 stopoutside=1");
			System.err.println("Usage (flanked-record):  fasta=<flanked.fa[.gz]> outstart=<start_table.tsv> outstop=<stop_table.tsv> startinside=9 startoutside=0 stopinside=8 stopoutside=1");
			System.err.println("Usage (discrimination):  fasta=<flanked.fa[.gz]> ks=<comma-separated k values> outdiscrim=<report.tsv>");
			System.err.println("Exactly one of list= / fasta= is required (for the table-build usages).");
			System.err.println("startinside=/startoutside=/stopinside=/stopoutside= are ALL required together (no default -- "
				+"the whole point is an explicit, unambiguous window per boundary type; see TrnaBoundaryFeatures.BoundaryType).");
			System.exit(1);
		}
		assert(haveAnyCounts) : "unreachable given haveAllCounts check above";
		assert(startInside>=0 && startOutside>=0 && startInside+startOutside>=1 && startInside+startOutside<=15)
			: "start window size out of supported range: inside="+startInside+" outside="+startOutside;
		assert(stopInside>=0 && stopOutside>=0 && stopInside+stopOutside>=1 && stopInside+stopOutside<=15)
			: "stop window size out of supported range: inside="+stopInside+" outside="+stopOutside;

		if(flankedFastaPath!=null){
			runFlankedFastaMode(flankedFastaPath, outStartPath, outStopPath, startInside, startOutside, stopInside, stopOutside);
		}else{
			runGenomeWalkMode(listPath, outStartPath, outStopPath, startInside, startOutside, stopInside, stopOutside);
		}
	}

	/**
	 * Grid discrimination-sweep mode (Noire+Brian, Aug 21, supersedes the original
	 * single-argmax discrimination mode). For each k in ks and each offset
	 * out in {0..k}, builds the boundary tables at that (k,out) window THEN
	 * re-evaluates the SAME (already-loaded, single fasta read) records with a
	 * DECISIVENESS metric instead of plain argmax:
	 *
	 * <p>Brian's correction (Aug 21) to the earlier k-only sweep: varying k at a
	 * FIXED offset (out=5, the original design) confounds width with inside/outside
	 * balance -- growing k that way adds bases INSIDE the tRNA for a START boundary
	 * but OUTSIDE (into the flank) for STOP, since inside/outside flip direction
	 * depending on which end of the tRNA boundaryPos marks. So "START wants bigger
	 * k, STOP wants smaller k" from that sweep might have been an offset preference
	 * misread as a width preference. This mode varies BOTH dimensions independently.
	 *
	 * <p>Brian's second correction (Aug 21): plain argmax treats a decisive win
	 * (trueCount=972 vs neighbor=125) identically to a coin-flip tie (123 vs 125) --
	 * both just count as "+1". What the feature actually needs is DECISIVENESS: at
	 * threshold m, TPR(m) = fraction where the true position's count is >=m times
	 * every neighbor's (safe to act on); FPR(m) = fraction where some neighbor's
	 * count is >=m times the true position's (would decisively drag the boundary to
	 * the WRONG place -- the harmful case). These are NOT complementary; the gap is
	 * "ambiguous, feature abstains." Both use the SAME +1-smoothed counts as
	 * elsewhere, with leave-one-out on the true position only (see the class's
	 * general leave-one-out rationale, unchanged from the single-k mode this
	 * replaces: this record's own true kmer was counted into its own table during
	 * the k-build, so excluding it is required, not optional, or every locus
	 * trivially "recognizes itself").
	 *
	 * <p>Flanked-fasta only -- reads the corpus ONCE and reuses the in-memory Read
	 * list across the whole (k,out) grid (cheap compared to the one fasta read).
	 */
	private static void runGridDiscriminationMode(String fastaPath, int[] ks, double[] ms, String outPath){
		Shared.TRIM_READ_DESCRIPTION=false;
		Shared.TRIM_RNAME=true;
		Read.TO_UPPER_CASE=true;

		final ArrayList<Read> reads=ReadInputStream.toReads(fastaPath, FileFormat.FA, -1);
		assert(!reads.isEmpty()) : "Empty or unreadable flanked fasta: "+fastaPath;
		final int n=reads.size();
		System.err.println("Loaded "+n+" records from "+fastaPath);

		//Pre-parse flank values ONCE (shared across the whole (k,out) grid) -- header
		//parsing doesn't depend on k or out, only the kmer windowing does.
		final int[] trueStarts=new int[n];
		final int[] trueStops=new int[n];//-1 in either array marks an ineligible record
		long noFlankHeader=0, malformedFlank=0;
		for(int i=0; i<n; i++){
			final Read r=reads.get(i);
			final int lf=parseFlankValue(r.id, "lflank=");
			final int rf=parseFlankValue(r.id, "rflank=");
			if(lf<0 || rf<0){trueStarts[i]=-1; trueStops[i]=-1; noFlankHeader++; continue;}
			final int ts=lf, te=r.bases.length-rf-1;
			if(ts<0 || te>=r.bases.length || te<=ts){trueStarts[i]=-1; trueStops[i]=-1; malformedFlank++; continue;}
			trueStarts[i]=ts; trueStops[i]=te;
		}
		System.err.println("No-flank-header: "+noFlankHeader+"  Malformed: "+malformedFlank);

		try(PrintStream out=new PrintStream(new java.io.FileOutputStream(outPath))){
			out.print("k\tout\tboundary\teligible");
			for(final double m : ms){out.print("\tTPR("+m+")\tFPR("+m+")");}
			out.println();

			for(final int k : ks){
				for(int o=0; o<=k; o++){
					final int out_=o;//effectively-final for the lambda-free inner loops below
					final IntHashMap startCounts=new IntHashMap(1<<16);
					final IntHashMap stopCounts=new IntHashMap(1<<16);
					for(int i=0; i<n; i++){
						if(trueStarts[i]<0){continue;}
						final byte[] bases=reads.get(i).bases;
						final int sk=packKmer(bases, gridWindowStart(trueStarts[i], k, out_), k);
						if(sk>=0){startCounts.increment(sk);}
						final int ek=packKmer(bases, gridWindowStart(trueStops[i], k, out_), k);
						if(ek>=0){stopCounts.increment(ek);}
					}

					final long[] startEligible={0}, stopEligible={0};
					final long[] startTpr=new long[ms.length], startFpr=new long[ms.length];
					final long[] stopTpr=new long[ms.length], stopFpr=new long[ms.length];
					for(int i=0; i<n; i++){
						if(trueStarts[i]<0){continue;}
						final byte[] bases=reads.get(i).bases;
						accumulateDecisiveness(bases, trueStarts[i], k, out_, startCounts, ms, startEligible, startTpr, startFpr);
						accumulateDecisiveness(bases, trueStops[i], k, out_, stopCounts, ms, stopEligible, stopTpr, stopFpr);
					}

					writeGridRow(out, k, out_, "start", startEligible[0], ms, startTpr, startFpr);
					writeGridRow(out, k, out_, "stop", stopEligible[0], ms, stopTpr, stopFpr);
				}
				System.err.println("k="+k+" done (out=0.."+k+")");
			}
		}catch(java.io.IOException e){
			throw new RuntimeException("Failed writing "+outPath, e);
		}
	}

	/**
	 * RAW signed-offset window formula, PRIVATE and scoped to the grid-discrimination
	 * sweep only -- deliberately explores both flip directions as a raw offset (that's
	 * the whole point of the sweep: it doesn't yet know which composition is right).
	 * This is NOT the shared production formula (TrnaBoundaryFeatures.kmerWindowStart,
	 * which takes semantic insideCount/outsideCount+BoundaryType and can't have this
	 * bug) -- do not resurrect this as a public/shared method; that's exactly the
	 * ambiguity that got caught and killed at the source (Noire, Aug 21). window=
	 * [pos-out, pos-out+k).
	 */
	private static int gridWindowStart(int boundaryPos, int k, int out){
		return boundaryPos-out;
	}

	private static void writeGridRow(PrintStream out, int k, int o, String boundary, long eligible,
			double[] ms, long[] tprHits, long[] fprHits){
		out.print(k+"\t"+o+"\t"+boundary+"\t"+eligible);
		for(int i=0; i<ms.length; i++){
			out.print("\t"+pct(tprHits[i],eligible)+"\t"+pct(fprHits[i],eligible));
		}
		out.println();
	}

	/**
	 * -1-sentinel via eligible[0] increment guard: if the TRUE boundary's own kmer is
	 * out of bounds/ambiguous, the locus contributes nothing (not eligible). Otherwise
	 * increments eligible[0] and, for every threshold m, increments tprHits[i] if the
	 * true position's smoothed leave-one-out count is >=m times every neighbor's
	 * smoothed count, or fprHits[i] if some neighbor's smoothed count is >=m times the
	 * true position's. For m>1 these are mutually exclusive per locus per m (a locus
	 * can't satisfy both c0>=m*max(other) and max(other)>=m*c0 simultaneously unless
	 * m<=1, which is asserted against at the CLI level).
	 *
	 * <p>Smoothing (Brian, Aug 21): the RATIO uses (count+1) for every term, true and
	 * neighbors alike, so a 0-count neighbor can't make a ratio blow up or divide by
	 * zero. Leave-one-out is applied to the true position's RAW count first (rawCount-1,
	 * since this record contributed exactly +1 to its own true kmer during the k-build),
	 * then +1-smoothed -- net effect: smoothedTrue=rawCount (algebraically, (rawCount-1)+1).
	 */
	private static void accumulateDecisiveness(byte[] seq, int boundaryPos, int k, int out, IntHashMap table,
			double[] ms, long[] eligible, long[] tprHits, long[] fprHits){
		final int trueKmer=packKmer(seq, gridWindowStart(boundaryPos, k, out), k);
		if(trueKmer<0){return;}
		final int trueRawCount=table.get(trueKmer);
		assert(trueRawCount>=1) : "This record's own true kmer must have been counted during the k-build "
			+"(same boundaryPos/k/out) -- table.get returned "+trueRawCount+" for kmer "+trueKmer
			+", meaning the build pass and this pass disagree on this record's true-position kmer.";
		final double smoothedTrue=trueRawCount;//leave-one-out (rawCount-1) then +1-smoothed = rawCount

		double maxSmoothedNeighbor=0;
		for(final int offset : new int[]{-1, 1, -2, 2}){
			final int neighborKmer=packKmer(seq, gridWindowStart(boundaryPos+offset, k, out), k);
			final int c=(neighborKmer<0 ? -1 : table.get(neighborKmer));
			final double smoothedNeighbor=(c<0 ? 0 : c)+1.0;//absent/invalid -> count=0 -> smoothed=1
			if(smoothedNeighbor>maxSmoothedNeighbor){maxSmoothedNeighbor=smoothedNeighbor;}
		}

		eligible[0]++;
		for(int i=0; i<ms.length; i++){
			final double m=ms[i];
			if(smoothedTrue>=m*maxSmoothedNeighbor){tprHits[i]++;}
			if(maxSmoothedNeighbor>=m*smoothedTrue){fprHits[i]++;}
		}
	}

	/** Flanked-record mode (Brian's pivot, Aug 21): reads a cutgff flank= output fasta
	 * directly -- see class header for the header-parsing convention. Each boundary
	 * type gets its OWN insideCount/outsideCount (hardened API, Noire, Aug 21 -- see
	 * TrnaBoundaryFeatures.BoundaryType for why a single shared "k+offset" can't
	 * express this safely). */
	private static void runFlankedFastaMode(String fastaPath, String outStartPath, String outStopPath,
			int startInside, int startOutside, int stopInside, int stopOutside){
		//Deliberately NOT TRIM_READ_DESCRIPTION -- we need the full header (past the first
		//token) to read lflank=/rflank=. TRIM_RNAME is unrelated (SAM-only) but set for
		//consistency with the genome-walk mode.
		Shared.TRIM_READ_DESCRIPTION=false;
		Shared.TRIM_RNAME=true;
		Read.TO_UPPER_CASE=true;

		final int startK=startInside+startOutside, stopK=stopInside+stopOutside;
		final IntHashMap startCounts=new IntHashMap(1<<16);
		final IntHashMap stopCounts=new IntHashMap(1<<16);

		long recordsTotal=0, noFlankHeader=0, malformedFlank=0;
		long startCounted=0, stopCounted=0, edgeSkippedBoundary=0, ambigSkipped=0;

		final ArrayList<Read> reads=ReadInputStream.toReads(fastaPath, FileFormat.FA, -1);
		assert(!reads.isEmpty()) : "Empty or unreadable flanked fasta: "+fastaPath;

		for(Read r : reads){
			recordsTotal++;
			final int lflank=parseFlankValue(r.id, "lflank=");
			final int rflank=parseFlankValue(r.id, "rflank=");
			if(lflank<0 || rflank<0){
				//Loud, not silent: a record with no flank header can't be used for this table.
				noFlankHeader++;
				continue;
			}
			final byte[] bases=r.bases;
			final int trueStart=lflank;
			final int trueStop=bases.length-rflank-1;
			if(trueStart<0 || trueStop>=bases.length || trueStop<=trueStart){malformedFlank++; continue;}

			final int startPos=TrnaBoundaryFeatures.kmerWindowStart(trueStart,
				TrnaBoundaryFeatures.BoundaryType.START, startInside, startOutside);
			final int startKmer=packKmer(bases, startPos, startK);
			if(startKmer<0){
				if(startKmer==EDGE){edgeSkippedBoundary++;} else {assert(startKmer==AMBIG); ambigSkipped++;}
			}else{
				startCounts.increment(startKmer);
				startCounted++;
			}

			final int stopPos=TrnaBoundaryFeatures.kmerWindowStart(trueStop,
				TrnaBoundaryFeatures.BoundaryType.STOP, stopInside, stopOutside);
			final int stopKmer=packKmer(bases, stopPos, stopK);
			if(stopKmer<0){
				if(stopKmer==EDGE){edgeSkippedBoundary++;} else {assert(stopKmer==AMBIG); ambigSkipped++;}
			}else{
				stopCounts.increment(stopKmer);
				stopCounted++;
			}
		}

		assert(startCounted+stopCounted+2*noFlankHeader+2*malformedFlank+edgeSkippedBoundary+ambigSkipped==2*recordsTotal)
			: "Boundary observations don't exactly account for 2x every record (start+stop) -- a boundary was "
			+"silently dropped somewhere. recordsTotal="+recordsTotal+" startCounted="+startCounted+" stopCounted="+stopCounted
			+" noFlankHeader="+noFlankHeader+" malformedFlank="+malformedFlank
			+" edgeSkippedBoundary="+edgeSkippedBoundary+" ambigSkipped="+ambigSkipped;

		writeTable(outStartPath, startCounts, "start", startInside, startOutside);
		writeTable(outStopPath, stopCounts, "stop", stopInside, stopOutside);

		System.err.println("=== TrnaNinemerTableBuilder flanked-record mode ===");
		System.err.println("START window: inside="+startInside+" outside="+startOutside+" (k="+startK+")");
		System.err.println("STOP window:  inside="+stopInside+" outside="+stopOutside+" (k="+stopK+")");
		System.err.println("Records total:       "+recordsTotal);
		System.err.println("No-flank-header:     "+noFlankHeader);
		System.err.println("Malformed flank:     "+malformedFlank);
		System.err.println("Start observations:  "+startCounted+" ("+pct(startCounted,recordsTotal)+"%)");
		System.err.println("Stop observations:   "+stopCounted+" ("+pct(stopCounted,recordsTotal)+"%)");
		System.err.println("Edge-skipped:        "+edgeSkippedBoundary);
		System.err.println("Ambiguous-base-skip: "+ambigSkipped);
		printOccupancyStats("start", startCounts, startK);
		printOccupancyStats("stop", stopCounts, stopK);
	}

	/** Parses an integer value out of a "key=NNN" substring in a fasta header. Returns -1
	 * if the key isn't present or has no digits following it. Package-private (not
	 * private): TrnaBoundaryVectorGen's flanked-batch mode reuses this exact parser
	 * rather than duplicating it, so both readers can never drift on the header format. */
	static int parseFlankValue(String header, String key){
		final int idx=header.indexOf(key);
		if(idx<0){return -1;}
		final int start=idx+key.length();
		int end=start;
		while(end<header.length() && Character.isDigit(header.charAt(end))){end++;}
		if(end==start){return -1;}
		return Integer.parseInt(header.substring(start, end));
	}

	/** Genome-walk mode (original design): pairs each genome with its GFF and walks every
	 * tRNA locus directly. Preview-only for the boundary-precision NN as of Brian's Aug 21
	 * pivot to the flanked-record corpus (see class header) -- still useful standalone.
	 * Each boundary type gets its OWN insideCount/outsideCount, same hardened API as
	 * runFlankedFastaMode. */
	private static void runGenomeWalkMode(String listPath, String outStartPath, String outStopPath,
			int startInside, int startOutside, int stopInside, int stopOutside){
		//AnnotateAnticodon sets these before parsing genome FASTAs so Read.id is the first
		//header token only; without this, a multi-word FASTA header makes r.id the FULL
		//header (with description), which never matches a GFF seqid -- exactly the bug this
		//tool hit on its first cluster run (9,390/45,964 loci with a spurious contig-missing).
		Shared.TRIM_READ_DESCRIPTION=Shared.TRIM_RNAME=true;
		Read.TO_UPPER_CASE=true;

		final int startK=startInside+startOutside, stopK=stopInside+stopOutside;
		final List<String> fnaPaths=readList(listPath);
		assert(!fnaPaths.isEmpty()) : "Empty or unreadable genome list: "+listPath;

		//9-mer over {A,C,G,T} packs into 18 bits (2 bits/base); always in [0,262144), never
		//negative -- safe as an IntHashMap key with no reserved-sentinel collision risk.
		final IntHashMap startCounts=new IntHashMap(1<<14);
		final IntHashMap stopCounts=new IntHashMap(1<<14);

		long genomesOk=0, genomesNoGff=0, lociTotal=0;
		long startCounted=0, stopCounted=0;
		//edgeSkippedLocus: locus-level skip, before either boundary is evaluated -- costs BOTH
		//boundary-observations at once. edgeSkippedBoundary: packKmer's own EDGE return, one
		//boundary at a time. Kept separate (not merged into one edgeSkipped) so the accounting
		//assert below can be an exact equality, not just >=, per Noire's review (Aug 21).
		long edgeSkippedLocus=0, edgeSkippedBoundary=0, ambigSkipped=0, contigMissing=0, malformedLocus=0;

		for(String fna : fnaPaths){
			final String gff=inferGff(fna);
			if(gff==null){
				genomesNoGff++;
				System.err.println("SKIP (no paired gff found): "+fna);
				continue;
			}
			final HashMap<String,byte[]> genome=loadFasta(fna);
			final List<Locus> loci=readGffLoci(gff);
			genomesOk++;

			for(Locus loc : loci){
				lociTotal++;
				byte[] contig=genome.get(loc.contig);
				if(contig==null){
					//AnnotateAnticodon already resolved this exact corpus's seqids during Stage 1
					//(x4/bench5x genomes) -- a miss here on a genome that loaded is a real
					//data-integrity surprise, not a benign edge case, so it's counted loudly,
					//not silently dropped.
					contigMissing++;
					continue;
				}
				if(loc.start<0 || loc.stop<loc.start || loc.stop>=contig.length){
					malformedLocus++;
					continue;
				}

				final int PAD=10;//matches TrnaBoundaryVectorGen's PAD; covers the +-5 window margin
				final int winStart=Math.max(0, loc.start-PAD);
				final int winStop=Math.min(contig.length-1, loc.stop+PAD);
				if(winStop-winStart<20){edgeSkippedLocus++; continue;}
				byte[] window=java.util.Arrays.copyOfRange(contig, winStart, winStop+1);
				if(loc.strand=='-'){window=revcomp(window);}
				final int trueStart=(loc.strand=='-' ? (winStop-loc.stop) : (loc.start-winStart));
				final int trueStop=(loc.strand=='-' ? (winStop-loc.start) : (loc.stop-winStart));
				if(trueStart<0 || trueStop>=window.length || trueStop<=trueStart){edgeSkippedLocus++; continue;}

				final int startPos=TrnaBoundaryFeatures.kmerWindowStart(trueStart,
					TrnaBoundaryFeatures.BoundaryType.START, startInside, startOutside);
				final int startKmer=packKmer(window, startPos, startK);
				if(startKmer<0){
					if(startKmer==EDGE){edgeSkippedBoundary++;} else {assert(startKmer==AMBIG); ambigSkipped++;}
				}else{
					startCounts.increment(startKmer);
					startCounted++;
				}

				final int stopPos=TrnaBoundaryFeatures.kmerWindowStart(trueStop,
					TrnaBoundaryFeatures.BoundaryType.STOP, stopInside, stopOutside);
				final int stopKmer=packKmer(window, stopPos, stopK);
				if(stopKmer<0){
					if(stopKmer==EDGE){edgeSkippedBoundary++;} else {assert(stopKmer==AMBIG); ambigSkipped++;}
				}else{
					stopCounts.increment(stopKmer);
					stopCounted++;
				}
			}
		}

		final long edgeSkipped=edgeSkippedLocus+edgeSkippedBoundary;//for the summary report only
		assert(startCounted+stopCounted+2*edgeSkippedLocus+edgeSkippedBoundary+ambigSkipped+2*contigMissing+2*malformedLocus==2*lociTotal)
			: "Boundary observations don't exactly account for 2x every locus (start+stop) -- a boundary was "
			+"silently dropped somewhere. lociTotal="+lociTotal+" startCounted="+startCounted+" stopCounted="+stopCounted
			+" edgeSkippedLocus="+edgeSkippedLocus+" edgeSkippedBoundary="+edgeSkippedBoundary
			+" ambigSkipped="+ambigSkipped+" contigMissing="+contigMissing+" malformedLocus="+malformedLocus;

		writeTable(outStartPath, startCounts, "start", startInside, startOutside);
		writeTable(outStopPath, stopCounts, "stop", stopInside, stopOutside);

		System.err.println("=== TrnaNinemerTableBuilder genome-walk mode ===");
		System.err.println("START window: inside="+startInside+" outside="+startOutside+" (k="+startK+")");
		System.err.println("STOP window:  inside="+stopInside+" outside="+stopOutside+" (k="+stopK+")");
		System.err.println("Genomes OK:          "+genomesOk);
		System.err.println("Genomes no-gff:      "+genomesNoGff);
		System.err.println("Loci total:          "+lociTotal);
		System.err.println("Start observations:  "+startCounted+" ("+pct(startCounted,lociTotal)+"%)");
		System.err.println("Stop observations:   "+stopCounted+" ("+pct(stopCounted,lociTotal)+"%)");
		System.err.println("Edge-skipped:        "+edgeSkipped);
		System.err.println("Ambiguous-base-skip: "+ambigSkipped);
		System.err.println("Contig missing:      "+contigMissing);
		System.err.println("Malformed locus:     "+malformedLocus);
		printOccupancyStats("start", startCounts, startK);
		printOccupancyStats("stop", stopCounts, stopK);
	}

	private static double pct(long num, long denom){return denom==0 ? 0 : (100.0*num/denom);}

	/** Occupancy/sparsity report: distinct kmers observed, median/max count, fraction of
	 * boundary observations landing on a kmer seen >=5x (Noire's sparsity gate, Aug 21). */
	private static void printOccupancyStats(String label, IntHashMap counts, int k){
		final long possible=1L<<(2*k);
		final int[] keys=counts.toKeyArray();
		final int distinct=keys.length;
		long totalObs=0, obsOnFivePlus=0;
		final long[] countVals=new long[distinct];
		for(int i=0; i<distinct; i++){
			int c=counts.get(keys[i]);
			assert(c>0) : "Non-positive count for a present key -- IntHashMap invariant violated: "+c;
			countVals[i]=c;
			totalObs+=c;
			if(c>=5){obsOnFivePlus+=c;}
		}
		java.util.Arrays.sort(countVals);
		final long median=(distinct==0 ? 0 : countVals[distinct/2]);
		final long max=(distinct==0 ? 0 : countVals[distinct-1]);
		System.err.println("--- "+label+" table occupancy (k="+k+") ---");
		System.err.println("Distinct "+k+"-mers observed: "+distinct+" / "+possible+" possible ("+pct(distinct,possible)+"%)");
		System.err.println("Total observations:       "+totalObs);
		System.err.println("Median count (per distinct "+k+"-mer): "+median);
		System.err.println("Max count:                 "+max);
		System.err.println("Fraction of observations on a "+k+"-mer seen >=5x: "+pct(obsOnFivePlus,totalObs)+"%");
	}

	private static final int EDGE=-1, AMBIG=-2;

	/** Packs the k-mer at seq[pos,pos+k) into a 2-bits/base int (A=0 C=1 G=2 T=3).
	 * Returns EDGE if the window is out of bounds, AMBIG if any base isn't A/C/G/T. */
	private static int packKmer(byte[] seq, int pos, int k){
		if(pos<0 || pos+k>seq.length){return EDGE;}
		int kmer=0;
		for(int i=0; i<k; i++){
			byte b=AminoAcid.baseToNumber[seq[pos+i]];
			if(b<0){return AMBIG;}
			kmer=(kmer<<2)|b;
		}
		assert(kmer>=0 && kmer<(1<<(2*k))) : "Packed "+k+"-mer out of range: "+kmer;
		return kmer;
	}

	/** A table loaded back from a TSV written by writeTable, for real (non-smoke)
	 * training-vector generation. Exposes boundaryType/insideCount/outsideCount (not
	 * just k) so the CALLER can pass EXACTLY what built this table to
	 * TrnaBoundaryFeatures.ninemerFeatures -- the table itself is the single source of
	 * truth for its own window shape, so there is no separate value for a caller to
	 * mis-specify (Noire's "kill the ambiguity at the source" hardening, Aug 21). */
	public static final class LoadedTable {
		public final int k;
		public final TrnaBoundaryFeatures.BoundaryType type;
		public final int insideCount, outsideCount;
		public final TrnaBoundaryFeatures.NinemerTable table;
		LoadedTable(int k, TrnaBoundaryFeatures.BoundaryType type, int insideCount, int outsideCount,
				TrnaBoundaryFeatures.NinemerTable table){
			this.k=k; this.type=type; this.insideCount=insideCount; this.outsideCount=outsideCount; this.table=table;
		}
	}

	/** Loads a table TSV (writeTable's format: "#boundary/#insidecount/#outsidecount/#k"
	 * headers + "kmer\tcount" rows) back into a queryable NinemerTable. Absent kmers
	 * read as count=0 (log2(1)=0), matching the same convention used everywhere else
	 * in this class. */
	public static LoadedTable loadTable(String path){
		final IntHashMap counts=new IntHashMap(1<<16);
		int insideCount=-1, outsideCount=-1;
		TrnaBoundaryFeatures.BoundaryType type=null;
		final ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(path, FileFormat.TEXT, null, false, true));
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			final String s=new String(line);
			if(s.startsWith("#boundary\t")){
				final String label=s.substring(10).trim();
				type=(label.equals("start") ? TrnaBoundaryFeatures.BoundaryType.START
					: label.equals("stop") ? TrnaBoundaryFeatures.BoundaryType.STOP : null);
				assert(type!=null) : "Unrecognized #boundary label in "+path+": "+label;
				continue;
			}
			if(s.startsWith("#insidecount\t")){insideCount=Integer.parseInt(s.substring(13).trim()); continue;}
			if(s.startsWith("#outsidecount\t")){outsideCount=Integer.parseInt(s.substring(14).trim()); continue;}
			if(s.charAt(0)=='#'){continue;}
			final int tab=s.indexOf('\t');
			assert(tab>0) : "Malformed table line (no tab) in "+path+": "+s;
			assert(insideCount>=0 && outsideCount>=0) : "Table data line seen before its #insidecount/#outsidecount header in "+path+": "+s;
			final int k=insideCount+outsideCount;
			final String kmer=s.substring(0, tab);
			final int count=Integer.parseInt(s.substring(tab+1));
			assert(kmer.length()==k) : "Kmer length "+kmer.length()+" != declared insideCount+outsideCount="+k+" in "+path+": "+kmer;
			counts.set(packKmerString(kmer), count);
		}
		bf.close();
		assert(type!=null) : "No #boundary header found in "+path+" -- not a table written by this class's writeTable.";
		assert(insideCount>=0 && outsideCount>=0) : "No #insidecount/#outsidecount header found in "+path;
		final int K=insideCount+outsideCount;
		final TrnaBoundaryFeatures.NinemerTable nt=(seq, pos)->{
			final int kmer=packKmer(seq, pos, K);
			if(kmer<0){return 0;}//out-of-bounds/ambiguous window -> never-observed, count=0
			final int c=counts.get(kmer);
			return (c<0 ? 0 : c);//IntHashMap.get returns -1 for absent -> count=0
		};
		return new LoadedTable(K, type, insideCount, outsideCount, nt);
	}

	/** Packs a kmer STRING (already-decoded ACGT text, as written by unpackKmer/writeTable)
	 * back into the same 2-bits/base int packKmer produces from raw bytes. */
	private static int packKmerString(String kmer){
		int out=0;
		for(int i=0; i<kmer.length(); i++){
			final byte b=AminoAcid.baseToNumber[(byte)kmer.charAt(i)];
			assert(b>=0) : "Non-ACGT character in table kmer: "+kmer;
			out=(out<<2)|b;
		}
		return out;
	}

	private static void writeTable(String path, IntHashMap counts, String label, int insideCount, int outsideCount){
		final int k=insideCount+outsideCount;
		final int[] keys=counts.toKeyArray();
		java.util.Arrays.sort(keys);
		try(PrintStream out=new PrintStream(new java.io.FileOutputStream(path))){
			out.println("#boundary\t"+label);
			out.println("#insidecount\t"+insideCount);
			out.println("#outsidecount\t"+outsideCount);
			out.println("#k\t"+k);
			out.println("#kmer\tcount");
			for(int key : keys){
				out.println(unpackKmer(key, k)+"\t"+counts.get(key));
			}
		}catch(java.io.IOException e){
			throw new RuntimeException("Failed writing "+path, e);
		}
	}

	private static final byte[] NUM_TO_BASE={'A','C','G','T'};
	private static String unpackKmer(int key, int k){
		byte[] b=new byte[k];
		for(int i=k-1; i>=0; i--){b[i]=NUM_TO_BASE[key&3]; key>>=2;}
		return new String(b);
	}

	private static byte[] revcomp(byte[] seq){
		byte[] out=new byte[seq.length];
		for(int i=0; i<seq.length; i++){out[i]=AminoAcid.baseToComplementExtended[seq[seq.length-1-i]];}
		return out;
	}

	/** Infers a .gff/.gff.gz path from a .fna path by AnnotateAnticodon.inferGff's rule:
	 * same directory, same basename, extension swapped. Returns null if neither exists. */
	private static String inferGff(String fna){
		String prefix=ReadWrite.stripExtension(fna);
		String gff=prefix+".gff";
		if(new File(gff).exists()){return gff;}
		String gz=gff+".gz";
		if(new File(gz).exists()){return gz;}
		return null;
	}

	private static List<String> readList(String path){
		List<String> out=new ArrayList<>();
		ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(path, FileFormat.TEXT, null, false, true));
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			out.add(new String(line));
		}
		bf.close();
		return out;
	}

	private static final class Locus {
		String contig; int start; int stop; char strand;//start/stop are 0-based, inclusive
	}

	/** Reads tRNA loci from a GFF via byte-level LineParser1 (no per-line regex/String[] split). */
	private static List<Locus> readGffLoci(String path){
		List<Locus> out=new ArrayList<>();
		ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(path, FileFormat.GFF, null, true, true));
		LineParser1 lp=new LineParser1((byte)'\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			lp.set(line);
			if(lp.terms()<7 || !lp.termEquals("tRNA", 2)){continue;}
			Locus loc=new Locus();
			loc.contig=new String(lp.parseByteArray(0));
			loc.start=lp.parseInt(3)-1;//GFF is 1-based inclusive
			loc.stop=lp.parseInt(4)-1;
			loc.strand=(char)lp.parseByte(6, 0);
			out.add(loc);
		}
		bf.close();
		return out;
	}

	/** Loads a genome FASTA into a seqid->bases map, registering both the full header key and
	 * its _tid_&lt;taxid&gt;-stripped form (mirrors AnnotateAnticodon.loadFasta exactly, so
	 * seqid resolution matches the already-proven Stage-1 pipeline on this same corpus). */
	private static HashMap<String,byte[]> loadFasta(String fna){
		ArrayList<Read> reads=ReadInputStream.toReads(fna, FileFormat.FA, -1);
		HashMap<String,byte[]> map=new HashMap<>(reads.size()*4);
		for(Read r : reads){
			map.put(r.id, r.bases);
			String stripped=stripTid(r.id);
			if(!stripped.equals(r.id) && !map.containsKey(stripped)){map.put(stripped, r.bases);}
		}
		return map;
	}

	private static String stripTid(String s){
		int idx=s.lastIndexOf("_tid_");
		if(idx<0){return s;}
		final int digitStart=idx+5;
		if(digitStart>=s.length()){return s;}
		for(int i=digitStart; i<s.length(); i++){
			if(!Character.isDigit(s.charAt(i))){return s;}
		}
		return s.substring(0, idx);
	}
}
