package prot;

import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import structures.ByteBuilder;
import structures.IntList;

/**
 * Offline builder for the production family-shortlist sidecar consumed by ProteinSearcher's
 * F4 shortlist scoring (Brian, 2026-09-03, relayed by UMP45: "integrate the hybrid shortlist
 * into the assigner" -- the assigner is the current bottleneck at 10 genes/s per 32 threads).
 * Bundles, per family, exactly what a query needs scored against WITHOUT re-scanning the family
 * corpus at search time:
 * <ul>
 * <li>a unit-normalized aa20 2-mer composition vector ({@code float[400]}, the same recipe as
 *     {@link CompositionProfileAssay}'s v3/HybridShortlistAssay's config -- see {@link #variant}
 *     and {@link #rawPseudocount} for the exact parameters, confirmed against UMP45 before
 *     building);</li>
 * <li>the family's covering-set k-mers, re-emitted from an existing already-gated covering-set
 *     file (this tool does not recompute covering sets, only bundles them);</li>
 * <li>length stats: median member length, sd(ln member length).</li>
 * </ul>
 * <p>The header also carries the build-time background frequency vector ({@code float[400]}) --
 * a live query's own composition vector must be "centered" (freq-background) against this EXACT
 * background to be comparable to the family vectors, so the sidecar is self-contained: scoring a
 * query needs nothing from {@code familiesdir=}.</p>
 *
 * <p><b>No held-out exclusion here, unlike the evaluation assays.</b> A production query is
 * never a member of the training corpus that built these profiles, so there is nothing to
 * exclude -- the full in-sample family composition/length statistics are used directly. This is
 * simpler than {@link CompositionProfileAssay}/{@link HybridShortlistAssay}'s held-out machinery,
 * not a missing feature.</p>
 *
 * <p>The sidecar's header records the live {@code consensus=} file's own SHA-256 and the
 * covering-set file's SHA-256; a loader (ProteinSearcher's production code, not this tool) is
 * expected to re-hash both live inputs at startup and refuse to load a sidecar built from a
 * different consensus/covering-set version -- a stale sidecar must never silently score against
 * a since-rebuilt consensus.</p>
 *
 * <p>Family order/identity is taken from {@code consensus=}'s own {@code rank=N} header field
 * (not assumed to equal file order or any other family list) -- cross-checked against the
 * covering-set file's family set (exact match required, both directions) and the family corpus
 * directory's {@code {rank}.faa} files (required to exist for every rank 0..nFamilies-1, via
 * {@link CompositionProfileAssay#scanFamilyCorpus}, which already crashes loud on a missing
 * file).</p>
 *
 * <p>Usage: {@code java -ea prot.FamilyShortlistSidecarBuilder consensus=<consensus_reps_v6.fasta>
 * familiesdir=<family corpus dir, {rank}.faa per family> kmersets=<covering_set/sets.tsv>
 * out=<sidecar.tsv> [variant=centered] [rawpseudo=0] [background=families]}<br>
 * Self-test: {@code java -ea prot.FamilyShortlistSidecarBuilder selftest}</p>
 *
 * <p><b>Sidecar v3 (2026-09-04, plans/PER_FAMILY_THRESHOLDS_v1.md WP-A A2, entirely optional --
 * everything above is UNCHANGED and byte-identical when {@code thresholds=} is omitted):</b> add
 * {@code thresholds=<results/family_thresholds_v1.tsv> aaalignerfile=<AAAligner.java>
 * blosum62file=<Blosum62.java> alignercommit=<pinned BBTools commit hash> percentilep=<chosen p>
 * nmin=<calibration-set-size fallback cutoff>
 * [globalminid=0] [globalminratio=0] [globalminscore=0] [globalmincovq=0] [globalmincovt=0]
 * [globalmindimer=0] [globalminkmers=0] [globalminlenratio=0] [globalmaxlenratio=Infinity]} to
 * emit the per-family APPLIED thresholds (max/min of each family's own percentile and the global
 * floor, sec 2: "the stricter wins") as nine extra columns, plus the aligner binding (name +
 * sha256 of {@link AAAligner}/{@link Blosum62} at the pinned commit) and every global floor value
 * used, in the header. {@link FamilyShortlistSidecar#load} refuses to load a v3 sidecar whose
 * recorded aligner does not match the one about to run it.</p>
 *
 * @author Ady
 */
public final class FamilyShortlistSidecarBuilder {

	public static void main(final String[] args) throws Exception{
		if(args.length==1 && args[0].equalsIgnoreCase("selftest")){selftest(); return;}
		new FamilyShortlistSidecarBuilder(args).process();
	}

	FamilyShortlistSidecarBuilder(final String[] args){
		for(String arg : args){
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase();
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("consensus")){consensusFile=b;}
			else if(a.equals("familiesdir") || a.equals("families")){familiesDir=b;}
			else if(a.equals("kmersets")){kmerSetsFile=b;}
			else if(a.equals("out")){outFile=b;}
			else if(a.equals("variant")){variant=b.toLowerCase();}
			else if(a.equals("rawpseudo")){rawPseudocount=Double.parseDouble(b);}
			else if(a.equals("background")){backgroundMode=b.toLowerCase();}
			else if(a.equals("thresholds")){thresholdsFile=b;}
			else if(a.equals("aaalignerfile")){aaalignerFile=b;}
			else if(a.equals("blosum62file")){blosum62File=b;}
			else if(a.equals("alignercommit")){alignerCommit=b;}
			else if(a.equals("percentilep")){percentileP=Double.parseDouble(b);}
			else if(a.equals("nmin")){nMin=Integer.parseInt(b);}
			else if(a.equals("globalminid")){globalMinId=Double.parseDouble(b);}
			else if(a.equals("globalminratio")){globalMinRatio=Double.parseDouble(b);}
			else if(a.equals("globalminscore")){globalMinScore=Double.parseDouble(b);}
			else if(a.equals("globalmincovq")){globalMinCovQ=Double.parseDouble(b);}
			else if(a.equals("globalmincovt")){globalMinCovT=Double.parseDouble(b);}
			else if(a.equals("globalmindimer")){globalMinDimer=Double.parseDouble(b);}
			else if(a.equals("globalminkmers")){globalMinKmers=Double.parseDouble(b);}
			else if(a.equals("globalminlenratio")){globalMinLenRatio=Double.parseDouble(b);}
			else if(a.equals("globalmaxlenratio")){globalMaxLenRatio=Double.parseDouble(b);}
			else{throw new RuntimeException("Unknown argument: "+arg);}
		}
		if(consensusFile==null || familiesDir==null || kmerSetsFile==null || outFile==null){
			throw new RuntimeException("consensus=, familiesdir=, kmersets=, and out= are required.");
		}
		if(!variant.equals("raw") && !variant.equals("centered") && !variant.equals("logodds")){
			throw new RuntimeException("variant must be raw/centered/logodds: "+variant);
		}
		if(rawPseudocount<0){throw new RuntimeException("rawpseudo must be >=0: "+rawPseudocount);}
		if(!backgroundMode.equals("members") && !backgroundMode.equals("families")){
			throw new RuntimeException("background must be 'members' or 'families': "+backgroundMode);
		}
		//thresholds= (plans/PER_FAMILY_THRESHOLDS_v1.md WP-A A2) is entirely optional -- when unset,
		//every new field below is unused and this class's behavior/output is BYTE-IDENTICAL to
		//before A2 (the sidecar v3 columns/header fields only appear when thresholds= is given).
		if(thresholdsFile!=null){
			if(aaalignerFile==null || blosum62File==null){
				throw new RuntimeException("thresholds= requires aaalignerfile= and blosum62file= "+
					"(sidecar v3 records their sha256 so a stale aligner build can never silently score "+
					"against a table calibrated for a different one).");
			}
			if(alignerCommit==null){
				throw new RuntimeException("thresholds= requires alignercommit= (sec 0: \"at the pinned BBTools "+
					"commit\" -- the source sha256 pins CONTENT, the commit pins PROVENANCE; both are recorded).");
			}
			if(percentileP<0){throw new RuntimeException("thresholds= requires percentilep= (the chosen percentile, >=0).");}
			if(nMin<1){throw new RuntimeException("thresholds= requires nmin= (>=1, the calibration-set-size fallback cutoff).");}
		}
	}

	private void process(){
		System.err.println("Hashing "+consensusFile+"...");
		final String consensusSha256=DigestSuffix.file(consensusFile);
		System.err.println("consensus sha80="+consensusSha256);

		System.err.println("Loading family order/identity from "+consensusFile+"...");
		final String[] repIds=loadConsensusFamilyOrder(consensusFile);
		final int nFamilies=repIds.length;
		System.err.println("Loaded "+nFamilies+" families.");

		System.err.println("Hashing "+kmerSetsFile+"...");
		final String kmerSetsSha256=DigestSuffix.file(kmerSetsFile);
		System.err.println("Loading covering sets from "+kmerSetsFile+"...");
		final ReducedAlphabetSeedAssay.KmerSets sets=ReducedAlphabetSeedAssay.loadKmerSets(kmerSetsFile);
		//Exact family-set match, both directions -- a covering-set file built against a different
		//family list (stale run, wrong kd/alphabet/k) must never be silently bundled.
		final HashSet<String> consensusIds=new HashSet<String>(Arrays.asList(repIds));
		final ArrayList<String> missingFromCovering=new ArrayList<String>();
		for(String rep : repIds){if(!sets.families.containsKey(rep)){missingFromCovering.add(rep);}}
		if(!missingFromCovering.isEmpty()){
			throw new RuntimeException(missingFromCovering.size()+" consensus families have no covering-set "+
				"entry in "+kmerSetsFile+" (e.g. '"+missingFromCovering.get(0)+"') -- covering-set file does not "+
				"match this consensus.");
		}
		final ArrayList<String> extraInCovering=new ArrayList<String>();
		for(String rep : sets.families.keySet()){if(!consensusIds.contains(rep)){extraInCovering.add(rep);}}
		if(!extraInCovering.isEmpty()){
			throw new RuntimeException(extraInCovering.size()+" covering-set families are not in the consensus "+
				"(e.g. '"+extraInCovering.get(0)+"') -- covering-set file does not match this consensus.");
		}
		System.err.println("Covering-set family set matches consensus exactly ("+nFamilies+" families, "+
			sets.words+" total k-mers, alphabet="+sets.alphabet.name+" k="+sets.k+").");

		//---- Composition: aa20 2-mer, full in-sample (no held-out -- see class javadoc). ----
		final CompositionProfileAssay.AlphaK config=new CompositionProfileAssay.AlphaK("aa20", ReducedAlphabet.named("aa20"), 2);
		System.err.println("Scanning family corpus for composition profiles ("+config.label+")...");
		final long scanStart=System.nanoTime();
		final double[][] famCounts=new double[nFamilies][config.dims];
		final int[] famMemberCount=new int[nFamilies];
		CompositionProfileAssay.scanFamilyCorpus(familiesDir, nFamilies, config, famCounts, famMemberCount);
		System.err.println("Scanned in "+((System.nanoTime()-scanStart)/1e9)+"s.");

		final double[] bgFreq=CompositionProfileAssay.computeBackgroundFreq(famCounts, nFamilies, config.dims,
			rawPseudocount, backgroundMode);
		//A live query at scoring time must apply the SAME "centered" transform (freq-bgFreq) using
		//THIS exact background -- recomputing a fresh background from a different/smaller corpus at
		//query time would silently produce inconsistent scores. Persisted below in the header so the
		//sidecar is self-contained (a query needs nothing from familiesdir= to be scored).
		final float[] bgFreqOut=new float[config.dims];
		for(int d=0; d<config.dims; d++){bgFreqOut[d]=(float)bgFreq[d];}
		final float[][] famVec=new float[nFamilies][config.dims];
		for(int r=0; r<nFamilies; r++){
			final double[] freq=CompositionProfileAssay.toFreq(famCounts[r], config.dims, rawPseudocount);
			final double[] dst=new double[config.dims];
			CompositionProfileAssay.applyVariant(freq, bgFreq, variant, dst);
			CompositionProfileAssay.normalize(dst);
			for(int d=0; d<config.dims; d++){famVec[r][d]=(float)dst[d];}
		}

		//---- Length stats: median, sd(ln len) -- full in-sample, no held-out. ----
		System.err.println("Scanning family corpus for member lengths...");
		final double[] medLen=new double[nFamilies], sdLnLen=new double[nFamilies];
		scanFamilyLengths(familiesDir, nFamilies, medLen, sdLnLen);

		final AppliedThresholds thresholds=(thresholdsFile==null ? null :
			loadAndApplyThresholds(thresholdsFile, repIds));

		writeSidecar(repIds, sets, famVec, bgFreqOut, medLen, sdLnLen, config, consensusSha256, kmerSetsSha256,
			thresholds);
	}

	/** One family's calibration-set size, fallback flag, and the nine APPLIED (max/min of the
	 *  global floor and the family's own percentile -- sec 2: "the stricter wins") threshold
	 *  values, in family-rank order (matching {@code repIds}). A fallback family (fewer than
	 *  {@link #nMin} calibration members, sec 2: "percentiles over &lt;20 points are noise, not
	 *  statistics") gets the global floor directly rather than a max/min combination -- "gets the
	 *  global floors only", not "gets its own noisy percentile combined with the floor". */
	private static final class AppliedThresholds{
		final int[] n; final boolean[] fallback;
		final double[] rLo, rHi, minDimer, minKmers, minId, minRatio, minScore, minCovQ, minCovT;
		AppliedThresholds(final int nFamilies){
			n=new int[nFamilies]; fallback=new boolean[nFamilies];
			rLo=new double[nFamilies]; rHi=new double[nFamilies]; minDimer=new double[nFamilies];
			minKmers=new double[nFamilies]; minId=new double[nFamilies]; minRatio=new double[nFamilies];
			minScore=new double[nFamilies]; minCovQ=new double[nFamilies]; minCovT=new double[nFamilies];
		}
	}

	/** Loads {@code results/family_thresholds_v1.tsv} (WP-B's product; sec 3 A2's exact column
	 *  list), validates it matches {@code repIds} exactly (same family-set-match discipline as
	 *  the covering-set check above -- a table built for a different family list must never be
	 *  silently applied), refuses NaN/negative/non-finite values, and combines each family's own
	 *  percentile with the global floor per sec 2 ("applied = max(global floor, family value) for
	 *  lower cutoffs, min(global ceiling, family value) for r_hi"). */
	private AppliedThresholds loadAndApplyThresholds(final String path, final String[] repIds){
		final int nFamilies=repIds.length;
		final HashMap<String,Integer> rankOf=new HashMap<String,Integer>();
		for(int r=0; r<nFamilies; r++){rankOf.put(repIds[r], r);}

		final ByteFile bf=ByteFile.makeByteFile(path, false);
		final String EXPECTED_HEADER="family\tn\tfallback\tr_lo\tr_hi\tmindimer\tminkmers\tminid\tminratio\tminscore\tmincovq\tmincovt";
		String headerLine=null;
		final AppliedThresholds t=new AppliedThresholds(nFamilies);
		final boolean[] seen=new boolean[nFamilies];
		int dataRows=0;
		for(byte[] lineBytes=bf.nextLine(); lineBytes!=null; lineBytes=bf.nextLine()){
			if(lineBytes.length==0){continue;}
			final String line=new String(lineBytes);
			if(headerLine==null){
				headerLine=line;
				if(!headerLine.equals(EXPECTED_HEADER)){
					throw new RuntimeException("Threshold table header mismatch in "+path+
						":\n  expected: "+EXPECTED_HEADER+"\n  got:      "+headerLine);
				}
				continue;
			}
			final String[] c=line.split("\t", -1);
			if(c.length!=12){
				throw new RuntimeException("Threshold table "+path+" row has "+c.length+
					" columns, expected 12: "+line);
			}
			final String family=c[0];
			final Integer rank=rankOf.get(family);
			if(rank==null){
				throw new RuntimeException("Threshold table "+path+" has family '"+family+
					"' which is not in the consensus/covering-set family list.");
			}
			if(seen[rank]){throw new RuntimeException("Threshold table "+path+" has duplicate family '"+family+"'.");}
			seen[rank]=true;
			final int n=parseNonNegInt(c[1], "n", family, path);
			final boolean fallback=parseBoolFlag(c[2], "fallback", family, path);
			final double rLoFam=parseFinite(c[3], "r_lo", family, path);
			final double rHiFam=parseFinite(c[4], "r_hi", family, path);
			final double mdFam=parseFinite(c[5], "mindimer", family, path);
			final double mkFam=parseFinite(c[6], "minkmers", family, path);
			final double miFam=parseFinite(c[7], "minid", family, path);
			final double mrFam=parseFinite(c[8], "minratio", family, path);
			final double msFam=parseFinite(c[9], "minscore", family, path);
			final double mcqFam=parseFinite(c[10], "mincovq", family, path);
			final double mctFam=parseFinite(c[11], "mincovt", family, path);

			t.n[rank]=n; t.fallback[rank]=fallback;
			//sec 2: a fallback family "gets the global floors only" -- not a max/min combination
			//with its own (noisy, <n_min-sample) percentile.
			t.rLo[rank]=(fallback ? globalMinLenRatio : Math.max(globalMinLenRatio, rLoFam));
			t.rHi[rank]=(fallback ? globalMaxLenRatio : Math.min(globalMaxLenRatio, rHiFam));
			t.minDimer[rank]=(fallback ? globalMinDimer : Math.max(globalMinDimer, mdFam));
			t.minKmers[rank]=(fallback ? globalMinKmers : Math.max(globalMinKmers, mkFam));
			t.minId[rank]=(fallback ? globalMinId : Math.max(globalMinId, miFam));
			t.minRatio[rank]=(fallback ? globalMinRatio : Math.max(globalMinRatio, mrFam));
			t.minScore[rank]=(fallback ? globalMinScore : Math.max(globalMinScore, msFam));
			t.minCovQ[rank]=(fallback ? globalMinCovQ : Math.max(globalMinCovQ, mcqFam));
			t.minCovT[rank]=(fallback ? globalMinCovT : Math.max(globalMinCovT, mctFam));
			dataRows++;
		}
		bf.close();
		if(headerLine==null){throw new RuntimeException("Threshold table "+path+" is empty (no header line).");}
		if(dataRows!=nFamilies){
			throw new RuntimeException("Threshold table "+path+" has "+dataRows+" families, expected "+
				nFamilies+" (must cover every consensus family exactly once).");
		}
		return t;
	}

	private static int parseNonNegInt(final String s, final String field, final String family, final String path){
		final int v;
		try{v=Integer.parseInt(s);}catch(NumberFormatException e){
			throw new RuntimeException("Threshold table "+path+", family '"+family+"': non-integer "+field+"='"+s+"'.");
		}
		if(v<0){throw new RuntimeException("Threshold table "+path+", family '"+family+"': negative "+field+"="+v+".");}
		return v;
	}
	private static boolean parseBoolFlag(final String s, final String field, final String family, final String path){
		if(s.equals("0")){return false;}
		if(s.equals("1")){return true;}
		throw new RuntimeException("Threshold table "+path+", family '"+family+"': "+field+" must be 0 or 1, got '"+s+"'.");
	}
	/** Refuses NaN, negative, and non-finite (Infinity) values -- sec 3 A2's loader contract
	 *  ("refuses NaN/negative"). r_hi has no natural upper bound in the table itself (it is
	 *  combined against {@link #globalMaxLenRatio}, which MAY be infinite by default), so this
	 *  intentionally rejects a literal Infinity/NaN token in the FAMILY value while still allowing
	 *  the separately-supplied global ceiling to be infinite. */
	private static double parseFinite(final String s, final String field, final String family, final String path){
		final double v;
		try{v=Double.parseDouble(s);}catch(NumberFormatException e){
			throw new RuntimeException("Threshold table "+path+", family '"+family+"': non-numeric "+field+"='"+s+"'.");
		}
		if(Double.isNaN(v) || Double.isInfinite(v) || v<0){
			throw new RuntimeException("Threshold table "+path+", family '"+family+"': "+field+"="+v+
				" must be finite and non-negative.");
		}
		return v;
	}

	/** Loads family rank order/identity directly from the consensus FASTA's own header
	 *  fields ({@code >repId rank=N ...} or {@code >repId active_index=N ...}), cross-checking
	 *  file order matches the declared dense index
	 *  (self-consistency: a reordered or corrupted file must not silently produce a wrong
	 *  rank-to-repid mapping). */
	static String[] loadConsensusFamilyOrder(final String consensusFile){
		final ArrayList<String> reps=new ArrayList<String>();
		final ByteFile bf=ByteFile.makeByteFile(consensusFile, false);
		int expectedRank=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]!='>'){continue;}
			final String header=new String(line, 1, line.length-1);
			final String[] tokens=header.split("\\s+");
			final String repId=tokens[0];
			int rank=-1,activeIndex=-1;
			for(int i=1; i<tokens.length; i++){
				if(tokens[i].startsWith("rank=")){
					if(rank>=0){throw new RuntimeException("Consensus header has duplicate rank= fields: '"+header+"'");}
					rank=Integer.parseInt(tokens[i].substring(5));
				}else if(tokens[i].startsWith("active_index=")){
					if(activeIndex>=0){throw new RuntimeException("Consensus header has duplicate active_index= fields: '"+header+"'");}
					activeIndex=Integer.parseInt(tokens[i].substring(13));
				}
			}
			if(rank<0 && activeIndex<0){throw new RuntimeException("Consensus header has neither rank= nor active_index= field: '"+header+"'");}
			if(rank>=0 && activeIndex>=0 && rank!=activeIndex){throw new RuntimeException("Consensus header rank=/active_index= disagreement: '"+header+"'");}
			final int denseIndex=activeIndex>=0 ? activeIndex : rank;
			if(denseIndex!=expectedRank){
				throw new RuntimeException("Consensus family order mismatch: file position "+expectedRank+
					" has dense index="+denseIndex+" (rep '"+repId+"') -- expected dense indices in ascending file order starting at 0.");
			}
			reps.add(repId);
			expectedRank++;
		}
		bf.close();
		if(reps.isEmpty()){throw new RuntimeException("No records found in "+consensusFile);}
		return reps.toArray(new String[0]);
	}

	/** Per-family median member length and sd(ln member length) -- full in-sample (no held-out;
	 *  see class javadoc). Lighter than {@link HybridShortlistAssay#scanFamilyLengths} since no
	 *  per-query exclusion machinery is needed here. */
	static void scanFamilyLengths(final String familiesDir, final int nFamilies,
			final double[] medLenOut, final double[] sdLnLenOut){
		//TODO: Probable bug - raw lengths include a conventional terminal '*', while
		//ProteinSequence and ProteinSearcher F4/triage use encoded lengths without it.
		//Correct with no-star/starred/multiline parity tests; preserve accepted sidecars.
		for(int rank=0; rank<nFamilies; rank++){
			final String file=familiesDir+"/"+rank+".faa";
			final ByteFile bf=ByteFile.makeByteFile(file, false);
			final IntList lens=new IntList(256);
			int curLen=0; boolean inRecord=false;
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0){continue;}
				if(line[0]=='>'){
					if(inRecord){lens.add(curLen);}
					curLen=0; inRecord=true;
				}else{
					curLen+=line.length;
				}
			}
			if(inRecord){lens.add(curLen);}
			bf.close();
			if(lens.size<1){throw new RuntimeException("Family rank "+rank+" ("+file+") has no members for length scan.");}
			final int[] arr=lens.toArray();
			Arrays.sort(arr);
			final int n=arr.length;
			medLenOut[rank]=(n%2==1 ? arr[n/2] : (arr[n/2-1]+arr[n/2])/2.0);
			double sumLn=0, sumLnSq=0;
			for(int len : arr){final double ln=Math.log(Math.max(1, len)); sumLn+=ln; sumLnSq+=ln*ln;}
			final double meanLn=sumLn/n;
			sdLnLenOut[rank]=Math.sqrt(Math.max(0, sumLnSq/n-meanLn*meanLn));
		}
	}

	private void writeSidecar(final String[] repIds, final ReducedAlphabetSeedAssay.KmerSets sets,
			final float[][] famVec, final float[] bgFreq, final double[] medLen, final double[] sdLnLen,
			final CompositionProfileAssay.AlphaK config, final String consensusSha256, final String kmerSetsSha256,
			final AppliedThresholds thresholds){
		final int nFamilies=repIds.length;
		final FileFormat ff=FileFormat.testOutput(outFile, FileFormat.TXT, null, true, true, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff); bsw.start();
		try{
			final ByteBuilder bb=new ByteBuilder(1<<20);
			//schema_version 2 (thresholds present) is a pure ADDITION -- every schema_version 1
			//field/column stays in the exact same position, so a schema-1-only reader (unaware of
			//A2) still parses everything it understands correctly; it just doesn't see the tail.
			bb.append("#schema_version\t").append(thresholds==null ? 1 : 2).append('\n');
			bb.append("#kind\tfamily_shortlist_sidecar_v1\n");
			bb.append("#consensus_file\t").append(consensusFile).append('\n');
			bb.append("#consensus_sha80\t").append(consensusSha256).append('\n');
			bb.append("#kmersets_file\t").append(kmerSetsFile).append('\n');
			bb.append("#kmersets_sha80\t").append(kmerSetsSha256).append('\n');
			bb.append("#composition_alphabet\t").append(config.alphabetName).append('\n');
			bb.append("#composition_k\t").append(config.k).append('\n');
			bb.append("#composition_dims\t").append(config.dims).append('\n');
			bb.append("#composition_variant\t").append(variant).append('\n');
			bb.append("#composition_rawpseudo\t").append(rawPseudocount, 8).append('\n');
			bb.append("#composition_background\t").append(backgroundMode).append('\n');
			bb.append("#covering_alphabet\t").append(sets.alphabet.name).append('\n');
			bb.append("#covering_k\t").append(sets.k).append('\n');
			bb.append("#n_families\t").append(nFamilies).append('\n');
			if(thresholds!=null){
				//A2 (plans/PER_FAMILY_THRESHOLDS_v1.md sec 3): aligner binding + the global floors
				//used at build time, recorded as VALUES (not just hashes) so a reader can audit them
				//without re-deriving anything.
				bb.append("#aligner\tAAAligner\n");
				bb.append("#aligner_sha80\t").append(DigestSuffix.file(aaalignerFile)).append('\n');
				bb.append("#aligner_commit\t").append(alignerCommit).append('\n');
				bb.append("#blosum62_sha80\t").append(DigestSuffix.file(blosum62File)).append('\n');
				//No #thresholds_file (path) header (UMP45's review, 2026-09-04): a machine-local build
				//directory path would make the sidecar's own artifact bytes -- and therefore its
				//sha256, which downstream tools bind to -- depend on WHERE it was built, not just what
				//it was built from. #thresholds_sha256 alone is the binding; that's host-independent.
				bb.append("#thresholds_sha80\t").append(DigestSuffix.file(thresholdsFile)).append('\n');
				bb.append("#percentile_p\t").append(percentileP, 4).append('\n');
				bb.append("#n_min\t").append(nMin).append('\n');
				bb.append("#global_minid\t").append(globalMinId, 6).append('\n');
				bb.append("#global_minratio\t").append(globalMinRatio, 6).append('\n');
				bb.append("#global_minscore\t").append(globalMinScore, 6).append('\n');
				bb.append("#global_mincovq\t").append(globalMinCovQ, 6).append('\n');
				bb.append("#global_mincovt\t").append(globalMinCovT, 6).append('\n');
				bb.append("#global_mindimer\t").append(globalMinDimer, 6).append('\n');
				bb.append("#global_minkmers\t").append(globalMinKmers, 6).append('\n');
				bb.append("#global_minlenratio\t").append(globalMinLenRatio, 6).append('\n');
				bb.append("#global_maxlenratio\t"); appendCeiling(bb, globalMaxLenRatio); bb.nl();
			}
			//The background vector used at build time for every family's "centered" transform -- a
			//live query MUST apply this same background at scoring time (see process()'s comment).
			bb.append("#background_"+config.dims+"f\t");
			for(int d=0; d<config.dims; d++){
				if(d>0){bb.append(',');}
				bb.append(bgFreq[d], 8);
			}
			bb.nl();
			bb.append("rank\trep_id\tmedian_len\tsd_ln_len\tcomposition_"+config.dims+"f\tcovering_kmers_hex");
			if(thresholds!=null){
				bb.append("\tn\tfallback\tr_lo\tr_hi\tmindimer\tminkmers\tminid\tminratio\tminscore\tmincovq\tmincovt");
			}
			bb.nl();
			bsw.print(bb); bb.clear();
			for(int r=0; r<nFamilies; r++){
				bb.append(r).append('\t').append(repIds[r]).append('\t').append(medLen[r], 4).append('\t')
					.append(sdLnLen[r], 6).append('\t');
				for(int d=0; d<config.dims; d++){
					if(d>0){bb.append(',');}
					bb.append(famVec[r][d], 8);
				}
				bb.append('\t');
				final long[] kmers=sets.families.get(repIds[r]).toArray();
				Arrays.sort(kmers);
				for(int i=0; i<kmers.length; i++){
					if(i>0){bb.append(',');}
					bb.append(hex16(kmers[i]));
				}
				if(thresholds!=null){
					bb.append('\t').append(thresholds.n[r]).append('\t').append(thresholds.fallback[r] ? 1 : 0)
						.append('\t').append(thresholds.rLo[r], 6).append('\t'); appendCeiling(bb, thresholds.rHi[r]);
					bb.append('\t').append(thresholds.minDimer[r], 6).append('\t').append(thresholds.minKmers[r], 6)
						.append('\t').append(thresholds.minId[r], 6).append('\t').append(thresholds.minRatio[r], 6)
						.append('\t').append(thresholds.minScore[r], 6).append('\t').append(thresholds.minCovQ[r], 6)
						.append('\t').append(thresholds.minCovT[r], 6);
				}
				bb.nl();
				bsw.print(bb); bb.clear();
			}
		}finally{
			if(bsw.poisonAndWait()){throw new RuntimeException("I/O error writing "+outFile);}
		}
		final String outSha256=DigestSuffix.file(outFile);
		final FileFormat ffSha=FileFormat.testOutput(outFile+".sha80", FileFormat.TXT, null, false, true, false, false);
		final ByteStreamWriter w=new ByteStreamWriter(ffSha); w.start();
		w.print(new ByteBuilder().append(outSha256).append("  ").append(outFile).nl());
		if(w.poisonAndWait()){throw new RuntimeException("I/O error writing "+outFile+".sha80");}
		System.err.println("Wrote sidecar for "+nFamilies+" families to "+outFile+" (sha80="+outSha256+").");
	}

	/** Writes a length-ratio CEILING value, handling {@code Double.POSITIVE_INFINITY} explicitly
	 *  (sec 2's "off" default for the one ceiling metric) -- {@link ByteBuilder#append(double,int)}
	 *  is a fixed-point formatter and silently corrupts Infinity into a large finite garbage value
	 *  (empirically verified: renders as 9223372036854775807.775807, a {@code Long.MAX_VALUE}
	 *  overflow artifact, not a real number and not parseable back as Infinity) -- self-caught
	 *  2026-09-04 before it ever reached a real sidecar. */
	private static void appendCeiling(final ByteBuilder bb, final double v){
		//TODO: Probable bug in structures/ByteBuilder.append(double,int) -- it is a fixed-point
		//formatter and silently corrupts Double.POSITIVE_INFINITY into 9223372036854775807.775807
		//(a Long.MAX_VALUE overflow artifact, not parseable back as Infinity) instead of throwing or
		//emitting "Infinity". Confirmed empirically 2026-09-04 (Eru); relayed to Brian by UMP45.
		//Worked around here rather than fixed at the source -- ByteBuilder is shared BBTools
		//infrastructure, not mag-qc-owned.
		if(Double.isInfinite(v)){bb.append("Infinity");}else{bb.append(v, 6);}
	}

	static String hex16(final long value){
		final char[] c=new char[16];
		for(int i=15, s=0; i>=0; i--, s+=4){final int d=(int)((value>>>s)&15); c[i]=(char)(d<10 ? '0'+d : 'a'+d-10);}
		return new String(c);
	}

	// ---------------- self-test ----------------
	static void selftest() throws Exception{
		final String dir=System.getProperty("java.io.tmpdir")+"/sidecar_selftest_"+System.nanoTime();
		new java.io.File(dir).mkdirs();
		final String famDir=dir+"/fam";
		new java.io.File(famDir).mkdirs();

		//3 tiny families, distinct compositions and lengths.
		writeFile(famDir+"/0.faa", ">a1\nAAAAAAAAAA\n>a2\nAAAAAAAAAAAA\n");//fam0: pure A, lens 10,12
		writeFile(famDir+"/1.faa", ">b1\nCCCCCCCCCC\n>b2\nCCCCCCCCCCCC\n>b3\nCCCCCCCCCCCCCC\n");//fam1: pure C, lens 10,12,14
		writeFile(famDir+"/2.faa", ">c1\nACACACACAC\n");//fam2: alternating, len 10

		writeFile(dir+"/consensus.fasta", ">rep0 rank=0 artifact=x\nAAAAAAAAAA\n"+
			">rep1 rank=1 artifact=x\nCCCCCCCCCC\n"+">rep2 rank=2 artifact=x\nACACACACAC\n");

		final ByteBuilder sets=new ByteBuilder();
		sets.append("#schema_version\t1\n#kind\tcovering_set\n#bbtools_commit\ttest\n");
		sets.append("#alphabet\taa20\n#k\t5\n");
		sets.append("rep0\tAAAAA\n");
		sets.append("rep1\tCCCCC\n");
		sets.append("rep2\tACACA\nrep2\tCACAC\n");
		writeFile(dir+"/sets.tsv", sets.toString());

		//1. Family order loading + rank-mismatch crash-loud.
		final String[] repIds=loadConsensusFamilyOrder(dir+"/consensus.fasta");
		if(!(repIds.length==3 && repIds[0].equals("rep0") && repIds[1].equals("rep1") && repIds[2].equals("rep2"))){
			throw new RuntimeException("loadConsensusFamilyOrder FAILED: "+Arrays.toString(repIds));
		}
		System.err.println("loadConsensusFamilyOrder: PASS");
		writeFile(dir+"/consensus_active.fasta", ">rep0 active_index=0 family_id=7\nAAAAAAAAAA\n"+
			">rep1 active_index=1 family_id=19\nCCCCCCCCCC\n");
		final String[] activeIds=loadConsensusFamilyOrder(dir+"/consensus_active.fasta");
		if(!(activeIds.length==2 && activeIds[0].equals("rep0") && activeIds[1].equals("rep1"))){
			throw new RuntimeException("active_index family-order load FAILED: "+Arrays.toString(activeIds));
		}
		writeFile(dir+"/consensus_disagree.fasta", ">rep0 rank=0 active_index=1\nAAAAAAAAAA\n");
		boolean threwDisagreement=false;
		try{loadConsensusFamilyOrder(dir+"/consensus_disagree.fasta");}catch(RuntimeException e){threwDisagreement=true;}
		if(!threwDisagreement){throw new RuntimeException("rank/active_index disagreement crash-loud FAILED: did not throw.");}

		writeFile(dir+"/consensus_badorder.fasta", ">rep0 rank=1 artifact=x\nAAAAAAAAAA\n");
		boolean threwOrder=false;
		try{loadConsensusFamilyOrder(dir+"/consensus_badorder.fasta");}catch(RuntimeException e){threwOrder=true;}
		if(!threwOrder){throw new RuntimeException("Rank-mismatch crash-loud FAILED: did not throw.");}
		System.err.println("Rank-mismatch crash-loud: PASS");

		//2. Length stats: fam0 median of {10,12}=11; fam1 median of {10,12,14}=12.
		final double[] medLen=new double[3], sdLnLen=new double[3];
		scanFamilyLengths(famDir, 3, medLen, sdLnLen);
		if(medLen[0]!=11.0 || medLen[1]!=12.0 || medLen[2]!=10.0){
			throw new RuntimeException("scanFamilyLengths FAILED: "+Arrays.toString(medLen));
		}
		System.err.println("scanFamilyLengths: PASS (medians "+Arrays.toString(medLen)+")");

		//3. Full run: covering-set/consensus mismatch crash-loud (extra family in sets.tsv).
		final ByteBuilder setsExtra=new ByteBuilder();
		setsExtra.append("#schema_version\t1\n#kind\tcovering_set\n#bbtools_commit\ttest\n#alphabet\taa20\n#k\t5\n");
		setsExtra.append("rep0\tAAAAA\nrep1\tCCCCC\nrep2\tACACA\nrepX\tGGGGG\n");
		writeFile(dir+"/sets_extra.tsv", setsExtra.toString());
		boolean threwExtra=false;
		try{
			new FamilyShortlistSidecarBuilder(new String[]{"consensus="+dir+"/consensus.fasta",
				"familiesdir="+famDir, "kmersets="+dir+"/sets_extra.tsv", "out="+dir+"/sidecar_bad.tsv"}).process();
		}catch(RuntimeException e){threwExtra=true;}
		if(!threwExtra){throw new RuntimeException("Covering-set-extra-family crash-loud FAILED: did not throw.");}
		System.err.println("Covering-set-extra-family crash-loud: PASS");

		//4. Full run, real (matching) inputs -- verify output structure.
		new FamilyShortlistSidecarBuilder(new String[]{"consensus="+dir+"/consensus.fasta",
			"familiesdir="+famDir, "kmersets="+dir+"/sets.tsv", "out="+dir+"/sidecar.tsv"}).process();
		final String content=readWhole(dir+"/sidecar.tsv");
		final String[] lines=content.split("\n");
		int dataLines=0;
		String fam0Row=null, fam1Row=null;
		for(String line : lines){
			if(line.startsWith("#") || line.startsWith("rank\t")){continue;}
			dataLines++;
			if(line.startsWith("0\trep0\t")){fam0Row=line;}
			if(line.startsWith("1\trep1\t")){fam1Row=line;}
		}
		if(dataLines!=3){throw new RuntimeException("Sidecar row count FAILED: expected 3, got "+dataLines);}
		if(fam0Row==null || fam1Row==null){throw new RuntimeException("Sidecar missing expected family rows.");}
		//fam0 (pure A) and fam1 (pure C) are orthogonal under aa20 2-mer composition -- verify the
		//covering-set hex field for fam0 contains exactly one k-mer (rep0 -> AAAAA, one 5-mer).
		final String[] fam0Cols=fam0Row.split("\t");
		if(fam0Cols.length!=6){throw new RuntimeException("Sidecar column count FAILED: "+fam0Cols.length);}
		final String[] fam0Kmers=fam0Cols[5].split(",");
		if(fam0Kmers.length!=1){throw new RuntimeException("fam0 covering-set k-mer count FAILED: "+fam0Cols[5]);}
		final String[] fam0Comp=fam0Cols[4].split(",");
		if(fam0Comp.length!=400){throw new RuntimeException("fam0 composition dim count FAILED: "+fam0Comp.length);}
		System.err.println("Full-run sidecar structure (3 rows, 400-dim composition, covering k-mers): PASS");

		//4b. Background vector header line: present, correctly-sized, and non-degenerate (fam0/fam1
		//are pure-A/pure-C -- their difference alone guarantees a non-uniform background).
		String bgLine=null;
		for(String line : lines){if(line.startsWith("#background_400f\t")){bgLine=line; break;}}
		if(bgLine==null){throw new RuntimeException("Sidecar missing #background_400f header line.");}
		final String[] bgVals=bgLine.substring(bgLine.indexOf('\t')+1).split(",");
		if(bgVals.length!=400){throw new RuntimeException("Background vector dim count FAILED: "+bgVals.length);}
		System.err.println("Background vector header (400-dim): PASS");

		//5. Sha80 sidecar exists and is well-formed.
		final String shaContent=readWhole(dir+"/sidecar.tsv.sha80");
		if(!shaContent.trim().matches("[0-9a-f]{20}\\s+.*")){
			throw new RuntimeException("sha80 sidecar malformed: "+shaContent);
		}
		System.err.println("sha80 sidecar: PASS");

		//---- 6. Sidecar v3 (A2): thresholds= round-trip, fallback handling, Infinity ceiling,
		//aligner-binding check, and loader crash-loud paths. Same 3-family fixture as above. ----
		writeFile(dir+"/fake_aaaligner.txt", "aligner source stand-in for hashing\n");
		writeFile(dir+"/fake_blosum62.txt", "blosum62 source stand-in for hashing\n");
		final String THRESH_HEADER="family\tn\tfallback\tr_lo\tr_hi\tmindimer\tminkmers\tminid\tminratio\tminscore\tmincovq\tmincovt";
		writeFile(dir+"/thresholds.tsv", THRESH_HEADER+"\n"+
			"rep0\t25\t0\t0.8\t1.2\t0.1\t2\t40\t0.5\t20\t0.6\t0.6\n"+//rank0: n>=nMin, family wins some cutoffs
			"rep1\t5\t1\t0.9\t1.1\t0.2\t3\t50\t0.6\t25\t0.65\t0.65\n"+//rank1: fallback -- family values IGNORED
			"rep2\t30\t0\t0.5\t2.0\t0.02\t0\t20\t0.1\t5\t0.5\t0.5\n");//rank2: n>=nMin, global floor wins every cutoff

		new FamilyShortlistSidecarBuilder(new String[]{"consensus="+dir+"/consensus.fasta",
			"familiesdir="+famDir, "kmersets="+dir+"/sets.tsv", "out="+dir+"/sidecar_v3.tsv",
			"thresholds="+dir+"/thresholds.tsv", "aaalignerfile="+dir+"/fake_aaaligner.txt",
			"blosum62file="+dir+"/fake_blosum62.txt", "alignercommit=54b62510", "percentilep=2", "nmin=20",
			"globalminid=30", "globalminratio=0.3", "globalminscore=15", "globalmincovq=0.7",
			"globalmincovt=0.7", "globalmindimer=0.05", "globalminkmers=1", "globalminlenratio=0.75"
			//globalmaxlenratio deliberately left at its Infinity default -- exercises the
			//Infinity-ceiling round-trip for rep1 (fallback) below.
			}).process();
		final String v3Content=readWhole(dir+"/sidecar_v3.tsv");
		if(!v3Content.contains("#schema_version\t2\n")){throw new RuntimeException("v3 sidecar missing schema_version=2.");}
		if(!v3Content.contains("#aligner\tAAAligner\n")){throw new RuntimeException("v3 sidecar missing #aligner header.");}
		if(!v3Content.contains("#aligner_commit\t54b62510\n")){throw new RuntimeException("v3 sidecar missing #aligner_commit.");}
		if(v3Content.contains("#thresholds_file\t")){
			throw new RuntimeException("v3 sidecar has #thresholds_file (host-dependent path) -- should be dropped.");
		}
		if(!v3Content.contains("rank\trep_id\tmedian_len\tsd_ln_len\tcomposition_400f\tcovering_kmers_hex\t"+
			"n\tfallback\tr_lo\tr_hi\tmindimer\tminkmers\tminid\tminratio\tminscore\tmincovq\tmincovt\n")){
			throw new RuntimeException("v3 sidecar column header line malformed/wrong order.");
		}
		System.err.println("Sidecar v3 header (schema_version=2, aligner+commit binding, no path leak, exact column header): PASS");

		final FamilyShortlistSidecar loaded=FamilyShortlistSidecar.load(dir+"/sidecar_v3.tsv",
			dir+"/consensus.fasta", dir+"/sets.tsv", "AAAligner");
		if(loaded.schemaVersion!=2){throw new RuntimeException("Loaded schemaVersion FAILED: "+loaded.schemaVersion);}
		//rank0 (rep0, not fallback): family beats global on r_lo/mindimer/minkmers/minid/minratio/
		//minscore; global beats family on mincovq/mincovt (family gave weaker 0.6 < floor 0.7).
		checkClose("rep0 r_lo", 0.8, loaded.rLo[0]); checkClose("rep0 r_hi", 1.2, loaded.rHi[0]);
		checkClose("rep0 minid", 40, loaded.minId[0]); checkClose("rep0 mincovq", 0.7, loaded.minCovQ[0]);
		if(loaded.fallback[0]){throw new RuntimeException("rep0 fallback FAILED: expected false.");}
		//rank1 (rep1, fallback): applied == global floor directly, ignoring the table's own 0.9/1.1/
		//etc -- including r_hi == Infinity (globalmaxlenratio was never set).
		checkClose("rep1 (fallback) r_lo", 0.75, loaded.rLo[1]);
		if(!Double.isInfinite(loaded.rHi[1])){throw new RuntimeException("rep1 (fallback) r_hi FAILED: expected Infinity, got "+loaded.rHi[1]);}
		checkClose("rep1 (fallback) minid", 30, loaded.minId[1]);
		if(!loaded.fallback[1]){throw new RuntimeException("rep1 fallback FAILED: expected true.");}
		if(loaded.n[1]!=5){throw new RuntimeException("rep1 n FAILED: expected 5, got "+loaded.n[1]);}
		//rank2 (rep2, not fallback): every family value is weaker than the global floor -- global
		//wins everywhere except r_hi (a CEILING: family's 2.0 beats an unset/infinite global ceiling).
		checkClose("rep2 r_lo", 0.75, loaded.rLo[2]); checkClose("rep2 r_hi", 2.0, loaded.rHi[2]);
		checkClose("rep2 minid", 30, loaded.minId[2]); checkClose("rep2 minscore", 15, loaded.minScore[2]);
		System.err.println("Sidecar v3 applied-threshold round-trip (family-wins/global-wins/fallback/Infinity ceiling): PASS");

		boolean threwAlignerMismatch=false;
		try{FamilyShortlistSidecar.load(dir+"/sidecar_v3.tsv", dir+"/consensus.fasta", dir+"/sets.tsv", "GlocalAmino");}
		catch(RuntimeException e){threwAlignerMismatch=true;}
		if(!threwAlignerMismatch){throw new RuntimeException("Aligner-mismatch crash-loud FAILED: did not throw.");}
		System.err.println("Aligner-mismatch crash-loud: PASS");

		//A v2 sidecar (no thresholds) has no aligner binding at all -- must load fine regardless of
		//expectedAligner (nothing to bind yet), preserving every pre-A2 caller's behavior.
		FamilyShortlistSidecar.load(dir+"/sidecar.tsv", dir+"/consensus.fasta", dir+"/sets.tsv", "AAAligner");
		FamilyShortlistSidecar.load(dir+"/sidecar.tsv", dir+"/consensus.fasta", dir+"/sets.tsv", "GlocalAmino");
		System.err.println("v2 sidecar ignores expectedAligner (no binding exists yet): PASS");

		writeFile(dir+"/thresholds_bad.tsv", THRESH_HEADER+"\n"+
			"rep0\t25\t0\t0.8\t1.2\t0.1\t2\t-40\t0.5\t20\t0.6\t0.6\n"+//negative minid -- must be refused
			"rep1\t5\t1\t0.9\t1.1\t0.2\t3\t50\t0.6\t25\t0.65\t0.65\n"+
			"rep2\t30\t0\t0.5\t2.0\t0.02\t0\t20\t0.1\t5\t0.5\t0.5\n");
		boolean threwNegative=false;
		try{
			new FamilyShortlistSidecarBuilder(new String[]{"consensus="+dir+"/consensus.fasta",
				"familiesdir="+famDir, "kmersets="+dir+"/sets.tsv", "out="+dir+"/sidecar_bad2.tsv",
				"thresholds="+dir+"/thresholds_bad.tsv", "aaalignerfile="+dir+"/fake_aaaligner.txt",
				"blosum62file="+dir+"/fake_blosum62.txt", "alignercommit=54b62510", "percentilep=2", "nmin=20"}).process();
		}catch(RuntimeException e){threwNegative=true;}
		if(!threwNegative){throw new RuntimeException("Negative-threshold crash-loud FAILED: did not throw.");}
		System.err.println("Negative-threshold-value crash-loud: PASS");

		//7. Column-header-line validation by NAME, not just count (UMP45's review, 2026-09-04): a
		//swapped pair of same-named-count columns must be refused, not silently loaded.
		final String swappedHeaderSidecar=v3Content.replaceFirst(
			"rank\trep_id\tmedian_len\tsd_ln_len", "rep_id\trank\tmedian_len\tsd_ln_len");
		writeFile(dir+"/sidecar_swapped.tsv", swappedHeaderSidecar);
		boolean threwSwapped=false;
		try{FamilyShortlistSidecar.load(dir+"/sidecar_swapped.tsv", dir+"/consensus.fasta", dir+"/sets.tsv", "AAAligner");}
		catch(RuntimeException e){threwSwapped=true;}
		if(!threwSwapped){throw new RuntimeException("Swapped-column-header crash-loud FAILED: did not throw.");}
		System.err.println("Swapped-column-header-name crash-loud: PASS");

		//8. CLI-shaped hang regression (UMP45's review, 2026-09-04): a REAL SUBPROCESS loading a
		//malformed schema-2 sidecar must exit nonzero within a bounded time. An in-process
		//try/catch (like every check above) cannot detect this bug class at all -- the test JVM's
		//main thread catches the exception fine either way; only watching a SEPARATE process's
		//actual exit proves FamilyShortlistSidecar.load() doesn't leave ByteFile4's non-daemon
		//producer/worker threads dangling (the real bug found 2026-09-04 via jstack: the row-parse
		//loop threw before reaching bf.close(), so the JVM printed its stack trace but never
		//exited). Corrupts rep0's minid (column 12, 0-based) to -999 in the valid v3 fixture --
		//any of the loader's many crash-loud checks would do; this one is simplest to construct.
		final String[] v3Lines=v3Content.split("\n");
		int rep0LineIdx=-1;
		for(int i=0; i<v3Lines.length; i++){if(v3Lines[i].startsWith("0\trep0\t")){rep0LineIdx=i; break;}}
		if(rep0LineIdx<0){throw new RuntimeException("Test setup FAILED: could not find rep0's data row in v3Content.");}
		final String[] rep0Cols=v3Lines[rep0LineIdx].split("\t", -1);
		rep0Cols[12]="-999";//minid column
		v3Lines[rep0LineIdx]=String.join("\t", rep0Cols);
		final StringBuilder corrupted=new StringBuilder();
		for(String ln : v3Lines){corrupted.append(ln).append('\n');}
		writeFile(dir+"/sidecar_hangtest.tsv", corrupted.toString());

		//System.getProperty("java.home")+"/bin/java", NOT ProcessHandle.current().info().command()
		//-- the latter is a Java-9+ API and BBTools must compile on Java 8 (bbtools skill's
		//"no lazy sysadmin left behind" rule; UMP45's review, 2026-09-04). Same precedent as
		//src/prot/IdentityGroupOutputVerifierTest.java's run().
		final String javaBin=System.getProperty("java.home")+"/bin/java";
		final ProcessBuilder pb=new ProcessBuilder(javaBin, "-ea", "--add-modules", "jdk.incubator.vector",
			"-cp", System.getProperty("java.class.path"), "prot.FamilyShortlistSidecar",
			dir+"/sidecar_hangtest.tsv", dir+"/consensus.fasta", dir+"/sets.tsv");
		pb.redirectErrorStream(true);
		final Process proc=pb.start();
		final boolean exitedInTime=proc.waitFor(30, java.util.concurrent.TimeUnit.SECONDS);
		if(!exitedInTime){
			proc.destroyForcibly();
			throw new RuntimeException("CLI hang regression FAILED: subprocess loading a malformed sidecar did "+
				"not exit within 30s -- this IS the bug (ByteFile4's non-daemon threads never got the close() "+
				"shutdown signal).");
		}
		if(proc.exitValue()==0){
			throw new RuntimeException("CLI hang regression FAILED: subprocess exited 0 -- expected a "+
				"crash-loud nonzero exit on the malformed (minid=-999) sidecar.");
		}
		System.err.println("CLI hang regression (real subprocess, malformed sidecar exits nonzero within 30s): PASS");

		System.err.println("FamilyShortlistSidecarBuilder selftest: ALL PASS");
	}

	private static void checkClose(final String what, final double expected, final double actual){
		if(Math.abs(expected-actual)>1e-9){
			throw new RuntimeException(what+" FAILED: expected "+expected+", got "+actual);
		}
	}

	private static void writeFile(final String path, final String content){
		final FileFormat ff=FileFormat.testOutput(path, FileFormat.TXT, null, false, true, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff); bsw.start();
		bsw.print(new ByteBuilder().append(content));
		if(bsw.poisonAndWait()){throw new RuntimeException("I/O error writing "+path);}
	}
	private static String readWhole(final String path){
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		final StringBuilder sb=new StringBuilder();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){sb.append(new String(line)).append('\n');}
		bf.close();
		return sb.toString();
	}

	private String consensusFile, familiesDir, kmerSetsFile, outFile;
	private String variant="centered";
	private double rawPseudocount=0.0;
	private String backgroundMode="families";

	//---- A2 (plans/PER_FAMILY_THRESHOLDS_v1.md): sidecar v3, all optional -- see the constructor's
	//thresholds!=null check for the fields this makes required. ----
	private String thresholdsFile, aaalignerFile, blosum62File, alignerCommit;
	private double percentileP=-1;
	private int nMin=-1;
	//Global floors: defaults match sec 2's own convention ("Each defaults to off (0/Infinity) so
	//triage=f reproduces today's assigner byte-for-byte") -- recorded in the header as the values
	//actually used at build time, never re-derived from a caller's runtime flags later.
	private double globalMinId=0, globalMinRatio=0, globalMinScore=0, globalMinCovQ=0, globalMinCovT=0,
		globalMinDimer=0, globalMinKmers=0, globalMinLenRatio=0, globalMaxLenRatio=Double.POSITIVE_INFINITY;
}
