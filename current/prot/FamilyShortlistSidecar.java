package prot;

import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.HashMap;

import fileIO.ByteFile;
import map.LongLongListHashMap;
import parse.LineParser1;
import structures.LongList;

/**
 * Loaded, ready-to-score representation of a {@link FamilyShortlistSidecarBuilder} artifact:
 * per-family composition vectors (flat row-major matrix for {@link simd.Vector#dotRows}), the
 * background frequency vector a live query must apply the same "centered" transform against,
 * a covering-set k-mer index (k-mer -&gt; family indices, for {@link ProteinSearcher}'s F4
 * shortlist scoring), and per-family length stats.
 *
 * <p><b>Hash-bound by construction.</b> {@link #load} re-hashes the LIVE consensus and
 * covering-set files and refuses to load if either does not match the sidecar's recorded
 * SHA-256 -- a stale sidecar can never silently score against a since-rebuilt consensus. This
 * check is baked into the loader (not left to callers) so it cannot be skipped.</p>
 *
 * @author Ady
 */
public final class FamilyShortlistSidecar {

	public final String[] repIds;
	public final int nFamilies, dims;
	/** Flat row-major composition matrix: family r's vector is compositionMatrix[r*dims .. r*dims+dims). */
	public final float[] compositionMatrix;
	/** Stored as double[] (not float[]) so a live query's per-call "centered" transform (which
	 *  operates in double precision, see {@link CompositionProfileAssay#applyVariant}) never
	 *  needs a fresh float-to-double conversion buffer -- UMP45's review, 2026-09-03. */
	public final double[] background;
	public final double[] medLen, sdLnLen;
	/** k-mer (packed, {@link ReducedAlphabetSeedAssay#kmerSet}'s encoding) -&gt; family indices. */
	public final LongLongListHashMap coveringIndex;
	/** |C_f|: the number of covering-set k-mers of family r (its row's {@code covering_kmers_hex} token
	 *  count), and their sum over all families. Derived at load, no format change. Used by {@link
	 *  ProteinSearcher#scoreFamiliesF4}'s {@code kmernorm=sub} candidate (UMP45, 2026-09-09,
	 *  plans/SHORTLIST_KMER_NORMALIZATION_ASSAY_v1.md §7): the query's expected coincidental k-mer
	 *  count against family r is proportional to |C_f|. */
	public final int[] coveringSetSize;
	public final long coveringSetSizeSum;
	public final ReducedAlphabet compositionAlphabet;
	/** Package-private {@link ReducedAlphabetSeedAssay.Alphabet} (a distinct type from {@link
	 *  ReducedAlphabet}, used specifically by {@link ReducedAlphabetSeedAssay#kmerSet}) -- the
	 *  covering-set module's own alphabet representation, not interchangeable with the composition
	 *  module's {@link ReducedAlphabet}; both exist in this codebase for their respective callers. */
	public final ReducedAlphabetSeedAssay.Alphabet coveringAlphabet;
	public final int compositionK, coveringK;
	public final String variant;
	public final double rawPseudocount;
	public final String consensusSha256, kmerSetsSha256;

	/** Schema version recorded in the sidecar (1 = v2, no per-family thresholds; 2 = v3, thresholds
	 *  present -- plans/PER_FAMILY_THRESHOLDS_v1.md WP-A A2, added 2026-09-04). */
	public final int schemaVersion;
	/** Non-null only when {@link #schemaVersion}&gt;=2. Header binding: the aligner name (e.g.
	 *  "AAAligner") every applied threshold below was calibrated against, plus sha256 of the
	 *  aligner/Blosum62 source and the threshold table itself, at the pinned commit. */
	public final String aligner, alignerCommit, alignerSha256, blosum62Sha256, thresholdsSha256;
	/** Non-null only when {@link #schemaVersion}&gt;=2: the chosen calibration percentile and the
	 *  minimum calibration-set size below which a family falls back to the global floor alone. */
	public final Double percentileP; public final Integer nMin;
	/** Global floors actually used when this sidecar's thresholds were built (values, not just a
	 *  hash) -- non-null only when {@link #schemaVersion}&gt;=2. */
	public final Double globalMinId, globalMinRatio, globalMinScore, globalMinCovQ, globalMinCovT,
		globalMinDimer, globalMinKmers, globalMinLenRatio, globalMaxLenRatio;
	/** Per-family APPLIED thresholds (max/min of the family's own percentile and the global floor,
	 *  sec 2: "the stricter wins"; a fallback family gets the global floor directly) -- non-null
	 *  only when {@link #schemaVersion}&gt;=2, indexed by family rank matching {@link #repIds}. */
	public final int[] n; public final boolean[] fallback;
	public final double[] rLo, rHi, minDimer, minKmers, minId, minRatio, minScore, minCovQ, minCovT;

	private FamilyShortlistSidecar(final String[] repIds_, final int dims_, final float[] compositionMatrix_,
			final double[] background_, final double[] medLen_, final double[] sdLnLen_,
			final LongLongListHashMap coveringIndex_, final ReducedAlphabet compositionAlphabet_,
			final int compositionK_, final ReducedAlphabetSeedAssay.Alphabet coveringAlphabet_, final int coveringK_,
			final String variant_, final double rawPseudocount_, final String consensusSha256_,
			final String kmerSetsSha256_, final int schemaVersion_, final String aligner_, final String alignerCommit_,
			final String alignerSha256_, final String blosum62Sha256_, final String thresholdsSha256_,
			final Double percentileP_, final Integer nMin_, final Double globalMinId_, final Double globalMinRatio_,
			final Double globalMinScore_, final Double globalMinCovQ_, final Double globalMinCovT_,
			final Double globalMinDimer_, final Double globalMinKmers_, final Double globalMinLenRatio_,
			final Double globalMaxLenRatio_, final int[] n_, final boolean[] fallback_, final double[] rLo_,
			final double[] rHi_, final double[] minDimer_, final double[] minKmers_, final double[] minId_,
			final double[] minRatio_, final double[] minScore_, final double[] minCovQ_, final double[] minCovT_,
			final int[] coveringSetSize_){
		repIds=repIds_; nFamilies=repIds_.length; dims=dims_; compositionMatrix=compositionMatrix_;
		background=background_; medLen=medLen_; sdLnLen=sdLnLen_; coveringIndex=coveringIndex_;
		coveringSetSize=coveringSetSize_;
		long sum=0; for(int v : coveringSetSize_){sum+=v;}
		coveringSetSizeSum=sum;
		//No positivity assertion here: the loader has always accepted empty covering sets (a family row with an
		//empty covering_kmers_hex column), and kmernorm=none must keep that behavior. An all-empty sidecar
		//(coveringSetSizeSum==0) is handled by ProteinSearcher.scoreFamiliesF4 as the zero-background case.
		assert(coveringSetSize_.length==nFamilies) : "one covering-set size per family: "+coveringSetSize_.length+" vs "+nFamilies+" (load() counts cols[5] tokens per row)";
		compositionAlphabet=compositionAlphabet_; compositionK=compositionK_;
		coveringAlphabet=coveringAlphabet_; coveringK=coveringK_; variant=variant_;
		rawPseudocount=rawPseudocount_; consensusSha256=consensusSha256_; kmerSetsSha256=kmerSetsSha256_;
		schemaVersion=schemaVersion_; aligner=aligner_; alignerCommit=alignerCommit_; alignerSha256=alignerSha256_;
		blosum62Sha256=blosum62Sha256_; thresholdsSha256=thresholdsSha256_; percentileP=percentileP_;
		nMin=nMin_; globalMinId=globalMinId_; globalMinRatio=globalMinRatio_; globalMinScore=globalMinScore_;
		globalMinCovQ=globalMinCovQ_; globalMinCovT=globalMinCovT_; globalMinDimer=globalMinDimer_;
		globalMinKmers=globalMinKmers_; globalMinLenRatio=globalMinLenRatio_; globalMaxLenRatio=globalMaxLenRatio_;
		n=n_; fallback=fallback_; rLo=rLo_; rHi=rHi_; minDimer=minDimer_; minKmers=minKmers_; minId=minId_;
		minRatio=minRatio_; minScore=minScore_; minCovQ=minCovQ_; minCovT=minCovT_;
	}

	/**
	 * Loads a sidecar TSV, verifying its recorded consensus/covering-set SHA-256 against the
	 * LIVE files (crash loud on any mismatch -- a stale sidecar must never silently score
	 * against a since-rebuilt consensus). Equivalent to {@link #load(String, String, String,
	 * String)} with a null {@code expectedAligner} (no aligner-binding check -- unchanged
	 * behavior for every pre-A2 caller).
	 * @param sidecarPath The sidecar TSV path (from {@link FamilyShortlistSidecarBuilder}).
	 * @param liveConsensusFile The live consensus FASTA this run is actually searching against.
	 * @param liveKmerSetsFile The live covering-set file (informational cross-check; the sidecar
	 *        already bundles the k-mers themselves, so this only detects a stale/rebuilt file).
	 * @return The loaded, verified sidecar.
	 */
	public static FamilyShortlistSidecar load(final String sidecarPath, final String liveConsensusFile,
			final String liveKmerSetsFile){
		return load(sidecarPath, liveConsensusFile, liveKmerSetsFile, null);
	}

	/** CLI entry point for {@link #selftest}'s CLI-shaped hang-regression check ONLY -- calls
	 *  {@link #load} from a genuinely separate process so a caller watching for PROCESS EXIT (not
	 *  an in-process try/catch, which cannot detect a dangling non-daemon thread) can verify the
	 *  crash-loud contract actually terminates the JVM. Not a general user-facing tool.
	 * @param args {@code <sidecarPath> <consensusFile> <kmerSetsFile>}. */
	public static void main(final String[] args){
		if(args.length!=3){
			throw new RuntimeException("Usage: prot.FamilyShortlistSidecar <sidecar> <consensus> <kmersets> "+
				"(selftest-only entry point for the CLI hang-regression check).");
		}
		load(args[0], args[1], args[2]);
		System.err.println("Loaded OK (unexpected for the hang-regression test's malformed fixture).");
	}

	/**
	 * As {@link #load(String, String, String)}, plus (2026-09-04, plans/PER_FAMILY_THRESHOLDS_v1.md
	 * WP-A A2) an aligner-binding check: when {@code expectedAligner} is non-null, a v3 (schema
	 * &gt;=2) sidecar whose recorded {@code #aligner} does not equal it is refused loud -- per-family
	 * thresholds calibrated against one aligner's score distribution must never be silently applied
	 * to hits from a different one. A v2 (schema 1) sidecar has no aligner binding at all and is
	 * accepted regardless of {@code expectedAligner} (thresholds don't exist yet to bind).
	 * @param expectedAligner The aligner name the caller is about to run hits through (e.g.
	 *        "AAAligner"), or null to skip the check entirely.
	 * @return The loaded, verified sidecar.
	 */
	public static FamilyShortlistSidecar load(final String sidecarPath, final String liveConsensusFile,
			final String liveKmerSetsFile, final String expectedAligner){
		final HashMap<String, String> header=new HashMap<String, String>();
		final ByteFile bf=ByteFile.makeByteFile(sidecarPath, false);
		String columnHeaderLine=null;
		int dims=-1;
		int schemaVersion=1;//pre-A2 sidecars never wrote this key; default preserves old behavior
		final java.util.ArrayList<String> repIdList=new java.util.ArrayList<String>();
		final java.util.ArrayList<Double> medLenList=new java.util.ArrayList<Double>();
		final java.util.ArrayList<Double> sdLnLenList=new java.util.ArrayList<Double>();
		final java.util.ArrayList<float[]> compRows=new java.util.ArrayList<float[]>();
		final java.util.ArrayList<Integer> setSizeList=new java.util.ArrayList<Integer>();
		final LongLongListHashMap coveringIndex=new LongLongListHashMap();
		double[] background=null;
		int expectedRank=0;
		final java.util.ArrayList<Integer> nList=new java.util.ArrayList<Integer>();
		final java.util.ArrayList<Boolean> fallbackList=new java.util.ArrayList<Boolean>();
		final java.util.ArrayList<Double> rLoList=new java.util.ArrayList<Double>(), rHiList=new java.util.ArrayList<Double>(),
			minDimerList=new java.util.ArrayList<Double>(), minKmersList=new java.util.ArrayList<Double>(),
			minIdList=new java.util.ArrayList<Double>(), minRatioList=new java.util.ArrayList<Double>(),
			minScoreList=new java.util.ArrayList<Double>(), minCovQList=new java.util.ArrayList<Double>(),
			minCovTList=new java.util.ArrayList<Double>();
		//try/finally around the whole parse loop (2026-09-04, found running A3's gate 2/4: a
		//RuntimeException thrown mid-loop -- by design, every parseThreshold* method above throws
		//on bad data -- used to skip this bf.close() entirely. ByteFile's underlying ByteFile4
		//reader runs its own non-daemon producer/worker thread pool (confirmed via jstack on a
		//real hung process: "ByteFile4-Producer"/"ByteFile4-Worker-*", parked, not daemon); with
		//bf never closed, those threads never got the shutdown signal, and the JVM never exited --
		//a crash printed its trace to stderr but the PROCESS HUNG FOREVER instead of terminating.
		//This is the exact "crash loud, never hang" violation the assertDie idiom exists to
		//prevent, except here it was a plain uncaught exception, not an assertion, and every
		//crash-loud path already in this loader (aligner mismatch, malformed row, negative
		//threshold, swapped header, family-order mismatch, ...) was equally exposed -- this fix
		//covers all of them at once, not just the mindimer case that surfaced it.
		try{
			for(byte[] lineBytes=bf.nextLine(); lineBytes!=null; lineBytes=bf.nextLine()){
				if(lineBytes.length==0){continue;}
				final String line=new String(lineBytes);
				if(line.charAt(0)=='#'){
					final int tab=line.indexOf('\t');
					if(tab<0){continue;}
					final String key=line.substring(1, tab);
					String value=line.substring(tab+1);
					if(key.endsWith("_sha256") || key.endsWith("_sha80")){
						if(key.endsWith("_sha256") && value.length()!=64){throw new RuntimeException("Legacy sidecar digest width is invalid: "+key);}
						if(key.endsWith("_sha80") && value.length()!=DigestSuffix.HEX_LENGTH){throw new RuntimeException("Sidecar sha80 width is invalid: "+key);}
						value=DigestSuffix.normalizeRecorded(value,key);
					}
					if(key.startsWith("background_")){
						final String[] parts=value.split(",");
						background=new double[parts.length];
						for(int i=0; i<parts.length; i++){background[i]=Double.parseDouble(parts[i]);}
					}else{
						header.put(key, value);
						if(key.equals("schema_version")){schemaVersion=Integer.parseInt(value);}
					}
					continue;
				}
				if(columnHeaderLine==null){
					columnHeaderLine=line;
					validateColumnHeader(columnHeaderLine, schemaVersion, sidecarPath);
					continue;
				}
				final String[] cols=line.split("\t", -1);
				final int expectedCols=(schemaVersion>=2 ? 17 : 6);
				if(cols.length!=expectedCols){
					throw new RuntimeException("Malformed sidecar data row (expected "+expectedCols+
						" columns for schema_version="+schemaVersion+", got "+cols.length+"): "+line);
				}
				final int rank=Integer.parseInt(cols[0]);
				if(rank!=expectedRank){
					throw new RuntimeException("Sidecar family order mismatch: expected rank "+expectedRank+
						", got "+rank+".");
				}
				final String repId=cols[1];
				final String[] compTokens=cols[4].split(",");
				if(dims<0){dims=compTokens.length;}
				else if(dims!=compTokens.length){
					throw new RuntimeException("Sidecar composition dimension mismatch at rank "+rank+": "+
						compTokens.length+" vs expected "+dims+".");
				}
				final float[] vec=new float[dims];
				for(int d=0; d<dims; d++){vec[d]=Float.parseFloat(compTokens[d]);}
				repIdList.add(repId);
				medLenList.add(Double.parseDouble(cols[2]));
				sdLnLenList.add(Double.parseDouble(cols[3]));
				compRows.add(vec);
				if(cols[5].length()>0){
					final String[] kmerTokens=cols[5].split(",");
					for(String hex : kmerTokens){coveringIndex.put(Long.parseUnsignedLong(hex, 16), (long)rank);}
					setSizeList.add(kmerTokens.length);
				}else{
					setSizeList.add(0);
				}
				if(schemaVersion>=2){
					nList.add(parseThresholdInt(cols[6], "n", repId, sidecarPath));
					fallbackList.add(parseThresholdBool(cols[7], "fallback", repId, sidecarPath));
					rLoList.add(parseThresholdDouble(cols[8], "r_lo", repId, sidecarPath));
					rHiList.add(parseThresholdCeiling(cols[9], "r_hi", repId, sidecarPath));
					//DO NOT "fix" this back to 0.0 -- mindimer scores a composition-vector COSINE
					//SIMILARITY, whose real range is [-1,1] (unlike every other metric column
					//here, which is a count/ratio/percentage bounded below by 0). -1 is that
					//range's actual floor, and it is the value the sidecar's OWN header records
					//as the global "off" default (#global_mindimer, parsed below into
					//globalMinDimer) -- a genuinely no-op sidecar (plans/PER_FAMILY_THRESHOLDS_v1.md
					//sec 2's "off" default) sets #global_mindimer=-1, and this loader must accept
					//that value or no no-op sidecar can ever load. Confirmed 2026-09-04 building
					//A3's gate 2: the OLD non-negative-only parseThresholdDouble refused this
					//legitimate value -- a real A2 loader gap (UMP45 confirmed, 2026-09-04), not a
					//gate-script bug.
					minDimerList.add(parseThresholdDouble(cols[10], "mindimer", repId, sidecarPath, -1.0));
					minKmersList.add(parseThresholdDouble(cols[11], "minkmers", repId, sidecarPath));
					minIdList.add(parseThresholdDouble(cols[12], "minid", repId, sidecarPath));
					minRatioList.add(parseThresholdDouble(cols[13], "minratio", repId, sidecarPath));
					minScoreList.add(parseThresholdDouble(cols[14], "minscore", repId, sidecarPath));
					minCovQList.add(parseThresholdDouble(cols[15], "mincovq", repId, sidecarPath));
					minCovTList.add(parseThresholdDouble(cols[16], "mincovt", repId, sidecarPath));
				}
				expectedRank++;
			}
		}finally{
			bf.close();
		}
		if(repIdList.isEmpty()){throw new RuntimeException("No data rows found in sidecar: "+sidecarPath);}
		if(background==null){
			throw new RuntimeException("Sidecar has no #background_Nf header line: "+sidecarPath+
				" -- this sidecar predates the background-vector fix and cannot be used for scoring.");
		}
		if(background.length!=dims){
			throw new RuntimeException("Sidecar background vector dims ("+background.length+
				") != composition dims ("+dims+").");
		}

		final int nFamilies=repIdList.size();
		final String[] repIds=repIdList.toArray(new String[0]);
		final double[] medLen=new double[nFamilies], sdLnLen=new double[nFamilies];
		final float[] compositionMatrix=new float[nFamilies*dims];
		final int[] coveringSetSize=new int[nFamilies];
		for(int r=0; r<nFamilies; r++){
			medLen[r]=medLenList.get(r); sdLnLen[r]=sdLnLenList.get(r);
			System.arraycopy(compRows.get(r), 0, compositionMatrix, r*dims, dims);
			coveringSetSize[r]=setSizeList.get(r);
		}

		final String recordedConsensusSha256=DigestSuffix.requiredCompatibleHeader(header,"consensus",sidecarPath);
		final String recordedKmerSetsSha256=DigestSuffix.requiredCompatibleHeader(header,"kmersets",sidecarPath);
		final String liveConsensusSha256=MagQCTextResource.sha80(liveConsensusFile);
		if(!liveConsensusSha256.equals(recordedConsensusSha256)){
			throw new RuntimeException("Sidecar consensus SHA-256 mismatch: sidecar was built from a "+
				"consensus with sha256="+recordedConsensusSha256+", but the live file "+liveConsensusFile+
				" has sha256="+liveConsensusSha256+" -- the consensus has changed since this sidecar was "+
				"built; rebuild the sidecar before using it.");
		}
		final String liveKmerSetsSha256=MagQCTextResource.sha80(liveKmerSetsFile);
		if(!liveKmerSetsSha256.equals(recordedKmerSetsSha256)){
			throw new RuntimeException("Sidecar covering-set SHA-256 mismatch: sidecar was built from "+
				"kmersets with sha256="+recordedKmerSetsSha256+", but the live file "+liveKmerSetsFile+
				" has sha256="+liveKmerSetsSha256+" -- rebuild the sidecar before using it.");
		}

		final ReducedAlphabet compositionAlphabet=ReducedAlphabet.named(required(header, "composition_alphabet", sidecarPath));
		final int compositionK=Integer.parseInt(required(header, "composition_k", sidecarPath));
		final ReducedAlphabetSeedAssay.Alphabet coveringAlphabet=
			ReducedAlphabetSeedAssay.Alphabet.named(required(header, "covering_alphabet", sidecarPath));
		final int coveringK=Integer.parseInt(required(header, "covering_k", sidecarPath));
		final String variant=required(header, "composition_variant", sidecarPath);
		final double rawPseudocount=Double.parseDouble(required(header, "composition_rawpseudo", sidecarPath));

		String aligner=null, alignerCommit=null, alignerSha256=null, blosum62Sha256=null, thresholdsSha256=null;
		Double percentileP=null, globalMinId=null, globalMinRatio=null, globalMinScore=null, globalMinCovQ=null,
			globalMinCovT=null, globalMinDimer=null, globalMinKmers=null, globalMinLenRatio=null, globalMaxLenRatio=null;
		Integer nMin=null;
		int[] n=null; boolean[] fallback=null;
		double[] rLo=null, rHi=null, minDimer=null, minKmers=null, minId=null, minRatio=null, minScore=null,
			minCovQ=null, minCovT=null;
		if(schemaVersion>=2){
			aligner=required(header, "aligner", sidecarPath);
			alignerCommit=required(header, "aligner_commit", sidecarPath);
			alignerSha256=DigestSuffix.requiredCompatibleHeader(header,"aligner",sidecarPath);
			blosum62Sha256=DigestSuffix.requiredCompatibleHeader(header,"blosum62",sidecarPath);
			thresholdsSha256=DigestSuffix.requiredCompatibleHeader(header,"thresholds",sidecarPath);
			percentileP=Double.parseDouble(required(header, "percentile_p", sidecarPath));
			nMin=Integer.parseInt(required(header, "n_min", sidecarPath));
			globalMinId=Double.parseDouble(required(header, "global_minid", sidecarPath));
			globalMinRatio=Double.parseDouble(required(header, "global_minratio", sidecarPath));
			globalMinScore=Double.parseDouble(required(header, "global_minscore", sidecarPath));
			globalMinCovQ=Double.parseDouble(required(header, "global_mincovq", sidecarPath));
			globalMinCovT=Double.parseDouble(required(header, "global_mincovt", sidecarPath));
			globalMinDimer=Double.parseDouble(required(header, "global_mindimer", sidecarPath));
			globalMinKmers=Double.parseDouble(required(header, "global_minkmers", sidecarPath));
			globalMinLenRatio=Double.parseDouble(required(header, "global_minlenratio", sidecarPath));
			globalMaxLenRatio=Double.parseDouble(required(header, "global_maxlenratio", sidecarPath));
			if(expectedAligner!=null && !expectedAligner.equals(aligner)){
				throw new RuntimeException("Sidecar aligner mismatch: thresholds in "+sidecarPath+
					" were calibrated against aligner='"+aligner+"', but the caller is running aligner='"+
					expectedAligner+"' -- per-family thresholds from one aligner's score distribution "+
					"must never be applied to another's hits.");
			}
			n=new int[nFamilies]; fallback=new boolean[nFamilies];
			rLo=new double[nFamilies]; rHi=new double[nFamilies]; minDimer=new double[nFamilies];
			minKmers=new double[nFamilies]; minId=new double[nFamilies]; minRatio=new double[nFamilies];
			minScore=new double[nFamilies]; minCovQ=new double[nFamilies]; minCovT=new double[nFamilies];
			for(int r=0; r<nFamilies; r++){
				n[r]=nList.get(r); fallback[r]=fallbackList.get(r); rLo[r]=rLoList.get(r); rHi[r]=rHiList.get(r);
				minDimer[r]=minDimerList.get(r); minKmers[r]=minKmersList.get(r); minId[r]=minIdList.get(r);
				minRatio[r]=minRatioList.get(r); minScore[r]=minScoreList.get(r); minCovQ[r]=minCovQList.get(r);
				minCovT[r]=minCovTList.get(r);
			}
		}

		return new FamilyShortlistSidecar(repIds, dims, compositionMatrix, background, medLen, sdLnLen,
			coveringIndex, compositionAlphabet, compositionK, coveringAlphabet, coveringK, variant,
			rawPseudocount, recordedConsensusSha256, recordedKmerSetsSha256, schemaVersion, aligner, alignerCommit,
			alignerSha256, blosum62Sha256, thresholdsSha256, percentileP, nMin, globalMinId, globalMinRatio,
			globalMinScore, globalMinCovQ, globalMinCovT, globalMinDimer, globalMinKmers, globalMinLenRatio,
			globalMaxLenRatio, n, fallback, rLo, rHi, minDimer, minKmers, minId, minRatio, minScore, minCovQ, minCovT,
			coveringSetSize);
	}

	/** Refuses a non-integer or negative {@code n} (calibration-set size) -- sec 3 A2's loader
	 *  contract ("refuses NaN/negative"). */
	private static int parseThresholdInt(final String s, final String field, final String repId, final String path){
		final int v;
		try{v=Integer.parseInt(s);}catch(NumberFormatException e){
			throw new RuntimeException("Sidecar "+path+", family '"+repId+"': non-integer "+field+"='"+s+"'.");
		}
		if(v<0){throw new RuntimeException("Sidecar "+path+", family '"+repId+"': negative "+field+"="+v+".");}
		return v;
	}
	private static boolean parseThresholdBool(final String s, final String field, final String repId, final String path){
		if(s.equals("0")){return false;}
		if(s.equals("1")){return true;}
		throw new RuntimeException("Sidecar "+path+", family '"+repId+"': "+field+" must be 0 or 1, got '"+s+"'.");
	}
	/** Refuses NaN, non-finite (Infinity), or below-floor applied thresholds -- sec 3 A2's loader
	 *  contract. An applied per-family FLOOR is always a concrete finite number for a non-fallback
	 *  family (max'd against a finite family percentile); Infinity here means a corrupted or
	 *  hand-edited table, not a legitimate "no cutoff" state. ({@code r_hi}, the one CEILING
	 *  metric, uses {@link #parseThresholdCeiling} instead -- it legitimately stays Infinity for a
	 *  fallback family when {@code globalmaxlenratio=} was left at its "off" default.)
	 *  Equivalent to {@link #parseThresholdDouble(String, String, String, String, double)} with
	 *  {@code minAllowed=0} -- every field except {@code mindimer} has a naturally non-negative
	 *  range (length ratio, k-mer count, identity, score ratio, raw score, coverage fraction). */
	private static double parseThresholdDouble(final String s, final String field, final String repId, final String path){
		return parseThresholdDouble(s, field, repId, path, 0.0);
	}
	/** As the 4-arg overload, but with an explicit lower bound -- {@code mindimer} (2026-09-04,
	 *  found building A3's gate 2 no-op sidecar) is a composition-cosine score, whose natural
	 *  range includes negatives; its genuinely-off floor is -1, not 0. */
	private static double parseThresholdDouble(final String s, final String field, final String repId, final String path, final double minAllowed){
		final double v;
		try{v=Double.parseDouble(s);}catch(NumberFormatException e){
			throw new RuntimeException("Sidecar "+path+", family '"+repId+"': non-numeric "+field+"='"+s+"'.");
		}
		if(Double.isNaN(v) || Double.isInfinite(v) || v<minAllowed){
			throw new RuntimeException("Sidecar "+path+", family '"+repId+"': "+field+"="+v+
				" must be finite and >= "+minAllowed+".");
		}
		return v;
	}
	/** As {@link #parseThresholdDouble}, but allows {@code +Infinity} for {@code r_hi} -- sec 2's
	 *  own "off" default for the length-ratio CEILING is {@code Infinity} (no upper bound), and a
	 *  fallback family (sec 2: "gets the global floors only") can genuinely carry that default
	 *  through into the sidecar when no real ceiling was ever set. Negative or {@code -Infinity}
	 *  still refused -- neither is ever a legitimate r_hi. */
	private static double parseThresholdCeiling(final String s, final String field, final String repId, final String path){
		final double v;
		try{v=Double.parseDouble(s);}catch(NumberFormatException e){
			throw new RuntimeException("Sidecar "+path+", family '"+repId+"': non-numeric "+field+"='"+s+"'.");
		}
		if(Double.isNaN(v) || v<0 || v==Double.NEGATIVE_INFINITY){
			throw new RuntimeException("Sidecar "+path+", family '"+repId+"': "+field+"="+v+
				" must be non-negative (Infinity allowed, meaning no upper bound).");
		}
		return v;
	}

	/** Validates the column header line BY NAME, not just by count (UMP45's review, 2026-09-04:
	 *  sec 3 A2 says the loader "validates column count/order" -- a count-only check would load a
	 *  swapped pair of same-count columns silently). The composition column's width is embedded in
	 *  its own name ({@code composition_<N>f}) and is not fixed, so that one field is checked by
	 *  pattern rather than exact match; every other field name is fixed and checked exactly. */
	private static void validateColumnHeader(final String line, final int schemaVersion, final String path){
		final String[] f=line.split("\t", -1);
		final int expectedCount=(schemaVersion>=2 ? 17 : 6);
		if(f.length!=expectedCount){
			throw new RuntimeException("Sidecar column header has "+f.length+" fields, expected "+expectedCount+
				" for schema_version="+schemaVersion+": "+line);
		}
		final String[] fixed={"rank","rep_id","median_len","sd_ln_len"};
		for(int i=0; i<fixed.length; i++){
			if(!f[i].equals(fixed[i])){
				throw new RuntimeException("Sidecar column header field "+i+" is '"+f[i]+"', expected '"+
					fixed[i]+"': "+path);
			}
		}
		if(!f[4].matches("composition_\\d+f")){
			throw new RuntimeException("Sidecar column header field 4 is '"+f[4]+
				"', expected 'composition_<N>f': "+path);
		}
		if(!f[5].equals("covering_kmers_hex")){
			throw new RuntimeException("Sidecar column header field 5 is '"+f[5]+
				"', expected 'covering_kmers_hex': "+path);
		}
		if(schemaVersion>=2){
			final String[] thresh={"n","fallback","r_lo","r_hi","mindimer","minkmers","minid","minratio",
				"minscore","mincovq","mincovt"};
			for(int i=0; i<thresh.length; i++){
				if(!f[6+i].equals(thresh[i])){
					throw new RuntimeException("Sidecar column header field "+(6+i)+" is '"+f[6+i]+
						"', expected '"+thresh[i]+"': "+path);
				}
			}
		}
	}

	private static String required(final HashMap<String, String> header, final String key, final String path){
		final String v=header.get(key);
		if(v==null){throw new RuntimeException("Sidecar missing required header key '"+key+"': "+path);}
		return v;
	}

}
