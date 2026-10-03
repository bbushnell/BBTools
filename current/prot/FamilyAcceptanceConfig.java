package prot;

import java.io.FileInputStream;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.security.DigestInputStream;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import fileIO.ByteFile;
import structures.ByteBuilder;

/**
 * In-memory per-family schema-4 D107 plus raw-HBM acceptance cutoffs for the canonical assigner ({@link
 * ProteinSearcher#assignFamily}, Patch 2). Every per-family array is indexed by family rank,
 * aligned 1:1 with the roster this artifact was built from (and therefore with any {@link
 * FamilyShortlistSidecar}/target list built from the same roster). Loads the {@code
 * d55_d107_hbm_family_thresholds} artifact format ({@code
 * /mnt/c/playground/UMP45/plans/d55_d107_threshold_artifact_schema_ump45_20260914.md}, draft v2),
 * produced by {@code calibrationreview.RuntimeThresholdExporter}.
 *
 * <p><b>Immutability.</b> Unlike {@link FamilyShortlistSidecar}'s public-array convention, every
 * per-family array here is PRIVATE and defensively copied at construction; callers read one
 * family's values through the {@code xxx(rank)} accessors, never through a shared mutable array
 * reference (Yoimiya's review, 2026-09-14, correcting UMP45's original Patch 2 draft, which had
 * these as public arrays). Scalar {@code String} fields stay public final -- Strings are
 * immutable, so there is no aliasing hazard to hide.</p>
 *
 * <p><b>Hash-bound by construction, like the sidecar.</b> {@link #load} re-hashes the LIVE family
 * list and consensus reference files and refuses to load if either does not match the artifact's
 * recorded SHA-256; it also re-derives the family RANK ORDER from the live family list (never
 * trusting the artifact's own row order in isolation) and cross-checks every data row's own
 * {@code rank}/{@code rep_id} against it, and recomputes {@code thresholds_sha256} over the data
 * rows actually read, refusing a hand-edited or corrupted artifact. {@code aligner} must be
 * literally {@code D55} (this artifact type is D55-specific, not a generic multi-aligner
 * container) and must match the caller's {@code expectedAligner}; {@code aligner_sha256}/{@code
 * blosum62_sha256} must match the caller-supplied values (this loader has no independent way to
 * re-hash "the aligner's source" from inside a running JVM, so the caller is responsible for
 * supplying it). <b>Source hash and compiled-class hash are NOT interchangeable</b> (Yoimiya's
 * correction, 2026-09-14: earlier wording here wrongly suggested a compiled-class hash would do
 * -- the schema's own {@code #aligner_sha256} key is explicitly the ALIGNER SOURCE hash; the
 * caller must supply that exact provenance type, never a class-file hash of the same aligner).</p>
 *
 * <p><b>Schema v4 requires HBM.</b> The loader accepts only {@code #hbm_enabled 1}, binds the
 * live MQHB bundle and its 16-row semantic-provenance file by SHA-256, and exposes exactly one
 * canonical per-family {@code min_hbm_raw_99} cutoff. Raw HBM is an inclusive post-base-filter
 * gate; it is not a cross-family winner metric (D118). Older HBM-off schemas fail at load time.</p>
 *
 * <p><b>Schema v2 (Yoimiya, 2026-09-14): kmer/covering-set binding.</b> {@code min_kmer_count}/
 * {@code min_kmer_qdensity} are only meaningful under the exact covering-set k-mer configuration
 * they were calibrated under. {@link #load} therefore ALSO reads the LIVE {@code
 * FamilyShortlistSidecar}'s own header (never trusting the artifact's claims about it in
 * isolation) and cross-checks: the sidecar's {@code consensus_sha256} against the same live
 * consensus hash already verified above; a re-hash of the live covering-set k-mer file against
 * both the sidecar's own {@code kmersets_sha256} AND the artifact's {@code
 * #covering_sets_sha256}; a re-hash of the live sidecar FILE itself against the artifact's {@code
 * #shortlist_sidecar_sha256}; and the sidecar's {@code covering_alphabet}/{@code covering_k}
 * against the artifact's {@code #covering_alphabet}/{@code #covering_k}. {@code #kmer_boundary}
 * must be literally {@code NONE} (the only boundary mode production code uses). {@code
 * #kmer_definition} has no live file to re-derive from (it names an algorithm, not a resource) and
 * is therefore caller-supplied via {@code expectedKmerDefinition} -- an explicit seam, since Sayu
 * is still defining its canonical stable form as of 2026-09-14; never a guessed value. Schema v1
 * artifacts (pre-dating this binding) are refused outright. These schema-v2 checks remain part of
 * the current schema-v3 loader; v1/v2 outputs already produced are preserved on disk as historical
 * records, not deleted, just no longer loadable by this class.</p>
 *
 * <p><b>Schema v3 (Yoimiya, 2026-09-16): role binding.</b> The artifact records
 * {@code #role_manifest_sha256}; {@link #load} re-hashes the live role manifest,
 * loads it through {@link FamilyRoleManifest}, and requires exact rank/rep-ID
 * agreement with the already-hash-bound roster. Schema v2 is now historical and
 * refused: otherwise a bait rank could be mistaken for a tracked network feature.
 * A tracked-only baseline uses the same schema-v3 role manifest with every row
 * declared {@code tracked,true}; there is no role-free production path.</p>
 *
 * @author Eru, 2026-09-14
 */
public final class FamilyAcceptanceConfig {

	/** Family count; every per-family accessor's {@code rank} argument must be in {@code [0,nFamilies)}. */
	public final int nFamilies;
	/** Runtime acceptance schema: 7 for the sealed production table, 0 for legacy schema 4/5/6 profiles. */
	public final int acceptanceSchema;
	private final String[] repIds;
	private final int[] familyIds;
	private final boolean[] acceptanceEnabled;
	private final int[] basePositiveN, hbmUniqueN;
	private final int[] lqLo, lqHi;
	private final int[] minRawScore;
	private final double[] minId, minOverlap;
	private final double[] minR, minCoreMutualOverlap;
	private final int[] minKmerCount;
	private final double[] minKmerQdensity;
	private final FamilyRoleManifest roles;
	/** Strict schema 4 always enables the bound raw-HBM acceptance gate. */
	public final boolean hbmEnabled=true;
	private final float[] minHbmRaw99;
	private final float[] minHbmPathRelative;
	public final String aligner, familyListSha256, consensusRefSha256, roleManifestSha256, thresholdsSha256;
	/** Hash of the complete sealed threshold artifact, or {@code NA} for legacy profiles. */
	public final String profileArtifactSha80;
	public final String coreCoordinatesSha80;
	/** Semantic predecessor identity bound by schema-7 core coordinates. */
	public final String coreFamilyListSha80;
	public final String hbmMetric, hbmCutoffRule, hbmScoreContract, hbmBundleSha256,
		hbmSemanticProvenanceSha256, hbmLandmarksSha256, hbmLandmarksProvenanceSha256,
		hbmPass4ManifestSha256, hbmFamilyHashesSha256;
	public final int hbmTargetBp, hbmFamilyCount, hbmMaxQueryLength;
	public final long hbmMemberTotal;
	/** Validated resource metadata, retained so a caller (e.g. UMP45's runtime bound-setup) can
	 *  compare this config to its actual live sidecar/contract without re-parsing the artifact
	 *  file itself (Yoimiya's review, 2026-09-14). All already cross-checked during {@link #load}. */
	public final String alignerSha256, blosum62Sha256, coveringSetsSha256, shortlistSidecarSha256,
		coveringAlphabet, kmerBoundary, kmerDefinition;
	public final int coveringK;

	private FamilyAcceptanceConfig(final String[] repIds_, final int[] basePositiveN_, final int[] lqLo_, final int[] lqHi_,
			final int[] minRawScore_, final double[] minId_, final double[] minOverlap_, final int[] minKmerCount_,
			final double[] minKmerQdensity_, final int[] hbmUniqueN_, final float[] minHbmRaw99_, final FamilyRoleManifest roles_,
			final String aligner_, final String familyListSha256_, final String consensusRefSha256_,
			final String roleManifestSha256_, final String thresholdsSha256_,
			final String alignerSha256_, final String blosum62Sha256_, final String coveringSetsSha256_,
			final String shortlistSidecarSha256_, final String coveringAlphabet_, final int coveringK_,
			final String kmerBoundary_, final String kmerDefinition_, final String hbmMetric_, final int hbmTargetBp_,
			final String hbmCutoffRule_, final String hbmScoreContract_, final String hbmBundleSha256_,
			final String hbmSemanticProvenanceSha256_, final String hbmLandmarksSha256_,
			final String hbmLandmarksProvenanceSha256_, final String hbmPass4ManifestSha256_,
			final String hbmFamilyHashesSha256_, final int hbmFamilyCount_, final long hbmMemberTotal_,
			final int hbmMaxQueryLength_){
		nFamilies=repIds_.length;
		acceptanceSchema=0;
		if(basePositiveN_.length!=nFamilies || lqLo_.length!=nFamilies || lqHi_.length!=nFamilies || minRawScore_.length!=nFamilies
				|| minId_.length!=nFamilies || minOverlap_.length!=nFamilies || minKmerCount_.length!=nFamilies
				|| minKmerQdensity_.length!=nFamilies || hbmUniqueN_.length!=nFamilies || minHbmRaw99_.length!=nFamilies){
			throw new IllegalArgumentException("FamilyAcceptanceConfig array length mismatch vs nFamilies="+nFamilies);
		}
		if(roles_==null || roles_.size()!=nFamilies){
			throw new IllegalArgumentException("FamilyAcceptanceConfig role-manifest size mismatch vs nFamilies="+nFamilies);
		}
		for(int r=0; r<nFamilies; r++){
			final String ctx="rank "+r+" ('"+repIds_[r]+"')";
			if(repIds_[r]==null || repIds_[r].isEmpty()){throw new IllegalArgumentException("Blank rep_id at "+ctx);}
			if(!repIds_[r].equals(roles_.repId(r))){
				throw new IllegalArgumentException("Role-manifest rep_id mismatch at "+ctx+": '"+roles_.repId(r)+"'.");
			}
			if(basePositiveN_[r]<1){throw new IllegalArgumentException("base_positive_n<1 at "+ctx+": "+basePositiveN_[r]);}
			if(hbmUniqueN_[r]<1){throw new IllegalArgumentException("hbm_unique_n<1 at "+ctx+": "+hbmUniqueN_[r]);}
			if(lqLo_[r]<1){throw new IllegalArgumentException("lq_lo<1 at "+ctx+": "+lqLo_[r]);}
			if(lqHi_[r]<lqLo_[r]){throw new IllegalArgumentException("lq_hi<lq_lo at "+ctx+": lo="+lqLo_[r]+" hi="+lqHi_[r]);}
			if(minRawScore_[r]==Integer.MIN_VALUE || minRawScore_[r]==Integer.MAX_VALUE){
				throw new IllegalArgumentException("min_raw_score at int-sentinel extreme at "+ctx+" (likely corrupt): "+minRawScore_[r]);
			}
			if(!Double.isFinite(minId_[r]) || minId_[r]<0 || minId_[r]>100){throw new IllegalArgumentException("min_identity out of [0,100] or non-finite at "+ctx+": "+minId_[r]);}
			if(!Double.isFinite(minOverlap_[r]) || minOverlap_[r]<0 || minOverlap_[r]>1){throw new IllegalArgumentException("min_overlap out of [0,1] or non-finite at "+ctx+": "+minOverlap_[r]);}
			if(minKmerCount_[r]<0){throw new IllegalArgumentException("min_kmer_count negative at "+ctx+": "+minKmerCount_[r]);}
			if(!Double.isFinite(minKmerQdensity_[r]) || minKmerQdensity_[r]<0 || minKmerQdensity_[r]>1){throw new IllegalArgumentException("min_kmer_qdensity out of [0,1] or non-finite at "+ctx+": "+minKmerQdensity_[r]);}
			if(!Float.isFinite(minHbmRaw99_[r])){throw new IllegalArgumentException("min_hbm_raw_99 non-finite at "+ctx+": "+minHbmRaw99_[r]);}
		}
		repIds=Arrays.copyOf(repIds_, nFamilies); basePositiveN=Arrays.copyOf(basePositiveN_, nFamilies);
		familyIds=new int[nFamilies]; acceptanceEnabled=new boolean[nFamilies];
		minR=new double[nFamilies]; minCoreMutualOverlap=Arrays.copyOf(minOverlap_,nFamilies);
		minHbmPathRelative=new float[nFamilies];
		Arrays.fill(acceptanceEnabled,true); Arrays.fill(minR,Double.NEGATIVE_INFINITY);
		Arrays.fill(minHbmPathRelative,Float.NEGATIVE_INFINITY);
		for(int i=0; i<nFamilies; i++){familyIds[i]=i;}
		hbmUniqueN=Arrays.copyOf(hbmUniqueN_,nFamilies); minHbmRaw99=Arrays.copyOf(minHbmRaw99_,nFamilies);
		lqLo=Arrays.copyOf(lqLo_, nFamilies); lqHi=Arrays.copyOf(lqHi_, nFamilies);
		minRawScore=Arrays.copyOf(minRawScore_, nFamilies);
		minId=Arrays.copyOf(minId_, nFamilies); minOverlap=Arrays.copyOf(minOverlap_, nFamilies);
		minKmerCount=Arrays.copyOf(minKmerCount_, nFamilies); minKmerQdensity=Arrays.copyOf(minKmerQdensity_, nFamilies);
		roles=roles_;
		aligner=aligner_; familyListSha256=familyListSha256_; consensusRefSha256=consensusRefSha256_;
		roleManifestSha256=roleManifestSha256_; thresholdsSha256=thresholdsSha256_;
		profileArtifactSha80="NA";
		coreCoordinatesSha80="NA";
		coreFamilyListSha80=familyListSha256_;
		alignerSha256=alignerSha256_; blosum62Sha256=blosum62Sha256_; coveringSetsSha256=coveringSetsSha256_;
		shortlistSidecarSha256=shortlistSidecarSha256_; coveringAlphabet=coveringAlphabet_; coveringK=coveringK_;
		kmerBoundary=kmerBoundary_; kmerDefinition=kmerDefinition_;
		hbmMetric=hbmMetric_; hbmTargetBp=hbmTargetBp_; hbmCutoffRule=hbmCutoffRule_; hbmScoreContract=hbmScoreContract_;
		hbmBundleSha256=hbmBundleSha256_; hbmSemanticProvenanceSha256=hbmSemanticProvenanceSha256_;
		hbmLandmarksSha256=hbmLandmarksSha256_; hbmLandmarksProvenanceSha256=hbmLandmarksProvenanceSha256_;
		hbmPass4ManifestSha256=hbmPass4ManifestSha256_; hbmFamilyHashesSha256=hbmFamilyHashesSha256_;
		hbmFamilyCount=hbmFamilyCount_; hbmMemberTotal=hbmMemberTotal_; hbmMaxQueryLength=hbmMaxQueryLength_;
	}

	/** Constructs the production config from one fully verified schema-7 table. */
	private FamilyAcceptanceConfig(final Schema7FamilyAcceptanceLoader.Loaded loaded){
		nFamilies=loaded.repIds.length; acceptanceSchema=7;
		if(nFamilies<1 || loaded.familyIds.length!=nFamilies || loaded.n.length!=nFamilies ||
				loaded.enabled.length!=nFamilies || loaded.lengthLo.length!=nFamilies ||
				loaded.lengthHi.length!=nFamilies || loaded.minRaw.length!=nFamilies ||
				loaded.minIdentity.length!=nFamilies || loaded.minR.length!=nFamilies ||
				loaded.minMutual.length!=nFamilies || loaded.minKmer.length!=nFamilies ||
				loaded.minKmerDensity.length!=nFamilies || loaded.minHbmPath.length!=nFamilies ||
				loaded.roles==null || loaded.roles.size()!=nFamilies){
			throw new IllegalArgumentException("Schema-7 FamilyAcceptanceConfig array length mismatch.");
		}
		for(int rank=0; rank<nFamilies; rank++){
			final String rep=loaded.repIds[rank];
			if(rep==null || rep.length()==0 || !rep.equals(loaded.roles.repId(rank)) ||
					loaded.familyIds[rank]<0 || loaded.n[rank]<1 || loaded.lengthLo[rank]<1 ||
					loaded.lengthHi[rank]<loaded.lengthLo[rank] || loaded.minKmer[rank]<0 ||
					!Double.isFinite(loaded.minIdentity[rank]) || loaded.minIdentity[rank]<0 ||
					loaded.minIdentity[rank]>100 || !Double.isFinite(loaded.minR[rank]) ||
					!Double.isFinite(loaded.minMutual[rank]) || loaded.minMutual[rank]<0 ||
					loaded.minMutual[rank]>1 || !Double.isFinite(loaded.minKmerDensity[rank]) ||
					loaded.minKmerDensity[rank]<0 || loaded.minKmerDensity[rank]>1 ||
					!Float.isFinite(loaded.minHbmPath[rank]) || loaded.minHbmPath[rank]<0 ||
					loaded.minHbmPath[rank]>1){
				throw new IllegalArgumentException("Invalid schema-7 threshold values at rank "+rank+".");
			}
		}
		repIds=Arrays.copyOf(loaded.repIds,nFamilies);
		familyIds=Arrays.copyOf(loaded.familyIds,nFamilies);
		acceptanceEnabled=Arrays.copyOf(loaded.enabled,nFamilies);
		basePositiveN=Arrays.copyOf(loaded.n,nFamilies);
		hbmUniqueN=Arrays.copyOf(loaded.n,nFamilies);
		lqLo=Arrays.copyOf(loaded.lengthLo,nFamilies);
		lqHi=Arrays.copyOf(loaded.lengthHi,nFamilies);
		minRawScore=Arrays.copyOf(loaded.minRaw,nFamilies);
		minId=Arrays.copyOf(loaded.minIdentity,nFamilies);
		minR=Arrays.copyOf(loaded.minR,nFamilies);
		minOverlap=Arrays.copyOf(loaded.minMutual,nFamilies);
		minCoreMutualOverlap=Arrays.copyOf(loaded.minMutual,nFamilies);
		minKmerCount=Arrays.copyOf(loaded.minKmer,nFamilies);
		minKmerQdensity=Arrays.copyOf(loaded.minKmerDensity,nFamilies);
		minHbmPathRelative=Arrays.copyOf(loaded.minHbmPath,nFamilies);
		minHbmRaw99=Arrays.copyOf(loaded.minHbmPath,nFamilies);
		roles=loaded.roles;
		aligner="D55"; familyListSha256=loaded.rosterSha80;
		consensusRefSha256=loaded.consensusSha80;
		roleManifestSha256=loaded.roleManifestSha80;
		thresholdsSha256=loaded.tableSha80;
		profileArtifactSha80=loaded.artifactSha80;
		coreCoordinatesSha80=loaded.coreSha80;
		coreFamilyListSha80=loaded.coreFamilyListSha80;
		alignerSha256="BOUND_BY_HBM_PROVENANCE";
		blosum62Sha256="BOUND_BY_HBM_PROVENANCE";
		coveringSetsSha256=loaded.coveringSetsSha80;
		shortlistSidecarSha256=loaded.sidecarSha80;
		coveringAlphabet=loaded.coveringAlphabet; coveringK=loaded.coveringK;
		kmerBoundary="NONE"; kmerDefinition="SIDECAR_BOUND";
		hbmMetric="hbm_path_relative"; hbmTargetBp=9900;
		hbmCutoffRule="schema7_effective_inclusive";
		hbmScoreContract=HbmScoreContract.VALUE;
		hbmBundleSha256=loaded.hbmBundleSha80;
		hbmSemanticProvenanceSha256=loaded.hbmProvenanceSha80;
		hbmLandmarksSha256="NA"; hbmLandmarksProvenanceSha256="NA";
		hbmPass4ManifestSha256="NA"; hbmFamilyHashesSha256="NA";
		hbmFamilyCount=nFamilies; hbmMemberTotal=loaded.members; hbmMaxQueryLength=0;
	}

	private void checkRank(final int rank){
		if(rank<0 || rank>=nFamilies){throw new IndexOutOfBoundsException("rank "+rank+" out of [0,"+nFamilies+")");}
	}
	public String repId(final int rank){checkRank(rank); return repIds[rank];}
	public int familyId(final int rank){checkRank(rank); return familyIds[rank];}
	public boolean acceptanceEnabled(final int rank){checkRank(rank); return acceptanceEnabled[rank];}
	public boolean isSchema7(){return acceptanceSchema==7;}
	/** Historical accessor retained for callers; schema 4 names this population explicitly. */
	public int n(final int rank){return basePositiveN(rank);}
	public int basePositiveN(final int rank){checkRank(rank); return basePositiveN[rank];}
	public int hbmUniqueN(final int rank){checkRank(rank); return hbmUniqueN[rank];}
	public int lqLo(final int rank){checkRank(rank); return lqLo[rank];}
	public int lqHi(final int rank){checkRank(rank); return lqHi[rank];}
	public int minRawScore(final int rank){checkRank(rank); return minRawScore[rank];}
	public double minId(final int rank){checkRank(rank); return minId[rank];}
	public double minR(final int rank){checkRank(rank); return minR[rank];}
	public double minOverlap(final int rank){checkRank(rank); return minOverlap[rank];}
	public double minCoreMutualOverlap(final int rank){checkRank(rank); return minCoreMutualOverlap[rank];}
	public int minKmerCount(final int rank){checkRank(rank); return minKmerCount[rank];}
	public double minKmerQdensity(final int rank){checkRank(rank); return minKmerQdensity[rank];}
	public FamilyRoleManifest.Role role(final int rank){checkRank(rank); return roles.role(rank);}
	public boolean isTracked(final int rank){checkRank(rank); return roles.isTracked(rank);}
	public boolean isBait(final int rank){checkRank(rank); return roles.isBait(rank);}
	public boolean isNetworkFeature(final int rank){checkRank(rank); return roles.isNetworkFeature(rank);}
	public float minHbmRaw99(final int rank){checkRank(rank); return minHbmRaw99[rank];}
	public float minHbmPathRelative(final int rank){checkRank(rank); return minHbmPathRelative[rank];}
	/** Compatibility name for the sole strict schema-4 HBM cutoff. */
	public double hbmMin(final int rank){return minHbmRaw99(rank);}

	/**
	 * Loads the canonical schema-7 production profile and binds every live
	 * assignment resource named by it.
	 */
	public static FamilyAcceptanceConfig loadSchema7(final String artifactPath,
			final String expectedArtifactSha80,
			final String rosterPath, final String consensusPath,
			final String roleManifestPath, final String coreCoordinatesPath,
			final String coveringSetsPath, final String sidecarPath,
			final String hbmBundlePath, final String hbmProvenancePath){
		return new FamilyAcceptanceConfig(Schema7FamilyAcceptanceLoader.load(
			artifactPath,expectedArtifactSha80,rosterPath,consensusPath,roleManifestPath,
			coreCoordinatesPath,coveringSetsPath,sidecarPath,hbmBundlePath,
			hbmProvenancePath));
	}

	/**
	 * Loads and fully verifies a {@code d55_d107_hbm_family_thresholds} artifact.
	 * @param artifactPath The exporter's output file.
	 * @param liveFamilyListPath The LIVE family-list/roster file ({@code #rank\trep_id\tocc_total}
	 *        schema) -- re-hashed against the artifact's {@code #family_list_sha256} and used to
	 *        derive the authoritative rank order every data row is checked against.
	 * @param liveConsensusRefPath The LIVE consensus reference FASTA -- re-hashed against {@code
	 *        #consensus_ref_sha256}.
	 * @param liveRoleManifestPath The LIVE per-rank role manifest -- re-hashed against
	 *        {@code #role_manifest_sha256} and loaded against the roster's exact rep-ID order.
	 * @param expectedAligner The aligner the caller is about to run hits through; must equal both
	 *        this and the artifact's own recorded {@code #aligner} (which must itself be {@code
	 *        D55} -- this artifact type is not a generic multi-aligner container).
	 * @param expectedAlignerSha256 Caller-computed hash of the live D55 aligner (this loader has
	 *        no independent way to re-derive it from inside a running JVM).
	 * @param expectedBlosum62Sha256 Caller-computed hash of the live BLOSUM62 source.
	 * @param liveCoveringSetsPath The LIVE covering-set k-mer file the sidecar was built from --
	 *        re-hashed against both the sidecar's own {@code kmersets_sha256} and the artifact's
	 *        {@code #covering_sets_sha256}.
	 * @param liveSidecarPath The LIVE {@code FamilyShortlistSidecar} TSV -- its header is read
	 *        (not the full composition matrix) for {@code consensus_sha256}/{@code kmersets_sha256}/
	 *        {@code covering_alphabet}/{@code covering_k}, and the whole file is re-hashed against
	 *        the artifact's {@code #shortlist_sidecar_sha256}.
	 * @param expectedKmerDefinition Caller-supplied expected {@code #kmer_definition} value (see
	 *        class javadoc -- no live file represents this, so it cannot be self-verified).
	 * @return The loaded, verified config.
	 */
	public static FamilyAcceptanceConfig load(final String artifactPath, final String liveFamilyListPath,
			final String liveConsensusRefPath, final String liveRoleManifestPath, final String expectedAligner, final String expectedAlignerSha256,
			final String expectedBlosum62Sha256, final String liveCoveringSetsPath, final String liveSidecarPath,
			final String expectedKmerDefinition, final String liveHbmBundlePath, final String liveHbmProvenancePath){
		return loadInternal(artifactPath,liveFamilyListPath,liveConsensusRefPath,liveRoleManifestPath,expectedAligner,
			expectedAlignerSha256,expectedBlosum62Sha256,liveCoveringSetsPath,liveSidecarPath,expectedKmerDefinition,
			liveHbmBundlePath,liveHbmProvenancePath,null);
	}

	/** Package-private deterministic trust-window fixture seam; production callers use {@link #load}. */
	static FamilyAcceptanceConfig loadWithBeforeFinalRehashForTest(final String artifactPath, final String liveFamilyListPath,
			final String liveConsensusRefPath, final String liveRoleManifestPath, final String expectedAligner, final String expectedAlignerSha256,
			final String expectedBlosum62Sha256, final String liveCoveringSetsPath, final String liveSidecarPath,
			final String expectedKmerDefinition, final String liveHbmBundlePath, final String liveHbmProvenancePath,
			final Runnable beforeFinalRehash){
		if(beforeFinalRehash==null){throw new IllegalArgumentException("Test trust-window hook must not be null.");}
		return loadInternal(artifactPath,liveFamilyListPath,liveConsensusRefPath,liveRoleManifestPath,expectedAligner,
			expectedAlignerSha256,expectedBlosum62Sha256,liveCoveringSetsPath,liveSidecarPath,expectedKmerDefinition,
			liveHbmBundlePath,liveHbmProvenancePath,beforeFinalRehash);
	}

	private static FamilyAcceptanceConfig loadInternal(final String artifactPath, final String liveFamilyListPath,
			final String liveConsensusRefPath, final String liveRoleManifestPath, final String expectedAligner, final String expectedAlignerSha256,
			final String expectedBlosum62Sha256, final String liveCoveringSetsPath, final String liveSidecarPath,
			final String expectedKmerDefinition, final String liveHbmBundlePath, final String liveHbmProvenancePath,
			final Runnable beforeFinalRehash){
		DigestSuffix.requireSuffix(expectedAlignerSha256,"expected aligner sha80");
		DigestSuffix.requireSuffix(expectedBlosum62Sha256,"expected BLOSUM62 sha80");
		final String liveFamilyListSha256=sha80File(liveFamilyListPath);
		final String[] rosterRepIds=readRosterOrdered(liveFamilyListPath);
		requireUnchanged(liveFamilyListPath,liveFamilyListSha256,"family list");
		final String liveConsensusSha256=sha80File(liveConsensusRefPath);
		final String liveRoleManifestSha256=sha80File(liveRoleManifestPath);
		final String liveHbmBundleSha256=sha80File(liveHbmBundlePath);
		final String liveHbmProvenanceSha256=sha80File(liveHbmProvenancePath);
		final FamilyRoleManifest roles=FamilyRoleManifest.load(liveRoleManifestPath,rosterRepIds);
		requireUnchanged(liveRoleManifestPath,liveRoleManifestSha256,"role manifest");
		final byte[][] hbmProvenance=HbmBundleLoader.loadSemanticProvenance(liveHbmProvenancePath);
		requireUnchanged(liveHbmProvenancePath,liveHbmProvenanceSha256,"HBM semantic provenance");
		if(!liveFamilyListSha256.equals(DigestSuffix.fromDigest(hbmProvenance[10]))){
			throw new IllegalArgumentException("HBM semantic provenance roster hash does not match the live family list: "+liveHbmProvenancePath);
		}
		if(!liveConsensusSha256.equals(DigestSuffix.fromDigest(hbmProvenance[11]))){
			throw new IllegalArgumentException("HBM semantic provenance consensus hash does not match the live consensus: "+liveHbmProvenancePath);
		}

		final String liveSidecarSha256=sha80File(liveSidecarPath);
		final SidecarHeaderFields sidecarHeader=readSidecarHeaderFields(liveSidecarPath);
		requireUnchanged(liveSidecarPath,liveSidecarSha256,"shortlist sidecar");
		if(!liveConsensusSha256.equals(sidecarHeader.consensusSha256)){
			throw new IllegalArgumentException("Live sidecar "+liveSidecarPath+"'s own consensus_sha256='"+
				sidecarHeader.consensusSha256+"' != live consensus hash "+liveConsensusSha256+
				" -- the sidecar and the live consensus reference are not the same build.");
		}
		final String liveCoveringSetsSha256=sha80File(liveCoveringSetsPath);
		if(!liveCoveringSetsSha256.equals(sidecarHeader.kmerSetsSha256)){
			throw new IllegalArgumentException("Live coveringsets file "+liveCoveringSetsPath+" hashes to "+
				liveCoveringSetsSha256+", but live sidecar "+liveSidecarPath+"'s own kmersets_sha256='"+
				sidecarHeader.kmerSetsSha256+"' -- mismatched covering-set file.");
		}
		final HashMap<String,String> header=new HashMap<String,String>();
		String columnHeaderLine=null;
		final ByteBuilder dataBytes=new ByteBuilder();
		final List<String[]> rows=new ArrayList<String[]>();
		final ByteFile bf=ByteFile.makeByteFile(artifactPath, false);
		try{
			for(byte[] lineBytes=bf.nextLine(); lineBytes!=null; lineBytes=bf.nextLine()){
				if(lineBytes.length==0){continue;}
				final String line=new String(lineBytes, StandardCharsets.UTF_8);
				if(line.charAt(0)=='#'){
					final int tab=line.indexOf('\t');
					if(tab<0){throw new IllegalArgumentException("Malformed header line (no tab) in "+artifactPath+": "+line);}
					final String hKey=line.substring(1, tab);
					String hValue=line.substring(tab+1);
					if(hKey.endsWith("_sha256") || hKey.endsWith("_sha80")){
						if(hKey.endsWith("_sha256") && hValue.length()!=64){throw new IllegalArgumentException("Legacy digest header must contain 64 lowercase hex characters: "+hKey);}
						if(hKey.endsWith("_sha80") && hValue.length()!=DigestSuffix.HEX_LENGTH){throw new IllegalArgumentException("sha80 header must contain 20 lowercase hex characters: "+hKey);}
						hValue=DigestSuffix.normalizeRecorded(hValue,hKey);
					}
					if(header.put(hKey, hValue)!=null){
						throw new IllegalArgumentException("Duplicate header key '"+hKey+"' in "+artifactPath);
					}
					continue;
				}
				if(columnHeaderLine==null){
					columnHeaderLine=line;
					if(!columnHeaderLine.equals("rank\trep_id\tbase_positive_n\tlq_lo\tlq_hi\tmin_raw_score\tmin_identity\tmin_overlap\tmin_kmer_count\tmin_kmer_qdensity\thbm_unique_n\tmin_hbm_raw_99")){
						throw new IllegalArgumentException("Unexpected column header in "+artifactPath+": "+columnHeaderLine);
					}
					continue;
				}
				dataBytes.append(line).append('\n');
				final String[] cols=line.split("\t", -1);
				if(cols.length!=12){throw new IllegalArgumentException("Data row width!=12 in "+artifactPath+": "+line);}
				rows.add(cols);
			}
		} finally { bf.close(); }

		final String schemaVersion=required(header,"schema_version",artifactPath);
		if(!("4".equals(schemaVersion) || "5".equals(schemaVersion) || "6".equals(schemaVersion))){
			throw new IllegalArgumentException("Unsupported schema_version in "+artifactPath+" (this loader requires "+
				"4 legacy, 5 sha80, or 6 purpose-bound sha80, all strict HBM-bound schemas; older HBM-off artifacts are refused): "+
				header.get("schema_version"));
		}
		if(!"d55_d107_hbm_family_thresholds".equals(required(header, "artifact_type", artifactPath))){
			throw new IllegalArgumentException("Unexpected artifact_type in "+artifactPath+": "+header.get("artifact_type"));
		}
		if("6".equals(schemaVersion) && !"consensus_refinement".equals(required(header,"profile_purpose",artifactPath))){
			throw new IllegalArgumentException("Schema 6 profile_purpose must be consensus_refinement: "+artifactPath);
		}
		final String aligner=required(header, "aligner", artifactPath);
		if(!"D55".equals(aligner)){
			throw new IllegalArgumentException("Artifact aligner is not D55 (this artifact type is D55-specific): "+aligner);
		}
		if(!aligner.equals(expectedAligner)){
			throw new IllegalArgumentException("Artifact aligner='"+aligner+"' != caller's expectedAligner='"+expectedAligner+"'.");
		}
		final String alignerSha256=requiredDigest(header, "aligner", artifactPath);
		if(!alignerSha256.equals(expectedAlignerSha256)){
			throw new IllegalArgumentException("Artifact aligner_sha256 mismatch: recorded="+alignerSha256+" live="+expectedAlignerSha256);
		}
		final String blosum62Sha256=requiredDigest(header, "blosum62", artifactPath);
		if(!blosum62Sha256.equals(expectedBlosum62Sha256)){
			throw new IllegalArgumentException("Artifact blosum62_sha256 mismatch: recorded="+blosum62Sha256+" live="+expectedBlosum62Sha256);
		}
		final String familyListSha256=requiredDigest(header, "family_list", artifactPath);
		if(!familyListSha256.equals(liveFamilyListSha256)){
			throw new IllegalArgumentException("Artifact family_list_sha256 mismatch: recorded="+familyListSha256+
				" but the live family list "+liveFamilyListPath+" hashes to "+liveFamilyListSha256+
				" -- the roster has changed since this artifact was built; rebuild it before using it.");
		}
		final String consensusRefSha256=requiredDigest(header, "consensus_ref", artifactPath);
		if(!consensusRefSha256.equals(liveConsensusSha256)){
			throw new IllegalArgumentException("Artifact consensus_ref_sha256 mismatch: recorded="+consensusRefSha256+
				" but the live consensus "+liveConsensusRefPath+" hashes to "+liveConsensusSha256+
				" -- rebuild the artifact before using it.");
		}
		final String roleManifestSha256=requiredDigest(header, "role_manifest", artifactPath);
		if(!roleManifestSha256.equals(liveRoleManifestSha256)){
			throw new IllegalArgumentException("Artifact role_manifest_sha256 mismatch: recorded="+roleManifestSha256+
				" but the live role manifest "+liveRoleManifestPath+" hashes to "+liveRoleManifestSha256+
				" -- role assignments changed since this artifact was built; rebuild it before use.");
		}
		final int nFamiliesHeader=Integer.parseInt(required(header, "n_families", artifactPath));
		if(nFamiliesHeader!=rosterRepIds.length){
			throw new IllegalArgumentException("Artifact n_families="+nFamiliesHeader+" != live family list size "+rosterRepIds.length);
		}
		final String coveringSetsSha256=requiredDigest(header, "covering_sets", artifactPath);
		if(!coveringSetsSha256.equals(liveCoveringSetsSha256)){
			throw new IllegalArgumentException("Artifact covering_sets_sha256 mismatch: recorded="+coveringSetsSha256+
				" but the live covering-set file "+liveCoveringSetsPath+" hashes to "+liveCoveringSetsSha256+
				" -- rebuild the artifact before using it.");
		}
		final String shortlistSidecarSha256=requiredDigest(header, "shortlist_sidecar", artifactPath);
		if(!shortlistSidecarSha256.equals(liveSidecarSha256)){
			throw new IllegalArgumentException("Artifact shortlist_sidecar_sha256 mismatch: recorded="+shortlistSidecarSha256+
				" but the live sidecar "+liveSidecarPath+" hashes to "+liveSidecarSha256+
				" -- the sidecar has changed since this artifact was built; rebuild it before using it.");
		}
		final String coveringAlphabet=required(header, "covering_alphabet", artifactPath);
		if(!coveringAlphabet.equals(sidecarHeader.coveringAlphabet)){
			throw new IllegalArgumentException("Artifact covering_alphabet='"+coveringAlphabet+
				"' != live sidecar's own covering_alphabet='"+sidecarHeader.coveringAlphabet+"'.");
		}
		final int coveringK=Integer.parseInt(required(header, "covering_k", artifactPath));
		if(coveringK<=0){throw new IllegalArgumentException("Artifact covering_k must be positive, got: "+coveringK);}
		if(coveringK!=sidecarHeader.coveringK){
			throw new IllegalArgumentException("Artifact covering_k="+coveringK+
				" != live sidecar's own covering_k="+sidecarHeader.coveringK+".");
		}
		final String kmerStatGating=required(header, "kmer_stat_gating", artifactPath);
		if(!"s_k_AND_s_k_over_Qsize".equals(kmerStatGating)){
			throw new IllegalArgumentException("Artifact kmer_stat_gating must be exactly s_k_AND_s_k_over_Qsize "+
				"(D105: count AND query-density, combined by AND), got: "+kmerStatGating);
		}
		final String kmerBoundary=required(header, "kmer_boundary", artifactPath);
		if(!"NONE".equals(kmerBoundary)){
			throw new IllegalArgumentException("Artifact kmer_boundary must be NONE (the only boundary mode "+
				"production code uses), got: "+kmerBoundary);
		}
		final String kmerDefinition=required(header, "kmer_definition", artifactPath);
		if(!kmerDefinition.equals(expectedKmerDefinition)){
			throw new IllegalArgumentException("Artifact kmer_definition='"+kmerDefinition+
				"' != caller's expectedKmerDefinition='"+expectedKmerDefinition+"'.");
		}
		final String hbmEnabledStr=required(header, "hbm_enabled", artifactPath);
		if(!"1".equals(hbmEnabledStr)){throw new IllegalArgumentException("Strict schema 4 requires hbm_enabled=1, got '"+hbmEnabledStr+"': "+artifactPath);}
		final String hbmMetric=required(header,"hbm_metric",artifactPath);
		if(!"hbm_raw".equals(hbmMetric)){throw new IllegalArgumentException("Strict schema 4 requires hbm_metric=hbm_raw, got: "+hbmMetric);}
		final int hbmTargetBp=parsePositiveHeaderInt(header,"hbm_target_bp",artifactPath);
		if(hbmTargetBp!=9900){throw new IllegalArgumentException("Strict schema 4 requires hbm_target_bp=9900, got: "+hbmTargetBp);}
		final String hbmCutoffRule=required(header,"hbm_cutoff_rule",artifactPath);
		if(!"highest inclusive cutoff retaining at least ceil(hbm_unique_n*target) positives".equals(hbmCutoffRule)){
			throw new IllegalArgumentException("Unexpected hbm_cutoff_rule: "+hbmCutoffRule);
		}
		final String hbmScoreContract=required(header,"hbm_score_contract",artifactPath);
		if(!"HbmBundleLoader.Loaded.score_raw_index0_v1".equals(hbmScoreContract)){
			throw new IllegalArgumentException("Unexpected hbm_score_contract: "+hbmScoreContract);
		}
		final String hbmBundleSha256=requiredHashHeader(header,"hbm_bundle_sha256",artifactPath);
		if(!hbmBundleSha256.equals(liveHbmBundleSha256)){throw new IllegalArgumentException("HBM bundle hash mismatch: recorded="+hbmBundleSha256+" live="+liveHbmBundleSha256);}
		final String hbmSemanticProvenanceSha256=requiredHashHeader(header,"hbm_semantic_provenance_sha256",artifactPath);
		if(!hbmSemanticProvenanceSha256.equals(liveHbmProvenanceSha256)){throw new IllegalArgumentException("HBM semantic provenance hash mismatch: recorded="+hbmSemanticProvenanceSha256+" live="+liveHbmProvenanceSha256);}
		final String hbmLandmarksSha256=requiredHashHeader(header,"hbm_landmarks_sha256",artifactPath);
		final String hbmLandmarksProvenanceSha256=requiredHashHeader(header,"hbm_landmarks_provenance_sha256",artifactPath);
		final String hbmPass4ManifestSha256=requiredHashHeader(header,"hbm_pass4_manifest_sha256",artifactPath);
		final String hbmFamilyHashesSha256=requiredHashHeader(header,"hbm_family_hashes_sha256",artifactPath);
		final int hbmFamilyCount=parsePositiveHeaderInt(header,"hbm_family_count",artifactPath);
		if(hbmFamilyCount!=rosterRepIds.length){throw new IllegalArgumentException("hbm_family_count="+hbmFamilyCount+" != roster count "+rosterRepIds.length);}
		final long hbmMemberTotal=parsePositiveHeaderLong(header,"hbm_member_total",artifactPath);
		final int hbmMaxQueryLength=parsePositiveHeaderInt(header,"hbm_max_query_length",artifactPath);
		if("4".equals(schemaVersion) || "5".equals(schemaVersion)){
			requiredHashHeader(header,"base_profile_sha256",artifactPath);
		}else{
			//Schema 6 is assembled directly from the independently accepted current-family
			//base/length tables.  Requiring a fictional intermediate "base profile" would
			//weaken provenance by naming an artifact that never existed.
			requiredDigest(header,"base_landmarks",artifactPath);
			requiredDigest(header,"length_bounds",artifactPath);
			requiredDigest(header,"base_landmarks_provenance",artifactPath);
		}
		final String recordedThresholdsSha256=requiredDigest(header, "thresholds", artifactPath);
		final String liveThresholdsSha256=DigestSuffix.bytes(dataBytes.toBytes());
		if(!recordedThresholdsSha256.equals(liveThresholdsSha256)){
			throw new IllegalArgumentException("Artifact thresholds_sha256 mismatch: recorded="+recordedThresholdsSha256+
				" but the actual data rows hash to "+liveThresholdsSha256+" -- the artifact was hand-edited or "+
				"corrupted after export. See "+artifactPath);
		}

		if(rows.size()!=rosterRepIds.length){
			throw new IllegalArgumentException("Artifact has "+rows.size()+" data rows but the live family list has "+
				rosterRepIds.length+" families: "+artifactPath);
		}
		final int nFam=rosterRepIds.length;
		final String[] repIds=new String[nFam]; final int[] basePositiveN=new int[nFam], lqLo=new int[nFam], lqHi=new int[nFam],
			minRawScore=new int[nFam], minKmerCount=new int[nFam], hbmUniqueN=new int[nFam];
		final double[] minId=new double[nFam], minOverlap=new double[nFam], minKmerQdensity=new double[nFam];
		final float[] minHbmRaw99=new float[nFam]; long observedHbmMembers=0;
		for(int i=0; i<nFam; i++){
			final String[] cols=rows.get(i);
			final int rank=parseIntField(cols[0], "rank", i, artifactPath);
			if(rank!=i){throw new IllegalArgumentException("Artifact row "+i+" has rank="+rank+", expected "+i+
				" (rows must be in the live family list's rank order): "+artifactPath);}
			final String repId=cols[1];
			if(!repId.equals(rosterRepIds[i])){
				throw new IllegalArgumentException("Artifact row "+i+" rep_id='"+repId+"' != live family list's rep_id '"+
					rosterRepIds[i]+"' at rank "+i+" -- artifact row order does not match the live roster: "+artifactPath);
			}
			repIds[i]=repId;
			basePositiveN[i]=parseIntField(cols[2], "base_positive_n", i, artifactPath);
			lqLo[i]=parseIntField(cols[3], "lq_lo", i, artifactPath);
			lqHi[i]=parseIntField(cols[4], "lq_hi", i, artifactPath);
			minRawScore[i]=parseIntField(cols[5], "min_raw_score", i, artifactPath);
			minId[i]=parseDoubleField(cols[6], "min_identity", i, artifactPath);
			minOverlap[i]=parseDoubleField(cols[7], "min_overlap", i, artifactPath);
			minKmerCount[i]=parseIntField(cols[8], "min_kmer_count", i, artifactPath);
			minKmerQdensity[i]=parseDoubleField(cols[9], "min_kmer_qdensity", i, artifactPath);
			hbmUniqueN[i]=parseIntField(cols[10], "hbm_unique_n", i, artifactPath);
			minHbmRaw99[i]=parseFloatField(cols[11], "min_hbm_raw_99", i, artifactPath);
			observedHbmMembers=Math.addExact(observedHbmMembers,hbmUniqueN[i]);
		}
		if(observedHbmMembers!=hbmMemberTotal){throw new IllegalArgumentException("Sum hbm_unique_n="+observedHbmMembers+" != hbm_member_total="+hbmMemberTotal);}
		if(beforeFinalRehash!=null){beforeFinalRehash.run();}
		requireUnchanged(liveFamilyListPath,liveFamilyListSha256,"family list");
		requireUnchanged(liveConsensusRefPath,liveConsensusSha256,"consensus reference");
		requireUnchanged(liveRoleManifestPath,liveRoleManifestSha256,"role manifest");
		requireUnchanged(liveCoveringSetsPath,liveCoveringSetsSha256,"covering sets");
		requireUnchanged(liveSidecarPath,liveSidecarSha256,"shortlist sidecar");
		requireUnchanged(liveHbmBundlePath,liveHbmBundleSha256,"HBM bundle");
		requireUnchanged(liveHbmProvenancePath,liveHbmProvenanceSha256,"HBM semantic provenance");
		return new FamilyAcceptanceConfig(repIds, basePositiveN, lqLo, lqHi, minRawScore, minId, minOverlap, minKmerCount,
			minKmerQdensity, hbmUniqueN, minHbmRaw99, roles, aligner, familyListSha256, consensusRefSha256, roleManifestSha256, recordedThresholdsSha256,
			alignerSha256, blosum62Sha256, coveringSetsSha256, shortlistSidecarSha256, coveringAlphabet, coveringK,
			kmerBoundary, kmerDefinition, hbmMetric, hbmTargetBp, hbmCutoffRule, hbmScoreContract, hbmBundleSha256,
			hbmSemanticProvenanceSha256, hbmLandmarksSha256, hbmLandmarksProvenanceSha256, hbmPass4ManifestSha256,
			hbmFamilyHashesSha256, hbmFamilyCount, hbmMemberTotal, hbmMaxQueryLength);
	}

	private static int parseIntField(final String s, final String field, final int rowIdx, final String path){
		try{return Integer.parseInt(s);}
		catch(NumberFormatException e){throw new IllegalArgumentException("Non-integer "+field+" at row "+rowIdx+" in "+path+": '"+s+"'.");}
	}
	private static double parseDoubleField(final String s, final String field, final int rowIdx, final String path){
		final double v;
		try{v=Double.parseDouble(s);}
		catch(NumberFormatException e){throw new IllegalArgumentException("Non-numeric "+field+" at row "+rowIdx+" in "+path+": '"+s+"'.");}
		if(!Double.isFinite(v)){throw new IllegalArgumentException("Non-finite "+field+" at row "+rowIdx+" in "+path+": "+v);}
		return v;
	}
	private static float parseFloatField(final String s, final String field, final int rowIdx, final String path){
		final float v;
		try{v=Float.parseFloat(s);}
		catch(NumberFormatException e){throw new IllegalArgumentException("Non-numeric "+field+" at row "+rowIdx+" in "+path+": '"+s+"'.");}
		if(!Float.isFinite(v) || !Float.toString(v).equals(s)){
			throw new IllegalArgumentException("Non-finite or noncanonical float "+field+" at row "+rowIdx+" in "+path+": '"+s+"'.");
		}
		return v;
	}
	private static String requiredHashHeader(final HashMap<String,String> header, final String key, final String path){
		final String suffix="_sha256";
		final String base=key.endsWith(suffix) ? key.substring(0,key.length()-suffix.length()) : key;
		return requiredDigest(header,base,path);
	}
	private static String requiredDigest(final HashMap<String,String> header, final String base, final String path){
		return DigestSuffix.requiredProfileHeader(header,base,path);
	}
	private static int parsePositiveHeaderInt(final HashMap<String,String> header, final String key, final String path){
		final long value=parsePositiveHeaderLong(header,key,path);
		if(value>Integer.MAX_VALUE){throw new IllegalArgumentException("Header "+key+" overflows int: "+value+" in "+path);}
		return (int)value;
	}
	private static long parsePositiveHeaderLong(final HashMap<String,String> header, final String key, final String path){
		final long value=HbmCanonicalParse.parseCanonicalNonNegLong(required(header,key,path),key);
		if(value<1){throw new IllegalArgumentException("Header "+key+" must be positive in "+path);}
		return value;
	}

	/** Reads the live family-list/roster ({@code #rank\trep_id\tocc_total}), returning rep_ids
	 *  ordered by rank 0..N-1 contiguous (fail-closed on gap/dup rank or dup rep_id). */
	private static String[] readRosterOrdered(final String path){
		final HashMap<Integer,String> byRank=new HashMap<Integer,String>();
		final java.util.HashSet<String> seenRepIds=new java.util.HashSet<String>();
		final ByteFile file=ByteFile.makeByteFile(path, false);
		try{
			final byte[] header=file.nextLine();
			if(header==null || !new String(header, StandardCharsets.UTF_8).equals("#rank\trep_id\tocc_total")){
				throw new IllegalArgumentException("Unexpected family-list header: "+path);
			}
			byte[] line;
			while((line=file.nextLine())!=null){
				if(line.length==0){continue;}
				final String[] f=new String(line, StandardCharsets.UTF_8).split("\t", -1);
				if(f.length!=3){throw new IllegalArgumentException("Family-list row width!=3: "+new String(line, StandardCharsets.UTF_8));}
				final int rank=Integer.parseInt(f[0]);
				final String repId=f[1];
				if(repId.isEmpty()){throw new IllegalArgumentException("Blank rep_id in family list");}
				if(rank<0){throw new IllegalArgumentException("Negative rank in family list: "+rank);}
				if(byRank.put(rank, repId)!=null){throw new IllegalArgumentException("Duplicate family-list rank: "+rank);}
				if(!seenRepIds.add(repId)){throw new IllegalArgumentException("Duplicate family-list rep_id: "+repId);}
			}
		} finally { file.close(); }
		if(byRank.isEmpty()){throw new IllegalArgumentException("Empty family list: "+path);}
		final String[] out=new String[byRank.size()];
		for(int r=0; r<out.length; r++){
			final String repId=byRank.get(r);
			if(repId==null){throw new IllegalArgumentException("Family-list rank sequence has a gap at rank "+r);}
			out[r]=repId;
		}
		return out;
	}

	static final class SidecarHeaderFields{
		final String consensusSha256, kmerSetsSha256, coveringAlphabet; final int coveringK;
		SidecarHeaderFields(String c,String k,String a,int ck){consensusSha256=c;kmerSetsSha256=k;coveringAlphabet=a;coveringK=ck;}
	}
	/** Reads ONLY the {@code #key\tvalue} header block of a {@code FamilyShortlistSidecar} TSV
	 *  (stopping at the first non-header line) -- a deliberately cheap partial parse, avoiding the
	 *  full composition-matrix load {@link FamilyShortlistSidecar#load} would do, since only 4
	 *  scalar header fields are needed here. */
	static SidecarHeaderFields readSidecarHeaderFields(final String path){
		final HashMap<String,String> header=new HashMap<String,String>();
		final ByteFile file=ByteFile.makeByteFile(path, false);
		try{
			byte[] line;
			while((line=file.nextLine())!=null){
				if(line.length==0){continue;}
				if(line[0]!='#'){break;}
				final String s=new String(line, StandardCharsets.UTF_8);
				final int tab=s.indexOf('\t');
				if(tab<0){continue;}
				final String key=s.substring(1,tab);
				String value=s.substring(tab+1);
				if(key.endsWith("_sha256") || key.endsWith("_sha80")){
					if(key.endsWith("_sha256") && value.length()!=64){throw new IllegalArgumentException("Legacy sidecar digest must contain 64 lowercase hex characters: "+key);}
					if(key.endsWith("_sha80") && value.length()!=DigestSuffix.HEX_LENGTH){throw new IllegalArgumentException("Sidecar sha80 digest must contain 20 lowercase hex characters: "+key);}
					value=DigestSuffix.normalizeRecorded(value,key);
				}
				header.put(key,value);
			}
		} finally { file.close(); }
		final String consensusSha256=DigestSuffix.requiredCompatibleHeader(header, "consensus", path);
		final String kmerSetsSha256=DigestSuffix.requiredCompatibleHeader(header, "kmersets", path);
		final String coveringAlphabet=required(header, "covering_alphabet", path);
		final int coveringK;
		try{coveringK=Integer.parseInt(required(header, "covering_k", path));}
		catch(NumberFormatException e){throw new IllegalArgumentException("Sidecar covering_k is not an integer: "+path);}
		if(coveringK<=0){throw new IllegalArgumentException("Sidecar covering_k must be positive, got "+coveringK+": "+path);}
		return new SidecarHeaderFields(consensusSha256, kmerSetsSha256, coveringAlphabet, coveringK);
	}

	private static String required(final HashMap<String,String> header, final String key, final String path){
		final String v=header.get(key);
		if(v==null){throw new IllegalArgumentException("File missing required header key '"+key+"': "+path);}
		return v;
	}

	private static String sha80File(final String path){return DigestSuffix.file(path);}
	private static void requireUnchanged(final String path, final String expectedSha256, final String label){
		final String observedSha256=sha80File(path);
		if(!expectedSha256.equals(observedSha256)){
			throw new IllegalArgumentException("Live "+label+" changed during FamilyAcceptanceConfig load: "+path+
				" (before="+expectedSha256+", after="+observedSha256+").");
		}
	}
}
