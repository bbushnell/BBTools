package var2;

import dna.Data;

/**
 * Static registry mapping sequencing platform and ploidy to a bundled network resource.
 *
 * <p>Each registered network declares which platforms it supports (a bitmask over
 * {@link VectorUMP45}'s {@code PLATFORM_*} indices) and the ploidy range it covers
 * ({@code minPloidy}..{@code maxPloidy}, inclusive). {@link #choose(int, int)} scans the
 * registry and resolves the first matching entry. It does not rank networks by measured
 * accuracy or inspect network contents. A matched but missing resource returns null;
 * there is no retry with another entry or the default in that case. When no entry matches,
 * assertions fail under {@code -ea}; with assertions disabled it resolves the default.</p>
 *
 * <p>This is invoked from the variant tools (CallVariants, CallVariants2, GradeVCF) when
 * NN scoring is requested ({@code nn}) but no explicit network was given ({@code net=}).
 * Adding a platform- or ploidy-specific network later is a one-line registry addition.</p>
 * Resource comments retain the original model-selection history; selection here performs
 * no training or validation of those historical results.
 *
 * @author UMP45
 * @contributor Sucrose
 */
public final class NNChooser{

	/*--------------------------------------------------------------*/
	/*----------------          Registry            ----------------*/
	/*--------------------------------------------------------------*/

	/** One registered network and the conditions it serves. */
	private static final class NetEntry{
		/** OR of {@code (1<<PLATFORM_*)} for every platform this net supports. */
		final int platformMask;
		/** Inclusive ploidy range this net is valid for. */
		final int minPloidy, maxPloidy;
		/** {@link Data#findPath(String, boolean)} token, e.g. {@code "?callvars_illumina_hap_dip.bbnet"}. */
		final String resource;

		/** Stores the platform mask, inclusive ploidy bounds and unresolved resource token. */
		NetEntry(int platformMask, int minPloidy, int maxPloidy, String resource){
			assert(minPloidy>=1 && maxPloidy>=minPloidy) : "Bad ploidy range "+minPloidy+".."+maxPloidy;
			assert(resource!=null && resource.startsWith("?")) : "Resource must be a findPath token: "+resource;
			this.platformMask=platformMask;
			this.minPloidy=minPloidy;
			this.maxPloidy=maxPloidy;
			this.resource=resource;
		}

		/** Tests the platform bit; callers supply a declared platform index. */
		boolean supportsPlatform(int platform){return (platformMask&(1<<platform))!=0;}

		/** True if this net covers the (platform, ploidy) combination. */
		boolean supports(int platform, int ploidy){
			return supportsPlatform(platform) && ploidy>=minPloidy && ploidy<=maxPloidy;
		}
	}

	/** The default/fallback network token, used when no registry entry matches and assertions are disabled. */
	static final String DEFAULT_RESOURCE="?callvars_illumina_hap_dip.bbnet";

	/** The Illumina polyploid (ploidy&gt;=3) network token.  Trained on synthetic tetraploid (ploidy=4) data;
	 *  serves ploidy 3 and up, generalizing across ploidy via the inverse-ploidy dims (dim0=1/ploidy, dim32=1). */
	static final String POLYPLOID_RESOURCE="?callvars_illumina_polyploid.bbnet";

	/** The Roche SBX network token.  Raw seed5 from the s_roche hap+dip net set (2026-06-25), trained on combined
	 *  haploid+diploid Roche calls (60/30/15x); serves all ploidies via the inverse-ploidy dims. */
	static final String ROCHE_RESOURCE="?callvars_roche.bbnet";

	/** The PacBio HiFi/CCS network token.  s_pe_swarm seed46 (2026-07-02), trained on Shred+BBMap (pe) HG001 PacBio
	 *  calls; best-on-pbmm2 of the 48-seed swarm under the platform-correct eval (pbmm2 30x diploid CROSSOVER_NN 55,938,
	 *  the lowest of all 48). Serves ploidy 1 (haploid) and 2 (diploid); ploidy 3+ is served by
	 *  {@link #PACBIO_POLYPLOID_RESOURCE}. (Replaced seed18, which was selected under a wrong-platform eval.) */
	static final String PACBIO_RESOURCE="?callvars_pacbio.bbnet";

	/** The PacBio HiFi/CCS polyploid (ploidy&gt;=3) network token.  s_pe_swarm seed38 (2026-07-02) -- a DIPLOID-trained
	 *  net (not a dedicated tetraploid net; the dedicated tet-jumbo nets all lost). Best ploidy-4 extrapolator of the
	 *  48-seed swarm on the pbmm2 tetraploid eval (CROSSOVER_NN 84,280 vs no-NN 244,555 = -66%); the ploidy encoding
	 *  (dim0=1/ploidy, dim32) lets a diploid net generalize to higher ploidy. Note diploid accuracy does NOT predict
	 *  tetraploid extrapolation: seed46 (best diploid) is worst on tet, seed38 (mediocre diploid) is best on tet.
	 *  Serves ploidy 3 and up. (Replaced s_pe_tetjumbo seed5, which was selected under a wrong-platform eval.) */
	static final String PACBIO_POLYPLOID_RESOURCE="?callvars_pacbio_polyploid.bbnet";

	/**
	 * Registered networks in selection order; the first match wins if entries overlap.
	 * Illumina and PacBio each have separate ploidy1–2 and ploidy3+ entries; Roche
	 * has one entry for all positive ploidies. No Nanopore entry is registered here.
	 * (Declared after the *_RESOURCE tokens: Java initializes statics in textual order.)
	 */
	private static final NetEntry[] REGISTRY={
		new NetEntry(1<<VectorUMP45.PLATFORM_ILLUMINA, 1, 2, DEFAULT_RESOURCE),
		new NetEntry(1<<VectorUMP45.PLATFORM_ILLUMINA, 3, Integer.MAX_VALUE, POLYPLOID_RESOURCE),
		new NetEntry(1<<VectorUMP45.PLATFORM_ROCHE, 1, Integer.MAX_VALUE, ROCHE_RESOURCE),
		new NetEntry(1<<VectorUMP45.PLATFORM_PACBIO, 1, 2, PACBIO_RESOURCE),
		new NetEntry(1<<VectorUMP45.PLATFORM_PACBIO, 3, Integer.MAX_VALUE, PACBIO_POLYPLOID_RESOURCE),
	};

	/*--------------------------------------------------------------*/
	/*----------------         Construction         ----------------*/
	/*--------------------------------------------------------------*/

	/** Static-only utility; never instantiated. */
	private NNChooser(){}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Resolves the first matching registry token with Data.findPath warnings disabled.
	 * Returns its result directly, including null for an unavailable matched resource.
	 * Only an unmatched platform/ploidy reaches the default fallback, and only when
	 * assertions are disabled. Does not load the network or validate its feature schema.
	 *
	 * @param platform a {@link VectorUMP45} {@code PLATFORM_*} index (e.g. {@code PLATFORM_ILLUMINA}).
	 * @param ploidy   sample ploidy ({@code >=1}).
	 * @return Location from {@link Data#findPath(String, boolean)}, or null if unresolved
	 * @throws AssertionError With assertions enabled, if ploidy is below one or no entry matches
	 */
	public static String choose(int platform, int ploidy){
		assert(ploidy>=1) : "Ploidy must be >=1; got "+ploidy;

		// Match platform first, then ploidy; return the first network that covers both.
		boolean platformKnown=false;
		for(NetEntry e : REGISTRY){
			if(e.supportsPlatform(platform)){
				platformKnown=true;
				if(e.supports(platform, ploidy)){return Data.findPath(e.resource, false);}
			}
		}

		// No registered network covers this case. Assert (so -ea crashes loud and debuggable);
		// with -da, fall through and use the default network anyway.
		assert(platformKnown) :
			"No neural network registered for platform="+platform
			+". Specify a network with net=, or run with -da to use the default ("+DEFAULT_RESOURCE+").";
		assert(false) :
			"No neural network registered for platform="+platform+" at ploidy="+ploidy
			+" (registered nets do not cover this ploidy)."
			+" Specify a network with net=, or run with -da to use the default ("+DEFAULT_RESOURCE+").";

		return Data.findPath(DEFAULT_RESOURCE, false);
	}

	/**
	 * Convenience overload that reads the platform from the static {@link VectorUMP45#platform} flag.
	 *
	 * @param ploidy sample ploidy ({@code >=1}).
	 * @return the resolved network file path, or {@code null} if unresolved.
	 */
	public static String choose(int ploidy){return choose(VectorUMP45.platform, ploidy);}

}
