package var2;

import dna.Data;

/**
 * Static registry that selects the best bundled neural network for variant scoring,
 * based on sequencing platform and ploidy.
 *
 * <p>Each registered network declares which platforms it supports (a bitmask over
 * {@link VectorUMP45}'s {@code PLATFORM_*} indices) and the ploidy range it covers
 * ({@code minPloidy}..{@code maxPloidy}, inclusive). {@link #choose(int, int)} scans the
 * registry, matching platform first and then ploidy, and returns the file path of the
 * best-matching network. When nothing matches it asserts (so a run with {@code -ea}
 * fails loudly and debuggably) and otherwise falls back to the default network
 * (so a run with {@code -da} proceeds using the default).</p>
 *
 * <p>This is invoked from the variant tools (CallVariants, CallVariants2, GradeVCF) when
 * NN scoring is requested ({@code nn}) but no explicit network was given ({@code net=}).
 * Adding a platform- or ploidy-specific network later is a one-line registry addition.</p>
 *
 * @author UMP45
 */
public final class NNChooser {

	/*--------------------------------------------------------------*/
	/*----------------          Registry            ----------------*/
	/*--------------------------------------------------------------*/

	/** One registered network and the conditions it serves. */
	private static final class NetEntry {
		/** OR of {@code (1<<PLATFORM_*)} for every platform this net supports. */
		final int platformMask;
		/** Inclusive ploidy range this net is valid for. */
		final int minPloidy, maxPloidy;
		/** {@link Data#findPath(String, boolean)} token, e.g. {@code "?callvars_illumina_hap_dip.bbnet"}. */
		final String resource;

		NetEntry(int platformMask, int minPloidy, int maxPloidy, String resource){
			assert minPloidy>=1 && maxPloidy>=minPloidy : "Bad ploidy range "+minPloidy+".."+maxPloidy;
			assert resource!=null && resource.startsWith("?") : "Resource must be a findPath token: "+resource;
			this.platformMask=platformMask;
			this.minPloidy=minPloidy;
			this.maxPloidy=maxPloidy;
			this.resource=resource;
		}

		/** True if this net's mask includes the platform bit. */
		boolean supportsPlatform(int platform){return (platformMask & (1<<platform))!=0;}

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

	/** The PacBio HiFi/CCS network token.  s_pe_swarm seed18 (2026-06-28), trained on Shred+BBMap (pe) HG001 PacBio
	 *  calls; best-on-pbmm2 of a 48-seed swarm, and won concordance across pbmm2/BBMap/MapPacBio test sets.
	 *  Serves ploidy 1 (haploid) and 2 (diploid); ploidy 3+ is served by {@link #PACBIO_POLYPLOID_RESOURCE}. */
	static final String PACBIO_RESOURCE="?callvars_pacbio.bbnet";

	/** The PacBio HiFi/CCS polyploid (ploidy&gt;=3) network token.  s_pe_tetjumbo seed5 (2026-07-02), trained on a
	 *  75%-tetraploid + 12%-diploid + 12%-haploid PacBio jumbo (HG001+HG002 pe-shred, ploidy-encoded so it
	 *  generalizes across ploidy via dim0=1/ploidy and dim32). Best of a 16-seed swarm on the tetraploid eval
	 *  (CROSSOVER_NN 116,378 vs the diploid net's 170,372 = -31.7%). Serves ploidy 3 and up. */
	static final String PACBIO_POLYPLOID_RESOURCE="?callvars_pacbio_polyploid.bbnet";

	/**
	 * Registered networks, most-specific first (choose() returns the first match).
	 * Two Illumina nets covering disjoint ploidy ranges:
	 *   - hap+dip combined net (DEFAULT_RESOURCE): ploidy 1 (haploid, dim32=0) and 2 (diploid, dim32=1).
	 *   - polyploid net (POLYPLOID_RESOURCE): ploidy 3 and up (synthetic-tetraploid-trained).
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
	 * Resolve the file path of the best network for a platform and ploidy.
	 *
	 * @param platform a {@link VectorUMP45} {@code PLATFORM_*} index (e.g. {@code PLATFORM_ILLUMINA}).
	 * @param ploidy   sample ploidy ({@code >=1}).
	 * @return the resolved network file path (via {@link Data#findPath(String, boolean)}),
	 *         or {@code null} if even the default cannot be located on disk.
	 */
	public static String choose(int platform, int ploidy){
		assert ploidy>=1 : "Ploidy must be >=1; got "+ploidy;

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
		assert platformKnown :
			"No neural network registered for platform="+platform
			+". Specify a network with net=, or run with -da to use the default ("+DEFAULT_RESOURCE+").";
		assert false :
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
