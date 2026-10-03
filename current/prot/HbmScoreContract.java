package prot;

/** Stable semantic contract for the two-value HBM scorer result. */
public final class HbmScoreContract {

	private HbmScoreContract(){}

	/** Index 0 is raw graph score; index 1 is path-relative score. */
	public static final String VALUE=
		"HbmBundleLoader.Loaded.score_raw_index0_path_relative_index1_v1";
}
