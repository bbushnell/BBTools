package idaligner;

/**
 * Exact score-only affine scorer with a candidate-centered band prepass.
 * The band supplies an achievable incumbent to the exact active-frontier
 * hybrid; a better path outside the band is therefore recovered rather than
 * silently excluded. One instance is intended per worker and is not thread-safe.
 * Scoring enforces the conservative numeric/dimension contract of {@link AffineScore32}.
 *
 * <p>If the banded Vector backend is absent or has an unsupported class
 * version, the caller's achievable fallback is used directly. The requested
 * mode still controls the exact scorer. Other backend defects propagate.</p>
 *
 * <p>The delegates use reusable linear-memory arrays. This first composition
 * does not yet share their arrays between phases.</p>
 */
public final class AffineScore32SeededHybrid {

	public AffineScore32SeededHybrid(final int match, final int substitution,
			final int ambiguity, final int insertionOpen,
			final int insertionExtend, final int deletionOpen,
			final int deletionExtend, final AffineScore32Hybrid.Mode mode,
			final int simdQueryThreshold) {
		banded=bandBackendPresent() ? new AffineScore32BandedSIMD(match, substitution, ambiguity,
				insertionOpen, insertionExtend, deletionOpen, deletionExtend) : null;
		exact=new AffineScore32Hybrid(match, substitution, ambiguity,
				insertionOpen, insertionExtend, deletionOpen, deletionExtend,
				mode, simdQueryThreshold);
	}

	public static AffineScore32SeededHybrid msaLike() {
		return msaLike(AffineScore32Hybrid.Mode.AUTO);
	}

	public static AffineScore32SeededHybrid msaLike(
			final AffineScore32Hybrid.Mode mode) {
		return new AffineScore32SeededHybrid(100, -127, 0, -395, -39, -472, -33,
				mode, AffineScore32Hybrid.DEFAULT_SIMD_QUERY_THRESHOLD);
	}

	/**
	 * Scores one nonempty inclusive reference window exactly. The supplied
	 * fallback must be an achievable complete-query score on this same window.
	 * {@code centerDiagonal} is relative to {@code refStart}.
	 */
	public int score(final byte[] query, final byte[] ref, final int refStart,
			final int refEnd, final int centerDiagonal, final int halfWidth,
			final int fallbackIncumbent) {
		if(banded==null) {
			if(query==null) {throw new NullPointerException("A query is required.");}
			if(ref==null) {throw new NullPointerException("A reference is required.");}
			if(refStart<0 || refEnd<refStart || refEnd>=ref.length) {
				throw new IndexOutOfBoundsException("Invalid inclusive reference range "+
						refStart+".."+refEnd+" for length "+ref.length);
			}
			if(halfWidth<0) {throw new IllegalArgumentException("Negative half-width: "+halfWidth);}
			lastBandScore=AffineScore32.BAD;
			lastBandUsable=false;
			lastIncumbent=fallbackIncumbent;
			lastBandCells=lastBandVectorCells=0;
		} else {
			lastBandScore=banded.score(query, ref, refStart, refEnd,
					centerDiagonal, halfWidth);
			lastBandUsable=lastBandScore!=AffineScore32.BAD;
			lastIncumbent=lastBandUsable ?
					Math.max(fallbackIncumbent, lastBandScore) : fallbackIncumbent;
			lastBandCells=banded.lastVisitedCells();
			lastBandVectorCells=banded.lastVectorCells();
		}
		final int answer=exact.score(query, ref, refStart, refEnd, lastIncumbent);
		lastExactCells=exact.lastProcessedCells();
		lastExactVectorCells=exact.lastVectorCells();
		return answer;
	}

	public int lastBandScore() {return lastBandScore;}
	public boolean lastBandUsable() {return lastBandUsable;}
	public int lastIncumbent() {return lastIncumbent;}
	public long lastBandCells() {return lastBandCells;}
	public long lastBandVectorCells() {return lastBandVectorCells;}
	public long lastDenseCells() {return exact.lastDenseCells();}
	public long lastExactCells() {return lastExactCells;}
	public long lastExactVectorCells() {return lastExactVectorCells;}
	public long lastTotalCells() {return lastBandCells+lastExactCells;}
	public AffineScore32Hybrid.Route lastExactRoute() {return exact.lastRoute();}

	private static boolean bandBackendPresent() {
		try {
			final ClassLoader loader=AffineScore32SeededHybrid.class.getClassLoader();
			Class.forName("jdk.incubator.vector.IntVector", false, loader);
			Class.forName("idaligner.AffineScore32BandedSIMD", false, loader);
			return true;
		} catch(ClassNotFoundException | UnsupportedClassVersionError unavailable) {
			return false;
		}
	}

	private final AffineScore32BandedSIMD banded;
	private final AffineScore32Hybrid exact;
	private int lastBandScore, lastIncumbent;
	private boolean lastBandUsable;
	private long lastBandCells, lastBandVectorCells;
	private long lastExactCells, lastExactVectorCells;
}
