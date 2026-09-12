package idaligner;

/**
 * Exact score-only dispatcher for the row-active and anti-diagonal SIMD
 * affine scorers. One instance is intended per worker thread and is not
 * thread-safe.
 * Scoring enforces the conservative numeric/dimension contract of {@link AffineScore32}.
 *
 * <p>AUTO falls back to row scoring when the Vector API or backend class is
 * absent or the backend class version is unsupported. FORCE_SIMD retains
 * loading failures. Other backend defects are not suppressed.</p>
 *
 * <p>The default query-length boundary is a provisional result from the
 * exclusive-node synthetic benchmark. Callers can force either implementation
 * for validation and deployment control. This class does not claim MSA score
 * equivalence or provide traceback.</p>
 */
public final class AffineScore32Hybrid {

	/** Dispatch policy selected when the instance is constructed. */
	public enum Mode {AUTO, FORCE_ROW, FORCE_SIMD}

	/** Concrete implementation selected by the most recent valid call. */
	public enum Route {NONE, ROW, SIMD}

	public AffineScore32Hybrid(final int match, final int substitution,
			final int ambiguity, final int insertionOpen,
			final int insertionExtend, final int deletionOpen,
			final int deletionExtend, final Mode mode_,
			final int simdQueryThreshold_) {
		if(mode_==null) {throw new NullPointerException("A dispatch mode is required.");}
		if(simdQueryThreshold_<0) {
			throw new IllegalArgumentException("SIMD query threshold must be nonnegative.");
		}
		row=new AffineScore32ActiveWords(match, substitution, ambiguity,
				insertionOpen, insertionExtend, deletionOpen, deletionExtend);
		final boolean useSimd=mode_==Mode.FORCE_SIMD ||
				(mode_==Mode.AUTO && vectorBackendPresent());
		simd=useSimd ? new AffineScore32DiagonalActiveSIMDWords(match, substitution, ambiguity,
				insertionOpen, insertionExtend, deletionOpen, deletionExtend) : null;
		mode=mode_;
		simdQueryThreshold=simdQueryThreshold_;
	}

	/** Returns the provisional automatic MSA-like scorer. */
	public static AffineScore32Hybrid msaLike() {
		return msaLike(Mode.AUTO);
	}

	/** Returns an MSA-like scorer with an explicit automatic or forced policy. */
	public static AffineScore32Hybrid msaLike(final Mode mode) {
		return new AffineScore32Hybrid(100, -127, 0, -395, -39, -472, -33,
				mode, DEFAULT_SIMD_QUERY_THRESHOLD);
	}

	/**
	 * Returns the exact full affine score under the same achievable-incumbent
	 * contract as both underlying scorers.
	 */
	public int score(final byte[] query, final byte[] ref, final int incumbent) {
		if(query==null) {throw new NullPointerException("A query is required.");}
		if(ref==null) {throw new NullPointerException("A reference is required.");}
		final Route route=selectRoute(query.length);
		lastRoute=route;
		if(route==Route.SIMD) {
			simdCalls++;
			return simd.score(query, ref, incumbent);
		}
		rowCalls++;
		return row.score(query, ref, incumbent);
	}

	/** Returns the exact score in an inclusive, nonempty reference window. */
	public int score(final byte[] query, final byte[] ref, final int refStart,
			final int refEnd, final int incumbent) {
		if(query==null) {throw new NullPointerException("A query is required.");}
		if(ref==null) {throw new NullPointerException("A reference is required.");}
		if(refStart<0 || refEnd<refStart || refEnd>=ref.length) {
			throw new IndexOutOfBoundsException("Invalid inclusive reference range "+
					refStart+".."+refEnd+" for length "+ref.length);
		}
		final Route route=selectRoute(query.length);
		lastRoute=route;
		if(route==Route.SIMD) {
			simdCalls++;
			return simd.score(query, ref, refStart, refEnd, incumbent);
		}
		rowCalls++;
		return row.score(query, ref, refStart, refEnd, incumbent);
	}

	private static boolean vectorBackendPresent() {
		try {
			final ClassLoader loader=AffineScore32Hybrid.class.getClassLoader();
			Class.forName("jdk.incubator.vector.IntVector", false, loader);
			Class.forName("idaligner.AffineScore32DiagonalActiveSIMDWords", false, loader);
			return true;
		} catch(ClassNotFoundException | UnsupportedClassVersionError unavailable) {
			return false;
		}
	}

	private Route selectRoute(final int queryLength) {
		if(mode==Mode.FORCE_ROW || simd==null) {return Route.ROW;}
		if(mode==Mode.FORCE_SIMD) {return Route.SIMD;}
		return queryLength>=simdQueryThreshold ? Route.SIMD : Route.ROW;
	}

	public Mode mode() {return mode;}
	public int simdQueryThreshold() {return simdQueryThreshold;}
	public Route lastRoute() {return lastRoute;}
	public long rowCalls() {return rowCalls;}
	public long simdCalls() {return simdCalls;}

	/** Nominal cells for the most recently selected scorer invocation. */
	public long lastDenseCells() {
		return lastRoute==Route.SIMD ? simd.lastDenseCells() :
			(lastRoute==Route.ROW ? row.lastDenseCells() : 0);
	}

	/** Cells processed by the most recently selected scorer invocation. */
	public long lastProcessedCells() {
		return lastRoute==Route.SIMD ? simd.lastProcessedCells() :
			(lastRoute==Route.ROW ? row.lastProcessedCells() : 0);
	}

	/** Vector cells for the most recent call, or zero when row-active ran. */
	public long lastVectorCells() {
		return lastRoute==Route.SIMD ? simd.lastVectorCells() : 0;
	}

	public static final int DEFAULT_SIMD_QUERY_THRESHOLD=250;
	private final AffineScore32ActiveWords row;
	private final AffineScore32DiagonalActiveSIMDWords simd;
	private final Mode mode;
	private final int simdQueryThreshold;
	private Route lastRoute=Route.NONE;
	private long rowCalls, simdCalls;
}
