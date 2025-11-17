package aligner;

import java.util.Arrays;

class Eh {
	int h, e;
}

/**
 * Java implementation of KSW2 global alignment algorithm with affine gap penalties.
 * Global sequence alignment using dynamic programming with affine gap penalty model
 * for biological realism. Uses fixed scoring matrix with standard parameters.
 * @author Brian Bushnell
 */
public class KswGgJava implements IDAligner {

    /** Main() passes the args and class to Test to avoid redundant code */
    public static <C extends IDAligner> void main(String[] args) throws Exception {
        StackTraceElement[] stackTrace = Thread.currentThread().getStackTrace();
        @SuppressWarnings("unchecked")
        Class<C> c=(Class<C>)Class.forName(stackTrace[(stackTrace.length<3 ? 1 : 2)].getClassName());
        Test.testAndPrint(c, args);
    }

	// Scoring constants
	/** Match score for identical bases */
	static final int MATCH = 1;
	/** Mismatch penalty for differing bases */
	static final int MISMATCH = -1;
	/** Insertion penalty for gaps in reference */
	static final int INS = -1;
	/** Deletion penalty for gaps in query */
	static final int DEL = -1;
	/**
	 * Sentinel value representing negative infinity to prevent invalid alignment paths
	 */
	static final int KSW_NEG_INF = Integer.MIN_VALUE;
	/** Gap opening penalty for starting new indel */
	static final int GAPO = 1;
	/** Gap extension penalty for continuing existing indel */
	static final int GAPE = 1;

	/** Counter tracking computational loops for performance analysis */
	private long loops = 0;

	@Override
	public String name() {
		return "KswGgJava";
	}

	@Override
	public float align(byte[] q, byte[] r) {
		return align(q, r, null);
	}

	@Override
	public float align(byte[] q, byte[] r, int[] posVector) {
		return align(q, r, posVector, 0, 0);
	}

	@Override
	public float align(byte[] q, byte[] r, int[] posVector, int rStart, int rStop) {
		return align(q, r, posVector, Integer.MIN_VALUE, rStart, rStop);
	}

	@Override
	public float align(byte[] q, byte[] r, int[] posVector, int minScore) {
		return align(q, r, posVector, minScore, 0, 0);
	}

	/**
	 * Core alignment implementation using KSW2 global alignment algorithm.
	 * Performs global sequence alignment with affine gap penalties using dynamic
	 * programming. Uses Eh cell structure storing horizontal and gap extension scores.
	 *
	 * @param q Query sequence as byte array
	 * @param r Reference sequence as byte array
	 * @param posVector Array to store alignment positions (unused in this implementation)
	 * @param minScore Minimum score threshold (unused in this implementation)
	 * @param rStart Start position in reference (unused in this implementation)
	 * @param rStop Stop position in reference (unused in this implementation)
	 * @return Float alignment score from final DP cell
	 */
	private float align(byte[] q, byte[] r, int[] posVector, int minScore, int rStart, int rStop) {
		int qlen = q.length;
		int tlen = r.length;
		int w = 2; // Example bandwidth.  Should be set dynamically.
		Eh[] eh = new Eh[qlen + 1];
		for (int i = 0; i < eh.length; i++) {
			eh[i] = new Eh();
		}
		Arrays.fill(eh, new Eh());

		// fill the first row
		eh[0].h = 0;
		eh[0].e = -(GAPO + GAPO + GAPE);
		for (int j = 1; j <= qlen && j <= w; ++j) {
			eh[j].h = -(GAPO + GAPE * (j - 1));
			eh[j].e = -(GAPO + GAPO + GAPE * j);
		}
		for (int j = qlen; j <= qlen; ++j) {
			eh[j].h = eh[j].e = KSW_NEG_INF;
		}

		// DP loop
		for (int i = 0; i < tlen; ++i) {
			int f = KSW_NEG_INF;
			int h1 = 0;
			int st = (i > w) ? i - w : 0;
			int en = (i + w + 1 < qlen) ? i + w + 1 : qlen;
			h1 = (st > 0) ? KSW_NEG_INF : -(GAPO + GAPE * i);
			f = (st > 0) ? KSW_NEG_INF : -(GAPO + GAPO + GAPE * i);
			for (int j = st; j < en; ++j) {
				Eh p = eh[j];
				int h = p.h;
				int e = p.e;
				p.h = h1;
				h += (q[j] == r[i] ? MATCH : MISMATCH); // Ternary operator for efficiency
				h = (h >= e) ? h : e;
				h = (h >= f) ? h : f;
				h1 = h;
				h -= (GAPO + GAPE);
				e -= GAPE;
				e = (e > h) ? e : h;
				p.e = e;
				f -= GAPE;
				f = (f > h) ? f : h;
			}
			eh[en].h = h1;
			eh[en].e = KSW_NEG_INF;
		}

		return eh[qlen].h;
	}


	@Override
	public long loops() {
		return loops;
	}

	@Override
	public void setLoops(long i) {
		loops = i;
	}
}