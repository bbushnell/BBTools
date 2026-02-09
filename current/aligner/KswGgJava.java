package aligner;

import java.util.Arrays;

import idaligner.IDAligner;
import idaligner.Test;

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

    /**
     * Program entry point that delegates to Test class for alignment validation.
     * Uses reflection to determine calling class and passes it to Test.testAndPrint().
     * @param args Command-line arguments passed to Test framework
     * @throws Exception If reflection fails or Test execution encounters errors
     */
    public static <C extends IDAligner> void main(String[] args) throws Exception {
        StackTraceElement[] stackTrace = Thread.currentThread().getStackTrace();
        @SuppressWarnings("unchecked")
        Class<C> c=(Class<C>)Class.forName(stackTrace[(stackTrace.length<3 ? 1 : 2)].getClassName());
        Test.testAndPrint(c, args);
    }

	// Scoring constants
	static final int MATCH = 1;
	static final int MISMATCH = -1;
	static final int INS = -1;
	static final int DEL = -1;
	static final int KSW_NEG_INF = Integer.MIN_VALUE;
	static final int GAPO = 1;
	static final int GAPE = 1;

	private long loops = 0;

	/** Returns the algorithm name for identification.
	 * @return String "KswGgJava" identifying this alignment implementation */
	@Override
	public String name() {
		return "KswGgJava";
	}

	/**
	 * Aligns two sequences and returns alignment identity.
	 * Delegates to align(q, r, null) without position vector output.
	 *
	 * @param q Query sequence as byte array
	 * @param r Reference sequence as byte array
	 * @return Float identity score from alignment
	 */
	@Override
	public float align(byte[] q, byte[] r) {
		return align(q, r, (int[])null);
	}

	/**
	 * Aligns two sequences with position vector output.
	 * Delegates to align(q, r, posVector, 0, 0) with full reference length.
	 *
	 * @param q Query sequence as byte array
	 * @param r Reference sequence as byte array
	 * @param posVector Array to store alignment start/stop positions (may be null)
	 * @return Float identity score from alignment
	 */
	@Override
	public float align(byte[] q, byte[] r, int[] posVector) {
		return align(q, r, posVector, 0, 0);
	}

	/**
	 * Aligns two sequences within reference window.
	 * Delegates to align(q, r, posVector, Integer.MIN_VALUE, rStart, rStop).
	 *
	 * @param q Query sequence as byte array
	 * @param r Reference sequence as byte array
	 * @param posVector Array to store alignment positions (may be null)
	 * @param rStart Start position in reference sequence
	 * @param rStop Stop position in reference sequence
	 * @return Float identity score from alignment
	 */
	@Override
	public float align(byte[] q, byte[] r, int[] posVector, int rStart, int rStop) {
		return align(q, r, posVector, Integer.MIN_VALUE, rStart, rStop);
	}

	/**
	 * Aligns two sequences with minimum score threshold.
	 * Delegates to align(q, r, posVector, minScore, 0, 0) with full reference.
	 *
	 * @param q Query sequence as byte array
	 * @param r Reference sequence as byte array
	 * @param posVector Array to store alignment positions (may be null)
	 * @param minScore Minimum alignment score threshold
	 * @return Float identity score from alignment
	 */
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


	/** Returns the current loop count for performance tracking.
	 * @return Long value representing number of computational loops performed */
	@Override
	public long loops() {
		return loops;
	}

	/**
	 * Sets the loop counter to specified value.
	 * Used for performance analysis and benchmarking.
	 * @param i New loop count value
	 */
	@Override
	public void setLoops(long i) {
		loops = i;
	}
}