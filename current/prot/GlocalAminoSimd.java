package prot;

import shared.Shared;

/**
 * SIMD glocal (query-global, reference-local) protein scorer: BLOSUM62 substitution, linear
 * gap, score only, no traceback. A drop-in for {@link GlocalAminoScoreOnly#scoreOnly} with
 * IDENTICAL results (same integer score, same first-maximum reference end) on every input;
 * the vector kernel ({@code simd.SIMDAminoGlocal}, striped along the query, 16 short lanes)
 * is used when {@code Shared.SIMD} is true AND the pair's scores provably fit 16-bit lanes;
 * otherwise the scalar scorer runs. Callers never see which path ran, except through
 * {@link #simdCalls}/{@link #scalarCalls} (diagnostics for the gate).
 *
 * <p>Usage shape for the family assigner: build one {@link Profile} per query (the striped
 * BLOSUM rows for that query, ~{@code 22*ceil(qLen/16)*16} shorts), then score every
 * shortlisted candidate against it; {@link #scoreOnly(byte[], byte[], int[])} builds and
 * discards a profile per call for tests and one-off use.</p>
 *
 * <p>Java-8 safety: this class imports nothing from {@code jdk.incubator.vector}; the kernel
 * class is linked lazily at the first call, which is guarded by {@code Shared.SIMD}
 * (false without AVX2 or without {@code --add-modules jdk.incubator.vector}), the same
 * pattern as {@code simd.Vector} → {@code simd.SIMD}.</p>
 *
 * <p>Design record: mag-qc {@code plans/AA_SIMD_ALIGNER_DESIGN_v1.md}; Brian's go and
 * answers 2026-09-09 06:00 PDT are recorded in that document's §6.</p>
 *
 * @author UMP45
 */
public final class GlocalAminoSimd {

	private GlocalAminoSimd(){}

	/** Linear per-residue gap cost; same value as GlocalAminoScoreOnly so scores are identical. */
	public static final int GAP=GlocalAminoScoreOnly.GAP;
	/** Profile rows: residue codes 0-19 and X_CODE (21); code 20 is never emitted by Blosum62.encode. */
	private static final int ALPHA=Blosum62.X_CODE+1;
	/** Largest BLOSUM62 entry (W/W=11); bounds the per-cell increase. */
	private static final int MAX_SUB=11;
	/** Safety margin inside the 16-bit range; every intermediate must stay within +-(32767-MARGIN). */
	private static final int MARGIN=700;
	/** Lane count of the kernel; the profile stride is built on it. Read via reflection-free constant to
	 *  keep this class free of Vector-API references: SIMDAminoGlocal.W is 16 for ShortVector.SPECIES_256. */
	static final int W=16;

	/** Diagnostic counters (not synchronized; approximate under threads, exact single-threaded). */
	public static long simdCalls=0, scalarCalls=0;

	/** A query's striped BLOSUM62 profile plus the numbers the kernel needs. Immutable; share across threads. */
	public static final class Profile {
		/** Encoded query, an OWNED snapshot (kept for the scalar fallback). Private: nothing outside this
		 *  outer class may reach the arrays, so a Profile cannot be mutated after construction. */
		private final byte[] qEnc;
		/** Striped segment length L=ceil(qLen/W) and stride n=L*W. */
		private final int segLen, n;
		/** Row a at offset a*n: S(a, q[v+k*L]) at v*W+k, 0 for padding. */
		private final short[] profile;
		/** Largest target length whose worst-case scores still fit the lanes (see fits()). */
		public final int maxRLen;

		Profile(final byte[] qEnc0){
			//Own snapshot: the profile rows are built from these residues NOW, and the scalar fallback reads qEnc
			//later; if the caller mutated its array after profile() the two paths would diverge (Yoimiya's review,
			//2026-09-09). Both must always see the same residues.
			this.qEnc=qEnc0.clone();
			final int qLen=qEnc.length;
			segLen=(qLen+W-1)/W;
			n=segLen*W;
			profile=new short[ALPHA*n];
			for(int a=0; a<ALPHA; a++){
				if(a==20){continue;}//never a target code; row stays 0
				final int base=a*n;
				for(int v=0; v<segLen; v++){
					for(int k=0; k<W; k++){
						final int p=v+k*segLen;
						profile[base+v*W+k]=(p<qLen ? (short)Blosum62.score((byte)a, qEnc[p]) : 0);
					}
				}
			}
			//Lower bound of any cell (real or padding): -(n+rLen)*GAP; upper bound: MAX_SUB*min(qLen,rLen).
			//Solve for the largest rLen keeping both inside the margin.
			final int limit=Short.MAX_VALUE-MARGIN;
			final int byLow=limit/GAP-n-1;
			final int byHigh=limit/MAX_SUB;
			maxRLen=Math.min(byLow, byHigh);
		}

		public int queryLength(){return qEnc.length;}
		/** True when the kernel's 16-bit arithmetic is provably exact for this target length. */
		public boolean fits(final int rLen){return rLen>=1 && rLen<=maxRLen && qEnc.length*MAX_SUB<=Short.MAX_VALUE-MARGIN;}
	}

	/** Builds the striped profile for one PRE-ENCODED query (Blosum62.encode output). */
	public static Profile profile(final byte[] qEnc){
		assert(qEnc!=null && qEnc.length>0) : "empty query";
		return new Profile(qEnc);
	}

	/** Per-thread scratch: two striped rows, grown never shrunk. */
	private static final class Scratch {
		short[] a=new short[0], b=new short[0];
		void ensureCapacity(final int n){if(a.length<n){a=new short[n]; b=new short[n];}}
	}
	private static final ThreadLocal<Scratch> SCRATCH=ThreadLocal.withInitial(Scratch::new);

	/**
	 * Scores a glocal alignment of the profiled query against a PRE-ENCODED reference.
	 * Result is identical to {@link GlocalAminoScoreOnly#scoreOnly}(profile's query, rEnc, out).
	 * @param out Receives {score, refEndPos (0-based inclusive, first maximum on the last row)}; length >= 2.
	 */
	public static void scoreOnly(final Profile p, final byte[] rEnc, final int[] out){
		assert(rEnc.length>0) : "empty reference";
		assert(out.length>=2) : "out needs 2 slots";
		if(Shared.SIMD && p.fits(rEnc.length)){
			final Scratch s=SCRATCH.get();
			s.ensureCapacity(p.n);
			simd.SIMDAminoGlocal.scoreStriped(p.profile, p.segLen, p.qEnc.length, rEnc, GAP, s.a, s.b, out);
			simdCalls++;
		}else{
			GlocalAminoScoreOnly.scoreOnly(p.qEnc, rEnc, out);
			scalarCalls++;
		}
		assert(out[0]>=-(p.qEnc.length)*GAP-(rEnc.length)*GAP) : "score below the all-gap floor: "+out[0];
	}

	/** Convenience: builds a profile for {@code qEnc} and scores once (allocates; tests and one-offs). */
	public static void scoreOnly(final byte[] qEnc, final byte[] rEnc, final int[] out){
		scoreOnly(profile(qEnc), rEnc, out);
	}
}
