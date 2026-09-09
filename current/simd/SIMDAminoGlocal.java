package simd;

import jdk.incubator.vector.ShortVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * Vector-API body of the striped glocal protein scorer (Farrar 2007 layout, linear gap,
 * 16-bit lanes). Reproduces {@code prot.GlocalAminoScoreOnly.scoreOnly} EXACTLY (same integer
 * score, same first-maximum reference end) for every pair whose scores fit the caller's
 * 16-bit bound; the caller ({@code prot.GlocalAminoSimd}) checks that bound and
 * {@code Shared.SIMD} before calling here, so this class is never linked on a JVM without
 * the incubator module (BBTools' Java-8 rule: lazy linking behind the Shared.SIMD gate,
 * exactly as {@code simd.Vector} → {@code simd.SIMD}).
 *
 * <p>Why striped along the QUERY: with 21 residue codes the substitution term is a table
 * lookup; on an anti-diagonal both indices change lane to lane (a gather, which we do not
 * have), whereas with one target residue per column and query positions in the lanes the
 * term is one contiguous load from a per-query profile row (mag-qc
 * {@code plans/AA_SIMD_ALIGNER_DESIGN_v1.md} §2-3). Brian 2026-09-09: "whether it is anti
 * diagonal or ref-per-lane doesn't really matter to me as long as it doesn't add too much
 * complexity and cause a time-consuming pipeline redesign" — this is a drop-in for the
 * scalar scorer, same signature shape, no pipeline change.</p>
 *
 * <p>Recurrence (0-based query position p = i-1, target column j = 1..rLen):
 * {@code H[i][j]=max(H[i-1][j-1]+S(q_i,r_j), H[i][j-1]-GAP, H[i-1][j]-GAP)},
 * {@code H[0][j]=0} (free reference start), {@code H[i][0]=-i*GAP} (query may begin before
 * the reference), score = first maximum over the last query row, exactly the scalar contract
 * ({@code GlocalAminoScoreOnly.java:17-21, :93-94}).</p>
 *
 * <p>Striped layout: segment length {@code L=ceil(qLen/W)}, W=16 lanes; vector {@code v}
 * lane {@code k} holds query position {@code p=v+k*L}. The vertical (up) dependency inside
 * one lane's run of positions {@code k*L..k*L+L-1} is carried through the v loop; the wrap
 * from position {@code k*L+L-1} to {@code (k+1)*L} (vector L-1 lane k → vector 0 lane k+1)
 * is Farrar's lazy-F correction: shift the last vector's up-term by one lane and re-sweep
 * until no lane improves (it terminates because the term strictly decreases by GAP per
 * vector and every cell is bounded below). Padding positions {@code p>=qLen} carry profile
 * 0 and can only feed positions below themselves, never a real cell.</p>
 *
 * <p>Exactness argument for 16-bit lanes: every intermediate is {@code H+S}, {@code H-GAP}
 * or a max of those; the caller admits a pair only when
 * {@code (L*W+rLen)*GAP + maxS} and {@code maxS*min(qLen,rLen)} both stay inside
 * {@code [-32000, 32000]}, so no lane ever wraps. The one sentinel, {@code Short.MIN_VALUE}
 * in the unknown lanes of the first-pass up-term, is only ever compared (max), never
 * decremented.</p>
 *
 * @author UMP45
 * @date 2026-09-09
 */
public final class SIMDAminoGlocal {

	private SIMDAminoGlocal(){}

	private static final VectorSpecies<Short> SS=ShortVector.SPECIES_256;
	/** Lane count; the striped layout and the caller's profile stride are built on it. */
	public static final int W=SS.length();

	/**
	 * Scores one target against a striped query profile.
	 * @param profile {@code short[ALPHA*n]}: row {@code a} at offset {@code a*n} holds, for vector v lane k,
	 *   {@code S(a, q[v+k*L])} or 0 for padding positions; {@code n=L*W}.
	 * @param segLen L, the striped segment length (vectors per column).
	 * @param qLen Real query length (1..L*W).
	 * @param rEnc Encoded target residues (each a valid profile row index), length >= 1.
	 * @param gap Linear per-residue gap cost (positive).
	 * @param hPrev Scratch, length >= n; overwritten.
	 * @param hCur Scratch, length >= n; overwritten.
	 * @param out Receives {score, refEndPos (0-based inclusive)}; length >= 2.
	 */
	public static void scoreStriped(final short[] profile, final int segLen, final int qLen, final byte[] rEnc,
			final int gap, final short[] hPrev0, final short[] hCur0, final int[] out){
		final int L=segLen, n=L*W, rLen=rEnc.length;
		assert(qLen>=1 && qLen<=n && rLen>=1 && gap>0) : "qLen="+qLen+" n="+n+" rLen="+rLen+" gap="+gap;
		assert(hPrev0.length>=n && hCur0.length>=n) : "scratch too small";
		short[] hPrev=hPrev0, hCur=hCur0;
		//Column 0: H[i][0]=-i*GAP for every striped position (padding included; harmless, see class javadoc).
		for(int v=0; v<L; v++){
			for(int k=0; k<W; k++){hPrev[v*W+k]=(short)(-(v+k*L+1)*gap);}
		}
		final ShortVector vGap=ShortVector.broadcast(SS, (short)gap);
		//First-pass up-term for vector 0: lane 0 sees row 0 (H[0][j]-GAP=-GAP); lanes k>=1 depend on vector L-1
		//of THIS column, unknown yet -> sentinel that never wins a max and is never decremented.
		final ShortVector vF0=ShortVector.broadcast(SS, Short.MIN_VALUE).withLane(0, (short)(-gap));
		final int lastOff=(qLen-1)%L*W+(qLen-1)/L;//striped index of the last real query position
		int best=Integer.MIN_VALUE, bestJ=0;
		//Lane shift by one (lane 0 <- fill, lane k <- v[k-1]) is done through memory: store at offset 1, reload at
		//offset 0, two ops per column. Vector.slice() is not reliably intrinsified for 16- or 32-bit lanes on this
		//JIT/AVX2 (measured 2026-09-09: a slice-per-step vertical prefix-max kernel ran at scalar speed, 1.5-2.1
		//ns/cell), so it is avoided on every path.
		//Measured cost breakdown (laptop i7-1265U, JDK 25, 300x300 random): main loop alone ~0.22 ns/cell (~11 cycles
		//per 16-lane vector for 9 vector ops -- the JIT's code, not the op count, is the floor); lazy-F adds ~0.15-0.2
		//because under the query-global boundary the below-diagonal cells are long vertical gap chains that must
		//cross every lane wrap (measured 0.55 lazy vector revisits per main visit on random pairs, 1.6 on 10%-diverged
		//homologs). Seeding the first pass from the previous column cannot help: any such seed is dominated by the
		//left term the first pass already sees. A candidate-parallel layout (one target per lane) has no lazy-F and
		//is the next step if more speed is wanted (design §3 B); it needs a batched call shape.
		final short[] shift=new short[W+1];

		for(int j=1; j<=rLen; j++){
			final int pOff=rEnc[j-1]*n;
			//Diagonal for vector 0: previous column shifted one query position: lane 0 <- H[0][j-1]=0,
			//lane k <- previous column's vector L-1 lane k-1.
			shift[0]=0;
			ShortVector.fromArray(SS, hPrev, (L-1)*W).intoArray(shift, 1);
			ShortVector diag=ShortVector.fromArray(SS, shift, 0);
			ShortVector vF=vF0;
			for(int v=0; v<L; v++){
				final ShortVector hp=ShortVector.fromArray(SS, hPrev, v*W);
				final ShortVector h=diag.add(ShortVector.fromArray(SS, profile, pOff+v*W)).max(hp.sub(vGap)).max(vF);
				h.intoArray(hCur, v*W);
				vF=h.sub(vGap);
				diag=hp;
			}
			//Lazy-F: carry the last vector's up-term across the lane wrap and re-sweep until nothing improves.
			shift[0]=(short)(-gap);//lane 0 <- -GAP (row-0 boundary), lane k <- hCur[L-1][k-1]-GAP
			vF.intoArray(shift, 1);
			vF=ShortVector.fromArray(SS, shift, 0);
			lazy:
			while(true){
				for(int v=0; v<L; v++){
					ShortVector h=ShortVector.fromArray(SS, hCur, v*W);
					final VectorMask<Short> improves=vF.compare(VectorOperators.GT, h);
					if(!improves.anyTrue()){break lazy;}//no lane improves here, so none can further down: h[v+1]>=h[v]-GAP>=vF-GAP
					h=h.max(vF);
					h.intoArray(hCur, v*W);
					vF=h.sub(vGap);
				}
				vF.intoArray(shift, 1);//wrapped again (shift[0] is still -GAP): another pass
				vF=ShortVector.fromArray(SS, shift, 0);
			}
			final int cell=hCur[lastOff];
			if(cell>best){best=cell; bestJ=j-1;}
			final short[] t=hPrev; hPrev=hCur; hCur=t;
		}
		out[0]=best;
		out[1]=bestJ;
	}
}
