package test.simd;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import align2.QuantumRanker;
import idaligner.AffineScore32;
import idaligner.AffineScore32Hybrid;
import idaligner.AffineScore32SeededHybrid;
import idaligner.BandedByteAligner;
import idaligner.BandedPlusAligner;
import idaligner.BandedPlusAligner2;
import idaligner.BandedPlusAligner3;
import idaligner.BandedPlusAlignerInt;
import idaligner.DriftingPlusAligner2;
import idaligner.DriftingPlusAligner3;
import idaligner.GlocalPlusAligner;
import idaligner.GlocalPlusAligner3;
import idaligner.GlocalPlusAligner4;
import idaligner.GlocalPlusAligner5;
import idaligner.GlocalPlusAligner6;
import idaligner.IDAligner;
import idaligner.QuantumPlusAligner3;
import idaligner.QuantumPlusAligner4;
import idaligner.WobblePlusAligner2;
import idaligner.WobblePlusAligner3;
import idaligner.WobblePlusAligner4;
import idaligner.WobblePlusAligner5;
import shared.Shared;
import simd.Vector;

/** Verifies fixed-256 fallbacks stay scalar when only narrower SIMD is available.
 * @author Jean */
public final class VectorWidthFallbackTest{

	/** Runs bounded width-fallback checks.
	 * @param args No arguments are expected */
	public static void main(final String[] args){
		if(args.length>0){throw new IllegalArgumentException("No arguments expected");}
		final boolean oldSimd=Shared.SIMD;
		int cases=0, skipped=0, failures=0;
		try{
			if(Vector.simd256){
				skipped+=ALIGNERS.length*ALIGN_CASES.length+
						QUANTUM_ALIGNERS.length*QUANTUM_ALIGN_CASES.length;
			}
			else{
				for(final AlignerFactory factory : ALIGNERS){
					for(final AlignCase alignCase : ALIGN_CASES){
						failures+=checkIdAlignerFallback(factory, alignCase);
						cases++;
					}
				}
				for(final AlignerFactory factory : QUANTUM_ALIGNERS){
					for(final AlignCase alignCase : QUANTUM_ALIGN_CASES){
						failures+=checkIdAlignerFallback(factory, alignCase);
						if(alignCase.expectedPos!=null){
							failures+=checkQuantumExactResult(factory, alignCase);
						}
						cases++;
					}
				}
			}

			failures+=checkWobble4Rejects(false, "scalar");
			cases++;
			if(Vector.simd256){skipped++;}
			else{
				failures+=checkWobble4Rejects(true, "no256");
				cases++;
			}

			if(Vector.simd256){skipped+=4;}
			else{
				failures+=checkQuantumRankerFallback();
				failures+=checkAffineHybridFallback();
				failures+=checkAffineSeededFallback();
				cases+=4;
			}

			System.out.println("VECTOR_WIDTH_FALLBACK_RESULT simd_bits="+
					Vector.SIMD_BITS+" simd128="+Vector.simd128+" simd256="+
					Vector.simd256+" cases="+cases+" skipped="+skipped+
					" failures="+failures);
			if(failures>0){throw new AssertionError("Vector width fallback failures="+failures);}
		}finally{
			Shared.SIMD=oldSimd;
		}
	}

	/** Verifies a guarded IDAligner uses the same scalar path when fixed-256 is absent.
	 * @param factory Aligner factory
	 * @param alignCase Bounded query/reference case
	 * @return One on mismatch, otherwise zero */
	private static int checkIdAlignerFallback(final AlignerFactory factory,
			final AlignCase alignCase){
		final int[] scalarPos=new int[4], requestedPos=new int[4];
		Shared.SIMD=false;
		final float scalar=factory.make().align(alignCase.query(), alignCase.ref(), scalarPos);
		Shared.SIMD=true;
		final float requested=factory.make().align(alignCase.query(), alignCase.ref(), requestedPos);
		if(Float.compare(scalar, requested)==0 && Arrays.equals(scalarPos, requestedPos)){
			return 0;
		}
		System.out.println("FAIL idaligner label="+alignCase.label+" aligner="+
				factory.name()+" scalar="+scalar+" requested="+requested+
				" scalarPos="+Arrays.toString(scalarPos)+
				" requestedPos="+Arrays.toString(requestedPos));
		return 1;
	}

	/** Verifies a repaired QuantumPlus exact case returns its ordinary coordinates.
	 * @param factory Aligner factory
	 * @param alignCase Exact alignment case with expected coordinates
	 * @return One on mismatch, otherwise zero */
	private static int checkQuantumExactResult(final AlignerFactory factory,
			final AlignCase alignCase){
		final int[] pos=new int[4];
		Shared.SIMD=false;
		final float identity=factory.make().align(alignCase.query(), alignCase.ref(), pos);
		if(Float.compare(identity, 1f)==0 && Arrays.equals(pos, alignCase.expectedPos)){
			return 0;
		}
		System.out.println("FAIL quantum-exact label="+alignCase.label+
				" aligner="+factory.name()+" identity="+identity+
				" expectedPos="+Arrays.toString(alignCase.expectedPos)+
				" pos="+Arrays.toString(pos));
		return 1;
	}

	/** Verifies WobblePlusAligner4 rejects a width it cannot support.
	 * @param simd Shared.SIMD setting to test
	 * @param label Diagnostic label
	 * @return One on mismatch, otherwise zero */
	private static int checkWobble4Rejects(final boolean simd, final String label){
		Shared.SIMD=simd;
		try{
			new WobblePlusAligner4().align(bytes("ACGTTGCA"), bytes("TTACGTTGCAAA"),
					new int[4]);
		}catch(UnsupportedOperationException expected){
			if(WOBBLE4_MESSAGE.equals(expected.getMessage())){return 0;}
			System.out.println("FAIL wobble4 label="+label+" unexpected message='"+
					expected.getMessage()+"'");
			return 1;
		}
		System.out.println("FAIL wobble4 label="+label+" accepted unsupported width");
		return 1;
	}

	/** Verifies QuantumRanker saturates a SIMD request to scalar without 256-bit lanes.
	 * @return One on mismatch, otherwise zero */
	private static int checkQuantumRankerFallback(){
		int failures=0;
		failures+=checkQuantumRankerCase("exact", bytes("ACGTTGCAAGTC"),
				bytes("TTACGTTGCAAGTCTT"), 2);
		failures+=checkQuantumRankerCase("substitution", bytes("ACGTTGCAAGTC"),
				bytes("TTACGTTGCCAGTCTT"), 2);
		return failures;
	}

	/** Verifies one QuantumRanker fallback case returns a supported scalar answer.
	 * @param label Diagnostic case label
	 * @param query Query bases
	 * @param ref Reference bases
	 * @param expectedStart Correct zero-based reference start hint
	 * @return One on mismatch, otherwise zero */
	private static int checkQuantumRankerCase(final String label,
			final byte[] query, final byte[] ref, final int expectedStart){
		final QuantumRanker.Result scalar=new QuantumRanker().align(query, ref, 0,
				ref.length-1, expectedStart, false);
		final QuantumRanker.Result requested=new QuantumRanker().align(query, ref, 0,
				ref.length-1, expectedStart, true);
		if(scalar.supported && requested.supported &&
				sameQuantumResult(scalar, requested) && !requested.vectorized){return 0;}
		System.out.println("FAIL quantum label="+label+" scalar="+scalar+
				" requested="+requested+" scalarSupported="+scalar.supported+
				" requestedSupported="+requested.supported+
				" requestedVectorized="+requested.vectorized);
		return 1;
	}

	/** Verifies AffineScore32Hybrid AUTO uses its row scorer without 256-bit lanes.
	 * @return One on mismatch, otherwise zero */
	private static int checkAffineHybridFallback(){
		final byte[] query=bytes("ACGTTGCAAGTC");
		final byte[] ref=bytes("TTACGTTGCAAGTCTT");
		final int incumbent=AffineScore32.msaLike().score(query, ref, 2, 13);
		final AffineScore32Hybrid hybrid=new AffineScore32Hybrid(100, -127, 0,
				-395, -39, -472, -33, AffineScore32Hybrid.Mode.AUTO, 0);
		final int score=hybrid.score(query, ref, 2, 13, incumbent);
		if(score==incumbent && hybrid.lastRoute()==AffineScore32Hybrid.Route.ROW &&
				hybrid.rowCalls()==1 && hybrid.simdCalls()==0 &&
				hybrid.lastVectorCells()==0){return 0;}
		System.out.println("FAIL affineHybrid score="+score+" incumbent="+incumbent+
				" route="+hybrid.lastRoute()+" rowCalls="+hybrid.rowCalls()+
				" simdCalls="+hybrid.simdCalls()+" vectorCells="+
				hybrid.lastVectorCells());
		return 1;
	}

	/** Verifies AffineScore32SeededHybrid AUTO skips its band and uses row exact scoring.
	 * @return One on mismatch, otherwise zero */
	private static int checkAffineSeededFallback(){
		final byte[] query=bytes("ACGTTGCAAGTC");
		final byte[] ref=bytes("TTACGTTGCAAGTCTT");
		final int incumbent=AffineScore32.msaLike().score(query, ref, 2, 13);
		final AffineScore32SeededHybrid seeded=new AffineScore32SeededHybrid(100,
				-127, 0, -395, -39, -472, -33,
				AffineScore32Hybrid.Mode.AUTO, 0);
		final int score=seeded.score(query, ref, 2, 13, 0, 2, incumbent);
		if(score==incumbent && seeded.lastExactRoute()==AffineScore32Hybrid.Route.ROW &&
				!seeded.lastBandUsable() && seeded.lastBandCells()==0 &&
				seeded.lastBandVectorCells()==0 && seeded.lastExactVectorCells()==0){return 0;}
		System.out.println("FAIL affineSeeded score="+score+" incumbent="+incumbent+
				" route="+seeded.lastExactRoute()+" bandUsable="+
				seeded.lastBandUsable()+" bandCells="+seeded.lastBandCells()+
				" bandVectorCells="+seeded.lastBandVectorCells()+
				" exactVectorCells="+seeded.lastExactVectorCells());
		return 1;
	}

	/** Compares the stable fields of two QuantumRanker results.
	 * @param a First result
	 * @param b Second result
	 * @return True when their score and coordinates match */
	private static boolean sameQuantumResult(final QuantumRanker.Result a,
			final QuantumRanker.Result b){
		return a.supported==b.supported && a.uncertain==b.uncertain &&
				a.ambiguous==b.ambiguous &&
				a.editBudgetExceeded==b.editBudgetExceeded &&
				a.prealignmentOnly==b.prealignmentOnly && a.score==b.score &&
				a.rStart==b.rStart && a.rStop==b.rStop &&
				a.matches==b.matches && a.substitutions==b.substitutions &&
				a.insertions==b.insertions && a.deletions==b.deletions &&
				a.bandwidth==b.bandwidth && a.editBudget==b.editBudget &&
				a.failureCode==b.failureCode && Float.compare(a.identity, b.identity)==0;
	}

	/** Converts a literal DNA sequence to ASCII bytes.
	 * @param bases DNA letters
	 * @return ASCII byte array */
	private static byte[] bytes(final String bases){
		return bases.getBytes(StandardCharsets.US_ASCII);
	}

	/** Builds one aligner instance. */
	private interface AlignerFactory{
		/** Builds a new aligner.
		 * @return New aligner */
		IDAligner make();

		/** Reports the aligner label.
		 * @return Aligner label */
		String name();
	}

	/** One bounded nonempty alignment case. */
	private static final class AlignCase{
		/** Constructs one alignment case without an exact-position oracle.
		 * @param label_ Diagnostic label
		 * @param query_ Query DNA letters
		 * @param ref_ Reference DNA letters */
		AlignCase(final String label_, final String query_, final String ref_){
			this(label_, query_, ref_, null);
		}

		/** Constructs one alignment case.
		 * @param label_ Diagnostic label
		 * @param query_ Query DNA letters
		 * @param ref_ Reference DNA letters
		 * @param expectedPos_ Expected position array for exact-oracle cases */
		AlignCase(final String label_, final String query_, final String ref_,
				final int[] expectedPos_){
			label=label_;
			query=query_;
			ref=ref_;
			expectedPos=expectedPos_;
		}

		/** Converts this case's query to ASCII bytes.
		 * @return Query ASCII bytes */
		byte[] query(){return bytes(query);}

		/** Converts this case's reference to ASCII bytes.
		 * @return Reference ASCII bytes */
		byte[] ref(){return bytes(ref);}

		/** Diagnostic label. */
		final String label;

		/** Query DNA letters. */
		final String query;

		/** Reference DNA letters. */
		final String ref;

		/** Expected position array for exact-oracle cases, or null. */
		final int[] expectedPos;
	}

	/** Prevents construction of this standalone test. */
	private VectorWidthFallbackTest(){}

	/** Exact exception message from WobblePlusAligner4's unsupported-width guard. */
	private static final String WOBBLE4_MESSAGE=
			"WobblePlusAligner4 requires enabled 256-bit Vector support";

	/** IDAligner implementations whose fixed-256 SIMD branches must fall back cleanly. */
	private static final AlignerFactory[] ALIGNERS={
		new AlignerFactory(){
			@Override public IDAligner make(){return new WobblePlusAligner2();}
			@Override public String name(){return "WobblePlusAligner2";}
		},
		new AlignerFactory(){
			@Override public IDAligner make(){return new WobblePlusAligner3();}
			@Override public String name(){return "WobblePlusAligner3";}
		},
		new AlignerFactory(){
			@Override public IDAligner make(){return new WobblePlusAligner5();}
			@Override public String name(){return "WobblePlusAligner5";}
		},
		new AlignerFactory(){
			@Override public IDAligner make(){return new BandedPlusAligner();}
			@Override public String name(){return "BandedPlusAligner";}
		},
		new AlignerFactory(){
			@Override public IDAligner make(){return new BandedPlusAligner2();}
			@Override public String name(){return "BandedPlusAligner2";}
		},
		new AlignerFactory(){
			@Override public IDAligner make(){return new BandedPlusAligner3();}
			@Override public String name(){return "BandedPlusAligner3";}
		},
		new AlignerFactory(){
			@Override public IDAligner make(){return new BandedByteAligner();}
			@Override public String name(){return "BandedByteAligner";}
		},
		new AlignerFactory(){
			@Override public IDAligner make(){return new BandedPlusAlignerInt();}
			@Override public String name(){return "BandedPlusAlignerInt";}
		},
		new AlignerFactory(){
			@Override public IDAligner make(){return new GlocalPlusAligner();}
			@Override public String name(){return "GlocalPlusAligner";}
		},
		new AlignerFactory(){
			@Override public IDAligner make(){return new GlocalPlusAligner3();}
			@Override public String name(){return "GlocalPlusAligner3";}
		},
		new AlignerFactory(){
			@Override public IDAligner make(){return new GlocalPlusAligner4();}
			@Override public String name(){return "GlocalPlusAligner4";}
		},
		new AlignerFactory(){
			@Override public IDAligner make(){return new GlocalPlusAligner5();}
			@Override public String name(){return "GlocalPlusAligner5";}
		},
		new AlignerFactory(){
			@Override public IDAligner make(){return new GlocalPlusAligner6();}
			@Override public String name(){return "GlocalPlusAligner6";}
		},
		new AlignerFactory(){
			@Override public IDAligner make(){return new DriftingPlusAligner2();}
			@Override public String name(){return "DriftingPlusAligner2";}
		},
		new AlignerFactory(){
			@Override public IDAligner make(){return new DriftingPlusAligner3();}
			@Override public String name(){return "DriftingPlusAligner3";}
		},
	};

	/** Quantum aligners checked on offset and same-span windows. */
	private static final AlignerFactory[] QUANTUM_ALIGNERS={
		new AlignerFactory(){
			@Override public IDAligner make(){return new QuantumPlusAligner3();}
			@Override public String name(){return "QuantumPlusAligner3";}
		},
		new AlignerFactory(){
			@Override public IDAligner make(){return new QuantumPlusAligner4();}
			@Override public String name(){return "QuantumPlusAligner4";}
		},
	};

	/** Paired alignment cases that stay finite while covering exact and mismatch paths. */
	private static final AlignCase[] ALIGN_CASES={
		new AlignCase("exact", "ACGTTGCAACGT", "TTACGTTGCAACGTAA"),
		new AlignCase("substitution", "ACGTTGCAACGT", "TTACGTTGCCACGTAA"),
	};

	/** Paired QuantumPlus cases covering offset and same-span dense-row guards. */
	private static final AlignCase[] QUANTUM_ALIGN_CASES={
		new AlignCase("quantum-offset-exact", "ACGTTGCAACGT",
				"TTACGTTGCAACGTAA", new int[] {2, 13, 12, 0}),
		new AlignCase("quantum-offset-substitution", "ACGTTGCAACGT",
				"TTACGTTGCCACGTAA"),
		new AlignCase("quantum-exact", "ACGTTGCAACGT", "ACGTTGCAACGT",
				new int[] {0, 11, 12, 0}),
		new AlignCase("quantum-substitution", "ACGTTGCAACGT", "ACGTTGCCACGT"),
	};
}
