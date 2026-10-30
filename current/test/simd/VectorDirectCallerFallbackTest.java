package test.simd;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import idaligner.CrossCutAligner;
import ifa.IndelFreeAligner;
import ifa.IndelFreeAligner2;
import ifa.IndelFreeAligner3;
import ifa.IndelFreeAligner4;
import prot.Blosum62;
import prot.GlocalAminoSimd;
import shared.Shared;
import simd.Vector;
import structures.IntList;

/** Verifies direct fixed-256 callers fall back under 128-bit or absent Vector support.
 * @author Sucrose */
public final class VectorDirectCallerFallbackTest{

	/** Runs bounded direct-caller fallback checks.
	 * @param args No arguments are expected */
	public static void main(final String[] args){
		if(args.length>0){throw new IllegalArgumentException("No arguments expected");}
		final boolean oldSimd=Shared.SIMD;
		int cases=0, skipped=0, failures=0;
		try{
			if(Vector.simd256){
				skipped+=6;
			}else{
				failures+=checkIfaAlignAllPositions(); cases++;
				failures+=checkIfa2AlignAllPositions(); cases++;
				failures+=checkIfa3AlignAllPositions(); cases++;
				failures+=checkIfa4AlignAllPositions(); cases++;
				failures+=checkGlocalAmino(); cases++;
				failures+=checkCrossCut(); cases++;
			}
			System.out.println("VECTOR_DIRECT_CALLER_FALLBACK_RESULT simd_bits="+
					Vector.SIMD_BITS+" simd128="+Vector.simd128+" simd256="+
					Vector.simd256+" cases="+cases+" skipped="+skipped+
					" failures="+failures);
			if(failures>0){throw new AssertionError("Direct-caller fallback failures="+failures);}
		}finally{
			Shared.SIMD=oldSimd;
		}
	}

	/** Verifies IndelFreeAligner falls back to scalar offset enumeration.
	 * @return One on mismatch, otherwise zero */
	private static int checkIfaAlignAllPositions(){
		return checkIfaAlignAllPositions("IndelFreeAligner", new IfaAligner(){
			@Override
			public IntList alignAllPositions(final byte[] query, final byte[] ref,
					final int maxSubs, final int maxClips){
				return IndelFreeAligner.alignAllPositions(query, ref, maxSubs, maxClips);
			}
		});
	}

	/** Verifies IndelFreeAligner2 falls back to scalar offset enumeration.
	 * @return One on mismatch, otherwise zero */
	private static int checkIfa2AlignAllPositions(){
		return checkIfaAlignAllPositions("IndelFreeAligner2", new IfaAligner(){
			@Override
			public IntList alignAllPositions(final byte[] query, final byte[] ref,
					final int maxSubs, final int maxClips){
				return IndelFreeAligner2.alignAllPositions(query, ref, maxSubs, maxClips);
			}
		});
	}

	/** Verifies IndelFreeAligner3 falls back to scalar offset enumeration.
	 * @return One on mismatch, otherwise zero */
	private static int checkIfa3AlignAllPositions(){
		return checkIfaAlignAllPositions("IndelFreeAligner3", new IfaAligner(){
			@Override
			public IntList alignAllPositions(final byte[] query, final byte[] ref,
					final int maxSubs, final int maxClips){
				return IndelFreeAligner3.alignAllPositions(query, ref, maxSubs, maxClips);
			}
		});
	}

	/** Verifies IndelFreeAligner4 falls back to scalar offset enumeration.
	 * @return One on mismatch, otherwise zero */
	private static int checkIfa4AlignAllPositions(){
		return checkIfaAlignAllPositions("IndelFreeAligner4", new IfaAligner(){
			@Override
			public IntList alignAllPositions(final byte[] query, final byte[] ref,
					final int maxSubs, final int maxClips){
				return IndelFreeAligner4.alignAllPositions(query, ref, maxSubs, maxClips);
			}
		});
	}

	/** Verifies one IFA class preserves scalar offsets after a SIMD request.
	 * @param label Diagnostic class label
	 * @param aligner Class-specific alignAllPositions caller
	 * @return One on mismatch or exception, otherwise zero */
	private static int checkIfaAlignAllPositions(final String label,
			final IfaAligner aligner){
		final byte[] query=repeat('A', 40);
		final byte[] ref=bytes("CCCCCCCCCCCC"+new String(query, StandardCharsets.US_ASCII)+
				"GGGGGGGGGG");
		try{
			Shared.SIMD=false;
			final IntList scalar=aligner.alignAllPositions(query, ref, 0, 0);
			Shared.SIMD=true;
			final IntList requested=aligner.alignAllPositions(query, ref, 0, 0);
			if("12".equals(toCsv(scalar)) && sameList(scalar, requested)){return 0;}
			System.out.println("FAIL ifa label="+label+" scalar="+toCsv(scalar)+
					" requested="+toCsv(requested)+" expected=12");
		}catch(Throwable t){
			System.out.println("FAIL ifa label="+label+" exception="+t);
			t.printStackTrace(System.out);
		}
		return 1;
	}

	/** Verifies GlocalAminoSimd falls back to the scalar scorer.
	 * @return One on mismatch or exception, otherwise zero */
	private static int checkGlocalAmino(){
		final byte[] query=encode("ACDEFGHIKLMNPQRSTVWY", "q");
		final byte[] ref=encode("GGACDEFGHIKLMNPQRSTVWYGG", "r");
		final int[] scalar=new int[2], requested=new int[2];
		final long oldScalarCalls=GlocalAminoSimd.scalarCalls;
		final long oldSimdCalls=GlocalAminoSimd.simdCalls;
		try{
			Shared.SIMD=false;
			GlocalAminoSimd.scoreOnly(query, ref, scalar);
			final long scalarCallsAfterScalar=GlocalAminoSimd.scalarCalls;
			final long simdCallsAfterScalar=GlocalAminoSimd.simdCalls;
			Shared.SIMD=true;
			GlocalAminoSimd.scoreOnly(query, ref, requested);
			if(Arrays.equals(scalar, requested) && scalar[0]==116 && scalar[1]==21 &&
					GlocalAminoSimd.scalarCalls==scalarCallsAfterScalar+1 &&
					GlocalAminoSimd.simdCalls==simdCallsAfterScalar){return 0;}
			System.out.println("FAIL glocal scalar="+Arrays.toString(scalar)+
					" requested="+Arrays.toString(requested)+" expected=[116, 21]");
		}catch(Throwable t){
			System.out.println("FAIL glocal exception="+t);
			t.printStackTrace(System.out);
		}finally{
			GlocalAminoSimd.scalarCalls=oldScalarCalls;
			GlocalAminoSimd.simdCalls=oldSimdCalls;
		}
		return 1;
	}

	/** Verifies CrossCutAligner falls back to scalar diagonal processing.
	 * @return One on mismatch or exception, otherwise zero */
	private static int checkCrossCut(){
		final byte[] query=bytes("ACGTTGCAACGT");
		final byte[] ref=bytes("TTACGTTGCAACGTAA");
		final int[] scalarPos=new int[2], requestedPos=new int[2];
		try{
			Shared.SIMD=false;
			final float scalar=new CrossCutAligner().align(query, ref, scalarPos);
			Shared.SIMD=true;
			final float requested=new CrossCutAligner().align(query, ref, requestedPos);
			if(Float.compare(scalar, 1f)==0 && Float.compare(requested, scalar)==0 &&
					Arrays.equals(scalarPos, requestedPos) &&
					Arrays.equals(scalarPos, new int[] {2, 13})){return 0;}
			System.out.println("FAIL crossCut scalar="+scalar+" requested="+requested+
					" scalarPos="+Arrays.toString(scalarPos)+
					" requestedPos="+Arrays.toString(requestedPos)+
					" expectedPos=[2, 13]");
		}catch(Throwable t){
			System.out.println("FAIL crossCut exception="+t);
			t.printStackTrace(System.out);
		}
		return 1;
	}

	/** Compares two IntList values.
	 * @param a First list
	 * @param b Second list
	 * @return True if both lists have the same values in the same order */
	private static boolean sameList(final IntList a, final IntList b){
		if(a==null || b==null){return a==b;}
		if(a.size!=b.size){return false;}
		for(int i=0; i<a.size; i++){
			if(a.get(i)!=b.get(i)){return false;}
		}
		return true;
	}

	/** Converts an IntList to a compact comma-separated diagnostic string.
	 * @param list List to convert, possibly null
	 * @return Comma-separated values, or an empty string for null/empty */
	private static String toCsv(final IntList list){
		if(list==null || list.isEmpty()){return "";}
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<list.size; i++){
			if(i>0){sb.append(',');}
			sb.append(list.get(i));
		}
		return sb.toString();
	}

	/** Encodes amino acid letters to Blosum62 residue codes.
	 * @param aa Amino acid letters
	 * @param id Sequence identifier for diagnostics
	 * @return Encoded residues */
	private static byte[] encode(final String aa, final String id){
		return Blosum62.encode(bytes(aa), id);
	}

	/** Converts a literal ASCII sequence to bytes.
	 * @param bases Sequence letters
	 * @return ASCII byte array */
	private static byte[] bytes(final String bases){
		return bases.getBytes(StandardCharsets.US_ASCII);
	}

	/** Builds a repeated ASCII sequence.
	 * @param base Base to repeat
	 * @param length Number of copies
	 * @return ASCII byte array */
	private static byte[] repeat(final char base, final int length){
		final byte[] out=new byte[length];
		Arrays.fill(out, (byte)base);
		return out;
	}

	/** One IFA alignAllPositions call. */
	private interface IfaAligner{
		/** Calls the class-specific alignAllPositions method.
		 * @param query Query bases
		 * @param ref Reference bases
		 * @param maxSubs Maximum substitution or paid-clipping cost
		 * @param maxClips Total free clipping allowance
		 * @return Accepted offsets, or null */
		IntList alignAllPositions(byte[] query, byte[] ref, int maxSubs, int maxClips);
	}

	/** Prevents construction of this standalone test. */
	private VectorDirectCallerFallbackTest(){}
}
