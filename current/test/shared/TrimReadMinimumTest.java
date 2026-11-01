package test.shared;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import shared.TrimRead;
import stream.Read;
import stream.SamLine;

/** Finite minimum-length and SAM-position checks with literal expected outputs.
 * Each case owns its arrays and optional SAM record; no files or workers are used.
 * Capped mapped cases preserve the existing right-array-trim fallback policy
 * and check reference displacement on both strands. Minimum-length cases retain
 * the complete read when no removal is allowed.
 * @author Sucrose
 */
public final class TrimReadMinimumTest{

	/** Runs every case, or one case selected with case=name; fails after reporting all results. */
	public static void main(final String[] args){
		String selected="all";
		if(args.length>1 || (args.length==1 && !args[0].startsWith("case="))){
			throw new IllegalArgumentException("Expected zero arguments or case=name");
		}
		if(args.length==1){selected=args[0].substring(5);}
		final String[] names={"minimum", "unmapped-control", "mapped-cap", "mapped-left-control", "mapped-right-control",
				"reverse-cap", "reverse-left-control", "reverse-right-control", "mapped-minimum", "reverse-minimum"};
		int cases=0, failures=0;
		for(final String name : names){
			if(!selected.equals("all") && !selected.equals(name)){continue;}
			cases++;
			final boolean ok;
			if(name.equals("minimum")){
				ok=check(name, "ACGT", 1, 0, 4, false, "ACGT", new byte[] {30, 31, 32, 33}, 0, 0);
			}else if(name.equals("unmapped-control")){
				ok=check(name, "ACGT", 1, 0, 3, false, "CGT", new byte[] {31, 32, 33}, 1, 0);
			}else if(name.equals("mapped-cap")){
				ok=check(name, "ACGATGCA", 3, 0, 7, true, "ACGATGC", new byte[] {30, 31, 32, 33, 34, 35, 36}, 1, 100);
			}else if(name.equals("mapped-left-control")){
				ok=check(name, "ACGATGCA", 1, 0, 7, true, "CGATGCA", new byte[] {31, 32, 33, 34, 35, 36, 37}, 1, 101);
			}else if(name.equals("mapped-right-control")){
				ok=check(name, "ACGATGCA", 0, 1, 7, true, "ACGATGC", new byte[] {30, 31, 32, 33, 34, 35, 36}, 1, 100);
			}else if(name.equals("reverse-cap")){
				ok=check(name, "ACGATGCA", 3, 0, 7, true, "ACGATGC", new byte[] {30, 31, 32, 33, 34, 35, 36}, 1, 101, 16);
			}else if(name.equals("reverse-left-control")){
				ok=check(name, "ACGATGCA", 1, 0, 7, true, "ACGATGC", new byte[] {30, 31, 32, 33, 34, 35, 36}, 1, 101, 16);
			}else if(name.equals("reverse-right-control")){
				ok=check(name, "ACGATGCA", 0, 1, 7, true, "CGATGCA", new byte[] {31, 32, 33, 34, 35, 36, 37}, 1, 100, 16);
			}else if(name.equals("mapped-minimum")){
				ok=check(name, "ACGT", 1, 0, 4, true, "ACGT", new byte[] {30, 31, 32, 33}, 0, 100);
			}else{
				ok=check(name, "ACGT", 1, 0, 4, true, "ACGT", new byte[] {30, 31, 32, 33}, 0, 100, 16);
			}
			if(!ok){failures++;}
		}
		if(cases==0){throw new IllegalArgumentException("Unknown case: "+selected);}
		System.out.println("TrimReadMinimumTest cases="+cases+" failures="+failures);
		if(failures!=0){throw new AssertionError("Minimum-length or retained-base SAM coordinate contract failed in "+failures+" cases");}
	}

	/** Runs a forward-strand or unmapped case. */
	private static boolean check(final String name, final String input, final int left, final int right,
			final int minimum, final boolean mapped, final String expectedBases, final byte[] expectedQual,
			final int expectedTrimmed, final int expectedPosition){
		return check(name, input, left, right, minimum, mapped, expectedBases, expectedQual, expectedTrimmed, expectedPosition, 0);
	}

	/** Calls the public dispatcher and compares bases, qualities, count, and mapped SAM state.
	 * Expected qualities are supplied independently, not derived from trimming code.
	 * SAM position is one-based; all alignments begin at position 100. Arrays are in
	 * original read orientation, matching SamLine's normal FLIP_ON_LOAD representation.
	 * On the reverse strand, removing the array's rightmost base removes the
	 * alignment's leftmost base and advances the reference start by one.
	 */
	private static boolean check(final String name, final String input, final int left, final int right,
			final int minimum, final boolean mapped, final String expectedBases, final byte[] expectedQual,
			final int expectedTrimmed, final int expectedPosition, final int flag){
		assert(expectedBases.length()==expectedQual.length) : "Literal base and quality oracle lengths must agree for "+name;
		final byte[] bases=input.getBytes(StandardCharsets.US_ASCII), quals=new byte[bases.length];
		for(int i=0; i<quals.length; i++){quals[i]=(byte)(30+i);}
		final Read read=new Read(bases, quals, name, 0, false);
		final SamLine sam;
		if(mapped){
			sam=new SamLine();
			sam.flag=flag;
			read.setStrand(sam.strand());
			sam.pos=100;
			sam.setCigar(bases.length+"M");
			sam.setSeq(bases);
			sam.setQual(quals);
			read.samline=sam;
		}else{sam=null;}
		final int trimmed=TrimRead.trimByAmount(read, left, right, minimum, false);
		final byte[] expected=expectedBases.getBytes(StandardCharsets.US_ASCII);
		boolean ok=trimmed==expectedTrimmed && Arrays.equals(read.bases, expected) && Arrays.equals(read.quality, expectedQual);
		if(sam!=null){
			final String cigarM=expected.length+"M", cigarEq=expected.length+"=";
			ok&=read.samline==sam && sam.pos==expectedPosition && (cigarM.equals(sam.cigar) || cigarEq.equals(sam.cigar))
					&& Arrays.equals(sam.seq, expected) && Arrays.equals(sam.qual, expectedQual);
		}
		System.out.println(name+" "+(ok ? "PASS" : "FAIL")+" expected="+expectedBases+"/"+Arrays.toString(expectedQual)
				+" trimmed="+expectedTrimmed+" pos="+expectedPosition+" actual="+new String(read.bases, StandardCharsets.US_ASCII)
				+"/"+Arrays.toString(read.quality)+" trimmed="+trimmed+" pos="+(sam==null ? 0 : sam.pos)
				+" cigar="+(sam==null ? "-" : sam.cigar));
		return ok;
	}
}
