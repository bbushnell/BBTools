package test.barcode;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;

import barcode.PCRMatrixHDist;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Checks full-distance candidates and preserves the existing second-best distance margin.
 * Literal expectations exercise the protected byte-array search, including packed distances.
 * @author Shinobu
 */
public final class MatrixHammingBoundaryTest{

	/** Requires mode=baseline|fixed; baseline qualifies the full-distance rejection. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, MatrixHammingBoundaryTest.class, false);
		try{
			String mode=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(!Parser.parseStatic(arg, key, value)){
					throw new IllegalArgumentException("Unknown argument: "+arg);
				}
			}
			check("baseline".equals(mode) || "fixed".equals(mode), "Set mode=baseline or mode=fixed.");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Boundary checks require -ea to exercise the search preconditions.");
			final int fullIndex="baseline".equals(mode) ? -1 : 0;
			checkCase("sole-full", new String[] {"TTTT"}, 4, 0, fullIndex, 4, 4);
			checkCase("sole-full-threshold", new String[] {"TTTT"}, 3, 0, -1, 0, 0);
			checkCase("sole-full-margin", new String[] {"TTTT"}, 4, 1, -1, 0, 0);
			checkCase("full-tie-zero-margin", new String[] {"TTTT", "CCCC"}, 4, 0, fullIndex, 4, 4);
			checkCase("full-tie-positive-margin", new String[] {"TTTT", "CCCC"}, 4, 1, -1, 0, 0);
			checkCase("sole-near-margin-three", new String[] {"AAAT"}, 1, 3, 0, 1, 4);
			checkCase("sole-near-margin-four", new String[] {"AAAT"}, 1, 4, -1, 0, 0);
			checkCase("nearest-last", new String[] {"TTTT", "AAAT"}, 1, 3, 1, 1, 4);
			checkCase("nearest-first", new String[] {"AAAT", "TTTT"}, 1, 3, 0, 1, 4);
			checkCase("near-tie-zero-margin", new String[] {"AAAT", "AATA"}, 1, 0, 0, 1, 1);
			checkCase("near-tie-positive-margin", new String[] {"AAAT", "AATA"}, 1, 1, -1, 0, 0);
			checkCase("sole-exact-margin-four", new String[] {"AAAA"}, 0, 4, 0, 0, 4);
			checkCase("sole-exact-margin-five", new String[] {"AAAA"}, 0, 5, -1, 0, 0);
			checkCase("exact-neighbor-margin-one", new String[] {"AAAA", "AAAT"}, 0, 1, 0, 0, 1);
			checkCase("exact-neighbor-margin-two", new String[] {"AAAA", "AAAT"}, 0, 2, -1, 0, 0);
			System.out.println("MATRIX_HAMMING_BOUNDARY mode="+mode+" cases=15 failed=0");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Compares literal expected index/distances; rejected results must be the exact -1 sentinel. */
	private static void checkCase(final String label, final String[] refs, final int maximum,
			final int margin, final int expectedIndex, final int expectedDistance, final int expectedSecond){
		final Matrix matrix=new Matrix(refs);
		final long packed=matrix.search(maximum, margin);
		if(expectedIndex<0){
			check(packed==-1, label+": expected rejection, packed="+packed);
		}else{
			final int index=(int)packed;
			final int distance=(int)((packed>>>32)&0xFFFFL);
			final int second=(int)((packed>>>48)&0xFFFFL);
			check(index==expectedIndex && distance==expectedDistance && second==expectedSecond,
					label+": expected index/distance/second="+expectedIndex+"/"+expectedDistance+"/"+
					expectedSecond+", observed="+index+"/"+distance+"/"+second);
		}
		System.out.println(label+" maximum="+maximum+" margin="+margin+" index="+(int)packed);
	}

	/** Exposes the existing protected search without changing production visibility or static options. */
	private static final class Matrix extends PCRMatrixHDist{
		Matrix(final String[] names){
			super(4, 0, 0, false);
			check(names.length>0, "Fixture needs expected names for the search precondition.");
			populateExpected(Arrays.asList(names));
			refs=new byte[names.length][];
			for(int i=0; i<names.length; i++){
				check(names[i].length()==4, "Fixture references must have four bases.");
				refs[i]=names[i].getBytes(StandardCharsets.US_ASCII);
			}
		}
		long search(final int maximum, final int margin){
			return findClosestHDist(QUERY, refs, maximum, margin);
		}
		private final byte[][] refs;
	}

	/** Keeps fixture failures visible even when downstream assertions are disabled accidentally. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
	private MatrixHammingBoundaryTest(){}
	private static final byte[] QUERY="AAAA".getBytes(StandardCharsets.US_ASCII);
}
