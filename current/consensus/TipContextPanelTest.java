package consensus;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

/** Tests coordinate-preserving context construction and conservative edit combination. @author Fischl */
public final class TipContextPanelTest {

	/** Exhausts three-allele combinations and both context orientations. */
	public static void main(final String[] args){
		int checks=0;
		for(byte old : DNA){
			for(byte a : DNA){
				for(byte b : DNA){
					final byte[] original={old}, left={a}, right={b};
					final byte expected=a==b ? a : a==old ? b : b==old ? a : old;
					final byte[] observed=TipContextPanel.combine(original, left, right);
					check(observed.length==1 && observed[0]==expected, "Disagreeing substitutions must abstain.");
					check(original[0]==old && left[0]==a && right[0]==b, "Merge must not mutate its evidence.");
					checks+=2;
				}
			}
		}
		final byte[] target="ACTGAGTCAA".getBytes(StandardCharsets.US_ASCII);
		final TipContextGuard.Candidate donor=new TipContextGuard.Candidate("CAATTTGC".getBytes(StandardCharsets.US_ASCII), 3);
		final byte[] expected="ACTGAGTCAATTTGC".getBytes(StandardCharsets.US_ASCII);
		check(Arrays.equals(TipContextPanel.context(target, donor, false), expected), "Right context must preserve complete target.");
		check(Arrays.equals(TipContextPanel.context(TipContextGuard.reverse(target), donor, true), TipContextGuard.reverse(expected)),
				"Left context is the exact reverse complement of the same oriented construction.");
		checks+=2;
		System.out.println("TIP_CONTEXT_PANEL_TEST_PASS checks="+checks);
	}

	/** Test expectations remain fatal without assertions. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final byte[] DNA="ACGT".getBytes(StandardCharsets.US_ASCII);
}
