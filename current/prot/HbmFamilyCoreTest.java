package prot;

import java.util.Arrays;

/** Exact paired-column and half-member endpoint oracles, independent of the DP implementation. */
public final class HbmFamilyCoreTest {
	public static void main(String[] args){
		final long[] depth=new long[4];
		HbmFamilyCore.accumulate(new HbmPositionModel.Result(0, 0, 3, new byte[]{'m', 'I', 'm', 'D', 'm'}), 4, depth);
		check(Arrays.equals(depth, new long[]{1, 1, 0, 1}), "Insertion/deletion columns were counted as paired coverage");
		check(Arrays.equals(HbmFamilyCore.findCore(new long[]{0, 3, 5, 2, 3, 0}, 5), new int[]{1, 4}), "Odd-member threshold or internal low-depth handling differs");
		check(Arrays.equals(HbmFamilyCore.findCore(new long[]{1, 2, 1}, 4), new int[]{1, 1}), "Even-member threshold differs");
		reject(()->HbmFamilyCore.findCore(new long[]{2, 2}, 5), "No reference column reaches");
		reject(()->HbmFamilyCore.findCore(new long[]{6}, 5), "Paired depth cannot exceed");
		reject(()->HbmFamilyCore.accumulate(new HbmPositionModel.Result(0, 0, 0, new byte[]{'m'}), 2, new long[1]), "does not consume");
		reject(()->HbmFamilyCore.accumulate(new HbmPositionModel.Result(0, 0, 0, new byte[]{'Z'}), 1, new long[1]), "Unknown core path");
		System.err.println("HBM_FAMILY_CORE_TEST_PASS paired_only=true threshold=true internal_holes=true malformed_rejected=true");
	}
	private static void reject(Runnable action, String text){
		try{action.run();}catch(IllegalArgumentException expected){check(expected.getMessage().contains(text), "Wrong rejection: "+expected); return;}
		throw new AssertionError("Malformed core fixture was accepted: "+text);
	}
	private static void check(boolean ok, String message){if(!ok){throw new AssertionError(message);}}
}
