package prot;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/** Literal membership differences, including equal-count swaps. @author Keqing */
public final class HbmMembershipCompareTest {
	public static void main(String[] args){
		check(list("b", "CCC", "a", "AAA*"), list("a", "aaa", "b", "CCC"), 0, 0, 0);
		check(list("a", "AAA", "b", "CCC"), list("a", "AAA", "c", "DDD"), 1, 1, 0);
		check(list("a", "AAA"), list("a", "CCC"), 0, 0, 1);
		check(list("a", "AAA", "b", "CCC"), list(), 0, 2, 0);
		check(list(), list("a", "AAA"), 1, 0, 0);
		for(boolean left : new boolean[]{false, true}){
			boolean rejected=false;
			try{HbmMembershipCompare.compare(left ? list("a", "AAA", "a", "AAA") : list(), left ? list() : list("a", "AAA", "a", "CCC"));}
			catch(IllegalArgumentException e){if(!e.getMessage().contains("Duplicate member query ID")){throw e;} rejected=true;}
			if(!rejected){throw new AssertionError("Duplicate membership was accepted");}
		}
		System.err.println("HBM_MEMBERSHIP_COMPARE_TEST_PASS reorder=true equal_count_swap=true changed_sequence=true empty=true duplicate_rejected=true");
	}
	private static List<ProteinSequence> list(String... values){
		if(values.length%2!=0){throw new AssertionError("Fixture needs ID/sequence pairs");}
		final ArrayList<ProteinSequence> out=new ArrayList<ProteinSequence>();
		for(int i=0; i<values.length; i+=2){out.add(new ProteinSequence(values[i], values[i+1]));} return out;
	}
	private static void check(List<ProteinSequence> a, List<ProteinSequence> b, long added, long removed, long changed){
		final long[] actual=HbmMembershipCompare.compare(a, b), expected={added, removed, changed};
		if(!Arrays.equals(actual, expected)){throw new AssertionError("Membership difference mismatch: "+Arrays.toString(actual));}
	}
}
