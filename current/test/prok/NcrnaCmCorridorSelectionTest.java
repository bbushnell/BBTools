package test.prok;

import java.util.Arrays;
import prok.NcrnaCmCorridorSelection;
import structures.IntList;

/** Literal support, alias-group, trimming-order and held-out-fold decisions.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCmCorridorSelectionTest{
	public static void main(final String[] args){
		check(args.length==0, "Whole-path selection fixtures use only declared literal populations");
		final NcrnaCmCorridorSelection selector=new NcrnaCmCorridorSelection(100, 10, 5);
		final IntList calls=new IntList();
		final int[] groups=new int[130], departures=new int[130];
		final byte[] retained=new byte[130];
		for(int i=0; i<101; i++){
			calls.add(i);
			groups[i]=i==0 ? 0 : 1+(i-1)%9;
		}
		departures[0]=99;
		check(selector.select(calls, groups, departures, -1, retained, 1)==101 && (retained[0]&1)!=0,
			"Dropping the largest departure would erase the tenth group, so trimming must stop without skipping to an easier trajectory");
		groups[100]=10;
		check(selector.select(calls, groups, departures, -1, retained, 2)==100 && (retained[0]&2)==0,
			"With eleven groups, one whole trajectory can be removed while retaining100calls from10groups");
		Arrays.fill(departures, 0);
		check(selector.select(calls, groups, departures, -1, retained, 4)==100 && (retained[0]&4)==0 && (retained[1]&4)!=0,
			"Equal departures deterministically remove the lower ordinal first, preserving the other fit masks");
		Arrays.fill(groups, 0);
		check(selector.select(calls, groups, departures, -1, retained, 8)==0,
			"Different formal call IDs from the same resolved assembly group cannot satisfy ten independent groups");
		calls.clear();
		for(int i=0; i<130; i++){
			calls.add(i);
			groups[i]=i%20;
		}
		departures[1]=-1;
		check(selector.select(calls, groups, departures, 0, retained, 16)==102,
			"Fold0excludes26calls before fitting; one unavailable anchor leaves103, and one allowed whole-path removal leaves102");
		for(int i=0; i<130; i++){
			if(groups[i]%5==0 || i==1){check((retained[i]&16)==0, "Held-out groups and unavailable anchors never support learned bounds");}
		}
		percentages();
		System.err.println("NCRNA_CM_CORRIDOR_SELECTION_FIXTURE_PASS");
	}

	private static void percentages(){
		final IntList calls=new IntList();
		final int[] groups=new int[201], departures=new int[201];
		for(int i=0; i<201; i++){
			calls.add(i);
			groups[i]=i%20;
			departures[i]=i;
		}
		final byte[] omitted=new byte[201], explicit=new byte[201];
		new NcrnaCmCorridorSelection(100, 10, 5).select(calls, groups, departures, -1, omitted, 1);
		new NcrnaCmCorridorSelection(100, 10, 5, 1).select(calls, groups, departures, -1, explicit, 1);
		check(Arrays.equals(omitted, explicit), "Omitted trimming must select exactly the explicit one-percent population");
		for(final int percent:new int[]{0, 5, 10, 20}){
			final byte[] retained=new byte[201];
			final int removed=201*percent/100;
			final int kept=new NcrnaCmCorridorSelection(100, 10, 5, percent).select(calls, groups, departures, -1, retained, 1);
			check(kept==201-removed, "Whole-path percentage rounds down on eligible paths before removal");
			for(int i=0; i<201; i++){
				check((retained[i]!=0)==(i<kept), "Only the greatest departures may be removed within the declared budget");
			}
		}
		check(new NcrnaCmCorridorSelection(100, 10, 5, 99).select(calls, groups, departures, -1, new byte[201], 1)==100,
			"A large budget must stop at retained-call support instead of forcing the percentage");
		check(new NcrnaCmCorridorSelection(100, 10, 5, 20).select(calls, groups, departures, 0, new byte[201], 1)==128,
			"Held-out groups are removed before applying twenty percent to the160eligible training trajectories");
		for(final int invalid:new int[]{-1, 100}){
			boolean rejected=false;
			try{new NcrnaCmCorridorSelection(100, 10, 5, invalid);}
			catch(IllegalArgumentException expected){rejected=true;}
			check(rejected, "Invalid trimming percentages cannot reach trajectory selection");
		}
	}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
}
