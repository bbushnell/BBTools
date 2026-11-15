package test.icecream;

import java.util.Arrays;

import icecream.ZMW;
import stream.Read;

/** Literal read-selection cases, including wells whose interior reads were discarded.
 * @author Brian Bushnell, Sucrose
 */
public final class ZMWSelectionTest{

	/** Runs the eighteen baseline observations and two mutable-ID/positive-pass controls.
	 * @param args No arguments are accepted
	 */
	public static void main(String[] args){
		require(args.length==0, "The selection fixture takes no arguments");
		require(ZMW.class.desiredAssertionStatus(), "Production assertions must be enabled");
		System.out.println("case\tinclude_discarded\texpected_index\tobserved_index\tresult");
		int failures=0;
		failures+=observe("empty", new int[]{}, new boolean[]{}, false, -1);
		failures+=observe("one_kept", new int[]{4}, new boolean[]{false}, false, 0);
		failures+=observe("one_discarded", new int[]{4}, new boolean[]{true}, false, -1);
		failures+=observe("two_one_kept", new int[]{9, 4}, new boolean[]{true, false}, false, 1);
		failures+=observe("normal_middle", new int[]{9, 3, 7, 5, 11}, new boolean[5], false, 3);
		failures+=observe("upper_middle", new int[]{9, 3, 5, 11}, new boolean[4], false, 2);
		failures+=observe("missing_middle_tie", new int[]{6, 6, 6}, new boolean[]{false, true, false}, false, 0);
		failures+=observe("missing_middle_last_longest", new int[]{4, 5, 8}, new boolean[]{false, true, false}, false, 2);
		failures+=observe("only_last_survives", new int[]{4, 5, 8}, new boolean[]{true, true, false}, false, 2);
		failures+=observe("only_first_survives", new int[]{8, 5, 4}, new boolean[]{false, true, true}, false, 0);
		failures+=observe("all_discarded", new int[]{4, 5, 8}, new boolean[]{true, true, true}, false, -1);
		failures+=observe("include_all_discarded", new int[]{4, 5, 8}, new boolean[]{true, true, true}, true, 1);
		failures+=observe("zero_length_middle", new int[]{4, 0, 5}, new boolean[3], false, 1);
		failures+=observe("one_eligible_middle", new int[]{20, 3, 7, 5, 22}, new boolean[]{false, true, false, true, false}, false, 2);
		failures+=observe("first_middle_tie", new int[]{20, 5, 5, 5, 22}, new boolean[5], false, 1);
		System.out.println("id_case\texpected_id\tobserved_id\tresult");
		failures+=observeID("empty_id", null, -1);
		failures+=observeID("coordinate_header_id", "m/42/0_4", 42);
		failures+=observeInvalidPasses();
		failures+=observeMutatedID();
		failures+=observePositivePasses();
		require(failures==0, failures+" selections differ from the literal eligible-read oracle");
	}

	/** Uses literal expected positions and checks that selection preserves the original list and flags. */
	private static int observe(String name, int[] lengths, boolean[] discarded, boolean includeDiscarded, int expected){
		require(lengths.length==discarded.length, "Each fixture read needs one explicit discard flag");
		final ZMW well=new ZMW(lengths.length);
		final Read[] original=new Read[lengths.length];
		for(int i=0; i<lengths.length; i++){
			final byte[] bases=new byte[lengths[i]];
			Arrays.fill(bases, (byte)'A');
			final Read read=new Read(bases, null, "movie/11/0_"+lengths[i], i);
			read.setDiscarded(discarded[i]);
			original[i]=read;
			well.add(read);
		}
		final Read selected=well.medianRead(includeDiscarded);
		int observed=-1;
		for(int i=0; i<original.length; i++){
			if(selected==original[i]){observed=i;}
		}
		boolean pass=observed==expected && selected==(expected<0 ? null : original[expected]);
		pass&=well.size()==original.length;
		for(int i=0; i<original.length && i<well.size(); i++){
			pass&=well.get(i)==original[i] && original[i].discarded()==discarded[i] && original[i].length()==lengths[i];
		}
		System.out.println(name+"\t"+includeDiscarded+"\t"+expected+"\t"+observed+"\t"+(pass ? "PASS" : "FAIL"));
		return pass ? 0 : 1;
	}

	/** Checks literal header IDs independently of list ordinals and Read.numericID. */
	private static int observeID(String name, String header, int expected){
		final ZMW well=new ZMW();
		if(header!=null){well.add(new Read(new byte[]{'A', 'C', 'G', 'T'}, null, header, 7));}
		final int observed=well.zid();
		final boolean pass=observed==expected;
		System.out.println(name+"\t"+expected+"\t"+observed+"\t"+(pass ? "PASS" : "FAIL"));
		return pass ? 0 : 1;
	}

	/** A zero internal median cannot produce a meaningful endpoint length ratio. */
	private static int observeInvalidPasses(){
		final ZMW well=new ZMW();
		well.add(new Read(new byte[]{'A', 'C', 'G', 'T'}, null, "m/42/0_4", 0));
		well.add(new Read(new byte[0], null, "m/42/4_4", 1));
		well.add(new Read(new byte[]{'A', 'C', 'G', 'T', 'A'}, null, "m/42/4_9", 2));
		String observed;
		try{observed=Float.toString(well.estimatePasses());}
		catch(IllegalStateException expected){observed="IllegalStateException";}
		final boolean pass="IllegalStateException".equals(observed);
		System.out.println("pass_case\texpected\tobserved\tresult");
		System.out.println("zero_internal_median\tIllegalStateException\t"+observed+"\t"+(pass ? "PASS" : "FAIL"));
		return pass ? 0 : 1;
	}

	/** Looks up IDs after replacement, clearing and repopulation of the same inherited list. */
	private static int observeMutatedID(){
		final ZMW well=new ZMW();
		well.add(new Read(new byte[]{'A', 'C', 'G', 'T'}, null, "m/42/0_4", 7));
		final int initial=well.zid();
		well.set(0, new Read(new byte[]{'A', 'C', 'G', 'T'}, null, "m/73/0_4", 7));
		final int replaced=well.zid();
		well.clear();
		final int cleared=well.zid();
		well.add(new Read(new byte[]{'A', 'C', 'G', 'T'}, null, "m/9/0_4", 7));
		final String observed=initial+":"+replaced+":"+cleared+":"+well.zid();
		final boolean pass="42:73:-1:9".equals(observed);
		System.out.println("mutation_case\texpected\tobserved\tresult");
		System.out.println("current_first_header\t42:73:-1:9\t"+observed+"\t"+(pass ? "PASS" : "FAIL"));
		return pass ? 0 : 1;
	}

	/** Empty endpoints with a positive internal median contribute exactly zero partial passes. */
	private static int observePositivePasses(){
		final ZMW well=new ZMW();
		well.add(new Read(new byte[0], null, "m/42/0_0", 0));
		well.add(new Read(new byte[]{'A', 'C', 'G', 'T'}, null, "m/42/0_4", 1));
		well.add(new Read(new byte[0], null, "m/42/4_4", 2));
		final float observed=well.estimatePasses();
		final boolean pass=observed==1.0f;
		System.out.println("pass_case\texpected\tobserved\tresult");
		System.out.println("positive_internal_median\t1.0\t"+observed+"\t"+(pass ? "PASS" : "FAIL"));
		return pass ? 0 : 1;
	}

	/** Keeps the fixture oracle active independently of assertion mode. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Static-only fixture. */
	private ZMWSelectionTest(){}
}
