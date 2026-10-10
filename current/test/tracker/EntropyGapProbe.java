package test.tracker;

import java.nio.charset.StandardCharsets;

import tracker.EntropyTracker;

/** Checks JT006 against the contiguous-block oracle agreed by the polish lead.
 * @author Brian Bushnell, Jean */
public final class EntropyGapProbe {

	public static void main(String[] args){
		assert(args.length==0 || (args.length==1 && (args[0].equals("mode=observe") || args[0].equals("mode=verify")))) :
			"Use mode=observe or mode=verify; literal fixtures are fixed.";
		verify=args.length==1 && args[0].equals("mode=verify");
		System.out.println("case\tsequence\tallowNs\tk\twindow\tcutoff\tmonomer_min\tobserved\tproposed_contiguous\tproposal_matches");
		probe("excluded_gap_four", "AAAAAANNNNAAAAAA", false, 6);
		probe("excluded_gap_one", "AAAAAANAAAAAA", false, 6);
		probe("unequal_islands", "AAAAAAAANNNNAAAAAA", false, 8);
		probe("continuous_control", "AAAAAAAAAAAAAAAA", false, 16);
		probe("high_entropy_gap_control", "AAAAAACGTCGAAAAAA", false, 6);
		probe("included_N_control", "AAAAAANNNNAAAAAA", true, 16);
		probe("all_N_control", "NNNNNN", false, 0);
		probe("edge_N_control", "NNNAAAAAANNN", false, 6);
		probe("partial_window_control", "AAA", false, 3);
		trace("AAAAAANNNNAAAAAA");
		System.err.println("checks=9 failures="+failures);
		if(verify && failures>0){throw new AssertionError("JT006 contiguous-block failures: "+failures);}
	}

	/** Expected lengths are literal lengths of eligible contiguous islands, agreed under T007. */
	private static void probe(String name, String sequence, boolean allowNs, int proposed){
		assert(proposed>=0 && proposed<=sequence.length()) : "A contiguous block cannot exceed its source sequence: "+name;
		EntropyTracker tracker=new EntropyTracker(2, 4, false, 0.5f, true);
		final int observed=tracker.longestLowEntropyBlock(bytes(sequence), allowNs, 0.75f);
		if(observed!=proposed){failures++;}
		System.out.println(name+"\t"+sequence+"\t"+allowNs+"\t2\t4\t0.5\t0.75\t"+observed+"\t"+proposed+"\t"+(observed==proposed));
	}

	/** Observes each full window independently of the method's internal streak counter. */
	private static void trace(String sequence){
		assert(sequence.length()>=4) : "Window trace requires at least one four-base window.";
		EntropyTracker tracker=new EntropyTracker(2, 4, false, 0.5f, true);
		byte[] bases=bytes(sequence);
		System.err.println("start\tend\twindow\tundefined\tentropy\tmonomer\teligible_low");
		for(int i=0; i<bases.length; i++){
			tracker.add(bases[i]);
			if(i>=3){
				final float entropy=tracker.calcEntropy(), monomer=tracker.calcMaxMonomerFraction();
				final boolean eligible=tracker.ns()==0 && entropy<0.5f && monomer>=0.75f;
				System.err.println((i-3)+"\t"+i+"\t"+sequence.substring(i-3, i+1)+"\t"+
					tracker.ns()+"\t"+entropy+"\t"+monomer+"\t"+eligible);
			}
		}
	}

	private static byte[] bytes(String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static boolean verify;
	private static int failures;
}
