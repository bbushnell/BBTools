package test.tracker;

import java.nio.charset.StandardCharsets;
import java.util.function.Supplier;

import repeat.Palindrome;
import tracker.EntropyTracker;
import tracker.PalindromeTracker;
import tracker.PolymerTracker;
import tracker.ReadStats;

/** Literal regression oracles for tracker lifecycle, coordinates, and window counts.
 * @author Brian Bushnell, Jean */
public final class TrackerRegressionTest {

	/** mode=observe prints baseline discrepancies; mode=verify requires all oracles. */
	public static void main(String[] args){
		assert(args.length==1 && (args[0].equals("mode=observe") || args[0].equals("mode=verify"))) :
			"Use mode=observe or mode=verify; the mode controls whether known baseline defects fail the run.";
		final boolean verify=args[0].equals("mode=verify");
		System.out.println("case\texpected\tactual\tpass");
		check("polymer_add_after_query", "2", () -> {
			PolymerTracker p=new PolymerTracker();
			p.add(bytes("AAAA"));
			p.getCountCumulative((byte)'A', 4);
			p.add(bytes("AAAA"));
			return ""+p.getCountCumulative((byte)'A', 4);
		});
		check("polymer_reset_after_query", "0", () -> {
			PolymerTracker p=new PolymerTracker();
			p.add(bytes("AAAA"));
			p.accumulate();
			p.reset();
			return ""+p.getCountCumulative((byte)'A', 4);
		});
		check("polymer_merge_after_query", "2", () -> {
			PolymerTracker p=new PolymerTracker(), q=new PolymerTracker();
			p.add(bytes("AAAA"));
			q.add(bytes("AAAA"));
			p.accumulate();
			p.add(q);
			return ""+p.getCountCumulative((byte)'A', 4);
		});
		check("polymer_longest_per_sequence_control", "0,1,1", () -> {
			PolymerTracker p=new PolymerTracker();
			p.add(bytes("AATAAA"));
			return p.getCount((byte)'A', 2)+","+p.getCount((byte)'A', 3)+","+p.getCountCumulative((byte)'A', 2);
		});
		check("polymer_per_run_control", "1,1,2", () -> {
			PolymerTracker.PER_SEQUENCE=false;
			try{
				PolymerTracker p=new PolymerTracker();
				p.add(bytes("AATAAA"));
				return p.getCount((byte)'A', 2)+","+p.getCount((byte)'A', 3)+","+p.getCountCumulative((byte)'A', 2);
			}finally{PolymerTracker.PER_SEQUENCE=true;}
		});
		check("palindrome_asymmetric_tails", "1,1,0,1", () -> {
			PalindromeTracker p=new PalindromeTracker();
			p.add(new Palindrome(3, 8, 3, 0), 0, 12);
			return p.tailList.get(3)+","+p.tailList.get(4)+","+p.tailList.get(9)+","+p.tailDifList.get(1);
		});
		check("palindrome_symmetric_tails", "2,1", () -> {
			PalindromeTracker p=new PalindromeTracker();
			p.add(new Palindrome(5, 10, 3, 0), 2, 13);
			return p.tailList.get(3)+","+p.tailDifList.get(0);
		});
		check("palindrome_arm_loop_control", "1,1,1", () -> {
			PalindromeTracker p=new PalindromeTracker();
			p.add(new Palindrome(3, 8, 3, 0), 0, 12);
			return p.plenList.get(3)+","+p.loopList.get(0)+","+p.found;
		});
		check("entropy_monomer_eviction", "1.0", () -> {
			EntropyTracker p=new EntropyTracker(2, 4, false);
			p.clear();
			for(byte b : bytes("ACGTTTT")){p.add(b);}
			return ""+p.calcMaxMonomerFraction();
		});
		check("entropy_amino_symbols", "1.0", () -> {
			EntropyTracker p=new EntropyTracker(1, 4, true);
			p.clear();
			for(byte b : bytes("WWWW")){p.add(b);}
			return ""+p.calcMaxMonomerFraction();
		});
		check("entropy_fresh_constructor", "1.0", () -> {
			EntropyTracker p=new EntropyTracker(2, 4, false);
			for(byte b : bytes("ACGT")){p.add(b);}
			return ""+p.calcEntropy();
		});
		check("entropy_clear_entropy_control", "1.0", () -> {
			EntropyTracker p=new EntropyTracker(2, 4, false);
			p.clear();
			for(byte b : bytes("ACGT")){p.add(b);}
			return ""+p.calcEntropy();
		});
		check("quality_short_high", "3", () -> quality(new byte[]{40, 40, 40}));
		check("quality_empty", "0", () -> quality(new byte[0]));
		check("quality_low_control", "2", () -> quality(new byte[]{40, 35, 20}));
		System.err.println("checks="+checks+" failures="+failures);
		if(verify && failures>0){throw new AssertionError("Tracker regression oracles failed: "+failures);}
	}

	/** Independent prefix-length oracle uses short literal quality arrays and eight available bins. */
	private static String quality(byte[] quals){
		final int oldMax=ReadStats.MAXLEN;
		ReadStats.clear();
		ReadStats.MAXLEN=8;
		ReadStats.COLLECT_QUALITY_STATS=true;
		try{
			ReadStats p=new ReadStats(false);
			p.addToQualityHistogram(quals, 0);
			return ""+p.q30(0);
		}finally{
			ReadStats.clear();
			ReadStats.MAXLEN=oldMax;
		}
	}

	private static byte[] bytes(String s){return s.getBytes(StandardCharsets.US_ASCII);}

	/** Continue after individual failures so baseline evidence includes every independent case. */
	private static void check(String name, String expected, Supplier<String> operation){
		assert(name!=null && expected!=null && operation!=null) : "Each regression needs an identity, literal oracle and operation.";
		String actual;
		try{actual=operation.get();}
		catch(AssertionError | RuntimeException e){actual=e.getClass().getSimpleName();}
		final boolean pass=expected.equals(actual);
		checks++;
		if(!pass){failures++;}
		System.out.println(name+"\t"+expected+"\t"+actual+"\t"+pass);
	}

	private static int checks=0, failures=0;
}
