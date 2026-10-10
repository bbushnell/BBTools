package test.tracker;

import java.nio.charset.StandardCharsets;

import tracker.EntropyTracker;

/** Independent fresh-window count oracle and matched checks for JT012.
 * @author Brian Bushnell, Jean */
public final class StrandednessWindowProbe {

	public static void main(String[] args){
		boolean verify=false;
		for(String arg : args){
			if(arg.equals("mode=verify")){verify=true;}
			else if(!arg.equals("mode=observe")){throw new IllegalArgumentException("Unknown JT012 argument: "+arg);}
		}
		System.out.println("case\tk\twindow\tsequence\toracle\tobserved\tagrees");
		for(int k=2; k<=3; k++){
			final String boundary=(k==2 ? "AATT" : "AAATTT");
			final int window=k+1;
			checkLiteral(boundary, k, window, 0);
			checkLiteral(boundary, k, boundary.length(), k==2 ? 0.5 : 1.0/3);
			checkLiteral("AAAAAA", k, 6, 0);
			checkLiteral("NNNNNN", k, 4, 0);
			checkLiteral(boundary, k, k, 0);
			checkLiteral("A", k, k+1, 0);
			for(int entry=0; entry<(k==2 ? 2 : 1); entry++){
				final boolean direct=entry==1;
				final String prefix=(direct ? "directK2_" : "general_");
				observe(prefix+"boundary_null", boundary, k, window, null, direct);
				observe(prefix+"boundary_fresh", boundary, k, window, new int[1<<(2*k)], direct);
				final int[] dirty=new int[1<<(2*k)];
				dirty[dirty.length-1]=4;//Four all-T kmers from an earlier sequence must not count here.
				observe(prefix+"dirty_scratch", "AAAAAA", k, 6, dirty, direct);
				final int[] reused=new int[1<<(2*k)];
				observe(prefix+"repeat_first", boundary, k, window, reused, direct);
				observe(prefix+"repeat_second", boundary, k, window, reused, direct);
				final int[] prior=new int[1<<(2*k)];
				call("AAAAAA", k, 6, prior, direct);
				observe(prefix+"previous_sequence", "TTTTTT", k, 6, prior, direct);
				observe(prefix+"whole_control", boundary, k, boundary.length(), null, direct);
				observe(prefix+"short_control", boundary, k, boundary.length()+2, null, direct);
				observe(prefix+"single_strand_control", "AAAAAA", k, 4, null, direct);
				observe(prefix+"undefined_control", "NNNNNN", k, 4, null, direct);
				observe(prefix+"gap", "AAANNTTT", k, 4, null, direct);
				observe(prefix+"window_equals_k", boundary, k, k, null, direct);
				observe(prefix+"shorter_than_k_control", "A", k, k+1, null, direct);
			}
		}
		if(verify && failures>0){throw new AssertionError("JT012 oracle mismatches: "+failures);}
	}

	/** Independently enumerate every wholly contained k-mer; never call product counting/score methods. */
	private static double oracle(String sequence, int k, int window){
		assert(window>=k && k>=2 && k<=3) : "Finite JT012 oracle enumerates only valid k2/k3 windows.";
		final int width=Math.min(window, sequence.length());
		final int windows=Math.max(1, sequence.length()-width+1);
		double sum=0;
		for(int start=0; start<windows; start++){
			final int[] counts=new int[1<<(2*k)];
			for(int pos=start; pos+k<=start+width; pos++){
				int code=0;
				boolean valid=true;
				for(int j=0; j<k; j++){
					final int base="ACGT".indexOf(sequence.charAt(pos+j));
					if(base<0){valid=false; break;}
					code=4*code+base;
				}
				if(valid){counts[code]++;}
			}
			long lower=0, upper=0;
			for(int code=0; code<counts.length/2; code++){
				//Basewise complement, intentionally not reverse-complement.
				final int a=counts[code], b=counts[counts.length-1-code];
				lower+=Math.min(a, b);
				upper+=Math.max(a, b);
			}
			sum+=lower/(double)Math.max(1, upper);
		}
		return sum/windows;
	}

	private static void checkLiteral(String sequence, int k, int window, double expected){
		final double actual=oracle(sequence, k, window);
		assert(Math.abs(actual-expected)<1e-7) : "Independent oracle disagrees with hand-counted fixture: "+
			sequence+", k="+k+", window="+window+", expected="+expected+", actual="+actual;
	}

	private static float call(String sequence, int k, int window, int[] scratch, boolean direct){
		final byte[] bases=sequence.getBytes(StandardCharsets.US_ASCII);
		return direct ? EntropyTracker.strandednessWindowedK2(bases, scratch, window) :
			EntropyTracker.strandednessWindowed(bases, scratch, k, window);
	}

	private static void observe(String name, String sequence, int k, int window, int[] scratch, boolean direct){
		final double expected=oracle(sequence, k, window);
		final float observed=call(sequence, k, window, scratch, direct);
		final boolean agrees=Math.abs(expected-observed)<1e-7;
		System.out.println(name+"\t"+k+"\t"+window+"\t"+sequence+"\t"+expected+"\t"+observed+
			"\t"+agrees);
		if(!agrees){failures++;}
	}
	private static int failures=0;
}
