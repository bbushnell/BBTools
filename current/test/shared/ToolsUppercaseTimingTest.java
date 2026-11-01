package test.shared;

import shared.Tools;

/** Allocation-free timing of the scalar byte-uppercase helper on valid ASCII data.
 * The independent checksum uses the ASCII alphabet bounds, not the helper table.
 * @author Shinobu
 */
public final class ToolsUppercaseTimingTest{

	/** Runs one prebuilt corpus; rounds and warmup are positive repetition counts. */
	public static void main(String[] args){
		long rounds=100000, warmup=10000;
		String mode="ascii";
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			if(split.length!=2){throw new IllegalArgumentException("Expected key=value: "+arg);}
			if(split[0].equals("rounds")){rounds=Long.parseLong(split[1]);}
			else if(split[0].equals("warmup")){warmup=Long.parseLong(split[1]);}
			else if(split[0].equals("mode")){mode=split[1];}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(rounds<1 || warmup<1 || (!mode.equals("ascii") && !mode.equals("dna"))){
			throw new IllegalArgumentException("Positive rounds/warmup and mode=ascii|dna required");
		}
		final byte[] data=new byte[4096];
		final String dna="ACGTNacgtn";
		long oneRound=0;
		for(int i=0; i<data.length; i++){
			final int value=mode.equals("ascii") ? ((i*37+11)&127) : dna.charAt((i*7+3)%dna.length());
			data[i]=(byte)value;
			oneRound+=(value>='a' && value<='z' ? value-32 : value);
		}
		if(rounds>Long.MAX_VALUE/oneRound || warmup>Long.MAX_VALUE/oneRound){
			throw new IllegalArgumentException("Repetition count would overflow checksum");
		}
		require(loop(data, warmup)==oneRound*warmup, "warmup checksum");
		final long start=System.nanoTime();
		final long sum=loop(data, rounds);
		final long elapsed=System.nanoTime()-start;
		require(sum==oneRound*rounds, "measured checksum");
		System.out.println("TOOLS_UPPERCASE_TIMING mode="+mode+" rounds="+rounds+" warmup="+warmup+
				" bytes="+(rounds*data.length)+" nanos="+elapsed+" checksum="+sum);
	}

	/** Traverses the same preallocated input with no allocation or per-element dispatch. */
	private static long loop(byte[] data, long rounds){
		long sum=0;
		for(long r=0; r<rounds; r++){
			for(byte value : data){sum+=Tools.toUpperCase(value);}
		}
		return sum;
	}

	/** Verifies the oracle even when assertions are disabled. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
