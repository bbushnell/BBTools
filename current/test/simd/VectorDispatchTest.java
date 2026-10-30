package test.simd;

import parse.Parser;
import shared.Shared;
import simd.Vector;
import structures.ByteBuilder;

/** Checks small exact integer results and records SIMD initialization order.
 * Runs with or without the optional Vector module; this is not a throughput test
 * or a substitute for execution on an actual Java8 runtime.
 * @author Shinobu */
public final class VectorDispatchTest{

	/** Tests default or forced-scalar dispatch in a fresh JVM.
	 * @param args load=shared|vector, module=present|absent, scalar=t|f, optional handles=t|f
	 * @throws Exception If the requested first class cannot be initialized */
	public static void main(final String[] args) throws Exception{
		String load="shared", module="present";
		boolean scalar=false;
		Boolean expectedHandles=null;
		for(String arg : args){
			if(arg.equals("load=shared")){load="shared";}
			else if(arg.equals("load=vector")){load="vector";}
			else if(arg.equals("module=present")){module="present";}
			else if(arg.equals("module=absent")){module="absent";}
			else if(arg.equals("scalar=t")){scalar=true;}
			else if(arg.equals("scalar=f")){scalar=false;}
			else if(arg.equals("handles=t")){expectedHandles=Boolean.TRUE;}
			else if(arg.equals("handles=f")){expectedHandles=Boolean.FALSE;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		Class.forName(load.equals("shared") ? "shared.Shared" : "simd.Vector");
		final boolean oldSimd=Shared.SIMD;
		try{
			final boolean present=modulePresent();
			check(present==module.equals("present"), "Unexpected Vector module availability: "+present);
			check(oldSimd==Vector.simd128,
					"Shared.SIMD default must follow 128-bit Vector availability; default="+
					oldSimd+" simd128="+Vector.simd128);
			if(!present){
				check(!oldSimd && Vector.maxSimdWidth==0,
						"Absent Vector module must leave default dispatch scalar: SIMD="+
						oldSimd+" width="+Vector.maxSimdWidth);
			}
			checkSimdParserRoundTrip(oldSimd);
			if(scalar){Shared.SIMD=false;}
			final int[] values=new int[65];
			for(int i=0; i<values.length; i++){values[i]=i+1;}
			final long sum=Vector.sum(values), range=Vector.sum(values, 7, 31);
			final int max=Vector.max(values);
			check(sum==2145 && range==500 && max==65,
					"Expected sum2145/range500/max65, got "+sum+"/"+range+"/"+max);
			checkStringAppend();
			System.out.println("VECTOR_DISPATCH_RESULT load="+load+" module="+module+
					" scalar="+scalar+" default_simd="+oldSimd+" active_simd="+Shared.SIMD+
					" width="+Vector.maxSimdWidth+" simd_bits="+Vector.SIMD_BITS+
					" int_lanes="+Vector.INT_LANES+" var_handles="+Vector.varHandles+
					" sum="+sum+" range="+range+" max="+max);
			final int expectedSimdBits=(Vector.simd256 ? 256 : (Vector.simd128 ? 128 : 0));
			check(Vector.SIMD_BITS==expectedSimdBits,
					"SIMD_BITS must report the selected public width; width="+
					Vector.maxSimdWidth+" expected="+expectedSimdBits+
					" actual="+Vector.SIMD_BITS);
			final int expectedIntLanes=(Vector.simd256 ? 8 : (Vector.simd128 ? 4 : 1));
			check(Vector.INT_LANES==expectedIntLanes,
					"Integer lanes must follow the selected vector width; width="+
					Vector.maxSimdWidth+" expected="+expectedIntLanes+
					" actual="+Vector.INT_LANES);
			check(Vector.simd256==(Vector.maxSimdWidth>=256) &&
					Vector.simd128==(Vector.maxSimdWidth>=128) &&
					Vector.simd64==(Vector.maxSimdWidth>=64),
					"SIMD capability flags must follow maxSimdWidth="+Vector.maxSimdWidth+
					" flags="+Vector.simd256+"/"+Vector.simd128+"/"+Vector.simd64);
			if(expectedHandles!=null){
				check(Vector.varHandles==expectedHandles.booleanValue(),
						"Unexpected optional String access: expected="+expectedHandles+" actual="+Vector.varHandles);
			}
		}finally{
			Shared.SIMD=oldSimd;
		}
	}

	/** Checks that the common parser can force SIMD and restore width-based auto.
	 * @param oldSimd Startup SIMD flag to restore after parser mutation */
	private static void checkSimdParserRoundTrip(final boolean oldSimd){
		try{
			check(Parser.parseStatic("simd=f", "simd", "f"), "Parser rejected simd=f");
			check(!Shared.SIMD, "simd=f must disable Shared.SIMD");
			check(Parser.parseStatic("simd=t", "simd", "t"), "Parser rejected simd=t");
			check(Shared.SIMD, "simd=t must enable Shared.SIMD");
			check(Parser.parseStatic("simd=auto", "simd", "auto"), "Parser rejected simd=auto");
			check(Shared.SIMD==Vector.simd128,
					"simd=auto must follow 128-bit Vector availability; SIMD="+
					Shared.SIMD+" simd128="+Vector.simd128);
		}finally{
			Shared.SIMD=oldSimd;
		}
	}

	/** Checks growth, Latin-1 copying and low-byte UTF-16 conversion. */
	private static void checkStringAppend(){
		final ByteBuilder bb=new ByteBuilder(1);
		bb.append('!');
		Vector.append(bb, "A\u00e9");
		Vector.append(bb, "\u0100B");
		check(bb.length==5 && bb.array[0]=='!' && bb.array[1]=='A' &&
				bb.array[2]==(byte)233 && bb.array[3]==0 && bb.array[4]=='B',
				"String append must preserve the prefix and the low byte of every UTF-16 code unit");
	}

	/** Reports whether the current JVM can load the optional Vector API.
	 * @return True when ByteVector is available */
	private static boolean modulePresent(){
		try{
			Class.forName("jdk.incubator.vector.ByteVector");
			return true;
		}catch(ClassNotFoundException e){return false;}
	}

	/** Keeps correctness checks active independently of JVM assertion flags.
	 * @param condition Required condition
	 * @param message Failure explanation */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Prevents construction of this standalone test. */
	private VectorDispatchTest(){}
}
