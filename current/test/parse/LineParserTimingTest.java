package test.parse;

import java.nio.charset.StandardCharsets;

import parse.LineParser1;
import parse.LineParserS2;
import parse.LineParserS3;
import shared.Shared;
import simd.Vector;

/** Single-threaded, in-memory timing driver with independently checked checksums.
 * Inputs and parsers are allocated before warmup; output is emitted after timing.
 * The team recipe must calibrate rounds to at least ten seconds of compute and
 * retain repeated trials, host/JDK/heap, wall time and CPU time for both revisions.
 * Short runs are pilots only. ASCII inputs give identical results before repairs.
 * @author Brian Bushnell, Sucrose */
public final class LineParserTimingTest {

	/** Parses flag=value options and runs a selected workload.
	 * @param args mode=lp1|s2char|s3, rounds=100000, warmup=10000,
	 * length=256, spacing=16, simd=f */
	public static void main(String[] args){
		String mode="lp1";
		long rounds=100000, warmup=10000;
		int length=256, spacing=16;
		boolean simd=false;
		for(String arg : args){
			final int equals=arg.indexOf('=');
			if(equals<1){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String key=arg.substring(0, equals), value=arg.substring(equals+1);
			if(key.equalsIgnoreCase("mode")){mode=value;}
			else if(key.equalsIgnoreCase("rounds")){rounds=Long.parseLong(value);}
			else if(key.equalsIgnoreCase("warmup")){warmup=Long.parseLong(value);}
			else if(key.equalsIgnoreCase("length")){length=Integer.parseInt(value);}
			else if(key.equalsIgnoreCase("spacing")){spacing=Integer.parseInt(value);}
			else if(key.equalsIgnoreCase("simd")){
				if(value.equalsIgnoreCase("t")){simd=true;}
				else if(value.equalsIgnoreCase("f")){simd=false;}
				else{throw new IllegalArgumentException("simd expects t or f");}
			}else{throw new IllegalArgumentException("Unknown flag: "+key);}
		}
		if(!mode.equals("lp1") && !mode.equals("s2char") && !mode.equals("s3")){
			throw new IllegalArgumentException("Unknown mode: "+mode);
		}
		if(length<8 || length>1048576 || spacing<2 || spacing>length-2){
			throw new IllegalArgumentException("Require length in [8,1048576], spacing in [2,length-2]");
		}
		if(rounds<1 || warmup<1 || rounds>Long.MAX_VALUE/64/length || warmup>Long.MAX_VALUE/64/length){
			throw new IllegalArgumentException("Positive rounds/warmup must fit processed-byte counters");
		}
		if(simd && !Vector.simd128){throw new IllegalStateException("SIMD requested but unavailable");}
		final boolean saved=Shared.SIMD;
		Shared.SIMD=simd;
		try{
			final LineParserTimingTest test=new LineParserTimingTest(length, spacing);
			test.verify(mode, warmup, test.measure(mode, warmup));
			final long start=System.nanoTime();
			final long checksum=test.measure(mode, rounds);
			final long nanos=System.nanoTime()-start;
			test.verify(mode, rounds, checksum);
			System.out.println("LINE_PARSER_TIMING mode="+mode+" rounds="+rounds+
				" warmup="+warmup+" length="+length+" spacing="+spacing+" simd="+simd+
				" vectorBits="+Vector.SIMD_BITS+" lines="+(rounds*64)+
				" inputBytes="+(rounds*64*length)+" nanos="+nanos+" checksum="+checksum);
		}finally{Shared.SIMD=saved;}
	}

	/** Builds 64 reusable ASCII lines and computes expectations without a production parser. */
	private LineParserTimingTest(int length, int spacing){
		for(int row=0; row<bytes.length; row++){
			final byte[] data=new byte[length];
			for(int i=0; i<length; i++){data[i]=(byte)('a'+(row+i)%26);}
			for(int i=spacing-1; i<length-1; i+=spacing){data[i]=',';}
			data[length-2]=' ';
			bytes[row]=data;
			strings[row]=new String(data, StandardCharsets.US_ASCII);
			int commas=0, last=-1, first=-1, whitespaceIndex=0;
			for(int i=0; i<length; i++){
				if(data[i]==','){
					commas++;
					last=i;
					if(first<0){first=i;}
					if(i<length-2){whitespaceIndex++;}
				}
			}
			assert(first>=0 && first+1<length) : "spacing bounds guarantee a nonempty second comma field";
			expectedLp1+=commas+1+length-last-1;
			expectedS2+=data[first+1];
			expectedS3+=commas+2+whitespaceIndex;
		}
	}

	/** Reuses parser instances; the selected loop only parses and accumulates observable results. */
	private long measure(String mode, long rounds){
		long checksum=0;
		if(mode.equals("lp1")){
			for(long r=0; r<rounds; r++){
				for(byte[] line : bytes){
					lp1.set(line);
					final int count=lp1.terms();
					checksum+=count+lp1.length(count-1);
				}
			}
		}else if(mode.equals("s2char")){
			for(long r=0; r<rounds; r++){
				for(String line : strings){s2.set(line); checksum+=s2.parseChar(1, 0);}
			}
		}else{
			for(long r=0; r<rounds; r++){
				for(String line : strings){s3.set(line); checksum+=s3.terms()+s3.indexOfWhitespace();}
			}
		}
		return checksum;
	}

	/** Checks the measured checksum against the independent corpus oracle, even under -da.
	 * Java long wraparound is identical on both sides for very large round counts. */
	private void verify(String mode, long rounds, long actual){
		final long expected=rounds*(mode.equals("lp1") ? expectedLp1 : mode.equals("s2char") ? expectedS2 : expectedS3);
		if(actual!=expected){throw new AssertionError("Checksum mismatch for "+mode+": "+actual+" != "+expected);}
	}

	/** Prebuilt corpus; UTF-16 string and ASCII byte versions contain identical characters. */
	private final byte[][] bytes=new byte[64][];
	private final String[] strings=new String[64];
	/** Reused parsers, owned exclusively by this single-threaded driver. */
	private final LineParser1 lp1=new LineParser1(',');
	private final LineParserS2 s2=new LineParserS2(',');
	private final LineParserS3 s3=new LineParserS3(',');
	/** Expected checksums for one pass through all 64 lines. */
	private long expectedLp1, expectedS2, expectedS3;
}
