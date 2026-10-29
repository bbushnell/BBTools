package test.tax;

import java.nio.charset.StandardCharsets;

import parse.Parse;
import tax.GiToTaxid;
import tax.GiToTaxidInt;

/** Measures explicit-ID header lookup without allocating in the timed loop.
 * Corpus construction and literal result checks precede warmup and measurement.
 * No GI table is needed because every header contains an explicit TaxID.
 * @author Shinobu */
public final class GiHeaderBenchmark{

	/** Measures one String or byte-array lookup implementation.
	 * @param args kind=array-string|array-bytes|dense-string|dense-bytes,
	 *        batch=65536 and rounds=16384 */
	public static void main(String[] args){
		int batch=65536, rounds=16384;
		String kind="array-string";
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			if(split.length!=2){throw new IllegalArgumentException("Expected key=value: "+arg);}
			if(split[0].equals("batch")){batch=Parse.parseIntKMG(split[1]);}
			else if(split[0].equals("rounds")){rounds=Parse.parseIntKMG(split[1]);}
			else if(split[0].equals("kind")){kind=split[1];}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		if(batch<1 || rounds<1){throw new IllegalArgumentException("batch and rounds must be positive");}
		final int mode=mode(kind);
		final String[] strings=new String[batch];
		final byte[][] bytes=new byte[batch][];
		long expected=0;
		for(int i=0; i<batch; i++){
			final int id=1+(int)(((long)i*31)%2000000);
			switch(i&3){
				case 0: strings[i]="tid|"+id+"|sequence"; break;
				case 1: strings[i]="ncbi|"+id+"|sequence"; break;
				case 2: strings[i]="tid_"+id+"_sequence"; break;
				default: strings[i]="ncbi_"+id+"_sequence"; break;
			}
			bytes[i]=strings[i].getBytes(StandardCharsets.US_ASCII);
			if(GiToTaxid.getID(strings[i])!=id || GiToTaxid.getID(bytes[i])!=id ||
					GiToTaxidInt.getID(strings[i])!=id || GiToTaxidInt.getID(bytes[i])!=id){
				throw new AssertionError("Explicit-ID oracle mismatch for "+strings[i]+", expected "+id);
			}
			expected+=id;
		}
		if(run(strings, bytes, 8, mode)!=expected*8){throw new AssertionError("Warmup checksum mismatch");}
		final long start=System.nanoTime();
		final long sum=run(strings, bytes, rounds, mode);
		final long elapsed=System.nanoTime()-start;
		if(sum!=expected*rounds){throw new AssertionError("Timed checksum mismatch");}
		System.out.println("kind="+kind+" checksum="+sum+" expected="+(expected*rounds)+
				" batch="+batch+" rounds="+rounds+" compute_ns="+elapsed);
	}

	/** Prevents construction of this standalone benchmark. */
	private GiHeaderBenchmark(){}

	/** Converts the CLI choice into a dispatch code outside the timed loop.
	 * @param kind Requested implementation
	 * @return Code from zero through three
	 * @throws IllegalArgumentException For an unknown implementation */
	private static int mode(final String kind){
		if(kind.equals("array-string")){return 0;}
		if(kind.equals("array-bytes")){return 1;}
		if(kind.equals("dense-string")){return 2;}
		if(kind.equals("dense-bytes")){return 3;}
		throw new IllegalArgumentException("Unknown kind: "+kind);
	}

	/** Parses prebuilt headers; dispatch occurs once per invocation, outside loops.
	 * @param strings Explicit-ID String corpus
	 * @param bytes Matching ASCII byte corpus
	 * @param rounds Number of corpus passes
	 * @param mode Implementation code from {@link #mode(String)}
	 * @return Sum of all returned TaxIDs */
	private static long run(final String[] strings, final byte[][] bytes,
			final int rounds, final int mode){
		assert(strings.length==bytes.length && rounds>0 && mode>=0 && mode<4) :
				"The benchmark dispatch requires matching corpora, positive rounds and a valid mode.";
		long sum=0;
		switch(mode){
			case 0:
				for(int r=0; r<rounds; r++){for(String s : strings){sum+=GiToTaxid.getID(s);}}
				break;
			case 1:
				for(int r=0; r<rounds; r++){for(byte[] b : bytes){sum+=GiToTaxid.getID(b);}}
				break;
			case 2:
				for(int r=0; r<rounds; r++){for(String s : strings){sum+=GiToTaxidInt.getID(s);}}
				break;
			case 3:
				for(int r=0; r<rounds; r++){for(byte[] b : bytes){sum+=GiToTaxidInt.getID(b);}}
				break;
			default: throw new IllegalArgumentException("Unknown benchmark mode: "+mode);
		}
		return sum;
	}
}
