package test.simd;

import java.util.Arrays;

import shared.Shared;
import simd.Vector;
import structures.IntList;

/** Checks exact line and FASTA boundary scans through the Vector byte facade.
 * Cases span empty ranges, nonzero starts, scalar tails and 64/256-bit chunks.
 * @author Sucrose */
public final class VectorByteScanTest{

	/** Runs finite scan fixtures.
	 * @param args simd=t or simd=f
	 */
	public static void main(final String[] args){
		final boolean simd=parseSimd(args);
		if(simd && !Vector.simd128){
			throw new IllegalStateException("SIMD test requires 128/256-bit Vector support");
		}
		final boolean oldSimd=Shared.SIMD;
		int cases=0, failures=0;
		try{
			Shared.SIMD=simd;
			for(final int length : LENGTHS){
				for(final int from : starts(length)){
					for(final int to : ends(length, from)){
						failures+=checkSymbolScan(makeSymbolBuffer(length, from, to), from, to); cases++;
						failures+=checkFastaScan(makeFastaBuffer(length, from, to), from, to); cases++;
					}
				}
				for(final int limit : limits(length)){
					failures+=checkLastScan(makeLastBuffer(length), limit); cases++;
				}
			}
			System.out.println("BYTE_SCAN_RESULT simd="+simd+" cases="+cases+" failures="+failures);
			if(failures!=0){throw new AssertionError("Byte scan failures="+failures);}
		}finally{
			Shared.SIMD=oldSimd;
		}
	}

	/** Parses the SIMD selector.
	 * @param args Command-line arguments
	 * @return True for simd=t
	 */
	private static boolean parseSimd(final String[] args){
		boolean simd=false;
		for(final String arg : args){
			if(arg.equals("simd=t")){simd=true;}
			else if(arg.equals("simd=f")){simd=false;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		return simd;
	}

	/** Builds a symbol-scan buffer with in-range and just-outside-range hits.
	 * @param length Buffer length
	 * @param from Inclusive start
	 * @param to Exclusive end
	 * @return New buffer
	 */
	private static byte[] makeSymbolBuffer(final int length, final int from, final int to){
		final byte[] buffer=new byte[length];
		Arrays.fill(buffer, (byte)'A');
		if(from<to){buffer[from]='\n';}
		if(from+7<to){buffer[from+7]='\n';}
		if(from+15<to){buffer[from+15]='\n';}
		if(from+31<to){buffer[from+31]='\n';}
		if(to>from){buffer[to-1]='\n';}
		if(to<length){buffer[to]='\n';}
		return buffer;
	}

	/** Builds a FASTA-boundary buffer.
	 * @param length Buffer length
	 * @param from Inclusive start
	 * @param to Exclusive end
	 * @return New buffer
	 */
	private static byte[] makeFastaBuffer(final int length, final int from, final int to){
		final byte[] buffer=new byte[length];
		Arrays.fill(buffer, (byte)'A');
		placeBoundary(buffer, from);
		placeBoundary(buffer, from+7);
		placeBoundary(buffer, from+15);
		placeBoundary(buffer, from+31);
		placeBoundary(buffer, to-2);
		if(to<length){buffer[to]='>';}
		if(length>0){buffer[0]='>';}
		return buffer;
	}

	/** Builds a buffer with hits in scalar tail and vector chunks.
	 * @param length Buffer length
	 * @return New buffer
	 */
	private static byte[] makeLastBuffer(final int length){
		final byte[] buffer=new byte[length];
		Arrays.fill(buffer, (byte)'A');
		if(length>0){buffer[0]='\n';}
		if(length>8){buffer[8]='\n';}
		if(length>16){buffer[16]='\n';}
		if(length>32){buffer[32]='\n';}
		if(length>0){buffer[length-1]='\n';}
		return buffer;
	}

	/** Places one newline/greater-than boundary if it fits.
	 * @param buffer Target buffer
	 * @param pos Newline position
	 */
	private static void placeBoundary(final byte[] buffer, final int pos){
		if(pos>=0 && pos+1<buffer.length){
			buffer[pos]='\n';
			buffer[pos+1]='>';
		}
	}

	/** Checks findSymbols and countSymbols.
	 * @param buffer Input buffer
	 * @param from Inclusive start
	 * @param to Exclusive end
	 * @return One on mismatch, otherwise zero
	 */
	private static int checkSymbolScan(final byte[] buffer, final int from, final int to){
		final IntList expected=new IntList(), actual=new IntList();
		expected.add(SENTINEL);
		actual.add(SENTINEL);
		for(int i=from; i<to; i++){
			if(buffer[i]=='\n'){expected.add(i);}
		}
		final int returned=Vector.findSymbols(buffer, from, to, (byte)'\n', actual);
		final int counted=Vector.countSymbols(buffer, from, to, (byte)'\n');
		if(returned==expected.size() && counted==expected.size()-1 && equals(expected, actual)){
			return 0;
		}
		System.out.println("FAIL symbol length="+buffer.length+" from="+from+" to="+to+
				" expected="+listString(expected)+" actual="+listString(actual)+
				" returned="+returned+" counted="+counted);
		return 1;
	}

	/** Checks findFastaHeaders.
	 * @param buffer Input buffer
	 * @param from Inclusive start
	 * @param to Exclusive end
	 * @return One on mismatch, otherwise zero
	 */
	private static int checkFastaScan(final byte[] buffer, final int from, final int to){
		final IntList expected=new IntList(), actual=new IntList();
		expected.add(SENTINEL);
		actual.add(SENTINEL);
		for(int i=from+1; i<to; i++){
			if(buffer[i]=='>' && buffer[i-1]=='\n'){expected.add(i-1);}
		}
		final int returned=Vector.findFastaHeaders(buffer, from, to, actual);
		if(returned==expected.size() && equals(expected, actual)){return 0;}
		System.out.println("FAIL fasta length="+buffer.length+" from="+from+" to="+to+
				" expected="+listString(expected)+" actual="+listString(actual)+
				" returned="+returned);
		return 1;
	}

	/** Checks findLastSymbol.
	 * @param buffer Input buffer
	 * @param limit Exclusive scan limit
	 * @return One on mismatch, otherwise zero
	 */
	private static int checkLastScan(final byte[] buffer, final int limit){
		int expected=-1;
		for(int i=limit-1; i>=0 && expected<0; i--){
			if(buffer[i]=='\n'){expected=i;}
		}
		final int actual=Vector.findLastSymbol(buffer, limit, (byte)'\n');
		if(actual==expected){return 0;}
		System.out.println("FAIL last length="+buffer.length+" limit="+limit+
				" expected="+expected+" actual="+actual);
		return 1;
	}

	/** Returns selected starts for one length.
	 * @param length Buffer length
	 * @return Start positions
	 */
	private static int[] starts(final int length){
		return new int[] {0, Math.min(3, length), Math.min(9, length)};
	}

	/** Returns selected exclusive ends for one length and start.
	 * @param length Buffer length
	 * @param from Inclusive start
	 * @return End positions
	 */
	private static int[] ends(final int length, final int from){
		return new int[] {from, Math.max(from, length-1), length};
	}

	/** Returns selected findLastSymbol limits.
	 * @param length Buffer length
	 * @return Exclusive limits
	 */
	private static int[] limits(final int length){
		return new int[] {0, Math.min(1, length), Math.min(7, length), Math.min(8, length),
				Math.min(9, length), Math.min(15, length), Math.min(16, length),
				Math.min(17, length), Math.min(31, length), Math.min(32, length),
				Math.min(33, length), length};
	}

	/** Compares IntList live prefixes.
	 * @param a Expected list
	 * @param b Actual list
	 * @return True when equal
	 */
	private static boolean equals(final IntList a, final IntList b){
		if(a.size()!=b.size()){return false;}
		for(int i=0; i<a.size(); i++){
			if(a.get(i)!=b.get(i)){return false;}
		}
		return true;
	}

	/** Formats an IntList live prefix.
	 * @param list List to format
	 * @return Text form
	 */
	private static String listString(final IntList list){
		final int[] array=new int[list.size()];
		for(int i=0; i<array.length; i++){array[i]=list.get(i);}
		return Arrays.toString(array);
	}

	/** Prevents construction of this standalone test. */
	private VectorByteScanTest(){}

	/** Lengths spanning 8-byte, 16-byte, 32-byte and tail thresholds. */
	private static final int[] LENGTHS={0, 1, 7, 8, 9, 15, 16, 17, 23, 24, 25, 31,
			32, 33, 39, 40, 41, 47, 48, 49, 65, 79, 80, 81};
	/** Sentinel proving scan methods append without clearing caller state. */
	private static final int SENTINEL=-1234567;
}
