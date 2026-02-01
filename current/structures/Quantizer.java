package structures;

import java.util.ArrayList;

import parse.Parse;
import shared.KillSwitch;
import shared.Tools;
import stream.Read;

/**
 * Quality score quantization utility for reducing quality score granularity.
 * Provides configurable quantization of Phred quality scores to reduce storage requirements
 * and noise in downstream analysis. Supports both uniform interval and custom bin quantization
 * with optional "sticky" mode that preserves adjacent identical scores.
 *
 * @author Brian Bushnell
 * @date 2025
 */
public class Quantizer {
	
	/**
	 * Parses command-line arguments related to quantization settings.
	 * Handles quantize, quantizesticky parameters and delegates to appropriate handlers.
	 *
	 * @param arg The full argument string
	 * @param a The parameter name (e.g., "quantize", "quantizesticky")
	 * @param b The parameter value
	 * @return true if argument was recognized and processed, false otherwise
	 */
	public static boolean parse(String arg, String a, String b){
		if(a.equals("quantize")){
			if(b!=null && b.equalsIgnoreCase("sticky")){
				STICKY=true;
				return true;
			}
		}else if(a.equals("quantizesticky")){
			if(b!=null && (b.charAt(0)=='/' || Tools.isDigit(b.charAt(0)))) {
				STICKY=true;
			}else {
				STICKY=Parse.parseBoolean(b);
				return true;
			}
		}
		
		if(b==null || b.length()<1 || Character.isLetter(b.charAt(0))){
			return Parse.parseBoolean(b);
		}
		return setArray(b);
	}
	
	/**
	 * Configures quantization array from string specification.
	 * Supports both uniform interval format ("/N") and comma-separated bin list.
	 * For uniform intervals, creates bins starting at 0 with step size N.
	 *
	 * @param s String specification: "/N" for uniform intervals or "v1,v2,v3" for custom bins
	 * @return true if array was successfully configured, false if disabled (quant=1)
	 */
	private static boolean setArray(String s){
		final byte[] array;
		if(s.indexOf(',')<0){
			if(s.charAt(0)=='/') {s=s.substring(1);}
			int quant=Integer.parseInt(s);
			assert(quant>0 && quant<128);
			if(quant==1){return false;}
			ByteBuilder bb=new ByteBuilder();
			for(int i=0, max=Read.MAX_CALLED_QUALITY(); i<=max; i+=quant){
				bb.append((byte)i);
			}
			array=bb.toBytes();
		}else{
			array=Parse.parseByteArray(s, ",");
		}
		setArray(array);
		return true;
	}
	
	/** Sets the quantization array and rebuilds the quality remapping table.
	 * @param a Array of quantization bin values */
	private static void setArray(byte[] a){
		quantizeArray=a;
		qualityRemapArray=makeQualityRemapArray(quantizeArray);
	}
	
	/**
	 * Quantizes quality scores for all reads in a list.
	 * Processes both primary reads and their mates if present.
	 * @param list List of reads to quantize (may be null)
	 */
	public static void quantize(ArrayList<Read> list) {
		if(list==null) {return;}
		for(Read r : list) {
			if(r!=null) {
				quantize(r);
				quantize(r.mate);
			}
		}
	}
	
	/**
	 * Quantizes quality scores for a pair of reads.
	 * @param r1 First read to quantize (may be null)
	 * @param r2 Second read to quantize (may be null)
	 */
	public static void quantize(Read r1, Read r2){
		quantize(r1);
		quantize(r2);
	}
	
	/** Quantizes quality scores for a single read.
	 * @param r Read to quantize (may be null) */
	public static void quantize(Read r){
		if(r!=null) {quantize(r.quality);}
	}
	
	/**
	 * Quantizes a quality score array in-place using configured mapping.
	 * In sticky mode, adjacent positions with similar original scores may retain
	 * the same quantized value to preserve local quality patterns.
	 * @param quals Quality score array to quantize (may be null)
	 */
	public static void quantize(byte[] quals){
		if(quals==null){return;}
		byte prev=0;
		for(int i=0; i<quals.length; i++){
			final byte qOld=quals[i];
			final byte q0=qualityRemapArray[qOld];
			byte q=q0;
			if(STICKY && q!=prev && prev>0 && q>0 && Tools.absdif(qOld, prev)<=Tools.absdif(qOld, q)){q=prev;}
//			assert(q==q0) : STICKY+", "+qOld+" -> "+q0+" -> "+q+", prev="+prev;
			quals[i]=q;
			prev=q;
		}
	}
	
	/**
	 * Creates a lookup table mapping original quality scores to quantized values.
	 * For each possible quality score (0-127), finds the closest quantization bin value.
	 * @param quantizeArray Array of quantization bin values
	 * @return Remapping array where index is original quality, value is quantized quality
	 */
	private static final byte[] makeQualityRemapArray(byte[] quantizeArray) {
		byte[] array=KillSwitch.allocByte1D(128);
		for(int i=0; i<array.length; i++){
			byte q=0;
			for(byte x : quantizeArray){
				if((i>0 && q==0 && x>0) || Tools.absdif(x, i)<=Tools.absdif(q, i)){q=x;}
			}
			array[i]=q;
		}
		return array;
	}
	
//	private static byte[] quantizeArray={0, 8, 13, 22, 27, 32, 37}; //Old
	private static byte[] quantizeArray={0, 14, 21, 27, 32, 36};
	/** Lookup table mapping original quality scores to quantized values */
	private static byte[] qualityRemapArray=makeQualityRemapArray(quantizeArray);
	/**
	 * Whether to use sticky mode that preserves adjacent identical quantized scores
	 */
	private static boolean STICKY=true;
	
}
