package structures;

import java.util.ArrayList;

import parse.Parse;
import shared.KillSwitch;
import shared.Tools;
import stream.Read;

/**
 * Static, in-place remapping of numeric quality scores through configurable bins.
 * Configuration is shared and unsynchronized; configure before concurrent use.
 * Sticky mode can reuse the preceding positive output when it is at least as
 * close to the original score as the selected positive bin. It starts afresh
 * for each quality array. Input quality bytes must be nonnegative table indices.
 *
 * @author Brian Bushnell
 * @date 2025
 */
public class Quantizer{

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Applies a quantize or quantizesticky option already selected by the caller.
	 * A boolean sticky setting returns true to enable quantization independently
	 * of the sticky value; ordinary boolean values return their parsed value.
	 * Numeric specifications replace bins, except interval 1 returns false without
	 * replacing them. This method does not validate an arbitrary option name.
	 * @param arg Full argument string; unused here.
	 * @param a Preselected option name, normally lower case.
	 * @param b Boolean, sticky keyword, interval, or comma-separated bin values.
	 * @return Whether the caller should enable quantization, not name recognition.
	 */
	public static boolean parse(String arg, String a, String b){
		if(a.equals("quantize")){
			if(b!=null && b.equalsIgnoreCase("sticky")){
				STICKY=true;
				return true;
			}
		}else if(a.equals("quantizesticky")){
			if(b!=null && (b.charAt(0)=='/' || Tools.isDigit(b.charAt(0)))){
				STICKY=true;
			}else{
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
	 * A single integer, optionally prefixed by '/', creates bins from 0 through
	 * Read.MAX_CALLED_QUALITY() in that interval. Comma-separated signed byte
	 * values are used in their supplied order; this method does not sort them.
	 *
	 * @param s String specification: "/N" for uniform intervals or "v1,v2,v3" for custom bins
	 * @return True after configuring bins, or false for interval 1 without changing bins.
	 */
	private static boolean setArray(String s){
		final byte[] array;
		if(s.indexOf(',')<0){
			if(s.charAt(0)=='/'){s=s.substring(1);}
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

	/** Retains the supplied bins and rebuilds the shared remapping table.
	 * @param a Bin values in tie-breaking iteration order; borrowed, not copied. */
	private static void setArray(byte[] a){
		quantizeArray=a;
		qualityRemapArray=makeQualityRemapArray(quantizeArray);
	}

	/**
	 * Quantizes quality scores for all reads in a list.
	 * Processes each nonnull list entry and its mate, in list order.
	 * The list is not deduplicated; listing both mates quantizes each twice.
	 * @param list List of reads to quantize (may be null)
	 */
	public static void quantize(ArrayList<Read> list){
		if(list==null){return;}
		for(Read r : list){
			if(r!=null){
				quantize(r);
				quantize(r.mate);
			}
		}
	}

	/**
	 * Quantizes two supplied reads independently; identical references are processed twice.
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
		if(r!=null){quantize(r.quality);}
	}

	/**
	 * Quantizes a quality score array in-place using configured mapping.
	 * Each original byte indexes the 128-entry table. Sticky mode retains the
	 * previous positive output if both candidates are positive and the previous
	 * output is no farther from the original score. Previous output starts at zero.
	 * @param quals Borrowed quality array to modify, with values 0..127; may be null.
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
	 * For each score 0..127, starts at zero and visits bins in supplied order.
	 * A bin replaces the candidate when it is at least as close, so later ties win.
	 * For a positive score, a positive bin also replaces a zero candidate regardless
	 * of distance. Zero input remains zero even if no zero bin was supplied.
	 * @param quantizeArray Bin values in iteration order.
	 * @return Remapping array where index is original quality, value is quantized quality
	 */
	private static final byte[] makeQualityRemapArray(byte[] quantizeArray){
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

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

//	private static byte[] quantizeArray={0, 8, 13, 22, 27, 32, 37}; //Old
	/** Current bins; replacement preserves the supplied iteration order. */
	private static byte[] quantizeArray={0, 14, 21, 27, 32, 36};
	/** Lookup table mapping original quality scores to quantized values */
	private static byte[] qualityRemapArray=makeQualityRemapArray(quantizeArray);
	/**
	 * Whether to reuse the preceding positive output when it is at least as close.
	 */
	private static boolean STICKY=true;

}
