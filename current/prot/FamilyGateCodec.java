package prot;

import parse.LineParser1;
import structures.ByteBuilder;

/** Directed float32 cutoff rounding and canonical A48 raw-bit encoding.
 * Cascade-v2 computed metrics and decoded gates are floats; raw calibration stays double.
 * @author Sayu
 */
public final class FamilyGateCodec{

	private FamilyGateCodec(){}

	/** Greatest finite float not greater than the original lower threshold. */
	public static float lower(final double value){
		representable(value);
		float result=(float)value;
		if((double)result>value){result=Math.nextDown(result);}
		assert(Float.isFinite(result) && (double)result<=value) : "Lower rounding must not lose previously retained scores";
		return result;
	}

	/** Identity acceptance uses original percent measurements and rounds the final bound upward once. */
	public static float upper(final double originalLower,final double originalUpper){
		return upper(originalLower,originalUpper,DEFAULT_IDENTITY_UPPER_MULTIPLIER,DEFAULT_IDENTITY_UPPER_ADDITIVE);
	}

	/** Identity acceptance with an explicit proportional and additive upper margin. */
	public static float upper(final double originalLower,final double originalUpper,
			final double multiplier,final double additive){
		representable(originalLower);
		representable(originalUpper);
		if(originalLower<0 || originalLower>100 || originalUpper<0 || originalUpper>100){
			throw new IllegalArgumentException("Raw identity gates must be percentages in [0,100]");
		}
		if(!Double.isFinite(multiplier) || !Double.isFinite(additive) || multiplier<1 || additive<0){
			throw new IllegalArgumentException("Identity margins require a finite multiplier >= 1 and a finite additive >= 0");
		}
		return roundUpper(Math.max(originalLower,Math.max(30,Math.max(multiplier*originalUpper,originalUpper+additive))));
	}

	/** Least finite float not below a computed bound; refuses overflow rather than emitting infinity. */
	public static float roundUpper(final double value){
		representable(value);
		float result=(float)value;
		if((double)result<value){result=Math.nextUp(result);}
		assert(Float.isFinite(result) && (double)result>=value) : "Strict acceptance must never round its threshold downward";
		return result;
	}

	/** Final cosine decision uses score>max(lower,0.72,1.1*upper,upper+0.05), with no outer cap or disabled sentinel. */
	public static float cosineAccept(final double originalLower,final double originalUpper){
		representable(originalLower);
		representable(originalUpper);
		if(originalLower<0 || originalLower>1 || originalUpper<0 || originalUpper>1){
			throw new IllegalArgumentException("Original cosine gates must be in [0,1]: "+originalLower+", "+originalUpper);
		}
		final double bound=Math.max(originalLower,Math.max(0.72,Math.max(originalUpper*1.1,originalUpper+0.05)));
		float result=(float)bound;
		if((double)result<bound){result=Math.nextUp(result);}
		assert(Float.isFinite(result) && (double)result>=bound && (double)result>=originalLower) :
			"Strict cosine acceptance must exceed the raw floor/additive/proportional formula without rounding it downward";
		return result;
	}

	/** Final BLOSUM decision uses score>this finite bound, never the disabled shortcut sentinel.
	 * Apply max(lower,0,upper*1.1,upper+0.005) to ORIGINAL doubles, then upward-round once.
	 */
	public static float finalBlosum(final double originalLower,final double originalUpper){
		representable(originalLower);
		representable(originalUpper);
		final double bound=Math.max(originalLower,Math.max(0,Math.max(originalUpper*1.1,originalUpper+0.005)));
		final float result=roundUpper(bound);
		assert(Float.isFinite(result) && (double)result>=bound && (double)result>=originalLower) :
			"Strict final acceptance must exceed the uncapped raw margin formula without rounding it downward";
		return result;
	}

	private static void representable(final double value){
		if(!Double.isFinite(value) || value>Float.MAX_VALUE || value<-Float.MAX_VALUE){
			throw new IllegalArgumentException("Original gate is nonfinite or outside finite float32 range: "+value);
		}
	}

	public static ByteBuilder appendLower(final ByteBuilder out,final float value){
		if(!Float.isFinite(value)){throw new IllegalArgumentException("Lower gate must be finite");}
		return out.appendFloatA48(value);
	}

	public static ByteBuilder appendUpper(final ByteBuilder out,final float value){
		if(!Float.isFinite(value)){throw new IllegalArgumentException("Upper gate must be finite");}
		return out.appendFloatA48(value);
	}

	public static ByteBuilder appendCosineAccept(final ByteBuilder out,final float value){
		if(!Float.isFinite(value)){throw new IllegalArgumentException("Final cosine gate must be finite; no disabled sentinel");}
		return out.appendFloatA48(value);
	}

	public static ByteBuilder appendFinalBlosum(final ByteBuilder out,final float value){
		if(!Float.isFinite(value)){throw new IllegalArgumentException("Final BLOSUM gate must be finite; no disabled sentinel");}
		return out.appendFloatA48(value);
	}

	public static float readLower(final LineParser1 line,final int term){
		final float value=decode(line,term);
		if(!Float.isFinite(value)){throw new IllegalArgumentException("Nonfinite encoded lower gate");}
		return value;
	}

	public static float readUpper(final LineParser1 line,final int term){
		final float value=decode(line,term);
		if(!Float.isFinite(value)){throw new IllegalArgumentException("Nonfinite encoded upper gate; re-export legacy gates");}
		return value;
	}

	public static float readCosineAccept(final LineParser1 line,final int term){
		final float value=decode(line,term);
		if(!Float.isFinite(value)){throw new IllegalArgumentException("Final cosine gate must be finite; re-export legacy sentinel from raw cutoffs");}
		return value;
	}

	public static float readFinalBlosum(final LineParser1 line,final int term){
		final float value=decode(line,term);
		if(!Float.isFinite(value)){throw new IllegalArgumentException("Final BLOSUM gate must be finite; re-export legacy sentinel from raw cutoffs");}
		return value;
	}

	/** Identity percentages cannot strictly exceed this finite endpoint. */
	public static boolean disabled(final float value){return Float.isFinite(value) && value>=DISABLED_UPPER;}

	private static float decode(final LineParser1 line,final int term){
		if(term<0 || term>=line.terms()){throw new IllegalArgumentException("Missing A48 gate field");}
		final int length=line.setBounds(term), start=line.a();
		final byte[] bytes=line.line();
		if(length<1 || length>6 || (length>1 && bytes[start]=='0')){
			throw new IllegalArgumentException("Noncanonical A48 gate length or leading zero");
		}
		for(int i=start; i<line.b(); i++){
			if(bytes[i]<48 || bytes[i]>111){throw new IllegalArgumentException("Invalid A48 gate alphabet");}
		}
		// ParseLongA48 itself does not validate alphabet or overflow; validate before narrowing to int.
		final long bits=line.parseLongA48(term);
		if(bits<0 || bits>0xffffffffL){throw new IllegalArgumentException("A48 gate exceeds unsigned32 bits");}
		return Float.intBitsToFloat((int)bits);
	}

	public static final float DISABLED_UPPER=100f;
	public static final int DISABLED_BITS=0x42c80000;
	public static final double DEFAULT_IDENTITY_UPPER_MULTIPLIER=1.1;
	public static final double DEFAULT_IDENTITY_UPPER_ADDITIVE=2;
	public static final String SCHEMA="family_gate_a48_f32_v6";
	public static final String LOWER_HEADER="roster_rank\tfamily_id\tmembers\tstatus\tcosine_cutoff_a48f32_v1\tcosine_retained_min\tidentity_pct_cutoff_a48f32_v1\tidentity_retained_min\tblosum_a_sq_cutoff_a48f32_v1\tblosum_retained_min";
	public static final String UPPER_HEADER="roster_rank\tfamily_id\tnegative_count\tcosine_accept_a48f32_v2\taligned_count\tidentity_pct_accept_a48f32_v3\tblosum_a_sq_final_a48f32_v4";
}
