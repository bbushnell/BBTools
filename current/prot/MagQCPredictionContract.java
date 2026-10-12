package prot;

import java.util.Arrays;
import java.util.HashSet;

import json.JsonObject;
import json.JsonParser;

/** Interprets model-declared value heads and blending without assuming their indices.
 * Loading allocates the contract; per-bin interpretation reuses the caller's output array.
 * @author Brian Bushnell, Yoimiya */
public final class MagQCPredictionContract {

	/** Validates model metadata against the independently parsed raw output names. */
	public static MagQCPredictionContract parse(final String text, final String[] rawNames){
		require(text!=null && rawNames!=null && rawNames.length>0, "Prediction metadata and raw names are required");
		final JsonParser validator=new JsonParser(text);
		require(validator.validate(), "Malformed prediction metadata JSON");
		final JsonObject root=new JsonParser(text).parseJsonObject();
		require("magqc_prediction_contract_v1".equals(string(root, "schema")), "Unsupported prediction metadata schema");
		final Object value=root.omap==null ? null : root.omap.get("heads");
		require(value instanceof Object[], "Prediction metadata requires a head array");
		final Object[] heads=(Object[])value;
		require(heads.length==rawNames.length, "Prediction metadata must describe every raw model output");
		final int[] linear={-1, -1}, magnified={-1, -1};
		final double[] minimum=new double[2], maximum=new double[2], scale=new double[2], offset=new double[2];
		final HashSet<String> seen=new HashSet<String>();
		for(int i=0; i<heads.length; i++){
			require(heads[i] instanceof JsonObject, "Every head declaration must be an object");
			final JsonObject head=(JsonObject)heads[i];
			final String name=string(head, "name"), kind=string(head, "kind"), quantity=string(head, "quantity");
			require(name.equals(rawNames[i]) && seen.add(name), "Head names/order differ from the raw model contract");
			require(KINDS.contains(kind) && QUANTITIES.contains(quantity), "Unknown head kind or quantity: "+name);
			if(!kind.equals("magnified_value")){
				require("identity".equals(string(head, "transform")), "Non-magnified heads require the identity reporting transform");
			}
			final int axis=quantity.equals("gene_completeness") ? 0 : quantity.equals("gene_contamination") ? 1 : -1;
			if(kind.equals("linear_value") && axis>=0){
				require(name.equals(quantity), "Reported gene value heads must use their quantity names");
				require(linear[axis]<0, "Duplicate linear value head for "+quantity);
				linear[axis]=i;
			}else if(kind.equals("magnified_value")){
				require(axis>=0 && magnified[axis]<0, "Duplicate or unsupported magnified quantity: "+quantity);
				require("log10".equals(string(head, "transform")) && "direct".equals(string(head, "orientation")),
					"Magnified heads require an explicitly declared direct log10 transform");
				magnified[axis]=i;
				minimum[axis]=number(head, "deviation_min_percent");
				maximum[axis]=number(head, "deviation_max_percent");
				scale[axis]=number(head, "log_scale"); offset[axis]=number(head, "log_offset");
				require(minimum[axis]>0 && maximum[axis]>minimum[axis] && maximum[axis]<=100 && scale[axis]>0,
					"Magnified transform has invalid percentage bounds or scale");
			}
		}
		final JsonObject ramp=root.getJson("ramp");
		require(ramp!=null && "own_magnified_deviation".equals(string(ramp, "key")), "Unsupported blend key");
		final double lo=number(ramp, "lo_percent"), hi=number(ramp, "hi_percent");
		require(lo>=0 && hi>lo && hi<=100, "Blend percentage bounds must be ordered and finite");
		for(int axis=0; axis<2; axis++){
			require(linear[axis]>=0 && magnified[axis]>=0, "Both gene quantities require linear and magnified value heads");
			require(lo>=minimum[axis] && hi<=maximum[axis], "Blend window lies outside the declared inverse range");
		}
		return new MagQCPredictionContract(rawNames, linear, magnified, minimum, maximum, scale, offset, lo, hi);
	}

	/** Retains only validated immutable interpretation data, never the mutable JSON tree. */
	private MagQCPredictionContract(final String[] rawNames, final int[] linear_, final int[] magnified_,
			final double[] minimum_, final double[] maximum_, final double[] scale_, final double[] offset_,
			final double lo_, final double hi_){
		rawCount=rawNames.length;
		names=Arrays.copyOf(rawNames, rawCount+DIAGNOSTICS.length);
		final HashSet<String> unique=new HashSet<String>(Arrays.asList(rawNames));
		for(int i=0; i<DIAGNOSTICS.length; i++){
			require(unique.add(DIAGNOSTICS[i]), "Diagnostic name collides with a model head");
			names[rawCount+i]=DIAGNOSTICS[i];
		}
		linear=linear_.clone(); magnified=magnified_.clone();
		minimum=minimum_.clone(); maximum=maximum_.clone(); scale=scale_.clone(); offset=offset_.clone();
		lo=lo_; hi=hi_;
	}

	/** Raw network width is distinct from the report width, which includes four diagnostics. */
	public int rawOutputs(){return rawCount;}
	/** Returns private names in the same order as the interpreted result array. */
	public String[] outputNames(){return names.clone();}

	/** Replaces only the two gene value slots; raw magnified, BP, error and binary slots remain unchanged. */
	public void apply(final float[] output){
		require(output!=null && output.length==names.length, "Prediction scratch must include diagnostic slots");
		for(int i=0; i<rawCount; i++){
			if(!Float.isFinite(output[i])){throw new IllegalArgumentException("Nonfinite raw prediction at head "+i);}
		}
		for(int axis=0; axis<2; axis++){
			final double original=output[linear[axis]], f=output[magnified[axis]];
			final double d=Math.max(minimum[axis], Math.min(maximum[axis], Math.pow(10, scale[axis]*f-offset[axis])));
			final double decoded=(axis==0 ? 100-d : d)/100;
			final double weight=Math.max(0, Math.min(1, (hi-d)/(hi-lo)));
			final double blended=weight==0 ? original : weight==1 ? decoded : weight*decoded+(1-weight)*original;
			output[rawCount+axis]=(float)original;
			output[rawCount+2+axis]=(float)decoded;
			output[linear[axis]]=(float)blended;
			assert(Float.isFinite(output[linear[axis]])) : "Finite source predictions and bounded convex weights must produce a finite value";
		}
	}

	/** Reads typed metadata without relying on assertion-only JsonObject casts. */
	private static String string(final JsonObject object, final String key){
		final Object value=object==null || object.omap==null ? null : object.omap.get(key);
		require(value instanceof String && !((String)value).isEmpty(), "Required string metadata: "+key);
		return (String)value;
	}
	/** Metadata numbers must be finite before they reach inverse or blend arithmetic. */
	private static double number(final JsonObject object, final String key){
		final Object value=object==null || object.omap==null ? null : object.omap.get(key);
		require(value instanceof Number, "Required numeric metadata: "+key);
		final double result=((Number)value).doubleValue();
		require(Double.isFinite(result), "Nonfinite metadata: "+key);
		return result;
	}
	/** Invalid model semantics fail before inference or output publication. */
	private static void require(final boolean condition, final String why){
		if(!condition){throw new IllegalArgumentException(why);}
	}

	private final int rawCount;
	private final String[] names;
	private final int[] linear, magnified;
	private final double[] minimum, maximum, scale, offset;
	private final double lo, hi;
	private static final String[] DIAGNOSTICS={"linear_gene_completeness", "linear_gene_contamination",
		"decoded_magnified_gene_completeness", "decoded_magnified_gene_contamination"};
	private static final HashSet<String> KINDS=new HashSet<String>(Arrays.asList("linear_value", "magnified_value", "error", "binary"));
	private static final HashSet<String> QUANTITIES=new HashSet<String>(Arrays.asList("gene_completeness", "gene_contamination",
		"bp_completeness", "bp_contamination"));
}
