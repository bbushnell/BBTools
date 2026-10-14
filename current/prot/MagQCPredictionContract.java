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
		final String schema=string(root, "schema");
		final boolean version2="magqc_prediction_contract_v2".equals(schema);
		require(version2 || "magqc_prediction_contract_v1".equals(schema), "Unsupported prediction metadata schema");
		final Object value=root.omap==null ? null : root.omap.get("heads");
		require(value instanceof Object[], "Prediction metadata requires a head array");
		final Object[] heads=(Object[])value;
		require(heads.length==rawNames.length, "Prediction metadata must describe every raw model output");
		final int[] linear={-1, -1}, magnified={-1, -1};
		final double[] minimum=new double[2], maximum=new double[2], scale=new double[2], offset=new double[2];
		final boolean[] logErrors=new boolean[rawNames.length];
		final String[] errorQuantities=new String[rawNames.length];
		final HashSet<String> servedErrors=new HashSet<String>();
		final HashSet<String> linearQuantities=new HashSet<String>();
		final HashSet<String> seen=new HashSet<String>();
		for(int i=0; i<heads.length; i++){
			require(heads[i] instanceof JsonObject, "Every head declaration must be an object");
			final JsonObject head=(JsonObject)heads[i];
			final String name=string(head, "name"), kind=string(head, "kind"), quantity=string(head, "quantity");
			require(name.equals(rawNames[i]) && seen.add(name), "Head names/order differ from the raw model contract");
			require(KINDS.contains(kind) && QUANTITIES.contains(quantity), "Unknown head kind or quantity: "+name);
			if(version2 && kind.equals("error") && "log10".equals(string(head, "transform"))){
				require("served".equals(string(head, "reference")) && quantity.equals(string(head, "reference_output")),
					"Log-error heads must describe the served prediction of their own quantity");
				require("log_target".equals(string(head, "unit")) &&
					"fraction_absolute_error".equals(string(head, "reported_unit")), "Log-error reporting units must be explicit");
				require(number(head, "deviation_min_percent")==.1 && number(head, "deviation_max_percent")==100 &&
					number(head, "log_scale")==3 && number(head, "log_offset")==1,
					"Log-error heads must use the declared v1.5 0.1–100 percentage-point transform");
				logErrors[i]=true;
				errorQuantities[i]=quantity;
				require(servedErrors.add(quantity), "Duplicate log-error head for "+quantity);
			}else if(!kind.equals("magnified_value")){
				require("identity".equals(string(head, "transform")), "Non-magnified heads require the identity reporting transform");
			}
			final int axis=quantity.equals("gene_completeness") ? 0 : quantity.equals("gene_contamination") ? 1 : -1;
			if(kind.equals("linear_value")){linearQuantities.add(quantity);}
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
		final boolean blend=magnified[0]>=0 || magnified[1]>=0;
		require(version2 || blend, "Legacy prediction metadata requires magnified values");
		require(!version2 || !servedErrors.isEmpty(), "Version2 metadata requires at least one log-error head");
		// JsonObject keeps nested objects in jmap, separately from scalar/array omap.
		require(!root.omap.containsKey("ramp"), "Blend ramp must be an object");
		final JsonObject ramp=root.getJson("ramp");
		double lo=0, hi=0;
		if(blend){
			require(ramp!=null && "own_magnified_deviation".equals(string(ramp, "key")), "Unsupported blend key");
			lo=number(ramp, "lo_percent"); hi=number(ramp, "hi_percent");
			require(lo>=0 && hi>lo && hi<=100, "Blend percentage bounds must be ordered and finite");
		}else{require(ramp==null, "Unblended predictions cannot declare a blend ramp");}
		for(int axis=0; axis<2; axis++){
			require(linear[axis]>=0 && (!blend || magnified[axis]>=0), "Both gene quantities require their declared value heads");
			if(blend){require(lo>=minimum[axis] && hi<=maximum[axis], "Blend window lies outside the declared inverse range");}
		}
		require(linearQuantities.containsAll(servedErrors), "Log-error reference lacks a linear value head");
		return new MagQCPredictionContract(rawNames, linear, magnified, minimum, maximum, scale, offset, lo, hi,
			blend, logErrors, errorQuantities, servedErrors);
	}

	/** Retains only validated immutable interpretation data, never the mutable JSON tree. */
	private MagQCPredictionContract(final String[] rawNames, final int[] linear_, final int[] magnified_,
			final double[] minimum_, final double[] maximum_, final double[] scale_, final double[] offset_,
			final double lo_, final double hi_, final boolean blend_, final boolean[] logErrors_,
			final String[] errorQuantities_, final HashSet<String> servedErrors_){
		rawCount=rawNames.length;
		blend=blend_; logErrors=logErrors_.clone(); servedErrors=new HashSet<String>(servedErrors_);
		errorQuantities=errorQuantities_.clone();
		names=Arrays.copyOf(rawNames, rawCount+(blend ? DIAGNOSTICS.length : 0));
		final HashSet<String> unique=new HashSet<String>(Arrays.asList(rawNames));
		for(int i=0; blend && i<DIAGNOSTICS.length; i++){
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
	/** Error calibration and report provenance use the quantity the head actually modeled. */
	public String errorReferenceName(final String quantity){
		return blend && !servedErrors.contains(quantity) && quantity.startsWith("gene_") ? "linear_"+quantity : quantity;
	}
	/** Indicates whether reports include blended gene values and linear diagnostics. */
	public boolean hasBlend(){return blend;}
	/** Indicates that declared log errors are decoded into the existing fractional API. */
	public boolean hasLogErrors(){return !servedErrors.isEmpty();}

	/** Raw target metadata makes inverse decoding mandatory, never an optional report decoration. */
	static boolean requiresLogDecoding(final String rawContract){
		return rawContract!=null && rawContract.contains("\"log_abs_error\"");
	}

	/** Cross-binds decoded slots to their training targets so missing or swapped metadata fails loudly. */
	void validateRawContract(final String text){
		if(!hasLogErrors() && !requiresLogDecoding(text)){return;}
		require(text!=null && new JsonParser(text).validate(), "Log errors require valid raw output metadata");
		final JsonObject root=new JsonParser(text).parseJsonObject();
		require("explicit".equals(string(root, "mode")) && "explicit_targets_v1".equals(string(root, "schema")),
			"Log errors require explicit raw targets");
		final Object rows=root.omap.get("heads");
		require(rows instanceof Object[] && ((Object[])rows).length==rawCount, "Raw error head count differs from reporting metadata");
		final Object[] heads=(Object[])rows;
		for(int i=0; i<rawCount; i++){
			require(heads[i] instanceof JsonObject, "Raw output head must be an object");
			final JsonObject head=(JsonObject)heads[i];
			require(names[i].equals(string(head, "name")), "Raw and reported head order differs");
			require(logErrors[i]=="log_abs_error".equals(string(head, "kind")), "Log-error decoding slots differ from training targets");
			if(logErrors[i]){
				final int column=index(head, "column", QUANTITY_ORDER.length);
				require(errorQuantities[i].equals(QUANTITY_ORDER[column]), "Raw and reported error quantities differ");
				final int reference=index(head, "prediction", rawCount);
				require(names[reference].equals(errorQuantities[i]) && heads[reference] instanceof JsonObject,
					"Raw error target references another quantity");
				final JsonObject valueHead=(JsonObject)heads[reference];
				require("column".equals(string(valueHead, "kind")) && index(valueHead, "column", QUANTITY_ORDER.length)==column,
					"Raw error reference must be the linear value of its declared quantity");
				if(blend && column<2){
					require(index(head, "magnified_prediction", rawCount)==magnified[column], "Raw error target uses a different magnified head");
					final Object declared=head.omap.get("ramp");
					require(declared instanceof Object[] && ((Object[])declared).length==2, "Raw blend ramp requires two bounds");
					final Object[] bounds=(Object[])declared;
					require(bounds[0] instanceof Number && bounds[1] instanceof Number &&
						((Number)bounds[0]).doubleValue()==lo && ((Number)bounds[1]).doubleValue()==hi,
						"Error target and serving blend windows differ");
					require(heads[magnified[column]] instanceof JsonObject, "Raw magnified head must be an object");
					final JsonObject mag=(JsonObject)heads[magnified[column]];
					require("endpoint".equals(string(mag, "kind")) && "log10".equals(string(mag, "transform")) &&
						Boolean.FALSE.equals(mag.omap.get("complement")) &&
						(column==0 ? "one_minus" : "identity").equals(string(mag, "deviation")) &&
						index(mag, "column", QUANTITY_ORDER.length)==column, "Raw magnified target differs from serving inverse");
				}else{require(!head.omap.containsKey("magnified_prediction") && !head.omap.containsKey("ramp"),
					"Unblended error target cannot declare a blend");}
			}
		}
	}

	/** Typed integer indices prevent rounded or out-of-range metadata from selecting another head. */
	private static int index(final JsonObject head, final String key, final int limit){
		final Object value=head.omap.get(key);
		require(value instanceof Long, "Raw target index must be an integer: "+key);
		final long n=((Long)value).longValue();
		require(n>=0 && n<limit, "Raw target index outside its declared range: "+key);
		return (int)n;
	}

	/** Applies declared blends and log-error inverses; legacy auxiliary slots remain unchanged. */
	public void apply(final float[] output){
		require(output!=null && output.length==names.length, "Prediction scratch must include diagnostic slots");
		for(int i=0; i<rawCount; i++){
			if(!Float.isFinite(output[i])){throw new IllegalArgumentException("Nonfinite raw prediction at head "+i);}
		}
		for(int axis=0; blend && axis<2; axis++){
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
		for(int i=0; i<rawCount; i++){
			if(logErrors[i]){output[i]=(float)(Math.max(.1, Math.min(100, Math.pow(10, 3.0*output[i]-1)))/100);}
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
	private final boolean blend;
	private final boolean[] logErrors;
	private final String[] errorQuantities;
	private final HashSet<String> servedErrors;
	private static final String[] DIAGNOSTICS={"linear_gene_completeness", "linear_gene_contamination",
		"decoded_magnified_gene_completeness", "decoded_magnified_gene_contamination"};
	private static final HashSet<String> KINDS=new HashSet<String>(Arrays.asList("linear_value", "magnified_value", "error", "binary"));
	private static final String[] QUANTITY_ORDER={"gene_completeness", "gene_contamination", "bp_completeness", "bp_contamination"};
	private static final HashSet<String> QUANTITIES=new HashSet<String>(Arrays.asList("gene_completeness", "gene_contamination",
		"bp_completeness", "bp_contamination"));
}
