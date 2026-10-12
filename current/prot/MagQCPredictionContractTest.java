package prot;

import java.util.Arrays;

/** Regression checks for model-declared gene predictions and unchanged auxiliary heads.
 * Synthetic predictions test interpretation, not biological accuracy.
 * @author Brian Bushnell, Yoimiya */
public final class MagQCPredictionContractTest {

	/** Tests both release widths, permutation, inverse limits, blend boundaries and malformed metadata. */
	public static void main(String[] args){
		for(int count:new int[]{10, 12}){
			final String[] names=names(count);
			checkLayout(names);
			final String temp=names[0]; names[0]=names[count-1]; names[count-1]=temp;
			checkLayout(names);
		}
		final String[] names=names(12);
		final String metadata=metadata(names, 6, 10);
		for(String invalid:new String[]{metadata.replace("log10", "unsupported"),
			metadata.replace("own_magnified_deviation", "linear_deviation"),
			metadata.replace("\"lo_percent\":6.0", "\"lo_percent\":10.0"),
			metadata.replace("\"log_scale\":3", "\"log_scale\":0"),
			metadata.replace("\"deviation_min_percent\":0.1", "\"deviation_min_percent\":0"),
			metadata.replace("\"gene_completeness\",\"kind\":\"linear_value\"", "\"wrong_name\",\"kind\":\"linear_value\""),
			metadata.substring(0, metadata.length()-1)}){
			boolean failed=false;
			try{MagQCPredictionContract.parse(invalid, names);}catch(IllegalArgumentException expected){failed=true;}
			check(failed, "Malformed declaration must fail before producing predictions");
		}
		final MagQCPredictionContract contract=MagQCPredictionContract.parse(metadata, names);
		final float[] invalid=new float[16]; invalid[0]=Float.NaN;
		boolean failed=false;
		try{contract.apply(invalid);}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Nonfinite predictions must fail before blending");
		System.out.println("MagQCPredictionContractTest PASS: 10/12 heads, reordered semantics, independent quantities, metadata windows, preserved auxiliary bits, invalid metadata");
	}

	/** Expected values are hand-selected examples; completeness and contamination use different deviations. */
	private static void checkLayout(final String[] names){
		final MagQCPredictionContract contract=MagQCPredictionContract.parse(metadata(names, 6, 10), names);
		check(contract.rawOutputs()==names.length && contract.outputNames().length==names.length+4,
			"Raw and diagnostic widths must remain distinct");
		final String[] changed=contract.outputNames(); changed[0]="changed";
		check(!contract.outputNames()[0].equals("changed"), "Returned names must not mutate the loaded contract");
		// Below window: wholly magnified. Inside: equal blend at 8%. Above: wholly linear.
		checkPrediction(contract, names, 1, 20, .99, .20);
		checkPrediction(contract, names, 8, 8, .86, .14);
		checkPrediction(contract, names, 20, 1, .80, .01);
		checkPrediction(contract, names, 6, 10, .94, .20);
		checkPrediction(contract, names, .001, 1000, .999, .20);
		// A different declared window proves inference is metadata-driven, not hardcoded to 6/10.
		final MagQCPredictionContract alternative=MagQCPredictionContract.parse(metadata(names, 2, 4), names);
		checkPrediction(alternative, names, 3, 3, .885, .115);
	}

	/** Checks unchanged raw heads bitwise and the blended values against independent expected numbers. */
	private static void checkPrediction(final MagQCPredictionContract contract, final String[] names,
			final double compDeviation, final double contamDeviation, final double compExpected, final double contamExpected){
		final float[] values=new float[names.length+4];
		for(int i=0; i<names.length; i++){values[i]=(i+1)*.03125f;}
		final int comp=index(names, "gene_completeness"), contam=index(names, "gene_contamination");
		values[comp]=.8f; values[contam]=.2f;
		values[index(names, "magnified_completeness")]=(float)((Math.log10(compDeviation)+1)/3);
		values[index(names, "magnified_contamination")]=(float)((Math.log10(contamDeviation)+1)/3);
		final float[] before=values.clone();
		contract.apply(values);
		near(values[comp], compExpected); near(values[contam], contamExpected);
		for(int i=0; i<names.length; i++){
			if(i!=comp && i!=contam){check(Float.floatToRawIntBits(values[i])==Float.floatToRawIntBits(before[i]),
				"BP, error, binary and raw magnified heads must remain bit-identical");}
		}
		check(values[names.length]==before[comp] && values[names.length+1]==before[contam],
			"Diagnostic linear values must retain the pre-blend predictions");
		near(values[names.length+2], 1-Math.max(.1, Math.min(100, compDeviation))/100);
		near(values[names.length+3], Math.max(.1, Math.min(100, contamDeviation))/100);
	}

	/** Builds a complete declaration from names, independently of production output indices. */
	static String metadata(final String[] names, final double lo, final double hi){
		final StringBuilder text=new StringBuilder("{\"schema\":\"magqc_prediction_contract_v1\",\"heads\":[");
		for(int i=0; i<names.length; i++){
			if(i>0){text.append(',');}
			final String name=names[i];
			final boolean magnified=name.startsWith("magnified_");
			final String kind=magnified ? "magnified_value" : name.startsWith("error_") ? "error" :
				name.startsWith("binary_") ? "binary" : "linear_value";
			final String quantity=(name.startsWith("bp_") ? "bp_" : "gene_")+
				(name.endsWith("completeness") ? "completeness" : "contamination");
			text.append("{\"name\":\"").append(name).append("\",\"kind\":\"").append(kind)
				.append("\",\"quantity\":\"").append(quantity).append("\",\"transform\":\"")
				.append(magnified ? "log10" : "identity").append('"');
			if(magnified){text.append(",\"orientation\":\"direct\",\"deviation_min_percent\":0.1,\"deviation_max_percent\":100,\"log_scale\":3,\"log_offset\":1");}
			text.append('}');
		}
		return text.append("],\"ramp\":{\"key\":\"own_magnified_deviation\",\"lo_percent\":")
			.append(lo).append(",\"hi_percent\":").append(hi).append("}}").toString();
	}

	/** Uses distinct auxiliary names so accidental field overwrites are observable. */
	private static String[] names(final int count){
		final String[] all={"gene_completeness", "gene_contamination", "magnified_completeness", "magnified_contamination",
			"binary_gene_completeness", "binary_gene_contamination", "bp_completeness", "bp_contamination",
			"error_bp_completeness", "error_bp_contamination", "error_gene_completeness", "error_gene_contamination"};
		return Arrays.copyOf(all, count);
	}
	/** Missing fixture names indicate an invalid test, not an inference fallback. */
	private static int index(final String[] names, final String name){
		for(int i=0; i<names.length; i++){if(names[i].equals(name)){return i;}}
		throw new AssertionError("Missing fixture output: "+name);
	}
	/** Allows only rounding of float log predictions and the final float result. */
	private static void near(final double actual, final double expected){
		check(Math.abs(actual-expected)<=2e-7, "Expected "+expected+", observed "+actual);
	}
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
