package prot;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

/** Explicit schema3 size transforms and feature bindings; never inferred from network dimensions.
 * @author Brian Bushnell, Nilou
 */
final class MagQCV14Contract {

	/** Parses the released metadata and retains the actual declared arithmetic constants. */
	MagQCV14Contract(final Map<String, String> metadata) throws IOException{
		assert(metadata!=null) : "Size arithmetic must be read from the authenticated bundle metadata";
		final String mode=required(metadata, "v14_fixture");
		if(!mode.equals("true") && !mode.equals("false")){throw new IOException("v14_fixture must be literal true or false");}
		fixture=mode.equals("true");
		ordinaryExtra=integer(metadata, "v14_ordinary_input_extra");
		if(ordinaryExtra!=10 && !(fixture && ordinaryExtra==0)){throw new IOException("Production ordinary inputs require ten additional features");}
		requirednessPin=DigestSuffix.requireSuffix(required(metadata, "v14_requiredness_sha80"), "v1.4 requiredness");
		markersPin=DigestSuffix.requireSuffix(required(metadata, "v14_size_markers_sha80"), "v1.4 size markers");
		bpOffset=number(metadata, "size_log10_bp_offset"); bpScale=number(metadata, "size_log10_bp_scale");
		cdsOffset=number(metadata, "size_log10_cds_offset"); cdsScale=number(metadata, "size_log10_cds_scale");
		linearScale=number(metadata, "size_linear_bp_scale");
		if(bpOffset!=5 || bpScale!=3 || cdsOffset!=1 || cdsScale!=4 || linearScale!=10000000){
			throw new IOException("Size transforms differ from the released offset5/3, offset1/4 and linear1e7 contract");
		}
		if(!required(metadata, "size_log_target_clamp").equals("0,1")){
			throw new IOException("Both logarithmic size targets require the declared [0,1] clamp");
		}
		activations=required(metadata, "size_output_activations").split(",", -1);
		if(!Arrays.equals(activations, ACTIVATIONS)){throw new IOException("Size output activation order differs");}
		if(!required(metadata, "size_output_names").equals(String.join(",", NAMES))){throw new IOException("Size output names differ");}
	}

	/** Uses the raw linear head and its recorded unit scale, before six-decimal transport rounding. */
	double logRatio(final long binBp, final float rawLinear){
		if(binBp<=0 || !Float.isFinite(rawLinear)){throw new IllegalArgumentException("Ratio needs positive observed bp and finite raw size output");}
		final double predictedLog=rawLinear<=1/linearScale ? 0 : Math.log(rawLinear)+Math.log(linearScale);
		return (Math.log(binBp)-predictedLog)/Math.log(2);
	}

	/** Metadata comes from one pinned TSV; no duplicate or unknown parameters silently survive. */
	static Map<String, String> parse(final byte[] bytes) throws IOException{
		assert(bytes!=null) : "The packer must supply the exact pinned metadata bytes before interpreting fields";
		final LinkedHashMap<String, String> values=new LinkedHashMap<String, String>();
		for(String line:new String(bytes, StandardCharsets.UTF_8).split("\\r?\\n")){
			if(line.isEmpty() || line.startsWith("#")){continue;}
			final String[] fields=line.split("\t", -1);
			if(fields.length!=2 || !Arrays.asList(KEYS).contains(fields[0]) || values.put(fields[0], fields[1])!=null){
				throw new IOException("Unknown, repeated or malformed v1.4 metadata field");
			}
		}
		new MagQCV14Contract(values);
		return values;
	}

	/** Checks the activation token on each actual serialized output cell, including sparse networks. */
	void requireActivations(final byte[] bytes, final int[] dims) throws IOException{
		assert(bytes!=null) : "Activation validation must inspect the serialized network actually placed in the bundle";
		if(dims==null || dims.length<2 || dims[dims.length-1]!=10){throw new IOException("Size activation validation requires ten output cells");}
		int first=1; for(int i=0; i<dims.length-1; i++){first=Math.addExact(first, dims[i]);}
		final boolean[] seen=new boolean[10];
		for(int start=0; start<bytes.length; ){
			int end=start;while(end<bytes.length && bytes[end]!='\n'){end++;}
			if(end>start && (bytes[start]=='C' || bytes[start]=='W')){
				int pos=start+1;while(pos<end && bytes[pos]>='0' && bytes[pos]<='9'){pos++;}
				final int id=parse.Parse.parseInt(bytes, start+1, pos);
				if(id>=first){
					final int index=id-first;
					if(index>=10 || seen[index]){throw new IOException("Repeated or extra size output cell");}
					while(pos<end && bytes[pos]==' '){pos++;} final int word=pos;
					while(pos<end && bytes[pos]!=' ' && bytes[pos]!='\r'){pos++;}
					if(!activations[index].equals(new String(bytes, word, pos-word, StandardCharsets.US_ASCII))){
						throw new IOException("Size output activation disagrees with metadata at index "+index);
					}
					seen[index]=true;
				}
			}
			start=end+1;
		}
		for(boolean value:seen){if(!value){throw new IOException("Missing size output cell");}}
	}

	private static String required(final Map<String, String> md, final String key) throws IOException{
		final String value=md.get(key);if(value==null || value.isEmpty()){throw new IOException("Missing v1.4 metadata: "+key);}return value;
	}
	private static double number(final Map<String, String> md, final String key) throws IOException{
		try{final double value=Double.parseDouble(required(md, key));if(Double.isFinite(value)){return value;}}
		catch(NumberFormatException e){throw new IOException("Malformed v1.4 constant: "+key, e);}
		throw new IOException("Nonfinite v1.4 constant: "+key);
	}
	private static int integer(final Map<String, String> md, final String key) throws IOException{
		try{return Integer.parseInt(required(md, key));}catch(NumberFormatException e){throw new IOException("Malformed v1.4 integer: "+key, e);}
	}

	final boolean fixture;
	final int ordinaryExtra;
	final String requirednessPin, markersPin;
	final double bpOffset, bpScale, cdsOffset, cdsScale, linearScale;
	private final String[] activations;
	static final String[] NAMES={"size01", "size_lin", "cds01", "gene_completeness", "gene_contamination",
		"size01_error", "size_lin_error", "cds01_error", "gene_completeness_error", "gene_contamination_error"};
	static final String[] UNITS={"clamped affine log10 bp", "bp / size_linear_bp_scale", "clamped affine log10 CDS", "gene fraction", "gene fraction",
		"absolute size01 error", "absolute linear-size error", "absolute cds01 error", "absolute gene-completeness error", "absolute gene-contamination error"};
	private static final String[] ACTIVATIONS={"SIG", "LINEAR", "SIG", "SIG", "SIG", "SIG", "LINEAR", "SIG", "SIG", "SIG"};
	static final String[] KEYS={"v14_fixture", "v14_ordinary_input_extra", "v14_requiredness_sha80", "v14_size_markers_sha80",
		"size_log10_bp_offset", "size_log10_bp_scale", "size_log10_cds_offset", "size_log10_cds_scale", "size_linear_bp_scale",
		"size_log_target_clamp", "size_output_names", "size_output_activations"};
}
