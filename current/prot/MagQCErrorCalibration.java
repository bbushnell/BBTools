package prot;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import fileIO.ByteFile;
import parse.LineParser1;
import structures.ByteBuilder;

/**
 * Immutable post-composite error scaling. Optional tables contain explicit
 * rectangular knots and two nonnegative multiplier surfaces, never fitted values
 * inferred by the reader. Lookup clamps coordinates, not reported predictions.
 * All units are fractions. This does not turn error estimates into confidence
 * intervals or alter the raw six-head composite API.
 * @author Collei
 */
final class MagQCErrorCalibration {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads and validates an optional table before workers or outputs are opened. */
	static MagQCErrorCalibration fromOptions(final HashMap<String,String> options){
		final double comp=positive(required(options, "comperrormultiplier"));
		final double contam=positive(required(options, "contamerrormultiplier"));
		final boolean path=options.containsKey("errorcalibration"), pin=options.containsKey("errorcalibrationsha80");
		if(path!=pin){throw new IllegalArgumentException("errorcalibration and errorcalibrationsha80 are required together");}
		if(!path){return new MagQCErrorCalibration(comp, contam);}
		if(comp!=1 || contam!=1){throw new IllegalArgumentException("A calibration table cannot also use nonidentity scalar error multipliers");}
		final MagQCErrorCalibration table=load(required(options, "errorcalibration"),
			required(options, "errorcalibrationsha80"), options);
		// Config defaults are placeholders, not permission to mislabel a loaded table.
		for(int i=0; i<FIT_KEYS.length; i++){
			final String key=FIT_KEYS[i], prior=options.get(key), actual=table.metadata.get(key);
			if(prior!=null && !prior.equals(FIT_DEFAULTS[i]) && !prior.equals(actual)){
				throw new IllegalArgumentException("Calibration metadata conflicts with configured "+key);
			}
		}
		for(String key:FIT_KEYS){options.put(key, table.metadata.get(key));}
		return table;
	}

	/** The legacy case retains precisely one multiplication per error head. */
	private MagQCErrorCalibration(final double comp, final double contam){
		compKnots=contamKnots=null; compFactors=new double[]{comp}; contamFactors=new double[]{contam};
		metadata=null; pin=null;
	}

	/** Takes ownership only of freshly parsed arrays and metadata; none are exposed. */
	private MagQCErrorCalibration(final double[] c, final double[] t, final double[] cm, final double[] tm,
			final HashMap<String,String> headers, final String pin_){
		assert(c.length>0 && t.length>0 && cm.length==c.length*t.length && tm.length==cm.length) :
			"A complete rectangular grid is required by the flattened interpolation indexes";
		compKnots=c; contamKnots=t; compFactors=cm; contamFactors=tm; metadata=headers; pin=pin_;
	}

	/** Strict metadata followed by one row per node; scientific pins cover decoded text bytes. */
	private static MagQCErrorCalibration load(final String path, final String pin, final Map<String,String> options){
		checkPin(pin, "errorcalibrationsha80");
		if(!pin.equals(MagQCTextResource.sha80(path))){throw new IllegalArgumentException("Calibration table sha80 mismatch: "+path);}
		final HashMap<String,String> headers=new HashMap<String,String>();
		final LineParser1 parser=new LineParser1('\t');
		final ByteFile input=ByteFile.makeByteFile(path, true);
		double[] c=null, t=null, cm=null, tm=null;
		boolean[] seen=null;
		int rows=0, lineNumber=0;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lineNumber++;
				if(line.length==0){continue;}
				parser.set(line);
				if(line[0]=='#'){
					if(cm!=null){throw new IllegalArgumentException("Calibration headers must precede node rows");}
					if(parser.termEquals("#columns", 0)){
						if(parser.terms()!=5 || !parser.termEquals("comp_index", 1) || !parser.termEquals("contam_index", 2) ||
							!parser.termEquals("comp_multiplier", 3) || !parser.termEquals("contam_multiplier", 4)){
							throw new IllegalArgumentException("Unexpected calibration node columns");
						}
						validateHeaders(headers, options);
						c=knots(required(headers, "comp_knots")); t=knots(required(headers, "contam_knots"));
						final long size=(long)c.length*t.length;
						if(size>Integer.MAX_VALUE-8){throw new IllegalArgumentException("Calibration grid exceeds array capacity");}
						cm=new double[(int)size]; tm=new double[(int)size]; seen=new boolean[(int)size];
					}else{
						if(parser.terms()!=2){throw new IllegalArgumentException("Calibration metadata requires key and value");}
						final String key=parser.parseString(0).substring(1), value=parser.parseString(1);
						if(!HEADER_KEYS.contains(key) || value.isEmpty() || headers.put(key, value)!=null){
							throw new IllegalArgumentException("Unknown, empty or repeated calibration metadata: "+key);
						}
						for(int i=0; i<value.length(); i++){
							if(Character.isISOControl(value.charAt(i))){throw new IllegalArgumentException("Control character in calibration metadata: "+key);}
						}
					}
				}else{
					if(cm==null || parser.terms()!=4){throw new IllegalArgumentException("Calibration node requires declared columns and four fields");}
					final int ci=Integer.parseInt(parser.parseString(0)), ti=Integer.parseInt(parser.parseString(1));
					if(ci<0 || ci>=c.length || ti<0 || ti>=t.length){throw new IllegalArgumentException("Calibration node index outside knot grid");}
					final int index=ci*t.length+ti;
					if(seen[index]){throw new IllegalArgumentException("Duplicate calibration node "+ci+","+ti);}
					cm[index]=nodeFactor(parser.parseString(2)); tm[index]=nodeFactor(parser.parseString(3));
					seen[index]=true; rows++;
				}
			}
		}catch(RuntimeException e){throw new IllegalArgumentException("Invalid calibration table "+path+" at line "+lineNumber, e);}
		finally{if(input.close()){throw new IllegalArgumentException("Error reading calibration table: "+path);}}
		if(cm==null || rows!=cm.length){throw new IllegalArgumentException("Calibration table must supply every rectangular node: "+path);}
		return new MagQCErrorCalibration(c, t, cm, tm, headers, pin);
	}

	/** Prevents attaching a surface fitted for other weights, features or units. */
	private static void validateHeaders(final HashMap<String,String> headers, final Map<String,String> options){
		for(String key:HEADER_KEYS){if(!key.equals("fitobjective")){required(headers, key);}}
		match(headers, "schema", SCHEMA); match(headers, "units", "fraction");
		final String truth=required(headers,"truth");
		if(!truth.equals("gene") && !truth.equals("bp")){throw new IllegalArgumentException("Calibration truth must explicitly be gene or bp");}
		match(headers, "interpolation", "bilinear_multiplier");
		match(headers, "layout", MagQCVectorMaker.JOINT_COMPOSITE_LAYOUT);
		match(headers, "output_names", OUTPUT_NAMES);
		for(String key:BINDING_KEYS){
			final String value=required(headers, key);
			checkPin(value, key); match(headers, key, required(options, key));
		}
		checkPin(required(headers, "truthmanifestsha80"), "truthmanifestsha80");
		checkPin(required(headers, "splitmanifestsha80"), "splitmanifestsha80");
		if(required(headers, "errorfitset").equals("UNCALIBRATED") || required(headers, "errorfitset").equals("NA")){
			throw new IllegalArgumentException("A table must identify its fitted set; identity defaults are not a fitted calibration");
		}
		final String date=required(headers, "errorfitdate");
		if(!date.matches("[0-9]{4}-[0-9]{2}-[0-9]{2}")){throw new IllegalArgumentException("Calibration fit date must be YYYY-MM-DD");}
		java.time.LocalDate.parse(date);
		final String objective=headers.get("fitobjective");
		if(objective!=null && !Arrays.asList("ratio_means","mean_ratio","quantile95","least_squares","mae_min").contains(objective)){
			throw new IllegalArgumentException("Unknown calibration fit objective");
		}
		final String coverageText=required(headers,"errorcoverage");
		if(coverageText.equals("NA")){
			if(!Arrays.asList("ratio_means","mean_ratio","least_squares","mae_min").contains(objective)){
				throw new IllegalArgumentException("NA nominal coverage requires an explicit point-error fit objective");
			}
		}else{
			final double coverage=positive(coverageText);
			if(coverage>1){throw new IllegalArgumentException("Per-head nominal calibration coverage must be in (0,1]");}
		}
	}

	/** Explicit knots may be nonuniform; one knot is a constant along that axis. */
	private static double[] knots(final String value){
		final String[] terms=value.split(",", -1);
		final double[] result=new double[terms.length];
		for(int i=0; i<result.length; i++){
			result[i]=Double.parseDouble(terms[i]);
			if(!Double.isFinite(result[i]) || result[i]<0 || result[i]>1 || (i>0 && result[i]<=result[i-1])){
				throw new IllegalArgumentException("Calibration knots must be finite, strictly increasing and in [0,1]");
			}
		}
		return result;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Lookup             ----------------*/
	/*--------------------------------------------------------------*/

	/** Scales only completeness error; predictions remain owned by the raw scorer. */
	double scaleComp(final double comp, final double contam, final double error){return scale(comp, contam, error, compFactors);}
	/** Scales only contamination error using the same immutable coordinate grid. */
	double scaleContam(final double comp, final double contam, final double error){return scale(comp, contam, error, contamFactors);}

	/** Rejects malformed outputs even with assertions disabled, before a report can publish. */
	private double scale(final double comp, final double contam, final double error, final double[] factors){
		if(!Double.isFinite(comp) || !Double.isFinite(contam) || !Double.isFinite(error) || error<0){
			throw new IllegalArgumentException("Error calibration requires finite predictions and a finite nonnegative raw error");
		}
		final double factor=metadata==null ? factors[0] : interpolate(comp, contam, factors);
		final double result=error*factor;
		if(!Double.isFinite(result)){throw new IllegalArgumentException("Calibrated error overflow");}
		return result;
	}

	/** Two coordinate-weighted linear interpolations per axis; no per-bin allocations. */
	private double interpolate(final double comp, final double contam, final double[] factors){
		return interpolateSurface(comp,contam,factors,compKnots,contamKnots);
	}

	/** Shared consumer arithmetic for evaluating a fitted surface before its final metadata is emitted. */
	static double interpolateSurface(final double comp,final double contam,final double[] factors,
		final double[] compKnots,final double[] contamKnots){
		assert(compKnots.length>0 && contamKnots.length>0 && factors.length==compKnots.length*contamKnots.length) :
			"Validated rectangular node arrays are required by the shared consumer interpolation";
		final double c=Math.max(compKnots[0], Math.min(compKnots[compKnots.length-1], comp));
		final double t=Math.max(contamKnots[0], Math.min(contamKnots[contamKnots.length-1], contam));
		final int ci=lower(compKnots, c), ti=lower(contamKnots, t);
		final int cj=Math.min(ci+1, compKnots.length-1), tj=Math.min(ti+1, contamKnots.length-1);
		final double cw=ci==cj ? 0 : (c-compKnots[ci])/(compKnots[cj]-compKnots[ci]);
		final double tw=ti==tj ? 0 : (t-contamKnots[ti])/(contamKnots[tj]-contamKnots[ti]);
		final int width=contamKnots.length;
		final double a=lerp(factors[ci*width+ti], factors[ci*width+tj], tw);
		final double b=lerp(factors[cj*width+ti], factors[cj*width+tj], tw);
		return lerp(a, b, cw);
	}

	/** Brackets an already clamped coordinate, retaining the final knot exactly. */
	private static int lower(final double[] knots, final double value){
		assert(value>=knots[0] && value<=knots[knots.length-1]) : "Interpolation clamps the lookup coordinate before bracketing";
		final int found=Arrays.binarySearch(knots, value);
		return found>=0 ? found : -found-2;
	}

	/** Endpoint cases avoid cancellation; positive finite endpoints keep subtraction bounded. */
	private static double lerp(final double a, final double b, final double weight){
		assert(weight>=0 && weight<=1) : "A bracketed coordinate gives a convex interpolation weight";
		return weight==0 ? a : weight==1 ? b : a+(b-a)*weight;
	}

	/** Scalar mode emits nothing, preserving legacy TSV metadata byte for byte. */
	void appendProvenance(final ByteBuilder out){
		if(metadata==null){return;}
		out.append("#errorcalibration_schema\t").append(SCHEMA).nl();
		out.append("#errorcalibrationsha80\t").append(pin).nl();
		for(String key:PROVENANCE_KEYS){out.append("#errorcalibration_").append(key).tab().append(metadata.get(key)).nl();}
		if(metadata.containsKey("fitobjective")){out.append("#errorcalibration_fitobjective\t").append(metadata.get("fitobjective")).nl();}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Validation           ----------------*/
	/*--------------------------------------------------------------*/

	private static String required(final Map<String,String> map, final String key){
		final String value=map.get(key);
		if(value==null || value.isEmpty()){throw new IllegalArgumentException("Missing calibration setting: "+key);}
		return value;
	}
	private static double positive(final String value){
		final double number=Double.parseDouble(value);
		if(!Double.isFinite(number) || number<=0){throw new IllegalArgumentException("Calibration multipliers/coverage must be positive and finite: "+value);}
		return number;
	}
	/** A valid zero-residual cell may have factor zero; scalar/coverage contracts stay positive. */
	private static double nodeFactor(final String value){
		final double number=Double.parseDouble(value);
		if(!Double.isFinite(number) || number<0){throw new IllegalArgumentException("Calibration node multiplier must be finite and nonnegative: "+value);}
		return number;
	}
	private static void match(final Map<String,String> map, final String key, final String expected){
		if(!expected.equals(required(map, key))){throw new IllegalArgumentException("Calibration binding mismatch: "+key);}
	}
	private static void checkPin(final String value, final String key){
		if(!value.matches("[0-9a-f]{20}")){throw new IllegalArgumentException("Calibration "+key+" requires a lowercase 20-hex sha80");}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	private final double[] compKnots, contamKnots, compFactors, contamFactors;
	private final HashMap<String,String> metadata;
	private final String pin;
	static final String SCHEMA="prokcc_error_calibration_v1";
	static final String OUTPUT_NAMES="gene_completeness,gene_contamination,abs_residual_gene_completeness,"+
		"abs_residual_gene_contamination,contamination_ge_0.01,completeness_gt_0.99";
	static final String[] BINDING_KEYS={"netsha80", "bundlesha80", "familylistsha80", "subnetmanifestsha80",
		"expectedcopytablesha80", "subnetpopulationssha80", "profilesha80"};
	private static final String[] FIT_KEYS={"errorfitset", "errorfitdate", "errorcoverage"};
	private static final String[] FIT_DEFAULTS={"UNCALIBRATED", "NA", "NA"};
	private static final String[] PROVENANCE_KEYS={"units", "truth", "truthpolicy", "truthmanifestsha80",
		"splitmanifestsha80", "fitrecipe", "smoothing", "interpolation"};
	private static final HashSet<String> HEADER_KEYS=new HashSet<String>();
	static{
		HEADER_KEYS.addAll(Arrays.asList(BINDING_KEYS)); HEADER_KEYS.addAll(Arrays.asList(FIT_KEYS));
		HEADER_KEYS.addAll(Arrays.asList(PROVENANCE_KEYS));
		HEADER_KEYS.addAll(Arrays.asList("schema", "layout", "output_names", "comp_knots", "contam_knots"));
		HEADER_KEYS.add("fitobjective");// Optional for historical tables with a numeric coverage target.
	}
}
