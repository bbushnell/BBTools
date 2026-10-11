package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Locale;

import fileIO.ByteFile;
import parse.LineParser1;
import structures.ByteBuilder;

/** Independent numerical and malformed-input fixtures for the optional error grid.
 * @author Collei */
public final class MagQCErrorCalibrationTest {

	/** Tests the reader, writes a fixture surface, or grades already-produced integration reports. */
	public static void main(String[] args) throws Exception{
		if((args.length==3 || args.length==4) && args[0].startsWith("scalar=") && args[1].startsWith("table=") && args[2].startsWith("screen=")){
			if(args.length==4 && !args[3].startsWith("calibrationconfig=")){throw new IllegalArgumentException("Expected calibrationconfig=<fitted release.config>");}
			compareReports(args[0].substring(7), args[1].substring(6), args[2].substring(7), args.length==4 ? args[3].substring(18) : null);
			return;
		}
		if(args.length!=0){
			if(args.length!=2 || !args[0].startsWith("fixtureconfig=") || !args[1].startsWith("fixtureout=")){
				throw new IllegalArgumentException("Use no arguments, or fixtureconfig=<release.config> fixtureout=<new.tsv>");
			}
			final HashMap<String,String> options=MagQCAssemblyBatch.parseOptions(new String[]{"config="+args[0].substring(14)});
			final Path out=Paths.get(args[1].substring(11));
			Files.write(out, surface(options, "0", "0", "0\t0\t2\t3\n").getBytes(StandardCharsets.US_ASCII),
				java.nio.file.StandardOpenOption.CREATE_NEW);
			System.out.println("FIXTURE_ONLY_TABLE "+MagQCTextResource.sha80(out.toString()));
			return;
		}
		final HashMap<String,String> options=defaults();
		final MagQCErrorCalibration scalar=MagQCErrorCalibration.fromOptions(options);
		check(scalar.scaleComp(.5, .2, .03)==.03 && scalar.scaleContam(.5, .2, .04)==.04, "Identity scalar changed");
		options.put("comperrormultiplier", "2"); options.put("contamerrormultiplier", "3");
		final MagQCErrorCalibration scaled=MagQCErrorCalibration.fromOptions(options);
		check(scaled.scaleComp(.5, .2, .03)==.03*2 && scaled.scaleContam(.5, .2, .04)==.04*3, "Legacy scalar arithmetic changed");
		final ByteBuilder metadata=new ByteBuilder(); scaled.appendProvenance(metadata);
		check(metadata.length()==0, "Scalar mode must not add table metadata");
		for(double value:new double[]{-1, Double.NaN, Double.POSITIVE_INFINITY}){
			rejectScale(scaled, .5, .2, value);
		}
		rejectScale(scaled, Double.NaN, .2, .1); rejectScale(scaled, .5, Double.NEGATIVE_INFINITY, .1);
		rejectScale(scaled, .5, .2, Double.MAX_VALUE);
		final Path root=Files.createTempDirectory("prokcc-error-calibration-");
		final Path table=root.resolve("table.tsv"), config=root.resolve("release.config");
		try{
			final String text=surface(defaults(), "0,1", "0,0.25,1", NODES);
			final HashMap<String,String> active=write(table, text);
			final MagQCErrorCalibration grid=MagQCErrorCalibration.fromOptions(active);
			// Nonuniform asymmetric values are an independent oracle for axis order and weights.
			close(grid.scaleComp(.25, .625, 1), 11.5); close(grid.scaleContam(.25, .625, 1), 8.5);
			close(grid.scaleComp(.25, .625, .02), .23); close(grid.scaleContam(.25, .625, .03), .255);
			close(grid.scaleComp(0, .125, 1), 3); close(grid.scaleContam(0, .125, 1), 3);
			close(grid.scaleComp(1, 1, 1), 20); close(grid.scaleContam(1, 1, 1), 19);
			close(grid.scaleComp(-.2, 1.4, 1), 16); close(grid.scaleContam(-.2, 1.4, 1), 9);
			close(grid.scaleComp(1, 0, 1), 6); close(grid.scaleContam(1, .25, 1), 7);
			check(grid.scaleComp(.5, .5, 0)==0, "Multiplicative calibration must preserve zero raw error without an invented floor");
			check(active.get("errorfitset").equals("FIXTURE_ONLY") && active.get("errorcoverage").equals("0.95"), "Table metadata did not replace placeholders");
			grid.appendProvenance(metadata);
			check(metadata.toString().contains("#errorcalibrationsha80\t"+active.get("errorcalibrationsha80")), "Missing table provenance");
			check(metadata.toString().contains("#errorcalibration_truth\tgene"), "Missing truth semantics");
			final String pin=active.get("errorcalibrationsha80");
			Files.write(config, ("errorcalibration=table.tsv\nerrorcalibrationsha80="+pin+"\n").getBytes(StandardCharsets.US_ASCII));
			final HashMap<String,String> relocated=MagQCAssemblyBatch.parseOptions(new String[]{"config="+config});
			check(relocated.get("errorcalibration").equals(table.toString()), "Table path lost config-relative resolution");
			for(String key:MagQCErrorCalibration.BINDING_KEYS){relocated.put(key, PIN);}
			close(MagQCErrorCalibration.fromOptions(relocated).scaleComp(.25, .625, 1), 11.5);
			final HashMap<String,String> missing=defaults(); missing.put("errorcalibration", table.toString()); reject(missing);
			missing.remove("errorcalibration"); missing.put("errorcalibrationsha80", pin); reject(missing);
			final HashMap<String,String> conflict=new HashMap<String,String>(active);
			conflict.put("comperrormultiplier", "2"); reject(conflict);
			conflict.put("comperrormultiplier", "1"); conflict.put("errorfitset", "wrong-fit"); reject(conflict);
			conflict.put("errorfitset", "FIXTURE_ONLY"); conflict.put("errorcoverage", ".5"); reject(conflict);
			conflict.put("errorcoverage", "0.95"); conflict.put("errorcalibrationsha80", "00000000000000000000"); reject(conflict);
			for(String key:MagQCErrorCalibration.BINDING_KEYS){
				final HashMap<String,String> wrong=new HashMap<String,String>(active);
				wrong.put(key, "00000000000000000000"); reject(wrong);
			}
			final String[][] substitutions={
				{"#schema\t"+MagQCErrorCalibration.SCHEMA, "#schema\tunknown"},
				{"#units\tfraction", "#units\tpercent"}, {"#truth\tgene", "#truth\tunknown"},
				{"#comp_knots\t0,1", "#comp_knots\t0,0"}, {"#comp_knots\t0,1", "#comp_knots\t1,0"},
				{"#comp_knots\t0,1", "#comp_knots\t-1,1"}, {"#comp_knots\t0,1", "#comp_knots\t0,NaN"},
				{"#contam_knots\t0,0.25,1", "#contam_knots\t0,0.25,2"},
				{"#errorfitdate\t2026-10-03", "#errorfitdate\t2026-02-30"},
				{"#errorcoverage\t0.95", "#errorcoverage\t1.1"}, {"#errorcoverage\t0.95", "#errorcoverage\t0"},
				{"#errorfitset\tFIXTURE_ONLY", "#errorfitset\tUNCALIBRATED"},
				{"#layout\t"+MagQCVectorMaker.JOINT_COMPOSITE_LAYOUT, "#layout\twrong"},
				{"#output_names\t"+MagQCErrorCalibration.OUTPUT_NAMES, "#output_names\tgene_contamination,gene_completeness"},
				{"#smoothing\tnone", "#smoothing\tnone\n#smoothing\tnone"},
				{"#smoothing\tnone\n", ""}, {"#smoothing\tnone", "#unexpected\tnone"},
				{"1\t2\t20\t19\n", ""}, {"1\t2\t20\t19", "1\t1\t20\t19"},
				{"1\t2\t20\t19", "2\t2\t20\t19"}, {"1\t2\t20\t19", "-1\t2\t20\t19"},
				{"1\t2\t20\t19", "1x\t2\t20\t19"}, {"1\t2\t20\t19", "1\t2\t20"},
				{"1\t2\t20\t19", "1\t2\tNaN\t19"}, {"1\t2\t20\t19", "1\t2\t20\tInfinity"},
				{"1\t2\t20\t19", "1\t2\t-1\t19"}, {"1\t2\t20\t19", "1\t2\t20\t-1"}
			};
			for(String[] replacement:substitutions){
				check(text.contains(replacement[0]), "Malformed fixture must actually change its intended field");
				reject(write(table, text.replace(replacement[0], replacement[1])));
			}
			reject(write(table, text+"#smoothing\tlate\n"));
			reject(write(table, ""));
			final MagQCErrorCalibration bp=MagQCErrorCalibration.fromOptions(write(table,text.replace("#truth\tgene","#truth\tbp")));
			close(bp.scaleComp(.25,.625,1),11.5);
			final ByteBuilder bpMetadata=new ByteBuilder();bp.appendProvenance(bpMetadata);
			check(bpMetadata.toString().contains("#errorcalibration_truth\tbp"),"BP truth must remain explicit in reports");
			final String meanText=text.replace("#errorcoverage\t0.95","#fitobjective\tratio_means\n#errorcoverage\tNA");
			final HashMap<String,String> meanOptions=write(table,meanText);
			final MagQCErrorCalibration mean=MagQCErrorCalibration.fromOptions(meanOptions);close(mean.scaleComp(.25,.625,1),11.5);
			check(meanOptions.get("errorcoverage").equals("NA"),"Mean calibration must not invent a nominal coverage target");
			final ByteBuilder meanMetadata=new ByteBuilder();mean.appendProvenance(meanMetadata);
			check(meanMetadata.toString().contains("#errorcalibration_fitobjective\tratio_means"),"Mean fit objective provenance lost");
			reject(write(table,meanText.replace("ratio_means","quantile95")));
			reject(write(table,meanText.replace("ratio_means","unknown")));
			for(String objective:new String[]{"least_squares","mae_min"}){
				close(MagQCErrorCalibration.fromOptions(write(table,meanText.replace("ratio_means",objective))).scaleComp(.25,.625,1),11.5);
			}
			reject(write(table,text.replace("#errorcoverage\t0.95","#errorcoverage\tNA")));
			checkDegenerate(table, "0", "0", "0\t0\t2\t3\n", 2, 3);
			checkDegenerate(table, "0", "0", "0\t0\t0\t0\n", 0, 0);
			checkDegenerate(table, "0,1", "0", "0\t0\t0\t4\n1\t0\t2\t0\n", 1, 2);
			final HashMap<String,String> zeroScalar=defaults();zeroScalar.put("comperrormultiplier","0");reject(zeroScalar);
			checkDegenerate(table, "0", "0,1", "0\t0\t2\t3\n0\t1\t6\t11\n", 4, 7);
			checkDegenerate(table, "0,1", "0", "0\t0\t2\t3\n1\t0\t6\t11\n", 4, 7);
			checkDegenerate(table, "0.3,0.7", "0.2", "0\t0\t2\t3\n1\t0\t6\t11\n", 4, 7);
			System.out.println("MagQCErrorCalibrationTest PASS: scalar parity, bilinear/degenerate grids, boundaries, metadata and malformed inputs");
		}finally{Files.deleteIfExists(table); Files.deleteIfExists(config); Files.deleteIfExists(root);}
	}

	/** Checks the 2x/3x synthetic surface without rerunning gene calling or model inference. */
	private static void compareReports(final String scalarPath, final String tablePath, final String screenPath,final String calibrationConfig){
		final String[] scalar=readRow(scalarPath), table=readRow(tablePath);
		check(scalar.length==table.length, "Scalar/table column counts differ");
		for(int i=0; i<scalar.length-1; i++){
			if(i!=8 && i!=9){check(scalar[i].equals(table[i]), "Calibration changed biological column "+i);}
		}
		// Scaled columns are doubles; raw Float.toString heads are shorter decimals.
		// Use the scalar double values so text rounding cannot invent a discrepancy.
		final double comp=Double.parseDouble(table[8]), contam=Double.parseDouble(table[9]);
		final MagQCErrorCalibration fitted=calibrationConfig==null ? null : MagQCErrorCalibration.fromOptions(
			MagQCAssemblyBatch.parseOptions(new String[]{"config="+calibrationConfig}));
		// Value-head text is Float.toString: parse it back to the original float coordinates.
		final double expectedComp=fitted==null ? Double.parseDouble(scalar[8])*2 : fitted.scaleComp(
			Float.parseFloat(scalar[2]),Float.parseFloat(scalar[3]),Double.parseDouble(scalar[8]));
		final double expectedContam=fitted==null ? Double.parseDouble(scalar[9])*3 : fitted.scaleContam(
			Float.parseFloat(scalar[2]),Float.parseFloat(scalar[3]),Double.parseDouble(scalar[9]));
		close(comp, expectedComp); close(contam, expectedContam);
		final String expectedC=String.format(Locale.ROOT, "Completeness: %.2f +-%.2f",
			100.0*Float.parseFloat(table[2]), 100*expectedComp);
		final String expectedT=String.format(Locale.ROOT, "Contamination: %.2f +-%.2f",
			100.0*Float.parseFloat(table[3]), 100*expectedContam);
		final ByteFile input=ByteFile.makeByteFile(screenPath, true);
		int c=0, t=0;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				final String text=new String(line, StandardCharsets.UTF_8).replaceAll(" +", " ");
				if(text.startsWith("Completeness:")){check(text.equals(expectedC), "Completeness screen/TSV mismatch"); c++;}
				if(text.startsWith("Contamination:")){check(text.equals(expectedT), "Contamination screen/TSV mismatch"); t++;}
			}
		}finally{check(!input.close(), "Could not close screen report");}
		check(c==1 && t==1, "Single-bin fixture requires exactly one screen estimate for each target");
		System.out.println("MagQCErrorCalibrationReportComparison PASS: unchanged raw heads, "+(fitted==null ? "2x/3x" : "fitted-table")+" errors and screen/TSV parity");
	}

	/** Requires one v4 row and explicit error-column positions before interpreting fixture numbers. */
	private static String[] readRow(final String path){
		final ByteFile input=ByteFile.makeByteFile(path, true);
		final LineParser1 parser=new LineParser1('\t');
		String[] row=null;
		int columns=-1;
		boolean schema=false;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				parser.set(line);
				if(parser.termEquals("#schema_version", 0)){
					check(!schema && parser.terms()==2 && parser.termEquals("prokcc_assembly_scores_v4", 1), "Unexpected report schema");
					schema=true;
				}else if(parser.termEquals("#columns", 0)){
					check(columns<0 && parser.terms()>11, "Missing or repeated report columns");
					check(parser.termEquals("gene_completeness", 3) && parser.termEquals("gene_contamination", 4)
						&& parser.termEquals("scaled_error_gene_completeness", 9)
						&& parser.termEquals("scaled_error_gene_contamination", 10)
						&& parser.termEquals("bin_worker_wall_seconds", parser.terms()-1), "Report column layout changed");
					columns=parser.terms()-1;
				}else if(line.length>0 && line[0]!='#'){
					check(schema && row==null && parser.terms()==columns, "Fixture requires one complete report row after its header");
					row=new String[columns];
					for(int i=0; i<columns; i++){row[i]=parser.parseString(i);}
				}
			}
		}finally{check(!input.close(), "Could not close report: "+path);}
		check(row!=null, "Missing fixture report row: "+path);
		return row;
	}

	/** A synthetic constant or one-axis grid, never a biological fit. */
	private static void checkDegenerate(Path table, String comp, String contam, String nodes, double c, double t) throws Exception{
		final MagQCErrorCalibration grid=MagQCErrorCalibration.fromOptions(write(table, surface(defaults(), comp, contam, nodes)));
		close(grid.scaleComp(.5, .5, 1), c); close(grid.scaleContam(.5, .5, 1), t);
	}
	/** Uses production option defaults so table/scalar conflict handling is actually exercised. */
	private static HashMap<String,String> defaults(){
		final HashMap<String,String> options=MagQCAssemblyBatch.parseOptions(new String[0]);
		for(String key:MagQCErrorCalibration.BINDING_KEYS){options.put(key, PIN);}
		return options;
	}
	/** Writes new fixture bytes and binds precisely those bytes, including malformed-content tests. */
	private static HashMap<String,String> write(Path path, String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII));
		final HashMap<String,String> options=defaults();
		options.put("errorcalibration", path.toString()); options.put("errorcalibrationsha80", MagQCTextResource.sha80(path.toString()));
		return options;
	}
	/** Small hand-built schema, explicitly marked FIXTURE_ONLY and not eligible for release. */
	private static String surface(HashMap<String,String> options, String comp, String contam, String nodes){
		final StringBuilder out=new StringBuilder("#schema\t"+MagQCErrorCalibration.SCHEMA+"\n");
		for(String key:MagQCErrorCalibration.BINDING_KEYS){out.append('#').append(key).append('\t').append(options.get(key)).append('\n');}
		out.append("#layout\t").append(MagQCVectorMaker.JOINT_COMPOSITE_LAYOUT).append('\n');
		out.append("#output_names\t").append(MagQCErrorCalibration.OUTPUT_NAMES).append('\n');
		out.append("#units\tfraction\n#truth\tgene\n#truthpolicy\tfixture\n#interpolation\tbilinear_multiplier\n");
		out.append("#truthmanifestsha80\t").append(PIN).append("\n#splitmanifestsha80\t").append(PIN).append('\n');
		out.append("#fitrecipe\tfixture\n#smoothing\tnone\n#errorfitset\tFIXTURE_ONLY\n#errorfitdate\t2026-10-03\n#errorcoverage\t0.95\n");
		out.append("#comp_knots\t").append(comp).append("\n#contam_knots\t").append(contam).append('\n');
		out.append("#columns\tcomp_index\tcontam_index\tcomp_multiplier\tcontam_multiplier\n").append(nodes);
		return out.toString();
	}
	private static void reject(HashMap<String,String> options){
		boolean failed=false;
		try{MagQCErrorCalibration.fromOptions(options);}catch(RuntimeException e){failed=true;}
		check(failed, "Malformed calibration configuration was accepted");
	}
	private static void rejectScale(MagQCErrorCalibration calibration, double c, double t, double error){
		boolean failed=false;
		try{calibration.scaleComp(c, t, error);}catch(IllegalArgumentException e){failed=true;}
		check(failed, "Invalid/overflowing error estimate was accepted");
	}
	private static void close(double actual, double expected){
		check(Double.isFinite(actual) && Math.abs(actual-expected)<=1e-12, "Expected "+expected+" but got "+actual);
	}
	/** Fixture checks remain loud without JVM assertions. */
	private static void check(boolean pass, String message){if(!pass){throw new AssertionError(message);}}
	private static final String PIN="12345678901234567890";
	private static final String NODES="0\t0\t2\t1\n0\t1\t4\t5\n0\t2\t16\t9\n1\t0\t6\t3\n1\t1\t12\t7\n1\t2\t20\t19\n";
}
