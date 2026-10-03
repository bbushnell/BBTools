package prot;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;

import fileIO.ByteFile;
import parse.LineParser1;
import structures.ByteBuilder;

/**
 * Independently verifies every header, row, gate, and digest in one sealed
 * schema-7 family-threshold artifact.
 *
 * @author Yoimiya
 */
public final class Schema7ThresholdVerifier {

	private Schema7ThresholdVerifier(){}

	/** Parses the strict command line and verifies the complete artifact. */
	public static void main(final String[] args){
		final Result result=run(parseArgs(args));
		System.out.println("[schema7-threshold-verifier] COMPLETE in="+result.in+
			" families="+result.families+" members="+result.members+
			" table_sha80="+result.tableSha80);
	}

	/** Verifies the input before and after parsing and returns its totals. */
	static Result run(final Config config){
		checkHash(config.in, config.inSha80, "INPUT");
		final HashMap<String,String> header=new HashMap<String,String>(64);
		final HashSet<Integer> familyIds=new HashSet<Integer>(config.expectedFamilies*2);
		final HashSet<String> repIds=new HashSet<String>(config.expectedFamilies*2);
		final ByteBuilder table=new ByteBuilder(1<<20);
		final ByteFile input=ByteFile.makeByteFile(config.in, false);
		final LineParser1 parser=new LineParser1((byte)'\t');
		boolean columns=false;
		int rank=0;
		long members=0;
		Throwable failure=null;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				if(line.length==0){fail("EMPTY_LINE", rank);}
				if(line[0]=='#'){
					if(columns){fail("LATE_HEADER", rank);}
					parseHeader(line, header);
					continue;
				}
				if(!columns){
					requireColumns(line);
					table.append(line).nl();
					columns=true;
					continue;
				}
				table.append(line).nl();
				parser.set(line);
				if(parser.terms()!=33){fail("ROW_FIELDS", parser.terms());}
				final int observedRank=parseInt(parser.parseString(0), "rank", rank);
				final int familyId=parseInt(parser.parseString(1), "family_id", rank);
				final String repId=parser.parseString(2);
				final int n=parseInt(parser.parseString(3), "n", rank);
				final String status=parser.parseString(4);
				final boolean enabled=parseBoolean(parser.parseString(5),
					"acceptance_enabled", rank);
				if(observedRank!=rank || familyId<0 || repId.length()==0 || n<1 ||
						!familyIds.add(familyId) || !repIds.add(repId) ||
						enabled!=acceptanceEnabled(status, n, rank)){
					fail("ROW_IDENTITY", rank);
				}
				checkIntGate(parser, 6, "raw_score", rank, false);
				checkFloatGate(parser, 9, "identity", rank, 0, 100,
					IDENTITY_FLOOR);
				checkFloatGate(parser, 12, "R", rank,
					Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, R_FLOOR);
				checkFloatGate(parser, 15, "mutual_overlap", rank, 0, 1,
					MUTUAL_FLOOR);
				final int lengthLo=checkIntGate(parser, 18, "length_lo", rank, true);
				final int lengthHi=checkIntGate(parser, 21, "length_hi", rank, true);
				if(lengthLo<1 || lengthHi<lengthLo){fail("LENGTH_ORDER", rank);}
				checkIntGate(parser, 24, "kmer_count", rank, true);
				checkDoubleGate(parser, 27, "kmer_density", rank, 0, 1);
				checkFloatGate(parser, 30, "hbm_path", rank, 0, 1,
					HBM_PATH_FLOOR);
				members=Math.addExact(members, n);
				rank++;
			}
		}catch(RuntimeException e){failure=e; throw e;}
		catch(Error e){failure=e; throw e;}
		finally{closeInput(input, failure);}
		if(!columns || rank!=config.expectedFamilies ||
				members!=config.expectedMembers){
			throw new RuntimeException("SCHEMA7_VERIFY_TOTALS: columns="+columns+
				" families="+rank+"/"+config.expectedFamilies+" members="+
				members+"/"+config.expectedMembers);
		}
		verifyHeader(header, config);
		final String tableSha80=DigestSuffix.bytes(table.toBytes());
		if(!config.tableSha80.equals(tableSha80) ||
				!tableSha80.equals(header.get("threshold_table_sha80"))){
			throw new RuntimeException("SCHEMA7_VERIFY_TABLE_HASH: "+tableSha80);
		}
		checkHash(config.in, config.inSha80, "INPUT_CHANGED");
		return new Result(config.in, rank, members, tableSha80);
	}

	/** Parses one unique comment header. */
	private static void parseHeader(final byte[] line,
			final HashMap<String,String> header){
		int tab=-1;
		for(int i=1; i<line.length; i++){
			if(line[i]=='\t'){tab=i; break;}
		}
		if(tab<2 || tab==line.length-1){fail("HEADER_FORMAT", header.size());}
		final String key=new String(line, 1, tab-1, StandardCharsets.US_ASCII);
		final String value=new String(line, tab+1, line.length-tab-1,
			StandardCharsets.US_ASCII);
		if(header.put(key, value)!=null){
			throw new RuntimeException("SCHEMA7_VERIFY_DUPLICATE_HEADER: "+key);
		}
	}

	/** Requires the exact schema-7 column order. */
	private static void requireColumns(final byte[] line){
		if(!COLUMN_HEADER.equals(new String(line, StandardCharsets.US_ASCII))){
			throw new RuntimeException("SCHEMA7_VERIFY_COLUMN_HEADER");
		}
	}

	/** Verifies the complete, closed schema-7 metadata header. */
	private static void verifyHeader(final HashMap<String,String> header,
			final Config config){
		if(header.size()!=32){
			throw new RuntimeException("SCHEMA7_VERIFY_HEADER_COUNT: "+header.size());
		}
		requireHeader(header, "artifact_type", "schema7_family_thresholds");
		requireHeader(header, "schema", "7");
		requireHeader(header, "families", Integer.toString(config.expectedFamilies));
		requireHeader(header, "members", Long.toString(config.expectedMembers));
		requireHeader(header, "roster_sha80", config.rosterSha80);
		requireHeader(header, "base_landmarks_sha80", config.baseLandmarksSha80);
		requireHeader(header, "base_provenance_sha80", config.baseProvenanceSha80);
		requireHeader(header, "length_bounds_sha80", config.lengthBoundsSha80);
		requireHeader(header, "final_empirical_sha80", config.finalEmpiricalSha80);
		requireHeader(header, "consensus_sha80", config.consensusSha80);
		requireHeader(header, "core_sha80", config.coreSha80);
		requireHeader(header, "covering_sets_sha80", config.coveringSetsSha80);
		requireHeader(header, "sidecar_sha80", config.sidecarSha80);
		requireHeader(header, "role_manifest_sha80", config.roleManifestSha80);
		requireHeader(header, "hbm_bundle_sha80", config.hbmBundleSha80);
		requireHeader(header, "hbm_provenance_sha80", config.hbmProvenanceSha80);
		requireHeader(header, "hbm_score_contract", HbmScoreContract.VALUE);
		requireHeader(header, "acceptance_enabled_rule",
			"all_roster_rows_compete_v1");
		requireHeader(header, "d125_output_manifest_sha80", config.d125OutputSha80);
		requireHeader(header, "d125_hbm_path_curve_sha80", config.d125CurveSha80);
		requireHeader(header, "d125_floor_summary_sha80", config.d125SummarySha80);
		requireHeader(header, "d127_identity_curve_sha80",
			config.d127IdentityCurveSha80);
		requireHeader(header, "d127_r_curve_sha80", config.d127RCurveSha80);
		requireHeader(header, "d127_comparison_manifest_sha80",
			config.d127ComparisonSha80);
		requireHeader(header, "identity_global_floor", IDENTITY_FLOOR_TEXT);
		requireHeader(header, "R_global_floor", R_FLOOR_TEXT);
		requireHeader(header, "mutual_overlap_global_floor", MUTUAL_FLOOR_TEXT);
		requireHeader(header, "hbm_path_global_floor", HBM_PATH_FLOOR_TEXT);
		requireHeader(header, "empirical_targets", EMPIRICAL_TARGETS);
		requireHeader(header, "effective_rule",
			"max_empirical_global_where_defined");
		requireHeader(header, "threshold_table_sha80_contract",
			"column_header_lf_plus_rows_each_lf");
		requireHeader(header, "threshold_table_sha80", config.tableSha80);
	}

	/** Requires one exact metadata value. */
	private static void requireHeader(final HashMap<String,String> header,
			final String key, final String expected){
		final String observed=header.get(key);
		if(!expected.equals(observed)){
			throw new RuntimeException("SCHEMA7_VERIFY_HEADER_"+key+": "+observed+
				" != "+expected);
		}
	}

	/** Verifies one integer empirical/effective/no-floor triple. */
	private static int checkIntGate(final LineParser1 parser, final int start,
			final String label, final int rank, final boolean nonnegative){
		final int empirical=parseInt(parser.parseString(start), label, rank);
		final int effective=parseInt(parser.parseString(start+1), label, rank);
		final boolean bound=parseBoolean(parser.parseString(start+2), label, rank);
		if(empirical!=effective || bound || (nonnegative && empirical<0)){
			fail("INT_GATE_"+label, rank);
		}
		return empirical;
	}

	/** Verifies one float empirical/effective/global-floor triple. */
	private static void checkFloatGate(final LineParser1 parser, final int start,
			final String label, final int rank, final double min, final double max,
			final float floor){
		final float empirical=parseFloat(parser.parseString(start), label, rank);
		final float effective=parseFloat(parser.parseString(start+1), label, rank);
		final boolean bound=parseBoolean(parser.parseString(start+2), label, rank);
		if(empirical<min || empirical>max){fail("RANGE_"+label, rank);}
		final float expected=Math.max(empirical, floor);
		if(Float.compare(effective, expected)!=0 || bound!=(expected>empirical)){
			fail("FLOAT_GATE_"+label, rank);
		}
	}

	/** Verifies one bounded double empirical/effective/no-floor triple. */
	private static void checkDoubleGate(final LineParser1 parser, final int start,
			final String label, final int rank, final double min, final double max){
		final double empirical=parseDouble(parser.parseString(start), label, rank);
		final double effective=parseDouble(parser.parseString(start+1), label, rank);
		final boolean bound=parseBoolean(parser.parseString(start+2), label, rank);
		if(empirical<min || empirical>max ||
				Double.compare(empirical, effective)!=0 || bound){
			fail("DOUBLE_GATE_"+label, rank);
		}
	}

	/** Validates the D122 vocabulary; every roster row remains a competitor. */
	private static boolean acceptanceEnabled(final String status, final int n,
			final int rank){
		if("KEEP".equals(status) && n>=10){return true;}
		if("KEEP_FLAG_LOW_N".equals(status) && n>=2 && n<10){return true;}
		if("EXCLUDE_LOW_N".equals(status) && n>=2 && n<=3){return true;}
		if("EXCLUDE_SINGLETON_STANDING_RULE".equals(status) && n==1){return true;}
		fail("STATUS", rank);
		return false;
	}

	/** Parses one exact boolean token. */
	private static boolean parseBoolean(final String text, final String label,
			final int rank){
		if("true".equals(text)){return true;}
		if("false".equals(text)){return false;}
		throw new RuntimeException("SCHEMA7_VERIFY_BOOLEAN_"+label+": rank="+
			rank+" value="+text);
	}

	/** Parses one finite float. */
	private static float parseFloat(final String text, final String label,
			final int rank){
		final float value;
		try{value=Float.parseFloat(text);}
		catch(NumberFormatException e){
			throw new RuntimeException("SCHEMA7_VERIFY_FLOAT_"+label+": rank="+
				rank+" value="+text, e);
		}
		if(!Float.isFinite(value)){fail("NONFINITE_"+label, rank);}
		if(!text.equals(Float.toString(value))){fail("NONCANONICAL_"+label, rank);}
		return value;
	}

	/** Parses one finite double. */
	private static double parseDouble(final String text, final String label,
			final int rank){
		final double value;
		try{value=Double.parseDouble(text);}
		catch(NumberFormatException e){
			throw new RuntimeException("SCHEMA7_VERIFY_DOUBLE_"+label+": rank="+
				rank+" value="+text, e);
		}
		if(!Double.isFinite(value)){fail("NONFINITE_"+label, rank);}
		if(!text.equals(Double.toString(value))){fail("NONCANONICAL_"+label, rank);}
		return value;
	}

	/** Parses one canonical decimal integer independently of assertion mode. */
	private static int parseInt(final String text, final String label,
			final int rank){
		final int value;
		try{value=Integer.parseInt(text);}
		catch(NumberFormatException e){
			throw new RuntimeException("SCHEMA7_VERIFY_INTEGER_"+label+": rank="+
				rank+" value="+text, e);
		}
		if(!text.equals(Integer.toString(value))){
			fail("NONCANONICAL_"+label, rank);
		}
		return value;
	}

	/** Preserves an in-flight parse failure while reporting deferred I/O errors. */
	private static void closeInput(final ByteFile input,
			final Throwable priorFailure){
		if(!input.close()){return;}
		final RuntimeException closeFailure=new RuntimeException(
			"SCHEMA7_VERIFY_INPUT_IO: deferred input error");
		if(priorFailure==null){throw closeFailure;}
		priorFailure.addSuppressed(closeFailure);
	}

	/** Requires one file to match its sha80. */
	private static void checkHash(final String path, final String expected,
			final String label){
		final String observed=MagQCTextResource.sha80(path);
		if(!expected.equals(observed)){
			throw new RuntimeException("SCHEMA7_VERIFY_"+label+"_HASH: "+observed+
				" != "+expected);
		}
	}

	/** Throws one stable verifier diagnostic. */
	private static void fail(final String label, final int rank){
		throw new RuntimeException("SCHEMA7_VERIFY_"+label+": rank="+rank);
	}

	/** Parses the strict BBTools-style verifier command line. */
	static Config parseArgs(final String[] args){
		final HashSet<String> known=new HashSet<String>(Arrays.asList(
			"in", "insha80", "expectedfamilies", "expectedmembers", "tablesha80",
			"rostersha80", "baselandmarkssha80", "baseprovenancesha80",
			"lengthboundssha80", "finalempiricalsha80", "consensussha80",
			"coresha80", "coveringsetssha80", "sidecarsha80", "hbmbundlesha80",
			"rolemanifestsha80", "hbmprovenancesha80", "d125outputsha80", "d125curvesha80",
			"d125summarysha80", "d127identitycurvesha80", "d127rcurvesha80",
			"d127comparisonsha80"));
		final HashMap<String,String> map=new HashMap<String,String>();
		for(final String arg : args){
			final int equals=arg.indexOf('=');
			if(equals<1){throw new IllegalArgumentException("SCHEMA7_VERIFY_ARG: "+arg);}
			final String key=arg.substring(0, equals).toLowerCase(Locale.ROOT);
			if(!known.contains(key)){
				throw new IllegalArgumentException("SCHEMA7_VERIFY_UNKNOWN_ARG: "+key);
			}
			if(map.put(key, arg.substring(equals+1))!=null){
				throw new IllegalArgumentException("SCHEMA7_VERIFY_DUPLICATE_ARG: "+key);
			}
		}
		final Config config=new Config();
		config.in=require(map, "in");
		config.inSha80=hash(map, "insha80");
		config.expectedFamilies=positiveInt(map, "expectedfamilies");
		config.expectedMembers=positiveLong(map, "expectedmembers");
		config.tableSha80=hash(map, "tablesha80");
		config.rosterSha80=hash(map, "rostersha80");
		config.baseLandmarksSha80=hash(map, "baselandmarkssha80");
		config.baseProvenanceSha80=hash(map, "baseprovenancesha80");
		config.lengthBoundsSha80=hash(map, "lengthboundssha80");
		config.finalEmpiricalSha80=hash(map, "finalempiricalsha80");
		config.consensusSha80=hash(map, "consensussha80");
		config.coreSha80=hash(map, "coresha80");
		config.coveringSetsSha80=hash(map, "coveringsetssha80");
		config.sidecarSha80=hash(map, "sidecarsha80");
		config.roleManifestSha80=hash(map, "rolemanifestsha80");
		config.hbmBundleSha80=hash(map, "hbmbundlesha80");
		config.hbmProvenanceSha80=hash(map, "hbmprovenancesha80");
		config.d125OutputSha80=hash(map, "d125outputsha80");
		config.d125CurveSha80=hash(map, "d125curvesha80");
		config.d125SummarySha80=hash(map, "d125summarysha80");
		config.d127IdentityCurveSha80=hash(map, "d127identitycurvesha80");
		config.d127RCurveSha80=hash(map, "d127rcurvesha80");
		config.d127ComparisonSha80=hash(map, "d127comparisonsha80");
		return config;
	}

	/** Returns one required nonempty argument. */
	private static String require(final HashMap<String,String> map,
			final String key){
		final String value=map.get(key);
		if(value==null || value.length()==0){
			throw new IllegalArgumentException("SCHEMA7_VERIFY_MISSING_ARG: "+key);
		}
		return value;
	}

	/** Returns one required sha80 argument. */
	private static String hash(final HashMap<String,String> map,
			final String key){
		return DigestSuffix.requireSuffix(require(map, key), key);
	}

	/** Parses one positive integer argument. */
	private static int positiveInt(final HashMap<String,String> map,
			final String key){
		final int value=Integer.parseInt(require(map, key));
		if(value<1){throw new IllegalArgumentException("SCHEMA7_VERIFY_POSITIVE_"+key);}
		return value;
	}

	/** Parses one positive long argument. */
	private static long positiveLong(final HashMap<String,String> map,
			final String key){
		final long value=Long.parseLong(require(map, key));
		if(value<1){throw new IllegalArgumentException("SCHEMA7_VERIFY_POSITIVE_"+key);}
		return value;
	}

	/** Complete verifier configuration. */
	static final class Config {
		String in, inSha80, tableSha80, rosterSha80, baseLandmarksSha80;
		String baseProvenanceSha80, lengthBoundsSha80, finalEmpiricalSha80;
		String consensusSha80, coreSha80, coveringSetsSha80, sidecarSha80;
		String roleManifestSha80;
		String hbmBundleSha80, hbmProvenanceSha80, d125OutputSha80;
		String d125CurveSha80, d125SummarySha80, d127IdentityCurveSha80;
		String d127RCurveSha80, d127ComparisonSha80;
		int expectedFamilies;
		long expectedMembers;
	}

	/** Immutable verification receipt. */
	static final class Result {
		final String in, tableSha80;
		final int families;
		final long members;
		Result(final String in_, final int families_, final long members_,
				final String tableSha80_){
			in=in_;
			families=families_;
			members=members_;
			tableSha80=tableSha80_;
		}
	}

	private static final float IDENTITY_FLOOR=36.177475f;
	private static final float R_FLOOR=0.28840125f;
	private static final float MUTUAL_FLOOR=0.80f;
	private static final float HBM_PATH_FLOOR=0.51504076f;
	private static final String IDENTITY_FLOOR_TEXT=Float.toString(IDENTITY_FLOOR);
	private static final String R_FLOOR_TEXT=Float.toString(R_FLOOR);
	private static final String MUTUAL_FLOOR_TEXT=Float.toString(MUTUAL_FLOOR);
	private static final String HBM_PATH_FLOOR_TEXT=Float.toString(HBM_PATH_FLOOR);
	private static final String EMPIRICAL_TARGETS=
		"raw=0.99;identity=0.99;R=0.99;mutual_overlap=0.99;"+
		"length=middle_0.998;kmer_count=0.995;kmer_density=0.995;hbm_path=0.99";
	private static final String COLUMN_HEADER=
		"rank\tfamily_id\trep_id\tn\tcalibration_status\tacceptance_enabled"+
		"\traw_score_empirical\traw_score_effective\traw_score_floor_bound"+
		"\tidentity_empirical\tidentity_effective\tidentity_floor_bound"+
		"\tR_empirical\tR_effective\tR_floor_bound"+
		"\tmutual_overlap_empirical\tmutual_overlap_effective"+
		"\tmutual_overlap_floor_bound\tlength_lo_empirical"+
		"\tlength_lo_effective\tlength_lo_floor_bound\tlength_hi_empirical"+
		"\tlength_hi_effective\tlength_hi_floor_bound\tkmer_count_empirical"+
		"\tkmer_count_effective\tkmer_count_floor_bound"+
		"\tkmer_density_empirical\tkmer_density_effective"+
		"\tkmer_density_floor_bound\thbm_path_empirical\thbm_path_effective"+
		"\thbm_path_floor_bound";
}
