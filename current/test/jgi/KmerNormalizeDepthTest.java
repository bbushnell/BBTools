package test.jgi;

import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Locale;

import jgi.KmerNormalize;
import parse.Parse;
import shared.Shared;

/**
 * Public parseDepth fixture for KmerNormalize depth/error field placement.
 * @author Jean
 */
public final class KmerNormalizeDepthTest{

	/** Requires a new {@code outdir=...}; {@code phase=baseline|fixed} is mandatory. */
	public static void main(final String[] args) throws Exception{
		String outdir=null, phase=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			final String key=pair[0].toLowerCase(Locale.ROOT);
			if(key.equals("outdir")){outdir=pair[1];}
			else if(key.equals("phase")){phase=pair[1].toLowerCase(Locale.ROOT);}
			else{throw new IllegalArgumentException("Unknown fixture option: "+arg);}
		}
		require(outdir!=null && !outdir.isEmpty(), "Set outdir=<new case directory>.");
		require(phase!=null && !phase.isEmpty(), "Set phase=baseline or phase=fixed.");
		require(phase.equals("baseline") || phase.equals("fixed"), "Set phase=baseline or phase=fixed.");

		if(!Shared.EA()){throw new AssertionError("Assertions must be enabled.");}
		Shared.setThreads(1);

		new KmerNormalizeDepthTest(Paths.get(outdir).toAbsolutePath().normalize(), phase).run();
	}

	/** Runs all six public-helper probes and writes exact observations. */
	private void run() throws Exception{
		Files.createDirectory(root);
		report=root.resolve("observations.tsv");
		allowedEntries.add(root.toAbsolutePath().normalize());
		allowedEntries.add(report.toAbsolutePath().normalize());
		write(report, "case\tphase\tinput\tdesired_fields\texpected_fields\tactual_fields\t"+
			"expected_exception\tactual_exception\tactual_message\texpected_top\tactual_top\tmatch\n");

		boolean all=true;
		for(final CaseSpec spec : CASES){
			final Observation observation=observe(spec);
			final String expectedFailure=expectedFailure(spec);
			final String expectedTop=expectedFailure.isEmpty() ? "" : EXPECTED_TOP;
			final boolean match=matches(spec, observation, expectedFailure, expectedTop);
			all&=match;
			append(report, spec.name+'\t'+phase+'\t'+scalar(spec.input)+'\t'+array(spec.desired)+'\t'+
				array(expectedFields(spec))+'\t'+observation.fields()+'\t'+
				scalar(expectedFailure)+'\t'+scalar(type(observation.failure))+'\t'+
				scalar(message(observation.failure))+'\t'+scalar(expectedTop)+'\t'+
				scalar(top(observation.failure))+'\t'+match+'\n');
		}

		final ArrayList<String> unexpectedEntries=unexpectedEntries();
		final boolean sideEffectMatch=unexpectedEntries.isEmpty();
		all&=sideEffectMatch;
		append(report, "<side-effects>\t"+phase+"\t\t[]\t[]\t"+unexpectedEntries+
			"\t\t\t\t\t\t"+sideEffectMatch+'\n');
		if(!all){throw new AssertionError("Unexpected parseDepth observation; inspect "+report);}
		System.out.println("KMER_NORMALIZE_DEPTH_PASS phase="+phase);
	}

	/** Calls only the public parseDepth helper with a newly allocated output array. */
	private static Observation observe(final CaseSpec spec){
		try{
			return new Observation(KmerNormalize.parseDepth(spec.input, null), null);
		}catch(final Throwable caught){
			return new Observation(null, caught);
		}
	}

	/** Returns whether the actual helper outcome matches the phase-specific expected outcome. */
	private boolean matches(final CaseSpec spec, final Observation observation,
			final String expectedFailure, final String expectedTop){
		if(expectedFailure.isEmpty()){
			return observation.failure==null && Arrays.equals(expectedFields(spec), observation.fields);
		}
		return expectedFailure.equals(type(observation.failure)) &&
			expectedTop.equals(top(observation.failure));
	}

	/** Returns the phase-specific fields expected from the current public helper. */
	private int[] expectedFields(final CaseSpec spec){
		return phase.equals("baseline") ? spec.baselineFields : spec.desired;
	}

	/** Returns the phase-specific exception expected from the current public helper. */
	private String expectedFailure(final CaseSpec spec){
		return phase.equals("baseline") && spec.baselineFailure!=null ? spec.baselineFailure : "";
	}

	/** Writes text to a newly created report. */
	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Appends one observation to the report. */
	private static void append(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
	}

	/** Returns unexpected files or directories below the required new output directory. */
	private ArrayList<String> unexpectedEntries() throws Exception{
		final ArrayList<String> extra=new ArrayList<String>();
		collectUnexpected(root, extra);
		Collections.sort(extra);
		return extra;
	}

	/** Recursively records files this fixture did not create intentionally. */
	private void collectUnexpected(final Path path, final ArrayList<String> extra) throws Exception{
		final Path absolute=path.toAbsolutePath().normalize();
		if(!allowedEntries.contains(absolute)){extra.add(root.relativize(absolute).toString());}
		if(Files.isDirectory(path)){
			final DirectoryStream<Path> stream=Files.newDirectoryStream(path);
			try{
				for(final Path child : stream){collectUnexpected(child, extra);}
			}finally{
				stream.close();
			}
		}
	}

	/** Formats an int array for a single TSV field. */
	private static String array(final int[] fields){
		return fields==null ? "null" : Arrays.toString(fields);
	}

	/** Formats nullable scalar values for a single TSV field. */
	private static String scalar(final String value){
		return value==null ? "<java-null>" : value.replace('\t', ' ').replace('\n', ' ').replace('\r', ' ');
	}

	/** Returns the exact throwable type name, if any. */
	private static String type(final Throwable failure){
		return failure==null ? "" : failure.getClass().getName();
	}

	/** Returns the exact throwable message, if any. */
	private static String message(final Throwable failure){
		return failure==null ? "" : failure.getMessage();
	}

	/** Returns the top stack frame class and method for an exception, if any. */
	private static String top(final Throwable failure){
		if(failure==null || failure.getStackTrace().length<1){return "";}
		final StackTraceElement top=failure.getStackTrace()[0];
		return top.getClassName()+'.'+top.getMethodName();
	}

	/** Returns four int fields. */
	private static int[] fields(final int a, final int b, final int c, final int d){
		return new int[] {a, b, c, d};
	}

	/** Requires a condition to be true. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** One parseDepth probe and its baseline expectation. */
	private static final class CaseSpec{

		CaseSpec(final String name_, final String input_, final int[] desired_,
				final int[] baselineFields_, final String baselineFailure_){
			name=name_;
			input=input_;
			desired=desired_;
			baselineFields=baselineFields_;
			baselineFailure=baselineFailure_;
		}

		final String name;
		final String input;
		final int[] desired;
		final int[] baselineFields;
		final String baselineFailure;
	}

	/** One observed parseDepth result or exception. */
	private static final class Observation{

		Observation(final int[] fields_, final Throwable failure_){
			fields=fields_;
			failure=failure_;
		}

		String fields(){
			return failure==null ? array(fields) : "<exception>";
		}

		final int[] fields;
		final Throwable failure;
	}

	private static final String EXPECTED_TOP="jgi.KmerNormalize.parseDepth";

	private static final CaseSpec[] CASES={
		new CaseSpec("null", null, null, null, null),
		new CaseSpec("non_id", "not-an-id", null, null, null),
		new CaseSpec("d1", "id=7,d1=11", fields(11, -1, -1, -1), fields(11, -1, -1, -1), null),
		new CaseSpec("d1_e1", "id=7,d1=11,e1=2", fields(11, -1, 2, -1), fields(11, -1, -1, 2), null),
		new CaseSpec("d1_d2", "id=7,d1=11,d2=13 /1", fields(11, 13, -1, -1), fields(11, -1, 13, -1), null),
		new CaseSpec("d1_d2_e1_e2", "id=7,d1=11,d2=13,e1=2,e2=3 /1", fields(11, 13, 2, 3),
			null, "java.lang.ArrayIndexOutOfBoundsException"),
	};

	/** Stores the report directory and expected comparison phase. */
	private KmerNormalizeDepthTest(final Path root_, final String phase_){
		root=root_;
		phase=phase_;
	}

	private final Path root;
	private final String phase;
	private final HashSet<Path> allowedEntries=new HashSet<Path>();
	private Path report;
}
