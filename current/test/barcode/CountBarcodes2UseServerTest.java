package test.barcode;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.Locale;

import barcode.CountBarcodes2;
import barcode.PCRMatrix;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import server.ServerTools;
import shared.Shared;

/** Verifies that explicit useserver=f survives the unloaded-PROB fallback. */
public final class CountBarcodes2UseServerTest{

	/** Requires mode=baseline|fixed, case=..., and a new dir=path. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, CountBarcodes2UseServerTest.class, false);
		try{
			String mode=null, testCase=null, directory=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(key.equals("case")){testCase=value;}
				else if(key.equals("dir")){directory=value;}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			if(!"baseline".equals(mode) && !"fixed".equals(mode)){
				throw new IllegalArgumentException("Set mode=baseline or mode=fixed.");
			}
			if(!"explicit-f".equals(testCase) && !"omitted".equals(testCase) &&
					!"explicit-t".equals(testCase) && !"hdist-f".equals(testCase)){
				throw new IllegalArgumentException("Set case=explicit-f, omitted, explicit-t, or hdist-f.");
			}
			if(directory==null){throw new IllegalArgumentException("Set dir to a nonexistent output directory.");}

			final Path dir=Files.createDirectory(Paths.get(directory));
			Files.write(dir.resolve("counts.tsv"), counts(testCase).getBytes(StandardCharsets.US_ASCII),
					StandardOpenOption.CREATE_NEW);
			Files.write(dir.resolve("expected.txt"), "ACGT\n".getBytes(StandardCharsets.US_ASCII),
					StandardOpenOption.CREATE_NEW);

			check(!probLoaded(), "The use-server fixture expects probability matrices to start unloaded.");
			final RunResult result=runCommand(dir, testCase);
			Files.write(dir.resolve("stderr.txt"), result.stderr.getBytes(StandardCharsets.UTF_8),
					StandardOpenOption.CREATE_NEW);

			if("hdist-f".equals(testCase)){
				check(result.thrown==null, "Hamming useserver=f should finish locally; thrown="+result.thrown+
						"\nstderr:\n"+result.stderr);
				check(!result.client, "Hamming useserver=f must not enter DemuxClient; stderr:\n"+
						result.stderr);
				final byte[] observed=Files.readAllBytes(dir.resolve("map.tsv"));
				final byte[] expected="ACGT\tACGT\t2\t0\n".getBytes(StandardCharsets.US_ASCII);
				check(Arrays.equals(expected, observed), "Unexpected Hamming map: "+
						new String(observed, StandardCharsets.US_ASCII));
			}else if("fixed".equals(mode) && "explicit-f".equals(testCase)){
				check(!result.client, "Fixed CountBarcodes2 must reject explicit useserver=f before " +
						"entering DemuxClient; stderr:\n"+result.stderr);
				check(result.explicitReject(), "Fixed CountBarcodes2 must reject unloaded probability " +
						"mode locally; thrown="+result.thrown);
				check(!Files.exists(dir.resolve("map.tsv")), "Rejected run must not create mapout.");
			}else{
				check(result.client, testCase+" must enter DemuxClient; stderr:\n"+result.stderr);
				check(result.nullResponseInGetMap(), testCase+" must die at DemuxClient.getMap " +
						"after the malformed URL returns no server response; thrown="+result.thrown+
						"\nstderr:\n"+result.stderr);
				check(!Files.exists(dir.resolve("map.tsv")), "Failed client run must not create mapout.");
			}

			System.out.println("COUNT_BARCODES2_USESERVER mode="+mode+" case="+testCase+
					" client="+result.client+" thrown="+
					(result.thrown==null ? "none" : result.thrown.getClass().getName()));
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** No instances; each fixture call owns one scratch directory and one JVM. */
	private CountBarcodes2UseServerTest(){}

	private static String[] command(final Path dir, final String testCase){
		final String[] prefix=new String[] {
				"countsin="+dir.resolve("counts.tsv"),
				"expected="+dir.resolve("expected.txt"),
				"mapout="+dir.resolve("map.tsv"),
				"type="+("hdist-f".equals(testCase) ? "hdist" : "prob"),
				"ow=f",
		};
		final String useServerFlag=useServerFlag(testCase);
		if(useServerFlag==null){return prefix;}
		return new String[] {prefix[0], prefix[1], prefix[2], prefix[3], useServerFlag, prefix[4]};
	}

	private static String counts(final String testCase){
		return "hdist-f".equals(testCase) ? "ACGT\t2\n" : "ACGT\t2\nTGCA\t1\n";
	}

	private static String useServerFlag(final String testCase){
		if("omitted".equals(testCase)){return null;}
		return "explicit-t".equals(testCase) ? "useserver=t" : "useserver=f";
	}

	private static RunResult runCommand(final Path dir, final String testCase) throws Exception{
		final Field demuxServer=Shared.class.getDeclaredField("demuxServer");
		demuxServer.setAccessible(true);
		final Object oldDemuxServer=demuxServer.get(null);
		final boolean oldSuppressErrors=ServerTools.suppressErrors;
		final PrintStream oldErr=System.err;
		final ByteArrayOutputStream stderr=new ByteArrayOutputStream();
		final PrintStream capture=new PrintStream(stderr, true, "UTF-8");
		Throwable thrown=null;
		try{
			demuxServer.set(null, "not-a-url/");
			ServerTools.suppressErrors=true;
			System.setErr(capture);
			CountBarcodes2.main(command(dir, testCase));
		}catch(final Throwable t){
			thrown=t;
		}finally{
			System.setErr(oldErr);
			ServerTools.suppressErrors=oldSuppressErrors;
			demuxServer.set(null, oldDemuxServer);
			capture.close();
		}
		return new RunResult(thrown, stderr.toString("UTF-8"));
	}

	private static boolean probLoaded() throws Exception{
		final Field field=PCRMatrix.class.getDeclaredField("probLoaded");
		field.setAccessible(true);
		return field.getBoolean(null);
	}

	/** Always enforce fixture oracles independent of assertion configuration. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Captured process result for one CountBarcodes2 main invocation. */
	private static final class RunResult{

		RunResult(final Throwable thrown_, final String stderr_){
			thrown=thrown_;
			stderr=stderr_;
			client=stderr.contains("Using client-server mode for barcode analysis.");
		}

		boolean explicitReject(){
			return thrown instanceof IllegalStateException &&
					thrown.getMessage()!=null &&
					thrown.getMessage().contains("useserver=f") &&
					thrown.getMessage().contains("probability");
		}

		boolean nullResponseInGetMap(){
			if(!(thrown instanceof NullPointerException)){return false;}
			for(final StackTraceElement ste : thrown.getStackTrace()){
				if("barcode.DemuxClient".equals(ste.getClassName()) &&
						"getMap".equals(ste.getMethodName())){
					return true;
				}
			}
			return false;
		}

		final Throwable thrown;
		final String stderr;
		final boolean client;
	}
}
