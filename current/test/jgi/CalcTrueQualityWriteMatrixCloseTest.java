package test.jgi;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Locale;

import jgi.CalcTrueQuality;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Finite fixture for the five CalcTrueQuality.writeMatrix overloads.
 * Normal cases pin the exact emitted row shape for 1D through 5D matrices;
 * /dev/full cases capture return from the current ignored poisonAndWait path.
 * This fixture calls the static matrix writers directly and does not initialize
 * CalcTrueQuality's global calibration matrices or worker pipeline.
 * @author Jean
 */
public final class CalcTrueQualityWriteMatrixCloseTest{

	/** Runs the normal file controls and error-output checks. */
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CalcTrueQualityWriteMatrixCloseTest.class, false);
		try{
			final Args parsed=new Args(pp.args);
			new CalcTrueQualityWriteMatrixCloseTest(Paths.get(parsed.outdir), parsed.fixed).process();
		}finally{
			Shared.closeStream(pp.outstream);
		}
	}

	/** Stores validated arguments. */
	private CalcTrueQualityWriteMatrixCloseTest(final Path outdir_, final boolean fixed_){
		outdir=outdir_;
		fixed=fixed_;
	}

	/** Runs every bounded case. */
	private void process(){
		createDirectory(outdir);
		writeNormalOutputs();
		verifyDevFull();
		writeErrorOutputs();
		System.out.println("CALCTRUEQUALITY_WRITEMATRIX_CLOSE_PASS fixed="+fixed);
	}

	/** Verifies exact normal output literals for all matrix ranks. */
	private void writeNormalOutputs(){
		write1D(outdir.resolve("matrix-1d.tsv"));
		checkFile("1D", "matrix-1d.tsv", "1\t7\t2\n");
		write2D(outdir.resolve("matrix-2d.tsv"));
		checkFile("2D", "matrix-2d.tsv", "1\t1\t7\t2\n");
		write3D(outdir.resolve("matrix-3d.tsv"));
		checkFile("3D", "matrix-3d.tsv", "1\t1\t1\t7\t2\n");
		write4D(outdir.resolve("matrix-4d.tsv"));
		checkFile("4D", "matrix-4d.tsv", "1\t1\t1\t1\t7\t2\n");
		write5D(outdir.resolve("matrix-5d.tsv"));
		checkFile("5D", "matrix-5d.tsv", "1\t1\t1\t1\t1\t7\t2\n");
	}

	/** Captures return or the intended future error from all ranks. */
	private void writeErrorOutputs(){
		writeError("1D");
		writeError("2D");
		writeError("3D");
		writeError("4D");
		writeError("5D");
	}

	/** Runs one /dev/full call and applies baseline or fixed expectations. */
	private void writeError(final String rank){
		try{
			if("1D".equals(rank)){write1D(DEV_FULL);}
			else if("2D".equals(rank)){write2D(DEV_FULL);}
			else if("3D".equals(rank)){write3D(DEV_FULL);}
			else if("4D".equals(rank)){write4D(DEV_FULL);}
			else if("5D".equals(rank)){write5D(DEV_FULL);}
			else{throw new AssertionError("Unknown matrix rank "+rank);}
		}catch(Throwable t){
			if(!fixed){throw asRuntime(t);}
			check(isExpectedOutputError(t), "Unexpected fixed "+rank+
					" writeMatrix failure: "+stackSummary(t));
			System.out.println("CAUGHT_WRITEMATRIX_OUTPUT_ERROR rank="+rank);
			return;
		}
		if(fixed){
			throw fail("Fixed "+rank+" writeMatrix unexpectedly returned from "+DEV_FULL);
		}
		System.out.println("RETURNED_FROM_WRITEMATRIX rank="+rank);
	}

	/** Writes the 1D nonzero-cell matrix. */
	private static void write1D(final Path path){
		final long[] good=new long[2], bad=new long[2];
		good[1]=5;
		bad[1]=2;
		CalcTrueQuality.writeMatrix(path.toString(), good, bad, true, false, 1);
	}

	/** Writes the 2D nonzero-cell matrix. */
	private static void write2D(final Path path){
		final long[][] good=new long[2][2], bad=new long[2][2];
		good[1][1]=5;
		bad[1][1]=2;
		CalcTrueQuality.writeMatrix(path.toString(), good, bad, true, false, 1);
	}

	/** Writes the 3D nonzero-cell matrix. */
	private static void write3D(final Path path){
		final long[][][] good=new long[2][2][2], bad=new long[2][2][2];
		good[1][1][1]=5;
		bad[1][1][1]=2;
		CalcTrueQuality.writeMatrix(path.toString(), good, bad, true, false, 1);
	}

	/** Writes the 4D nonzero-cell matrix. */
	private static void write4D(final Path path){
		final long[][][][] good=new long[2][2][2][2], bad=new long[2][2][2][2];
		good[1][1][1][1]=5;
		bad[1][1][1][1]=2;
		CalcTrueQuality.writeMatrix(path.toString(), good, bad, true, false, 1);
	}

	/** Writes the 5D nonzero-cell matrix. */
	private static void write5D(final Path path){
		final long[][][][][] good=new long[2][2][2][2][2], bad=new long[2][2][2][2][2];
		good[1][1][1][1][1]=5;
		bad[1][1][1][1][1]=2;
		CalcTrueQuality.writeMatrix(path.toString(), good, bad, true, false, 1);
	}

	/** Requires the finite output-error sink used by output-close baselines. */
	private static void verifyDevFull(){
		try{
			final int mode=((Number)Files.getAttribute(DEV_FULL, "unix:mode")).intValue();
			check((mode&0170000)==0020000,
					DEV_FULL+" must be the existing character device for this fixture");
		}catch(IOException e){
			throw new RuntimeException("Could not stat "+DEV_FULL, e);
		}
	}

	/** Compares one normal matrix output to its complete expected bytes. */
	private void checkFile(final String rank, final String name, final String expected){
		try{
			final Path path=outdir.resolve(name);
			final String observed=new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
			check(expected.equals(observed), rank+" writeMatrix output differs for "+path+
					" expected="+escape(expected)+" observed="+escape(observed));
		}catch(IOException e){
			throw new RuntimeException("Could not read "+name, e);
		}
	}

	/** Creates a single-use output directory. */
	private static void createDirectory(final Path path){
		try{
			Files.createDirectory(path);
		}catch(IOException e){
			throw new RuntimeException("Could not create new output directory "+path, e);
		}
	}

	/** True only for the intended future writer-close propagation. */
	private static boolean isExpectedOutputError(final Throwable t){
		return t.getClass()==RuntimeException.class &&
				EXPECTED_OUTPUT_ERROR.equals(t.getMessage());
	}

	/** Builds a compact Throwable-chain string for matching and diagnostics. */
	private static String stackSummary(final Throwable t){
		final StringBuilder sb=new StringBuilder();
		for(Throwable x=t; x!=null; x=x.getCause()){
			if(sb.length()>0){sb.append(" caused-by ");}
			sb.append(x.getClass().getName());
			if(x.getMessage()!=null){sb.append(": ").append(x.getMessage());}
		}
		return sb.toString();
	}

	/** Escapes control characters in exact-oracle diagnostics. */
	private static String escape(final String s){
		return s.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
	}

	/** Throws when a condition is false. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw fail(message);}
	}

	/** Converts a Throwable to an unchecked failure. */
	private static RuntimeException asRuntime(final Throwable t){
		if(t instanceof RuntimeException){return (RuntimeException)t;}
		if(t instanceof Error){throw (Error)t;}
		return new RuntimeException(t);
	}

	/** Returns a RuntimeException with a diagnostic message. */
	private static RuntimeException fail(final String message){
		return new RuntimeException(message);
	}

	/** Parses and stores CalcTrueQualityWriteMatrixCloseTest arguments. */
	private static final class Args{
		Args(final String[] args0){
			final String[] args=(args0==null ? new String[0] : args0);
			final Parser parser=new Parser();
			for(final String arg : args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String key=split[0].toLowerCase(Locale.ROOT);
				final String value=(split.length>1 ? split[1] : null);
				if("outdir".equals(key)){
					outdir=value;
				}else if("fixed".equals(key)){
					fixed=Parse.parseBoolean(value);
				}else if("phase".equals(key)){
					fixed="fixed".equalsIgnoreCase(value);
					check(fixed || "baseline".equalsIgnoreCase(value),
							"phase must be baseline or fixed: "+value);
				}else if(parser.parse(arg, key, value)){
					// Consumed by BBTools common parsing.
				}else{
					throw fail("Unknown argument: "+key);
				}
			}
			check(outdir!=null && outdir.length()>0,
					"Usage: outdir=<new directory> [phase=baseline|fixed]");
		}

		/** New output directory for normal matrix files. */
		String outdir;
		/** Whether output-close errors should be propagated by the product code. */
		boolean fixed;
	}

	/** New output directory for normal files. */
	private final Path outdir;
	/** Whether to require propagated /dev/full failures. */
	private final boolean fixed;

	/** Exact future error message for an observed writeMatrix close-status failure. */
	private static final String EXPECTED_OUTPUT_ERROR="Error writing calibration matrix: /dev/full";

	/** Finite output-error sink. */
	private static final Path DEV_FULL=Paths.get("/dev/full");
}
