package test.hiseq;

import java.io.File;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import hiseq.AnalyzeFlowCell;
import hiseq.FlowCell;
import stream.Streamer;
import stream.StreamerFactory;

/** Finite real-backend checks for optional SAM failure propagation, separate from report fitting.
 * @author Brian Bushnell, Sucrose */
public class AnalyzeFlowCellSamFailureTest{

	/** Runs phase=baseline or phase=fixed in a dedicated JVM with a fresh outdir. */
	public static void main(String[] args) throws Exception{
		File out=null;
		boolean fixed=false;
		for(String arg : args){
			if(arg.startsWith("outdir=")){out=new File(arg.substring(7));}
			else if(arg.equals("phase=fixed")){fixed=true;}
			else if(!arg.equals("phase=baseline")){throw new IllegalArgumentException(arg);}
		}
		require(out!=null && out.mkdir(), "Use a fresh fixture directory");
		final File primary=new File(out, "empty-primary.fq");
		Files.write(primary.toPath(), new byte[0]);
		final String prefix="INST:1:FC:1:1101:100:100\t";
		final String suffix="\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";
		final File good=new File(out, "valid.sam"), bad=new File(out, "invalid-flag.sam");
		Files.write(good.toPath(), (prefix+"4"+suffix).getBytes(StandardCharsets.US_ASCII));
		Files.write(bad.toPath(), (prefix+"x"+suffix).getBytes(StandardCharsets.US_ASCII));
		System.out.println("case\tobserved\texpected\tresult");
		//Fresh tool instances isolate the caller's sticky error flag and aggregate counters.
		runCaller(primary, good, "valid-sam", false);
		runCaller(primary, bad, "invalid-sam", fixed);
		probeBackend(bad);
	}

	/** Exercises the actual optional-reader creation/consumption/cleanup path with no primary observations. */
	private static void runCaller(File primary, File sam, String name, boolean expectedError) throws Exception{
		final AnalyzeFlowCell tool=new AnalyzeFlowCell(new String[] {"in="+primary, "sam="+sam,
			"threads=1", "loadkmers=f", "align=f", "simd=f", "overwrite=t", "qin=33"});
		final FlowCell flowcell=new FlowCell(31);
		final Field field=AnalyzeFlowCell.class.getDeclaredField("flowcell");
		field.setAccessible(true); field.set(tool, flowcell);
		final Method fill=AnalyzeFlowCell.class.getDeclaredMethod("fillTiles");
		fill.setAccessible(true);
		try{fill.invoke(tool);}
		catch(InvocationTargetException e){throw new AssertionError("fillTiles threw instead of returning its error flag: "+name, e.getCause());}
		check(name+"-error", tool.errorState ? 1 : 0, expectedError ? 1 : 0);
		check(name+"-primary-reads", tool.readsProcessed, 0);
		check(name+"-primary-bases", tool.basesProcessed, 0);
		check(name+"-mapped-tiles", flowcell.toList().size(), 0);
	}

	/** Confirms the real one-thread backend reports the malformed flag, then explicitly closes that probe. */
	private static void probeBackend(File bad){
		final FileFormat ff=FileFormat.testInput(bad.getPath(), FileFormat.SAM, null, true, false);
		final Streamer stream=StreamerFactory.makeSamOrBamStreamer(ff, 1, false, false, -1, false);
		boolean closeError=false;
		try{
			require("stream.SamStreamerST".equals(stream.getClass().getName()),
				"threads=1 must select the investigated backend: "+stream.getClass().getName());
			System.out.println("OBSERVED\tbackend="+stream.getClass().getName());
			stream.start();
			check("backend-first-batch-is-terminal", stream.nextLines()==null ? 1 : 0, 1);
			check("backend-error-before-close", stream.errorState() ? 1 : 0, 1);
		}finally{closeError=ReadWrite.closeStream(stream);}
		check("backend-close-error", closeError ? 1 : 0, 1);
	}

	/** Reports literal counters/flags rather than relying on the fixture's own JVM assertion setting. */
	private static void check(String name, long observed, long expected){
		System.out.println(name+"\t"+observed+"\t"+expected+"\t"+(observed==expected ? "PASS" : "FAIL"));
		require(observed==expected, name+": observed="+observed+", expected="+expected);
	}

	/** Fails loudly when an independent fixture expectation is not met. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
