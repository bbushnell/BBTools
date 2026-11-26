package test.jgi;

import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.regex.Pattern;

import jgi.PhylipToFasta;
import parse.Parse;

/**
 * Observes one diagnostic-stream ownership case per fresh JVM.
 * Crosses file/stdout/stderr diagnostics with normal, constructor-error and process-error exits.
 * Linux descriptor probes and observable System stream closes avoid reflection or forced GC.
 * PreParser-internal failures, auto-exit and diagnostic write failures are outside this fixture.
 * @author Shinobu
 */
public final class PhylipDiagnosticOwnershipTest{

	/** Requires phase=baseline|fixed, a literal case and a fresh or existing-empty outdir. */
	public static void main(final String[] args) throws Exception{
		String phase=null, caseName=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			final String key=pair[0].toLowerCase(Locale.ROOT);
			if(key.equals("phase")){phase=pair[1];}
			else if(key.equals("case")){caseName=pair[1];}
			else if(key.equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown diagnostic fixture option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline or phase=fixed");
		boolean known=false;
		for(final String name : CASES){known|=name.equals(caseName);}
		require(known, "Unknown diagnostic case: "+caseName);
		require(directory!=null && !directory.isEmpty(), "Set a fresh outdir");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run dev/phylipdiagnosticownershiptest.sh with assertions enabled");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		if(!Files.exists(root)){Files.createDirectory(root);}
		require(Files.isDirectory(root), "Output path must be a directory: "+root);
		try(final DirectoryStream<Path> existing=Files.newDirectoryStream(root)){
			require(!existing.iterator().hasNext(), "Output directory must be empty: "+root);
		}
		runCase(root, phase, caseName);
	}

	/** Saves actual bytes and ownership observations before asserting the phase-specific prediction. */
	private static void runCase(final Path root, final String phase, final String caseName) throws Exception{
		final String sink=caseName.substring(0, caseName.indexOf('-'));
		final boolean file=sink.equals("file");
		final boolean constructorError=caseName.endsWith("-constructor-error");
		final boolean processError=caseName.endsWith("-process-error");
		final boolean normal=!constructorError && !processError;
		final Path input=root.resolve("input.phylip"), output=root.resolve("output.fa");
		final Path diagnostic=root.resolve("diagnostic.log");
		writeText(input, processError ? MALFORMED : NORMAL);
		Files.createFile(diagnostic);
		final int probeBefore=countDescriptors(diagnostic);
		int probeDuring;
		try(final FileInputStream probe=new FileInputStream(diagnostic.toFile())){
			assert(probe.getFD().valid()) : "The held probe must own a live descriptor for "+diagnostic;
			probeDuring=countDescriptors(diagnostic);
		}
		final int probeAfter=countDescriptors(diagnostic);
		require(probeBefore==0 && probeDuring==1 && probeAfter==0, "Expected diagnostic descriptor probe 0/1/0");
		final int fdBefore=countDescriptors(diagnostic);
		final ArrayList<String> arguments=new ArrayList<String>();
		arguments.add("outstream="+(file ? diagnostic.toString() : sink));
		if(!constructorError){arguments.add("in="+input);}
		arguments.add("out="+output);
		arguments.add("threads=1");
		arguments.add("simd=f");
		arguments.add("pigz=f");
		arguments.add("unpigz=f");
		final String[] invocation=arguments.toArray(new String[0]);
		writeText(root.resolve("arguments.txt"), String.join("\n", invocation)+"\n");

		final PrintStream oldOut=System.out, oldErr=System.err;
		final ByteArrayOutputStream stdout=new ByteArrayOutputStream(), stderr=new ByteArrayOutputStream();
		final TrackedPrintStream captureOut=new TrackedPrintStream(stdout), captureErr=new TrackedPrintStream(stderr);
		Throwable failure=null;
		try{
			System.setOut(captureOut);
			System.setErr(captureErr);
			try{PhylipToFasta.main(invocation);}
			catch(final Throwable caught){failure=caught;}
		}finally{
			captureOut.flush();
			captureErr.flush();
			System.setOut(oldOut);
			System.setErr(oldErr);
		}
		final int fdAfter=countDescriptors(diagnostic), inputFdAfter=countDescriptors(input);
		final int stdoutCloses=captureOut.closeCalls, stderrCloses=captureErr.closeCalls;
		final boolean captureError=captureOut.checkError() || captureErr.checkError();
		captureOut.close();
		captureErr.close();//Fixture-owned cleanup follows the ownership observation.
		final String diagnosticText=file ? new String(Files.readAllBytes(diagnostic), StandardCharsets.UTF_8) :
			(sink.equals("stdout") ? stdout : stderr).toString("UTF-8");
		final String exception=failure==null ? "none" : failure.getClass().getName();
		final String expectedException=constructorError ? RuntimeException.class.getName() :
			processError ? ArithmeticException.class.getName() : "none";
		final String expectedMessage=constructorError ? "Error - at least one input file is required." :
			processError ? "/ by zero" : null;
		final String top=failure==null || failure.getStackTrace().length==0 ? "none" :
			failure.getStackTrace()[0].getClassName()+"."+failure.getStackTrace()[0].getMethodName();
		final String expectedTop=constructorError ? "jgi.PhylipToFasta.<init>" :
			processError ? "jgi.PhylipToFasta.process" : "none";
		final boolean messageMatched=failure==null ? expectedMessage==null : expectedMessage!=null &&
			expectedMessage.equals(failure.getMessage());
		final boolean countsMatched=normal ?
			Pattern.compile("Reads Processed:\\s+2\\s").matcher(diagnosticText).find() &&
			Pattern.compile("Bases Processed:\\s+8\\s").matcher(diagnosticText).find() :
			!diagnosticText.contains("Reads Processed:");
		final boolean fastaMatched=normal ? Files.exists(output) &&
			Arrays.equals(FASTA.getBytes(StandardCharsets.UTF_8), Files.readAllBytes(output)) : !Files.exists(output);
		final int expectedFd=file && !normal && phase.equals("baseline") ? 1 : 0;
		final boolean streamsPreserved=stdoutCloses==0 && stderrCloses==0;
		final boolean desired=fdAfter==0 && inputFdAfter==0 && streamsPreserved;
		final boolean prediction=fdBefore==0 && fdAfter==expectedFd && inputFdAfter==0 &&
			streamsPreserved && !captureError && exception.equals(expectedException) &&
			top.equals(expectedTop) && messageMatched && countsMatched && fastaMatched &&
			diagnosticText.contains("Executing jgi.PhylipToFasta");
		Files.write(root.resolve("stdout.txt"), stdout.toByteArray(), StandardOpenOption.CREATE_NEW);
		Files.write(root.resolve("stderr.txt"), stderr.toByteArray(), StandardOpenOption.CREATE_NEW);
		final StringWriter stack=new StringWriter();
		if(failure==null){stack.write("none\n");}
		else{failure.printStackTrace(new PrintWriter(stack));}
		writeText(root.resolve("exception.txt"), stack.toString());
		final String row=caseName+"\t"+phase+"\t"+sink+"\t"+exception+"\t"+expectedException+"\t"+top+"\t"+
			messageMatched+"\t"+probeBefore+"\t"+probeDuring+"\t"+probeAfter+"\t"+fdBefore+"\t"+fdAfter+"\t"+
			expectedFd+"\t"+inputFdAfter+"\t"+stdoutCloses+"\t"+stderrCloses+"\t"+countsMatched+"\t"+
			fastaMatched+"\t"+desired+"\t"+prediction+"\n";
		writeText(root.resolve("ownership.tsv"), HEADER+"\n"+row);
		require(prediction, "Diagnostic ownership prediction failed; inspect "+root.resolve("ownership.tsv"));
		System.out.print(row);
	}

	/** Counts descriptors targeting one canonical fixture path, excluding JVM infrastructure. */
	private static int countDescriptors(final Path path) throws IOException{
		final Path directory=Paths.get("/proc/self/fd");
		require(Files.isDirectory(directory) && Files.isReadable(directory), "Linux descriptor probe unavailable");
		final Path real=path.toRealPath();
		int count=0;
		try(final DirectoryStream<Path> descriptors=Files.newDirectoryStream(directory)){
			for(final Path fd : descriptors){
				try{
					final Path target=Files.readSymbolicLink(fd);
					final Path absolute=target.isAbsolute() ? target : fd.getParent().resolve(target);
					if(real.equals(absolute.normalize())){count++;}
				}catch(final NoSuchFileException closed){
					// Descriptors can close during enumeration; other inspection errors propagate.
				}
			}
		}
		return count;
	}

	/** Creates evidence without replacing an earlier capture. */
	private static void writeText(final Path path, final String text) throws IOException{
		Files.write(path, text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Keeps fixture checks enabled independently of the program's assertion settings. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Observes attempted closure even when the backing byte-array stream ignores close. */
	private static final class TrackedPrintStream extends PrintStream{
		TrackedPrintStream(final ByteArrayOutputStream bytes) throws IOException{super(bytes, true, "UTF-8");}
		@Override
		public void close(){closeCalls++; super.close();}
		int closeCalls;
	}

	private static final String NORMAL="2 4\nalpha     AC\nbeta      TG\n\n          GT\n          CA\n";
	private static final String MALFORMED="1 2\n        AC\n";
	private static final String FASTA=">alpha\nACGT\n>beta\nTGCA\n";
	private static final String[] CASES={"file-normal", "file-constructor-error", "file-process-error",
		"stdout-normal", "stdout-constructor-error", "stdout-process-error",
		"stderr-normal", "stderr-constructor-error", "stderr-process-error"};
	private static final String HEADER="case\tphase\tsink\texception\texpectedException\ttopFrame\tmessageMatched"+
		"\tprobeBefore\tprobeDuring\tprobeAfter\tfdBefore\tfdAfter\texpectedFdAfter\tinputFdAfter"+
		"\tstdoutCloseCalls\tstderrCloseCalls\tcountsMatched\tfastaMatched\tdesiredCleanup\tpredictionMatched";
}
