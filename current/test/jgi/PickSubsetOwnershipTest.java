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

import fileIO.ByteStreamWriter;
import jgi.PickSubset;
import parse.Parse;

/**
 * Observes one PickSubset exceptional-ownership case per fresh JVM.
 * Crosses file/stdout/stderr diagnostics with constructor and post-writer failures.
 * Linux descriptor probes, live writer-thread counts and observable System stream closes avoid reflection or forced GC.
 * @author Jean
 */
public final class PickSubsetOwnershipTest{

	/** Parses one literal case, then exits explicitly so baseline orphan writers cannot keep the JVM alive. */
	public static void main(final String[] args){
		int exitCode=1;
		Path root=null;
		try{
			String phase=null, caseName=null, directory=null;
			for(final String arg : args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				require(pair.length==2, "Expected flag=value: "+arg);
				final String key=pair[0].toLowerCase(Locale.ROOT);
				if(key.equals("phase")){phase=pair[1];}
				else if(key.equals("case")){caseName=pair[1];}
				else if(key.equals("outdir")){directory=pair[1];}
				else{throw new IllegalArgumentException("Unknown PickSubset ownership fixture option: "+arg);}
			}
			require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline or phase=fixed");
			boolean known=false;
			for(final String name : CASES){known|=name.equals(caseName);}
			require(known, "Unknown ownership case: "+caseName);
			require(directory!=null && !directory.isEmpty(), "Set a fresh outdir");
			boolean assertions=false;
			assert(assertions=true);
			require(assertions, "Run dev/picksubsetownershiptest.sh with assertions enabled");
			root=Paths.get(directory).toAbsolutePath().normalize();
			if(!Files.exists(root)){Files.createDirectory(root);}
			require(Files.isDirectory(root), "Output path must be a directory: "+root);
			try(final DirectoryStream<Path> existing=Files.newDirectoryStream(root)){
				require(!existing.iterator().hasNext(), "Output directory must be empty: "+root);
			}
			runCase(root, phase, caseName);
			exitCode=0;
		}catch(final Throwable failure){
			failure.printStackTrace();
			if(root!=null){
				try{
					final StringWriter stack=new StringWriter();
					failure.printStackTrace(new PrintWriter(stack));
					writeText(root.resolve("fixture-failure.txt"), stack.toString());
				}catch(final Throwable ignored){
					// The stack printed above is the fallback report when evidence writing itself failed.
				}
			}
		}finally{
			System.exit(exitCode);
		}
	}

	/** Saves input, output, exception and ownership observations before checking phase-specific predictions. */
	private static void runCase(final Path root, final String phase, final String caseName) throws Exception{
		final boolean baseline=phase.equals("baseline");
		final boolean fileDiagnostic=caseName.startsWith("file-");
		final boolean constructorError=caseName.endsWith("constructor-error");
		final boolean eofError=caseName.equals("file-eof-error");
		final boolean preEofError=caseName.endsWith("preeof-error");
		final boolean normal=caseName.endsWith("normal");
		final String sink=fileDiagnostic ? "file" : caseName.substring(0, caseName.indexOf('-'));
		final Path input=root.resolve("input.tsv"), kept=root.resolve("kept.txt"),
			discarded=root.resolve("discarded.txt"), diagnostic=root.resolve("diagnostic.log");
		writeText(kept, "");
		writeText(discarded, "");
		if(!constructorError){writeText(input, eofError ? EOF_INCOMPLETE : COMPLETE);}
		if(fileDiagnostic){writeText(diagnostic, "");}
		final Path probe=fileDiagnostic ? diagnostic : kept;
		final int probeBefore=countDescriptors(probe);
		int probeDuring;
		try(final FileInputStream held=new FileInputStream(probe.toFile())){
			assert(held.getFD().valid()) : "The descriptor probe must hold a live handle for "+probe;
			probeDuring=countDescriptors(probe);
		}
		final int probeAfter=countDescriptors(probe);
		require(probeBefore==0 && probeDuring==1 && probeAfter==0, "Descriptor probe did not observe 0/1/0 for "+probe);
		final int diagnosticFdBefore=fileDiagnostic ? countDescriptors(diagnostic) : 0;
		final int inputFdBefore=countDescriptors(input), keptFdBefore=countDescriptors(kept),
			discardedFdBefore=countDescriptors(discarded), writersBefore=liveWriters();

		final ArrayList<String> arguments=new ArrayList<String>();
		arguments.add("outstream="+(fileDiagnostic ? diagnostic.toString() : sink));
		arguments.add("out="+kept);
		arguments.add("invalid="+discarded);
		if(!constructorError){arguments.add("in="+input);}
		arguments.add(normal || constructorError ? "files=1" : "files=0");
		arguments.add("maxani=0.1");
		if(preEofError){arguments.add("lines=2");}
		arguments.add("bf1=t");
		arguments.add("ow=t");
		arguments.add("threads=1");
		arguments.add("simd=f");
		arguments.add("pigz=f");
		arguments.add("unpigz=f");
		final String[] invocation=arguments.toArray(new String[0]);
		writeText(root.resolve("arguments.txt"), String.join("\n", invocation)+"\n");
		if(constructorError){writeText(root.resolve("input.bytes"), "");}
		else{Files.write(root.resolve("input.bytes"), Files.readAllBytes(input), StandardOpenOption.CREATE_NEW);}

		final PrintStream oldOut=System.out, oldErr=System.err;
		final ByteArrayOutputStream stdout=new ByteArrayOutputStream(), stderr=new ByteArrayOutputStream();
		final TrackedPrintStream captureOut=new TrackedPrintStream(stdout), captureErr=new TrackedPrintStream(stderr);
		Throwable failure=null;
		try{
			System.setOut(captureOut);
			System.setErr(captureErr);
			try{PickSubset.main(invocation);}
			catch(final Throwable caught){failure=caught;}
		}finally{
			captureOut.flush();
			captureErr.flush();
			System.setOut(oldOut);
			System.setErr(oldErr);
		}

		final int diagnosticFdAfter=fileDiagnostic ? countDescriptors(diagnostic) : 0;
		final int inputFdAfter=countDescriptors(input), keptFdAfter=countDescriptors(kept),
			discardedFdAfter=countDescriptors(discarded), writersAfter=liveWriters();
		final int stdoutCloses=captureOut.closeCalls, stderrCloses=captureErr.closeCalls;
		final boolean captureError=captureOut.checkError() || captureErr.checkError();
		captureOut.close();
		captureErr.close();
		final String diagnosticText=fileDiagnostic ? readIfPresent(diagnostic) :
			(sink.equals("stdout") ? stdout : stderr).toString("UTF-8");
		final byte[] keptBytes=readBytesIfPresent(kept), discardedBytes=readBytesIfPresent(discarded);
		Files.write(root.resolve("stdout.txt"), stdout.toByteArray(), StandardOpenOption.CREATE_NEW);
		Files.write(root.resolve("stderr.txt"), stderr.toByteArray(), StandardOpenOption.CREATE_NEW);
		Files.write(root.resolve("kept.bytes"), keptBytes, StandardOpenOption.CREATE_NEW);
		Files.write(root.resolve("discarded.bytes"), discardedBytes, StandardOpenOption.CREATE_NEW);
		writeText(root.resolve("diagnostic.txt"), diagnosticText);
		final StringWriter stack=new StringWriter();
		if(failure==null){stack.write("none\n");}
		else{failure.printStackTrace(new PrintWriter(stack));}
		writeText(root.resolve("exception.txt"), stack.toString());

		final String exception=failure==null ? "none" : failure.getClass().getName();
		final String top=failure==null || failure.getStackTrace().length==0 ? "none" :
			failure.getStackTrace()[0].getClassName()+"."+failure.getStackTrace()[0].getMethodName();
		final String expectedException=constructorError ? RuntimeException.class.getName() :
			(normal ? "none" : (baseline ? NullPointerException.class.getName() :
			IllegalArgumentException.class.getName()));
		final String expectedTop=constructorError ? "jgi.PickSubset.fixExtensions" :
			(normal ? "none" : (baseline ? "jgi.PickSubset.pickNodes" :
			"jgi.PickSubset.referenceNode"));
		final int expectedDiagnosticFd=baseline && fileDiagnostic && !normal ? 1 : 0;
		final int expectedInputFd=baseline && preEofError ? 1 : 0;
		final int expectedKeptFd=baseline && !constructorError && !normal ? 1 : 0;
		final int expectedDiscardedFd=expectedKeptFd;
		final int expectedWriters=baseline && !constructorError && !normal ? 2 : 0;
		final boolean normalOutputs=normal && Arrays.equals(ALPHA.getBytes(StandardCharsets.UTF_8), keptBytes) &&
			Arrays.equals(BETA.getBytes(StandardCharsets.UTF_8), discardedBytes);
		final boolean failureOutputs=!normal && keptBytes.length==0 && discardedBytes.length==0;
		final boolean summaryMatched=normal ?
			Pattern.compile("Nodes In:\\s+2\\s").matcher(diagnosticText).find() &&
			Pattern.compile("Nodes Out:\\s+1\\s").matcher(diagnosticText).find() :
			!diagnosticText.contains("Nodes In:") && !diagnosticText.contains("Processed:");
		final boolean borrowedPreserved=stdoutCloses==0 && stderrCloses==0;
		final boolean descriptorsMatched=diagnosticFdBefore==0 && inputFdBefore==0 && keptFdBefore==0 &&
			discardedFdBefore==0 && diagnosticFdAfter==expectedDiagnosticFd && inputFdAfter==expectedInputFd &&
			keptFdAfter==expectedKeptFd && discardedFdAfter==expectedDiscardedFd;
		final boolean prediction=writersBefore==0 && writersAfter==expectedWriters && descriptorsMatched &&
			borrowedPreserved && !captureError && exception.equals(expectedException) && top.equals(expectedTop) &&
			summaryMatched && (normal ? normalOutputs : failureOutputs);
		final String row=caseName+"\t"+phase+"\t"+sink+"\t"+exception+"\t"+expectedException+"\t"+top+"\t"+
			probeBefore+"\t"+probeDuring+"\t"+probeAfter+"\t"+diagnosticFdBefore+"\t"+diagnosticFdAfter+"\t"+
			expectedDiagnosticFd+"\t"+inputFdBefore+"\t"+inputFdAfter+"\t"+expectedInputFd+"\t"+
			keptFdBefore+"\t"+keptFdAfter+"\t"+expectedKeptFd+"\t"+discardedFdBefore+"\t"+discardedFdAfter+"\t"+
			expectedDiscardedFd+"\t"+writersBefore+"\t"+writersAfter+"\t"+expectedWriters+"\t"+
			stdoutCloses+"\t"+stderrCloses+"\t"+summaryMatched+"\t"+normalOutputs+"\t"+prediction+"\n";
		writeText(root.resolve("ownership.tsv"), HEADER+"\n"+row);
		require(prediction, "PickSubset ownership prediction failed; inspect "+root.resolve("ownership.tsv"));
		System.out.print(row);
	}

	/** Counts descriptors targeting one canonical fixture path, excluding JVM infrastructure. */
	private static int countDescriptors(final Path path) throws IOException{
		if(!Files.exists(path)){return 0;}
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

	/** Counts live ByteStreamWriter workers visible to this fresh fixture JVM. */
	private static int liveWriters(){
		int count=0;
		for(final Thread thread : Thread.getAllStackTraces().keySet()){
			if(thread instanceof ByteStreamWriter && thread.isAlive()){count++;}
		}
		return count;
	}

	/** Creates evidence without replacing an earlier capture. */
	private static void writeText(final Path path, final String text) throws IOException{
		Files.write(path, text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Reads an output that may not be created when constructor validation fails. */
	private static byte[] readBytesIfPresent(final Path path) throws IOException{
		return Files.exists(path) ? Files.readAllBytes(path) : new byte[0];
	}

	/** Decodes a diagnostic that may be borrowed stdout/stderr rather than a file. */
	private static String readIfPresent(final Path path) throws IOException{
		return Files.exists(path) ? new String(Files.readAllBytes(path), StandardCharsets.UTF_8) : "";
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

	private static final String HEADER="case\tphase\tsink\texception\texpectedException\ttopFrame"+
		"\tprobeBefore\tprobeDuring\tprobeAfter\tdiagnosticFdBefore\tdiagnosticFdAfter\texpectedDiagnosticFd"+
		"\tinputFdBefore\tinputFdAfter\texpectedInputFd\tkeptFdBefore\tkeptFdAfter\texpectedKeptFd"+
		"\tdiscardedFdBefore\tdiscardedFdAfter\texpectedDiscardedFd\twritersBefore\twritersAfter"+
		"\texpectedWriters\tstdoutCloseCalls\tstderrCloseCalls\tsummaryMatched\tnormalOutputs\tpredictionMatched";
	private static final String ALPHA="alpha\n";
	private static final String BETA="beta\n";
	private static final String COMPLETE="#Query\tRef\tANI\nalpha\tbeta\t0.9\nbeta\talpha\t0.9\n"+
		"alpha\talpha\t1.0\nbeta\tbeta\t1.0\n";
	private static final String EOF_INCOMPLETE="#Query\tRef\tANI\nalpha\tbeta\t0.9\n";
	private static final String[] CASES={"file-normal", "file-constructor-error", "file-eof-error",
		"file-preeof-error", "stdout-preeof-error", "stderr-preeof-error"};
}
