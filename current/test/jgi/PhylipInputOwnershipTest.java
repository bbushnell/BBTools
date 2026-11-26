package test.jgi;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
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
import java.util.Locale;
import java.util.regex.Pattern;
import java.util.zip.GZIPOutputStream;

import jgi.PhylipToFasta;
import parse.Parse;

/**
 * Captures one direct-main input-ownership case per fresh JVM.
 * Linux descriptor observations are calibrated with a separately held stream.
 * Borrowed raw/gzip stdin uses an observable close counter, never reflection.
 * File-option and format policy are outside this fixture's scope.
 * @author Shinobu
 */
public final class PhylipInputOwnershipTest{

	/** Requires phase=baseline|fixed, a literal case, and a fresh or existing-empty output directory. */
	public static void main(final String[] args) throws Exception{
		String phase=null, caseName=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			final String key=pair[0].toLowerCase(Locale.ROOT);
			if(key.equals("phase")){phase=pair[1];}
			else if(key.equals("case")){caseName=pair[1];}
			else if(key.equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown ownership fixture option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline or phase=fixed");
		boolean known=false;
		for(final String name : CASES){known|=name.equals(caseName);}
		require(known, "Unknown ownership case: "+caseName);
		require(directory!=null && !directory.isEmpty(), "Set a fresh outdir");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run dev/phylipinputownershiptest.sh with assertions enabled");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		if(!Files.exists(root)){Files.createDirectory(root);}
		require(Files.isDirectory(root), "Output path must be a directory: "+root);
		try(final DirectoryStream<Path> existing=Files.newDirectoryStream(root)){
			require(!existing.iterator().hasNext(), "Output directory must be empty: "+root);
		}
		runCase(root, phase, caseName);
	}

	/** Records bytes, descriptors, exceptions and borrowed-input closure before checking predictions. */
	private static void runCase(final Path root, final String phase, final String caseName) throws Exception{
		final boolean borrowed=caseName.startsWith("stdin-");
		final boolean gzip=caseName.contains("-gzip-");
		final boolean parseFailure=caseName.endsWith("-error");
		final byte[] decoded=(parseFailure ? MALFORMED : NORMAL).getBytes(StandardCharsets.UTF_8);
		final byte[] encoded=gzip ? gzip(decoded) : decoded;
		final String leaf=borrowed ? "stdin-seed.bin" : caseName.equals("file-stdin-name") ?
			"stdin.phylip" : gzip ? "input.phylip.gz" : "input.phylip";
		final Path input=root.resolve(leaf);
		final boolean stdinNamedFile=caseName.equals("file-stdin-name");
		if(stdinNamedFile){
			require(root.equals(Paths.get("").toAbsolutePath().normalize()),
				"Run file-stdin-name with its fresh output directory as the working directory");
		}
		final String inputArg=borrowed ? (gzip ? "stdin.phylip.gz" : "stdin.phylip") :
			stdinNamedFile ? "stdin.phylip" : input.toString();
		if(borrowed){require(!Files.exists(Paths.get(inputArg)), "A real file would shadow borrowed stdin: "+inputArg);}
		Files.write(input, encoded, StandardOpenOption.CREATE_NEW);

		final int probeBefore=countDescriptors(input);
		int probeDuring;
		try(final FileInputStream probe=new FileInputStream(input.toFile())){
			assert(probe.getFD().valid()) : "The held probe must own a live descriptor for "+input;
			probeDuring=countDescriptors(input);
		}
		final int probeAfter=countDescriptors(input);
		final boolean probeMatched=probeBefore==0 && probeDuring==1 && probeAfter==0;
		require(probeMatched, "Descriptor probe did not observe 0/1/0 for "+input);
		final int fdBefore=countDescriptors(input);

		final InputStream oldIn=System.in;
		final PrintStream oldOut=System.out, oldErr=System.err;
		final ByteArrayOutputStream stdout=new ByteArrayOutputStream(), stderr=new ByteArrayOutputStream();
		final PrintStream captureOut=new PrintStream(stdout, true, "UTF-8");
		final PrintStream captureErr=new PrintStream(stderr, true, "UTF-8");
		final TrackedInput stdin=new TrackedInput(borrowed ? encoded : "Unused stdin\n".getBytes(StandardCharsets.UTF_8));
		final String[] invocation={"in="+inputArg, "out=null", "threads=1", "simd=f", "pigz=f", "unpigz=f"};
		writeText(root.resolve("arguments.txt"), String.join("\n", invocation)+"\n");
		Throwable failure=null;
		try{
			System.setOut(captureOut);
			System.setErr(captureErr);
			System.setIn(stdin);
			try{PhylipToFasta.main(invocation);}
			catch(final Throwable caught){failure=caught;}
		}finally{
			captureOut.flush();
			captureErr.flush();
			System.setIn(oldIn);
			System.setOut(oldOut);
			System.setErr(oldErr);
		}
		final int fdAfter=countDescriptors(input);
		final int closeCalls=stdin.closeCalls, readCalls=stdin.readCalls;
		final boolean outputError=captureOut.checkError() || captureErr.checkError();
		captureOut.close();
		captureErr.close();
		stdin.close();//Fixture-owned cleanup occurs after the observation.

		final String stderrText=stderr.toString("UTF-8");
		final String exception=failure==null ? "none" : failure.getClass().getName();
		final String expectedException=parseFailure ? ArithmeticException.class.getName() : "none";
		final String top=failure==null || failure.getStackTrace().length==0 ? "none" :
			failure.getStackTrace()[0].getClassName()+"."+failure.getStackTrace()[0].getMethodName();
		final int expectedFd="baseline".equals(phase) && !borrowed ? 1 : 0;
		final boolean countsMatched=parseFailure ? !stderrText.contains("Reads Processed:") :
			Pattern.compile("Reads Processed:\\s+2\\s").matcher(stderrText).find() &&
			Pattern.compile("Bases Processed:\\s+8\\s").matcher(stderrText).find();
		final boolean exceptionMatched=exception.equals(expectedException) &&
			(!parseFailure || top.equals("jgi.PhylipToFasta.process"));
		final boolean stdinPreserved=closeCalls==0 && (borrowed ? readCalls>0 : readCalls==0);
		final boolean desired=fdAfter==0 && stdinPreserved;
		final boolean prediction=fdBefore==0 && fdAfter==expectedFd && exceptionMatched &&
			stdinPreserved && stdout.size()==0 && !outputError && countsMatched;
		Files.write(root.resolve("stdout.txt"), stdout.toByteArray(), StandardOpenOption.CREATE_NEW);
		Files.write(root.resolve("stderr.txt"), stderr.toByteArray(), StandardOpenOption.CREATE_NEW);
		final StringWriter stack=new StringWriter();
		if(failure==null){stack.write("none\n");}
		else{failure.printStackTrace(new PrintWriter(stack));}
		writeText(root.resolve("exception.txt"), stack.toString());
		final String row=caseName+"\t"+phase+"\t"+borrowed+"\t"+gzip+"\t"+exception+"\t"+
			expectedException+"\t"+top+"\t"+probeBefore+"\t"+probeDuring+"\t"+probeAfter+"\t"+
			fdBefore+"\t"+fdAfter+"\t"+expectedFd+"\t"+closeCalls+"\t"+readCalls+"\t"+stdout.size()+"\t"+
			countsMatched+"\t"+desired+"\t"+prediction+"\n";
		writeText(root.resolve("ownership.tsv"), HEADER+"\n"+row);
		require(prediction, "Ownership prediction failed; inspect "+root.resolve("ownership.tsv"));
		System.out.print(row);
	}

	/** Creates a standard Java gzip member without external codecs or availability gates. */
	private static byte[] gzip(final byte[] data) throws IOException{
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		try(final GZIPOutputStream gzip=new GZIPOutputStream(bytes)){gzip.write(data);}
		return bytes.toByteArray();
	}

	/** Counts only descriptors whose symlink targets identify this input, not JVM infrastructure. */
	private static int countDescriptors(final Path input) throws IOException{
		final Path directory=Paths.get("/proc/self/fd");
		require(Files.isDirectory(directory) && Files.isReadable(directory), "Linux descriptor probe unavailable");
		final Path real=input.toRealPath();
		int count=0;
		try(final DirectoryStream<Path> descriptors=Files.newDirectoryStream(directory)){
			for(final Path fd : descriptors){
				try{
					final Path target=Files.readSymbolicLink(fd);
					final Path absolute=target.isAbsolute() ? target : fd.getParent().resolve(target);
					if(real.equals(absolute.normalize())){count++;}
				}catch(final NoSuchFileException closed){
					// A descriptor may close during enumeration; other inspection errors propagate.
				}
			}
		}
		return count;
	}

	/** Writes fresh evidence; it must not silently replace an earlier capture. */
	private static void writeText(final Path path, final String text) throws IOException{
		Files.write(path, text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Keeps fixture requirements active independently of the converter's assertions. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Makes accidental closure of borrowed stdin observable even though byte-array close is harmless. */
	private static final class TrackedInput extends ByteArrayInputStream{
		TrackedInput(final byte[] data){super(data);}
		@Override
		public int read(){readCalls++; return super.read();}
		@Override
		public int read(final byte[] data, final int off, final int len){readCalls++; return super.read(data, off, len);}
		@Override
		public void close() throws IOException{closeCalls++; super.close();}
		int closeCalls, readCalls;
	}

	private static final String NORMAL="2 4\n \t\nalpha     AC\nbeta      TG\n\n\t \t\n          GT\n          CA\n";
	private static final String MALFORMED="1 2\n        AC\n";
	private static final String[] CASES={"file-normal", "file-error", "file-gzip-normal", "file-gzip-error",
		"file-stdin-name", "stdin-normal", "stdin-error", "stdin-gzip-normal", "stdin-gzip-error"};
	private static final String HEADER="case\tphase\tborrowed\tgzip\texception\texpectedException\ttopFrame"+
		"\tprobeBefore\tprobeDuring\tprobeAfter\tfdBefore\tfdAfter\texpectedFdAfter\tstdinCloseCalls"+
		"\tstdinReadCalls\tstdoutBytes\tcountsMatched\tdesiredCleanup\tpredictionMatched";
}
