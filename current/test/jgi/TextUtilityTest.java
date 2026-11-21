package test.jgi;

import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.OutputStream;
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
import java.util.HashSet;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import jgi.Difference;
import jgi.FindString;
import parse.Parse;

/**
 * Direct-main baseline and fixed-behavior capture for Difference and FindString.
 * Each JVM handles one literal case, writes tiny inputs through ByteStreamWriter,
 * counts open descriptors to those inputs, captures stdout/stderr, records the
 * complete thrown Throwable, and preserves current physical-line versus
 * logical-line behavior.
 * @author Jean, Shinobu
 */
public final class TextUtilityTest{

	/** Requires a phase, a named case and a fresh output directory. */
	public static void main(final String[] args) throws Exception{
		String phase=null, caseName=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			final String key=pair[0].toLowerCase(Locale.ROOT);
			if(key.equals("phase")){phase=pair[1];}
			else if(key.equals("case")){caseName=pair[1];}
			else if(key.equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown text utility fixture option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase),
			"Set phase=baseline or phase=fixed.");
		require(caseName!=null && !caseName.isEmpty(), "Set case=<literal>.");
		final CaseLiteral literal=CaseLiteral.byName(caseName);
		require(directory!=null && !directory.isEmpty(), "Set a fresh outdir.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/textutilitytest.sh with assertions enabled.");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		Files.createDirectory(root);
		new TextUtilityTest(root, phase).row(literal);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates the stable one-row report for this process. */
	private TextUtilityTest(final Path root_, final String phase_) throws IOException{
		root=root_;
		phase=phase_;
		Files.write(root.resolve("text-utility-cases.tsv"),
			(HEADER+"\n").getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes inputs, invokes one selected main directly and records the outcome. */
	private void row(final CaseLiteral literal) throws Exception{
		final ArrayList<Path> inputs=new ArrayList<Path>();
		final Path probe=root.resolve("descriptor-probe.txt");
		writeInput(probe, "probe\n");
		final DescriptorProbe probeResult=verifyDescriptorProbe(probe);
		final String[] mainArgs=literal.writeInputs(root, inputs);
		final int fdBefore=countOpenDescriptors(inputs);
		final MainResult result=invoke(literal, mainArgs);
		final int fdAfter=countOpenDescriptors(inputs);
		final int fdDelta=fdAfter-fdBefore;

		Files.write(root.resolve(literal.name+".stdout"),
			result.stdout.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		Files.write(root.resolve(literal.name+".stderr"),
			result.stderr.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		Files.write(root.resolve(literal.name+".exception.txt"),
			throwableText(result.failure).getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);

		final Expected expected="baseline".equals(phase) ? literal.old : literal.fixed;
		final boolean stdoutMatched=expected.stdout.equals(result.stdout);
		final boolean stderrMatched=expected.stderr.equals(result.stderr);
		final boolean exceptionMatched=expected.exception.equals(type(result.failure));
		final boolean messageMatched=messageMatches(expected, result.failure);
		final boolean stdoutErrorMatched=expected.stdoutError==result.stdoutError;
		final boolean stderrErrorMatched=expected.stderrError==result.stderrError;
		final boolean fdBeforeMatched=fdBefore==0;
		final boolean fdMatched=expected.fdDelta==fdDelta;
		final boolean probeMatched=probeResult.matches();
		final boolean predictionMatched=stdoutMatched && stderrMatched &&
			exceptionMatched && messageMatched && stdoutErrorMatched &&
			stderrErrorMatched && fdBeforeMatched && fdMatched && probeMatched;
		final boolean desiredCleanup=fdBefore==0 && fdAfter==0;
		final boolean desiredDiagnostic=!literal.needsDiagnostic ||
			(literal.fixed.exception.equals(type(result.failure)) &&
				messageMatches(literal.fixed, result.failure));
		final boolean outputErrorUnreported=literal.suspectedOutput &&
			result.failure==null && (result.stdoutError || result.stderrError);
		final boolean legacyFailure=!desiredCleanup || !desiredDiagnostic ||
			outputErrorUnreported;

		final String row=literal.name+'\t'+phase+'\t'+literal.tool+'\t'+
			literal.description+'\t'+literal.physicalLines+'\t'+literal.logicalLines+
			'\t'+probeResult.before+'\t'+probeResult.during+'\t'+probeResult.after+
			'\t'+probeMatched+'\t'+fdBefore+'\t'+fdAfter+'\t'+fdDelta+
			'\t'+expected.fdDelta+'\t'+type(result.failure)+'\t'+
			clean(result.failure==null ? "" : String.valueOf(result.failure.getMessage()))+
			'\t'+topFrame(result.failure)+'\t'+result.stdoutError+
			'\t'+result.stderrError+'\t'+result.stdout.length()+
			'\t'+result.stderr.length()+'\t'+stdoutMatched+'\t'+stderrMatched+
			'\t'+exceptionMatched+'\t'+messageMatched+'\t'+stdoutErrorMatched+
			'\t'+stderrErrorMatched+'\t'+fdBeforeMatched+'\t'+fdMatched+
			'\t'+desiredCleanup+'\t'+desiredDiagnostic+'\t'+outputErrorUnreported+
			'\t'+legacyFailure+'\t'+predictionMatched+'\n';
		Files.write(root.resolve("text-utility-cases.tsv"),
			row.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);

		final String summary="TEXT_UTILITY_CAPTURE_"+
			(predictionMatched ? "PASS" : "FAIL")+" phase="+phase+
			" case="+literal.name+" utility="+literal.tool+
			" predictionMatched="+predictionMatched+
			" fdBefore="+fdBefore+" fdAfter="+fdAfter+
			" directException="+type(result.failure)+
			" descriptorProbeMatched="+probeMatched+
			" desiredCleanupMatched="+desiredCleanup+
			" desiredDiagnosticMatched="+desiredDiagnostic+
			" outputErrorUnreported="+outputErrorUnreported+
			" legacyFailureObserved="+legacyFailure+"\n";
		Files.write(root.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		System.out.print(summary);
		if(!predictionMatched){throw new AssertionError(summary+
			"Inspect saved text utility row in "+root);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Tests the supplied expectation, independently of the selected reporting phase. */
	private static boolean messageMatches(final Expected expected, final Throwable failure){
		return expected.messageToken.length()<1 ||
			(failure!=null && String.valueOf(failure.getMessage()).contains(expected.messageToken));
	}

	/** Directly invokes one selected utility main with stdout/stderr restored. */
	private static MainResult invoke(final CaseLiteral literal, final String[] mainArgs)
			throws Exception{
		final PrintStream oldOut=System.out, oldErr=System.err;
		final ByteArrayOutputStream stdout=new ByteArrayOutputStream();
		final ByteArrayOutputStream stderr=new ByteArrayOutputStream();
		final PrintStream newOut=new PrintStream(literal.failStdout ?
			new ThrowingOutputStream("stdout") : stdout, true, "UTF-8");
		final PrintStream newErr=new PrintStream(literal.failStderr ?
			new ThrowingOutputStream("stderr") : stderr, true, "UTF-8");
		Throwable failure=null;
		boolean stdoutError=false, stderrError=false;
		try{
			System.setOut(newOut);
			System.setErr(newErr);
			if(literal.tool.equals("Difference")){Difference.main(mainArgs);}
			else if(literal.tool.equals("FindString")){FindString.main(mainArgs);}
			else{throw new AssertionError("Unknown utility "+literal.tool);}
		}catch(final Throwable caught){
			failure=caught;
		}finally{
			newOut.flush();
			newErr.flush();
			stdoutError=newOut.checkError();
			stderrError=newErr.checkError();
			System.setOut(oldOut);
			System.setErr(oldErr);
			newOut.close();
			newErr.close();
		}
		return new MainResult(
			new String(stdout.toByteArray(), StandardCharsets.UTF_8),
			new String(stderr.toByteArray(), StandardCharsets.UTF_8),
			stdoutError, stderrError, failure);
	}

	/** Writes one tiny literal file with a BBTools byte writer. */
	private static void writeInput(final Path path, final String text) throws IOException{
		final ByteStreamWriter bsw=new ByteStreamWriter(path.toString(), true, false, false);
		bsw.start();
		bsw.print(text);
		if(bsw.poisonAndWait()){throw new IOException("Error writing "+path);}
	}

	/** Verifies that the fd probe observes one independently held descriptor. */
	private static DescriptorProbe verifyDescriptorProbe(final Path path) throws IOException{
		final ArrayList<Path> one=new ArrayList<Path>();
		one.add(path);
		final int before=countOpenDescriptors(one);
		int during;
		try(final FileInputStream fis=new FileInputStream(path.toFile())){
			assert(fis.getFD().valid()) : "probe stream should be live for "+path;
			during=countOpenDescriptors(one);
		}
		final int after=countOpenDescriptors(one);
		return new DescriptorProbe(before, during, after);
	}

	/** Counts live descriptors pointing at this case's input files. */
	private static int countOpenDescriptors(final ArrayList<Path> inputs) throws IOException{
		final Path fdRoot=Paths.get("/proc/self/fd");
		if(!Files.isDirectory(fdRoot) || !Files.isReadable(fdRoot)){
			throw new IOException("Required descriptor probe is unavailable: "+fdRoot);
		}
		final HashSet<String> realInputs=new HashSet<String>();
		for(final Path input : inputs){realInputs.add(input.toRealPath().toString());}
		int count=0;
		try(final DirectoryStream<Path> fds=Files.newDirectoryStream(fdRoot)){
			for(final Path fd : fds){
				try{
					final Path target=Files.readSymbolicLink(fd);
					final Path absolute=target.isAbsolute() ? target :
						fd.getParent().resolve(target);
					if(realInputs.contains(absolute.normalize().toString())){count++;}
				}catch(final NoSuchFileException nsfe){
					// Ignore descriptors that closed while /proc/self/fd was enumerated.
				}catch(final IOException ioe){
					throw new IOException("Could not inspect "+fd, ioe);
				}
			}
		}
		return count;
	}

	/** Writes the complete Throwable chain instead of a flattened TSV message. */
	private static String throwableText(final Throwable error){
		if(error==null){return "none\n";}
		final StringWriter sw=new StringWriter();
		error.printStackTrace(new PrintWriter(sw));
		return sw.toString();
	}

	/** Makes one value fit in a TSV cell. */
	private static String clean(final String value){
		return value.replace('\t', ' ').replace('\r', ' ').replace('\n', ' ');
	}

	/** Records absent exceptions explicitly. */
	private static String type(final Throwable error){
		return error==null ? "none" : error.getClass().getName();
	}

	/** Returns the top application frame without pinning a line number or JDK message. */
	private static String topFrame(final Throwable error){
		if(error==null){return "none";}
		final StackTraceElement[] trace=error.getStackTrace();
		if(trace.length<1){return "no-stack";}
		final StackTraceElement top=trace[0];
		return top.getClassName()+"."+top.getMethodName();
	}

	/** Rejects fixture misuse. */
	private static void require(final boolean correct, final String message){
		if(!correct){throw new IllegalArgumentException(message);}
	}

	/** Builds an expectation with no message token and no reported stream errors. */
	private static Expected e(final String out, final String err, final String exc,
			final int fd){
		return new Expected(out, err, exc, "", false, false, fd);
	}

	/** Builds an expectation with explicit message-token and stream-error requirements. */
	private static Expected e(final String out, final String err, final String exc,
			final String token, final boolean outErr, final boolean errErr,
			final int fd){
		return new Expected(out, err, exc, token, outErr, errErr, fd);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	/** Output stream for observing whether PrintStream errors are reported. */
	private static final class ThrowingOutputStream extends OutputStream{
		ThrowingOutputStream(final String name_){name=name_;}

		@Override
		public void write(final int b) throws IOException{
			throw new IOException("Injected "+name+" failure");
		}

		private final String name;
	}

	/** Held descriptor probe result. */
	private static final class DescriptorProbe{
		DescriptorProbe(final int before_, final int during_, final int after_){
			before=before_;
			during=during_;
			after=after_;
		}

		boolean matches(){return before==0 && during==1 && after==0;}

		final int before, during, after;
	}

	/** Captured direct-main outcome. */
	private static final class MainResult{
		MainResult(final String stdout_, final String stderr_,
				final boolean stdoutError_, final boolean stderrError_,
				final Throwable failure_){
			stdout=stdout_;
			stderr=stderr_;
			stdoutError=stdoutError_;
			stderrError=stderrError_;
			failure=failure_;
		}

		final String stdout, stderr;
		final boolean stdoutError, stderrError;
		final Throwable failure;
	}

	/** Expected direct-main result for one phase. */
	private static final class Expected{
		Expected(final String stdout_, final String stderr_, final String exception_,
				final String messageToken_, final boolean stdoutError_,
				final boolean stderrError_, final int fdDelta_){
			stdout=stdout_;
			stderr=stderr_;
			exception=exception_;
			messageToken=messageToken_;
			stdoutError=stdoutError_;
			stderrError=stderrError_;
			fdDelta=fdDelta_;
		}

		final String stdout, stderr, exception, messageToken;
		final boolean stdoutError, stderrError;
		final int fdDelta;
	}

	/** One literal case with current and desired one-row expectations. */
	private static final class CaseLiteral{
		static CaseLiteral byName(final String name){
			for(final CaseLiteral literal : CASES){
				if(literal.name.equals(name)){return literal;}
			}
			throw new IllegalArgumentException("Unknown case "+name);
		}

		CaseLiteral(final String name_, final String tool_, final String left_,
				final String right_, final boolean directoryRight_,
				final int differenceArgs_, final String[] terms_,
				final boolean failStdout_, final boolean failStderr_,
				final Expected old_, final Expected fixed_, final String description_,
				final boolean physicalLines_, final boolean logicalLines_,
				final boolean needsDiagnostic_, final boolean suspectedOutput_){
			name=name_;
			tool=tool_;
			left=left_;
			right=right_;
			directoryRight=directoryRight_;
			differenceArgs=differenceArgs_;
			terms=terms_;
			failStdout=failStdout_;
			failStderr=failStderr_;
			old=old_;
			fixed=fixed_;
			description=description_;
			physicalLines=physicalLines_;
			logicalLines=logicalLines_;
			needsDiagnostic=needsDiagnostic_;
			suspectedOutput=suspectedOutput_;
		}

		String[] writeInputs(final Path root, final ArrayList<Path> inputs)
				throws IOException{
			if(tool.equals("Difference") && differenceArgs==0){return new String[0];}
			if(tool.equals("FindString") && left==null){return new String[0];}
			final Path leftPath=root.resolve(name+".left.txt");
			writeInput(leftPath, left);
			inputs.add(leftPath);
			if(tool.equals("Difference")){
				if(differenceArgs==1){return new String[] {leftPath.toString()};}
				final Path rightPath=root.resolve(name+".right.txt");
				if(directoryRight){
					Files.createDirectory(rightPath);
				}else{
					writeInput(rightPath, right);
					inputs.add(rightPath);
				}
				return differenceArgs==3 ? new String[] {leftPath.toString(),
					rightPath.toString(), "ignored"} :
					new String[] {leftPath.toString(), rightPath.toString()};
			}
			final String[] args=new String[terms.length+1];
			args[0]=leftPath.toString();
			System.arraycopy(terms, 0, args, 1, terms.length);
			return args;
		}

		final String name, tool, left, right, description;
		final String[] terms;
		final Expected old, fixed;
		final int differenceArgs;
		final boolean directoryRight, failStdout, failStderr;
		final boolean physicalLines, logicalLines, needsDiagnostic, suspectedOutput;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final Path root;
	private final String phase;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	private static final String AIOOBE="java.lang.ArrayIndexOutOfBoundsException";
	private static final String ASSERTION="java.lang.AssertionError";
	private static final String IAE="java.lang.IllegalArgumentException";
	private static final String RUNTIME="java.lang.RuntimeException";
	private static final String NONE="none";
	private static final String LINE1="Line 1:\na\nb\n\n";
	private static final String LINE2="Line 2:\nc\nd\n\n";
	private static final String LINE3="Line 3:\ne\nf\n\n";
	private static final String LINE4="Line 4:\ng\nh\n\n";
	private static final String LINE5="Line 5:\ni\nj\n\n";
	private static final String FOUR_DIFF=LINE1+LINE2+LINE3+LINE4;
	private static final String FIVE_DIFF=FOUR_DIFF+LINE5;
	private static final String DIFF_USAGE="Usage: Difference";
	private static final String FIND_USAGE="Usage: FindString";

	private static final String HEADER="case\tphase\tutility\tdescription"+
		"\tcomparesPhysicalLines\tcomparesLogicalLines"+
		"\tdescriptorProbeBefore\tdescriptorProbeDuring"+
		"\tdescriptorProbeAfter\tdescriptorProbeMatched"+
		"\tfdBefore\tfdAfter\tfdDelta\tpredictedFdDelta"+
		"\tdirectExceptionClass\tdirectExceptionMessage"+
		"\tdirectExceptionTopFrame\tstdoutCheckError\tstderrCheckError"+
		"\tstdoutChars\tstderrChars\tstdoutMatched\tstderrMatched"+
		"\texceptionMatched\tmessageMatched\tstdoutErrorMatched"+
		"\tstderrErrorMatched\tfdBeforeMatched\tfdDeltaMatched"+
		"\tdesiredCleanupMatched\tdesiredDiagnosticMatched"+
		"\toutputErrorUnreported\tlegacyFailureObserved\tpredictionMatched";

	private static final CaseLiteral[] CASES={
		new CaseLiteral("difference-equal", "Difference", "same\n", "same\n",
			false, 2, null, false, false, e("", "", NONE, 0),
			e("", "", NONE, 0), "equal physical lines", true, false, false,
			false),
		new CaseLiteral("difference-one-mismatch", "Difference", "same\na\n",
			"same\nb\n", false, 2, null, false, false,
			e("", "Line 2:\na\nb\n\n", NONE, 0),
			e("", "Line 2:\na\nb\n\n", NONE, 0),
			"one physical mismatch", true, false, false, false),
		new CaseLiteral("difference-four-mismatches", "Difference",
			"a\nc\ne\ng\n", "b\nd\nf\nh\n", false, 2, null,
			false, false, e("", FOUR_DIFF, NONE, 0),
			e("", FOUR_DIFF, NONE, 0), "four mismatches return normally",
			true, false, false, false),
		new CaseLiteral("difference-fifth-mismatch", "Difference",
			"a\nc\ne\ng\ni\n", "b\nd\nf\nh\nj\n", false, 2, null,
			false, false, e("", FIVE_DIFF, ASSERTION, 2),
			e("", FIVE_DIFF, ASSERTION, 0),
			"fifth mismatch asserts after printing", true, false, false, false),
		new CaseLiteral("difference-unequal-lengths", "Difference", "same\na\n",
			"same\n", false, 2, null, false, false, e("", "", ASSERTION, 2),
			e("", "", ASSERTION, 0), "unequal physical line counts assert",
			true, false, false, false),
		new CaseLiteral("difference-missing-filename", "Difference", null, null,
			false, 0, null, false, false, e("", "", AIOOBE, 0),
			e("", "", IAE, DIFF_USAGE, false, false, 0),
			"zero supplied filenames", true, false, true, false),
		new CaseLiteral("difference-one-filename", "Difference", "open\n", null,
			false, 1, null, false, false, e("", "", AIOOBE, 1),
			e("", "", IAE, DIFF_USAGE, false, false, 0),
			"one supplied filename", true, false,
			true, false),
		new CaseLiteral("difference-second-open-failure", "Difference", "open\n",
			null, true, 2, null, false, false, e("", "", RUNTIME, 1),
			e("", "", RUNTIME, 0),
			"second input is an existing directory", true, false, false,
			false),
		new CaseLiteral("difference-physical-blank", "Difference", "a\n\nb\n",
			"a\nx\nb\n", false, 2, null, false, false,
			e("", "Line 2:\n\nx\n\n", NONE, 0),
			e("", "Line 2:\n\nx\n\n", NONE, 0),
			"blank physical line is compared", true, false, false, false),
		new CaseLiteral("difference-ignored-extra", "Difference", "same\n",
			"same\n", false, 3, null, false, false, e("", "", NONE, 0),
			e("", "", NONE, 0), "third argument is ignored", true, false,
			false, false),
		new CaseLiteral("findstring-literal-hit", "FindString",
			"alpha\nprefix need.le suffix\nprefix needle suffix\n", null, false, 0,
			new String[] {"need."}, false, false,
			e("prefix need.le suffix\n", "", NONE, 0),
			e("prefix need.le suffix\n", "", NONE, 0), "literal substring in line interior",
			false, true, false, false),
		new CaseLiteral("findstring-overlap-or", "FindString",
			"first\nsecond\nfirst second\nnone\n", null, false, 0,
			new String[] {"first", "second"}, false, false,
			e("first\nsecond\nfirst second\n", "", NONE, 0),
			e("first\nsecond\nfirst second\n", "", NONE, 0),
			"OR terms print first-only, second-only and both-term lines once",
			false, true, false, false),
		new CaseLiteral("findstring-no-hit", "FindString", "Alpha\n", null,
			false, 0, new String[] {"alpha"}, false, false,
			e("", "", NONE, 0), e("", "", NONE, 0),
			"literal search is case-sensitive", false, true, false, false),
		new CaseLiteral("findstring-zero-terms", "FindString", "alpha\n", null,
			false, 0, new String[] {}, false, false, e("", "", NONE, 0),
			e("", "", NONE, 0), "zero terms emit no output", false, true,
			false, false),
		new CaseLiteral("findstring-empty-term-skips-blanks", "FindString",
			"\n alpha \n \n\t\nbeta\t\n", null, false, 0,
			new String[] {""}, false, false, e(" alpha \nbeta\t\n", "", NONE, 0),
			e(" alpha \nbeta\t\n", "", NONE, 0),
			"empty term hits retained logical lines and preserves spaces",
			false, true, false, false),
		new CaseLiteral("findstring-missing-filename", "FindString", null, null,
			false, 0, null, false, false, e("", "", AIOOBE, 0),
			e("", "", IAE, FIND_USAGE, false, false, 0),
			"missing filename", false, true, true,
			false),
		new CaseLiteral("difference-failing-stderr", "Difference", "a\n", "b\n",
			false, 2, null, false, true,
			e("", "", NONE, "", false, true, 0),
			e("", "", RUNTIME, "output", false, true, 0),
			"injected stderr write failure", true, false, false,
			true),
		new CaseLiteral("findstring-failing-stdout", "FindString", "needle\n",
			null, false, 0, new String[] {"needle"}, true, false,
			e("", "", NONE, "", true, false, 0),
			e("", "", RUNTIME, "output", true, false, 0),
			"injected stdout write failure", false, true, false,
			true)
	};
}
