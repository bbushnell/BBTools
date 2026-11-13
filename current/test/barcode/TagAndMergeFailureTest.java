package test.barcode;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Locale;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPOutputStream;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Bounded process-level checks for TagAndMerge overwrite, append, and close
 * error propagation. Inputs are generated locally and no real service is used.
 */
public final class TagAndMergeFailureTest{

	/** Requires mode=baseline|fixed, case=..., {@code dir=<new directory>}, {@code installation=<BBTools root>}. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, TagAndMergeFailureTest.class, false);
		try{
			String mode=null, testCase=null, directory=null, installation=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT);
				final String value=pair.length>1 ? pair[1] : null;
				if("mode".equals(key)){mode=value;}
				else if("case".equals(key)){testCase=value;}
				else if("dir".equals(key)){directory=value;}
				else if("installation".equals(key)){installation=value;}
				else if(!Parser.parseStatic(arg, key, value)){
					throw new IllegalArgumentException("Unknown argument: "+arg);
				}
			}
			check("baseline".equals(mode) || "fixed".equals(mode),
					"Set mode=baseline or mode=fixed.");
			check(isCase(testCase), "Unknown case: "+testCase);
			check(directory!=null, "Set dir to a nonexistent output directory.");
			check(installation!=null, "Set installation to the BBTools root.");

			final Path dir=Files.createDirectory(Paths.get(directory));
			final CaseResult result=runCase(mode, testCase, dir, Paths.get(installation));
			System.out.println("TAGANDMERGE_FAILURE mode="+mode+" case="+testCase+
					" exit="+result.exit+" sentinel="+result.sentinelState+
					" stderrBytes="+result.stderr.length()+" failed=0");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private TagAndMergeFailureTest(){}

	private static CaseResult runCase(final String mode, final String testCase,
			final Path dir, final Path installation) throws Exception{
		if("normal-fresh".equals(testCase)){return normalFresh(mode, dir, installation);}
		if("valid-gzip".equals(testCase)){return validGzip(mode, dir, installation);}
		if("overwrite-ignored".equals(testCase)){return overwriteIgnored(mode, dir, installation);}
		if("append-ignored".equals(testCase)){return appendIgnored(mode, dir, installation);}
		if("output-dev-full".equals(testCase)){return outputDevFull(mode, dir, installation);}
		if("barcodes-dev-full".equals(testCase)){return barcodesDevFull(mode, dir, installation);}
		if("truncated-gzip".equals(testCase)){return truncatedGzip(mode, dir, installation);}
		throw new AssertionError(testCase);
	}

	private static CaseResult normalFresh(final String mode, final Path dir,
			final Path installation) throws Exception{
		final Path out=dir.resolve("merged.no_extension");
		final CaseResult result=runTagAndMerge(dir, installation, fastqInput(dir),
				"out="+out);
		checkSuccess(result);
		checkEquals(EXPECTED_FASTQ, read(out), "normal fresh output");
		return result.withSentinel("none");
	}

	private static CaseResult validGzip(final String mode, final Path dir,
			final Path installation) throws Exception{
		final Path out=dir.resolve("merged.fq");
		final CaseResult result=runTagAndMerge(dir, installation, gzipInput(dir),
				"out="+out);
		checkSuccess(result);
		checkEquals(EXPECTED_FASTQ, read(out), "valid gzip output");
		return result.withSentinel("none");
	}

	private static CaseResult overwriteIgnored(final String mode, final Path dir,
			final Path installation) throws Exception{
		final Path out=dir.resolve("preexisting.fq");
		write(out, SENTINEL);
		final CaseResult result=runTagAndMerge(dir, installation, fastqInput(dir),
				"out="+out, "ow=f");
		final String contents=read(out);
		if("baseline".equals(mode)){
			checkSuccess(result);
			check(contents.equals(EXPECTED_FASTQ),
					"baseline should overwrite despite ow=f; observed:\n"+contents);
			return result.withSentinel("overwritten");
		}
		check(result.exit!=0, "fixed ow=f should reject an existing output file.");
		checkEquals(SENTINEL, contents, "fixed ow=f must leave sentinel output intact");
		return result.withSentinel("preserved");
	}

	private static CaseResult appendIgnored(final String mode, final Path dir,
			final Path installation) throws Exception{
		final Path out=dir.resolve("append.fq");
		write(out, SENTINEL);
		final CaseResult result=runTagAndMerge(dir, installation, fastqInput(dir),
				"out="+out, "append=t");
		final String contents=read(out);
		if("baseline".equals(mode)){
			checkSuccess(result);
			checkEquals(EXPECTED_FASTQ, contents,
					"baseline should overwrite the sentinel despite append=t");
			return result.withSentinel("discarded");
		}
		checkSuccess(result);
		checkEquals(SENTINEL+EXPECTED_FASTQ, contents,
				"fixed append=t should retain the sentinel and append reads");
		return result.withSentinel("appended");
	}

	private static CaseResult outputDevFull(final String mode, final Path dir,
			final Path installation) throws Exception{
		requireDevFull();
		final CaseResult result=runTagAndMerge(dir, installation, fastqInput(dir),
				"out=/dev/full");
		checkDevFullStillDevice();
		checkMergedOutputWriteFailure(result.stderr);
		return expectIoFailure(mode, result, "merged output close");
	}

	private static CaseResult barcodesDevFull(final String mode, final Path dir,
			final Path installation) throws Exception{
		requireDevFull();
		final Path merged=dir.resolve("merged.fq");
		final CaseResult result=runTagAndMerge(dir, installation, fastqInput(dir),
				"out="+merged, "barcodes=/dev/full");
		checkDevFullStillDevice();
		checkBarcodeListWriteFailure(result.stderr);
		checkEquals(EXPECTED_FASTQ, read(merged),
				"barcode-list failure must preserve regular merged output");
		return expectIoFailure(mode, result, "barcode-list close");
	}

	private static CaseResult truncatedGzip(final String mode, final Path dir,
			final Path installation) throws Exception{
		final Path out=dir.resolve("merged.fq");
		final CaseResult result=runTagAndMerge(dir, installation, truncatedGzipInput(dir),
				"out="+out);
		checkTruncatedGzipReadFailure(result.stderr);
		return expectIoFailure(mode, result, "truncated gzip input close");
	}

	private static CaseResult expectIoFailure(final String mode, final CaseResult result,
			final String label){
		if("baseline".equals(mode)){
			checkSuccess(result);
		}else{
			check(result.exit!=0, "fixed "+label+" should fail; stderr:\n"+result.stderr);
			check(result.stderr.contains("barcode.TagAndMerge terminated in an error state"),
					"fixed "+label+" should propagate TagAndMerge I/O error; stderr:\n"+
					result.stderr);
		}
		return result.withSentinel("none");
	}

	private static CaseResult runTagAndMerge(final Path dir, final Path installation,
			final Path in, final String... extraArgs) throws Exception{
		final Path stdout=dir.resolve("tagandmerge.stdout");
		final Path stderr=dir.resolve("tagandmerge.stderr");
		final ArrayList<String> command=new ArrayList<String>();
		command.add(installation.resolve("dev").resolve("tagandmergefailuretest.sh").toString());
		command.add("--tag-child");
		command.add("in="+in);
		command.add("t=1");
		command.add("simd=f");
		command.add("interleaved=f");
		command.add("qin=33");
		command.add("qout=33");
		command.add("pigz=f");
		command.add("unpigz=f");
		for(final String arg : extraArgs){command.add(arg);}
		Files.write(dir.resolve("tagandmerge.command"),
				join(command).getBytes(StandardCharsets.UTF_8),
				StandardOpenOption.CREATE_NEW);
		final ProcessBuilder pb=new ProcessBuilder(command);
		pb.redirectOutput(stdout.toFile());
		pb.redirectError(stderr.toFile());
		final Process process=pb.start();
		try{
			final boolean finished=process.waitFor(TIMEOUT_SECONDS, TimeUnit.SECONDS);
			if(!finished){
				throw new AssertionError("TagAndMerge timed out after "+TIMEOUT_SECONDS+
						"s; timeout is not accepted as a close-error reproduction.");
			}
			final int exitValue=process.exitValue();
			write(dir.resolve("tagandmerge.exit"), Integer.toString(exitValue)+"\n");
			return new CaseResult(exitValue, read(stdout), read(stderr));
		}finally{
			if(process.isAlive()){
				process.destroyForcibly();
				check(process.waitFor(2, TimeUnit.SECONDS),
						"TagAndMerge child did not terminate after forced cleanup.");
			}
			check(!process.isAlive(), "TagAndMerge child was not reaped.");
		}
	}

	private static Path fastqInput(final Path dir) throws IOException{
		final Path path=dir.resolve(INPUT_FASTQ_NAME);
		write(path, INPUT_FASTQ);
		return path;
	}

	private static Path gzipInput(final Path dir) throws IOException{
		final Path path=dir.resolve(INPUT_GZIP_NAME);
		Files.write(path, gzip(INPUT_FASTQ.getBytes(StandardCharsets.UTF_8)),
				StandardOpenOption.CREATE_NEW);
		return path;
	}

	private static Path truncatedGzipInput(final Path dir) throws IOException{
		final byte[] full=gzip(INPUT_FASTQ.getBytes(StandardCharsets.UTF_8));
		final int length=Math.max(1, full.length-8);
		final byte[] truncated=new byte[length];
		System.arraycopy(full, 0, truncated, 0, length);
		final Path path=dir.resolve(INPUT_GZIP_NAME);
		Files.write(path, truncated, StandardOpenOption.CREATE_NEW);
		return path;
	}

	private static byte[] gzip(final byte[] data) throws IOException{
		final ByteArrayOutputStream baos=new ByteArrayOutputStream();
		final GZIPOutputStream gzip=new GZIPOutputStream(baos);
		gzip.write(data);
		gzip.close();
		return baos.toByteArray();
	}

	private static void requireDevFull() throws IOException{
		checkDevFullCharacterDevice("/dev/full must exist as a character device.");
	}

	private static void checkDevFullStillDevice() throws IOException{
		checkDevFullCharacterDevice("/dev/full was unexpectedly replaced.");
	}

	private static void checkDevFullCharacterDevice(final String label) throws IOException{
		final Path devFull=Paths.get("/dev/full");
		check(Files.exists(devFull), label);
		final int mode=((Integer)Files.getAttribute(devFull, "unix:mode")).intValue();
		check((mode&0170000)==0020000, label+" unix:mode="+Integer.toOctalString(mode));
	}

	private static void checkMergedOutputWriteFailure(final String stderr){
		check(stderr.contains("No space left") || stderr.contains("ENOSPC"),
				"merged output should emit a raw ENOSPC write/close diagnostic; stderr:\n"+
				stderr);
	}

	private static void checkBarcodeListWriteFailure(final String stderr){
		check(stderr.contains("No space left") || stderr.contains("ENOSPC") ||
				stderr.contains("IOException"),
				"barcode-list close should emit raw ENOSPC or its PrintWriter " +
				"IOException wrapper; stderr:\n"+stderr);
	}

	private static void checkTruncatedGzipReadFailure(final String stderr){
		check(stderr.contains("Error: truncated or corrupt input for '") &&
				stderr.contains("'; data may be incomplete."),
				"truncated gzip input should emit the distinctive ByteFile read " +
				"diagnostic; stderr:\n"+
				stderr);
	}

	private static void checkSuccess(final CaseResult result){
		check(result.exit==0, "Expected success, observed exit "+result.exit+
				"; stdout:\n"+result.stdout+"\nstderr:\n"+result.stderr);
	}

	private static boolean isCase(final String name){
		return "normal-fresh".equals(name) || "valid-gzip".equals(name) ||
				"overwrite-ignored".equals(name) || "append-ignored".equals(name) ||
				"output-dev-full".equals(name) || "barcodes-dev-full".equals(name) ||
				"truncated-gzip".equals(name);
	}

	private static String read(final Path path) throws IOException{
		return new String(Files.readAllBytes(path), StandardCharsets.UTF_8);
	}

	private static void write(final Path path, final String text) throws IOException{
		Files.write(path, text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	private static String join(final ArrayList<String> command){
		final StringBuilder sb=new StringBuilder();
		for(final String word : command){
			if(sb.length()>0){sb.append('\n');}
			sb.append(word);
		}
		sb.append('\n');
		return sb.toString();
	}

	private static void checkEquals(final String expected, final String observed,
			final String label){
		check(expected.equals(observed), label+" mismatch; expected:\n"+expected+
				"\nobserved:\n"+observed);
	}

	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final class CaseResult{
		CaseResult(final int exit_, final String stdout_, final String stderr_){
			this(exit_, stdout_, stderr_, "unset");
		}

		CaseResult(final int exit_, final String stdout_, final String stderr_,
				final String sentinelState_){
			exit=exit_;
			stdout=stdout_;
			stderr=stderr_;
			sentinelState=sentinelState_;
		}

		CaseResult withSentinel(final String sentinelState_){
			return new CaseResult(exit, stdout, stderr, sentinelState_);
		}

		final int exit;
		final String stdout;
		final String stderr;
		final String sentinelState;
	}

	private static final long TIMEOUT_SECONDS=10;
	private static final String TAG="GAGGCCGCCA+TTATCTAGCT";
	private static final String INPUT_FASTQ_NAME=
			"52866.4.475040.GAGGCCGCCA-TTATCTAGCT.filter-DNA.fastq";
	private static final String INPUT_GZIP_NAME=INPUT_FASTQ_NAME+".gz";
	private static final String INPUT_FASTQ="@SEQ0\nACGT\n+\nIIII\n";
	private static final String EXPECTED_FASTQ="@SEQ0\t"+TAG+"\nACGT\n+\nIIII\n";
	private static final String SENTINEL="@SENTINEL\nNNNN\n+\n####\n";
}
