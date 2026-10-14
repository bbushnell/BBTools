package test.template;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/** Literal A_SampleByteFile maxLines and empty-line limit cases.
 * @author Sucrose
 */
public final class TemplateByteFileMaxLinesCase{

	/** Runs one maxLines case through the patched ByteFile template. */
	public static void main(final String[] args) throws Exception{
		Path work=null;
		String mode=null, expect=null, classpath=null;
		for(String arg : args){
			final String[] pair=arg.split("=", 2);
			check(pair.length==2, "Arguments require key=value: "+arg);
			if(pair[0].equals("classpath")){classpath=pair[1];}
			else if(pair[0].equals("work")){work=Paths.get(pair[1]).toAbsolutePath().normalize();}
			else if(pair[0].equals("mode")){mode=pair[1];}
			else if(pair[0].equals("expect")){expect=pair[1];}
			else{throw new IllegalArgumentException(arg);}
		}
		check(classpath!=null, "classpath= is required");
		check(work!=null, "work= is required");
		check("baseline".equals(expect) || "fixed".equals(expect),
				"expect=baseline or expect=fixed is required");
		check(!Files.exists(work), "Fresh work directory required: "+work);
		Files.createDirectories(work);

		final CaseSpec spec=CaseSpec.make(mode, "fixed".equals(expect));
		final Path dir=work.resolve(mode);
		Files.createDirectories(dir);
		final Path input=dir.resolve("input.tsv");
		final Path valid=dir.resolve("valid.tsv");
		final Path invalid=dir.resolve("invalid.tsv");
		Files.write(input, spec.input);

		final ArrayList<String> sampleArgs=new ArrayList<String>();
		sampleArgs.add("in="+input);
		sampleArgs.add("out="+valid);
		sampleArgs.add("invalid="+invalid);
		sampleArgs.add("lines="+spec.lines);
		sampleArgs.add("t=1");
		sampleArgs.add("ow=t");
		final RunResult result=runSample(classpath, dir, sampleArgs);
		check(result.exit==0, mode+" should pass: "+result.output);
		check(Arrays.equals(spec.valid, Files.readAllBytes(valid)), mode+" valid bytes");
		check(Arrays.equals(spec.invalid, Files.readAllBytes(invalid)), mode+" invalid bytes");
		check(count(result.output, "Valid Lines")==spec.validCount,
				mode+" valid count output: "+result.output);
		check(count(result.output, "Invalid Lines")==spec.invalidCount,
				mode+" invalid count output: "+result.output);
		System.out.println("PASS\t"+mode+"\t"+expect+"\t"+result.exit);
	}

	/** Runs A_SampleByteFile in a subprocess so each case has clean statics. */
	private static RunResult runSample(final String classpath, final Path dir,
			final ArrayList<String> args) throws Exception{
		final ArrayList<String> command=new ArrayList<String>();
		command.add(Paths.get(System.getProperty("java.home"), "bin", "java").toString());
		command.add("-ea");
		command.add("--add-modules");
		command.add("jdk.incubator.vector");
		final Path classLog=dir.resolve("child-class-load.txt");
		command.add("-Xlog:class+load=info:file="+classLog);
		command.add("-Xmx256m");
		command.add("-XX:ActiveProcessorCount=1");
		command.add("-cp");
		command.add(classpath);
		command.add("template.A_SampleByteFile");
		command.addAll(args);
		final Process process=new ProcessBuilder(command).directory(dir.toFile()).
				redirectErrorStream(true).start();
		final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		final IOException[] readError=new IOException[1];
		final Thread reader=new Thread(new Runnable(){
			@Override
			public void run(){
				try{
					final byte[] temp=new byte[4096];
					for(int len=process.getInputStream().read(temp); len>=0;
							len=process.getInputStream().read(temp)){
						buffer.write(temp, 0, len);
					}
				}catch(IOException e){readError[0]=e;}
			}
		});
		reader.start();
		Files.write(dir.resolve("child-command.txt"), command, StandardCharsets.UTF_8);
		final boolean timeout=!process.waitFor(20, TimeUnit.SECONDS);
		if(timeout){process.destroyForcibly(); process.waitFor();}
		reader.join(5000);
		check(!reader.isAlive(), "Child output reader did not finish after process exit");
		if(readError[0]!=null && !timeout){throw readError[0];}
		final RunResult result=new RunResult(timeout ? -1 : process.exitValue(),
				new String(buffer.toByteArray(), StandardCharsets.UTF_8), timeout);
		Files.write(dir.resolve("child-output.txt"),
				result.output.getBytes(StandardCharsets.UTF_8));
		Files.write(dir.resolve("child-status.tsv"),
				("exit\ttimeout\n"+result.exit+'\t'+result.timeout+'\n').
				getBytes(StandardCharsets.UTF_8));
		checkClassLoadLog(classpath, classLog);
		return result;
	}

	/** Verifies that the child loaded the sample from the prepared classes root. */
	private static void checkClassLoadLog(final String classpath, final Path classLog)
			throws IOException{
		final String[] entries=classpath.split(Pattern.quote(File.pathSeparator), 2);
		final Path classes=Paths.get(entries[0]).toAbsolutePath().normalize();
		final String log=new String(Files.readAllBytes(classLog), StandardCharsets.UTF_8);
		check(log.contains("template.A_SampleByteFile source: file:"+
				classes.toUri().getPath()), "A_SampleByteFile did not load from "+
				classes+" in "+classLog);
	}

	/** Parses one summary count from the sample's stdout/stderr stream. */
	private static long count(final String text, final String label){
		final Matcher matcher=Pattern.compile(label+":\\s+(\\d+)").matcher(text);
		check(matcher.find(), "Missing "+label+" in "+text);
		return Long.parseLong(matcher.group(1));
	}

	/** Fixture check independent of the tested template's assertions. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Input with two valid and two invalid nonempty lines. */
	private static final byte[] MIXED=("alpha\tignored\n#bad\n"+
			"beta\tignored\n#ugly\n").getBytes(StandardCharsets.US_ASCII);
	/** Input whose empty lines should not consume the nonempty-line budget. */
	private static final byte[] EMPTY_SPACED=("\nalpha\tignored\n\n#bad\n\n"+
			"beta\tignored\n").getBytes(StandardCharsets.US_ASCII);
	/** Both valid prefixes from MIXED. */
	private static final byte[] MIXED_VALID="alpha\nbeta\n".getBytes(StandardCharsets.US_ASCII);
	/** Both invalid lines from MIXED. */
	private static final byte[] MIXED_INVALID="#bad\n#ugly\n".getBytes(StandardCharsets.US_ASCII);
	/** Only the first valid prefix from MIXED. */
	private static final byte[] FIRST_VALID="alpha\n".getBytes(StandardCharsets.US_ASCII);
	/** Only the invalid line after the first valid line in EMPTY_SPACED. */
	private static final byte[] ONE_INVALID="#bad\n".getBytes(StandardCharsets.US_ASCII);
	/** Empty writer output. */
	private static final byte[] EMPTY=new byte[0];

	/** Literal oracle for one case. */
	private static final class CaseSpec{
		CaseSpec(final byte[] input_, final String lines_, final byte[] valid_,
				final byte[] invalid_, final long validCount_, final long invalidCount_){
			input=input_;
			lines=lines_;
			valid=valid_;
			invalid=invalid_;
			validCount=validCount_;
			invalidCount=invalidCount_;
		}

		static CaseSpec make(final String mode, final boolean fixed){
			if("zero".equals(mode)){
				return fixed ? new CaseSpec(MIXED, "0", EMPTY, EMPTY, 0, 0) :
						new CaseSpec(MIXED, "0", MIXED_VALID, MIXED_INVALID, 2, 2);
			}else if("one".equals(mode)){
				return new CaseSpec(MIXED, "1", FIRST_VALID, EMPTY, 1, 0);
			}else if("unlimited".equals(mode)){
				return new CaseSpec(MIXED, "-1", MIXED_VALID, MIXED_INVALID, 2, 2);
			}else if("empty".equals(mode)){
				return new CaseSpec(EMPTY_SPACED, "2", FIRST_VALID, ONE_INVALID, 1, 1);
			}
			throw new IllegalArgumentException("Unknown mode: "+mode);
		}

		final byte[] input;
		final String lines;
		final byte[] valid;
		final byte[] invalid;
		final long validCount;
		final long invalidCount;
	}

	/** Captured child process state. */
	private static final class RunResult{
		RunResult(final int exit_, final String output_, final boolean timeout_){
			exit=exit_;
			output=output_;
			timeout=timeout_;
		}
		final int exit;
		final String output;
		final boolean timeout;
	}
}
