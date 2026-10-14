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

/** Finite A_SampleByteFile cases with its customization hook filled locally.
 * @author Sucrose
 */
public final class TemplateByteFileCase{

	/** Runs one patched-template case and checks observed output. */
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
		check("noout".equals(mode) || "ordinary".equals(mode) ||
			"preexisting".equals(mode) || "alias".equals(mode), "Unknown mode: "+mode);
		check("pass".equals(expect) || "npehang".equals(expect) ||
			"reject".equals(expect) || "allow".equals(expect), "Unknown expect: "+expect);
		check(!Files.exists(work), "Fresh work directory required: "+work);
		Files.createDirectories(work);

		final CaseSpec spec=writeCase(work.resolve("case"), mode);
		final RunResult result=runSample(classpath, work.resolve("case"), spec.args);
		checkResult(mode, expect, spec, result);
		System.out.println("PASS\t"+mode+"\t"+expect+"\t"+result.exit);
	}

	/** Writes the input files and returns A_SampleByteFile arguments for one case. */
	private static CaseSpec writeCase(final Path dir, final String mode) throws IOException{
		Files.createDirectories(dir);
		final Path input=dir.resolve("mixed.tsv");
		final Path valid=dir.resolve("valid.tsv");
		final Path invalid=dir.resolve("invalid.tsv");
		final ArrayList<String> args=new ArrayList<String>();
		if(mode.equals("alias")){
			Files.write(input, INVALID_ONLY);
			args.add("in="+input);
			args.add("invalid="+input);
			args.add("ow=t");
		}else{
			Files.write(input, MIXED_INPUT);
			args.add("in="+input);
			if(mode.equals("ordinary")){
				args.add("out="+valid);
				args.add("invalid="+invalid);
				args.add("ow=t");
			}else if(mode.equals("preexisting")){
				Files.write(invalid, PREEXISTING);
				args.add("invalid="+invalid);
				args.add("ow=f");
			}else{
				args.add("invalid="+invalid);
				args.add("ow=t");
			}
		}
		return new CaseSpec(input, valid, invalid, args);
	}

	/** Runs the patched template in a subprocess so each case has clean statics. */
	private static RunResult runSample(final String classpath, final Path dir,
			final ArrayList<String> args) throws Exception{
		final ArrayList<String> command=new ArrayList<String>();
		command.add(Paths.get(System.getProperty("java.home"), "bin", "java").toString());
		command.add("-ea");
		command.add("--add-modules");
		command.add("jdk.incubator.vector");
		command.add("-Xmx256m");
		command.add("-XX:ActiveProcessorCount=1");
		command.add("-cp");
		command.add(classpath);
		command.add("template.A_SampleByteFile");
		command.addAll(args);
		final Process process=new ProcessBuilder(command).directory(dir.toFile()).redirectErrorStream(true).start();
		final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		final IOException[] readError=new IOException[1];
		final Thread reader=new Thread(new Runnable(){
			@Override
			public void run(){
				try{
					final byte[] temp=new byte[4096];
					for(int len=process.getInputStream().read(temp); len>=0; len=process.getInputStream().read(temp)){
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
		Files.write(dir.resolve("child-output.txt"), result.output.getBytes(StandardCharsets.UTF_8));
		Files.write(dir.resolve("child-status.tsv"),
			("exit\ttimeout\n"+result.exit+'\t'+result.timeout+'\n').getBytes(StandardCharsets.UTF_8));
		return result;
	}

	/** Checks the case-specific exit and byte oracle. */
	private static void checkResult(final String mode, final String expect, final CaseSpec spec,
			final RunResult result) throws IOException{
		if(mode.equals("ordinary")){
			check("pass".equals(expect), "ordinary requires expect=pass");
			check(result.exit==0, "ordinary should pass: "+result.output);
			check(Arrays.equals(VALID_OUTPUT, Files.readAllBytes(spec.valid)), "valid output bytes");
			check(Arrays.equals(INVALID_OUTPUT, Files.readAllBytes(spec.invalid)), "invalid output bytes");
		}else if(mode.equals("noout")){
			if(expect.equals("npehang")){
				check(result.timeout, "noout baseline should leave the invalid writer unterminated");
				check(result.output.contains("NullPointerException"), "noout failure should be the valid-line NPE: "+result.output);
			}else{
				check("pass".equals(expect), "noout expects pass or npehang");
				check(result.exit==0, "noout should pass after repair: "+result.output);
				check(Arrays.equals(INVALID_OUTPUT, Files.readAllBytes(spec.invalid)), "noout invalid output bytes");
			}
		}else if(mode.equals("preexisting")){
			check("reject".equals(expect), "preexisting requires expect=reject");
			check(!result.timeout, "ow=f preexisting invalid output should reject without timing out");
			check(result.exit!=0, "ow=f preexisting invalid output should reject");
			check(Arrays.equals(PREEXISTING, Files.readAllBytes(spec.invalid)), "preexisting invalid output must be preserved");
		}else{
			if(expect.equals("allow")){
				check(result.exit==0, "baseline alias case should run until processing exposes missing validation: "+result.output);
			}else{
				check("reject".equals(expect), "alias expects allow or reject");
				check(!result.timeout, "alias of input and invalid output should reject without timing out");
				check(result.exit!=0, "alias of input and invalid output should reject");
				check(Arrays.equals(INVALID_ONLY, Files.readAllBytes(spec.input)), "alias input must be preserved");
			}
		}
	}

	/** Fixture check independent of the tested template's assertions. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final byte[] MIXED_INPUT="alpha\tignored\n#bad\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] INVALID_ONLY="#bad\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] VALID_OUTPUT="alpha\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] INVALID_OUTPUT="#bad\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] PREEXISTING="keep\n".getBytes(StandardCharsets.US_ASCII);

	private static final class CaseSpec{
		CaseSpec(final Path input_, final Path valid_, final Path invalid_, final ArrayList<String> args_){
			input=input_;
			valid=valid_;
			invalid=invalid_;
			args=args_;
		}
		final Path input;
		final Path valid;
		final Path invalid;
		final ArrayList<String> args;
	}

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
