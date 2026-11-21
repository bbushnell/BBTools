package test.jgi;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import parse.Parse;

/** CountGC FQZ external-decoder registration controls.
 * @author Jean
 */
public final class CountGCFqzRegistrationTest{

	public static void main(final String[] args) throws Exception{
		for(final String arg : args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String key=split[0].toLowerCase();
			final String value=split.length>1 ? split[1] : null;
			if("phase".equals(key)){
				phase=value;
			}else if("outdir".equals(key)){
				outdir=value;
			}else if("scope".equals(key)){
				scope=value;
			}else{
				throw new AssertionError("Unknown argument: "+arg);
			}
		}

		check(phase!=null, "phase is required");
		check(outdir!=null, "outdir is required");
		final boolean baseline="baseline".equals(phase);
		check(baseline || "fixed".equals(phase), "phase must be baseline or fixed: "+phase);
		if(scope==null){scope="registration";}
		final boolean finalDecoder="finaldecoder".equals(scope);
		if(!finalDecoder && !"registration".equals(scope)){
			throw new AssertionError("scope must be registration or finaldecoder: "+scope);
		}
		final Path scratch=Files.createDirectory(Paths.get(outdir).toAbsolutePath());
		checkRawFastqControl(scratch);
		checkFqzRegistration(scratch, baseline);
		if(finalDecoder){checkFinalDecoderExitBaseline(scratch, baseline);}
		final String scopeText=finalDecoder ? " scope="+scope : "";
		System.out.println("PASS CountGCFqzRegistrationTest phase="+phase
				+scopeText+" checks="+checks);
	}

	/** Proves countgc.sh succeeds for the same one-read FASTQ without FQZ subprocesses. */
	private static void checkRawFastqControl(final Path scratch) throws Exception{
		final Path launcher=Paths.get("countgc.sh").toAbsolutePath().normalize();
		check(Files.isRegularFile(launcher), "countgc.sh not found at "+launcher);

		final Path evidence=Files.createDirectory(scratch.resolve("raw-fastq-control"));
		final Path input=evidence.resolve("input.fq");
		final Path output=evidence.resolve("countgc.out");
		final Path stdout=evidence.resolve("countgc.stdout");
		final Path stderr=evidence.resolve("countgc.stderr");
		Files.write(input, "@raw\nACGTN\n+\nIIIII\n".getBytes(StandardCharsets.US_ASCII));

		final int exit=run(evidence, "raw-fastq-control", stdout, stderr,
				"timeout", "--kill-after=5s", "30s",
				"bash", launcher.toString(), "-Xmx128m", "-Xms32m",
				"in="+input.toString(), "out="+output.toString(),
				"format=1", "ow=t", "threads=2");
		Files.write(evidence.resolve("countgc.exit"), Arrays.asList(Integer.toString(exit)),
				StandardCharsets.US_ASCII);

		final String stderrText=readString(stderr);
		check(exit==0, "raw FASTQ countgc exit="+exit
				+"\nstdout:\n"+readString(stdout)+"\nstderr:\n"+stderrText);
		check("".equals(readString(stdout)),
				"raw FASTQ countgc.sh wrote unexpected stdout: "+readString(stdout));
		check(stderrText.contains("jgi.CountGC") && stderrText.contains("-Xmx128m")
				&& stderrText.contains("in="+input.toString()),
				"raw FASTQ stderr lacked echoed Java command: "+stderrText);
		check(!stderrText.contains("Duplicate process for file"),
				"raw FASTQ stderr hit FQZ duplicate path: "+stderrText);
		check(rawRow().equals(readString(output)), "raw FASTQ output differed: "
				+readString(output));
	}

	/** Captures the default .fqz FileFormat probe and parser decoder ownership boundary. */
	private static void checkFqzRegistration(final Path scratch, final boolean baseline)
			throws Exception{
		final Path launcher=Paths.get("countgc.sh").toAbsolutePath().normalize();
		check(Files.isRegularFile(launcher), "countgc.sh not found at "+launcher);

		final Path evidence=Files.createDirectory(scratch.resolve("fqz-registration"));
		final Path bin=Files.createDirectory(evidence.resolve("bin"));
		final Path diag=Files.createDirectory(evidence.resolve("diag"));
		final Path marker=evidence.resolve("fqz.marker");
		final Path fakeFqz=bin.resolve("fqz_comp");
		final Path input=evidence.resolve("input.fqz");
		final Path output=evidence.resolve("countgc.out");
		final Path stdout=evidence.resolve("countgc.stdout");
		final Path stderr=evidence.resolve("countgc.stderr");
		Files.write(input, "not real fqz; fake decoder writes FASTQ\n"
				.getBytes(StandardCharsets.US_ASCII));
		writeFakeFqz(fakeFqz, 0);

		final int exit=runWithFakeFqz(evidence, "fqz-registration", stdout, stderr,
				bin, marker, diag,
				"timeout", "--kill-after=5s", "30s",
				"bash", launcher.toString(), "-Xmx128m", "-Xms32m",
				"in="+input.toString(), "out="+output.toString(),
				"format=1", "ow=t", "threads=2");
		Files.write(evidence.resolve("countgc.exit"), Arrays.asList(Integer.toString(exit)),
				StandardCharsets.US_ASCII);
		writeDiagnosticSnapshot(evidence.resolve("fqz.diag"), diag);
		validateFqzMarkers(marker, input, baseline ? -2 : 2, "registration");

		final String stderrText=readString(stderr);
		final String stdoutText=readString(stdout);
		final String outputText=readStringIfExists(output);
		check(exit==(baseline ? 1 : 0), "FQZ registration exit="+exit
				+"\nstdout:\n"+readString(stdout)+"\nstderr:\n"+stderrText
				+"\nmarker:\n"+readStringIfExists(marker));
		check("".equals(stdoutText),
				"FQZ countgc.sh wrote unexpected stdout: "+stdoutText);
		check(stderrText.contains("jgi.CountGC") && stderrText.contains("-Xmx128m")
				&& stderrText.contains("in="+input.toString()),
				"FQZ stderr lacked echoed Java command: "+stderrText);
		if(baseline){
			check(stderrText.contains("java.lang.Exception: Duplicate process for file "
					+input.toString()),
					"FQZ stderr lacked duplicate-process message: "+stderrText);
			check(stderrText.contains("fileIO.ReadWrite.addProcess"),
					"FQZ stderr lacked addProcess stack frame: "+stderrText);
			check(stderrText.contains("fileIO.FileFormat.getFirstOctet"),
					"FQZ stderr lacked getFirstOctet stack frame: "+stderrText);
			check(!outputText.contains("fqz\t"),
					"FQZ output unexpectedly contained parser row: "+outputText);
			check(!stderrText.contains(overallRow()),
					"FQZ stderr reached Overall row: "+stderrText);
			check(!stderrText.contains("\nTime: \t"),
					"FQZ stderr reached Time report: "+stderrText);
		}else{
			check(!stderrText.contains("Duplicate process for file"),
					"fixed FQZ stderr hit duplicate-process path: "+stderrText);
			check(fqzRow().equals(outputText), "fixed FQZ output differed: "+outputText);
			check(stderrText.contains(overallRow()),
					"fixed FQZ stderr lacked Overall row: "+stderrText);
			check(stderrText.contains("\nTime: \t"),
					"fixed FQZ stderr lacked Time report: "+stderrText);
		}
	}

	/** Captures the current final parser path ignoring a decoder that exits 7. */
	private static void checkFinalDecoderExitBaseline(final Path scratch,
			final boolean baseline) throws Exception{
		final Path launcher=Paths.get("countgc.sh").toAbsolutePath().normalize();
		check(Files.isRegularFile(launcher), "countgc.sh not found at "+launcher);

		final Path evidence=Files.createDirectory(scratch.resolve("fqz-final-exit7"));
		final Path bin=Files.createDirectory(evidence.resolve("bin"));
		final Path diag=Files.createDirectory(evidence.resolve("diag"));
		final Path marker=evidence.resolve("fqz.marker");
		final Path fakeFqz=bin.resolve("fqz_comp");
		final Path input=evidence.resolve("input.fqz");
		final Path output=evidence.resolve("countgc.out");
		final Path stdout=evidence.resolve("countgc.stdout");
		final Path stderr=evidence.resolve("countgc.stderr");
		Files.write(input, "not real fqz; fake decoder writes FASTQ then exits 7\n"
				.getBytes(StandardCharsets.US_ASCII));
		writeFakeFqz(fakeFqz, 7);

		final int exit=runWithFakeFqz(evidence, "fqz-final-exit7", stdout, stderr,
				bin, marker, diag,
				"timeout", "--kill-after=5s", "30s",
				"bash", launcher.toString(), "-Xmx128m", "-Xms32m",
				"in="+input.toString(), "out="+output.toString(),
				"format=1", "ow=t", "threads=2");
		Files.write(evidence.resolve("countgc.exit"), Arrays.asList(Integer.toString(exit)),
				StandardCharsets.US_ASCII);
		writeDiagnosticSnapshot(evidence.resolve("fqz.diag"), diag);
		validateFqzMarkers(marker, input, 2, "final exit7");
		validateFinalExitDiagnostics(diag);

		final String stderrText=readString(stderr);
		final String stdoutText=readString(stdout);
		final String outputText=readStringIfExists(output);
		check(exit==(baseline ? 0 : 1),
				"FQZ final decoder exit7 exit differed; exit="
				+exit+"\nstdout:\n"+stdoutText+"\nstderr:\n"+stderrText
				+"\nmarker:\n"+readStringIfExists(marker));
		check("".equals(stdoutText),
				"FQZ final decoder countgc.sh wrote unexpected stdout: "+stdoutText);
		if(baseline){
			check(fqzRow().equals(outputText),
					"FQZ final decoder baseline output differed: "+outputText);
			check(stderrText.contains(overallRow()),
					"FQZ final decoder stderr lacked Overall row: "+stderrText);
			check(stderrText.contains("\nTime: \t"),
					"FQZ final decoder stderr lacked Time report: "+stderrText);
			check(!stderrText.contains("decoder failed"),
					"FQZ final decoder baseline observed completion failure: "+stderrText);
		}else{
			check(fqzRow().equals(outputText),
					"FQZ final decoder fixed output differed: "+outputText);
			check(!stderrText.contains(overallRow()),
					"FQZ final decoder fixed stderr reached Overall row: "+stderrText);
			check(!stderrText.contains("\nTime: \t"),
					"FQZ final decoder fixed stderr reached Time report: "+stderrText);
			check(stderrText.contains("java.io.IOException: FQZ decoder failed for "
					+input.toString()),
					"FQZ final decoder fixed stderr lacked completion failure: "
					+stderrText);
		}
	}

	/** Fake FQZ decoder: finite, nonblocking, and records best-effort PID diagnostics. */
	private static void writeFakeFqz(final Path fakeFqz, final int successExit)
			throws IOException{
		final String script="#!/usr/bin/env bash\n"
				+"set -eo pipefail\n"
				+"marker=${FAKE_FQZ_MARKER:?}\n"
				+"diag=${FAKE_FQZ_DIAG_DIR:?}\n"
				+"pid=${BASHPID:-$$}\n"
				+"tab=$'\\t'\n"
				+"line=\"$pid${tab}$#\"\n"
				+"for arg in \"$@\"; do line=\"$line${tab}$arg\"; done\n"
				+"printf '%s\\n' \"$line\" >> \"$marker\"\n"
				+"printf 'pid=%s\\nargc=%s\\n' \"$pid\" \"$#\" > \"$diag/$pid.start\"\n"
				+"if [[ $# -eq 3 && $1 = -d && $2 = \"\" && -f $3 ]]; then\n"
				+"  role=final\n"
				+"  exitcode="+successExit+"\n"
				+"  if mkdir \"$diag/probe.lock\" 2>/dev/null; then\n"
				+"    role=probe\n"
				+"    exitcode=0\n"
				+"  fi\n"
				+"  printf 'role=%s\\n' \"$role\" >> \"$diag/$pid.start\"\n"
				+"  printf '@fqz\\nACGTN\\n+\\nIIIII\\n'\n"
				+"  printf 'role=%s\\nexit=%s\\n' \"$role\" \"$exitcode\" "
				+"> \"$diag/$pid.exit\"\n"
				+"  exit \"$exitcode\"\n"
				+"fi\n"
				+"printf 'unexpected fqz_comp argv' >&2\n"
				+"printf 'exit=64\\n' > \"$diag/$pid.exit\"\n"
				+"exit 64\n";
		Files.write(fakeFqz, script.getBytes(StandardCharsets.US_ASCII));
		check(fakeFqz.toFile().setExecutable(true, false),
				"Could not make fake fqz_comp executable at "+fakeFqz);
	}

	private static void validateFqzMarkers(final Path marker, final Path input,
			final int expected, final String label) throws IOException{
		final List<String> lines=Files.exists(marker)
				? Files.readAllLines(marker, StandardCharsets.US_ASCII)
				: Collections.<String>emptyList();
		if(expected<0){
			check(lines.size()<=2, "FQZ marker contained more than two records: "+lines);
		}else{
			check(lines.size()==expected, label+" FQZ marker count differed: "+lines);
		}
		for(final String line : lines){
			check(line.length()>0, "FQZ marker contained a blank record: "+lines);
			final String[] split=line.split("\t", -1);
			check(split.length==5, "FQZ marker field count differed: "+line);
			check(split[0].matches("[1-9][0-9]*"), "FQZ marker PID differed: "+line);
			check("3".equals(split[1]), "FQZ marker argc differed: "+line);
			check("-d".equals(split[2]), "FQZ marker first arg differed: "+line);
			check("".equals(split[3]), "FQZ marker empty second arg missing: "+line);
			check(input.toString().equals(split[4]),
					"FQZ marker input differed: "+line);
		}
	}

	private static void validateFinalExitDiagnostics(final Path diag) throws IOException{
		final File[] files=diag.toFile().listFiles();
		int probe=0;
		int parser=0;
		if(files!=null){
			for(final File file : files){
				if(file.getName().endsWith(".exit")){
					final String text=readString(file.toPath());
					if("role=probe\nexit=0\n".equals(text)){
						probe++;
					}else if("role=final\nexit=7\n".equals(text)){
						parser++;
					}else{
						throw new AssertionError("FQZ exit diagnostic differed in "
								+file+": "+text);
					}
				}
			}
		}
		check(probe==1, "FQZ probe exit diagnostic count differed: "+probe);
		check(parser==1, "FQZ final exit diagnostic count differed: "+parser);
	}

	private static int run(final Path scratch, final String stem,
			final Path stdout, final Path stderr, final String... command) throws Exception{
		Files.write(scratch.resolve(stem+".args"), Arrays.asList(command),
				StandardCharsets.US_ASCII);
		final ProcessBuilder pb=new ProcessBuilder(command);
		pb.directory(Paths.get(".").toAbsolutePath().normalize().toFile());
		pb.redirectOutput(stdout.toFile());
		pb.redirectError(stderr.toFile());
		final Process process=pb.start();
		Throwable primary=null;
		try{
			process.getOutputStream().close();
			if(!process.waitFor(40, TimeUnit.SECONDS)){
				throw new AssertionError("Owned countgc.sh supervisor exceeded 40 seconds; "
						+"inspect "+stderr);
			}
			return process.exitValue();
		}catch(Exception e){
			primary=e;
			throw e;
		}catch(Error e){
			primary=e;
			throw e;
		}finally{
			try{finishOwnedChild(process);}
			catch(RuntimeException e){
				if(primary==null){throw e;}
				else if(primary!=e){primary.addSuppressed(e);}
			}
			catch(Error e){
				if(primary==null){throw e;}
				else if(primary!=e){primary.addSuppressed(e);}
			}
		}
	}

	private static int runWithFakeFqz(final Path scratch, final String stem,
			final Path stdout, final Path stderr, final Path pathFirst, final Path marker,
			final Path diag, final String... command) throws Exception{
		Files.write(scratch.resolve(stem+".args"), Arrays.asList(command),
				StandardCharsets.US_ASCII);
		final ProcessBuilder pb=new ProcessBuilder(command);
		pb.directory(Paths.get(".").toAbsolutePath().normalize().toFile());
		pb.redirectOutput(stdout.toFile());
		pb.redirectError(stderr.toFile());
		final Map<String, String> env=pb.environment();
		final String path=env.get("PATH");
		env.put("PATH", pathFirst.toString()+File.pathSeparator+(path==null ? "" : path));
		env.put("FAKE_FQZ_MARKER", marker.toString());
		env.put("FAKE_FQZ_DIAG_DIR", diag.toString());
		final Process process=pb.start();
		Throwable primary=null;
		try{
			process.getOutputStream().close();
			if(!process.waitFor(40, TimeUnit.SECONDS)){
				throw new AssertionError("Owned countgc.sh supervisor exceeded 40 seconds; "
						+"inspect "+stderr);
			}
			return process.exitValue();
		}catch(Exception e){
			primary=e;
			throw e;
		}catch(Error e){
			primary=e;
			throw e;
		}finally{
			try{finishOwnedChild(process);}
			catch(RuntimeException e){
				if(primary==null){throw e;}
				else if(primary!=e){primary.addSuppressed(e);}
			}
			catch(Error e){
				if(primary==null){throw e;}
				else if(primary!=e){primary.addSuppressed(e);}
			}
		}
	}

	private static void finishOwnedChild(final Process process){
		if(!process.isAlive()){return;}
		process.destroy();
		if(waitOwnedChild(process, 8)){return;}
		process.destroyForcibly();
		if(!waitOwnedChild(process, 5)){
			throw new AssertionError("Owned countgc.sh supervisor did not terminate "
					+"after forced cleanup.");
		}
	}

	private static boolean waitOwnedChild(final Process process, final int seconds){
		final long deadline=System.nanoTime()+TimeUnit.SECONDS.toNanos(seconds);
		boolean interrupted=false;
		try{
			while(process.isAlive()){
				final long remaining=deadline-System.nanoTime();
				if(remaining<=0){return false;}
				try{return process.waitFor(remaining, TimeUnit.NANOSECONDS);}
				catch(InterruptedException e){interrupted=true;}
			}
			return true;
		}finally{
			if(interrupted){Thread.currentThread().interrupt();}
		}
	}

	private static void writeDiagnosticSnapshot(final Path out, final Path diag)
			throws IOException{
		final List<String> lines=new ArrayList<String>();
		final File[] files=diag.toFile().listFiles();
		if(files!=null){
			Arrays.sort(files);
			for(final File file : files){
				if(!file.isFile()){continue;}// The decoder's probe claim is a directory.
				final Path path=file.toPath();
				lines.add("## "+file.getName());
				lines.add(readString(path));
			}
		}
		Files.write(out, lines, StandardCharsets.US_ASCII);
	}

	private static String readString(final Path path) throws IOException{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}

	private static String readStringIfExists(final Path path){
		try{
			return Files.exists(path) ? readString(path) : "<missing>";
		}catch(IOException ioe){
			return "<read failure: "+ioe+">";
		}
	}

	private static String rawRow(){
		return "raw\t5\t0.25000\t0.25000\t0.25000\t0.25000\t0.20000\t0.50000\n";
	}

	private static String fqzRow(){
		return "fqz\t5\t0.25000\t0.25000\t0.25000\t0.25000\t0.20000\t0.50000\n";
	}

	private static String overallRow(){
		return "Overall\t5\t0.25000\t0.25000\t0.25000\t0.25000";
	}

	private static void check(final boolean condition, final String message){
		checks++;
		if(!condition){throw new AssertionError(message);}
	}

	private CountGCFqzRegistrationTest(){}

	private static String outdir;
	private static String phase;
	private static String scope;
	private static int checks;
}
