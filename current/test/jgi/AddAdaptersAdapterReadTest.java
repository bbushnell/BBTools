package test.jgi;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipException;

import parse.Parse;

/**
 * Two finite AddAdapters CLI cases for ignored adapter-reader gzip failures.
 * The malformed file has a valid empty gzip member and one flipped CRC byte so
 * construction opens the stream but reading reports the failure. The literal
 * adapter keeps baseline output generation from hiding the ignored read error
 * behind an empty-adapter rejection.
 * @author Jean
 */
public final class AddAdaptersAdapterReadTest{

	/** Requires a canonical {@code tree=...}, phase=baseline|fixed and a new {@code outdir=...}. */
	public static void main(final String[] args) throws Exception{
		String sourceTree=null, phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			if(pair[0].equals("tree")){sourceTree=pair[1];}
			else if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown adapter-reader fixture argument: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(sourceTree!=null && directory!=null, "Set tree and a new outdir.");
		final Path tree=Paths.get(sourceTree).toAbsolutePath().normalize();
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		require(safePath(tree) && safePath(root), "Paths must survive AddAdapters equals splitting and fit one evidence row.");
		require(Files.isRegularFile(tree.resolve("addadapters.sh")), "Missing actual AddAdapters launcher: "+tree);
		require(Files.isExecutable(Paths.get("/usr/bin/timeout")), "Owned child supervision requires GNU timeout.");
		Files.createDirectory(root);
		final AddAdaptersAdapterReadTest test=new AddAdaptersAdapterReadTest(tree, root, phase.equals("baseline"));
		test.scenario("valid-empty-gzip", false);
		test.scenario("corrupt-empty-gzip", true);
		test.check(test.cases==2 && test.children==2, "Frozen adapter-read matrix requires2 cases/2 terminal children: "+test.cases+"/"+test.children);
		if(test.failures>0){throw new AssertionError("AddAdapters adapter-read failed "+test.failures+" of "+test.checks+" checks; inspect "+root);}
		System.out.println("ADD_ADAPTERS_ADAPTER_READ_PASS phase="+phase+" cases="+test.cases+" children="+test.children+" checks="+test.checks);
	}

	private AddAdaptersAdapterReadTest(final Path tree_, final Path root_, final boolean baseline_){
		tree=tree_; root=root_; baseline=baseline_;
	}

	/** Launches AddAdapters on an owned one-read input and one generated gzip adapter file. */
	private void scenario(final String name, final boolean corrupt) throws Exception{
		final Path work=Files.createDirectory(root.resolve(name));
		final Path input=work.resolve("input.fq"), output=work.resolve("generated.fq"), adapter=work.resolve("adapter.fa.gz");
		final byte[] gzip=emptyGzip(corrupt);
		write(input, INPUT); write(work.resolve("input-before.fq"), INPUT);
		Files.write(adapter, gzip, StandardOpenOption.CREATE_NEW);
		write(work.resolve("adapter-before.hex"), hex(gzip)+"\n");
		oracle(work, gzip, corrupt);
		check(Files.isRegularFile(work.resolve("gzip-oracle.txt")), name+": gzip oracle must be recorded before launch.");
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("-ea", "-Xmx256m", "-Xms32m",
			"in="+input, "out="+output, "adapters="+adapter, "literal=GG", "mode=write",
			"reads=-1", "int=f", "qin=33", "qout=33", "simd=f", "t=1", "ow=t", "app=f",
			"arc=f", "adderrors=f", "rate=1", "right=t"));
		final int exit=run(args, work);
		final boolean sameInput=saveBytes(work, input, "input-after.fq", INPUT_BYTES);
		final boolean sameAdapter=saveBytes(work, adapter, "adapter-after.gz", gzip);
		final String diagnostics=read(work.resolve("stderr"));
		final boolean expectSuccess=!corrupt || baseline;
		write(work.resolve("observed.tsv"), "exit\tinputUnchanged\tadapterUnchanged\texpectedSuccess\n"+
			exit+"\t"+sameInput+"\t"+sameAdapter+"\t"+expectSuccess+"\n");
		if(Files.exists(work.resolve("timeout.txt"))){throw new AssertionError("Timeout/forced cleanup is never accepted: "+work);}
		check(sameInput, name+": original input bytes must survive the child unchanged.");
		check(sameAdapter, name+": adapter gzip bytes must survive the child unchanged.");
		check(Files.size(work.resolve("stdout"))==0, name+": file/statistics mode unexpectedly wrote stdout.");
		check(exit==(expectSuccess ? 0 : 1), name+": wrong Java exit, observed="+exit+", expectedSuccess="+expectSuccess);
		if(expectSuccess){
			check(matches(output, GENERATED), name+": successful generation must emit the literal-backed one-base read.");
			check(hasSuccessSummary(diagnostics), name+": successful generation must report adapter statistics.");
			if(corrupt){
				check(diagnostics.contains("Error reading FASTA input "+adapter) &&
					diagnostics.contains("Corrupt GZIP trailer"), name+": baseline corrupt case must expose the adapter read error.");
			}else{
				check(!diagnostics.contains("Error reading FASTA input") && !diagnostics.contains("Exception") &&
					!diagnostics.contains("AssertionError"), name+": valid empty gzip control must have clean diagnostics.");
			}
		}else{
			final Pattern primary=Pattern.compile("(?m)^(?:Exception in thread \\\"main\\\" )?java\\.lang\\.RuntimeException: "+
				Pattern.quote("Error reading adapter file: "+adapter)+"$");
			check(!Files.exists(output), name+": fixed corrupt rejection created data output.");
			check(!hasSuccessSummary(diagnostics), name+": fixed corrupt rejection printed a partial success summary.");
			check(primary.matcher(diagnostics).find(), name+": fixed corrupt rejection did not report the primary adapter-read error.");
		}
		cases++;
	}

	/** Generates one gzip member with no payload; the corrupt variant flips the first CRC trailer byte. */
	private static byte[] emptyGzip(final boolean corrupt) throws Exception{
		final ByteArrayOutputStream baos=new ByteArrayOutputStream();
		final GZIPOutputStream gzip=new GZIPOutputStream(baos);
		gzip.finish();
		gzip.close();
		final byte[] bytes=baos.toByteArray();
		if(corrupt){bytes[bytes.length-8]^=1;}
		return bytes;
	}

	/** Confirms, inside this fixture, that the bytes have the intended local gzip behavior. */
	private void oracle(final Path work, final byte[] bytes, final boolean corrupt) throws Exception{
		ZipException failure=null;
		byte[] decoded=null;
		try(final GZIPInputStream in=new GZIPInputStream(new ByteArrayInputStream(bytes))){
			decoded=readAll(in);
		}
		catch(ZipException e){failure=e;}
		if(corrupt){
			require(failure!=null && failure.getMessage()!=null && failure.getMessage().contains("Corrupt GZIP trailer"),
				work.getFileName()+": corrupt member must open but fail with a trailer ZipException while reading.");
		}else{
			require(failure==null && decoded!=null && decoded.length==0,
				work.getFileName()+": valid empty member must decode to zero bytes.");
		}
		write(work.resolve("gzip-oracle.txt"), failure==null ? "decoded-empty\n" : failure.toString()+"\n");
	}

	private static byte[] readAll(final GZIPInputStream in) throws Exception{
		final ByteArrayOutputStream out=new ByteArrayOutputStream();
		final byte[] buffer=new byte[64];
		for(int read=in.read(buffer); read>=0; read=in.read(buffer)){
			out.write(buffer, 0, read);
		}
		in.close();
		return out.toByteArray();
	}

	/** Uses only the supervisor handle returned by this launch; GNU timeout owns its isolated child group. */
	private int run(final List<String> args, final Path work) throws Exception{
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("/usr/bin/timeout", "--verbose",
			"--signal=TERM", "--kill-after=2s", "30s", "bash", tree.resolve("addadapters.sh").toString()));
		command.addAll(args);
		Files.write(work.resolve("command.argv"), command, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		final ProcessBuilder builder=new ProcessBuilder(command).directory(tree.toFile())
			.redirectOutput(work.resolve("stdout").toFile()).redirectError(work.resolve("stderr").toFile());
		final String inherited=builder.environment().get("JAVA_TOOL_OPTIONS");
		builder.environment().put("JAVA_TOOL_OPTIONS", (inherited==null ? "" : inherited+" ")+"-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		final Process process=builder.start();
		Throwable primary=null;
		boolean exceeded=false;
		try{
			process.getOutputStream().close();
			exceeded=!process.waitFor(45, TimeUnit.SECONDS);
		}catch(Exception e){primary=e; throw e;}
		catch(Error e){primary=e; throw e;}
		finally{
			try{finishOwnedChild(process);}
			catch(RuntimeException e){if(primary==null){throw e;}else{primary.addSuppressed(e);}}
			catch(Error e){if(primary==null){throw e;}else{primary.addSuppressed(e);}}
		}
		final int exit=process.exitValue();
		write(work.resolve("exit.txt"), exit+"\n");
		if(exceeded || exit==124 || exit==137 || read(work.resolve("stderr")).contains("timeout: sending signal")){
			write(work.resolve("timeout.txt"), "Owned CLI child exceeded its deadline or was killed.\n");
		}
		children++;
		return exit;
	}

	private static void finishOwnedChild(final Process process){
		if(!process.isAlive()){return;}
		process.destroy();
		if(waitOwnedChild(process, 5)){return;}
		process.destroyForcibly();
		if(!waitOwnedChild(process, 5)){throw new AssertionError("Owned AddAdapters CLI supervisor survived forced cleanup.");}
	}

	/** Interrupts cannot cancel the cleanup deadline; restores the flag after the owned process is reaped. */
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
		}finally{if(interrupted){Thread.currentThread().interrupt();}}
	}

	private static boolean saveBytes(final Path work, final Path path, final String name, final byte[] expected) throws Exception{
		if(!Files.isRegularFile(path)){write(work.resolve(name+".missing"), "Generated file is no longer a regular file.\n"); return false;}
		final byte[] actual=Files.readAllBytes(path);
		Files.write(work.resolve(name), actual, StandardOpenOption.CREATE_NEW);
		return Arrays.equals(actual, expected);
	}

	private static boolean matches(final Path path, final String expected) throws Exception{
		return Files.isRegularFile(path) && Arrays.equals(Files.readAllBytes(path), expected.getBytes(StandardCharsets.US_ASCII));
	}
	private static boolean hasSuccessSummary(final String text){
		return text.contains("Adapters Added:") || text.contains("Valid Output:");
	}
	private static String read(final Path path) throws Exception{return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);}
	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}
	private static boolean safePath(final Path path){
		final String text=path.toString();
		return text.indexOf('=')<0 && text.indexOf('\t')<0 && text.indexOf('\n')<0 && text.indexOf('\r')<0;
	}
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}
	private void check(final boolean condition, final String message){checks++; if(!condition){failures++; System.err.println("CHECK_FAILED "+message);}}
	private static String hex(final byte[] bytes){
		final char[] out=new char[bytes.length*2];
		for(int i=0; i<bytes.length; i++){
			final int value=bytes[i]&0xff;
			out[2*i]=HEX[value>>>4];
			out[2*i+1]=HEX[value&15];
		}
		return new String(out);
	}

	private static final char[] HEX="0123456789abcdef".toCharArray();
	private static final String INPUT="@ordinary\nG\n+\nI\n";
	private static final byte[] INPUT_BYTES=INPUT.getBytes(StandardCharsets.US_ASCII);
	private static final String GENERATED="@0_1_0\nG\n+\nI\n";
	private final Path tree, root;
	private final boolean baseline;
	private int checks, failures, cases, children;
}
