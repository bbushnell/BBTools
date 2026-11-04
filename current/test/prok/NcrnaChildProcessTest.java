package test.prok;

import java.io.File;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.nio.charset.StandardCharsets;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicReference;
import fileIO.ByteFile;
import parse.PreParser;
import prok.NcrnaChildProcess;
import prok.NcrnaTextIO.Output;
import fileIO.FileFormat;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** Real process tests, including startup races, blocked stdin and concurrent shutdown.
 * Preserved controls for NcrnaChildProcess, which remains necessary for the
 * standalone fixed-library regeneration mode. New public tools use NcrnaProcess;
 * its tests do not cover this older owner. Keep these controls while the
 * standalone mode retains that owner, including any future direct-library use.
 * @author Raiden */
public final class NcrnaChildProcessTest {

	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaChildProcessTest.class, false);
		shared.Shared.setThreads(1);fileIO.ReadWrite.setZipThreads(1);
		String root=null, installation=null, mode="test";
		for(String arg:pp.args){
			if(arg.startsWith("outdir=")){root=arg.substring(7);}
			else if(arg.startsWith("installation=")){installation=arg.substring(13);}
			else if(arg.startsWith("mode=")){mode=arg.substring(5);}
			else{throw new IllegalArgumentException(arg);}
		}
		check(root!=null && installation!=null && new File(root).mkdir(), "Use a fresh fixture directory and the sole installed program tree");
		final NcrnaChildProcessTest test=new NcrnaChildProcessTest(root, installation);
		if(mode.equals("shutdown")){test.shutdownChild();return;}
		check(mode.equals("test"), "Unknown fixture mode: "+mode);
		test.basic();test.cancel(false);test.cancel(true);test.heldInput();test.shutdown();
		System.out.println("NCRNA_CHILD_PROCESS_PASS exact_argv native_streaming exit_propagation startup_cancellation repeated_interruption held_stdin concurrent_shutdown");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private NcrnaChildProcessTest(String root_, String installation_){
		root=root_;installation=installation_;runner=new NcrnaChildProcess(new File(installation));
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private void basic() throws Exception{
		new NcrnaChildProcess(new File(installation), new File(root)).run(FIXTURE, new String[]{"mode=cwd"}, null, new File(root, "cwd"), 1, 64);
		check(read(root+"/cwd.stdout").equals(new File(root).getCanonicalPath()+"\n"),
			"Child resource resolution must use the declared working directory");
		check(!new File(root+"/cwd.cwd.txt").exists() && !new File(root+"/cwd.argv.tsv").exists(),
			"Child launch must not create removed argument or working-directory metadata");
		write(root+"/one.txt", "one\n");write(root+"/two.txt.gz", "two\n");
		final NcrnaChildProcess.Result ok=runner.run(FIXTURE, new String[]{"mode=echo", "payload=literal space; $(not-a-command)"},
			new String[]{root+"/one.txt", root+"/two.txt.gz"}, new File(root, "echo"), 1, 64);
		check(ok.exitCode==0 && ok.stdinBytes==8, "Child must consume both decompressed shards without truncation");
		check(read(root+"/echo.stdout").equals("literal space; $(not-a-command)\none\ntwo\n"), "Shell-looking payload must remain literal");
		check(read(root+"/echo.stderr").contains("stderr captured"), "Separate stderr must reach its file");
		boolean failed=false;
		try{runner.run(FIXTURE, new String[]{"mode=fail"}, null, new File(root, "failure"), 1, 64);}
		catch(IOException expected){failed=expected.getMessage().contains("exited 7");}
		check(failed, "A nonzero child exit must throw its actual exit7 IOException");
		check(!new File(root+"/failure.exit.tsv").exists() && !new File(root+"/echo.exit.tsv").exists(),
			"Neither failed nor successful children may write removed completion logs");
	}

	/** Stop once before ownership adoption, once while a descendant ignores TERM. */
	private void cancel(boolean duringStartup) throws Exception{
		final String label=duringStartup ? "startup" : "cancelled", prefix=root+"/"+label, descendant=prefix+".descendant.pid";
		final CountDownLatch entered=new CountDownLatch(1), hold=new CountDownLatch(1);
		final AtomicReference<Throwable> result=new AtomicReference<Throwable>();
		final NcrnaChildProcess.StartupAction observer=duringStartup ? ()->{entered.countDown();hold.await();} : null;
		final Thread task=new Thread(()->{
			try{runner.run(FIXTURE, new String[]{"mode=ignoreterm", "childfile="+descendant}, null, new File(prefix), 1, 64, observer);
				result.set(new AssertionError("Cancelled child returned success"));
			}catch(Throwable failure){result.set(failure);}
		}, label+"-fixture");
		task.setDaemon(true);task.start();
		try{
			if(duringStartup){check(entered.await(30, TimeUnit.SECONDS), "Parent did not reach startup observer");}
			awaitFile(duringStartup ? prefix+".pid" : descendant);
			task.interrupt();
			for(int i=0; i<20 && task.isAlive(); i++){Thread.sleep(25);task.interrupt();}
			task.join(15000);
			check(!task.isAlive() && result.get() instanceof InterruptedException, "Cancellation must finish and propagate interruption: "+result.get());
			check(!new File(prefix+".exit.tsv").exists(), "Cancellation must not create a removed completion log");
			if(duringStartup){check(!new File(descendant).exists(), "Handshake must prevent workload before parent ownership");}
			else{assertGone(read(descendant).trim(), prefix);}
		}finally{task.interrupt();hold.countDown();task.join(15000);}
	}

	private void heldInput() throws Exception{
		final String prefix=root+"/held", descendant=prefix+".descendant.pid";
		try(Output out=new Output(FileFormat.testOutput(root+"/large.txt.gz", FileFormat.TEXT, null, true, false, false, false))){
			final ByteBuilder block=new ByteBuilder();for(int i=0; i<65536; i++){block.append('a');}
			for(int i=0; i<32; i++){out.writer.print(block);}
		}
		final AtomicReference<Throwable> result=new AtomicReference<Throwable>();
		final Thread task=new Thread(()->{
			try{runner.run(FIXTURE, new String[]{"mode=holdstdin", "childfile="+descendant}, new String[]{root+"/large.txt.gz"}, new File(prefix), 1, 64);
				result.set(new AssertionError("Early exit with incomplete stdin returned success"));
			}catch(Throwable failure){result.set(failure);}
		}, "held-input-fixture");
		task.setDaemon(true);task.start();
		try{
			task.join(15000);check(!task.isAlive(), "Exited launcher with inherited stdin must not hang the input pump");
			check(result.get() instanceof IOException, "Truncated input must fail loudly: "+result.get());
			assertGone(read(descendant).trim(), prefix);
		}finally{task.interrupt();task.join(15000);}
	}

	/** Separate JVM receives TERM while its cleanup hook and interrupted main race. */
	private void shutdown() throws Exception{
		final String dir=root+"/shutdown";
		final Process child=new ProcessBuilder("bash", installation+"/dev/testncrnachild.sh", "mode=shutdown", "outdir="+dir)
			.redirectOutput(new File(root, "shutdown.stdout")).redirectError(new File(root, "shutdown.stderr")).start();
		try{
			awaitFile(dir+"/descendant.pid");awaitFile(dir+"/jvm.pid");
			final String pid=read(dir+"/jvm.pid").trim();check(pid.matches("[1-9][0-9]*"), "Signal target must be the fixture JVM PID");
			final Process signal=new ProcessBuilder("kill", "-TERM", "--", pid).start();
			check(signal.waitFor(5, TimeUnit.SECONDS) && signal.exitValue()==0, "Fixture JVM must receive TERM");
			check(child.waitFor(15, TimeUnit.SECONDS), "Concurrent cleanup callers must join rather than abandon cleanup or deadlock");
			assertGone(read(dir+"/descendant.pid").trim(), dir+"/child");
		}finally{
			if(child.isAlive()){
				child.destroy();
				if(!child.waitFor(15, TimeUnit.SECONDS)){child.destroyForcibly();check(child.waitFor(5, TimeUnit.SECONDS), "Fixture JVM survived forced cleanup");}
			}
		}
	}

	private void shutdownChild() throws Exception{
		final Thread main=Thread.currentThread();
		Runtime.getRuntime().addShutdownHook(new Thread(()->{
			for(int i=0; i<20; i++){main.interrupt();try{Thread.sleep(25);}catch(InterruptedException ignored){}}
		}, "fixture-concurrent-interrupt"));
		final String jvm=ManagementFactory.getRuntimeMXBean().getName();
		write(root+"/jvm.pid", jvm.substring(0, jvm.indexOf('@'))+"\n");
		runner.run(FIXTURE, new String[]{"mode=ignoreterm", "childfile="+root+"/descendant.pid"}, null, new File(root, "child"), 1, 64);
		throw new AssertionError("Shutdown fixture workload must never complete normally");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private static void awaitFile(String path) throws InterruptedException{
		final long end=System.nanoTime()+TimeUnit.SECONDS.toNanos(30);
		while(!new File(path).isFile() && System.nanoTime()<end){Thread.sleep(20);}
		check(new File(path).isFile(), "Fixture did not publish required identity: "+path);
	}
	private static void assertGone(String pid, String prefix) throws Exception{
		check(pid.matches("[1-9][0-9]*"), "Expected the owned descendant's numeric PID");
		String state="";
		for(int i=0; i<50; i++){
			final Process ps=new ProcessBuilder("ps", "-p", pid, "-o", "stat=").redirectOutput(new File(prefix+".stat")).redirectError(new File(prefix+".ps.err")).start();
			check(ps.waitFor(5, TimeUnit.SECONDS), "Process status check exceeded its fixture bound");
			state=read(prefix+".stat").trim();if(state.isEmpty() || state.startsWith("Z")){return;}Thread.sleep(20);
		}
		throw new AssertionError("A TERM-ignoring descendant survived group cleanup: "+state);
	}
	private static void write(String path, String content){
		try(Output out=new Output(FileFormat.testOutput(path, FileFormat.TEXT, null, true, false, false, false))){out.writer.print(new ByteBuilder(content));}
	}
	private static String read(String path){
		final ByteFile input=RrnaResourceIO.open(path);final ByteBuilder data=new ByteBuilder();
		try{for(byte[] line; (line=input.nextLine())!=null;){data.append(line).nl();}}
		finally{check(!input.close(), "Fixture output read failed: "+path);}return new String(data.toBytes(), StandardCharsets.UTF_8);
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final String root, installation;
	private final NcrnaChildProcess runner;
	private static final String FIXTURE="dev/ncrnachildfixture.sh";
}
