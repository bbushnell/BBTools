package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Locale;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;

import dna.Data;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaProcess;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Cancels a native writer while its parent JVM's stdout pipe is deliberately unread.
 * Only a finite synthetic byte stream is attempted; no bulk output is saved.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaProcessBackpressureTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaProcessBackpressureTest.class, false);
		try(AutoCloseable log=()->Shared.closeStream(pp.outstream)){
			final Parser parser=new Parser();
			boolean child=false;
			for(String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(a.equals("child")){child=Parse.parseBoolean(b);}
				else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown backpressure option: "+arg);}
			}
			check(parser.out1!=null, "Specify a fresh out= directory for the backpressure fixture");
			final Path root=new File(parser.out1).toPath().toAbsolutePath();
			if(child){blockedParent(root);}
			else{
				Files.createDirectory(root);
				check(root);
				System.out.println("Backpressure passed: cancelled child and returned while parent stdout remained unread");
			}
		}
	}

	/** The outer process never drains the inner JVM's stdout, even during cancellation. */
	private static void check(final Path root) throws Exception{
		final String launcher=Data.findPath("?dev/testncrnaprocessbackpressure.sh", false);
		check(launcher!=null, "The fixture must use its installed native launcher");
		final Process parent=new ProcessBuilder("bash", launcher, "child=t", "out="+root)
			.redirectError(root.resolve("parent.stderr").toFile()).start();
		try{
			check(parent.waitFor(20, TimeUnit.SECONDS), "Parent did not finish cancellation with unread stdout");
			check(parent.exitValue()==0, "Blocked-output parent failed; inspect "+root.resolve("parent.stderr"));
			final String result=new String(Files.readAllBytes(root.resolve("result.txt")), StandardCharsets.US_ASCII);
			check(result.equals("cancelled\n"), "Inner JVM must finish and report the cancellation, not merely exit");
		}finally{
			if(parent.isAlive()){
				parent.destroy();
				if(!parent.waitFor(5, TimeUnit.SECONDS)){
					parent.destroyForcibly();
					check(parent.waitFor(5, TimeUnit.SECONDS), "Fixture parent survived forced cleanup");
				}
			}
			parent.getInputStream().close();
		}
	}

	/** Four MiB exceeds the ordinary Linux ProcessBuilder pipe capacity; no reader drains it. */
	private static void blockedParent(final Path root) throws Exception{
		check(Files.isDirectory(root), "Outer fixture owns the existing child output directory");
		final Path ready=root.resolve("writer_ready"), finished=root.resolve("writer_finished");
		final Path script=root.resolve("writer.sh");
		try(Output output=new Output(script.toString())){
			output.writer.print("#!/usr/bin/env bash\nset -e\nprintf 'ready\\n' > "+quote(ready.toString())+
				"\ndd if=/dev/zero bs=65536 count=64 2>/dev/null\nprintf 'finished\\n' > "+quote(finished.toString())+"\n");
		}
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(script.toString());
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final AtomicBoolean retained=new AtomicBoolean();
		final Thread worker=new Thread(()->{
			try{NcrnaProcess.run(command, null, "Blocked-output fixture");}
			catch(Throwable caught){failure.set(caught);}
			finally{retained.set(Thread.currentThread().isInterrupted());}
		}, "backpressure-caller");
		// A failed assertion must still let this test JVM exit and run its cleanup hook.
		worker.setDaemon(true);
		worker.start();
		try{
			final long deadline=System.nanoTime()+TimeUnit.SECONDS.toNanos(10);
			while(!Files.exists(ready) && worker.isAlive() && System.nanoTime()<deadline){Thread.sleep(10);}
			check(Files.exists(ready), "The native writer must start before cancellation");
			Thread.sleep(250);
			check(worker.isAlive() && !Files.exists(finished), "The unread pipe must prevent the finite writer from completing");
			worker.interrupt();
			worker.join(5000);
			check(!worker.isAlive() && failure.get()!=null && retained.get(),
				"Cancellation must return with interruption preserved while parent stdout is still unread");
			try(Output output=new Output(root.resolve("result.txt").toString())){output.writer.print("cancelled\n");}
		}finally{
			worker.interrupt();
			worker.join(1000);
		}
	}

	private static String quote(final String text){return "'"+text.replace("'", "'\\''")+"'";}
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
