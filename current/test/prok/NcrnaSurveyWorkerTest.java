package test.prok;

import java.io.File;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.PreParser;
import prok.NcrnaWorkers;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;
import template.ThreadWaiter;

/** Controlled worker startup and native output lifetime. The legacy route is a
 * helper-level negative control, not a historical full-survey reproduction.
 * Gates hold producers across startup failure; observations precede fixture cleanup.
 * @author Raiden */
public final class NcrnaSurveyWorkerTest{

	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaSurveyWorkerTest.class, false);
		try{runFixture(pp.args);}finally{Shared.closeStream(pp.outstream);}
	}
	private static void runFixture(String[] args) throws Exception{
		String root=null, expect=null;
		for(String arg:args){
			if(arg.startsWith("outdir=")){root=arg.substring(7);}
			else if(arg.startsWith("expect=")){expect=arg.substring(7);}
			else{throw new IllegalArgumentException("Unknown worker fixture option: "+arg);}
		}
		check(root!=null && ("legacy".equals(expect) || "repaired".equals(expect)), "Declare outdir and expect=legacy|repaired");
		check(new File(root).mkdir(), "Worker fixture requires a fresh output directory");
		final boolean repaired=expect.equals("repaired");
		final ByteBuilder report=new ByteBuilder("case\tstarted\tlate_writes\trows\tprimary\tsuppressed\tfinal_interrupt\n");
		one(root, "success", repaired, 0, false, false, report);
		one(root, "start_before", repaired, 1, false, false, report);
		one(root, "start_after", repaired, 2, false, false, report);
		one(root, "start_and_close", repaired, 2, true, false, report);
		one(root, "interrupted_join", repaired, 0, true, true, report);
		try(Output out=new Output(FileFormat.testOutput(root+"/observations.tsv", FileFormat.TEXT, null, true, false, false, false))){
			out.writer.print(report);
		}
		System.out.println("NCRNA_SURVEY_WORKER_PASS cases=5 expectation="+expect);
	}

	/** The driver's release determines exactly which writes must precede closure. */
	private static void one(String root, String label, boolean repaired, int startMode,
			boolean closeFailure, boolean interrupt, ByteBuilder report) throws Exception{
		final Scenario scenario=new Scenario(root+"/"+label+".tsv", repaired, startMode, closeFailure);
		try{
			scenario.start();
			await(scenario.entered, "All deliberately started producers must reach their write gate");
			if(startMode>0 && !repaired){
				join(scenario);
				check(scenario.writer.closed, "Legacy helper must expose its premature owner closure before release");
			}else{
				awaitJoin(scenario);
				if(interrupt){scenario.interrupt();awaitJoin(scenario);}
				check(!scenario.writer.closing, "A joined producer still owns a borrow until its gated write finishes");
			}
			scenario.release.countDown();
			join(scenario);
			for(Producer producer:scenario.workers){if(producer.getState()!=Thread.State.NEW){join(producer);}}
			check(scenario.writer.closed && !scenario.writer.isAlive(), "Native writer must terminate before observations");
			final int started=startMode==0 ? 3 : startMode==1 ? 1 : 2;
			final int expectedLate=startMode>0 && !repaired ? started : 0;
			int late=0;
			for(int i=0; i<3; i++){
				final Producer producer=scenario.workers.get(i);
				check(producer.problem==null, "Controlled producer failed: "+producer.problem);
				check(producer.getState()==(i<started ? Thread.State.TERMINATED : Thread.State.NEW), "Only started workers may terminate");
				if(producer.late){late++;}
			}
			check(late==expectedLate, "Premature closure must differ between legacy and repaired helpers: "+label);
			final Throwable expected=startMode>0 ? scenario.startFailure : closeFailure ? scenario.writer.cleanup : null;
			check(scenario.problem==expected, "Original startup/cleanup failure identity must survive: "+label+" "+scenario.problem);
			if(startMode>0){
				final Throwable[] suppressed=scenario.problem.getSuppressed();
				check(suppressed.length==(closeFailure ? 1 : 0), "Cleanup must remain underneath the original startup error");
				if(closeFailure){check(suppressed[0]==scenario.writer.cleanup, "The actual native close exception must be retained");}
			}
			check(scenario.finalInterrupt==(repaired && interrupt), "The owned cohort preserves interruption until native output cleanup finishes");
			check(descriptors(scenario.path)==0, "The owned native file must have no remaining descriptor");
			final int rows=rows(scenario.path, expectedLate>0 ? 0 : (1<<started)-1);
			report.append(label).tab().append(started).tab().append(late).tab().append(rows).tab()
				.append(expected==null ? "none" : expected.getClass().getName()).tab()
				.append(expected==null ? 0 : expected.getSuppressed().length).tab().append(scenario.finalInterrupt).nl();
		}finally{
			// Release even after an assertion. Daemon controllers cannot trap the test on a NEW-thread spin.
			scenario.release.countDown();
			for(Producer producer:scenario.workers){if(producer.getState()!=Thread.State.NEW){join(producer);}}
			scenario.join(WAIT_MS);
			if(scenario.writer!=null && !scenario.writer.closed){
				try{scenario.writer.poisonAndWait();}catch(IllegalStateException failure){check(failure==scenario.writer.cleanup, "Only the declared close fault is expected during teardown");}
			}
			check(!scenario.isAlive(), "Controller must finish without spinning on unstarted workers");
		}
	}
	private static void await(CountDownLatch latch, String why) throws InterruptedException{
		check(latch.await(WAIT_MS, TimeUnit.MILLISECONDS), why);
	}
	private static void join(Thread thread) throws InterruptedException{
		assert(thread.getState()!=Thread.State.NEW):"Fixture joins only actually started threads";
		thread.join(WAIT_MS);check(!thread.isAlive(), "Controlled thread failed to terminate: "+thread.getName());
	}
	/** WAITING before close, with producers gated, identifies the join boundary. */
	private static void awaitJoin(Scenario scenario){
		final long deadline=System.nanoTime()+TimeUnit.MILLISECONDS.toNanos(WAIT_MS);
		while(System.nanoTime()<deadline){
			check(scenario.isAlive() && !scenario.writer.closing, "Output closed before gated producers finished");
			if(scenario.getState()==Thread.State.WAITING && !scenario.isInterrupted()){return;}
			Thread.yield();
		}
		throw new AssertionError("Controller never reached its worker join");
	}
	private static int rows(String path, int expected){
		final ByteFile in=ByteFile.makeByteFile(FileFormat.testInput(path, FileFormat.TEXT, null, false, true));
		int mask=0, count=0;
		try{
			for(byte[] row; (row=in.nextLine())!=null;){
				check(row.length==1 && row[0]>='1' && row[0]<='3', "Native output must contain literal worker IDs");
				final int bit=1<<(row[0]-'1');check((mask&bit)==0, "Each producer writes exactly once");mask|=bit;count++;
			}
		}finally{check(!in.close(), "Native fixture input must close without error");}
		check(mask==expected, "Native output IDs must match the independently declared started cohort");return count;
	}
	private static int descriptors(String path) throws Exception{
		int count=0;
		try(DirectoryStream<Path> entries=Files.newDirectoryStream(new File("/proc/self/fd").toPath())){
			for(Path entry:entries){try{if(Files.readSymbolicLink(entry).toString().equals(path)){count++;}}
			catch(NoSuchFileException closedConcurrently){}}
		}return count;
	}

	/** Owns the real writer and invokes either the unchanged helper or the survey repair. */
	private static final class Scenario extends Thread{
		Scenario(String path_, boolean repaired_, int mode_, boolean close_) throws Exception{
			path=new File(path_).getAbsolutePath();repaired=repaired_;mode=mode_;closeFailure=close_;
			entered=new CountDownLatch(mode==0 ? 3 : mode==1 ? 1 : 2);
			for(int i=1; i<=3; i++){workers.add(new Producer(this, i));}
			setDaemon(true);
		}
		@Override public void run(){
			try{
				try(NcrnaWorkers cohort=repaired ? new NcrnaWorkers() : null;
						Output owned=openOutput()){
					check(owned.writer==writer, "Fixture uses the actual native owner and its borrowed writer");
					if(repaired){
						cohort.startAndJoin(workers, null);
					}else{ThreadWaiter.startAndWait(workers);}
				}
			}catch(Throwable failure){problem=failure;}
			finally{finalInterrupt=isInterrupted();}
		}
		/** Transfers the unstarted test writer directly to the shared output owner. */
		private Output openOutput(){
			writer=new TrackedWriter(path, closeFailure);
			return new Output(writer);
		}
		final String path;final boolean repaired, closeFailure;final int mode;
		final ArrayList<Producer> workers=new ArrayList<Producer>();
		final CountDownLatch entered, release=new CountDownLatch(1);
		final AssertionError startFailure=new AssertionError("declared worker startup failure");
		volatile TrackedWriter writer;volatile Throwable problem;volatile boolean finalInterrupt;
	}
	private static final class Producer extends Thread{
		Producer(Scenario scenario_, int id_){scenario=scenario_;id=id_;setDaemon(true);}
		@Override public synchronized void start(){
			if(id==2 && scenario.mode==1){throw scenario.startFailure;}
			super.start();
			if(id==2 && scenario.mode==2){throw scenario.startFailure;}
		}
		@Override public void run(){
			try{
				scenario.entered.countDown();scenario.release.await();
				synchronized(scenario.writer){
					late=scenario.writer.closing;
					// Never enqueue to the poisoned queue in the intentionally broken control.
					if(!late){scenario.writer.print(new ByteBuilder().append(id).nl());}
				}
			}catch(Throwable failure){problem=failure;}
		}
		final Scenario scenario;final int id;boolean late;Throwable problem;
	}
	private static final class TrackedWriter extends ByteStreamWriter{
		TrackedWriter(String path, boolean failure){super(FileFormat.testOutput(path, FileFormat.TEXT, null, true, false, false, false));failClose=failure;}
		@Override public boolean poisonAndWait(){
			closing=true;final boolean error=super.poisonAndWait();closed=true;
			if(failClose){throw cleanup;}return error;
		}
		final boolean failClose;volatile boolean closing, closed;
		final IllegalStateException cleanup=new IllegalStateException("declared native close failure");
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private static final long WAIT_MS=10000;
}
