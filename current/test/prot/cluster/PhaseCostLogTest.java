package test.prot.cluster;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.AtomicReference;

import prot.cluster.PhaseCostLog;

/** Exact clock oracles for inclusive accounting, unavailable CPU and isolation.
 * The cake fixture also calls verifyReport on its actual instrumented runs.
 * @author Brian Bushnell, Fischl
 */
public final class PhaseCostLogTest{

	/** Runs only in the bounded Dori fixture suite; never overwrites evidence. */
	public static void main(final String[] args)throws Exception{
		check(args.length==1,"Supply one new phase-cost fixture directory");
		final Path root=Paths.get(args[0]);
		check(!Files.exists(root),"Phase-cost fixtures must preserve earlier results");
		Files.createDirectories(root);
		exactIntervals(root.resolve("nested.tsv"));
		unavailableCpu(root.resolve("unavailable.tsv"));
		closeGuards(root.resolve("guards.tsv"),root.resolve("forbidden.tsv"));
		try(PhaseCostLog.Span ignored=PhaseCostLog.begin("outside_context")){
			// Ordinary callers and worker threads have no reporting side effects.
		}
		try(PhaseCostLog log=new PhaseCostLog(root.resolve("real_clock.tsv").toString());
				PhaseCostLog.Span phase=PhaseCostLog.begin("jvm_clock")){
			check(!Files.exists(root.resolve("forbidden.tsv")),"A rejected nested recorder must not create its output");
		}
		verifyReport(root.resolve("real_clock.tsv"),"jvm_clock");
		System.out.println("PHASE_COST_LOG_TEST_PASS checks="+checks);
	}

	/** Child CPU can exceed wall time; parent CPU already includes that child. */
	private static void exactIntervals(final Path path)throws Exception{
		final FakeClock clock=new FakeClock();
		try(PhaseCostLog log=new PhaseCostLog(path.toString(),clock);
				PhaseCostLog.Span outer=PhaseCostLog.begin("outer")){
			clock.advance(5,20);
			try(PhaseCostLog.Span inner=PhaseCostLog.begin("inner")){
				final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
				final Thread worker=new Thread(()->{
					try(PhaseCostLog.Span ignored=PhaseCostLog.begin("worker_not_inherited")){}
					catch(Throwable problem){failure.set(problem);}
				});
				worker.start();worker.join();
				check(failure.get()==null,"A worker without the coordinator context must remain uninstrumented");
				clock.advance(7,40);
			}
			clock.advance(3,2);
		}
		final List<String> rows=Files.readAllLines(path,StandardCharsets.US_ASCII);
		check(rows.size()==3,"Only the coordinator's two intervals are emitted");
		check(rows.get(1).equals("1\t0\tinner\t7\t40\tmeasured"),"The child records its own parallel CPU interval");
		check(rows.get(2).equals("0\t-1\touter\t15\t62\tmeasured"),"The parent is inclusive, not the sum of reported intervals");
		verifyReport(path,"outer","inner");
	}

	/** Both permanent lack of support and a disappearing counter remain explicit NA. */
	private static void unavailableCpu(final Path path)throws Exception{
		final FakeClock clock=new FakeClock();
		clock.cpu=-1;
		try(PhaseCostLog log=new PhaseCostLog(path.toString(),clock);
				PhaseCostLog.Span outer=PhaseCostLog.begin("unsupported")){
			clock.wall+=9;
			clock.cpu=100;
			try(PhaseCostLog.Span child=PhaseCostLog.begin("became_unavailable")){
				clock.wall+=4;
				clock.cpu=-1;
			}
		}
		final List<String> rows=Files.readAllLines(path,StandardCharsets.US_ASCII);
		check(rows.get(1).equals("1\t0\tbecame_unavailable\t4\tNA\tunavailable"),"Lost CPU support must not become zero or a negative duration");
		check(rows.get(2).equals("0\t-1\tunsupported\t13\tNA\tunavailable"),"Unsupported CPU stays distinct from measured zero");
		verifyReport(path,"unsupported","became_unavailable");
	}

	/** Incorrect closure order and foreign-thread ownership fail before corrupting the tree. */
	private static void closeGuards(final Path path, final Path forbidden)throws Exception{
		final FakeClock clock=new FakeClock();
		try(PhaseCostLog log=new PhaseCostLog(path.toString(),clock)){
			final PhaseCostLog.Span outer=PhaseCostLog.begin("outer");
			final PhaseCostLog.Span inner=PhaseCostLog.begin("inner");
			reject(()->outer.close(),"An active child must close before its parent");
			reject(()->log.close(),"An active interval must close before its writer");
			reject(()->PhaseCostLog.begin("invalid\tname"),"Phase names cannot inject columns");
			reject(()->new PhaseCostLog(forbidden.toString(),clock),"Nested recorders cannot steal an active context");
			final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
			final Thread worker=new Thread(()->{
				try{inner.close();failure.set(new AssertionError("Foreign close unexpectedly succeeded"));}
				catch(IllegalStateException expected){}
				catch(Throwable problem){failure.set(problem);}
			});
			worker.start();worker.join();
			check(failure.get()==null,"A worker must not close the coordinator's interval");
			inner.close();outer.close();outer.close();
		}
		check(Files.readAllLines(path).size()==3,"Rejected operations and repeated close must not duplicate rows");
		verifyReport(path,"outer","inner");
	}

	/** Independently checks the report graph and direct-child conservation. All
	 * timings are inclusive, so summing parents and children would double-count.
	 */
	public static void verifyReport(final Path path, final String... required)throws Exception{
		final List<String> lines=Files.readAllLines(path,StandardCharsets.US_ASCII);
		check(lines.size()>1 && lines.get(0).equals(PhaseCostLog.HEADER),"A phase report needs the declared schema and real intervals");
		final HashMap<Integer,String[]> rows=new HashMap<Integer,String[]>();
		final HashMap<Integer,long[]> childSums=new HashMap<Integer,long[]>();
		final Set<String> names=new HashSet<String>();
		int roots=0;
		for(int i=1;i<lines.size();i++){
			final String[] row=lines.get(i).split("\t",-1);
			check(row.length==6,"Every phase row has six explicit fields");
			final int id=Integer.parseInt(row[0]),parent=Integer.parseInt(row[1]);
			final long wall=Long.parseLong(row[3]);
			check(id>=0 && parent>=-1 && parent<id && wall>=0 && rows.put(id,row)==null,"Timing IDs are unique and parents precede children");
			final boolean available=row[5].equals("measured");
			check(available ? Long.parseLong(row[4])>=0 : row[5].equals("unavailable") && row[4].equals("NA"),
				"CPU is a measured nonnegative counter difference or explicit NA");
			names.add(row[2]);
			if(parent<0){roots++;}
			else{
				long[] sums=childSums.get(parent);
				if(sums==null){sums=new long[2];childSums.put(parent,sums);}
				sums[0]=Math.addExact(sums[0],wall);
				if(available){sums[1]=Math.addExact(sums[1],Long.parseLong(row[4]));}
			}
		}
		check(roots==1,"Each fixture/arm report has one inclusive root interval");
		for(String[] row:rows.values()){
			final int id=Integer.parseInt(row[0]),parent=Integer.parseInt(row[1]);
			check(parent<0 || rows.containsKey(parent),"Every child must name a recorded parent");
			final long[] sums=childSums.get(id);
			if(sums!=null){
				check(sums[0]<=Long.parseLong(row[3]),"Disjoint direct-child wall intervals fit inside their inclusive parent");
				check(!row[5].equals("measured") || sums[1]<=Long.parseLong(row[4]),"Direct-child process CPU is already contained in the parent");
			}
		}
		for(String name:required){check(names.contains(name),"Actual execution must record phase "+name);}
	}

	/** Runtime misuse is rejected independently of whether Java assertions are enabled. */
	private static void reject(final Runnable action, final String reason){
		boolean rejected=false;
		try{action.run();}catch(IllegalArgumentException | IllegalStateException expected){rejected=true;}
		check(rejected,reason);
	}

	private static final class FakeClock implements PhaseCostLog.Clock{
		@Override
		public long wallNanos(){return wall;}
		@Override
		public long processCpuNanos(){return cpu;}
		void advance(final long elapsed, final long work){wall+=elapsed;cpu+=work;}
		long wall=100,cpu=1000;
	}

	private static void check(final boolean condition, final String reason){
		checks++;
		if(!condition){throw new AssertionError(reason);}
	}
	private static int checks;
	private PhaseCostLogTest(){}
}
