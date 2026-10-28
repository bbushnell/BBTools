package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import fileIO.ByteFile;
import prok.CallGenesFileCost;
import prok.NcrnaCpuProfile;
import prok.RrnaResourceIO;

/** Deterministic clocks exercise file identity, joins, zero work and conservation.
 * Actual caller parity is exercised by dev/callgenes_file_cost_test.sh.
 * @author Ganyu */
public final class CallGenesFileCostTest {
	public static void main(String[] args){
		check(args.length==1 && args[0].startsWith("outdir="), "Supply a fresh outdir= fixture directory");
		check(CallGenesFileCost.start()==null, "Absent diagnostic property must not open outputs or clocks");
		final String root=args[0].substring(7);check(new File(root).mkdir(), "Fixture root must be fresh");
		final Clock clock=new Clock();final CallGenesFileCost report=new CallGenesFileCost(root+"/good", clock);
		try{
			final CallGenesFileCost.Span load=report.beginStartup("B12_RESOURCES");
			clock.process+=6;clock.wall+=8;report.endStartup(load, true);
			report.skippedStartup("DISABLED_EXAMPLE");
			report.prepare(new String[]{"same.fa", "same.fa", "empty.fa"}, new String[]{"A", "B"}, true);
			clock.wall=1000;clock.cpu=100;
			final CallGenesFileCost.Worker a=report.worker(0);a.start();
			final long begin=a.beginCalling();clock.wall+=10;clock.cpu+=7;a.endCalling(begin);
			final NcrnaCpuProfile init=new NcrnaCpuProfile(2, clock);
			init.begin(0, NcrnaCpuProfile.INIT);clock.wall+=2;clock.cpu+=2;init.end();a.stop(true);
			report.accept(a, 2, 40, new long[]{2, 1}, new long[]{10, 4}, ZERO, ZERO, new long[]{4, 1}, new long[]{2, 1}, init);
			boolean duplicate=false;try{report.accept(a, 2, 40, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, null);}catch(IllegalArgumentException e){duplicate=true;}
			check(duplicate, "Joining one worker twice must not inflate counts");
			final CallGenesFileCost.Worker b=report.worker(0);b.start();
			final long next=b.beginCalling();clock.wall+=5;clock.cpu+=3;b.endCalling(next);b.stop(true);
			report.accept(b, 1, 20, new long[]{3, 0}, new long[]{20, 0}, ZERO, ZERO, new long[]{6, 0}, new long[]{3, 0}, null);
			final CallGenesFileCost.Worker c=report.worker(1);c.start();clock.wall+=4;clock.cpu++;c.stop(true);
			report.accept(c, 0, 0, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, null);
			final CallGenesFileCost.Worker empty=report.worker(2);empty.start();empty.stop(true);
			boolean missing=false;try{report.finish(new long[]{5, 1});}catch(IllegalArgumentException e){missing=true;}
			check(missing, "An unjoined empty input must not silently disappear");
			report.accept(empty, 0, 0, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, null);
			boolean mismatch=false;try{report.finish(new long[]{6, 1});}catch(IllegalArgumentException e){mismatch=true;}
			check(mismatch, "File alignment totals must match the existing global caller counters");
			report.finish(new long[]{5, 1});
		}finally{report.close();}
		final ArrayList<String> rows=rows(root+"/good.inputs.tsv");
		check(rows.size()==4 && rows.get(1).equals("0\tsame.fa\t2\t3\t60\t12\t17\t15")
			&& rows.get(2).equals("1\tsame.fa\t1\t0\t0\t1\t4\t0")
			&& rows.get(3).equals("2\tempty.fa\t1\t0\t0\t0\t0\t0"),
			"Repeated paths retain distinct input ordinals; multi-worker CPU sums while input latency spans earliest to latest worker");
		final ArrayList<String> activity=rows(root+"/good.families.tsv");
		check(activity.size()==7 && activity.get(1).equals("0\tsame.fa\tA\t5\t30\t0\t0\t10\t5\t2")
			&& activity.get(2).equals("0\tsame.fa\tB\t1\t4\t0\t0\t1\t1\t0")
			&& activity.get(6).equals("2\tempty.fa\tB\t0\t0\t0\t0\t0\t0\t0"),
			"Actual counters and lazy INIT CPU conserve per-family identity, including empty-file zero rows");
		final ArrayList<String> startup=rows(root+"/good.startup.tsv");
		check(startup.size()==3 && startup.get(1).equals("B12_RESOURCES\tCOMPLETE\t8\t6")
			&& startup.get(2).equals("DISABLED_EXAMPLE\tSKIPPED_DISABLED\t0\t0"), "Shared-load process CPU is separate from worker CPU");
		boolean collision=false;try{new CallGenesFileCost(root+"/good", clock);}catch(IllegalArgumentException e){collision=true;}
		check(collision, "A diagnostic prefix cannot overwrite accepted evidence");
		final CallGenesFileCost unknown=new CallGenesFileCost(root+"/unknown", clock);
		try{
			unknown.prepare(new String[]{"empty.fa"}, new String[]{"A"}, false);
			final CallGenesFileCost.Worker w=unknown.worker(0);w.start();w.stop(true);
			unknown.accept(w, 0, 0, new long[1], new long[1], new long[1], new long[1], new long[1], new long[1], null);
			unknown.finish(new long[1]);
		}finally{unknown.close();}
		check(rows(root+"/unknown.families.tsv").get(1).endsWith("\tNA"), "Unmeasured lazy INIT must not be reported as zero CPU");
		System.out.println("CALLGENES_FILE_COST_TEST_PASS");
	}
	private static ArrayList<String> rows(String path){
		final ArrayList<String> result=new ArrayList<String>();final ByteFile in=RrnaResourceIO.open(path);
		try{for(byte[] line; (line=in.nextLine())!=null; ){if(line.length>0 && line[0]!='#'){result.add(new String(line, StandardCharsets.UTF_8));}}}
		finally{check(!in.close(), "Fixture report read failed");}return result;
	}
	private static final class Clock implements CallGenesFileCost.Clock, NcrnaCpuProfile.Clock{
		long wall=10, cpu=20, process=30;
		public long wall(){return wall;}public long threadCpu(){return cpu;}public long processCpu(){return process;}public long cpu(){return cpu;}
	}
	private static final long[] ZERO={0, 0};
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
