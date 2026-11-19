package test.prok;

import prok.*;

import java.io.File;
import java.nio.charset.StandardCharsets;
import fileIO.ByteFile;

/** Disjoint CPU accounting, exact merge and diagnostic-output regression checks.
 * @author Ganyu, Raiden */
public final class NcrnaCpuProfileTest {
	public static void main(String[] args){
		final parse.PreParser pp=new parse.PreParser(args, NcrnaCpuProfileTest.class, false);
		try{run(pp.args);}finally{shared.Shared.closeStream(pp.outstream);}
	}
	private static void run(String[] args){
		check(args.length==1 && args[0].startsWith("outdir="), "Expected fresh outdir");
		check(!NcrnaCpuProfile.enabled() && NcrnaCpuProfile.create(2)==null, "Absent property must leave CPU observation disabled");
		final String root=args[0].substring(7);check(new File(root).mkdir(), "Fresh fixture directory required");
		final FakeClock clock=new FakeClock();
		final NcrnaCpuProfile first=new NcrnaCpuProfile(2, clock), second=new NcrnaCpuProfile(2, clock), total=new NcrnaCpuProfile(2, clock);
		first.begin(0, NcrnaCpuProfile.INIT);clock.add(8);first.end();
		first.begin(0, NcrnaCpuProfile.SCAVENGE);clock.add(21);first.end();
		first.begin(2, NcrnaCpuProfile.SHARED_SCAN);clock.add(3);first.end();
		first.begin(2, NcrnaCpuProfile.TRNA);
		clock.add(5);
		first.end();
		first.activity(0, 2, 100, 1, 44, 7, 3);
		first.verification(0, 14, 4, 2, 10, 3, 2);
		second.begin(1, NcrnaCpuProfile.INIT);clock.add(9);second.end();
		second.begin(1, NcrnaCpuProfile.SCAVENGE);clock.add(11);second.end();
		second.begin(2, NcrnaCpuProfile.SHARED_SCAN);clock.add(4);second.end();
		second.begin(2, NcrnaCpuProfile.SSU_APPLY);clock.add(6);second.end();
		second.begin(2, NcrnaCpuProfile.PATH);clock.add(12);second.end();
		second.activity(1, 5, 200, 2, 88, 9, 6);
		second.verification(1, 7, 2, 1, 8, 1, 1);
		total.merge(first);total.merge(second);
		final String path=root+"/cpu.tsv";total.write(path, new String[]{"familyA", "familyB"});
		final ByteFile in=RrnaResourceIO.open(path);int row=0;
		final String[] expected={
			"familyA\tINIT\t1\t8\tNA\tNA\tNA\tNA\tNA\tNA",
			"familyA\tSCAVENGE\t1\t21\t2\t100\t1\t44\t7\t3",
			"familyB\tINIT\t1\t9\tNA\tNA\tNA\tNA\tNA\tNA",
			"familyB\tSCAVENGE\t1\t11\t5\t200\t2\t88\t9\t6",
			"GLOBAL\tSHARED_SEED_SCAN\t2\t7\tNA\tNA\tNA\tNA\tNA\tNA",
			"GLOBAL\tSSU_APPLY\t1\t6\tNA\tNA\tNA\tNA\tNA\tNA",
			"GLOBAL\tSSU_FINISH\t0\t0\tNA\tNA\tNA\tNA\tNA\tNA",
			"GLOBAL\tFIND_PATH\t1\t12\tNA\tNA\tNA\tNA\tNA\tNA",
			"GLOBAL\tTRNA\t1\t5\tNA\tNA\tNA\tNA\tNA\tNA"};
		try{
			for(byte[] line; (line=in.nextLine())!=null;){
				final String text=new String(line, StandardCharsets.UTF_8);
				if(text.startsWith("#") || text.startsWith("family\t")){continue;}
				check(row<expected.length && text.equals(expected[row++]), "CPU deltas, invocations and activity must keep family/stage identity");
			}
		}finally{check(!in.close(), "Fixture output read failed");}
		check(row==9, "Each family has two stages; shared scan, SSU phases, path and independent tRNA appear once each");
		final ByteFile details=RrnaResourceIO.open(path+".verification.tsv");int detailRows=0;
		final String[] expectedDetails={"familyA\t14\t4\t2\t10\t3\t2", "familyB\t7\t2\t1\t8\t1\t1"};
		try{for(byte[] line; (line=details.nextLine())!=null; ){
			final String text=new String(line, StandardCharsets.UTF_8);if(text.startsWith("#") || text.startsWith("family\t")){continue;}
			check(detailRows<2 && text.equals(expectedDetails[detailRows++]), "Verifier subdivisions and considered/length-skip counts must merge without changing outer CPU spans");
		}}finally{check(!details.close(), "Verifier fixture output read failed");}
		check(detailRows==2 && NcrnaCpuProfile.cpuNow()==-1, "Every family has one detail row; disabled instrumentation must not sample a CPU clock");
		boolean duplicate=false;try{total.merge(first);}catch(IllegalArgumentException e){duplicate=true;}
		check(duplicate, "One worker cannot be merged twice");
		boolean overwrite=false;try{total.write(path, new String[]{"familyA", "familyB"});}catch(IllegalArgumentException e){overwrite=true;}
		check(overwrite, "Diagnostic output must fail loudly instead of overwriting prior evidence");
		windowHistograms(root);
		System.out.println("PASS default-off, disjoint per-family/shared CPU accounting, exact counters, post-join merge, duplicate merge and output collision rejection");
	}

	/** Unequal windows plus cache reuse distinguish per-window work from attempts
	 * and physical kernels. One deliberate outside request must remain visible. */
	private static void windowHistograms(final String root){
		final NcrnaAlignmentWork first=new NcrnaAlignmentWork(), second=new NcrnaAlignmentWork();
		first.window(0, 0);
		first.completed(20, 30, 2);
		first.completed(20, 30, 3);
		first.reused(20, 30);
		first.window(0, 0);
		first.completed(20, 30, 4);
		first.rankingQuery("ARN".getBytes(StandardCharsets.US_ASCII));
		second.reused(20, 30);
		second.window(0, 0);
		second.rankingQuery("acgtn".getBytes(StandardCharsets.US_ASCII));
		first.proposal(0, true);
		first.proposal(1, true);
		first.proposal(2, true);
		first.proposal(1, false);
		second.proposal(-1, true);
		second.proposal(2, false);
		final NcrnaAlignmentWork total=new NcrnaAlignmentWork();
		total.add(first);
		total.add(second);
		check(total.windows==3 && total.logical==5 && total.physical==3
			&& total.windowLogical==4 && total.windowPhysical==2,
			"Joined windows must retain logical3/physical2, logical1/physical0, zero work and one outside kernel");
		check(total.windowCount(true, 0)==1 && total.windowCount(true, 1)==1 && total.windowCount(true, 3)==1
			&& total.windowCount(false, 0)==2 && total.windowCount(false, 2)==1,
			"Primitive histogram bins must merge counts, not append worker bins as windows");
		check(total.windowPercentile(true, 50)==1 && total.windowPercentile(true, 95)==3
			&& total.windowPercentile(false, 50)==0 && total.windowPercentile(false, 95)==2,
			"Nearest-rank percentiles include the real zero-work window and exclude outside-window requests");
		final NcrnaCpuProfile report=new NcrnaCpuProfile(1, new FakeClock());
		report.alignmentWork(0).add(total);
		final String path=root+"/window_cpu.tsv";
		report.write(path, new String[]{"windowFamily"});
		final parse.LineParser1 p=new parse.LineParser1('\t');
		final ByteFile in=RrnaResourceIO.open(path+".windows.tsv");
		boolean found=false;
		try{
			for(byte[] row=in.nextLine(); row!=null; row=in.nextLine()){
				p.set(row);
				if(!p.termEquals("windowFamily", 0)){continue;}
				check(p.terms()==16 && p.parseLong(1)==3 && p.parseLong(2)==4 && p.parseLong(3)==2
					&& p.parseLong(4)==1 && p.parseLong(5)==1 && p.parseInt(7)==1 && p.parseInt(8)==3
					&& p.parseInt(9)==3 && p.parseInt(11)==0 && p.parseInt(12)==2 && p.parseInt(13)==2
					&& p.parseLong(14)==2 && p.parseLong(15)==1,
					"Written per-window quantiles and outside counts must preserve the literal merged fixture");
				check(Math.abs(Double.parseDouble(p.parseString(6))-4.0/3)<1e-12
					&& Math.abs(Double.parseDouble(p.parseString(10))-2.0/3)<1e-12,
					"Window means divide by the three actual windows, not five requests or worker count");
				found=true;
			}
		}finally{check(!in.close(), "Window fixture output read failed");}
		check(found, "Window work summary must retain its family row");
		final String[] proposals={
			"windowFamily\tnative\tBELOW_SHIPPED_FLOOR_2\t1",
			"windowFamily\tnative\tAT_LEAST_SHIPPED_FLOOR_2\t1",
			"windowFamily\tsketch\tBELOW_SHIPPED_FLOOR_2\t2",
			"windowFamily\tsketch\tAT_LEAST_SHIPPED_FLOOR_2\t1",
			"windowFamily\tsketch_only\tNA\t1"};
		final ByteFile accepts=RrnaResourceIO.open(path+".proposals.tsv");
		int acceptedRows=0;
		try{
			for(byte[] row; (row=accepts.nextLine())!=null;){
				if(row.length==0 || row[0]=='#'){continue;}
				p.set(row);
				if(p.termEquals("family", 0)){continue;}
				check(acceptedRows<proposals.length && proposals[acceptedRows].equals(new String(row, StandardCharsets.UTF_8)),
					"Merged accepted-event counts must preserve both sides of floor2 and sketch-only NA, without counting final calls");
				acceptedRows++;
			}
		}finally{check(!accepts.close(), "Proposal fixture output read failed");}
		check(acceptedRows==5, "All five observed stream/support classes must survive worker merge and report output");
	}
	static final class FakeClock implements NcrnaCpuProfile.Clock{
		long cpu=100;
		public long cpu(){return cpu;}
		void add(long c){cpu+=c;}
	}
	static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
