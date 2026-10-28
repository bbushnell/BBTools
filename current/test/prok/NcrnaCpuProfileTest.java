package test.prok;

import prok.*;

import java.io.File;
import java.nio.charset.StandardCharsets;
import fileIO.ByteFile;

/** Disjoint CPU accounting, exact merge and diagnostic-output regression checks.
 * @author Ganyu */
public final class NcrnaCpuProfileTest {
	public static void main(String[] args){
		check(args.length==1 && args[0].startsWith("outdir="), "Expected fresh outdir");
		check(!NcrnaCpuProfile.enabled() && NcrnaCpuProfile.create(2)==null, "Absent property must leave CPU observation disabled");
		final String root=args[0].substring(7);check(new File(root).mkdir(), "Fresh fixture directory required");
		final FakeClock clock=new FakeClock();
		final NcrnaCpuProfile first=new NcrnaCpuProfile(2, clock), second=new NcrnaCpuProfile(2, clock), total=new NcrnaCpuProfile(2, clock);
		first.begin(0, NcrnaCpuProfile.INIT);clock.add(8);first.end();
		first.begin(0, NcrnaCpuProfile.SCAVENGE);clock.add(21);first.end();
		first.begin(2, NcrnaCpuProfile.SHARED_SCAN);clock.add(3);first.end();
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
			"GLOBAL\tFIND_PATH\t1\t12\tNA\tNA\tNA\tNA\tNA\tNA"};
		try{
			for(byte[] line; (line=in.nextLine())!=null;){
				final String text=new String(line, StandardCharsets.UTF_8);
				if(text.startsWith("#") || text.startsWith("family\t")){continue;}
				check(row<expected.length && text.equals(expected[row++]), "CPU deltas, invocations and activity must keep family/stage identity");
			}
		}finally{check(!in.close(), "Fixture output read failed");}
		check(row==8, "Each family has two disjoint stages; global shared scan, SSU phases and path appear once each");
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
		System.out.println("PASS default-off, disjoint per-family/shared CPU accounting, exact counters, post-join merge, duplicate merge and output collision rejection");
	}
	static final class FakeClock implements NcrnaCpuProfile.Clock{
		long cpu=100;
		public long cpu(){return cpu;}
		void add(long c){cpu+=c;}
	}
	static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
