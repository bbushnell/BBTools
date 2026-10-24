package test.prok;

import java.io.File;
import prok.MitoSsuTierLog;
import prok.NcrnaMitoWindowReport;

/** Saved-window projection conserves accepts and counts real verifier calls.
 * @author Brian, Ganyu */
public final class NcrnaMitoWindowReportTest {
	public static void main(String[] args)throws Exception{
		check(args.length==2 && args[0].startsWith("outdir=") && args[1].startsWith("sha80="), "Require fresh outdir and canonical sha80");
		final String root=args[0].substring(7), sha=args[1].substring(6);check(new File(root).mkdir(), "Fixture directory must be fresh");
		NcrnaMitoCrossLibraryTest.write(root+"/time.txt", "wall_seconds=2 user_seconds=3 system_seconds=1\n");
		final String a=row(85, 60, 40, 1, 1, 0, 7, 1, 1, 85, "ACCEPT", 1);
		final String b=row(120, 100, 90, 2, 3, 1, 8, 3, 2, 90, "ACCEPT", 3);
		final String c=row(150, 130, 100, 3, 5, 1, 9, 4, 3, 85, "ACCEPT", 4);
		final String d=row(300, 200, 100, 4, 6, 2, -1, -1, 0, 0, "REJECT_VERIFY", 3);
		final String e=row(3, 2, 1, 0, 0, 0, -1, -1, 0, 0, "REJECT_ENTRY", 0);
		write(root, a+b+c+d+e);run(root, bind(root, sha), "good");
		check(read(root+"/good.summary.tsv").contains("TEST\t5\t3\t2\t2\t1\t10\t15\t11\t4\t4\t2\t1\t85\t85\t"),
			"Cap two must preserve first/second attempts, lose late accepts, and sum actual prefix calls rather than assume two calls per model");
		final String histogram=read(root+"/good.top_shared.tsv");
		check(histogram.contains("TEST\t85\t1\t1\t") && histogram.contains("TEST\t300\t1\t0\t") && histogram.contains("TEST\t3\t1\t0\t"),
			"Exact top-shared histogram must distinguish all entered window events from accepted events, including zero-attempt windows");
		final String[] fields=read(root+"/good.summary.tsv").split("\n")[1].split("\t");
		check(fields.length==20 && Double.parseDouble(fields[15])==2 && Double.parseDouble(fields[18])==4, "Report actual wall and total CPU independently of the cap-two projection");
		write(root, "");run(root, bind(root, sha), "empty");
		check(read(root+"/empty.summary.tsv").contains("TEST\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tNA\tNA\t"), "No observed windows gives measured zero counts and unavailable score extrema");
		write(root, row(120, 100, 90, 2, 3, 0, 8, 2, 2, 90, "ACCEPT", 4));bad(root, bind(root, sha), "prefix_overflow");
		write(root, row(120, 100, 90, 2, 3, 0, 8, 2, 2, 90, "ACCEPT", 2));bad(root, bind(root, sha), "prefix_not_total");
		write(root, row(120, 100, 90, 2, 3, 0, 8, 3, 3, 90, "ACCEPT", 3));bad(root, bind(root, sha), "wrong_attempt");
		write(root, row(120, 100, 90, 0, 0, 0, 8, 1, 1, 90, "REJECT_ENTRY", 0));bad(root, bind(root, sha), "false_accept");
		write(root, row(120, 100, 90, 0, 0, 0, -1, -1, 0, 0, "UNKNOWN", 0));bad(root, bind(root, sha), "unknown_outcome");
		write(root, a.replace("\t1\t0\t7\t", "\t9223372036854775808\t0\t7\t"));bad(root, bind(root, sha), "long_overflow");
		write(root, a.replace("\t1000\t", "\t-1\t"));bad(root, bind(root, sha), "negative_count");
		write(root, a.trim()+"\textra\n");bad(root, bind(root, sha), "extra_field");
		write(root, a);NcrnaMitoCrossLibraryTest.write(root+"/time.txt", "wall_seconds=2 user_seconds=NaN system_seconds=1\n");bad(root, bind(root, sha), "invalid_time");
		System.out.println("MITO_WINDOW_REPORT_TEST_PASS cap2_prefix late_accepts zero_windows score_extrema event_histogram malformed_conservation overflow");
	}
	private static String row(int top, int second, int third, int models, int calls, int skips, int model, int rank, int attempt, int shared, String outcome, int prefix){
		return "ctg\t0\t1\t0\t3000\t1000\t"+top+'\t'+second+'\t'+third+'\t'+models+'\t'+calls+'\t'+skips+'\t'+model+'\t'+rank+'\t'+attempt+'\t'+shared+'\t'+outcome+'\t'+prefix+'\n';
	}
	private static void write(String root, String rows){NcrnaMitoCrossLibraryTest.write(root+"/windows.tsv", MitoSsuTierLog.HEADER+'\n'+rows);}
	private static String bind(String root, String sha)throws Exception{return NcrnaMitoCrossLibraryTest.bind(root, sha, "windows.tsv", "time.txt");}
	private static String read(String path){return NcrnaMitoCrossLibraryTest.read(path);}
	private static void run(String root, String pins, String name){NcrnaMitoWindowReport.main(new String[]{"in="+root+"/windows.tsv", "time="+root+"/time.txt", "verified="+pins, "outprefix="+root+'/'+name, "label=TEST"});}
	private static void bad(String root, String pins, String name){boolean failed=false;try{run(root, pins, name);}catch(IllegalArgumentException e){failed=true;}
		check(failed && !new File(root+'/'+name+".summary.tsv").exists(), "Invalid observer evidence must fail before output: "+name);}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
