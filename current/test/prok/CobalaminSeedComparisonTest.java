package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.CobalaminSeedComparison;
import prok.RrnaResourceIO;

/** Workload ratios use actual window/genome denominators; absent totals fail.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminSeedComparisonTest {
	public static void main(String[] args)throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Supply a fresh fixture directory");}
		final Path root=Paths.get(args[0]);Files.createDirectory(root);
		final String header="view\tscope\tgenomes\ttruth\tcalls\tTP\tFP\tFN\tunscored\texact_duplicates\tprecision\trecall\tapparent_FP_per_genome\texact5\texact3\tboth_exact\n";
		final StringBuilder grades=new StringBuilder(header);
		for(String view:new String[]{"SOURCE", "FULL_CM"}){
			for(String scope:new String[]{"RF00174", "RF01482", "RF01689", "union"}){
				final boolean positive=scope.equals("RF00174") || scope.equals("union");
				grades.append(view).append('\t').append(scope).append("\t300\t").append(positive?"2\t2\t1\t1\t1":"0\t0\t0\t0\t0")
					.append("\t0\t0\t").append(positive?"0.5\t0.5\t0.00333333":"0\t0\t0").append("\t0\t0\t0\n");
			}
		}
		write(root.resolve("grade.tsv"), grades.toString());
		final String counters="ncRNA Windows (cobalamin_rf00174):\t 10\nncRNA Alignments (cobalamin_rf00174):\t 40\n"+
			"ncRNA Windows (cobalamin_rf01482):\t 5\nncRNA Alignments (cobalamin_rf01482):\t 10\n"+
			"ncRNA Windows (cobalamin_rf01689):\t 0\nncRNA Alignments (cobalamin_rf01689):\t 0\n";
		write(root.resolve("caller.stderr"), counters);write(root.resolve("caller.time"), "wall_seconds=12.50 user_seconds=20.25 system_seconds=0.75\n");
		manifest(root, "good", "caller.stderr");run(root, "good");
		final ByteFile in=RrnaResourceIO.open(root.resolve("good.tsv").toString());final LineParser1 p=new LineParser1('\t');int rows=0;
		try{in.nextLine();for(byte[] row; (row=in.nextLine())!=null; ){
			p.set(row);rows++;check(p.terms()==28, "Comparison schema must preserve every native grade field and workload unit");
			if(p.termEquals("union", 7)){
				check(p.parseLong(22)==15 && p.parseLong(23)==50 && Math.abs(p.parseDouble(24)-50d/15)<0.000001
					&& Math.abs(p.parseDouble(25)-50d/300)<0.000001, "Union sums workload once across subtypes, then divides by windows/genomes, not calls or TP");
			}else if(p.termEquals("RF01689", 7)){check(p.termEquals("NA", 24), "A zero-window ratio is unavailable, not fabricated zero");}
			check(p.parseDouble(26)==12.5 && p.parseDouble(27)==21, "Reported CPU includes user and system while wall stays separate");
		}}finally{check(!in.close(), "Fixture output read failed");}
		check(rows==8, "Both views and all scopes must survive collection");
		write(root.resolve("missing.stderr"), counters.replace("ncRNA Alignments (cobalamin_rf01482):\t 10\n", ""));manifest(root, "missing", "missing.stderr");
		boolean failed=false;try{run(root, "missing");}catch(IllegalArgumentException expected){failed=true;}
		check(failed && !Files.exists(root.resolve("missing.tsv")), "Absent counters must fail before any output, never become zero");
		write(root.resolve("duplicate.stderr"), counters+"ncRNA Windows (cobalamin_rf00174):\t 10\n");manifest(root, "duplicate", "duplicate.stderr");
		failed=false;try{run(root, "duplicate");}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Repeated caller totals must fail rather than inflate workload");
		System.out.println("COBALAMIN_SEED_COMPARISON_TEST_PASS");
	}
	private static void manifest(Path root, String name, String counters){
		write(root.resolve(name+".runs.tsv"), "run\tseed_minimum\tseed_keys\tmap_k\ttop_n\tid_pass\tgrade\tstderr\ttime\n"+
			"one_seed\t1\t1337\t9\t100\t0.64\t"+root.resolve("grade.tsv")+'\t'+root.resolve(counters)+'\t'+root.resolve("caller.time")+'\n');
	}
	private static void run(Path root, String out){CobalaminSeedComparison.main(new String[]{"in="+root.resolve(out+".runs.tsv"), "out="+root.resolve(out+".tsv"), "expected=1"});}
	private static void write(Path path, String value){
		final ByteStreamWriter out=RrnaResourceIO.writer(path.toString());try{out.print(value);}finally{RrnaResourceIO.finish(out);}
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
