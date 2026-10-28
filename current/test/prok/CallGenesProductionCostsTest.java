package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.CallGenesProductionCosts;
import prok.RrnaResourceIO;

/** Independent clock populations and signed file deltas must survive reduction.
 * @author Brian Bushnell, Ganyu
 */
public final class CallGenesProductionCostsTest {
	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply a fresh fixture directory");final Path root=Paths.get(args[0]);Files.createDirectory(root);
		final StringBuilder manifest=new StringBuilder(CallGenesProductionCosts.RUN_HEADER).append('\n');
		for(int t:new int[]{1,16}){for(boolean enabled:new boolean[]{false,true}){
			final String name="t"+t+(enabled?"_on":"_off");final Path prefix=root.resolve(name);
			final String input=CallGenesProductionCosts.INPUT_HEADER+"\n0\ta.fa\t1\t10\t1000\t"+(enabled?900000000:800000000)+"\t2000000000\t1000000000\n"+
				"1\tb.fa\t1\t20\t2000\t"+(enabled?400000000:600000000)+"\t1500000000\t900000000\n";
			write(root.resolve(name+".inputs.tsv"), "# preserved native comments\n"+input);
			final StringBuilder families=new StringBuilder(CallGenesProductionCosts.FAMILY_HEADER).append('\n');
			for(int i=0; i<2; i++){
				final String head=i+"\t"+(i==0?"a.fa":"b.fa")+'\t';
				families.append(head).append("rnasep\t1\t100\t0\t0\t2\t1\t1000000\n");
				if(enabled){
					families.append(head).append("cobalamin_rf00174\t10\t1000\t0\t0\t5\t3\t2000000\n");
					families.append(head).append("cobalamin_rf01482\t4\t400\t0\t0\t4\t2\t3000000\n");
					families.append(head).append("cobalamin_rf01689\t2\t200\t0\t0\t2\t1\t4000000\n");
				}
			}
			write(root.resolve(name+".families.tsv"), families.toString());
			write(root.resolve(name+".startup.tsv"), CallGenesProductionCosts.STARTUP_HEADER+"\n"+
				(enabled?"B12_RESOURCES\tCOMPLETE\t100000000\t40000000\n":"B12_RESOURCES\tSKIPPED_DISABLED\t0\t0\n")+
				"SHARED_SEED_INDEX\tCOMPLETE\t20000000\t10000000\n");
			write(root.resolve(name+".time"), "wall_seconds="+(t==1?(enabled?12:10):(enabled?4:3))+" user_seconds="+(enabled?9:8)+" system_seconds="+(enabled?3:2)+"\n");
			manifest.append(name).append('\t').append(t).append('\t').append(enabled).append('\t').append(prefix).append('\t').append(root.resolve(name+".time")).append('\n');
		}}
		write(root.resolve("runs.tsv"), manifest.toString());run(root, "runs.tsv", "result");
		final ByteFile in=RrnaResourceIO.open(root.resolve("result.runs.tsv").toString());final LineParser1 p=new LineParser1('\t');int rows=0;
		try{in.nextLine();for(byte[] row; (row=in.nextLine())!=null; ){
			p.set(row);rows++;final boolean enabled=p.termEquals("true", 2);
			check(p.parseInt(3)==2 && close(p.parseDouble(7), enabled?12:10) && close(p.parseDouble(8), enabled?1.3:1.4),
				"Whole-process CPU may rise while summed worker CPU falls; neither replaces the other");
			check(close(p.parseDouble(9), enabled?.05:.01) && close(p.parseDouble(11), enabled?.02:.002)
				&& close(p.parseDouble(12), enabled?.018:0) && p.parseInt(13)==(enabled?32:0) && p.parseInt(14)==(enabled?12:0),
				"Selected startup CPU and lazy worker initialization are separate; all three B12 family events sum exactly");
		}}finally{check(!in.close(), "Fixture summary read failed");}
		check(rows==4, "Exactly four arms must survive reduction");
		final ByteFile deltas=RrnaResourceIO.open(root.resolve("result.deltas.tsv").toString());int negative=0;
		try{deltas.nextLine();for(byte[] row; (row=deltas.nextLine())!=null; ){
			p.set(row);if(p.parseInt(1)==1){check(close(p.parseDouble(7), -.2), "A negative observed input delta must not be clipped");negative++;}
		}}finally{check(!deltas.close(), "Fixture delta read failed");}check(negative==2, "Both thread pairs retain the signed delta");
		// A missing diagnostic is not free initialization; preserve an explicit failure.
		final String source=root.resolve("t1_on").toString(), bad=root.resolve("missing").toString();
		Files.copy(root.resolve("t1_on.inputs.tsv"), root.resolve("missing.inputs.tsv"));
		Files.copy(root.resolve("t1_on.startup.tsv"), root.resolve("missing.startup.tsv"));
		write(root.resolve("missing.families.tsv"), CallGenesProductionCosts.FAMILY_HEADER+"\n0\ta.fa\trnasep\t1\t100\t0\t0\t2\t1\tNA\n");
		write(root.resolve("missing_runs.tsv"), manifest.toString().replace("\t"+source+"\t", "\t"+bad+"\t"));
		boolean failed=false;try{run(root, "missing_runs.tsv", "missing_result");}catch(IllegalArgumentException e){failed=e.getMessage().contains("missing measurements are not zeros");}
		check(failed && !Files.exists(root.resolve("missing_result.runs.tsv")), "Unavailable required measurements fail before publishing cost summaries");
		System.out.println("CALLGENES_PRODUCTION_COSTS_TEST_PASS");
	}
	private static void run(Path root, String manifest, String out){CallGenesProductionCosts.main(new String[]{"in="+root.resolve(manifest), "outprefix="+root.resolve(out), "expected=2"});}
	private static boolean close(double a, double b){return Math.abs(a-b)<1e-8;}
	private static void write(Path path, String text){final ByteStreamWriter w=RrnaResourceIO.writer(path.toString());try{w.print(text);}finally{RrnaResourceIO.finish(w);}}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
