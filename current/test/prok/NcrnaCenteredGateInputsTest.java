package test.prok;

import java.io.File;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.EvaluateNcrnaRealEndpoints;
import prok.NcrnaCenteredGateFeatures;
import prok.NcrnaCenteredGateInputs;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Literal finite population and leakage controls; no networks or resources.
 * @author Brian, Ganyu
 */
public final class NcrnaCenteredGateInputsTest{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaCenteredGateInputsTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, pair[0].toLowerCase(Locale.ROOT), pair.length>1 ? pair[1] : null), "Unknown finite gate-input parameter");
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null, "Finite controls require only out=fresh_directory");
			final File root=new File(parser.out1);
			check(root.mkdir(), "Fixture root must be fresh");
			run(root);
			pp.outstream.println("NCRNA_CENTERED_GATE_INPUTS_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final File root){
		final ByteBuilder manifest=new ByteBuilder(NcrnaCenteredGateInputs.MANIFEST_HEADER).nl();
		for(int fold=0; fold<5; fold++){
			final boolean legacy=fold==0;
			final String prefix=path(root, "fold"+fold);
			final ByteBuilder audit=new ByteBuilder("#NCRNA_CENTERED_REAL_GATE_V1\n#schema\t").append(NcrnaCenteredGateFeatures.SCHEMA)
				.append("\n#panel\tdata").append(fold).append("\n#member\tRF00174\n#fold\t").append(fold)
				.append("\n#geometry\t/fixture/cache\n#folds\t/fixture/folds\n#offsets\t/fixture/offsets\n#network\t/fixture/net\n#seed\t58\n#end\t5prime\n")
				.append(legacy ? EvaluateNcrnaRealEndpoints.LEGACY_GATE_HEADER : EvaluateNcrnaRealEndpoints.GATE_HEADER).nl();
			final ByteBuilder matrix=new ByteBuilder(NcrnaCenteredGateInputs.MATRIX_HEADER).nl();
			for(int trial=0; trial<2; trial++){
				audit.append(trial+1).append("\tinput").append(fold).append("\tassembly").append(fold).append("\tref").append(trial)
					.append("\tRF00174\tmodel\tcontig\t+\t300\t0\t100\t220\t").append(trial==0 ? 103 : 100)
					.append("\t220\t109\t").append(legacy ? "2147483647" : "").append("\t102\t220\t10\tVALID\t5prime\tPROPOSAL\t").append(trial+1)
					.append("\t102\t0\t2\t").append(trial==0 ? 1 : 0).append("\tVALID\n");
				for(int i=0; i<35; i++){matrix.append(i==32 || i==33 ? 2 : 0).tab();}
				matrix.append(trial==0 ? 1 : 0).nl();
			}
			audit.append("3\tinput").append(fold).append("\tassembly").append(fold)
				.append("\tfixed\tRF00174\tmodel\tcontig2\t+\t300\t1\t0\t119\t0\t119\t109\t")
				.append(legacy ? "2147483647" : "").append("\t0\t119\t10\tVALID\t5prime\tFIXED\t0\t0\t0\t0\t0\tVALID\n#COMPLETE\t3\t2\n");
			write(prefix+".rows.tsv.gz", audit);
			write(prefix+".vectors.tsv.gz", matrix);
			manifest.append(fold).append("\tdata").append(fold).tab().append(prefix).append("\t3\t2\n");
		}
		write(path(root, "manifest.tsv"), manifest);
		new NcrnaCenteredGateInputs(arguments(root, "manifest.tsv", "result")).process();
		checkMatrix(path(root, "result/train.tsv.gz"), 8);
		checkMatrix(path(root, "result/monitor.tsv.gz"), 2);
		try(Input in=new Input(path(root, "result/summary.tsv"))){
			in.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			for(int fold=0; fold<5; fold++){
				p.set(in.reader.nextLine());
				check(p.terms()==11 && p.parseInt(2)==fold && p.parseInt(4)==3 && p.parseInt(5)==2 && p.parseInt(6)==1
					&& p.parseInt(7)==1 && p.parseInt(8)==35, "All-call, proposal, positive and hold counts remain separate for every fold");
			}
			check(in.reader.nextLine()==null, "Exactly five explicit fold summaries are required");
		}
		write(path(root, "outer.tsv"), new ByteBuilder(manifest.toString().replace("4\tdata4\t", "-1\tdata4\t")));
		boolean failed=false;
		try{new NcrnaCenteredGateInputs(arguments(root, "outer.tsv", "bad_outer")).process();}
		catch(IllegalArgumentException expected){failed=expected.getMessage().contains("outer data are forbidden");}
		check(failed && !new File(path(root, "bad_outer/summary.tsv")).exists(), "Outer data cannot acquire an accepted training result, even with valid vector shapes");
		write(path(root, "duplicate.tsv"), new ByteBuilder(manifest.toString()).append("0\tdata0\t").append(path(root, "fold0")).append("\t3\t2\n"));
		failed=false;
		try{new NcrnaCenteredGateInputs(arguments(root, "duplicate.tsv", "bad_duplicate")).process();}
		catch(IllegalArgumentException expected){failed=expected.getMessage().contains("once only");}
		check(failed && !new File(path(root, "bad_duplicate/summary.tsv")).exists(), "Repeated fold datasets cannot silently multiply teaching weight");
		final String changed=path(root, "cross_fold");
		final ByteBuilder crossRows=new ByteBuilder();
		try(Input in=new Input(path(root, "fold1.rows.tsv.gz"))){
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				crossRows.append(new String(row, java.nio.charset.StandardCharsets.UTF_8)
					.replace("\tassembly1\t", "\tassembly0\t").replace("\tref0\t", "\tdifferent0\t")
					.replace("\tref1\t", "\tdifferent1\t").replace("\tfixed\t", "\tdifferentFixed\t")).nl();
			}
		}
		write(changed+".rows.tsv.gz", crossRows);
		final ByteBuilder crossMatrix=new ByteBuilder();
		try(Input in=new Input(path(root, "fold1.vectors.tsv.gz"))){
			for(byte[] row; (row=in.reader.nextLine())!=null;){crossMatrix.append(row).nl();}
		}
		write(changed+".vectors.tsv.gz", crossMatrix);
		write(path(root, "cross_fold.tsv"), new ByteBuilder(manifest.toString().replace(path(root, "fold1"), changed)));
		failed=false;
		try{new NcrnaCenteredGateInputs(arguments(root, "cross_fold.tsv", "bad_cross_fold")).process();}
		catch(IllegalArgumentException expected){failed=expected.getMessage().contains("One assembly cannot enter multiple held folds");}
		check(failed && !new File(path(root, "bad_cross_fold/summary.tsv")).exists(),
			"Disjoint reference IDs from one assembly must not leak into a different gate fold");
	}

	private static void checkMatrix(final String file, final int expected){
		try(Input in=new Input(file)){
			in.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			for(int row=0; row<expected; row++){
				p.set(in.reader.nextLine());
				check(p.terms()==36 && p.parseInt(35)==(row%2==0 ? 1 : 0)
					&& p.parseFloat(32)==2 && p.parseFloat(33)==2 && p.parseFloat(34)==0,
					"Unchanged whole-move targets/inputs preserve fold order; monitor rows cannot enter teaching");
			}
			check(in.reader.nextLine()==null, "Merged matrix must contain exactly the accepted proposal count");
		}
	}
	private static String[] arguments(final File root, final String manifest, final String result){
		return new String[]{"in="+path(root, manifest), "out="+path(root, result), "member=RF00174", "end=5prime", "t=1", "unpigz=f", "pigz=f"};
	}
	private static String path(final File root, final String name){return new File(root, name).getPath();}
	private static void write(final String path, final ByteBuilder rows){try(Output out=new Output(path)){out.writer.print(rows);}}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
