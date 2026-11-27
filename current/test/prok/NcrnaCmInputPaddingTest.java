package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaCmInputPadding;
import prok.WriteNcrnaCmInputPadding;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Literal input-padding resource contract, exact axes and sparse fallback.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCmInputPaddingTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaCmInputPaddingTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, split[0].toLowerCase(Locale.ROOT), split.length>1 ? split[1] : null),
					"Unknown input-padding fixture option");
			}
			check(parser.out1!=null, "Bind a fresh data-only fixture directory with out=");
			run(new File(parser.out1));
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fixture           ----------------*/
	/*--------------------------------------------------------------*/

	private static void run(final File root){
		check(root.mkdir(), "Input-padding fixtures require a fresh directory");
		final String[] names={"a n=10", "b", "new"};
		final byte[][] bases={bytes("aCgTn"), bytes("GGA"), bytes("TTT")};
		final String literal="#BBTools CM input padding v1\n#policy\tP98_PLUS10\t5\t2\n"
			+"model_ordinal0\tmodel\tconsensus_sequence\tpad5\tpad3\ttraining_finite\ttraining_groups\tlearned\n"
			+"0\ta\tACGTN\t3\t9\t10\t2\ttrue\n"
			+"1\tb\tGGA\t20\t20\t4\t2\tfalse\n"
			+"2\tnew\tTTT\t20\t20\t0\t0\tfalse\n";
		final File source=write(root, "literal.tsv", literal);
		final NcrnaCmInputPadding.Data data=NcrnaCmInputPadding.read(source.getPath(), names, bases);
		check(!data.fixed && data.size()==3 && data.rule.equals("P98_PLUS10") && data.minCalls==5 && data.minGroups==2,
			"The literal policy must retain its declared rule, support minima and complete three-model axis");
		check(data.pad5(0)==3 && data.pad3(0)==9 && data.finiteCalls(0)==10 && data.sourceGroups(0)==2,
			"The supported model retains distinct sense-oriented end pads and finite/source counts");
		check(data.pad5(1)==20 && data.pad3(1)==20 && data.finiteCalls(1)==4
			&& data.pad5(2)==20 && data.pad3(2)==20 && data.finiteCalls(2)==0 && data.sourceGroups(2)==0,
			"Insufficient and absent training support both retain20/20 without inventing observed errors");
		final File exported=new File(root, "exported.tsv");
		NcrnaCmInputPadding.write(exported.getPath(), "P98_PLUS10", 5, 2, names, bases,
			new int[]{3, 20, 20}, new int[]{9, 20, 20}, new int[]{10, 4, 0}, new int[]{2, 2, 0});
		compare(source, exported);
		reject(root, "changed_sequence.tsv", literal.replace("ACGTN", "ACGTA"), names, bases);
		reject(root, "reordered.tsv", literal.replace("0\ta\t", "1\ta\t"), names, bases);
		reject(root, "renamed.tsv", literal.replace("0\ta\t", "0\tother\t"), names, bases);
		reject(root, "truncated.tsv", literal.substring(0, literal.indexOf("2\tnew\t")), names, bases);
		reject(root, "fallback_end.tsv", literal.replace("GGA\t20\t20", "GGA\t19\t20"), names, bases);
		reject(root, "false_support.tsv", literal.replace("4\t2\tfalse", "4\t2\ttrue"), names, bases);
		reject(root, "invalid_groups.tsv", literal.replace("10\t2\ttrue", "10\t0\ttrue"), names, bases);
		reject(root, "invalid_minimum.tsv", literal.replace("P98_PLUS10\t5\t2", "P98_PLUS10\t0\t2"), names, bases);
		reject(root, "unknown_rule.tsv", literal.replace("P98_PLUS10", "UNKNOWN"), names, bases);
		reject(root, "learned_fixed_rule.tsv", literal.replace("P98_PLUS10", "FIXED"), names, bases);
		reject(root, "learned_overflow.tsv", literal.replace("3\t9\t10", "4294967296\t9\t10"), names, bases);
		reject(root, "repeated.tsv", literal+"2\tnew\tTTT\t20\t20\t0\t0\tfalse\n", names, bases);
		boolean failed=false;
		final File badExport=new File(root, "bad_export.tsv");
		try{
			NcrnaCmInputPadding.write(badExport.getPath(), "P98_PLUS10", 5, 2, names, bases,
				new int[]{3, 19, 20}, new int[]{9, 20, 20}, new int[]{10, 4, 0}, new int[]{2, 2, 0});
		}catch(final IllegalArgumentException expected){failed=true;}
		check(failed && !badExport.exists(), "Invalid fallback export must fail before creating a partial resource");
		fixed(root, names, bases);
		selected(root, names, bases, data);
		System.err.println("NCRNA_CM_INPUT_PADDING_RESOURCE_FIXTURE_PASS");
	}

	/** Thin percentile fits keep the selected base, rather than inheriting old20/20 fallback. */
	private static void selected(final File root, final String[] names, final byte[][] bases,
			final NcrnaCmInputPadding.Data percentile){
		final String literal="#BBTools CM input padding selected v1\n#selected\tP98_PLUS10\t5\t2\t25\t10\n"
			+"model_ordinal0\tmodel\tconsensus_sequence\tpad5\tpad3\ttraining_finite\ttraining_groups\tselected5\tselected3\n"
			+"0\ta\tACGTN\t3\t10\t10\t2\ttrue\tfalse\n"
			+"1\tb\tGGA\t25\t10\t4\t2\tfalse\tfalse\n"
			+"2\tnew\tTTT\t25\t10\t0\t0\tfalse\tfalse\n";
		final File source=write(root, "selected_literal.tsv", literal);
		final NcrnaCmInputPadding.Data data=NcrnaCmInputPadding.read(source.getPath(), names, bases);
		check(data.selected && !data.fixed && data.base5==25 && data.base3==10 && data.rule.equals("P98_PLUS10"),
			"A selected mask declares its own base and source rule without pretending to be a newly fitted policy");
		check(data.pad5(0)==3 && data.pad3(0)==10 && data.pad5(1)==25 && data.pad3(1)==10
			&& data.pad5(2)==25 && data.pad3(2)==10 && data.finiteCalls(1)==4 && data.finiteCalls(2)==0,
			"End masks are independent; thin and absent fits keep25/10 while original support remains literal");
		final File exported=new File(root, "selected_exported.tsv");
		NcrnaCmInputPadding.writeSelected(exported.getPath(), 25, 10, percentile,
			new boolean[]{true, false, false}, new boolean[]{false, false, false});
		compare(source, exported);
		final File symmetric=new File(root, "selected_25_25.tsv");
		NcrnaCmInputPadding.writeSelected(symmetric.getPath(), 25, 25, percentile,
			new boolean[]{false, false, false}, new boolean[]{true, false, false});
		final NcrnaCmInputPadding.Data other=NcrnaCmInputPadding.read(symmetric.getPath(), names, bases);
		check(other.pad5(0)==25 && other.pad3(0)==9 && other.pad3(1)==25 && other.pad3(2)==25,
			"The second declared base retains25at unselected three-prime ends, independent of the25/10 mask");
		reject(root, "selected_changed_base.tsv", literal.replace("#selected\tP98_PLUS10\t5\t2\t25\t10", "#selected\tP98_PLUS10\t5\t2\t20\t10"), names, bases);
		reject(root, "selected_bad_mask.tsv", literal.replace("true\tfalse", "maybe\tfalse"), names, bases);
		reject(root, "selected_unsupported.tsv", literal.replace("25\t10\t4\t2\tfalse\tfalse", "20\t10\t4\t2\ttrue\tfalse"), names, bases);
		reject(root, "selected_tied.tsv", literal.replace("ACGTN\t3\t10", "ACGTN\t25\t10"), names, bases);
		reject(root, "selected_wrong_axis.tsv", literal.replace("ACGTN", "ACGTA"), names, bases);
		reject(root, "selected_truncated.tsv", literal.substring(0, literal.indexOf("2\tnew\t")), names, bases);
		final File invalid=new File(root, "selected_bad_export.tsv");
		boolean failed=false;
		try{
			NcrnaCmInputPadding.writeSelected(invalid.getPath(), 25, 10, percentile,
				new boolean[]{false, true, false}, new boolean[]{false, false, false});
		}catch(final IllegalArgumentException expected){failed=true;}
		check(failed && !invalid.exists(), "Unsupported percentile fallback cannot be selected or create a partial resource");
	}

	/** Literal fixed records have no support fields, including the previously unsupported model. */
	private static void fixed(final File root, final String[] names, final byte[][] bases){
		final String literal="#BBTools CM input padding fixed v1\n#fixed\t25\t10\n"
			+"model_ordinal0\tmodel\tconsensus_sequence\n0\ta\tACGTN\n1\tb\tGGA\n2\tnew\tTTT\n";
		final File source=write(root, "fixed_literal.tsv", literal);
		final NcrnaCmInputPadding.Data data=NcrnaCmInputPadding.read(source.getPath(), names, bases);
		check(data.fixed && data.rule.equals("FIXED") && data.size()==3 && data.minCalls==0 && data.minGroups==0,
			"Explicit fixed policy retains the whole three-model axis and declares no learned minima");
		for(int m=0; m<names.length; m++){
			check(data.pad5(m)==25 && data.pad3(m)==10 && data.finiteCalls(m)==0 && data.sourceGroups(m)==0,
				"Fixed constants apply even to unsupported models without fabricating fitted observations");
		}
		final File exported=new File(root, "fixed_exported.tsv");
		NcrnaCmInputPadding.writeFixed(exported.getPath(), 25, 10, names, bases);
		compare(source, exported);
		reject(root, "fixed_negative5.tsv", literal.replace("#fixed\t25\t10", "#fixed\t-1\t10"), names, bases);
		reject(root, "fixed_negative3.tsv", literal.replace("#fixed\t25\t10", "#fixed\t25\t-1"), names, bases);
		reject(root, "fixed_overflow.tsv", literal.replace("#fixed\t25\t10", "#fixed\t4294967296\t10"), names, bases);
		reject(root, "fixed_empty.tsv", literal.replace("#fixed\t25\t10", "#fixed\t\t10"), names, bases);
		reject(root, "fixed_nondigit.tsv", literal.replace("#fixed\t25\t10", "#fixed\tx\t10"), names, bases);
		reject(root, "fixed_support.tsv", literal.replace("#fixed\t25\t10", "#fixed\t25\t10\t5"), names, bases);
		reject(root, "fixed_policy.tsv", literal.replace("#fixed", "#policy"), names, bases);
		reject(root, "fixed_sequence.tsv", literal.replace("ACGTN", "ACGTA"), names, bases);
		reject(root, "fixed_name.tsv", literal.replace("0\ta\t", "0\tother\t"), names, bases);
		reject(root, "fixed_order.tsv", literal.replace("0\ta\t", "1\ta\t"), names, bases);
		reject(root, "fixed_truncated.tsv", literal.substring(0, literal.indexOf("2\tnew\t")), names, bases);
		reject(root, "fixed_repeated.tsv", literal+"2\tnew\tTTT\n", names, bases);
		final File fasta=write(root, "fixed.fa", ">a n=10\naCgTn\n>b\nGGA\n>new\nTTT\n");
		final File directory=new File(root, "fixed_cli");
		WriteNcrnaCmInputPadding.main(new String[]{"in="+fasta, "out="+directory,
			"pad5=25,20,0", "pad3=10,20,0", "t=1", "ow=f"});
		compare(source, new File(directory, "padding_25_10.tsv"));
		for(final int pad:new int[]{0, 20}){
			final NcrnaCmInputPadding.Data constant=NcrnaCmInputPadding.read(
				new File(directory, "padding_"+pad+"_"+pad+".tsv").getPath(), names, bases);
			check(constant.fixed && constant.pad5(2)==pad && constant.pad3(2)==pad,
				"Zero and matched20 remain explicit fixed policies for unsupported models");
		}
		final String[][] bad={{"pad5=20,25", "pad3=10"}, {"pad5=-1", "pad3=10"},
			{"pad5=20,", "pad3=10,20"}, {"pad5=25,25", "pad3=10,10"},
			{"pad5=25", "pad3=10", "reads=1"}, {"pad5=25", "pad3=10", "ow=t"},
			{"pad5=4294967296", "pad3=10"}};
		for(int i=0; i<bad.length; i++){
			final File rejected=new File(root, "fixed_cli_rejected"+i);
			final String[] args=new String[bad[i].length+2];
			args[0]="in="+fasta;
			args[1]="out="+rejected;
			System.arraycopy(bad[i], 0, args, 2, bad[i].length);
			boolean failed=false;
			try{WriteNcrnaCmInputPadding.main(args);}catch(final IllegalArgumentException expected){failed=true;}
			check(failed && !rejected.exists(), "Invalid fixed export options must fail before creating resources; case="+i);
		}
		System.err.println("NCRNA_CM_FIXED_PADDING_RESOURCE_AND_EXPORT_FIXTURE_PASS");
	}

	private static void reject(final File root, final String name, final String text,
			final String[] names, final byte[][] bases){
		final File file=write(root, name, text);
		boolean failed=false;
		try{NcrnaCmInputPadding.read(file.getPath(), names, bases);}catch(final IllegalArgumentException expected){failed=true;}
		check(failed, "Malformed or misbound input padding must fail loudly: "+name);
	}

	private static void compare(final File expected, final File actual){
		try(Input a=new Input(expected.getPath()); Input b=new Input(actual.getPath())){
			for(byte[] line; (line=a.reader.nextLine())!=null;){
				check(Arrays.equals(line, b.reader.nextLine()), "Resource export must match the independently written literal contract");
			}
			check(b.reader.nextLine()==null, "Resource export must not append extra descriptors");
		}
	}

	private static File write(final File root, final String name, final String text){
		final File file=new File(root, name);
		check(!file.exists(), "Fixture writes must preserve earlier evidence: "+file);
		try(Output output=new Output(file.getPath())){output.writer.print(text);}
		return file;
	}

	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}

	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
