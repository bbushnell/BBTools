package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.PrepareNcrnaRosterBaseline;
import shared.Shared;
import structures.ByteBuilder;

/** Keeps older, zero-gene and unavailable source frames distinct in the full
 * original-gene population. This is metadata correctness, not caller recall.
 * Adapted from template.A_SampleBasic2's output-only parser lifecycle.
 * @author Brian Bushnell, Raiden
 */
public final class PrepareNcrnaRosterBaselineTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, PrepareNcrnaRosterBaselineTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown fixture argument: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append,
				"Use out=fresh_fixture_directory without input or overwrite");
			new PrepareNcrnaRosterBaselineTest(new File(parser.out1).getAbsoluteFile()).run();
			pp.outstream.println("NCRNA_ROSTER_BASELINE_PREPARATION_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private PrepareNcrnaRosterBaselineTest(final File root_){
		root=root_;
		check(root.mkdir(), "A baseline preparation fixture needs a fresh directory");
	}

	private void run(){
		write("older.fa", new ByteBuilder(">same\nACGTACGT\n"));
		write("live.fa", new ByteBuilder(">same\nAAAACCCC\n"));
		write("zero.fa", new ByteBuilder(">same\nGGGGTTTT\n"));
		final ByteBuilder frames=new ByteBuilder(PrepareNcrnaRosterBaseline.FRAME_HEADER).nl();
		frame(frames, 7, "older.fa", "live.fa", 1);
		frame(frames, 7, "live.fa", "live.fa", 1);
		frame(frames, 8, "zero.fa", "zero.fa", 0);
		frame(frames, 9, "missing.fa", "missing.fa", 1);
		write("frames.tsv", frames);
		final ByteBuilder genes=new ByteBuilder(PrepareNcrnaRosterBaseline.ORIGINAL_HEADER).nl();
		gene(genes, 40, 7, "older.fa", false);
		gene(genes, 41, 7, "live.fa", true);
		gene(genes, 42, 9, "missing.fa", false);
		write("genes.tsv", genes);
		prepare("accepted", "genes.tsv");
		final LineParser1 p=new LineParser1('\t');
		int rows=0, unavailable=0, retainedTruncated=0;
		try(Input in=new Input(path("accepted.original_calls.tsv.gz"))){
			in.reader.nextLine();
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				check(p.terms()==31, "Each original row retains all29fields plus exact input ID and source state");
				if(p.termEquals("UNRESOLVED_SOURCE_UNAVAILABLE", 30)){unavailable++;}
				if(p.parseInt(17)>0){retainedTruncated++;}
				rows++;
			}
		}
		check(rows==3 && unavailable==1 && retainedTruncated==1,
			"Unavailable and construction-excluded original calls remain in the recall ledger");
		try(Input in=new Input(path("accepted.inputs.config"))){
			final String list=new String(in.reader.nextLine(), StandardCharsets.UTF_8);
			check(list.equals("in="+path("older.fa")+","+path("live.fa")+","+path("zero.fa")),
				"Use both exact older/live coordinate frames and the zero-gene frame; never replace older with its live alias");
			check(in.reader.nextLine()==null, "The input configuration contains one complete population list");
		}
		try(Input in=new Input(path("accepted.frames.tsv"))){
			in.reader.nextLine();
			rows=0;
			while(in.reader.nextLine()!=null){rows++;}
			check(rows==4, "The frame ledger includes available, zero-call and unavailable inputs");
		}
		final ByteBuilder duplicate=new ByteBuilder().append(genes.toBytes());
		gene(duplicate, 40, 7, "older.fa", false);
		write("duplicate.tsv", duplicate);
		reject("duplicate", "duplicate.tsv", "unique ID");
		check(new File(path("older.fa")).renameTo(new File(path("older.before_size_change.fa"))),
			"Preserve the initial fixture sequence before testing a changed source size");
		write("older.fa", new ByteBuilder(">same\nACGTACGTA\n"));
		reject("changed", "genes.tsv", "Source size changed");
	}

	private void frame(final ByteBuilder out, final int group, final String name, final String alias, final int genes){
		final File file=new File(path(name));
		out.append(group).tab().append(file.getPath()).tab().append(path(alias)).tab()
			.append(file.exists() ? file.length() : 15).append("\t8\t1\t1\t")
			.append(group==7 ? 2 : 1).tab().append(genes).tab().append(genes>0).tab()
			.append(group==7).tab().append(file.exists()).append("\tFIXTURE\n");
	}

	private void gene(final ByteBuilder out, final int id, final int group, final String name, final boolean truncated){
		out.append(id).tab().append(group).tab().append(path(name))
			.append("\tsame\t1\t8\t+\ttrue\tfalse\t1\t0\t0\tQUALIFIED\t")
			.append(id).tab().append(truncated ? -1 : id).tab()
			.append(truncated ? "EXCLUDED_NATIVE_TRUNCATION" : "REPRESENTATIVE")
			.append("\t1\t").append(truncated ? 1 : 0)
			.append("\t1\t8\t8\tfalse\t").append(truncated ? "NATIVE_TRUNCATION" : "ELIGIBLE_NATIVE_HIT")
			.append("\t0\t0\t0\t0\tCONTIG_END_REACH\t50.0\n");
	}

	private void prepare(final String prefix, final String genes){
		PrepareNcrnaRosterBaseline.main(new String[]{"in="+path("frames.tsv"), "originals="+path(genes),
			"out="+path(prefix), "expectedgenes=3", "expectedframes=4"});
	}

	private void reject(final String prefix, final String genes, final String reason){
		boolean rejected=false;
		try{prepare(prefix, genes);}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains(reason);}
		check(rejected, "Malformed or changed source metadata must fail loudly: "+reason);
	}

	private String path(final String name){return new File(root, name).getPath();}
	private void write(final String name, final ByteBuilder bytes){try(Output out=new Output(path(name))){out.writer.print(bytes);}}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private final File root;
}
