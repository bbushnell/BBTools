package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.CaptureNcrnaWindows;
import prok.MeasureNcrnaRejectedWindows;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.PrepareNcrnaRosterBaseline;
import prok.SummarizeNcrnaRosterBaseline;
import structures.ByteBuilder;
import shared.Shared;

/** Two real primary models: an exact80base fragment versus a lower-identity
 *120base candidate. Checks separate span-conditioned winners, original-forward
 *minus coordinates, original physical-edge flags and missing own-model state.
 * @author Brian Bushnell, Raiden
 */
public final class MeasureNcrnaRejectedWindowsTest{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, MeasureNcrnaRejectedWindowsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, split[0].toLowerCase(java.util.Locale.ROOT), split.length>1 ? split[1] : null), "Unknown fixture parameter: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append, "Use out=fresh_fixture_directory");
			new MeasureNcrnaRejectedWindowsTest(new File(parser.out1)).run();
			pp.outstream.println("NCRNA_REJECTED_WINDOW_MEASUREMENT_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private MeasureNcrnaRejectedWindowsTest(final File root_){
		root=root_.getAbsoluteFile();
		check(root.mkdir(), "The integration fixture must preserve prior attempts by requiring a fresh directory");
	}
	private void run(){
		final byte[] query=Prok16sPartialRescueTest.random(120, 20261115), variant=query.clone();
		for(int i=5; i<variant.length; i+=10){variant[i]=variant[i]=='A' ? (byte)'C' : (byte)'A';}
		final byte[] window=Prok16sPartialRescueTest.random(320, 20261116);
		System.arraycopy(query, 0, window, 100, query.length);
		write("models.fa", new ByteBuilder().append('>').append(A).nl().append(Arrays.copyOf(query, 80)).nl()
			.append('>').append(B).nl().append(variant).nl());
		write("capture.windows.fa.gz", new ByteBuilder(">i0e0\n").append(window).append("\n>i0e1\n").append(window).nl());
		write("capture.windows.tsv", new ByteBuilder(CaptureNcrnaWindows.WINDOW_HEADER).nl()
			.append(("i0e0\t0\tframe_0\t/unused.fa\t0\t"+CONTIG+"%09description\t0\t1000\t1\t100\t419\t100\t419\t0\t320\n").getBytes(StandardCharsets.UTF_8))
			.append(("i0e1\t0\tframe_0\t/unused.fa\t0\t"+CONTIG+" description\t1\t1000\t1\t0\t319\t680\t999\t0\t320\n").getBytes(StandardCharsets.UTF_8)));
		final ByteBuilder population=new ByteBuilder(PrepareNcrnaRosterBaseline.ORIGINAL_HEADER).append("\tbaseline_input_id\tbaseline_input_state\n");
		original(population, 10, 201, 320, false, false);
		original(population, 11, 781, 900, true, true);
		write("population.original_calls.tsv.gz", population);
		final ByteBuilder genes=new ByteBuilder(SummarizeNcrnaRosterBaseline.GENES_HEADER).nl();
		grade(genes, 10, 201, 320, false, A);
		grade(genes, 11, 781, 900, true, "NA");
		write("grade.genes.tsv.gz", genes);
		write("grade.gene_windows.tsv.gz", new ByteBuilder(SummarizeNcrnaRosterBaseline.ASSOCIATION_HEADER).nl()
			.append("10\tframe_0\ti0s0\tshard_0:i0e0\ttrue\ttrue\tfalse\t2\ttrue\n")
			.append("11\tframe_0\ti0s1\tshard_0:i0e1\ttrue\ttrue\tfalse\t2\ttrue\n"));
		write("misses.tsv", new ByteBuilder("library\toriginal_gene\tnative_truncated\tincomplete_coverage\toutcome\tcaptured_overlap\n")
			.append("fixture\t10\tfalse\tfalse\tCAPTURED_NO_GENE_ACCEPTANCE\ttrue\n")
			.append("fixture\t11\tfalse\ttrue\tCAPTURED_NO_GENE_ACCEPTANCE\ttrue\n"));
		MeasureNcrnaRejectedWindows.main(new String[]{"library=fixture", "captureid=shard_0", "capture="+path("capture"),
			"population="+path("population"), "grade="+path("grade"), "in="+path("misses.tsv"), "consensus="+path("models.fa"),
			"out="+path("result"), "expectedmodels=2", "expectedclass2=2", "expectedgenes=2", "expectedincomplete=1",
			"expectedtruncated=0", "expectedshortboth=1", "t=2"});
		final LineParser1 p=new LineParser1('\t');
		boolean shortWinner=false, longWinner=false, reverse=false, absentOwn=false;
		try(Input in=new Input(path("result.best.tsv.gz"))){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(!p.termEquals("PRIMARY_RAW", 5)){continue;}
				if(p.parseInt(2)==10 && p.termEquals("ID064_NO_SPAN", 6)){
					shortWinner=p.termEquals(A, 8) && p.parseInt(9)==80;
				}
				if(p.parseInt(2)==10 && p.termEquals("ID064_SPAN109", 6)){
					longWinner=p.termEquals(B, 8) && p.parseInt(9)>=109 && p.parseFloat(7)<1;
				}
				if(p.parseInt(2)==11 && p.termEquals("BEST_POSITIVE_OVERLAP", 6)){
					reverse=p.parseInt(10)==820 && p.parseInt(11)==899 && p.termEquals("true", 3) && p.termEquals("true", 4);
				}
				if(p.parseInt(2)==11 && p.termEquals("OWN_MODEL_ANY_OVERLAP_STATE", 6)){
					absentOwn=p.termEquals("NA", 7) && p.termEquals("NA", 17);
				}
			}
		}
		check(shortWinner && longWinner, "A higher-identity short fragment must not hide the lower-identity candidate passing the span gate");
		check(reverse && absentOwn, "Minus coordinates, native subsets and a missing construction-model binding must stay explicit");
		boolean work=false;
		try(Input in=new Input(path("result.summary.tsv"))){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(p.termEquals("physical_edge_model_alignments", 0)){work=p.parseInt(1)==2;}
			}
		}
		check(work, "Only the original physical-edge window may request clipped kernels; interior saved-window ends are not contig ends");
	}
	private static void original(ByteBuilder out, int id, int start, int stop, boolean reverse, boolean incomplete){
		out.append(id).append("\t7\t/unused.fa\t").append(CONTIG.getBytes(StandardCharsets.UTF_8)).tab().append(start).tab().append(stop).tab().append(reverse ? '-' : '+')
			.append("\ttrue\tfalse\t1\t0\t0\tQUALIFIED\t").append(id).tab().append(id).append("\tCONSTRUCTION_REPRESENTATIVE\t1\t0\t")
			.append(incomplete ? "45\t80\t100\ttrue" : "1\t100\t100\tfalse").append("\tFIXTURE\t0\t0\t0\t0\tFIXTURE\t50\tframe_0\tAVAILABLE_RECORDED_FRAME\n");
	}
	private static void grade(ByteBuilder out, int id, int start, int stop, boolean reverse, String own){
		out.append(id).append("\tframe_0\t").append(CONTIG.getBytes(StandardCharsets.UTF_8)).tab().append(start).tab().append(stop).tab().append(reverse ? '-' : '+').append("\t0\t").append(id).tab()
			.append(own).append("\t1\t1\tWINDOW_WITHOUT_ACCEPTED_MODEL\t-1\t-1\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tfalse\tNA\t2\tNA\tCONSTRUCTION_REPRESENTATIVE\n");
	}
	private void write(String name, ByteBuilder bytes){try(Output out=new Output(path(name))){out.writer.print(bytes);}}
	private String path(String name){return new File(root, name).getPath();}
	private static void check(boolean condition, String why){if(!condition){throw new AssertionError(why);}}
	private final File root;
	private static final String A="cobalamin_RF00174_A", B="cobalamin_RF00174_B";
	private static final String CONTIG="ctgΩ";
}
