package test.prok;

import java.io.File;
import java.util.Arrays;
import java.util.Locale;

import parse.Parse;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaBoundaryInvocationSummary;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Literal counts from the already-closed real-seam fixture; no caller repetition.
 * OFF retains only one fixed-end invocation; ON includes rejected candidates
 * while both fixed ends never invoke. Wrong headers/declared populations fail.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaBoundaryInvocationSummaryTest{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaBoundaryInvocationSummaryTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown invocation-summary fixture argument: "+arg);
			}
			check(parser.in1!=null && parser.in2!=null && parser.out1!=null && !parser.overwrite && !parser.append
				&& new File(parser.out1).mkdir(), "Supply in1=closed_OFF_trace in2=closed_ON_trace out=fresh_directory");
			final String off=new File(parser.out1, "off.tsv").getPath(), on=new File(parser.out1, "on.tsv").getPath();
			NcrnaBoundaryInvocationSummary.summarize(parser.in1, off, 1);
			NcrnaBoundaryInvocationSummary.summarize(parser.in2, on, 1);
			row(off, "sam1\t5\t5\t1\t1\t1\t3\t0");
			row(on, "sam1\t5\t5\t4\t3\t2\t3\t0");
			for(int version=1; version<=2; version++){
				final String old=new File(parser.out1, "version"+version+".tsv").getPath();
				maximumEncoding(parser.in1, old, version, "2147483647");
				final String summary=new File(parser.out1, "version"+version+"_summary.tsv").getPath();
				NcrnaBoundaryInvocationSummary.summarize(old, summary, 1);
				row(summary, "sam1\t5\t5\t1\t1\t1\t3\t0");
			}
			final String badMaximum=new File(parser.out1, "bad_maximum.tsv").getPath();
			maximumEncoding(parser.in1, badMaximum, 3, "4294967496");
			fails(badMaximum, new File(parser.out1, "bad_maximum_summary.tsv").getPath(), 1);
			final String wrong=new File(parser.out1, "wrong_header.tsv").getPath();
			try(Output output=new Output(wrong)){
				output.writer.println("#ncrna_boundary_trace\t1");
			}
			fails(wrong, new File(parser.out1, "bad_header_out.tsv").getPath(), 1);
			fails(parser.in1, new File(parser.out1, "bad_input_count_out.tsv").getPath(), 2);
			final String unknown=new File(parser.out1, "unknown_trace.tsv").getPath();
			mutate(parser.in1, unknown, true);
			final String unknownSummary=new File(parser.out1, "unknown_summary.tsv").getPath();
			NcrnaBoundaryInvocationSummary.summarize(unknown, unknownSummary, 1);
			row(unknownSummary, "sam1\t5\t4\t1\t1\t1\t3\t1");
			final String mismatch=new File(parser.out1, "mismatched_total.tsv").getPath();
			mutate(parser.in1, mismatch, false);
			fails(mismatch, new File(parser.out1, "bad_total_out.tsv").getPath(), 1);
			pp.outstream.println("NCRNA_BOUNDARY_INVOCATION_SUMMARY_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Independent literal old-format and malformed maximum controls. */
	private static void maximumEncoding(final String in, final String out, final int version, final String maximum){
		check(version>=1 && version<=3, "Invocation fixture version must be explicitly supported");
		final LineParser1 parser=new LineParser1('\t');
		final ByteBuilder text=new ByteBuilder();
		try(Input input=new Input(in); Output output=new Output(out)){
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				parser.set(line);
				text.clear();
				if(parser.termEquals("#ncrna_boundary_invocations", 0)){
					text.append("#ncrna_boundary_invocations\t").append(version);
				}else if(parser.termEquals("C", 0)){
					parser.setBounds(22);
					text.append(line, 0, parser.a()).append(maximum).append(line, parser.b(), line.length-parser.b());
				}else{text.append(line);}
				output.writer.print(text.nl());
			}
		}
	}
	private static void row(final String path, final String expected){
		try(Input input=new Input(path)){
			check(Arrays.equals(input.reader.nextLine(), NcrnaBoundaryInvocationSummary.HEADER.getBytes()), "Exact summary schema");
			check(Arrays.equals(input.reader.nextLine(), expected.getBytes()) && input.reader.nextLine()==null,
				"Literal candidate/returned/retained counts differ: "+expected);
		}
	}
	private static void fails(final String in, final String out, final int expected){
		boolean rejected=false;
		try{NcrnaBoundaryInvocationSummary.summarize(in, out, expected);}
		catch(IllegalArgumentException exception){rejected=true;}
		check(rejected && !new File(out).exists(), "Invalid provenance must fail before summary output is opened");
	}
	/** Changes literal fixture bytes only: the first uninvoked selected row becomes
	 * unknown, or the contig declares six candidates while the stream contains five. */
	private static void mutate(final String in, final String out, final boolean unknown){
		boolean candidate=false, total=false;
		try(Input input=new Input(in); Output output=new Output(out)){
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				if(line.length>1 && line[0]=='C' && unknown && !candidate){
					replace(line, 14, '1', '0');
					candidate=true;
				}else if(line.length>1 && line[0]=='T'){
					check(!total, "The literal input has exactly one contig total");
					if(unknown){replace(line, 9, '0', '1');}
					else{replace(line, 6, '5', '6');}
					total=true;
				}
				output.writer.println(line);
			}
		}
		check(total && (!unknown || candidate), "The requested literal mutation must have occurred");
	}
	private static void replace(final byte[] line, final int column, final char before, final char after){
		int field=0, start=0;
		for(int i=0; i<=line.length; i++){
			if(i==line.length || line[i]=='\t'){
				if(field==column){
					check(i==start+1 && line[start]==before, "The finite mutation needs its independently specified one-byte value");
					line[start]=(byte)after;
					return;
				}
				field++;
				start=i+1;
			}
		}
		throw new AssertionError("Missing literal fixture field "+column);
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
