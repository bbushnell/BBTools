package test.prok;

import java.io.File;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.SummarizeNcrnaRejectedWindows;
import shared.Shared;
import structures.ByteBuilder;

/** Reuses the previously passed two-gene raw measurement; no new alignment.
 * Checks global curves/own-model outcomes and rejects duplicate or incomplete
 * shard collections instead of treating their smaller output as a denominator.
 * @author Brian Bushnell, Raiden
 */
public final class SummarizeNcrnaRejectedWindowsTest{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SummarizeNcrnaRejectedWindowsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String measurement=null, population=null, misses=null;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(a.equals("measurement")){measurement=b;}
				else if(a.equals("population")){population=b;}
				else if(a.equals("misses")){misses=b;}
				else{check(parser.parse(arg,a,b), "Unknown fixture argument: "+arg);}
			}
			check(parser.out1!=null && measurement!=null && population!=null && misses!=null && !parser.overwrite && !parser.append,
				"Supply out=fresh measurement=passed_two_gene_prefix population=prefix misses=file without repeating the measurement");
			new SummarizeNcrnaRejectedWindowsTest(new File(parser.out1), measurement, population, misses).run();
			pp.outstream.println("NCRNA_REJECTED_WINDOW_SUMMARY_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private SummarizeNcrnaRejectedWindowsTest(File root_, String measurement_, String population_, String misses_){
		root=root_.getAbsoluteFile();
		measurement=measurement_;
		population=population_;
		misses=misses_;
		check(root.mkdir(), "The summary fixture requires a fresh directory preserving prior attempts");
	}
	private void run(){
		final String header="library\tcapture\tresult_prefix\n", row="fixture\tshard_0\t"+measurement+"\n";
		write("complete.tsv", new ByteBuilder(header).append(row));
		summarize("complete", 1);
		final LineParser1 p=new LineParser1('\t');
		boolean curve=false, own=false;
		try(Input in=new Input(path("complete.curves.tsv"))){
			in.reader.nextLine();
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				if(p.termEquals("CLASS2_ALL",1) && p.termEquals("PRIMARY_RAW",2) && p.termEquals("0.64",4) && p.parseInt(5)==109){
					curve=p.parseInt(3)==2 && p.parseInt(6)==2;
				}
			}
		}
		try(Input in=new Input(path("complete.own_models.tsv"))){
			in.reader.nextLine();
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				if(p.termEquals("CLASS2_ALL",1) && p.termEquals("PRIMARY_RAW",2)){
					own=p.parseInt(3)==2 && p.parseInt(4)==1 && p.parseInt(5)==0 && p.parseInt(6)==1;
				}
			}
		}
		check(curve && own, "The full summary must preserve a lower-identity span-eligible winner and distinguish missing from ineligible construction models");
		write("duplicate.tsv", new ByteBuilder(header).append(row).append(row));
		reject("duplicate", 2, "namespaces must be unique");
		try(Input in=new Input(measurement+".best.tsv.gz"); Output out=new Output(path("partial_source.best.tsv.gz"))){
			final ByteBuilder bytes=new ByteBuilder();
			out.writer.print(bytes.append(in.reader.nextLine()).nl());
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				if(p.parseInt(2)==10){out.writer.print(bytes.clear().append(line).nl());}
			}
		}
		try(Input in=new Input(measurement+".summary.tsv"); Output out=new Output(path("partial_source.summary.tsv"))){
			final ByteBuilder bytes=new ByteBuilder();
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){out.writer.print(bytes.clear().append(line).nl());}
		}
		write("partial.tsv", new ByteBuilder(header).append("fixture\tshard_0\t").append(path("partial_source")).nl());
		reject("partial", 1, "separately counted shard summary");
		malformedSummary("missing_counter", "selected_saved_windows", null, false);
		reject("missing_counter", 1, "counter must occur exactly once");
		malformedSummary("duplicate_counter", "primary_model_alignments", null, true);
		reject("duplicate_counter", 1, "cannot repeat");
		malformedSummary("negative_counter", "primary_model_alignments", "-1", false);
		reject("negative_counter", 1, "cannot be negative");
		malformedSummary("wrong_population", "full_library_class2_genes", "1", false);
		reject("wrong_population", 1, "full-library class2 counter");
	}
	private void malformedSummary(String name, String field, String replacement, boolean duplicate){
		final String source=path(name+"_source");
		final ByteBuilder bytes=new ByteBuilder();
		try(Input in=new Input(measurement+".best.tsv.gz"); Output out=new Output(source+".best.tsv.gz")){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){out.writer.print(bytes.clear().append(row).nl());}
		}
		try(Input in=new Input(measurement+".summary.tsv"); Output out=new Output(source+".summary.tsv")){
			final LineParser1 p=new LineParser1('\t');
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(!p.termEquals(field,0)){out.writer.print(bytes.clear().append(row).nl());}
				else if(duplicate){
					out.writer.print(bytes.clear().append(row).nl().append(row).nl());
				}else if(replacement!=null){out.writer.print(bytes.clear().append(field).tab().append(replacement).nl());}
			}
		}
		write(name+".tsv", new ByteBuilder("library\tcapture\tresult_prefix\nfixture\tshard_0\t").append(source).nl());
	}
	private void summarize(String name, int captures){
		SummarizeNcrnaRejectedWindows.main(new String[]{"in="+path(name+".tsv"), "out="+path(name), "misses="+misses,
			"population="+population, "expectedlibraries=1", "expectedcaptures="+captures, "expectedrows=2", "expectedoriginals=2"});
	}
	private void reject(String name, int captures, String why){
		boolean rejected=false;
		try{summarize(name,captures);}catch(IllegalArgumentException expected){rejected=expected.getMessage().contains(why);}
		check(rejected, "The intended incomplete/duplicate collection must fail loudly: "+name);
	}
	private void write(String name, ByteBuilder bytes){try(Output out=new Output(path(name))){out.writer.print(bytes);}}
	private String path(String name){return new File(root,name).getPath();}
	private static void check(boolean valid, String why){if(!valid){throw new AssertionError(why);}}
	private final File root;
	private final String measurement, population, misses;
}
