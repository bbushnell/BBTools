package test.prok;

import java.io.File;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.ReplayNcrnaWindows;
import prok.SummarizeNcrnaRosterBaseline;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Reuses the independent saved-baseline literal fixture, adding acceptance-only
 * schemas, real gene-set parity failures and explicit genomic-background counts.
 * No alignments or new biological measurements are made by this report test.
 * @author Brian Bushnell, Raiden
 */
public final class SummarizeNcrnaAcceptanceTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SummarizeNcrnaAcceptanceTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown acceptance-summary fixture parameter: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append,
				"Use out=fresh_fixture_directory without input or overwrite");
			new SummarizeNcrnaAcceptanceTest(new File(parser.out1).getAbsoluteFile()).run();
			pp.outstream.println("NCRNA_ACCEPTANCE_SUMMARY_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private SummarizeNcrnaAcceptanceTest(final File root_){
		root=root_;
		check(root.mkdir(), "The acceptance report fixture requires a fresh directory");
	}
	private void run(){
		SummarizeNcrnaRosterBaselineTest.main(new String[]{"out="+path("literal")});
		copy("literal/population.frames.tsv", "population.frames.tsv");
		try(Input in=new Input(path("literal/population.original_calls.tsv.gz")); Output out=new Output(path("population.original_calls.tsv.gz"))){
			out.writer.print(row.clear().append(in.reader.nextLine()).nl());
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				row.clear();
				final boolean partial=p.parseInt(0)==10;
				for(int column=0; column<31; column++){
					if(column>0){row.tab();}
					if(partial && column==18){row.append(45);}
					else if(partial && column==19){row.append(80);}
					else if(partial && column==21){row.append("true");}
					else if(p.parseInt(0)==11 && column==18){row.append(2);}
					else{field(row, line, column);}
				}
				out.writer.print(row.nl());
			}
		}
		try(Input in=new Input(path("literal/native.inputs.tsv")); Output out=new Output(path("replay.inputs.tsv"))){
			in.reader.nextLine();
			out.writer.println(ReplayNcrnaWindows.INPUT_HEADER);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				final boolean first=p.parseInt(0)==0;
				row.clear();
				for(int column=0; column<9; column++){
					if(column>0){row.tab();}
					if(first && (column==5 || column==6 || column==8)){row.append(6);}
					else if(first && column==7){row.append(29);}
					else{field(row, line, column);}
				}
				out.writer.print(row.nl());
			}
		}
		try(Input in=new Input(path("literal/native.results.tsv")); Output out=new Output(path("replay.results.tsv"))){
			in.reader.nextLine();
			out.writer.println(ReplayNcrnaWindows.RESULT_HEADER);
			boolean background=false;
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				if(!background && p.parseInt(0)!=0){
					background(out);
					background=true;
					p.set(line);
				}
				if(model==null){model=p.parseString(18);}
				row.clear();
				for(int column=0; column<26; column++){
					if(column>=14 && column<=16){continue;}
					if(row.length()>0){row.tab();}
					field(row, line, column);
				}
				out.writer.print(row.nl());
			}
			check(background, "The original literal must contain both original source frames");
		}
		write("replay.replay.tsv", new ByteBuilder("metric\tvalue\nidentity\t0.64\n"));
		summarize("accepted", "literal/summary.genes.tsv.gz");
		metric("unique_best_members", 3);
		metric("captured_window_events", 7);
		metric("saved_window_bases", 1197);
		metric("recorded_input_bases", 5000);
		metric("actual_alignment_events", 35);
		metric("accepted_events_without_any_strand_roster_overlap_unadjudicated", 1);
		metric("accepted_events_without_same_strand_roster_overlap_unadjudicated", 2);
		metric("parity_original_genes_checked", 8);
		check(!new File(path("accepted.calls.tsv.gz")).exists() && !new File(path("accepted.model_gc.tsv")).exists(),
			"Acceptance-only reporting must not invent final-call metrics or refit GC models");
		try(Input in=new Input(path("accepted.genes.tsv.gz"))){
			check(new String(in.reader.nextLine(), java.nio.charset.StandardCharsets.UTF_8).equals(SummarizeNcrnaRosterBaseline.ACCEPTANCE_GENES_HEADER),
				"The explicit acceptance-only gene schema must omit final path fields");
		}
		try(Input in=new Input(path("accepted.cohorts.tsv"))){
			in.reader.nextLine();
			final long[][] expected={{8,3},{6,1},{1,1},{1,1},{2,0},{7,2}};
			for(final long[] counts:expected){
				p.set(in.reader.nextLine());
				check(p.parseLong(1)==counts[0] && p.parseLong(2)==counts[1], "Coverage/truncation subset counts must come from original metadata and actual attributed acceptance");
			}
			check(in.reader.nextLine()==null, "Report the coverage partition and distinct exact-span, truncation and joint-short subsets");
		}
		badParity();
	}
	private void background(final Output out){
		assert(model!=null) : "Background events use an actual construction model from the independent literal fixture";
		out.writer.print(row.clear().append("0\tframe_0\t").append(path("literal/older.fa"))
			.append("\t0\tsame%25name description\t1\t1000\t0\t779\t919\t80\t220\ti0e4\tACCEPTED_PROPOSAL\t8\t")
			.append(model).append("\t0.8\t1\t1\t100\t199\t20\t100\n"));
		out.writer.print(row.clear().append("0\tframe_0\t").append(path("literal/older.fa"))
			.append("\t1\tbackground\t0\t1000\t0\t700\t850\t700\t850\ti0e5\tACCEPTED_PROPOSAL\t7\t")
			.append(model).append("\t0.8\t1\t1\t720\t819\t20\t100\n"));
	}
	private void badParity(){
		for(int mode=0; mode<3; mode++){
			final String file="bad_parity_"+mode+".tsv.gz";
			try(Input in=new Input(path("literal/summary.genes.tsv.gz")); Output out=new Output(path(file))){
				out.writer.print(row.clear().append(in.reader.nextLine()).nl());
				for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
					p.set(line);
					if(mode==2 && p.parseInt(0)==20){continue;}
					if(mode==0 && p.parseInt(0)==10){
						p.setBounds(11);
						out.writer.print(row.clear().append(line, 0, p.a()).append("WINDOW_WITHOUT_ACCEPTED_MODEL")
							.append(line, p.b(), line.length-p.b()).nl());
					}else{out.writer.print(row.clear().append(line).nl());}
					if(mode==1 && p.parseInt(0)==10){out.writer.print(row.clear().append(line).nl());}
				}
			}
			boolean failed=false;
			try{summarize("rejected_"+mode, file);}catch(IllegalArgumentException | IllegalStateException expected){failed=true;}
			check(failed, "Changed, duplicated or missing original-gene parity rows must fail instead of accepting equal aggregate counts");
		}
	}
	private void summarize(final String output, final String parity){
		SummarizeNcrnaRosterBaseline.main(new String[]{"in="+path("replay"), "population="+path("population"),
			"assignments="+path("literal/assignments.tsv"), "out="+path(output), "expectedgenes=8", "expectedframes=4",
			"expectedmembers=6", "acceptanceonly=t", "identity=0.64", "parity="+path(parity)});
	}
	private void metric(final String key, final long expected){
		boolean found=false;
		try(Input in=new Input(path("accepted.summary.tsv"))){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				if(!p.termEquals(key, 0)){continue;}
				check(!found && p.parseLong(1)==expected, "The finite acceptance summary must preserve its independently expected count");
				found=true;
			}
		}
		check(found, "A required acceptance-summary metric is absent");
	}
	private void field(final ByteBuilder out, final byte[] line, final int column){p.setBounds(column);out.append(line, p.a(), p.b()-p.a());}
	private void copy(final String from, final String to){
		try(Input in=new Input(path(from)); Output out=new Output(path(to))){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){out.writer.print(row.clear().append(line).nl());}
		}
	}
	private String path(final String name){return new File(root, name).getPath();}
	private void write(final String name, final ByteBuilder bytes){try(Output out=new Output(path(name))){out.writer.print(bytes);}}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private final File root;
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder row=new ByteBuilder();
	private String model;
}
