package test.prok;

import java.io.File;
import java.util.Arrays;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.AnalyzeCobalaminInputPadding;
import prok.NcrnaBoundaryTrace;
import prok.NcrnaCmInputPadding;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Four literal inputs: shared alias half, both strands, a model switch,
 * a squared-error tie and one lost matched reference. No caller/search runs.
 * @author Brian Bushnell, Raiden
 */
public final class AnalyzeCobalaminInputPaddingTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, AnalyzeCobalaminInputPaddingTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, split[0], split.length>1 ? split[1] : null), "Use out=fresh_fixture_directory");
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append, "Finite fixtures need a fresh output directory");
			run(parser.out1);
			pp.outstream.println("COBALAMIN_INPUT_PADDING_FIXTURES_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}
	public static void run(final String root){
		check(new File(root).mkdir(), "Do not overwrite a prior finite fixture");
		make(root, 0);
		invoke(root);
		checkOutputs(root);
		selectedBases(root);
		pooledFit(root);
		for(int fault=1; fault<=6; fault++){
			final String variant=root+"/bad"+fault;
			check(new File(variant).mkdir(), "Each invalid-input test keeps its own literal inputs");
			make(variant, fault);
			boolean failed=false;
			try{invoke(variant);}catch(IllegalArgumentException expected){failed=true;}
			check(failed, "The reader must reject bad alias/offset/terminal/profile/sequence/reference evidence: "+fault);
		}
	}
	private static void make(final String root, final int fault){
		try(Output library=new Output(root+"/consensus.fa"); Output pads=new Output(root+"/padding.tsv");
			Output panel=new Output(root+"/panel.tsv"); Output audit=new Output(root+"/audit.tsv")){
			library.writer.println(">A\nAAAA\n>B\nCCCC\n>C\nGGGG");
			pads.writer.println("#BBTools CM input padding v1\n#policy\tP98_PLUS10\t5\t2");
			pads.writer.println("model_ordinal0\tmodel\tconsensus_sequence\tpad5\tpad3\ttraining_finite\ttraining_groups\tlearned");
			pads.writer.println(fault==5 ? "0\tA\tAAAC\t40\t30\t5\t2\ttrue" : "0\tA\tAAAA\t40\t30\t5\t2\ttrue");
			pads.writer.println("1\tB\tCCCC\t10\t20\t5\t2\ttrue\n2\tC\tGGGG\t20\t20\t0\t0\tfalse");
			panel.writer.println("input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff");
			audit.writer.println("assembly_id\tinput_id\tcanonical_fasta\thistorical_alias_group\tresolved_group\theld_out\tduplicate_original_group\tsampled");
			final ByteBuilder b=new ByteBuilder();
			for(int i=0; i<4; i++){
				panel.writer.print(b.clear().append(i).append("\tassembly_").append(i).tab().append(i)
					.append("\t/literal_").append(i).append(".fa\t/literal_").append(i).append(".gff\n"));
				audit.writer.print(b.clear().append(fault==1 && i==2 ? "WRONG_" : "assembly_").append(i).tab().append(i)
					.append("\t/literal_").append(i).append(".fa\t").append(i).tab().append(fault==8 && i==2 ? 2 : GROUPS[i])
					.append("\tfalse\tfalse\ttrue\n"));
			}
		}
		for(int arm=0; arm<2; arm++){
			final String name=arm==0 ? "baseline" : "candidate";
			try(Output trace=new Output(root+'/'+name+".trace.tsv"); Output profile=new Output(root+'/'+name+".profile.tsv");
				Output witnesses=new Output(root+'/'+name+".witnesses.tsv"); Output genomes=new Output(root+'/'+name+".genomes.tsv")){
				trace.writer.print(NcrnaBoundaryTrace.CM_HEADER);
				profile.writer.println("member\tcalls");
				profile.writer.println(fault==4 && arm==1 ? "cobalamin_rf00174\t3" : "cobalamin_rf00174\t4");
				witnesses.writer.println(WITNESSES);
				genomes.writer.println("input_id\tcalls\tmatched\tmasked\tunmatched");
				final ByteBuilder b=new ByteBuilder();
				for(int i=0; i<4; i++){
					final boolean selected=arm==0 || i!=3;
					final int strand=i==1 ? 1 : 0;
					final int start=arm==0 ? OLD_START[i] : NEW_START[i];
					final int stop=arm==0 ? OLD_STOP[i] : NEW_STOP[i];
					final String model=arm==1 && i==2 && fault<7 ? "B" : "A";
					final int from=i<2 ? 5 : i==2 ? 102 : 100;
					final int to=i<2 ? 114 : i==2 ? 203 : 199;
					final int refinedStart=selected ? (strand==0 ? start-1 : 1000-stop) : from;
					final int refinedStop=selected ? (strand==0 ? stop-1 : 1000-start) : to;
					trace.writer.print(b.clear().append("I\t").append(i).append("\t/literal_").append(i).append(".fa\n"));
					final String[] fields=new String[37];
					Arrays.fill(fields, "0");
					fields[0]="C";
					fields[1]=Integer.toString(i);
					fields[5]="contig";
					fields[6]=Integer.toString(strand);
					fields[7]="cobalamin_rf00174";
					fields[8]=model;
					fields[9]="1000";
					fields[10]=Integer.toString(from);
					fields[11]=Integer.toString(to);
					fields[14]=fields[15]=fields[16]=fields[20]=fields[21]=fields[31]="1";
					fields[22]="1000";
					fields[23]=Integer.toString(refinedStart);
					fields[24]=Integer.toString(refinedStop);
					fields[25]=selected ? Integer.toString(start-1) : "";
					fields[26]=selected ? Integer.toString(stop-1) : "";
					fields[27]=selected ? "0" : "";
					fields[28]="CM";
					fields[32]=selected ? "1" : "0";
					fields[33]=selected ? "SELECTED" : "PATH_REJECTED";
					fields[35]="999";
					fields[36]="4";
					b.clear();
					for(int field=0; field<fields.length; field++){
						if(field>0){b.tab();}
						b.append(fields[field]);
					}
					trace.writer.print(b.nl());
					trace.writer.print(b.clear().append("M\t").append(i).append("\t0\t0\t").append(from).tab().append(to)
						.tab().append(to-from+1).tab().append(refinedStart).tab().append(refinedStop)
						.tab().append(refinedStop-refinedStart+1).append("\t4\t4\t")
						.append(Float.floatToIntBits(40)).tab().append(Float.floatToIntBits(30))
						.append("\t-2\t0\t0\t3\t0\t0\tACCEPTED\n"));
					if(!(fault==3 && arm==1 && i==3)){
						trace.writer.print(b.clear().append("T\t").append(i).append("\t0\t0\tcontig\t1000\t1\t1\t")
							.append(selected ? 1 : 0).append("\t0\t0\n"));
					}
					if(selected){
						final long d5=strand==0 ? start-REF_START[i] : REF_STOP[i]-stop;
						final long d3=strand==0 ? stop-REF_STOP[i] : REF_START[i]-start;
						b.clear().append(i).append("\tr\tRF00174\t").append(model).append("\tcontig\t").append(start).tab()
							.append(stop).tab().append(REF_START[i]).tab().append(REF_STOP[i]).tab().append(strand==0 ? '+' : '-')
							.tab().append(d5+(fault==2 && i==1 && arm==1 ? 1 : 0)).tab().append(d3).nl();
						witnesses.writer.print(b);
						if(fault==6 && i==0 && arm==1){witnesses.writer.print(b);}
					}
					genomes.writer.print(b.clear().append(i).tab().append(selected ? "1\t1\t0\t0\n" : "0\t0\t0\t0\n"));
				}
			}
		}
	}
	private static void invoke(final String root){
		invoke(root, "result", new String[0]);
	}
	private static void invoke(final String root, final String prefix, final String[] extras){
		final String[] ordinary={
			"in="+root+"/baseline,"+root+"/candidate",
			"names=baseline,candidate",
			"traces="+root+"/baseline.trace.tsv,"+root+"/candidate.trace.tsv",
			"profiles="+root+"/baseline.profile.tsv,"+root+"/candidate.profile.tsv",
			"bindings="+root+"/panel.tsv", "audit="+root+"/audit.tsv",
			"library="+root+"/consensus.fa", "padding="+root+"/padding.tsv",
			"expectedinputs=4", "expectedharm5=2", "cmcpu=100", "out="+root+'/'+prefix};
		final String[] args=Arrays.copyOf(ordinary, ordinary.length+extras.length);
		System.arraycopy(extras, 0, args, ordinary.length, extras.length);
		AnalyzeCobalaminInputPadding.main(args);
	}
	private static void selectedBases(final String root){
		for(final int three:new int[]{10, 25}){
			final String base=root+"/base25_"+three+".tsv", prefix="selected25_"+three;
			base(base, three);
			invoke(root, prefix, new String[]{"basepadding="+base, "minselectioncalls=1", "minselectiongroups=1", "exportselected=t"});
			final NcrnaCmInputPadding.Data selected=readSelected(root+'/'+prefix+".selected.candidate.tsv");
			check(selected.selected && selected.base5==25 && selected.base3==three
				&& selected.pad5(0)==25 && selected.pad3(0)==30 && selected.pad5(1)==25 && selected.pad3(1)==20
				&& selected.pad5(2)==25 && selected.pad3(2)==three,
				"Each declared base keeps tied/worse five-prime and unsupported ends while exporting its improving three-prime overrides");
			final LineParser1 p=new LineParser1('\t');
			try(Input input=new Input(root+'/'+prefix+".costs.tsv")){
				input.reader.nextLine();
				for(byte[] line; (line=input.reader.nextLine())!=null;){
					p.set(line);
					if(p.termEquals("KEEP_BASE_SPLICE", 1)){
						check(p.parseDouble(5)==(three==10 ? 137532 : 170772),
							"Literal base work uses25five-prime and the declared10or25three-prime flank, including contig clipping");
					}else if(p.termEquals("SELECTION_HALF_MODEL_END", 1)){
						check(p.parseDouble(5)==176552, "The selected literal work uses25five-prime, not the old20fallback");
					}
				}
			}
		}
		for(int variant=7; variant<=8; variant++){
			final String directory=root+"/selection_groups"+variant;
			check(new File(directory).mkdir(), "Group fixtures use separate data-only directories");
			make(directory, variant);
			final String base=directory+"/base.tsv";
			base(base, 10);
			invoke(directory, "groups", new String[]{"basepadding="+base, "minselectioncalls=2", "minselectiongroups=2", "exportselected=t"});
			final NcrnaCmInputPadding.Data selected=readSelected(directory+"/groups.selected.candidate.tsv");
			check(selected.pad5(0)==25 && selected.pad3(0)==(variant==7 ? 10 : 30),
				"Two aliases in group0 do not satisfy two groups; distinct groups0and2 may select the improving end");
			invoke(directory, "thin", new String[]{"basepadding="+base, "minselectioncalls=3", "minselectiongroups=2", "exportselected=t"});
			check(readSelected(directory+"/thin.selected.candidate.tsv").pad3(0)==10,
				"Two selection loci cannot meet a three-locus floor even when their groups are distinct");
		}
	}
	private static void base(final String path, final int three){
		NcrnaCmInputPadding.writeFixed(path, 25, three, new String[]{"A", "B", "C"}, bases());
	}
	private static void pooledFit(final String root){
		invoke(root, "pooled", new String[]{"basepadding="+root+"/base25_10.tsv",
			"minselectioncalls=2", "minselectiongroups=2", "exportselected=t", "fitall=t"});
		final NcrnaCmInputPadding.Data selected=readSelected(root+"/pooled.selected.candidate.tsv");
		check(selected.pad5(0)==25 && selected.pad3(0)==30 && selected.pad3(1)==10,
			"Pooling includes A's groups0and1 to meet two-locus/two-group support; B's single group cannot acquire an override");
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(root+"/pooled.coverage.tsv")){
			final String[] names=columns(input, p);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				check(value(p, names, "fitting_common")==3 && value(p, names, "independent_judgment_common")==0,
					"All three common references are in-sample fitting data; none remains an independent judging observation");
			}
		}
		try(Input input=new Input(root+"/pooled.model_ends.tsv")){
			final String[] names=columns(input, p);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				if(p.termEquals("A", 1)){
					check(value(p, names, "selection_n")==2 && value(p, names, "selection_groups")==2,
						"The pooled fit retains two original resolved groups rather than replacing both by one fake selection ID");
				}
			}
		}
		try(Input input=new Input(root+"/pooled.estimates.tsv")){
			final String[] names=columns(input, p);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				check(p.termEquals("ALL_ELIGIBLE_FIXED_REFERENCE_ENDPOINT_SPLICE_IN_SAMPLE_NOT_COUPLED_CALLER", column(names, "qualification")),
					"Pooled endpoint splices cannot be labelled held-out or actual coupled caller measurements");
				if(p.termEquals("ALL", 2)){
					check(value(p, names, "in_sample_n")==3, "The in-sample estimate includes both former parity halves");
				}
			}
		}
	}
	private static NcrnaCmInputPadding.Data readSelected(final String path){
		return NcrnaCmInputPadding.read(path, new String[]{"A", "B", "C"}, bases());
	}
	private static byte[][] bases(){
		return new byte[][]{"AAAA".getBytes(java.nio.charset.StandardCharsets.US_ASCII),
			"CCCC".getBytes(java.nio.charset.StandardCharsets.US_ASCII), "GGGG".getBytes(java.nio.charset.StandardCharsets.US_ASCII)};
	}
	private static void checkOutputs(final String root){
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(root+"/result.coverage.tsv")){
			input.reader.nextLine();
			input.reader.nextLine();
			p.set(input.reader.nextLine());
			check(p.termEquals("candidate", 0) && p.parseInt(1)==3 && p.parseInt(3)==3 && p.parseInt(4)==2
				&& p.parseInt(5)==3 && p.parseInt(6)==2 && p.parseInt(7)==1,
				"Opposite-parity inputs1/2 share selection group0; input0 is judged instead, with one reference lost from the common cohort");
		}
		int masks=0;
		try(Input input=new Input(root+"/result.model_ends.tsv")){
			final String[] columns=columns(input, p);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				check(p.terms()==columns.length, "Model/end report fields must agree with the declared schema");
				final boolean end5=p.termEquals("5prime", 2);
				final boolean chosen=p.termEquals("true", column(columns, "use_candidate_pad"));
				if(p.termEquals("A", 1)){
					check(chosen!=end5 && value(p, columns, "selection_n")==1, "Only A's3prime error improves on the selection half");
					if(end5){
						check(value(p, columns, "selection_baseline_sse")==0 && value(p, columns, "selection_candidate_sse")==25
							&& value(p, columns, "worse")==2 && value(p, columns, "crossing3")==1
							&& value(p, columns, "worse_clipped_equal20")==2,
							"A must learn from input1's error5/SSE25, not input0's error1/SSE1; both flanks are clipped to five bases, with only one harm crossing3");
					}else{
						check(value(p, columns, "selection_baseline_sse")==100 && value(p, columns, "selection_candidate_sse")==0,
							"The selected3prime mask uses the literal selection-side squared errors");
					}
				}else if(p.termEquals("B", 1)){
					check(chosen!=end5, "B's equal squared5prime error retains20 despite a signed movement;3prime improves");
				}else{
					check(!chosen && value(p, columns, "selection_n")==0, "An unsupported unseen model retains20 at both ends");
				}
				masks++;
			}
		}
		check(masks==6, "The policy mask includes every model and both sense ends");
		int judged=0;
		try(Input input=new Input(root+"/result.estimates.tsv")){
			final String[] columns=columns(input, p);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				check(p.terms()==columns.length, "Every endpoint estimate uses the declared schema");
				if(!p.termEquals("ALL", 2)){continue;}
				check(value(p, columns, "judgment_n")==1, "Selection rows cannot leak into judgment metrics");
				if(p.termEquals("FULL_CANDIDATE", 1)){
					check(value(p, columns, "both_exact")==0 && value(p, columns, "prior_close_harms5")==1
						&& value(p, columns, "crossing3_5")==0 && value(p, columns, "within3_5")==1,
						"The judgment-side error1 is worse but stays within3; ordinal parity would incorrectly judge the error5 crossing instead");
				}else if(!p.termEquals("KEEP20_SPLICE", 1)){
					check(value(p, columns, "both_exact")==1 && value(p, columns, "prior_close_harms5")==0,
						"Both split-end estimates retain the exact5prime end and recover the3prime end on the judgment input");
				}
				judged++;
			}
		}
		check(judged==4, "All four declared endpoint policies report their same judgment cohort");
		int harms=0, staysWithin=0, crossings=0;
		try(Input input=new Input(root+"/result.harms.tsv")){
			final String[] columns=columns(input, p);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				check(p.terms()==columns.length && p.parseInt(13)==40 && p.parseInt(15)==5 && p.parseInt(16)==5,
					"Contig clipping removes the configured padding increase in both literal harms");
				if(p.termEquals("0", 1)){
					check(p.parseLong(10)==1 && p.termEquals("false", 11), "Input0 moves one base and remains within3");
					staysWithin++;
				}else{
					check(p.termEquals("1", 1) && p.parseLong(10)==5 && p.termEquals("true", 11),
						"Input1 moves five sense bases and crosses3");
					crossings++;
				}
				harms++;
			}
		}
		check(harms==2 && staysWithin==1 && crossings==1, "Two prior-close harms split into one within3 and one crossing3");
		try(Input input=new Input(root+"/result.costs.tsv")){
			final String[] columns=columns(input, p);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				check(p.terms()==columns.length && value(p, columns, "attempted_requests")==4,
					"Unselected searched requests remain in the work projection");
				if(p.termEquals("SELECTION_HALF_MODEL_END", 1)){
					check(p.parseDouble(4)==177932 && p.parseDouble(5)==170592,
						"Literal four-window dense-work totals are177932 full and170592 under the selected end mask");
				}
			}
		}
	}
	private static String[] columns(final Input input, final LineParser1 p){
		p.set(input.reader.nextLine());
		final String[] names=new String[p.terms()];
		for(int i=0; i<names.length; i++){names[i]=p.parseString(i);}
		return names;
	}
	private static int column(final String[] names, final String name){
		for(int i=0; i<names.length; i++){if(names[i].equals(name)){return i;}}
		throw new AssertionError("Missing declared output column: "+name);
	}
	private static long value(final LineParser1 p, final String[] columns, final String name){return p.parseLong(column(columns, name));}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private static final int[] OLD_START={6, 886, 103, 101}, OLD_STOP={115, 995, 204, 200};
	private static final int[] NEW_START={7, 896, 99, 101}, NEW_STOP={105, 990, 201, 200};
	private static final int[] GROUPS={1, 0, 0, 1};
	private static final int[] REF_START={6, 896, 101, 101}, REF_STOP={105, 995, 200, 200};
	private static final String WITNESSES="input_id\treference_id\temitted_member\tmodel\tseqid\tcall_start\tcall_stop\tref_start\tref_stop\tstrand\td5\td3";
}
