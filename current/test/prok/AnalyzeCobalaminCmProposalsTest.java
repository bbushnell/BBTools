package test.prok;

import java.io.File;
import java.util.ArrayList;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.AnalyzeCobalaminCmProposals;
import prok.NcrnaBoundaryTrace;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaReferenceGff;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Independent seven-input oracle: two improvements, two bad shrinks, an
 * ambiguous reference, an unassociated call and an unavailable CM search.
 * One improving invocation is not selected. Repeated contig names and minus
 * coordinates must retain exact source/candidate/query linkage.
 * @author Brian Bushnell, Raiden
 */
public final class AnalyzeCobalaminCmProposalsTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, AnalyzeCobalaminCmProposalsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0], b=pair.length>1 ? pair[1] : null;
				check(parser.parse(arg, a, b), "Use out=fresh_fixture_directory");
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append, "The literal proposal fixtures require a fresh output directory");
			run(new File(parser.out1).getAbsolutePath());
			pp.outstream.println("COBALAMIN_CM_PROPOSAL_FIXTURES_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	public static void run(final String root){
		check(new File(root).mkdir(), "Proposal-analysis fixture directory must be fresh");
		final ByteBuilder row=new ByteBuilder();
		try(Output bindings=new Output(root+"/bindings.tsv"); Output trace=new Output(root+"/trace.tsv");
			Output refs=new Output(root+"/grade.references.tsv"); Output queries=new Output(root+"/grade.queries.tsv");
			Output genomes=new Output(root+"/grade.genomes.tsv"); Output profile=new Output(root+"/profile.tsv")){
			profile.writer.println("member\tcalls\ncobalamin_rf00174\t7\ncobalamin_rf01482\t0\ncobalamin_rf01689\t0");
			bindings.writer.println("input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff");
			trace.writer.print(NcrnaBoundaryTrace.CM_HEADER);
			refs.writer.println("input_id\treference_id\tseqid\tstart\tstop\tstrand\tscored");
			queries.writer.println("input_id\tseqid\tstart\tstop\tstrand\temitted_member\tmodel\toutcome\tscored_eligible\tmatched_reference_id");
			genomes.writer.println("input_id\tcalls\tmatched\tmasked\tunmatched");
			for(int input=0; input<7; input++){
				bindings.writer.print(row.clear().append(input).tab().append("assembly_").append(input).tab().append(input)
					.append("\t/literal_").append(input).append(".fa\t/literal_").append(input).append(".gff\n"));
				trace.writer.print(row.clear().append("I\t").append(input).append("\t/literal_").append(input).append(".fa\n"));
				final boolean selected=input!=5, available=input!=6;
				final int strand=input==2 ? 1 : 0;
				final int proposedStart=input==1 || input==2 ? 119 : input==0 || input==5 ? 109 : 99;
				final int proposedStop=input==1 || input==2 ? 168 : input==0 || input==5 ? 208 : 198;
				final int effectiveStart=available ? proposedStart : 99, effectiveStop=available ? proposedStop : 198;
				final int finalStart=strand==0 ? effectiveStart+1 : 1000-effectiveStop;
				final int finalStop=strand==0 ? effectiveStop+1 : 1000-effectiveStart;
				row.clear().append("C\t").append(input).append("\t0\t0\t0\tC\t").append(strand)
					.append("\tcobalamin_rf00174\tcobalamin_RF00174_toy\t1000\t99\t198\t0\t0\t1\t1\t1\t0\t0\t0\t1\t1\t1000\t")
					.append(effectiveStart).tab().append(effectiveStop).tab();
				if(selected){row.append(finalStart-1);}
				row.tab();
				if(selected){row.append(finalStop-1);}
				row.tab();
				if(selected){row.append(0);}
				trace.writer.print(row.append("\tCM\t0\t0\t1\t").append(selected ? 1 : 0).tab().append(selected ? "SELECTED" : "PATH_REJECTED")
					.append("\t0\t999\t100\n"));
				row.clear().append("M\t").append(input).append("\t0\t0\t99\t198\t100\t");
				if(available){row.append(proposedStart);}
				row.tab();
				if(available){row.append(proposedStop);}
				row.tab();
				if(available){row.append(proposedStop-proposedStart+1);}
				row.append("\t100\t100\t");
				if(available){row.append(Float.floatToIntBits(input==1 || input==2 ? 10 : 50));}
				row.tab().append(Float.floatToIntBits(30)).tab().append(available ? 32 : 0)
					.append("\t0\t0\t").append(available ? 2 : 0).tab();
				if(available){row.append(0);}
				trace.writer.print(row.tab().append(available ? 0 : 1).tab().append(!available ? "UNAVAILABLE" : input==3 || input==4 ? "UNCHANGED" : "APPLIED").nl());
				trace.writer.print(row.clear().append("T\t").append(input).append("\t0\t0\tC\t1000\t1\t1\t").append(selected ? 1 : 0).append("\t0\t0\n"));
				final int refStart=input==2 ? 802 : input==0 || input==5 ? 110 : 100;
				if(input!=4){
					refs.writer.print(row.clear().append(input).append("\tr\tC\t").append(refStart).tab().append(refStart+99).tab().append(strand==0 ? '+' : '-').append("\ttrue\n"));
				}
				if(input==3){refs.writer.println("3\tr2\tC\t110\t209\t+\ttrue");}
				if(selected){
					queries.writer.print(row.clear().append(input).append("\tC\t").append(finalStart).tab().append(finalStop).tab().append(strand==0 ? '+' : '-')
						.append("\tRF00174\tcobalamin_RF00174_toy n=1\t").append(input==4 ? "UNMATCHED_OTHER\tfalse\t." : "TP\ttrue\tr").nl());
				}
				genomes.writer.print(row.clear().append(input).tab().append(selected ? 1 : 0).tab().append(selected && input!=4 ? 1 : 0)
					.append("\t0\t").append(input==4 ? 1 : 0).nl());
			}
		}
		invoke(root, "trace.tsv", "grade", "result");
		final LineParser1 p=new LineParser1('\t');
		boolean shrink=false, improving=false, minus=false, absent=false, ambiguous=false, unselected=false;
		try(Input in=new Input(root+"/result.calls.tsv")){
			in.reader.nextLine();
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				check(p.terms()==48, "Every call row retains explicit empty fields and endpoint categories");
				final int input=p.parseInt(2);
				if(input==0){
					improving=p.parseLong(38)==-10 && p.parseLong(39)==0 && p.parseLong(40)==-10 && p.parseLong(41)==0 && p.termEquals("EXACT", 45);
				}else if(input==1){shrink=p.parseDouble(22)==.5 && p.termEquals("true", 42) && p.termEquals("true", 43);}
				else if(input==2){minus=p.parseInt(12)==802 && p.parseInt(15)==832 && p.parseLong(39)==20 && p.parseLong(41)==-30;}
				else if(input==3){ambiguous=p.termEquals("AMBIGUOUS_ORIGINAL_ORACLE", 35) && p.termEquals("", 38);}
				else if(input==5){unselected=p.termEquals("NOT_EMITTED", 10) && p.termEquals("true", 37);}
				else if(input==6){absent=p.termEquals("", 25) && p.termEquals("", 15) && p.termEquals("UNAVAILABLE", 34);}
			}
		}
		check(improving && shrink && minus && absent && ambiguous && unselected,
			"Literal endpoint oracles, both strands, absent score and unselected proposal must survive exact grade linking");
		boolean separation=false;
		try(Input in=new Input(root+"/result.tests.tsv")){
			in.reader.nextLine();
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				if(p.termEquals("ALL", 0) && p.termEquals("proposed_original_ratio_below_0.6", 1)){
					separation=p.parseInt(2)==4 && p.parseInt(3)==2 && p.parseInt(4)==0 && p.parseInt(5)==0 && p.parseInt(6)==2
						&& p.parseInt(7)==4 && p.parseInt(8)==0 && p.parseInt(9)==4 && p.parseDouble(10)==1;
				}
			}
		}
		check(separation, "Two literal bad shrinks and two improvements yield perfect separation with four exact proposal ends kept and four original ends restored");
		try(Input in=new Input(root+"/trace.tsv"); Output dropped=new Output(root+"/dropped.tsv"); Output omitted=new Output(root+"/omitted.tsv")){
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				final boolean pair=line.length>0 && line[0]!='#' && (p.termEquals("C", 0) || p.termEquals("M", 0)) && p.termEquals("5", 1);
				final boolean terminal=line.length>0 && line[0]!='#' && p.termEquals("T", 0) && p.termEquals("5", 1);
				if(!pair){dropped.writer.println(line);}
				if(!pair && !terminal){omitted.writer.println(line);}
			}
		}
		reject(root, "dropped.tsv", "grade", "dropped_result", "An omitted unselected C/M pair must violate the unchanged contig terminal count");
		reject(root, "omitted.tsv", "grade", "omitted_result", "An omitted whole unselected contig must violate the independent CM invocation profile");
		for(final String mode:new String[]{"strand", "outcome", "counts"}){
			for(final String suffix:new String[]{"queries", "references", "genomes"}){
				try(Input in=new Input(root+"/grade."+suffix+".tsv"); Output bad=new Output(root+"/"+mode+"."+suffix+".tsv")){
					for(byte[] line; (line=in.reader.nextLine())!=null;){
						String text=new String(line, java.nio.charset.StandardCharsets.UTF_8);
						if(suffix.equals("queries") && mode.equals("strand")){text=text.replace("\t+\tRF00174", "\t?\tRF00174");}
						if(suffix.equals("queries") && mode.equals("outcome")){text=text.replace("\tTP\t", "\tINVENTED\t");}
						if(suffix.equals("genomes") && mode.equals("counts") && text.equals("0\t1\t1\t0\t0")){text="0\t1\t0\t1\t0";}
						bad.writer.println(text);
					}
				}
			}
			reject(root, "trace.tsv", mode, mode+"_result", "Malformed strand/outcome or changed matched/masked counts must fail: "+mode);
		}
		scoreGuard(root);
		final File codecRoot=new File(root, "codec");
		NcrnaReferenceGffTest.run(codecRoot);
		acceptance(root);
		capture(root);
	}
	/** Two minor events share bounds, worker and read but have distinct native
	 * candidate IDs. One measured impossible score and one missing score must
	 * remain different; endpoint counts alone cannot account for either worker. */
	private static void capture(final String root){
		final LineParser1 p=new LineParser1('\t');
		final ByteBuilder b=new ByteBuilder();
		try(Input input=new Input(root+"/acceptance_trace.tsv"); Output trace=new Output(root+"/capture_trace.tsv");
			Output scores=new Output(root+"/capture_scores.tsv")){
			scores.writer.println(prok.CobalaminCmAcceptance.HEADER);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				if(line.length>0 && line[0]!='#' && p.terms()>1 && p.termEquals("5", 1)){
					if(p.termEquals("M", 0)){continue;}
					b.clear();
					for(int field=0; field<p.terms(); field++){
						if(field>0){b.tab();}
						if(p.termEquals("C", 0) && (field==14 || field==16)){b.append(0);}
						else if(p.termEquals("T", 0) && field==9){b.append(1);}
						else{b.append(p.parseString(field));}
					}
					trace.writer.print(b.nl());
				}else{trace.writer.println(line);}
				if(line.length==0 || line[0]=='#' || !p.termEquals("C", 0)){continue;}
				final int ordinal=p.parseInt(1), strand=p.parseInt(6), start=p.parseInt(23), stop=p.parseInt(24);
				b.clear().append(ordinal).append("\tC\t").append(strand==0 ? '+' : '-').append("\tRF00174\tcobalamin_RF00174_toy\t1000\t")
					.append(strand==0 ? start : 999-stop).tab().append(strand==0 ? stop : 999-start).tab().append(p.parseInt(18)).tab();
				if(ordinal!=6){b.append(Float.floatToIntBits(ordinal==1 || ordinal==2 ? 10 : 50));}
				scores.writer.print(b.tab().append(Float.floatToIntBits(30))
					.append(ordinal==5 ? "\tSCORE_ONLY" : "\tREUSED_ENDPOINT").append("\ttrue\t0\t999\t0\t0\n"));
			}
			for(int member=1; member<3; member++){
				final String rf=member==1 ? "RF01482" : "RF01689";
				trace.writer.print(b.clear().append("C\t0\t9\t").append(16+member).append("\t1\tC\t0\tcobalamin_")
					.append(rf.toLowerCase(java.util.Locale.ROOT)).append("\tcobalamin_").append(rf)
					.append("_minor\t1000\t99\t198\t0\t0\t0\t0\t0\t0\t0\t0\t1\t1\t1000\t99\t198\t99\t198\t0\tNONE\t0\t0\t1\t1\tSELECTED\t0\t999\t100\n"));
				scores.writer.print(b.clear().append("0\tC\t+\t").append(rf).append("\tcobalamin_").append(rf)
					.append("_minor\t1000\t99\t198\t").append(Float.floatToIntBits(.62f)).tab()
					.append(member==1 ? Integer.toString(Float.floatToIntBits(Float.NEGATIVE_INFINITY)) : "")
					.tab().append(Float.floatToIntBits(30)).append("\tSCORE_ONLY\ttrue\t")
					.append(member==1 ? "0\t999" : "\t").append("\t9\t").append(16+member).nl());
			}
			trace.writer.println("C\t0\t9\t19\t1\tC\t0\tcobalamin_rf01482\tcobalamin_RF01482_minor\t1000\t99\t198\t0\t0\t1\t1\t1\t0\t0\t0\t1\t1\t1000\t99\t198\t\t\t\tNN\t0\t0\t0\t0\tLENGTH_REJECTED\t0\t999\t100");
			trace.writer.println("T\t0\t9\t1\tC\t1000\t3\t2\t2\t2\t0");
		}
		try(Input input=new Input(root+"/acceptance_grade.queries.tsv"); Output output=new Output(root+"/capture_grade.queries.tsv");
			Input raw=new Input(root+"/acceptance_raw.tsv"); Output rawOut=new Output(root+"/capture_raw.tsv")){
			for(byte[] line; (line=input.reader.nextLine())!=null;){output.writer.println(line);}
			for(byte[] line; (line=raw.reader.nextLine())!=null;){rawOut.writer.println(line);}
			for(final String rf:new String[]{"RF01482", "RF01689"}){
				output.writer.println("0\tC\t100\t199\t+\t"+rf+"\tcobalamin_"+rf+"_minor\tUNMATCHED_OTHER\tfalse\t.");
				rawOut.writer.println("0\tQUERY\tcobalamin\tC\t100\t199\t+\tFP\t.\t99\t.\t.\t.\tfalse\tfalse\tcobalamin_"+rf+"_minor\t.\t.");
			}
		}
		for(final String suffix:new String[]{"references", "genomes"}){
			try(Input input=new Input(root+"/acceptance_grade."+suffix+".tsv"); Output output=new Output(root+"/capture_grade."+suffix+".tsv")){
				for(byte[] line; (line=input.reader.nextLine())!=null;){
					p.set(line);
					if(suffix.equals("genomes") && p.termEquals("0", 0)){output.writer.println("0\t3\t1\t0\t2");}
					else{output.writer.println(line);}
				}
			}
		}
		try(Output costs=new Output(root+"/capture_costs.tsv")){
			costs.writer.println(prok.CobalaminCmAcceptance.COST_HEADER);
			for(int ordinal=0; ordinal<7; ordinal++){
				costs.writer.print(b.clear().append(ordinal).append("\tRF00174\t1\t1\t").append(ordinal==6 ? 0 : 1)
					.append(ordinal==5 ? "\t0\t1\t2\t10\t5\t3\n" : "\t1\t0\t0\t0\t0\t0\n"));
				for(final String rf:new String[]{"RF01482", "RF01689"}){
					costs.writer.print(b.clear().append(ordinal).tab().append(rf).append(ordinal==0
						? rf.equals("RF01482") ? "\t1\t1\t1\t0\t1\t2\t10\t5\t3\n" : "\t1\t1\t0\t0\t1\t0\t4\t0\t2\n"
						: "\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"));
				}
			}
		}
		captureInvoke(root, "capture_scores.tsv", "capture_costs.tsv", "capture");
		int rows=0;
		boolean impossible=false, missing=false;
		try(Input input=new Input(root+"/capture.acceptance_calls.tsv")){
			input.reader.nextLine();
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				rows++;
				if(p.termEquals("RF01482", 8)){impossible=p.termEquals("-Infinity", 15) && p.parseFloat(13)==.62f;}
				if(p.termEquals("RF01689", 8)){missing=p.termEquals("NA", 15) && p.parseFloat(13)==.62f;}
			}
		}
		check(rows==8 && impossible && missing, "The two native minor IDs must join distinct impossible/missing scores without inventing endpoint M records");
		for(final String fault:new String[]{"identity", "traceid", "reuse", "missing"}){
			try(Input input=new Input(root+"/capture_scores.tsv"); Output output=new Output(root+"/capture_"+fault+".tsv")){
				for(byte[] line; (line=input.reader.nextLine())!=null;){
					p.set(line);
					if(p.terms()!=17 || !p.termEquals("0", 0)){
						output.writer.println(line);
						continue;
					}
					if(fault.equals("missing") && p.termEquals("RF01689", 3)){continue;}
					b.clear();
					for(int field=0; field<17; field++){
						if(field>0){b.tab();}
						if(fault.equals("identity") && field==8){b.append(Float.floatToIntBits(.99f));}
						else if(fault.equals("traceid") && field==16){b.append(999);}
						else if(fault.equals("reuse") && field==9){b.append(Float.floatToIntBits(0));}
						else{b.append(p.parseString(field));}
					}
					output.writer.print(b.nl());
				}
			}
			boolean rejected=false;
			try{captureInvoke(root, "capture_"+fault+".tsv", "capture_costs.tsv", "capture_bad_"+fault);}
			catch(IllegalArgumentException expected){rejected=true;}
			check(rejected, "Exact native identities, reuse bits and complete selected captures cannot be replaced by a coordinate-only match: "+fault);
		}
		try(Input input=new Input(root+"/capture_costs.tsv"); Output output=new Output(root+"/capture_bad_costs.tsv")){
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				if(p.terms()!=11 || !p.termEquals("0", 0) || !p.termEquals("RF00174", 1)){
					output.writer.println(line);
					continue;
				}
				b.clear();
				for(int field=0; field<11; field++){
					if(field>0){b.tab();}
					if(field==2){b.append(2);}else{b.append(p.parseString(field));}
				}
				output.writer.print(b.nl());
			}
		}
		boolean badCosts=false;
		try{captureInvoke(root, "capture_scores.tsv", "capture_bad_costs.tsv", "capture_bad_costs");}
		catch(IllegalArgumentException expected){badCosts=true;}
		check(badCosts, "Independent member counts must detect a missing or invented candidate event even when native endpoint totals still match");
	}
	private static void captureInvoke(final String root, final String scores, final String costs, final String output){
		AnalyzeCobalaminCmProposals.main(new String[]{"mode=acceptance", "nativecm=f", "in="+root+"/capture_trace.tsv",
			"in2="+root+"/capture_grade", "bindings="+root+"/acceptance_bindings.tsv", "audit="+root+"/acceptance_audit.tsv",
			"profiles="+root+"/acceptance_profile.tsv", "rawgrades="+root+"/capture_raw.tsv", "names=literal", "expectedinputs=7",
			"acceptancescores="+root+'/'+scores, "acceptancecosts="+root+'/'+costs, "out="+root+'/'+output});
	}
	/** Independent literal curves: ties cross both source splits, an ambiguous
	 * endpoint oracle keeps its recorded TP, and only the explicitly measured
	 * matching CM supplies a native score. Reference assignments never rematch.
	 */
	private static void acceptance(final String root){
		final LineParser1 p=new LineParser1('\t');
		final ByteBuilder b=new ByteBuilder();
		final float[] identity={.65f, .65f, .7f, .75f, .8f, .9f, .95f};
		try(Input in=new Input(root+"/trace.tsv"); Output out=new Output(root+"/acceptance_trace.tsv")){
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				if(line.length==0 || line[0]=='#' || !p.termEquals("C", 0)){
					out.writer.println(line);
					continue;
				}
				b.clear();
				for(int field=0; field<p.terms(); field++){
					if(field>0){b.tab();}
					if(field==18){b.append(Float.floatToIntBits(identity[p.parseInt(1)]));}
					else if(field==19){b.append(Float.floatToIntBits(.91f));}
					else{b.append(p.parseString(field));}
				}
				out.writer.print(b.nl());
			}
		}
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		final ArrayList<ArrayList<NcrnaReferenceGff.Row>> nativeRows=new ArrayList<ArrayList<NcrnaReferenceGff.Row>>();
		for(int i=0; i<7; i++){nativeRows.add(new ArrayList<NcrnaReferenceGff.Row>());}
		final String[][] ids=new String[7][2];
		try(Input in=new Input(root+"/grade.references.tsv"); Output out=new Output(root+"/acceptance_grade.references.tsv")){
			out.writer.println(in.reader.nextLine());
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				final int input=p.parseInt(0), which=p.termEquals("r2", 1) ? 1 : 0;
				final NcrnaReferenceGff.Scores scores=input==0
					? codec.scores("RF01482,RF00174", new long[]{9990, 450}, new String[]{"30", "30"})
					: input==2 ? codec.scores("RF00174", new long[]{Long.MIN_VALUE}, new String[]{"30"})
					: input==3 ? codec.scores("RF01482", new long[]{1000}, new String[]{"30"}) : null;
				final NcrnaReferenceGff.Row ref=codec.create(Integer.toString(input),
					new Locus("C", "cobalamin", input==2 ? 850 : p.parseInt(3), input==2 ? 859 : p.parseInt(4), p.termEquals("-", 5) ? '-' : '+'),
					input==0, input==0 ? "historical_mask" : null, "RF01482", "literal_evidence", scores);
				ids[input][which]=ref.id;
				if(input!=6){nativeRows.get(input).add(ref);}
				out.writer.print(b.clear().append(input).tab().append(ref.id).append("\tC\t").append(ref.locus.start)
					.tab().append(ref.locus.stop).tab().append(ref.locus.strand).append("\ttrue\n"));
			}
		}
		try(Input in=new Input(root+"/grade.queries.tsv"); Output out=new Output(root+"/acceptance_grade.queries.tsv")){
			out.writer.println(in.reader.nextLine());
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				b.clear();
				for(int field=0; field<p.terms(); field++){
					if(field>0){b.tab();}
					if(field==9 && p.termEquals("TP", 7)){b.append(ids[p.parseInt(0)][0]);}
					else{b.append(p.parseString(field));}
				}
				out.writer.print(b.nl());
			}
		}
		try(Input in=new Input(root+"/grade.genomes.tsv"); Output out=new Output(root+"/acceptance_grade.genomes.tsv")){
			for(byte[] line; (line=in.reader.nextLine())!=null;){out.writer.println(line);}
		}
		try(Output bindings=new Output(root+"/acceptance_bindings.tsv"); Output audit=new Output(root+"/acceptance_audit.tsv")){
			bindings.writer.println("input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff");
			audit.writer.println("assembly_id\tinput_id\tcanonical_fasta\thistorical_alias_group\tresolved_group\theld_out\tduplicate_original_group\tsampled");
			for(int input=0; input<7; input++){
				final String path=root+"/native_"+input+".gff";
				nativeRows.get(input).add(codec.create(Integer.toString(input), new Locus("C", "tRNA", 1, 70, '+'),
					false, null, "RF00005", "unrelated_family"));
				codec.write(path, Integer.toString(input), nativeRows.get(input));
				check(codec.readFamily(path, Integer.toString(input), "cobalamin").size()==nativeRows.get(input).size()-1,
					"The family-scoped codec retains every requested row and omits unrelated ncRNAs");
				bindings.writer.print(b.clear().append(input).tab().append("assembly_").append(input).tab().append(input)
					.append("\t/literal_").append(input).append(".fa\t").append(path).nl());
				audit.writer.print(b.clear().append("assembly_").append(input).tab().append(input)
					.append("\t/historical_alias_").append(input).append(".fa\t").append(input).tab().append(input).append("\tfalse\tfalse\ttrue\n"));
			}
		}
		try(Output profile=new Output(root+"/acceptance_profile.tsv")){
			profile.writer.println("member\tcalls\ncobalamin_rf00174\t7");
		}
		try(Input in=new Input(root+"/acceptance_grade.queries.tsv"); Output raw=new Output(root+"/acceptance_raw.tsv")){
			raw.writer.println(prok.GradeNcrnaReference.RAW_HEADER);
			in.reader.nextLine();
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				final int input=p.parseInt(0);
				final boolean matched=p.termEquals("TP", 7);
				raw.writer.print(b.clear().append(input).append("\tQUERY\tcobalamin\tC\t").append(p.parseInt(2)).tab().append(p.parseInt(3))
					.tab().append(p.parseString(4)).tab().append(input==2 ? "TP_INEXACT" : matched ? "TP" : "FP")
					.append("\t.\t").append(input+1).tab().append(p.parseString(9)).append("\t.\t.\tfalse\t")
					.append(matched).tab().append(p.parseString(6)).append("\t.\t.\n"));
			}
		}
		AnalyzeCobalaminCmProposals.main(new String[]{"mode=acceptance", "members=RF00174", "in="+root+"/acceptance_trace.tsv",
			"in2="+root+"/acceptance_grade", "bindings="+root+"/acceptance_bindings.tsv", "audit="+root+"/acceptance_audit.tsv",
			"profiles="+root+"/acceptance_profile.tsv", "rawgrades="+root+"/acceptance_raw.tsv", "names=literal", "expectedinputs=7", "out="+root+"/acceptance"});
		int rows=0;
		boolean nativeMatch=false, noAttributes=false, noHit=false, differentModel=false, overlay=false, strictDifferent=false, residual=false;
		try(Input in=new Input(root+"/acceptance.acceptance_calls.tsv")){
			in.reader.nextLine();
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				rows++;
				final int input=p.parseInt(1);
				check(p.terms()==27 && p.parseFloat(13)==identity[input] && p.parseFloat(14)==.91f
					&& p.termEquals("NA", 24), "Acceptance identity is distinct from endpoint identity, and absent HBM is not zero");
				if(input==0){nativeMatch=p.parseDouble(21)==45 && p.parseDouble(23)==5;}
				if(input==1){noAttributes=p.termEquals("NO_SCORE_ATTRIBUTES", 20);}
				if(input==2){noHit=p.termEquals("NO_CURRENT_NATIVE_HIT", 20);}
				if(input==1){strictDifferent=p.termEquals("false", 11) && p.termEquals("true", 26);}
				if(input==2){residual=p.termEquals("true", 11) && p.termEquals("true", 26);}
				if(input==3){differentModel=p.termEquals("NO_SAME_CM_MEASUREMENT", 20) && p.termEquals("TP", 10);}
				if(input==6){overlay=p.termEquals("REFERENCE_NOT_IN_BOUND_GFF", 20) && p.termEquals("NA", 15);}
			}
		}
		check(rows==6 && nativeMatch && noAttributes && noHit && differentModel && overlay && strictDifferent && residual,
			"Only selected calls enter curves; native scores use the measured CM and retain the current grade label across an older mask disposition");
		boolean selection=false, judging=false, missing=false;
		try(Input in=new Input(root+"/acceptance.acceptance_curves.tsv")){
			in.reader.nextLine();
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				if(!p.termEquals("FAMILY", 1)){continue;}
				check(!p.termEquals("accepted_identity", 3) || p.parseDouble(4)!=(double).75f,
					"A judging-only identity value cannot enter the selection cutoff grid");
				if(p.termEquals("accepted_identity", 3) && p.parseDouble(4)==(double).7f){
					if(p.termEquals("SELECTION", 5)){
						selection=p.parseInt(6)==2 && p.parseInt(7)==1 && p.parseInt(8)==1 && p.parseInt(10)==1;
					}else{judging=p.parseInt(6)==1 && p.parseInt(7)==1 && p.parseInt(10)==0 && p.parseInt(11)==0;}
				}
				if(p.termEquals("cm_minus_ga", 3) && p.termEquals("SELECTION", 5)){
					missing=p.parseInt(14)==1;
				}
			}
		}
		check(selection && judging && missing, "Inclusive tied cutoffs conserve source-split matched/inexact counts and keep missing CM scores separate");
		boolean offset=false;
		try(Input in=new Input(root+"/acceptance.acceptance_native.tsv")){
			in.reader.nextLine();
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				if(p.termEquals("SELECTION", 1) && p.termEquals("MATCHED", 2)){
					offset=p.parseInt(3)==1 && p.parseInt(4)==0 && p.parseDouble(5)==5 && p.termEquals("NA", 6)
						&& p.parseDouble(7)==5 && p.parseDouble(11)==5;
				}
			}
		}
		check(offset, "The measured matching-CM pair has offset five bits; one observation has no fabricated sample standard deviation");
	}
	/** Original spans are 100..199 (+) and 802..901 (-). Two low-score
	 * proposals restore both true ends, two spoil exact ends; one oracle is
	 * ambiguous and one is masked. An untraced minor-member TP stays in totals.
	 */
	private static void scoreGuard(final String root){
		final LineParser1 p=new LineParser1('\t');
		final ByteBuilder b=new ByteBuilder();
		try(Input input=new Input(root+"/trace.tsv"); Output output=new Output(root+"/guard_trace.tsv")){
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				if(line.length==0 || line[0]=='#'){
					output.writer.println(line);
					continue;
				}
				p.set(line);
				final int inputID=p.parseInt(1);
				if(inputID==6 || (!p.termEquals("C", 0) && !p.termEquals("M", 0))){
					output.writer.println(line);
					continue;
				}
				b.clear();
				for(int field=0; field<p.terms(); field++){
					if(field>0){b.tab();}
					if(p.termEquals("C", 0) && (field==23 || field==24)){b.append(field==23 ? 99 : 198);}
					else if(p.termEquals("C", 0) && inputID!=5 && (field==25 || field==26)){
						b.append(inputID==2 ? field==25 ? 801 : 900 : field==25 ? 99 : 198);
					}else if(p.termEquals("M", 0) && field==12){b.append(Float.floatToIntBits(10));}
					else if(p.termEquals("M", 0) && field==20){b.append("SCORE_REJECTED");}
					else{b.append(p.parseString(field));}
				}
				output.writer.print(b.nl());
			}
			guardExtraTrace(output, 7, 139, 158);
			guardExtraTrace(output, 8, 111, 210);
			guardExtraTrace(output, 9, 109, 208);
		}
		try(Input input=new Input(root+"/grade.queries.tsv"); Output output=new Output(root+"/guard_grade.queries.tsv")){
			output.writer.println(input.reader.nextLine());
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				final int inputID=p.parseInt(0);
				b.clear();
				for(int field=0; field<p.terms(); field++){
					if(field>0){b.tab();}
					if(field==2 || field==3){b.append(inputID==2 ? field==2 ? 802 : 901 : field==2 ? 100 : 199);}
					else if(inputID==4 && field==7){b.append("MASKED");}
					else{b.append(p.parseString(field));}
				}
				output.writer.print(b.nl());
			}
			output.writer.println("0\tC\t500\t599\t+\tRF01482\tcobalamin_RF01482_minor\tTP\ttrue\tminor");
			for(int inputID=7; inputID<10; inputID++){
				output.writer.print(b.clear().append(inputID).append("\tC\t100\t199\t+\tRF00174\tcobalamin_RF00174_toy n=1\tTP\ttrue\tr\n"));
			}
		}
		try(Input input=new Input(root+"/grade.references.tsv"); Output output=new Output(root+"/guard_grade.references.tsv")){
			for(byte[] line; (line=input.reader.nextLine())!=null;){output.writer.println(line);}
			output.writer.println("0\tminor\tC\t500\t599\t+\ttrue");
			output.writer.println("4\tmasked\tC\t100\t199\t+\tfalse");
			output.writer.println("7\tr\tC\t140\t159\t+\ttrue");
			output.writer.println("8\tr\tC\t110\t209\t+\ttrue");
			output.writer.println("9\tr\tC\t120\t219\t+\ttrue");
		}
		try(Input input=new Input(root+"/grade.genomes.tsv"); Output output=new Output(root+"/guard_grade.genomes.tsv")){
			output.writer.println(input.reader.nextLine());
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				if(p.termEquals("0", 0)){output.writer.println("0\t2\t2\t0\t0");}
				else if(p.termEquals("4", 0)){output.writer.println("4\t1\t0\t1\t0");}
				else{output.writer.println(line);}
			}
			for(int inputID=7; inputID<10; inputID++){output.writer.print(b.clear().append(inputID).append("\t1\t1\t0\t0\n"));}
		}
		try(Input input=new Input(root+"/bindings.tsv"); Output output=new Output(root+"/guard_bindings.tsv")){
			for(byte[] line; (line=input.reader.nextLine())!=null;){output.writer.println(line);}
			for(int inputID=7; inputID<10; inputID++){
				output.writer.print(b.clear().append(inputID).tab().append("assembly_").append(inputID).tab().append(inputID)
					.append("\t/literal_").append(inputID).append(".fa\t/literal_").append(inputID).append(".gff\n"));
			}
		}
		try(Output output=new Output(root+"/guard_profile.tsv")){
			output.writer.println("member\tcalls\trejected_score\ncobalamin_rf00174\t10\t9");
		}
		try(Input input=new Input(root+"/guard_trace.tsv"); Output output=new Output(root+"/guard_omitted.tsv")){
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				final boolean missing=line.length>0 && line[0]!='#' && p.terms()>1 && p.termEquals("5", 1)
					&& (p.termEquals("C", 0) || p.termEquals("M", 0) || p.termEquals("T", 0));
				if(!missing){output.writer.println(line);}
			}
		}
		guardInvoke(root, "guard_trace.tsv", "guard_profile.tsv", "RF00174", "guard");
		int summaries=0;
		try(Input input=new Input(root+"/guard.score_guard_summary.tsv")){
			p.set(input.reader.nextLine());
			final String[] fields=new String[p.terms()];
			for(int i=0; i<fields.length; i++){fields[i]=p.parseString(i);}
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				summaries++;
				if(p.termEquals("ALL_HELD_ORIGINAL_ORACLE", 1)){
					check(value(p, fields, "held")==9 && value(p, fields, "scorable")==6
						&& value(p, fields, "no_scored_reference")==2 && value(p, fields, "ambiguous_scored_reference")==1
						&& value(p, fields, "no_scored_reference_masked_overlap")==1 && value(p, fields, "original_both_exact")==2
						&& value(p, fields, "original_inexact")==4 && value(p, fields, "inexact_to_both_exact")==2
						&& value(p, fields, "original_exact_spoiled")==2 && value(p, fields, "inexact_to_both_within3_nonexact")==1
						&& value(p, fields, "inexact_to_any_end_beyond3")==1 && value(p, fields, "inexact_any_end_worse")==0,
						"Nine holds retain six unique oracles, one ambiguity and two unscorable originals; a still-beyond-three improvement must not become a worse proposal");
				}else if(p.termEquals("SELECTED_HELD_RECORDED_TP_REFERENCE", 1)){
					check(value(p, fields, "held")==7 && value(p, fields, "scorable")==7
						&& value(p, fields, "original_both_exact")==3 && value(p, fields, "inexact_to_both_exact")==2,
						"Recorded TP assignments retain both ambiguous-original and below-half inclusive matches without weakening the original oracle");
				}else if(p.termEquals("UNSELECTED_HELD_ORIGINAL_ORACLE", 1)){
					check(value(p, fields, "held")==1 && value(p, fields, "not_emitted")==1 && value(p, fields, "inexact_to_both_exact")==1,
						"An exact unselected proposal remains a request counterfactual, not an emitted recovery");
				}
			}
		}
		check(summaries==4, "Every arm reports all, selected, unselected and recorded-TP reference scopes exactly once");
		boolean inclusive=false, near=false, improvedBeyond=false;
		try(Input input=new Input(root+"/guard.score_guard_calls.tsv")){
			input.reader.nextLine();
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				if(p.termEquals("7", 2)){inclusive=p.termEquals("TP", 10) && p.termEquals("NO_ORIGINAL_ORACLE", 35) && p.termEquals("", 38);}
				else if(p.termEquals("8", 2)){near=p.parseLong(38)==-10 && p.parseLong(39)==2 && p.termEquals("WITHIN3", 45);}
				else if(p.termEquals("9", 2)){improvedBeyond=p.parseLong(38)==-20 && p.parseLong(39)==-10
					&& p.termEquals("BEYOND3", 45) && p.termEquals("false", 42) && p.termEquals("false", 43);}
			}
		}
		check(inclusive && near && improvedBeyond, "Inclusive assignment, within-three nonexact accuracy and beyond-three improvement are distinct literal outcomes");
		try(Output output=new Output(root+"/guard_bad_profile.tsv")){
			output.writer.println("member\tcalls\trejected_score\ncobalamin_rf00174\t10\t8");
		}
		boolean rejected=false;
		try{guardInvoke(root, "guard_trace.tsv", "guard_bad_profile.tsv", "RF00174", "guard_bad_count");}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Independent rejected-score totals catch a missing or misclassified held request");
		rejected=false;
		try{guardInvoke(root, "guard_trace.tsv", "guard_profile.tsv", "RF01482", "guard_bad_scope");}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "A main-only trace must not masquerade as another member's declared scope");
		try(Output output=new Output(root+"/guard_missing_profile.tsv")){
			output.writer.println("member\tcalls\trejected_score");
		}
		rejected=false;
		try{guardInvoke(root, "guard_trace.tsv", "guard_missing_profile.tsv", "RF00174", "guard_missing_enabled_member");}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "An enabled member cannot disappear from the independent profile even when all disabled members are absent");
		rejected=false;
		try{guardInvoke(root, "guard_omitted.tsv", "guard_profile.tsv", "RF00174", "guard_missing_enabled_contig");}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "An omitted whole unselected enabled-member contig must fail under the explicit main-only scope too");
		for(final String fault:new String[]{"disjoint", "wrong_strand", "missing_id"}){
			final String grade="guard_"+fault;
			for(final String suffix:new String[]{"queries", "references", "genomes"}){
				try(Input input=new Input(root+"/guard_grade."+suffix+".tsv"); Output output=new Output(root+'/'+grade+'.'+suffix+".tsv")){
					for(byte[] line; (line=input.reader.nextLine())!=null;){
						p.set(line);
						if(!suffix.equals("references") || !p.termEquals("7", 0)){
							output.writer.println(line);
							continue;
						}
						b.clear();
						for(int field=0; field<p.terms(); field++){
							if(field>0){b.tab();}
							if(fault.equals("disjoint") && field==3){b.append(400);}
							else if(fault.equals("disjoint") && field==4){b.append(419);}
							else if(fault.equals("wrong_strand") && field==5){b.append('-');}
							else if(fault.equals("missing_id") && field==1){b.append("absent");}
							else{b.append(p.parseString(field));}
						}
						output.writer.print(b.nl());
					}
				}
			}
			rejected=false;
			try{guardInvoke(root, "guard_trace.tsv", grade, "guard_profile.tsv", "RF00174", grade+"_result");}
			catch(IllegalArgumentException expected){rejected=true;}
			check(rejected, "Inclusive recorded TPs still reject disjoint, wrong-strand and missing-ID assignments: "+fault);
		}
	}
	private static void guardExtraTrace(final Output output, final int inputID, final int proposalStart, final int proposalStop){
		check(inputID>=7 && inputID<=9, "Three additional literal requests have fixed source IDs and positive-overlap assignments");
		final ByteBuilder b=new ByteBuilder();
		output.writer.print(b.append("I\t").append(inputID).append("\t/literal_").append(inputID).append(".fa\n"));
		output.writer.print(b.clear().append("C\t").append(inputID)
			.append("\t0\t0\t0\tC\t0\tcobalamin_rf00174\tcobalamin_RF00174_toy\t1000\t99\t198\t0\t0\t1\t1\t1\t0\t0\t0\t1\t1\t1000\t99\t198\t99\t198\t0\tCM\t0\t0\t1\t1\tSELECTED\t0\t999\t100\n"));
		output.writer.print(b.clear().append("M\t").append(inputID).append("\t0\t0\t99\t198\t100\t").append(proposalStart)
			.tab().append(proposalStop).tab().append(proposalStop-proposalStart+1).append("\t100\t100\t")
			.append(Float.floatToIntBits(10)).tab().append(Float.floatToIntBits(30)).append("\t32\t0\t0\t2\t0\t0\tSCORE_REJECTED\n"));
		output.writer.print(b.clear().append("T\t").append(inputID).append("\t0\t0\tC\t1000\t1\t1\t1\t0\t0\n"));
	}
	private static long value(final LineParser1 p, final String[] fields, final String name){
		for(int i=0; i<fields.length; i++){if(fields[i].equals(name)){return p.parseLong(i);}}
		throw new AssertionError("Missing score-guard summary field: "+name);
	}
	private static void guardInvoke(final String root, final String trace, final String profile, final String members, final String result){
		guardInvoke(root, trace, "guard_grade", profile, members, result);
	}
	private static void guardInvoke(final String root, final String trace, final String grade, final String profile, final String members, final String result){
		AnalyzeCobalaminCmProposals.main(new String[]{"mode=scoreguard", "members="+members, "in="+root+'/'+trace,
			"in2="+root+'/'+grade, "bindings="+root+"/guard_bindings.tsv", "profiles="+root+'/'+profile,
			"names=literal", "expectedinputs=10", "out="+root+'/'+result});
	}
	private static void invoke(final String root, final String trace, final String grade, final String result){
		AnalyzeCobalaminCmProposals.main(new String[]{"in="+root+'/'+trace, "in2="+root+'/'+grade, "bindings="+root+"/bindings.tsv",
			"profiles="+root+"/profile.tsv", "names=literal", "expectedinputs=7", "out="+root+'/'+result});
	}
	private static void reject(final String root, final String trace, final String grade, final String result, final String why){
		boolean failed=false;
		try{invoke(root, trace, grade, result);}catch(IllegalArgumentException expected){failed=true;}
		check(failed, why);
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private AnalyzeCobalaminCmProposalsTest(){}
}
