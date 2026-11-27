package test.prok;

import java.io.File;
import java.nio.file.Files;
import java.util.Arrays;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.AnalyzeCobalaminCmProposals;
import prok.CobalaminCmAcceptance;
import prok.ScoreCobalaminClanCandidates;
import prok.SelectCobalaminClanAcceptance;
import prok.SelectCobalaminClanAcceptance.Grid;
import prok.SummarizeCobalaminMemberAssignments;
import prok.SummarizeCobalaminMemberEndpoints;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.RrnaResourceIO;
import shared.Shared;
import structures.ByteBuilder;

/** Prefix counts are checked against an independent call-by-call oracle for
* every grid cell. The file fixture distinguishes native events with identical
* coordinates and preserves source-qualified references and unavailable scores.
* @author Brian, Ganyu
*/
public final class SelectCobalaminClanAcceptanceTest{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, SelectCobalaminClanAcceptanceTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, pair[0], pair.length>1 ? pair[1] : null), "Use out=fresh_directory");
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append, "Fresh fixture directory required");
			final File root=new File(parser.out1);
			check(root.mkdir(), "Fixture root already exists");
			kernel(false);
			kernel(true);
			adapter(new File(root, "adapter"));
			final Fixture good=new Fixture(root, "good");
			good.run();
			checkGood(good);
			final Fixture explicit=new Fixture(root, "explicit_defaults");
			explicit.run("refit=f", "upperidentity=.778");
			for(String suffix:REPORTS){sameRows(good.path("result"+suffix), explicit.path("result"+suffix));}
			final Fixture refit=new Fixture(root, "refit");
			refit.run("refit=t", "upperidentity=.75");
			checkRefit(refit);
			final Fixture oldTail=new Fixture(root, "upper_legacy");
			oldTail.run();
			checkUnmatched(oldTail, 0);
			final Fixture newTail=new Fixture(root, "upper_refit");
			newTail.run("refit=t", "upperidentity=.75");
			checkUnmatched(newTail, 2);
			final Fixture below=new Fixture(root, "below_upper");
			below.run("refit=t", "upperidentity=.75");
			checkUnmatched(below, 1);
			for(String fault:new String[]{"duplicate_target", "missing_target", "wrong_event", "wrong_trace", "missing_assignment", "missing_negative", "wrong_support"}){
				final Fixture bad=new Fixture(root, fault);
				boolean rejected=false;
				try{
					bad.run();
				}catch(IllegalArgumentException expected){
					rejected=true;
				}
				check(rejected && !new File(bad.path("result.choice.tsv")).exists(), "Bad binding must fail before grid publication: "+fault);
			}
			for(String fault:new String[]{"wrong_group_refit", "missing_tail_target"}){
				final Fixture bad=new Fixture(root, fault);
				reject(bad, "refit=t", "upperidentity=.75");
			}
			for(String value:new String[]{"NaN", "Infinity", ".59", ".70", ".80"}){
				reject(new Fixture(root, "invalid_upper_"+value), "upperidentity="+value);
			}
			final String fixed="fixedpolicy="+good.path("result.policy.tsv");
			good.write("gains.tsv", "input_id\treference_id\tseqid\tstart\tstop\tstrand\tchange\tpresence_class\tbaseline_member\tbaseline_model\tbaseline_model_in_subset\tcandidate_member\tcandidate_model\n");
			final String gains="gains="+good.path("gains.tsv"); // Valid empty roster: no unrelated schema error may satisfy the mode-rejection checks.
			reject(new Fixture(root, "fixed_refit"), "mode=fixed", fixed, gains, "expectedgains=0", "refit=t");
			reject(new Fixture(root, "fixed_upper"), "mode=fixed", fixed, gains, "expectedgains=0", "upperidentity=.75");
			pp.outstream.println("COBALAMIN_CLAN_ACCEPTANCE_FIXTURES_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Independent Oracle     ----------------*/
	/*--------------------------------------------------------------*/

	private static void kernel(boolean refit){
		final double[][] values={{0, 2, 4}, {1, 1, 3}, {2, 0, 2}, {3, 3, 1}, {4, 4, 0},
			{Double.NaN, Double.NEGATIVE_INFINITY, Double.NaN}, {-5, 2, 2}, {2, 2, 2}, {99, 99, 99}, {-99, -99, -99}, {3, 0, 1}, {1, 4, 3}};
		final int[] split={0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1};
		final int[] kind={0, 2, 0, 2, 0, 2, 0, 2, 0, 1, 1, 2};
		final boolean[] band={true, true, true, true, true, true, true, true, false, true, true, true};
		final long[][] baseline={{2, 1, 3}, {2, 1, 2}};
		final int[] originalSplit=split.clone();
		final Grid g=refit ? new Grid(values, split, kind, band, baseline, 20, 5, 4, 1.5, true)
			: new Grid(values, split, kind, band, baseline, 20, 5, 4, 1.5);
		check(Arrays.equals(split, originalSplit), "Pooling must not rewrite the original reporting partitions");
		check(g.axes[0][0]==(refit ? -5 : 0), "Only refit may use the finite judging minimum; masks and bypass scores never enter quantiles");
		double minimum=Double.NEGATIVE_INFINITY;
		double difference=Double.POSITIVE_INFINITY;
		final int[] best=new int[3];
		for(int a=0; a<g.axes[0].length; a++){
			for(int b=0; b<g.axes[1].length; b++){
				for(int c=0; c<g.axes[2].length; c++){
					final long[][] kept=new long[2][3];
					for(int i=0; i<values.length; i++){
						boolean pass=!band[i];
						for(int m=0; m<3; m++){pass|=Double.isFinite(values[i][m]) && values[i][m]>=g.axes[m][m==0 ? a : m==1 ? b : c];}
						if(pass){kept[split[i]][kind[i]]++;}
					}
					for(int s=0; s<3; s++){
						for(int k=0; k<3; k++){
							check(g.kept(a, b, c, s, k)==(s==2 ? kept[0][k]+kept[1][k] : kept[s][k]), "3D prefix differs from independent OR oracle");
						}
					}
					final double gain=(refit ? kept[0][0]+kept[1][0]-4 : kept[0][0]-2)/20.0;
					final double drop=(refit ? 5-kept[0][2]-kept[1][2] : 3-kept[0][2])/5.0;
					final double lo=Math.min(gain, drop);
					final double gap=Math.abs(gain-drop);
					if(lo>minimum || lo==minimum && gap<difference){
						minimum=lo;
						difference=gap;
						best[0]=a;
						best[1]=b;
						best[2]=c;
					}
				}
			}
		}
		check(Arrays.equals(best, g.best) && minimum==g.bestMinimum, "Fitted-population objective or deterministic tie order differs from independent enumeration");
		final double[][] altered=new double[values.length][];
		for(int i=0; i<values.length; i++){
			altered[i]=values[i].clone();
			if(!refit && split[i]==1 || kind[i]==1){Arrays.fill(altered[i], 123456.25);}
		}
		final Grid other=new Grid(altered, split, kind, band, baseline, 20, 5, 4, 1.5, refit);
		for(int m=0; m<3; m++){
			check(Arrays.equals(g.axes[m], other.axes[m]), "Ineligible values must not define grid cutoffs in either fitting mode");
			check(g.axes[m].length<=6, "Four quantiles plus common/reject controls cannot silently expand the configured grid");
		}
		check(Arrays.equals(g.best, other.best), "Ineligible outcomes must not select the cutoff triple");
	}

	private static void sameRows(String a, String b){
		try(Input first=new Input(a); Input second=new Input(b)){
			for(byte[] row; (row=first.reader.nextLine())!=null; ){
				check(Arrays.equals(row, second.reader.nextLine()), "Explicit legacy options must preserve every existing report row: "+a);
			}
			check(second.reader.nextLine()==null, "Explicit legacy options added report rows");
		}
	}
	private static void reject(Fixture f, String... options){
		boolean rejected=false;
		try{f.run(options);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && !new File(f.path("result.choice.tsv")).exists(), "Invalid refit/band binding must fail before fitting: "+f.fault);
	}
	private static void checkUnmatched(Fixture f, long expected){
		final LineParser1 p=new LineParser1('\t');
		boolean seen=false;
		try(Input input=new Input(f.path("result.comparison.tsv"))){
			input.reader.nextLine();
			for(byte[] line; (line=input.reader.nextLine())!=null; ){
				p.set(line);
				if(p.termStartsWith("BEST", 0) && p.termEquals("ALL", 1)){
					check(p.parseLong(3)==expected, "Exact .75 must bypass missing/impossible CM scores; just below it must not");
					seen=true;
				}
			}
		}
		check(seen, "Upper-band comparison must contain its complete population");
	}
	private static void checkRefit(Fixture f){
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(f.path("result.policy.tsv"))){
			input.reader.nextLine();
			for(int m=0; m<3; m++){
				final byte[] line=input.reader.nextLine();
				check(line!=null, "Positive pooled fixture needs a three-member policy");
				p.set(line);
				check(p.parseFloat(2)==.60f && p.parseFloat(3)==.75f, "Export must bind the requested product identity band");
			}
			check(input.reader.nextLine()==null, "Pooled policy cannot add target members");
		}
		int partitions=0;
		try(Input input=new Input(f.path("result.choice.tsv"))){
			input.reader.nextLine();
			for(byte[] line; (line=input.reader.nextLine())!=null; ){
				p.set(line);
				check(p.termStartsWith("ALL_ELIGIBLE_SAVED_REFIT", 10), "Refit choice must not claim held-out judging performance");
				partitions|=p.termEquals("SELECTION", 3) ? 1 : p.termEquals("JUDGING", 3) ? 2 : p.termEquals("ALL", 3) ? 4 : 0;
			}
		}
		check(partitions==7, "Refit must retain both original partitions and their pooled report");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Adapter Regression     ----------------*/
	/*--------------------------------------------------------------*/

	private static void adapter(File root) throws Exception{
		// This existing suite exercises repeated candidate coordinates, missing and
		// impossible own scores, strict/inexact assignments and native trace binding.
		AnalyzeCobalaminCmProposalsTest.run(root.getPath());
		final String r=root.getPath();
		AnalyzeCobalaminCmProposals.main(new String[]{"mode=acceptance", "nativecm=f", "eventassignments=t",
				"in="+r+"/capture_trace.tsv", "in2="+r+"/capture_grade", "bindings="+r+"/acceptance_bindings.tsv",
				"audit="+r+"/acceptance_audit.tsv", "profiles="+r+"/acceptance_profile.tsv", "rawgrades="+r+"/capture_raw.tsv",
				"names=literal", "expectedinputs=7", "acceptancescores="+r+"/capture_scores.tsv",
				"acceptancecosts="+r+"/capture_costs.tsv", "out="+r+"/companion"});
		for(String suffix:new String[]{"acceptance_calls", "acceptance_curves", "acceptance_coverage", "acceptance_native"}){
			check(Arrays.equals(Files.readAllBytes(new File(root, "capture."+suffix+".tsv").toPath()),
				Files.readAllBytes(new File(root, "companion."+suffix+".tsv").toPath())), "Optional event export changed an existing report");
		}
		int count=0;
		boolean impossible=false;
		boolean missing=false;
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(r+"/companion.event_assignments.tsv")){
			input.reader.nextLine();
			for(byte[] line; (line=input.reader.nextLine())!=null; ){
				p.set(line);
				count++;
				check(p.terms()==17, "Companion event schema changed");
				if(p.termEquals("RF01482", 10)){impossible=p.parseInt(4)==9 && p.parseLong(5)==17 && Float.intBitsToFloat(p.parseInt(15))==Float.NEGATIVE_INFINITY;}
				if(p.termEquals("RF01689", 10)){missing=p.parseInt(4)==9 && p.parseLong(5)==18 && p.termEquals("", 15);}
			}
		}
		check(count==8 && impossible && missing, "Native event IDs must distinguish same-coordinate minor calls and score states");
	}

	private static void checkGood(Fixture f){
		final LineParser1 p=new LineParser1('\t');
		boolean best=false;
		boolean common=false;
		try(Input input=new Input(f.path("result.comparison.tsv"))){
			input.reader.nextLine();
			for(byte[] line; (line=input.reader.nextLine())!=null; ){
				p.set(line);
				if(!p.termEquals("ALL", 1)){continue;}
				if(p.termStartsWith("BEST", 0)){
					best=p.parseLong(2)==3 && p.parseLong(3)==0 && p.parseLong(5)==2 && p.parseLong(6)==0 && p.parseLong(7)==1;
				}
				if(p.termStartsWith("COMMON", 0)){common=p.parseLong(2)==1 && p.parseLong(3)==1;}
			}
		}
		check(best && common, "Literal OR rescue, ordinary turnover and common missing-score fallback changed");
		boolean rejected=false;
		try(Input input=new Input(f.path("result.policy.tsv"))){
			// This is the existing own-CM loader's exact header check. An OR policy
			// must not be silently loadable with those different acceptance semantics.
			RrnaResourceIO.header(input.reader, "member\tupper_identity\tcm_minus_ga", f.path("result.policy.tsv"));
		}catch(IllegalArgumentException expected){
			rejected=true;
		}
		check(rejected, "The old own-CM filter must reject the explicit OR policy schema");
		try(Input input=new Input(f.path("result.policy.tsv"))){
			RrnaResourceIO.header(input.reader, "acceptance_policy\tmember\tlower_identity\tupper_identity\tcm_minus_ga\tunavailable_score", f.path("result.policy.tsv"));
			for(int m=0; m<3; m++){
				final byte[] line=input.reader.nextLine();
				check(line!=null, "The positive literal choice needs all three target policies");
				p.set(line);
				check(p.terms()==6 && p.termEquals("CLAN_OR_ANY_FINITE_PASS", 0) && p.termEquals(MEMBERS[m], 1)
					&& p.termEquals("FAIL_THIS_CM", 5), "OR semantics and unavailable-score behavior must be explicit");
			}
			check(input.reader.nextLine()==null, "Exactly three OR policy rows required");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          File Fixture        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Fixture{
		Fixture(File parent, String fault_){
			root=new File(parent, fault_);
			fault=fault_;
			check(root.mkdir(), "Fixture must be fresh");
			write("sources.tsv", "input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff\n0\ta0\tg0\t/no0.fa\t/no0.gff\n1\ta1\tg1\t/no1.fa\t/no1.gff\n2\tan\tneg\t/non.fa\t/non.gff\n");
			write("audit.tsv", "assembly_id\tinput_id\tcanonical_fasta\thistorical_alias_group\tresolved_group\theld_out\tduplicate_original_group\tsampled\na0\tg0\t/no0.fa\t0\t10\tfalse\tfalse\ttrue\na1\tg1\t/no1.fa\t1\t21\tfalse\tfalse\ttrue\nan\tneg\t/non.fa\t2\t42\tfalse\tfalse\ttrue\n");
			final String reference="input_id\treference_id\tseqid\tstart\tstop\tstrand\tscored\ng0\tr0\tC\t100\t199\t+\ttrue\ng1\tr1\tC\t100\t199\t+\ttrue\ng0\tr2\tC\t300\t399\t+\ttrue\ng1\tr3\tC\t300\t399\t+\ttrue\ng0\tmask\tC\t500\t599\t+\tfalse\n";
			write("before.references.tsv", reference);
			write("low.references.tsv", reference);
			try(Output output=new Output(path("annotations.tsv"))){
				output.writer.println(SummarizeCobalaminMemberEndpoints.ANNOTATIONS);
				for(int i=0; i<5; i++){
					final String input=i==1 || i==3 ? "g1" : "g0";
					final String support=i==2 ? "RF00174,RF01689" : i==3 ? "RF01482" : i==4 ? "NA" : "RF00174";
					final String label=i==2 ? "AMBIGUOUS" : i==4 ? "UNKNOWN" : support;
					output.writer.println(input+'\t'+(i==4 ? "mask" : "r"+i)+'\t'+(i!=4)+"\tRF00174\tRF00174\t"+(fault.equals("wrong_support") && i==0 ? "UNKNOWN" : label)+'\t'+support+"\tNA\tfalse\tNA\t/closed");
				}
			}
			try(Output before=new Output(path("before.queries.tsv")); Output low=new Output(path("low.queries.tsv"))){
				before.writer.println(SummarizeCobalaminMemberAssignments.QUERY_HEADER);
				low.writer.println(SummarizeCobalaminMemberAssignments.QUERY_HEADER);
				query(before, 0);
				query(before, 2);
				query(before, 3);
				query(before, 4);
				for(int i=0; i<6; i++){query(low, i);}
			}
			final String header=SummarizeCobalaminMemberAssignments.GENOME_HEADER+'\n';
			write("before.genomes.tsv", header+"g0\t2\t1\t0\t1\ng1\t2\t1\t0\t1\nneg\t0\t0\t0\t0\n");
			write("low.genomes.tsv", header+"g0\t4\t2\t1\t1\ng1\t2\t1\t0\t1\n"+(fault.equals("missing_negative") ? "" : "neg\t0\t0\t0\t0\n"));
			try(Output raw=new Output(path("raw.tsv")); Output assigned=new Output(path("assign.tsv")); Output score=new Output(path("scores.tsv"))){
				raw.writer.println(CobalaminCmAcceptance.HEADER);
				assigned.writer.println(AnalyzeCobalaminCmProposals.EVENT_ASSIGNMENT_HEADER);
				score.writer.println(ScoreCobalaminClanCandidates.SCORE_HEADER);
				for(int i=0; i<7; i++){
					final int j=i==6 ? 0 : i; // An unselected same-coordinate/model candidate has its own trace ID.
					final int ordinal=j==2 || j==4 ? 1 : 0;
					final String input=ordinal==0 ? "g0" : "g1";
					final int start=start(j);
					final String member=MEMBERS[member(j)];
					final boolean tail=fault.startsWith("upper_") || fault.equals("below_upper") || fault.equals("missing_tail_target");
					final float identity=j==2 ? .8f : tail && j==3 ? fault.equals("below_upper") ? Math.nextDown(.75f) : .75f
						: tail && j==4 ? .777f : .7f;
					final float own=j==3 ? Float.NaN : j==4 ? Float.NEGATIVE_INFINITY : j==5 ? 30 : 10;
					final String ownBits=Float.isNaN(own) ? "" : Integer.toString(Float.floatToIntBits(own));
					final String event=ordinal+"\tC\t+\t"+member+"\tcobalamin_"+member+"_toy\t1000\t"+(start-1)+'\t'+(start+98)+'\t'+Float.floatToIntBits(identity)+'\t'+ownBits+'\t'+Float.floatToIntBits(30)+"\tSCORE_ONLY\ttrue\t0\t999\t0\t"+i;
					raw.writer.println(event);
					if(i<6 && !(fault.equals("missing_assignment") && i==0)){
						final int group=ordinal==0 ? 10 : 21;
						assigned.writer.println(i+"\t"+input+'\t'+(fault.equals("wrong_group_refit") && i==0 ? group+2 : group)+'\t'+ordinal+"\t0\t"+(fault.equals("wrong_trace") && i==0 ? 99 : i)+"\tC\t+\t"+start+'\t'+(start+99)+'\t'+member+"\tcobalamin_"+member+"_toy\t"+kind(j)+'\t'+ref(j)+'\t'+Float.floatToIntBits(identity)+'\t'+ownBits+'\t'+Float.floatToIntBits(30));
					}
					if(identity>=.778f){continue;}
					for(int target=0; target<3; target++){
						if(fault.equals("missing_target") && i==0 && target==2){continue;}
						if(fault.equals("missing_tail_target") && i==3 && target==2){continue;}
						final float value=target==member(j) ? own : j==3 ? Float.NaN : j==4 ? Float.NEGATIVE_INFINITY : target==0 ? 50 : 10;
						final String targetBits=Float.isNaN(value) ? "" : Integer.toString(Float.floatToIntBits(value));
						final String line=(fault.equals("wrong_event") && i==0 ? 6 : i)+"\t"+event+'\t'+MEMBERS[target]+"\tcobalamin_"+MEMBERS[target]+"_toy\t"+targetBits+'\t'+Float.floatToIntBits(30)+"\tLITERAL\t0\t0\t0\t0\t0";
						score.writer.println(line);
						if(fault.equals("duplicate_target") && i==0 && target==0){score.writer.println(line);}
					}
				}
			}
		}
		void run(String... options){
			final String[] base={"in="+path("scores.tsv"), "in2="+path("assign.tsv"), "raw="+path("raw.tsv"), "baseline="+path("before"), "capture="+path("low"), "sources="+path("sources.tsv"), "audit="+path("audit.tsv"), "annotations="+path("annotations.tsv"), "commoncutoff=0", "expectedinputs=3", "expectedreferences=4", "expectedmasked=1", "baselineunmatched=2", "quantiles=4", "out="+path("result")};
			final String[] args=Arrays.copyOf(base, base.length+options.length);
			System.arraycopy(options, 0, args, base.length, options.length);
			SelectCobalaminClanAcceptance.main(args);
		}
		void query(Output out, int j){
			final String member=MEMBERS[member(j)];
			out.writer.println((j==2 || j==4 ? "g1" : "g0")+"\tC\t"+start(j)+'\t'+(start(j)+99)+"\t+\t"+member+"\tcobalamin_"+member+"_toy\t"+kind(j)+'\t'+(j<3)+'\t'+ref(j));
		}
		String path(String name){return new File(root, name).getPath();}
		void write(String name, String value){try(Output output=new Output(path(name))){output.writer.print(value);}}
		final File root;
		final String fault;
	}
	private static int start(int j){return j==1 ? 300 : j==3 || j==4 ? 700 : j==5 ? 500 : 100;}
	private static int member(int j){return j==1 ? 2 : j==2 || j==5 ? 0 : 1;}
	private static String kind(int j){return j<3 ? "TP" : j==5 ? "MASKED" : "UNMATCHED_OTHER";}
	private static String ref(int j){return j==0 ? "r0" : j==1 ? "r2" : j==2 ? "r1" : ".";}
	private static void check(boolean ok, String why){
		if(!ok){
			throw new AssertionError(why);
		}
	}
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final String[] REPORTS={".axes.tsv", ".surface.tsv.gz", ".choice.tsv", ".policy.tsv", ".comparison.tsv", ".references.tsv", ".genomes.tsv", ".breakdown.tsv"};
}
