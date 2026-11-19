package prok;

import java.io.File;
import java.util.Arrays;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Summarizes saved regression gate cutoffs by median across three seeds.
 * Choices use only the internal monitor tables. Outer rows confirm the identical
 * chosen cutoffs; no network, score cache, fit or caller is evaluated.
 * Fixed rules have no numerical cutoff: historical NA and current empty input
 * cells are accepted only for those policies, and output omits their cutoff.
 * @author Brian, Ganyu, Raiden
 */
public final class NcrnaGateCutoffSummary{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaGateCutoffSummary.class, false);
		try{
			new NcrnaGateCutoffSummary(pp.args).process();
			pp.outstream.println("NCRNA_GATE_CUTOFF_SUMMARY_COMPLETE_MONITOR_CHOICES_OUTER_CONFIRMATION_ONLY");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Bind the original manifest and completed outer comparison; no implicit file discovery. */
	public NcrnaGateCutoffSummary(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("outer")){outer=b;}
			else{require(parser.parse(arg, a, b), "Unknown saved cutoff argument: "+arg);}
		}
		input=parser.in1;
		output=parser.out1;
		require(input!=null && outer!=null && output!=null && !new File(output).exists()
			&& !parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1, "Use in=original_cells.tsv outer=completed_outer_root out=fresh_directory");
		require(Tools.testInputFiles(false, true, input), "The original complete fit manifest must be readable");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Saved Comparisons      ----------------*/
	/*--------------------------------------------------------------*/

	/** Read every declared regression seed, then select independently of outer values. */
	public void process(){
		final File root=new File(input).getAbsoluteFile().getParentFile();
		require(one(new File(outer, "complete.txt")).equals("THIRTY_SIX_OUTER_GATE_COMPARISONS_COMPLETE_NO_SELECTION_NO_CALLER"),
			"Outer confirmation requires the completed predeclared comparison");
		final double[][][][][] values=new double[2][8][3][][];
		final LineParser1 p=new LineParser1('\t');
		int count=0;
		try(Input manifest=new Input(input)){
			RrnaResourceIO.header(manifest.reader, "cell\tarm\tend\tseed\toutput\trate\trole\tjob", input);
			for(byte[] line=manifest.reader.nextLine(); line!=null; line=manifest.reader.nextLine()){
				p.set(line);
				require(p.terms()==8 && (p.termEquals("R", 1) || p.termEquals("T", 1)), "Malformed original fit manifest");
				if(p.termEquals("T", 1)){continue;}
				final String cell=p.parseString(0);
				final int end=index(p.parseString(2), ENDS);
				final int activation=index(p.parseString(4), ACTIVATIONS);
				final int rate=index(p.parseString(5), RATES);
				final int seed=p.parseInt(3)-58;
				final int group=end*4+activation*2+rate;
				require(cell.matches("[A-Za-z0-9_.]+") && seed>=0 && seed<3 && values[0][group][seed]==null,
					"Each end/output/rate requires exactly the three distinct seeds58,59,60");
				final File fit=new File(root, "runs/"+cell);
				final File confirmation=new File(outer, "runs/"+cell);
				require(one(new File(fit, "exit.txt")).equals("0") && one(new File(confirmation, "exit.txt")).equals("0")
					&& one(new File(fit, "complete.txt")).equals("GATE_CONVERGENCE_FIT_AND_SAVED_REPLAY_COMPLETE_OUTER_UNPERFORMED")
					&& one(new File(confirmation, "complete.txt")).equals("OUTER_GATE_SCORED_AND_REPLAYED_NO_SELECTION_NO_CALLER"),
					"Every input must have its own successful complete receipt: "+cell);
				require(one(new File(fit, "source.txt")).equals(one(new File(confirmation, "fit.source.txt")))
					&& one(new File(fit, "build.seal")).equals(one(new File(confirmation, "fit.build.seal")))
					&& one(new File(fit, "job.txt")).equals(one(new File(confirmation, "fit.job.txt")))
					&& one(new File(confirmation, "selection.txt")).equals("monitor_mse"), "Outer rows must confirm the same selected regression fit");
				values[0][group][seed]=read(new File(fit, "replay/summary.tsv"), end, false);
				values[1][group][seed]=read(new File(confirmation, "replay/summary.tsv"), end, true);
				count++;
			}
		}
		require(count==24, "The summary requires all24regression fits, never a selected outer winner");
		require(new File(output).mkdir(), "Summary output must be fresh");
		try(Output grid=new Output(path("monitor_grid.tsv")); Output choices=new Output(path("choices.tsv"))){
			grid.writer.println(OUT_HEADER);
			choices.writer.println(OUT_HEADER);
			for(int group=0; group<8; group++){
				final double[][] monitor=medians(values[0][group]);
				// Outer is never an argument to the selection predicates below.
				final double[][] confirmation=medians(values[1][group]);
				for(int policy=0; policy<15; policy++){write(grid, group, "monitor", "measured_grid", policy, monitor);}
				for(int policy=0; policy<4; policy++){
					write(choices, group, "monitor", "fixed_rule", policy, monitor);
					write(choices, group, "outer", "fixed_rule", policy, confirmation);
				}
				for(int objective=0; objective<3; objective++){
					final int metric=new int[]{4, 11, 7}[objective];
					double best=monitor[4][metric];
					for(int policy=5; policy<15; policy++){
						best=objective==2 ? Math.min(best, monitor[policy][metric]) : Math.max(best, monitor[policy][metric]);
					}
					for(int policy=4; policy<15; policy++){
						if(monitor[policy][metric]!=best){continue;}
						write(choices, group, "monitor", OBJECTIVES[objective], policy, monitor);
						write(choices, group, "outer", OBJECTIVES[objective], policy, confirmation);
					}
				}
			}
		}
		try(Output summary=new Output(path("summary.tsv"))){
			summary.writer.println("models\tgroups\tseeds_per_group\tcutoffs_per_group\tqualification");
			summary.writer.println("24\t8\t3\t11\tCOMPONENTWISE_MEDIANS_NATIVE_TABLE_PRECISION_ALL_MONITOR_OPTIMUM_TIES_OUTER_CONFIRMATION_ONLY");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Table Validation      ----------------*/
	/*--------------------------------------------------------------*/

	/** Fixed schemas and complete policy grids prevent silently combining different assays. */
	private static double[][] read(final File file, final int end, final boolean outer){
		final double[][] result=new double[15][];
		final boolean[] seen=new boolean[26];
		final LineParser1 p=new LineParser1('\t');
		final long trials=outer ? (end==0 ? 13401 : 13406) : (end==0 ? 88966 : 88979);
		final long proposals=outer ? (end==0 ? 7616 : 9347) : (end==0 ? 50183 : 61986);
		int count=0;
		try(Input input=new Input(file.getPath())){
			RrnaResourceIO.header(input.reader, TABLE_HEADER, file.getPath());
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				require(p.terms()==40 && p.termEquals(ENDS[end], 2) && p.parseLong(3)==trials && p.parseLong(4)==proposals,
					"Saved cutoff row has the wrong schema/end/population: "+file);
				final String policy=p.parseString(0);
				int slot=-1;
				for(int i=0; i<4; i++){if(policy.equals(FIXED[i])){slot=i;}}
				if(slot>=0){require(p.length(1)==0 || p.termEquals("NA", 1), "Fixed rules have no gate cutoff");}
				else{
					require(policy.equals("gate") || policy.equals("gate_hard_oor_veto"), "Unknown policy must not silently enter the gate comparison");
					final double threshold=p.parseDouble(1);
					final int tenth=(int)Math.round(threshold*10);
					require(Double.isFinite(threshold) && tenth>=0 && tenth<=10 && Math.abs(threshold-tenth/10.0)<1e-7,
						"Only the predeclared measured coarse cutoff grid0..1 is compared");
					slot=(policy.equals("gate") ? 4 : 15)+tenth;
				}
				require(!seen[slot], "Duplicate policy/cutoff in "+file);
				seen[slot]=true;
				count++;
				if(slot>=15){continue;}
				final double[] row=new double[12];
				for(int i=0; i<COLUMNS.length; i++){
					row[i]=p.parseDouble(COLUMNS[i]);
					require(Double.isFinite(row[i]), "All compared counts/distances/rates must be finite");
				}
				require(row[0]>=0 && row[0]<=proposals && row[1]>=0 && row[2]>=0 && row[3]>=0
					&& row[1]+row[2]+row[3]==row[0] && row[4]>=0 && row[4]<=trials, "Movement counts must conserve closer/worse/tied moves");
				row[11]=row[1]-row[2];
				result[slot]=row;
			}
		}
		require(count==26, "The complete four fixed rules and22gate/veto rows must be present: "+file);
		return result;
	}

	/** Each displayed value is the middle of three seeds, not a synthetic pooled model. */
	private static double[][] medians(final double[][][] seeds){
		require(seeds.length==3 && seeds[0]!=null && seeds[1]!=null && seeds[2]!=null, "All three seeds are required before taking a median");
		final double[][] result=new double[15][12];
		final double[] scratch=new double[3];
		for(int policy=0; policy<15; policy++){
			if(policy<4){require(Arrays.equals(seeds[0][policy], seeds[1][policy]) && Arrays.equals(seeds[0][policy], seeds[2][policy]),
				"Same-population fixed rules must be identical across the three models");}
			for(int metric=0; metric<12; metric++){
				for(int seed=0; seed<3; seed++){scratch[seed]=seeds[seed][policy][metric];}
				Arrays.sort(scratch);
				result[policy][metric]=scratch[1];
			}
		}
		return result;
	}

	/** Compare componentwise medians with fixed-rule extrema at the saved table precision. */
	private static void write(final Output out, final int group, final String cohort, final String objective, final int policy, final double[][] rows){
		final double[] row=rows[policy];
		double exact=rows[0][4], rms=rows[0][7], net=rows[0][11];
		for(int i=1; i<4; i++){
			exact=Math.max(exact, rows[i][4]);
			rms=Math.min(rms, rows[i][7]);
			net=Math.max(net, rows[i][11]);
		}
		final boolean dominates=row[4]>=exact && row[7]<=rms;
		final ByteBuilder text=new ByteBuilder(ENDS[group/4]).tab().append(ACTIVATIONS[(group%4)/2]).tab().append(RATES[group%2])
			.tab().append(cohort).tab().append(objective).tab().append(policy<4 ? FIXED[policy] : "gate").tab();
		if(policy>=4){text.appendSlow((policy-4)/10.0);}
		for(final double value:row){text.tab().appendSlow(value);}
		text.tab().append(dominates ? 1 : 0).tab().append(dominates && row[11]>=net ? 1 : 0)
			.tab().append(dominates && (row[4]>exact || row[7]<rms) ? 1 : 0)
			.tab().append(Arrays.equals(row, rows[3]) ? 1 : 0).nl();
		out.writer.print(text);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Metadata            ----------------*/
	/*--------------------------------------------------------------*/

	/** Reject unknown settings instead of silently dropping a comparison group. */
	private static int index(final String value, final String[] options){
		for(int i=0; i<options.length; i++){if(options[i].equals(value)){return i;}}
		throw new IllegalArgumentException("Unexpected declared regression setting: "+value);
	}
	/** Require exactly one receipt row; presence alone is not completion. */
	private static String one(final File file){
		try(Input in=new Input(file.getPath())){
			final byte[] line=in.reader.nextLine();
			require(line!=null && in.reader.nextLine()==null, "Expected exactly one complete receipt: "+file);
			return new String(line, java.nio.charset.StandardCharsets.UTF_8);
		}
	}
	private String path(final String name){return new File(output, name).getPath();}
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}

	/*--------------------------------------------------------------*/
	/*----------------       Fields and Schemas     ----------------*/
	/*--------------------------------------------------------------*/

	/** Original fit manifest, fresh report directory and completed outer run root. */
	private String input, output, outer;
	private static final String[] ENDS={"5prime", "3prime"}, ACTIVATIONS={"sigmoid", "rslog"}, RATES={"0.003", "0.0003"};
	private static final String[] FIXED={"do_nothing", "margin0", "margin0.1", "all_proposals"};
	private static final String[] OBJECTIVES={"max_exact", "max_net_improvements", "min_rms"};
	/** Native summary columns; the additional net-count metric is closer minus worse. */
	private static final int[] COLUMNS={5,14,15,36,6,9,11,12,13,39,33};
	public static final String TABLE_HEADER="policy\tthreshold\tend\ttrials\tproposals\tmoved\texact\twithin1\twithin3\texact_pct\twithin3_pct"
		+"\trms_before\trms_after\trms_delta\tbetter\tworse\ttied\tbetter_fraction\tcorrect_before\tcorrect_harmed\treachable_wrong\treachable_fixed"
		+"\toor\toor_improved\toor_worsened\toor_tied\toor_held\tcrossed_proposals\taccepted_crossed\tshort_proposals\taccepted_short\texact_delta\texact_delta_pp"
		+"\tmoved_pct_trials\tbetter_pct_moves\tworse_pct_moves\ttied_moves\ttied_pct_moves\tmissed_gains\tnet_better_minus_worse_per100_trials";
	private static final String OUT_HEADER="end\toutput_activation\tlearning_rate\tcohort\tobjective\tpolicy\tcutoff"
		+"\tmedian_moves\tmedian_closer\tmedian_worse\tmedian_tied_moves\tmedian_exact\tmedian_exact_pct\tmedian_rms_before\tmedian_rms_after"
		+"\tmedian_rms_change\tmedian_net_per100_trials\tmedian_move_pct\tmedian_net_count"
		+"\tnonworse_than_all_fixed_exact_rms\tnonworse_than_all_fixed_exact_rms_net\tstrict_exact_or_rms_gain_over_fixed_extrema\tsame_reported_metrics_as_all_proposals";
}
