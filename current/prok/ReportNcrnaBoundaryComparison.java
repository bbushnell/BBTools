package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Presents saved endpoint and detection tables without calling or matching again.
 * Adapted from template.A_SampleByteFile's argument, ByteFile and owned-output
 * lifecycle. A small explicit manifest selects measured and unmeasured rows.
 * Retained row objects are report results; one thread owns all state.
 * @author Brian Bushnell, Raiden
 */
public final class ReportNcrnaBoundaryComparison{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, ReportNcrnaBoundaryComparison.class, false);
		try{
			new ReportNcrnaBoundaryComparison(pp.args).process();
			pp.outstream.println("SAVED_BOUNDARY_COMPARISON_REPORT_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private ReportNcrnaBoundaryComparison(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("family")){family=b;}
			else if(a.equals("mainmember")){mainMember=b;}
			else if(a.equals("historical")){
				require("t".equals(b) || "f".equals(b), "historical=t|f selects the schema explicitly, never by fallback");
				historical="t".equals(b);
			}
			else{require(parser.parse(arg, a, b), "Unknown saved-report argument: "+arg);}
		}
		manifest=parser.in1;
		root=parser.out1==null ? null : new File(parser.out1).getAbsoluteFile();
		require(manifest!=null && root!=null && family!=null && mainMember!=null
			&& family.matches("[A-Za-z0-9_.:-]+") && mainMember.matches("[A-Za-z0-9_.:-]+")
			&& new File(manifest).isFile() && !root.exists() && root.getParentFile().isDirectory()
			&& parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append
			&& parser.maxReads<0 && parser.samplerate==1,
			"Supply in=manifest.tsv out=fresh_directory family=label mainmember=label without filtering");
		Shared.setThreads(1);
	}

	private void process(){
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		for(final String[] entry:read(manifest, MANIFEST)){
			require(seen.add(entry[0]+"\t"+entry[1]), "Each configuration/panel needs one explicit manifest row");
			final long panel=number(entry[1]);
			require(panel>0 && (entry[2].equals("MEASURED") || entry[2].equals("NOT_RUN")), "Declare measured or explicitly unmeasured panel rows");
			if(entry[2].equals("NOT_RUN")){
				require(entry[3].equals(".") && entry[4].equals("."), "Unmeasured rows must not pretend to have saved output paths");
				rows.add(new Row(entry[0], panel, family));
				rows.add(new Row(entry[0], panel, mainMember));
				continue;
			}
			final String prefix=entry[3], comparison=entry[4];
			final String[] metrics=single(read(prefix+".metrics.tsv", METRICS));
			require(number(metrics[0])==panel && number(metrics[8])==0, "Saved panel size must match and duplicate calls remain explicit failures");
			final ArrayList<String[]> members=read(prefix+".members.tsv", MEMBERS);
			final long[][] losses=new long[2][3];
			if(!comparison.equals(".")){
				final ArrayList<String[]> changes=historical ? read(comparison+".summary.tsv", CHANGES) :
					read(comparison+".summary.tsv", CompareCobalaminBenchmarkGrades.SUMMARY_HEADER, CompareCobalaminBenchmarkGrades.SUMMARY_VERSION);
				for(final String[] locus:read(comparison+".loci.tsv", LOCI)){
					if(locus[6].equals("LOST_GRADED_REFERENCE")){
						final boolean absent=locus[7].equals("NO_OVERLAPPING_CLAN_CALL");
						require(absent || locus[7].equals("OVERLAPPING_CLAN_CALL_PRESENT"), "Lost reference needs its saved physical-presence class");
						losses[0][0]++;
						if(absent){losses[0][1]++;}
						if(locus[8].equals(mainMember)){
							losses[1][0]++;
							if(absent){losses[1][1]++;}
						}
					}else{
						require(locus[6].equals("GAINED_GRADED_REFERENCE"), "Saved comparison locus needs a declared lost/gained state");
						losses[0][2]++;
						if(locus[11].equals(mainMember)){losses[1][2]++;}
					}
				}
				require(losses[0][0]==change(changes, "graded_scored_references", "LOST")
					&& losses[0][2]==change(changes, "graded_scored_references", "GAINED")
					&& losses[0][1]==change(changes, "lost_references_without_any_candidate_overlap", "LOST"),
					"Saved locus categories must conserve the accepted aggregate loss/absence/gain counts");
			}
			for(int side=0; side<2; side++){
				final String cohort=side==0 ? family : mainMember;
				final String[] summary=single(read(prefix+'_'+cohort+".summary.tsv", SUMMARY));
				final String[] near=unique(read(prefix+'_'+cohort+".radii.tsv", RADII), 2, "3");
				require(summary[0].equals(cohort) && near[0].equals(cohort) && summary[1].equals(near[1])
					&& number(summary[2])==number(near[3]), "Endpoint summaries and within-three rows must bind one matching cohort");
				final Row row=new Row(entry[0], panel, cohort);
				row.measured=true;
				row.matched=number(summary[2]);
				row.exact5=number(summary[3]);
				row.exact3=number(summary[4]);
				row.both=number(summary[5]);
				row.near5=number(near[4]);
				row.near3=number(near[5]);
				if(side==0){
					row.calls=number(metrics[3]);
					row.unmatched=number(metrics[5]);
					require(row.matched==number(metrics[4]), "Whole-clan endpoint denominator must equal accepted matched count");
				}else{
					final String[] member=unique(members, 0, mainMember);
					row.calls=number(member[1]);
					row.unmatched=Math.addExact(number(member[4]), number(member[5]));
					require(row.matched==number(member[2]), "Main emitted-member endpoint denominator must equal its saved normalized match count");
				}
				row.lost=losses[side][0];
				row.absent=losses[side][1];
				row.gained=losses[side][2];
				require(row.both<=row.exact5 && row.both<=row.exact3 && row.exact5<=row.near5
					&& row.exact3<=row.near3 && row.near5<=row.matched && row.near3<=row.matched
					&& row.matched+row.unmatched<=row.calls && row.absent<=row.lost,
					"Endpoint and presence counts must conserve their declared matched/call populations");
				rows.add(row);
			}
		}
		require(!rows.isEmpty() && root.mkdir(), "Report needs explicit rows and a fresh output directory");
		try(Output table=new Output(new File(root, "comparison.tsv").getPath());
			Output states=new Output(new File(root, "population.tsv").getPath());
			Output html=new Output(new File(root, "report.html").getPath())){
			table.writer.println(historical ? "#HISTORICAL_NCRNA_BOUNDARY_COMPARISON_V2" : "#NCRNA_BOUNDARY_COMPARISON_V2");
			table.writer.println(OUTPUT);
			states.writer.println("#NCRNA_BOUNDARY_COMPARISON_POPULATION_V2");
			states.writer.println("configuration\tpanel\tcohort\tstate");
			html.writer.print("<!doctype html><html lang='en'><meta charset='utf-8'><title>HELD OUT - endpoint comparison</title><style>body{font:15px system-ui;color:#182435;margin:24px}table{border-collapse:collapse;font-size:13px}th,td{padding:7px;border:1px solid #d6dce4;text-align:right;white-space:nowrap}th{background:#e9eff6}td:first-child,td:nth-child(3){text-align:left}.scroll{overflow:auto}.missing{color:#657184}p{max-width:1250px;line-height:1.45}</style><body><h1>HELD OUT</h1><p>Saved cobalamin endpoint comparisons. Each rate uses that row's matched calls; the main cohort is the currently emitted main member, not independent subtype truth. All measured networks use margin 0.1.</p><div class='scroll'><table><thead><tr>");
			for(final String name:DISPLAY){html.writer.print("<th>"+escape(name)+"</th>");}
			html.writer.print("</tr></thead><tbody>");
			for(final Row row:rows){
				states.writer.print(buffer.clear().append(row.label).tab().append(row.panel).tab().append(row.cohort)
					.tab().append(row.measured ? "MEASURED" : "NOT_RUN").nl());
				if(!row.measured){continue;}// There are no counts for a run that never happened.
				final String[] values=row.values();
				buffer.clear();
				for(int i=0; i<values.length; i++){
					if(i>0){buffer.tab();}
					buffer.append(values[i]);
				}
				table.writer.print(buffer.nl());
				if(row.matched==0){continue;}// Counts stay in the TSV; percentages have no denominator.
				html.writer.print(row.measured ? "<tr>" : "<tr class='missing'>");
				final String[] shown={row.label, ""+row.panel, row.cohort.equals(mainMember) ? "Main member ("+mainMember+")" : "Whole clan", row.measured ? "Measured" : "NOT RUN",
					row.rate(row.both), row.rate(row.exact5), row.rate(row.exact3), row.rate(row.near5), row.rate(row.near3),
					row.count(row.calls), row.count(row.matched), row.count(row.unmatched), row.count(row.absent), row.count(row.lost-row.absent)};
				for(final String value:shown){html.writer.print("<td>"+escape(value)+"</td>");}
				html.writer.print("</tr>");
			}
			html.writer.print("</tbody></table></div>");
			if(historical){html.writer.print("<p><strong>Historical input format selected explicitly.</strong> Saved measurements were not regenerated.</p>");}
			for(final Row row:rows){
				if(!row.measured || row.matched==0){
					html.writer.print("<p>"+escape(row.label)+", "+row.panel+" genomes, "+escape(row.cohort)+": "
						+(row.measured ? row.calls+" calls, zero matched references; endpoint percentages do not exist." : "not run; no measurements exist.")+"</p>");
				}
			}
			html.writer.print("<p><strong>Limits:</strong> panel genomes were excluded from tables and training only; the library and method choices saw these panels; these one-hot networks were trained without leave-one-out.</p><p>Absent and overlap-present losses refer to baseline-detected scored references lost under the accepted grader; absence means no overlapping clan call on the same input, contig and strand. Main-cohort losses use the OFF emitted label; gains use the ON label. Physical overlap is descriptive, not proof of candidate ancestry. B main-only controls were not run. Unmatched calls are not established biological false positives. The historical 300-genome panel uses its stated legacy reference view.</p><p><a href='comparison.tsv'>Measured counts</a>; <a href='population.tsv'>all requested configurations and measurement states</a>. No caller, matcher, grader or neural fit was run for this page.</p></body></html>");
		}
	}

	/** Native byte parsing; retain only the small table cells needed by this report. */
	private static ArrayList<String[]> read(final String path, final String header){
		return read(path, header, null);
	}
	private static ArrayList<String[]> read(final String path, final String header, final String version){
		final ArrayList<String[]> rows=new ArrayList<String[]>();
		final LineParser1 parser=new LineParser1('\t');
		try(Input input=new Input(path)){
			if(version!=null){RrnaResourceIO.header(input.reader, version, path+" (old tables require explicit historical=t)");}
			require(Arrays.equals(input.reader.nextLine(), header.getBytes(StandardCharsets.UTF_8)), "Unexpected saved table schema: "+path);
			parser.set(header.getBytes(StandardCharsets.UTF_8));
			final int width=parser.terms();
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				parser.set(line);
				require(parser.terms()==width, "Saved report row has the wrong column count: "+path);
				final String[] fields=new String[width];
				for(int i=0; i<width; i++){fields[i]=parser.parseString(i);}
				rows.add(fields);
			}
		}
		return rows;
	}
	private long change(final ArrayList<String[]> rows, final String measure, final String state){
		if(historical){return number(unique(rows, 0, measure)[state.equals("LOST") ? 2 : 3]);}
		String[] found=null;
		for(final String[] row:rows){
			if(row[0].equals(measure) && row[1].equals(state)){
				require(found==null, "Duplicate comparison measure/state: "+measure+"/"+state);
				found=row;
			}
		}
		require(found!=null, "Missing comparison measure/state: "+measure+"/"+state);
		return number(found[2]);
	}
	private static String[] single(final ArrayList<String[]> rows){
		require(rows.size()==1, "Saved metric/summary must have one row");
		return rows.get(0);
	}
	private static String[] unique(final ArrayList<String[]> rows, final int column, final String value){
		String[] found=null;
		for(final String[] row:rows){
			if(row[column].equals(value)){
				require(found==null, "Duplicate saved row: "+value);
				found=row;
			}
		}
		require(found!=null, "Missing saved row: "+value);
		return found;
	}
	private static long number(final String value){
		final long number=Long.parseLong(value);
		require(number>=0, "Saved report counts cannot be negative");
		return number;
	}
	private static String escape(final String value){return value.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace("\"", "&quot;");}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}

	private static final class Row{
		Row(final String label_, final long panel_, final String cohort_){
			label=label_;
			panel=panel_;
			cohort=cohort_;
		}
		String count(final long value){
			require(measured, "Unmeasured rows have no numeric counts and belong only in the population table");
			return ""+value;
		}
		String rate(final long value){
			require(measured && matched>0, "Endpoint rates require a measured nonempty matched cohort");
			return String.format(Locale.ROOT, "%.2f%% (%d/%d)", 100.0*value/matched, value, matched);
		}
		String[] values(){return new String[]{label, ""+panel, cohort, measured ? "MEASURED" : "NOT_RUN", count(calls), count(matched), count(unmatched),
			count(both), count(exact5), count(exact3), count(near5), count(near3), count(absent), count(lost-absent), count(lost), count(gained)};}
		final String label, cohort;
		final long panel;
		boolean measured;
		long calls, matched, unmatched, both, exact5, exact3, near5, near3, absent, lost, gained;
	}

	private String family, mainMember, manifest;
	private boolean historical;
	private File root;
	private final ArrayList<Row> rows=new ArrayList<Row>();
	private final ByteBuilder buffer=new ByteBuilder();
	private static final String MANIFEST="configuration\tpanel\tstate\tprefix\tcomparison";
	private static final String METRICS="inputs\tscored_references\tmasked_references\tcalls\tmatched\tunmatched\tmissed\tmasked_calls\tduplicate_queries\texact5\texact3\tboth_exact\tprecision\trecall\tqualification\tcpu_user_seconds\tcpu_system_seconds\twall_seconds\talign_rf00174\talign_rf01482\talign_rf01689";
	private static final String MEMBERS="emitted_member\tcalls\tmatched\tmasked\tunmatched_duplicate\tunmatched_other\tqualification";
	private static final String SUMMARY="family\toverlap\tmatched\texact5\texact3\tboth_exact";
	private static final String RADII="family\toverlap\tradius_bp\tmatched\treachable5\treachable3\treachable_both";
	private static final String CHANGES="measure\tretained\tlost\tgained\tqualification";
	private static final String LOCI="input_id\treference_id\tseqid\tstart\tstop\tstrand\tchange\tpresence_class\tbaseline_member\tbaseline_model\tbaseline_model_in_subset\tcandidate_member\tcandidate_model";
	private static final String OUTPUT="configuration\tpanel\tcohort\tstate\tcalls\tmatched\tunmatched\tboth_exact\texact5\texact3\twithin3_5\twithin3_3\tabsent_baseline_matches\toverlap_present_grading_losses\tlost_baseline_matches\tgained_baseline_matches";
	private static final String[] DISPLAY={"Configuration", "Genomes", "Cohort", "State", "Both exact", "5-prime exact", "3-prime exact", "5-prime within3", "3-prime within3", "Calls", "Matched", "Unmatched", "Absent baseline matches", "Overlap-present losses"};
}
