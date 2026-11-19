package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Presents the completed six-configuration comparison through the existing
 * benchmark comparison launcher. Adapted from template.A_SampleByteFile's Parser,
 * byte-reader and owned-writer pattern. Small parsed table fields are retained
 * report results, not transient per-record parsing objects. No caller, grader,
 * fit, spatial join or configuration selection runs here. One thread owns all
 * report state; input schemas and requested unique rows are checked before output.
 * @author Brian Bushnell, Raiden
 */
final class CobalaminBoundaryReport{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	CobalaminBoundaryReport(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("mode")){require("report".equals(b), "Report mode must be explicit");}
			else if(a.equals("historical")){
				require("t".equals(b) || "f".equals(b), "historical=t|f must be explicit; old schemas are never inferred");
				historical="t".equals(b);
			}
			else{require(parser.parse(arg, a, b), "Unknown boundary-report option: "+arg);}
		}
		tables=parser.in1;
		replay=parser.in2;
		out=parser.out1;
		require(tables!=null && out!=null && new File(tables).isDirectory()
			&& !new File(out).exists() && !parser.overwrite
			&& !parser.append && parser.maxReads<0 && parser.samplerate==1,
			"Supply mode=report in=completed_tables out=fresh_directory without filtering");
		require(historical ? replay!=null && new File(replay).isDirectory() : replay==null,
			"Only explicit historical=t accepts and requires in2=original_replay; current reports never read historical replay tables");
		Shared.setThreads(1);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Report             ----------------*/
	/*--------------------------------------------------------------*/

	void process(){
		final Cell[] cells=new Cell[12];
		final String[][] headlines=new String[12][];
		int index=0;
		for(final String kind:new String[]{"sliding", "positional"}){
			for(final int radius:new int[]{6, 10, 16}){
				for(final int panel:new int[]{1000, 300}){cells[index++]=new Cell(kind, radius, panel);}
			}
		}
		for(int i=0; i<cells.length; i++){headlines[i]=cells[i].headline();}
		final Table[] replays=new Table[historical ? 6 : 0];
		for(int i=0; i<replays.length; i++){
			final Cell cell=cells[2*i];
			replays[i]=new Table(replay+"/"+cell.kind+"_r"+cell.radius+"/comparison.replay_endpoints.tsv", REPLAY);
			require(replays[i].rows.size()==32, "Original replay must retain four margins, four cohorts and both ends");
			for(final String margin:new String[]{"0", "0.005", "0.02", "0.1"}){
				for(final String cohort:COHORTS){
					for(final String end:ENDS){replays[i].unique(new int[]{2, 3, 4}, new String[]{margin, cohort, end});}
				}
			}
		}
		require(index==cells.length && new File(out).mkdir(), "Could not create a complete fresh report directory");
		writeTables(cells, headlines);
		try(Output output=new Output(out+"/report.html")){
			output.writer.print("<!doctype html><html lang='en'><meta charset='utf-8'><title>Cobalamin endpoint comparison</title><style>body{font:16px system-ui;max-width:1500px;margin:30px auto;padding:0 20px;color:#172033}p{line-height:1.5;max-width:1150px}.scroll{overflow:auto}table{border-collapse:collapse;font-size:13px;margin:18px 0}td,th{padding:7px;border-bottom:1px solid #ccd;text-align:left;white-space:nowrap}th{background:#eef2f6}img{max-width:100%;height:auto}summary{padding:12px;cursor:pointer}</style><body><h1>Cobalamin: small-radius endpoint comparison</h1>");
			output.writer.print("<p>All six configurations completed on the same 1,000-genome and historical 300-genome panels. Confidence was selected on 1,000 only and applied unchanged to old 300. All eligible original loci were retained without a reserved holdout; 14 third-member loci had no model-dependent vectors and remain accounted for. These are coverage comparisons, not independent validation.</p><p>Absent clan calls and failure of the matching rule are separate. Sliding six bases at margin 0.1 has zero absent losses and unchanged aggregate counts on both panels, despite one lost and one gained grading match on each. No 10,000-genome result, adoption, or main-member-only experiment is claimed.</p><p>Arrows mean OFF to ON. Exact/within3/within10 use each run's matched denominator. Full paired RMS and better/worse/tied counts use the same reference IDs in both runs, partitioned by the OFF emitted member. Unmatched calls are not established biological false positives.</p>");
			if(historical){output.writer.print("<p><strong>Historical format selected explicitly.</strong> These retained version-1 tables may contain undefined-value placeholders. They are not version-2 measurements.</p>");}
			else{output.writer.print("<p>Version-2 statistics contain only defined measurements. Empty populations remain in the population table; their nonexistent RMS statistics are omitted. Better, worse and tied counts replace the derived better fraction.</p>");}
			begin(output, HEADLINE);
			for(final String[] headline:headlines){row(output, headline);}
			end(output);
			for(final int panel:new int[]{1000, 300}){
				output.writer.print("<h2>"+panel+" genomes</h2><img src='endpoints_"+panel+".png' alt='Exact and near-exact endpoint counts by member and end'>");
			}
			output.writer.print("<h2>Unmatched calls and physical counterparts</h2><p>ON-only uses exact coordinates, ignoring labels/models. Physical overlap means at least one shared base on the same input, contig and strand. Shared OFF means multiple ON calls overlap one OFF call; it is geometric evidence, not splitting ancestry. Full witnesses remain on Dori. The retained first cell has no additional sharing headline: explicitly unavailable, never zero.</p>");
			begin(output, "Configuration\tGenomes\t"+SPATIAL);
			for(final Cell cell:cells){
				if(cell.spatial!=null){for(final String[] record:cell.spatial.rows){row(output, cell.prefix(record));}}
			}
			end(output);
			interpretation(output, cells);
			output.writer.print("<details><summary>All member counts</summary>");
			begin(output, "Configuration\tGenomes\t"+MEMBERS);
			for(final Cell cell:cells){for(final String[] record:cell.members.rows){row(output, cell.prefix(record));}}
			end(output);
			output.writer.print("</details><details><summary>Complete paired endpoint gains and harm</summary>");
			begin(output, "Configuration\tGenomes\t"+pairedHeader());
			for(final Cell cell:cells){for(final String[] record:cell.paired.rows){row(output, cell.prefix(record));}}
			end(output);
			if(!historical){
				output.writer.print("</details><details><summary>Paired populations, including zero</summary>");
				begin(output, "Configuration\tGenomes\t"+CompareCobalaminBenchmarkGrades.POPULATION_HEADER);
				for(final Cell cell:cells){for(final String[] record:cell.population.rows){row(output, cell.prefix(record));}}
				end(output);
				output.writer.print("</details><details><summary>Previously-near subset statistics</summary>");
				begin(output, "Configuration\tGenomes\t"+CompareCobalaminBenchmarkGrades.NEAR_HEADER);
				for(final Cell cell:cells){for(final String[] record:cell.near.rows){row(output, cell.prefix(record));}}
				end(output);
			}
			output.writer.print("</details>");
			if(historical){
				output.writer.print("<h2>Historical saved-score diagnostics</h2><p>These original version-1 replay tables retain their historical definitions and undefined-value notation. They are not new version-2 caller measurements. Fixed captured cohorts; diagnostic only. Margins 0, 0.005, 0.02 and 0.1 were predeclared. The later requested 0.2/0.4 extension is post hoc and is not included here.</p>");
			}
			for(int i=0; i<replays.length; i++){
				output.writer.print("<details><summary>"+cells[2*i].label+" full replay table</summary>");
				begin(output, REPLAY);
				for(final String[] record:replays[i].rows){row(output, record);}
				end(output);
				output.writer.print("</details>");
			}
			output.writer.print("<p>The Java tables accompany this report. Saved-call ledgers remain in cobalamin_gated_saved_audit_20261116_01: the retained first cell under run/, the rest under run_legacy_compatible/. No caller or grader was rerun for presentation. K/placement and three-table features are separate work.</p></body></html>");
		}
	}

	/** Describes the observed mixture without equating physical sharing with ancestry. */
	private void interpretation(final Output output, final Cell[] cells){
		for(final int index:new int[]{7, 9, 10, 11}){
			final Cell cell=cells[index];
			final String[] third=cell.members.unique(0, "RF01689");
			final String[] spatial=cell.spatial.unique(0, "RF01689");
			text(output, cell.label+", "+cell.panel+" genomes: third-member total calls "+arrow(third[1], third[2])
				+", unmatched "+arrow(third[7], third[8])+". Of "+spatial[3]+" ON-only unmatched third-member calls, "
				+spatial[4]+" have no OFF overlap, "+spatial[5]+" overlap an OFF call, and "+spatial[6]
				+" share an OFF span with another ON call.", "p");
		}
		for(final int index:new int[]{6, 8}){
			final Cell cell=cells[index];
			final String[] main=cell.members.unique(0, "RF00174");
			text(output, cell.label+", "+cell.panel+" genomes also changes main-member unmatched calls "
				+arrow(main[7], main[8])+"; the observed increases do not arise exclusively in the third member.", "p");
		}
		text(output, "The positional configurations show mixed geometry: the six- and ten-base old-panel results keep third-member totals unchanged, while the sixteen-base results include shared existing spans and, on 1,000 genomes, new places. These observations establish neither uniform splitting nor a proven main-member-only remedy. No such caller experiment was run.", "p");
	}

	/** Retains one compact table per result type instead of twelve copied OFF tables. */
	private void writeTables(final Cell[] cells, final String[][] headlines){
		try(Output summary=new Output(out+"/comparison.tsv"); Output pairs=new Output(out+"/paired_endpoints.tsv");
				Output near=historical ? null : new Output(out+"/paired_near_endpoints.tsv");
				Output population=historical ? null : new Output(out+"/paired_population.tsv");
				Output members=new Output(out+"/members.tsv"); Output spatial=new Output(out+"/unmatched_spatial.tsv")){
			summary.writer.println(historical ? "#HISTORICAL_COBALAMIN_BOUNDARY_REPORT_V1" : "#COBALAMIN_BOUNDARY_REPORT_V2");
			summary.writer.println(HEADLINE);
			pairs.writer.println(historical ? "#HISTORICAL_COBALAMIN_PAIRED_ENDPOINTS_V1" : CompareCobalaminBenchmarkGrades.PAIRED_VERSION);
			pairs.writer.println("Configuration\tGenomes\t"+pairedHeader());
			if(!historical){
				near.writer.println(CompareCobalaminBenchmarkGrades.NEAR_VERSION);
				near.writer.println("Configuration\tGenomes\t"+CompareCobalaminBenchmarkGrades.NEAR_HEADER);
				population.writer.println(CompareCobalaminBenchmarkGrades.POPULATION_VERSION);
				population.writer.println("Configuration\tGenomes\t"+CompareCobalaminBenchmarkGrades.POPULATION_HEADER);
			}
			members.writer.println("Configuration\tGenomes\t"+MEMBERS);
			spatial.writer.println("Configuration\tGenomes\t"+SPATIAL);
			for(final String[] headline:headlines){tsv(summary, headline);}
			for(final Cell cell:cells){
				for(final String[] record:cell.paired.rows){tsv(pairs, cell.prefix(record));}
				if(!historical){
					for(final String[] record:cell.near.rows){tsv(near, cell.prefix(record));}
					for(final String[] record:cell.population.rows){tsv(population, cell.prefix(record));}
				}
				for(final String[] record:cell.members.rows){tsv(members, cell.prefix(record));}
				if(cell.spatial!=null){for(final String[] record:cell.spatial.rows){tsv(spatial, cell.prefix(record));}}
			}
		}
	}
	private void tsv(final Output output, final String[] fields){
		assert(fields.length>0) : "Every saved report row needs its configuration and measured fields";
		buffer.clear();
		for(int i=0; i<fields.length; i++){
			if(i>0){buffer.tab();}
			buffer.append(fields[i]);
		}
		output.writer.print(buffer.nl());
	}

	/** One completed cell; the constructor validates every requested unique row. */
	private final class Cell{
		Cell(final String kind_, final int radius_, final int panel_){
			kind=kind_;
			radius=radius_;
			panel=panel_;
			label=(kind.equals("sliding") ? "Sliding" : "One-hot")+" / "+radius+" bases";
			final String root=tables+"/"+kind+"_r"+radius+"_n"+panel+"/";
			off=new Table(root+"off.metrics.tsv", METRICS).single();
			on=new Table(root+"on.metrics.tsv", METRICS).single();
			require(number(off[0])==panel && number(on[0])==panel, "Metric input count differs from requested cell");
			changes=historical ? new Table(root+"changes.tsv", CHANGES) :
				new Table(root+"changes.tsv", CompareCobalaminBenchmarkGrades.SUMMARY_HEADER, CompareCobalaminBenchmarkGrades.SUMMARY_VERSION);
			members=new Table(root+"members.tsv", MEMBERS);
			paired=new Table(root+"paired_endpoints.tsv", pairedHeader(), historical ? null : CompareCobalaminBenchmarkGrades.PAIRED_VERSION);
			near=historical ? null : new Table(root+"paired_near_endpoints.tsv", CompareCobalaminBenchmarkGrades.NEAR_HEADER, CompareCobalaminBenchmarkGrades.NEAR_VERSION);
			population=historical ? null : new Table(root+"paired_population.tsv", CompareCobalaminBenchmarkGrades.POPULATION_HEADER, CompareCobalaminBenchmarkGrades.POPULATION_VERSION);
			require(members.rows.size()==4, "Expected complete three-member/union counts");
			for(final String member:new String[]{"RF00174", "RF01482", "RF01689", "ALL"}){members.unique(0, member);}
			if(historical){
				require(paired.rows.size()==8, "Historical paired tables require all eight rows");
				for(final String cohort:COHORTS){for(final String end:ENDS){paired.unique(new int[]{0, 1}, new String[]{cohort, end});}}
			}else{validatePairs();}
			for(int side=0; side<2; side++){
				radii[side]=new Table(root+(side==0 ? "off" : "on")+".cobalamin.radii.tsv", RADII);
				for(final String r:new String[]{"0", "3", "10", "100"}){radii[side].unique(2, r);}
			}
			final File path=new File(root+"unmatched_spatial.tsv");
			spatial=path.isFile() ? new Table(path.getPath(), SPATIAL) : null;
			require(spatial!=null || (kind.equals("sliding") && radius==6 && panel==1000), "Missing sharing statistics outside the explicitly retained first cell");
			if(spatial!=null){
				require(spatial.rows.size()==4, "Expected complete spatial member counts");
				for(final String member:new String[]{"RF00174", "RF01482", "RF01689", "ALL"}){spatial.unique(0, member);}
			}
		}
		private void validatePairs(){
			require(population.rows.size()==8, "Version-2 paired populations require all four cohorts and both ends, including zero");
			int fullRows=0;
			int nearRows=0;
			for(final String cohort:COHORTS){
				for(final String end:ENDS){
					final int[] keys={0, 1};
					final String[] values={cohort, end};
					final String[] counts=population.unique(keys, values);
					final long n=number(counts[3]);
					final long close=number(counts[4]);
					require(number(counts[2])==radius && close<=n, "Paired population radius/subset differs from its cell");
					if(n>0){
						final String[] row=paired.unique(keys, values);
						finite(row, 2, 12);
						require(number(row[2])==radius && number(row[3])==n && number(row[6])+number(row[7])+number(row[8])==n
							&& number(row[9])+number(row[10])==number(row[7]), "Paired error changes must conserve their population");
						fullRows++;
					}
					if(close>0){
						final String[] row=near.unique(keys, values);
						finite(row, 2, 8);
						require(number(row[2])==radius && number(row[3])==n && number(row[4])==close,
							"Near statistics must bind the same full and OFF-selected populations");
						require(number(row[7])<=close && number(row[7])==number(paired.unique(keys, values)[9]),
							"Near-subset worsening must agree with the full table and cannot exceed its subset");
						nearRows++;
					}
				}
			}
			require(paired.rows.size()==fullRows && near.rows.size()==nearRows, "Empty or unknown paired cohorts must not acquire statistics rows");
		}
		private long change(final String measure, final String state){
			if(!historical){return number(changes.unique(new int[]{0, 1}, new String[]{measure, state})[2]);}
			return number(changes.unique(0, measure)[state.equals("RETAINED") ? 1 : state.equals("LOST") ? 2 : 3]);
		}
		String[] headline(){
			final ArrayList<String> result=new ArrayList<String>();
			result.add(label);
			result.add(""+panel);
			result.add(kind.equals("positional") && radius!=6 ? "0.005" : "0.1");
			for(final int field:new int[]{3, 4, 5, 7}){result.add(arrow(off[field], on[field]));}
			final long lost=change("graded_scored_references", "LOST"), absent=change("lost_references_without_any_candidate_overlap", "LOST");
			require(absent<=lost, "Absent losses must be a subset of graded losses before subtraction");
			result.add(""+lost);
			result.add(""+absent);
			result.add(""+(lost-absent));
			result.add(""+change("graded_scored_references", "GAINED"));
			result.add(arrow(off[11], on[11]));
			for(final int field:new int[]{4, 5}){
				final String[] values=new String[2];
				for(int side=0; side<2; side++){
					values[side]=radii[side].unique(2, "0")[field]+" / "+radii[side].unique(2, "3")[field]+" / "+radii[side].unique(2, "10")[field];
				}
				result.add(arrow(values[0], values[1]));
			}
			final String[] third=members.unique(0, "RF01689");
			result.add(arrow(third[1], third[2]));
			return result.toArray(new String[result.size()]);
		}
		String[] prefix(final String[] row){
			final String[] result=new String[row.length+2];
			result[0]=label;
			result[1]=""+panel;
			System.arraycopy(row, 0, result, 2, row.length);
			return result;
		}
		final String kind, label;
		final int radius, panel;
		final String[] off, on;
		final Table changes, members, paired, near, population, spatial;
		final Table[] radii=new Table[2];
	}

	/*--------------------------------------------------------------*/
	/*----------------       Small Text Tables      ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Table{
		Table(final String path, final String expected){
			this(path, expected, null);
		}
		Table(final String path, final String expected, final String version){
			final int width=expected.split("\t", -1).length;
			final LineParser1 parser=new LineParser1('\t');
			try(Input input=new Input(path)){
				if(version!=null){RrnaResourceIO.header(input.reader, version, path+" (old formats require explicit historical=t)");}
				RrnaResourceIO.header(input.reader, expected, path);
				for(byte[] line; (line=input.reader.nextLine())!=null;){
					parser.set(line);
					require(parser.terms()==width, "Changed report row width: "+path);
					final String[] record=new String[width];
					for(int column=0; column<width; column++){record[column]=parser.parseString(column);}
					rows.add(record);
				}
			}
		}
		String[] single(){
			require(rows.size()==1, "Expected exactly one metrics row");
			return rows.get(0);
		}
		String[] unique(final int field, final String key){return unique(new int[]{field}, new String[]{key});}
		String[] unique(final int[] fields, final String[] keys){
			assert(fields.length==keys.length) : "Each required report key needs its schema column";
			String[] found=null;
			for(final String[] row:rows){
				boolean matches=true;
				for(int i=0; i<fields.length; i++){matches&=row[fields[i]].equals(keys[i]);}
				if(matches){
					require(found==null, "Duplicate requested report row: "+keys[0]);
					found=row;
				}
			}
			require(found!=null, "Missing requested report row: "+keys[0]);
			return found;
		}
		final ArrayList<String[]> rows=new ArrayList<String[]>();
	}
	private void begin(final Output output, final String header){
		output.writer.print("<div class='scroll'><table><thead><tr>");
		for(final String field:header.split("\t", -1)){text(output, field, "th");}
		output.writer.print("</tr></thead><tbody>");
	}
	private void row(final Output output, final String[] fields){
		output.writer.print("<tr>");
		for(final String field:fields){text(output, field, "td");}
		output.writer.print("</tr>");
	}
	private void text(final Output output, final String text, final String tag){
		buffer.clear().append('<').append(tag).append('>');
		for(int i=0; i<text.length(); i++){
			final char c=text.charAt(i);
			if(c=='&'){buffer.append("&amp;");}
			else if(c=='<'){buffer.append("&lt;");}
			else if(c=='>'){buffer.append("&gt;");}
			else if(c=='\"'){buffer.append("&quot;");}
			else if(c>127){buffer.append("&#").append((int)c).append(';');}
			else{buffer.append(c);}
		}
		output.writer.print(buffer.append("</").append(tag).append('>'));
	}
	private static void end(final Output output){output.writer.print("</tbody></table></div>");}
	private static String arrow(final String before, final String after){return before+" -> "+after;}
	private static long number(final String value){
		final long n=Long.parseLong(value);
		require(n>=0, "Negative report count");
		return n;
	}
	private static void finite(final String[] row, final int from, final int to){
		for(int i=from; i<to; i++){
			final double value=Double.parseDouble(row[i]);
			require(Double.isFinite(value) && value>=0, "Defined endpoint statistics must be finite and nonnegative");
		}
	}
	private String pairedHeader(){return historical ? PAIRED : CompareCobalaminBenchmarkGrades.PAIRED_HEADER;}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	private final String tables, replay, out;
	private boolean historical;
	private final ByteBuilder buffer=new ByteBuilder();
	private static final String[] COHORTS={"cobalamin", "RF00174", "RF01482", "RF01689"}, ENDS={"5prime", "3prime"};
	private static final String METRICS="inputs\tscored_references\tmasked_references\tcalls\tmatched\tunmatched\tmissed\tmasked_calls\tduplicate_queries\texact5\texact3\tboth_exact\tprecision\trecall\tqualification\tcpu_user_seconds\tcpu_system_seconds\twall_seconds\talign_rf00174\talign_rf01482\talign_rf01689";
	private static final String CHANGES="measure\tretained\tlost\tgained\tqualification";
	private static final String MEMBERS="member\toff_calls\ton_calls\toff_matched\ton_matched\toff_masked\ton_masked\toff_unmatched\ton_unmatched\tretained_unmatched_off\tretained_unmatched_on\ton_only_unmatched\toff_only_unmatched\tlabel_switch_in\tlabel_switch_out\tunmatched_net_change";
	private static final String PAIRED="cohort\tend\tradius_bp\tpaired_n\toff_rms_nt\ton_rms_nt\tbetter\tworse\tunchanged\tbetter_fraction\tworse_off_within_radius\tworse_off_beyond_radius\toff_within_radius_n\toff_within_radius_rms_nt\ton_same_subset_rms_nt\tlabel_switches\tqualification";
	private static final String RADII="family\toverlap\tradius_bp\tmatched\treachable5\treachable3\treachable_both";
	private static final String SPATIAL="member\ton_unmatched\tretained_exact_unmatched\ton_only_unmatched\ton_only_no_off_overlap\ton_only_off_overlap\ton_only_shares_off_with_another_on";
	private static final String REPLAY="kind\tradius_bp\tmargin\toff_member_cohort\tend\tpaired_total\treplayable_n\tunknown_n\toff_rms_nt\treplay_rms_nt\tbetter\tworse\tunchanged\tbetter_fraction\texact\twithin3\twithin10\tmoved_from_anchor\tunmoved_from_anchor\tworse_off_within_radius\tworse_off_beyond_radius\toff_within_radius_n\toff_within_radius_rms_nt\treplay_same_subset_rms_nt\tqualification";
	private static final String HEADLINE="Configuration\tGenomes\tMargin\tCalls OFF -> ON\tMatched\tUnmatched\tMasked\tGraded lost\tAbsent lost\tOverlap-present lost\tGained\tBoth exact\t5prime exact / within3 / within10\t3prime exact / within3 / within10\tThird-member calls";
}
