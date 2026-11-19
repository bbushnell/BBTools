package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.nio.charset.StandardCharsets;
import java.util.Locale;

import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Pools saved conditional gate summaries with the unchanged third member.
 * Exact integer squared sums are added before RMS is computed. No calling,
 * matching, model inference or cutoff selection occurs. Every compatible tied
 * member-cutoff combination survives with its full identity.
 * @author Brian, Ganyu
 */
public final class NcrnaCenteredGatePool{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaCenteredGatePool.class, false);
		try{
			new NcrnaCenteredGatePool(pp.args).process();
			pp.outstream.println("NCRNA_CENTERED_GATE_POOL_COMPLETE_EXACT_SUMS_NO_INFERENCE");
		}finally{Shared.closeStream(pp.outstream);}
	}

	public NcrnaCenteredGatePool(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			require(parser.parse(arg, a, b), "Unknown saved gate-pooling option: "+arg);
		}
		manifest=parser.in1;
		witnessManifest=parser.in2;
		output=parser.out1;
		require(manifest!=null && witnessManifest!=null && output!=null && parser.out2==null
			&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1
			&& !new File(output).exists() && Tools.testInputFiles(false, true, manifest, witnessManifest),
			"Bind complete summary and raw-witness manifests with in=/in2= and a fresh out= directory");
	}

	/*--------------------------------------------------------------*/
	/*----------------          Clan Report         ----------------*/
	/*--------------------------------------------------------------*/

	/** Validates all saved inputs before pooling each compatible cutoff tuple.
	 * No alternative is selected or collapsed because its metrics tie another. */
	public void process(){
		readWitnesses();
		readSummaries();
		require(new File(output).mkdir(), "Cannot create fresh clan-report directory");
		final ByteBuilder row=new ByteBuilder();
		try(Output clan=new Output(path("clan.tsv")); Output accounting=new Output(path("accounting.tsv"));
				Output raw=new Output(path("third_raw.tsv")); Output missing=new Output(path("unavailable.tsv"))){
			clan.writer.println(POOL_HEADER);
			accounting.writer.println("panel\tseed\tmain_rows\tsecond_rows\tclan_rows\tmain_n\tsecond_n\tthird_n\tclan_n");
			raw.writer.println("panel\tmember\tend\tmatched\texact\twithin3\twithin10\tboth_exact\tboth_within3\tboth_within10\traw_sum_squared_bp2\traw_near10_n\tnear_raw_sum_squared_bp2\tqualification");
			missing.writer.println("panel\tmember\tseed\tpolicy\tend\tside\tmonitor_target_moves\tqualification");
			for(final Panel panel:panels){
				for(int end=0; end<2; end++){
					final Stats third=panel.raw[2][end];
					row.clear().append(panel.name).append("\tRF01689\t").append(ENDS[end]);
					for(int column=C_MATCHED; column<=C_BOTH_WITHIN10; column++){
						row.tab().append(third.counts[column]);
					}
					raw.writer.print(row.tab().append(third.sums[S_RAW]).tab().append(third.counts[C_NEAR]).tab().append(third.sums[S_NEAR_RAW])
						.append("\tUNCHANGED_SAVED_OFF_MEMBER_NO_NEW_REFINEMENT\n"));
				}
				for(int seed=58; seed<=60; seed++){
					final Cell cell=panel.cells[seed-58];
					require(cell.present[0] && cell.present[1], "Every matching seed needs both complete member summaries");
					final int[] used=new int[cell.second.size()];
					long combinations=0;
					for(final Row a:cell.main){
						int found=0;
						for(int i=0; i<cell.second.size(); i++){
							final Row b=cell.second.get(i);
							if(!a.sameGroup(b)){continue;}
							final Stats pooled=Stats.pool(a.stats, b.stats, panel.raw[2][a.end], a.joint);
							row.clear().append(panel.name).tab().append(seed).tab().append(a.meta[3]).tab().append(a.meta[8]).tab().append(ENDS[a.end]);
							for(int k=4; k<=7; k++){
								row.tab().append(a.meta[k]);
							}
							for(int k=4; k<=7; k++){
								row.tab().append(b.meta[k]);
							}
							pooled.write(row, a.joint);
							row.tab().append(a.meta[10]).tab().append(a.meta[11]).tab().append(b.meta[10]).tab().append(b.meta[11]);
							clan.writer.print(row.nl());
							used[i]++;
							found++;
							combinations++;
						}
						require(found>0, "A main policy/scope/end has no second-member counterpart; cannot silently omit it");
					}
					for(final int count:used){require(count>0, "Every tied second-member row must enter the clan report");}
					accounting.writer.print(row.clear().append(panel.name).tab().append(seed).tab().append(cell.main.size()).tab().append(cell.second.size())
						.tab().append(combinations).tab().append(panel.expected[0]).tab().append(panel.expected[1]).tab().append(panel.expected[2])
						.tab().append(panel.expected[0]+panel.expected[1]+panel.expected[2]).nl());
					for(final String absent:cell.unavailable){missing.writer.println(panel.name+'\t'+absent);}
				}
			}
		}
		try(Output columns=new Output(path("columns.tsv"))){
			columns.writer.println("column_group\tscope\tmeaning");
			columns.writer.println("matched/exact/within/both/better/worse/tied/moved\tALL_THREE_MEMBERS\tCounts include every unchanged third-member locus; single-end rows have no both-end metrics.");
			columns.writer.println("raw_sum_squared_bp2/sum_squared_bp2\tALL_THREE_MEMBERS\tUnrounded integer error squares are added before RMS; rounded member RMS values are never averaged.");
			columns.writer.println("raw_near10_n/near_*\tFIXED_RAW_SUBSET_ALL_THREE_MEMBERS\tThe original absolute error is at most ten bases; the same loci remain in before/after sums.");
			columns.writer.println("refined_members_*\tRF00174_AND_RF01482_ONLY\tSupport, fixed-end, feature and refinement-geometry diagnostics do not assess the unchanged third member, whose witness file lacks those fields.");
			columns.writer.println("main_*/second_*\tFROZEN_MEMBER_CHOICES\tBoth members' cutoff sides and values are retained for every tied combination; no outer selection is performed.");
			columns.writer.println("qualification\tCONDITIONAL_ENDPOINT_REPLAY\tSaved matched calls do not measure changed detection, path selection or release adoption.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Saved Input Checks     ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads the two explicit source panels and independently counts all three
	 * emitted-member populations, including unchanged third-member errors. */
	private void readWitnesses(){
		try(Input input=new Input(witnessManifest)){
			RrnaResourceIO.header(input.reader, WITNESS_MANIFEST_HEADER, witnessManifest);
			final LineParser1 p=new LineParser1('\t');
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==5, "Raw witness manifest requires panel/path and all three accepted member counts");
				final int panelIndex=panelIndex(p.parseString(0));
				require(panels[panelIndex]==null && p.termStartsWith("/", 1), "Exactly one absolute raw witness source per panel is required");
				final Panel panel=new Panel(PANELS[panelIndex]);
				for(int m=0; m<3; m++){
					panel.expected[m]=p.parseLong(2+m);
					require(panel.expected[m]>0, "All three raw member counts must be explicit and nonempty");
				}
				readWitnessFile(panel, p.parseString(1));
				panels[panelIndex]=panel;
			}
		}
		for(final Panel panel:panels){require(panel!=null, "Both original panels need raw witness evidence");}
	}

	/** Preserves the original one-to-one witness and strand contract. The error
	 * columns already have gene-sense orientation; coordinates only verify them. */
	private void readWitnessFile(final Panel panel, final String file){
		final ObjectSet<String> references=new ObjectSet<String>(String.class);
		final ObjectSet<String> calls=new ObjectSet<String>(String.class);
		final ByteBuilder key=new ByteBuilder();
		try(Input in=new Input(file)){
			RrnaResourceIO.header(in.reader, WITNESS_HEADER, file);
			final LineParser1 p=new LineParser1('\t');
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==12, "Saved witnesses require explicit member, coordinates, strand and both signed errors");
				final int m=memberIndex(p.parseString(2));
				require(p.termStartsWith("cobalamin_"+MEMBERS[m]+"_", 3), "Partition by emitted member, never a potentially misleading reference-name prefix");
				final String source=p.parseString(0);
				require(references.add(source+'\t'+p.parseString(1)), "One source-qualified reference cannot be counted twice");
				final long cs=p.parseInt(5), ce=p.parseInt(6), rs=p.parseInt(7), re=p.parseInt(8);
				require(cs>=1 && ce>=cs && rs>=1 && re>=rs && (p.termEquals("+", 9) || p.termEquals("-", 9)), "Witnesses retain original one-based valid coordinate pairs and explicit strand");
				require(calls.add(key.clear().append(source).tab().append(p.parseString(4)).tab().append(p.parseString(9))
					.tab().append(cs).tab().append(ce).toString()), "A physical saved call cannot consume two references");
				final long d5=p.parseInt(10), d3=p.parseInt(11);
				require(p.termEquals("+", 9) ? d5==cs-rs && d3==ce-re : d5==re-ce && d3==rs-cs,
					"Saved d5/d3 are already gene-sense errors; they must agree with strand and coordinates without a second remapping");
				final long overlap=Math.max(0, Math.min(ce, re)-Math.max(cs, rs)+1);
				require(2*overlap>=Math.max(ce-cs+1, re-rs+1), "The original witness contract uses reciprocal-half one-to-one matches (NcrnaModelEndpointOffsets.readRow)");
				panel.raw[m][0].addRaw(d5, d3);
				panel.raw[m][1].addRaw(d3, d5);
			}
		}
		for(int m=0; m<3; m++){
			require(panel.raw[m][0].counts[C_MATCHED]==panel.expected[m], "Raw witness count differs from the accepted emitted-member cohort");
		}
	}

	/** Binds each panel/member/initialization once; monitor data cannot enter. */
	private void readSummaries(){
		try(Input in=new Input(manifest)){
			RrnaResourceIO.header(in.reader, MANIFEST_HEADER, manifest);
			final LineParser1 p=new LineParser1('\t');
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==5 && p.termStartsWith("/", 3) && p.termStartsWith("/", 4), "Every saved summary binds panel, seed, member and absolute summary/unavailable paths");
				final Panel panel=panels[panelIndex(p.parseString(0))];
				final int seed=p.parseInt(1), member=memberIndex(p.parseString(2));
				require(seed>=58 && seed<=60 && member<2 && !panel.cells[seed-58].present[member], "Each matching initialization/member summary enters once only");
				final Cell cell=panel.cells[seed-58];
				cell.present[member]=true;
				readSummary(panel, cell, seed, member, p.parseString(3));
				readUnavailable(cell, seed, member, p.parseString(4));
			}
		}
	}

	/** Checks every source row against the independently read raw cohort and
	 * retains all scope/end/cutoff alternatives for the matching member pair. */
	private void readSummary(final Panel panel, final Cell cell, final int seed, final int member, final String file){
		final ArrayList<Row> rows=member==0 ? cell.main : cell.second;
		final ObjectSet<String> identities=new ObjectSet<String>(String.class);
		final boolean[][][] coverage=new boolean[POLICIES.length][2][2];
		try(Input in=new Input(file)){
			final byte[] header=in.reader.nextLine();
			final boolean legacy=Arrays.equals(header, NcrnaCenteredGateReplay.LEGACY_RESULT_HEADER.getBytes(StandardCharsets.US_ASCII));
			require(legacy || Arrays.equals(header, NcrnaCenteredGateReplay.RESULT_HEADER.getBytes(StandardCharsets.US_ASCII)),
				"Saved summary requires its complete legacy or explicit-mode schema");
			final LineParser1 p=new LineParser1('\t');
			final ByteBuilder key=new ByteBuilder();
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==(legacy ? 46 : SUMMARY_COLUMNS) && p.termEquals(MEMBERS[member], 0) && p.parseInt(1)==seed
					&& p.termEquals(MEMBERS[member]+"_"+panel.name, 2), "Summary identity must match its declared original panel/member/seed");
				final Row row=new Row();
				key.clear();
				for(int i=0; i<10; i++){
					row.meta[i]=legacy && i>=4 && i<=7 && p.termEquals("NA", i) ? "" : p.parseString(i);
				}
				final int policy=policyIndex(row.meta[3]);
				require(p.termEquals("SINGLE_OTHER_RAW", 8) || p.termEquals("PAIRED_JOINT", 8), "Do not mix single-end and paired evaluation scopes");
				row.joint=p.termEquals("PAIRED_JOINT", 8);
				row.end=endIndex(row.meta[9]);
				require(p.termEquals("CONDITIONAL_REAL_CALLS_NOT_CALLER_OR_NESTED_MONITOR", QUALIFICATION_COLUMN), "Saved comparisons must retain their conditional, non-caller qualification");
				for(int end=0; end<2; end++){
					if(policy<6 || !row.joint && end!=row.end){
						require(row.meta[4+end].isEmpty() && row.meta[6+end].isEmpty() && (legacy || p.length(46+end)==0),
							"Controls and inactive ends must not acquire invented gate choices");
						row.meta[10+end]="";
					}else{
						final String side=row.meta[4+end];
						require(policy==6 ? side.equals("ALL_TIED_BEST") : side.equals("EXACT") || side.equals("LOW") || side.equals("HIGH"), "Preserve each actual frozen cutoff side without forcing an unavailable bracket");
						final boolean never=legacy ? row.meta[6+end].equals("Infinity") : p.termEquals("NEVER", 46+end);
						row.meta[10+end]=never ? "NEVER" : "THRESHOLD";
						if(never){
							require(legacy || row.meta[6+end].isEmpty(), "Never-move has no numeric cutoff");
							row.meta[6+end]="";
						}else{
							require((legacy || p.termEquals("THRESHOLD", 46+end)) && Float.isFinite(Float.parseFloat(row.meta[6+end])),
								"Threshold mode requires a finite numeric cutoff");
						}
					}
				}
				for(final String field:row.meta){key.append(field).tab();}
				require(identities.add(key.toString()), "Duplicate policy/cutoff/scope/end rows cannot multiply the report");
				row.stats=Stats.read(p, row.joint, legacy);
				final Stats raw=panel.raw[member][row.end];
				require(row.stats.counts[C_MATCHED]==raw.counts[C_MATCHED] && row.stats.counts[C_NEAR]==raw.counts[C_NEAR]
					&& row.stats.sums[S_RAW]==raw.sums[S_RAW] && row.stats.sums[S_NEAR_RAW]==raw.sums[S_NEAR_RAW], "Every policy must retain the exact original raw cohort and raw-near square sums");
				if(policy==0){
					for(final int c:new int[]{11,12,13}){require(row.stats.counts[c]==raw.counts[c], "Raw exact/tolerance counts must agree with the original witnesses");}
					if(row.joint){
						for(int c=C_BOTH_EXACT; c<=C_BOTH_WITHIN10; c++){
							require(row.stats.counts[c]==raw.counts[c], "Raw both-end counts must agree with original error pairs");
						}
					}
					require(row.stats.sums[S_AFTER]==raw.sums[S_RAW] && row.stats.sums[S_NEAR_AFTER]==raw.sums[S_NEAR_RAW] && row.stats.counts[C_MOVED]==0, "The RAW control must not change any endpoint");
				}
				coverage[policy][row.joint ? 1 : 0][row.end]=true;
				rows.add(row);
			}
		}
		for(final boolean[][] scopes:coverage){
			for(final boolean[] ends:scopes){
				for(final boolean seen:ends){
					require(seen, "Every declared policy needs both ends and both scopes");
				}
			}
		}
		for(final Row row:rows){
			if(!row.joint){continue;}
			boolean found=false;
			for(final Row other:rows){
				if(!other.joint || row.end==other.end || !row.samePair(other)){continue;}
				for(int c=C_BOTH_EXACT; c<=C_BOTH_WITHIN10; c++){
					require(row.stats.counts[c]==other.stats.counts[c], "A paired outcome has identical both-end counts in its two end rows");
				}
				found=true;
			}
			require(found, "Each tied joint cutoff pair needs both endpoint rows");
		}
	}

	/** Carries saved unavailable brackets into the report; never invents a
	 * cutoff or substitutes outer movement for the original monitor target. */
	private void readUnavailable(final Cell cell, final int seed, final int member, final String file){
		try(Input in=new Input(file)){
			RrnaResourceIO.header(in.reader, "member\tseed\tpolicy\tend\tside\tmonitor_target_moves\tqualification", file);
			final LineParser1 p=new LineParser1('\t');
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==7 && p.termEquals(MEMBERS[member], 0) && p.parseInt(1)==seed
					&& (p.termEquals("LOW_UNAVAILABLE", 4) || p.termEquals("HIGH_UNAVAILABLE", 4)), "Unavailable movement brackets must retain their member/seed identity");
				require((p.termEquals("RATE_MARGIN0", 2) || p.termEquals("RATE_MARGIN01", 2))
					&& (p.termEquals("5prime", 3) || p.termEquals("3prime", 3)) && p.parseLong(5)>=0
					&& p.termEquals("MONITOR_UNATTAINABLE_NOT_RETUNED_ON_OUTER", 6), "Unavailable rows retain monitor targets and never claim outer retuning");
				cell.unavailable.add(new String(line, java.nio.charset.StandardCharsets.UTF_8));
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Exact Accounting      ----------------*/
	/*--------------------------------------------------------------*/

	/** Sparse count indices match NcrnaCenteredGateReplay.RESULT_HEADER. The
	 * four square-sum slots are separate, exact integer accumulators. */
	private static final class Stats{
		/** Adds one unchanged raw endpoint; the other error supplies pair counts. */
		void addRaw(final long error, final long other){
			final long absolute=Math.abs(error), second=Math.abs(other);
			counts[C_MATCHED]++;
			counts[C_TIED]++;
			if(absolute==0){counts[C_EXACT]++;}
			if(absolute<=3){counts[C_WITHIN3]++;}
			if(absolute<=10){counts[C_WITHIN10]++;}
			if(absolute==0 && second==0){counts[C_BOTH_EXACT]++;}
			if(absolute<=3 && second<=3){counts[C_BOTH_WITHIN3]++;}
			if(absolute<=10 && second<=10){counts[C_BOTH_WITHIN10]++;}
			final long square=Math.multiplyExact(error, error);
			sums[S_RAW]=Math.addExact(sums[S_RAW], square);
			sums[S_AFTER]=sums[S_RAW];
			if(absolute<=10){
				counts[C_NEAR]++;
				sums[S_NEAR_RAW]=Math.addExact(sums[S_NEAR_RAW], square);
				sums[S_NEAR_AFTER]=sums[S_NEAR_RAW];
			}
		}
		/** Imports exact counts/sums and verifies the rounded display fields. */
		static Stats read(final LineParser1 p, final boolean joint, final boolean legacy){
			final Stats s=new Stats();
			for(final int column:COUNT_COLUMNS){
				if(!joint && column>=C_BOTH_EXACT && column<=C_BOTH_WITHIN10){
					require(absent(p, column, legacy), "Single-end rows cannot silently supply paired counts");
					s.counts[column]=-1;
				}else{
					s.counts[column]=p.parseLong(column);
					require(s.counts[column]>=0, "Scientific counts must be nonnegative");
				}
			}
			for(int i=0; i<SUM_SLOTS; i++){
				final double value=p.parseDouble(FIRST_SUM_COLUMN+i);
				require(Double.isFinite(value) && value>=0 && value<0x1p53 && value==Math.rint(value), "Saved integer-error square sums must be exactly represented; rounded RMS is never a pooling input");
				s.sums[i]=(long)value;
			}
			s.validate(joint);
			checkDecimal(p, RAW_RMS_COLUMN, Math.sqrt((double)s.sums[S_RAW]/s.counts[C_MATCHED]));
			checkDecimal(p, AFTER_RMS_COLUMN, Math.sqrt((double)s.sums[S_AFTER]/s.counts[C_MATCHED]));
			checkDecimal(p, RMS_DELTA_COLUMN, Math.sqrt((double)s.sums[S_AFTER]/s.counts[C_MATCHED])-Math.sqrt((double)s.sums[S_RAW]/s.counts[C_MATCHED]));
			if(s.counts[C_NEAR]==0){
				require(absent(p, NEAR_RAW_RMS_COLUMN, legacy) && absent(p, NEAR_AFTER_RMS_COLUMN, legacy), "Empty raw-near subsets have no RMS");
			}else{
				checkDecimal(p, NEAR_RAW_RMS_COLUMN, Math.sqrt((double)s.sums[S_NEAR_RAW]/s.counts[C_NEAR]));
				checkDecimal(p, NEAR_AFTER_RMS_COLUMN, Math.sqrt((double)s.sums[S_NEAR_AFTER]/s.counts[C_NEAR]));
			}
			if(s.counts[C_BETTER]+s.counts[C_WORSE]==0){
				require(absent(p, BETTER_FRACTION_COLUMN, legacy), "A no-change cohort has no better fraction");
			}else{
				checkDecimal(p, BETTER_FRACTION_COLUMN, (double)s.counts[C_BETTER]/(s.counts[C_BETTER]+s.counts[C_WORSE]));
			}
			return s;
		}
		private static void checkDecimal(final LineParser1 p, final int column, final double expected){
			final double actual=p.parseDouble(column);
			require(Double.isFinite(actual) && Math.abs(actual-expected)<=1e-7,
				"Saved scalar fields are rounded to eight decimal places and must agree with their unrounded counts/sums: column="+column);
		}
		/** Adds both refined members and the unchanged third member before RMS. */
		static Stats pool(final Stats a, final Stats b, final Stats third, final boolean joint){
			final Stats result=new Stats();
			for(final int c:COUNT_COLUMNS){
				result.counts[c]=!joint && c>=C_BOTH_EXACT && c<=C_BOTH_WITHIN10 ? -1 : Math.addExact(Math.addExact(a.counts[c], b.counts[c]), third.counts[c]);
			}
			for(int i=0; i<SUM_SLOTS; i++){
				result.sums[i]=Math.addExact(Math.addExact(a.sums[i], b.sums[i]), third.sums[i]);
			}
			result.validate(joint);
			return result;
		}
		/** Enforces complete outcome partitions and nested fixed-cohort subsets. */
		void validate(final boolean joint){
			final long n=counts[C_MATCHED];
			for(final int c:COUNT_COLUMNS){require(counts[c]<=n, "An endpoint count cannot exceed its full cohort");}
			require(n>0 && counts[C_BETTER]+counts[C_WORSE]+counts[C_TIED]==n && counts[C_BETTER]+counts[C_WORSE]<=counts[C_MOVED] && counts[C_MOVED]<=n,
				"Better/worse/tied partition the full cohort; moved also includes equal-distance changes");
			require(counts[C_EXACT]<=counts[C_WITHIN3] && counts[C_WITHIN3]<=counts[C_WITHIN10] && counts[C_WITHIN10]<=n && counts[C_NEAR]<=n
				&& counts[C_EXACT_HARMED]<=counts[C_WORSE] && counts[C_NEAR_HARMED]+counts[C_DISTANT_HARMED]==counts[C_WORSE] && counts[C_DISTANT_IMPROVED]<=counts[C_BETTER], "Exactness, tolerances and raw-near/distant harms must be consistent");
			require(!joint || counts[C_BOTH_EXACT]<=counts[C_BOTH_WITHIN3] && counts[C_BOTH_WITHIN3]<=counts[C_BOTH_WITHIN10] && counts[C_BOTH_WITHIN10]<=n, "Both-end tolerances must be nested");
			require(sums[S_NEAR_RAW]<=sums[S_RAW] && sums[S_NEAR_AFTER]<=sums[S_AFTER] && (counts[C_NEAR]>0 || sums[S_NEAR_RAW]==0 && sums[S_NEAR_AFTER]==0), "Raw-near sums describe a fixed subset of the full cohort");
			require(counts[C_INVALID]==counts[C_CROSSED]+counts[C_OUTSIDE]+counts[C_BELOW_MINIMUM]+counts[C_ABOVE_MAXIMUM], "Invalid-pair reasons must conserve their reported total");
		}
		/** Serializes the replay metric suffix, computing RMS only from exact sums. */
		void write(final ByteBuilder row, final boolean joint){
			for(int c=C_MATCHED; c<=QUALIFICATION_COLUMN; c++){
				row.tab();
				if(c==RAW_RMS_COLUMN || c==AFTER_RMS_COLUMN){
					row.append(Math.sqrt((double)sums[c==RAW_RMS_COLUMN ? S_RAW : S_AFTER]/counts[C_MATCHED]), 8);
				}else if(c==RMS_DELTA_COLUMN){
					row.append(Math.sqrt((double)sums[S_AFTER]/counts[C_MATCHED])-Math.sqrt((double)sums[S_RAW]/counts[C_MATCHED]), 8);
				}else if(c==NEAR_RAW_RMS_COLUMN || c==NEAR_AFTER_RMS_COLUMN){
					if(counts[C_NEAR]>0){row.append(Math.sqrt((double)sums[c==NEAR_RAW_RMS_COLUMN ? S_NEAR_RAW : S_NEAR_AFTER]/counts[C_NEAR]), 8);}
				}else if(c==BETTER_FRACTION_COLUMN){
					if(counts[C_BETTER]+counts[C_WORSE]>0){row.append((double)counts[C_BETTER]/(counts[C_BETTER]+counts[C_WORSE]), 8);}
				}else if(c>=FIRST_SUM_COLUMN && c<FIRST_SUM_COLUMN+SUM_SLOTS){
					row.append(sums[c-FIRST_SUM_COLUMN]);
				}else if(c==QUALIFICATION_COLUMN){
					row.append("CONDITIONAL_CLAN_INCLUDES_UNCHANGED_THIRD_REFINE_DIAGNOSTICS_TWO_MEMBERS_ONLY");
				}else if(joint || c<C_BOTH_EXACT || c>C_BOTH_WITHIN10){row.append(counts[c]);}
			}
		}
		final long[] counts=new long[SUMMARY_COLUMNS];
		final long[] sums=new long[SUM_SLOTS];
	}
	private static final class Row{
		boolean sameGroup(final Row other){return end==other.end && joint==other.joint && meta[3].equals(other.meta[3]);}
		/** Joint partners differ only by endpoint, not policy or cutoff tuple. */
		boolean samePair(final Row other){
			for(int i=3; i<=8; i++){
				if(!meta[i].equals(other.meta[i])){return false;}
			}
			return meta[10].equals(other.meta[10]) && meta[11].equals(other.meta[11]);
		}
		/** Original columns0..9, followed by each end's normalized cutoff mode. */
		final String[] meta=new String[12];
		Stats stats;
		int end;
		boolean joint;
	}
	private static final class Cell{
		final boolean[] present=new boolean[2];
		final ArrayList<Row> main=new ArrayList<Row>(), second=new ArrayList<Row>();
		final ArrayList<String> unavailable=new ArrayList<String>();
	}
	private static final class Panel{
		Panel(final String name_){
			name=name_;
			for(int m=0; m<3; m++){
				for(int e=0; e<2; e++){
					raw[m][e]=new Stats();
				}
			}
		}
		final String name;
		final long[] expected=new long[3];
		final Stats[][] raw=new Stats[3][2];
		final Cell[] cells={new Cell(), new Cell(), new Cell()};
	}
	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	/** Only the old complete summary header permits its historical NA cells. */
	private static boolean absent(final LineParser1 p, final int field, final boolean legacy){
		return legacy ? p.termEquals("NA", field) : p.length(field)==0;
	}

	private static int index(final String value, final String[] options){
		for(int i=0; i<options.length; i++){
			if(options[i].equals(value)){return i;}
		}
		throw new IllegalArgumentException("Unknown reporting identity: "+value);
	}
	private static int panelIndex(final String value){return index(value, PANELS);}
	private static int memberIndex(final String value){return index(value, MEMBERS);}
	private static int endIndex(final String value){return index(value, ENDS);}
	private static int policyIndex(final String value){return index(value, POLICIES);}
	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}
	private String path(final String name){return new File(output, name).getPath();}

	/*--------------------------------------------------------------*/
	/*----------------       Fields and Schema      ----------------*/
	/*--------------------------------------------------------------*/
	private final String manifest, witnessManifest, output;
	private final Panel[] panels=new Panel[2];
	private static final String[] PANELS={"new1000", "old300"};
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final String[] ENDS={"5prime", "3prime"};
	private static final String[] POLICIES={"RAW", "CENTER_ONLY", "CENTER_MARGIN0", "CENTER_MARGIN01", "ALL_WHOLE_PROPOSALS", "IGNORE_HEAD_DIAGNOSTIC", "BEST", "RATE_MARGIN0", "RATE_MARGIN01"};
	/** Sparse indices in NcrnaCenteredGateReplay.RESULT_HEADER. Count slots
	 * retain source column numbers; RMS/fraction/text columns are not counts. */
	private static final int C_MATCHED=10;
	private static final int C_EXACT=11;
	private static final int C_WITHIN3=12;
	private static final int C_WITHIN10=13;
	private static final int C_BOTH_EXACT=14;
	private static final int C_BOTH_WITHIN3=15;
	private static final int C_BOTH_WITHIN10=16;
	private static final int C_BETTER=20;
	private static final int C_WORSE=21;
	private static final int C_TIED=22;
	private static final int C_MOVED=23;
	private static final int C_EXACT_HARMED=24;
	private static final int C_NEAR=25;
	private static final int C_NEAR_HARMED=28;
	private static final int C_DISTANT_IMPROVED=29;
	private static final int C_DISTANT_HARMED=30;
	private static final int C_FIXED=31;
	private static final int C_UNSUPPORTED=32;
	private static final int C_MISSING=33;
	private static final int C_CENTER_HELD=34;
	private static final int C_INVALID=35;
	private static final int C_CROSSED=36;
	private static final int C_OUTSIDE=37;
	private static final int C_BELOW_MINIMUM=38;
	private static final int C_ABOVE_MAXIMUM=39;
	private static final int[] COUNT_COLUMNS={C_MATCHED,C_EXACT,C_WITHIN3,C_WITHIN10,C_BOTH_EXACT,C_BOTH_WITHIN3,C_BOTH_WITHIN10,
		C_BETTER,C_WORSE,C_TIED,C_MOVED,C_EXACT_HARMED,C_NEAR,C_NEAR_HARMED,C_DISTANT_IMPROVED,C_DISTANT_HARMED,
		C_FIXED,C_UNSUPPORTED,C_MISSING,C_CENTER_HELD,C_INVALID,C_CROSSED,C_OUTSIDE,C_BELOW_MINIMUM,C_ABOVE_MAXIMUM};
	/** Exact integer sums, imported from replay columns41..44 in this order. */
	private static final int S_RAW=0;
	private static final int S_AFTER=1;
	private static final int S_NEAR_RAW=2;
	private static final int S_NEAR_AFTER=3;
	private static final int SUM_SLOTS=4;
	private static final int RAW_RMS_COLUMN=17;
	private static final int AFTER_RMS_COLUMN=18;
	private static final int RMS_DELTA_COLUMN=19;
	private static final int NEAR_RAW_RMS_COLUMN=26;
	private static final int NEAR_AFTER_RMS_COLUMN=27;
	private static final int BETTER_FRACTION_COLUMN=40;
	private static final int FIRST_SUM_COLUMN=41;
	private static final int QUALIFICATION_COLUMN=45;
	private static final int SUMMARY_COLUMNS=48;
	public static final String MANIFEST_HEADER="panel\tseed\tmember\tsummary\tunavailable";
	public static final String WITNESS_MANIFEST_HEADER="panel\twitnesses\texpected_main\texpected_second\texpected_third";
	public static final String WITNESS_HEADER="input_id\treference_id\temitted_member\tmodel\tseqid\tcall_start\tcall_stop\tref_start\tref_stop\tstrand\td5\td3";
	public static final String POOL_HEADER="panel\tseed\tpolicy\tscope\tend\tmain_side5\tmain_side3\tmain_cutoff5\tmain_cutoff3\tsecond_side5\tsecond_side3\tsecond_cutoff5\tsecond_cutoff3\t"+
		NcrnaCenteredGateReplay.RESULT_HEADER.substring(NcrnaCenteredGateReplay.RESULT_HEADER.indexOf("matched\t"),
			NcrnaCenteredGateReplay.RESULT_HEADER.indexOf("\tcutoff5_mode"))
		.replace("\tfixed\t", "\trefined_members_fixed\t").replace("\tunsupported\t", "\trefined_members_unsupported\t")
		.replace("\tmissing_features\t", "\trefined_members_missing_features\t").replace("\tcenter_held\t", "\trefined_members_center_held\t")
		.replace("\tinvalid_joint\t", "\trefined_members_invalid_joint\t").replace("\tcrossed\t", "\trefined_members_crossed\t")
		.replace("\toutside_contig\t", "\trefined_members_outside_contig\t").replace("\tbelow_minimum\t", "\trefined_members_below_minimum\t")
		.replace("\tabove_maximum\t", "\trefined_members_above_maximum\t")
		+"\tmain_cutoff5_mode\tmain_cutoff3_mode\tsecond_cutoff5_mode\tsecond_cutoff3_mode";
}
