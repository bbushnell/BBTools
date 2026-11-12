package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.time.Duration;
import java.time.Instant;
import java.util.ArrayList;

import map.IntLongHashMap;
import map.ObjectMap;
import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Reports saved native search costs; never executes searches or changes truth.
 * Linear projections use each final sample's nucleotide count, including the same
 * two-strand work on both sides. They do not measure per-genome launch/I/O overhead,
 * queue wait, composition effects or parallel efficiency. Expanded timing tables
 * retain failed, short and unstarted points without extrapolating them; a partial
 * roster yields only a qualified-subset projection. One thread owns state.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceCost{

	static void process(final String input, final String census, final String previous, final String prefix){
		final NcrnaReferenceCost report=new NcrnaReferenceCost(input);
		report.load(census);
		if(previous!=null){report.compare(previous);}
		report.write(prefix);
	}
	private NcrnaReferenceCost(final String root_){root=new File(root_);}

	/*--------------------------------------------------------------*/
	/*----------------             Input            ----------------*/
	/*--------------------------------------------------------------*/

	private void load(final String census){
		try(Input in=new Input(path("costs.tsv"))){
			broad=new String(required(in, "costs.tsv"), StandardCharsets.UTF_8).equals(BROAD_HEADER);
		}
		if(broad){
			final int exit=Integer.parseInt(line(path("exit.txt")));
			require(exit>=0 && line(path("report.exit")).equals("0"), "Expanded timing table requires a terminal pilot and completed telemetry report");
			pilotComplete=exit==0;
			if(pilotComplete){
				require(!new File(root, "STOP.txt").exists() && line(path("result.txt")).equals("BROAD_TIMING_TASKS_FINISHED_CHECK_UNPRICED_ROWS"),
					"Successful expanded pilot requires its completion marker and no stop receipt");
			}
		}else{
			require(line(path("exit.txt")).equals("0") && !new File(root, "STOP.txt").exists(), "Pilot must have completed without a budget stop");
			require(line(path("result.txt")).equals("ADDITIONAL_MODEL_COST_PILOT_COMPLETE"), "Expected completed timing pilot");
		}
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(census)){
			RrnaResourceIO.header(in.reader, "genomes\tinput_file_bytes\tnucleotide_bases\tfasta_records\tboth_strands_search_space_Mb", census);
			p.set(required(in, census));
			require(p.terms()==5 && in.reader.nextLine()==null, "Full census summary requires exactly one five-field row");
			genomes=p.parseLong(0);
			bases=p.parseLong(2);
			require(genomes>0 && bases>0, "Full input population must be positive");
		}
		try(Input in=new Input(path("population.tsv"))){
			RrnaResourceIO.header(in.reader, "stage\tgenomes\tnucleotide_bases\tfasta_records", "population.tsv");
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==4, "Sample population requires four fields");
				final int stage=p.parseInt(0);
				final long sample=p.parseLong(2);
				require(stage>0 && p.parseInt(1)==stage && sample>0 && p.parseLong(3)>0 && !sampleBases.contains(stage),
					"Sample stages must be distinct and retain positive measured dimensions");
				sampleBases.put(stage, sample);
			}
		}
		if(broad){
			loadExpanded();
			return;
		}
		try(Input in=new Input(path("candidates.tsv"))){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==7, "Candidate metadata requires seven fields");
				final String accession=p.parseString(0), name=p.parseString(1);
				require(accession.matches("RF[0-9]{5}") && !name.isEmpty() && !byModel.contains(accession), "Candidate accessions must be unique Rfam IDs");
				final Cost cost=new Cost(accession, name);
				byModel.put(accession, cost);
				costs.add(cost);
			}
		}
		try(Input in=new Input(path("costs.tsv"))){
			RrnaResourceIO.header(in.reader, "model\tstage\tuser_CPU_seconds\tsystem_CPU_seconds\twall_seconds\tmax_RSS_KiB\ttiming_status", "costs.tsv");
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==7, "Native cost rows require seven fields");
				final Cost cost=byModel.get(p.parseString(0));
				require(cost!=null && !cost.seen, "Every candidate must have exactly one final timing row");
				cost.stage=p.parseInt(1);
				require(sampleBases.contains(cost.stage), "Unknown timing sample stage");
				cost.bases=sampleBases.get(cost.stage);
				cost.cpu=nonnegative(p, 2)+nonnegative(p, 3);
				cost.wall=nonnegative(p, 4);
				require(p.termEquals("MEASURED", 6) && cost.cpu>=10 && cost.wall>0, "Do not extrapolate short or incomplete final timings");
				final String point="search_"+cost.stage+'/'+cost.accession+'/';
				require(line(path(point+"exit.txt")).equals("0"), "Final timed native point must have exited successfully");
				cost.bytes=Math.addExact(bytes(point+"hits.tbl"), bytes(point+"search.out"));
				cost.projectedCpu=cost.cpu*bases/cost.bases;
				cost.projectedBytes=cost.bytes*(double)bases/cost.bases;
				cost.seen=true;
			}
		}
		require(!costs.isEmpty(), "Cannot report an empty candidate roster");
		for(final Cost cost:costs){require(cost.seen, "Candidate lacks a completed measured cost: "+cost.accession);}
		for(int i=0; i<PHASES.length; i++){
			p.set(line(path(PHASES[i]+".time")).getBytes(StandardCharsets.US_ASCII));
			require(p.terms()==4, "Phase time requires user/system/wall/RSS fields");
			final double cpu=nonnegative(p, 0)+nonnegative(p, 1);
			pilotCpu+=cpu;
			if(i<2){preparationCpu+=cpu;}
		}
		elapsed=Duration.between(Instant.parse(line(path("start.txt"))), Instant.parse(line(path("end.txt")))).getSeconds();
		require(elapsed>=0, "Pilot end precedes start");
	}

	/** Reads saved final-point telemetry, including explicitly unpriced rows.
	 * The producer's table supplies measurements; this reader does not rerun native
	 * searches or infer completion from file sizes. Reused points are not pilot CPU.
	 * Earlier sample points remain in points.tsv; projections use costs.tsv only.
	 */
	private void loadExpanded(){
		assert(broad) : "Expanded parsing must follow exact eleven-column header recognition";
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path("model_assets.tsv"))){
			RrnaResourceIO.header(in.reader, "accession\trelease\tCLEN\tGA_bits\tsha80\tpath", "model_assets.tsv");
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==6, "Model asset rows require six fields");
				final String accession=p.parseString(0);
				require(accession.matches("RF[0-9]{5}") && !byModel.contains(accession), "Model roster requires unique Rfam accessions");
				require(p.parseInt(2)>0 && nonnegative(p, 3)>0, "Model span and original GA must be positive");
				final Cost cost=new Cost(accession, accession);
				byModel.put(accession, cost);
				costs.add(cost);
			}
		}
		try(Input in=new Input(path("costs.tsv"))){
			RrnaResourceIO.header(in.reader, BROAD_HEADER, "costs.tsv");
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==11, "Expanded cost rows require eleven fields");
				final Cost cost=byModel.get(p.parseString(0));
				require(cost!=null && !cost.seen, "Every model must have exactly one selected timing row");
				cost.stage=p.parseInt(1);
				require(sampleBases.contains(cost.stage), "Unknown timing sample stage");
				cost.bases=sampleBases.get(cost.stage);
				cost.status=p.parseString(6);
				cost.outputStatus=p.parseString(10);
				cost.qualified=p.termEquals("MEASURED", 6) || p.termEquals("MEASURED_REUSED", 6);
				final boolean shortPoint=p.termEquals("UNPRICED_TOO_SHORT", 6);
				cost.hasTiming=cost.qualified || shortPoint;
				require(cost.hasTiming || p.termEquals("UNPRICED_NOT_COMPLETED", 6), "Unknown timing qualification for "+cost.accession);
				cost.hitsBytes=optionalBytes(p, 8);
				cost.reportBytes=optionalBytes(p, 9);
				if(cost.hasTiming){
					cost.cpu=nonnegative(p, 2)+nonnegative(p, 3);
					cost.wall=nonnegative(p, 4);
					cost.rss=p.parseLong(5);
					require(cost.wall>0 && cost.rss>0 && (cost.qualified ? cost.cpu>=10 : cost.cpu<10), "Timing qualification contradicts duration or native RSS for "+cost.accession);
					require(p.termEquals("COMPLETE_NATIVE_OUTPUTS", 10) && cost.hitsBytes>=0 && cost.reportBytes>=0,
						"Successful native points require both measured output sizes");
					cost.bytes=Math.addExact(cost.hitsBytes, cost.reportBytes);
				}else{
					for(int i=2; i<=5; i++){require(p.termEquals("NA", i), "Incomplete timings must remain explicitly unavailable");}
					require(p.termEquals("PARTIAL_NATIVE_OUTPUTS", 10) || p.termEquals("NOT_STARTED", 10), "Incomplete native outputs cannot be labeled complete");
					if(p.termEquals("NOT_STARTED", 10)){
						require(cost.hitsBytes<0 && cost.reportBytes<0, "An unstarted point cannot claim measured file sizes");
					}
				}
				if(cost.qualified){
					cost.projectedCpu=cost.cpu*bases/cost.bases;
					cost.projectedBytes=cost.bytes*(double)bases/cost.bases;
					require(Double.isFinite(cost.projectedCpu) && Double.isFinite(cost.projectedBytes), "Sample projection overflow for "+cost.accession);
				}
				cost.seen=true;
			}
		}
		require(!costs.isEmpty(), "Cannot report an empty model roster");
		for(final Cost cost:costs){require(cost.seen, "Model lacks its selected timing row: "+cost.accession);}
		if(pilotComplete){
			for(final String phase:new String[]{"prepare", "startup", "searches"}){
				p.set(line(path(phase+".time")).getBytes(StandardCharsets.US_ASCII));
				require(p.terms()==4, "Phase time requires user/system/wall/RSS fields");
				final double cpu=nonnegative(p, 0)+nonnegative(p, 1);
				pilotCpu+=cpu;
				if(phase.equals("prepare")){preparationCpu=cpu;}
			}
		}
		elapsed=Duration.between(Instant.parse(line(path("start.txt"))), Instant.parse(line(path("end.txt")))).getSeconds();
		require(elapsed>=0, "Pilot end precedes start");
	}
	/** Missing output is unknown, distinct from an existing zero-byte file. */
	private static long optionalBytes(final LineParser1 p, final int field){
		if(p.termEquals("NA", field)){return -1;}
		final long value=p.parseLong(field);
		require(value>=0, "Recorded native output bytes must be nonnegative or NA");
		return value;
	}

	/** These are logical bytes of the native reports, not filesystem allocation or genomes. */
	private long bytes(final String relative){
		final File file=new File(root, relative);
		require(file.isFile() && file.canRead(), "Missing native output: "+file);
		return file.length();
	}
	private String path(final String relative){return new File(root, relative).getPath();}
	private static byte[] required(final Input in, final String path){
		final byte[] row=in.reader.nextLine();
		require(row!=null, "Missing data row: "+path);
		return row;
	}
	private static String line(final String path){
		try(Input in=new Input(path)){
			final String value=new String(required(in, path), StandardCharsets.UTF_8);
			require(in.reader.nextLine()==null, "Expected exactly one line: "+path);
			return value;
		}
	}
	private static double nonnegative(final LineParser1 p, final int field){
		final double value=p.parseDouble(field);
		require(Double.isFinite(value) && value>=0, "Timing values must be finite and nonnegative");
		return value;
	}

	/** Compares the independently generated table numerically at its declared precision.
	 * The previous report is comparison evidence, never the source of current values.
	 */
	private void compare(final String path){
		require(!broad, "Legacy comparison tables do not describe expanded or partial timing rows");
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, "accession\tname\tsample_genomes\tsample_bases\tCPU_seconds\twall_seconds\tCPU_seconds_per_Gbp\tprojected_full_CPU_hours\tstatus", path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==9, "Previous cost report requires nine fields");
				final Cost cost=byModel.get(p.parseString(0));
				require(cost!=null && !cost.compared && p.termEquals(cost.name, 1)
					&& p.parseInt(2)==cost.stage && p.parseLong(3)==cost.bases
					&& p.termEquals("LINEAR_SAMPLE_PROJECTION_NOT_FULL_SCAN", 8), "Previous report identity or population changed");
				require(equal(p, 4, cost.cpu) && equal(p, 5, cost.wall)
					&& equal(p, 6, cost.cpu*1e9/cost.bases) && equal(p, 7, cost.projectedCpu/3600),
					"Previous report numeric values disagree for "+cost.accession);
				cost.compared=true;
				compared++;
			}
		}
		require(compared==costs.size(), "Previous report must cover every current candidate exactly once");
	}
	private static boolean equal(final LineParser1 p, final int column, final double expected){
		final double actual=p.parseDouble(column);
		return Double.isFinite(actual) && Math.abs(actual-expected)<=0.00000051;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(final String prefix){
		if(broad){
			writeExpanded(prefix);
			return;
		}
		final ByteBuilder b=new ByteBuilder();
		double total=0, viral=0, tbox=0, projectedBytes=0;
		long nativeBytes=0;
		try(Output out=new Output(prefix+SUFFIXES[0])){
			out.writer.println("accession\tname\tsample_genomes\tsample_bases\tCPU_seconds\twall_seconds\tCPU_seconds_per_Gbp\tprojected_full_CPU_hours\tstatus");
			for(final Cost cost:costs){
				b.clear().append(cost.accession).tab().append(cost.name).tab().append(cost.stage).tab().append(cost.bases).tab();
				b.append(cost.cpu, 2).tab().append(cost.wall, 2).tab().append(cost.cpu*1e9/cost.bases, 6).tab();
				out.writer.print(b.append(cost.projectedCpu/3600, 6).append("\tLINEAR_SAMPLE_PROJECTION_NOT_FULL_SCAN\n"));
				total+=cost.projectedCpu;
				if(cost.accession.equals("RF00498") || cost.accession.equals("RF01102")){viral+=cost.projectedCpu;}
				if(cost.accession.equals("RF00230")){tbox=cost.projectedCpu;}
				nativeBytes=Math.addExact(nativeBytes, cost.bytes);
				projectedBytes+=cost.projectedBytes;
			}
		}
		try(Output out=new Output(prefix+SUFFIXES[1])){
			out.writer.println("metric\tvalue\tunits\tinterpretation");
			metric(out, b, "pilot_elapsed", elapsed, "seconds", "MEASURED");
			metric(out, b, "pilot_CPU", pilotCpu, "CPU_seconds", "SUM_OF_NONOVERLAPPING_PHASE_PROCESS_TIMES");
			metric(out, b, "preparation_and_extraction_CPU", preparationCpu, "CPU_seconds", "MEASURED_INCLUDED_IN_PILOT_CPU_ONLY");
			metric(out, b, "all_new_models", costs.size(), "models", "TIMING_CANDIDATES_NOT_APPROVED_ROSTER");
			if(compared>0){metric(out, b, "previous_report_model_rows", compared, "models", "NUMERIC_PARITY_AT_SIX_DECIMAL_PRECISION");}
			metric(out, b, "projected_all_new_models", total/3600, "CPU_hours", "LINEAR_SAMPLE_PROJECTION_NOT_FULL_SCAN");
			metric(out, b, "projected_without_two_viral_candidates", (total-viral)/3600, "CPU_hours", "UNAPPROVED_VIRAL_EXCLUSION_ALTERNATIVE");
			if(byModel.contains("RF00230")){metric(out, b, "projected_T_box", tbox/3600, "CPU_hours", "LINEAR_SAMPLE_PROJECTION_NOT_FULL_SCAN");}
			metric(out, b, "final_sample_native_reports", nativeBytes, "bytes", "MEASURED_HITS_TABLE_PLUS_SEARCH_REPORT");
			metric(out, b, "projected_native_report_bytes", projectedBytes, "bytes", "LINEAR_BASE_SCALING_INCLUDES_SAMPLE_HEADER_BYTES");
			metric(out, b, "projected_native_reports", projectedBytes/GIB, "GiB", "LINEAR_BASE_SCALING_INCLUDES_SAMPLE_HEADER_BYTES");
			// A deliberately generous allocation scenario: give each genome the entire
			// measured multi-genome report for each model. Not a statistical bound.
			metric(out, b, "per_genome_full_sample_report_allowance", nativeBytes*(double)genomes/GIB,
				"GiB", "PLANNING_ALLOWANCE_NOT_BOUND_EACH_GENOME_GETS_WHOLE_SAMPLE_REPORT");
		}
	}

	/** Emits subset-qualified projections without filling unknown model costs. */
	private void writeExpanded(final String prefix){
		assert(!costs.isEmpty()) : "Roster conservation is checked before any cost output opens";
		final ByteBuilder b=new ByteBuilder();
		int measured=0;
		double total=0, projectedBytes=0;
		long nativeBytes=0, maximumRss=0;
		try(Output out=new Output(prefix+SUFFIXES[0])){
			out.writer.println("accession\tname\tsample_genomes\tsample_bases\tCPU_seconds\twall_seconds\tmax_RSS_KiB\thits_bytes\tsearch_report_bytes\tprojected_full_CPU_hours\tprojected_native_report_bytes\tstatus\toutput_status");
			for(final Cost cost:costs){
				b.clear().append(cost.accession).tab().append(cost.name).tab().append(cost.stage).tab().append(cost.bases).tab();
				if(cost.hasTiming){
					b.append(cost.cpu, 2).tab().append(cost.wall, 2).tab().append(cost.rss).tab();
				}else{b.append("NA\tNA\tNA\t");}
				if(cost.hitsBytes<0){b.append("NA");}else{b.append(cost.hitsBytes);}
				b.tab();
				if(cost.reportBytes<0){b.append("NA");}else{b.append(cost.reportBytes);}
				b.tab();
				if(cost.qualified){
					b.append(cost.projectedCpu/3600, 6).tab().append(cost.projectedBytes, 6).tab();
					measured++;
					total+=cost.projectedCpu;
					projectedBytes+=cost.projectedBytes;
					nativeBytes=Math.addExact(nativeBytes, cost.bytes);
				}else{b.append("NA\tNA\t");}
				maximumRss=Math.max(maximumRss, cost.rss);
				out.writer.print(b.append(cost.status).tab().append(cost.outputStatus).nl());
			}
		}
		final boolean allMeasured=measured==costs.size();
		final String qualification=allMeasured ? "LINEAR_SAMPLE_PROJECTION_NOT_FULL_SCAN" : "MEASURED_SUBSET_ONLY_NOT_WHOLE_ROSTER_COST";
		try(Output out=new Output(prefix+SUFFIXES[1])){
			out.writer.println("metric\tvalue\tunits\tinterpretation");
			metric(out, b, "pilot_elapsed", elapsed, "seconds", "MEASURED_EXCLUDES_QUEUE_WAIT");
			if(pilotComplete){
				metric(out, b, "pilot_CPU", pilotCpu, "CPU_seconds", "NONOVERLAPPING_PHASES_EXCLUDES_REUSED_POINT");
				metric(out, b, "preparation_CPU", preparationCpu, "CPU_seconds", "INCLUDED_IN_PILOT_CPU");
			}else{out.writer.println("pilot_CPU\tNA\tCPU_seconds\tUNAVAILABLE_AFTER_FAILED_OR_STOPPED_PILOT");}
			metric(out, b, "roster_models", costs.size(), "models", "DECLARED_MODEL_ROSTER");
			metric(out, b, "measured_models", measured, "models", "COMPLETE_POINTS_AT_LEAST_TEN_CPU_SECONDS");
			metric(out, b, "unpriced_models", costs.size()-measured, "models", "NO_EXTRAPOLATION");
			metric(out, b, allMeasured ? "projected_all_models" : "projected_measured_subset", total/3600, "CPU_hours", qualification);
			metric(out, b, "qualified_sample_native_reports", nativeBytes, "bytes", "MEASURED_HITS_TABLE_PLUS_SEARCH_REPORT");
			metric(out, b, "projected_native_report_bytes", projectedBytes, "bytes", qualification+"_INCLUDES_SAMPLE_HEADERS");
			metric(out, b, "maximum_final_point_RSS", maximumRss, "KiB", "NATIVE_FINAL_POINTS_ONLY_NOT_JVM_HEAP_REQUIREMENT");
		}
	}
	private static void metric(final Output out, final ByteBuilder b, final String key,
		final double value, final String units, final String interpretation){
		assert(Double.isFinite(value) && value>=0) : "Completed finite timings and positive populations must yield finite nonnegative accounting";
		out.writer.print(b.clear().append(key).tab().append(value, 6).tab().append(units).tab().append(interpretation).nl());
	}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}
	private static final class Cost{
		Cost(final String accession_, final String name_){
			accession=accession_;
			name=name_;
		}
		final String accession, name;
		int stage;
		long bases, bytes, rss, hitsBytes=-1, reportBytes=-1;
		double cpu, wall, projectedCpu, projectedBytes;
		boolean seen, compared, qualified, hasTiming;
		String status, outputStatus;
	}
	private final File root;
	private final IntLongHashMap sampleBases=new IntLongHashMap();
	private final ObjectMap<String,Cost> byModel=new ObjectMap<String,Cost>(String.class, Cost.class);
	private final ArrayList<Cost> costs=new ArrayList<Cost>();
	private long genomes, bases, elapsed;
	private int compared;
	private double pilotCpu, preparationCpu;
	private boolean broad, pilotComplete;
	private static final String BROAD_HEADER="model\tstage\tuser_CPU_seconds\tsystem_CPU_seconds\twall_seconds\tmax_RSS_KiB\ttiming_status\tresult_directory\thits_bytes\tsearch_report_bytes\toutput_status";
	private static final double GIB=1073741824.0;
	private static final String[] PHASES={"prepare", "extract", "startup", "models64", "models256"};
	static final String[] SUFFIXES={".model_costs.tsv", ".cost_summary.tsv"};
}
