package prok;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import map.ObjectIntMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Compares completed native Infernal timing points with their original GA controls.
 * Derived from template.A_SampleByteFile. Streams one model at a time and retains
 * only distinct included control-hit keys. A key is target, alignment mode,
 * model/target bounds, strand and truncation; duplicate multiplicity is preserved.
 * This compares native reports, not biological loci or complete reference truth.
 * Repeated target names in concatenated genomes cannot identify individual genomes.
 * All state belongs to the calling thread; output owners join their writers.
 * @author Brian Bushnell, Raiden
 */
public final class CompareCmSearchCosts{

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final CompareCmSearchCosts tool=new CompareCmSearchCosts(args);
		try{
			tool.process();
			timer.stop();
			tool.outstream.println("Models compared: "+tool.models+"; time: "+timer);
		}finally{
			Shared.closeStream(tool.outstream);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses native flag=value arguments and requires fresh output files. */
	public CompareCmSearchCosts(String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		args=pp.args;
		outstream=pp.outstream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase();
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("census")){
				census=b;
			}else if(a.equals("population")){
				population=b;
			}else if(a.equals("stage")){
				stage=Integer.parseInt(b);
			}else if(parser.parse(arg, a, b)){
				// Standard BBTools argument.
			}else{
				throw new IllegalArgumentException("Unknown parameter: "+arg);
			}
		}
		root=parser.in1;
		prefix=parser.out1;
		require(root!=null && prefix!=null && census!=null && population!=null && stage>0,
			"Required: in=completed_pilot census=full_summary population=sample_table out=fresh_prefix stage=64");
		require(!parser.overwrite && !parser.append, "Comparison outputs must be fresh; overwrite/append are unsupported");
		require(new File(root).isDirectory(), "Pilot directory not found: "+root);
		require(Tools.testInputFiles(false, true, census, population, path("bindings.tsv")), "Unreadable comparison metadata");
		for(final String suffix:SUFFIXES){
			require(Tools.testOutputFiles(false, false, false, prefix+suffix), "Output already exists or cannot be written: "+prefix+suffix);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Report Input         ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads completed points before writing a final summary; no searches are executed. */
	private void process(){
		require(line(path("exit.txt")).equals("0") && !new File(root, "STOP.txt").exists(), "Pilot failed or stopped at a resource limit");
		require(line(path("result.txt")).equals("SUBGA_NATIVE_TIMING_COMPLETE_NOT_REFERENCE_OR_ADMISSION"), "Pilot is not complete");
		loadPopulation();
		final LineParser1 p=new LineParser1('\t');
		final ObjectIntMap<String> seen=new ObjectIntMap<String>(String.class);
		final ByteBuilder b=new ByteBuilder();
		try(Input in=new Input(path("bindings.tsv"));
			Output modelOut=new Output(prefix+SUFFIXES[0]);
			Output differences=new Output(prefix+SUFFIXES[1]);
			Output rounding=new Output(prefix+SUFFIXES[2])){
			RrnaResourceIO.header(in.reader, BINDING_HEADER, "bindings.tsv");
			modelOut.writer.println(MODEL_HEADER);
			differences.writer.println("accession\tchange\tcount\tnative_table\twitness_line\ttarget\tmode\tmodel_from\tmodel_to\tseq_from\tseq_to\tstrand\ttruncation");
			rounding.writer.println("accession\tnative_table\twitness_line\tprinted_bits\tGA_bits\tincluded\tinterpretation");
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==8, "Expected eight binding fields");
				final String model=p.parseString(0);
				require(model.matches("RF[0-9]{5}") && !seen.contains(model), "Duplicate or malformed model: "+model);
				seen.put(model, models++);
				final double ga=p.parseDouble(1), report=p.parseDouble(2);
				require(Double.isFinite(ga) && Double.isFinite(report) && ga>0 && report>0 && report<ga, "Expected positive sub-GA report threshold");
				final String modelPath=p.parseString(4), oldCommand=p.parseString(5);
				final String oldHits=p.parseString(6), oldReport=p.parseString(7);
				final File control=new File(oldCommand).getParentFile();
				final String point=path("search_"+stage+'/'+model)+File.separator;
				require(line(new File(control, "exit.txt").getPath()).equals("0") && line(point+"exit.txt").equals("0"), "Incomplete native point for "+model);
				final String oldProfile=checkCommand(line(oldCommand), modelPath, ga, report, false);
				final String newProfile=checkCommand(line(point+"command.txt"), modelPath, ga, report, true);
				require(oldProfile.equals(newProfile), "Control and candidate search profiles or inputs differ for "+model);
				final double[] before=timing(new File(control, "time.tsv").getPath());
				final double[] after=timing(point+"time.tsv");
				final long oldBytes=Math.addExact(bytes(oldHits), bytes(oldReport));
				final long newBytes=Math.addExact(bytes(point+"hits.tbl"), bytes(point+"search.out"));
				final Counts counts=compareHits(oldHits, point+"hits.tbl", model, ga, differences, rounding);
				final boolean eligible=before[0]>=10 && after[0]>=10;
				if(!eligible){shortTimings++;}
				oldCpu+=before[0];
				newCpu+=after[0];
				nativeBytes=Math.addExact(nativeBytes, newBytes);
				b.clear().append(model).tab().append(ga, 6).tab().append(report, 6).tab();
				b.append(before[0], 2).tab().append(after[0], 2).tab().append(before[1], 2).tab().append(after[1], 2).tab();
				b.append(oldBytes).tab().append(newBytes).tab().append(counts.oldHits).tab().append(counts.newHits).tab();
				b.append(counts.oldIncluded).tab().append(counts.newIncluded).tab().append(counts.retained).tab();
				b.append(counts.lost).tab().append(counts.gained).tab().append(counts.roundedDisagreements).tab();
				if(eligible){
					b.append(after[0]*1e9/sampleBases, 6).tab().append(after[0]*fullBases/sampleBases/3600, 6).tab();
				}else{
					b.append("NA\tNA\t");
				}
				modelOut.writer.print(b.append(eligible ? "LINEAR_PROJECTION" : "SHORT_TIMING_NO_PROJECTION").nl());
			}
		}
		require(models>0, "No model bindings found");
		writeSummary();
	}

	/** Reads an explicit full census and the selected whole-genome sample dimensions. */
	private void loadPopulation(){
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(census)){
			RrnaResourceIO.header(in.reader, "genomes\tinput_file_bytes\tnucleotide_bases\tfasta_records\tboth_strands_search_space_Mb", census);
			p.set(required(in, census));
			require(p.terms()==5 && in.reader.nextLine()==null, "Expected one full population summary row");
			fullGenomes=p.parseLong(0);
			fullBases=p.parseLong(2);
		}
		try(Input in=new Input(population)){
			RrnaResourceIO.header(in.reader, "stage\tgenomes\tnucleotide_bases\tfasta_records", population);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==4, "Expected four sample population fields");
				if(p.parseInt(0)!=stage){continue;}
				require(sampleBases==0, "Repeated sample stage: "+stage);
				sampleGenomes=p.parseLong(1);
				sampleBases=p.parseLong(2);
			}
		}
		require(sampleBases>0 && sampleGenomes>0 && fullBases>=sampleBases && fullGenomes>=sampleGenomes,
			"Projection needs positive sample/full populations, with sample contained in full population");
	}

	/** Binds inclusion semantics to a deliberately narrow native command profile.
	 * Only the recorded unescaped, single-space argv form is accepted; this is not
	 * a shell evaluator. The returned profile omits output paths and the declared
	 * threshold changes; all other argv tokens must agree between paired points.
	 */
	public static String checkCommand(final String command, final String model, final double ga,
		final double report, final boolean candidate){
		require(command.indexOf('\\')<0 && command.indexOf('\'')<0 && command.indexOf('"')<0,
			"Escaped shell commands require an argv manifest; cannot infer inclusion binding");
		final String[] argv=command.trim().split(" +");
		require(argv.length>4 && argv[argv.length-2].equals(model), "Command model path differs from binding");
		int cut=0, inc=0, threshold=0;
		final ByteBuilder profile=new ByteBuilder();
		for(int i=0; i<argv.length-2; i++){
			final String arg=argv[i];
			if(arg.equals("--cut_ga")){
				cut++;
			}else if(arg.equals("--incT")){
				require(i+1<argv.length-2 && Double.parseDouble(argv[++i])==ga, "Native inclusion threshold must equal original GA");
				inc++;
			}else if(arg.equals("-T")){
				require(i+1<argv.length-2 && Double.parseDouble(argv[++i])==report, "Reporting threshold differs from binding");
				threshold++;
			}else if(arg.equals("--default")){
				// Explicit default matches the control's implicit default filters.
			}else if(arg.equals("--tblout") || arg.equals("-o")){
				require(++i<argv.length-2, "Missing native output path");
			}else{
				if(arg.startsWith("-")){
					require(arg.equals("--cpu") || arg.equals("-Z") || arg.equals("--tformat")
						|| arg.equals("--nohmmonly") || arg.equals("--noali") || arg.equals("--notextw"),
						"Unsupported search option; inclusion/filter binding is ambiguous: "+arg);
				}
				profile.append(arg).tab();
			}
		}
		require(candidate ? cut==0 && inc==1 && threshold==1 : cut==1 && inc==0 && threshold==0,
			"Candidate requires exactly -T/--incT; control requires exactly --cut_ga");
		return profile.append(model).tab().append(argv[argv.length-1]).toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------       Hit Comparison         ----------------*/
	/*--------------------------------------------------------------*/

	/** Counts native included rows, never calling that structural or biological confidence. */
	public static Counts compareHits(final String control, final String candidate, final String model,
		final double ga, final Output differences, final Output rounding){
		final Counts counts=new Counts();
		final ObjectIntMap<Key> remaining=new ObjectIntMap<Key>(Key.class);
		readHits(control, model, ga, true, remaining, counts, differences, rounding);
		readHits(candidate, model, ga, false, remaining, counts, differences, rounding);
		final Key[] keys=remaining.keys();
		final int[] values=remaining.values();
		final ByteBuilder b=new ByteBuilder();
		for(int i=0; i<keys.length; i++){
			if(keys[i]==null || values[i]==0){continue;}
			counts.lost+=values[i];
			witness(differences, b, model, "LOST", values[i], control, keys[i].line, keys[i]);
		}
		assert(counts.oldIncluded==counts.retained+counts.lost && counts.newIncluded==counts.retained+counts.gained) :
			"Each included report row must be consumed once or retained as a multiset difference";
		return counts;
	}

	/** Streams byte tokens; only previously unseen control keys allocate retained bytes. */
	private static void readHits(final String path, final String model, final double ga, final boolean control,
		final ObjectIntMap<Key> remaining, final Counts counts, final Output differences, final Output rounding){
		final RrnaCmTruthBuilder.Fields f=new RrnaCmTruthBuilder.Fields();
		final ByteBuilder keyBytes=new ByteBuilder(), b=new ByteBuilder();
		final Key probe=new Key();
		long line=0;
		boolean complete=false;
		try(Input in=new Input(path)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				line++;
				if(row.length==0){continue;}
				if(row[0]=='#'){
					if(Arrays.equals(row, TABLE_END)){complete=true;}
					continue;
				}
				if(complete){throw new IllegalArgumentException("Native data follows completion marker: "+path);}
				f.set(row);
				if(!f.eq(3, model) || !(f.eq(16, "!") || f.eq(16, "?"))){
					throw new IllegalArgumentException("Wrong model or inclusion flag: "+path+':'+line);
				}
				final double bits=Parse.parseDouble(row, f.start[14], f.end[14]);
				if(!Double.isFinite(bits)){throw new IllegalArgumentException("Nonfinite native score: "+path+':'+line);}
				final boolean included=f.eq(16, "!");
				if(control){counts.oldHits++;}else{counts.newHits++;}
				if(included!=(bits>=ga)){
					counts.roundedDisagreements++;
					if(rounding!=null){
						b.clear().append(model).tab().append(path).tab().append(line).tab().append(bits, 1).tab().append(ga, 6).tab();
						rounding.writer.print(b.append(included ? '!' : '?').append("\tPRINTED_SCORE_DIFFERS_FROM_BOUND_NATIVE_GA_DECISION\n"));
					}
				}
				if(!included){continue;}
				keyBytes.clear();
				for(final int field:KEY_FIELDS){
					if(keyBytes.length()>0){keyBytes.tab();}
					keyBytes.append(row, f.start[field], f.end[field]-f.start[field]);
				}
				probe.set(keyBytes.array, keyBytes.length(), line);
				final int old=remaining.get(probe);
				if(control){
					counts.oldIncluded++;
					require(old<Integer.MAX_VALUE, "Control duplicate count exceeds integer capacity");
					if(old<0){remaining.put(probe.copy(), 1);}else{remaining.put(probe, old+1);}
				}else{
					counts.newIncluded++;
					if(old>0){
						remaining.put(probe, old-1);
						counts.retained++;
					}else{
						counts.gained++;
						witness(differences, b, model, "GAINED", 1, path, line, probe);
					}
				}
			}
		}
		require(complete, "Missing native completion marker: "+path);
	}

	/** Writes one coordinate witness; gained duplicates may occupy separate rows. */
	private static void witness(final Output out, final ByteBuilder b, final String model, final String change,
		final int count, final String path, final long line, final Key key){
		if(out==null){return;}
		out.writer.print(b.clear().append(model).tab().append(change).tab().append(count).tab().append(path).tab()
			.append(line).tab().append(key.bytes, 0, key.length).nl());
	}

	/*--------------------------------------------------------------*/
	/*----------------        Cost Accounting       ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes total new-pass cost, separate from growth over already-spent control CPU. */
	private void writeSummary(){
		final double scale=(double)fullBases/sampleBases;
		final ByteBuilder b=new ByteBuilder();
		try(Output out=new Output(prefix+SUFFIXES[3])){
			out.writer.println("metric\tvalue\tunits\tinterpretation");
			metric(out, b, "models", models, "models", "COMPLETED_NATIVE_PAIRS");
			metric(out, b, "short_timings", shortTimings, "models", "BELOW_TEN_CPU_SECONDS_NO_TIMING_PROJECTION");
			metric(out, b, "sample_bases", sampleBases, "bases", "MEASURED_POPULATION");
			metric(out, b, "full_bases", fullBases, "bases", "MEASURED_POPULATION");
			metric(out, b, "control_CPU", oldCpu, "CPU_seconds", "HISTORICAL_NONCONCURRENT_CONTROLS");
			metric(out, b, "candidate_CPU", newCpu, "CPU_seconds", "MEASURED_SUM_OF_POINT_CPU");
			metric(out, b, "native_reports", nativeBytes, "bytes", "HITS_TABLE_PLUS_SEARCH_REPORT");
			if(shortTimings==0){
				metric(out, b, "projected_new_pass_CPU", newCpu*scale/3600, "CPU_hours", "LINEAR_BASE_SCALING_NOT_MEASURED_FULL_SCAN");
				metric(out, b, "projected_control_CPU", oldCpu*scale/3600, "CPU_hours", "HISTORICAL_CONTROL_SAME_SAMPLE_SCALING");
				metric(out, b, "candidate_to_control_CPU_ratio", newCpu/oldCpu, "ratio", "NONCONCURRENT_TIMINGS");
				metric(out, b, "ideal_8worker_wall", newCpu*scale/3600/8, "hours", "IDEAL_BALANCE_NO_QUEUE_STARTUP_IO_OR_CONTENTION");
			}
			metric(out, b, "projected_reports", nativeBytes*scale/GIB, "GiB", "LINEAR_BASE_SCALING_INCLUDES_SAMPLE_HEADERS");
			metric(out, b, "per_genome_sample_report_allowance", nativeBytes*(double)fullGenomes/GIB, "GiB",
				"PLANNING_ALLOWANCE_NOT_BOUND_EACH_GENOME_GETS_ENTIRE_SAMPLE_REPORT");
			metric(out, b, "two_report_files_per_pair", 2.0*models*fullGenomes, "files", "IF_ONE_NATIVE_SEARCH_PER_GENOME_MODEL_EXCLUDES_METADATA_FILES");
		}
	}

	private static void metric(final Output out, final ByteBuilder b, final String name,
		final double value, final String unit, final String meaning){
		require(Double.isFinite(value) && value>=0, "Nonfinite or negative cost metric: "+name);
		out.writer.print(b.clear().append(name).tab().append(value, 6).tab().append(unit).tab().append(meaning).nl());
	}
	private static double[] timing(final String path){
		final LineParser1 p=new LineParser1('\t');
		p.set(line(path).getBytes(StandardCharsets.US_ASCII));
		require(p.terms()==4, "Native timing needs user/system/wall/RSS: "+path);
		final double user=p.parseDouble(0), system=p.parseDouble(1), wall=p.parseDouble(2);
		require(Double.isFinite(user) && Double.isFinite(system) && Double.isFinite(wall)
			&& user>=0 && system>=0 && wall>0, "Invalid native timing: "+path);
		return new double[]{user+system, wall};
	}
	private static long bytes(final String path){
		final File file=new File(path);
		require(file.isFile() && file.canRead(), "Missing native report: "+path);
		return file.length();
	}
	private String path(final String relative){return new File(root, relative).getPath();}
	private static byte[] required(final Input in, final String path){
		final byte[] row=in.reader.nextLine();
		require(row!=null, "Empty required file: "+path);
		return row;
	}
	private static String line(final String path){
		try(Input in=new Input(path)){
			final String value=new String(required(in, path), StandardCharsets.UTF_8);
			require(in.reader.nextLine()==null, "Expected one line: "+path);
			return value;
		}
	}
	private static void require(final boolean condition, final String reason){
		if(!condition){throw new IllegalArgumentException(reason);}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Retained Objects       ----------------*/
	/*--------------------------------------------------------------*/

	/** Event counts; retained/lost/gained refer only to the bound native inclusion cutoff. */
	public static final class Counts{
		public long oldHits, newHits, oldIncluded, newIncluded, retained, lost, gained, roundedDisagreements;
	}
	/** Borrowed byte probe until copy(); stored keys own immutable byte arrays. */
	private static final class Key{
		void set(final byte[] bytes_, final int length_, final long line_){
			assert(bytes_!=null && length_>=0 && length_<=bytes_.length) : "Borrowed key range must fit its normalization buffer";
			bytes=bytes_;
			length=length_;
			line=line_;
			hash=1;
			for(int i=0; i<length; i++){hash=31*hash+bytes[i];}
		}
		Key copy(){
			final Key key=new Key();
			key.set(Arrays.copyOf(bytes, length), length, line);
			return key;
		}
		@Override
		public int hashCode(){return hash;}
		@Override
		public boolean equals(final Object object){
			if(!(object instanceof Key)){return false;}
			final Key other=(Key)object;
			if(length!=other.length || hash!=other.hash){return false;}
			for(int i=0; i<length; i++){if(bytes[i]!=other.bytes[i]){return false;}}
			return true;
		}
		byte[] bytes;
		int length, hash;
		long line;
	}

	private String root, prefix, census, population;
	private int stage=64, models, shortTimings;
	private long sampleBases, sampleGenomes, fullBases, fullGenomes, nativeBytes;
	private double oldCpu, newCpu;
	private PrintStream outstream=System.err;
	private static final int[] KEY_FIELDS={0, 4, 5, 6, 7, 8, 9, 10};
	private static final byte[] TABLE_END="# [ok]".getBytes(StandardCharsets.US_ASCII);
	private static final double GIB=1073741824.0;
	private static final String BINDING_HEADER="accession\tGA_bits\treport_bits\tmodel_sha80\tmodel_path\tcontrol_command\tcontrol_hits\tcontrol_report";
	private static final String MODEL_HEADER="accession\tGA_bits\treport_bits\tcontrol_CPU_seconds\tcandidate_CPU_seconds\tcontrol_wall_seconds\tcandidate_wall_seconds\tcontrol_bytes\tcandidate_bytes\tcontrol_hits\tcandidate_hits\tcontrol_included\tcandidate_included\tretained\tlost\tgained\tprinted_score_flag_disagreements\tCPU_seconds_per_Gbp\tprojected_full_CPU_hours\tstatus";
	private static final String[] SUFFIXES={".model_costs.tsv", ".hit_differences.tsv", ".rounding.tsv", ".cost_summary.tsv"};
}
