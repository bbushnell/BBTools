package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import map.LongHashSet;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Compares completed GradeNcrnaReference reports without rematching calls.
 * Adapted from template.A_SampleByteFile: native byte input/output, standard
 * argument parsing and checked stream closure. Input populations, reference
 * identities/dispositions and grading settings must agree; caller settings may
 * differ. This is a fixed-reference caller comparison, not a reference-swap tool.
 *
 * Reference keys are scoped by input and family. A reusable borrowed byte probe
 * performs lookups; only retained keys and records allocate. Original summaries
 * are checked against detailed rows before any output opens. Declared family-map
 * equality is not independent proof of how historical reports were produced.
 * Overlap equality uses the grader's six-decimal report, not an unavailable
 * historical float. The owning recipe retains the original exact configurations.
 * Memory is proportional to retained reference/TP rows and the total query-line
 * set capacity across inputs, not genome sequence size. Clearing a set reuses
 * its capacity. All computation is confined to one thread.
 * @author Brian Bushnell, Raiden
 */
public final class CompareNcrnaGrades{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final CompareNcrnaGrades tool=new CompareNcrnaGrades(args);
		try{tool.process();}
		finally{Shared.closeStream(tool.outstream);}
		timer.stop();
		System.err.println("Compared completed grades without rematching. Time: "+timer);
	}

	private CompareNcrnaGrades(final String[] args){
		final PreParser pp=new PreParser(args, CompareNcrnaGrades.class, false);
		outstream=pp.outstream;
		Shared.setThreads(1);
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("familymap1")){maps[0]=b;}
			else if(a.equals("familymap2")){maps[1]=b;}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		prefixes[0]=parser.in1;
		prefixes[1]=parser.in2;
		out=parser.out1;
		require(out!=null && prefixes[0]!=null && prefixes[1]!=null && maps[0]!=null && maps[1]!=null,
			"Supply in1=baseline_prefix in2=candidate_prefix familymap1=file familymap2=file out=fresh_prefix");
		require(!parser.overwrite && !parser.append, "Comparison reports must be fresh; overwrite/append are unsupported");
		for(final String prefix:prefixes){
			require(Tools.testInputFiles(false, true, prefix+".raw.tsv", prefix+".inputs.tsv", prefix+".summary.tsv"),
				"Both grades require completed raw, input and summary reports");
		}
		require(Tools.testInputFiles(true, true, maps), "Both original grading family maps must be readable; a shared map is expected");
		require(Tools.testOutputFiles(false, false, false, out+".loci.tsv", out+".counts.tsv", out+".summary.tsv"),
			"Comparison outputs must be fresh and writable");
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads each saved report once; no genome, caller or matching operation runs. */
	private void process(){
		compareMaps();
		for(int side=0; side<2; side++){
			readInputs(side);
			readSummary(side);
			readRows(side);
			validate(side);
		}
		writeReports();
		outstream.println("Inputs: "+inputs.size()+"; families: "+families.size());
	}

	/** The recipe supplies the original maps; equal declarations do not certify their history. */
	private void compareMaps(){
		try(Input left=new Input(maps[0]); Input right=new Input(maps[1])){
			RrnaResourceIO.header(left.reader, GradeNcrnaReference.RULE_HEADER, maps[0]);
			RrnaResourceIO.header(right.reader, GradeNcrnaReference.RULE_HEADER, maps[1]);
			while(true){
				final byte[] a=left.reader.nextLine(), b=right.reader.nextLine();
				require(Arrays.equals(a, b), "Grading family maps differ; fixed-reference comparisons require identical rules");
				if(a==null){break;}
			}
		}
	}

	/** Retains input IDs and paths, including inputs with no query features. */
	private void readInputs(final int side){
		final String path=prefixes[side]+".inputs.tsv";
		int count=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, GradeNcrnaReference.INPUT_HEADER, path);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==10, "Input report requires ten fields: "+path);
				probe.set(p, 0);
				InputState state=byInput.get(probe);
				if(side==0){
					require(state==null, "Duplicate baseline input ID");
					state=new InputState(line);
					byInput.put(probe.copy(), state);
					inputs.add(state);
				}else{
					require(state!=null && !state.candidateSeen, "Candidate input ID is absent from baseline or duplicated");
					q.set(state.metadata);
					require(equal(1, 1) && equal(3, 3), "FASTA or reference GFF paths differ for an input");
					state.candidateSeen=true;
					state.queryLines.clear();
				}
				state.expectedIncluded=p.parseLong(8);
				state.expectedExcluded=p.parseLong(9);
				require(state.expectedIncluded>=0 && state.expectedExcluded>=0, "Query counts cannot be negative");
				state.included=state.excluded=0;
				count++;
			}
		}
		require(count==inputs.size(), "Input populations differ, including zero-call inputs");
	}

	/** Summary presence is the producer's completion marker; raw counts corroborate it. */
	private void readSummary(final int side){
		final String path=prefixes[side]+".summary.tsv";
		int count=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, GradeNcrnaReference.SUMMARY_HEADER, path);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==17 && p.parseLong(1)==inputs.size(), "Summary schema or input population mismatch");
				probe.set(p, 0);
				Stats stats=byFamily.get(probe);
				if(side==0){
					require(stats==null, "Duplicate baseline summary family");
					stats=new Stats();
					byFamily.put(probe.copy(), stats);
					families.add(stats);
				}else{
					require(stats!=null && stats.summary[1]==null, "Candidate family is missing or duplicated in the baseline");
					q.set(stats.summary[0]);
					require(equal(15, 15) && equal(16, 16), "Overlap settings or grading qualification differ");
				}
				stats.summary[side]=line;
				for(int i=0; i<11; i++){
					stats.expected[side][i]=p.parseLong(i+2);
					require(stats.expected[side][i]>=0, "Summary counts cannot be negative");
				}
				count++;
			}
		}
		require(count==families.size() && count>0, "Summary family populations differ or are empty");
	}

	/** Processes byte rows with retained reference state and primitive query counters. */
	private void readRows(final int side){
		final String path=prefixes[side]+".raw.tsv";
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, GradeNcrnaReference.RAW_HEADER, path);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==18, "Detailed grade report requires eighteen fields");
				probe.set(p, 0);
				final InputState input=byInput.get(probe);
				require(input!=null, "Detailed grade row names an undeclared input");
				probe.set(p, 2);
				final Stats stats=byFamily.get(probe);
				require(stats!=null, "Detailed grade row names an undeclared family");
				if(p.termEquals("REFERENCE", 1)){reference(input, stats, side, line);}
				else if(p.termEquals("QUERY", 1)){query(input, stats, side);}
				else{throw new IllegalArgumentException("Unknown detailed grade role");}
			}
		}
	}

	/** Reference identity/geometry/model/disposition must remain fixed across callers. */
	private void reference(final InputState input, final Stats stats, final int side, final byte[] line){
		final Entry entry=entry(input, stats, 8, side);
		final Observation observation=entry.observation(side);
		require(!observation.seen, "Duplicate reference ID within an input and family");
		observation.seen=true;
		final int outcome=p.termEquals("TP", 7) ? TP : p.termEquals("FN", 7) ? FN : p.termEquals("DONT_CARE", 7) ? MASK : -1;
		require(outcome>=0, "Unknown reference outcome");
		observation.outcome=outcome;
		observation.matchedLine=p.termEquals(".", 11) ? -1 : p.parseLong(11);
		if(side==0){entry.reference=line;}
		else{
			q.set(entry.reference);
			for(final int column:REFERENCE_COLUMNS){
				require(equal(column, column), "Reference identity, geometry or model changed between caller grades");
			}
			require((outcome==MASK)==(entry.baseline.outcome==MASK), "Reference disposition changed; reference-swap comparisons are unsupported");
		}
		final long[] counts=stats.actual[side];
		if(outcome==MASK){counts[4]++;}
		else{
			counts[3]++;
			if(outcome==FN){counts[2]++;}
		}
	}

	/** Attaches TP offsets through matched_reference_id, never the query's empty reference_id. */
	private void query(final InputState input, final Stats stats, final int side){
		final long line=p.parseLong(9);
		require(line>0 && input.queryLines.add(line), "Query line must be positive and unique within its input");
		final long[] counts=stats.actual[side];
		if(p.termEquals("EXCLUDED", 7)){
			input.excluded++;
			counts[7]++;
			return;
		}
		input.included++;
		require(p.termEquals("true", 13) || p.termEquals("false", 13), "Evaluated duplicate status must be explicit");
		if(p.termEquals("true", 13)){counts[6]++;}
		if(p.termEquals("TP", 7)){
			final Entry entry=entry(input, stats, 10, side);
			final Observation observation=entry.observation(side);
			require(observation.queryLine<0, "Two query TPs consume the same reference");
			observation.queryLine=line;
			observation.query=p.line();
			observation.d5=p.parseInt(16);
			observation.d3=p.parseInt(17);
			counts[0]++;
			if(observation.d5==0){counts[8]++;}
			if(observation.d3==0){counts[9]++;}
			if(observation.d5==0 && observation.d3==0){counts[10]++;}
		}else if(p.termEquals("MASKED", 7)){
			entry(input, stats, 12, side).observation(side).maskCalls++;
			counts[5]++;
		}else if(p.termEquals("FP", 7)){counts[1]++;}
		else{throw new IllegalArgumentException("Unknown query outcome");}
	}

	/** Borrowed family/tab/reference key; one owned copy is stored per retained locus. */
	private Entry entry(final InputState input, final Stats stats, final int idColumn, final int side){
		require(!p.termEquals(".", idColumn), "A reference or witness must name an explicit reference ID");
		key.clear();
		appendField(key, p, 2).tab();
		appendField(key, p, idColumn);
		probe.set(key.array, 0, key.length());
		Entry entry=input.byReference.get(probe);
		if(entry==null){
			require(side==0, "Candidate grade contains a reference identity absent from the baseline");
			entry=new Entry(stats);
			input.byReference.put(probe.copy(), entry);
			input.references.add(entry);
		}
		return entry;
	}

	/** Detects missing/duplicated detailed records before publishing comparison results. */
	private void validate(final int side){
		for(final InputState input:inputs){
			require(input.included==input.expectedIncluded && input.excluded==input.expectedExcluded,
				"Detailed query counts disagree with the input report; incomplete or mixed grade files");
			for(final Entry entry:input.references){
				final Observation observation=entry.observation(side);
				require(observation.seen, "Missing reference row or unresolved query witness");
				if(observation.outcome==TP){
					require(observation.queryLine>0 && observation.queryLine==observation.matchedLine,
						"TP reference/query witnesses do not agree on the original query line");
					validateWitness(entry, observation);
				}else{
					require(observation.queryLine<0 && observation.matchedLine<0, "Unmatched/masked reference has a TP witness");
				}
				require(observation.maskCalls==0 || observation.outcome==MASK, "Masked query points to a scored reference");
			}
		}
		for(final Stats stats:families){
			require(Arrays.equals(stats.actual[side], stats.expected[side]), "Detailed family counts disagree with the completed summary");
		}
	}

	/** Checks the saved pairing's geometry and transcript-direction offsets.
	 * Does not search for a different partner or repeat the overlap matcher. */
	private void validateWitness(final Entry entry, final Observation observation){
		assert(entry.reference!=null && observation.query!=null) : "A validated TP line link must retain both rows for geometry checking";
		p.set(observation.query);
		q.set(entry.reference);
		require(equal(3, 3) && equal(6, 6), "TP witness contig or strand differs from its paired reference");
		require(p.termEquals("+", 6) || p.termEquals("-", 6), "TP witness requires an explicit strand");
		final int queryStart=p.parseInt(4), queryStop=p.parseInt(5);
		final int referenceStart=q.parseInt(4), referenceStop=q.parseInt(5);
		require(queryStart>0 && queryStop>=queryStart && referenceStart>0 && referenceStop>=referenceStart,
			"TP witness coordinates must be positive ordered intervals");
		final boolean minus=p.termEquals("-", 6);
		final int d5=minus ? referenceStop-queryStop : queryStart-referenceStart;
		final int d3=minus ? referenceStart-queryStart : queryStop-referenceStop;
		require(observation.d5==d5 && observation.d3==d3, "TP witness signed offsets disagree with its paired geometry");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Reports           ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes all loci and the unchanged native count rows; summary closes last. */
	private void writeReports(){
		try(Output loci=new Output(out+".loci.tsv"); Output counts=new Output(out+".counts.tsv")){
			loci.writer.println("input_id\tfamily\treference_id\tseqid\tstart\tstop\tstrand\tbaseline\tcandidate\tbaseline_d5\tbaseline_d3\tcandidate_d5\tcandidate_d3\tendpoint_changed\tlost_exact5\tlost_exact3");
			counts.writer.println("run\t"+GradeNcrnaReference.SUMMARY_HEADER);
			for(final Stats stats:families){
				for(int side=0; side<2; side++){
					buffer.clear().append(side==0 ? "baseline\t" : "candidate\t").append(stats.summary[side]).nl();
					counts.writer.print(buffer);
				}
			}
			for(final InputState input:inputs){
				for(final Entry entry:input.references){writeLocus(entry, loci);}
			}
		}
		try(Output summary=new Output(out+".summary.tsv")){
			summary.writer.println("family\tretained_TP\tlost_TP\tgained_TP\tshared_FN\tmasked_references\tchanged_endpoint_loci\tlost_exact5\tlost_exact3\tlost_either_exact_end");
			for(final Stats stats:families){
				p.set(stats.summary[0]);
				buffer.clear();
				appendField(buffer, p, 0);
				for(final long value:stats.changes){buffer.tab().append(value);}
				summary.writer.print(buffer.nl());
			}
		}
	}

	/** Endpoint losses count only retained TPs; detection losses have their own column. */
	private void writeLocus(final Entry entry, final Output out){
		final Observation baseline=entry.baseline, candidate=entry.candidate;
		assert(baseline.seen && candidate.seen) : "Validation requires each reference in both completed grades before reporting";
		final int before=baseline.outcome, after=candidate.outcome;
		final boolean retained=before==TP && after==TP;
		final boolean changed=retained && (baseline.d5!=candidate.d5 || baseline.d3!=candidate.d3);
		final boolean lost5=retained && baseline.d5==0 && candidate.d5!=0;
		final boolean lost3=retained && baseline.d3==0 && candidate.d3!=0;
		final long[] changes=entry.stats.changes;
		if(before==MASK){changes[4]++;}
		else if(retained){changes[0]++;}
		else if(before==TP){changes[1]++;}
		else if(after==TP){changes[2]++;}
		else{changes[3]++;}
		if(changed){changes[5]++;}
		if(lost5){changes[6]++;}
		if(lost3){changes[7]++;}
		if(lost5 || lost3){changes[8]++;}
		p.set(entry.reference);
		buffer.clear();
		for(final int column:LOCUS_COLUMNS){appendField(buffer, p, column).tab();}
		buffer.append(OUTCOMES[before]).tab().append(OUTCOMES[after]);
		for(int side=0; side<2; side++){
			final Observation observation=entry.observation(side);
			if(observation.outcome==TP){buffer.tab().append(observation.d5).tab().append(observation.d3);}
			else{buffer.append("\t.\t.");}
		}
		out.writer.print(buffer.tab().append(changed).tab().append(lost5).tab().append(lost3).nl());
	}

	/** Compares borrowed fields without decoding or allocating transient Strings. */
	private boolean equal(final int left, final int right){
		p.setBounds(left);
		q.setBounds(right);
		final int length=p.b()-p.a();
		if(length!=q.b()-q.a()){return false;}
		for(int i=0; i<length; i++){if(p.line()[p.a()+i]!=q.line()[q.a()+i]){return false;}}
		return true;
	}
	private static ByteBuilder appendField(final ByteBuilder buffer, final LineParser1 parser, final int column){
		parser.setBounds(column);
		return buffer.append(parser.line(), parser.a(), parser.b()-parser.a());
	}
	private static void require(final boolean valid, final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Types            ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained metadata and keys for one genome; query-line sets reset between runs. */
	private static final class InputState{
		InputState(final byte[] metadata_){metadata=metadata_;}
		final byte[] metadata;
		final ObjectMap<Target, Entry> byReference=new ObjectMap<Target, Entry>(16, Target.class, Entry.class);
		final ArrayList<Entry> references=new ArrayList<Entry>();
		final LongHashSet queryLines=new LongHashSet(16);
		long expectedIncluded, expectedExcluded, included, excluded;
		boolean candidateSeen;
	}
	/** One reference and its already assigned query witnesses in the two runs. */
	private static final class Entry{
		Entry(final Stats stats_){stats=stats_;}
		Observation observation(final int side){return side==0 ? baseline : candidate;}
		final Stats stats;
		byte[] reference;
		final Observation baseline=new Observation(), candidate=new Observation();
	}
	/** Primitive state retained for one locus/run; no per-field wrappers or arrays. */
	private static final class Observation{
		byte[] query;
		boolean seen;
		int outcome, d5, d3;
		long queryLine=-1, matchedLine=-1, maskCalls;
	}
	/** Counts follow native summary columns2..12; original rows preserve precision/recall. */
	private static final class Stats{
		final byte[][] summary=new byte[2][];
		final long[][] expected=new long[2][11], actual=new long[2][11];
		final long[] changes=new long[9];
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final PrintStream outstream;
	private final String[] prefixes=new String[2], maps=new String[2];
	private String out;
	private final LineParser1 p=new LineParser1('\t'), q=new LineParser1('\t');
	private final ByteBuilder buffer=new ByteBuilder(), key=new ByteBuilder();
	private final Target probe=new Target();
	private final ObjectMap<Target, InputState> byInput=new ObjectMap<Target, InputState>(16, Target.class, InputState.class);
	private final ObjectMap<Target, Stats> byFamily=new ObjectMap<Target, Stats>(16, Target.class, Stats.class);
	private final ArrayList<InputState> inputs=new ArrayList<InputState>();
	private final ArrayList<Stats> families=new ArrayList<Stats>();
	private static final int TP=0, FN=1, MASK=2;
	private static final String[] OUTCOMES={"TP", "FN", "DONT_CARE"};
	private static final int[] REFERENCE_COLUMNS={2, 3, 4, 5, 6, 8, 15};
	private static final int[] LOCUS_COLUMNS={0, 2, 8, 3, 4, 5, 6};
}
