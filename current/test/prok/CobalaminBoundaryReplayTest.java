package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.CompareCobalaminBenchmarkGrades;
import prok.NcrnaBoundaryTrace;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Literal final calls and score windows, independent of a trained network.
 * Both strands share contig names; one end is physically fixed. A scalar first
 * move is harmful, while gating it exposes a helpful held-second context.
 * Positional windows test the independent-end route. Unknowns cannot masquerade
 * as no movement, and wrong bindings/observed endpoints must fail loudly.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminBoundaryReplayTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, CobalaminBoundaryReplayTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown replay fixture argument: "+arg);
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append && parser.in1==null && parser.in2==null
				&& parser.maxReads<0 && parser.samplerate==1 && new File(parser.out1).mkdir(), "Supply out=fresh_fixture_directory");
			new CobalaminBoundaryReplayTest(parser.out1).run();
			pp.outstream.println("COBALAMIN_BOUNDARY_REPLAY_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private CobalaminBoundaryReplayTest(final String directory){root=new File(directory);}

	private void run() throws Exception{
		write("sources.tsv", "input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff\n"
			+"0\ta\ta\t/a.fa\t/a.gff\n1\tb\tb\t/b.fa\t/b.gff\n2\tc\tc\t/c.fa\t/c.gff\n3\tempty\tempty\t/empty.fa\t/empty.gff\n");
		views(false);
		write("scalar.trace.tsv", trace(false, false));
		invoke("scalar.trace.tsv", "sources.tsv", "scalar");
		final String unbounded=trace(false, false).replace("\t16\t1000\t", "\t16\t\t");
		write("unbounded.trace.tsv", unbounded);
		invoke("unbounded.trace.tsv", "sources.tsv", "unbounded");
		write("numeric_max_v2.trace.tsv", unbounded.replace("#ncrna_boundary_trace\t3", "#ncrna_boundary_trace\t2")
			.replace("\t16\t\t", "\t16\t2147483647\t"));
		invoke("numeric_max_v2.trace.tsv", "sources.tsv", "numeric_max_v2");
		check(read("unbounded.replay_endpoints.tsv").equals(read("numeric_max_v2.replay_endpoints.tsv")),
			"Version-three absence and the exact version-two unbounded control preserve complete replay results");
		for(final String bad:new String[]{"-1", "2147483647", "2147483648", "4294967496", "NA"}){
			final String name="bad_max_"+bad;
			write(name+".tsv", unbounded.replace("\t16\t\t", "\t16\t"+bad+"\t"));
			fails(name+".tsv", "sources.tsv", name);
		}
		LineParser1 row=row("scalar.replay_endpoints.tsv", "cobalamin", "5prime", 0);
		check(row.parseLong(5)==3 && row.parseLong(6)==3 && row.parseLong(7)==0, "All three paired references must stay in every replay denominator");
		close(row.parseDouble(8), 0);
		close(row.parseDouble(9), Math.sqrt(2.0/3));
		check(row.parseLong(10)==0 && row.parseLong(11)==2 && row.parseLong(12)==1, "Ungated first-end harm must remain visible on both strands");
		row=row("scalar.replay_endpoints.tsv", "cobalamin", "3prime", .1f);
		close(row.parseDouble(8), 1);
		close(row.parseDouble(9), 0);
		check(row.parseLong(10)==3 && row.parseLong(11)==0 && row.parseLong(12)==0 && row.parseLong(17)==3,
			"Gating the first move must select the HELD context and improve three stop ends");
		check(read("scalar.replay_selection.tsv").contains("\t0.1\tSELECTED_DIAGNOSTIC_MARGIN\t"), "Equal full squared-error optima choose the largest declared positive margin");
		invoke("scalar.trace.tsv", "sources.tsv", "extended", "0,0.005,0.02,0.1,0.2,0.4");
		check(read("extended.replay_selection.tsv").contains("\t0.4\tSELECTED_DIAGNOSTIC_MARGIN\t")
			&& read("scalar.replay_selection.tsv").contains("\t0.1\tSELECTED_DIAGNOSTIC_MARGIN\t"),
			"An explicitly extended diagnostic list must be evaluated without changing the original selection");
		row=row("extended.replay_endpoints.tsv", "cobalamin", "3prime", .4f);
		close(row.parseDouble(9), 0);
		check(row.parseLong(5)==3 && row.parseLong(6)==3, "Added margins keep the same fixed reference cohort");
		invoke("scalar.trace.tsv", "sources.tsv", "precise", "0,0.0001,0.1");
		row=row("precise.replay_endpoints.tsv", "cobalamin", "5prime", .0001f);
		close(row.parseDouble(9), Math.sqrt(2.0/3));
		int badIndex=0;
		for(final String margins:new String[]{"0,0.1,0.005", "0,0.1,0.1", "0,NaN", "0.1,0.2"}){
			fails("scalar.trace.tsv", "sources.tsv", "bad_margins"+(badIndex++), margins);
		}
		row=row("scalar.replay_endpoints.tsv", "RF01482", "5prime", .1f);
		check(row.parseLong(5)==0 && row.length(8)==0 && row.length(13)==0, "ON label switches cannot populate an OFF-member cohort");
		check(read("scalar.replay_loci.tsv").contains("b\tr\tRF00174\tRF01482\t0.1\tREPLAYED\t0\t0\t-1\t0\t"),
			"Reverse coordinates must produce literal zero errors after the positive-margin correction");
		check(read("scalar.replay_accounting.tsv").contains("a\t3\t3\t2\t0\t0\n"), "Rejected and unmatched calls must survive conservation outside paired error metrics");
		CompareCobalaminBenchmarkGrades.main(new String[]{"in1="+path("on"), "in2="+path("on"), "out="+path("unchanged"),
			"radius=1", "trace="+path("scalar.trace.tsv"), "bindings="+path("sources.tsv")});
		check(read("unchanged.replay_selection.tsv").contains("NONE_QUALIFIES"), "An ungated no-change comparison cannot select a demonstrated correction");
		final String strong=Float.floatToIntBits(.1f)+","+Float.floatToIntBits(.2f)+","+Float.floatToIntBits(.9f);
		final String weak=Float.floatToIntBits(.1f)+","+Float.floatToIntBits(.2f)+","+Float.floatToIntBits(.203f);
		write("no_move.trace.tsv", trace(false, false).replace(strong, weak));
		invoke("no_move.trace.tsv", "sources.tsv", "no_move");
		check(read("no_move.replay_selection.tsv").contains("NONE_QUALIFIES"), "A no-move positive margin with undefined better fraction cannot qualify");
		row=row("no_move.replay_endpoints.tsv", "cobalamin", "3prime", .1f);
		check(row.length(13)==0 && row.parseLong(17)==0 && row.parseLong(18)==3, "No-move ends remain explicit, with no better fraction");
		CompareCobalaminBenchmarkGrades.main(new String[]{"in1="+path("off"), "in2="+path("on"), "out="+path("plain"), "radius=1"});
		for(final String suffix:new String[]{"summary", "genomes", "loci", "overlaps", "lengths", "paired_endpoints", "paired_loci"}){
			check(read("scalar."+suffix+".tsv").equals(read("plain."+suffix+".tsv")), "Replay must not alter an existing comparison output: "+suffix);
		}
		write("unknown.trace.tsv", trace(false, true));
		invoke("unknown.trace.tsv", "sources.tsv", "unknown");
		check(read("unknown.replay_selection.tsv").contains("UNAVAILABLE_INCOMPLETE_CAPTURE")
			&& read("unknown.replay_loci.tsv").contains("UNKNOWN_ELIGIBILITY"), "Unknown eligibility cannot be assigned a no-move prediction or a selected margin");
		final String clipped="S\t2\t2\t0\t1\t0\t0\t59\t";
		final String present=clipped+Float.floatToIntBits(.1f)+",";
		final String absent=clipped+",";
		final String optional=trace(false, true).replace(present, absent);
		check(!optional.equals(trace(false, true)), "The compatibility case must actually omit the clipped score position");
		write("optional.trace.tsv", optional);
		invoke("optional.trace.tsv", "sources.tsv", "optional");
		write("legacy.trace.tsv", legacyTrace(optional));
		check(read("optional.trace.tsv").contains("SELECTED\t\t\t\n")
			&& read("legacy.trace.tsv").contains("SELECTED\t-1\t-1\t-1\n"),
			"Persisted unknown-window metadata must use the literal current and legacy encodings");
		invoke("legacy.trace.tsv", "sources.tsv", "legacy");
		for(final String suffix:new String[]{"endpoints", "loci", "selection", "accounting"}){
			check(read("optional.replay_"+suffix+".tsv").equals(read("legacy.replay_"+suffix+".tsv")),
				"Legacy absence must normalize to the complete current replay output: "+suffix);
		}
		for(final float invalid:new float[]{Float.NaN, Float.POSITIVE_INFINITY, Float.NEGATIVE_INFINITY}){
			final String name="invalid_score_"+Float.floatToIntBits(invalid);
			write(name+".tsv", optional.replace(absent, clipped+Float.floatToIntBits(invalid)+","));
			fails(name+".tsv", "sources.tsv", name);
		}
		write("legacy_infinite.tsv", legacyTrace(optional).replace(clipped+Float.floatToIntBits(Float.NaN)+",",
			clipped+Float.floatToIntBits(Float.POSITIVE_INFINITY)+","));
		fails("legacy_infinite.tsv", "sources.tsv", "legacy_infinite");
		write("wrong.sources.tsv", read("sources.tsv").replace("/a.fa", "/wrong.fa"));
		fails("scalar.trace.tsv", "wrong.sources.tsv", "wrong_source");
		write("wrong.trace.tsv", trace(false, false).replace("\t79\t139\t79\t139\t0\tSLIDING", "\t80\t139\t79\t139\t0\tSLIDING"));
		fails("wrong.trace.tsv", "sources.tsv", "wrong_zero");
		final CobalaminBoundaryReplayTest positional=freshCase("positional_case");
		positional.views(true);
		positional.write("trace.tsv", trace(true, false));
		positional.invoke("trace.tsv", "sources.tsv", "result");
		check(positional.read("result.replay_selection.tsv").contains("\t0.1\tSELECTED_DIAGNOSTIC_MARGIN\t"), "Independent positional ends use the same confidence and tie selection rule");
		row=positional.row("result.replay_endpoints.tsv", "cobalamin", "5prime", .1f);
		close(row.parseDouble(9), 0);
		row=positional.row("result.replay_endpoints.tsv", "cobalamin", "3prime", .1f);
		close(row.parseDouble(9), 0);
		final CobalaminBoundaryReplayTest stopFirst=freshCase("stop_first_case");
		stopFirst.views(false, true);
		stopFirst.write("trace.tsv", trace(false, false, true));
		stopFirst.invoke("trace.tsv", "sources.tsv", "result");
		check(stopFirst.read("result.replay_selection.tsv").contains("\t0.1\tSELECTED_DIAGNOSTIC_MARGIN\t"), "Stop-first replay must gate its first end before choosing the held second-end features");
		row=stopFirst.row("result.replay_endpoints.tsv", "cobalamin", "5prime", .1f);
		check(row.parseLong(10)==2 && row.parseLong(11)==0 && row.parseLong(12)==1, "The stop-first held context corrects two first ends and preserves the fixed end");
		close(row.parseDouble(9), 0);
		row=stopFirst.row("result.replay_endpoints.tsv", "cobalamin", "3prime", .1f);
		close(row.parseDouble(9), 0);
		check(stopFirst.read("result.replay_loci.tsv").contains("b\tr\tRF00174\tRF01482\t0.1\tREPLAYED\t1\t0\t0\t0\t"),
			"Reverse stop-first signed errors require the opposite held/moved dependency from start-first");
	}
	private CobalaminBoundaryReplayTest freshCase(final String name) throws Exception{
		final File directory=new File(root, name);
		check(directory.mkdir(), "Each fixture case must have a fresh namespace without overwriting evidence");
		final CobalaminBoundaryReplayTest result=new CobalaminBoundaryReplayTest(directory.getPath());
		result.write("sources.tsv", read("sources.tsv"));
		return result;
	}

	/** Source-qualified normalized views are explicit records, not a new matcher. */
	private void views(final boolean positional){views(positional, false);}
	private void views(final boolean positional, final boolean stopFirst){
		final String header="input_id\tcalls\tmatched\tmasked\tunmatched\n";
		write("off.genomes.tsv", header+"a\t1\t1\t0\t0\nb\t1\t1\t0\t0\nc\t1\t1\t0\t0\nempty\t0\t0\t0\t0\n");
		write("on.genomes.tsv", header+"a\t2\t1\t0\t1\nb\t1\t1\t0\t0\nc\t1\t1\t0\t0\nempty\t0\t0\t0\t0\n");
		final int aStart=stopFirst ? 80 : 81, aStop=stopFirst ? 140 : 141;
		final int bStart=stopFirst ? 101 : 100, bStop=stopFirst ? 161 : 160;
		final int off5=stopFirst ? 1 : 0, off3=stopFirst ? 0 : -1;
		final int on5=stopFirst ? 1 : -1, on3=stopFirst ? 1 : positional ? 0 : -1;
		final String refs="input_id\treference_id\tseqid\tstart\tstop\tstrand\tscored\n"
			+"a\tr\tshared\t"+aStart+'\t'+aStop+"\t+\ttrue\nb\tr\tshared\t"+bStart+'\t'+bStop+"\t-\ttrue\nc\tr\tshared\t1\t61\t+\ttrue\n";
		write("off.references.tsv", refs);
		write("on.references.tsv", refs);
		final String q="input_id\tseqid\tstart\tstop\tstrand\temitted_member\tmodel\toutcome\tscored_eligible\tmatched_reference_id\n";
		write("off.queries.tsv", q+"a\tshared\t81\t140\t+\tRF00174\ta\tTP\ttrue\tr\n"
			+"b\tshared\t101\t160\t-\tRF00174\tb\tTP\ttrue\tr\nc\tshared\t1\t60\t+\tRF01689\tc\tTP\ttrue\tr\n");
		final int plusStart=stopFirst ? 81 : 80, plusStop=positional || stopFirst ? 141 : 140;
		final int minusStart=positional || stopFirst ? 100 : 101, minusStop=stopFirst ? 160 : 161;
		write("on.queries.tsv", q+"a\tshared\t"+plusStart+'\t'+plusStop+"\t+\tRF00174\ta\tTP\ttrue\tr\n"
			+"a\tshared\t171\t200\t+\tRF00174\tu\tUNMATCHED_OTHER\tfalse\t.\n"
			+"b\tshared\t"+minusStart+'\t'+minusStop+"\t-\tRF01482\tb\tTP\ttrue\tr\nc\tshared\t1\t61\t+\tRF01689\tc\tTP\ttrue\tr\n");
		final String w="input_id\treference_id\temitted_member\tmodel\tseqid\tcall_start\tcall_stop\tref_start\tref_stop\tstrand\td5\td3\n";
		write("off.witnesses.tsv", w+"a\tr\tRF00174\ta\tshared\t81\t140\t"+aStart+'\t'+aStop+"\t+\t"+off5+'\t'+off3+"\n"
			+"b\tr\tRF00174\tb\tshared\t101\t160\t"+bStart+'\t'+bStop+"\t-\t"+off5+'\t'+off3+"\nc\tr\tRF01689\tc\tshared\t1\t60\t1\t61\t+\t0\t-1\n");
		write("on.witnesses.tsv", w+"a\tr\tRF00174\ta\tshared\t"+plusStart+'\t'+plusStop+'\t'+aStart+'\t'+aStop+"\t+\t"+on5+'\t'+on3+"\n"
			+"b\tr\tRF01482\tb\tshared\t"+minusStart+'\t'+minusStop+'\t'+bStart+'\t'+bStop+"\t-\t"+on5+'\t'+on3+"\n"
			+"c\tr\tRF01689\tc\tshared\t1\t61\t1\t61\t+\t0\t0\n");
	}

	private String trace(final boolean positional, final boolean unknown){return trace(positional, unknown, false);}
	private String trace(final boolean positional, final boolean unknown, final boolean stopFirst){
		final ByteBuilder b=new ByteBuilder(NcrnaBoundaryTrace.HEADER);
		for(int input=0; input<4; input++){b.append("I\t").append(input).tab().append(input==3 ? "/empty.fa" : input==0 ? "/a.fa" : input==1 ? "/b.fa" : "/c.fa").nl();}
		for(int input=0; input<3; input++){
			candidate(b, input, 0, true, unknown && input==0, positional, false, stopFirst);
			if(input==0){
				candidate(b, 0, 1, false, false, positional, false, stopFirst);
				candidate(b, 0, 2, true, false, positional, true, stopFirst);
			}
			b.append("T\t").append(input).tab().append(input).append("\t0\tshared\t240\t")
				.append(input==0 ? 3 : 1).tab().append(input==0 ? 3 : 1).tab().append(input==0 ? 2 : 1)
				.tab().append(unknown && input==0 ? 1 : 0).append("\t0\n");
		}
		return b.append("T\t3\t3\t0\tempty\t0\t0\t0\t0\t0\t0\n").toString();
	}
	private static void candidate(final ByteBuilder b, final int input, final int id, final boolean selected,
			final boolean unknown, final boolean positional, final boolean unmatched, final boolean stopFirst){
		final boolean edge=input==2;
		final int start=unmatched ? 170 : edge ? 0 : unknown ? 79 : 80, stop=unmatched ? 199 : edge ? 59 : 139;
		final int observedStart=unmatched ? 170 : edge ? 0 : stopFirst ? 80 : 79;
		final int observedStop=unmatched ? 199 : edge ? 60 : positional || stopFirst ? 140 : 139;
		final int finalStart=input==1 ? 239-observedStop : observedStart, finalStop=input==1 ? 239-observedStart : observedStop;
		final String family=input==0 ? "cobalamin_rf00174" : input==1 ? "cobalamin_rf01482" : "cobalamin_rf01689";
		b.append("C\t").append(input).tab().append(input).tab().append(id).append("\t0\tshared\t").append(input==1 ? 1 : 0)
			.tab().append(family).tab().append(unmatched ? "u" : input==0 ? "a" : input==1 ? "b" : "c")
			.append("\t240\t").append(start).tab().append(stop).tab();
		if(!unknown){b.append(edge ? 1 : 0);}
		b.tab().append(edge ? 1 : 0)
			.tab().append(unknown ? 0 : 1).tab().append(unknown || unmatched ? 0 : 1).tab().append(unknown || unmatched ? 0 : 1)
			.append(unknown ? "\t\t\t\t\t16\t1000\t" : "\t0\t1061997773\t1061997773\t1\t16\t1000\t")
			.append(observedStart).tab().append(observedStop).tab();
		if(selected){b.append(finalStart);}
		b.tab();
		if(selected){b.append(finalStop);}
		b.tab();
		if(selected){b.append(0);}
		b.tab().append(positional ? "1-HOT" : "SLIDING").append("\t1\t").append(unknown || unmatched ? 0 : edge ? 2 : 3)
			.append("\t1\t").append(selected ? 1 : 0).tab().append(selected ? "SELECTED" : "PATH_REJECTED")
			.append(unknown ? "\t\t\t\n" : "\t0\t239\t60\n");
		if(unknown || unmatched){return;}
		if(edge){window(b, input, id, 1, positional ? 3 : 0, start, stop, .1f, .2f, .9f);}
		else if(positional){
			window(b, input, id, 0, 3, start, stop, .9f, .2f, .8995f);
			window(b, input, id, 1, 3, start, stop, .1f, .2f, .9f);
		}else{
			if(stopFirst){
				window(b, input, id, 1, 0, start, stop, .8995f, .2f, .9f);
				window(b, input, id, 0, 2, start, stop+1, 1, 1, 1);
				window(b, input, id, 0, 1, start, stop, .9f, .2f, .1f);
			}else{
				window(b, input, id, 0, 0, start, stop, .9f, .2f, .8995f);
				window(b, input, id, 1, 2, start-1, stop, 1, 1, 1);
				window(b, input, id, 1, 1, start, stop, .1f, .2f, .9f);
			}
		}
	}
	/** Builds the frozen version-one encoding independently of the current writer. */
	private static String legacyTrace(final String current){
		final ByteBuilder b=new ByteBuilder();
		final ByteBuilder window=new ByteBuilder();
		for(final String line:current.split("\n")){
			final String[] fields=line.split("\t", -1);
			if(line.equals("#ncrna_boundary_trace\t3")){fields[1]="1";}
			else if(fields[0].equals("C")){
				if(fields[22].isEmpty()){fields[22]="2147483647";}
				if(fields[14].equals("0")){
					check(fields[34].isEmpty() && fields[35].isEmpty() && fields[36].isEmpty(),
						"Unknown current traces must omit all window metadata");
					fields[12]="-1";
					for(int i=17; i<=19; i++){fields[i]=Integer.toString(Float.floatToIntBits(Float.NaN));}
					fields[20]="0";
					for(int i=34; i<=36; i++){fields[i]="-1";}
				}
				if(fields[32].equals("0")){
					for(int i=25; i<=27; i++){fields[i]="-1";}
				}
			}else if(fields[0].equals("S")){
				final String[] scores=fields[8].split(",", -1);
				window.clear();
				for(int i=0; i<scores.length; i++){
					if(i>0){window.append(',');}
					window.append(scores[i].isEmpty() ? Integer.toString(Float.floatToIntBits(Float.NaN)) : scores[i]);
				}
				fields[8]=window.toString();
			}
			for(int i=0; i<fields.length; i++){
				if(i>0){b.tab();}
				b.append(fields[i]);
			}
			b.nl();
		}
		return b.toString();
	}
	private static void window(final ByteBuilder b, final int input, final int id, final int end, final int context,
			final int start, final int stop, final float left, final float current, final float right){
		b.append("S\t").append(input).tab().append(input).tab().append(id).tab().append(end).tab().append(context).tab().append(start).tab().append(stop)
			.tab().append(Float.floatToIntBits(left)).append(',').append(Float.floatToIntBits(current)).append(',').append(Float.floatToIntBits(right)).nl();
	}
	private void invoke(final String trace, final String sources, final String output){invoke(trace, sources, output, null);}
	private void invoke(final String trace, final String sources, final String output, final String margins){
		final String[] args=new String[margins==null ? 6 : 7];
		args[0]="in1="+path("off");
		args[1]="in2="+path("on");
		args[2]="out="+path(output);
		args[3]="radius=1";
		args[4]="trace="+path(trace);
		args[5]="bindings="+path(sources);
		if(margins!=null){args[6]="replaymargins="+margins;}
		CompareCobalaminBenchmarkGrades.main(args);
	}
	private void fails(final String trace, final String sources, final String output){fails(trace, sources, output, null);}
	private void fails(final String trace, final String sources, final String output, final String margins){
		boolean rejected=false;
		try{invoke(trace, sources, output, margins);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && !new File(path(output+".replay_selection.tsv")).exists(), "Bad provenance or zero reconstruction must fail before selection exists");
	}
	private LineParser1 row(final String name, final String cohort, final String end, final float margin) throws Exception{
		final LineParser1 p=new LineParser1('\t');
		for(final String line:read(name).split("\n")){
			p.set(line.getBytes(StandardCharsets.UTF_8));
			// Float.toString can emit exponents; the native decimal-only parser cannot read those labels.
			if(p.termEquals(cohort, 3) && p.termEquals(end, 4) && Float.parseFloat(p.parseString(2))==margin){return p;}
		}
		throw new AssertionError("Missing replay margin/cohort/end row");
	}
	private void write(final String name, final String value){try(Output output=new Output(path(name))){output.writer.print(value);}}
	private String read(final String name) throws Exception{return new String(Files.readAllBytes(new File(path(name)).toPath()), StandardCharsets.UTF_8);}
	private String path(final String name){return new File(root, name).getPath();}
	private static void close(final double actual, final double expected){check(Math.abs(actual-expected)<1e-6, "Literal RMS differs: "+actual+" expected="+expected);}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private final File root;
}
