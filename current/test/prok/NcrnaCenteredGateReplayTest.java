package test.prok;

import java.io.File;
import java.util.Arrays;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaCenteredGateReplay;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Finite paired-cache proof: internal objective rejects a harm-balanced gate,
 * outer truth would favor it but cannot choose again, and a centered control's
 * movement count exceeds the whole-move gate's attainable count.
 * @author Brian, Ganyu
 */
public final class NcrnaCenteredGateReplayTest{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaCenteredGateReplayTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, pair[0].toLowerCase(Locale.ROOT), pair.length>1 ? pair[1] : null), "Unknown finite replay-test option");
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null, "Finite replay requires only out=fresh_directory");
			final File root=new File(parser.out1);
			check(root.mkdir(), "Fresh finite replay root required");
			pure();
			run(root);
			pp.outstream.println("NCRNA_CENTERED_GATE_REPLAY_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void pure(){
		final int[] result=new int[2];
		check(NcrnaCenteredGateReplay.apply(30, 149, 40, 139, true, false, 200, 109, 200, result)==0
			&& Arrays.equals(result, new int[]{40, 149}), "A single inward move is legal with the other RAW endpoint");
		check(NcrnaCenteredGateReplay.apply(30, 149, 40, 139, true, true, 200, 109, 200, result)==3
			&& Arrays.equals(result, new int[]{30, 149}), "Two individually valid inward moves can jointly fail and must restore BOTH raw ends");
		check(NcrnaCenteredGateReplay.apply(30, 149, -1, 139, true, false, 200, 109, 200, result)==2
			&& Arrays.equals(result, new int[]{30, 149}), "An out-of-contig proposal cannot be clipped into an apparent improvement");
		check(NcrnaCenteredGateReplay.apply(30, 100, 40, 139, false, false, 200, 109, 200, result)==0
			&& Arrays.equals(result, new int[]{30, 100}), "An unchanged inherited short raw span remains counted without inventing a new move failure");
	}

	private static void run(final File root){
		final File scored=new File(root, "scorer_fixture");
		ScoreNcrnaCenteredGateTest.main(new String[]{"out="+scored.getPath(), "t=1"});
		final String original=new File(scored, "result/scores.tsv.gz").getPath();
		for(final boolean outer:new boolean[]{false, true}){
			for(int end=0; end<2; end++){rewrite(original, path(root, (outer ? "outer" : "internal")+end+".tsv.gz"), end, outer);}
		}
		new NcrnaCenteredGateReplay(arguments(root, false)).process();
		boolean bestNever=false, unavailable=false;
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(path(root, "selected/choices.tsv"))){
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				if(row[0]=='#'){continue;}
				p.set(row);
				if(p.termEquals("BEST", 0) && p.termEquals("5prime", 1)){
					check(p.parseInt(3)==1 && p.length(4)==0 && p.length(5)==0 && p.termEquals("NEVER", 10) && p.parseLong(6)==0,
						"The internal better-minus-worse tie must prefer zero harms, not outer improvements or lower RMS");
					bestNever=true;
				}
				if(p.termEquals("RATE_MARGIN0", 0) && p.termEquals("5prime", 1) && p.termEquals("HIGH_UNAVAILABLE", 2)){
					check(p.parseInt(3)==0 && p.length(4)==0 && p.length(5)==0 && p.length(6)==0 && p.length(8)==0 && p.length(9)==0 && p.length(10)==0 && p.parseLong(7)==3,
						"Three centered control moves cannot be forced from only two whole-move gate candidates");
					unavailable=true;
				}
			}
		}
		check(bestNever && unavailable, "Finite controls must actually exercise selection and unavailable matched-rate brackets");
		checkSummary(path(root, "selected/summary.tsv"), 4, 1, 1);
		checkCurves(path(root, "selected/curves.tsv"), false);
		new NcrnaCenteredGateReplay(arguments(root, true)).process();
		checkSummary(path(root, "outer_result/summary.tsv"), 16, 10, 0);
		checkCurves(path(root, "outer_result/curves.tsv"), true);
		for(final String directory:new String[]{"selected", "outer_result"}){
			for(final String name:new String[]{"summary.tsv", "curves.tsv"}){checkDefined(path(root, directory+"/"+name));}
		}
		checkDefined(path(root, "selected/choices.tsv"));
		legacy(root);
		versionTwo(root);
		invalidScores(root);
		zeroCutoff(root);
	}

	/** Saved version-one placeholders are accepted only under their original version. */
	private static void legacy(final File root){
		for(final boolean outer:new boolean[]{false, true}){
			for(int end=0; end<2; end++){
				final String name=(outer ? "outer" : "internal")+end+".tsv.gz";
				recodeCache(path(root, name), path(root, "legacy_"+name), true, null);
			}
		}
		final String[] selected=arguments(root, false);
		selected[1]="in="+path(root, "legacy_internal0.tsv.gz");
		selected[2]="in2="+path(root, "legacy_internal1.tsv.gz");
		selected[3]="out="+path(root, "legacy_selected");
		new NcrnaCenteredGateReplay(selected).process();
		same(path(root, "selected/summary.tsv"), path(root, "legacy_selected/summary.tsv"));
		final ByteBuilder text=new ByteBuilder();
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path(root, "selected/choices.tsv")); Output out=new Output(path(root, "legacy_choices.tsv"))){
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(p.termEquals("#NCRNA_CENTERED_GATE_CHOICES_V2", 0)){
					out.writer.println("#NCRNA_CENTERED_GATE_CHOICES_V1");
				}else if(row[0]=='#'){out.writer.print(text.clear().append(row).nl());}
				else if(p.termEquals("policy", 0)){
					out.writer.println("policy\tend\tside\tavailable\tcutoff_bits\tcutoff\tmonitor_moved\ttarget_moved\tmonitor_better\tmonitor_worse");
				}else{
					check(p.terms()==11, "The new choice schema must declare its distinct never-move mode");
					final String[] fields=new String[10];
					for(int i=0; i<10; i++){fields[i]=p.length(i)==0 ? "NA" : p.parseString(i);}
					if(p.termEquals("NEVER", 10)){
						fields[4]="2139095040";
						fields[5]="Infinity";
					}
					writeRow(out, text, fields);
				}
			}
		}
		final String[] outer=arguments(root, true);
		outer[1]="in="+path(root, "legacy_outer0.tsv.gz");
		outer[2]="in2="+path(root, "legacy_outer1.tsv.gz");
		outer[3]="out="+path(root, "legacy_outer_result");
		outer[outer.length-1]="choices="+path(root, "legacy_choices.tsv");
		new NcrnaCenteredGateReplay(outer).process();
		same(path(root, "outer_result/summary.tsv"), path(root, "legacy_outer_result/summary.tsv"));
		checkDefined(path(root, "legacy_outer_result/summary.tsv"));
	}

	/** The preceding cache has optional scores but a required numeric maximum. */
	private static void versionTwo(final File root){
		for(final boolean outer:new boolean[]{false, true}){
			for(int end=0; end<2; end++){
				final String name=(outer ? "outer" : "internal")+end+".tsv.gz";
				recodeCache(path(root, name), path(root, "v2_"+name), 2, null);
			}
			final String[] args=arguments(root, outer);
			args[1]="in="+path(root, "v2_"+(outer ? "outer" : "internal")+"0.tsv.gz");
			args[2]="in2="+path(root, "v2_"+(outer ? "outer" : "internal")+"1.tsv.gz");
			args[3]="out="+path(root, outer ? "v2_outer_result" : "v2_selected");
			new NcrnaCenteredGateReplay(args).process();
			final String original=outer ? "outer_result" : "selected";
			for(final String name:new String[]{"summary.tsv", "curves.tsv", "unavailable.tsv"}){
				same(path(root, original+"/"+name), path(root, "v2_"+original+"/"+name));
			}
			if(!outer){same(path(root, "selected/choices.tsv"), path(root, "v2_selected/choices.tsv"));}
		}
	}

	/** Missing and computed-nonfinite proposal evidence both fail before any result is created. */
	private static void invalidScores(final File root){
		final String[] values={"", "NA", "2143289344", "2139095040"};
		for(int i=0; i<values.length; i++){
			final String cache=path(root, "bad_score"+i+".tsv.gz");
			recodeCache(path(root, "internal0.tsv.gz"), cache, false, values[i]);
			final String[] args=arguments(root, false);
			args[1]="in="+cache;
			args[3]="out="+path(root, "bad_result"+i);
			boolean rejected=false;
			try{new NcrnaCenteredGateReplay(args).process();}
			catch(IllegalArgumentException expected){rejected=true;}
			catch(AssertionError expected){
				rejected=values[i].equals("NA") && expected.getMessage()!=null
					&& expected.getMessage().startsWith("Non-digit character with ASCII code ");
			}
			check(rejected && !new File(path(root, "bad_result"+i)).exists(),
				"An absent/NA/NaN/infinite gate score cannot be silently converted to a hold");
		}
	}

	/** Numeric zero remains an inclusive threshold, distinct from both absence and never-move. */
	private static void zeroCutoff(final File root){
		for(int end=0; end<2; end++){
			recodeCache(path(root, "internal"+end+".tsv.gz"), path(root, "zero"+end+".tsv.gz"), false, "0");
		}
		final String[] args=arguments(root, false);
		args[1]="in="+path(root, "zero0.tsv.gz");
		args[2]="in2="+path(root, "zero1.tsv.gz");
		args[3]="out="+path(root, "zero_result");
		new NcrnaCenteredGateReplay(args).process();
		checkSummary(path(root, "zero_result/summary.tsv"), 4, 1, 1);
		checkCurves(path(root, "zero_result/curves.tsv"), false, 0);
	}

	/** A finite fixture transformer, never an inference or scientific cache generator. */
	private static void recodeCache(final String input, final String output, final boolean legacy, final String score){
		recodeCache(input, output, legacy ? 1 : 3, score);
	}
	private static void recodeCache(final String input, final String output, final int version, final String score){
		check(version>=1 && version<=3, "The finite transformer covers exactly the three declared cache versions");
		final boolean legacy=version==1, numericMaximum=version<3;
		final ByteBuilder text=new ByteBuilder();
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(input); Output out=new Output(output)){
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(numericMaximum && p.termEquals("#NCRNA_CENTERED_GATE_SCORES_V3", 0)){
					out.writer.println("#NCRNA_CENTERED_GATE_SCORES_V"+version);
				}else if(legacy && p.termStartsWith("#Optional scores:", 0)){continue;}
				else if(numericMaximum && p.termEquals("dataset", 0)){out.writer.print(prok.ScoreNcrnaCenteredGate.cacheHeader(true).nl());}
				else if(row[0]=='#' || p.termEquals("dataset", 0)){out.writer.print(text.clear().append(row).nl());}
				else{
					check(p.terms()==53, "Only complete literal score rows may be recoded");
					final String[] fields=new String[53];
					for(int i=0; i<fields.length; i++){
						fields[i]=legacy && i>=30 && p.length(i)==0 ? "NA" : p.parseString(i);
					}
					if(numericMaximum && fields[17].isEmpty()){fields[17]="2147483647";}
					if(score!=null && p.termEquals("PROPOSAL", 23)){fields[30]=score;}
					writeRow(out, text, fields);
				}
			}
		}
	}

	private static void same(final String first, final String second){
		try(Input a=new Input(first); Input b=new Input(second)){
			for(byte[] row; (row=a.reader.nextLine())!=null;){
				check(Arrays.equals(row, b.reader.nextLine()), "Legacy normalization must preserve the complete finite result");
			}
			check(b.reader.nextLine()==null, "Legacy normalization cannot add result rows");
		}
	}

	private static void checkDefined(final String path){
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path)){
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				for(int i=0; i<p.terms(); i++){
					check(!p.termEquals("NA", i) && !p.termEquals("NaN", i) && !p.termEquals("Infinity", i)
						&& !p.termEquals("-Infinity", i), "New result fields must omit nonexistent numbers and never emit nonfinite measurements");
				}
			}
		}
	}

	/** The literal internal fixture has two raw-exact starts, one of which
	 * the gate harms. Outer truth leaves only its five fixed starts exact. */
	private static void checkCurves(final String path, final boolean outer){
		checkCurves(path, outer, -.75f);
	}
	private static void checkCurves(final String path, final boolean outer, final float cutoff){
		boolean never=false;
		boolean always=false;
		try(Input in=new Input(path)){
			in.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				check(p.terms()==20 && p.parseInt(5)==(outer ? 16 : 4)
					&& p.termEquals(outer ? "outer_panel" : "fold4", 17) && p.termEquals("SINGLE_OTHER_RAW", 18),
					"Every curve point retains the full frozen cohort, panel, scope and presentation schema");
				check(p.termEquals(outer ? "OUTER_DESCRIPTIVE_NO_SELECTION" : "INTERNAL_MONITOR_SELECTION", 15),
					"Outer tradeoff curves must never masquerade as cutoff selection");
				if(!p.termEquals("5prime", 2)){continue;}
				if(p.termEquals("NEVER", 19)){
					check(p.parseInt(9)==0 && p.parseFloat(14)==0 && p.parseInt(13)==(outer ? 5 : 2)
						&& p.parseString(16).contains("BEST:ALL_TIED_BEST"), "Never-move retains raw exactness and the monitor-selected marker on both populations");
					never=true;
				}else{
					check(p.parseFloat(4)==cutoff && p.termEquals("THRESHOLD", 19) && p.parseInt(9)==(outer ? 10 : 2)
						&& p.parseFloat(14)==(outer ? .625f : .5f) && p.parseInt(13)==(outer ? 5 : 1),
						"Always-move curve reports actual physical movement and independently specified exact counts");
					always=true;
				}
			}
		}
		check(never && always, "Both ends of the cutoff tradeoff must be present, even though only never-move was chosen internally");
	}

	/** Creates paired literal caches from the separately tested scorer fixture.
	 * The extra fourth record in fold4 has C=raw+2 and a neural reversal -2:
	 * the gate holds raw, while the ordinary head veto retains the shifted C. */
	private static void rewrite(final String input, final String output, final int end, final boolean outer){
		final ByteBuilder text=new ByteBuilder();
		final LineParser1 p=new LineParser1('\t');
		String[] extra=null;
		try(Input in=new Input(input); Output out=new Output(output)){
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(p.termEquals("#COMPLETE", 0)){break;}
				if(row[0]=='#'){
					if(p.termEquals("#end", 0)){out.writer.println("#end\t"+(end==0 ? "5prime" : "3prime"));}
					else if(p.termEquals("#partition", 0)){out.writer.println("#partition\t"+(outer ? "outer" : "internal"));}
					else if(end==1 && p.termEquals("#network", 0)){out.writer.println("#network\t/fixture/stop_gate.net");}
					else{out.writer.print(text.clear().append(row).nl());}
					continue;
				}
				if(p.termEquals("dataset", 0)){out.writer.print(text.clear().append(row).nl());continue;}
				check(p.terms()==53, "The preceding finite scorer supplies the complete real cache layout");
				final String[] fields=new String[53];
				for(int i=0; i<53; i++){fields[i]=p.parseString(i);}
				if(p.parseInt(1)==4 && p.parseInt(2)==1){extra=fields.clone();}
				transform(fields, end, outer);
				writeRow(out, text, fields);
			}
			check(extra!=null, "Cancellation fixture must be based on the actual fold4 cohort");
			extra[2]="4";
			extra[5]="cancelled";
			extra[23]="SAME_AS_RAW_HOLD";
			extra[24]="0";
			extra[25]="100";
			extra[26]="-2";
			extra[28]="0";
			extra[30]="";
			for(int i=31; i<53; i++){extra[i]="0";}
			extra[39]=Integer.toString(Float.floatToRawIntBits(1));
			extra[52]=Integer.toString(Float.floatToRawIntBits(2));
			transform(extra, end, outer);
			writeRow(out, text, extra);
			out.writer.print(text.clear().append("#COMPLETE\t16\t").append(end==0 ? 10 : 0).nl());
		}
	}

	private static void transform(final String[] fields, final int end, final boolean outer){
		if(fields[23].equals("PROPOSAL")){fields[30]=Integer.toString(Float.floatToRawIntBits(-.75f));}
		if(outer){
			fields[0]="outer_panel";
			fields[1]="-1";
			if(!fields[11].equals("1")){
				fields[14]="103";
				if(fields[23].equals("PROPOSAL")){fields[28]="1";}
			}
		}
		if(end==1){
			fields[22]="3prime";
			fields[23]="SAME_AS_RAW_HOLD";
			fields[24]="0";
			fields[25]=fields[13];
			fields[26]="0";
			fields[27]="0";
			fields[28]="0";
			fields[30]="";
			for(int i=31; i<53; i++){fields[i]="0";}
		}
	}
	private static void writeRow(final Output out, final ByteBuilder text, final String[] fields){
		text.clear();
		for(int i=0; i<fields.length; i++){
			if(i>0){text.tab();}
			text.append(fields[i]);
		}
		out.writer.print(text.nl());
	}
	private static void checkSummary(final String path, final int expected, final int allBetter, final int allWorse){
		boolean best=false, all=false;
		try(Input in=new Input(path)){
			in.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				check(p.terms()==48 && p.parseInt(10)==expected && p.parseInt(20)+p.parseInt(21)+p.parseInt(22)==expected,
					"Every policy retains the exact matched-call denominator and complete metric schema");
				if(!p.termEquals("5prime", 9) || !p.termEquals("PAIRED_JOINT", 8)){continue;}
				if(p.termEquals("BEST", 3)){
					check(p.parseInt(20)==0 && p.parseInt(21)==0 && p.parseInt(23)==0, "The frozen never-move choice must remain unchanged on outer truth that favors movement");
					best=true;
				}
				if(p.termEquals("ALL_WHOLE_PROPOSALS", 3)){
					check(p.parseInt(20)==allBetter && p.parseInt(21)==allWorse && p.parseInt(23)==allBetter+allWorse,
						"The same proposals have explicit different outcomes internally and externally; the outer labels cannot select a new cutoff");
					all=true;
				}
			}
		}
		check(best && all, "Finite summary must actually contain both selected and all-proposal controls");
	}
	private static String[] arguments(final File root, final boolean outer){
		final String[] common={"mode="+(outer ? "evaluate" : "select"), "in="+path(root, (outer ? "outer" : "internal")+"0.tsv.gz"),
			"in2="+path(root, (outer ? "outer" : "internal")+"1.tsv.gz"), "out="+path(root, outer ? "outer_result" : "selected"),
			"member=RF00174", "seed=59", "expected=16", "t=1", "unpigz=f", "pigz=f"};
		if(!outer){return common;}
		final String[] args=Arrays.copyOf(common, common.length+1);
		args[common.length]="choices="+path(root, "selected/choices.tsv");
		return args;
	}
	private static String path(final File root, final String name){return new File(root, name).getPath();}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
