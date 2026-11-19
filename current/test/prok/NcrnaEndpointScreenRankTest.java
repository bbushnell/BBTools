package test.prok;

import java.io.File;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaEndpointScreen;
import prok.NcrnaTextIO.Input;
import shared.Shared;
import structures.ByteBuilder;

/** Literal saved-count ranking oracle: validation first, then ratio3, smallerK,
 * retained placement ties, unequal denominators, unavailable and malformed rows.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaEndpointScreenRankTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaEndpointScreenRankTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown rank fixture argument: "+arg);
			}
			check(parser.out1!=null && new File(parser.out1).mkdir(), "Supply out=fresh_fixture_directory");
			new NcrnaEndpointScreenRankTest(parser.out1).run();
			pp.outstream.println("NCRNA_ENDPOINT_SCREEN_RANK_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private NcrnaEndpointScreenRankTest(final String root_){root=root_;}

	private void run(){
		final ByteBuilder data=new ByteBuilder(HEADER);
		geometry(data, "START", 5, 0, 10, 6, 1);
		geometry(data, "START", 7, 0, 20, 12, 4);
		geometry(data, "START", 6, 0, 20, 12, 4);
		geometry(data, "START", 6, 1, 20, 12, 4);
		geometry(data, "START", 5, 2, 0, 0, 0);
		geometry(data, "STOP", 5, 0, 10, 7, 2);
		geometry(data, "STOP", 7, 0, 20, 12, 12);
		write("curves.tsv", data.toString());
		invoke("curves.tsv", "result");
		int choices=0, starts=0, stops=0, startPlacements=0;
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(path("result.choices.tsv"))){
			check(input.reader.nextLine()!=null, "Choice schema must exist");
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				check(p.terms()==20 && p.parseInt(1)==1 && p.termEquals("TOP_OR_TIED_COUNT_SCREEN", 16), "Only complete valid top/tied rows belong in choices");
				if(p.termEquals("START", 0)){
					final int outside=p.parseInt(4);
					check(p.parseInt(2)==6 && outside>=0 && outside<=1 && (startPlacements&(1<<outside))==0 && p.parseLong(6)==20 && p.parseLong(7)==12
						&& p.parseLong(9)==4 && p.parseLong(11)==40 && p.parseLong(12)==8,
						"Ratio3 beats smallerK5, then K6 beats tied K7; both K6 placements survive despite inferior training fractions");
					startPlacements|=1<<outside;
					check(Math.abs(p.parseDouble(18)+.4)<1e-7 && Math.abs(p.parseDouble(19)+.15)<1e-7,
						"Gap is training minus validation with separate denominators:8/40-12/20 and2/40-4/20");
					starts++;
				}else{
					check(p.termEquals("STOP", 0) && p.parseInt(2)==5 && p.parseLong(6)==10 && p.parseLong(7)==7,
						"Validation ratio2 fraction0.7 outranks0.6 even when the latter has a stronger ratio3 count");
					stops++;
				}
				choices++;
			}
		}
		check(choices==3 && starts==2 && stops==1 && startPlacements==3, "The independent oracle has exactly two distinct tied START placements and one STOP choice");
		int perK=0, tiedAtSix=0;
		try(Input input=new Input(path("result.byk.tsv"))){
			input.reader.nextLine();
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				check(p.terms()==20 && p.parseInt(1)==1 && p.termEquals("TOP_OR_TIED_AT_K_COUNT_SCREEN", 16),
					"Per-K output uses the same native counts, signed gaps and full ordering rule");
				if(p.termEquals("START", 0) && p.parseInt(2)==6){tiedAtSix++;}
				perK++;
			}
		}
		check(perK==6 && tiedAtSix==2, "Every K keeps its best placements, including both tied K6 START rows");
		final ByteBuilder empty=new ByteBuilder(HEADER);
		geometry(empty, "START", 5, 0, 0, 0, 0);
		geometry(empty, "STOP", 5, 0, 0, 0, 0);
		write("empty.tsv", empty.toString());
		invoke("empty.tsv", "empty");
		try(Input input=new Input(path("empty.choices.tsv"))){
			input.reader.nextLine();
			for(int end=0; end<2; end++){
				final byte[] line=input.reader.nextLine();
				check(line!=null, "Both unavailable endpoints remain explicit");
				p.set(line);
				check(p.termEquals("NA", 1) && p.termEquals("NA", 2) && p.termEquals("UNAVAILABLE_NO_VALIDATION", 16)
					&& p.termEquals("NA", 18) && p.termEquals("NA", 19),
					"Zero validation cannot invent a rank or choose a K");
			}
			check(input.reader.nextLine()==null, "Unavailable accounting emits one row per end");
		}
		write("duplicate.tsv", data.toString()+data.toString().substring(HEADER.length()));
		fails("duplicate.tsv", "bad_duplicate");
		write("mismatched.tsv", data.toString().replace("START\tVALIDATION\t16\t3.0\t10\t1\t0\t9\t", "START\tVALIDATION\t16\t3.0\t11\t1\t0\t10\t"));
		fails("mismatched.tsv", "bad_cohort");
		write("missing_validation.tsv", removeLine(data.toString(), "5\t5\t0\tSTART\tVALIDATION\t16\t3.0\t"));
		fails("missing_validation.tsv", "bad_missing_validation");
		write("nonranking_only.tsv", data.toString()
			+"8\t8\t0\tSTART\tTRAIN_LOO\t16\t1.5\t40\t2\t0\t38\t0.05\t0\n"
			+"8\t8\t0\tSTART\tVALIDATION\t16\t1.5\t20\t1\t0\t19\t0.05\t0\n");
		fails("nonranking_only.tsv", "bad_missing_target_geometry");
	}
	private static String removeLine(final String text, final String prefix){
		final int from=text.indexOf(prefix), to=from<0 ? -1 : text.indexOf('\n', from);
		check(from>=0 && to>from, "The missing-row oracle must remove its literal validation record");
		return text.substring(0, from)+text.substring(to+1);
	}

	/** Counts and expected winners above are hand chosen, not generated by the ranker. */
	private static void geometry(final ByteBuilder b, final String end, final int k, final int outside,
			final int n, final int two, final int three){
		assert(three<=two && two<=n) : "Literal ratio3 decisions must be a subset of ratio2 decisions";
		for(int phase=0; phase<2; phase++){
			for(int ratio=0; ratio<2; ratio++){
				final int eligible=phase==0 ? 40 : n;
				final int correct=phase==0 ? ratio==1 ? 2 : k==5 ? 30 : 8 : ratio==0 ? two : three;
				b.append(k).tab().append(k-outside).tab().append(outside).tab().append(end)
					.tab().append(phase==0 ? "TRAIN_LOO" : "VALIDATION").append("\t16\t").append(ratio==0 ? "2" : "3.0")
					.tab().append(eligible).tab().append(correct).append("\t0\t").append(eligible-correct).tab();
				if(eligible==0){b.append("NA");}else{b.append(correct/(double)eligible, 8);}
				b.append("\t0\n");
			}
		}
	}
	private void invoke(final String input, final String output){
		NcrnaEndpointScreen.main(new String[]{"mode=rank", "in="+path(input), "out="+path(output), "t=1"});
	}
	private void fails(final String input, final String output){
		boolean rejected=false;
		try{invoke(input, output);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && !new File(path(output+".choices.tsv")).exists(), "Malformed native rows must fail before choosing geometry");
	}
	private void write(final String file, final String text){NcrnaLearningPopulationTest.write(path(file), text);}
	private String path(final String name){return new File(root, name).getPath();}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private final String root;
	private static final String HEADER="k\tinside\toutside\tboundary\tpartition\tradius\tratio\teligible_loci\tdecisive_true\tdecisive_other\tabstain\ttrue_fraction\tother_fraction\n";
}
