package test.prok;

import java.io.File;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlHoldout;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.PrepareNcrnaCmTraining;
import prok.TraceNcrnaCmTraining;
import shared.Shared;

/** Literal paired AC grammar: full paths, reverse strand, empty END state,
 * unexpanded glocal state anchors, disjoint pilot/remainder ranges and plan-only.
 * Expected intervals are written independently of the scorer's path export.
 * @author Brian Bushnell, Raiden
 */
public final class TraceNcrnaCmTrainingTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, TraceNcrnaCmTrainingTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, split[0], split.length>1 ? split[1] : null), "Use out=fresh_directory");
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append && new File(parser.out1).mkdir(), "Training trace fixture owns fresh data files");
			// The exported state anchors also feed existing band layouts. Retain
			// their independent exhaustive linear/paired/branching comparisons.
			test.aligner.CovarianceModelCykTest.main(new String[0]);
			run(new File(parser.out1).getAbsolutePath());
			System.err.println("NCRNA_CM_TRAINING_PATH_FIXTURE_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fixtures           ----------------*/
	/*--------------------------------------------------------------*/

	private static void run(final String root){
		// A literal %25 must survive one decode, not recursively become '%'.
		write(root+"/a.fa", ">C%25 original contig\nTTACGG\n");
		write(root+"/b.fa", ">C%25\nTTGTGG\n");
		write(root+"/sources.tsv", NcrnaDdlHoldout.PANEL_HEADER+'\n'
			+"0\tassembly_a\ta\t"+root+"/a.fa\t/unused_a.gff\n"
			+"1\tassembly_b\tb\t"+root+"/b.fa\t/unused_b.gff\n");
		write(root+"/calls.tsv", PrepareNcrnaCmTraining.CALL_HEADER+'\n'
			+"0\tassembly_a\t0\ta\tC%2525\t2\t5\t+\tRF00174\tcobalamin_RF00174_literal\tTP\ttrue\tr_a\n"
			+"1\tassembly_b\t1\tb\tC%2525\t2\t5\t-\tRF01482\tcobalamin_RF01482_literal\tUNMATCHED\ttrue\t.\n"
			+"2\tassembly_a\t2\ta\tC%2525\t2\t5\t+\tRF01689\tcobalamin_RF01689_literal\tTP\ttrue\tr_c\n");
		for(final String member:MEMBERS){
			write(root+'/'+member+".fa", ">cobalamin_"+member+"_literal\nAC\n");
			if(member.equals("RF01689")){continue;}
			// Exact two-residue AC grammar from the independently enumerated
			// endpoint fixtures. Other pairs are impossible, not low-score ties.
			write(root+'/'+member+".cm", "INFERNAL1/a\nNAME LiteralAC\nACC "+member
				+"\nSTATES 3\nNODES 3\nCLEN 2\nW 6\nALPH RNA\nNULL 0 0 0 0\nELSELF -1\nEFP7GF 0 1\nGA 100\nQDBBETA1 1e-7\nQDBBETA2 1e-15\nCM\n"
				+"[ ROOT 0 ] - - - - - -\nS 0 -1 0 1 1 0 0 6 6 0\n"
				+"[ MATP 1 ] - - - - - -\nMP 1 0 1 2 1 0 0 6 6 0 * 0 * * * * * * * * * * * * * *\n"
				+"[ END 2 ] - - - - - -\nE 2 1 1 -1 0 0 0 0 0\n//\n");
			write(root+'/'+member+".sto", "# STOCKHOLM 1.0\ncobalamin_"+member+"_literal AC\n#=GC RF xx\n//\n");
		}
		invoke(root, "whole", 0, 3, 1, false);
		checkCalls(root+"/whole/calls.tsv.gz", new int[]{0, 1});
		checkPaths(root+"/whole/paths.tsv.gz", 2);
		invoke(root, "pilot", 0, 1, 1, false);
		checkCalls(root+"/pilot/calls.tsv.gz", new int[]{0});
		checkPaths(root+"/pilot/paths.tsv.gz", 1);
		invoke(root, "remainder", 1, 2, 2, false);
		checkCalls(root+"/remainder/calls.tsv.gz", new int[]{1});
		checkPaths(root+"/remainder/paths.tsv.gz", 1);
		invoke(root, "plan", 0, 3, 1, true);
		check(new File(root+"/plan/plan.tsv").isFile() && !new File(root+"/plan/paths.tsv.gz").exists(), "Planning must allocate no exact decks or path outputs");
	}

	private static void checkCalls(final String path, final int[] expected){
		final LineParser1 p=new LineParser1('\t');
		final boolean[] seen=new boolean[3];
		int count=0;
		try(Input input=new Input(path)){
			check(new String(input.reader.nextLine(), java.nio.charset.StandardCharsets.UTF_8).equals(TraceNcrnaCmTraining.CALL_HEADER), "Trace call schema must retain every range/source and timing field");
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				check(p.terms()==23, "Call path summary must have exactly23 declared fields");
				final int index=p.parseInt(0);
				check(index>=0 && index<2 && !seen[index], "Third member is not searched and exact call identities cannot repeat");
				seen[index]=true;
				count++;
				check(p.termEquals("C%2525", 4), "Literal percent-bearing FASTA ID must round-trip the saved TSV framing exactly once");
				check(p.termEquals(index==0 ? "+" : "-", 7) && p.parseInt(8)==1 && p.parseInt(9)==6 && p.parseInt(10)==6,
					"Both literal strands use the exact full six-base oriented window");
				check(p.termEquals("OK", 11) && p.parseFloat(12)==4 && p.parseInt(13)==3 && p.parseInt(14)==4,
					"The unique AC pair has probability1 versus null1/16, hence4bits at oriented bases3..4, belowGA100");
				check(p.termEquals("true", 15) && p.parseInt(16)==0 && p.parseLong(17)>=0 && p.parseLong(18)>=0
					&& p.parseLong(20)==147, "Literal path lies in its zero-radius glocal anchors and exact matrix has3x7x7 cells");
				check(p.termEquals(index==0 ? "TP" : "UNMATCHED", 21), "Unmatched accepted calls must retain path evidence");
			}
		}
		check(count==expected.length, "Every selected nonthird-member call needs exactly one outcome");
		for(final int index:expected){check(seen[index], "Missing literal selected call "+index);}
	}

	private static void checkPaths(final String path, final int calls){
		final LineParser1 p=new LineParser1('\t');
		int rows=0;
		try(Input input=new Input(path)){
			input.reader.nextLine();
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				check(p.terms()==13, "Every trace step must preserve state/children and four raw anchor bounds");
				final int state=p.parseInt(3), left=state==2 ? 3 : 2, right=state==2 ? 3 : 4;
				check(state>=0 && state<3 && p.parseInt(4)==left && p.parseInt(5)==right,
					"Root/paired interval is[2,4); END consumes the empty[3,3) interval");
				check(p.parseInt(9)==left && p.parseInt(10)==left && p.parseInt(11)==right && p.parseInt(12)==right,
					"Exported zero-radius state anchors must equal the independent literal consumed-prefix bounds");
				rows++;
			}
		}
		check(rows==3*calls, "Each finite paired grammar has root, paired-emitter and empty-END steps");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Helpers           ----------------*/
	/*--------------------------------------------------------------*/

	private static void invoke(final String root, final String label, final int first, final int limit, final int threads, final boolean plan){
		assert(first>=0 && limit>0 && first+limit<=3) : "Fixture ranges must lie on the same declared three-call axis";
		TraceNcrnaCmTraining.main(new String[]{"in="+root+"/calls.tsv", "bindings="+root+"/sources.tsv",
			"refs="+root+"/RF00174.fa,"+root+"/RF01482.fa,"+root+"/RF01689.fa",
			"maps="+root+"/RF00174.sto,"+root+"/RF01482.sto", "models="+root+"/RF00174.cm,"+root+"/RF01482.cm",
			"expectedcalls=3", "expectedinputs=2", "first="+first, "limit="+limit, "padding=1", "maxcells=1000",
			"plan="+plan, "out="+root+'/'+label, "t="+threads, "trd=f", "itn=f", "utot=f", "dotdashxton=f", "lctn=f"});
	}
	private static void write(final String path, final String text){
		assert(path!=null && text!=null) : "Literal fixture files require explicit paths and complete owned content";
		try(Output output=new Output(path)){output.writer.print(text);}
	}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
}
