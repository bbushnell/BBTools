package test.prok;

import java.io.File;
import java.util.ArrayList;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import parse.LineParser1;
import prok.NcrnaBoundaryComparison;
import prok.NcrnaBoundaryComparison.Change;
import prok.NcrnaBoundaryComparison.Result;
import prok.NcrnaLocusComparison.Locus;
import prok.NcrnaNativeTrainer.Method;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.ReduceNcrnaModels;
import shared.Shared;

/** Literal detection, strand, eligibility and tie-order oracles; no biological data. */
public final class NcrnaBoundaryComparisonTest {

	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, null, false);
		try(AutoCloseable log=()->Shared.closeStream(pp.outstream)){
			run(pp.args);
		}
	}

	/** Main owns the selected log through parsing and literal comparison checks. */
	private static void run(final String[] args){
		final Parser parser=new Parser();
		for(String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			if(!parser.parse(arg, split[0].toLowerCase(), split.length>1 ? split[1] : null)){
				throw new IllegalArgumentException("Unknown comparison fixture argument: "+arg);
			}
		}
		if(parser.out1==null){throw new IllegalArgumentException("Specify a fresh out= directory");}
		final File root=new File(parser.out1);
		if(!root.mkdir()){throw new IllegalArgumentException("Fixture directory must be fresh");}
		final ArrayList<Locus> truth=loci(new int[][]{{100,199,0},{400,499,1},{700,799,0},{1000,1099,0}});
		final Result off=grade(truth, new int[][]{{98,199,0},{400,502,1},{700,799,0}}, null);
		assert(off.matched==3 && off.unmatched==0 && off.missed==1 && off.exactStart==1
			&& off.exactStop==3 && off.exactBoth==1 && off.precision()==1 && off.recall()==.75)
			: "Literal baseline has three detections, one both-exact locus and one miss";
		assert(NcrnaBoundaryComparison.offset(off.calls.get(0), truth.get(0), true)==-2
			&& NcrnaBoundaryComparison.offset(off.calls.get(1), truth.get(1), true)==-3)
			: "Both strands use upstream-negative transcription-oriented offsets";
		final Result classifier=grade(truth, new int[][]{{100,199,0},{400,499,1},{700,799,0}}, off);
		final Result regression=grade(truth, new int[][]{{100,199,0},{401,499,1},{700,799,0}}, off);
		assert(NcrnaBoundaryComparison.offset(regression.calls.get(1), truth.get(1), false)==-1)
			: "A minus-strand genomic-left shift must use the transcriptional 3-prime sign";
		assert(new Change(off, classifier).exactGain==2 && new Change(off, classifier).lost==0)
			: "Correcting two endpoints changes exactness without changing detections";
		assert(NcrnaBoundaryComparison.choose(off, classifier, regression, 0, .95, Method.REGRESSION)==Method.CLASSIFIER)
			: "Greater raw both-exact gain precedes tie preference";
		assert(NcrnaBoundaryComparison.choose(off, classifier, classifier, 0, .95, Method.REGRESSION)==Method.REGRESSION
			&& NcrnaBoundaryComparison.choose(off, classifier, classifier, 0, .95, Method.CLASSIFIER)==Method.CLASSIFIER)
			: "Declared tie preference decides otherwise equal eligible candidates";
		final Result loss=grade(truth, new int[][]{{100,199,0},{400,499,1},{1000,1099,0}}, off);
		final Change changed=new Change(off, loss);
		assert(changed.lost==1 && changed.gained==1 && changed.exactGain==2
			&& !changed.eligible(0, .95) && changed.eligible(1, .95))
			: "Equal TP counts do not erase a protected loss or a distinct gain";
		assert(NcrnaBoundaryComparison.choose(off, classifier, loss, 1, .95, Method.REGRESSION)==Method.CLASSIFIER)
			: "Fewer lost references precedes tie preference at equal exact gain";
		final Result more=grade(truth, new int[][]{{100,199,0},{400,499,1},{700,799,0},{999,1099,0}}, off);
		assert(NcrnaBoundaryComparison.choose(off, classifier, more, 0, .95, Method.CLASSIFIER)==Method.REGRESSION)
			: "More detections precedes tie preference at equal exact gain and losses";
		final Result duplicates=grade(truth, new int[][]{{100,199,0},{400,499,1},{700,799,0},{100,199,0}}, off);
		assert(duplicates.matched==3 && duplicates.unmatched==1 && !new Change(off, duplicates).eligible(0, .95))
			: "Duplicate predictions consume no second reference and lower precision";
		final Change lowerPrecision=new Change(off, duplicates);
		assert(lowerPrecision.eligible(0, 0) && lowerPrecision.eligible(0, .75)
			&& !lowerPrecision.eligible(0, .750001))
			: "Explicit lower precision limits apply exactly to the measured three-of-four candidate";
		assert(!changed.eligible(0, 0) && !new Change(off, off).eligible(0, 0))
			: "A zero precision floor must retain the detection-loss and positive endpoint-gain requirements";
		for(final double invalid:new double[]{-.01, 1.01, Double.NaN, Double.POSITIVE_INFINITY}){
			boolean rejected=false;
			try{lowerPrecision.eligible(0, invalid);}
			catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("minprecision");}
			assert(rejected) : "A selected precision floor must be finite and in [0,1]: "+invalid;
		}
		final Result empty=grade(truth, new int[0][], off);
		assert(Double.isNaN(empty.precision()) && NcrnaBoundaryComparison.choose(off, empty, off, 0, .95, Method.REGRESSION)==null)
			: "Undefined precision and zero exact gain cannot displace OFF";
		precisionBoundary();
		boolean wrongAxis=false;
		try{NcrnaBoundaryComparison.grade(new ArrayList<Locus>(truth), classifier.calls, off);}
		catch(IllegalArgumentException expected){wrongAxis=true;}
		assert(wrongAxis) : "A different reference axis must not be compared by coincidental row counts";
		NcrnaBoundaryComparison.writeSummary(new File(root, "summary.tsv").getPath(), off, classifier, regression);
		// Literal three-of-four, empty and complete calls exercise fractional/zero/unit report rates.
		NcrnaBoundaryComparison.writeSummary(new File(root, "format.tsv").getPath(), off, empty, more);
		NcrnaBoundaryComparison.writeTransitions(new File(root, "transitions.tsv").getPath(), off, loss);
		savedComparison(root);
		System.out.println("Boundary comparison passed: detection conservation, both strands, losses/gains, duplicates, precision boundary, OFF and tie selection");
	}

	/** The public saved mode must expose the existing rule without needing any sequence/model files. */
	private static void savedComparison(final File root){
		write(new File(root, "saved.truth.tsv"), "seqid\tfamily\tstart\tstop\tstrand\n"
			+"chr\tfixture\t100\t199\t+\nchr\tfixture\t400\t499\t-\n"
			+"chr\tfixture\t700\t799\t+\nchr\tfixture\t1000\t1099\t+\n");
		write(new File(root, "saved.before.gff"), savedCall(98, 199, '+')
			+savedCall(400, 502, '-')+savedCall(700, 799, '+'));
		write(new File(root, "saved.after.gff"), savedCall(100, 399, '+')
			+savedCall(400, 499, '-')+savedCall(1000, 1099, '+')+savedCall(1000, 1099, '+'));
		final File compared=new File(root, "saved_comparison");
		compareSaved(root, "saved.after.gff", compared);
		checkSavedSummary(compared, new int[]{3, 1, 1, 3, 2, 2, 1, 1, 1});
		final File same=new File(root, "saved_self");
		compareSaved(root, "saved.before.gff", same);
		checkSavedSummary(same, new int[]{3, 0, 1, 1, 3, 1, 0, 0, 0});
		final LineParser1 lp=new LineParser1('\t');
		try(Input input=new Input(new File(compared, "half_shorter.transitions.tsv").getPath())){
			lp.set(input.reader.nextLine());
			assert(lp.termEquals("after_3prime_offset", 12)) : "The saved comparison must retain the native transition schema";
			lp.set(input.reader.nextLine());
			assert(lp.parseInt(9)==-2 && lp.parseInt(11)==0 && lp.parseInt(12)==200 && lp.termEquals("retained", 13))
				: "A long candidate matches half of the shorter span, with its actual 200-base 3-prime error";
			int rows=1, lost=0, gained=0;
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				lp.set(line);
				if(lp.termEquals("lost", 13)){lost++;}
				if(lp.termEquals("gained", 13)){gained++;}
				rows++;
			}
			assert(rows==4 && lost==1 && gained==1) : "Per-reference output must expose one loss and one gain despite equal TP totals";
		}
		for(final String option:new String[]{"target=0.98", "callerverify=f", "maxmodelloss=0", "ref=absent.fa", "rcomp=f"}){
			final File rejected=new File(root, "saved_rejected");
			boolean failed=false;
			try{compareSaved(root, "saved.after.gff", rejected, option);}
			catch(IllegalArgumentException expected){failed=expected.getMessage().contains("mode=compare accepts saved calls and truth only");}
			assert(failed && !rejected.exists()) : "Saved mode must reject explicit unrelated work before creating output: "+option;
		}
		assert(compared.list().length==2 && same.list().length==2)
			: "Read-only comparison produces only two reports, never a caller directory or model library";
	}

	/** Independent literal expected counts include a duplicate FP and one gain replacing a loss. */
	private static void checkSavedSummary(final File root, final int[] expected){
		final LineParser1 lp=new LineParser1('\t');
		try(Input input=new Input(new File(root, "half_shorter.summary.tsv").getPath())){
			lp.set(input.reader.nextLine());
			assert(lp.termStartsWith("# overlap_rule=HALF_SHORTER", 0)) : "Saved reports must explicitly identify the native overlap rule";
			input.reader.nextLine();
			lp.set(input.reader.nextLine());
			assert(lp.termEquals("method", 0) && lp.termEquals("both_exact_gain", 11)) : "Summary has the documented native count schema";
			lp.set(input.reader.nextLine());
			assert(lp.termEquals("baseline", 0) && lp.parseInt(1)==3 && lp.parseInt(2)==0 && lp.parseInt(3)==1)
				: "Literal saved baseline detects three of four references without extras";
			lp.set(input.reader.nextLine());
			assert(lp.termEquals("candidate", 0)) : "Exactly one candidate row follows baseline";
			final int[] columns={1, 2, 3, 6, 7, 8, 9, 10, 11};
			for(int i=0; i<columns.length; i++){
				assert(lp.parseInt(columns[i])==expected[i]) : "Saved-call literal oracle differs at column "+columns[i]+" in "+root;
			}
			assert(input.reader.nextLine()==null) : "Only baseline and candidate belong in a two-GFF comparison";
		}
	}

	private static void compareSaved(final File root, final String candidate, final File out, final String... extra){
		final ArrayList<String> args=new ArrayList<String>();
		args.add("mode=compare");
		args.add("baseline="+new File(root, "saved.before.gff"));
		args.add("calls="+new File(root, candidate));
		args.add("truth="+new File(root, "saved.truth.tsv"));
		args.add("family=example");
		args.add("truthfamily=fixture");
		args.add("out="+out);
		for(final String arg:extra){args.add(arg);}
		ReduceNcrnaModels.main(args.toArray(new String[args.size()]));
	}

	private static String savedCall(final int start, final int stop, final char strand){
		return "chr\tfixture\tRNA\t"+start+"\t"+stop+"\t.\t"+strand+"\t.\tfamily:example;model:fixture\n";
	}

	private static void write(final File path, final String text){
		try(Output output=new Output(path.getPath())){output.writer.print(text);}
	}

	private static void precisionBoundary(){
		final ArrayList<Locus> truth=new ArrayList<Locus>();
		final ArrayList<Locus> baseline=new ArrayList<Locus>();
		final ArrayList<Locus> candidate=new ArrayList<Locus>();
		for(int i=0; i<20; i++){
			truth.add(new Locus("chr", 1000*i+10, 1000*i+109, 0, null));
			baseline.add(new Locus("chr", 1000*i+9, 1000*i+109, 0, "model"));
			if(i<19){candidate.add(new Locus("chr", 1000*i+10, 1000*i+109, 0, "model"));}
		}
		candidate.add(new Locus("other", 10, 109, 0, "model"));
		final Result off=NcrnaBoundaryComparison.grade(truth, baseline, null);
		final Change change=new Change(off, NcrnaBoundaryComparison.grade(truth, candidate, off));
		assert(change.after.matched==19 && change.after.unmatched==1 && change.eligible(1, .95)
			&& !change.eligible(1, .96)) : "Nineteen correct calls among twenty is eligible exactly at precision0.95";
		final ArrayList<Locus> all=new ArrayList<Locus>(truth);
		final Result clean=NcrnaBoundaryComparison.grade(truth, all, off);
		final ArrayList<Locus> withExtra=new ArrayList<Locus>(all);
		withExtra.add(new Locus("other", 10, 109, 0, "model"));
		final Result noisy=NcrnaBoundaryComparison.grade(truth, withExtra, off);
		assert(NcrnaBoundaryComparison.choose(off, clean, noisy, 0, .95, Method.REGRESSION)==Method.CLASSIFIER)
			: "Fewer unmatched predictions precedes tie preference when both candidates pass precision";
	}

	private static Result grade(ArrayList<Locus> truth, int[][] calls, Result off){
		return NcrnaBoundaryComparison.grade(truth, loci(calls), off);
	}
	private static ArrayList<Locus> loci(int[][] positions){
		final ArrayList<Locus> loci=new ArrayList<Locus>();
		for(int[] p:positions){loci.add(new Locus("chr", p[0], p[1], p[2], "model"));}
		return loci;
	}
}
