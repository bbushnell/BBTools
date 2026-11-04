package test.prok;

import java.io.File;
import java.util.ArrayList;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaBoundaryComparison;
import prok.NcrnaBoundaryComparison.Change;
import prok.NcrnaBoundaryComparison.Result;
import prok.NcrnaLocusComparison.Locus;
import prok.NcrnaNativeTrainer.Method;
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
		System.out.println("Boundary comparison passed: detection conservation, both strands, losses/gains, duplicates, precision boundary, OFF and tie selection");
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
