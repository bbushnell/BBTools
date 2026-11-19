package test.prok;

import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlMoments;
import shared.Shared;

/** Literal weighted population with unequal event/model counts, plus reuse.
 * @author Brian, Raiden
 */
public final class NcrnaDdlMomentsTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaDdlMomentsTest.class, false);
		try{
			final Parser parser=new Parser();
			boolean run=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(a.equals("run")){run=Parse.parseBoolean(b);}
				else{check(parser.parse(arg, a, b), "Unknown moment fixture parameter: "+arg);}
			}
			boolean assertions=false;
			assert(assertions=true);
			check(run && assertions && parser.in1==null && parser.in2==null && parser.out1==null && parser.out2==null,
				"Use run=t and -ea; weighted-moment fixtures take no external inputs or outputs");
			verify();
			pp.outstream.println("NCRNA_DDL_MOMENTS_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void verify(){
		final NcrnaDdlMoments population=new NcrnaDdlMoments();
		check(Double.isNaN(population.mean()) && Double.isNaN(population.standardDeviation()),
			"An empty fit population is unavailable, not a mean or variance of zero");
		// Locus A: events [0,4] and [2]. Locus B: one event [10,10].
		// Model weights: 1/8,1/8,1/4,1/4,1/4; mean6,second53,variance17.
		population.addLocus(4, 12, 2, 3);
		population.addLocus(10, 100, 1, 2);
		check(population.mean()==6 && population.standardDeviation()==Math.sqrt(17),
			"Equal locus/event weights must retain within-event variation: variance17, not variance16 of the event averages");
		check(population.loci()==2 && population.events()==3 && population.models()==5,
			"Raw model and event observations remain separate from two locus weights");
		population.clear();
		check(population.loci()==0 && population.events()==0 && population.models()==0 && Double.isNaN(population.mean()),
			"Reusing a cell cannot retain the prior fit population");
		population.addLocus(6, 18, 2, 10);
		check(population.mean()==3 && population.standardDeviation()==0,
			"A constant population has zero variance; no positive floor may make an empirical weight available");
		invalid(population, 0, 0, 0, 0);
		invalid(population, 2, 1, 1, 1);
		invalid(population, Double.NaN, 1, 1, 1);
		invalid(population, 1, Double.POSITIVE_INFINITY, 1, 1);
		check(population.loci()==1 && population.mean()==3,
			"A rejected empty/invalid class cannot change the previously accepted moments");
		population.clear();
		population.addLocus(1, 1-1e-14, 1, 1);
		check(population.standardDeviation()==0, "Tiny negative cancellation is reported as zero variance");
		population.clear();
		population.addLocus(1, 1+1e-14, 1, 1);
		check(population.standardDeviation()>0, "Tiny positive variance is retained without a positive floor");
	}

	private static void invalid(final NcrnaDdlMoments population, final double first,
			final double second, final int events, final long models){
		try{population.addLocus(first, second, events, models);}
		catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Invalid model/event moments were accepted");
	}
	private static void check(final boolean valid, final String why){
		if(!valid){throw new AssertionError(why);}
	}
}
