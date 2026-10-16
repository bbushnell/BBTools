package assemble;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Selects exact attainable curve rows under fixed false-rescue budgets.
 * Never retunes a model or accesses TEST. Ties retain the higher cutoff.
 * @author Fischl
 */
public final class LiveBranchCurveSummary {

	/** Summarizes one model's exact validation curve through standard argument parsing. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, LiveBranchCurveSummary.class, false);
		try{
			final Parser parser=new Parser();
			long expected=-1;
			for(String arg : pp.args){
				final String[] split=arg.split("=", 2);
				String name=split[0].toLowerCase(Locale.ROOT);
				while(name.startsWith("-")){name=name.substring(1);}
				final String value=split.length>1 ? split[1] : null;
				if(name.equals("expected")){expected=Long.parseLong(value);}
				else if(!parser.parse(arg, name, value)){throw new IllegalArgumentException("Unknown curve-summary option: "+arg);}
			}
			if(parser.in1==null || parser.out1==null || expected<1 || parser.in2!=null || parser.out2!=null || new File(parser.out1).exists()){
				throw new IllegalArgumentException("Required in=grey.tsv out=NEW.tsv expected=validation_rows.");
			}
			summarize(parser.in1, parser.out1, expected);
			pp.outstream.println("LIVE_BRANCH_CURVE_SUMMARY_PASS rows="+expected+" budgets="+BUDGETS.length);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Validates count conservation and keeps the first maximal-TP row under each FP budget. */
	static void summarize(final String input, final String output, final long expected){
		assert(expected>0) : "False-rescue curves must describe a nonempty validation population.";
		final ByteBuilder[] best=new ByteBuilder[BUDGETS.length];
		final long[] bestTP=new long[BUDGETS.length];
		Arrays.fill(bestTP, -1);
		for(int i=0; i<best.length; i++){best[i]=new ByteBuilder();}
		final ByteFile file=ByteFile.makeByteFile(input, true);
		final LineParser1 p=new LineParser1('\t');
		long positive=-1, negative=-1, previousTP=0, previousFP=0, count=0;
		float previousCutoff=Float.POSITIVE_INFINITY;
		try{
			if(!Arrays.equals(file.nextLine(), HEADER.getBytes(StandardCharsets.UTF_8))){throw new IllegalArgumentException("Wrong exact-curve header.");}
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				p.set(line);
				if(p.terms()!=12){throw new IllegalArgumentException("Wrong curve width.");}
				final int bits=p.parseInt(0);
				final float cutoff=Float.intBitsToFloat(bits);
				final long gp=p.parseLong(2), gn=p.parseLong(3), tp=p.parseLong(4), fp=p.parseLong(5), tn=p.parseLong(6), fn=p.parseLong(7);
				// The emitter uses Float.toString; small thresholds legitimately contain exponents.
				final float printedCutoff=Float.parseFloat(p.parseString(1));
				if(!Float.isFinite(cutoff) || cutoff<0 || cutoff>Math.nextUp(1.0f) || printedCutoff!=cutoff ||
						!(cutoff<previousCutoff) || gp<0 || gn<0 || tp<0 || fp<0 || tn<0 || fn<0){
					throw new IllegalArgumentException("Invalid threshold or negative curve counts.");
				}
				if(count==0){
					positive=fn; negative=tn;
					if(tp!=0 || fp!=0 || gp!=0 || gn!=0 || positive+negative!=expected){
						throw new IllegalArgumentException("Missing reject-all state or wrong population size.");
					}
				}else if(gp+gn==0 || tp-previousTP!=gp || fp-previousFP!=gn){
					throw new IllegalArgumentException("Tied-group increments disagree with cumulative counts.");
				}
				if(tp+fn!=positive || fp+tn!=negative){throw new IllegalArgumentException("Curve population changed between thresholds.");}
				checkRate(p, 8, tp, tp+fp); checkRate(p, 9, tp, positive);
				checkRate(p, 10, fp, negative); checkRate(p, 11, fn, positive);
				for(int i=0; i<BUDGETS.length; i++){
					if(fp<=BUDGETS[i] && tp>bestTP[i]){bestTP[i]=tp; best[i].clear().append(line);}
				}
				previousCutoff=cutoff; previousTP=tp; previousFP=fp; count++;
			}
			if(count<2 || previousTP!=positive || previousFP!=negative){throw new IllegalArgumentException("Incomplete exact curve.");}
		}finally{if(file.close()){throw new IllegalStateException("Curve read failed.");}}
		final ByteStreamWriter writer=new ByteStreamWriter(output, false, false, false);
		writer.start();
		try{
			writer.print("fp_budget\t"+HEADER+"\n");
			final ByteBuilder out=new ByteBuilder();
			for(int i=0; i<BUDGETS.length; i++){
				assert(bestTP[i]>=0 && best[i].length()>0) : "Reject-all is feasible for every nonnegative false-rescue budget.";
				writer.print(out.clear().append(BUDGETS[i]).tab().append(best[i]).nl());
			}
		}finally{if(writer.poisonAndWait()){throw new IllegalStateException("Curve summary write failed.");}}
	}

	/** Verifies existing rate columns without replacing undefined rates by zero. */
	private static void checkRate(final LineParser1 p, final int field, final long numerator, final long denominator){
		assert(numerator>=0 && numerator<=denominator) : "A class-conditional numerator cannot exceed its population.";
		if(denominator==0){
			if(!p.termEquals("NA", field)){throw new IllegalArgumentException("Undefined rate must be NA.");}
		}else{
			final double rate=p.parseDouble(field);
			if(!Double.isFinite(rate) || Math.abs(rate-numerator/(double)denominator)>1e-9){
				throw new IllegalArgumentException("Curve rate field "+field+" is "+rate+", expected "+numerator+"/"+denominator+" from its confusion counts.");
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	static final String HEADER="threshold_bits\tthreshold\tgroup_positive\tgroup_negative\ttp\tfp\ttn\tfn\tprecision\trecall\tfpr\tfnr";
	static final long[] BUDGETS={0, 1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000, Long.MAX_VALUE};
}
