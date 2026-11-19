package test.cardinality;

import java.util.Arrays;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.RrnaResourceIO;
import shared.Shared;
import structures.ByteBuilder;

/** Validates the literal ten-run alternating loglog comparison and reports
 * medians and ranges. This is a bounded no-penalty gate, not a significance test.
 * @author Ganyu
 */
public final class HashSmallTimingSummary{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, HashSmallTimingSummary.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				require(parser.parse(arg, a, b), "Unknown timing summary argument: "+arg);
			}
			require(parser.in1!=null && parser.out1!=null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use in=complete_ten_run_ledger.tsv out=fresh_summary.tsv without sampling or overwrite");
			RrnaResourceIO.fresh(parser.out1, "");
			final double[][] cpu=new double[2][5];
			final double[][] elapsed=new double[2][5];
			final LineParser1 fields=new LineParser1('\t');
			try(Input in=new Input(parser.in1)){
				final byte[] header=in.reader.nextLine();
				require(header!=null && Arrays.equals(header, "arm\trepeat\tuser_seconds\tsystem_seconds\telapsed_seconds".getBytes(java.nio.charset.StandardCharsets.US_ASCII)),
					"Timing ledger must have the exact five-column schema");
				for(int row=0; row<10; row++){
					final byte[] line=in.reader.nextLine();
					require(line!=null, "Ten completed alternating runs are required; missing row="+row);
					fields.set(line);
					final int arm=row&1;
					final int repeat=row/2;
					require(fields.terms()==5 && fields.termEquals(arm==0 ? "old" : "new", 0)
						&& fields.parseInt(1)==repeat+1, "Runs must alternate old/new five times, without omitted observations");
					final double user=fields.parseDouble(2);
					final double system=fields.parseDouble(3);
					final double wall=fields.parseDouble(4);
					require(Double.isFinite(user) && Double.isFinite(system) && Double.isFinite(wall)
						&& user>=0 && system>=0 && user+system>=10 && wall>=10,
						"Each real-file one-worker measurement must exceed ten CPU and elapsed seconds");
					cpu[arm][repeat]=user+system;
					elapsed[arm][repeat]=wall;
				}
				require(in.reader.nextLine()==null, "Only the predeclared ten observations belong in this comparison");
			}
			for(int arm=0; arm<2; arm++){
				Arrays.sort(cpu[arm]);
				Arrays.sort(elapsed[arm]);
			}
			final boolean pass=cpu[1][2]<=cpu[0][4];
			try(Output out=new Output(parser.out1)){
				out.writer.println("arm\truns\tcpu_median_seconds\tcpu_min_seconds\tcpu_max_seconds\telapsed_median_seconds\telapsed_min_seconds\telapsed_max_seconds");
				final ByteBuilder row=new ByteBuilder();
				for(int arm=0; arm<2; arm++){
					row.clear().append(arm==0 ? "old" : "new").tab().append(5).tab()
						.append(cpu[arm][2], 3).tab().append(cpu[arm][0], 3).tab().append(cpu[arm][4], 3).tab()
						.append(elapsed[arm][2], 3).tab().append(elapsed[arm][0], 3).tab().append(elapsed[arm][4], 3).nl();
					out.writer.print(row);
				}
			}
			pp.outstream.println("TIMING_RULE\tnew CPU median <= old CPU maximum; all observations retained");
			pp.outstream.println(pass ? "HASH_SMALL_TIMING_WITHIN_BASELINE_RANGE" : "HASH_SMALL_TIMING_PENALTY_NOT_ACCEPTED");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}
}
