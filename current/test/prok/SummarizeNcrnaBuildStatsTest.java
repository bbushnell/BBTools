package test.prok;

import java.io.File;

import parse.Parse;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.SummarizeNcrnaBuildStats;
import shared.Shared;

/** Literal rounded-bin boundaries and count-conservation cases; no alignment.
 * @author Brian Bushnell, Raiden
 */
public final class SummarizeNcrnaBuildStatsTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SummarizeNcrnaBuildStatsTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(java.util.Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown fixture parameter: "+arg);}
			}
			if(parser.out1==null){throw new IllegalArgumentException("Supply out=fresh absolute fixture directory");}
			final File root=new File(parser.out1);
			if(!root.isAbsolute() || !root.mkdir()){throw new IllegalArgumentException("Fixture directory must be fresh and absolute");}
			run(root);
			pp.outstream.println("NCRNA_BUILD_IDENTITY_HISTOGRAM_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final File root){
		final String input=new File(root, "eight.tsv").getPath();
		try(Output out=new Output(input)){
			out.writer.println("kind\tvalue\tcount");
			out.writer.println("input_members\tall\t8");
			out.writer.println("identity_nearest_thousandth\t638\t2");
			out.writer.println("identity_nearest_thousandth\t640\t1");
			out.writer.println("identity_nearest_thousandth\t700\t2");
			out.writer.println("identity_nearest_thousandth\t751\t3");
		}
		final String output=new File(root, "summary.tsv").getPath();
		SummarizeNcrnaBuildStats.main(new String[]{"in="+input, "labels=literal", "out="+output, "expectedmembers=8"});
		final int[] lower={2,3,5}, upper={3,5,5};
		final float[] cuts={.640f,.700f,.750f};
		final double[] lowPercent={25,37.5,62.5}, highPercent={37.5,62.5,62.5};
		final LineParser1 fields=new LineParser1('\t');
		int rows=0;
		try(Input in=new Input(output)){
			in.reader.nextLine();
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				fields.set(line);
				check(rows<3 && fields.terms()==10 && fields.termEquals("literal", 0) && fields.parseInt(1)==8,
					"Every literal cutoff needs one complete ten-field row on the eight-member population");
				check(fields.parseFloat(2)==.638f && fields.parseFloat(3)==.700f && fields.parseFloat(4)==.751f,
					"Nearest-rank quartiles must use members, not nonempty histogram bins");
				check(fields.parseFloat(5)==cuts[rows] && fields.parseInt(6)==lower[rows] && fields.parseInt(7)==upper[rows]
					&& fields.parseDouble(8)==lowPercent[rows] && fields.parseDouble(9)==highPercent[rows],
					"Threshold-bin uncertainty must retain literal lower/upper counts and percent units, independent of trailing zero formatting");
				rows++;
			}
		}
		check(rows==3, "Both uncertain threshold bins and the exact empty-bin case must be reported");
		boolean rejected=false;
		try{SummarizeNcrnaBuildStats.readHistogram(input, 9);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "A denominator mismatch must fail, not silently change reported fractions");
		final String duplicate=new File(root, "duplicate.tsv").getPath();
		try(Output out=new Output(duplicate)){
			out.writer.println("kind\tvalue\tcount");
			out.writer.println("identity_nearest_thousandth\t640\t4");
			out.writer.println("identity_nearest_thousandth\t640\t4");
		}
		rejected=false;
		try{SummarizeNcrnaBuildStats.readHistogram(duplicate, 8);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "A repeated bin must not masquerade as a conserved complete histogram");
	}

	private static void check(final boolean valid, final String why){
		if(!valid){throw new AssertionError(why);}
	}
}
