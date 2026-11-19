package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;

import cardinality.DynamicDemiLog;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.CompareNcrnaHashBoundary;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Literal query-versus-profile attribution and end-to-end saved-cache conservation.
 * @author Ganyu
 */
public final class CompareNcrnaHashBoundaryTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, CompareNcrnaHashBoundaryTest.class, false);
		final int exponent=DynamicDemiLog.exponentBits();
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				require(parser.parse(arg, a, b), "Unknown comparison fixture argument: "+arg);
			}
			require(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use out=fresh_fixture_directory without input, sampling or overwrite");
			final File root=new File(parser.out1).getAbsoluteFile();
			require(root.mkdir(), "Do not overwrite previous comparison evidence");
			DynamicDemiLog.setExponent(4);
			try(Output out=new Output(root+"/unicode_profile.tsv")){
				out.writer.print(CompareNcrnaHashBoundary.profileRow(new ByteBuilder(), 0, 0, bytes("model_Ω"), 1, 0, 0));
			}
			try(Input in=new Input(root+"/unicode_profile.tsv")){
				require(Arrays.equals(in.reader.nextLine(), bytes("0\t10\t2048\t0\t1\t0\tmodel_Ω\t1\t0\t0")) && in.reader.nextLine()==null,
					"The actual profile formatter must preserve the complete non-ASCII model label as UTF-8 bytes");
			}
			final int[] query=CompareNcrnaHashBoundary.literalRanks(new byte[][]{bytes("AAAA"),bytes("GAAA")}, bytes("AAANAAA"), 4, 32);
			require(Arrays.equals(query, new int[]{1,2,2}), "Removing the spurious query word must change only query-fixed rank first");
			final int[] members=CompareNcrnaHashBoundary.literalRanks(new byte[][]{bytes("AAANAAA"),bytes("GAAA")}, bytes("GAAA"), 4, 32);
			require(Arrays.equals(members, new int[]{2,2,1}), "N-free query rank can improve only when its N-containing member profile is corrected");
			HybridSavedNcrnaWindowsTest.main(new String[]{"out="+root+"/cache", "t=1"});
			final String fixture=root+"/cache/good";
			CompareNcrnaHashBoundary.run(config(fixture, root+"/one", 1));
			CompareNcrnaHashBoundary.run(config(fixture, root+"/two", 2));
			DynamicDemiLog.setExponent(6);
			CompareNcrnaHashBoundary.main(new String[]{"panels="+fixture+"/panels.tsv", "selection="+fixture+"/selection.tsv",
				"baseline="+fixture+"/one.loci.tsv.gz", "members="+fixture+"/members.fa", "ledger="+fixture+"/ledger.tsv",
				"assignments="+fixture+"/assignments.tsv", "census="+fixture+"/models.tsv", "summary="+fixture+"/members.summary.tsv",
				"ref="+fixture+"/ref.fa", "out="+root+"/cli", "expectedmodels=2", "expectedgroups=3", "t=1"});
			require(DynamicDemiLog.exponentBits()==6, "The standalone entry point must restore the caller's exponent after using exponent4");
			DynamicDemiLog.setExponent(4);
			for(final String suffix:new String[]{".profiles.tsv", ".changes.tsv.gz", ".loci.tsv.gz", ".summary.tsv"}){
				equal(root+"/one"+suffix, root+"/two"+suffix);
				equal(root+"/one"+suffix, root+"/cli"+suffix);
			}
			try(Input in=new Input(root+"/one.changes.tsv.gz")){
				require(in.reader.nextLine()!=null && in.reader.nextLine()==null, "N-free cache must produce no invented changed tuple");
			}
			final LineParser1 p=new LineParser1('\t');
			try(Input in=new Input(root+"/one.summary.tsv")){
				require(in.reader.nextLine()!=null, "Comparison summary requires a header");
				int rows=0;
				for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
					p.set(row);
					require(p.terms()==12 && p.parseLong(6)==0 && p.parseLong(7)==0 && p.parseLong(8)==11 && p.parseLong(9)==0
						&& p.parseLong(10)==0 && p.parseLong(11)==0,
						"Literal no-N cache has eleven windows, no changed tuple/locus, and keeps every reporting stratum");
					require(p.parseLong(3)==p.parseLong(4) && p.parseLong(4)==p.parseLong(5), "All three rank sums must match on the N-free cache");
					rows++;
				}
				require(rows==57, "Three strata times nineteen methods must remain explicit");
			}
			pp.outstream.println("NCRNA_HASH_BOUNDARY_COMPARISON_TEST_PASS");
		}finally{
			DynamicDemiLog.setExponent(exponent);
			Shared.closeStream(pp.outstream);
		}
	}

	private static CompareNcrnaHashBoundary.Config config(final String fixture, final String output, final int workers){
		final CompareNcrnaHashBoundary.Config c=new CompareNcrnaHashBoundary.Config();
		c.panels=fixture+"/panels.tsv";
		c.selection=fixture+"/selection.tsv";
		c.baseline=fixture+"/one.loci.tsv.gz";
		c.members=fixture+"/members.fa";
		c.ledger=fixture+"/ledger.tsv";
		c.assignments=fixture+"/assignments.tsv";
		c.census=fixture+"/models.tsv";
		c.summary=fixture+"/members.summary.tsv";
		c.references=fixture+"/ref.fa";
		c.models=2;
		c.groups=3;
		c.threads=workers;
		c.output=output;
		return c;
	}
	private static byte[] bytes(final String s){return s.getBytes(StandardCharsets.UTF_8);}
	private static void equal(final String a, final String b){
		try(Input left=new Input(a); Input right=new Input(b)){
			for(byte[] row=left.reader.nextLine(); row!=null; row=left.reader.nextLine()){
				require(Arrays.equals(row, right.reader.nextLine()), "Worker count cannot change a literal comparison result: "+a);
			}
			require(right.reader.nextLine()==null, "Parallel comparison cannot append rows");
		}
	}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}
}
